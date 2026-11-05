#!/usr/bin/env python3
"""Generate golden/stats_golden.json for the fits/stats package.

Every expected value comes from an independent source, never from the Go
code under test:

  numpy 2.x      mean, std(ddof=1), min/max, numpy.histogram counts/edges,
                 the CDF of those counts
  fractions      exact (rational) mean and sample stdev, correctly rounded
                 to float64 -- the truth for MeanStdevTwoPass -- and the
                 exact histogram-percentile definition
  plain Python   the one-pass formula MeanStdev documents, evaluated in
  floats         IEEE double in the same order (pins its documented
                 cancellation behaviour bit-for-bit)
  scipy          scipy.stats.rv_histogram ppf: the histogram quantile,
                 independent of Histogram.Percentile's own definition
  astropy        astropy.stats.SigmaClip (sigma_clip) for SigmaClip; the median centre
                 is Siril's histogram median, run from the verbatim Siril
                 source via ./sref/sref median-f32|median-f64 (build it
                 first: make -C sref)

Usage (from stats/testdata, with the project venv, see PROVENANCE.md):
    make -C sref
    python gen_golden.py [output_dir]      # default: golden/

Re-running reproduces the committed file byte-for-byte.

JSON number encoding: finite floats are JSON numbers (Python repr, which
round-trips float64 exactly); NaN and +-Inf are the strings "nan", "inf",
"-inf". Inputs: float32 as "input_f32", a space-separated list of hex
IEEE bit patterns; integer types as "input_int", space-separated decimal
integers (exact for int64/uint64); float64 as "input", JSON numbers, or as
"input_runs", [value, count] pairs, when runs of equal values make that
under a quarter of the length, or as "input_cycle", [pattern, repeats],
when a short pattern repeats. Histogram edges and CDF are {"sha256",
"sample"}: the sha256 of the float64 values as big-endian IEEE bits, and
a few explicit values (see digest_f64).
Large synthetic inputs are not committed: "formula" names a deterministic
xorshift32 recurrence (formula_f32 / formula_f64 below) that the Go test
(formulaValues / formulaValues64 in siril_test.go; uniform24, astro_u16
and uniform53 are shared with the Siril harness) regenerates bit-for-bit,
and "input_sha256" (sha256 of the big-endian IEEE bits) proves it did.
A case with "gotype" is run in Go as that element type (namedf32: a named
type whose underlying type is float32, which takes the float64 path).
"""

import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import warnings
from decimal import Decimal, getcontext
from fractions import Fraction

import numpy as np
import scipy.stats
from astropy.stats import SigmaClip

HERE = os.path.dirname(os.path.abspath(__file__))
SREF = os.environ.get("SREF_BIN") or os.path.join(HERE, "sref", "sref")
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "golden")

getcontext().prec = 80


def jf(x):
    """JSON-safe float: finite as number, else "nan"/"inf"/"-inf"."""
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    return x


def xorshift32(x):
    x ^= (x << 13) & 0xFFFFFFFF
    x ^= x >> 17
    x ^= (x << 5) & 0xFFFFFFFF
    return x


def formula_f32(kind, seed, n):
    """formula() of sref/sref.cc and formulaData of siril_test.go. Only
    exactly rounded IEEE operations: an integer times 2^-24, or one
    float32 division of two integers."""
    x = seed
    ints = []
    for _ in range(n):
        x = xorshift32(x)
        b = x
        if kind == "uniform24":
            ints.append(b >> 8)
        elif kind == "astro_u16":
            adu = 2900 + (b & 0xFF) + ((b >> 8) & 0xFF) + ((b >> 16) & 0xFF) + (b >> 24)
            if b % 509 == 0:
                adu += (b >> 7) % 60000
            ints.append(min(adu, 65535))
        elif kind in ("sky_hot", "adu_f32"):
            k = (b & 0xFF) + ((b >> 8) & 0xFF) + ((b >> 16) & 0xFF) + (b >> 24) - 510
            hot = kind == "sky_hot" and b % 100 == 7
            ints.append((k, hot, (b >> 16) & 0xFFF))
        else:
            raise SystemExit(f"unknown formula {kind}")
    if kind == "sky_hot":
        # Sky 0.05 + k*2^-17 (k a sum of four bytes minus 510: sigma ~0.0011),
        # 1% hot pixels 0.9 + j*2^-16. One float32 rounding per value.
        k = np.array([v[0] for v in ints], dtype=np.float32)
        j = np.array([v[2] for v in ints], dtype=np.float32)
        hot = np.array([v[1] for v in ints])
        sky = np.float32(0.05) + k * np.float32(2.0 ** -17)
        return np.where(hot, np.float32(0.9) + j * np.float32(2.0 ** -16), sky).astype(np.float32)
    if kind == "adu_f32":
        # 30000 + k/256: exact in float32 (sigma ~0.58 ADU).
        k = np.array([v[0] for v in ints], dtype=np.float32)
        return (np.float32(30000) + k * np.float32(2.0 ** -8)).astype(np.float32)
    a = np.array(ints, dtype=np.float32)
    if kind == "uniform24":
        return a * np.float32(5.9604644775390625e-8)
    return a / np.float32(65535)


def formula_f64(kind, seed, n):
    """float64 formula inputs (formulaValues64 in siril_test.go):
    "uniform53_1e9" is 1e9 + u, u = 53 bits from two draws times 2^-53
    (exact), the sum rounded once."""
    x = seed
    out = []
    for _ in range(n):
        x = xorshift32(x)
        a = x
        if kind == "wide_i32":
            # int32(b) * 2^467: up to ~1.4e150, exact.
            out.append(float(a - (1 << 32) if a >= 1 << 31 else a) * 2.0 ** 467)
            continue
        x = xorshift32(x)
        m = ((a << 32) | x) >> 11
        if kind == "uniform53_1e9":
            out.append(1e9 + m * 2.0 ** -53)
        elif kind == "uniform53":
            out.append(m * 2.0 ** -53)
        else:
            raise SystemExit(f"unknown formula64 {kind}")
    return np.array(out, dtype=np.float64)


# Cases whose input is a formula, not committed (name -> (kind, seed)).
FORMULA_CASES = {
    "float32_large": ("astro_u16", 21), "float64_offset_noise": ("uniform53_1e9", 23),
    "float64_random": ("uniform53", 24), "float32_adu_offset": ("adu_f32", 25),
    "float32_outliers": ("sky_hot", 26), "float64_wide": ("wide_i32", 27),
}
# Cases run in Go as a named element type.
GOTYPE = {"namedf32_random": "namedf32"}


# ---------------------------------------------------------------------
# golden/mean_rounding.json: MeanStdevTwoPass's mean is correctly rounded.
# The inputs are not committed: mean_rounding_sets() is an xorshift32
# recurrence that TestMeanStdevTwoPass_MeanCorrectlyRounded regenerates
# (meanRoundingSets in stats_test.go). Each set holds 2..64 integers
# +-(2^e + k), e in 40..62, |k| < 2^23: beyond 2^53 they are rounded to
# float64 first (round to nearest even, as Go's float64(int64)), and the
# truth is the exact mean of those float64 values, correctly rounded
# (Python int / int true division). A quotient rounded once more than
# necessary, or a remainder taken from a rounded product, is visibly
# wrong on a large share of these sets.
# ---------------------------------------------------------------------

MEAN_ROUNDING_SEED = 0x2545F491
MEAN_ROUNDING_SETS = 20000


def mean_rounding_sets():
    x = MEAN_ROUNDING_SEED
    for _ in range(MEAN_ROUNDING_SETS):
        x = xorshift32(x)
        n = 2 + x % 63
        x = xorshift32(x)
        e = 40 + x % 23
        x = xorshift32(x)
        neg = x & 1
        vals = []
        for _ in range(n):
            x = xorshift32(x)
            v = (1 << e) + (x >> 8) - (1 << 23)
            vals.append(-v if neg else v)
        yield vals


def mean_rounding_golden():
    digest = hashlib.sha256()
    samples = {}
    inexact = 0
    for i, vals in enumerate(mean_rounding_sets()):
        fl = [float(v) for v in vals]  # correctly rounded int -> float64
        s = sum(int(f) for f in fl)  # every value is an integer
        n = len(fl)
        mean = s / n  # int / int: correctly rounded
        # How many sets have an exact mean that is not a float64 (so a
        # rounding decision is made at all).
        if mean * n != s:
            inexact += 1
        bits = np.array([mean], dtype=">f8").tobytes()
        digest.update(bits)
        if i < 16 or i % 1000 == 0:
            samples[str(i)] = bits.hex()
    return {
        "about": "MeanStdevTwoPass mean of xorshift32 integer sets, correctly rounded; see gen_golden.py mean_rounding_sets",
        "seed": MEAN_ROUNDING_SEED,
        "sets": MEAN_ROUNDING_SETS,
        "inexact": inexact,
        "sha256_be_f64": digest.hexdigest(),
        "samples": samples,
    }


# ---------------------------------------------------------------------
# golden/histogram_ranges.json: numpy.histogram on ranges at the limits of
# float64 (the |min| >= 2^52 padding, overflowing max-min, subnormal and
# 2^60 ranges just wide or just too narrow for nbins edges). numpy either
# bins (counts and edges, as float64 bit patterns) or raises; its exception
# is recorded and BuildHistogram must return ErrHistogramRange there.
# ---------------------------------------------------------------------

HIST_RANGES = [
    ([2.0 ** 52] * 2, 1), ([2.0 ** 52] * 2, 2), ([-(2.0 ** 52)] * 2, 1),
    ([2.0 ** 53] * 2, 1), ([2.0 ** 53] * 2, 2),
    ([-1e308, 1e308], 1), ([-1e308, 1e308], 2), ([1.7e308, -1.7e308, 0.0], 1),
    ([1e308] * 3, 1),
    ([0.0, 5e-324], 1), ([0.0, 5e-324], 2), ([0.0, 5e-324], 3),
    ([0.0, 1e-323], 2), ([0.0, 1e-322], 7),
    ([1.0], 1), ([2.0 ** 60, 2.0 ** 60 + 768], 3), ([2.0 ** 60, 2.0 ** 60 + 768], 4),
]


def f64hex(x):
    return np.array([x], dtype=">f8").tobytes().hex()


def histogram_ranges_golden():
    out = []
    for values, nbins in HIST_RANGES:
        r = {"values": [f64hex(v) for v in values], "nbins": nbins}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                counts, edges = np.histogram(np.array(values, dtype=np.float64), bins=nbins)
            except (ValueError, IndexError) as e:
                r["numpy_error"] = f"{type(e).__name__}: {e}"
            else:
                r["counts"] = [int(c) for c in counts]
                r["edges"] = [f64hex(e) for e in edges]
        out.append(r)
    return out


# ---------------------------------------------------------------------
# golden/histogram_edges.json: values on and next to numpy's bin edges,
# where its computed index ((x-first)/(last-first))*nbins is off by one
# and numpy corrects it (x < edges[i]: decrement; x >= edges[i+1] outside
# the last bin: increment). Each "edges" case's input is not committed:
# it is [edges[0], edges[nbins]] then, for every edge e, nextafter(e,
# -inf) (not for edges[0]), e (interior edges only) and nextafter(e, +inf)
# (not for edges[nbins]), in that order -- built by the Go test from the
# committed numpy edges with math.Nextafter. "values" cases commit their
# input. "decrements" / "increments" count how often numpy's correction
# fires, from numpy's own index computation replayed here; the Go test
# requires both to be exercised.
# ---------------------------------------------------------------------

HIST_EDGE_SEED = 5
HIST_EDGE_RANDOM = 24
HIST_EDGE_FIXED = [(0.0, 1.0, 10), (9.943696218777372, 103.20010657904987, 39)]
HIST_EDGE_VALUES = [([0.0, 0.3, 0.6, 0.7, 1.0], 10)]


def numpy_corrections(x, nbins):
    """numpy.histogram's uniform-bin index computation (numpy 2.x
    _histograms_impl.histogram), replayed to count its corrections; the
    resulting counts must equal np.histogram's."""
    counts, edges = np.histogram(x, bins=nbins)
    first, last = edges[0], edges[-1]
    idx = (((x - first) / (last - first)) * nbins).astype(np.intp)
    idx[idx == nbins] -= 1
    dec = x < edges[idx]
    idx[dec] -= 1
    inc = (x >= edges[idx + 1]) & (idx != nbins - 1)
    idx[inc] += 1
    if not (np.bincount(idx, minlength=nbins) == counts).all():
        raise SystemExit("numpy index replay disagrees with np.histogram")
    return counts, edges, int(dec.sum()), int(inc.sum())


def edge_input(edges):
    """The "edges" case input (see above)."""
    n = len(edges) - 1
    v = [edges[0], edges[n]]
    for i, e in enumerate(edges):
        if i > 0:
            v.append(np.nextafter(e, -np.inf))
        if 0 < i < n:
            v.append(e)
        if i < n:
            v.append(np.nextafter(e, np.inf))
    return np.array(v, dtype=np.float64)


def histogram_edges_golden():
    rng = np.random.RandomState(HIST_EDGE_SEED)
    configs = list(HIST_EDGE_FIXED)
    for _ in range(HIST_EDGE_RANDOM):
        lo = rng.uniform(-100, 100)
        hi = lo + 10 ** rng.uniform(-3, 4)
        configs.append((float(lo), float(hi), int(rng.randint(3, 80))))
    out = []
    for lo, hi, nbins in configs:
        _, edges = np.histogram(np.array([lo, hi]), bins=nbins)
        x = edge_input(edges)
        counts, e2, dec, inc = numpy_corrections(x, nbins)
        if not np.array_equal(e2.view(np.uint64), edges.view(np.uint64)):
            raise SystemExit("edge input changed numpy's edges")
        out.append({"edges": [f64hex(e) for e in edges], "counts": [int(c) for c in counts],
                    "decrements": dec, "increments": inc})
    for values, nbins in HIST_EDGE_VALUES:
        x = np.array(values, dtype=np.float64)
        counts, edges, dec, inc = numpy_corrections(x, nbins)
        out.append({"values": [f64hex(v) for v in values], "nbins": nbins,
                    "edges": [f64hex(e) for e in edges], "counts": [int(c) for c in counts],
                    "decrements": dec, "increments": inc})
    return out


def gen_test_arrays():
    """Test arrays. The RandomState(42) draws keep their original order so
    the arrays committed before this generator gained new cases are
    unchanged."""
    rng = np.random.RandomState(42)
    cases = {}

    cases["float32_random"] = rng.uniform(0, 1, 500).astype(np.float32)
    cases["float32_gradient"] = np.linspace(0.01, 0.3, 200).astype(np.float32)
    cases["float32_with_nan"] = np.array(
        [0.1, 0.2, float("nan"), 0.4, 0.5, float("nan"), 0.3, 0.15],
        dtype=np.float32,
    )
    cases["float32_with_zeros"] = np.array(
        [0.0, 0.1, 0.0, 0.2, 0.0, 0.3, 0.15, 0.25, 0.0, 0.05],
        dtype=np.float32,
    )
    cases["float32_all_same"] = np.full(100, 0.42, dtype=np.float32)
    cases["float32_single"] = np.array([3.14], dtype=np.float32)
    cases["float32_empty"] = np.array([], dtype=np.float32)
    cases["float32_all_nan"] = np.array([float("nan")] * 5, dtype=np.float32)
    # float32_large was rng.normal(100, 15, 100000); it is now a formula
    # input (not committed). The draw is kept so the arrays after it are
    # unchanged.
    rng.normal(100, 15, 100000)
    cases["float32_large"] = formula_f32(*FORMULA_CASES["float32_large"], 100000)
    cases["float32_negative"] = rng.uniform(-1, 1, 300).astype(np.float32)
    cases["int16_random"] = rng.randint(0, 30000, 500).astype(np.int16)
    cases["uint8_random"] = rng.randint(0, 255, 500).astype(np.uint8)
    cases["int32_negative"] = rng.randint(-1000, 1000, 200).astype(np.int32)

    rng = np.random.RandomState(7)
    # float64_random, float32_adu_offset and float32_outliers were
    # RandomState(7) draws (uniform, 30000 + normal(0, 0.5), a normal sky
    # with 50 hot pixels); they are now formula inputs (not committed). The
    # draws are kept so the arrays after them are unchanged.
    rng.uniform(0, 1, 1000)
    cases["float64_random"] = formula_f64(*FORMULA_CASES["float64_random"], 1000)
    # MeanStdev's documented cancellation case: true sample stdev
    # sqrt(250/999) = 0.50025...
    cases["float64_offset_1e9"] = np.array([1e9, 1e9 + 1] * 500, dtype=np.float64)
    rng.normal(0, 0.5, 2000)
    cases["float32_adu_offset"] = formula_f32(*FORMULA_CASES["float32_adu_offset"], 2000)
    sky = rng.normal(0.05, 0.003, 5000)
    rng.choice(5000, 50, replace=False)
    rng.uniform(0, 0.05, 50)
    cases["float32_outliers"] = formula_f32(*FORMULA_CASES["float32_outliers"], 5000)
    ints = np.round(rng.normal(1000, 30, 500))
    ints[rng.choice(500, 10, replace=False)] = 30000
    cases["int16_hot"] = ints.astype(np.int16)
    withinf = rng.uniform(0, 1, 50).astype(np.float32)
    cases["float32_with_inf"] = np.concatenate(
        [withinf, np.array([np.inf, -np.inf], dtype=np.float32), withinf[:5]]
    )
    cases["float64_overflow"] = np.array([1e308, 1e308, 1e308])
    cases["float64_with_posinf"] = np.array([1.0, 2.0, np.inf, 4.0])
    cases["float64_square_overflow"] = np.array([1e200, -1e200, 1e200, -1e200])
    # More MeanStdevTwoPass accuracy cases (exact rational truth).
    # float64_offset_noise was 1e9 + rng.uniform(0, 1, 3000); it is now a
    # formula input (not committed). The draw is kept so the arrays after
    # it are unchanged.
    rng.uniform(0, 1, 3000)
    cases["float64_offset_noise"] = formula_f64(*FORMULA_CASES["float64_offset_noise"], 3000)
    cases["float32_tiny"] = rng.uniform(0, 1e-30, 500).astype(np.float32)
    rng.normal(0, 1e150, 1000)  # was float64_wide; now a formula input
    cases["float64_wide"] = formula_f64(*FORMULA_CASES["float64_wide"], 1000)
    # A value larger than the running sum arriving after many small ones:
    # the branch of the Neumaier sum that keeps the running sum's lost low
    # bits (|x| > |s|). Without it the mean of the first is ~60 ULP off and
    # the second gives 0.25 instead of 0.5.
    cases["float64_big_after_small"] = np.array([0.1] * 1000 + [1e16, -1e16])
    cases["float64_neumaier_classic"] = np.array([1.0, 1e100, 1.0, -1e100])

    # Element types with their own float64-conversion failure modes, and
    # sigma-clip cases for asymmetric bounds and maxiters=None.
    rng = np.random.RandomState(8)
    cases["int8_random"] = rng.randint(-128, 128, 600).astype(np.int8)
    cases["uint32_big"] = rng.randint(0, 2**32, 600, dtype=np.uint64).astype(np.uint32)
    cases["uint64_big"] = rng.randint(0, 2**64 - 1, 600, dtype=np.uint64)
    cases["int64_big"] = rng.randint(-(2**62), 2**62, 600, dtype=np.int64)
    # Clustered far beyond 2^53: float64 spacing is 256 (int64) and 2048
    # (uint64) there, so converting to float64 visibly changes the stdev.
    cases["int64_offset_2p60"] = np.int64(2**60) + rng.randint(0, 1000, 500).astype(np.int64)
    cases["uint64_offset_2p63"] = np.uint64(2**63) + rng.randint(0, 10000, 500).astype(np.uint64)
    cases["namedf32_random"] = rng.uniform(0, 1, 500).astype(np.float32)
    # A geometric tail that sigma clipping peels off one value per
    # iteration: astropy needs more than 100 iterations with maxiters=None.
    cases["float64_geometric_tail"] = np.concatenate([rng.normal(0, 1, 20), 2.0 ** np.arange(1, 301)])
    # Sigma clipping that rejects every value (SIGMA_EXTRA): with k = 0.5
    # the first iteration drops 100, the second drops both clusters (each
    # is exactly one sigma from the mean), the third runs on an empty set.
    cases["float64_all_rejected"] = np.array([0.0] * 5 + [1.0] * 5 + [100.0])
    cases["float32_all_rejected"] = np.array([0.0] * 5 + [1.0] * 5 + [100.0], dtype=np.float32)
    # max - min overflows float64: numpy raises, BuildHistogram errors.
    cases["float64_range_overflow"] = np.array([-1e308, 1e308, 5.0])
    return cases


def exact_sqrt(fr):
    """Correctly rounded float64 sqrt of a non-negative Fraction (80
    significant digits before the final rounding)."""
    if fr == 0:
        return 0.0
    return float((Decimal(fr.numerator) / Decimal(fr.denominator)).sqrt())


def fma(a, b, c):
    """a*b + c with one rounding (IEEE fusedMultiplyAdd), via exact
    rationals; non-finite operands follow IEEE (as a*b + c does)."""
    if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(c)):
        return a * b + c
    exact = Fraction(a) * Fraction(b) + Fraction(c)
    try:
        return float(exact)
    except OverflowError:
        return math.inf if exact > 0 else -math.inf


def onepass_meanstdev(values, fuse_acc=False, fuse_sub=False):
    """The formula MeanStdev documents, in IEEE double, same order.

    MeanStdev is frozen and has no explicit float64() conversions, so the
    Go compiler may fuse its multiply-adds where the target has FMA
    instructions (Go spec, "Floating-point operators"). fuse_acc models
    `sum2 += x*x` fused (fma(x, x, sum2)); fuse_sub models
    `sum2 - float64(n)*mean*mean` fused (one rounding of
    sum2 - (n*mean)*mean). Which build fuses what is measured by
    TestMeanStdev_OnePassFormula (stats/fusion_*_test.go)."""
    s = 0.0
    s2 = 0.0
    n = 0
    for x in values:
        x = float(x)
        s += x
        s2 = fma(x, x, s2) if fuse_acc else s2 + x * x
        n += 1
    if n == 0:
        return 0.0, 0.0
    mean = s / float(n)
    if n < 2:
        return mean, 0.0
    nm = float(n) * mean
    num = fma(-nm, mean, s2) if fuse_sub else s2 - nm * mean
    var = num / float(n - 1)
    if var < 0:
        var = 0.0
    return mean, math.sqrt(var)


def sref_median(values, f32):
    """Siril's histogram median (float32) or its float64 compile."""
    if f32:
        arr = np.asarray(values, dtype=np.float64).astype(np.float32)
        if not np.array_equal(arr.astype(np.float64), np.asarray(values, dtype=np.float64)):
            raise SystemExit("sref_median: values are not float32")
        hexs = " ".join(f"{u:08x}" for u in arr.view(np.uint32))
        mode = "median-f32"
    else:
        arr = np.asarray(values, dtype=np.float64)
        hexs = " ".join(f"{u:016x}" for u in arr.view(np.uint64))
        mode = "median-f64"
    if arr.size == 0:
        return float("nan")
    out = subprocess.run([SREF, mode], input=hexs, capture_output=True, text=True, check=True).stdout.strip()
    if f32:
        return float(np.array([int(out, 16)], dtype=np.uint32).view(np.float32)[0])
    return float(np.array([int(out, 16)], dtype=np.uint64).view(np.float64)[0])


# (sigma_lower, sigma_upper, maxiters, center); maxiters 0 is astropy's
# maxiters=None (SigmaClip's maxIter <= 0: iterate until nothing changes).
SIGMA_VARIANTS = [
    (2.0, 2.0, 10, "mean"), (3.0, 3.0, 10, "mean"), (5.0, 5.0, 10, "mean"),
    (2.0, 2.0, 10, "median"), (3.0, 3.0, 10, "median"), (5.0, 5.0, 10, "median"),
    (3.0, 3.0, 0, "mean"), (3.0, 3.0, 0, "median"),
    (2.0, 4.0, 0, "mean"), (2.0, 4.0, 0, "median"),
    (4.0, 1.5, 0, "mean"), (4.0, 1.5, 0, "median"),
    (3.0, 3.0, 1, "mean"),
]

# Extra variants for the all-rejected cases: maxiters None, 10, 2 and 1
# with both centres, so the iteration count of the empty-set branch is
# pinned (astropy counts the iteration that runs on the empty set).
SIGMA_EXTRA = {
    name: [(0.5, 0.5, m, c) for c in ("mean", "median") for m in (0, 10, 2, 1)]
    for name in ("float64_all_rejected", "float32_all_rejected")
}


def exact_mean_popstd(values):
    """Correctly rounded mean and population (ddof=0) standard deviation of
    float64 values, in exact integer arithmetic: every value is V_i * 2^k
    with integer V_i and one common k. NaN for empty input."""
    vals = [float(v) for v in values]
    n = len(vals)
    if n == 0:
        return float("nan"), float("nan")
    parts = []
    for v in vals:
        num, den = v.as_integer_ratio()  # den is a power of two
        parts.append((num, den.bit_length() - 1))
    kmax = max(e for _, e in parts)
    V = [num << (kmax - e) for num, e in parts]  # value = V / 2^kmax
    S = sum(V)
    scale = Fraction(1, 1 << kmax)
    mean = float(Fraction(S, n) * scale)
    m2 = sum((n * x - S) ** 2 for x in V)  # = n^2 * sum (x - S/n)^2
    var = Fraction(m2, n ** 3) * scale * scale
    return mean, exact_sqrt(var)


def compute_sigma_clip(finite, f32, klow, khigh, maxiters, center):
    """astropy.stats.SigmaClip with the documented SigmaClip semantics;
    also returns astropy's iteration count (_niterations: the iterations
    run, including the last one that rejects nothing).

    The centre and sigma astropy is given are exact: the correctly rounded
    mean and population std (exact_mean_popstd), or Siril's histogram
    median. numpy's np.mean/np.std would put numpy's own rounding into the
    truth, and it is far from exact for data with a large offset (for
    uint64_offset_2p63 np.std is 3119.42, the exact value 2981.4)."""
    if center == "median":
        def cen(a, axis=None):
            return sref_median(a, f32)
    else:
        def cen(a, axis=None):
            return exact_mean_popstd(a)[0]

    def std(a, axis=None):
        return exact_mean_popstd(a)[1]

    sc = SigmaClip(sigma_lower=klow, sigma_upper=khigh, maxiters=maxiters or None,
                   cenfunc=cen, stdfunc=std)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        kept = sc(finite, masked=False, axis=None, copy=True)
    kept = np.asarray(kept, dtype=np.float64)
    kept = kept[np.isfinite(kept)]
    r = {"klow": klow, "khigh": khigh, "maxiters": maxiters, "center": center}
    if kept.size == 0:
        r.update({"mean": 0.0, "stdev": 0.0, "median": 0.0, "ngood": 0})
    else:
        km, ks = exact_mean_popstd(kept)
        r.update({
            "mean": jf(km),
            "stdev": jf(ks),
            "median": jf(sref_median(kept, f32)),
            "ngood": int(kept.size),
        })
    r["niter"] = int(sc._niterations)
    return r


def hist_percentile_exact(counts, edges, p):
    """Histogram.Percentile's documented definition, exact arithmetic
    (t = p*Total in float64, as documented)."""
    total = int(sum(counts))
    if total == 0:
        return 0.0
    p = min(max(p, 0.0), 1.0)
    t = Fraction(p * float(total))
    cum = 0
    for i, c in enumerate(counts):
        cum += int(c)
        if cum >= t and c > 0:
            prev = cum - int(c)
            e0, e1 = Fraction(edges[i]), Fraction(edges[i + 1])
            return float(e0 + (t - prev) / int(c) * (e1 - e0))
    return float(edges[-1])


def scipy_hist_ppf(counts, edges):
    """The quantile function of the histogram's piecewise-linear
    distribution from scipy.stats.rv_histogram (its ppf interpolates the
    CDF linearly inside each bin): an independent evaluation of
    Histogram.Percentile for 0 < p < 1. Left out where the CDF is flat at
    t = p*Total (t equals the cumulative count before an empty bin): every
    x on the flat part has CDF p, Percentile documents the left end and
    scipy's ppf takes the right end (or either end, by its own rounding of
    the cumulative sums)."""
    total = int(sum(counts))
    if total == 0:
        return {}
    cum = np.cumsum(counts)
    flat = {int(cum[i]) for i in range(len(counts) - 1) if counts[i + 1] == 0}
    rv = scipy.stats.rv_histogram((np.asarray(counts), np.asarray(edges)), density=False)
    out = {}
    for p in HIST_PS:
        if not 0 < p < 1:
            continue
        t = Fraction(p * float(total))
        if t.denominator == 1 and int(t) in flat:
            continue
        out[repr(p)] = jf(rv.ppf(p))
    return out


HIST_PS = [-0.5, 0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0, 1.5]


PCT_PS = [0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]


def order_stats(clean, f32path):
    """numpy order statistics at the rank Percentile's histogram walk
    targets: j = max(1, ceil(t)) with t = p*n computed as the code does
    (float32(float32(p)*float32(n)) on the []float32 path, p*n in float64
    otherwise). Percentile's doc bounds its distance from this value."""
    srt = np.sort(clean)
    n = len(srt)
    out = {}
    for p in PCT_PS:
        if f32path:
            t = Fraction(float(np.float32(np.float32(p) * np.float32(n))))
        else:
            t = Fraction(p * float(n))
        j = max(1, math.ceil(t))
        out[repr(p)] = jf(srt[j - 1])
    return out


def shortest_cycle(vals):
    """Length of the shortest pattern whose repetition is vals (len(vals)
    if there is none)."""
    n = len(vals)
    for k in range(1, n // 2 + 1):
        if n % k == 0 and all(vals[i] == vals[i % k] for i in range(n)):
            return k
    return n


def digest_f64(values):
    """sha256 of float64 values as big-endian IEEE bits, with a sample of
    explicit values (first two, middle, last two) for debugging."""
    a = np.asarray(values, dtype=">f8")
    n = len(a)
    idx = sorted({i for i in (0, 1, n // 2, n - 2, n - 1) if 0 <= i < n})
    return {"sha256": hashlib.sha256(a.tobytes()).hexdigest(),
            "sample": {str(i): jf(a[i]) for i in idx}}


def compute_case(name, data):
    r = {"dtype": str(data.dtype), "len": int(len(data))}
    if name in GOTYPE:
        r["gotype"] = GOTYPE[name]
    # The []float32 Siril path (a named float32 type takes the float64 path).
    f32 = data.dtype == np.float32 and name not in GOTYPE
    if name in FORMULA_CASES:
        kind, seed = FORMULA_CASES[name]
        r["formula"] = {"kind": kind, "seed": seed}
        # The Go test checks it regenerated these bits before anything else.
        be = ">f4" if data.dtype == np.float32 else ">f8"
        r["input_sha256"] = hashlib.sha256(data.astype(be).tobytes()).hexdigest()
    elif data.dtype == np.float32:
        # float32 inputs as IEEE bit patterns (compact and exact).
        r["input_f32"] = " ".join(f"{u:08x}" for u in data.view(np.uint32))
    elif data.dtype.kind in "iu":
        r["input_int"] = " ".join(str(int(v)) for v in data.flat)
    else:
        vals = [jf(x) for x in data.flat]
        runs = []
        for v in vals:
            if runs and runs[-1][0] == v:
                runs[-1][1] += 1
            else:
                runs.append([v, 1])
        cycle = shortest_cycle(vals)
        if len(runs) * 4 < len(vals):
            # Mostly repeated values: run-length encoded [value, count].
            r["input_runs"] = runs
        elif cycle * 4 < len(vals):
            # A short pattern repeated: [pattern, repeats].
            r["input_cycle"] = [vals[:cycle], len(vals) // cycle]
        else:
            r["input"] = vals
    x64 = data.astype(np.float64)
    if data.dtype.kind in "iu":
        # The Go code converts with float64(v), round to nearest even;
        # numpy's cast must agree (Python's float(int) is correctly rounded).
        if any(float(int(v)) != float(f) for v, f in zip(data.flat, x64)):
            raise SystemExit("numpy integer->float64 cast is not correctly rounded")
    clean = x64[~np.isnan(x64)]  # what the Go functions see (NaN skipped)
    finite = clean[np.isfinite(clean)]
    n = len(clean)

    om, os_ = onepass_meanstdev(clean)
    r["onepass_mean"] = jf(om)
    r["onepass_stdev"] = jf(os_)
    # The same formula with the multiply-adds fused (see onepass_meanstdev).
    for key, acc, sub in (("fma_acc", True, False), ("fma_all", True, True)):
        fm, fs = onepass_meanstdev(clean, acc, sub)
        if not (fm == om or (math.isnan(fm) and math.isnan(om))):
            raise SystemExit("fusion changed the mean, which has no multiply-add")
        r["onepass_stdev_" + key] = jf(fs)

    if n == 0:
        r["min"] = r["max"] = r["mean"] = r["stdev"] = None
    else:
        r["min"] = jf(clean.min())
        r["max"] = jf(clean.max())
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r["mean"] = jf(np.mean(clean))
            r["stdev"] = jf(np.std(clean, ddof=1)) if n >= 2 else None
        if len(finite) == n:
            fr = [Fraction(float(v)) for v in clean]
            m = sum(fr) / n
            r["exact_mean"] = float(m)
            if n >= 2:
                var = sum((v - m) ** 2 for v in fr) / (n - 1)
                r["exact_stdev"] = exact_sqrt(var)
            if data.dtype.kind in "iu" and n >= 2 and any(int(v) != int(float(v)) for v in data.flat):
                # Exact statistics of the integers themselves, before the
                # float64 conversion MeanStdevTwoPass documents.
                fi = [Fraction(int(v)) for v in data.flat]
                mi = sum(fi) / n
                r["exact_mean_int"] = float(mi)
                r["exact_stdev_int"] = exact_sqrt(sum((v - mi) ** 2 for v in fi) / (n - 1))
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if f32:
                    rng_ok = np.isfinite(np.float32(clean.max()) - np.float32(clean.min()))
                else:
                    rng_ok = np.isfinite(clean.max() - clean.min())
            # Percentile is NaN where max-min overflows (no order-statistic
            # bound; TestPercentile_RangeOverflow).
            if rng_ok:
                r["order_stats"] = order_stats(clean, f32)

    if len(data) >= 10 and len(finite) > 0:
        r["sigma_clip"] = [compute_sigma_clip(finite, f32, *v) for v in SIGMA_VARIANTS]
        if name in SIGMA_EXTRA:
            r["sigma_clip_extra"] = [compute_sigma_clip(finite, f32, *v) for v in SIGMA_EXTRA[name]]

    if len(data) > 0:
        nbins = min(256, max(10, len(data) // 5))
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                counts, edges = np.histogram(finite, bins=nbins)
        except ValueError as e:
            # numpy cannot bin this range; BuildHistogram documents its
            # own fallback for it, which the Go test checks.
            r["histogram"] = {"nbins": nbins, "numpy_error": str(e)}
            return r
        total = int(counts.sum())
        cdf = (np.cumsum(counts) / total) if total else np.zeros(nbins)
        r["histogram"] = {
            "nbins": nbins,
            "counts": [int(c) for c in counts],
            "edges": digest_f64(edges),
            "cdf": digest_f64(cdf),
            "percentiles": {repr(p): jf(hist_percentile_exact(counts, edges, p)) for p in HIST_PS},
            "scipy_ppf": scipy_hist_ppf(counts, edges),
        }
    return r


def main():
    if not os.access(SREF, os.X_OK):
        raise SystemExit(f"{SREF} not built: run `make -C {os.path.join(HERE, 'sref')}` first")
    os.makedirs(OUTDIR, exist_ok=True)
    print("Generating stats golden fixtures")
    out = {}
    for name, data in gen_test_arrays().items():
        print(f"  {name}: len={len(data)} dtype={data.dtype}")
        out[name] = compute_case(name, data)
    # One case per line: compact, but diffs stay per case. Committed
    # gzipped (stats_golden.json.gz: level 9, mtime 0); the Go tests read
    # it through testutil.ReadGolden and TestGoldensReproduce compares the
    # decompressed bytes.
    path = os.path.join(OUTDIR, "stats_golden.json.gz")
    names = list(out)
    text = "{\n" + "".join(
        f"{json.dumps(name)}: {json.dumps(out[name], allow_nan=False)}{',' if i + 1 < len(names) else ''}\n"
        for i, name in enumerate(names)) + "}\n"
    with open(path, "wb") as f:
        f.write(gzip.compress(text.encode(), compresslevel=9, mtime=0))
    print(f"Written to {path}")
    path = os.path.join(OUTDIR, "histogram_ranges.json")
    with open(path, "w") as f:
        f.write("[\n")
        rows = histogram_ranges_golden()
        for i, row in enumerate(rows):
            f.write(json.dumps(row, separators=(",", ":")) + ("," if i + 1 < len(rows) else "") + "\n")
        f.write("]\n")
    print(f"Written to {path}")
    path = os.path.join(OUTDIR, "histogram_edges.json")
    with open(path, "w") as f:
        f.write("[\n")
        rows = histogram_edges_golden()
        for i, row in enumerate(rows):
            f.write(json.dumps(row, separators=(",", ":")) + ("," if i + 1 < len(rows) else "") + "\n")
        f.write("]\n")
    print(f"Written to {path}")
    path = os.path.join(OUTDIR, "mean_rounding.json")
    with open(path, "w") as f:
        json.dump(mean_rounding_golden(), f, sort_keys=True, separators=(",", ":"))
        f.write("\n")
    print(f"Written to {path}")


if __name__ == "__main__":
    main()
