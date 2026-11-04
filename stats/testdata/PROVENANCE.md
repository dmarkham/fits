# stats/testdata provenance

## golden/stats_golden.json (committed as stats_golden.json.gz)

Generator: `gen_golden.py` (this directory). It writes the JSON gzipped
(`gzip.compress`, level 9, mtime 0) to keep the repository inside its size
budget; the tests read it through `testutil.ReadGolden`. Re-running it
reproduces the committed JSON byte-for-byte (`TestGoldensReproduce/gen_golden.py`
compares the decompressed bytes).

    make -C sref                       # the Siril median used for SigmaClip
    /path/to/venv/bin/python gen_golden.py

Tool versions used for the committed file:

- Python 3.11.11 (pyenv), venv at `fits-proof-audit/venv`
- numpy 2.4.6 — mean, `std(ddof=1)`, min/max, `numpy.histogram`, CDF,
  order statistics (`np.sort`) for the Percentile one-bin bound
- scipy 1.17.1 — `scipy.stats.rv_histogram(...).ppf` (`scipy_ppf`), an
  independent evaluation of `Histogram.Percentile` for 0 < p < 1 where
  the histogram CDF is not flat
- astropy 8.0.2 — `astropy.stats.SigmaClip` (SigmaClip truth: survivors
  and `_niterations`), run with the exactly rounded mean and population
  std as cenfunc/stdfunc (or Siril's median), for 13 variants per case:
  symmetric k = 2, 3, 5 with maxiters 10, maxiters=None with symmetric and
  asymmetric (2/4, 4/1.5) bounds, and maxiters 1
- Python `fractions` / `decimal` (stdlib) and exact integer arithmetic —
  exact mean and sample stdev (MeanStdevTwoPass truth), exact mean and
  population std of the sigma-clip survivors, and the exact
  Histogram.Percentile definition
- `sref/sref median-f32|median-f64` — Siril's histogram median compiled
  from the Siril 1.4.0 source (see `sref/PROVENANCE.md`), the centre for
  `SigmaClip(..., CenterMedian)`

The arrays keep their original seeds (numpy `RandomState(42)` for the
cases that existed before, `RandomState(7)` for the next ones,
`RandomState(8)` for the element-type and sigma-clip cases added in round
2: int8, uint32, uint64 and int64 beyond 2^53 and 2^63, a named float32
type, and a geometric tail that needs more than 100 clipping iterations).

Size: `float32_large` (100000 values) is not committed. It used to be
`RandomState(42).normal(100, 15, 100000)` stored as 909857 bytes of hex;
it is now the `astro_u16` xorshift32 formula (seed 21), implemented in
`gen_golden.py` (`formula_f32`) and in Go (`formulaValues`, shared with the
Siril harness), so both sides generate the same float32 values. The numpy
draw is still made and discarded so every later array is unchanged. The
file went from 1227271 to under 480000 bytes.

History: the committed file was originally produced by this script with
`ddof=0`, then hand-edited to sample stdev (`ddof=1`) values in 9db8930,
so re-running the script broke `TestMeanStdev`. With `ddof=1` the old
script reproduced the old file byte-for-byte except `float32_single`,
where numpy's `std(ddof=1)` is NaN and the hand edit had written 0.0
(MeanStdev's documented n=1 convention). The script now writes `null`
there and the Go test checks the n<2 convention directly.

Round 2 (second pass):

- `float64_offset_noise` (3000 values) is the `uniform53_1e9` formula
  (seed 23): 1e9 + a 53-bit uniform value, rounded once, in
  `gen_golden.py` (`formula_f64`) and Go (`formulaValues64`). It was
  `1e9 + RandomState(7).uniform(0, 1, 3000)`, 73 KB of JSON; the draw is
  still made and discarded.
- float64 inputs that are mostly runs of one value are written as
  `input_runs` ([value, count] pairs): `float64_big_after_small`.
- New cases: `float64_all_rejected` and `float32_all_rejected`
  ({0 x5, 1 x5, 100}) with `sigma_clip_extra`, astropy runs at k = 0.5
  with maxiters None / 10 / 2 / 1 that reject every value (astropy's
  `_niterations` counts the iteration on the empty set: 3, 3, 2);
  `float64_range_overflow` ({-1e308, 1e308, 5}: max-min overflows, numpy
  raises).
- The file went from 472648 to 416694 bytes.

Round 2 (third pass, size):

- Four more inputs are formula inputs (not committed), each with
  `input_sha256` (sha256 of its big-endian IEEE bits) that the Go loader
  checks before any test runs: `float64_random` (`uniform53`, seed 24,
  shared with the Siril harness), `float32_adu_offset` (`adu_f32`, seed
  25: 30000 + k/256, k a sum of four bytes minus 510, exact in float32),
  `float32_outliers` (`sky_hot`, seed 26: sky 0.05 + k*2^-17 with 1%
  hot pixels 0.9 + j*2^-16; sigma clipping at k = 3 rejects the 49 hot
  ones), `float64_wide` (`wide_i32`, seed 27: int32 * 2^467, up to
  ~1.4e150). They were `RandomState(7)` uniform / normal draws (the draws
  are still made and discarded, so later arrays are unchanged). Go:
  `formulaValues` / `formulaValues64` in `siril_test.go`.
- `float64_offset_1e9` ({1e9, 1e9+1} x 500) is `input_cycle`
  ([pattern, repeats]).
- Histogram edges and CDF are `{"sha256", "sample"}` (`digest_f64`): the
  sha256 of the numpy float64 values as big-endian bits plus five
  explicit values; the counts and percentiles stay explicit.
- New: `scipy_ppf` per histogram (scipy 1.17.1).
- The file went from 422992 (with `scipy_ppf`) to 218148 bytes.

## golden/mean_rounding.json

Also written by `gen_golden.py`: the correctly rounded mean (Python exact
integer arithmetic; `int / int` true division is correctly rounded) of
20000 xorshift32 sets of 2..64 integers ±(2^e + k), e in 40..62,
|k| < 2^23, each value first rounded to float64. The inputs are not
committed (`mean_rounding_sets`, mirrored by `meanRoundingSets` in
`mean_rounding_test.go`); the file holds the sha256 of all 20000 expected
means as big-endian float64 bits, plus 35 explicit samples. 1137 bytes.

## golden/histogram_ranges.json

Also written by `gen_golden.py` (`HIST_RANGES`): numpy 2.4.6
`numpy.histogram` on 17 inputs at the limits of float64, with counts and
edges as float64 bit patterns, or numpy's exception (ValueError "Too many
bins for data range", or IndexError for an overflowing range with one
bin). 2626 bytes.

## golden/histogram_edges.json

Also written by `gen_golden.py` (`histogram_edges_golden`): numpy 2.4.6
`numpy.histogram` on values on its bin edges and one float64 either side
(`np.nextafter`), for 26 (lo, hi, nbins) ranges (two fixed, 24 from
`RandomState(5)`), plus the explicit input {0, 0.3, 0.6, 0.7, 1} with 10
bins. These are the inputs where numpy's computed index is one off and
numpy corrects it; the golden records how often each correction fired
(numpy's index computation replayed in the script, checked against
`np.histogram`): 106 decrements, 238 increments. The edge-case inputs are
not committed; the Go test builds them from the committed numpy edges
with `math.Nextafter`. 19556 bytes.

## sref/golden/siril_stats_golden.json

See `sref/PROVENANCE.md`.
