#!/usr/bin/env python3
"""Generate every Siril-derived stretch golden.

Truth source: the VERBATIM Siril 1.4.0 C/C++ functions compiled in sref/
(`make -C sref`), driven exactly like the corresponding siril-cli commands.
This script only builds inputs, runs the oracle, and serialises its float32
output in the compact exact format of goldenfmt.py; no expected value is
computed here or in Go.

It writes (under --out, default golden/):

  inputs/*.f32           the committed input images (little-endian float32)
  stretch_<preset>.json  the 18 DefaultPresets commands on inputs/synth.f32;
                         each digest must equal the digest of the original
                         siril-cli output (golden/siril_cli_sha256.json,
                         siril_cli_digests.py) or this script fails
  sref/<case>.json       every other oracle case, listed in sref/MANIFEST

    make -C sref goldens     # rewrite golden/
    make -C sref verify      # regenerate into a scratch dir, diff -r
    python gen_siril_goldens.py [--out DIR]

Re-running reproduces every file byte for byte (TestSirilGoldensReproduce).
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

from goldenfmt import Inputs, arr32, canon32, const_input, dump, f32hex, f64hex, hx32, sha32

HERE = os.path.dirname(os.path.abspath(__file__))
SREF = os.environ.get("SREF", os.path.join(HERE, "sref", "sref"))
SREF_CLAHE = os.environ.get("SREF_CLAHE", os.path.join(HERE, "sref", "sref_clahe"))
# ght.c built with double-rounded logf/powf/expf: locates libm-sensitive
# pixels only (libm_record), never a source of expected values.
SREF_LIBM = os.environ.get("SREF_LIBM", os.path.join(HERE, "sref", "sref_libm"))
SCRATCH = tempfile.mkdtemp(prefix="sref-")
F32 = np.float32


# ---------------------------------------------------------------- inputs --

def make_synthetic(nx=64, ny=48):
    """gen_golden.py:make_synthetic_fits minus the FITS write: the input the
    siril-cli goldens were made from (checked against its siril-cli digest)."""
    rng = np.random.RandomState(42)
    data = np.zeros((3, ny, nx), dtype=np.float32)
    for ch in range(3):
        gradient = np.linspace(0.001, 0.3, nx).reshape(1, nx)
        gradient = np.broadcast_to(gradient, (ny, nx)).copy()
        gradient += 0.02 * (ch - 1)
        noise = rng.normal(0, 0.01, (ny, nx)).astype(np.float32)
        data[ch] = np.clip(gradient + noise, 0, 1)
    data[:, 0, 0] = 0
    data[:, ny - 1, nx - 1] = 0
    return data


def make_bright(nx=48, ny=40):
    """Bright background (every channel median > 0.5) so Siril's inverted
    autostretch branch (mtf.c:299 / :414) is taken; some exact zeros (excluded
    from the stats) and saturated pixels."""
    rng = np.random.RandomState(7)
    data = np.zeros((3, ny, nx), dtype=np.float32)
    for ch in range(3):
        grad = np.linspace(0.62, 0.86, nx).reshape(1, nx) + 0.015 * ch
        data[ch] = np.clip(grad + rng.normal(0, 0.02, (ny, nx)), 0, 1).astype(np.float32)
    data[:, 3, 5] = 0
    data[:, 10, 40] = 1
    data[0, 20, 7] = 0
    return data


def make_mixed(nx=40, ny=32):
    """R bright, G and B dark: only 1 of 3 channels has median > 0.5, so
    Siril keeps the normal branch for ALL channels (mtf.c:394)."""
    rng = np.random.RandomState(11)
    d = np.zeros((3, ny, nx), dtype=np.float32)
    x = np.linspace(0, 1, nx).reshape(1, nx)
    d[0] = np.clip(0.7 + 0.1 * x + rng.normal(0, 0.02, (ny, nx)), 0, 1)
    d[1] = np.clip(0.05 + 0.1 * x + rng.normal(0, 0.01, (ny, nx)), 0, 1)
    d[2] = np.clip(0.08 + 0.05 * x + rng.normal(0, 0.01, (ny, nx)), 0, 1)
    return d.astype(np.float32)


def make_clahe_clip(nch, nx, ny, seed):
    """Smooth sky (gradient + low noise) with a few stars: in every tile most
    pixels fall in a few hundred 16-bit bins, so histogram bins exceed the
    clip limit and the clip limit decides the output."""
    rng = np.random.RandomState(seed)
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
    d = np.zeros((nch, ny, nx), dtype=np.float32)
    stars = rng.uniform(0, 1, (6, 2)) * [nx, ny]
    for ch in range(nch):
        sky = 0.08 + 0.04 * xx / nx + 0.02 * yy / ny + 0.01 * ch
        img = sky + rng.normal(0, 0.002, (ny, nx))
        for sx, sy in stars:
            img += 0.6 * np.exp(-((xx - sx) ** 2 + (yy - sy) ** 2) / 8.0)
        d[ch] = np.clip(img, 0, 1)
    return d


def median_image(target, seed):
    """Mono 33x33 image whose median (of the non-zero pixels) is exactly
    `target`: 400 pixels below, 289 equal to target, 400 above."""
    rng = np.random.RandomState(seed)
    v = np.concatenate([rng.uniform(0.2, 0.45, 400), np.full(289, target, np.float64),
                        rng.uniform(0.56, 0.8, 400)]).astype(np.float32)
    v[400:689] = np.float32(target)
    rng.shuffle(v)
    return v.reshape(1, 33, 33)


def clahe_tile_total(nx, ny, tile):
    """OpenCV CLAHE_Impl::apply (clahe.cpp): the image is padded to a
    multiple of the grid, tileSize = padded size / grid."""
    px = nx if nx % tile == 0 else nx + tile - nx % tile
    py = ny if ny % tile == 0 else ny + tile - ny % tile
    return (px // tile) * (py // tile)


# ---------------------------------------------------------------- oracle --

def run(prog, op, data, args):
    nch, ny, nx = data.shape
    fin = os.path.join(SCRATCH, "in.f32")
    fout = os.path.join(SCRATCH, "out.f32")
    np.ascontiguousarray(data, dtype="<f4").tofile(fin)
    if prog in (SREF, SREF_LIBM):
        cmd = [prog, op, str(nch), str(nx), str(ny), fin, fout] + [repr(a) if isinstance(a, float) else str(a) for a in args]
    else:
        cmd = [prog, str(nch), str(nx), str(ny), fin, fout] + [repr(a) if isinstance(a, float) else str(a) for a in args]
    subprocess.run(cmd, check=True)
    if op == "asparams":
        rows = []
        for line in open(fout):
            r, s, m, h = line.split()
            # %.9g round-trips float32 exactly (and prints nan / -nan)
            rows.append([int(r)] + [f32hex(F32(float(v))) for v in (s, m, h)])
        return rows
    if op == "ghtcoef":
        return open(fout).read().split()
    return np.fromfile(fout, dtype="<f4").reshape(data.shape)


def libm_record(data, args, got):
    """Pixels where glibc's float libm decides Siril's GHT output: those
    where the verbatim oracle and the same code with double-rounded
    logf/powf/expf (sref_libm) disagree. The golden lists them with
    Siril's value, so the Go test can accept exactly there a value that a
    1-ULP different logf/powf/expf result explains (ght_test.go)."""
    alt = run(SREF_LIBM, "ght", data, args)
    a, b = canon32(got.reshape(-1)), canon32(alt.reshape(-1))
    idx = [int(i) for i in np.where(a != b)[0]]
    return {"idx": idx, "bits": hx32(a[idx])}


def opencv_version():
    return subprocess.run([SREF_CLAHE, "--version"], capture_output=True, text=True, check=True).stdout.strip()


# DefaultPresets (stretch/pipeline.go): name -> siril command. The Go test
# requires each golden's command to equal Preset.Command, so a drift between
# this table and the Go presets fails there.
PRESETS = [
    ("kstars1", "autostretch -2.8 0.25"), ("kstars2", "autostretch -1.5 0.25"),
    ("kstars3", "autostretch -1.0 0.125"), ("kstars4", "autostretch -0.5 0.125"),
    ("kstars5", "autostretch -2.0 0.125"), ("kstars6", "autostretch -4.0 0.125"),
    ("kstars7", "autostretch -5.5 0.25"), ("custom_dark", "autostretch -6.5 0.03"),
    ("custom_extreme", "autostretch -7.5 0.02"), ("custom_ultra", "autostretch -9.0 0.01"),
    ("custom_mask", "autostretch -12.0 0.005"), ("auto_linked", "autostretch -linked"),
    ("asinh_light", "asinh 1 0"), ("asinh_medium", "asinh 50 0"), ("asinh_strong", "asinh 500 0"),
    ("ght_mild", "ght -D=2.0 -B=0.0 -SP=0.5"), ("ght_strong", "ght -D=5.0 -B=2.0 -SP=0.3"),
    ("linear", "linstretch -BP=0.05"),
]


def command_to_oracle(cmd):
    """siril command -> (oracle op, args), with siril-cli's defaults
    (command.c process_autostretch: -2.8 0.25; process_ght_args: LP 0, HP 1)."""
    w = cmd.split()
    if w[0] == "autostretch":
        if w[1] == "-linked":
            return "autostretch", [1] + ([float(w[2]), float(w[3])] if len(w) == 4 else [-2.8, 0.25])
        return "autostretch", [0, float(w[1]), float(w[2])]
    if w[0] == "asinh":
        return "asinh", [float(w[1]), float(w[2]), 0]
    if w[0] == "ght":
        v = {"-LP": 0.0, "-HP": 1.0}
        for kv in w[1:]:
            k, x = kv.split("=")
            v[k] = float(x)
        return "ght", [v["-D"], v["-B"], v["-LP"], v["-SP"], v["-HP"]]
    if w[0] == "linstretch":
        return "linstretch", [float(w[1].split("=")[1])]
    raise ValueError(cmd)


class Writer:
    def __init__(self, out):
        self.out = out
        self.sref = os.path.join(out, "sref")
        os.makedirs(self.sref, exist_ok=True)
        self.inputs = Inputs(out)
        self.manifest = []

    def case(self, name, inp, op, args, command, note, prog=None, **extra):
        """Run the oracle on input record/array `inp` and write sref/<name>.json."""
        rec_in, data = inp
        got = run(prog or SREF, op, data, args)
        rec = {"name": name, "op": op, "args": args, "command": command, "note": note, "input": rec_in}
        if any(isinstance(a, float) and not np.isfinite(a) for a in args):
            # JSON has no NaN/Inf: such args are null, and args_f64 holds
            # every argument as big-endian float64 bits in hex.
            rec["args"] = [None if isinstance(a, float) and not np.isfinite(a) else a for a in args]
            rec["args_f64"] = [f64hex(float(a)) for a in args]
        rec.update(extra)
        rec["out"] = arr32(got)
        if op == "ght":
            rec["libm"] = libm_record(data, args, got)
        if op == "autostretch":
            rec["params"] = run(SREF, "asparams", data, args)
        if prog == SREF_CLAHE:
            rec["opencv"] = opencv_version()
        dump(os.path.join(self.sref, name + ".json"), rec)
        self.manifest.append(name)
        print(f"  {name}: {got.shape} nan-out {int(np.isnan(got).sum())}")
        return got, rec

    def finish(self):
        with open(os.path.join(self.sref, "MANIFEST"), "w") as f:
            f.write("\n".join(self.manifest) + "\n")
        print(f"Done: {len(self.manifest)} oracle goldens in {self.sref}")


# ----------------------------------------------------------------- cases --

def siril_cli_cases(w, synth):
    """The 18 preset goldens. The oracle must reproduce the original siril-cli
    output (digests in golden/siril_cli_sha256.json) bit for bit."""
    cli = json.load(open(os.path.join(HERE, "golden", "siril_cli_sha256.json")))
    if sha32(synth) != cli["input"]["sha256"]:
        sys.exit("make_synthetic() no longer reproduces the siril-cli input")
    ref, _ = w.inputs.ref("synth.f32")
    bad = 0
    for name, command in PRESETS:
        op, args = command_to_oracle(command)
        got = run(SREF, op, synth, args)
        want = cli["presets"][name]
        if command != want["command"] or sha32(got) != want["sha256"]:
            print(f"  oracle vs siril-cli {name}: DIFFERENT")
            bad += 1
            continue
        print(f"  oracle vs siril-cli {name:15s} identical")
        rec = {"name": name, "command": command, "op": op, "args": args,
               "source": "Siril 1.4.0 oracle (sref); sha256 equals the original siril-cli output, "
                         "siril_cli_sha256.json blob " + want["blob"],
               "input": ref, "out": arr32(got)}
        if op == "ght":
            rec["libm"] = libm_record(synth, args, got)
        dump(os.path.join(w.out, f"stretch_{name}.json"), rec)
    if bad:
        sys.exit(f"oracle does not reproduce {bad} siril-cli goldens")


def basic_cases(w, ins):
    """Per-operation cases beyond the presets."""
    R = w.inputs.ref
    for name, src, chans, op, args, command, note in [
        ("asinh_mono_b50", "synth.f32", [0], "asinh", [50, 0, 0], "asinh 50", "mono image"),
        ("asinh_mono_b50_off", "synth.f32", [0], "asinh", [50, 0.02, 0], "asinh 50 0.02", "mono, offset"),
        ("asinh_rgb_b20_off", "synth.f32", None, "asinh", [20, 0.015, 0], "asinh 20 0.015", "offset != 0"),
        ("asinh_rgb_b20_human", "synth.f32", None, "asinh", [20, 0, 1], "asinh -human 20", "human luminance"),
        ("asinh_rgb_b200_human_off", "synth.f32", None, "asinh", [200, 0.01, 1], "asinh -human 200 0.01", "human + offset"),
        ("asinh_bright_b500", "bright.f32", None, "asinh", [500, 0, 0], "asinh 500", "channels clip at 1"),
        ("autostretch_bright", "bright.f32", None, "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "inverted branch, unlinked"),
        ("autostretch_bright_linked", "bright.f32", None, "autostretch", [1, -2.8, 0.25], "autostretch -linked -2.8 0.25", "inverted branch, linked"),
        ("autostretch_bright_mono", "bright.f32", [1], "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "inverted branch, mono"),
        ("autostretch_bright_mono_strong", "bright.f32", [2], "autostretch", [0, -0.5, 0.125], "autostretch -0.5 0.125", "inverted, c1 not clipped"),
        ("autostretch_mixed", "mixed.f32", None, "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "1 of 3 channels bright: normal branch for all"),
        ("autostretch_mixed_linked", "mixed.f32", None, "autostretch", [1, -2.8, 0.25], "autostretch -linked -2.8 0.25", "1 of 3 bright, linked"),
        ("autostretch_bright_c1_clip", "bright.f32", [0], "autostretch", [0, -100, 0.25], "autostretch -100 0.25", "inverted branch, c1 clipped to 1"),
        ("autostretch_m2_le_0", "synth.f32", [0], "autostretch", [0, 5.0, 0.25], "autostretch 5 0.25", "c0 > median: m2 <= 0 -> MTF(m2)=0 -> midtones 0"),
        ("mtf_hl_08", "synth.f32", None, "mtf", [0.01, 0.2, 0.8], "mtf 0.01 0.2 0.8", "MTF highlights != 1 (function level)"),
        ("mtf_hl_03", "synth.f32", None, "mtf", [0.0, 0.7, 0.3], "mtf 0 0.7 0.3", "MTF highlights=0.3, midtones>0.5"),
        ("ght_b_neg1", "synth.f32", None, "ght", [3.0, -1.0, 0.0, 0.1, 1.0], None, "B=-1 (log branch); siril-cli rejects B<0, GUI allows"),
        ("ght_b_neg05", "synth.f32", None, "ght", [3.0, -0.5, 0.0, 0.1, 1.0], None, "B=-0.5; siril-cli rejects B<0, GUI allows"),
        ("ght_lp_hp", "synth.f32", None, "ght", [4.0, 1.5, 0.05, 0.12, 0.25], "ght -D=4 -B=1.5 -LP=0.05 -SP=0.12 -HP=0.25", "LP>0 and HP<1"),
        ("ght_b0_lp_hp", "synth.f32", None, "ght", [2.5, 0.0, 0.02, 0.1, 0.2], "ght -D=2.5 -B=0 -LP=0.02 -SP=0.1 -HP=0.2", "B=0 with LP>0, HP<1"),
        ("ght_bneg1_lp_hp", "synth.f32", None, "ght", [2.0, -1.0, 0.03, 0.1, 0.22], None, "B=-1 with LP>0, HP<1"),
        ("ght_bright_hp", "bright.f32", None, "ght", [1.5, 3.0, 0.1, 0.5, 0.9], "ght -D=1.5 -B=3 -LP=0.1 -SP=0.5 -HP=0.9", "pixels above HP and at 1"),
        ("linstretch_bright", "bright.f32", None, "linstretch", [0.6], "linstretch -BP=0.6", "large black point"),
    ]:
        w.case(name, R(src, [["channels", chans]] if chans else None), op, args, command, note)

    # No good pixel in ONE channel of an RGB image: Siril's
    # compute_all_channels_statistics_single_image returns -1, and
    # find_unlinked_midtones_balance (mtf.c:386-396) / find_linked_ (mtf.c:
    # 276-281) set (0, 0.2, 1) for EVERY channel, which process_autostretch
    # (command.c:3349-3360) applies unconditionally, so the good channels
    # change too.
    zero = f32hex(0)
    for name, src, c, args, command, note in [
        ("autostretch_zero_g", "synth.f32", 1, [0, -2.8, 0.25], "autostretch -2.8 0.25", "G all zero, unlinked: fallback for all channels"),
        ("autostretch_zero_g_linked", "synth.f32", 1, [1, -2.8, 0.25], "autostretch -linked -2.8 0.25", "G all zero, linked: fallback"),
        ("autostretch_bright_zero_b", "bright.f32", 2, [0, -2.8, 0.25], "autostretch -2.8 0.25", "bright, B all zero, unlinked"),
        ("autostretch_bright_zero_b_linked", "bright.f32", 2, [1, -1.0, 0.125], "autostretch -linked -1.0 0.125", "bright, B all zero, linked"),
    ]:
        _, rec = w.case(name, R(src, [["fill", c, zero]]), "autostretch", args, command, note)
        if any(p[0] == 0 for p in rec["params"]):
            sys.exit(f"{name}: oracle statistics did not fail")

    # Every DefaultPreset on the bright and mixed inputs (synth is covered by
    # the siril-cli goldens). On mixed, Siril's unlinked autostretch differs
    # from per-channel autostretch (one bright channel).
    for src, tag in [("bright.f32", "bright"), ("mixed.f32", "mixed")]:
        for pname, command in PRESETS:
            op, args = command_to_oracle(command)
            w.case(f"preset_{pname}_{tag}", R(src), op, args, command, f"DefaultPresets {pname} on {src}", preset=pname)

    # Degenerate statistics: constant image (MAD = 0 -> 0.001 guard) and an
    # all-zero image (no good pixels -> Siril falls back to 0 / 0.2 / 1).
    for name, value in [("autostretch_const", 0.3), ("autostretch_allzero", 0.0)]:
        w.case(name, const_input(value, (1, 16, 16)), "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "degenerate stats")


def clahe_cases(w, ins):
    """cvClahe_float against the system OpenCV."""
    R = w.inputs.ref
    for name, inp, args, command, note in [
        ("clahe_rgb_3_8", R("synth.f32"), [3.0, 8], "clahe 3 8", "fitview CLAHE preset"),
        ("clahe_rgb_2_4", R("synth.f32"), [2.0, 4], "clahe 2 4", "other clip/tile"),
        ("clahe_mono_3_8", R("synth.f32", [["channels", [1]]]), [3.0, 8], "clahe 3 8", "mono path"),
        ("clahe_rgb_odd", R("synth.f32", [["crop", 3, 32, 5, 42]]), [2.5, 8], "clahe 2.5 8", "37x29: tiles padded (BORDER_REFLECT_101), row tail"),
        ("clahe_bright_rgb", R("bright.f32"), [4.0, 6], "clahe 4 6", "bright image, 48x40"),
        ("clahe_mono_odd", R("bright.f32", [["channels", [0]], ["crop", 0, 40, 1, 46]]), [1.5, 5], "clahe 1.5 5", "mono 45x40"),
    ]:
        w.case(name, inp, "clahe", args, command, note, SREF_CLAHE)

    # Clip limit actually in effect. OpenCV's integer clip limit is
    # int(clip * tileTotal / 65536), at least 1; in the cases above it is
    # < 1 and clamps to 1, so they cannot tell clip limits apart. These
    # cases have an effective limit >= 2, and the oracle run at clip values
    # giving limit - 1 and limit + 1 must change the output (checked here),
    # so the goldens pin the clip-limit arithmetic.
    for name, fname, clip, tile, note in [
        ("clahe_clip_mono", "clahe_clip_mono.f32", 70.0, 2, "mono 128x96, 2x2 tiles of 3072 px"),
        ("clahe_clip_rgb", "clahe_clip_rgb.f32", 70.0, 2, "RGB 128x96, 2x2 tiles of 3072 px"),
        ("clahe_clip_rgb_pad", "clahe_clip_rgb_pad.f32", 300.0, 3, "RGB 90x70 padded to 90x72, 3x3 tiles of 720 px"),
    ]:
        rec_in, data = R(fname)
        nch, ny, nx = data.shape
        tt = clahe_tile_total(nx, ny, tile)
        limit = int(clip * tt / 65536)
        if limit < 2:
            sys.exit(f"{name}: effective clip limit {limit} < 2")
        got, _ = w.case(name, (rec_in, data), "clahe", [clip, tile], f"clahe {clip:g} {tile}",
                        note + f"; effective clip limit {limit}", SREF_CLAHE, clip_limit_int=limit)
        for k in (limit - 1, limit + 1):
            other = run(SREF_CLAHE, "clahe", data, [(k + 0.5) * 65536 / tt, tile])
            if np.array_equal(other.view(np.uint32), got.view(np.uint32)):
                sys.exit(f"{name}: output does not depend on the clip limit ({k} vs {limit})")


def edge_cases(w, ins):
    """NaN input, GHT output > 1, small D, the |B| < 0.001 threshold, B < 0
    coefficient rounding, the median = 0.5 boundary, asinh beta = 0 and
    inputs > 1. Each case checks here that the input really exercises what
    it is for."""
    R = w.inputs.ref
    nan = [["nanpattern"]]

    # NaN input. Siril excludes NaN from the autostretch statistics, its
    # max() macro propagates NaN in asinh, GHT/linstretch clamp NaN to 0
    # (fminf/fmaxf), and the CLAHE Lab conversions follow OpenCV's SSE
    # min/max (NaN -> the second operand).
    for name, src, ops, op, args, command, note, prog in [
        ("nan_autostretch", "synth.f32", nan, "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "NaN input, unlinked", None),
        ("nan_autostretch_linked", "synth.f32", nan, "autostretch", [1, -2.8, 0.25], "autostretch -linked -2.8 0.25", "NaN input, linked", None),
        ("nan_autostretch_mono", "synth.f32", nan + [["channels", [0]]], "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "NaN input, mono", None),
        ("nan_autostretch_bright", "bright.f32", nan, "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "NaN input, inverted branch", None),
        ("nan_asinh_rgb", "synth.f32", nan, "asinh", [50, 0.01, 0], "asinh 50 0.01", "NaN input, RGB", None),
        ("nan_asinh_rgb_human", "synth.f32", nan, "asinh", [20, 0, 1], "asinh -human 20", "NaN input, human weights", None),
        ("nan_asinh_mono", "synth.f32", nan + [["channels", [0]]], "asinh", [50, 0, 0], "asinh 50", "NaN input, mono", None),
        ("nan_linstretch", "synth.f32", nan, "linstretch", [0.05], "linstretch -BP=0.05", "NaN input", None),
        ("nan_mtf", "synth.f32", nan, "mtf", [0.01, 0.2, 1.0], None, "NaN input, MTF (function level)", None),
        ("nan_ght", "synth.f32", nan, "ght", [2.0, 0.0, 0.0, 0.5, 1.0], "ght -D=2.0 -B=0.0 -SP=0.5", "NaN input", None),
        ("nan_clahe_rgb", "synth.f32", nan, "clahe", [3.0, 8], "clahe 3 8", "NaN input, RGB", SREF_CLAHE),
        ("nan_clahe_mono", "synth.f32", nan + [["channels", [1]]], "clahe", [3.0, 8], "clahe 3 8", "NaN input, mono", SREF_CLAHE),
    ]:
        w.case(name, R(src, ops), op, args, command, note, prog)

    # A channel that is entirely NaN has no good pixel, like an all-zero
    # one: Siril's statistics fail and every channel gets (0, 0.2, 1).
    alln = [["fill", 1, f32hex(np.nan)]]
    for name, ops, args, command in [
        ("nan_channel_g", alln, [0, -2.8, 0.25], "autostretch -2.8 0.25"),
        ("nan_channel_g_linked", alln, [1, -2.8, 0.25], "autostretch -linked -2.8 0.25"),
        ("nan_allnan_mono", alln + [["channels", [1]]], [0, -2.8, 0.25], "autostretch -2.8 0.25"),
    ]:
        _, rec = w.case(name, R("synth.f32", ops), "autostretch", args, command, "one channel all NaN: no good pixel, fallback")
        if any(p[0] == 0 for p in rec["params"]):
            sys.exit(f"{name}: oracle statistics did not fail")

    # GHT output above 1: Siril does not clamp GHT's output.
    near1 = [["put", -1, list(range(8)), [f32hex(v) for v in
                                           [1, 1 - 2 ** -24, 1 - 2 ** -23, 0.999, 0.99, 0.97, 0.96, 0.95]]]]
    for name, args, command in [
        ("ght_out_above_1", [0.5, 2.0, 0.0, 0.1, 0.95], "ght -D=0.5 -B=2 -LP=0 -SP=0.1 -HP=0.95"),
        ("ght_out_above_1_b0", [1.0, 0.0, 0.0, 0.1, 0.9], "ght -D=1 -B=0 -LP=0 -SP=0.1 -HP=0.9"),
    ]:
        got, _ = w.case(name, R("bright.f32", near1), "ght", args, command, "output > 1 (not clamped)")
        if not (got > 1).any():
            sys.exit(f"{name}: no output above 1")

    # Small D: float32(expm1(D)) != float32(exp(D)) - 1 for these D.
    for name, args, command in [
        ("ght_small_d_005", [0.05, 1.0, 0.0, 0.2, 1.0], "ght -D=0.05 -B=1 -SP=0.2"),
        ("ght_small_d_001", [0.01, 0.0, 0.0, 0.1, 1.0], "ght -D=0.01 -B=0 -SP=0.1"),
    ]:
        D = args[0]
        if F32(np.expm1(D)) == F32(F32(np.exp(D)) - F32(1)):
            sys.exit(f"{name}: expm1 and exp-1 agree at D={D}")
        w.case(name, R("synth.f32"), "ght", args, command, "small D (expm1)")

    # |B| < 0.001 is treated as 0 (process_ght_args); B = 0.0009 and
    # -0.0009 take the exponential branch, B = 0.0011 the power branch.
    # (siril-cli refuses B < 0, so the negative ones have no command.)
    for name, B in [("ght_b_tiny_pos", 0.0009), ("ght_b_tiny_neg", -0.0009), ("ght_b_above_thresh", 0.0011)]:
        w.case(name, R("synth.f32"), "ght", [3.0, B, 0.0, 0.1, 1.0], f"ght -D=3 -B={B!r} -SP=0.1" if B > 0 else None,
               f"B={B}: threshold |B| < 0.001")

    # B < 0, not -1: (B-1)/B and 1-1/B round differently in float32 for |B|=0.35.
    b = F32(0.35)
    if F32((b - F32(1)) / b) == F32(F32(1) - F32(1) / b):
        sys.exit("B=-0.35: (B-1)/B == 1-1/B in float32")
    w.case("ght_bneg035_lp_hp", R("synth.f32"), "ght", [2.5, -0.35, 0.04, 0.12, 0.25], None,
           "B=-0.35 with LP>0, HP<1; siril-cli rejects B<0, GUI allows")
    # B = -1 with many pixels below LP (the b1 * in segment, no libm call).
    # At D=1, LP=0.08, SP=0.2 the float32 coefficient
    # (1 + D*(SP-LP)) / (D*q) is one ULP away from ((1 + D*(SP-LP)) / D) / q,
    # so the golden pins Siril's grouping.
    w.case("ght_bneg1_big_lp", R("synth.f32"), "ght", [1.0, -1.0, 0.08, 0.2, 1.0], None,
           "B=-1, LP=0.08: many pixels below LP; siril-cli rejects B<0")

    # Inverted-branch boundary (mtf.c:414, median > 0.5 compared in
    # double). Siril's median is the histogram median of
    # findMinMaxPercentile, not an order statistic, so on these images a
    # plateau at 0.4999 already reads as above 0.5 and one at 0.499 does
    # not; the cases sit on both sides of that edge. Branch check: with
    # target 0.25 and a median near 0.5, the normal branch gives midtones
    # MTF(~0.5, 0.25) ~ 0.75, the inverted one 1 - MTF(~0.5, 0.25) ~ 0.25.
    for name, target, inverted in [("autostretch_median_049", F32(0.49), False),
                                   ("autostretch_median_0495", F32(0.495), False),
                                   ("autostretch_median_0499", F32(0.499), False),
                                   ("autostretch_median_04999", F32(0.4999), True),
                                   ("autostretch_median_05", F32(0.5), True)]:
        fname = name.replace("autostretch_", "") + ".f32"
        w.inputs.add(fname, median_image(target, 31))
        _, rec = w.case(name, R(fname), "autostretch", [0, -2.8, 0.25],
                        "autostretch -2.8 0.25", f"plateau median {float(target)!r}")
        mid = np.array([int(rec["params"][0][2], 16)], np.uint32).view(np.float32)[0]
        if (mid < 0.5) != inverted:
            sys.exit(f"{name}: oracle inverted={mid < 0.5}, want {inverted}")

    # asinh with beta = 0 (k = 1). siril-cli rejects beta < 1
    # (command.c:2655); the oracle calls command_asinh directly.
    w.case("asinh_beta0_mono", R("synth.f32", [["channels", [0]]]), "asinh", [0, 0.01, 0], None, "beta=0, mono; siril-cli rejects beta<1")
    w.case("asinh_beta0_rgb", R("synth.f32"), "asinh", [0, 0.01, 0], None, "beta=0, RGB; siril-cli rejects beta<1")
    # mono asinh on input above 1: the min(1, .) clamp.
    hi = R("synth.f32", [["channels", [0]], ["mul", f32hex(4)]])
    if not (hi[1] > 1).any():
        sys.exit("asinh_mono_above_1: no input above 1")
    w.case("asinh_mono_above_1", hi, "asinh", [50, 0, 0], "asinh 50", "mono input up to 1.2: output clamped to 1")


def round2_cases(w, ins):
    """Inputs no earlier golden had: negative pixels, -0, +-Inf, a median of
    exactly 0.5, B = +-0.001 exactly."""
    R = w.inputs.ref
    hx = f32hex
    neg0 = ["put", -1, [5, 6, 7, 70], [hx(-0.0)] * 4]          # -0 pixels in every channel
    neg = [["add", hx(-0.08)], neg0]                          # synth - 0.08: about a quarter negative
    rec_in, data = R("synth.f32", neg)
    frac = float((data < 0).mean())
    if frac < 0.2:
        sys.exit(f"negative input: only {frac:.0%} negative pixels")
    print(f"  negative input: {frac:.1%} of pixels < 0")
    # bright with negative pixels: every channel median stays > 0.5
    negb = [["put", -1, list(range(0, 400, 10)), [hx(-0.3 + 0.007 * k) for k in range(40)]], neg0]
    for name, src, ops, op, args, command, note, prog in [
        ("neg_autostretch", "synth.f32", neg, "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "negative pixels count in the statistics (only 0 and NaN are excluded)", None),
        ("neg_autostretch_linked", "synth.f32", neg, "autostretch", [1, -2.8, 0.25], "autostretch -linked -2.8 0.25", "negative pixels, linked", None),
        ("neg_autostretch_mono", "synth.f32", neg + [["channels", [1]]], "autostretch", [0, -1.0, 0.125], "autostretch -1 0.125", "negative pixels, mono", None),
        ("neg_autostretch_bright", "bright.f32", negb, "autostretch", [0, -2.8, 0.25], "autostretch -2.8 0.25", "negative pixels, inverted branch", None),
        ("neg_autostretch_bright_linked", "bright.f32", negb, "autostretch", [1, -2.8, 0.25], "autostretch -linked -2.8 0.25", "negative pixels, inverted branch, linked", None),
        ("neg_asinh_rgb", "synth.f32", neg, "asinh", [50, 0, 0], "asinh 50", "negative pixels, RGB", None),
        ("neg_asinh_rgb_off", "synth.f32", neg, "asinh", [20, -0.05, 0], "asinh 20 -0.05", "negative pixels, negative offset (siril-cli takes any offset)", None),
        ("neg_asinh_human_off", "synth.f32", neg, "asinh", [100, 0.01, 1], "asinh -human 100 0.01", "negative pixels, human, offset", None),
        ("neg_asinh_mono", "synth.f32", neg + [["channels", [0]]], "asinh", [50, 0.005, 0], "asinh 50 0.005", "negative pixels, mono", None),
        ("neg_linstretch", "synth.f32", neg, "linstretch", [0.02], "linstretch -BP=0.02", "negative pixels", None),
        ("neg_linstretch_bp0", "synth.f32", neg, "linstretch", [0.0], "linstretch -BP=0", "negative pixels and -0, BP = 0: Siril returns +0", None),
        ("neg_mtf", "synth.f32", neg, "mtf", [-0.05, 0.3, 0.9], None, "negative pixels, negative shadows (function level)", None),
        ("neg_ght", "synth.f32", neg, "ght", [3.0, 1.0, 0.0, 0.05, 1.0], "ght -D=3 -B=1 -SP=0.05", "negative pixels and -0 clamp to +0", None),
        ("neg_clahe_rgb", "synth.f32", neg, "clahe", [3.0, 8], "clahe 3 8", "negative pixels, RGB", SREF_CLAHE),
        ("neg_clahe_mono", "synth.f32", neg + [["channels", [2]]], "clahe", [2.0, 4], "clahe 2 4", "negative pixels, mono", SREF_CLAHE),
    ]:
        got, rec = w.case(name, R(src, ops), op, args, command, note, prog)
        if name == "neg_linstretch_bp0":
            u = canon32(got.reshape(3, -1)[:, [5, 6, 7, 70]].reshape(-1))
            if (u != 0).any():
                sys.exit(f"{name}: Siril did not return +0 for the -0 pixels: {u}")

    # A median of exactly 0.5: Siril compares median > 0.5 (mtf.c:299 /
    # :414), so a constant 0.5 image takes the normal branch and one ULP
    # above 0.5 the inverted one. (findMinMaxPercentile's fast exit returns
    # the constant itself.)
    for name, value, shape, linked, inverted in [
        ("autostretch_const_05", 0.5, (1, 10, 10), 0, False),
        ("autostretch_const_05_rgb", 0.5, (3, 10, 10), 0, False),
        ("autostretch_const_05_linked", 0.5, (3, 10, 10), 1, False),
        ("autostretch_const_05_up", float(np.nextafter(F32(0.5), F32(1))), (1, 10, 10), 0, True),
        ("autostretch_const_05_up_linked", float(np.nextafter(F32(0.5), F32(1))), (3, 10, 10), 1, True),
    ]:
        args = [linked, -2.8, 0.25]
        _, rec = w.case(name, const_input(value, shape), "autostretch", args,
                        "autostretch " + ("-linked " if linked else "") + "-2.8 0.25", f"constant {value!r}")
        h = [np.array([int(p[3], 16)], np.uint32).view(np.float32)[0] for p in rec["params"]]
        if any((x < 1) != inverted for x in h):
            sys.exit(f"{name}: oracle branch (highlights {h}) is not inverted={inverted}")

    # B = +-0.001 exactly: process_ght_args compares fabs(B) < 1.e-3f, the
    # float literal promoted to double (0.0010000000475), so the double
    # 0.001 becomes 0 (exponential branch). The GUI compares in float
    # (gui/histogram.c:1723) and keeps 0.001.
    for name, B in [("ght_b_0001", 0.001), ("ght_b_neg0001", -0.001)]:
        if not abs(B) < float(F32(1e-3)):
            sys.exit("threshold arithmetic")
        w.case(name, R("synth.f32"), "ght", [3.0, B, 0.0, 0.1, 1.0], f"ght -D=3 -B={B!r} -SP=0.1" if B > 0 else None,
               f"B={B}: |B| < (double)1e-3f, so B = 0")

    # +-Inf input. Siril's statistics use findMinMaxPercentile, whose
    # histogram index static_cast<uint16_t>(scale * (x - min)) is undefined
    # behaviour for an infinite min or max (NaN or Inf argument), so a
    # channel with an Inf pixel gets a NaN or meaningless median. The Go
    # port excludes +-Inf from the statistics like 0 and NaN instead. Truth
    # for that rule: Siril's parameters for the same image with the Inf
    # pixels set to 0 (which Siril excludes), and Siril's MTF applied with
    # those parameters to the image WITH the Inf pixels (op mtf, per channel
    # for unlinked). Siril's own parameters for the Inf image are recorded
    # as "siril_inf_params" to document the divergence.
    inf, ninf = hx(np.inf), hx(-np.inf)
    infops = [["put", 0, [100, 1500], [inf, inf]], ["put", 1, [200], [ninf]], ["put", 2, [300, 301], [inf, ninf]]]
    zops = [["put", 0, [100, 1500], [hx(0)] * 2], ["put", 1, [200], [hx(0)]], ["put", 2, [300, 301], [hx(0)] * 2]]
    for name, src, ops, z, args, note in [
        ("inf_autostretch", "synth.f32", infops, zops, [0, -2.8, 0.25], "+-Inf pixels, unlinked"),
        ("inf_autostretch_linked", "synth.f32", infops, zops, [1, -2.8, 0.25], "+-Inf pixels, linked"),
        ("inf_autostretch_mono", "synth.f32", infops + [["channels", [0]]], zops + [["channels", [0]]], [0, -2.8, 0.25], "+Inf pixels, mono"),
        ("inf_autostretch_mono_neg", "synth.f32", infops + [["channels", [1]]], zops + [["channels", [1]]], [0, -2.8, 0.25], "-Inf pixel, mono"),
        ("inf_autostretch_bright", "bright.f32", infops, zops, [0, -2.8, 0.25], "+-Inf pixels, inverted branch"),
    ]:
        rec_in, data = R(src, ops)
        _, zdata = R(src, z)
        siril_inf = run(SREF, "asparams", data, args)
        params = run(SREF, "asparams", zdata, args)
        dec = lambda p: [float(np.array([int(h, 16)], np.uint32).view(np.float32)[0]) for h in p[1:]]
        if args[0]:
            got = run(SREF, "mtf", data, dec(params[0]))
        else:
            got = np.concatenate([run(SREF, "mtf", data[[c]], dec(params[c])) for c in range(data.shape[0])])
        if siril_inf == params:
            sys.exit(f"{name}: Siril's own parameters for the Inf image equal the Inf-excluded ones; revisit the Inf rule")
        rec = {"name": name, "op": "autostretch", "args": args, "command": None,
               "note": note + "; truth: Siril's parameters with the Inf pixels set to 0, then Siril's MTF on the Inf image",
               "input": rec_in, "truth_input_ops": z, "params": params, "siril_inf_params": siril_inf,
               "out": arr32(got)}
        dump(os.path.join(w.sref, name + ".json"), rec)
        w.manifest.append(name)
        print(f"  {name}: Siril on the Inf image {[dec(p) for p in siril_inf]}, Inf excluded {[dec(p) for p in params]}")
    # The other operations on +-Inf input, straight from the oracle.
    for name, op, args, command, note, prog in [
        ("inf_asinh_rgb", "asinh", [50, 0, 0], "asinh 50", "+-Inf pixels", None),
        ("inf_asinh_rgb_off", "asinh", [20, 0.01, 1], "asinh -human 20 0.01", "+-Inf pixels, human, offset", None),
        ("inf_linstretch", "linstretch", [0.05], "linstretch -BP=0.05", "+-Inf pixels", None),
        ("inf_mtf", "mtf", [0.01, 0.2, 0.9], None, "+-Inf pixels (function level)", None),
        ("inf_ght", "ght", [2.0, 0.5, 0.0, 0.2, 1.0], "ght -D=2 -B=0.5 -SP=0.2", "+-Inf pixels clamp to [0, 1]", None),
        ("inf_clahe_rgb", "clahe", [3.0, 8], "clahe 3 8", "+-Inf pixels", SREF_CLAHE),
        ("inf_clahe_mono", "clahe", [3.0, 8], "clahe 3 8", "+-Inf pixels, mono", SREF_CLAHE),
    ]:
        ops = infops + ([["channels", [0]]] if name.endswith("mono") else [])
        w.case(name, R("synth.f32", ops), op, args, command, note, prog)
        if name == "inf_asinh_rgb":
            w.case("inf_asinh_mono", R("synth.f32", infops + [["channels", [2]]]), op, args, command, "+-Inf pixels, mono")


def round3_cases(w, ins):
    """Inputs the round-2 goldens missed: B exactly at the threshold
    float32(1e-3), HP = 0 (siril-cli accepts it with SP = LP = 0) and
    linstretch's BP = 1 edge."""
    R = w.inputs.ref
    # process_ght_args keeps B unless fabs(B) < (double)1.e-3f; the
    # threshold itself, 0.0010000000474974513 (what a float32 0.001 widens
    # to), is not below it, so B stays (power branch, not exponential).
    thr = float(F32(1e-3))
    if abs(thr) < float(F32(1e-3)) or not abs(np.nextafter(thr, 0)) < float(F32(1e-3)):
        sys.exit("threshold arithmetic")
    for name, B, command in [
        ("ght_b_f32thr", thr, f"ght -D=3 -B={thr!r} -SP=0.1"),
        ("ght_b_negf32thr", -thr, None),
    ]:
        w.case(name, R("synth.f32"), "ght", [3.0, B, 0.0, 0.1, 1.0], command,
               f"B={B!r} = (double)1e-3f: not below the threshold, B kept")
    # HP = 0: valid in siril-cli only with LP = SP = 0 (SP <= HP <= 1). The
    # Go zero value of HP means Siril's default 1; GHTParams.HPZero asks for 0.
    for name, D, B, command in [
        ("ght_hp0", 1.0, 0.0, "ght -D=1 -B=0 -SP=0 -HP=0"),
        ("ght_hp0_bpos", 2.0, 2.0, "ght -D=2 -B=2 -SP=0 -HP=0"),
        ("ght_hp0_bneg", 2.0, -0.5, None),
        ("ght_hp0_bm1", 2.0, -1.0, None),
    ]:
        w.case(name, R("synth.f32"), "ght", [D, B, 0.0, 0.0, 0.0], command, "HP = 0 with LP = SP = 0")
    # linstretch -BP=1, the upper edge siril-cli accepts: (x - 1) / 0.
    one = [["put", -1, [0, 1, 2, 3], [f32hex(1.0), f32hex(1.0), f32hex(0.5), f32hex(-0.0)]]]
    w.case("linstretch_bp1", R("synth.f32", one), "linstretch", [1.0], "linstretch -BP=1", "BP = 1: (x - 1) / 0")
    # Tiny and subnormal pixels for asinh: asinhf(beta*x) of a tiny x is x
    # itself (log(x + sqrt(x*x + 1)) would give 0), and with beta < 1.17 the
    # smallest ones make 0/0 (asinh(beta*x) and x*asinh(beta) both round to
    # 0), so that pixel is NaN in Siril.
    tiny = [["put", -1, list(range(8)), [f32hex(v) for v in (1e-30, 1e-38, 1.4e-45, 1e-20, 1e-10, -1e-30, 3e-39, 1e-44)]]]
    for name, chans, args, command in [
        ("asinh_tiny_mono", [["channels", [0]]], [1.0, 0.0, 0], "asinh 1"),
        ("asinh_tiny_mono_b05", [["channels", [0]]], [0.5, 0.0, 0], None),
        ("asinh_tiny_rgb", [], [50.0, 0.0, 0], "asinh 50"),
        ("asinh_tiny_rgb_b05", [], [0.5, 0.0, 0], None),
    ]:
        got, _ = w.case(name, R("synth.f32", tiny + chans), "asinh", args, command, f"tiny and subnormal pixels, beta {args[0]}")
        if name.endswith("_b05") and not np.isnan(got[0].reshape(-1)[2]):
            sys.exit(f"{name}: expected the 0/0 NaN at the smallest subnormal")
    # CLAHE on images smaller than one tile: OpenCV pads to a multiple of
    # the grid with BORDER_REFLECT_101, and when the padding is wider than
    # the image the reflection runs past the start (negative index), and a
    # width of 1 reflects onto itself.
    for name, ops, args in [
        ("clahe_tiny_mono", [["channels", [0]], ["crop", 0, 3, 0, 5]], [2.0, 8]),
        ("clahe_tiny_rgb", [["crop", 0, 2, 0, 3]], [2.0, 4]),
        ("clahe_1col_mono", [["channels", [1]], ["crop", 0, 6, 0, 1]], [2.0, 4]),
    ]:
        w.case(name, R("synth.f32", ops), "clahe", args, f"clahe {args[0]:g} {args[1]}", "image smaller than the tile grid", SREF_CLAHE)


def round4_cases(w, ins):
    """CLAHE's clip-limit conversion, MTF and asinh parameters at the edges
    of what the Go functions accept, and GHT D = 0."""
    R = w.inputs.ref
    hx = f32hex
    # OpenCV clahe.cpp: clipLimit = static_cast<int>(clipLimit_ *
    # tileSizeTotal / histSize), then max(clipLimit, 1). On x86-64 the
    # conversion is cvttsd2si: a value >= 2^31 (or +Inf) gives INT_MIN, so
    # the limit becomes 1, the strongest. A 32x32 image with tile 1 has
    # tileSizeTotal 1024, so clip = 2^37 is exactly 2^31 and 2^37 - 64 is
    # 2^31 - 1 (limit 2147483647: no bin reaches it). A NaN clip fails
    # OpenCV's clipLimit_ > 0 test: no limit. siril-cli accepts all of
    # these (it refuses only clip <= 0).
    crop = [["crop", 0, 32, 0, 32]]
    mono = R("synth.f32", [["channels", [0]]] + crop)
    rgb = R("synth.f32", crop)
    if clahe_tile_total(32, 32, 1) != 1024:
        sys.exit("tile arithmetic")
    below, at = float(2 ** 37 - 64), float(2 ** 37)
    if below * 1024 / 65536 != 2 ** 31 - 1 or at * 1024 / 65536 != 2 ** 31:
        sys.exit("clip arithmetic")
    lim1 = run(SREF_CLAHE, "clahe", mono[1], [1.0, 1])     # int(1*1024/65536) = 0 -> 1
    nolim = run(SREF_CLAHE, "clahe", mono[1], [0.0, 1])    # clip 0: no limit
    if np.array_equal(lim1.view(np.uint32), nolim.view(np.uint32)):
        sys.exit("clahe clip conversion: limit 1 and no limit give the same output")
    for name, clip, want in [
        ("clahe_clipconv_below", below, nolim), ("clahe_clipconv_at", at, lim1),
        ("clahe_clipconv_1e12", 1e12, lim1), ("clahe_clipconv_1e30", 1e30, lim1),
        ("clahe_clipconv_inf", float("inf"), lim1), ("clahe_clipconv_nan", float("nan"), nolim),
    ]:
        got, _ = w.case(name, mono, "clahe", [clip, 1], f"clahe {clip!r} 1" if np.isfinite(clip) else f"clahe {clip} 1",
                        f"clip {clip!r} on 32x32, tile 1: clip*1024/65536 converted by cvttsd2si", SREF_CLAHE)
        if not np.array_equal(got.view(np.uint32), want.view(np.uint32)):
            sys.exit(f"{name}: not the expected {'limit 1' if want is lim1 else 'unlimited'} output")
    w.case("clahe_clipconv_rgb_at", rgb, "clahe", [at, 1], f"clahe {at!r} 1", "RGB, clip 2^37: limit 1", SREF_CLAHE)

    # MTF with finite parameters outside siril-cli's mtf ranges. Siril's
    # autostretch hands MTF such values (find_unlinked_midtones_balance:
    # midtones 0 when c0 >= median, 1 in the inverted branch when c1 <=
    # median, shadows above 1, highlights below 0), so the Go MTF takes
    # every finite shadows/highlights and 0 <= midtones <= 1.
    edge = [["put", -1, list(range(10)), [hx(v) for v in (0.0, -0.0, 1.0, 0.3, 0.4, 0.5, -0.2, 1.5, 0.29999998, 0.30000001)]]]
    for name, args, command in [
        ("mtf_m0", [0.1, 0.0, 1.0], None), ("mtf_m1", [0.1, 1.0, 1.0], None),
        ("mtf_lo_gt_hi", [0.5, 0.3, 0.4], "mtf 0.5 0.3 0.4"), ("mtf_lo_eq_hi", [0.3, 0.6, 0.3], "mtf 0.3 0.6 0.3"),
        ("mtf_lo_ge_1", [1.2, 0.3, 1.0], None), ("mtf_hi_le_0", [0.0, 0.3, -0.1], None),
        ("mtf_wide", [-0.5, 0.25, 1.5], None),
    ]:
        w.case(name, R("synth.f32", edge), "mtf", args, command, f"MTF{tuple(args)}: finite parameters at the edge of the Go range")

    # Asinh at the edge of the Go range: offset just below 1, a large beta
    # (beta*x overflows to +Inf for input above 1), beta in (0, 1) (Siril's
    # GUI range is 0..1000; siril-cli refuses beta < 1).
    big = [["put", -1, list(range(6)), [hx(v) for v in (0.0, -0.0, 1.0, 1.5, 3.0, 0.999)]]]
    for name, chans, args, command in [
        ("asinh_off099_mono", [["channels", [0]]], [10.0, 0.99, 0], "asinh 10 0.99"),
        ("asinh_off099_rgb", [], [10.0, 0.99, 0], "asinh 10 0.99"),
        ("asinh_bigbeta_mono", [["channels", [0]]], [3e38, 0.0, 0], "asinh 3e+38"),
        ("asinh_bigbeta_rgb", [], [3e38, 0.0, 0], "asinh 3e+38"),
        ("asinh_b03_rgb_human", [], [0.3, 0.0, 1], None),
    ]:
        w.case(name, R("synth.f32", big + chans), "asinh", args, command, f"asinh beta {args[0]!r}, offset {args[1]!r}")

    # GHT D = 0: siril-cli refuses it, Siril's GUI accepts it (the D
    # slider starts at 0) and GHT() returns the clamped input there.
    px = [["put", -1, list(range(8)), [hx(v) for v in (0.0, -0.0, -0.25, 1.0, 1.5, 0.1, 0.5, 0.9)]]]
    for name, args in [
        ("ght_d0_b0", [0.0, 0.0, 0.0, 0.3, 1.0]), ("ght_d0_bpos", [0.0, 3.0, 0.05, 0.2, 0.8]),
        ("ght_d0_bneg1", [0.0, -1.0, 0.0, 0.5, 1.0]), ("ght_d0_bneg", [0.0, -2.5, 0.1, 0.1, 0.9]),
    ]:
        w.case(name, R("synth.f32", [["nanpattern"]] + px), "ght", args, None, "D = 0 (Siril's GUI): the clamp alone")


GHT_COEF = ["b1", "a2", "b2", "c2", "d2", "e2", "a3", "b3", "c3", "d3", "e3", "a4", "b4"]


def ght_sweep(w):
    """GHTsetup coefficients and GHT pixels over random parameter sets, one
    file per B family (B = -1, B < 0, B = 0, B > 0): 256 sets each, all
    within the range the Go GHT accepts (0 <= D <= 10, -5 <= B <= 15,
    0 <= LP <= SP <= HP <= 1), including D = 10, D = 1e-9 and 4e-9 (where
    expm1(D) and exp(D) - 1 round to different float32), LP = 0, LP = SP,
    SP = HP and HP = 1. Coefficients come from the verbatim GHTsetup (sref ghtcoef);
    the first 64 sets of each family also carry pixels: 0, -0, a negative
    value, NaN, LP, SP and HP exactly and one float32 step either side,
    1, 1.5 and six random values, through the verbatim GHT (sref ght).

    Purpose: a coefficient rounding that differs from Siril's float/double
    promotion on some fraction p of parameter sets survives all 256 with
    probability (1-p)^256 (round-2 verifier: p = 4.3% to 15.8% for the
    rounding mutants found, so < 2e-5); the exact-LP/SP/HP pixels pin the
    `<` branch comparisons, and 0 / -0 the "0 stays 0" rule."""
    rng = np.random.RandomState(2026)
    out_dir = os.path.join(w.out, "sweep")
    os.makedirs(out_dir, exist_ok=True)
    one = np.zeros((1, 1, 1), np.float32)
    for fam in ["bneg1", "bneg", "b0", "bpos"]:
        cases = []
        for k in range(256):
            D = float(rng.choice([rng.uniform(0.01, 3.0), rng.uniform(3.0, 10.0)]))
            if k == 0:
                D = 10.0
            if k in (2, 3):
                # tiny D where float32(expm1(D)) != float32(exp(D) - 1) even
                # with exp(D) - 1 in double: pins expm1 (process_ght_args)
                D = [1e-9, 4e-9][k - 2]
                if F32(np.expm1(D)) == F32(np.exp(D) - 1.0):
                    sys.exit(f"D={D}: expm1 and exp-1 agree in float32")
            SP = float(rng.uniform(0, 1))
            r = rng.uniform()
            LP = 0.0 if r < 0.35 else SP if r < 0.45 else float(SP * rng.uniform())
            r = rng.uniform()
            HP = 1.0 if r < 0.3 else SP if r < 0.4 else float(SP + (1 - SP) * rng.uniform())
            B = {"bneg1": -1.0, "bneg": float(rng.uniform(-5, -0.002)), "b0": 0.0,
                 "bpos": float(rng.uniform(0.002, 15))}[fam]
            if fam == "bneg" and k == 1:
                B = -5.0
            if fam == "bpos" and k == 1:
                B = 15.0
            if fam == "b0" and k % 16 == 1:
                B = float(rng.uniform(-0.00099, 0.00099))   # |B| < 0.001 counts as 0
            args = [D, B, LP, SP, HP]
            case = {"args": args, "coef": "".join(run(SREF, "ghtcoef", one, args))}
            if k < 64:
                f = F32
                px = [0.0, -0.0, -0.25, np.nan]
                for v in (LP, SP, HP):
                    v = f(v)
                    px += [v, np.nextafter(v, f(2)), np.nextafter(v, f(-1))]
                px += [1.0, 1.5] + list(rng.uniform(0, 1, 6))
                data = np.array(px, np.float32).reshape(1, 1, -1)
                got = run(SREF, "ght", data, args)
                case["input"] = {"shape": list(data.shape), "sha256": sha32(data), "bits": hx32(canon32(data))}
                case["out"] = arr32(got)
                case["libm"] = libm_record(data, args, got)
            cases.append(case)
        dump(os.path.join(out_dir, f"ght_{fam}.json"), {
            "name": f"ght_{fam}", "op": "ght", "coef_fields": GHT_COEF,
            "note": "random GHT parameter sets; coef = verbatim GHTsetup fields as float32 hex, "
                    "out = verbatim GHT on the input pixels", "cases": cases}, gz=True)
        print(f"  sweep ght_{fam}: {len(cases)} parameter sets")


def autostretch_sweep(w):
    """Autostretch parameters (unlinked and linked) on 96 RGB images made by
    scaling and offsetting the three base inputs, with random shadows clip
    and target. Pins the per-channel statistics and the linked averaging
    (c0 / n and c1 / n in float32, mtf.c:299-325) on both branches with
    c0 > 0 and c1 < 1, where the clamps do not hide them."""
    rng = np.random.RandomState(77)
    R = w.inputs.ref
    cases, nrm, inv = [], 0, 0
    for k in range(96):
        src = ["synth.f32", "bright.f32", "mixed.f32"][k % 3]
        if src == "bright.f32":   # keep every channel median > 0.5 (inverted branch)
            mul, add = F32(rng.uniform(0.75, 1.15)), F32(rng.uniform(-0.05, 0.1))
        else:
            mul, add = F32(rng.uniform(0.3, 1.6)), F32(rng.uniform(-0.15, 0.35))
        clip = float(F32(rng.uniform(-4, 1.5)))
        target = float(F32(rng.uniform(0.02, 0.5)))
        rec_in, data = R(src, [["mul", f32hex(mul)], ["add", f32hex(add)]])
        pu = run(SREF, "asparams", data, [0, clip, target])
        pl = run(SREF, "asparams", data, [1, clip, target])
        s = np.array([int(pl[0][1], 16)], np.uint32).view(np.float32)[0]
        h = np.array([int(pl[0][3], 16)], np.uint32).view(np.float32)[0]
        nrm += s > 0
        inv += h < 1
        cases.append({"input": rec_in, "clip": clip, "target": target, "unlinked": pu, "linked": pl})
    if nrm < 16 or inv < 16:
        sys.exit(f"autostretch sweep: linked shadows > 0 in {nrm}, highlights < 1 in {inv} cases (want >= 16 each)")
    dump(os.path.join(w.out, "sweep", "autostretch.json"), {
        "name": "autostretch", "op": "asparams",
        "note": "Siril's unlinked and linked autostretch parameters [retval, s, m, h]",
        "cases": cases}, gz=True)
    print(f"  sweep autostretch: {len(cases)} images, linked shadows > 0 in {nrm}, highlights < 1 in {inv}")


def main():
    out = os.path.join(HERE, "golden")
    if "--out" in sys.argv:
        out = os.path.abspath(sys.argv[sys.argv.index("--out") + 1])
    w = Writer(out)
    ins = w.inputs
    synth = ins.add("synth.f32", make_synthetic())
    ins.add("bright.f32", make_bright())
    ins.add("mixed.f32", make_mixed())
    ins.add("clahe_clip_mono.f32", make_clahe_clip(1, 128, 96, 21))
    ins.add("clahe_clip_rgb.f32", make_clahe_clip(3, 128, 96, 22))
    ins.add("clahe_clip_rgb_pad.f32", make_clahe_clip(3, 90, 70, 23))

    siril_cli_cases(w, synth)
    basic_cases(w, ins)
    clahe_cases(w, ins)
    edge_cases(w, ins)
    round2_cases(w, ins)
    round3_cases(w, ins)
    round4_cases(w, ins)
    ght_sweep(w)
    autostretch_sweep(w)
    w.finish()


if __name__ == "__main__":
    try:
        main()
    finally:
        shutil.rmtree(SCRATCH, ignore_errors=True)
