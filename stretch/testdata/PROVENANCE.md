# stretch test data: provenance

Every expected value used by the `stretch` tests comes from one of the three
truth sources below. No golden is computed by this package's Go code, and a
missing golden fails the test (it is never skipped).

## Golden format

All goldens are compact JSON written by `goldenfmt.py` (mirrored by
`../golden_test.go`):

* An array (image, vector set) is `{"shape", "sha256", "bits"}` with every
  value when it has at most 1024 values, else `{"shape", "sha256", "idx",
  "bits"}`. `sha256` covers the canonical encoding of **every** element in
  row-major (channel-major for images) order: float32 as big-endian IEEE
  bits with every NaN written as `0x7fc00000` (so -0 and +0 differ, NaN
  payloads do not), uint16 as big-endian. `bits` holds the canonical words
  in hex for every element, or for the sorted sample `idx`: first and last
  8, the four corners of every image plane, the first minimum and maximum,
  the first 8 NaNs and 64 evenly strided elements. A test fails on any hash
  mismatch; the sample only names the pixels in the failure message.
* Inputs are committed raw little-endian files in `golden/inputs/`
  (`*.f32` float32, `*.u16` uint16, channel-major), referenced as
  `{"file", "file_shape", "ops", "shape", "sha256"}`. `ops` are
  deterministic transformations applied by both the generator and the Go
  test (`goldenfmt.derive` / `golden_test.go derive`: channel selection,
  fill a channel, the NaN pattern, crop, multiply/add a float32, put
  values at pixels), and `sha256` is the digest of the result, so the Go
  test proves it rebuilt the same input. Small inputs are inline; a
  constant image is `{"const": "<hex>"}`.
* Autostretch parameters are `[retval, shadows, midtones, highlights]`,
  the floats as hex bits (they may be NaN).
* A GHT golden also has `"libm": {"idx", "bits"}`: the pixels where the
  verbatim oracle and the same code built with double-rounded
  logf/powf/expf (`sref/sref_libm`, `sref/shim/libm_double.h`) disagree,
  with Siril's values. `sref_libm` only locates these pixels, it is never a
  source of expected values. Only there may the Go value differ, and only
  if a 1-ULP different logf/powf/expf result reproduces Siril's value
  (`../ght_test.go`); the sha256 then covers every pixel.

This replaced pretty decimal pixel lists (13.5 MB for `golden/`) with
digests (1.4 MB with the round-2 cases and sweeps, about half of it the raw input images; largest file 147 KB). The cross-check of the
conversion: the sha256 of every old pixel array (Siril, siril-cli and cv2
outputs and inputs) equals the new digest, and every old autostretch
parameter set equals the new one bit for bit.

## 1. siril-cli goldens — `golden/stretch_*.json` (except `stretch_clahe_cv2.json`)

* Originally written by `gen_golden.py` (commit c17580b, 2026-04-05)
  running a real `siril-cli` on the synthetic input (now
  `golden/inputs/synth.f32`) for the 18 `DefaultPresets` commands, as full
  pixel lists.
* **The siril-cli build used then was not recorded.** The package doc
  said "siril-cli 1.5.0", `gen_golden.py` writes `requires 1.3.6` (a Siril
  script directive giving the *minimum* version), and the local build at
  `~/git/dmarkham/Siril` is git master `1.4.0-rc1-235-g551111560`
  (`meson.build` version `1.5.0`), rebuilt on 2026-07-12, after the goldens,
  and no longer runnable (`libopencv_calib3d.so.413` and FFmpeg/libjxl
  sonames missing). So the original build cannot be named; what follows
  establishes that the goldens are what a known siril-cli produces.
* `siril_cli_digests.py` reads the original siril-cli pixel files from git
  (`git show c17580b:...`, blob ids recorded) and writes their canonical
  digests to `golden/siril_cli_sha256.json` (`--check` recomputes and
  compares; `TestSirilCLIDigestsFromHistory`). `gen_siril_goldens.py`
  writes `stretch_<preset>.json` from the verbatim Siril 1.4.0 oracle
  (section 2) and **refuses to write them unless every output digest
  equals the siril-cli digest**; `TestSirilCLIDigests` checks the
  committed files against `siril_cli_sha256.json` without any tool.
* **A real siril-cli 1.4.0** is built from the release tarball by
  `build_siril_cli.sh` (2026-11-01) under
  `/var/lib/aurbuild/pkgbuild/siril`, nothing installed system-wide:
  Siril 1.4.0 needs OpenCV 4 and the system has OpenCV 5.0.0, so the Arch
  package `opencv-4.13.0-9-x86_64.pkg.tar.zst` (sha256
  `69ce38e8db804ea58f187996e4771b9a6245c762cd22c047a0ee7cb93c107635`, from
  archive.archlinux.org) is unpacked into `deps/opencv4` and linked through
  a local pkg-config file; meson 1.12.1, ninja 1.13.2, gcc 16.2.1,
  glibc 2.44, buildtype release, ffmpeg/libjxl/libheif/criterion off,
  libcurl on (Siril 1.4.0 does not compile without it). The wrapper
  `siril-cli-1.4.0` sets `LD_LIBRARY_PATH` and prints `siril 1.4.0`.
  `gen_golden.py <siril-cli>` (and `TestSirilCLIReproduce`, which runs it
  with `$FITS_SIRIL_CLI` or that wrapper) checks, with that siril-cli:
  - the 18 preset outputs: all 18 digests equal `siril_cli_sha256.json`,
    i.e. siril-cli 1.4.0 reproduces the original goldens bit for bit;
  - every oracle golden in `golden/sref` that names a siril command
    (146): all come out of siril-cli with the golden's digest, CLAHE
    (OpenCV 4.13 in siril-cli, 5.0.0 in the oracle) included. siril-cli's
    FITS loader divides float data whose maximum exceeds 10 by that
    maximum (`io/image_format_fits.c:705`, `:794`), so the 7 `inf_*`
    inputs are written with the markers +7 / -7 at the ±Inf pixels (the
    script checks the input has no such value and nothing above 10) and
    the script runs `pm "iif($T==7, 1e39, iif($T==-7, -1e39, $T))"`
    before the command: pixel math evaluates in double and stores
    (float)1e39 = +Inf (`pixelMath/pixel_math_runner.c:804`); the script
    first saves the image after pm alone and requires it to equal the
    ±Inf input bit for bit;
  - `golden/siril_cli_args.json`, written by `gen_golden.py <siril-cli>
    --write-args`: 102 commands at the edges of linstretch's `-BP=`,
    autostretch's shadows clip and target background, every GHT range
    check, asinh's beta and offset, mtf's three parameters and clahe's
    clip and tile (NaN and ±Inf included; clahe clip 2^42 and 2^42 - 2048,
    where clip*32/65536 is 2^31 and 2^31 - 1 on the 8x4 image with tile 1,
    and the script requires 2^42, 1e30 and inf to give the output of clip
    2 (limit 1) and 2^42 - 2048, 1e5 and nan the unlimited one) on a small
    mono image (0, 1, -0, four equal pixels so that CLAHE's limit matters,
    and a spread); for each,
    whether siril-cli runs it (exit status 0 and an output file), Siril's
    error message if not, and the output pixels if so. Each accepted
    command is run twice: `ght -B=nan` reads GHTsetup coefficients it
    never set and its pixels differ run to run (`out_varies`). The table
    is re-written identically on every run. `TestSirilCLIArgs` requires
    the Go functions to refuse exactly the refused commands (except the
    divergences it lists) and to give siril-cli's pixels for the rest.
* `golden/inputs/synth.f32` is reproduced byte for byte by
  `make_synthetic` (numpy `RandomState(42)`, numpy 2.4.6); its digest equals
  the digest of the original `input.json`.

## 2. Verbatim Siril 1.4.0 oracle — `sref/`, goldens in `golden/sref/`

* Source: `siril-1.4.0.tar.bz2` from https://free-astro.org/download/
  (sha256 `439def7c40ad783afb82e87f3c656d85b449c701f67ab0a9b97ca372ea2d73c9`,
  the hash pinned by the Arch AUR PKGBUILD), unpacked at
  `/var/lib/aurbuild/pkgbuild/siril/src/siril-1.4.0`. License GPL-3.0-or-later.
* `sref/extract.sh` copies, without edits, `filters/mtf.c`, `filters/mtf.h`,
  `filters/ght.c`, `filters/ght.h`, `core/arithm.h`, `rt/rt_math.h`, and cuts
  whole functions by line from `filters/asinh.c` (asinhlut_ushort,
  asinhlut_float, asinhlut, command_asinh), `algos/quantize.c`
  (FnMeanSigma_float), `algos/statistics_float.c`
  (reassign_to_non_null_data_float, siril_stats_float_mad), `algos/sorting.c`
  (histogram_median_float), `rt/rt_algo.cc` (findMinMaxPercentile) and
  `opencv/opencv.cpp` (fits_to_bgrbgr_float, cvClahe_float). Each excerpt
  names its source file and line. `make -C sref verify-excerpt` re-extracts
  and diffs.
* Hand-written parts, kept minimal and documented in each file:
  `sref/shim/` (headers: Siril's macros/enums copied with line references,
  the `fits`/`imstats` structs reduced to the members used, logging as
  no-ops), `sref/glue.c` (calls the verbatim statistics functions in the order
  `statistics_internal_float` does), `sref/sref.c` and `sref/sref_clahe.cpp`
  (drivers mirroring `process_asinh`, `process_autostretch`,
  `process_ght_args`, `process_linstretch`, `process_clahe` from
  `core/command.c`; the GHT driver skips the CLI's range checks so B < 0 —
  allowed by the GUI, `gui/uifiles/histogram_dialog.ui:77`, range -5..15 —
  can be generated; such goldens say so in their `note`).
  `sref/shim/libm_double.h` builds `sref_libm` (see Golden format).
* Build: `make -C sref` (gcc/g++ 16.2.1, `-O2 -ffp-contract=off`, no
  `-march`, no OpenMP; glibc 2.44; glib2 2.88.3; system OpenCV 5.0.0, Arch
  `opencv 5.0.0-12`, baseline SSE/SSE2, so no FMA in the OpenCV code paths
  used). Binaries are gitignored.
* Goldens: `make -C sref goldens` runs `gen_siril_goldens.py` with
  `$FITS_ASTROPY_PYTHON` if set, else the repository's `venv/bin/python`,
  else `python3` (it needs numpy; override with `PY=<venv python>`). It
  writes the inputs (`golden/inputs/*.f32` except `cv2_*`), the 18 preset
  goldens (section 1), every `golden/sref/*.json` and
  `golden/sref/MANIFEST`. `make -C sref verify` (and
  `TestSirilGoldensReproduce`) regenerates into a scratch directory and
  byte-compares every file. Each golden records the oracle op, arguments,
  the equivalent siril command (or null when siril-cli cannot express it)
  and a note. Besides the per-operation cases they include:
  `autostretch_*zero*` (one channel of an RGB image all zero, unlinked and
  linked: Siril's no-good-pixel fallback 0/0.2/1 applied to every
  channel); `preset_<name>_{bright,mixed}` (every `DefaultPresets` command
  on `inputs/bright.f32`, every channel median > 0.5, and
  `inputs/mixed.f32`, one bright channel of three); and `clahe_clip_*`
  (images whose CLAHE tiles are large enough that OpenCV's integer clip
  limit `int(clip*tileTotal/65536)` is 3 rather than the clamped minimum 1;
  the generator aborts unless the oracle output changes at limits 2 and 4).
  Edge cases (each checked by the generator to exercise what it is for):
  `nan_*` (NaN at different pixels in each channel plus one pixel NaN in
  all channels, for every operation; and `nan_channel_g*` /
  `nan_allnan_mono`, a channel entirely NaN, where Siril's statistics
  fail); `ght_out_above_1*` (Siril's GHT output exceeds 1);
  `ght_small_d_*` (D where float32(expm1(D)) differs from
  float32(exp(D))-1); `ght_b_tiny_*` / `ght_b_above_thresh` (B on both
  sides of `|B| < 0.001`); `ght_bneg035_lp_hp` (B = -0.35, where (B-1)/B
  and 1-1/B round differently in float32); `ght_bneg1_big_lp`;
  `autostretch_median_*` (mono images with a plateau at 0.49 .. 0.5:
  Siril's histogram median puts 0.4999 and 0.5 on the inverted side and
  0.499 on the normal side); `asinh_beta0_*` (beta = 0, which siril-cli
  rejects; the oracle calls `command_asinh` directly) and
  `asinh_mono_above_1`.
* Round-2 cases (`round2_cases`): `neg_*` (synth - 0.08, 26.8% of pixels
  negative, plus -0 pixels; bright with 40 negative pixels for the
  inverted branch; every operation; `neg_linstretch_bp0` checks Siril
  returns +0 for -0); `autostretch_const_05*` (a median of exactly 0.5:
  normal branch; one float32 step above: inverted; checked against the
  branch the oracle took); `ght_b_0001` / `ght_b_neg0001` (B = ±0.001
  exactly, which siril-cli's double comparison sets to 0); `inf_*` (±Inf
  pixels). For `inf_autostretch*` the expected values are Siril's
  parameters for the same image with the Inf pixels set to 0
  (`truth_input_ops`) and Siril's MTF with those parameters applied to the
  Inf image, per channel for unlinked; Siril's own parameters for the Inf
  image are kept as `siril_inf_params` (NaN or meaningless: its
  findMinMaxPercentile casts an infinite histogram position to uint16,
  undefined behaviour) and the generator requires them to differ. The
  other `inf_*` cases are the oracle on the Inf image directly.
* Round-3 cases (`round3_cases`): `ght_b_f32thr` / `ght_b_negf32thr` (B =
  ±0.0010000000474974513 = (double)1e-3f, the threshold itself, which
  process_ght_args keeps: not below it); `ght_hp0*` (HP = 0 with LP = SP =
  0, accepted by siril-cli; the output is the clamp alone for B = 0, 2,
  -0.5 and -1); `linstretch_bp1` (BP = 1, every pixel +0, with 1, 0.5 and
  -0 pixels added); `asinh_tiny_*` (tiny and subnormal pixels, beta 1, 50
  and 0.5; at 0.5 the smallest subnormal is 0/0 = NaN in Siril, which the
  generator checks); `clahe_tiny_*` / `clahe_1col_mono` (images smaller
  than the tile grid, padded by more than their own size, and 1 pixel
  wide).
* Round-4 cases (`round4_cases`): `clahe_clipconv_*` (a 32x32 crop with
  tile 1, tileSizeTotal 1024: clip 2^37 - 64, 2^37, 1e12, 1e30, +Inf and
  NaN, and RGB at 2^37; OpenCV's `static_cast<int>(clip*1024/65536)` is
  cvttsd2si on x86-64, INT_MIN from 2^31 up, so the limit is 1; the
  generator requires each output to equal the oracle at clip 1 (limit 1)
  or clip 0 (no limit) as stated); `mtf_*` (midtones 0 and 1, shadows
  above highlights or equal, shadows >= 1, highlights <= 0, shadows < 0
  with highlights > 1: the finite values Siril's autostretch can hand
  MTF); `asinh_off099_*`, `asinh_bigbeta_*` (beta 3e38, input above 1),
  `asinh_b03_rgb_human`; `ght_d0_*` (D = 0, the GUI's clamp, in each B
  family, NaN/-0/negative/above-1 pixels). Arguments that are NaN or
  ±Inf are `null` in `args`, and `args_f64` lists every argument as
  big-endian float64 bits in hex.
* Sweeps, `golden/sweep/`: `ght_{bneg1,bneg,b0,bpos}.json`, 256 random GHT
  parameter sets per B family (`RandomState(2026)`, D in [0.01, 10] with
  D = 10 and D = 1e-9, 4e-9 (where float32(expm1(D)) != float32(exp(D) - 1)), B in [-5, 15] with B = -5 and 15, |B| < 0.001 in the B = 0
  family, LP = 0, LP = SP, SP = HP and HP = 1 among them): `coef` is the
  verbatim GHTsetup struct (`sref ghtcoef`, fields `b1 a2 b2 c2 d2 e2 a3 b3
  c3 d3 e3 a4 b4`), and the first 64 per family carry 20 pixels (0, -0,
  -0.25, NaN, LP, SP, HP exactly and one float32 step either side, 1, 1.5
  and 6 random) through the verbatim GHT with their libm list.
  `autostretch.json`: Siril's unlinked and linked parameters on 96 images
  (`RandomState(77)`: synth, bright or mixed times a float32 in [0.3, 1.6]
  plus one in [-0.15, 0.35], bright kept above 0.5), random clip in
  [-4, 1.5] and target in [0.02, 0.5]; the generator requires the linked
  shadows > 0 in at least 16 cases (48) and highlights < 1 in at least 16
  (23), where the c0/n and c1/n divisions are visible.

## 3. OpenCV (opencv-python) — `gen_clahe_cv2.py`

* Python 3.11.11 venv (pyenv), numpy 2.4.6, opencv-python-headless 5.0.0.93
  (`cv2.__version__` 5.0.0, baseline SSE/SSE2/SSE3).
* `golden/clahe_mono_u8.json`: cv2 8-bit CLAHE (clip 3, 8×8) on a
  deterministic image (`RandomState(3)`).
* `golden/clahe_mono_u16.json`: cv2 16-bit CLAHE, 6 cases (padded tiles,
  non-square grids, clip 0 = off, large clip).
* `golden/clahe_clip_cv2.json`: cv2 8- and 16-bit CLAHE, 6 cases whose
  integer clip limit is 3 to 11 (the cases above all clamp it to 1, so they
  cannot tell clip limits apart); the script aborts unless cv2 at limit - 1
  and limit + 1 gives a different output.
* `golden/lab_cv2.json`: cv2 float `COLOR_BGR2Lab` on all 33³ nodes of
  OpenCV's internal RGB→Lab LUT plus random/out-of-range vectors, and
  `COLOR_Lab2BGR` on random Lab vectors, in rows of width 1, 13 and 64
  (OpenCV's SIMD block of 8 pixels and its scalar tail round differently).
* `golden/nan_cv2.json`: cv2 float `COLOR_BGR2Lab` and `COLOR_Lab2BGR` on
  vectors with NaN in every combination of components (rows of width 1, 13
  and 64), and `cv2.add(x, 0, dtype=CV_32S)` (saturate_cast<int>(float),
  i.e. cvRound) on NaN, ±inf, out-of-range and half-way values.
* `golden/stretch_clahe_cv2.json`: Siril's float CLAHE pipeline rebuilt from
  cv2 calls on `golden/inputs/synth.f32` (clahe 3 8). The script asserts it
  equals the verbatim-Siril golden `golden/sref/clahe_rgb_3_8.json` bit for
  bit, so the pip wheel and the system OpenCV agree.
* `golden/jpeg_dqt_cv2.json`: the quantization tables (DQT segments, id,
  precision, 64 values in zigzag order) that `cv2.imencode('.jpg',
  IMWRITE_JPEG_QUALITY=q)` writes for q in 1, 5, 10, 25, 49, 50, 51, 75, 90,
  95, 100 (cv2 bundles libjpeg-turbo 3.1.2); the script requires every
  quality to give different tables.
* Inputs: `golden/inputs/cv2_*`. Re-running (`--out DIR`) reproduces every
  file byte for byte (`TestCLAHECv2GoldensReproduce`).

### Replaced goldens

`stretch_clahe_cv2.json` and `clahe_mono_u8.json` existed before (commit
26912ff) with no generator and no test loading them. The old
`stretch_clahe_cv2.json` held multiples of 1/255, i.e. an 8-bit pipeline, which
is not what Siril runs on float images (its float path quantizes L* to 16
bits), and no cv2 5.0.0 pipeline tried reproduced it (closest: uint8
truncation → BGR2Lab → CLAHE → Lab2BGR, 4699 of 9216 values different, up
to 16/255). The old `clahe_mono_u8.json` output was reproducible by cv2 from
its stored input, but that input had no recorded origin. Both were replaced
by the generator output above.

## glibc float libm (GHT tolerance)

`make -C sref libm-probe` (`sref/libm_probe.c`): of 2,000,000 deterministic
arguments, glibc 2.44 `powf`, `expf`, `logf` differ from the correctly rounded
float32 for 1290, 1253 and 18910 of them; `asinhf` for none. Siril's GHT calls
`logf`/`powf`/`expf`; the committed GHT goldens have 0 to 58 libm-sensitive
pixels each (their `libm` lists), and every other operation is compared
bit for bit.

## Storage of the sweep goldens

`golden/sweep/*.json` are written gzipped as `golden/sweep/<name>.json.gz`
(`goldenfmt.dump(..., gz=True)`: Python `gzip.compress`, level 9, mtime 0).
The Go tests read them through `testutil.ReadGolden`, and
`TestSirilGoldensReproduce` compares the regenerated JSON with the
decompressed committed bytes, byte for byte.
