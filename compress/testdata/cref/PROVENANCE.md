# compress/testdata/cref: cfitsio truth for the compression code

`cref.c` is linked against the system libcfitsio. Everything in `golden/`
is output of the real cfitsio library run by this harness; nothing in it
was produced by the Go code under test. Inputs are generated inside
`cref.c` by a fixed-seed LCG, so a regeneration is deterministic.

The inputs are **not committed**. Each case's JSON records the generator
(`fill`, or `kind` + `gen_seed`) and `in_sha256`, the sha256 of the exact
input bytes cfitsio was given; `internal/crefgen` is a line-for-line Go
port of the generators, and every test regenerates its input and fails
unless the hash matches (no tolerance: one differing bit fails). Large
cfitsio outputs (quantized integers, the dither table, decodes of lossy
images) are stored as exact digests in the format of
`testdata/cfref/cfref.c` (`internal/testutil.Values`): the element count,
sha256 over a canonical big-endian encoding of every element (flag byte +
8-byte value; NaN canonical), every null as a run, and explicit values at
the first/last 8 indices, the plane corners, the first min/max and 64
strided indices. A hash mismatch is a failure; the sample says where.
Arrays of at most 64 values stay explicit. Codec streams (`.rice`, `.hc`,
`.plio`) are committed as they are. The cfitsio-written FITS files
(`<case>.fits`, `<case>.src.fits`) are committed gzipped as
`<case>.fits.gz` / `<case>.src.fits.gz` (`gzip -n -9`: no file name, no
timestamp); the tests inflate them before reading (`goldenFITS` in
`write_compressed_cfitsio_test.go`), and `make verify` compares the
inflated bytes with a fresh cfitsio run, so the check is on the exact bytes
cfitsio wrote whatever gzip made the `.gz`. The digests cut the directory
from 2,502,607 to about 1.2 MB; round 2 added the img64 and imgw goldens
(about 1.7 MB); round 3 gzipped the 124 FITS files (1,218,240 bytes ->
278,673), taking the tracked `compress/` tree from 1,746,053 to 940,619
bytes; the largest file is `cref.c` (about 93 KB).

## Tool versions used for the committed goldens

| tool       | version |
|------------|---------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`, `/usr/lib/libcfitsio.so.4.7.0`, sha256 `946cc8aef45d158d559ea722da95e8242ed7a9007ab9770911960dd02bbab65d`) |
| zlib       | 1.3.2 (Arch `zlib 1:1.3.2-3`; only affects the deflate bytes inside the GZIP goldens, which the Go tests compare after inflating) |
| gzip       | 1.15 (Arch `gzip`; makes the committed `.fits.gz` files; `make verify` inflates them, so another gzip version cannot fail the check) |
| compiler   | gcc 16.2.1 (Arch `gcc 16.2.1+r23+gd564253eb6c8-1`), `-O2 -Wall -Wextra -Wno-missing-field-initializers -std=c99` |
| host       | Linux x86_64 (LP64) |

Live interop tests that use this harness (`cref decode`, `cref
decodenull`) require cfitsio 4.7.0 and check the version it reports; they
build their own copy of the harness from `cref.c` and the `Makefile`. The interop matrix in the
root package also runs `fitsverify` and `funpack` (both from the same
Arch `cfitsio 1:4.7.0-1` package; funpack reports "1.7.0 (Dec 2013)
CFITSIO version 4.070") and astropy 8.0.2 / numpy 2.4.6 (repository venv,
pyenv Python 3.11.11).

## Commands

```sh
make -C compress/testdata/cref            # build ./cref (gitignored)
make -C compress/testdata/cref goldens    # regenerate every golden in place
                                          # (then gzip -n -9 every .fits)
make -C compress/testdata/cref verify     # regenerate into verify.out/ and
                                          # diff -r against golden/ (its
                                          # .fits.gz inflated in verify.ref/)
```

`make clean` removes only the binary, `verify.out/` and `verify.ref/`; no target deletes
committed goldens. `make verify` reproduced every committed golden
byte-for-byte with the versions above.

## What is in golden/

| files | produced by | used by |
|-------|-------------|---------|
| `<case>.out.json` (36 cases, incl. `runs3_*`: runs of exactly 3 equal pixels, which make FnNoise5 take the 2nd-order median over `nvals` entries, stale ones included; `runs3_stale_200x2_q4_dither1` makes the result depend on the exact order quick_select_float leaves. Round 2: `range_fsub_256x1` and `center_fsum_64x1` pin that cfitsio computes `maxval - minval` and `minval + maxval` in float, not double; `margin_64x1` (4294967292.8 steps: inside the `2*2147483647 - N_RESERVED_VALUES` margin, not quantized; round 3 `margin_far_64x1`, 4294967286.75 steps, whose centered extremes, about -+2147483643, are clear of the sentinels, so the margin stays pinned now that the Go quantizer also refuses a pixel on a sentinel) and `halfmargin_64x1` (2147483638.4 steps: inside the `2147483647 - N_RESERVED_VALUES` margin, centered) reach both margins with a positive qlevel, found by searching qlevel and one extreme pixel with `fits_img_stats_float` (with an absolute qlevel the ratio of two floats can never fall inside either margin); `wide_center_64x1` the centered zero point; `gauss_64x64_q0_dither1` qlevel 0; `tiny_2x8`, `tiny_4x5_abs`, `tiny_2x2_abs` the FnNoise1 nx < 3 and FnNoise3 nx < 5 rules. Round 2, literal inputs that each pin one branch: `n5_one_9x3` / `n5_two_10x3` / `n5_zero_nulls_9x3` (FnNoise5 rows with nvals 1 and nvals2 0, nvals 2 and nvals2 1, nvals 0), `n1_short_nulls_5x3` (a FnNoise1 row with one difference), `clip_tie_27x1` (the 5-sigma clip on an exact tie: |26 - 1| == 5 * 5), `one_good_3x3` (FnMeanSigma with one good pixel), `abs_nulls_10x4[_dither1]` (absolute qlevel with nulls; the minimum is among FnNoise3's first 4 pixels of a row), `all_null_*` / `all_nan_10x10_q16_dither2` (every pixel null: the ngood == 0 dummy minval/maxval/stdev, qlevel 4/16/-0.5, no dither / dither 1 / dither 2)) | `cref gen`: `fits_img_stats_float`, `fits_quantize_float` | `compress/cref_test.go` |
| `fits_rand_value.json` | `cref gen`: `fits_init_randoms` -> `fits_rand_value[0..9999]` (digest) | `compress/randoms_cref_test.go` |
| `rice/` (24 cases) | `cref encode`: `fits_rcomp`, `fits_rcomp_short`, `fits_rcomp_byte` on full-range-jump, extreme, mixed-block, tail and block-16 inputs | `compress/encode_cref_test.go` |
| `hcomp/` (20 cases) | `cref encode`: `fits_hcompress` (BITPIX 8/16) and `fits_hcompress64` (BITPIX 32), lossless and scaled, odd/minimum sizes, full-range int32; round 2: 8x10, 10x8, 16x18, 18x16, 8x9, 9x8 and i32 16x18 (a quadrant side that is a power of two next to a larger one: `qtree_encode`'s nqmax / log2n) | `compress/encode_cref_test.go` |
| `plio/` (9 cases) | `cref encode`: `pl_p2li` (long runs, values up to 2^24-1, isolated pixels, a 36003-word list that sets header word 5, and `plio_dv_boundary`: value steps of exactly +-4095/4096/4097, the 12-bit/24-bit set-high boundary) | `compress/encode_cref_test.go` |
| `img/` (96 cases: `<case>.json` with cfitsio's decode digest for lossy cases, `<case>.fits.gz` unless cfitsio refused) | `cref encode`: whole images compressed by `fits_img_compress` (the fpack path) for u8/s8/i16/u16/i32/u32/f32 x RICE_1/GZIP_1/GZIP_2/HCOMPRESS_1/PLIO_1 x tilings x NO_DITHER/SUBTRACTIVE_DITHER_1/2 x ZDITHER0 1/2/3/4/5/7/9000/10000, NaNs, constant and unquantizable tiles, 1-D and 3-D (incl. HCOMPRESS/PLIO cubes), HCOMPRESS default tile heights 30/16/24/20/28 (NAXIS2 30, 31, 33, 49, 241), in-tile dither-table wrap (`rice_f32_wrapin_*`, ZDITHER0 10000 on one 10000-pixel tile), caller keywords (`rice_i16_hdr`, whose cfitsio-written source image is kept as `rice_i16_hdr.src.fits`), float IEEE specials (`rice_f32_inf_*`, `gzip2_f32_inf_fallback_t40x5`: +-Inf, NaN, denormals, -0, +-FLT_MIN, also in a GZIP_COMPRESSED_DATA fallback tile), HCOMPRESS tile limits (`hcomp_i16_t3x8` refused, `hcomp_i16_11x8_t8x8` enlarged by cfitsio, `hcomp_i16_12x24_t4x12` and `hcomp_i16_12x9_t8x5` accepted at 4 pixels) and tiles longer than the axis (`rice_i16_11x8_t100x3`, `gzip1_f32_11x8_t100x3`: ZTILE1 = 100 kept), more than 3 axes (`rice_i16_4d_t10x6x2x1`, `gzip2_i16_6d`, `rice_i16_7d` refused). Round 2: a NaN and a real -1e30 in one tile (`rice_f32_nan_realnull_t40x5`) and in most whole-row tiles (`rice_f32_nan_realnull_rows`), a real FLOATNULLVALUE beside a NaN (`rice_f32_floatnullvalue_t40x5`: cfitsio nulls it), the quantizer range margins through the writer (`rice_f32_margin_64x2`, `rice_f32_center_64`, `rice_f32_halfmargin_64`), the dither re-seed in a tile with nulls (`rice_f32_wrap_nan_z9000`) and a float source with BLANK = 5 (`rice_f32_blank_hdr`, source kept as `.src.fits`) | `write_compressed_cfitsio_test.go`, `write_compressed_astropy_test.go`, `write_compressed_nulls_test.go`, `compressed_cref_read_test.go` |
| `imgw/` (5 cases, same JSON as `img/`, input seed 7000 + index) | `cref encode`: images compressed while they are written (`fits_write_imgnull` into a compressing HDU, null value -1e30, the generator's NaNs passed as -1e30) rather than by `fits_img_compress`. Each 40 x 10 image has, in both 40 x 5 tiles, a NaN and a real FLOATNULLVALUE (-9.11912e-36), so cfitsio quantizes that pixel as a value (`fits_img_compress` would null it: it uses FLOATNULLVALUE as its own null): `rice_f32_fnv_quant_t40x5` (SUBTRACTIVE_DITHER_1), `gzip2_f32_fnv_quant_nodither` (NO_DITHER, qlevel 16). Round 3, a caller null value (`"nullvalue"` in the JSON; the NaNs passed as it, and real pixels equal to it nulled by cfitsio): `rice_f32_nullval_m999_t40x5` and `rice_f64_nullval_m999_t40x5` (-999, float32 and float64 via TDOUBLE), `rice_f32_nullval_denorm_t40x5` (the denormal 1.04488e-40: `fits_write_imgnull` compares the raw value, so those pixels are nulls, while another denormal, -3e-42, stays a value) | `write_compressed_cfitsio_test.go` (`TestCompressedWriterFloatNullValueQuantizedMatchesCfitsio`), `write_compressed_quantlimits_test.go` (`TestCompressedWriterNullValueMatchesCfitsio`) |
| `img64/` (20 cases, same JSON as `img/`, input seed 8000 + index; input hashed as float64 little-endian; `dec` is cfitsio's TDOUBLE decode, not rounded) | `cref encode`: BITPIX -64 images compressed by `fits_img_compress` (each tile quantized by `fits_quantize_double`): RICE_1/GZIP_1/GZIP_2, SUBTRACTIVE_DITHER_1/2, NO_DITHER, absolute and default qlevel, constant and huge rows (float64 GZIP_COMPRESSED_DATA fallback), IEEE specials (float64 denormals 1e-310/-4e-320, -0, +-DBL_MIN, +-Inf, NaN), NaN beside a real -1e30, a real DOUBLENULLVALUE beside a NaN (`rice_f64_dnv_t40x5`: cfitsio nulls it), a 1e8 background with 0.01 noise (`rice_f64_precise_1e8`), every pixel above FLT_MAX (`rice_f64_beyond_flt_max[_abs]`), a cube, ZDITHER0 9000 on a 100x100 tile; round 3: SUBTRACTIVE_DITHER_2 with NaNs and zeros in the same tiles (`rice_f64_sd2_nan_zeros_t16x8`) and the int64 zero-point overflow (`rice_f64_zpt_overflow_abs`, see below). The float64 generator's noise is a sum of 4 uniforms, not `lcg_gauss`: libm's log/cos and Go's differ in the last bit (float32 rounding hides that, float64 does not) | `write_compressed_cfitsio_test.go` (`TestCompressedWriterFloat64MatchesCfitsio`) |
| `hcomp_default_tiles.txt` | `cref encode`: cfitsio's default HCOMPRESS_1 ZTILE1/ZTILE2 (imcomp_init_table via `fits_create_img` on a compressing `mem://` file) for an 8 x NAXIS2 image, NAXIS2 = 4..2000, 1681, 8401, 13443, 166323 (the first heights giving 26, 22, 18, 17 rows) | `write_compressed_cfitsio_test.go` (`TestHcompressDefaultTileMatchesCfitsio`) |
| `img2/` (5 cases, same JSON as `img/`, input seed 9000 + index; kept out of `img/` so the set the reader tests pin is unchanged) | `cref encode`, round 3: SUBTRACTIVE_DITHER_2 with NaNs and zeros in the same tiles (`rice_f32_sd2_nan_zeros_t16x8`), a source image whose SIMPLE card is patched to F (`rice_i16_simple_f`, source kept as `.src.fits.gz`: cfitsio writes ZSIMPLE = F), and three cfitsio quantizer defects (below): `rice_f32_zpt_overflow_abs`, `rice_f32_center_overflow`, `rice_f32_null_sentinel_hit` | `write_compressed_cfitsio_test.go`, `write_compressed_quantlimits_test.go` |
| `MANIFEST` | `cref encode`: one line per encode case (kinds rice, hcomp, plio, img, imgw, img64, img2) | both test files fail when a listed golden is missing or the counts change |

`img/rice_i16_hdr.src.fits` carries `CHECKSUM = 'AAAAAAAAAAAAAAAA'` and
`DATASUM = '1234567890'` as fixed strings, not real checksums
(`fits_write_chksum` stamps the current date into the comment, which would
make the golden unreproducible); cfitsio copies them opaquely as
ZHECKSUM/ZDATASUM, which is what the test checks. fitsverify will
therefore report checksum warnings on that one source file.

Img cases are appended, never inserted: each case's input seed is
5000 + its index, so inserting a case would change every later golden.

Five `img` cases record a cfitsio refusal instead of a file: with
`"status": 413` `plio_i16_negative`, `plio_i32_adversarial`, `plio_u32`
("data out of range for PLIO compression (0 - 2**24)") and
`hcomp_i16_t3x8` ("Hcompress minimum tile dimension is 4 pixels"); with
`"status": 212` `rice_i16_7d` (7 axes; cfitsio compresses at most
MAX_COMPRESS_DIM = 6). The Go writer must refuse them too.

cfitsio reads every float tile it compresses with a null value
(`imcomp_compress_image` -> `ffgsve` -> `fffr4r4`), so its `fnan` test
applies: an all-ones exponent (NaN, +Inf, -Inf) is a null, an all-zero
exponent (+-0, denormals) becomes +0 -- in quantized tiles and in the
lossless GZIP_COMPRESSED_DATA fallback alike (`gzip2_f32_inf_fallback_t40x5`:
1e-40 decodes as 0). The Go writer does the same.

## What cfitsio itself gets wrong (observed, recorded here, not hidden)

- **HCOMPRESS output buffer overflow.** `imcomp_calc_max_elem` sizes the
  HCOMPRESS buffer at 2.2 (8/16-bit) or 4.4 (32-bit) bytes/pixel + 26 and
  `fits_hcompress*` does not check it; full-range random tiles need more
  (`hc_i16_jumps_16x16_s0`: 593 bytes vs 589, `hc_i16_min_4x4_s0`: 95 vs
  61, `hc_i32_jumps_33x17_s0`: 2517 vs 2494). Compressing such an image
  through `fits_img_compress` corrupted the heap ("corrupted size vs.
  prev_size"). The `img` HCOMPRESS cases therefore use noise-like content;
  full-range data is covered by `hcomp/`, which calls the codec directly
  with a large buffer.
- **Reading HCOMPRESS images with BZERO as TLONG crashes.** On cfitsio's
  own `img/hcomp_u16.fits`, `hcomp_s8.fits`, `hcomp_u32.fits`,
  `hcomp_u16_noise_s2.fits`, `fits_read_imgnull(..., TLONG, ...)`
  segfaults inside `imcomp_decompress_tile` (so `testdata/cfref dump`
  crashes on them). TDOUBLE reads and funpack work; `cref decode` reads
  TDOUBLE only.
- **funpack mis-decodes float images with nulls and tiles narrower than
  NAXIS1** (`img/rice_f32_nan_z1_t16x8.fits`: 888 of 1000 pixels differ
  from `fits_read_img` and astropy, which agree). Pinned by
  `TestFunpackMisdecodesNullTiledFloats`.
- **PLIO_1 uint16 is not interoperable.** cfitsio stores the raw unsigned
  value and its reader corrects BZERO for PLIO only; astropy 8.0.2 reads
  every pixel of cfitsio's own `img/plio_u16.fits` wrong (+32768, wrapping
  above 32767). Pinned by `TestPLIOUint16AstropyMisreadsCfitsio`; the Go
  writer refuses PLIO_1 for uint16.
- **A 4x4 HCOMPRESS tile through `fits_img_compress` fails** (status 413,
  heap damage from the undersized buffer above), so no `img` case uses a
  4x4 tile; `hcomp/hc_i16_min_4x4_s0` covers 4x4 with a large buffer and
  `img/hcomp_i16_12x24_t4x12` covers a 4-pixel tile dimension.
- **FnNoise5 2nd-order median.** `quantize.c` takes
  `quick_select_float(differences2, nvals)` where `nvals2` entries were
  written, so the median includes zeros or entries a previous row's
  `quick_select` left behind. The Go port reproduces it (and ports
  `quick_select_float` verbatim, since the permutation it leaves matters);
  goldens `runs3_*`, `img/rice_f32_inf_sd2_z3_t16x8`.
- **`fits_read_imgnull` keeps stale null flags across fallback tiles.**
  On a float image whose GZIP_COMPRESSED_DATA tiles hold NaNs, each such
  tile also gets the null flags of the earlier fallback tiles at the same
  offsets: on cfitsio's own `img/rice_f32_nan_realnull_rows.fits` 327 real
  pixels are flagged null besides the 26 NaNs (`cref decodenull`,
  `TestCfitsioImgnullStaleFallbackFlags`). `fits_read_img` with a null
  value is right, and it is what the `dec` digests and `cref decode` use.
- **A caller BLANK on a float image is copied.** `fits_img_compress`
  copies BLANK into the ZBITPIX -32 header; with no NaN there is no ZBLANK,
  and cfitsio's own decode then nulls the real pixels whose quantized
  integer equals BLANK (`img/rice_f32_blank_hdr`: 64 of 1000; astropy 8.0.2
  the same pixels, with "BLANK only applicable to integer data"). The FITS
  standard does not allow BLANK with BITPIX < 0; the Go writer refuses it.
- **A real FLOATNULLVALUE pixel becomes a null.** In a tile that has NaNs
  cfitsio reads the tile with FLOATNULLVALUE (-9.11912e-36) as the null
  value, so a real pixel of that value is stored as a null
  (`img/rice_f32_floatnullvalue_t40x5`, pixel 10). The Go writer keeps it;
  funpack 1.7.0 also returns NaN for it on files the Go writer wrote, since
  it writes its output with FLOATNULLVALUE as the null.
- **The quantizer overflows without a check, status 0.** `fits_quantize_float`
  / `_double` convert `zeropt/delta + 0.5` to the int64 `iqfactor` and each
  pixel to int32 without a range check. A constant tile of 100 at qlevel
  -1e-20 (`img2/rice_f32_zpt_overflow_abs`) decodes as -0.0922; 1e10 at
  -1e-12 (`img64/rice_f64_zpt_overflow_abs`) as -9223372.002; a
  4294966272-step ramp from 1e10 at qlevel -1 takes the centered branch,
  whose float sum `minval + maxval` rounds 1024 low, so its maximum
  quantizes to 2^31 and wraps (`img2/rice_f32_center_overflow`: 14294966272
  decodes as 9999998976); and a centered minimum can land on NULL_VALUE,
  read back as a null once the image has ZBLANK
  (`img2/rice_f32_null_sentinel_hit`, found by searching ramp endpoints and
  qlevel). The Go quantizer refuses such tiles and the writer stores them
  losslessly (`TestQuantizerOverflowStoredLosslessly`).
- `fits_rcomp` computes `pdiff = nextpix - lastpix` in `int`, which
  overflows for full-range int32 jumps (undefined behaviour in C). The
  goldens record what the compiled 4.7.0 library does (two's-complement
  wrap), which its decoder undoes; the Go encoder reproduces it.
