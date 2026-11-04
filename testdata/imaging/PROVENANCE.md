# testdata/imaging provenance

Samples and truth for `fits.ReadFloat32`, `ReadFrame`, `ReadFrameRGB`
(readfloat32_proof_test.go, frame_test.go) and `AsImage` / `image.Decode`
(stdimage_test.go). Nothing here is produced by the Go library.

## Tools

| tool | version |
|------|---------|
| Python | 3.11.11 (pyenv), shared venv `fits-proof-audit/venv` |
| numpy / astropy | 2.4.6 / 8.0.2 (`gen_expected.py` refuses other versions) |
| libcfitsio, fpack, fitsverify | 4.7.0 (Arch `cfitsio 1:4.7.0-1`), fpack 1.7.0 |
| Siril source | git 55111156048b (1.4.0-rc1-235-g551111560, 2026-01-10), `~/git/dmarkham/Siril` |
| gcc | 16.2.1 |

## Regenerate (all byte-for-byte reproducible)

```sh
cd testdata/imaging
<venv>/bin/python gen_samples.py            # every *.fits except the camera crop
<venv>/bin/python crop_camera.py <Light_001.fits>   # camera crop (source not committed)
./make_compressed.sh                        # comp_*.fits.fz via fpack -C + zblank_variants.py
make -C cfnull goldens                      # cfnull/golden/*.json (cfitsio's reading of every comp_*.fz)
make -C cfnull verify                       # regenerate + cmp
<venv>/bin/python gen_expected.py           # readfloat32_expected.json, stdimage_expected.json
make -C sirilref goldens                    # sirilref/golden/*.json
make -C sirilref verify                     # regenerate + cmp
make -C sirilref verbatim                   # VERBATIM functions == Siril source
```

`TestImagingSamplesReproduce` deletes every generated sample in a scratch
copy (everything but the camera crop) and re-runs `gen_samples.py` and
`make_compressed.sh`: the file set and every byte must come back.

All sample files pass `fitsverify` except `i16_bscale0.fits` and
`f32_bscale0_2d.fits` (BSCALE = 0, 1 warning each) and
`i16_bscale_str_2d.fits` / `i16_bscale_bad_2d.fits` /
`i16_bscale_qnum_2d.fits` (BSCALE not a number) and
`f32_bzero_qnum_2d.fits` / `f32_bzero_hex_2d.fits` (BZERO not a number),
1 error each — all deliberately invalid.

## Samples

`gen_samples.py` writes the files card by card and byte by byte (astropy
would rescale): BITPIX 8 (plain, all 256 values; BLANK; BZERO=-128;
BSCALE=2/BZERO=10), BITPIX 16 (all 65536 values with BZERO=0 and with
BZERO=32768; no BZERO card; BSCALE=-1; BSCALE=0.5/BZERO=-3; BLANK with
BZERO=32768; BSCALE=0), BITPIX 32 (extremes + seeded random; BZERO=2^31;
BLANK + BSCALE=2.5/BZERO=1e6; BSCALE=-3), BITPIX 64 (extremes + random;
BZERO=2^63 written as a real-valued card, and as an integer literal in
`u64_intcard.fits`; BLANK), BITPIX -32 (unit values
with NaN and denormal; wide values with ±Inf; BSCALE=2/BZERO=1; physical
overflow), BITPIX -64 (underflow, NaN, rounding; wide; overflow), a planar
RGB cube and an orientation fixture (pixel (x,y) holds 1000·(4y+x)).

Added by the second proof pass (appended, so earlier files are
byte-identical): `i64_round.fits` / `i64_round_neg.fits` (BSCALE=-1
mirror) hold BITPIX 64 values where rounding u/(2^64-1) through a double
gives a different float32 than the exact rounding (found by a
numpy/fractions search near float32 midpoints; includes raw 2^39+1);
`f32_bscale_inexact_2d.fits` (BITPIX -32, BSCALE 0.1: physical values
float32 cannot hold exactly) and `f32_bscale0_2d.fits` (BSCALE 0).

Added by the third proof pass (appended): `i32_round.fits` /
`i32_round_neg.fits` (BSCALE=-1 mirror) hold the only two BITPIX 32
values, u = raw + 2^31 = 4294966655 and 4294967167, where
float32(float64(u)/(2^32-1)) is not the float32 nearest u/(2^32-1), plus
the two extremes. `find_i32_round.py` finds them by an exhaustive search
of all 2^32 values (numpy filters quotients within 4 double ULP of a
float32 midpoint; fractions.Fraction decides each candidate exactly).

Added by the fourth proof pass (appended): `u64_intcard.fits`, the
`u64.fits` pixels with `BZERO = 9223372036854775808` written as an integer
literal, the form cfitsio and astropy write for uint64 (decoded as
`TypeBigInt`; its expected values equal `u64.fits`'s; Siril: unsupported,
BITPIX 64).

Added by the round-2 proof pass (appended):
`f32_bzero_2d.fits` / `f64_bzero_2d.fits` (float BITPIX with BZERO and no
BSCALE card: the offset must be applied; every physical value exact in
float32) and `f32_bzero_inexact_2d.fits` (BZERO 10 on values whose sum is
not a float32, so AsImage must refuse it); `i16_bscale_str_2d.fits`
(BSCALE = 'abc'), `i16_bscale_bad_2d.fits` (BSCALE = 2.0.0) and
`f32_bzero_huge_2d.fits` (BZERO = 1E400, a valid FITS real that overflows
a double): scaling cards that are present but not finite numbers, which
cfitsio 4.7.0 `fits_read_img` silently replaces by the default (status 0)
and ReadFloat32 / AsImage refuse (gen_expected.py classifies them from
the raw cards: rule "error", kind "reject"; astropy cannot open the
2.0.0 file); `frame_naxis1_primary.fits` (a 1-D primary in front of an
IMAGE extension holding the `u16_orient.fits` pixels, for ReadFrame's
"primary if NAXIS >= 2, else HDU 1" rule).

Added by the round-2 second pass (appended): `f32_bzero_qnum_2d.fits`
(BZERO = '10.0') and `i16_bscale_qnum_2d.fits` (BSCALE = '0.5'): scaling
cards whose value is a quoted numeric string, i.e. a character string
(FITS 4.0 §4.2.1), not a number; `f32_bzero_hex_2d.fits` (BZERO = 0x1p3,
a C99 hex float, not a FITS integer or real per §4.2.3/§4.2.4). fitsverify
4.7.0 reports 1 error for each; cfitsio 4.7.0 `fits_read_img` (TFLOAT)
returns status 0 with BZERO 0 for '10.0' (0.5 1.5 -2.25 nan 100 -0),
BSCALE 1 for '0.5' (32768 32769 32868 ...), and strtod's 8 for 0x1p3
(8.5 9.5 5.75 nan 108 8). gen_expected.py classifies all three as rule
"error" / kind "reject". `i16_dexp_2d.fits` (BSCALE = 5.0D-1,
BZERO = -3.0D0) and `f32_bzero_dexp_2d.fits` (BZERO = 1.0D1): valid
D-exponent reals (fitsverify OK); astropy reads 0.5 / -3 / 10 and so
does cfitsio (fits_read_img: -16387 -16386.5 -503 -3.5 -3 -2.5 and
10.5 11.5 7.75 nan 110 10).

Added by round-2 proof2 (appended): `card_*.fits`, BSCALE / BZERO /
BLANK card variants (gen_samples.py, cards written verbatim) for the card
rule of ReadFloat32 / AsImage (`TestReadFloat32_ScalingCardVariants`).
`make -C cfnull goldens` writes cfitsio 4.7.0's reading of each
(fits_read_pixnull TDOUBLE) to `cfnull/cards/`; gen_expected.py records it
and astropy 8.0.2's scaled `getdata` under `refs`, and asserts both equal
the documented physical values for every accepted sample. Readings
(cfitsio / astropy), on the stored floats 0.5 1.5 -2.25 NaN 100 -0:

| sample | card(s) | cfitsio | astropy | rule |
|---|---|---|---|---|
| `card_f32_bzero_dup_2d` | BZERO 10.0, BZERO 20.0 | +20 | +20 | last card: +20 |
| `card_f32_bzero_dupbadfirst_2d` | 'abc', 20.0 | +20 | +20 | last card: +20 |
| `card_f32_bzero_dupbadlast_2d` | 20.0, 'abc' | +0 | fails | error |
| `card_f32_bzero_noeqsp_2d` | `BZERO   =10.0` | +10 | +0 | error |
| `card_f32_bzero_col10eq_2d` | `=` in column 10 | +10 | +0 | error |
| `card_f32_bzero_col11eq_2d` | `=` in column 11 | +10 | +0 | error |
| `card_f32_bzero_spacein_2d` | `1 0.0` | +1 | fails | error |
| `card_f32_bzero_trail_2d` | `10.0 x` | +10 | fails | error |
| `card_f32_bzero_lowerd_2d` | `1.0d1` | +0 | +10 | error |
| `card_f32_bzero_lowere_2d` | `1.0e1` | +10 | +10 | +10 |
| `card_f32_bzero_lowkey_2d` | `bzero` | +0 | +10 | error |
| `card_f32_bzero_empty_2d` | no value | +0 | fails | error |
| `card_f32_bzero_{int,plus,dotlead,dottrail}_2d` | `10`, `+1.0E1`, `.5E1`, `10.` | +10/+10/+5/+10 | same | read |
| `card_i16_bscale_dup_2d` | BSCALE 2.0, BZERO 32768, BSCALE -1.0 | BSCALE -1 | BSCALE -1 | last card |
| `card_i16_blank_dup_2d` (BSCALE 2) | BLANK 1, BLANK 100 | null at 100 | NaN at 100 | last card |
| `card_i16_blank_float_2d` | BLANK 100.0 | ignored | ignored | error |
| `card_i16_blank_str_2d` | BLANK 'abc' | ignored | ignored | error |
| `card_i16_blank_noeq_2d` | `BLANK   =100` | null at 100 | ignored | error |

fitsverify 4.7.0 reports a warning (duplicated keyword) for the dup files
and 1 error for each malformed one, `1.0e1` included (lower-case exponent).
`comp_card_i16_bscale_dup_2d.fits.fz` (fpack -C of the BSCALE dup sample;
fpack copies both cards): on a tile-compressed image cfitsio
(imcomp_get_compressed_image_par, ffgky) and astropy apply the FIRST
BSCALE (2.0), so ReadFloat32 does too.

### Tile-compressed samples

`make_compressed.sh` fpacks (RICE_1, `-C`) `u16_blank`, `i32_blank_bscale`,
`u8_blank`, `asi533_u16_crop`, `i16_bscale_neg` (BSCALE -1),
`i16_bscale0` (BSCALE 0), `i8_signed` (BZERO -128) and `u16_rgb` (the
planar RGB cube, read by `ReadFrameRGB` through the compressed HDU-1
path). fpack copies the
BLANK, BSCALE and BZERO cards into the compressed header. Then
`zblank_variants.py` (stdlib only; rewrites header cards and the main
table, copies the compressed heap verbatim) derives the other null
conventions of the tiled-image convention as cfitsio reads them
(`imcomp_get_compressed_image_par`, cfitsio 4.7.0 imcompress.c:5534-5568:
a ZBLANK column, else the ZBLANK keyword, else BLANK; BSCALE/BZERO
keywords; no ZBSCALE/ZBZERO):

| file | change |
|------|--------|
| `comp_u16_zblank.fits.fz` | BLANK card renamed ZBLANK |
| `comp_i32_zblank_bscale.fits.fz` | same, BITPIX 32 with BSCALE 2.5 / BZERO 1e6 |
| `comp_u16_zblank_blank.fits.fz` | ZBLANK = -32768 and BLANK = -32767 (cfitsio uses ZBLANK; pixel 1 stores -32767 and is valid) |
| `comp_u16_zblank_col.fits.fz` | a per-tile ZBLANK column (TFORM 1J) instead of the keyword |
| `comp_u8_zbscale.fits.fz` | comp_u8_blank plus non-standard ZBSCALE = 2 / ZBZERO = 5 cards (cfitsio ignores them) |
| `comp_card_i16_bscale_str.fits.fz` | comp_i16_bscale_neg with `BSCALE  = 'abc'` (cfitsio uses BSCALE 1, astropy 8.0.2 fails: refused by the scaling-card rule) |
| `comp_card_i16_bscale_hierarch.fits.fz` | comp_i16_bscale_neg with its BSCALE card replaced by `HIERARCH BSCALE = -1.0` (cfitsio and astropy apply it on this tile-compressed header; astropy ignores it on a plain image, so the one scaling-card rule refuses it) |

Both readings of the three `comp_card_*` samples are recorded under
`refs` in `readfloat32_expected.json` (`gen_expected.py comp_card_refs`:
cfitsio's cfnull golden and astropy's `getdata`); the generator asserts
that the accepted one (`comp_card_i16_bscale_dup_2d`, first card) is read
alike by both and that `comp_card_i16_bscale_str` is not.

All pass `fitsverify` (`comp_card_i16_bscale_dup_2d` with its
duplicate-keyword warning). `testdata/cfref/cfref dump` segfaults inside
libcfitsio (`imcomp_decompress_tile` via `ffgpfj`, a TLONG read) on
`comp_i16_bscale_neg` and `comp_i16_bscale0`, so the compressed truth uses
the small `cfnull` harness below (TDOUBLE reads), which does not crash.

`u64.fits` uses the real-valued spelling `BZERO = 9.223372036854776E+18`;
`u64_intcard.fits` uses the integer literal `9223372036854775808` (what
cfitsio and astropy write). Both are in the expected matrix.

### Real camera frame

`asi533_u16_crop.fits` — REAL data: a 128x128 crop (origin x=1440, y=1440,
0-based) of a 60 s light frame of HD 16582 from the repository author's
ZWO ASI533MC Pro (Bayer RGGB, BITPIX=16, BZERO=32768, ROWORDER='TOP-DOWN',
written by INDI), source sha256
`2f4ae03a13ab74f7c6c6494ad8de9d187bd4e645b86f54ce566b555aeee5922d`. The data
rows are copied byte for byte (`crop_camera.py`); the header is kept except
SITELAT, SITELONG, OBJCTAZ, OBJCTALT and AIRMASS, which were removed
because with DATE-OBS they locate the observing site; NAXIS1/2 and
CRPIX1/2 adjusted for the crop; HISTORY records the crop. The 18 MB source
frame is not committed.

## Truth

- `readfloat32_expected.json` (`gen_expected.py`): pixels read with astropy
  `do_not_scale_image_data=True`, mapped by the documented contract with
  exact rationals (`fractions.Fraction`, correctly rounded to float32) or
  numpy float32 arithmetic for the Siril-constant rule; BLANK -> null (NaN).
  For the tile-compressed samples the pixels and null flags are cfitsio's
  own (`cfnull/golden`), and the scaling keywords those of the compressed
  header; the old compressed entries came out identical to the earlier
  astropy-based values.
- `cfnull/` : `cfnull.c` opens each comp_*.fz with libcfitsio 4.7.0 and
  writes every pixel of the first NAXIS>0 HDU as `fits_read_pixnull`
  (TDOUBLE) returns it, null where cfitsio sets its null flag. Note that
  cfitsio returns zeros, without an error, for the BSCALE = 0 file;
  ReadFloat32 rejects it.
- `sirilref/` : `sirilref.c` reproduces Siril's `readfits(..., force_float
  = TRUE)` on the first image HDU: `manage_bitpix`, `fit_stats`,
  `convert_floats` copied verbatim from Siril `src/io/image_format_fits.c`,
  constants from `src/core/siril.h:70-77`, linked against libcfitsio. Its
  JSON per sample is the bit pattern of every float Siril produces (or
  `unsupported` for BITPIX=64, or the cfitsio status when Siril's read
  fails). Siril is GPL-3.0-or-later; the harness is a test tool only.

### Encoding of large truth arrays

Arrays of more than 1024 values (`readfloat32_expected.json` "bits",
`stdimage_expected.json` "values", `sirilref/golden` "bits": the camera
crops and the two 65536-value files) are stored as an exact digest
`{"n", "sha256", "sample"}` instead of every value: sha256 over the
big-endian int64 of every value (a float32 bit pattern or pixel value;
null/NaN as -1 in readfloat32_expected; every NaN as 0x7fc00000 in a
float stdimage digest) and `[index, value]` samples at the first and last
4 indices, 32 evenly spaced ones and the first minimum and maximum.
`gen_expected.py` and `sirilref.c` (the SHA-256 of `testdata/cfref/cfref.c`)
write the same digest; the tests (`imaging_codec_test.go`) compare their
full result with it, so a single differing pixel fails. This took
`sirilref/golden` from 3.0 MB to 0.1 MB, `readfloat32_expected.json` from
386 KB to 29 KB and `stdimage_expected.json` from 480 KB to 6 KB. The
`cfnull/golden` dumps stay explicit (`gen_expected.py` reads cfitsio's
values from them element by element); the largest is 94 KB.
