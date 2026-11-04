# Tile-compression read corpus: provenance

This file covers the compressed-image samples in `testdata/`:

- `comp_fz_<case>.fits` - the compressed READ corpus (148 cases, plus 12
  ZBITPIX 64 cases), and `comp_fzsrc_<src>.fits`, the
  source images that a test compares against directly. The corpus has a
  sample for every algorithm x ZBITPIX x BZERO convention x ZQUANTIZ cell
  that cfitsio can write (`TestCompressedCorpusComplete` reads the cells
  from the goldens; the cells cfitsio cannot write are listed there with
  the reason: PLIO uint32, RICE/HCOMPRESS lossless floats, HCOMPRESS with
  SUBTRACTIVE_DITHER_2).
- `comp_fzun_<case>.fits`, the funpack twin of each case (astropy's decode
  for ZBITPIX 64), is NOT committed: `gen_twins.sh` regenerates every twin
  from the committed samples (into `cfref/twins.out`, run by
  `make -C testdata/cfref goldens` and by gen_fpack.sh), and what the tests
  use of it is its cfitsio golden, `cfref/golden/twins/comp_fzun_<case>
  .fits.json`. A twin byte-identical to an earlier one by name has no
  golden; `comp_fzun_same.txt` names the kept one.
- `comp_bad/u16_rice_lowercase.fits` - `comp_fz_u16_rice.fits` with the
  ZCMPTYPE card byte-patched to `'rice_1'` (`gen_fpack patch`; cfitsio will
  not write it). cfitsio compares ZCMPTYPE case-sensitively and refuses the
  HDU (funpack 1.7.0 (CFITSIO 4.7.0): status 414 "Unknown image compression type"); its
  golden `comp_bad/golden/u16_rice_lowercase.fits.json` (`cfref dump`)
  stops at the primary HDU. It lives in a subdirectory so the top-level
  sample loops, which compare HDU counts with cfitsio, do not see it.
- `comp_bad/u16_rice_ztile0.fits` (ZTILE1 byte-patched to 0) and
  `comp_bad/u16_rice_nocompdata.fits` (TTYPE1 'COMPRESSED_DATA' patched to
  'COMPRESSED_DATX'): cfitsio will not move to either HDU (goldens stop at
  the primary); the library refuses both (`TestCompressedBadSamplesRefused`).
- `comp_bad/i32_plio_bzero32768.fits`, `comp_bad/u8_plio_bzero32768.fits`
  (`setkey BZERO 32768` on the PLIO_1 ZBITPIX 32 / 8 corpus files):
  cfitsio's reader applies its uint16 shift at any ZBITPIX while its
  writer adds it only for 16-bit images; the goldens show cfitsio reading
  values, the library refuses.
- `comp_bad/ztable_tbl500.fits` - `fpack -table` of a 500-row cfitsio
  table (`gen_fpack.c` `tbl500`): tiled TABLE compression (ZTABLE = T),
  refused with `ErrTiledTableCompression` (`TestTiledTableCompressionRefused`;
  the golden shows the 1-row compressed-byte view a plain reader gets).
- `comp_bad/u16_gzip1_extrabyte.fits` (`gen_fpack rewidth ... 2 1`: every
  GZIP_1 tile one zero byte longer than 2 x tilelen): cfitsio reads the
  header but every pixel read fails (status 414, "uncompressed tile has
  wrong size"; funpack refuses it); the library refuses it.
- `comp_bad/i16_gzip1_w4_overflow.fits` (`gen_fpack rewidth ... 4 0 70000`:
  4-byte GZIP_1 tiles on the int16 BLANK image, 70000 added to the first
  value): cfitsio's TDOUBLE read (the golden) returns the out-of-range
  value, a read as the image type fails with status 412 and funpack
  refuses the file; the library refuses it for every target type.
- `comp_bad/u16_rice_znaxis0.fits`, `comp_bad/u16_rice_znaxism1.fits`
  (`ZNAXIS` byte-patched to 0 and -1) and `comp_bad/gzip1_7d.fits` (a
  7-axis GZIP_1 image written by astropy 8.0.2, `generate_compressed.py`):
  cfitsio will not move to any of them (BAD_NAXIS, imcompress.c:5370-5379;
  goldens stop at the primary); the library refuses them.
- `comp_<name>.fits` / `comp_raw_<name>.fits` - the astropy-written fixtures
  from `generate_compressed.py`.

Every one of these files has a cfitsio golden in `cfref/golden/` (see
`cfref/PROVENANCE.md`). The golden is what libcfitsio decodes from the file,
and it is the truth that the Go tests compare against bit-for-bit
(`TestCompressedCorpusRead`, `TestCompressedFixtures`).

## Tool versions

| tool | version |
|------|---------|
| libcfitsio | 4.7.0 (Arch `cfitsio 1:4.7.0-1`) |
| fpack / funpack | 1.7.0 (Dec 2013), linked to CFITSIO 4.070 |
| gcc | 16.2.1 (Arch `gcc 16.2.1+r23+gd564253eb6c8-1`), `-O2 -Wall -Wextra -std=c99` |
| Python | 3.11.11 (pyenv), the repo venv |
| astropy | 8.0.2 |
| numpy | 2.4.6 |

## Regenerating

```sh
FITS_ASTROPY_PYTHON=$PWD/venv/bin/python sh testdata/gen_fpack.sh   # comp_fz_*, comp_fzsrc_*, comp_bad/*, comp_fzun_same.txt
venv/bin/python testdata/generate_compressed.py   # comp_*, comp_raw_*
FITS_ASTROPY_PYTHON=$PWD/venv/bin/python make -C testdata/cfref goldens   # cfitsio truth: samples, twins (gen_twins.sh), comp_bad
```

`gen_fpack.sh` needs `cc`, libcfitsio, `fpack` and `funpack`. It also needs
a Python with astropy (set by `FITS_ASTROPY_PYTHON`, default `python3`), which
it uses only to find astropy's test-data directory for the real samples.
Running each generator twice gives byte-identical output (checked with
sha256 over every file), for these reasons:

- The source pixels come from a fixed 64-bit LCG in `gen_fpack.c`.
- Every quantization uses an explicit dither seed (`fpack -qN`,
  `fits_set_dither_seed`, astropy `dither_seed=`).
- `fpack`/`funpack` run with `-C`, because CHECKSUM comments carry a
  timestamp.
- `generate_compressed.py` pins the gzip header mtime to 0, because astropy's
  GZIP codecs otherwise stamp the current time into every tile.

## Corpus construction (gen_fpack.sh / gen_fpack.c)

- Source images (40x24, plus 20x12x4 cubes, a 10x6x4x2 4-D and a
  5x4x3x2x2x2 6-D int16 image) are written by libcfitsio. They
  cover u8, s8 (BZERO -128), i16, u16 (BZERO 32768), u32 (BZERO 2^31), i32,
  i16 with BLANK, i16 with BSCALE 0.5 / BZERO 100, a PLIO mask, f32, f32 with
  NaN, f32 with zeros, f32 with one constant row, f64, f64 with NaN, f64 with
  one constant row, and a uint16 PLIO image with BLANK = -32743
  (`u16blank_plio`: cfitsio's fits_read_img compares BLANK with the PLIO
  value, stored + 32768, and flags no pixel; funpack's plain twin nulls the
  30 pixels whose stored value is -32743).
- cfitsio 4.7.0's `fits_read_img` refuses compressed images of more than 3
  axes (status 414, "only 1D, 2D, or 3D images are currently supported");
  `cfref` reads them through `fits_read_subset` instead, which is right for
  4-D (`cube4d16_rice`) but returns shifted pixels and zeros for 6-D
  (`cube6d16_gzip2`), whose truth is therefore the funpack twin (equal to
  astropy 8.0.2's decode and the source).
- `fpack` makes most cases. The exact flags are listed in `gen_fpack.sh`.
- fpack 1.7.0 (CFITSIO 4.7.0) corrupts its heap and aborts (`malloc(): corrupted top
  size`, `free(): invalid next size`) on PLIO_1 images whose rows PLIO
  expands: a uint16 image with row tiles (the default), 40x8, 20x8, 20x24
  or 10x24 tiles, and an int32 image with pixel-to-pixel noise. The
  `u16_plio` case therefore uses `-t 40,6` (one 40x24 tile would hit the
  cfref single-tile + BZERO crash below), and the `i32mask` source uses
  runs of 5 equal pixels. Both were checked against funpack and the
  source (0 of 960 pixels differ).
- fpack cannot make some cases, so `gen_fpack.c` makes them with libcfitsio:
  - `hcomp`: HCOMPRESS with `fits_set_hcomp_smooth` (SMOOTH=1).
  - `zblankcol`: moves the null value into a per-tile ZBLANK column. For
    i16, tile r gets -32768 + 7(r-1).
  - `uncompcol`: writes tiles into a pre-2011 UNCOMPRESSED_DATA column of
    the ZBITPIX type (1PB/1PI/1PJ/1PE/1PD), holding the stored values read
    back by libcfitsio with BSCALE/BZERO off. With no rows given, it moves
    the GZIP_COMPRESSED_DATA fallback tiles of a float image. With rows, it
    moves those tiles out of COMPRESSED_DATA (descriptor set to length 0),
    which is how the integer cases are made.
  - `editkey`: removes ZDITHER0 or ZTILEn, or sets a floating-point keyword:
    BSCALE 2 / BZERO 10 on a quantized float image (fpack refuses to
    compress a float image that has BSCALE: "Implicit datatype conversion
    is not supported"), and ZSCALE 2 / ZZERO 10 keywords on a lossless one.
  - `setkey`: sets an integer or string keyword: ZBLANK 26 next to BLANK
    -32768, ZQUANTIZ 'FUTURE_METHOD' on a lossless image, swapped
    ZNAME1/ZNAME2 on RICE_1 and HCOMPRESS_1 images, the pre-2009
    ZNAME2 'NOISEBIT' with ZVAL2 16 (`f32_rice_noisebit`, cfitsio then uses
    BYTEPIX 4), and ZVAL1 2 / ZVAL2 32 (`u16_rice_zvalswap`, which cfitsio
    swaps back to BLOCKSIZE 32 / BYTEPIX 2).
  - `rewidth`: re-encodes every tile of an integer GZIP_1/GZIP_2/RICE_1
    file (row tiles) with tile integers W bytes wide, using libcfitsio's
    own coders (`compress2mem_from_mem`, `fits_rcomp`/`_short`/`_byte`) on
    the stored values cfitsio reads back with BSCALE/BZERO off: the layout
    old cfitsio wrote (16-bit images in 4-byte GZIP tiles; RICE_1 without
    ZVAL2, i.e. BYTEPIX 4). cfitsio takes the tile type from the inflated
    length (GZIP, imcompress.c:6372-6417) or BYTEPIX (RICE, :6335-6346).
    Cases `u8_gzip1_w4`, `u8_gzip2_w2`, `u16_gzip2_w4`, `i16blank_gzip1_w8`,
    `i32_gzip2_w8`, `u8_rice_bp4`; their funpack twins equal the original
    twins and astropy 8.0.2 decodes each one to the same values.
  - The BITPIX -64 twins of the scaled -32 cases (`editkey`):
    `f64scaled_rice` (quantized + BSCALE 2 / BZERO 10),
    `f64flat_rice_scaled` and `f64flat_uncomp_scaled` (GZIP_COMPRESSED_DATA
    and UNCOMPRESSED_DATA fallback tiles + BSCALE/BZERO) and
    `f64_lossless_zscalekw` (lossless + ZSCALE 2 / ZZERO 10 keywords).
  - In each case libcfitsio reads the result back to produce the golden, so
    the truth is still cfitsio's.
- The real samples come from astropy 8.0.2's test data (BSD-3-Clause;
  https://github.com/astropy/astropy, `astropy/io/fits/`). They are copied
  verbatim, except that the row-tiled M13 files (gzip, rice, plio) are cut
  to their first 40 of 300 tiles by `gen_fpack croprows`: the kept table
  rows and their compressed tile bytes are the original 2006 ones
  (checked: every kept COMPRESSED_DATA cell equal to the source's, the
  decode equal to the first 40 rows of m13_hcomp, fitsverify OK); ZNAXIS2
  becomes 40, the heap is compacted and CHECKSUM/DATASUM/ZDATASUM, which
  no longer match, are removed. The property they prove (old-cfitsio
  4-byte GZIP tiles, RICE_1 without ZVAL2, PLIO_1) survives. The sha256
  column is the astropy source file:

  | corpus file | astropy file | sha256 |
  |-------------|--------------|--------|
  | `comp_fz_real_with_nan.fits` | `hdu/compressed/tests/data/compressed_with_nan.fits` (fpack'd by CFITSIO in 2022; BITPIX -64, SUBTRACTIVE_DITHER_1, ZDITHER0 2961, 6x6 tiles, source `arange(144)` with NaN at [1,1]) | `d7022e3e092abb2ecab24cd58c14201332b0b3111c2079668e26d42d50f4a6a7` |
  | `comp_fz_real_image.fits` | `tests/data/compressed_image.fits` (10x10 int16 RICE_1) | `dc08d7fd26b8ebda5db7d9fbd3dc889999face95fbe42af6a6f418735d3e04fe` |
  | `comp_fz_real_m13_gzip.fits` | `hdu/compressed/tests/data/m13_gzip.fits` (M13, a 2006 SkyView 300x300 int16 image; GZIP_1 with 4-byte tiles: 1200 inflated bytes per 300-pixel row) | `5ecfcdac4ba05bfbdda1912b17b165eb9cd445c6f268f89fa1ca13467bb13914` |
  | `comp_fz_real_m13_rice.fits` | `hdu/compressed/tests/data/m13_rice.fits` (RICE_1 without ZVAL2, so BYTEPIX 4 on ZBITPIX 16) | `de6d2f69750e11fb053c024ebe4efddbc49c4e70b5e270c3594f096a9fcdbcaa` |
  | `comp_fz_real_m13_hcomp.fits` | `hdu/compressed/tests/data/m13_hcomp.fits` (HCOMPRESS_1, one 300x300 tile, scale 0) | `cfac7cefe161d12a4187a507eb82109fad1630654a033a941c5ffc42045d2e38` |
  | `comp_fz_real_m13_plio.fits` | `hdu/compressed/tests/data/m13_plio.fits` (PLIO_1, row tiles) | `9aa01f3237b11cbd851e232b75048e1dd17d11212c4515e3eb97ec39647369c3` |

  The four M13 source files decode to the same 90000 pixels in libcfitsio, funpack
  1.7.0 and astropy 8.0.2, equal to astropy's uncompressed `m13.fits`
  (sha256 `eb3e208edbe302cae0ea45d17ab618930d85847da3f5e6ffd53d9410ec0a5a45`,
  not committed: the funpack twins carry the same pixels).

  astropy's `tests/data/compressed_float_bzero.fits` was not used, because
  cfitsio/funpack cannot read it ("1st key not SIMPLE or XTENSION"), so there
  is no truth for it.

## What cfitsio 4.7.0 itself gets wrong (recorded, not hidden)

- funpack writes wrong pixels for quantized float images that have NaN and 2D
  tiles: `comp_fz_f32nan_rice_t` (840 of 960 pixels; NaN positions differ
  from the source) and the real file `comp_fz_real_with_nan` (30 pixels, NaN
  lost). The library read (`fits_read_img`, which is what the golden holds)
  matches the source. The test therefore uses the source as the arbiter and
  flags these twins `twinWrong`.
- The funpack twin of `comp_fz_i16_zblankcol` has no nulls. A single BLANK
  keyword cannot express per-tile null values, so the test compares that twin
  on the non-null pixels only.
- libcfitsio segfaults in `imcomp_decompress_tile` when an HCOMPRESS image
  with BZERO is re-read through the same handle after `fits_set_bscale`
  (its tile cache). The u16 and u32 HCOMPRESS cases (including the lossy
  `u32_hcomp_s4` and `u32_hcomp_s1024`) use `-t 40,8`, which gives 3
  tiles; `comp_fz_u16sat_hcomp_s2` crashed the harness with any tiling, so
  the cfref harness now reads the raw values of a compressed image through
  a second handle (cfref/PROVENANCE.md).
- Fallback tiles (GZIP_COMPRESSED_DATA, UNCOMPRESSED_DATA) in an image
  with BSCALE/BZERO or BLANK: `fits_read_img` returns them unscaled and
  without null checking (imcompress.c `imcomp_decompress_tile` reads them
  with `fffr4r4(..., 1., 0., ...)` or `fits_read_col`), while it scales and
  null-checks every other tile of the same image. funpack writes all tiles'
  stored values under the image's BSCALE/BZERO/BLANK, and astropy decodes
  the same way. The Go reader follows funpack/astropy. In these rows the
  test uses the funpack twin as truth and requires the cfitsio golden to
  still differ: `comp_fz_f32flat_rice_scaled`, `comp_fz_f32flat_uncomp_scaled`
  (row 5: cfitsio 7.25, funpack/astropy/ours 24.5), `comp_fz_u16_uncomp`
  (row 4: cfitsio returns the stored int16, without the 32768) and
  `comp_fz_i16blank_uncomp` (rows 1 and 4: cfitsio does not flag BLANK).
- The funpack twin of `comp_fz_i16blank_zblankkw` keeps BLANK -32768 and
  drops ZBLANK 26, so its null flags differ from the compressed file's. The
  test compares that twin only where both are valid.
- cfitsio cannot tile-compress BITPIX 64 (fpack: status 113) and funpack
  cannot decompress it ("unsupported datatype for uncompressing image").
  astropy 8.0.2 writes it with GZIP_1/GZIP_2/NOCOMPRESS (RICE_1,
  HCOMPRESS_1 and PLIO_1 data it converts to ZBITPIX 32), and libcfitsio's
  TLONG `fits_read_img` reads the GZIP ones; its uint64 (BZERO 2^63) TLONG
  read overflows (status 412) and its stored values (raw) are used. The
  ZBITPIX 64 samples (`comp_fz_i64_*`, `comp_fz_u64_*`, `comp_fz_i64small_*`)
  are made inside gen_fpack.sh by an inline astropy script; their twins are
  astropy's decode as plain images. RICE_1 (BYTEPIX 4), HCOMPRESS_1 and
  PLIO_1 at ZBITPIX 64 (`comp_fz_i64_rice`, `_hcomp`, `_plio`) are built by
  `gen_fpack mk64` from libcfitsio's own coders (fits_rcomp,
  fits_hcompress64, pl_p2li) into a hand-made table; cfitsio decodes them
  into C ints and widens them (its TLONG read is the truth, astropy 8.0.2
  decodes them to the same values: the twin).
- cfitsio's TDOUBLE `fits_read_img` of an integer image stored in 8-byte
  GZIP tiles (`comp_fz_i16blank_gzip1_w8`, `comp_fz_i32_gzip2_w8`) returns
  its tile buffer unwritten: imcomp_decompress_tile takes the 8-byte tile
  as a double array and has no branch for it in an integer image, so the
  pixels are uninitialized memory (status 0). Earlier goldens showed zeros
  only because the heap happened to be fresh; the cfref harness now runs
  with glibc's malloc fill (cfref/PROVENANCE.md), so the goldens hold the
  0xA5 byte pattern and the tests pin exactly that. Its TLONG read (the
  golden's `i64`), funpack and astropy give the source values, so the
  corpus test uses the twin as truth there (`twinTruth`).
- `comp_nocompress_*.fits` (NOCOMPRESS, u8/i16/i32/i64/f32/f64): cfitsio cannot
  decode it ("unknown compression algorithm"), so it is checked against its
  astropy raw twin only (for the floats, which astropy quantizes, the twin
  is astropy's own decode).
- 8-byte GZIP tiles on an 8-bit image (`comp_fz_u8_gzip2_w8`): cfitsio's
  TDOUBLE read returns unwritten memory and funpack zeros (no TDOUBLE-tile
  branch for a TBYTE output); its TLONG read, equal to astropy and the u8 twin, is the
  truth (`i64Truth`).
- ZQUANTIZ = 'NONE' with ZSCALE/ZZERO columns (`comp_fz_f32_rice_zqnone`):
  cfitsio and funpack take the tile integers as IEEE floats scaled by the
  tile's ZSCALE/ZZERO; astropy 8.0.2 dequantizes them instead. The library
  follows cfitsio. One tile integer is -1, a NaN pattern: `fits_read_img`
  returns that NaN as a valid value, because it does not null-check tiles
  of an image without ZBLANK; funpack's plain twin flags it (rawRows).
- A null value on tiles narrower than 4 bytes is cast to the tile type
  before the compare (`(unsigned char)` / `(short)` tnull): BLANK -32768
  nulls the zeros of 1-byte tiles (`comp_fz_i16small_*`), BLANK 65536 the
  zeros of 2-byte tiles. funpack does not null-check integer tiles, so its
  twins are compared where both are valid.

## Round-2 additions and size

- New cases (gen_fpack.sh "round 2" block): `f64_rice_wrap` (one 102x100
  -64 SD1 tile: the 10000-value dither sequence wraps inside it),
  `f32/f64flatnan_rice/_uncomp` (NaN in fallback tiles), `f32_rice_zqnone`,
  the tile-width cases (`gen_fpack rewidth`/`rewidthq` on new 40x6 sources
  `i16small`/`i32small` with values 0..255), `f64_rice_bp8` (setkey ZVAL2 8),
  ZBITPIX 64, and `comp_bad/u16_rice_emptytile`, `u16_rice_gzipcol`,
  `u16_rice_extrarow`, `u16_uncomp_short` (cfitsio's behaviour in their
  goldens). New sources are written after every earlier G() call, so the
  earlier sources and samples are unchanged.
- The astropy fixtures (`generate_compressed.py`) are 40x37 (were
  100x100) and the odd-size one 23x19 (was 103x97); multitile uses 16x16
  tiles (partial edge tiles kept). Byte-identical raw twins are written
  once (`raw_of`).
- Every comp_* / comp_bad golden reproduces from regenerated samples:
  `TestCompressedGoldensReproduce`.

## Round-3 additions

gen_fpack.sh "round 3" block, all made by libcfitsio 4.7.0 / fpack /
funpack 1.7.0 from existing corpus files (nothing earlier changes, except
that three byte-identical funpack twins are now kept under the new case's
name, `comp_fzun_same.txt`):

- `*_edge` (`gen_fpack rewidth ... i:v` sets tile-1 values, `croprows` to
  2 rows): the lowest and highest stored value of ZBITPIX 8/16/32 in 2-,
  4- and 8-byte tiles; `comp_bad/*_lo` / `*_hi` hold one past each edge
  (cfitsio TDOUBLE/TLONG return the value, funpack fails with 412).
- `comp_bad/i16_gzip1_w4_nullout`: stored 37232 = BLANK in 4-byte tiles of
  a 16-bit image (cfitsio flags it null; funpack fails 412; astropy 8.0.2
  flags nothing). A BLANK outside ZBITPIX is refused by the library's
  scaling-card rule (`ErrInvalidScaling`); the golden records cfitsio.
- `comp_bad/i32_rice_blank2p32`, `i64_gzip1_blank2p40`,
  `i16_rice_zblank2p33`, `i16_rice_zblankcol1k` (`zblankcol ... 1K v`):
  null values outside a C int. cfitsio flags nothing / falls back to BLANK
  / fails the read with 412; astropy 8.0.2 refuses the first two.
- `f32_rice_bp1_zblank406`: ZBLANK 406 on 1-byte quantized tiles
  (`(unsigned char)` 406 = 150 flags 6 pixels).
- `i16small_hcomp_blank65536`, `i16small_plio_blank65536`: BLANK 65536 is
  compared with 4-byte ints, no pixel null (cfitsio). These two,
  `i16small_rice_blankcast` (cfitsio casts 65536 to the short 0; astropy
  8.0.2 fails: out of bounds for int16) and `i16_zblankcol_blank2p32`
  carry a BLANK outside ZBITPIX 16, which the library's scaling-card rule
  refuses (`fzCase.blankRefused`); the goldens record cfitsio.
- `f32_gzip1_q_w1`, `f32_gzip2_q_w2`, `f64_gzip1_q_w2` (`rewidthq` now
  re-deflates quantized GZIP tiles at 1/2/8 bytes; astropy 8.0.2 agrees with
  cfitsio on all three) and `comp_bad/f32_gzip1_q_w8` (cfitsio returns
  unwritten memory, funpack and astropy zeros).
- `f64_rice_zqnone`: ZQUANTIZ NONE on a -64 RICE file (cfitsio fffr4r8).
- ZBITPIX 64 with BLANK 37: `i64small_gzip1_blank` (8-byte tiles),
  `_w2_blank`, `_w1_blank`; the twin is astropy's float64 decode (NaN at
  the null).

## Round-2 fix pass additions

gen_fpack.sh "round 2 fix pass" block (nothing earlier changes except
`comp_fz_f32_rice_noisebit`, whose ZVAL2 is now 1 instead of 16: with
ZNAME2 'NOISEBIT' cfitsio ignores ZVAL2 and decodes BYTEPIX 4, so the file
reads as `comp_fz_f32_rice_sd1` it was made from, while a reader taking
ZVAL2 would decode 1-byte Rice; ZVAL2 16 proved nothing, as the default rule
reads it as 4 anyway):

- `gen_fpack.c src` adds `u8sat`, `i16sat`, `u16sat`, `i32sat` (40x6, noise
  of +-20 / +-300 below the type maximum, a 16x4 block at the maximum; own
  LCG). fpack -h -s 2 -t 20,6: `comp_fz_u16sat_hcomp_s2`,
  `comp_fz_i16sat_hcomp_s2` (13 pixels reconstructed above the maximum:
  cfitsio's TDOUBLE read returns them unclipped, its TSHORT read and
  funpack clip them, its TUSHORT read of the uint16 file fails with status
  412; funpack's twin is the truth),
  `comp_bad/u8sat_hcomp_s2` (funpack fails 412), `comp_bad/i32sat_hcomp_s2`
  (cfitsio and funpack return the saturated pixels wrapped to about -2^31).
- `comp_bad/u16_plio_bzero0` (BZERO 0 on the uint16 PLIO file),
  `i32_plio_zbitpix8`, `i32_plio_zbitpix16_bzero32768` (i32_plio cut to 6
  rows, ZBITPIX / BZERO edited): PLIO values the image type cannot hold;
  cfitsio's TDOUBLE read returns them, funpack fails 412.
- `comp_fz_i32_gzip1_blankmin`, `_blankmax` (4-byte GZIP_1 tiles, pixel 0 =
  -2^31, pixel 1 = 2^31 - 1, 2 rows; BLANK at either end),
  `comp_fz_i32_gzip1_zblankcolmin` (1K ZBLANK column -2^31),
  `comp_bad/i32_gzip1_blank2p31`, `_blankm2p31m1` (cfitsio flags nothing),
  `comp_bad/i32_gzip1_zblankcol2p31`, `_zblankcolm2p31m1` (cfitsio fails 412).
- `comp_fz_i16_zblankcol_blank2p32`: BLANK 2^32 next to a ZBLANK column,
  which cfitsio uses.
- `comp_bad/i16_gzip1_w4_nullout2` (`rewidth ... 0:70000,2@0:70001`, per-tile
  ZBLANK 70000 / 70001, 2 rows): two null values int16 cannot hold.
- `comp_fz_i16small_hcomp_blank0`, `comp_fz_i16mask_plio_blank3`: BLANK
  values the decoded HCOMPRESS / PLIO ints hold.
- `comp_bad/u16_rice_zbitpix12`, `_zbitpix0` (cfitsio decodes the tiles
  anyway), `comp_bad/u16_rice_znaxisneg` (ZNAXIS1 -1, ZTILE1 1, ZNAXIS2 -24
  byte-patched: cfitsio refuses the HDU, 414).
- `gen_fpack mk64`: `comp_fz_i64_rice`, `_hcomp`, `_plio` (ZBITPIX 64,
  12x8, two 12x4 tiles) and `comp_bad/i64_hcomp_2p31`, `_m2p31m1` (pixel 5
  one past an int32 end; cfitsio's (int) cast wraps it).
