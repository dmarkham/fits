#!/bin/sh
# Regenerate the tile-compression READ corpus testdata/comp_fz_*.fits, the
# committed source images testdata/comp_fzsrc_*.fits, and (gen_twins.sh)
# their decompressed twins in testdata/cfref/twins.out (not committed: their
# cfitsio goldens testdata/cfref/golden/twins/ are) plus the twin map
# testdata/comp_fzun_same.txt.
#
#   sh testdata/gen_fpack.sh
#
# Truth chain (nothing here uses the Go library):
#   gen_fpack.c (libcfitsio) writes deterministic source images
#   -> /usr/bin/fpack (cfitsio) or gen_fpack.c (fits_img_compress, for the
#      HCOMPRESS SMOOTH and per-tile ZBLANK-column cases fpack cannot make)
#   -> /usr/bin/funpack (cfitsio) writes the decompressed twin (gen_twins.sh;
#      astropy for the ZBITPIX 64 samples, which funpack cannot decompress).
# After this script, refresh the cfitsio goldens for the new files (and the
# twins, which that target regenerates into testdata/cfref/twins.out) with
#   FITS_ASTROPY_PYTHON=... make -C testdata/cfref goldens
# (which also writes testdata/comp_bad/golden).
# Every float quantization uses a fixed dither seed (-qN), and fpack/funpack
# run with -C (no CHECKSUM keywords, which embed a timestamp), so re-running
# reproduces every byte. Tool versions: testdata/COMP_FZ_PROVENANCE.md.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
work=$(mktemp -d "${TMPDIR:-/tmp}/gen_fpack.XXXXXX")
trap 'rm -rf "$work"' EXIT

py=${FITS_ASTROPY_PYTHON:-python3}
cc -O2 -Wall -Wextra -std=c99 -o "$work/gen_fpack" gen_fpack.c -lcfitsio -lm
"$work/gen_fpack" src "$work"

rm -f comp_fz_*.fits comp_fzsrc_*.fits comp_bad/*.fits comp_fzun_same.txt

# Source images whose original pixels a test compares against directly
# (NaN positions, lossless round trips).
for s in f32nan f64nan f32flat f64flat; do
  cp "$work/$s.fits" "comp_fzsrc_$s.fits"
done

pk() { # pk <case> <source> <fpack flags...>
  c=$1 s=$2
  shift 2
  fpack -C "$@" -O "comp_fz_$c.fits" "$work/$s.fits"
}

# --- integer BITPIX, BZERO/BSCALE in the compressed header ---------------
pk u8_rice       u8        -r
pk u8_gzip1      u8        -g1
pk u8_hcomp      u8        -h
pk s8_rice       s8        -r            # BZERO = -128
pk s8_gzip2      s8        -g2
pk u16_rice      u16       -r            # BZERO = 32768
pk u16_gzip2     u16       -g2
pk u16_hcomp     u16       -h -t 40,8    # (single-tile + BZERO crashes cfitsio 4.7.0 tile cache on re-read; see PROVENANCE)
pk u32_rice      u32       -r            # BZERO = 2147483648
pk u32_hcomp     u32       -h -t 40,8
pk u32_hcomp_s4  u32       -h -s 4 -t 40,8    # lossy HCOMPRESS on wide-range 32-bit data:
pk u32_hcomp_s1024 u32     -h -s 1024 -t 40,8 # undigitize needs 64-bit products
pk i32_hcomp     i32       -h
pk i32_gzip2     i32       -g2
pk i16scaled_rice i16scaled -r           # BSCALE = 0.5, BZERO = 100
pk i16blank_rice i16blank  -r            # BLANK keyword -> masked read
pk i16blank_gzip1 i16blank -g1
pk i16mask_plio  i16mask   -p
pk i16_hcomp_s4  i16       -h -s 4       # lossy HCOMPRESS, scale 4 x noise
pk cube16_rice_t3 cube16   -r -t 20,12,2 # 3D tiles
pk cube16_gzip2_t3 cube16  -g2 -t 10,6,4
# Remaining algorithm x ZBITPIX x BZERO-convention cells of the claimed
# matrix (TestCompressedCorpusComplete lists the required cells).
pk s8_gzip1      s8        -g1
pk u16_gzip1     u16       -g1
pk i32_gzip1     i32       -g1
pk u32_gzip1     u32       -g1
pk u8_gzip2      u8        -g2
pk u32_gzip2     u32       -g2
pk s8_hcomp      s8        -h -t 40,8    # (single tile + BZERO: see u16_hcomp)
# PLIO_1 at every ZBITPIX cfitsio can write it for.
pk u8_plio       u8        -p
pk s8_plio       s8        -p            # BZERO = -128
pk i32_plio      i32mask   -p            # values 0 .. 2^24-1
# uint16 (BZERO = 32768): cfitsio adds 32768 before PLIO and takes it off
# on read. 40x6 tiles: with row tiles (the default), 40x8, 20x8, 20x24 or
# 10x24 fpack 1.7.0 (CFITSIO 4.7.0) corrupts its heap and aborts, and with one 40x24 tile
# the cfref read crashes (single tile + BZERO; see COMP_FZ_PROVENANCE.md).
pk u16_plio      u16       -p -t 40,6
# uint16 PLIO with BLANK: cfitsio compares BLANK with the PLIO value
# (stored + 32768), so the stored-value BLANK flags no pixel.
pk u16blank_plio u16pliob  -p -t 40,6
# 4-D and 6-D images (cfitsio's MAX_COMPRESS_DIM is 6).
pk cube4d16_rice cube4d16  -r -t 10,6,2,1
pk cube6d16_gzip2 cube6d16 -g2 -t 5,4,1,2,1,2

# --- HCOMPRESS SMOOTH (fpack has no flag; fits_set_hcomp_smooth) ---------
"$work/gen_fpack" hcomp "$work/i16.fits" comp_fz_i16_hcomp_smooth.fits 4 1
"$work/gen_fpack" hcomp "$work/i32.fits" comp_fz_i32_hcomp_smooth.fits 4 1
"$work/gen_fpack" hcomp "$work/f32.fits" comp_fz_f32_hcomp_smooth.fits 2 1 4 17

# --- float, quantized ----------------------------------------------------
pk f32_rice_sd1     f32     -q17 4        # SUBTRACTIVE_DITHER_1, ZDITHER0 17
pk f32_rice_sd2     f32zero -qz17 4       # SUBTRACTIVE_DITHER_2, zeros kept
pk f32_rice_nodith  f32     -q0 4         # NO_DITHER
pk f32_rice_abs     f32     -q23 -0.05    # absolute step 0.05
pk f32_gzip1_sd1    f32     -g1 -q17 4
pk f32_gzip2_sd1    f32     -g2 -q17 4
pk f32nan_rice_sd1  f32nan  -q17 4        # NaN -> ZBLANK
pk f32nan_rice_t    f32nan  -q17 4 -t 20,8
pk f32_hcomp_q      f32     -h -q17 4     # HCOMPRESS quantized float
pk f32nan_hcomp_q   f32nan  -h -q17 4
pk f32flat_rice     f32flat -q17 4        # constant row -> GZIP_COMPRESSED_DATA
pk cubef32_rice_t3  cubef32 -q17 4 -t 20,12,2
pk f64_rice_sd1     f64     -q17 4        # BITPIX -64: 4-byte quantized ints
pk f64_gzip1_sd1    f64     -g1 -q17 4
pk f64_gzip2_sd1    f64     -g2 -q17 4
pk f64_hcomp_q      f64     -h -q17 4
pk f64nan_rice_sd1  f64nan  -q17 4
pk f64flat_rice     f64flat -q17 4        # BITPIX -64 constant row -> GZIP_COMPRESSED_DATA
pk f64_rice_sd2     f64zero -qz17 4       # BITPIX -64 SUBTRACTIVE_DITHER_2, zeros kept
pk f64_rice_nodith  f64     -q0 4         # BITPIX -64 NO_DITHER
pk f64_gzip2_sd2    f64zero -g2 -qz17 4
pk f64_hcomp_nodith f64     -h -q0 4
# Every algorithm x float ZBITPIX x dither method.
pk f32_gzip1_sd2    f32zero -g1 -qz17 4
pk f32_gzip1_nodith f32     -g1 -q0 4
pk f64_gzip1_sd2    f64zero -g1 -qz17 4
pk f64_gzip1_nodith f64     -g1 -q0 4
pk f32_gzip2_sd2    f32zero -g2 -qz17 4
pk f32_gzip2_nodith f32     -g2 -q0 4
pk f64_gzip2_nodith f64     -g2 -q0 4
pk f32_hcomp_nodith f32     -h -q0 4      # (HCOMPRESS + SD2: cfitsio refuses, uses SD1)

# --- float, lossless (ZQUANTIZ = 'NONE', no ZSCALE/ZZERO) ----------------
pk f32_gzip1_lossless   f32    -g1 -q 0
pk f32_gzip2_lossless   f32    -g2 -q 0
pk f32nan_gzip2_lossless f32nan -g2 -q 0
pk f64_gzip1_lossless   f64    -g1 -q 0
pk f64_gzip2_lossless   f64    -g2 -q 0

# --- per-tile ZBLANK column ----------------------------------------------
fpack -C -r -O "$work/rb.fz" "$work/i16rowblank.fits"
"$work/gen_fpack" zblankcol "$work/rb.fz" comp_fz_i16_zblankcol.fits
"$work/gen_fpack" zblankcol comp_fz_f32nan_rice_sd1.fits comp_fz_f32nan_zblankcol.fits

# --- cases fpack cannot write directly ----------------------------------
fpack -C -q1 4 -O "$work/nozd0.fz" "$work/f32.fits"
"$work/gen_fpack" editkey "$work/nozd0.fz" comp_fz_f32_rice_nozd0.fits ZDITHER0
# f32_rice_nozd0: SD1 without ZDITHER0 (pre-2009)
"$work/gen_fpack" editkey comp_fz_f32_rice_sd1.fits "$work/sc.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/sc.fits" comp_fz_f32scaled_rice.fits BZERO 10
# f32scaled_rice: quantized float + BSCALE 2 / BZERO 10
"$work/gen_fpack" editkey comp_fz_u16_rice.fits "$work/nt.fits" ZTILE2
"$work/gen_fpack" editkey "$work/nt.fits" comp_fz_u16_rice_noztile.fits ZTILE1
# u16_rice_noztile: no ZTILEn: defaults ZNAXIS1 x 1 x ...
"$work/gen_fpack" uncompcol comp_fz_f32flat_rice.fits comp_fz_f32flat_uncomp.fits
# f32flat_uncomp: pre-2011 UNCOMPRESSED_DATA column
"$work/gen_fpack" uncompcol comp_fz_f64flat_rice.fits comp_fz_f64flat_uncomp.fits
# f64flat_uncomp: UNCOMPRESSED_DATA 1PD
# Fallback tiles in a quantized float image with BSCALE 2 / BZERO 10:
# libcfitsio's fits_read_img leaves those tiles unscaled (see the case).
"$work/gen_fpack" editkey comp_fz_f32flat_rice.fits "$work/fs.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/fs.fits" comp_fz_f32flat_rice_scaled.fits BZERO 10
"$work/gen_fpack" editkey comp_fz_f32flat_uncomp.fits "$work/fu.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/fu.fits" comp_fz_f32flat_uncomp_scaled.fits BZERO 10
# Integer tiles in UNCOMPRESSED_DATA (1PB/1PI/1PJ), chosen rows (1-based).
fpack -C -r -O "$work/i32r.fz" "$work/i32.fits"
"$work/gen_fpack" uncompcol "$work/i32r.fz" comp_fz_i32_uncomp.fits 4 11
"$work/gen_fpack" uncompcol comp_fz_u8_rice.fits comp_fz_u8_uncomp.fits 3
"$work/gen_fpack" uncompcol comp_fz_u16_rice.fits comp_fz_u16_uncomp.fits 4
# u16_uncomp: BZERO 32768
"$work/gen_fpack" uncompcol comp_fz_i16blank_rice.fits comp_fz_i16blank_uncomp.fits 1 4
# i16blank_uncomp: BLANK -32768 (row 1 has nulls)
# Keyword rules cfitsio applies (imcomp_get_compressed_image_par).
"$work/gen_fpack" setkey comp_fz_i16blank_rice.fits comp_fz_i16blank_zblankkw.fits ZBLANK int 26
# i16blank_zblankkw: ZBLANK keyword wins over BLANK
"$work/gen_fpack" editkey comp_fz_f32_gzip1_lossless.fits "$work/zs.fits" ZSCALE 2
"$work/gen_fpack" editkey "$work/zs.fits" comp_fz_f32_lossless_zscalekw.fits ZZERO 10
# f32_lossless_zscalekw: lossless + ZSCALE/ZZERO keywords
"$work/gen_fpack" setkey comp_fz_f32_gzip1_lossless.fits comp_fz_f32_lossless_futurezq.fits ZQUANTIZ str FUTURE_METHOD
# f32_lossless_futurezq: ZQUANTIZ ignored when unquantized
"$work/gen_fpack" setkey comp_fz_u16_rice.fits "$work/zn.fits" ZNAME1 str BYTEPIX
"$work/gen_fpack" setkey "$work/zn.fits" comp_fz_u16_rice_znameswap.fits ZNAME2 str BLOCKSIZE
# u16_rice_znameswap: RICE params are positional (ZVAL1, ZVAL2)
"$work/gen_fpack" setkey comp_fz_i16_hcomp_smooth.fits "$work/hn.fits" ZNAME1 str SMOOTH
"$work/gen_fpack" setkey "$work/hn.fits" comp_fz_i16_hcomp_smooth_znameswap.fits ZNAME2 str SCALE
# i16_hcomp_smooth_znameswap: HCOMPRESS smooth is ZVAL2 whatever ZNAME2 says
"$work/gen_fpack" setkey comp_fz_f32_rice_sd1.fits "$work/nb.fits" ZNAME2 str NOISEBIT
"$work/gen_fpack" setkey "$work/nb.fits" comp_fz_f32_rice_noisebit.fits ZVAL2 int 1
# f32_rice_noisebit: ZNAME2 'NOISEBIT' (pre-2009): BYTEPIX 4, ZVAL2 1 ignored
# BITPIX -64 twins of the four scaled -32 cases above.
"$work/gen_fpack" editkey comp_fz_f64_rice_sd1.fits "$work/sc64.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/sc64.fits" comp_fz_f64scaled_rice.fits BZERO 10
# f64scaled_rice: quantized -64 + BSCALE 2 / BZERO 10
"$work/gen_fpack" editkey comp_fz_f64flat_rice.fits "$work/fs64.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/fs64.fits" comp_fz_f64flat_rice_scaled.fits BZERO 10
# f64flat_rice_scaled: -64 GZIP_COMPRESSED_DATA fallback + BSCALE/BZERO
"$work/gen_fpack" editkey comp_fz_f64flat_uncomp.fits "$work/fu64.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/fu64.fits" comp_fz_f64flat_uncomp_scaled.fits BZERO 10
# f64flat_uncomp_scaled: -64 UNCOMPRESSED_DATA + BSCALE/BZERO
"$work/gen_fpack" editkey comp_fz_f64_gzip1_lossless.fits "$work/zs64.fits" ZSCALE 2
"$work/gen_fpack" editkey "$work/zs64.fits" comp_fz_f64_lossless_zscalekw.fits ZZERO 10
# f64_lossless_zscalekw: lossless -64 + ZSCALE/ZZERO keywords
# Tile integers wider or narrower than ZBITPIX (old cfitsio wrote 16-bit
# GZIP tiles as int32, and RICE_1 without ZVAL2 means BYTEPIX 4): cfitsio
# takes the tile type from the inflated length (GZIP) or BYTEPIX (RICE)
# and converts to ZBITPIX. Re-encoded by libcfitsio's own coders.
"$work/gen_fpack" rewidth comp_fz_u8_gzip1.fits comp_fz_u8_gzip1_w4.fits 4
# u8_gzip1_w4: ZBITPIX 8, 4-byte GZIP_1 tiles
"$work/gen_fpack" rewidth comp_fz_u8_gzip2.fits comp_fz_u8_gzip2_w2.fits 2
# u8_gzip2_w2: ZBITPIX 8, 2-byte shuffled GZIP_2 tiles
"$work/gen_fpack" rewidth comp_fz_u16_gzip2.fits comp_fz_u16_gzip2_w4.fits 4
# u16_gzip2_w4: uint16, 4-byte shuffled GZIP_2 tiles
"$work/gen_fpack" rewidth comp_fz_i16blank_gzip1.fits comp_fz_i16blank_gzip1_w8.fits 8
# i16blank_gzip1_w8: int16 + BLANK, 8-byte GZIP_1 tiles
"$work/gen_fpack" rewidth comp_fz_i32_gzip2.fits comp_fz_i32_gzip2_w8.fits 8
# i32_gzip2_w8: int32, 8-byte shuffled GZIP_2 tiles
"$work/gen_fpack" rewidth comp_fz_u8_rice.fits comp_fz_u8_rice_bp4.fits 4
# u8_rice_bp4: ZBITPIX 8, RICE_1 BYTEPIX 4
# ... and the layouts cfitsio refuses: a tile one byte longer than 2 x
# tilelen ("uncompressed tile has wrong size"), and 4-byte tiles holding a
# value int16 cannot hold (status 412, datatype conversion overflow).
"$work/gen_fpack" rewidth comp_fz_u16_gzip1.fits comp_bad/u16_gzip1_extrabyte.fits 2 1
"$work/gen_fpack" rewidth comp_fz_i16blank_gzip1.fits comp_bad/i16_gzip1_w4_overflow.fits 4 0 70000
# ZNAXIS outside cfitsio's 1..6 (BAD_NAXIS, imcompress.c:5370-5379).
"$work/gen_fpack" patch comp_fz_u16_rice.fits comp_bad/u16_rice_znaxis0.fits "ZNAXIS  =                    2" "ZNAXIS  =                    0"
"$work/gen_fpack" patch comp_fz_u16_rice.fits comp_bad/u16_rice_znaxism1.fits "ZNAXIS  =                    2" "ZNAXIS  =                   -1"
"$work/gen_fpack" setkey comp_fz_u16_rice.fits "$work/sw.fits" ZVAL1 int 2
"$work/gen_fpack" setkey "$work/sw.fits" comp_fz_u16_rice_zvalswap.fits ZVAL2 int 32
# u16_rice_zvalswap: ZVAL1 2 / ZVAL2 32: cfitsio swaps them back
"$work/gen_fpack" setkey comp_fz_u16_rice.fits comp_fz_u16_rice_riceone.fits ZCMPTYPE str RICE_ONE
# u16_rice_riceone: ZCMPTYPE 'RICE_ONE': cfitsio's alias for RICE_1
# ZCMPTYPE in lower case: cfitsio compares case-sensitively (FSTRCMP) and
# refuses it, so there is no funpack twin; the library refuses it too.
# (cfitsio will not write it either, so the card is byte-patched.)
"$work/gen_fpack" patch comp_fz_u16_rice.fits comp_bad/u16_rice_lowercase.fits "ZCMPTYPE= 'RICE_1  '" "ZCMPTYPE= 'rice_1  '"

# Layouts the library refuses (comp_bad/, no funpack twin; the cfitsio
# goldens in comp_bad/golden record what cfitsio makes of each).
"$work/gen_fpack" patch comp_fz_u16_rice.fits comp_bad/u16_rice_ztile0.fits "ZTILE1  =                   40" "ZTILE1  =                    0"
"$work/gen_fpack" patch comp_fz_u16_rice.fits comp_bad/u16_rice_nocompdata.fits "TTYPE1  = 'COMPRESSED_DATA'" "TTYPE1  = 'COMPRESSED_DATX'"
"$work/gen_fpack" setkey comp_fz_i32_plio.fits comp_bad/i32_plio_bzero32768.fits BZERO int 32768
"$work/gen_fpack" setkey comp_fz_u8_plio.fits comp_bad/u8_plio_bzero32768.fits BZERO int 32768
# Tiled TABLE compression (fpack -table, ZTABLE = T): not supported, refused.
fpack -C -table -O comp_bad/ztable_tbl500.fits "$work/tbl500.fits" 2>/dev/null

# --- round 2: cases the first corpus left unproven ---------------------
# BITPIX -64 SUBTRACTIVE_DITHER_1 in one 10200-pixel tile: the dither
# random sequence (10000 values) wraps inside the tile.
pk f64_rice_wrap    f64wrap    -q17 4 -t 102,100
# NaN inside the float fallback tiles (constant row 5 with two NaN), in
# GZIP_COMPRESSED_DATA and in pre-2011 UNCOMPRESSED_DATA, -32 and -64.
pk f32flatnan_rice  f32flatnan -q17 4
pk f64flatnan_rice  f64flatnan -q17 4
"$work/gen_fpack" uncompcol comp_fz_f32flatnan_rice.fits comp_fz_f32flatnan_uncomp.fits
"$work/gen_fpack" uncompcol comp_fz_f64flatnan_rice.fits comp_fz_f64flatnan_uncomp.fits
# ZQUANTIZ = 'NONE' on a file with ZSCALE/ZZERO columns: cfitsio sets
# NO_QUANTIZE (imcompress.c:5330) and takes the 4-byte tile integers as
# IEEE floats, ignoring the columns.
"$work/gen_fpack" setkey comp_fz_f32_rice_sd1.fits comp_fz_f32_rice_zqnone.fits ZQUANTIZ str NONE
# ... with BSCALE 2 / BZERO 10: cfitsio combines them with ZSCALE/ZZERO.
"$work/gen_fpack" editkey comp_fz_f32_rice_zqnone.fits "$work/zqs.fits" BSCALE 2
"$work/gen_fpack" editkey "$work/zqs.fits" comp_fz_f32_rice_zqnone_scaled.fits BZERO 10
# Tile integers of every width cfitsio converts (imcompress.c:6335-6417):
# RICE_1 BYTEPIX 1/2 on 8/16/32-bit images, GZIP 1-byte tiles on 16- and
# 32-bit, 2-byte on 32-bit and 8-byte on 8-bit images, on the 0..255
# sources (BLANK cast to the tile type: see gen_fpack.c i16small).
fpack -C -r  -O "$work/i16s_r.fz"  "$work/i16small.fits"
fpack -C -g1 -O "$work/i16s_g1.fz" "$work/i16small.fits"
fpack -C -r  -O "$work/i32s_r.fz"  "$work/i32small.fits"
fpack -C -g2 -O "$work/i32s_g2.fz" "$work/i32small.fits"
rw() { # rw <case> <in.fz> <W>
  "$work/gen_fpack" rewidth "$2" "comp_fz_$1.fits" "$3"
}
rw i16small_rice_bp1 "$work/i16s_r.fz"  1
rw i16small_gzip1_w1 "$work/i16s_g1.fz" 1
rw i32small_rice_bp1 "$work/i32s_r.fz"  1
rw i32small_rice_bp2 "$work/i32s_r.fz"  2
rw i32small_gzip2_w1 "$work/i32s_g2.fz" 1
rw i32small_gzip2_w2 "$work/i32s_g2.fz" 2
rw u8_rice_bp2       comp_fz_u8_rice.fits  2
rw u8_gzip2_w8       comp_fz_u8_gzip2.fits 8
# RICE_1 ZVAL2 = 8: cfitsio decodes every BYTEPIX other than 1 and 2 with
# the 4-byte coder, so this reads as the BYTEPIX 4 file it was made from.
"$work/gen_fpack" setkey comp_fz_f64_rice_sd1.fits comp_fz_f64_rice_bp8.fits ZVAL2 int 8
# Quantized float tiles in 2-byte Rice (cfitsio unquantize_i2r4/_i2r8).
"$work/gen_fpack" rewidthq comp_fz_f32_rice_sd1.fits comp_fz_f32_rice_bp2.fits 2
"$work/gen_fpack" rewidthq comp_fz_f64_rice_nodith.fits comp_fz_f64_rice_bp2.fits 2
# ... and in 1-byte Rice (NO_DITHER, step 0.4: the integers fit 0..255 and
# pass 127, so they must be read unsigned).
fpack -C -q0 -0.4 -O "$work/f32q2.fz" "$work/f32.fits"
"$work/gen_fpack" rewidthq "$work/f32q2.fz" comp_fz_f32_rice_bp1.fits 1
# BLANK 65536 on a 16-bit RICE image: cfitsio casts it to the 2-byte tile
# type ((short) 65536 = 0), so the zeros are null.
"$work/gen_fpack" setkey "$work/i16s_r.fz" comp_fz_i16small_rice_blankcast.fits BLANK int 65536
# Refused layouts cfitsio fails on too (goldens in comp_bad/golden): a
# tile with no payload, an integer tile in GZIP_COMPRESSED_DATA, more
# table rows than tiles, an UNCOMPRESSED_DATA tile one value short.
"$work/gen_fpack" emptytile comp_fz_u16_rice.fits comp_bad/u16_rice_emptytile.fits 3
"$work/gen_fpack" gzipcol comp_fz_u16_rice.fits comp_bad/u16_rice_gzipcol.fits 3
"$work/gen_fpack" patch comp_fz_u16_rice.fits comp_bad/u16_rice_extrarow.fits "ZNAXIS2 =                   24" "ZNAXIS2 =                   23"
"$work/gen_fpack" uncompcol comp_fz_u16_rice.fits comp_bad/u16_uncomp_short.fits -4

# --- ZBITPIX 64 ---------------------------------------------------------
# fpack and fits_img_compress refuse 64-bit images and funpack cannot
# decompress them ("unsupported datatype for uncompressing image"), but
# astropy 8.0.2 writes them, with GZIP_1 or GZIP_2 (it converts RICE_1,
# HCOMPRESS_1 and PLIO_1 data to ZBITPIX 32). gen_int64.py writes int64
# (full range) and uint64 (BZERO 2^63) images (their twins, astropy's own
# decode written as a plain image, come from gen_twins.sh: there is no
# funpack twin).
cat > "$work/gen_int64.py" <<'PY'
#!/usr/bin/env python3
"""ZBITPIX 64 samples of the compressed read corpus (run by gen_fpack.sh).

    gen_int64.py <workdir>          write comp_fz_{i64,u64}_gzip{1,2}.fits
                                    into the current directory, and
                                    <workdir>/i64small.fz

cfitsio cannot write a tile-compressed 64-bit image (fpack and
fits_img_compress refuse it) and funpack cannot decompress one, so astropy
8.0.2 writes them (GZIP_1/GZIP_2; gzip header mtime pinned to 0 for
reproducible bytes) and the twin is astropy's own decode. libcfitsio's
TLONG fits_read_img (the cfref golden's i64, or raw + BZERO for uint64,
whose TLONG read overflows) is the other truth.
"""
import functools
import gzip
import sys

import numpy as np
import astropy.io.fits as fits
import astropy.io.fits.hdu.compressed._codecs as _codecs

_codecs.gzip_compress = functools.partial(gzip.compress, mtime=0)


def comp(path, data, ct):
    hdu = fits.CompImageHDU(data=data, compression_type=ct)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path, overwrite=True, checksum=False)


def main():
    work = sys.argv[1]
    # int64 across the full range: steps of 2^40 (not exact in float64)
    # and both extremes.
    i64 = (np.arange(120, dtype=np.int64).reshape(10, 12) - 60) * (2**40) + 7
    i64.flat[1], i64.flat[2] = 2**63 - 1, -(2**63)
    # uint64 (stored with BZERO = 2^63): steps of 2^57, 0 and 2^64 - 1.
    u64 = np.arange(120, dtype=np.uint64).reshape(10, 12) * np.uint64(2**57) + np.uint64(3)
    u64.flat[1], u64.flat[2] = np.uint64(2**64 - 1), np.uint64(0)
    for name, data, ct in (("i64_gzip1", i64, "GZIP_1"), ("i64_gzip2", i64, "GZIP_2"),
                           ("u64_gzip1", u64, "GZIP_1"), ("u64_gzip2", u64, "GZIP_2")):
        comp(f"comp_fz_{name}.fits", data, ct)
    # 0..255 int64, re-encoded by gen_fpack rewidth with 1- and 4-byte
    # GZIP_1 tiles on the ZBITPIX 64 image.
    small = ((np.arange(120, dtype=np.int64) * 37) % 256).reshape(10, 12)
    comp(f"{work}/i64small.fz", small, "GZIP_1")


if __name__ == "__main__":
    main()
PY
"$py" "$work/gen_int64.py" "$work"
"$work/gen_fpack" rewidth "$work/i64small.fz" comp_fz_i64small_gzip1_w4.fits 4
"$work/gen_fpack" rewidth "$work/i64small.fz" comp_fz_i64small_gzip1_w1.fits 1
# ... with a BLANK keyword (37, a stored value): native 8-byte tiles and
# 2- and 1-byte ones. cfitsio flags the 37s (TLONG fits_read_imgnull);
# astropy's decode turns them into NaN (its twin is a float64 image).
"$work/gen_fpack" setkey "$work/i64small.fz" comp_fz_i64small_gzip1_blank.fits BLANK int 37
for w in 2 1; do
  "$work/gen_fpack" rewidth "$work/i64small.fz" "$work/i64s_w$w.fz" $w
  "$work/gen_fpack" setkey "$work/i64s_w$w.fz" comp_fz_i64small_gzip1_w${w}_blank.fits BLANK int 37
done

# --- round 3 ---------------------------------------------------------------
# The edges of the ZBITPIX range in wider tiles (tile 1 values 0/1 set to
# the lowest and highest value the image type holds; cut to 2 rows): read.
# One past either edge: cfitsio's read as the image type fails (status 412)
# and funpack refuses the file, so the library refuses it (comp_bad).
edge() { # edge <case> <in.fz> <W> <i:v,...>
  "$work/gen_fpack" rewidth "$2" "$work/edge.fz" "$3" 0 0 "$4"
  "$work/gen_fpack" croprows "$work/edge.fz" "$1" 2
}
edge comp_fz_u8_gzip2_w2_edge.fits  comp_fz_u8_gzip2.fits  2 0:0,1:255
edge comp_fz_u16_gzip1_w4_edge.fits comp_fz_u16_gzip1.fits 4 0:-32768,1:32767
edge comp_fz_i32_gzip1_w8_edge.fits comp_fz_i32_gzip1.fits 8 0:-2147483648,1:2147483647
edge comp_bad/u8_gzip2_w2_lo.fits   comp_fz_u8_gzip2.fits  2 0:-1
edge comp_bad/u8_gzip2_w2_hi.fits   comp_fz_u8_gzip2.fits  2 0:256
edge comp_bad/u16_gzip1_w4_lo.fits  comp_fz_u16_gzip1.fits 4 0:-32769
edge comp_bad/u16_gzip1_w4_hi.fits  comp_fz_u16_gzip1.fits 4 0:32768
edge comp_bad/i32_gzip1_w8_lo.fits  comp_fz_i32_gzip1.fits 8 0:-2147483649
edge comp_bad/i32_gzip1_w8_hi.fits  comp_fz_i32_gzip1.fits 8 0:2147483648
# A null value the image type cannot hold: 4-byte tiles on a 16-bit image
# with stored 37232 in pixel 0 and BLANK = 37232. cfitsio flags it null
# (fits_read_imgnull), funpack and a read as the image type fail (412).
"$work/gen_fpack" rewidth comp_fz_i16blank_gzip1.fits "$work/no.fz" 4 0 70000
"$work/gen_fpack" setkey "$work/no.fz" "$work/no2.fz" BLANK int 37232
"$work/gen_fpack" croprows "$work/no2.fz" comp_bad/i16_gzip1_w4_nullout.fits 2
# BLANK / ZBLANK outside the C int cfitsio reads them as (ffgky TINT
# overflows): cfitsio then sets no null value (or falls back from ZBLANK to
# BLANK), astropy 8.0.2 refuses the file; refused (comp_bad).
"$work/gen_fpack" setkey "$work/i32s_r.fz" "$work/b32.fz" BLANK int 4294967296
"$work/gen_fpack" croprows "$work/b32.fz" comp_bad/i32_rice_blank2p32.fits 2
"$work/gen_fpack" setkey comp_fz_i64_gzip1.fits comp_bad/i64_gzip1_blank2p40.fits BLANK int 1099511627783
"$work/gen_fpack" setkey comp_fz_i16blank_rice.fits "$work/zb.fz" ZBLANK int 8589934592
"$work/gen_fpack" croprows "$work/zb.fz" comp_bad/i16_rice_zblank2p33.fits 2
# ... and in a per-tile 1K ZBLANK column: cfitsio fails the read (412).
"$work/gen_fpack" croprows comp_fz_i16blank_rice.fits "$work/zc.fz" 2
"$work/gen_fpack" zblankcol "$work/zc.fz" comp_bad/i16_rice_zblankcol1k.fits 1K 4294967296
# A ZBLANK on 1-byte RICE tiles of a quantized float: cfitsio casts it to
# the unsigned tile type ((unsigned char) 406 = 150; signed it would be
# -106), so the 150s are null.
"$work/gen_fpack" setkey comp_fz_f32_rice_bp1.fits comp_fz_f32_rice_bp1_zblank406.fits ZBLANK int 406
# HCOMPRESS_1 / PLIO_1 decode to 4-byte ints, which cfitsio compares with
# the int null value: BLANK 65536 on a 16-bit image flags nothing (cast to
# the 16-bit image type it would be 0, and the zeros would be null).
fpack -C -h -O "$work/i16s_h.fz" "$work/i16small.fits"
"$work/gen_fpack" setkey "$work/i16s_h.fz" comp_fz_i16small_hcomp_blank65536.fits BLANK int 65536
fpack -C -p -O "$work/i16s_p.fz" "$work/i16small.fits"
"$work/gen_fpack" setkey "$work/i16s_p.fz" comp_fz_i16small_plio_blank65536.fits BLANK int 65536
# Quantized float GZIP tiles of 1 and 2 bytes (cfitsio unquantize_i1r4 /
# unquantize_i2r4 / _i2r8 from the inflated length), and of 8 bytes, which
# cfitsio, funpack and astropy all return as zeros (refused, comp_bad).
fpack -C -g1 -q0 -0.4 -O "$work/f32q2g1.fz" "$work/f32.fits"
"$work/gen_fpack" rewidthq "$work/f32q2g1.fz" comp_fz_f32_gzip1_q_w1.fits 1
"$work/gen_fpack" rewidthq comp_fz_f32_gzip2_sd1.fits comp_fz_f32_gzip2_q_w2.fits 2
"$work/gen_fpack" rewidthq comp_fz_f64_gzip1_sd1.fits comp_fz_f64_gzip1_q_w2.fits 2
"$work/gen_fpack" rewidthq comp_fz_f32_gzip1_sd1.fits comp_bad/f32_gzip1_q_w8.fits 8
# ZQUANTIZ = 'NONE' on a -64 RICE_1 file with ZSCALE/ZZERO columns: cfitsio
# takes the 4-byte tile integers as IEEE float32 (fffr4r8) and scales them.
"$work/gen_fpack" setkey comp_fz_f64_rice_sd1.fits comp_fz_f64_rice_zqnone.fits ZQUANTIZ str NONE

# --- round 2 fix pass ------------------------------------------------------
# Lossy HCOMPRESS of images saturated at the type maximum (gen_fpack.c
# u8sat/i16sat/u16sat/i32sat, two 20x6 tiles: one 40x6 tile + BZERO crashes
# the cfref re-read): some reconstructed pixels lie above the maximum.
# 16-bit: cfitsio's read as the image type and funpack clip them (the twin
# is the truth; cfitsio's TDOUBLE read returns them unclipped). 8-bit:
# funpack fails (status 412), refused. 32-bit: fits_hdecompress64's (int)
# cast wraps them (cfitsio and funpack return the wrapped values), refused.
pk u16sat_hcomp_s2 u16sat -h -s 2 -t 20,6
pk i16sat_hcomp_s2 i16sat -h -s 2 -t 20,6
fpack -C -h -s 2 -t 20,6 -O comp_bad/u8sat_hcomp_s2.fits "$work/u8sat.fits"
fpack -C -h -s 2 -t 20,6 -O comp_bad/i32sat_hcomp_s2.fits "$work/i32sat.fits"
# PLIO_1 values the ZBITPIX type cannot hold: the uint16 PLIO file with
# BZERO 0 (values 32768 above int16), the 0..2^24-1 values of i32_plio as
# ZBITPIX 8, and as ZBITPIX 16 with BZERO 32768 (value - 32768 above int16).
# cfitsio's TDOUBLE read returns the PLIO values, its read as the image
# type and funpack fail (412): refused.
"$work/gen_fpack" setkey comp_fz_u16_plio.fits comp_bad/u16_plio_bzero0.fits BZERO int 0
"$work/gen_fpack" croprows comp_fz_i32_plio.fits "$work/p6.fz" 6
"$work/gen_fpack" setkey "$work/p6.fz" comp_bad/i32_plio_zbitpix8.fits ZBITPIX int 8
"$work/gen_fpack" setkey "$work/p6.fz" "$work/p16.fz" ZBITPIX int 16
"$work/gen_fpack" setkey "$work/p16.fz" comp_bad/i32_plio_zbitpix16_bzero32768.fits BZERO int 32768
# The int32 boundary of the null value. 4-byte GZIP_1 tiles of an int32
# image with pixel 0 = -2^31 and pixel 1 = 2^31 - 1 (2 rows): BLANK at
# either extreme is read (cfitsio flags that pixel); one past either is
# refused (cfitsio's C int read of the keyword overflows; astropy refuses).
edge "$work/i32e.fz" comp_fz_i32_gzip1.fits 4 0:-2147483648,1:2147483647
"$work/gen_fpack" setkey "$work/i32e.fz" comp_fz_i32_gzip1_blankmin.fits BLANK int -2147483648
"$work/gen_fpack" setkey "$work/i32e.fz" comp_fz_i32_gzip1_blankmax.fits BLANK int 2147483647
"$work/gen_fpack" setkey "$work/i32e.fz" comp_bad/i32_gzip1_blank2p31.fits BLANK int 2147483648
"$work/gen_fpack" setkey "$work/i32e.fz" comp_bad/i32_gzip1_blankm2p31m1.fits BLANK int -2147483649
# ... and in a per-tile 1K ZBLANK column: -2^31 is read, 2^31 and
# -2^31 - 1 are refused (cfitsio fails the read, 412).
"$work/gen_fpack" zblankcol comp_fz_i32_gzip1_blankmin.fits comp_fz_i32_gzip1_zblankcolmin.fits 1K -2147483648
"$work/gen_fpack" zblankcol comp_fz_i32_gzip1_blankmin.fits comp_bad/i32_gzip1_zblankcol2p31.fits 1K 2147483648
"$work/gen_fpack" zblankcol comp_fz_i32_gzip1_blankmin.fits comp_bad/i32_gzip1_zblankcolm2p31m1.fits 1K -2147483649
# A ZBLANK column wins over a BLANK keyword, which cfitsio then does not
# use: BLANK 2^32 next to the column is read.
"$work/gen_fpack" setkey comp_fz_i16_zblankcol.fits comp_fz_i16_zblankcol_blank2p32.fits BLANK int 4294967296
# Two null values the image type cannot hold, in two tiles (4-byte tiles
# of an int16 image, per-tile ZBLANK 70000 and 70001 at pixel 0 of tile 1
# and of tile 2): the unmasked read names the first.
"$work/gen_fpack" rewidth comp_fz_i16blank_gzip1.fits "$work/n2.fz" 4 0 0 0:70000,2@0:70001
"$work/gen_fpack" setkey "$work/n2.fz" "$work/n2b.fz" BLANK int 70000
"$work/gen_fpack" setkey "$work/n2b.fz" "$work/n2c.fz" ZBLKSTEP int 1
"$work/gen_fpack" croprows "$work/n2c.fz" "$work/n2d.fz" 2
"$work/gen_fpack" zblankcol "$work/n2d.fz" comp_bad/i16_gzip1_w4_nullout2.fits
# A null value that HCOMPRESS_1 / PLIO_1 tiles hold: BLANK 0 on the
# lossless HCOMPRESS i16small (its zeros) and BLANK 3 on the PLIO mask (its
# 3s); cfitsio flags them (compared with the decoded C ints).
"$work/gen_fpack" setkey "$work/i16s_h.fz" comp_fz_i16small_hcomp_blank0.fits BLANK int 0
"$work/gen_fpack" setkey comp_fz_i16mask_plio.fits comp_fz_i16mask_plio_blank3.fits BLANK int 3
# Header values the library refuses with an error: a ZBITPIX that is no
# BITPIX (cfitsio reads the tiles anyway) and negative ZNAXISn with a
# matching tile count (cfitsio: status 414).
"$work/gen_fpack" setkey comp_fz_u16_rice.fits comp_bad/u16_rice_zbitpix12.fits ZBITPIX int 12
"$work/gen_fpack" setkey comp_fz_u16_rice.fits comp_bad/u16_rice_zbitpix0.fits ZBITPIX int 0
# (cfitsio will not write the negative sizes: the cards are byte-patched)
"$work/gen_fpack" patch comp_fz_u16_rice.fits "$work/ng1.fz" "ZNAXIS1 =                   40" "ZNAXIS1 =                   -1"
"$work/gen_fpack" patch "$work/ng1.fz" "$work/ng2.fz" "ZTILE1  =                   40" "ZTILE1  =                    1"
"$work/gen_fpack" patch "$work/ng2.fz" comp_bad/u16_rice_znaxisneg.fits "ZNAXIS2 =                   24" "ZNAXIS2 =                  -24"
# ZBITPIX 64 in RICE_1, HCOMPRESS_1 and PLIO_1 (gen_fpack mk64; fpack and
# fits_img_compress refuse 64-bit images): cfitsio decodes them through
# C ints; the twin is astropy 8.0.2's decode. HCOMPRESS_1 with a pixel one
# past either end of the int32 range, which cfitsio's (int) cast wraps:
# refused.
for a in rice:RICE_1 hcomp:HCOMPRESS_1 plio:PLIO_1; do
  "$work/gen_fpack" mk64 "comp_fz_i64_${a%%:*}.fits" "${a#*:}"
done
"$work/gen_fpack" mk64 comp_bad/i64_hcomp_2p31.fits HCOMPRESS_1 2147483648
"$work/gen_fpack" mk64 comp_bad/i64_hcomp_m2p31m1.fits HCOMPRESS_1 -2147483649

# --- real files: astropy's test data (BSD-3-Clause), copied verbatim -----
apdir=$("$py" -c 'import astropy, os; print(os.path.dirname(astropy.__file__))')/io/fits
cp "$apdir/hdu/compressed/tests/data/compressed_with_nan.fits" comp_fz_real_with_nan.fits
# real_with_nan: fpack 2022: -64 SD1, 6x6 tiles, NaN
cp "$apdir/tests/data/compressed_image.fits" comp_fz_real_image.fits
# real_image: 10x10 int16 RICE_1 with ZTENSION

# M13 (a 2006 SkyView 300x300 int16 image) as astropy ships it in four
# compressions. m13_gzip holds its 16-bit pixels in 4-byte GZIP_1 tiles and
# m13_rice is RICE_1 without ZVAL2 (BYTEPIX 4), the layouts old cfitsio wrote.
# m13_hcomp (one 300x300 tile) is copied verbatim; the row-tiled gzip, rice
# and plio files are cut to their first 40 tiles (rows) byte for byte
# (gen_fpack croprows): the kept tiles are the original 2006 bytes.
cp "$apdir/hdu/compressed/tests/data/m13_hcomp.fits" comp_fz_real_m13_hcomp.fits
for m in gzip rice plio; do
  "$work/gen_fpack" croprows "$apdir/hdu/compressed/tests/data/m13_$m.fits" "comp_fz_real_m13_$m.fits" 40
done

# The decompressed twins (gen_twins.sh): funpack, or astropy at ZBITPIX 64.
# A twin byte-identical to an earlier one by name is listed in
# comp_fzun_same.txt as "<removed> <kept>"; the corpus test reads a case's
# twin golden through that list.
sh "$here/gen_twins.sh" "$here/cfref/twins.out"
cp "$here/cfref/twins.out/comp_fzun_same.txt" comp_fzun_same.txt

ls comp_fz_*.fits | wc -l | sed 's/^/comp_fz files: /'
