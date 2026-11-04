#!/bin/sh
# fpack (cfitsio 4.7.0 / fpack 1.7.0, default RICE_1 lossless for integers)
# the samples whose compressed form is tested. -C: no CHECKSUM keywords
# (they carry a timestamp and would make the output irreproducible).
# zblank_variants.py then derives the ZBLANK samples (stdlib python).
set -e
cd "$(dirname "$0")"
for f in u16_blank i32_blank_bscale u8_blank asi533_u16_crop i16_bscale_neg i16_bscale0 i8_signed u16_rgb card_i16_bscale_dup_2d; do
  rm -f "comp_$f.fits.fz"
  fpack -C -O "comp_$f.fits.fz" "$f.fits"
done
${PYTHON:-python3} zblank_variants.py
