#!/bin/sh
# Write the decompressed twin of every compressed read corpus sample:
#
#   sh testdata/gen_twins.sh <outdir>
#
# testdata/comp_fz_<case>.fits -> <outdir>/comp_fzun_<case>.fits, written
# by /usr/bin/funpack (cfitsio's fits_img_decompress, an independent decode
# path), or for a ZBITPIX 64 sample, which funpack cannot decompress, by
# astropy 8.0.2 (its decode written as a plain image). The twins are not
# committed: they are derived from the committed samples, and what the
# tests use of them is their cfitsio golden,
# testdata/cfref/golden/twins/comp_fzun_<case>.fits.json (written by
# `make -C testdata/cfref goldens`, which runs this script into
# testdata/cfref/twins.out). A twin byte-identical to an earlier one (by
# name) is removed and listed in <outdir>/comp_fzun_same.txt as
# "<removed> <kept>"; gen_fpack.sh commits that list as
# testdata/comp_fzun_same.txt. funpack runs with -C (no CHECKSUM), so the
# twins reproduce byte for byte. Tool versions: testdata/COMP_FZ_PROVENANCE.md.
set -eu

here=$(cd "$(dirname "$0")" && pwd)
out=$1
py=${FITS_ASTROPY_PYTHON:-python3}
mkdir -p "$out"
rm -f "$out"/comp_fzun_*.fits "$out"/comp_fzun_same.txt "$out"/.twinsums

for f in "$here"/comp_fz_*.fits; do
  c=$(basename "$f" .fits)
  c=${c#comp_fz_}
  if head -c 5760 "$f" | fold -w 80 | grep -q '^ZBITPIX = *64 '; then
    "$py" - "$f" "$out/comp_fzun_$c.fits" <<'PY'
import sys

import astropy.io.fits as fits

with fits.open(sys.argv[1]) as hl:
    d = hl[1].data.copy()
fits.HDUList([fits.PrimaryHDU(data=d)]).writeto(sys.argv[2], overwrite=True)
PY
  else
    funpack -C -O "$out/comp_fzun_$c.fits" "$f"
  fi
done

: > "$out/comp_fzun_same.txt"
for f in $(cd "$out" && ls comp_fzun_*.fits | LC_ALL=C sort); do
  h=$(sha256sum "$out/$f" | cut -c1-64)
  k=$(grep "^$h " "$out/.twinsums" 2>/dev/null | cut -d' ' -f2 || true)
  if [ -n "$k" ]; then
    echo "$f $k" >> "$out/comp_fzun_same.txt"
    rm -f "$out/$f"
  else
    echo "$h $f" >> "$out/.twinsums"
  fi
done
rm -f "$out/.twinsums"
