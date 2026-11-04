"""Write strarray_astropy.fits: string-array columns written by astropy.

Usage: python gen_strarray.py <out.fits>

astropy 8.0.2 writes a 2-D string column as TFORM rA with TDIM (w, n)
and pads each element with NUL bytes, not blanks: the cell of
['ab', 'cdefg'] is 'ab\\0\\0\\0cdefg'. A reader that treats the cell as one
string cut at its first NUL keeps only 'ab'. The values come from
cfitsio's reading of the file (golden/strarray_astropy.fits.json) and
astropy's own (golden_astropy/), not from this script.
"""

import sys

import numpy as np
from astropy.io import fits


def main(out):
    s = np.array([["ab", "cdefg"], ["hello", "world"], ["", "x"]], dtype="S5")
    t = np.array([[["a", "bb"], ["ccc", ""]], [["d", "e"], ["f", "ggg"]], [["", ""], ["", "z"]]], dtype="S3")
    one = np.array(["single", "", "x"], dtype="S6")
    cols = [
        fits.Column(name="S", format="10A", dim="(5,2)", array=s),
        fits.Column(name="T", format="12A", dim="(3,2,2)", array=t),
        fits.Column(name="ONE", format="6A", array=one),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="STRARRAY")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out, overwrite=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: gen_strarray.py <out.fits>")
    main(sys.argv[1])
