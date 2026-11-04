"""Write gaia_dr3_sample.fits: a small Gaia-DR3-style source table.

Usage: python gen_gaia.py <out.fits>

The table is written by astropy (version pinned in PROVENANCE.md) so that
the Go reader is checked against a file it did not write, in the layout a
real catalogue extract has: 64-bit source_id identifiers that do not fit a
float64 mantissa, float64 coordinates, NaN parallaxes, float32 photometry,
a masked (TNULL) int16, a logical flag, a string column and the unsigned
16/32-bit TZERO convention astropy uses for uint16/uint32 columns.

The values are synthetic (Gaia-style, not copied from the archive); the
source_id values use the real DR3 encoding (HEALPix level-12 index * 2^35 +
running number), so they span the full positive int64 range like real ones.
"""

import sys

import numpy as np
from astropy.io import fits
from astropy.table import MaskedColumn, Table


def main(out):
    healpix = np.array([0, 1, 12345, 98765, 201326591], dtype=np.int64)
    running = np.array([1, 2**35 - 1, 77, 4096, 123456789], dtype=np.int64)
    source_id = healpix * (2**35) + running
    # Extremes a reader must keep exact.
    source_id = np.concatenate([source_id, np.array([2**62 + 1, 2**63 - 1], dtype=np.int64)])
    n = len(source_id)

    t = Table()
    t["source_id"] = source_id
    t["ra"] = np.linspace(0.0, 359.9999999999, n)
    t["dec"] = np.linspace(-89.999999, 89.999999, n)
    plx = np.linspace(-1.5, 25.0, n)
    plx[[2, 5]] = np.nan
    t["parallax"] = plx
    t["phot_g_mean_mag"] = np.linspace(3.0, 21.0, n).astype(np.float32)
    t["phot_variable_flag"] = np.array(
        ["NOT_AVAILABLE", "VARIABLE", "", "NOT_AVAILABLE", "VARIABLE", "X", "NOT_AVAILABLE"]
    )
    t["has_xp"] = np.array([True, False, True, True, False, False, True])
    t["visibility_periods_used"] = MaskedColumn(
        np.array([10, 0, 25, 3, 7, 9, 31], dtype=np.int16),
        mask=[False, True, False, False, True, False, False],
    )
    t["nu_flags"] = np.array([0, 1, 32767, 32768, 40000, 65534, 65535], dtype=np.uint16)
    t["solution_id"] = np.array(
        [0, 1, 2**31 - 1, 2**31, 3000000000, 2**32 - 2, 2**32 - 1], dtype=np.uint32
    )
    hdu = fits.BinTableHDU(t, name="GAIA_DR3_STYLE")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out, overwrite=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: gen_gaia.py <out.fits>")
    main(sys.argv[1])
