"""Write astropy_special.fits: E, D, C and M cells (fixed and VLA) holding
the IEEE values cfitsio's null-checking reads flush or flag: +-Inf, -0,
the smallest subnormals of float32 and float64, NaN, the largest finite
values.

Usage: python gen_special.py <out.fits>

Written by astropy 8.0.2 itself (fits.BinTableHDU.from_columns), so the
writer is independent of the Go code. The sample lives outside
testdata/tables/*.fits because the shared cfitsio harness cannot dump a
complex cell that cfitsio leaves unset (its null-checking complex read
does that for Inf/NaN), so the truth is astropy's own read of the file
(golden_astropy/astropy_special.fits.json, `astropy_dump.py --golden`).
TestSpecialValuesSample also checks that this golden holds exactly the
bit patterns below, so the sample cannot lose them silently.
"""

import sys

import numpy as np
from astropy.io import fits

INF = float("inf")
NAN = float("nan")
SUB32 = float(np.float32(1.401298464324817e-45))  # smallest float32 subnormal
SUB64 = 5e-324  # smallest float64 subnormal
MAX32 = float(np.finfo(np.float32).max)
MAX64 = float(np.finfo(np.float64).max)

# One value per row; the same list for the float32 and float64 columns
# (SUB64 underflows to 0 in float32, so the float32 columns use SUB32).
F32 = [INF, -INF, -0.0, SUB32, -SUB32, NAN, 1.5, MAX32]
F64 = [INF, -INF, -0.0, SUB64, -SUB64, NAN, 1.5, MAX64]


def main(out):
    n = len(F32)
    e = np.array(F32, dtype=">f4")
    d = np.array(F64, dtype=">f8")
    e2 = np.array([[F32[r], F32[(r + 3) % n]] for r in range(n)], dtype=">f4")
    c = np.array([complex(F32[r], F32[(r + 2) % n]) for r in range(n)], dtype=">c8")
    m = np.array([complex(F64[r], F64[(r + 2) % n]) for r in range(n)], dtype=">c16")
    ve = [np.array(F32[: r % 4], dtype=">f4") for r in range(n)]
    vd = [np.array(F64[r:], dtype=">f8") for r in range(n)]
    vc = [np.array([complex(F32[k], F32[(k + 5) % n]) for k in range(r % 3)], dtype=">c8") for r in range(n)]
    vm = [np.array([complex(F64[k], F64[(k + 1) % n]) for k in range(r, min(n, r + 2))], dtype=">c16") for r in range(n)]
    cols = [
        fits.Column(name="E1", format="1E", array=e),
        fits.Column(name="D1", format="1D", array=d),
        fits.Column(name="E2", format="2E", array=e2),
        fits.Column(name="C1", format="1C", array=c),
        fits.Column(name="M1", format="1M", array=m),
        fits.Column(name="VE", format="PE(3)", array=ve),
        fits.Column(name="VD", format="QD(8)", array=vd),
        fits.Column(name="VC", format="PC(2)", array=vc),
        fits.Column(name="VM", format="QM(2)", array=vm),
    ]
    hdu = fits.BinTableHDU.from_columns(cols, name="SPECIAL")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(out, overwrite=True)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: gen_special.py <out.fits>")
    main(sys.argv[1])
