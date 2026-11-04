"""Write vla_theap_gap.fits: a 1PJ / 1QD VLA table whose heap does not start
right after the fixed rows (THEAP > NAXIS1*NAXIS2, PCOUNT = gap + heap).

Usage: python gen_theap.py <out.fits>

The bytes are laid out by hand (header cards via astropy.io.fits.Header,
data with numpy big-endian dtypes) because neither cfitsio nor astropy
writes a heap gap on its own. The file is checked by fitsverify (0
warnings, 0 errors) and its values come from cfitsio's reading of it
(golden/vla_theap_gap.fits.json), not from this script.

Heap layout (offsets relative to THEAP): the rows' payloads are stored in
REVERSE row order and row 3 re-uses row 1's payload, so a reader that
assumes "heap = rows in order right after the table" gets wrong values.
"""

import sys

import numpy as np
from astropy.io import fits

NROWS = 4
THEAP = 2880 + 16  # fixed rows are 4 * 28 = 112 bytes; the heap starts 2784 bytes later


def main(out):
    ids = np.array([10, 20, 30, 40], dtype=">i4")
    jrows = [[1, -1], [2, -2, 4], [], [7]]
    drows = [[0.5], [], [1.25, -2.5], [3.0, 4.0, 5.0]]

    heap = bytearray()
    jdesc = [None] * NROWS
    ddesc = [None] * NROWS
    for r in reversed(range(NROWS)):
        if r == 2:
            continue
        b = np.array(jrows[r], dtype=">i4").tobytes()
        jdesc[r] = (len(jrows[r]), len(heap))
        heap += b
    jdesc[2] = (0, 0)
    for r in reversed(range(NROWS)):
        if r == 1:
            ddesc[r] = (0, 0)
            continue
        b = np.array(drows[r], dtype=">f8").tobytes()
        ddesc[r] = (len(drows[r]), len(heap))
        heap += b
    # Row 3's J payload aliases row 1's (shared heap region).
    jrows_alias = list(jdesc)
    jrows_alias[2] = (1, jdesc[0][1])  # first element of row 0's payload
    jdesc = jrows_alias

    rowdt = np.dtype([("ID", ">i4"), ("JARR", ">i4", (2,)), ("DARR", ">i8", (2,))])
    rows = np.zeros(NROWS, dtype=rowdt)
    for r in range(NROWS):
        rows["ID"][r] = ids[r]
        rows["JARR"][r] = jdesc[r]
        rows["DARR"][r] = ddesc[r]
    table = rows.tobytes()
    naxis1 = rowdt.itemsize
    assert naxis1 == 28
    gap = THEAP - naxis1 * NROWS
    assert gap > 0
    data = table + b"\0" * gap + bytes(heap)
    pcount = gap + len(heap)

    h = fits.Header()
    h["XTENSION"] = ("BINTABLE", "binary table extension")
    h["BITPIX"] = 8
    h["NAXIS"] = 2
    h["NAXIS1"] = naxis1
    h["NAXIS2"] = NROWS
    h["PCOUNT"] = (pcount, "heap gap + heap size")
    h["GCOUNT"] = 1
    h["TFIELDS"] = 3
    h["TTYPE1"] = "ID"
    h["TFORM1"] = "1J"
    h["TTYPE2"] = "JARR"
    h["TFORM2"] = "1PJ(3)"
    h["TTYPE3"] = "DARR"
    h["TFORM3"] = "1QD(3)"
    h["THEAP"] = (THEAP, "heap starts after a gap")
    h["EXTNAME"] = "THEAP_GAP"

    prim = fits.Header()
    prim["SIMPLE"] = True
    prim["BITPIX"] = 8
    prim["NAXIS"] = 0
    prim["EXTEND"] = True

    def block(b, fill):
        return b + fill * ((2880 - len(b) % 2880) % 2880)

    with open(out, "wb") as f:
        f.write(prim.tostring().encode("ascii"))
        f.write(h.tostring().encode("ascii"))
        f.write(block(data, b"\0"))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: gen_theap.py <out.fits>")
    main(sys.argv[1])
