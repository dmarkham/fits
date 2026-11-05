"""Write the VLA heap boundary samples in this directory.

Usage: python gen_heap.py <outdir>

Each file is a primary HDU plus a one-column 1PJ(2) BINTABLE of one row
(NAXIS1 = 8) with an 8-byte heap holding the int32 values 7 and 9
(PCOUNT = 8), laid out by hand (plain Python, no libraries) so that the
descriptor or THEAP sits exactly on, or one byte past, a boundary of
FITS 4.0 7.3.5: the heap starts THEAP bytes into the data unit (default
NAXIS1*NAXIS2) and ends at NAXIS1*NAXIS2 + PCOUNT, so it holds
PCOUNT - (THEAP - NAXIS1*NAXIS2) bytes.

  end_exact            descriptor (2, 0): the payload ends at the heap end
  end_plus1            descriptor (2, 1): one byte past the heap end
  off_eq               descriptor (1, 8): offset equal to the heap size
  theap_end            THEAP = 16 = NAXIS1*NAXIS2 + PCOUNT, an empty row
  theap_end_plus1      THEAP = 17, one past the data unit, an empty row
  theap_fixed_minus1   THEAP = 7, inside the fixed rows
  theap_gap_exact      THEAP = 12 (a 4-byte gap), descriptor (1, 0): ends
                       at the heap end
  theap_gap_plus1      THEAP = 12, descriptor (1, 1): one byte past it
  neg_off              descriptor (1, -4): a negative offset (int32)
  neg_off_gap          THEAP = 12, descriptor (1, -4): a negative offset
                       that would land in the 4-byte gap before the heap
  theap_float          THEAP = 12.0 (a float literal, not an integer),
                       descriptor (1, 0)

fitsverify's verdict on each file is committed in fitsverify.txt
(`make heap-goldens`); TestVLAHeapBoundaries reads both.
"""

import os
import struct
import sys


def card(s):
    return s.ljust(80).encode("ascii")


def block(b, fill):
    return b + fill * ((-len(b)) % 2880)


def write(path, desc, theap=None):
    prim = b"".join(card(c) for c in (
        "SIMPLE  =                    T", "BITPIX  =                    8",
        "NAXIS   =                    0", "EXTEND  =                    T", "END"))
    ext = ["XTENSION= 'BINTABLE'", "BITPIX  =                    8",
           "NAXIS   =                    2", "NAXIS1  =                    8",
           "NAXIS2  =                    1", "PCOUNT  =                    8",
           "GCOUNT  =                    1", "TFIELDS =                    1",
           "TTYPE1  = 'V       '", "TFORM1  = '1PJ(2)  '"]
    if theap is not None:
        ext.append("THEAP   = %20s" % theap)
    ext.append("END")
    data = struct.pack(">ii", *desc) + struct.pack(">ii", 7, 9)
    with open(path, "wb") as f:
        f.write(block(prim, b" ") + block(b"".join(card(c) for c in ext), b" ") + block(data, b"\0"))


CASES = [
    ("end_exact", (2, 0), None),
    ("end_plus1", (2, 1), None),
    ("off_eq", (1, 8), None),
    ("theap_end", (0, 0), 16),
    ("theap_end_plus1", (0, 0), 17),
    ("theap_fixed_minus1", (1, 0), 7),
    ("theap_gap_exact", (1, 0), 12),
    ("theap_gap_plus1", (1, 1), 12),
    ("neg_off", (1, -4), None),
    ("neg_off_gap", (1, -4), 12),
    ("theap_float", (1, 0), "12.0"),
]

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: gen_heap.py <outdir>")
    for name, desc, theap in CASES:
        write(os.path.join(sys.argv[1], "heap_%s.fits" % name), desc, theap)
