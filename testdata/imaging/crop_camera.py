#!/usr/bin/env python3
"""Cut the real camera sample asi533_u16_crop.fits out of a full frame.

    <venv>/bin/python crop_camera.py /path/to/Light_001.fits

Source: a 60 s light frame from a ZWO ASI533MC Pro (Bayer RGGB, 3008x3008,
BITPIX=16 BZERO=32768) written by INDI, owned by the repository author (not
committed: 18 MB). The 128x128 crop centred on the target star is copied
byte for byte from the data unit (no re-encoding); the header keeps every
card except four that, together with DATE-OBS, locate the observing site
(SITELAT, SITELONG, OBJCTAZ, OBJCTALT, AIRMASS removed), NAXIS1/NAXIS2 are
set to 128 and CRPIX1/CRPIX2 are shifted by the crop origin; HISTORY cards
record the crop. Re-running on the same source reproduces the file.
"""
import hashlib
import os
import sys

SRC_SHA256 = "2f4ae03a13ab74f7c6c6494ad8de9d187bd4e645b86f54ce566b555aeee5922d"
X0, Y0, N = 1440, 1440, 128   # 0-based crop origin (column, row) and size
DROP = {"SITELAT", "SITELONG", "OBJCTAZ", "OBJCTALT", "AIRMASS"}

src = sys.argv[1]
raw = open(src, "rb").read()
digest = hashlib.sha256(raw).hexdigest()
if SRC_SHA256 != "SET_BELOW" and digest != SRC_SHA256:
    sys.exit(f"source sha256 {digest} != {SRC_SHA256}")

cards = []
off = 0
while True:
    c = raw[off:off + 80].decode("ascii")
    off += 80
    if c.startswith("END     ") or c.rstrip() == "END":
        break
    cards.append(c)
data_start = (off + 2879) // 2880 * 2880


def val(key):
    for c in cards:
        if c[:8].rstrip() == key:
            return c[10:].split("/")[0].strip()
    raise KeyError(key)


W, H = int(val("NAXIS1")), int(val("NAXIS2"))
assert val("BITPIX") == "16" and W == 3008 and H == 3008
out = []
for c in cards:
    k = c[:8].rstrip()
    if k in DROP:
        continue
    if k in ("NAXIS1", "NAXIS2"):
        c = f"{k:<8}= {N:>20d} / length of data axis {k[-1]}".ljust(80)
    elif k in ("CRPIX1", "CRPIX2"):
        shift = X0 if k == "CRPIX1" else Y0
        v = float(val(k)) - shift
        c = f"{k:<8}= {v:>20.10E} / {k} (shifted by crop origin)".ljust(80)
    out.append(c)
out.append(f"HISTORY cropped {N}x{N} at 0-based (x={X0}, y={Y0}) from a {W}x{H} frame".ljust(80))
out.append(f"HISTORY source sha256 {digest[:32]}".ljust(80))
out.append(f"HISTORY   {digest[32:]}".ljust(80))
out.append("HISTORY data copied verbatim; removed SITELAT SITELONG OBJCTAZ OBJCTALT AIRMASS".ljust(80))
out.append("END".ljust(80))
hdr = "".join(out).encode("ascii")
assert all(len(c) == 80 for c in out)
hdr += b" " * (-len(hdr) % 2880)
rows = []
for y in range(Y0, Y0 + N):
    a = data_start + 2 * (y * W + X0)
    rows.append(raw[a:a + 2 * N])
data = b"".join(rows)
data += b"\0" * (-len(data) % 2880)
dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), "asi533_u16_crop.fits")
with open(dst, "wb") as f:
    f.write(hdr + data)
print("wrote", dst, "source sha256", digest)
