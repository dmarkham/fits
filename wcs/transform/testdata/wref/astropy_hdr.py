#!/usr/bin/env python3
"""astropy.wcs's reading of a header case file (the hdr2sky case format).

Usage: astropy_hdr.py <case file>

Reads the case exactly as hdr2sky.c does (GRID / PIX directives before
HEADER, '#' lines skipped; then one card per line, padded to 80 columns,
stopping at END), builds astropy.wcs.WCS from the cards (warnings
recorded) and prints one JSON object:

    {"astropy": ..., "wcslib": ..., "numpy": ...,
     "warnings": [sorted distinct warning messages],
     "log": [sorted distinct astropy log messages, whitespace collapsed],
     "sip": true if astropy kept a SIP polynomial,
     "p2s": [[x, y, lng, lat], ...],   all_pix2world(x, y, 1)
     "s2p": [[lng, lat, x, y], ...]}   all_world2pix of each p2s output, 1

Pixels are 1-based FITS coordinates (origin 1), as in hdr2sky; every float
is Python's repr (the shortest decimal that reads back to the same
binary64). Used by TestHdrDeliberateDifferences (a recorded tool reading,
testutil.AstropyScript), with the venv interpreter of FITS_ASTROPY_PYTHON
(astropy 8.0.2, numpy 2.4.6).
"""

import contextlib
import io
import json
import sys
import warnings

import numpy as np
import astropy
import astropy.wcs
from astropy.io import fits


def lin(v0, v1, n, i):
    # hdr2sky.c lin(): the same operations in the same order.
    if n == 1:
        return v0
    return v0 + (v1 - v0) * i / (n - 1)


def read_case(path):
    xs, ys, cards = [], [], []
    in_hdr = False
    for line in open(path, encoding="latin-1").read().split("\n"):
        line = line.split("\r")[0]
        if in_hdr:
            if line == "":
                continue
            card = line.ljust(80)
            if card.startswith("END     "):
                break
            cards.append(card)
            continue
        p = line.lstrip(" \t")
        if p == "" or p[0] == "#":
            continue
        f = p.split()
        if p.startswith("HEADER"):
            in_hdr = True
        elif f[0] == "GRID":
            v = [float(x) for x in f[1:7]]
            na, nb = int(v[2]), int(v[5])
            for j in range(nb):
                for i in range(na):
                    xs.append(lin(v[0], v[1], na, i))
                    ys.append(lin(v[3], v[4], nb, j))
        elif f[0] == "PIX":
            xs.append(float(f[1]))
            ys.append(float(f[2]))
        elif f[0] in ("SKYGRID", "SKY", "ALT"):
            sys.exit("astropy_hdr.py: %s is not supported" % f[0])
    return np.array(xs), np.array(ys), cards


def main():
    xs, ys, cards = read_case(sys.argv[1])
    hdr = fits.Header.fromstring("".join(cards))
    # The log messages are recorded, not printed (astropy's console handler
    # writes INFO to stdout, which holds the JSON).
    with warnings.catch_warnings(record=True) as caught, astropy.log.log_to_list() as logged, \
            contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        warnings.simplefilter("always")
        w = astropy.wcs.WCS(hdr)
    msgs = sorted({str(c.message).strip() for c in caught})
    log = sorted({" ".join(r.getMessage().split()) for r in logged})
    lng, lat = w.all_pix2world(xs, ys, 1)
    px, py = w.all_world2pix(lng, lat, 1)
    out = {
        "astropy": astropy.__version__,
        "wcslib": astropy.wcs.WCSLIB_VERSION if hasattr(astropy.wcs, "WCSLIB_VERSION") else astropy.wcs._wcs.WCSLIB_VERSION,
        "numpy": np.__version__,
        "warnings": msgs,
        "log": log,
        "sip": w.sip is not None,
        "p2s": [[float(a), float(b), float(c), float(d)] for a, b, c, d in zip(xs, ys, lng, lat)],
        "s2p": [[float(a), float(b), float(c), float(d)] for a, b, c, d in zip(lng, lat, px, py)],
    }
    json.dump(out, sys.stdout, indent=0)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
