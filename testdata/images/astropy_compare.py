#!/usr/bin/env python3
"""Second-reader check used by TestWriteImageAstropyReadsBack.

stdin: JSON list of {"name", "go", "ref"}: a FITS file the Go library wrote
and the cfitsio-written reference for the same case. For each pair, astropy
(default uint=True, so BZERO-offset integers come back as unsigned / int8)
reads both primary images; they must have the same dtype, the same values
(NaN equal to NaN, -0.0 distinct from 0.0) and the same BITPIX/BSCALE/BZERO
/BLANK values, an absent BSCALE / BZERO counting as its FITS default 1 / 0
(cfitsio's writer always writes both once either is set; a caller header
may carry only one). stdout: JSON {name: {"equal": bool, "detail": str}}.
"""
import json
import sys
import warnings

import numpy as np
from astropy.io import fits


def load(path):
    with warnings.catch_warnings():
        warnings.simplefilter("error")  # any astropy complaint is a failure
        with fits.open(path) as hl:
            h = hl[0]
            hdr = {k: h.header.get(k) for k in ("BITPIX", "BSCALE", "BZERO", "BLANK")}
            hdr["BSCALE"] = 1.0 if hdr["BSCALE"] is None else hdr["BSCALE"]
            hdr["BZERO"] = 0.0 if hdr["BZERO"] is None else hdr["BZERO"]
            return np.array(h.data), hdr, len(hl)


def same(a, b):
    if a.dtype != b.dtype or a.shape != b.shape:
        return False
    if a.dtype.kind == "f":
        return bool(np.all((a == b) & (np.signbit(a) == np.signbit(b)) | (np.isnan(a) & np.isnan(b))))
    return bool(np.array_equal(a, b))


out = {}
for p in json.load(sys.stdin):
    try:
        ga, gh, gn = load(p["go"])
        ra, rh, rn = load(p["ref"])
        ok = same(ga, ra) and gh == rh and gn == rn == 1
        detail = "" if ok else f"go {ga.dtype}{ga.shape} {gh} {ga.ravel()[:8]} / cfitsio {ra.dtype}{ra.shape} {rh} {ra.ravel()[:8]}"
    except Exception as e:  # reported, never swallowed
        ok, detail = False, f"{type(e).__name__}: {e}"
    out[p["name"]] = {"equal": ok, "detail": detail}
json.dump(out, sys.stdout)
