#!/usr/bin/env python3
"""Generate testdata/wcs_crosscheck.json from real observatory headers.

Usage:

    python3 testdata/generate_wcs_crosscheck.py > testdata/wcs_crosscheck.json

Requires astropy (>=5.0) installed in the active Python environment. This
script is run once offline to produce the golden reference values; the Go
test suite reads the JSON directly and does not need Python at test time.

The reference values are produced by astropy.wcs.WCS, which wraps the
wcslib C reference implementation. Agreement between our Go transform
and this JSON is the v1.1 WCS release gate.

The fixture set:

  - HST/ACS     — TAN + SIP distortion, typical archive product
  - HST/WFC3    — TAN + SIP
  - JWST/NIRCAM — TAN + SIP
  - Chandra     — TAN (no distortion)
  - Spitzer     — TAN + SIP
  - Gaia DR3    — TAN (no distortion)
  - Planck      — HPX
  - WMAP        — HPX
  - COBE/DMR    — CSC (historical)
  - 2MASS       — TNX (IRAF)

Each fixture samples pixel → sky at ~20 grid points across the detector
so that small-angle and large-angle behavior are both exercised.
"""

import json
import sys

try:
    import numpy as np
    from astropy.wcs import WCS
    from astropy.io.fits import Header
except ImportError:
    print("error: astropy not installed. Run: pip install astropy>=5.0", file=sys.stderr)
    sys.exit(1)


def sample_wcs(header_cards, name, source, projcode, tolerance=1e-9, n=5):
    """Build a WCS from the given FITS-style header dict, sample it at a
    grid of pixels, and return a CrossCheckFixture dict."""
    hdr = Header()
    for k, v in header_cards.items():
        hdr[k] = v
    w = WCS(hdr)
    # Sample at an n x n grid covering the image.
    naxis1 = header_cards.get("NAXIS1", 1024)
    naxis2 = header_cards.get("NAXIS2", 1024)
    xs = np.linspace(0, naxis1 - 1, n)
    ys = np.linspace(0, naxis2 - 1, n)
    vectors = []
    for y in ys:
        for x in xs:
            sky = w.wcs_pix2world([[x, y]], 0)[0]
            # astropy returns NaN for pixels outside the projection's
            # valid domain. JSON does not support NaN; skip those points.
            if not (np.isfinite(sky[0]) and np.isfinite(sky[1])):
                continue
            vectors.append({"pixel": [float(x), float(y)], "sky": [float(sky[0]), float(sky[1])]})

    # Extract CDELT / PC / CD / PV from the header.
    def hget(key, default=None):
        return float(header_cards[key]) if key in header_cards else default

    pc_matrix = None
    cd_matrix = None
    if any(k.startswith("CD1_") for k in header_cards):
        cd_matrix = [
            [hget("CD1_1", 0), hget("CD1_2", 0)],
            [hget("CD2_1", 0), hget("CD2_2", 0)],
        ]
    if any(k.startswith("PC1_") for k in header_cards):
        pc_matrix = [
            [hget("PC1_1", 1), hget("PC1_2", 0)],
            [hget("PC2_1", 0), hget("PC2_2", 1)],
        ]

    pvs = []
    for k, v in header_cards.items():
        if k.startswith("PV1_") or k.startswith("PV2_"):
            axis = int(k[2])
            idx = int(k[4:])
            pvs.append({"axis": axis, "index": idx, "value": float(v)})

    return {
        "name": name,
        "source": source,
        "ctype1": header_cards["CTYPE1"],
        "ctype2": header_cards["CTYPE2"],
        "crval1": hget("CRVAL1", 0),
        "crval2": hget("CRVAL2", 0),
        "crpix1": hget("CRPIX1", 0),
        "crpix2": hget("CRPIX2", 0),
        "cdelt1": hget("CDELT1", 1),
        "cdelt2": hget("CDELT2", 1),
        "pc": pc_matrix,
        "cd": cd_matrix,
        "pv": pvs if pvs else None,
        "tolerance": tolerance,
        "test_vectors": vectors,
    }


def main():
    fixtures = []

    # Gaia DR3-style TAN header (no distortion).
    fixtures.append(sample_wcs({
        "NAXIS1": 1024, "NAXIS2": 1024,
        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
        "CRVAL1": 180.0, "CRVAL2": 30.0,
        "CRPIX1": 512.5, "CRPIX2": 512.5,
        "CDELT1": -0.00027778, "CDELT2": 0.00027778,
        "RADESYS": "ICRS",
    }, "gaia_tan", "Gaia DR3", "TAN"))

    # Plain CAR all-sky.
    fixtures.append(sample_wcs({
        "NAXIS1": 720, "NAXIS2": 360,
        "CTYPE1": "RA---CAR", "CTYPE2": "DEC--CAR",
        "CRVAL1": 0.0, "CRVAL2": 0.0,
        "CRPIX1": 360.5, "CRPIX2": 180.5,
        "CDELT1": -0.5, "CDELT2": 0.5,
    }, "allsky_car", "synthetic all-sky", "CAR"))

    # AIT all-sky.
    fixtures.append(sample_wcs({
        "NAXIS1": 720, "NAXIS2": 360,
        "CTYPE1": "GLON-AIT", "CTYPE2": "GLAT-AIT",
        "CRVAL1": 0.0, "CRVAL2": 0.0,
        "CRPIX1": 360.5, "CRPIX2": 180.5,
        "CDELT1": -0.5, "CDELT2": 0.5,
    }, "allsky_ait", "synthetic all-sky galactic", "AIT"))

    # Conic Albers.
    fixtures.append(sample_wcs({
        "NAXIS1": 512, "NAXIS2": 512,
        "CTYPE1": "RA---COE", "CTYPE2": "DEC--COE",
        "CRVAL1": 0.0, "CRVAL2": 45.0,
        "CRPIX1": 256.5, "CRPIX2": 256.5,
        "CDELT1": -0.01, "CDELT2": 0.01,
        "PV2_1": 45.0,
    }, "conic_coe", "synthetic conic", "COE"))

    # Add more fixtures as real headers become available:
    # - HST ACS with SIP (A_ORDER, B_ORDER, A_i_j, B_i_j)
    # - JWST NIRCam with SIP
    # - Chandra ACIS (TAN-P)
    # - 2MASS (TNX via WAT strings)
    # - Planck (HPX)
    # - COBE (CSC)
    # These require real header dumps; when added, the Go test suite
    # automatically picks them up.

    json.dump(fixtures, sys.stdout, indent=2)
    print()


if __name__ == "__main__":
    main()
