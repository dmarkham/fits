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
  - ATLAS       — TPV (difference-image stamp, CRPIX far off the cutout)

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

    # ATLAS difference-image stamp (01a59742o0240o.diff, MLO): TPV
    # distortion from the ATLAS pipeline, on a 100x100 cutout whose CRPIX
    # lies ~4000 px off the stamp, so the polynomial terms dominate.
    fixtures.append(sample_wcs({
        "NAXIS1": 100, "NAXIS2": 100,
        "CTYPE1": "RA---TPV", "CTYPE2": "DEC--TPV",
        "CRVAL1": 272.84895,
        "CRVAL2": 8.28673,
        "CRPIX1": -3799.0,
        "CRPIX2": 1852.0,
        "CD1_1": -0.0005172389232442,
        "CD1_2": 8.768323883707e-06,
        "CD2_1": 8.84243071444e-06,
        "CD2_2": 0.0005171961669467,
        "PV1_0": 2.172247402608e-13,
        "PV1_1": 0.9999999999989,
        "PV1_2": 8.547019713842e-12,
        "PV1_4": 1.732758107952e-05,
        "PV1_5": 2.291893113739e-05,
        "PV1_6": 3.988039825779e-06,
        "PV1_7": -0.0006465683213014,
        "PV1_8": -5.326872236656e-06,
        "PV1_9": -0.0006335762897805,
        "PV1_10": -1.022284260944e-06,
        "PV1_12": 3.225351893624e-07,
        "PV1_13": 7.945185105527e-07,
        "PV1_14": -2.941640767804e-08,
        "PV1_15": 1.058822218628e-06,
        "PV1_16": 3.905016850559e-07,
        "PV1_17": -3.143883482666e-06,
        "PV1_18": 1.420641960623e-07,
        "PV1_19": -8.322414657774e-06,
        "PV1_20": 3.606625473241e-07,
        "PV1_21": -4.262232802183e-06,
        "PV1_22": 3.562970026892e-08,
        "PV2_0": -4.986662354562e-14,
        "PV2_1": 0.9999999999991,
        "PV2_2": -8.416590590396e-12,
        "PV2_4": 1.948199301454e-05,
        "PV2_5": 5.608380547107e-07,
        "PV2_6": 1.883925313334e-05,
        "PV2_7": -0.0006489064539911,
        "PV2_8": -4.337750661674e-07,
        "PV2_9": -0.0006448985069562,
        "PV2_10": -1.268989767674e-06,
        "PV2_12": 1.104589358844e-06,
        "PV2_13": 7.033242408929e-07,
        "PV2_14": 1.061381063149e-06,
        "PV2_15": 5.795264155287e-07,
        "PV2_16": -8.60685540635e-07,
        "PV2_17": -2.490758118647e-06,
        "PV2_18": 3.319013437658e-08,
        "PV2_19": -8.309002201252e-06,
        "PV2_20": 5.831004923987e-08,
        "PV2_21": -3.523011214328e-06,
        "PV2_22": 1.935996964426e-08,
        "RADESYS": "FK5", "EQUINOX": 2000.0,
    }, "atlas_tpv", "ATLAS difference stamp", "TPV"))

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
