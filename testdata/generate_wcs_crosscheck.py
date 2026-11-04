#!/usr/bin/env python3
"""Generate testdata/wcs_crosscheck.json: astropy.wcs truth through the
FULL header path.

Usage (from the repository root, with the project venv):

    venv/bin/python testdata/generate_wcs_crosscheck.py > testdata/wcs_crosscheck.json

Re-running it with the pinned versions (astropy 8.0.2 / numpy 2.4.6, see
the "astropy"/"numpy" fields of the JSON) reproduces the committed file
byte for byte; TestCrossCheckFixturesReproduce checks that.

For every fixture the script builds an astropy.io.fits.Header from the
keyword dict below, and records:

  - "cards": the header exactly as 80-column card images
    (Header.tostring), which the Go test feeds through
    header.ParseCards -> wcs.Parse -> transform.New, so the Go parser is
    part of what is checked (no hand-built wcs.Header);
  - "p2s": pixel -> sky from astropy.wcs.WCS.all_pix2world (wcslib core
    plus astropy's SIP), 0-based pixels on a grid over the image; points
    astropy maps to NaN (outside the projection) are dropped;
  - "s2p": sky -> pixel from WCS.all_world2pix at each p2s sky position
    (tolerance 1e-12 px, at most 100 iterations; astropy raises if it
    does not converge).

These fixtures are synthetic or transcribed headers (see "source"). The
REAL archive headers (HST, JWST, Spitzer, ZTF, PTF, astrometry.net) are
checked against wcslib directly by the header goldens in
wcs/transform/testdata/wref/hdrcases.
"""

import json
import sys

try:
    import numpy as np
    import astropy
    import astropy.wcs
    from astropy.wcs import WCS
    from astropy.io.fits import Header
except ImportError:
    print("error: astropy not installed; use the project venv", file=sys.stderr)
    sys.exit(1)


def sample_wcs(header_cards, name, source, tol_arcsec=1e-8, tol_px=1e-9, n=5):
    """Build the WCS from header_cards (an ordered dict of keywords), sample
    pixel -> sky on an n x n grid and sky -> pixel at the results."""
    hdr = Header()
    for k, v in header_cards.items():
        hdr[k] = v
    w = WCS(hdr)
    naxis1 = header_cards.get("NAXIS1", 1024)
    naxis2 = header_cards.get("NAXIS2", 1024)
    xs = np.linspace(0, naxis1 - 1, n)
    ys = np.linspace(0, naxis2 - 1, n)
    p2s, s2p = [], []
    for y in ys:
        for x in xs:
            sky = w.all_pix2world([[x, y]], 0)[0]
            if not (np.isfinite(sky[0]) and np.isfinite(sky[1])):
                continue
            p2s.append({"pixel": [float(x), float(y)], "sky": [float(sky[0]), float(sky[1])]})
    for v in p2s:
        pix = w.all_world2pix([v["sky"]], 0, tolerance=1e-12, maxiter=100, quiet=False)[0]
        s2p.append({"sky": v["sky"], "pixel": [float(pix[0]), float(pix[1])]})
    text = hdr.tostring(sep="", endcard=False, padding=False)
    cards = [text[i:i + 80] for i in range(0, len(text), 80)]
    return {
        "name": name,
        "source": source,
        "cards": cards,
        "tol_arcsec": tol_arcsec,
        "tol_px": tol_px,
        "p2s": p2s,
        "s2p": s2p,
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
    }, "gaia_tan", "synthetic TAN (Gaia-like 1 arcsec/px tile)"))

    # Plain CAR all-sky.
    fixtures.append(sample_wcs({
        "NAXIS1": 720, "NAXIS2": 360,
        "CTYPE1": "RA---CAR", "CTYPE2": "DEC--CAR",
        "CRVAL1": 0.0, "CRVAL2": 0.0,
        "CRPIX1": 360.5, "CRPIX2": 180.5,
        "CDELT1": -0.5, "CDELT2": 0.5,
    }, "allsky_car", "synthetic all-sky CAR"))

    # AIT all-sky.
    fixtures.append(sample_wcs({
        "NAXIS1": 720, "NAXIS2": 360,
        "CTYPE1": "GLON-AIT", "CTYPE2": "GLAT-AIT",
        "CRVAL1": 0.0, "CRVAL2": 0.0,
        "CRPIX1": 360.5, "CRPIX2": 180.5,
        "CDELT1": -0.5, "CDELT2": 0.5,
    }, "allsky_ait", "synthetic all-sky galactic AIT"))

    # Conic Albers.
    fixtures.append(sample_wcs({
        "NAXIS1": 512, "NAXIS2": 512,
        "CTYPE1": "RA---COE", "CTYPE2": "DEC--COE",
        "CRVAL1": 0.0, "CRVAL2": 45.0,
        "CRPIX1": 256.5, "CRPIX2": 256.5,
        "CDELT1": -0.01, "CDELT2": 0.01,
        "PV2_1": 45.0,
    }, "conic_coe", "synthetic conic COE"))

    # ATLAS difference-image stamp (01a59742o0240o.diff, MLO): TPV
    # distortion from the ATLAS pipeline, on a 100x100 cutout whose CRPIX
    # lies ~4000 px off the stamp, so the polynomial terms dominate. The
    # stamp's real WCS cards (these TPV values plus the SIP copy of the
    # same correction the pipeline writes next to them) are committed
    # since PR #4: wcs/transform/testdata/wref/hdrcases_differ/
    # real_atlas_k2_tpv_sip.case.
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
    }, "atlas_tpv", "ATLAS difference stamp 01a59742o0240o.diff (values transcribed in PR #3; file not in repo)"))

    # UH88 TAN-SIP header (astrometry.net solution), full 5280x5280 frame:
    # SIP corrections reach ~2 px at the corners.
    fixtures.append(sample_wcs({
        "NAXIS": 2, "NAXIS1": 5280, "NAXIS2": 5280,
        "CTYPE1": "RA---TAN-SIP", "CTYPE2": "DEC--TAN-SIP",
        "CUNIT1": "deg", "CUNIT2": "deg",
        "CRVAL1": 101.531931895, "CRVAL2": 28.9319723058,
        "CRPIX1": 2640.5, "CRPIX2": 2640.5,
        "CD1_1": 4.56664458396e-05, "CD1_2": 1.65789570401e-07,
        "CD2_1": 1.70901177979e-07, "CD2_2": -4.56699743556e-05,
        "A_ORDER": 2,
        "A_0_2": -4.51136534482e-08, "A_1_1": -1.02357905543e-07, "A_2_0": -5.43669607604e-08,
        "B_ORDER": 2,
        "B_0_2": -7.143588012e-08, "B_1_1": -2.07643869017e-07, "B_2_0": -1.51807826685e-07,
        "AP_ORDER": 2,
        "AP_0_0": 5.53851868084e-07, "AP_0_1": 1.77638829679e-07, "AP_0_2": 4.51135313836e-08,
        "AP_1_0": 1.98630755541e-07, "AP_1_1": 1.0235777946e-07, "AP_2_0": 5.43668700797e-08,
        "BP_ORDER": 2,
        "BP_0_0": 1.1123140637e-06, "BP_0_1": 3.53518741358e-07, "BP_0_2": 7.143564447e-08,
        "BP_1_0": 4.11686097018e-07, "BP_1_1": 2.07643618938e-07, "BP_2_0": 1.51807636879e-07,
        "RADESYS": "FK5", "EQUINOX": 2000.0,
    }, "uh88_sip", "UH88 astrometry.net solution (values transcribed in PR #3; file not in repo)"))

    out = {
        "generator": "testdata/generate_wcs_crosscheck.py",
        "astropy": astropy.__version__,
        "wcslib": astropy.wcs._wcs.WCSLIB_VERSION,
        "numpy": np.__version__,
        "fixtures": fixtures,
    }
    json.dump(out, sys.stdout, indent=1)
    print()


if __name__ == "__main__":
    main()
