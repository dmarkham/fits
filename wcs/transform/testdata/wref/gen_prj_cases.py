#!/usr/bin/env python3
"""Write the per-projection full-pipeline case files hdrcases/prj_*.case.

Every projection code gets (at least) two synthetic headers run through
hdr2sky (wcslib wcspih -> wcsset -> wcsp2s/wcss2p):

  *_c0   CRVAL2 = 0 (default PV unless the projection requires one)
  *_cN   CRVAL2 != theta_0 (and != 0), usually with non-default PVi_m

Both use 1 deg/pixel and a pixel grid whose intermediate world
coordinates span well past the projection's boundary, so wcslib's
accept/reject decision for off-projection pixels is part of the golden.
The grid nodes and CRPIX are deliberately non-commensurate (CRPIX =
250.37/250.61, grid step ~35.5 px) so no grid point lands exactly on a
seam (x = +-180 for the cylindricals, HPX facet edges, cube face edges)
where the sky->pixel answer is two-valued. The sky grid likewise avoids
the celestial poles and CRVAL1 +- 180.

No LONPOLE/LATPOLE card is written: these cases prove the projections and
the default native-pole choice; explicit LONPOLE belongs to the pipeline
cases.

Usage (stdlib only):  python3 gen_prj_cases.py [OUTROOT]   (then: make hdr-goldens)
OUTROOT defaults to this directory (cases go to OUTROOT/hdrcases);
TestPrjCaseGeneratorReproduce runs it into a temp OUTROOT and requires
the committed case files back byte for byte.
The case files are committed; this script documents and reproduces them.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE, "hdrcases")

# (case suffix, code, CRVAL1, CRVAL2, {m: PV2_m}, note)
CASES = [
    # zenithal: theta_0 = 90
    ("azp_c0", "AZP", 37.3, 0.0, {}, "AZP defaults (mu = 0, gamma = 0)"),
    ("azp_c45_mu2_g30", "AZP", 37.3, 45.0, {1: 2.0, 2: 30.0}, "AZP mu = 2, tilted gamma = 30"),
    ("szp_c0", "SZP", 37.3, 0.0, {}, "SZP defaults"),
    ("szp_cm40_mu1p5", "SZP", 37.3, -40.0, {1: 1.5, 2: 90.0, 3: 45.0}, "SZP mu = 1.5, phi_c = 90, theta_c = 45"),
    ("tan_c0", "TAN", 37.3, 0.0, {}, "TAN"),
    ("tan_c60", "TAN", 37.3, 60.0, {}, "TAN"),
    ("stg_c0", "STG", 37.3, 0.0, {}, "STG"),
    ("stg_cm35", "STG", 37.3, -35.0, {}, "STG"),
    ("sin_c0", "SIN", 37.3, 0.0, {}, "SIN"),
    ("sin_c50_slant", "SIN", 37.3, 50.0, {1: 0.2, 2: -0.3}, "slant SIN xi = 0.2, eta = -0.3"),
    ("arc_c0", "ARC", 37.3, 0.0, {}, "ARC"),
    ("arc_cm35", "ARC", 37.3, -35.0, {}, "ARC"),
    ("zpn_c0_cubic", "ZPN", 37.3, 0.0, {1: 1.0, 3: -0.2}, "ZPN degree 3 with a stationary point"),
    ("zpn_c70_quintic", "ZPN", 37.3, 70.0, {1: 0.9, 3: 0.1, 5: -0.02}, "ZPN degree 5"),
    ("zpn_c0_deg12", "ZPN", 37.3, 0.0, {1: 1.0, 8: 2e-4, 10: -1e-4, 12: 1e-5}, "ZPN degree 12 (PV2_8, PV2_10, PV2_12)"),
    ("zea_c0", "ZEA", 37.3, 0.0, {}, "ZEA"),
    ("zea_cm35", "ZEA", 37.3, -35.0, {}, "ZEA"),
    ("air_c0", "AIR", 37.3, 0.0, {}, "AIR theta_b = 90 (default)"),
    ("air_c30_tb45", "AIR", 37.3, 30.0, {1: 45.0}, "AIR theta_b = 45"),
    # cylindrical: theta_0 = 0
    ("cyp_c0", "CYP", 37.3, 0.0, {}, "CYP defaults (mu = 1, lambda = 1)"),
    ("cyp_c30_gall", "CYP", 37.3, 30.0, {1: 1.0, 2: 0.70710678118654752}, "CYP Gall: mu = 1, lambda = sqrt(2)/2"),
    ("cea_c0", "CEA", 37.3, 0.0, {}, "CEA lambda = 1"),
    ("cea_c30_l0p5", "CEA", 37.3, 30.0, {1: 0.5}, "CEA lambda = 0.5"),
    ("car_c0", "CAR", 37.3, 0.0, {}, "CAR"),
    ("car_c30", "CAR", 37.3, 30.0, {}, "CAR"),
    ("mer_c0", "MER", 37.3, 0.0, {}, "MER"),
    ("mer_cm40", "MER", 37.3, -40.0, {}, "MER"),
    # pseudo-cylindrical: theta_0 = 0
    ("sfl_c0", "SFL", 37.3, 0.0, {}, "SFL"),
    ("sfl_c30", "SFL", 37.3, 30.0, {}, "SFL"),
    ("par_c0", "PAR", 37.3, 0.0, {}, "PAR"),
    ("par_cm40", "PAR", 37.3, -40.0, {}, "PAR"),
    ("mol_c0", "MOL", 37.3, 0.0, {}, "MOL"),
    ("mol_c30", "MOL", 37.3, 30.0, {}, "MOL"),
    ("ait_c0", "AIT", 37.3, 0.0, {}, "AIT"),
    ("ait_cm40", "AIT", 37.3, -40.0, {}, "AIT"),
    # conic: theta_0 = theta_a
    ("cop_c0_ta45", "COP", 37.3, 0.0, {1: 45.0}, "COP theta_a = 45"),
    ("cop_cm20_tam30_eta15", "COP", 37.3, -20.0, {1: -30.0, 2: 15.0}, "COP southern, theta_a = -30, eta = 15"),
    ("coe_c0_ta45", "COE", 37.3, 0.0, {1: 45.0}, "COE theta_a = 45"),
    ("coe_cm20_tam30_eta15", "COE", 37.3, -20.0, {1: -30.0, 2: 15.0}, "COE southern, theta_a = -30, eta = 15"),
    ("cod_c0_ta45", "COD", 37.3, 0.0, {1: 45.0}, "COD theta_a = 45"),
    ("cod_c60_ta45_eta10", "COD", 37.3, 60.0, {1: 45.0, 2: 10.0}, "COD theta_a = 45, eta = 10"),
    ("coo_c0_ta45", "COO", 37.3, 0.0, {1: 45.0}, "COO theta_a = 45"),
    ("coo_cm20_tam30_eta15", "COO", 37.3, -20.0, {1: -30.0, 2: 15.0}, "COO southern, theta_a = -30, eta = 15"),
    # polyconic: theta_0 = 0
    ("bon_c0_t30", "BON", 37.3, 0.0, {1: 30.0}, "BON theta_1 = 30"),
    ("bon_c50_tm45", "BON", 37.3, 50.0, {1: -45.0}, "BON theta_1 = -45"),
    ("pco_c0", "PCO", 37.3, 0.0, {}, "PCO"),
    ("pco_c30", "PCO", 37.3, 30.0, {}, "PCO"),
    # quad-cube: theta_0 = 0
    ("tsc_c0", "TSC", 37.3, 0.0, {}, "TSC"),
    ("tsc_c30", "TSC", 37.3, 30.0, {}, "TSC"),
    ("csc_c0", "CSC", 37.3, 0.0, {}, "CSC"),
    ("csc_c30", "CSC", 37.3, 30.0, {}, "CSC"),
    ("qsc_c0", "QSC", 37.3, 0.0, {}, "QSC"),
    ("qsc_c30", "QSC", 37.3, 30.0, {}, "QSC"),
    # HEALPix
    ("hpx_c0", "HPX", 37.3, 0.0, {}, "HPX H = 4, K = 3 (default)"),
    ("hpx_c30", "HPX", 37.3, 30.0, {}, "HPX H = 4, K = 3 (default)"),
    ("hpx_c0_h3k2", "HPX", 37.3, 0.0, {1: 3.0, 2: 2.0}, "HPX H = 3, K = 2 (even K: offset southern facets)"),
    ("hpx_c0_h6k1", "HPX", 37.3, 0.0, {1: 6.0, 2: 1.0}, "HPX H = 6, K = 1 (no equatorial zone)"),
    ("xph_c0", "XPH", 37.3, 0.0, {}, "XPH"),
    ("xph_cm50", "XPH", 37.3, -50.0, {}, "XPH"),
]

# Pixel grid: 1 deg/pixel; intermediate coordinates span roughly
# [-249, +249] deg in both axes (quad-cube: [-81, +417] in x, where the
# unfolded cube reaches 315 deg).
GRID = "GRID 1.13 498.71 15 1.29 498.53 15"


def case_text(name, code, crval1, crval2, pv, note):
    crpix1 = 250.37
    if code in ("TSC", "CSC", "QSC"):
        crpix1 = 418.37  # x = -(p - crpix1) spans about [-81, +417]
    lines = [
        "# Generated by gen_prj_cases.py: synthetic %s header for the" % code,
        "# full-pipeline projection proof (see that script for the design).",
        "NOTE projection %s: %s; CRVAL2 = %g" % (code, note, crval2),
        GRID,
        "SKYGRID 0.71 359.13 13 -89.37 89.21 13",
    ]
    if code == "COP":
        # COP diverges as theta -> theta_a - 90: sky points near there map
        # to |y| ~ 4e4 px, where 1e-9 px is below 3e-14 relative (~ULP
        # level). Measured worst there: 6.3e-9 px (1.6e-13 relative).
        lines += [
            "# Sky points near the COP divergence map to |pixel| ~ 4e4, so the",
            "# default 1e-9 px absolute budget is below the float64 relative",
            "# precision there; 5e-8 px is 1.2e-12 relative at 4e4 px.",
            "TOL_PX 5e-8",
        ]
    lines += [
        "HEADER",
        "NAXIS   = 2",
        "NAXIS1  = 500",
        "NAXIS2  = 500",
        "CTYPE1  = 'RA---%s'" % code,
        "CTYPE2  = 'DEC--%s'" % code,
        "CRVAL1  = %r" % crval1,
        "CRVAL2  = %r" % crval2,
        "CRPIX1  = %r" % crpix1,
        "CRPIX2  = 250.61",
        "CDELT1  = -1.0",
        "CDELT2  = 1.0",
    ]
    for m in sorted(pv):
        lines.append("PV2_%-4d= %r" % (m, pv[m]))
    return "\n".join(lines) + "\n"


# Native fiducial point given on the LONGITUDE axis (PV1_1 = phi_0,
# PV1_2 = theta_0; PV1_0 = 1 also offsets (x, y) so the fiducial point is
# at (0, 0)). wcslib applies them; these cases pin what a header that sets
# them must produce. (case suffix, code, CRVAL2, {m: PV2_m}, {m: PV1_m}, note)
FIDUCIAL = [
    ("car_pv1_t20", "CAR", 20.0, {}, {1: 0.0, 2: 20.0},
     "CAR with theta_0 = 20 via PV1_2 (no PV1_0 offset)"),
    ("car_pv1_off", "CAR", 20.0, {}, {0: 1.0, 1: 10.0, 2: 20.0},
     "CAR with PV1_0 = 1, phi_0 = 10, theta_0 = 20"),
    ("tan_pv1_off", "TAN", 50.0, {}, {0: 1.0, 1: 0.0, 2: 70.0},
     "TAN with PV1_0 = 1, theta_0 = 70"),
]


# Latitude on the FIRST axis (CTYPE1 = 'DEC--xxx', CTYPE2 = 'RA---xxx'):
# the projection parameters are then PV1_m, on the latitude axis, and a
# PV2_m card belongs to the longitude axis. Proves that the PVi_m are read
# from the latitude axis, whichever it is. (case suffix, code, CRVAL lat,
# {m: PV on the latitude axis}, note)
SWAPPED = [
    ("azp_mu2_g30", "AZP", 45.0, {1: 2.0, 2: 30.0}, "AZP mu = 2, gamma = 30, latitude first"),
    ("zpn_cubic", "ZPN", 0.0, {1: 1.0, 3: -0.2}, "ZPN degree 3, latitude first"),
    ("cea_l0p5", "CEA", 30.0, {1: 0.5}, "CEA lambda = 0.5, latitude first"),
    ("cop_tam30_eta15", "COP", -20.0, {1: -30.0, 2: 15.0}, "COP southern, latitude first"),
]


def swapped_text(name, code, crlat, pv, note):
    text = case_text(name, code, 37.3, crlat, {}, note)
    text = text.replace("# full-pipeline projection proof", "# latitude-first axis order proof")
    text = text.replace("; CRVAL2 = %g" % crlat, "; latitude CRVAL1 = %g" % crlat)
    text = text.replace("CTYPE1  = 'RA---%s'" % code, "CTYPE1  = 'DEC--%s'" % code)
    text = text.replace("CTYPE2  = 'DEC--%s'" % code, "CTYPE2  = 'RA---%s'" % code)
    text = text.replace("CRVAL1  = 37.3\nCRVAL2  = %r" % crlat, "CRVAL1  = %r\nCRVAL2  = 37.3" % crlat)
    text = text.replace("CDELT1  = -1.0\nCDELT2  = 1.0", "CDELT1  = 1.0\nCDELT2  = -1.0")
    text += "".join("PV1_%-4d= %r\n" % (m, pv[m]) for m in sorted(pv))
    return text


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, code, crlat, pv, note in SWAPPED:
        with open(os.path.join(OUT, "prjswap_%s.case" % name), "w") as f:
            f.write(swapped_text(name, code, crlat, pv, note))
    for name, code, c1, c2, pv, note in CASES:
        path = os.path.join(OUT, "prj_%s.case" % name)
        with open(path, "w") as f:
            f.write(case_text(name, code, c1, c2, pv, note))
    for name, code, c2, pv, pv1, note in FIDUCIAL:
        text = case_text(name, code, 37.3, c2, pv, note)
        text = text.replace("# full-pipeline projection proof", "# native-fiducial-point (PV1_m) proof")
        text += "".join("PV1_%-4d= %r\n" % (m, pv1[m]) for m in sorted(pv1))
        with open(os.path.join(OUT, "pvfid_%s.case" % name), "w") as f:
            f.write(text)
    codes = sorted({c[1] for c in CASES})
    print("wrote %d cases covering %d codes" % (len(CASES), len(codes)))


if __name__ == "__main__":
    main()
