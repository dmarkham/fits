#!/usr/bin/env python3
"""Write the synthetic high-order distortion cases hdrcases/hiorder_*.case.

The real SIP headers in hdrcases/ stop at order 5 and the real TPV
headers use PVi_0..16 (plus the radial slots in tpv_radial_cd), so these
cases set EVERY term the code accepts and nothing proves otherwise:

  hiorder_sip7, hiorder_sip9  TAN-SIP A/B order 7 and 9 (wcs.MaxSIPOrder),
                              every term 2 <= p+q <= order on both axes,
                              no AP/BP (both sides invert A/B numerically)
  hiorder_tpv_cd, _pc         TPV with every PVi_0..PVi_39 non-zero on both
                              axes (degree 0..7 and r, r^3, r^5, r^7), CD
                              form and PC + CDELT form

The coefficients follow a fixed formula (in each case file's comment) so
each term moves the chip edge by a visible amount: dropping a term or
giving a degree-6/7 slot the wrong powers changes the result by far more
than the tolerance. The truth is wcslib 8.10 (hdr2sky; make hdr-goldens).

Usage (stdlib only):  python3 gen_highorder_cases.py [OUTROOT]
OUTROOT defaults to this directory (cases go to OUTROOT/hdrcases);
TestHighOrderCaseGeneratorReproduce runs it into a temp OUTROOT and
requires the committed case files back byte for byte.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE, "hdrcases")

# TOL_PX: the s2p input is wcslib's p2s sky position, whose angular
# round-off (~2e-10 arcsec, the p2s budget measured on every case) is a
# few 1e-10 px at these plate scales; the Newton/fixed-point inverses
# converge to ~1e-12 px. Known SIP/TPV inverse bugs were >= 1e-4 px.
TOL = "TOL_PX 1e-8"


def card(k, v):
    return f"{k:<8}= {v}"


def sip_case(order, half):
    n = 2 * half
    L = [
        f"# TAN-SIP A/B order {order} with EVERY term 2 <= p+q <= {order} set on both",
        f"# axes (synthetic, gen_highorder_cases.py). Coefficient of u^p v^q:",
        f"# s * c / {half}^(p+q) with c = 0.8 px (A) or 0.6 px (B) and",
        "# s = (1 + 0.11 p + 0.07 q), sign (-1)^(p+q) for A and (-1)^p for B, so",
        f"# every term moves the chip corner ({half} px from CRPIX on both axes) by",
        "# 0.6-1.9 px: a dropped or swapped order-6+ term is visible. No AP/BP:",
        "# sky -> pixel inverts A/B numerically.",
        f"# {TOL[:6]}: see gen_highorder_cases.py.",
        TOL,
        f"NOTE pipeline: SIP order {order}, every A/B term set (high-order SIP)",
        f"GRID 1 {n} 9 1 {n} 9",
        f"PIX {half + 0.5} {half + 0.5}",
        "HEADER",
    ]
    hdr = [("NAXIS", "2"), ("NAXIS1", str(n)), ("NAXIS2", str(n)),
           ("CTYPE1", "'RA---TAN-SIP'"), ("CTYPE2", "'DEC--TAN-SIP'"),
           ("CRVAL1", "210.80227"), ("CRVAL2", "54.34895"),
           ("CRPIX1", f"{half + 0.5}"), ("CRPIX2", f"{half + 0.5}"),
           ("CD1_1", "-2.8E-04"), ("CD1_2", "1.5E-06"), ("CD2_1", "1.2E-06"), ("CD2_2", "2.8E-04")]
    for ax, c0 in (("A", 0.8), ("B", 0.6)):
        hdr.append((f"{ax}_ORDER", str(order)))
        for d in range(2, order + 1):
            for p in range(d, -1, -1):
                q = d - p
                sign = (-1) ** (p + q) if ax == "A" else (-1) ** p
                v = sign * (1 + 0.11 * p + 0.07 * q) * c0 / half ** d
                hdr.append((f"{ax}_{p}_{q}", "%.6E" % v))
    return L + [card(k, v) for k, v in hdr]


# TPV slot m -> (xi power, eta power, r power) for axis 1 (the
# registry/wcslib layout; axis 2 swaps the roles of xi and eta).
TPV_SLOTS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1),
             (2, 0, 0), (1, 1, 0), (0, 2, 0),
             (3, 0, 0), (2, 1, 0), (1, 2, 0), (0, 3, 0), (0, 0, 3),
             (4, 0, 0), (3, 1, 0), (2, 2, 0), (1, 3, 0), (0, 4, 0),
             (5, 0, 0), (4, 1, 0), (3, 2, 0), (2, 3, 0), (1, 4, 0), (0, 5, 0), (0, 0, 5),
             (6, 0, 0), (5, 1, 0), (4, 2, 0), (3, 3, 0), (2, 4, 0), (1, 5, 0), (0, 6, 0),
             (7, 0, 0), (6, 1, 0), (5, 2, 0), (4, 3, 0), (3, 4, 0), (2, 5, 0), (1, 6, 0), (0, 7, 0), (0, 0, 7)]
assert len(TPV_SLOTS) == 40


def tpv_case(pc):
    L = [
        "# TPV with EVERY coefficient PVi_0..PVi_39 non-zero on BOTH axes",
        "# (synthetic, gen_highorder_cases.py): degree 0-7 polynomial terms and",
        "# the radial r, r^3, r^5, r^7 terms. Slot m of degree d gets",
        ("# s_m * 1.37 / 2048^d (px; TPV acts on PC*(p - CRPIX) before CDELT)," if pc else
         "# s_m * 1e-4 / 0.15^d (deg),"),
        "# s_m = +-(1 + 0.13 (m mod 7)), except PVi_1 = 1 + 1e-3 s_1. At the chip",
        "# edge (|eta| ~ 0.15 deg, |xi| ~ 0.075 deg) a term is up to ~0.7 arcsec,",
        "# so a degree-6/7 slot with the wrong xi/eta powers is visible.",
        "# PC + CDELT form (the polynomial applies before CDELT, as in wcslib)." if pc else "# Rotated CD form.",
        f"# {TOL[:6]}: see gen_highorder_cases.py.",
        TOL,
        "NOTE pipeline: TPV every PVi_0..PVi_39 on both axes" + (" (PC + CDELT)" if pc else " (CD form)"),
        "GRID 1 2048 9 1 4096 9",
        "PIX 1024.5 2048.5",
        "HEADER",
    ]
    h = [("NAXIS", "2"), ("NAXIS1", "2048"), ("NAXIS2", "4096"),
         ("CTYPE1", "'RA---TPV'"), ("CTYPE2", "'DEC--TPV'"),
         ("CRVAL1", "35.5"), ("CRVAL2", "-4.75"), ("CRPIX1", "1024.5"), ("CRPIX2", "2048.5")]
    if pc:
        h += [("CDELT1", "7.3E-05"), ("CDELT2", "7.3E-05"), ("PC1_1", "0.99999"),
              ("PC1_2", "1.5E-03"), ("PC2_1", "-1.6E-03"), ("PC2_2", "0.99998")]
    else:
        h += [("CD1_1", "7.3E-05"), ("CD1_2", "1.1E-07"), ("CD2_1", "-1.2E-07"), ("CD2_2", "7.3E-05")]
    for ax in (1, 2):
        for m, (a, b, r) in enumerate(TPV_SLOTS):
            s = (1 + 0.13 * (m % 7)) * (-1 if (m + ax) % 2 else 1)
            # CD form: xi, eta in degrees (chip edge ~0.15 deg). PC + CDELT
            # form: wcslib applies TPV before CDELT, on PC*(p - CRPIX) in
            # pixels (edge ~2048 px), so the same corrections in pixels
            # are 1e-4 deg / 7.3e-5 deg/px = 1.37 px over 2048^d.
            amp, edge = (1.37, 2048.0) if pc else (1e-4, 0.15)
            v = 1 + 1e-3 * s if m == 1 else s * amp / edge ** (a + b + r)
            h.append((f"PV{ax}_{m}", "%.6E" % v))
    return L + [card(k, v) for k, v in h]


def main():
    os.makedirs(OUT, exist_ok=True)
    for name, lines in (("hiorder_sip7", sip_case(7, 1024)),
                        ("hiorder_sip9", sip_case(9, 1024)),
                        ("hiorder_tpv_cd", tpv_case(False)),
                        ("hiorder_tpv_pc", tpv_case(True))):
        with open(os.path.join(OUT, name + ".case"), "w") as f:
            f.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
