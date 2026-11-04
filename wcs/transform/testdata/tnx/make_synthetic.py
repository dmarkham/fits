#!/usr/bin/env python3
"""Write the synthetic TNX headers in headers/syn_*.hdr.

The real headers (noao_ct4m_tnx1, registry_tnxdoc_tnx3) cover Chebyshev
half-cross and polynomial half-cross surfaces with CD matrices. These
synthetic headers fill in the rest of the TNX surface space: Legendre,
no and full cross terms, unequal x/y orders, PC + CDELT instead of CD, a
reference point at the celestial pole (TNX default LONPOLE = 180), a
"longpole" WAT attribute, a LONPOLE card alone and beside a WAT longpole,
the WAT1/WAT2 precedence of longpole, lngcor and latcor, the
cross-axis lngcor/latcor fallback, galactic axes, a TNX header with no
WAT cards, and a WAT string at the WCSTools length limit. Truth for
all of them comes from WCSTools (see
gen_golden.py); nothing here is evaluated.

The WAT strings are split exactly as IRAF does: the attribute string is
cut every 68 characters with no regard for token boundaries, so chunk
boundaries fall inside numbers and after trailing blanks. Run from this
directory:

    python3 make_synthetic.py           # rewrite headers/syn_*.hdr
    python3 make_synthetic.py --check   # exit 1 if any would change, is
                                        # missing, or is not one of CASES
"""


def card(key, value=None, comment=None):
    if value is None:
        return key
    if isinstance(value, bool):
        v = ("T" if value else "F").rjust(20)
    elif isinstance(value, str):
        v = "'" + value.replace("'", "''") + "'"
        v = v.ljust(20)
    else:
        v = repr(value).upper().rjust(20)
    s = "%-8s= %s" % (key, v)
    if comment:
        s += " / " + comment
    assert len(s) <= 80, s
    return s


def wat(axis, text):
    chunks = [text[i:i + 68] for i in range(0, len(text), 68)]
    return [card("WAT%d_%03d" % (axis, n + 1), c) for n, c in enumerate(chunks)]


def surface(stype, xo, yo, xt, dom, coeffs):
    nums = [float(stype), float(xo), float(yo), float(xt)] + list(dom) + list(coeffs)
    out = []
    for v in nums:
        s = repr(v)
        if s.endswith(".0"):
            s = s[:-1]  # IRAF writes "3." not "3.0"
        out.append(s)
    return " ".join(out) + " "


def header(naxis, lin, crval, crpix, wat1, wat2, extra=(), ctype=("RA---TNX", "DEC--TNX")):
    cards = [card("SIMPLE", True), card("BITPIX", 16), card("NAXIS", 2),
             card("NAXIS1", naxis[0]), card("NAXIS2", naxis[1]),
             card("WCSDIM", 2), card("CTYPE1", ctype[0]), card("CTYPE2", ctype[1]),
             card("CRVAL1", crval[0]), card("CRVAL2", crval[1]),
             card("CRPIX1", crpix[0]), card("CRPIX2", crpix[1])]
    cards += [card(k, v) for k, v in lin]
    cards += list(extra)
    cards += [card("WAT0_001", "system=image")] + wat(1, wat1) + wat(2, wat2)
    cards.append("END")
    return "\n".join(c.rstrip() for c in cards) + "\n"


CASES = {}

# Legendre, full cross terms, 3x3, rotated CD.
CASES["syn_leg_full"] = header(
    (2048, 2048),
    [("CD1_1", -7.2e-05), ("CD1_2", 4.1e-06), ("CD2_1", 3.9e-06), ("CD2_2", 7.25e-05)],
    (150.11916667, 2.20583333), (1024.5, 1030.25),
    'wtype=tnx axtype=ra lngcor = "%s"' % surface(
        2, 3, 3, 1, (-0.085, 0.087, -0.084, 0.086),
        (1.2e-4, -3.4e-4, 2.2e-4, 5.6e-4, -7.1e-5, 3.3e-5, -2.9e-4, 1.1e-4, -6.6e-5)),
    'wtype=tnx axtype=dec latcor = "%s"' % surface(
        2, 3, 3, 1, (-0.085, 0.087, -0.084, 0.086),
        (-2.1e-4, 1.4e-4, -3.2e-5, -4.4e-4, 8.8e-5, -1.2e-5, 3.1e-4, -9.9e-5, 4.5e-5)))

# Chebyshev, no cross terms, xorder 4 / yorder 3, PC + CDELT.
CASES["syn_cheb_none_pc"] = header(
    (1500, 1200),
    [("CDELT1", -0.000135), ("CDELT2", 0.000135),
     ("PC1_1", 0.9961946981), ("PC1_2", -0.0871557427),
     ("PC2_1", 0.0871557427), ("PC2_2", 0.9961946981)],
    (30.25, -45.5), (750.0, 600.0),
    'wtype=tnx axtype=ra lngcor = "%s"' % surface(
        1, 4, 3, 0, (-0.11, 0.11, -0.09, 0.09),
        (2.5e-4, -1.5e-4, 7.5e-5, -3.0e-5, 1.8e-4, -6.0e-5)),
    'wtype=tnx axtype=dec latcor = "%s"' % surface(
        1, 4, 3, 0, (-0.11, 0.11, -0.09, 0.09),
        (-1.0e-4, 2.2e-4, -4.0e-5, 1.5e-5, -2.6e-4, 9.0e-5)))

# Polynomial, half cross terms, xorder 3 / yorder 4, reference point at the
# celestial pole with no LONPOLE (TNX default 180 deg).
CASES["syn_poly_half_pole"] = header(
    (1024, 1024),
    [("CD1_1", -0.0001), ("CD1_2", 0.0), ("CD2_1", 0.0), ("CD2_2", 0.0001)],
    (0.0, 90.0), (512.5, 512.5),
    'wtype=tnx axtype=ra lngcor = "%s"' % surface(
        3, 3, 4, 2, (-0.06, 0.06, -0.06, 0.06),
        (1.0e-5, 2.0e-3, -3.0e-2, -1.5e-3, 4.0e-2, -2.0e-2, 2.5e-2, 1.0e-2, -5.0e-1)),
    'wtype=tnx axtype=dec latcor = "%s"' % surface(
        3, 3, 4, 2, (-0.06, 0.06, -0.06, 0.06),
        (-2.0e-5, -1.0e-3, 2.0e-2, 3.0e-3, -1.0e-2, 3.0e-2, -4.0e-2, 2.0e-2, 3.0e-1)))

# Polynomial, full cross terms, 2x2, with a "longpole" WAT attribute and no
# LONPOLE card (WCSTools tnxinit reads longpole from the WAT strings).
CASES["syn_poly_full_watlongpole"] = header(
    (800, 800),
    [("CD1_1", -0.0002), ("CD1_2", 0.0), ("CD2_1", 0.0), ("CD2_2", 0.0002)],
    (210.5, 54.25), (400.5, 400.5),
    'wtype=tnx axtype=ra longpole=170. lngcor = "%s"' % surface(
        3, 2, 2, 1, (-0.09, 0.09, -0.09, 0.09), (3.0e-5, 1.2e-3, -8.0e-4, 2.0e-2)),
    'wtype=tnx axtype=dec latcor = "%s"' % surface(
        3, 2, 2, 1, (-0.09, 0.09, -0.09, 0.09), (-1.0e-5, 6.0e-4, 1.1e-3, -1.5e-2)))

# LONPOLE card (175) and no WAT longpole: the card must be used, not the
# TNX default 180. Same surfaces as syn_leg_full.
CASES["syn_leg_lonpole_card"] = header(
    (2048, 2048),
    [("CD1_1", -7.2e-05), ("CD1_2", 4.1e-06), ("CD2_1", 3.9e-06), ("CD2_2", 7.25e-05)],
    (150.11916667, 2.20583333), (1024.5, 1030.25),
    'wtype=tnx axtype=ra lngcor = "%s"' % surface(
        2, 3, 3, 1, (-0.085, 0.087, -0.084, 0.086),
        (1.2e-4, -3.4e-4, 2.2e-4, 5.6e-4, -7.1e-5, 3.3e-5, -2.9e-4, 1.1e-4, -6.6e-5)),
    'wtype=tnx axtype=dec latcor = "%s"' % surface(
        2, 3, 3, 1, (-0.085, 0.087, -0.084, 0.086),
        (-2.1e-4, 1.4e-4, -3.2e-5, -4.4e-4, 8.8e-5, -1.2e-5, 3.1e-4, -9.9e-5, 4.5e-5)),
    extra=[card("LONPOLE", 175.0)])

# LONPOLE card (175) AND a WAT longpole (170): the card wins (WCSTools
# tnxinit reads the WAT longpole only when LONPOLE is absent).
CASES["syn_poly_lonpole_card_and_wat"] = header(
    (800, 800),
    [("CD1_1", -0.0002), ("CD1_2", 0.0), ("CD2_1", 0.0), ("CD2_2", 0.0002)],
    (210.5, 54.25), (400.5, 400.5),
    'wtype=tnx axtype=ra longpole=170. lngcor = "%s"' % surface(
        3, 2, 2, 1, (-0.09, 0.09, -0.09, 0.09), (3.0e-5, 1.2e-3, -8.0e-4, 2.0e-2)),
    'wtype=tnx axtype=dec latcor = "%s"' % surface(
        3, 2, 2, 1, (-0.09, 0.09, -0.09, 0.09), (-1.0e-5, 6.0e-4, 1.1e-3, -1.5e-2)),
    extra=[card("LONPOLE", 175.0)])

# WAT precedence (tnxinit): longpole in both WAT1 (165) and WAT2 (195),
# WAT1 wins; lngcor AND latcor in both strings with different
# coefficients: lngcor is taken from WAT1, latcor from WAT2.
_dom = (-0.09, 0.09, -0.09, 0.09)
CASES["syn_poly_wat_precedence"] = header(
    (800, 800),
    [("CD1_1", -0.0002), ("CD1_2", 0.0), ("CD2_1", 0.0), ("CD2_2", 0.0002)],
    (210.5, 54.25), (400.5, 400.5),
    'wtype=tnx axtype=ra longpole=165. lngcor = "%s" latcor = "%s"' % (
        surface(3, 2, 2, 1, _dom, (3.0e-5, 1.2e-3, -8.0e-4, 2.0e-2)),
        surface(3, 2, 2, 1, _dom, (4.0e-5, -9.0e-4, 1.3e-3, 1.7e-2))),
    'wtype=tnx axtype=dec longpole=195. lngcor = "%s" latcor = "%s"' % (
        surface(3, 2, 2, 1, _dom, (-5.0e-5, 7.0e-4, -1.6e-3, -1.8e-2)),
        surface(3, 2, 2, 1, _dom, (-1.0e-5, 6.0e-4, 1.1e-3, -1.5e-2))))

# lngcor/latcor cross-axis fallback (tnxinit igets str1/str2): lngcor
# only in WAT2 and latcor only in WAT1, so each is read from the other
# axis's string.
_cd800 = [("CD1_1", -0.0002), ("CD1_2", 0.0), ("CD2_1", 0.0), ("CD2_2", 0.0002)]
_L = surface(3, 2, 2, 1, _dom, (3.0e-5, 1.2e-3, -8.0e-4, 2.0e-2))
_B = surface(3, 2, 2, 1, _dom, (-1.0e-5, 6.0e-4, 1.1e-3, -1.5e-2))
CASES["syn_poly_cor_fallback"] = header(
    (800, 800), _cd800, (210.5, 54.25), (400.5, 400.5),
    'wtype=tnx axtype=ra latcor = "%s"' % _B,
    'wtype=tnx axtype=dec lngcor = "%s"' % _L)

# Galactic axes, GLON-TNX / GLAT-TNX (xy2sky reports "galactic").
CASES["syn_poly_galactic"] = header(
    (800, 800), _cd800, (121.75, -33.5), (400.5, 400.5),
    'wtype=tnx axtype=glon lngcor = "%s"' % _L,
    'wtype=tnx axtype=glat latcor = "%s"' % _B,
    ctype=("GLON-TNX", "GLAT-TNX"))

# The other celestial pairs wcs.Parse recognises: ecliptic (xy2sky
# reports "ecliptic"), supergalactic and helioecliptic (no WCSTools
# system name; the TNX math is the same, done in the native frame).
# For a system it cannot name, xy2sky prints a fixed-width field that
# overflows ("****") at 12 decimals once a value reaches 10 deg, so those
# two reference points are kept below 10 deg.
for _name, _ct, _crval in (("syn_poly_ecliptic", ("ELON-TNX", "ELAT-TNX"), (45.125, 12.5)),
                           ("syn_poly_supergalactic", ("SLON-TNX", "SLAT-TNX"), (7.25, -3.5)),
                           ("syn_poly_helioecliptic", ("HLON-TNX", "HLAT-TNX"), (5.5, 2.75))):
    CASES[_name] = header(
        (800, 800), _cd800, _crval, (400.5, 400.5),
        'wtype=tnx lngcor = "%s"' % _L, 'wtype=tnx latcor = "%s"' % _B, ctype=_ct)

# RA---TNX / DEC--TNX with no WAT1/WAT2 cards: WCSTools falls back to TAN
# (tnxinit finds no correction), with the TNX default LONPOLE 180.
CASES["syn_tan_no_wat"] = header(
    (800, 800), _cd800, (210.5, 54.25), (400.5, 400.5), "", "")


# The WCSTools length limit: WCSTools joins WATn into a 2000-byte buffer
# (SZ_ATSTRING) and keeps a string of at most 1999 characters intact. The
# WAT1 string here is exactly 1999 characters and its LAST token is
# longpole=175.5, so a reader that cut even one character would see
# longpole=175. (0.5 deg of rotation). Legendre 9x9 full cross terms.
def _wat_1999():
    co = [((k * 37) % 19 - 9) * 1.3e-6 + k / 7.0e9 for k in range(81)]
    head = 'wtype=tnx axtype=ra lngcor = "%s"' % surface(2, 9, 9, 1, _dom, co)
    tail = " longpole=175.5"
    pad = 1999 - len(head) - len(tail)
    assert pad > 0, len(head)
    s = head + " " * pad + tail
    assert len(s) == 1999
    return s


CASES["syn_leg_wat_1999"] = header(
    (800, 800), _cd800, (210.5, 54.25), (400.5, 400.5),
    _wat_1999(),
    'wtype=tnx axtype=dec latcor = "%s"' % _B)


# --- Accepted header variants that wcslib and WCSTools could read
# differently; each golden proves NewTNX agrees with WCSTools. Base: the
# syn_poly_full_watlongpole surfaces (WAT longpole 170).
_LW = 'wtype=tnx axtype=ra longpole=170. lngcor = "%s"' % _L
_BW = 'wtype=tnx axtype=dec latcor = "%s"' % _B


def _base(extra=(), lin=None, wat1=_LW, wat2=_BW, ctype=("RA---TNX", "DEC--TNX")):
    return header((800, 800), _cd800 if lin is None else lin, (210.5, 54.25),
                  (400.5, 400.5), wat1, wat2, extra=extra, ctype=ctype)


# PV2_m on the latitude axis: TAN has no latitude-axis parameters in
# wcslib, and WCSTools ignores PVi_m on a TNX header, so both ignore them.
CASES["syn_poly_pv2_ignored"] = _base([card("PV2_1", 1.5), card("PV2_2", 0.001)])
# CUNIT1/2 = 'deg' (the only CUNIT NewTNX accepts besides blank).
CASES["syn_poly_cunit_deg"] = _base([card("CUNIT1", "deg"), card("CUNIT2", "deg")])
# LONPOLE = 360: WCSTools tnxinit uses a LONPOLE card only when it is
# <= 360, so 360 itself is the largest value it reads.
CASES["syn_poly_lonpole_360"] = _base([card("LONPOLE", 360.0)])
# A negative LONPOLE beyond -360 is used as given by WCSTools.
CASES["syn_poly_lonpole_neg400"] = _base([card("LONPOLE", -400.0)])
# Keys that start WATn_ but do not end in digits are not WAT cards
# (WCSTools hgetm looks up WATn_<number> only): they must be ignored, not
# counted as stray WAT cards.
CASES["syn_poly_wat_nondigit_keys"] = _base(
    [card("WAT1_A", 'lngcor = "3. 1. 1. 0. -1 1 -1 1 0.5"'), card("WAT2_1X", "junk")])
# PROJPn: wcs.Parse reads it as PVi_n of the LATITUDE axis (wcslib
# WCSHDR_PROJPn), which TAN ignores; WCSTools does not read PROJPn on a
# TNX header. Accepted, with and without TNX surfaces.
CASES["syn_poly_projp1"] = _base([card("PROJP1", 0.5)])
CASES["syn_tan_no_wat_projp1"] = _base([card("PROJP1", 0.5)], wat1="", wat2="")
# Linear-step forms WCSTools wcsinit and wcslib read the same way:
# PC00i00j with CDELT (wcsinit reads PC001001 when PC1_1 is absent),
# PC1_1 present with PC2_2 absent (both default it), CD with only the
# off-diagonal terms (both default CD1_1/CD2_2 to 0), and CD beside
# CDELT (both ignore CDELT).
_cdelt = [("CDELT1", -0.0002), ("CDELT2", 0.0002)]
CASES["syn_poly_pc001001"] = _base(lin=_cdelt + [
    ("PC001001", 0.9), ("PC001002", 0.1), ("PC002001", -0.1), ("PC002002", 0.9)])
CASES["syn_poly_pc_no_pc2_2"] = _base(lin=_cdelt + [
    ("PC1_1", 0.9), ("PC1_2", 0.1), ("PC2_1", -0.1)])
CASES["syn_poly_cd_offdiag_only"] = _base(lin=[("CD1_2", 0.0002), ("CD2_1", 0.0002)])
CASES["syn_poly_cd_and_cdelt"] = _base(lin=_cd800 + [("CDELT1", -0.0003), ("CDELT2", 0.0003)])
# No CD, CDELT or PC at all: 1 deg/pixel in both readers (WCSTools
# wcsinit's no-scaling branch, wcslib's CDELT default 1). A 20 x 20 image
# keeps the field within 15 deg of the reference point.
CASES["syn_poly_no_scale"] = header((20, 20), [], (210.5, 54.25), (10.5, 10.5), _LW, _BW)
# CD00i00j beside its CDi_j with the same value: WCSTools reads CDi_j,
# wcs.Parse (wcslib) both, equal; only a CD00i00j WITHOUT its CDi_j is
# refused (refused/cd001001_only).
CASES["syn_poly_cd_and_cd001001"] = _base(lin=_cd800 + [
    ("CD001001", -0.0002), ("CD001002", 0.0), ("CD002001", 0.0), ("CD002002", 0.0002)])
# CD beside CDELT1 alone: the CDELT1/CDELT2 refusal applies only without
# a CD matrix (both readers use CD here).
CASES["syn_poly_cd_and_cdelt1"] = _base(lin=_cd800 + [("CDELT1", -0.0003)])
# One correction surface (lngcor only) plus PV2_m: tnxinit succeeds, so
# WCSTools ignores the PVi_m, as wcslib TAN does; the no-surface PV2_m
# refusal (refused/nowat_pv2) does not apply.
CASES["syn_poly_lngcor_only_pv2"] = _base(
    [card("PV2_1", 1.0), card("PV2_2", 0.001), card("PV2_4", 0.01)],
    wat2="wtype=tnx axtype=dec")


# --- Refused headers: refused/<name>.hdr. Each is a header WCSTools reads
# one way and the wcslib conventions (wcs.Parse / New) another, or one
# where WCSTools ignores a card, so NewTNX must return an error. Their
# refused/<name>.json goldens (gen_golden.py) record the WCSTools answer
# NewTNX gave something else for before the refusal existed.
REFUSED = {}
# PV1_1 / PV1_2 (phi0, theta0) without PV1_0: wcslib moves the fiducial
# point, WCSTools ignores PVi_m (about 1 deg apart at the reference).
REFUSED["pv1_fiducial"] = _base([card("PV1_1", 0.0), card("PV1_2", 89.0)])
# PV1_3 (LONPOLE as a PV card) with no WAT longpole: wcslib uses 150,
# WCSTools the TNX default 180.
REFUSED["pv1_3_lonpole"] = _base([card("PV1_3", 150.0)],
                                 wat1='wtype=tnx axtype=ra lngcor = "%s"' % _L)
# Both CD and PC: WCSTools uses CD, wcslib PC (Paper I forbids both).
REFUSED["cd_and_pc"] = _base([card("CDELT1", -0.0002), card("CDELT2", 0.0002),
                              card("PC1_1", 0.9), card("PC1_2", 0.1),
                              card("PC2_1", -0.1), card("PC2_2", 0.9)])
# LONPOLE > 360: WCSTools ignores the card (its "unset" sentinel) and uses
# the WAT longpole 170; wcslib uses the card.
REFUSED["lonpole_361"] = _base([card("LONPOLE", 361.0)])
REFUSED["lonpole_400"] = _base([card("LONPOLE", 400.0)])
# CUNIT other than deg: wcslib scales CDELT to degrees, WCSTools ignores
# CUNIT (CDELT read as degrees); TNX surfaces are in degrees of xi/eta.
REFUSED["cunit_arcmin"] = _base([card("CUNIT1", "arcmin"), card("CUNIT2", "arcmin"),
                                 card("CDELT1", -0.012), card("CDELT2", 0.012)],
                                lin=[])
# The same WAT key twice: WCSTools reads the first card and the second
# card's text is lost.
REFUSED["wat_duplicate_key"] = _base().replace(
    "WAT1_002= '0.09 0.09 3e-05 0.0012 -0.0008 0.02 \"'\n",
    "WAT1_002= '0.09 0.09 3e-05 0.0012 -0.0008 0.02 \"'\n"
    "WAT1_002= '0.09 0.09 3e-05 0.0012 -0.0008 0.05 \"'\n")
assert REFUSED["wat_duplicate_key"].count("WAT1_002") == 2
# Linear step without a CD matrix, where WCSTools wcsinit and wcslib
# disagree. wcsinit reads CDELT scaling only when CDELT1 is present
# (else 1 deg/pixel, CDELT2 and PC ignored), sets CDELT2 = -CDELT1 when
# CDELT2 is absent, reads PC only when PC1_1 (or PC001001) is present,
# and reads CDi_j only (not the defunct CD00i00j). wcslib: CDELTi default
# 1, PC applied whenever any PCi_j is present, CD00i00j read as CDi_j.
REFUSED["cdelt1_only"] = _base(lin=[("CDELT1", -0.0002)])
REFUSED["cdelt2_only"] = _base(lin=[("CDELT2", 0.0002)])
REFUSED["pc_no_cdelt"] = _base(lin=[("PC1_1", 0.0002), ("PC1_2", 0.00002),
                                    ("PC2_1", -0.00002), ("PC2_2", 0.0002)])
REFUSED["pc_no_pc1_1"] = _base(lin=_cdelt + [("PC1_2", 0.1), ("PC2_1", -0.1)])
REFUSED["cd001001_only"] = _base(lin=_cdelt + [
    ("CD001001", -0.00018), ("CD001002", 0.00002), ("CD002001", 0.00002), ("CD002002", 0.00018)])
# RA---TNX with no correction surface (no WAT1/WAT2, or WAT strings with
# neither lngcor nor latcor) plus latitude-axis PV2_m: WCSTools tnxinit
# fails, wcsinit demotes the projection to TAN and its TPV block applies
# the PVi_m (PV1_m default 0); wcslib TAN ignores PV2_m.
_pv2 = [card("PV2_1", 1.0), card("PV2_2", 0.001), card("PV2_4", 0.01)]
REFUSED["nowat_pv2"] = _base(_pv2, wat1="", wat2="")
REFUSED["watnocor_pv2"] = _base(_pv2, wat1="wtype=tnx axtype=ra", wat2="wtype=tnx axtype=dec")
# SIP keywords on a TNX header. With the -SIP suffix WCSTools applies SIP
# to the pixel offsets before the linear step and TNX (and sky2xy, with no
# AP/BP, inverts without SIP); without the suffix WCSTools ignores the
# cards; wcslib 8.10 applies SIP with TAN and drops the WAT correction.
_sip = [card("A_ORDER", 2), card("A_2_0", 2e-6), card("A_0_2", 1e-6),
        card("B_ORDER", 2), card("B_1_1", 1.5e-6)]
REFUSED["tnx_sip_suffix"] = _base(_sip, ctype=("RA---TNX-SIP", "DEC--TNX-SIP"))
REFUSED["tnx_sip_nosuffix"] = _base(_sip)
# Lowercase CTYPE: WCSTools 3.9.7 does not recognise 'RA---tnx' (it reads
# the axes as LINEAR and xy2sky prints asterisks), so there is no truth
# to prove an answer against.
REFUSED["ctype_lowercase"] = _base(ctype=("RA---tnx", "DEC--tnx"))


def main(argv):
    import glob
    import os
    import sys
    check = "--check" in argv
    bad = []
    for d, cases, pat in (("headers", CASES, "syn_*.hdr"), ("refused", REFUSED, "*.hdr")):
        for name, text in sorted(cases.items()):
            path = "%s/%s.hdr" % (d, name)
            if check:
                old = open(path).read() if os.path.exists(path) else None
                if old != text:
                    bad.append(path + (" (missing)" if old is None else " (differs)"))
                continue
            os.makedirs(d, exist_ok=True)
            with open(path, "w") as f:
                f.write(text)
        for path in sorted(glob.glob("%s/%s" % (d, pat))):
            name = os.path.splitext(os.path.basename(path))[0]
            if name not in cases:
                bad.append(path + " (not written by make_synthetic.py)")
    if bad:
        print("make_synthetic.py: synthetic headers out of date: " + ", ".join(bad))
        sys.exit(1)


if __name__ == "__main__":
    import sys
    main(sys.argv[1:])
