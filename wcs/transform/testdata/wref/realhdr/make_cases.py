#!/usr/bin/env python3
"""Build hdrcases/real_*.case from REAL archive headers (header cards only).

For each source below this script downloads the file (or uses a local copy
given with --src DIR), checks its sha256 against the pinned value, walks the
FITS HDUs with a minimal pure-Python reader (no astropy, so the cards are
copied byte for byte), takes the cards of the named HDU, drops the cards
listed in `drop` (each drop is recorded in the case file), and writes

    hdrcases/<name>.case          (or <casedir>/<name>.case)

with the case's grid directives followed by the verbatim header cards
(trailing blanks stripped; hdr2sky re-pads them to 80 columns). A source
may also name cards to REPLACE (`replace`: regex -> new card, recorded in
the case as "cards replaced"), extra directive lines (`extra`, e.g. the
"# SAME-AS:" line of an hdrcases_differ case) and the directory it is
written to (`casedir`, default hdrcases).

Usage (from wcs/transform/testdata/wref):

    python3 realhdr/make_cases.py [--src DIR] [--only NAME ...]
    make hdr-goldens

Re-running it reproduces the committed case files byte for byte (checked
by hand against the pinned sources; no test downloads them).
TestRealHdrCasesIntact only checks that each case's cards still match the
"header cards sha256" line written next to them, which catches an
accidental edit of the cards, not a deliberate edit of both.
Python >= 3.8, standard library only. See PROVENANCE.md for licenses.
"""

import argparse
import gzip
import hashlib
import os
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))

# The real ATLAS header (PR #4): an 80-column text file of WCS cards, no
# END card, two COMMENT cards naming the stamp.
ATLAS = dict(
    url="https://github.com/user-attachments/files/34297216/atlas_k2_tpv_sip_wcs.txt",
    local="atlas_k2_tpv_sip_wcs.txt",
    sha256="3b537ea0fa4cb12a1f41539bef0abfc02ecfe022fbfd64534600c031cb8a7d4b",
    hdu=0,
    text_header=True,
    license="WCS cards only (no image data) of an ATLAS survey product, posted by Daryl Manning "
            "on dmarkham/fits PR #4 (comment of 2026-11-03) as a test sample",
)
ATLAS_GRID = ["GRID 1 100 12 1 100 12", "PIX 50.5 50.5"]
SIP_CARD = r"^(A|B|AP|BP)_(ORDER|[0-9]_[0-9])"

# Directive comments of the hdrcases_differ cases (read by
# TestHdrDeliberateDifferences; hdr2sky skips '#' lines).
ATLAS_DIFFER = [
    "# SAME-AS: real_atlas_k2_tpv",
    "# WCSLIB-DIFFERS-ARCSEC: 24.4 26.9",
    "# ASTROPY-TOL-ARCSEC: 1e-8",
    "# ASTROPY-TOL-PX: 1e-9",
]
PTF_DIFFER = [
    "# SAME-AS: real_ptf_tpv4",
    "# WCSLIB-DIFFERS-ARCSEC: 72.6 133.4",
    "# ASTROPY-TOL-ARCSEC: 1e-8",
    "# ASTROPY-TOL-PX: 1e-9",
]

SOURCES = [
    dict(
        name="real_hst_acs_wfc_sip4",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/j94f05bgq_flt.fits",
        local="j94f05bgq_flt.fits",
        sha256="900038e0d853828140a757e2656934cb268ff9f315c5c6f617de85a632ad526b",
        hdu=1,
        drop=[],
        what="HST ACS/WFC FLT j94f05bgq, SCI,1: TAN-SIP A/B order 4, no AP/BP, CD matrix; "
             "also carries an alternate WCS 'O' (OPUS) the primary parse must ignore",
        license="astropy test data, BSD-3-Clause (header of a public HST/MAST product)",
        tol_px=("0.05", "4e-09", "2.5e-9"),
        grid=["GRID 1 4096 9 1 2048 7", "PIX 2048 1024"],
    ),
    dict(
        name="real_hst_wfc3_uvis_sip4",
        url="https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/ic3t33svq_flt.fits",
        local="ic3t33svq.part",
        partial=300001,
        sha256="73994d4498b679d4c0140bbee70b3b59d87fc49004396ad181fbee79f478a380",
        hdu=1,
        # Lookup-table distortions live in WCSDVARR/D2IMARR extensions,
        # which are not part of a header-cards-only sample. They are removed
        # here; wcs/transform rejects the unmodified header loudly
        # (TestNewRefusesUnimplementedDistortions).
        drop=[r"^D2IM", r"^CPDIS", r"^DP[12] ", r"^CPERR", r"^NPOL"],
        what="HST WFC3/UVIS FLT ic3t33svq, SCI,1: TAN-SIP A/B order 4, no AP/BP; "
             "lookup-table (CPDIS/D2IM) records removed",
        license="public HST data from MAST (STScI); header cards only",
        tol_px=("0.04", "5e-09", "2.2e-9"),
        grid=["GRID 1 4096 9 1 2051 7", "PIX 2048 1026"],
    ),
    dict(
        name="real_jwst_nircam_sip5",
        url="https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/jw01345001001_02201_00001_nrca1_cal.fits",
        local="jw01345001001_02201_00001_nrca1_cal.part",
        partial=400001,
        sha256="a69f669dda4da9cb16c8afe9c98fc0edae1a909c78e43afe391c5ea16b89ae2f",
        hdu=1,
        drop=[],
        what="JWST NIRCam NRCA1 cal jw01345001001_02201_00001 (CEERS), SCI: TAN-SIP A/B order 5 "
             "with approximate AP/BP order 5 (the jwst pipeline's SIP approximation), PC+CDELT",
        license="public JWST data from MAST (STScI); header cards only",
        tol_px=("0.031", "6e-09", "2.7e-9"),
        grid=["GRID 1 2048 9 1 2048 9", "PIX 1024.5 1024.5"],
    ),
    dict(
        name="real_jwst_nircam_sip4",
        url="https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/jw02736001001_02105_00001_nrca1_cal.fits",
        local="jw02736001001_02105_00001_nrca1_cal.part",
        partial=400001,
        sha256="6c4ddd65ab390ae070a7ddb4aed4fe2f26e7438c90b9bad99801ec72a0a3eafc",
        hdu=1,
        drop=[],
        what="JWST NIRCam NRCA1 cal jw02736001001_02105_00001 (ERO SMACS 0723), SCI: TAN-SIP "
             "A/B order 4 with approximate AP/BP order 4",
        license="public JWST data from MAST (STScI); header cards only",
        tol_px=("0.031", "6e-09", "3.0e-9"),
        grid=["GRID 1 2048 9 1 2048 9", "PIX 1024.5 1024.5"],
    ),
    dict(
        name="real_spitzer_irac_cryo_sip2",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/irac_sip.hdr",
        local="irac_sip.hdr",
        sha256="e5ceeee30461dd50b2d674c8c923bf2e22ad3cb3f365b5cfee78a82afc7ddd2f",
        hdu=0,
        raw_header=True,
        drop=[],
        what="Spitzer IRAC ch1 BCD (cryogenic mission, S8.9.0 pipeline): TAN-SIP A/B order 2 "
             "with AP/BP order 2, CD matrix",
        license="astropy test data, BSD-3-Clause (header of a public Spitzer/IRSA product)",
        grid=["GRID 1 256 9 1 256 9", "PIX 128 128"],
    ),
    dict(
        name="real_spitzer_irac_warm_sip5",
        url="https://irsa.ipac.caltech.edu/data/SPITZER/SHA/archive/proc/IRAC041600/r53849088/ch1/bcd/SPITZER_I1_53849088_0128_0000_1_cbcd.fits",
        local="irac.fits",
        sha256="c171fe84e010d6a4c2e7c89dd777bdafcd13d7f63d21252c80b031121604d51d",
        hdu=0,
        drop=[],
        what="Spitzer IRAC ch1 CBCD (warm mission, AOR 53849088, 2016): TAN-SIP A/B order 5 "
             "with AP/BP order 5",
        license="public Spitzer Heritage Archive data (IRSA/IPAC); header cards only",
        grid=["GRID 1 256 9 1 256 9", "PIX 128 128"],
    ),
    dict(
        name="real_anet_sip3",
        url="https://nova.astrometry.net/wcs_file/11126704",
        local="j11126704.fits",
        sha256="87fb0ac3da96cfc424014caae15580a4fd8bd0e4bbcc74c13eaadc44b92a13b2",
        hdu=0,
        drop=[],
        what="astrometry.net (nova) job 11126704 wcs.fits: TAN-SIP A/B order 3 with AP/BP order 3, "
             "2044x2044 image",
        license="nova.astrometry.net public job output (WCS solution only; no image data)",
        grid=["GRID 1 2044 9 1 2044 9", "PIX 1022.5 1022.5"],
    ),
    dict(
        name="real_anet_sip4",
        url="https://nova.astrometry.net/wcs_file/11087109",
        local="j11087109.fits",
        sha256="e148bc064534badece6b3d84e213c24d9b318a7964ace29f0ea3368856638671",
        hdu=0,
        drop=[],
        what="astrometry.net (nova) job 11087109 wcs.fits: TAN-SIP A/B order 4 with AP/BP order 4, "
             "3192x2128 image",
        license="nova.astrometry.net public job output (WCS solution only; no image data)",
        grid=["GRID 1 3192 9 1 2128 9", "PIX 1596.5 1064.5"],
    ),
    dict(
        name="real_anet_sip3_wide",
        url="https://nova.astrometry.net/wcs_file/11522654",
        local="j11522654.fits",
        sha256="8c63705e10e3d7a59514a114970c30681faf04ce08bf86990f093dc1502b09ad",
        hdu=0,
        drop=[],
        what="astrometry.net (nova) job 11522654 wcs.fits: TAN-SIP order 3 with AP/BP order 3, "
             "3920x1016 image at 2.5 arcsec/px (2.7 deg wide)",
        license="nova.astrometry.net public job output (WCS solution only; no image data)",
        grid=["GRID 1 3920 11 1 1016 5", "PIX 1960.5 508.5"],
    ),
    dict(
        name="real_ztf_tpv4",
        url="https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/2020/0331/486227/ztf_20200331486227_000535_zr_c11_o_q3_sciimg.fits",
        local="ztf_sciimg.fits",
        sha256="c516d7fe4292388352f366c02a91f5ae507056696e5b723e1095e2c07e52c89f",
        hdu=0,
        drop=[],
        what="ZTF sciimg field 535 zr c11 q3 (2020-03-31): SCAMP TPV, PV1_0..PV1_16 / PV2_0..PV2_16 "
             "(degree 4, no radial terms), CD matrix",
        license="public ZTF data (IRSA/IPAC); header cards only",
        grid=["GRID 1 3072 9 1 3080 9", "PIX 1536.5 1540.5"],
    ),
    dict(
        name="real_ptf_tpv4",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/tpvonly.hdr",
        local="tpvonly.hdr",
        sha256="ca3f753c704fadc7afb74b2c085afd4539c494271afbd212c93476a92025899d",
        hdu=0,
        raw_header=True,
        drop=[],
        what="PTF/MOSAIC CCD 5 (2014-07-31): TPV degree 4 with CRPIX ~4000 px off the chip, CD matrix",
        license="astropy test data, BSD-3-Clause (header of a public PTF/IRSA product)",
        grid=["GRID 1 2048 9 1 4096 9"],
    ),
    dict(
        # The SAME PTF exposure as real_ptf_tpv4 carries a second, SIP
        # solution (astropy's siponly.hdr; the two files differ only in
        # CTYPE, the PV cards and the SIP cards). This case is that real
        # pair merged the way astropy's test_tpv_ctype_tpv merges it: the
        # TPV header with the SIP cards of siponly.hdr appended, CTYPE
        # still RA---TPV (no -SIP). Since PR #4 (the real ATLAS header
        # below) a RA---TPV header's SIP keywords are ignored, as astropy
        # does, so this is a deliberate difference from wcslib
        # (hdrcases_differ): the answer must be real_ptf_tpv4's.
        name="real_ptf_tpv_plus_sip",
        casedir="hdrcases_differ",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/tpvonly.hdr",
        local="tpvonly.hdr",
        sha256="ca3f753c704fadc7afb74b2c085afd4539c494271afbd212c93476a92025899d",
        hdu=0,
        raw_header=True,
        drop=[],
        append=dict(
            url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/siponly.hdr",
            local="siponly.hdr",
            sha256="ac786a124934a3cd00f11921c6f2647b470f05ef04a2f4033960b408635a34a9",
            hdu=0,
            raw_header=True,
            select=r"^(A|B|AP|BP)_(ORDER|DMAX|[0-9]+_[0-9]+)",
        ),
        what="PTF/MOSAIC CCD 5 (2014-07-31): RA---TPV (no -SIP) TPV degree 4 PLUS the SIP order-4 "
             "A/B/AP/BP solution of the same exposure appended (wcslib applies both; we and astropy "
             "ignore the SIP set)",
        license="astropy test data, BSD-3-Clause (headers of a public PTF/IRSA product)",
        extra=PTF_DIFFER,
        grid=["GRID 1 2048 9 1 4096 9"],
    ),
    dict(
        # A REAL ATLAS header with the TPV + redundant SIP form: the WCS
        # cards of the ATLAS difference stamp 01a59742o0240o.diff (C/2017
        # K2, MLO, MJD 59742.348, 100x100 cutout), posted by Daryl Manning
        # on PR #4 (comment of 2026-11-03, attachment below). RA---TPV with
        # a TPV solution (CRPIX ~4000 px off the stamp) PLUS an order-5
        # A/B/AP/BP SIP set that is a copy of the same correction, no -SIP
        # suffix. wcslib 8.10 applies both (13.3-14.6 px off the image);
        # astropy 8.0.2 drops the SIP set; so do we (hdrcases_differ).
        name="real_atlas_k2_tpv_sip",
        casedir="hdrcases_differ",
        **ATLAS,
        drop=[],
        what="ATLAS difference stamp 01a59742o0240o.diff (C/2017 K2, MLO, MJD 59742.348, 100x100): "
             "RA---TPV (no -SIP) TPV degree 5 PLUS an order-5 A/B/AP/BP SIP copy of the same "
             "correction (wcslib applies both; we and astropy ignore the SIP set)",
        extra=ATLAS_DIFFER,
        grid=ATLAS_GRID,
    ),
    dict(
        # The ATLAS header with its SIP cards removed: TPV only. The answer
        # real_atlas_k2_tpv_sip must give (as astropy gives it).
        name="real_atlas_k2_tpv",
        **ATLAS,
        drop=[SIP_CARD],
        what="ATLAS difference stamp 01a59742o0240o.diff (C/2017 K2, MLO, MJD 59742.348, 100x100): "
             "the TPV solution alone (the SIP cards removed), CRPIX ~4000 px off the stamp",
        grid=ATLAS_GRID,
    ),
    dict(
        # The ATLAS header's SIP set on its own: CTYPE made RA---TAN-SIP,
        # the PV cards removed, CRPIX/CRVAL/CD unchanged. wcslib's reading
        # of this and of real_atlas_k2_tpv agree to < 1e-6 arcsec on the
        # stamp: the SIP set is a copy of the TPV correction
        # (TestATLASSIPCopiesTPV).
        name="real_atlas_k2_tansip",
        **ATLAS,
        drop=[r"^PV[12]_"],
        replace=[(r"^CTYPE1 ", "CTYPE1  = 'RA---TAN-SIP'"), (r"^CTYPE2 ", "CTYPE2  = 'DEC--TAN-SIP'")],
        what="ATLAS difference stamp 01a59742o0240o.diff (C/2017 K2, MLO, MJD 59742.348, 100x100): "
             "the SIP set alone (order-5 A/B/AP/BP) as RA---TAN-SIP, the PV cards removed",
        grid=ATLAS_GRID,
    ),
    dict(
        # A Digitized Sky Survey (STScI) cutout: the plate solution in the
        # DSS keywords (CNPIX, PPO, X/YPIXELSZ, PLTRA*/PLTDEC*, AMDX/AMDY)
        # together with an approximate TAN WCS (CRPIX/CRVAL/CDELT/CROTA,
        # CD and the defunct PC00i00j). wcslib wcspih replaces the
        # approximate WCS with the plate solution (Paper IV sect. 5.2).
        name="real_dss_plate_solution",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/dss.14.29.56-62.41.05.fits.gz",
        local="dss.14.29.56-62.41.05.fits.gz",
        sha256="d9e8a9adced55200b418517d5fba5f61a403842aa5eb57df9281eb38beb8d68a",
        gzip=True,
        hdu=0,
        drop=[],
        what="DSS (UK Schmidt plate S134 0025, STScI cutout at 14 29 56 -62 41 05): DSS plate "
             "solution AMDX/AMDY 1..13 plus the approximate TAN WCS (CD, PC00i00j, CROTA) "
             "wcslib replaces",
        license="astropy test data, BSD-3-Clause (header of a public STScI DSS cutout)",
        grid=["GRID 1 100 9 1 100 9", "GRID -2000 2000 9 -2000 2000 9", "PIX 50.5 50.5"],
    ),
    dict(
        name="real_rosat_glon_ait_lp0",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/visualization/wcsaxes/tests/data/rosat_header",
        local="rosat_header",
        sha256="4b693fb7d43dc7ee78f5031bda925de219715b3c00bc2ff1b300c4033c6976c1",
        hdu=0,
        text_header=True,
        drop=[],
        what="ROSAT all-sky survey map: GLON-AIT / GLAT-AIT, LONPOLE = 0, LATPOLE = 90, "
             "CDELT 0.675 deg; grid runs off the all-sky ellipse",
        license="astropy test data, BSD-3-Clause (header of a public ROSAT survey map)",
        grid=["GRID -20 500 9 -10 250 9", "PIX 240.5 120.5"],
    ),
    dict(
        name="real_msx_glon_car_lp0",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/visualization/wcsaxes/tests/data/msx_header",
        local="msx_header",
        sha256="a9cccf6d7bf2365e24583bf6f10b7f7e75f4a731d35b9bade7ba6b6651445b4c",
        hdu=0,
        text_header=True,
        drop=[],
        what="MSX Galactic-plane image: GLON-CAR / GLAT-CAR, LONPOLE = 0, LATPOLE = 90, "
             "CDELT 0.00667 deg",
        license="astropy test data, BSD-3-Clause (header of a public MSX survey image)",
        grid=["GRID 1 150 9 1 150 9", "PIX 75.907 74.8485"],
    ),
    dict(
        name="real_ps1_pc00i00j",
        url="https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/defunct_keywords.hdr",
        local="defunct_keywords.hdr",
        sha256="30d1197ab65e77a2846f1117cf03143daa634488291a14e7ca6a3c9e435caf2b",
        hdu=0,
        raw_header=True,
        drop=[],
        what="Pan-STARRS GPC1 skycell 0680.071 deep stack: RA---TAN with the defunct "
             "PC00i00j matrix keywords (PC001001 = -1) and CDELT, CRPIX far off the image",
        license="astropy test data, BSD-3-Clause (header of a public Pan-STARRS1 product)",
        grid=["GRID 1 720 9 1 720 9", "PIX 360.5 360.5"],
    ),
]


def read_text_header(b):
    """Cards of a newline-separated header text file (astropy test data
    style), each padded to 80 columns; stops at END."""
    cards = []
    for line in b.decode("latin-1").split("\n"):
        line = line.rstrip("\r")
        if line.rstrip() == "END":
            break
        if len(line) > 80:
            sys.exit("card longer than 80 columns: %r" % line)
        cards.append(line.ljust(80))
    return [cards]


def read_headers(b, raw=False):
    """Yield the card lists (END excluded) of each HDU in b."""
    if raw and len(b) % 80:
        # A raw header block whose final END card is not padded to 80
        # columns (astropy's defunct_keywords.hdr ends in "END").
        b = b + b" " * (80 - len(b) % 80)
    pos, out = 0, []
    while pos + 80 <= len(b):
        cards = []
        end = False
        while pos + 80 <= len(b):
            c = b[pos:pos + 80].decode("latin-1")
            pos += 80
            if c.startswith("END     ") or c.rstrip() == "END":
                end = True
                break
            cards.append(c)
        if not end:
            break
        out.append(cards)
        if raw:
            break
        pos = (pos + 2879) // 2880 * 2880
        kv = {}
        for c in cards:
            if c[8:10] == "= ":
                kv[c[:8].strip()] = c[10:].split("/")[0].strip()
        naxis = int(kv.get("NAXIS", "0"))
        n = 0
        if naxis > 0:
            n = 1
            for k in range(1, naxis + 1):
                n *= int(kv["NAXIS%d" % k])
        nbytes = abs(int(kv.get("BITPIX", "8"))) // 8 * n * int(kv.get("GCOUNT", "1")) + int(kv.get("PCOUNT", "0"))
        pos += (nbytes + 2879) // 2880 * 2880
    return out


def fetch(s, srcdir):
    """The bytes of source s (local copy under srcdir if present), checked
    against its pinned sha256."""
    if srcdir and os.path.exists(os.path.join(srcdir, s["local"])):
        data = open(os.path.join(srcdir, s["local"]), "rb").read()
    else:
        req = urllib.request.Request(s["url"])
        if s.get("partial"):
            req.add_header("Range", "bytes=0-%d" % (s["partial"] - 1))
        data = urllib.request.urlopen(req).read()
        if s.get("partial"):
            data = data[: s["partial"]]
    got = hashlib.sha256(data).hexdigest()
    if got != s["sha256"]:
        sys.exit("%s: sha256 %s != pinned %s" % (s.get("name", s["local"]), got, s["sha256"]))
    if s.get("gzip"):
        data = gzip.decompress(data)
    return data


def hdu_cards(s, data):
    if s.get("text_header"):
        hdrs = read_text_header(data)
    else:
        hdrs = read_headers(data, raw=s.get("raw_header", False))
    return hdrs[s["hdu"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", help="directory holding already downloaded sources")
    ap.add_argument("--only", action="append", help="regenerate only this case (repeatable)")
    args = ap.parse_args()
    for s in SOURCES:
        if args.only and s["name"] not in args.only:
            continue
        cards = hdu_cards(s, fetch(s, args.src))
        if s.get("append"):
            a = s["append"]
            cards = cards + [c for c in hdu_cards(a, fetch(a, args.src)) if re.match(a["select"], c)]
        kept, dropped, replaced = [], [], []
        for c in cards:
            if any(re.match(p, c) for p in s["drop"]):
                dropped.append(c[:8].rstrip())
                continue
            for p, new in s.get("replace", []):
                if re.match(p, c):
                    replaced.append("%s -> %s" % (c.rstrip(), new))
                    c = new
            line = c.rstrip()
            if line == "":
                continue
            kept.append(line)
        body = "\n".join(kept) + "\n"
        lines = [
            "# REAL header (header cards only), generated by realhdr/make_cases.py; do not edit.",
            "# " + s["what"],
            "# source: " + s["url"],
            "# license: " + s["license"],
        ]
        lines.append("# source sha256: " + s["sha256"] + ", HDU %d" % s["hdu"])
        if s.get("append"):
            a = s["append"]
            lines.append("# appended: the cards matching %s of %s" % (a["select"], a["url"]))
            lines.append("# appended source sha256: " + a["sha256"] + ", HDU %d" % a["hdu"])
        if dropped:
            lines.append("# cards removed: " + " ".join(sorted(set(dropped))))
        for r in replaced:
            lines.append("# card replaced: " + r)
        lines.append("# header cards sha256: " + hashlib.sha256(body.encode("latin-1")).hexdigest())
        if s.get("tol_px"):
            sc, floor, meas = s["tol_px"]
            lines += [
                "# TOL_PX 1e-8: the s2p input is wcslib's p2s sky position, whose angular",
                "# round-off (~2e-10 arcsec) is %s px at this header's %s arcsec/px plate" % (floor, sc),
                "# scale; measured Go vs wcslib max %s px. Known SIP inverse bugs were >= 1e-4 px." % meas,
                "TOL_PX 1e-8",
            ]
        lines += s.get("extra", [])
        lines.append("NOTE real: " + s["name"])
        lines += s["grid"]
        lines.append("HEADER")
        text = "\n".join(lines) + "\n" + body
        casedir = os.path.join(HERE, "..", s.get("casedir", "hdrcases"))
        with open(os.path.join(casedir, s["name"] + ".case"), "w", encoding="latin-1", newline="\n") as f:
            f.write(text)
        print("wrote", s["name"], len(kept), "cards", ("dropped " + ",".join(sorted(set(dropped)))) if dropped else "")


if __name__ == "__main__":
    main()
