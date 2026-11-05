#!/usr/bin/env python3
"""Regenerate every sample and golden in testdata/images.

Run through the Makefile (``make -C testdata/images all``), which builds the
cfitsio harness ./imgref and ../cfref/cfref first. Requires the repository
venv (astropy 8.0.2, numpy 2.4.6; see PROVENANCE.md).

Nothing here uses the Go library: samples are written by cfitsio (imgref),
by astropy, copied from astropy's own test data, or assembled byte by byte
(framing cases); goldens are what cfitsio reports about each sample.

Every output is a pure function of the inputs, so re-running reproduces the
committed bytes (``make verify`` checks this).
"""
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys

import astropy
import numpy as np
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
IMGREF = os.path.join(HERE, "imgref")
CFREF = os.path.join(HERE, "..", "cfref", "cfref")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE  # verify writes elsewhere

ASTROPY_DATA = os.path.join(os.path.dirname(fits.__file__), "tests", "data")

# Real files copied verbatim from astropy's io/fits test data (BSD-3).
# Their sha256 is recorded in golden/manifest.json; `make verify` fails if
# the installed astropy ships different bytes.
ASTROPY_FILES = [
    "arange.fits",          # 3-D int32 cube [11,10,7]
    "blank.fits",           # BITPIX 64 with BLANK
    "checksum.fits",        # valid CHECKSUM/DATASUM, image + bintable
    "checksum_false.fits",  # wrong CHECKSUM/DATASUM
    "fixed-1890.fits",      # uint16 (BZERO 32768) 100x100 image
    "group.fits",           # random groups
    "random_groups.fits",   # random groups (real UV data layout)
    "o4sp040b0_raw.fits",   # HST/STIS raw: real detector uint16 frames
    "scale.fits",           # BITPIX 16 with non-integer BSCALE
    "test0.fits",           # multi-extension int16, BZERO primary
    "memtest.fits",         # CHECKSUM without DATASUM (+ a bintable with wrong sums)
]


def run(args, **kw):
    r = subprocess.run(args, capture_output=True, text=True, **kw)
    if r.returncode != 0:
        sys.exit(f"{' '.join(args)} failed ({r.returncode}): {r.stderr}")
    return r.stdout


def compact(text):
    """JSON text without the whitespace between tokens (strings, numbers
    and every other byte kept exactly), plus a final newline: the goldens
    are stored compact, never pretty-printed."""
    out, in_str, esc = [], False, False
    for ch in text:
        if in_str:
            out.append(ch)
            if esc:
                esc = False
            elif ch == "\\":
                esc = True
            elif ch == '"':
                in_str = False
        elif ch == '"':
            in_str = True
            out.append(ch)
        elif ch not in " \t\r\n":
            out.append(ch)
    return "".join(out) + "\n"


GZIP_OVER = 16384  # golden text larger than this is committed gzipped


def put_golden(path, text):
    """Write golden text to path, or, over GZIP_OVER bytes, to path + ".gz"
    (gzip level 9, mtime 0). The Go tests read both (testutil.ReadGolden);
    make verify compares the decompressed bytes (cmp_golden.py)."""
    data = text.encode()
    for p in (path, path + ".gz"):
        if os.path.exists(p):
            os.remove(p)
    if len(data) > GZIP_OVER:
        open(path + ".gz", "wb").write(gzip.compress(data, compresslevel=9, mtime=0))
    else:
        open(path, "wb").write(data)


def write_golden(path, text):
    """Write harness JSON output compactly (the parsed value is unchanged)."""
    c = compact(text)
    assert json.loads(c) == json.loads(text), path
    put_golden(path, c)


def dump_golden(obj, path, sort_keys=False):
    put_golden(path, json.dumps(obj, sort_keys=sort_keys, separators=(",", ":")) + "\n")


def outpath(*p):
    path = os.path.join(OUT, *p)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return path


def card(s):
    assert len(s) <= 80, s
    return s.ljust(80).encode("ascii")


def header_bytes(cards):
    b = b"".join(card(c) for c in cards) + card("END")
    return b + b" " * (-len(b) % 2880)


def data_bytes(raw):
    return raw + b"\0" * (-len(raw) % 2880)


def gen_cfwrite():
    cases = json.load(open(os.path.join(HERE, "cases.json")))["cases"]
    status = {}
    for c in cases:
        dst = outpath("cfwrite", c["name"] + ".fits")
        opt = lambda v: "-" if v is None else repr(v) if isinstance(v, float) else str(v)
        args = [IMGREF, "write", dst, str(c["cf_imgtype"]), str(c["cf_datatype"]),
                opt(c["bscale"]), opt(c["bzero"]), opt(c["blank"]),
                ",".join(str(n) for n in c["naxes"])] + c["values"]
        status[c["name"]] = json.loads(run(args))
    return status


def gen_astropy_written():
    """Integer conventions as astropy writes them (uint=True default)."""
    w = lambda name, arr, **kw: fits.PrimaryHDU(arr, **kw).writeto(outpath("astropy_written", name), overwrite=True)
    w("i8.fits", np.array([-128, -127, -1, 0, 1, 126, 127], dtype=np.int8))
    w("u8.fits", np.array([0, 1, 127, 128, 255], dtype=np.uint8))
    w("u16.fits", np.array([[0, 1, 32767], [32768, 40000, 65535]], dtype=np.uint16))
    w("u32.fits", np.array([0, 1, 2147483647, 2147483648, 4294967295], dtype=np.uint32))
    w("u64.fits", np.array([0, 1, 2**63 - 1, 2**63, 2**64 - 2, 2**64 - 1], dtype=np.uint64))
    w("i64.fits", np.array([-2**63, -2**53 - 1, -1, 0, 2**62 + 1, 2**63 - 1], dtype=np.int64))
    w("f32_nan.fits", np.array([1.5, np.nan, -np.inf, 3.4028235e38], dtype=np.float32))
    cube = (np.arange(5 * 4 * 3, dtype=np.uint64) * 307445734561825860).astype(np.uint64).reshape(3, 4, 5)
    w("u64_cube.fits", cube)
    # A uint16 frame in the layout a camera driver writes (BZERO 32768):
    # the full 62x44 HST/STIS raw science frame (real detector counts)
    # re-written by astropy as a primary image, so the cards are astropy's.
    stis = fits.getdata(os.path.join(ASTROPY_DATA, "o4sp040b0_raw.fits"), 1)
    assert stis.dtype == np.uint16, stis.dtype
    hdu = fits.PrimaryHDU(np.ascontiguousarray(stis))
    hdu.header["ORIGIN"] = "astropy rewrite of o4sp040b0_raw.fits[1]"
    hdu.writeto(outpath("astropy_written", "u16_stis_frame.fits"), overwrite=True)
    # BSCALE/BZERO written by astropy's scale() (int16, non-integer scale)
    sc = fits.PrimaryHDU(np.linspace(-3.0, 9.0, 24).reshape(4, 6))
    sc.scale("int16", bscale=0.01, bzero=3.0)
    sc.writeto(outpath("astropy_written", "i16_bscale.fits"), overwrite=True)


def gen_astropy_data():
    shas = {}
    for name in ASTROPY_FILES:
        dst = outpath("astropy_data", name)
        shutil.copyfile(os.path.join(ASTROPY_DATA, name), dst)
        shas[name] = hashlib.sha256(open(dst, "rb").read()).hexdigest()
    return shas


def gen_framing():
    """Byte-assembled framing / structure cases (truth: cfitsio's reading)."""
    img = ["SIMPLE  =                    T", "BITPIX  =                   16",
           "NAXIS   =                    1", "NAXIS1  =                   10"]
    raw = b"".join(int(v).to_bytes(2, "big", signed=True) for v in range(-5, 5))
    base = header_bytes(img) + data_bytes(raw)
    ext = header_bytes(["XTENSION= 'IMAGE   '", "BITPIX  =                   16",
                        "NAXIS   =                    1", "NAXIS1  =                    3",
                        "PCOUNT  =                    0", "GCOUNT  =                    1"]) + \
        data_bytes(b"\x00\x07\x00\x08\x00\x09")
    w = lambda name, b: open(outpath("framing", name), "wb").write(b)
    w("trail_zero_block.fits", base + b"\0" * 2880)
    w("trail_junk_bytes.fits", base + b"junk\n")
    w("trail_junk_block.fits", base + b"not a FITS extension " * 137 + b"xxx")
    w("missing_final_pad.fits", header_bytes(img) + raw)
    w("missing_final_pad_ext.fits", base + ext[:2880 + 6])
    # NAXIS3 card although NAXIS = 1: cfitsio ignores NAXISn with n > NAXIS.
    stray = header_bytes(img + ["NAXIS3  =                    3"]) + data_bytes(raw)
    w("stray_naxis3.fits", stray + ext)
    # The stray NAXIS3 twice (different values): still ignored, not a
    # duplicate axis (cfitsio reads only NAXIS1..NAXIS).
    w("stray_naxis3_twice.fits", header_bytes(img + ["NAXIS3  =                    3",
                                                     "NAXIS3  =                    5"]) + data_bytes(raw))
    # A malformed non-structural card (unquoted string value): cfitsio
    # opens the file and reads the pixels.
    w("malformed_card.fits", header_bytes(img + ["DATE-OBS= 2020-01-01T00:00:00"]) + data_bytes(raw))
    # Random groups primary followed by an IMAGE extension.
    rg = ["SIMPLE  =                    T", "BITPIX  =                  -32",
          "NAXIS   =                    3", "NAXIS1  =                    0",
          "NAXIS2  =                    3", "NAXIS3  =                    2",
          "GROUPS  =                    T", "PCOUNT  =                    2",
          "GCOUNT  =                    4"]
    rgdata = np.arange(4 * (2 + 3 * 2), dtype=">f4").tobytes()
    w("random_groups_ext.fits", header_bytes(rg) + data_bytes(rgdata) + ext)
    # Legacy XTENSION = 'IUEIMAGE', which cfitsio reads as IMAGE_HDU.
    iue = ext.replace(b"XTENSION= 'IMAGE   '", b"XTENSION= 'IUEIMAGE'", 1)
    assert iue != ext
    w("xtension_iueimage.fits", base + iue)
    # A conforming extension type the library does not decode (FOREIGN,
    # registered with the IAU FITS WG) between the primary and an IMAGE:
    # its size is BITPIX * GCOUNT * (PCOUNT + NAXIS1*...) like any
    # extension (FITS 4.0 eq. 2), so the IMAGE after it is found at the
    # same address cfitsio finds it. cfitsio reads the FOREIGN HDU as an
    # IMAGE_HDU ("unknown type of extension; press on anyway").
    foreign = header_bytes(["XTENSION= 'FOREIGN '", "BITPIX  =                    8",
                            "NAXIS   =                    1", "NAXIS1  =                 3000",
                            "PCOUNT  =                  100", "GCOUNT  =                    1",
                            "EXTNAME = 'BLOB    '"]) + data_bytes(bytes(range(256)) * 12 + bytes(28))
    w("xtension_foreign.fits", base + foreign + ext)
    # An IMAGE extension without PCOUNT or GCOUNT: cfitsio accepts it
    # (defaults 0 / 1), unlike a table extension (see gen_refused).
    for k in ("PCOUNT", "GCOUNT"):
        cards = ["XTENSION= 'IMAGE   '", "BITPIX  =                   16",
                 "NAXIS   =                    1", "NAXIS1  =                    3",
                 "PCOUNT  =                    0", "GCOUNT  =                    1"]
        w(f"img_no_{k.lower()}.fits", base + header_bytes([c for c in cards if not c.startswith(k)]) +
          data_bytes(b"\x00\x07\x00\x08\x00\x09"))


def gen_refused():
    """Byte-assembled files the library refuses with a typed error (truth:
    cfitsio's status for each, from `imgref open`). TestRefusedFramingMatchesCfitsio."""
    img = ["SIMPLE  =                    T", "BITPIX  =                   16",
           "NAXIS   =                    1", "NAXIS1  =                   10"]
    raw = b"".join(int(v).to_bytes(2, "big", signed=True) for v in range(-5, 5))
    base = header_bytes(img) + data_bytes(raw)
    ext = ["XTENSION= 'IMAGE   '", "BITPIX  =                   16",
           "NAXIS   =                    1", "NAXIS1  =                    3",
           "PCOUNT  =                    0", "GCOUNT  =                    1"]
    extdata = data_bytes(b"\x00\x07\x00\x08\x00\x09")
    w = lambda name, b: open(outpath("refused", name), "wb").write(b)
    # Declared data runs past EOF: 2000 int16 = 4000 bytes, 2880 present.
    w("truncated_data.fits", header_bytes(img[:3] + ["NAXIS1  =                 2000"]) + b"\x00\x01" * 1440)
    # The extension's 6 data bytes cut to 2.
    w("truncated_data_ext.fits", base + header_bytes(ext) + b"\x00\x07")
    # Mandatory keywords missing.
    w("missing_naxis2.fits", header_bytes(img[:2] + ["NAXIS   =                    2", img[3]]) + data_bytes(raw))
    w("missing_bitpix.fits", header_bytes([img[0]] + img[2:]) + data_bytes(raw))
    w("missing_naxis.fits", header_bytes(img[:2] + [img[3]]) + data_bytes(raw))
    w("ext_missing_naxis1.fits", base + header_bytes([c for c in ext if not c.startswith("NAXIS1")]) + extdata)
    w("ext_missing_bitpix.fits", base + header_bytes([c for c in ext if not c.startswith("BITPIX")]) + extdata)
    # XTENSION with a non-string value (cfitsio reads it as an IMAGE).
    w("ext_xtension_not_string.fits", base + header_bytes(["XTENSION=                    1"] + ext[1:]) + extdata)
    # A table extension without PCOUNT, GCOUNT or TFIELDS (cfitsio will not
    # move to it: 228 NO_PCOUNT, 229 NO_GCOUNT, 230 NO_TFIELDS). An IMAGE
    # extension without PCOUNT/GCOUNT is accepted by cfitsio and by us.
    tbl = ["XTENSION= 'BINTABLE'", "BITPIX  =                    8",
           "NAXIS   =                    2", "NAXIS1  =                    4",
           "NAXIS2  =                    2", "PCOUNT  =                    0",
           "GCOUNT  =                    1", "TFIELDS =                    1",
           "TFORM1  = '1J      '", "TTYPE1  = 'A       '"]
    tbldata = data_bytes(b"\x00\x00\x00\x07\x00\x00\x00\x09")
    for k in ("PCOUNT", "GCOUNT", "TFIELDS"):
        w(f"tbl_no_{k.lower()}.fits", base + header_bytes([c for c in tbl if not c.startswith(k)]) + tbldata)
    asc = ["XTENSION= 'TABLE   '", "BITPIX  =                    8",
           "NAXIS   =                    2", "NAXIS1  =                    3",
           "NAXIS2  =                    2", "PCOUNT  =                    0",
           "GCOUNT  =                    1", "TFIELDS =                    1",
           "TFORM1  = 'I3      '", "TBCOL1  =                    1", "TTYPE1  = 'A       '"]
    w("ascii_no_tfields.fits", base + header_bytes([c for c in asc if not c.startswith("TFIELDS")]) + b"  7  9" + b" " * 2874)


# Image scaling keywords as cfitsio reads them (ffgphd, getkey.c): every
# card is visited and the LAST BSCALE / BZERO / BLANK wins; the value token
# is cut by ffpsvc ('=' in column 9, else the first '=' on the card) and
# converted with ffc2dd (strtod after replacing the first uppercase 'D';
# junk after the number, ERANGE, Inf or NaN -> the default 1 / 0) or
# ffc2jj (strtoll base 10; failure -> no BLANK). Each sample is an int16
# (or int32) image whose stored values start 0, 1 so the BZERO and BSCALE
# cfitsio applies can be read off its pixels. Truth: cfref dump (TDOUBLE
# with null flags) and imgref pix (TFLOAT, and TDOUBLE without null checks).
KEYWORD_RAW16 = [0, 1, 2, 3, -4, 7]
KEYWORD_RAW32 = [0, 1, 2, 1234554321, -4, 7]
KEYWORD_CASES = [
    # name, extra cards, bitpix
    ("dexp_upper", ["BSCALE  =                1.5D0", "BZERO   =               2.0D+1"], 16),
    ("dexp_lower", ["BSCALE  =                1.5d0", "BZERO   =               2.0d+1"], 16),
    ("dexp_hex", ["BSCALE  =                 0x1D"], 16),
    # ffc2dd replaces only the FIRST 'D' (strchr): 0xDD reads as 0xED = 237
    # (replacing both would give 0xEE = 238)
    ("dexp_hex_two_d", ["BSCALE  =                 0xDD"], 16),
    ("dup_bscale", ["BSCALE  =                    2", "BSCALE  =                    3"], 16),
    ("dup_bzero", ["BZERO   =                   10", "BZERO   =                   20"], 16),
    ("dup_blank", ["BLANK   =                    1", "BLANK   =                    2"], 16),
    ("dup_bscale_bad_last", ["BSCALE  =                    2", "BSCALE  = 'x'"], 16),
    ("dup_bscale_bad_first", ["BSCALE  = 'x'", "BSCALE  =                    3"], 16),
    ("dup_blank_bad_last", ["BLANK   =                    1", "BLANK   =                  1.5"], 16),
    ("bscale_string", ["BSCALE  = '2.0'"], 16),
    ("bzero_string", ["BZERO   = '5'"], 16),
    ("bscale_inf", ["BSCALE  =                  inf"], 16),
    ("bscale_neg_infinity", ["BSCALE  =            -Infinity"], 16),
    ("bscale_nan", ["BSCALE  =                  NaN"], 16),
    ("bzero_nan_payload", ["BZERO   =               nan(1)"], 16),
    ("bscale_overflow", ["BSCALE  =               1.0E400"], 16),
    ("bscale_underflow", ["BSCALE  =              1.0E-400"], 16),
    ("bscale_subnormal", ["BSCALE  =              1.0E-310"], 16),
    ("bscale_subnormal_exact", ["BSCALE  =            0x1p-1074"], 16),
    ("bscale_noeq_space", ["BSCALE  =2.0"], 16),
    ("bscale_key_eq", ["BSCALE=2.0"], 16),
    # '= ' in columns 9-10 sets the value position even with an earlier
    # '=' (ffpsvc): the value is 2, not "1=" (from the first '=')
    ("bscale_eq_col9_after_eq", ["BSCALE=1= 2"], 16),
    ("bscale_hex", ["BSCALE  =                0x1p1"], 16),
    ("bscale_hex_noexp", ["BSCALE  =                  0x3"], 16),
    ("bscale_empty", ["BSCALE  =", "BZERO   =                    5"], 16),
    ("bscale_slash", ["BSCALE  = / no value", "BZERO   =                    5"], 16),
    ("bscale_no_equals", ["BSCALE    2.0", "BZERO   =                    5"], 16),
    ("bscale_comment_equals", ["BSCALE    2.0 / x=3"], 16),
    ("bscale_lower", ["bscale  =                  2.0", "bzero   =                  5.0"], 16),
    ("bscale_hierarch", ["HIERARCH BSCALE = 2.0", "HIERARCH BZERO = 5.0"], 16),
    ("bscale_plus", ["BSCALE  =                 +2.5"], 16),
    ("bscale_leading_dot", ["BSCALE  =                   .5"], 16),
    ("bscale_trailing_dot", ["BSCALE  =                   2."], 16),
    ("bscale_exponent_only", ["BSCALE  =                  1E1"], 16),
    ("bscale_trailing_e", ["BSCALE  =                 2.0E"], 16),
    ("bscale_logical", ["BSCALE  =                    T"], 16),
    ("bscale_complex", ["BSCALE  = (2.0, 0.0)"], 16),
    ("bscale_junk", ["BSCALE  =                 2.0x"], 16),
    ("bscale_zero", ["BSCALE  =                  0.0", "BZERO   =                  5.0"], 16),
    # glibc's underflow test is "tiny after rounding to 53 bits": the first
    # rounds up to 2^-1022 and is accepted, the second is below it (ERANGE)
    ("bscale_min_normal_rounded", ["BSCALE  = 2.2250738585072013e-308"], 16),
    ("bscale_below_min_normal", ["BSCALE  = 2.2250738585072012e-308"], 16),
    ("blank_plus", ["BLANK   =                   +2"], 16),
    ("blank_float", ["BLANK   =                  2.0"], 16),
    ("blank_empty", ["BLANK   ="], 16),
    ("blank_overflow", ["BLANK   = 99999999999999999999"], 16),
    ("blank_string", ["BLANK   = '2'"], 16),
    ("blank_hex", ["BLANK   =                  0x2"], 16),
    ("blank_wrap16", ["BLANK   =                65538"], 16),
    ("blank_wrap32", ["BLANK   =           4294967298"], 32),
    ("blank_magic", ["BLANK   =           1234554321"], 32),
]


def gen_keywords():
    for name, extra, bp in KEYWORD_CASES:
        raw = KEYWORD_RAW16 if bp == 16 else KEYWORD_RAW32
        cards = ["SIMPLE  =                    T", f"BITPIX  = {bp:20d}",
                 "NAXIS   =                    1", f"NAXIS1  = {len(raw):20d}"] + extra
        data = b"".join(int(v).to_bytes(bp // 8, "big", signed=True) for v in raw)
        open(outpath("keywords", name + ".fits"), "wb").write(header_bytes(cards) + data_bytes(data))


# OverwritePixels truth: cfitsio's fits_write_img over the primary image of
# each keywords/ sample, scaled with the BSCALE/BZERO/BLANK cfitsio read
# from that header (imgref overwrite). golden/keywords_overwrite.json: the
# values, cfitsio's status and the written data unit (hex, unpadded).
OVERWRITE_VALUES = ["5", "6", "-7", "100", "0", "1"]


def gen_keyword_overwrite():
    out = {}
    tmp = outpath("golden", "_overwrite.fits")
    for name, _, bp in KEYWORD_CASES:
        dt = 21 if bp == 16 else 31
        r = json.loads(run([IMGREF, "overwrite", outpath("keywords", name + ".fits"), tmp, str(dt)] + OVERWRITE_VALUES))
        b = open(tmp, "rb").read()
        assert (b"END" + b" " * 77) in [b[i:i + 80] for i in range(0, 2880, 80)], name  # one header block
        data = b[2880:2880 + len(OVERWRITE_VALUES) * bp // 8]
        out[name + ".fits"] = dict(datatype=dt, values=OVERWRITE_VALUES, status=r["status"], data=data.hex())
    os.remove(tmp)
    dump_golden(out, outpath("golden", "keywords_overwrite.json"), sort_keys=True)


# Duplicate structural keywords. cfitsio reads BITPIX, NAXIS and NAXISn by
# position (cards 2, 3, 4..) and PCOUNT/GCOUNT of an image in the ffgphd
# loop (last wins); the library refuses every such header instead
# (TestDuplicateStructuralKeywordsRefused). Truth recorded: cfref dump.
def gen_duplicate():
    S = "SIMPLE  =                    T"
    img = [S, "BITPIX  =                   16", "NAXIS   =                    1", "NAXIS1  =                    4"]
    raw = data_bytes(b"\x00\x01\x00\x02\x00\x03\x00\x04")
    w = lambda name, b: open(outpath("duplicate", name), "wb").write(b)
    w("dup_bitpix.fits", header_bytes(img + ["BITPIX  =                    8"]) + raw)
    w("dup_naxis.fits", header_bytes(img + ["NAXIS2  =                    1", "NAXIS   =                    2"]) + raw)
    w("dup_naxis1_smaller.fits", header_bytes(img + ["NAXIS1  =                    2"]) + raw)
    w("dup_naxis1_larger.fits", header_bytes(img[:3] + ["NAXIS1  =                    2", "NAXIS1  =                    4"]) + raw)
    w("dup_simple.fits", header_bytes(img + ["SIMPLE  =                    F"]) + raw)
    # a real duplicate axis (NAXIS1) followed by a duplicated stray NAXIS3:
    # the NAXIS1 duplicate must still be refused
    w("dup_naxis1_and_stray_naxis3.fits", header_bytes(img + ["NAXIS1  =                    2",
      "NAXIS3  =                    1", "NAXIS3  =                    1"]) + raw)
    base = header_bytes(img) + raw
    ext = ["XTENSION= 'IMAGE   '", "BITPIX  =                   16", "NAXIS   =                    1",
           "NAXIS1  =                    3", "PCOUNT  =                    0", "GCOUNT  =                    1"]
    extdata = data_bytes(b"\x00\x07\x00\x08\x00\x09")
    w("dup_pcount_ext.fits", base + header_bytes(ext + ["PCOUNT  =                    2880"]) + extdata)
    w("dup_gcount_ext.fits", base + header_bytes(ext + ["GCOUNT  =                    2"]) + extdata)
    w("dup_xtension_ext.fits", base + header_bytes(ext + ["XTENSION= 'BINTABLE'"]) + extdata)
    tbl = ["XTENSION= 'BINTABLE'", "BITPIX  =                    8",
           "NAXIS   =                    2", "NAXIS1  =                    4",
           "NAXIS2  =                    2", "PCOUNT  =                    0",
           "GCOUNT  =                    1", "TFIELDS =                    1",
           "TFORM1  = '1J      '", "TTYPE1  = 'A       '", "TFIELDS =                    2"]
    w("dup_tfields_ext.fits", base + header_bytes(tbl) + data_bytes(b"\x00\x00\x00\x07\x00\x00\x00\x09"))


# WriteImage refusals with a cfitsio counterpart: the scaling cfitsio's
# fits_set_bscale refuses (imgref write reports the status).
WRITE_REFUSED = [
    # name, imgtype, datatype, bscale, bzero
    ("i16_bscale0", 16, 21, "0", "5"),
    ("i16_bscale0_nobzero", 16, 21, "0", "0"),
    ("u8_bscale0", 8, 11, "0", "0"),
    ("i64_bscale0", 64, 81, "0", "0"),
    ("f32_bscale0", -32, 42, "0", "1"),
]


def gen_write_refused():
    out = {}
    tmp = outpath("golden", "_refused_write.fits")
    for name, it, dt, bscale, bzero in WRITE_REFUSED:
        r = json.loads(run([IMGREF, "write", tmp, str(it), str(dt), bscale, bzero, "-", "4", "1", "2", "3", "4"]))
        if os.path.exists(tmp):
            sys.exit(f"write_refused {name}: cfitsio left a file")
        out[name] = dict(imgtype=it, datatype=dt, bscale=bscale, bzero=bzero, **r)
    dump_golden(out, outpath("golden", "write_refused.json"), sort_keys=True)


# Values WriteImage refuses with ErrNumOverflow (writing nothing). Each
# value is written ALONE by cfitsio (imgref write: fits_create_img +
# BSCALE/BZERO + fits_write_img) with the case's scaling, and its status
# (412 NUM_OVERFLOW or 0) and stored bytes recorded, so the golden says,
# per value, whether cfitsio overflowed and what it stored. Integer cases
# are the values of the former clamp cases of cases.json (which now keep
# only their in-range values); the float cases scale finite values to
# infinity, which cfitsio stores without an error.
# golden/write_overflow.json. TestWriteImageOverflowMatchesCfitsio.
WRITE_OVERFLOW = [
    # name, gotype, imgtype, datatype, bscale, bzero, values
    ("u8_bzero_neg100", "uint8", 8, 11, "1", "-100",
     ["0", "1", "100", "154", "155", "156", "200", "255"]),
    ("u8_bzero_10", "uint8", 8, 11, "1", "10",
     ["0", "5", "9", "10", "11", "200", "255"]),
    ("i16_bzero1000", "int16", 16, 21, "1", "1000",
     ["-32768", "-31769", "-31768", "-1", "0", "1000", "32767"]),
    ("i16_bzero_only", "int16", 16, 21, None, "1000",
     ["1", "2", "3", "4", "-32768"]),
    ("i32_bscale_milli", "int32", 32, 31, "0.001", "0",
     ["-2147484", "-2147483", "-1", "0", "1", "2147483", "2147484", "1000000000"]),
    ("i64_bzero1000", "int64", 64, 81, "1", "1000",
     ["-9223372036854775808", "-9223372036854774808", "-1", "0", "1000", "123456789012"]),
    ("i64_bscale_half", "int64", 64, 81, "0.5", "0.25",
     ["-9223372036854775808", "-4611686018427387904", "4611686018427387904", "4611686018427387903",
      "9223372036854775807", "3", "-3", "0", "1", "-1", "9007199254740993", "-9007199254740995"]),
    ("i16_bscale_neg1", "int16", 16, 21, "-1", "0",
     ["-32768", "-32767", "-2", "-1", "0", "1", "2", "32767"]),
    ("u8_bzero_half", "uint8", 8, 11, "1", "0.5",
     ["0", "1", "2", "255"]),
    ("i16_bzero_half", "int16", 16, 21, "1", "0.5",
     ["-32768", "-32767", "0", "1", "32767"]),
    ("i32_bzero_half", "int32", 32, 31, "1", "0.5",
     ["-2147483648", "-2147483647", "0", "2147483647"]),
    ("u16_bscale_1p1", "uint16", 16, 20, "1.1", "0.25",
     ["0", "14", "25", "58", "102", "113", "36044", "65535"]),
    ("u32_bscale_1p1", "uint32", 32, 30, "1.1", "0.25",
     ["0", "14", "25", "2147483649", "2147483660", "4294967295"]),
    ("u64_bscale_1p1", "uint64", 64, 80, "1.1", "0.25",
     ["0", "14", "25", "58", "18446744073709551615"]),
    ("u64_bscale_neg1p5", "uint64", 64, 80, "-1.5", "0",
     ["0", "1", "2", "3", "18446744073709551615"]),
    ("i8_bscale_neg1", "int8", 8, 12, "-1", "0",
     ["-128", "-1", "0", "1", "127"]),
    ("f32_bscale_tiny", "float32", -32, 42, "1e-30", "0",
     ["0", "1", "300000000", "400000000", "-400000000", "inf", "nan"]),
    ("f64_bscale_tiny", "float64", -64, 82, "1e-300", "0",
     ["0", "1", "100000000", "1000000000", "-1000000000", "-inf"]),
]


def gen_write_overflow():
    out = {}
    tmp = outpath("golden", "_overflow.fits")
    for name, gotype, it, dt, bscale, bzero, values in WRITE_OVERFLOW:
        size = abs(it) // 8
        status, data = [], []
        for v in values:
            r = json.loads(run([IMGREF, "write", tmp, str(it), str(dt), bscale or "-", bzero or "-", "-", "1", v]))
            b = open(tmp, "rb").read()
            end = [i for i in range(0, 2880, 80) if b[i:i + 80] == b"END" + b" " * 77]
            assert end, name  # one header block
            status.append(r["status"])
            data.append(b[2880:2880 + size].hex())
        out[name] = dict(gotype=gotype, imgtype=it, datatype=dt, bscale=bscale, bzero=bzero,
                         values=values, status=status, data=data)
    os.remove(tmp)
    dump_golden(out, outpath("golden", "write_overflow.json"), sort_keys=True)


def gen_strtod():
    """strtod / strtoll truth (imgref strtod, i.e. the C library cfitsio's
    ffc2dd / ffc2jj call) for a fixed list of boundary inputs plus seeded
    random tokens and numbers. golden/strtod.txt: "input<TAB>strtod<TAB>
    strtoll" per line (inputs never contain a tab or newline)."""
    import random
    fixed = ["", ".", "+", "-", "1e", "1e+", "1E5", "+.5", "-0", "00012", "  12", "\x0b3", "1_0", "2.0x",
             "1.5d0", "1.5D0", "inf", "INFINITY", "infin", "-Inf", "nan", "NaN", "nan(123)", "nan()", "nan(_a)",
             "nan(a b)", "nanx", "-nan", "0x", "0x.", "0x.8", "0x.p1", "0x1.", "0x1p1", "0x2", "0X1P3", "0x1_0p0",
             "0x1p-1074", "0x1p-1075", "0x3p-1076", "0x1.8p-1074", "0x1p-1022", "0x1.0000000000001p-1075",
             "0x0.fffffffffffff8p-1022", "0x0.fffffffffffffcp-1022", "0x1.fffffffffffffp-1023", "0x1p1023",
             "0x1p1024", "0x0p-99999", "1e-310", "2e-324", "4.9406564584124654e-324", "2.4703282292062328e-324",
             "2.4703282292062327e-324", "1e-400", "0e-400", "0e400", "1e400", "-1e-320",
             "2.2250738585072011e-308", "2.2250738585072012e-308", "2.2250738585072013e-308",
             "2.2250738585072014e-308", "2.225073858507201e-308", "1.7976931348623157e308",
             "1.7976931348623158e308", "1.7976931348623159e308", "1.797693134862315807e308",
             "1e-99999999999999999999", "1e99999999999999999999",
             "0.000000000000000000000000000000000000000000000000000000000000001e-260",
             "9223372036854775807", "9223372036854775808", "-9223372036854775808", "-9223372036854775809",
             "99999999999999999999", "+2", "0x2", "2.0", "1234554321", "4294967298", "'2'", "T", "(2.0, 0.0)"]
    rnd = random.Random(20261101)
    alpha = "0123456789.eEpPxX+-dDinfaINFAtyNy() _"
    toks = [''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 12))) for _ in range(1500)]
    nums = []
    for _ in range(750):
        if rnd.random() < 0.5:
            f = np.frombuffer(rnd.getrandbits(64).to_bytes(8, "little"), dtype="<f8")[0]
        else:
            f = np.ldexp(rnd.random(), -1000 - rnd.randint(0, 80))
        if not np.isfinite(f):
            f = 1.0
        nums.append(f"{float(f):.{rnd.randint(0, 19)}e}")
        nums.append(float(f).hex())
    inputs = fixed + toks + nums
    gen_c2(fixed + toks)
    assert all("\t" not in x and "\n" not in x for x in inputs)
    r = subprocess.run([IMGREF, "strtod"], input="\n".join(inputs) + "\n", capture_output=True, text=True, check=True)
    outs = r.stdout.rstrip("\n").split("\n")
    assert len(outs) == len(inputs), (len(outs), len(inputs))
    put_golden(outpath("golden", "strtod.txt"), "".join(i + "\t" + o + "\n" for i, o in zip(inputs, outs)))


def gen_c2(base):
    """cfitsio's ffc2dd / ffc2jj called directly (imgref c2) on the strtod
    boundary inputs and random tokens plus value strings that exercise the
    wrappers themselves: a trailing blank or blank-separated junk
    (accepted), several 'D's (only the first is replaced), a 'D' string
    longer than 72 characters (refused), Inf/NaN (refused).
    golden/c2.txt: "input<TAB>ffc2dd status bits<TAB>ffc2jj status value"."""
    extra = ["2 ", "2 x", "2x", "-3 ", "0xDD", "1D1D", "1DD", "1.5D0 ", "1.5D0x", "D", "1D", "0x1D",
             "1" * 70 + "D0", "1" * 71 + "D0", "1" * 80, "+2 ", "7 8", "1234554321 ", "9223372036854775808 "]
    inputs = base + extra
    assert all("\t" not in x and "\n" not in x for x in inputs)
    r = subprocess.run([IMGREF, "c2"], input="\n".join(inputs) + "\n", capture_output=True, text=True, check=True)
    outs = r.stdout.rstrip("\n").split("\n")
    assert len(outs) == len(inputs), (len(outs), len(inputs))
    put_golden(outpath("golden", "c2.txt"), "".join(i + "\t" + o + "\n" for i, o in zip(inputs, outs)))


def gen_datasum():
    """DATASUM values cfitsio parses with atof + (unsigned long) (ffvcks):
    leading blanks, a fraction (truncated), a number followed by junk, an
    exponent, a hex float, a value 2^32 too large, a negative fraction on a
    header-only HDU (data sum 0, so -0.5 truncates to a match), -1 and
    -2^32 there. Truth: imgref chksum (fits_verify_chksum)."""
    src = outpath("cfwrite", "i16_plain.fits")
    good = open(src, "rb").read()
    ds = json.loads(run([IMGREF, "chksum", src]))["hdus"][0]["datasum"]
    hdr, data = good[:2880], good[2880:]
    cards = [hdr[i:i + 80].decode() for i in range(0, 2880, 80)]
    end = cards.index("END".ljust(80))
    def with_datasum(value):
        c = cards[:end] + [("DATASUM = " + value).ljust(80)]
        return header_bytes([x.rstrip() for x in c]) + data
    variants = {
        "datasum_leading_blanks": f"'   {ds}'",
        "datasum_fraction": f"'{ds}.9'",
        "datasum_junk": f"'{ds}abc'",
        "datasum_exponent": f"'{ds}e0'",
        "datasum_hex": f"'{hex(ds)}'",
        "datasum_plus_2_32": f"'{ds + 2**32}'",
        "datasum_unquoted": f"{ds}",
    }
    for name, v in variants.items():
        open(outpath("checksum", name + ".fits"), "wb").write(with_datasum(v))
    ho = ["SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0"]
    open(outpath("checksum", "datasum_neg_fraction_empty.fits"), "wb").write(header_bytes(ho + ["DATASUM = '-0.5'"]))
    open(outpath("checksum", "datasum_neg_one_empty.fits"), "wb").write(header_bytes(ho + ["DATASUM = '-1'"]))
    # -2^32: its low 32 bits are those of a zero sum, but cfitsio's cast
    # gives 2^64-2^32 (x86-64), which does not match
    open(outpath("checksum", "datasum_neg_2_32_empty.fits"), "wb").write(header_bytes(ho + ["DATASUM = '-4294967296'"]))
    # a DATASUM / CHECKSUM card with no value: fits_verify_chksum fails
    # (VALUE_UNDEFINED 204) instead of reporting a status
    open(outpath("checksum", "datasum_no_value_empty.fits"), "wb").write(header_bytes(ho + ["DATASUM ="]))
    open(outpath("checksum", "checksum_no_value_empty.fits"), "wb").write(header_bytes(ho + ["CHECKSUM="]))
    # 2^32 on a zero data sum: (unsigned long) 2^32 is not 0, so no match
    # (a 32-bit truncation would make it match).
    open(outpath("checksum", "datasum_2_32_empty.fits"), "wb").write(header_bytes(ho + ["DATASUM = '4294967296'"]))
    # Card order and name case, from cfitsio's own checksummed f64_plain
    # (cards ... DATASUM, CHECKSUM, END). ffvcks reads CHECKSUM first with
    # ffgcrd (the search starts at the top of the header after
    # fits_movabs_hdu) and then DATASUM, whose ffgcrd search starts after
    # the CHECKSUM card it found and wraps to the top; names are compared
    # in upper case.
    fp = open(outpath("checksum", "f64_plain.fits"), "rb").read()
    fcards = [fp[i:i + 80].decode().rstrip() for i in range(0, 2880, 80)]
    fend = fcards.index("END")
    real_ds = next(c for c in fcards[:fend] if c.startswith("DATASUM ="))
    real_ck = next(c for c in fcards[:fend] if c.startswith("CHECKSUM="))
    rest = [c for c in fcards[:fend] if c not in (real_ds, real_ck)]
    fdata = fp[2880:]
    def reorder(name, cards):
        open(outpath("checksum", name + ".fits"), "wb").write(header_bytes(rest + cards) + fdata)
    wrong = "DATASUM = '1'"
    reorder("datasum_dup_after_checksum", [real_ds, real_ck, wrong])     # the second one is read: wrong
    reorder("datasum_dup_before_checksum", [wrong, real_ck, real_ds])    # the second one is read: right
    reorder("datasum_dup_no_checksum", [real_ds, wrong])                 # no CHECKSUM: the first is read
    reorder("datasum_dup_wrap", [real_ds, wrong, real_ck])  # none after CHECKSUM: wraps, the first is read
    reorder("datasum_dup_undefined_after", [real_ds, real_ck, "DATASUM ="])  # the read card has no value: 204
    reorder("checksum_dup_empty_first", ["CHECKSUM= ''", real_ds, real_ck])  # the first CHECKSUM is blank
    reorder("datasum_lowercase_ok", [real_ds.replace("DATASUM", "datasum", 1)])
    reorder("datasum_lowercase_bad", ["datasum = '1'"])


def gen_checksum():
    srcs = {
        "u16_conv.fits": outpath("cfwrite", "u16_conv.fits"),
        "f64_plain.fits": outpath("cfwrite", "f64_plain.fits"),
        "stis_raw.fits": outpath("astropy_data", "o4sp040b0_raw.fits"),
        "testprog.fits": os.path.join(HERE, "..", "testprog.std"),
        "header_only.fits": None,
    }
    hdr_only = outpath("checksum", "_header_only_src.fits")
    open(hdr_only, "wb").write(header_bytes(["SIMPLE  =                    T", "BITPIX  =                    8",
                                             "NAXIS   =                    0"]))
    srcs["header_only.fits"] = hdr_only
    for name, src in srcs.items():
        run([IMGREF, "addchk", src, outpath("checksum", name)])
    os.remove(hdr_only)
    good = open(outpath("checksum", "u16_conv.fits"), "rb").read()
    # Flip one data byte: DATASUM and CHECKSUM both become wrong.
    bad = bytearray(good)
    bad[2880 + 3] ^= 0x01
    open(outpath("checksum", "corrupt_data.fits"), "wb").write(bytes(bad))
    # Change one header byte inside a comment: only CHECKSUM becomes wrong.
    bad = bytearray(good)
    i = bad.index(b"/ number of bits per data pixel")
    bad[i + 2] = ord("N")
    open(outpath("checksum", "corrupt_header.fits"), "wb").write(bytes(bad))


SUBSETS = [
    # file (relative to testdata/images), hdu, type, fpixel, lpixel, inc (1-based, inclusive)
    ("cfwrite/i32_cube.fits", 0, "i64", [2, 1, 1], [4, 3, 2], [1, 2, 1]),
    ("cfwrite/i32_cube.fits", 0, "f64", [1, 2, 2], [3, 3, 2], [2, 1, 1]),
    ("cfwrite/u16_cube.fits", 0, "i64", [1, 1, 1], [5, 4, 3], [2, 3, 2]),
    ("cfwrite/u16_cube.fits", 0, "f64", [2, 2, 1], [4, 4, 3], [1, 1, 1]),
    ("cfwrite/i64_cube.fits", 0, "i64", [1, 1, 1], [3, 3, 2], [1, 1, 1]),
    ("cfwrite/i64_cube.fits", 0, "i64", [2, 1, 2], [3, 3, 2], [1, 2, 1]),
    ("cfwrite/f64_4d.fits", 0, "f64", [1, 1, 1, 1], [4, 3, 2, 2], [3, 2, 1, 1]),
    ("cfwrite/i16_bscale2.fits", 0, "f64", [2], [8], [3]),
    ("astropy_written/u64_cube.fits", 0, "u64", [1, 2, 1], [5, 4, 3], [2, 1, 2]),
    ("astropy_data/arange.fits", 0, "i64", [3, 2, 1], [11, 10, 7], [4, 3, 3]),
    ("astropy_data/scale.fits", 0, "f64", [5, 5], [20, 21], [5, 4]),
    ("astropy_data/o4sp040b0_raw.fits", 1, "i64", [10, 7], [62, 44], [13, 9]),
    # float32 targets: BITPIX 64 rounded once (not through float64), a
    # BSCALE result, and a uint16 cube
    ("cfwrite/i64_f32_round.fits", 0, "f32", [1], [7], [1]),
    ("cfwrite/i64_cube.fits", 0, "f32", [1, 1, 1], [3, 3, 2], [2, 1, 1]),
    ("cfwrite/i16_bscale_1p1.fits", 0, "f32", [1], [9], [1]),
    ("astropy_data/scale.fits", 0, "f32", [5, 5], [20, 21], [5, 4]),
]


def gen_goldens(shas, status):
    gold = outpath("golden", "x")[:-1]
    samples = []
    for sub in ("cfwrite", "astropy_written", "astropy_data", "framing", "checksum"):
        for name in sorted(os.listdir(os.path.join(OUT, sub))):
            if name.endswith(".fits"):
                samples.append(f"{sub}/{name}")
    for s in samples:
        p = os.path.join(OUT, s)
        flat = s.replace("/", "__")
        write_golden(os.path.join(gold, flat + ".json"), run([CFREF, "dump", p]))
        write_golden(os.path.join(gold, flat + ".chksum.json"), run([IMGREF, "chksum", p]))
        write_golden(os.path.join(gold, flat + ".pix.json"), run([IMGREF, "pix", p]))
    # keywords/ and duplicate/ are read by their own tests, not the corpus
    # iteration (their cfitsio scaling differs from fits_read_key's on
    # purpose): the same goldens, kept out of "samples".
    for sub in ("keywords", "duplicate"):
        for name in sorted(os.listdir(os.path.join(OUT, sub))):
            p = os.path.join(OUT, sub, name)
            flat = f"{sub}__{name}"
            write_golden(os.path.join(gold, flat + ".json"), run([CFREF, "dump", p]))
            if sub == "keywords":
                write_golden(os.path.join(gold, flat + ".pix.json"), run([IMGREF, "pix", p]))
    # The top-level testdata samples (cfitsio's own test files and the
    # compression corpus) have cfref goldens in ../cfref/golden; add the
    # reads that dump does not make.
    top = sorted(n for n in os.listdir(os.path.join(HERE, ".."))
                 if n.endswith((".fits", ".fit", ".std")))
    for n in top:
        p = os.path.join(HERE, "..", n)
        write_golden(os.path.join(gold, "testdata__" + n + ".chksum.json"), run([IMGREF, "chksum", p]))
        # pixels only where the corpus tests read them: a file whose image
        # HDUs are all tile-compressed (cfref golden "compressed") has none
        # (its pixels are proven by the compressed-read corpus).
        cf = json.load(open(os.path.join(HERE, "..", "cfref", "golden", n + ".json")))
        if any(h.get("image") and not h["image"].get("compressed") and h["image"].get("npix", 0) > 0
               for h in cf["hdus"]):
            write_golden(os.path.join(gold, "testdata__" + n + ".pix.json"), run([IMGREF, "pix", p]))
    subsets = []
    for f, hdu, ty, fp, lp, inc in SUBSETS:
        out = run([IMGREF, "subset", os.path.join(OUT, f), str(hdu), ty,
                   ",".join(map(str, fp)), ",".join(map(str, lp)), ",".join(map(str, inc))])
        subsets.append(dict(file=f, hdu=hdu, type=ty, fpixel=fp, lpixel=lp, inc=inc,
                            values=json.loads(out)))
    dump_golden(dict(cfitsio=json.loads(run([IMGREF, "chksum", os.path.join(OUT, "cfwrite", "u8_plain.fits")]))["cfitsio"],
                     subsets=subsets), os.path.join(gold, "subsets.json"))
    for name in sorted(os.listdir(os.path.join(OUT, "refused"))):
        write_golden(os.path.join(gold, "refused__" + name + ".open.json"),
                     run([IMGREF, "open", os.path.join(OUT, "refused", name)]))
    vals = [0, 1, 2, 0x7FFFFFFF, 0x80000000, 0xFFFFFFFE, 0xFFFFFFFF, 0x12345678,
            0xDEADBEEF, 0x30303030, 0x0F0F0F0F, 0xA5A5A5A5, 2503531142, 1723431930]
    esum = [dict(value=v, complement=c, ascii=run([IMGREF, "esum", str(v), str(c)]).strip())
            for v in vals for c in (0, 1)]
    dump_golden(esum, os.path.join(gold, "esum.json"))
    # Second reader for the framing cases: what astropy makes of each file
    # (cfitsio cannot read the pixels of a file whose final pad is missing).
    framing = {}
    for name in sorted(os.listdir(os.path.join(OUT, "framing"))):
        try:
            with fits.open(os.path.join(OUT, "framing", name)) as hl:
                hdus = []
                for h in hl:
                    if isinstance(h, fits.GroupsHDU) or h.data is None:
                        hdus.append(None)
                    else:
                        hdus.append([int(v) for v in np.asarray(h.data).ravel()])
                framing[name] = dict(nhdu=len(hl), data=hdus)
        except Exception as e:  # recorded: astropy refuses some of these
            framing[name] = dict(error=f"{type(e).__name__}: {e}")
    dump_golden(framing, os.path.join(gold, "framing_astropy.json"), sort_keys=True)
    dump_golden(dict(astropy=astropy.__version__, numpy=np.__version__, astropy_data_sha256=shas,
                     cfwrite_status=status, samples=samples),
                os.path.join(gold, "manifest.json"), sort_keys=True)


# The scaling-card rule (fits scalingCards) on every keywords/ sample: what
# the library must do, and why. "accept": the cards are read (and cfitsio
# and astropy read the file alike, asserted below); any other code is an
# error the library returns on every read:
#   not_value_card  the card has no "= " in columns 9-10 or a lower-case
#                   keyword (cfitsio and astropy read it differently)
#   not_number      the value is not a FITS integer/real (BLANK: integer)
#   cfitsio_range   a FITS real that cfitsio's strtod rejects with ERANGE
#                   (it substitutes the default; astropy reads the value)
#   hierarch        HIERARCH BSCALE / BZERO (cfitsio applies, astropy not)
#   blank_range     BLANK outside the BITPIX range (cfitsio compares a
#                   BITPIX 32 BLANK as 32 bits, astropy as the integer)
#   blank_undefined BLANK = 1234554321 (cfitsio: no BLANK; astropy: null)
# The generator asserts that every sample cfitsio and astropy read
# differently is an error, so the library never accepts a card the
# references disagree on.
KEYWORD_RULE = {
    "dup_bscale": "accept", "dup_bzero": "accept", "dup_blank": "accept",
    "dup_bscale_bad_first": "accept", "dexp_upper": "accept",
    "bscale_plus": "accept", "bscale_leading_dot": "accept",
    "bscale_trailing_dot": "accept", "bscale_exponent_only": "accept",
    "bscale_zero": "accept", "bscale_min_normal_rounded": "accept",
    "blank_plus": "accept",
    "bscale_noeq_space": "not_value_card", "bscale_key_eq": "not_value_card",
    "bscale_eq_col9_after_eq": "not_value_card", "bscale_no_equals": "not_value_card",
    "bscale_comment_equals": "not_value_card", "bscale_lower": "not_value_card",
    "bscale_string": "not_number", "bzero_string": "not_number", "bscale_inf": "not_number",
    "bscale_neg_infinity": "not_number", "bscale_nan": "not_number",
    "bzero_nan_payload": "not_number", "bscale_overflow": "not_number",
    "bscale_hex": "not_number", "bscale_hex_noexp": "not_number",
    "bscale_subnormal_exact": "not_number", "bscale_trailing_e": "not_number",
    "bscale_logical": "not_number", "bscale_complex": "not_number",
    "bscale_junk": "not_number", "dexp_lower": "not_number", "dexp_hex": "not_number",
    "dexp_hex_two_d": "not_number", "bscale_empty": "not_number",
    "bscale_slash": "not_number", "dup_bscale_bad_last": "not_number",
    "blank_float": "not_number", "blank_empty": "not_number",
    "blank_overflow": "not_number", "blank_string": "not_number",
    "blank_hex": "not_number", "dup_blank_bad_last": "not_number",
    "bscale_underflow": "cfitsio_range", "bscale_subnormal": "cfitsio_range",
    "bscale_below_min_normal": "cfitsio_range",
    "bscale_hierarch": "hierarch",
    "blank_wrap16": "blank_range", "blank_wrap32": "blank_range",
    "blank_magic": "blank_undefined",
}


def astropy_reading(path):
    """astropy 8.0.2's scaled getdata of the primary image: dtype and
    values (NaN -> None), or the exception it raises."""
    import warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            with fits.open(path) as hl:
                d = np.array(hl[0].data)
    except Exception as e:  # astropy refuses the file
        return {"error": type(e).__name__}, None
    # repr strings: exact, and JSON has no Infinity
    vals = [None if (d.dtype.kind == "f" and np.isnan(x)) else repr(float(x)) for x in d.ravel()]
    return {"dtype": str(d.dtype), "values": vals}, d


def gen_keyword_rule():
    names = sorted(n[:-5] for n in os.listdir(os.path.join(OUT, "keywords")))
    assert names == sorted(KEYWORD_RULE), sorted(set(names) ^ set(KEYWORD_RULE))
    out = {}
    for n in names:
        cf = json.load(open(os.path.join(OUT, "golden", f"keywords__{n}.fits.json")))["hdus"][0]["image"]["phys"]
        ap, arr = astropy_reading(os.path.join(OUT, "keywords", n + ".fits"))
        agree = False
        if arr is not None:
            # cfitsio's TDOUBLE values in astropy's dtype (astropy scales
            # BITPIX 16 into float32): the same reading when they match.
            want = [None if v is None else repr(float(np.array(v, dtype=np.float64).astype(arr.dtype))) for v in cf]
            agree = want == ap["values"]
        rule = KEYWORD_RULE[n]
        assert agree or rule != "accept", (n, "accepted, but cfitsio and astropy differ", cf, ap)
        out[n + ".fits"] = {"rule": rule, "refs_agree": agree, "astropy": ap}
    dump_golden(out, os.path.join(OUT, "golden", "keywords_rule.json"), sort_keys=True)


def main():
    status = gen_cfwrite()
    gen_astropy_written()
    shas = gen_astropy_data()
    gen_framing()
    gen_refused()
    gen_checksum()
    gen_datasum()
    gen_keywords()
    gen_keyword_overwrite()
    gen_duplicate()
    gen_write_refused()
    gen_write_overflow()
    gen_strtod()
    gen_goldens(shas, status)
    gen_keyword_rule()


if __name__ == "__main__":
    main()
