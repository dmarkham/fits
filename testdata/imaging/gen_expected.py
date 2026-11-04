#!/usr/bin/env python3
"""Expected fits.ReadFloat32 output for every sample in this directory.

    <venv>/bin/python gen_expected.py      # writes readfloat32_expected.json

Pixels are read with astropy (do_not_scale_image_data=True: the stored
integers or floats, independent of cfitsio and of the Go code) and mapped
with the documented ReadFloat32 contract, evaluated here with exact
rational arithmetic (fractions.Fraction) or numpy float32:

  rule "siril8"  BITPIX=8, BSCALE=1, BZERO=0:
                 float32(raw) * float32(.0039215686)      (float32 product)
  rule "siril16" BITPIX=16, BSCALE=1, BZERO=0 or 32768:
                 float32(raw + 32768) * float32(.000015259022)
  rule "exact"   any other integer BITPIX/BSCALE/BZERO:
                 the float32 nearest to u / (rawmax - rawmin), with
                 u = raw - rawmin (BSCALE > 0) or rawmax - raw (BSCALE < 0)
  rule "float"   BITPIX -32/-64: float32(raw * BSCALE + BZERO) evaluated in
                 float64 (raw itself when BSCALE=1, BZERO=0)
  rule "error"   BSCALE = 0, a last BSCALE / BZERO / BLANK card that is
                 not a valid value card (scaling_cards), or a finite float value whose float32
                 rounding overflows (cfitsio reports NUM_OVERFLOW):
                 ReadFloat32 must return an error
  BLANK (integer BITPIX): pixels whose stored value equals BLANK -> NaN.

Tile-compressed samples (*.fz) use cfitsio's own decoded values and null
flags instead (cfnull/golden, see expected_compressed).

Each expected pixel is a float32 bit pattern (uint32), or null for NaN.
The two 65536-pixel all-value files are compared with the Siril harness
output instead (sirilref/golden), so they are listed with "bits": null.
"""
import glob
import hashlib
import json
import math
import os
import re
import sys
import warnings
from fractions import Fraction

import numpy as np
from astropy.io import fits
import astropy

if astropy.__version__ != "8.0.2" or np.__version__ != "2.4.6":
    sys.exit("pinned versions: astropy 8.0.2, numpy 2.4.6 (see PROVENANCE.md)")

HERE = os.path.dirname(os.path.abspath(__file__))
RANGE = {8: (0, 255), 16: (-2**15, 2**15 - 1), 32: (-2**31, 2**31 - 1), 64: (-2**63, 2**63 - 1)}
C8 = np.float32(.0039215686)
C16 = np.float32(.000015259022)
SIRIL_ONLY = {"i16_all.fits", "u16_all.fits"}


DIGEST_MIN = 1024


def digest(vals):
    """An explicit list when it has at most DIGEST_MIN values, else the
    exact digest decoded by imaging_codec_test.go: {"n", "sha256",
    "sample"}, sha256 over the big-endian int64 of every value (null as
    -1; float32 bit patterns and integer pixel values are >= 0), and
    [index, value] samples at the first and last 4 indices, 32 evenly
    spaced ones and the first minimum and maximum."""
    if len(vals) <= DIGEST_MIN:
        return vals
    v = [-1 if x is None else int(x) for x in vals]
    n = len(v)
    h = hashlib.sha256(b"".join(x.to_bytes(8, "big", signed=True) for x in v)).hexdigest()
    idx = set(range(4)) | set(range(n - 4, n)) | {k * (n - 1) // 31 for k in range(32)}
    idx |= {v.index(min(v)), v.index(max(v))}
    return {"n": n, "sha256": h, "sample": [[k, v[k]] for k in sorted(idx)]}


def bits(f):
    return int(np.array(f, dtype=np.float32).view(np.uint32))


def cr32(fr):
    """float32 nearest to the Fraction fr (ties to even), as bits."""
    if fr == 0:
        return 0
    f = np.float32(float(fr))
    best = None
    for c in (np.nextafter(f, np.float32(-np.inf)), f, np.nextafter(f, np.float32(np.inf))):
        key = (abs(Fraction(float(c)) - fr), bits(c) & 1)
        if best is None or key < best[0]:
            best = (key, bits(c))
    return best[1]


def expected_compressed(path):
    """Tile-compressed sample: the pixels, null flags and scaling are what
    cfitsio itself reads (cfnull/golden, fits_read_pixnull TDOUBLE); the
    scaling keywords are the compressed HDU's BSCALE/BZERO cards, which is
    what cfitsio applies. The documented rules are then evaluated on
    cfitsio's physical values p: siril8 float32(p) * C8; siril16
    float32(p - BZERO + 32768) * C16; exact the float32 nearest to
    (p - physical value of the range end that maps to 0) / (|BSCALE| *
    (rawmax - rawmin)); float float32(p)."""
    name = os.path.basename(path)
    g = json.load(open(os.path.join(HERE, "cfnull", "golden", name + ".json")))
    zbitpix = int(next(i[10:].split("/")[0] for i in raw_headers(path)[1] if i[:8] == "ZBITPIX "))
    bad = scaling_cards(raw_headers(path)[1], zbitpix, first_wins=True)[0]
    if bad:
        return {"rule": "error", "bits": None, "invalid_card": bad}
    with fits.open(path, disable_image_compression=True) as hdul:
        h = hdul[1].header
        bitpix = h["ZBITPIX"]
        bscale = float(h.get("BSCALE", 1.0))
        bzero = float(h.get("BZERO", 0.0))
        cols = [str(h.get(f"TTYPE{i}", "")).strip().upper() for i in range(1, h["TFIELDS"] + 1)]
        blank = h.get("ZBLANK", h.get("BLANK")) if bitpix > 0 else None
        shape = [h[f"ZNAXIS{i}"] for i in range(1, h["ZNAXIS"] + 1)]
    rec = {"bitpix": bitpix, "bscale": bscale, "bzero": bzero, "blank": blank,
           "shape": shape, "npix": g["npix"]}
    if bscale == 0:
        # BSCALE=0 is invalid (cfitsio silently returns zeros): error. A
        # per-tile ZBLANK column is applied like cfitsio (the null flags of
        # the golden).
        rec.update(rule="error", bits=None)
        return rec
    assert g["status"] == 0, g
    if bitpix < 0:
        rec["rule"] = "float"
        rec["bits"] = [None if p is None else bits(np.float32(p)) for p in g["phys"]]
        return rec
    if bitpix == 8 and bscale == 1.0 and bzero == 0:
        rule = "siril8"
    elif bitpix == 16 and bscale == 1.0 and bzero in (0, 32768):
        rule = "siril16"
    else:
        rule = "exact"
    rec["rule"] = rule
    lo, hi = RANGE[bitpix]
    phys0 = Fraction(lo) * Fraction(bscale) + Fraction(bzero) if bscale > 0 else \
        Fraction(hi) * Fraction(bscale) + Fraction(bzero)
    out = []
    for p in g["phys"]:
        if p is None:
            out.append(None)
        elif rule == "siril8":
            out.append(bits(np.float32(p) * C8))
        elif rule == "siril16":
            out.append(bits(np.float32(p - bzero + 32768) * C16))
        else:
            u = (Fraction(p) - phys0) / abs(Fraction(bscale))
            out.append(cr32(u / (hi - lo)))
    rec["bits"] = out
    return rec


# A FITS real or integer value (standard 4.2.4 / 4.2.3). Upper-case E and
# D exponents, and lower-case e: cfitsio 4.7.0 and astropy 8.0.2 both read
# 1.0e1 as 10 (card_f32_bzero_lowere_2d.fits). Lower-case d is refused:
# cfitsio silently drops such a card (BZERO 0) while astropy reads 10
# (card_f32_bzero_lowerd_2d.fits).
NUMBER = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([EeD][+-]?\d+)?$")
INTEGER = re.compile(r"^[+-]?\d+$")
SCALING_KEYS = ("BSCALE", "BZERO", "BLANK")


def ffgknm(img):
    """cfitsio's keyword name of a card: columns 1 up to the first blank
    or '=' (at most 74 characters)."""
    t = img.rstrip(" ")[:74]
    m = re.search(r"[ =]", t)
    return t[:m.start()] if m else t


def raw_headers(path):
    """Every HDU header of the file as a list of 80-character card images
    (read byte by byte: astropy cannot open some invalid samples)."""
    with open(path, "rb") as f:
        data = f.read()
    off, out = 0, []
    while off + 2880 <= len(data):
        imgs = []
        while True:
            img = data[off:off + 80].decode("ascii")
            off += 80
            imgs.append(img)
            if img[:8] == "END     ":
                break
        off += -off % 2880
        kv = {i[:8].strip(): i[10:].split("/")[0].strip() for i in imgs if i[8:10] == "= "}
        naxis = int(kv["NAXIS"])
        n = 1 if naxis else 0
        for i in range(1, naxis + 1):
            n *= int(kv[f"NAXIS{i}"])
        n = (n + int(kv.get("PCOUNT", 0))) * int(kv.get("GCOUNT", 1)) * abs(int(kv["BITPIX"])) // 8
        off += n + (-n % 2880)
        out.append(imgs)
    return out


def scaling_cards(imgs, bitpix, first_wins=False):
    """The documented scaling-card rule (fits scalingCards): for each of
    BSCALE, BZERO and BLANK (BLANK on integer BITPIX only) the LAST card
    whose keyword name (ffgknm, any case) is that keyword is the one used,
    as cfitsio 4.7.0 and astropy 8.0.2 both use it on a plain image (the
    FIRST with first_wins: a tile-compressed header, where both use the
    first); it must be a valid value card: the keyword upper case in
    columns 1-8, "= " in columns 9-10 and a value that is a FITS number
    (NUMBER; for BLANK an integer in the BITPIX range, not 1234554321),
    not a quoted string. A HIERARCH card named BSCALE / BZERO / BLANK is
    invalid (cfitsio applies it, astropy ignores it). Returns (the invalid
    keyword or None, {keyword: value})."""
    last = {}
    for img in imgs:
        if img[:8] == "END     ":
            break
        if img.startswith("HIERARCH "):
            name = img[9:].split("=")[0].strip().upper() if "=" in img else ""
            if name in SCALING_KEYS and (name != "BLANK" or bitpix > 0):
                return name, None
            continue
        k = ffgknm(img).upper()
        if k in SCALING_KEYS and not (first_wins and k in last):
            last[k] = img
    vals = {}
    for key in SCALING_KEYS:
        img = last.get(key)
        if img is None or (key == "BLANK" and bitpix < 0):
            continue
        if img[:8] != f"{key:<8}" or img[8:10] != "= ":
            return key, None
        v = img[10:].split("/")[0].strip()
        if key == "BLANK":
            if not INTEGER.match(v) or not RANGE[bitpix][0] <= int(v) <= RANGE[bitpix][1] or int(v) == 1234554321:
                return key, None
            vals[key] = int(v)
            continue
        if not NUMBER.match(v) or not math.isfinite(float(v.replace("D", "E"))):
            return key, None
        vals[key] = float(v.replace("D", "E"))
    return None, vals


def image_header(path):
    """The raw cards of the first HDU with a non-empty data array (the HDU
    the tests read)."""
    for imgs in raw_headers(path):
        kv = {i[:8].strip(): i[10:].split("/")[0].strip() for i in imgs if i[8:10] == "= "}
        naxis = int(kv["NAXIS"])
        if naxis and all(int(kv[f"NAXIS{i}"]) > 0 for i in range(1, naxis + 1)):
            return imgs, int(kv["BITPIX"])
    raise SystemExit(f"{path}: no image HDU")


def bad_scaling(path):
    """The BSCALE / BZERO / BLANK keyword whose last card is not a valid
    value card (scaling_cards), or None."""
    imgs, bitpix = image_header(path)
    return scaling_cards(imgs, bitpix)[0]


def card_refs(path, rec):
    """card_* samples: what cfitsio (cfnull/cards golden, fits_read_pixnull
    TDOUBLE) and astropy (fits.getdata, scaled) read, recorded next to the
    expected entry; for an accepted sample both must equal the documented
    physical values raw * BSCALE + BZERO (BLANK -> null)."""
    name = os.path.basename(path)
    g = json.load(open(os.path.join(HERE, "cfnull", "cards", name + ".json")))
    assert g["status"] == 0, g
    cf = [None if p is None else float(p) for p in g["phys"]]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a = fits.getdata(path).astype(np.float64).ravel()
        ap = [None if np.isnan(x) else float(x) for x in a]
    except Exception as e:  # astropy refuses the file
        ap = f"error: {type(e).__name__}"
    refs = {"cfitsio": cf, "astropy": ap}
    if rec["rule"] != "error":
        imgs, bitpix = image_header(path)
        vals = scaling_cards(imgs, bitpix)[1]
        with fits.open(path, do_not_scale_image_data=True) as hdul:
            raw = np.asarray(hdul[0].data).ravel().astype(np.float64)
        want = []
        for r in raw:
            if "BLANK" in vals and r == vals["BLANK"]:
                want.append(None)
            elif np.isnan(r):
                want.append(None)
            else:
                want.append(float(r * vals.get("BSCALE", 1.0) + vals.get("BZERO", 0.0)))
        norm = lambda xs: [x if x is None else x + 0.0 for x in xs]  # -0 == 0 (cfnull prints %.17g of -0 as -0)
        assert norm(cf) == norm(want), (name, cf, want)
        assert norm(ap) == norm(want), (name, ap, want)
    return refs


def comp_card_refs(path, rec):
    """comp_card_* samples (tile-compressed scaling-card variants): what
    cfitsio (cfnull/golden, fits_read_pixnull TDOUBLE) and astropy
    (getdata of the decompressed image, scaled) read. An accepted sample
    must be read alike by both; a refused one is recorded as they read
    it."""
    name = os.path.basename(path)
    g = json.load(open(os.path.join(HERE, "cfnull", "golden", name + ".json")))
    cf = [None if p is None else float(p) for p in g["phys"]] if g["status"] == 0 else f"status {g['status']}"
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            a = fits.getdata(path, 1).astype(np.float64).ravel()
        ap = [None if np.isnan(x) else float(x) for x in a]
    except Exception as e:  # astropy refuses the file
        ap = f"error: {type(e).__name__}"
    if rec["rule"] != "error":
        assert cf == ap, (name, cf, ap)
    elif not any(i.startswith("HIERARCH ") for i in raw_headers(path)[1]):
        assert cf != ap, (name, "refused, but cfitsio and astropy agree", cf)
    # else: on a tile-compressed header both apply HIERARCH BSCALE, but on
    # a plain image (e.g. the funpack'd file) astropy ignores it, so the
    # one rule refuses it everywhere (keywords/bscale_hierarch).
    return {"cfitsio": cf, "astropy": ap}


def expected(path):
    if path.endswith(".fz"):
        return expected_compressed(path)
    bad = bad_scaling(path)
    if bad:
        return {"rule": "error", "bits": None, "invalid_card": bad}
    hdul = fits.open(path, do_not_scale_image_data=True)
    hdu = next(h for h in hdul if h.is_image and h.data is not None and h.data.size)
    imgs, bitpix = image_header(path)
    vals = scaling_cards(imgs, bitpix)[1]
    bscale = vals.get("BSCALE", 1.0)
    bzero = vals.get("BZERO", 0)
    blank = vals.get("BLANK") if bitpix > 0 else None
    raw = np.asarray(hdu.data).ravel()  # FITS order, NAXIS1 fastest
    name = os.path.basename(path)
    rec = {"bitpix": bitpix, "bscale": bscale, "bzero": float(bzero),
           "blank": blank, "shape": list(reversed(hdu.data.shape)), "npix": int(raw.size)}
    if bscale == 0:
        rec.update(rule="error", bits=None)
        return rec
    if bitpix < 0:
        rec["rule"] = "float"
        phys = raw.astype(np.float64) * bscale + float(bzero)
        with np.errstate(over="ignore"):
            if bscale == 1.0 and bzero == 0:
                vals = raw.astype(np.float32)
            else:
                vals = phys.astype(np.float32)
        if np.any(np.isfinite(phys) & np.isinf(vals)):
            rec.update(rule="error", bits=None)
            return rec
        rec["bits"] = [None if np.isnan(v) else bits(v) for v in vals]
        return rec
    if bitpix == 8 and bscale == 1.0 and bzero == 0:
        rule = "siril8"
    elif bitpix == 16 and bscale == 1.0 and bzero in (0, 32768):
        rule = "siril16"
    else:
        rule = "exact"
    rec["rule"] = rule
    if name in SIRIL_ONLY:
        rec["bits"] = None
        return rec
    lo, hi = RANGE[bitpix]
    out = []
    for r in raw.tolist():
        r = int(r)
        if blank is not None and r == int(blank):
            out.append(None)
        elif rule == "siril8":
            out.append(bits(np.float32(r) * C8))
        elif rule == "siril16":
            out.append(bits(np.float32(r + 32768) * C16))
        else:
            u = r - lo if bscale > 0 else hi - r
            out.append(cr32(Fraction(u, hi - lo)))
    rec["bits"] = out
    return rec


files = sorted(glob.glob(os.path.join(HERE, "*.fits")) + glob.glob(os.path.join(HERE, "*.fz")))
golden = {os.path.basename(p): expected(p) for p in files}
for p in files:
    if os.path.basename(p).startswith("card_"):
        golden[os.path.basename(p)]["refs"] = card_refs(p, golden[os.path.basename(p)])
    if os.path.basename(p).startswith("comp_card_"):
        golden[os.path.basename(p)]["refs"] = comp_card_refs(p, golden[os.path.basename(p)])
for rec in golden.values():
    if rec.get("bits") is not None:
        rec["bits"] = digest(rec["bits"])
with open(os.path.join(HERE, "readfloat32_expected.json"), "w") as f:
    json.dump({"astropy": astropy.__version__, "numpy": np.__version__, "files": golden},
              f, separators=(",", ":"), sort_keys=True)
    f.write("\n")
print(f"{len(golden)} files:", ", ".join(f"{k}={v['rule']}" for k, v in golden.items()))


# --------------------------------------------------------------------------
# stdimage_expected.json: what AsImage / image.Decode must return for each
# uncompressed sample. The accept/kind decision is the documented AsImage
# contract; the pixel values come from astropy.
#   gray    2D BITPIX=8, BSCALE=1, BZERO=0, no pixel equal to BLANK:
#           image.Gray, values = stored bytes
#   gray16  2D BITPIX=16, BSCALE=1, BZERO=32768, no pixel equal to BLANK:
#           image.Gray16, values = physical uint16 (stored + 32768)
#   float   2D BITPIX=-32 unscaled: FloatGray, values = stored float32 bits
#           (lossless)
#   floatphys 2D BITPIX=-64, or -32 with BSCALE/BZERO, when every physical
#           value float32(raw * BSCALE + BZERO) (evaluated in float64) is
#           exactly representable in float32 (NaN allowed): FloatGray with
#           those values; otherwise reject (AsImage is lossless or fails)
#   reject  BSCALE=0, a BSCALE/BZERO card that is not a finite number,
#           and everything else: an error (ErrNotImageCompatible)
# --------------------------------------------------------------------------
def stdimage(path):
    bad = bad_scaling(path)
    if bad:
        return {"kind": "reject", "invalid_card": bad}
    with fits.open(path, do_not_scale_image_data=True) as hdul:
        hdu = hdul[0]
        h = hdu.header
        bitpix, naxis = h["BITPIX"], h["NAXIS"]
        imgs = raw_headers(path)[0]
        vals = scaling_cards(imgs, bitpix)[1]
        bscale, bzero = vals.get("BSCALE", 1.0), vals.get("BZERO", 0.0)
        blank = vals.get("BLANK")
        rec = {"bitpix": bitpix, "naxis": naxis}
        if naxis != 2:
            rec["kind"] = "reject"
            return rec
        raw = np.asarray(hdu.data)
        rec["w"], rec["h"] = int(raw.shape[1]), int(raw.shape[0])
        flat = raw.ravel()
        hit = blank is not None and bitpix > 0 and bool(np.any(flat.astype(np.int64) == int(blank)))
        if bscale == 0:
            rec["kind"] = "reject"  # BSCALE=0 is invalid
        elif bitpix == 8 and bscale == 1 and bzero == 0 and not hit:
            rec["kind"], rec["values"] = "gray", [int(v) for v in flat]
        elif bitpix == 16 and bscale == 1 and bzero == 32768 and not hit:
            rec["kind"], rec["values"] = "gray16", [int(v) + 32768 for v in flat]
        elif bitpix == -32 and bscale == 1 and bzero == 0:
            rec["kind"], rec["values"] = "float", [bits(v) for v in flat.astype(np.float32)]
        elif bitpix < 0:
            phys = flat.astype(np.float64)
            if bscale != 1 or bzero != 0:  # unscaled data is not touched (keeps -0)
                phys = phys * bscale + bzero
            with np.errstate(over="ignore"):
                f32 = phys.astype(np.float32)
            exact = np.isnan(phys) | (f32.astype(np.float64) == phys)
            if not np.all(exact):
                rec["kind"] = "reject"
            else:
                rec["kind"], rec["values"] = "floatphys", [bits(v) for v in f32]
        else:
            rec["kind"] = "reject"
        return rec


std = {os.path.basename(p): stdimage(p) for p in sorted(glob.glob(os.path.join(HERE, "*.fits")))}
for rec in std.values():
    if "values" in rec:
        vals = rec["values"]
        if rec["kind"] in ("float", "floatphys"):  # every NaN as 0x7fc00000 in a digest
            vals = [0x7FC00000 if (v & 0x7F800000) == 0x7F800000 and v & 0x7FFFFF else v for v in vals] \
                if len(vals) > DIGEST_MIN else vals
        rec["values"] = digest(vals)
with open(os.path.join(HERE, "stdimage_expected.json"), "w") as f:
    json.dump({"astropy": astropy.__version__, "files": std}, f, separators=(",", ":"), sort_keys=True)
    f.write("\n")
print("stdimage:", ", ".join(f"{k}={v['kind']}" for k, v in std.items()))
