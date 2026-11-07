#!/usr/bin/env python3
"""Cut real astrometry.net index files down to committed table samples.

The 5200-series index files (index-5205-*, index-5206-*; the LITE builds
KStars installs) carry two kdtree HDUs (AN_FILE 'CKDT' and 'SKDT') whose
headers were written by ESO-QFITS and repeat BITPIX, NAXIS, NAXIS1,
NAXIS2, PCOUNT, GCOUNT and TFIELDS, each time with the same value.
v1.2.2 refused every such header (DuplicateKeywordError); cfitsio opens
them. index-4107 (Tycho-2, 2013) has the same HDU layout without the
repeats and is the control.

For each source this script reads the file (a local copy from --src DIR,
or the pinned URL), checks its sha256, walks the HDUs with a minimal
pure-Python reader (no astropy, so every card is copied byte for byte)
and writes

    anet_<source without "index-">.fits

with every HDU header verbatim, except that a table of more than ROWS
rows is cut to its first ROWS rows: every NAXIS2 card of that header gets
the new row count in the same 20-column value field (comment kept), and
the data unit is the first ROWS rows, zero-padded to a block. The kdtree
headers that carry the repeats have NAXIS2 = 0 and are not touched at
all. PCOUNT must be 0 everywhere (no heap), so a cut row set is complete.

Usage (from testdata/tables):

    python3 extract_anet.py [--src DIR] [--only NAME ...]
    make goldens tfloat-goldens astropy-goldens

Re-running it reproduces the committed samples byte for byte (checked by
hand against the pinned sources; no test downloads them: the committed
sample is pinned by its own sha256 in PROVENANCE.md and by its goldens).
Python >= 3.8, standard library only.
"""

import argparse
import hashlib
import os
import re
import sys
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
BLOCK = 2880
CARD = 80
ROWS = 8

SOURCES = [
    dict(
        name="index-5205-00.fits",
        url="https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5205-00.fits",
        sha256="2b1732c0d80ff0fd871c83f102ea8c0ce7edea69e117c90c0864d50dd0e631bf",
    ),
    dict(
        name="index-5206-38.fits",
        url="https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/LITE/index-5206-38.fits",
        sha256="1b171a956636155e6fa89b607baa1db3e45dd1f3cde73caa706aba0120ce7222",
    ),
    dict(
        name="index-4107.fits",
        url="http://data.astrometry.net/4100/index-4107.fits",
        sha256="0bfe80012338923d1a7db3a0660a00f67a79c9e8ca51563619001c32d6cb94d7",
    ),
]


def out_name(src):
    return "anet_" + src["name"].replace("index-", "", 1)


def load(src, srcdir):
    if srcdir:
        with open(os.path.join(srcdir, src["name"]), "rb") as fh:
            data = fh.read()
    else:
        sys.stderr.write(f"downloading {src['url']}\n")
        with urllib.request.urlopen(src["url"]) as r:
            data = r.read()
    got = hashlib.sha256(data).hexdigest()
    if got != src["sha256"]:
        sys.exit(f"{src['name']}: sha256 {got}, pinned {src['sha256']}")
    return data


def cards_of(header):
    return [header[i:i + CARD] for i in range(0, len(header), CARD)]


def keyval(card):
    """Keyword and the 20-column integer value field of a fixed-format card."""
    key = card[:8].decode("ascii").rstrip()
    if card[8:10] != b"= ":
        return key, None
    try:
        return key, int(card[10:30])
    except ValueError:
        return key, None


def header_end(data, off):
    """Offset just past the header block that holds the END card."""
    pos = off
    while pos < len(data):
        block = data[pos:pos + BLOCK]
        if len(block) < BLOCK:
            sys.exit(f"truncated header at {off}")
        for i in range(0, BLOCK, CARD):
            if block[i:i + CARD] == b"END".ljust(CARD):
                return pos + BLOCK
        pos += BLOCK
    sys.exit(f"no END card in the header at {off}")


def structure(header):
    """BITPIX, NAXISn, PCOUNT, GCOUNT, XTENSION from the FIRST card of each
    (what cfitsio reads by position)."""
    first = {}
    for c in cards_of(header):
        key, val = keyval(c)
        if key and key not in first:
            first[key] = (val, c)
    bitpix = first["BITPIX"][0]
    naxis = first["NAXIS"][0]
    axes = [first[f"NAXIS{i}"][0] for i in range(1, naxis + 1)]
    pcount = first.get("PCOUNT", (0, None))[0]
    gcount = first.get("GCOUNT", (1, None))[0]
    xt = ""
    if "XTENSION" in first:
        # a free-format string ('BINTABLE' / comment): the quoted value
        m = re.match(rb"= *'([^']*)'", first["XTENSION"][1][8:])
        xt = m.group(1).strip().decode("ascii") if m else ""
    return bitpix, axes, pcount, gcount, xt


def padded(n):
    return (n + BLOCK - 1) // BLOCK * BLOCK


def cut_naxis2(header, rows):
    """Every NAXIS2 card (QFITS repeats them) with the value field rewritten
    to rows: columns 11-30 right-aligned, comment untouched."""
    out = bytearray(header)
    n = 0
    for i in range(0, len(header), CARD):
        key, _ = keyval(header[i:i + CARD])
        if key == "NAXIS2":
            out[i + 10:i + 30] = b"%20d" % rows
            n += 1
    if n == 0:
        sys.exit("table without NAXIS2")
    return bytes(out)


def extract(data):
    out = bytearray()
    off = 0
    nhdu = 0
    while off < len(data):
        hend = header_end(data, off)
        header = data[off:hend]
        bitpix, axes, pcount, gcount, xt = structure(header)
        n = 1 if axes else 0  # NAXIS = 0: no data
        for a in axes:
            n *= a
        size = abs(bitpix) // 8 * gcount * (pcount + n)
        dstart, dend = hend, hend + padded(size)
        if pcount != 0:
            sys.exit(f"HDU {nhdu}: PCOUNT {pcount}, a heap is not cut")
        if xt in ("BINTABLE", "TABLE") and len(axes) == 2 and axes[1] > ROWS:
            header = cut_naxis2(header, ROWS)
            rows = data[dstart:dstart + axes[0] * ROWS]
            out += header + rows + b"\0" * (padded(len(rows)) - len(rows))
        else:
            out += header + data[dstart:dend]
        off = dend
        nhdu += 1
    return bytes(out), nhdu


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", help="directory holding the source files (else download)")
    ap.add_argument("--only", nargs="*", help="source names to extract")
    args = ap.parse_args()
    for src in SOURCES:
        if args.only and src["name"] not in args.only:
            continue
        data = load(src, args.src)
        sample, nhdu = extract(data)
        path = os.path.join(HERE, out_name(src))
        with open(path, "wb") as fh:
            fh.write(sample)
        print(f"{out_name(src)}: {nhdu} HDUs, {len(sample)} bytes, sha256 {hashlib.sha256(sample).hexdigest()}")


if __name__ == "__main__":
    main()
