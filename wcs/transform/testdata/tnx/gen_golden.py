#!/usr/bin/env python3
"""Generate golden/<name>.json for every headers/<name>.hdr with WCSTools.

Truth source: WCSTools 3.9.7 xy2sky / sky2xy (Jessica Mink, SAO), whose
libwcs tnxpos.c is a C translation of IRAF mwcs/wftnx.x + wfgsurfit.x, the
reference TNX implementation. It reads the WAT strings with hgetm(), which
joins WATn_001, WATn_002, ... verbatim (no blanks trimmed), as the TNX
convention requires. Nothing here is computed by the Go code under test.

Each .hdr file is FITS header card text, one card per line (trailing
blanks may be stripped; every line is padded back to 80 columns). The
script writes it as a header-only FITS file in a temporary directory
(honours TMPDIR), then

  * p2s: xy2sky -d -n 12 <file> x1 y1 x2 y2 ...   on a pixel grid
         (1-based FITS pixels: a 7x7 grid over [1, NAXISn] plus
          off-grid fractional points and points 10% outside each edge);
  * s2p: sky2xy -n 10 <file> ra1 dec1 ...          on the p2s outputs
         (exact printed strings), giving WCSTools' TNX inverse.

The same is done for every refused/<name>.hdr (headers NewTNX must
refuse), writing refused/<name>.json: the WCSTools answer that a wcslib
reading of the header would contradict. A refused header that WCSTools
does not read as celestial at all (xy2sky prints asterisks for every
point) records that line as "wcstools_no_answer" with empty p2s/s2p.

Usage (from this directory, needs xy2sky/sky2xy on PATH or WCSTOOLS_BIN):

    python3 gen_golden.py            # rewrite every golden
    python3 gen_golden.py --check    # exit 1 if any golden would change

Output is deterministic for a given WCSTools build.
"""
import glob
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
WCSTOOLS_VERSION = "3.9.7"


def tool(name):
    d = os.environ.get("WCSTOOLS_BIN")
    p = os.path.join(d, name) if d else shutil.which(name)
    if not p or not os.path.exists(p):
        sys.exit("gen_golden.py: %s not found (set WCSTOOLS_BIN)" % name)
    return p


def wcstools_version(xy2sky):
    out = subprocess.run([xy2sky], capture_output=True, text=True).stderr
    first = (out.splitlines() or [""])[0]
    if ("WCSTools %s," % WCSTOOLS_VERSION) not in first:
        sys.exit("gen_golden.py: need WCSTools %s, got %r" % (WCSTOOLS_VERSION, first))
    return first.strip()


def read_cards(path):
    with open(path, "rb") as f:
        text = f.read().decode("ascii")
    cards = []
    for line in text.split("\n"):
        if line == "":
            continue
        if len(line) > 80:
            sys.exit("%s: card longer than 80 columns: %r" % (path, line))
        cards.append(line.ljust(80))
    if not cards or cards[-1].rstrip() != "END":
        sys.exit("%s: last card must be END" % path)
    return cards


def card_int(cards, key):
    for c in cards:
        if c[:8].rstrip() == key:
            return int(c[10:30])
    sys.exit("missing %s" % key)


def pixel_points(nx, ny):
    pts = []
    for j in range(7):
        for i in range(7):
            pts.append((1 + (nx - 1) * i / 6.0, 1 + (ny - 1) * j / 6.0))
    pts += [(17.25, 33.75), (nx * 0.37 + 0.125, ny * 0.81 - 0.375),
            (nx * 0.5 + 0.5, ny * 0.5 + 0.5),
            (1 - 0.1 * nx, 1 - 0.1 * ny), (1.1 * nx, 1 - 0.1 * ny),
            (1 - 0.1 * nx, 1.1 * ny), (1.1 * nx, 1.1 * ny)]
    return [(round(x, 6), round(y, 6)) for x, y in pts]


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("%s failed (%d): %s" % (cmd[0], r.returncode, r.stderr))
    return r.stdout


def golden_for(path, xy2sky, sky2xy, workdir, allow_no_answer=False):
    name = os.path.splitext(os.path.basename(path))[0]
    rel = os.path.relpath(path, HERE)
    cards = read_cards(path)
    blob = "".join(cards).encode("ascii")
    blob += b" " * ((-len(blob)) % 2880)
    fits = os.path.join(workdir, name + ".fits")
    with open(fits, "wb") as f:
        f.write(blob)
    nx, ny = card_int(cards, "NAXIS1"), card_int(cards, "NAXIS2")
    pix = pixel_points(nx, ny)

    args = [xy2sky, "-d", "-n", "12", fits]
    for x, y in pix:
        args += [repr(x), repr(y)]
    lines = [l for l in run(args).splitlines() if l.strip()]
    if len(lines) != len(pix):
        sys.exit("%s: xy2sky gave %d lines for %d points" % (name, len(lines), len(pix)))
    with open(path, "rb") as f:
        sha = hashlib.sha256(f.read()).hexdigest()
    if allow_no_answer and lines[0].startswith("*"):
        # WCSTools does not recognise the header as celestial (xy2sky
        # prints asterisks): there is no answer to prove one against. Only
        # a refused case may record this; every xy2sky line must agree.
        if not all(l.startswith("*") for l in lines):
            sys.exit("%s: xy2sky answers some points and not others" % name)
        return name, {"header": rel, "header_sha256": sha,
                      "wcstools_no_answer": lines[0].strip(), "p2s": [], "s2p": []}
    p2s, sky_strs = [], []
    for (x, y), l in zip(pix, lines):
        f = l.split()
        sky_strs.append((f[0], f[1]))
        p2s.append([x, y, float(f[0]), float(f[1])])

    args = [sky2xy, "-n", "10", fits]
    for ra, dec in sky_strs:
        args += [ra, dec]
    lines = [l for l in run(args).splitlines() if l.strip()]
    if len(lines) != len(sky_strs):
        sys.exit("%s: sky2xy gave %d lines for %d points" % (name, len(lines), len(sky_strs)))
    s2p = []
    for (ra, dec), l in zip(sky_strs, lines):
        f = l.split("->")[1].split()
        # sky2xy prints "nan" where its own inverse fails (seen only on
        # refused headers); JSON has no NaN, so it is written as null.
        s2p.append([float(ra), float(dec)] + [None if math.isnan(float(v)) else float(v) for v in f[:2]])

    return name, {"header": rel, "header_sha256": sha,
                  "pixel_convention": "1-based FITS pixels",
                  "p2s": p2s, "s2p": s2p}


def render(doc):
    head = {k: v for k, v in doc.items() if k not in ("p2s", "s2p")}
    out = json.dumps(head)[:-1]
    for k in ("p2s", "s2p"):
        if not doc[k]:
            out += ', "%s": []' % k
            continue
        out += ', "%s": [\n' % k + ",\n".join(json.dumps(r, allow_nan=False) for r in doc[k]) + "\n]"
    return out + "}\n"


def main():
    check = "--check" in sys.argv[1:]
    xy2sky, sky2xy = tool("xy2sky"), tool("sky2xy")
    banner = wcstools_version(xy2sky)
    workdir = tempfile.mkdtemp(prefix="tnxgold")
    changed = []
    try:
        jobs = [(p, "golden") for p in sorted(glob.glob(os.path.join(HERE, "headers", "*.hdr")))]
        jobs += [(p, "refused") for p in sorted(glob.glob(os.path.join(HERE, "refused", "*.hdr")))]
        for path, outdir in jobs:
            name, doc = golden_for(path, xy2sky, sky2xy, workdir, allow_no_answer=(outdir == "refused"))
            doc = {"wcstools": banner, **doc}
            text = render(doc)
            out = os.path.join(HERE, outdir, name + ".json")
            name = outdir + "/" + name
            old = open(out).read() if os.path.exists(out) else None
            if old != text:
                changed.append(name)
                if not check:
                    with open(out, "w") as f:
                        f.write(text)
    finally:
        shutil.rmtree(workdir)
    if check and changed:
        print("goldens differ from WCSTools output: %s" % ", ".join(changed))
        sys.exit(1)


if __name__ == "__main__":
    main()
