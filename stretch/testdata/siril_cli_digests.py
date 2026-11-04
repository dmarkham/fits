#!/usr/bin/env python3
"""Digest the 18 original siril-cli goldens (and their input) from git history.

The siril-cli pixel files golden/stretch_<preset>.json written by
gen_golden.py were committed in c17580b (2026-04-05) and replaced, in the
round-2 size reduction, by exact sha256 digests written by
gen_siril_goldens.py from the verbatim Siril 1.4.0 oracle. This script keeps
the link to the real siril-cli output: it reads the original files from git
(`git show c17580b:<path>`), computes the same canonical digest the Go tests
use (sha256 over the big-endian bits of every float32, channel-major, every
NaN as 0x7fc00000) and writes golden/siril_cli_sha256.json. The generator
refuses to run unless the oracle reproduces every one of these digests, and
TestSirilCLIDigests requires each committed stretch_<preset>.json to carry
the siril-cli digest.

    python siril_cli_digests.py           # rewrite golden/siril_cli_sha256.json
    python siril_cli_digests.py --check   # compare with the committed file

Needs git (with the repository history) and numpy.
"""

import hashlib
import json
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "golden", "siril_cli_sha256.json")
COMMIT = "c17580b"
REL = "stretch/testdata/golden"
PRESETS = ["kstars1", "kstars2", "kstars3", "kstars4", "kstars5", "kstars6", "kstars7",
           "custom_dark", "custom_extreme", "custom_ultra", "custom_mask", "auto_linked",
           "asinh_light", "asinh_medium", "asinh_strong", "ght_mild", "ght_strong", "linear"]


def canon_sha256(a):
    """sha256 of the canonical float32 encoding (gen_siril_goldens.py canon)."""
    u = np.ascontiguousarray(a, dtype=np.float32).reshape(-1).view(np.uint32).copy()
    u[np.isnan(np.asarray(a, dtype=np.float32).reshape(-1))] = 0x7FC00000
    return hashlib.sha256(u.astype(">u4").tobytes()).hexdigest()


def git(*args):
    return subprocess.run(["git", "-C", HERE] + list(args), check=True, capture_output=True).stdout


def load(name):
    path = f"{COMMIT}:{REL}/{name}"
    blob = git("rev-parse", path).decode().strip()
    g = json.loads(git("show", path))
    a = np.array(g["pixels"], dtype=np.float64)
    f = a.astype(np.float32)
    if not np.array_equal(f.astype(np.float64)[~np.isnan(a)], a[~np.isnan(a)]):
        sys.exit(f"{name}: values are not float32")
    return g, f.reshape(g["shape"]), blob


def main():
    inp, data, blob = load("input.json")
    rec = {"source": f"siril-cli output committed in {COMMIT} (gen_golden.py; siril-cli build not recorded)",
           "encoding": "sha256 over big-endian float32 bits, channel-major (nch, ny, nx), NaN as 0x7fc00000",
           "input": {"blob": blob, "shape": list(data.shape), "sha256": canon_sha256(data)},
           "presets": {}}
    for p in PRESETS:
        g, px, blob = load(f"stretch_{p}.json")
        rec["presets"][p] = {"command": g["command"], "blob": blob, "shape": list(px.shape),
                             "sha256": canon_sha256(px)}
    text = json.dumps(rec, indent=1, sort_keys=True) + "\n"
    if "--check" in sys.argv:
        if open(OUT).read() != text:
            sys.exit(f"{OUT} does not match the digests of the siril-cli files in {COMMIT}")
        print(f"{OUT}: matches {COMMIT}")
        return
    with open(OUT, "w") as f:
        f.write(text)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
