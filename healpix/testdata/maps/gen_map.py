#!/usr/bin/env python3
"""Write a small HEALPix map with healpy.write_map and a JSON golden that
the Go test (healpix/map_test.go) checks the map against.

    venv/bin/python gen_map.py      # in healpix/testdata/maps

healpy_map_nside16_{nest,ring}.fits: one float64 column whose value at
pixel p is the z component of pixel p's centre (healpy.pix2vec), written
by healpy.write_map with the matching ORDERING. healpy_map_nside16.json:
seeded-random directions, healpy's ang2pix for both orderings and the map
value at that pixel read back with healpy.read_map.
"""
import json
import os
import sys

import healpy as hp
import numpy as np

if hp.__version__ != "1.20.1" or np.__version__ != "2.4.6":
    sys.exit("pinned versions: healpy 1.20.1, numpy 2.4.6 (see ../PROVENANCE.md)")

NSIDE = 16
here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(16)
theta = np.arccos(rng.uniform(-1, 1, 500))
phi = rng.uniform(0, 2 * np.pi, 500)
out = {"nside": NSIDE, "theta": theta.tolist(), "phi": phi.tolist()}
for nest, name in ((True, "nest"), (False, "ring")):
    npix = hp.nside2npix(NSIDE)
    z = np.array(hp.pix2vec(NSIDE, np.arange(npix), nest=nest)[2], dtype=np.float64)
    path = os.path.join(here, f"healpy_map_nside16_{name}.fits")
    hp.write_map(path, z, nest=nest, dtype=np.float64, overwrite=True,
                 column_names=["ZCENTRE"], extra_header=[("ORIGIN", "healpy 1.20.1 gen_map.py")])
    back = hp.read_map(path, nest=nest, dtype=np.float64)
    assert np.array_equal(back, z)
    pix = hp.ang2pix(NSIDE, theta, phi, nest=nest)
    out[name] = {"pix": [int(p) for p in pix], "value": [float(back[p]) for p in pix]}
with open(os.path.join(here, "healpy_map_nside16.json"), "w") as f:
    json.dump(out, f, separators=(",", ":"))
    f.write("\n")
print("wrote healpy_map_nside16_{nest,ring}.fits and healpy_map_nside16.json")
