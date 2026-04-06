#!/usr/bin/env python3
"""Generate golden test fixtures for the Go HEALPix indexing package.

Calls healpy functions on a comprehensive set of inputs and writes JSON
files to golden/. The Go test suite loads these and diffs against the
pure-Go implementation — any disagreement is a test failure, zero
tolerance on integer pixel indices, bit-exact on float64 angles.

Usage:
    python gen_golden.py           # writes to golden/ in cwd
    python gen_golden.py <outdir>  # writes to specified directory
"""

import json
import math
import os
import sys
import numpy as np
import healpy as hp

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "golden"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Test point generators
# ---------------------------------------------------------------------------

def standard_directions():
    """50+ (theta, phi) pairs spanning the full sphere."""
    pts = []
    # Poles
    pts.append((0.0, 0.0))           # north pole
    pts.append((math.pi, 0.0))       # south pole
    pts.append((1e-10, 0.0))         # near north pole
    pts.append((math.pi - 1e-10, 0.0))  # near south pole
    # Equator sweep
    for phi in np.linspace(0, 2*math.pi, 17, endpoint=False):
        pts.append((math.pi/2, phi))
    # Mid-latitudes
    for theta in [math.pi/6, math.pi/3, 2*math.pi/3, 5*math.pi/6]:
        for phi in np.linspace(0, 2*math.pi, 9, endpoint=False):
            pts.append((theta, phi))
    # Face boundaries (theta where |cos(theta)| = 2/3)
    theta_boundary = math.acos(2.0/3.0)
    for phi in np.linspace(0, 2*math.pi, 13, endpoint=False):
        pts.append((theta_boundary, phi))
        pts.append((math.pi - theta_boundary, phi))
    return pts

def pixel_sample(nside, count=200):
    """Sample pixel indices spanning the full range."""
    npix = 12 * nside * nside
    if npix <= count:
        return list(range(npix))
    # Include first, last, and evenly-spaced interior
    indices = [0, npix - 1]
    step = max(1, npix // count)
    indices.extend(range(0, npix, step))
    # Add some near ring boundaries
    for ring in [1, nside, nside+1, 2*nside, 3*nside-1, 3*nside, 4*nside-2]:
        if ring >= 1:
            # First pixel in ring (approximate)
            if ring <= nside:
                idx = 2 * ring * (ring - 1)
            elif ring < 3 * nside:
                idx = 2 * nside * (nside - 1) + 4 * nside * (ring - nside)
            else:
                ri = 4 * nside - ring
                idx = 12 * nside * nside - 2 * ri * (ri + 1)
            if 0 <= idx < npix:
                indices.append(idx)
    return sorted(set(indices))

def neighbor_sample(nside, count=30):
    """Sample pixel indices for neighbor testing."""
    npix = 12 * nside * nside
    if npix <= count:
        return list(range(npix))
    indices = [0, npix - 1, npix // 2]
    step = max(1, npix // count)
    indices.extend(range(0, npix, step))
    return sorted(set(indices))[:count]

# ---------------------------------------------------------------------------
# Golden file writers
# ---------------------------------------------------------------------------

NSIDES_SMALL = [1, 2, 4, 8, 16]
NSIDES_LARGE = [64, 256, 1024, 8192]
NSIDES_ALL = NSIDES_SMALL + NSIDES_LARGE

def write_ang2pix():
    """ang2pix golden: for each nside, map standard directions to pixel
    indices in both RING and NESTED ordering."""
    directions = standard_directions()
    for nside in NSIDES_ALL:
        cases = []
        for theta, phi in directions:
            pix_ring = int(hp.ang2pix(nside, theta, phi, nest=False))
            pix_nest = int(hp.ang2pix(nside, theta, phi, nest=True))
            cases.append({
                "theta": theta,
                "phi": phi,
                "pix_ring": pix_ring,
                "pix_nest": pix_nest,
            })
        data = {"nside": nside, "cases": cases}
        path = os.path.join(OUTDIR, f"ang2pix_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    print(f"  ang2pix: {len(NSIDES_ALL)} nside values × {len(directions)} directions")

def write_pix2ang():
    """pix2ang golden: for each nside, map sample pixel indices back to
    (theta, phi) in both orderings."""
    for nside in NSIDES_ALL:
        indices = pixel_sample(nside)
        cases = []
        for pix in indices:
            theta_r, phi_r = hp.pix2ang(nside, pix, nest=False)
            theta_n, phi_n = hp.pix2ang(nside, pix, nest=True)
            cases.append({
                "pix": pix,
                "ring_theta": float(theta_r),
                "ring_phi": float(phi_r),
                "nest_theta": float(theta_n),
                "nest_phi": float(phi_n),
            })
        data = {"nside": nside, "cases": cases}
        path = os.path.join(OUTDIR, f"pix2ang_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    print(f"  pix2ang: {len(NSIDES_ALL)} nside values, {sum(len(pixel_sample(n)) for n in NSIDES_ALL)} total pixels")

def write_ordering():
    """ring2nest / nest2ring golden: exhaustive for small nside, sampled
    for large."""
    for nside in NSIDES_SMALL:
        npix = 12 * nside * nside
        ring2nest = []
        nest2ring = []
        for pix in range(npix):
            ring2nest.append(int(hp.ring2nest(nside, pix)))
            nest2ring.append(int(hp.nest2ring(nside, pix)))
        data = {
            "nside": nside,
            "npix": npix,
            "ring2nest": ring2nest,
            "nest2ring": nest2ring,
        }
        path = os.path.join(OUTDIR, f"ordering_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    # Sampled for large nside
    for nside in NSIDES_LARGE:
        indices = pixel_sample(nside, 500)
        cases = []
        for pix in indices:
            cases.append({
                "ring": pix,
                "nest": int(hp.ring2nest(nside, pix)),
            })
        data = {"nside": nside, "cases": cases}
        path = os.path.join(OUTDIR, f"ordering_sampled_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    print(f"  ordering: {len(NSIDES_SMALL)} exhaustive + {len(NSIDES_LARGE)} sampled")

def write_neighbors():
    """get_all_neighbours golden."""
    for nside in NSIDES_ALL:
        sample = neighbor_sample(nside)
        cases = []
        for pix in sample:
            # RING ordering
            nb_ring = hp.get_all_neighbours(nside, pix, nest=False)
            # NESTED ordering
            nb_nest = hp.get_all_neighbours(nside, pix, nest=True)
            cases.append({
                "pix": pix,
                "neighbors_ring": [int(x) for x in nb_ring],
                "neighbors_nest": [int(x) for x in nb_nest],
            })
        data = {"nside": nside, "cases": cases}
        path = os.path.join(OUTDIR, f"neighbors_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    print(f"  neighbors: {len(NSIDES_ALL)} nside values")

def write_utilities():
    """nside2npix, npix2nside, pixel area, max pixel radius."""
    cases = []
    for nside in NSIDES_ALL:
        npix = int(hp.nside2npix(nside))
        area = float(hp.nside2pixarea(nside))
        resol = float(hp.nside2resol(nside))
        # max_pixrad not in all healpy versions; compute analytically
        max_rad = float(hp.max_pixrad(nside))
        cases.append({
            "nside": nside,
            "npix": npix,
            "pixel_area_sr": area,
            "pixel_resol_rad": resol,
            "max_pixrad": max_rad,
        })
    data = {"cases": cases}
    path = os.path.join(OUTDIR, "utilities.json")
    with open(path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"  utilities: {len(cases)} nside values")

def write_xyf():
    """pix2xyf / xyf2pix golden for validating the internal XY scheme."""
    for nside in NSIDES_SMALL:
        npix = 12 * nside * nside
        cases = []
        for pix in range(npix):
            # RING ordering
            x_r, y_r, f_r = hp.pix2xyf(nside, pix, nest=False)
            # NESTED ordering
            x_n, y_n, f_n = hp.pix2xyf(nside, pix, nest=True)
            cases.append({
                "pix": pix,
                "ring_x": int(x_r), "ring_y": int(y_r), "ring_f": int(f_r),
                "nest_x": int(x_n), "nest_y": int(y_n), "nest_f": int(f_n),
            })
        data = {"nside": nside, "npix": npix, "cases": cases}
        path = os.path.join(OUTDIR, f"xyf_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, indent=2)
    print(f"  xyf: {len(NSIDES_SMALL)} nside values (exhaustive)")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"Generating HEALPix golden fixtures with healpy {hp.__version__}")
    print(f"Output: {os.path.abspath(OUTDIR)}/")
    write_ang2pix()
    write_pix2ang()
    write_ordering()
    write_neighbors()
    write_utilities()
    write_xyf()
    nfiles = len([f for f in os.listdir(OUTDIR) if f.endswith(".json")])
    print(f"Done: {nfiles} golden files written.")
