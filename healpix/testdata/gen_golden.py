#!/usr/bin/env python3
"""Generate golden test fixtures for the Go HEALPix indexing package.

Calls healpy functions on a comprehensive set of inputs and writes JSON
files to golden/. The Go test suite loads these and diffs against the
pure-Go implementation: integer results (pixel indices, neighbours, ring
<-> nest, isnsideok) must agree exactly; float results within the budgets
stated in healpix_test.go / healpix_ext_test.go (angles 1 ULP, unit-vector
components 2^-51 absolute, max_pixrad 2 ULP of the mpmath exact value).
Pinned: healpy 1.20.1, numpy 2.4.6, mpmath 1.4.1 (the script refuses
other versions). Re-running it reproduces every committed file
byte-for-byte; see PROVENANCE.md.

Usage:
    python gen_golden.py           # writes to golden/ in cwd
    python gen_golden.py <outdir>  # writes to specified directory
"""

import base64
import gzip
import hashlib
import json
import math
import os
import struct
import sys
import numpy as np
import healpy as hp
import mpmath

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "golden"
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Compact exact encodings (decoded by golden_codec_test.go)
#   f64(a): "f64:" + base64 of the big-endian IEEE float64 bytes of every
#           value (exact, about 10.7 characters per value).
#   i64(a): an explicit JSON list when it has at most DIGEST_MIN values,
#           else {"n", "sha256", "sample"}: sha256 over the big-endian
#           int64 bytes of every value, and [index, value] pairs at the
#           first and last 4 indices, 32 evenly spaced ones and the first
#           minimum and maximum. A test compares its own full result
#           against the sha256 (any differing value fails) and reports the
#           sampled values that differ.
# ---------------------------------------------------------------------------
DIGEST_MIN = 256


def f64(a):
    a = [float(v) for v in a]
    return "f64:" + base64.b64encode(struct.pack(f">{len(a)}d", *a)).decode("ascii")


def i64(a):
    a = [int(v) for v in a]
    n = len(a)
    if n <= DIGEST_MIN:
        return a
    h = hashlib.sha256(struct.pack(f">{n}q", *a)).hexdigest()
    idx = set(range(4)) | set(range(n - 4, n)) | {k * (n - 1) // 31 for k in range(32)}
    idx |= {a.index(min(a)), a.index(max(a))}
    return {"n": n, "sha256": h, "sample": [[k, a[k]] for k in sorted(idx)]}

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
            json.dump(data, f, separators=(",", ":"))
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
            json.dump(data, f, separators=(",", ":"))
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
            "ring2nest": i64(ring2nest),
            "nest2ring": i64(nest2ring),
        }
        path = os.path.join(OUTDIR, f"ordering_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, separators=(",", ":"))
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
            json.dump(data, f, separators=(",", ":"))
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
            json.dump(data, f, separators=(",", ":"))
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
        json.dump(data, f, separators=(",", ":"))
    print(f"  utilities: {len(cases)} nside values")

def write_xyf():
    """pix2xyf / xyf2pix golden for validating the internal XY scheme."""
    for nside in NSIDES_SMALL:
        npix = 12 * nside * nside
        p = np.arange(npix, dtype=np.int64)  # every pixel, in order
        x_r, y_r, f_r = hp.pix2xyf(nside, p, nest=False)
        x_n, y_n, f_n = hp.pix2xyf(nside, p, nest=True)
        data = {"nside": nside, "npix": npix,
                "ring_x": i64(x_r), "ring_y": i64(y_r), "ring_f": i64(f_r),
                "nest_x": i64(x_n), "nest_y": i64(y_n), "nest_f": i64(f_n)}
        path = os.path.join(OUTDIR, f"xyf_nside{nside}.json")
        with open(path, "w") as f:
            json.dump(data, f, separators=(",", ":"))
    print(f"  xyf: {len(NSIDES_SMALL)} nside values (exhaustive)")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"Generating HEALPix golden fixtures with healpy {hp.__version__} numpy {np.__version__}")
    if hp.__version__ != "1.20.1" or np.__version__ != "2.4.6" or mpmath.__version__ != "1.4.1":
        sys.exit("pinned versions: healpy 1.20.1, numpy 2.4.6, mpmath 1.4.1 (see PROVENANCE.md)")
    print(f"Output: {os.path.abspath(OUTDIR)}/")
    write_ang2pix()
    write_pix2ang()
    write_ordering()
    write_neighbors()
    write_utilities()
    write_xyf()

# ---------------------------------------------------------------------------
# Extended goldens (added by the proof audit). Everything above this line
# writes the original fixture set unchanged; everything below adds coverage
# for claims that previously had no golden: non-power-of-2 RING nside,
# nside up to 2^29, Vec2Pix/Pix2Vec, IsNsidePow2, ordered neighbours at
# face edges/corners, random directions, max_pixrad/nside2resol for every
# nside used.
# ---------------------------------------------------------------------------

NSIDES_NONPOW2 = [3, 5, 6, 7, 12, 100, 1000, 3**10, 3 * 2**27, 2**29 - 1]
NSIDES_HUGE = [2**20, 2**29]
SEED = 20261101


GZIP_OVER = 16384  # goldens larger than this are committed gzipped


def dump(name, data):
    """Write golden/<name>; one over GZIP_OVER bytes is written as
    <name>.gz (gzip level 9, no name, mtime 0), which the Go tests read
    through testutil.ReadGolden and TestHealpixGoldensReproduce compares
    decompressed."""
    path = os.path.join(OUTDIR, name)
    text = (json.dumps(data, separators=(",", ":")) + "\n").encode()
    for p in (path, path + ".gz"):
        if os.path.exists(p):
            os.remove(p)
    if len(text) > GZIP_OVER:
        with open(path + ".gz", "wb") as f:
            f.write(gzip.compress(text, compresslevel=9, mtime=0))
    else:
        with open(path, "wb") as f:
            f.write(text)


def huge_pixel_sample(nside, count=300):
    """Pixel sample for very large nside: ends, ring-zone boundaries,
    face boundaries and seeded-random interior pixels (exact Python ints)."""
    npix = 12 * nside * nside
    rng = np.random.default_rng(SEED + nside % 1000003)
    s = set([0, 1, 2, 3, 4, npix - 1, npix - 2, npix - 5, npix // 2])
    ncap = 2 * nside * (nside - 1)
    for b in [ncap - 1, ncap, ncap + 1, npix - ncap - 1, npix - ncap, npix - ncap + 1,
              ncap + 4 * nside * nside - 1, ncap + 4 * nside * nside]:
        if 0 <= b < npix:
            s.add(b)
    for f in range(12):
        base = f * nside * nside
        for off in (0, 1, nside * nside - 1, nside * nside // 2):
            s.add(base + off)
    # interior: random 63-bit values reduced mod npix (exact)
    for v in rng.integers(0, 2**62, size=count, dtype=np.int64):
        s.add(int(v) % npix)
    return sorted(s)


def random_directions(n, seed):
    rng = np.random.default_rng(seed)
    z = rng.uniform(-1.0, 1.0, n)
    theta = np.arccos(z)
    phi = rng.uniform(0.0, 2 * math.pi, n)
    # near-pole points (theta < 0.01 and > pi - 0.01) where healpix_cxx
    # switches to the sin(theta) branch.
    t2 = rng.uniform(0.0, 0.02, n // 4)
    theta = np.concatenate([theta, t2, math.pi - t2])
    phi = np.concatenate([phi, rng.uniform(0.0, 2 * math.pi, n // 2)])
    return theta, phi


def write_ang2pix_random():
    """ang2pix of seeded-random directions (uniform on the sphere plus
    near-pole points) for every nside, both orderings (RING only for
    non-power-of-2)."""
    theta, phi = random_directions(2000, SEED)
    results = []
    for nside in NSIDES_ALL + NSIDES_NONPOW2 + NSIDES_HUGE:
        rec = {"nside": nside, "pix_ring": i64(hp.ang2pix(nside, theta, phi, nest=False))}
        if hp.isnsideok(nside, nest=True):
            rec["pix_nest"] = i64(hp.ang2pix(nside, theta, phi, nest=True))
        results.append(rec)
    dump("ang2pix_random.json", {"theta": f64(theta), "phi": f64(phi), "results": results})
    print(f"  ang2pix_random: {len(NSIDES_ALL + NSIDES_NONPOW2 + NSIDES_HUGE)} nside values x {len(theta)} directions")


def write_ang2pix_phiwrap():
    """ang2pix with phi outside [0, 2pi): negative, > 2pi, large."""
    phis = [-0.1, -math.pi, -2 * math.pi, 2 * math.pi, 2 * math.pi + 0.1, 100.0, -100.0, 1e6, -1e6, 7 * math.pi / 2]
    thetas = [0.005, 0.3, math.acos(2.0 / 3.0), math.pi / 2, 2.5, math.pi - 0.005]
    th = np.array([t for t in thetas for _ in phis])
    ph = np.array([p for _ in thetas for p in phis])
    results = []
    for nside in [1, 4, 64, 1000, 8192, 2**29]:
        rec = {"nside": nside, "pix_ring": [int(x) for x in hp.ang2pix(nside, th, ph, nest=False)]}
        if hp.isnsideok(nside, nest=True):
            rec["pix_nest"] = [int(x) for x in hp.ang2pix(nside, th, ph, nest=True)]
        results.append(rec)
    dump("ang2pix_phiwrap.json", {"theta": th.tolist(), "phi": ph.tolist(), "results": results})
    print(f"  ang2pix_phiwrap: {len(th)} directions with phi outside [0, 2pi)")


def write_pixfuncs_ext():
    """pix2ang + pix2vec + ring2nest + neighbours for non-power-of-2 RING
    nside and for nside 2^20 / 2^29."""
    for nside in NSIDES_NONPOW2 + NSIDES_HUGE:
        npix = 12 * nside * nside
        if npix <= 1200:
            pix = list(range(npix))
        elif nside in NSIDES_HUGE:
            pix = huge_pixel_sample(nside)
        else:
            pix = huge_pixel_sample(nside, 400)
        p = np.array(pix, dtype=np.int64)
        nest_ok = bool(hp.isnsideok(nside, nest=True))
        rec = {"nside": nside, "pix": pix}

        def scheme(th, ph, vx, vy, vz, x, y, f, nb):
            return {"theta": f64(th), "phi": f64(ph), "vx": f64(vx), "vy": f64(vy), "vz": f64(vz),
                    "x": i64(x), "y": i64(y), "f": i64(f),
                    "neighbours": i64([int(nb[k][i]) for i in range(len(pix)) for k in range(8)])}
        th, ph = hp.pix2ang(nside, p, nest=False)
        vx, vy, vz = hp.pix2vec(nside, p, nest=False)
        x, y, f = hp.pix2xyf(nside, p, nest=False)
        nb = hp.get_all_neighbours(nside, p, nest=False)  # shape (8, n)
        rec["ring"] = scheme(th, ph, vx, vy, vz, x, y, f, nb)
        if nest_ok:
            th, ph = hp.pix2ang(nside, p, nest=True)
            vx, vy, vz = hp.pix2vec(nside, p, nest=True)
            x, y, f = hp.pix2xyf(nside, p, nest=True)
            nb = hp.get_all_neighbours(nside, p, nest=True)
            rec["nest"] = scheme(th, ph, vx, vy, vz, x, y, f, nb)
            rec["ring2nest"] = i64(hp.ring2nest(nside, p))
            rec["nest2ring"] = i64(hp.nest2ring(nside, p))
        dump(f"pixfuncs_nside{nside}.json", rec)
    print(f"  pixfuncs: {len(NSIDES_NONPOW2 + NSIDES_HUGE)} nside values (non-pow2 RING + 2^20, 2^29)")


def write_pix2vec():
    """pix2vec for the original nside set and the original pixel samples."""
    for nside in NSIDES_ALL:
        p = np.array(pixel_sample(nside), dtype=np.int64)
        rv = hp.pix2vec(nside, p, nest=False)
        nv = hp.pix2vec(nside, p, nest=True)
        dump(f"pix2vec_nside{nside}.json", {
            "nside": nside, "pix": [int(a) for a in p],
            "ring": [f64(a) for a in rv], "nest": [f64(a) for a in nv]})
    print(f"  pix2vec: {len(NSIDES_ALL)} nside values")


def write_vec2pix():
    """vec2pix of non-unit vectors: random directions scaled by assorted
    lengths, plus axis-aligned and exact-boundary vectors."""
    rng = np.random.default_rng(SEED + 7)
    theta, phi = random_directions(1000, SEED + 3)
    scale = rng.choice([1e-3, 0.5, 1.0, 3.0, 1e6], size=len(theta))
    vx = scale * np.sin(theta) * np.cos(phi)
    vy = scale * np.sin(theta) * np.sin(phi)
    vz = scale * np.cos(theta)
    extra = [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
             (1, 1, 0), (1, -1, 0), (2, 0, 1), (5, 0, 3), (3, 4, 12), (-3, -4, -12)]
    vx = np.concatenate([vx, [e[0] for e in extra]]).astype(float)
    vy = np.concatenate([vy, [e[1] for e in extra]]).astype(float)
    vz = np.concatenate([vz, [e[2] for e in extra]]).astype(float)
    results = []
    for nside in NSIDES_ALL + NSIDES_NONPOW2 + NSIDES_HUGE:
        rec = {"nside": nside, "pix_ring": i64(hp.vec2pix(nside, vx, vy, vz, nest=False))}
        if hp.isnsideok(nside, nest=True):
            rec["pix_nest"] = i64(hp.vec2pix(nside, vx, vy, vz, nest=True))
        results.append(rec)
    dump("vec2pix.json", {"vx": f64(vx), "vy": f64(vy), "vz": f64(vz), "results": results})
    print(f"  vec2pix: {len(NSIDES_ALL + NSIDES_NONPOW2 + NSIDES_HUGE)} nside values x {len(vx)} vectors")


def write_neighbours_edges():
    """Ordered get_all_neighbours for every pixel on a face edge or corner
    (where the face-crossing logic runs), both orderings."""
    for nside in [1, 2, 4, 8, 16, 64, 8192]:
        if nside <= 16:
            pix = list(range(12 * nside * nside))
        else:
            pix = set()
            edge = sorted(set([0, 1, 2, nside // 2, nside - 2, nside - 1]))
            for f in range(12):
                for a in edge:
                    for b in (0, nside - 1):
                        pix.add(int(hp.xyf2pix(nside, a, b, f, nest=True)))
                        pix.add(int(hp.xyf2pix(nside, b, a, f, nest=True)))
            pix = sorted(pix)
        p = np.array(pix, dtype=np.int64)
        nr = hp.get_all_neighbours(nside, p, nest=False)
        nn = hp.get_all_neighbours(nside, p, nest=True)
        dump(f"neighbours_ordered_nside{nside}.json", {
            "nside": nside, "pix": pix,
            "ring": i64([int(nr[k][i]) for i in range(len(pix)) for k in range(8)]),
            "nest": i64([int(nn[k][i]) for i in range(len(pix)) for k in range(8)])})
    print("  neighbours_ordered: 7 nside values (edges/corners, ordered SW,W,NW,N,NE,E,SE,S)")


def max_pixrad_exact(nside):
    """The healpix_cxx max_pixrad definition (angle between the point
    z=2/3, phi=pi/(4 nside) and the point z=1-(1-1/nside)^2/3, phi=0)
    evaluated in 60-digit arithmetic with mpmath and rounded once to
    double: an exact reference for the quantity, independent of libm."""
    with mpmath.workdps(60):
        n = mpmath.mpf(nside)
        za = mpmath.mpf(2) / 3
        pa = mpmath.pi / (4 * n)
        zb = 1 - (1 - 1 / n) ** 2 / 3
        sa = mpmath.sqrt(1 - za * za)
        sb = mpmath.sqrt(1 - zb * zb)
        ax, ay = sa * mpmath.cos(pa), sa * mpmath.sin(pa)
        bx = sb
        cx, cy, cz = ay * zb, za * bx - ax * zb, -ay * bx
        cross = mpmath.sqrt(cx * cx + cy * cy + cz * cz)
        dot = ax * bx + za * zb
        return float(mpmath.atan2(cross, dot))


def write_utilities_ext():
    """nside2npix, nside2pixarea, nside2resol, max_pixrad and isnsideok for
    every nside used anywhere, plus invalid values."""
    cases = []
    for nside in sorted(set(NSIDES_ALL + NSIDES_NONPOW2 + NSIDES_HUGE + [2**k for k in range(30)])):
        cases.append({
            "nside": nside,
            "npix": int(hp.nside2npix(nside)),
            "pixel_area_sr": float(hp.nside2pixarea(nside)),
            "pixel_resol_rad": float(hp.nside2resol(nside)),
            "max_pixrad": float(hp.max_pixrad(nside)),
            "max_pixrad_exact": max_pixrad_exact(nside),
        })
    ok = []
    for n in list(range(-3, 70)) + [2**k for k in range(31)] + [2**k + 1 for k in range(2, 31)] + [3 * 2**20, 2**29 - 1, 2**30, 2**31]:
        ok.append({"nside": n, "nest_ok": bool(hp.isnsideok(n, nest=True)),
                   "ring_ok": bool(hp.isnsideok(n, nest=False))})
    npix_cases = []
    for npix in [0, 1, 11, 12, 13, 24, 36, 47, 48, 108, 192, 12 * 9, 12 * 2**58, 12 * 2**58 + 12, 12 * (2**29 - 1)**2, 12 * (2**29)**2, 12 * (2**29 + 1)**2]:
        npix_cases.append({"npix": npix, "valid": bool(hp.isnpixok(npix)),
                           "nside": int(hp.npix2nside(npix)) if hp.isnpixok(npix) else 0})
    dump("utilities_ext.json", {"cases": cases, "isnsideok": ok, "npix2nside": npix_cases})
    print(f"  utilities_ext: {len(cases)} nside values, {len(ok)} isnsideok cases")


POLAR_NSIDES = [8192, 2**20, 2**29 - 1, 2**29]
HPREF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hpref", "hpref")


def cr(f, *args):
    """f(args) in 256-bit mpmath, rounded once to the nearest double.

    float(mpf) rounds to 53 bits first, then to the subnormal grid (double
    rounding), and 256 bits cannot resolve f(t) - t for subnormal t. So a
    result at or below 2^-1022 in magnitude is recomputed with enough bits
    to hold its exponent and rounded once on the 2^-1074 grid (no exact
    tie is possible for these transcendental values)."""
    with mpmath.workprec(256):
        v = f(*[mpmath.mpf(a) for a in args])
    if v == 0 or abs(v) > mpmath.ldexp(1, -1022):
        return float(v)
    with mpmath.workprec(256 + 3 * 1100):
        v = f(*[mpmath.mpf(a) for a in args])
        q = int(mpmath.nint(mpmath.ldexp(v, 1074)))
    return math.copysign(math.ldexp(q, -1074), float(v)) if q else math.copysign(0.0, float(v))


def hpref_run(lines):
    """Run the healpix_cxx harness (hpref/hpref.cc) on query lines."""
    import subprocess
    if not os.path.exists(HPREF):
        sys.exit(f"{HPREF} missing: build it with `make -C hpref` (see PROVENANCE.md)")
    out = subprocess.run([HPREF], input="".join(lines), capture_output=True, text=True, check=True).stdout
    vals = [int(v) for v in out.split()]
    assert len(vals) == len(lines)
    return vals


def polar_directions():
    """Directions on, and up to 2 ULP either side of, polar-cap ring and
    pixel boundaries. Each boundary is solved in 50-digit mpmath for the
    double tp = frac(phi * 2/pi) that healpix itself computes: the
    healpix_cxx cap condition tp * tmp = k (or (1 - tp) * tmp = k) with
    tmp = nside * sqrt(6) * sin(theta / 2). Rings k run log-uniformly
    from 1 to the ring at theta = 0.1415 (beyond acos(0.99), the vec2pix
    switch-over), for both poles and every POLAR_NSIDES value."""
    rng = np.random.default_rng(SEED + 99)
    inv_halfpi = 2.0 / math.pi
    th_all, ph_all = [], []
    with mpmath.workdps(50):
        s6 = mpmath.sqrt(6)
        for nside in POLAR_NSIDES:
            kmax = int(nside * math.sqrt(6) * math.sin(0.1415 / 2))
            for south in (False, True):
                # 2000 candidate boundaries; keep the 100 whose nearest
                # double lies closest to the boundary relative to the
                # distance from the pole (near pi the theta grid is coarse
                # relative to pi - theta, so most candidates sit far from
                # their boundary in sin(theta) terms).
                cands = []
                for _ in range(2000):
                    k = int(round(math.exp(rng.uniform(0.0, math.log(kmax)))))
                    k = max(1, min(k, kmax))
                    q = int(rng.integers(0, 4))
                    phi = (q + rng.uniform(0.02, 0.98)) * (math.pi / 2)
                    tt = math.fmod(phi * inv_halfpi, 4.0)
                    tp = tt - int(tt)
                    frac = mpmath.mpf(tp) if rng.uniform() < 0.5 else 1 - mpmath.mpf(tp)
                    arg = (mpmath.mpf(k) / frac) / (nside * s6)
                    if arg >= 1:
                        continue
                    t_exact = 2 * mpmath.asin(arg)
                    t_target = mpmath.pi - t_exact if south else t_exact
                    t0 = float(t_target)
                    rel = float(abs(mpmath.mpf(t0) - t_target) / t_exact)
                    cands.append((rel, len(cands), t0, phi))
                cands.sort()
                for _, _, t0, phi in cands[:100]:
                    for d in (-2, -1, 0, 1, 2):
                        t = t0
                        for _ in range(abs(d)):
                            t = float(np.nextafter(t, np.inf if d > 0 else -np.inf))
                        th_all.append(t)
                        ph_all.append(phi)
    return th_all, ph_all


def write_ang2pix_polar():
    """ang2pix / vec2pix at polar-cap boundaries (see polar_directions).

    Within an ULP of a boundary the pixel depends on how cos/sin/atan2
    round. healpy calls glibc, which is not correctly rounded for ~0.1% of
    inputs, and the Go package calls Go's math.Cos/Sin/Atan2, which differ
    from both in the last bit more often. So the expected value is the real
    healpix_cxx loc2pix (hpref harness) fed correctly rounded cos/sin/atan2
    (mpmath), healpy's own answer is recorded next to it, and the test
    accepts a different pixel only where the Go answer itself changes
    under a few-ULP nudge (boundary-ambiguous). The script checks that the harness reproduces healpy
    exactly when it calls glibc itself, and that every place where healpy
    differs from the expected value is one where glibc's cos, sin or atan2
    is not correctly rounded."""
    th_all, ph_all = polar_directions()
    n = len(th_all)
    with mpmath.workprec(256):
        vx = [cr(lambda a, b: mpmath.sin(a) * mpmath.cos(b), a, b) for a, b in zip(th_all, ph_all)]
        vy = [cr(lambda a, b: mpmath.sin(a) * mpmath.sin(b), a, b) for a, b in zip(th_all, ph_all)]
        vz = [cr(mpmath.cos, a) for a in th_all]
    cos_cr = [cr(mpmath.cos, a) for a in th_all]
    sin_cr = [cr(mpmath.sin, a) for a in th_all]
    atan_cr = [0.0 if (x == 0 and y == 0) else cr(mpmath.atan2, y, x) for x, y in zip(vx, vy)]
    glibc_bad_ang = [math.cos(a) != c or math.sin(a) != s for a, c, s in zip(th_all, cos_cr, sin_cr)]
    glibc_bad_vec = [(0.0 if (x == 0 and y == 0) else math.atan2(y, x)) != f for x, y, f in zip(vx, vy, atan_cr)]
    th, ph = np.array(th_all), np.array(ph_all)
    ax, ay, az = np.array(vx), np.array(vy), np.array(vz)
    results = []
    ndiff = 0
    for nside in POLAR_NSIDES:
        rec = {"nside": nside}
        for nest in ((0, 1) if hp.isnsideok(nside, nest=True) else (0,)):
            key = "nest" if nest else "ring"
            hp_ang = [int(x) for x in hp.ang2pix(nside, th, ph, nest=bool(nest))]
            hp_vec = [int(x) for x in hp.vec2pix(nside, ax, ay, az, nest=bool(nest))]
            q = []
            for i in range(n):
                t, p = th_all[i].hex(), ph_all[i].hex()
                x, y, z = vx[i].hex(), vy[i].hex(), vz[i].hex()
                q.append(f"a {nside} {nest} {t} {p}\n")
                q.append(f"v {nside} {nest} {x} {y} {z}\n")
                q.append(f"A {nside} {nest} {t} {p} {math.cos(th_all[i]).hex()} {math.sin(th_all[i]).hex()}\n")
                q.append(f"A {nside} {nest} {t} {p} {cos_cr[i].hex()} {sin_cr[i].hex()}\n")
                q.append(f"V {nside} {nest} {x} {y} {z} {atan_cr[i].hex()}\n")
            r = hpref_run(q)
            h_ang, h_vec, h_ang_glibc, exp_ang, exp_vec = r[0::5], r[1::5], r[2::5], r[3::5], r[4::5]
            assert h_ang == hp_ang, f"hpref ang2pix != healpy at nside {nside} {key}"
            assert h_vec == hp_vec, f"hpref vec2pix != healpy at nside {nside} {key}"
            assert h_ang_glibc == hp_ang, f"hpref ang2pix body != healpy at nside {nside} {key}"
            diffs = []
            for i in range(n):
                if exp_ang[i] != hp_ang[i]:
                    assert glibc_bad_ang[i], f"ang {i}: healpy differs although glibc is correctly rounded"
                    diffs.append(["ang", i, hp_ang[i]])
                if exp_vec[i] != hp_vec[i]:
                    assert glibc_bad_vec[i], f"vec {i}: healpy differs although glibc is correctly rounded"
                    diffs.append(["vec", i, hp_vec[i]])
            ndiff += len(diffs)
            # Explicit lists, not digests: the test needs every expected
            # value to tell a boundary-ambiguous mismatch from a real one.
            rec["ang_" + key] = exp_ang
            rec["vec_" + key] = exp_vec
            rec["healpy_differs_" + key] = diffs
        results.append(rec)
    dump("ang2pix_polar.json", {"theta": f64(th_all), "phi": f64(ph_all), "vx": f64(vx), "vy": f64(vy),
                                "vz": f64(vz), "cos_cr": f64(cos_cr), "sin_cr": f64(sin_cr),
                                "atan2_cr": f64(atan_cr), "results": results})
    print(f"  ang2pix_polar: {n} polar boundary directions x {len(POLAR_NSIDES)} Nside; "
          f"glibc cos/sin not correctly rounded for {sum(glibc_bad_ang)}, atan2 for {sum(glibc_bad_vec)}; "
          f"healpy differs from correctly rounded healpix_cxx in {ndiff} answers")


def write_vec2pix_extreme():
    """Vectors whose squared length overflows or underflows double:
    unit-scale vectors multiplied by 2^k (exact, so the direction is
    exactly the same). The expected pixel is healpy's answer at unit
    scale. healpy's own answer for the scaled vector is recorded too; it
    is wrong whenever the squared length underflows or overflows."""
    import warnings
    base = [(1.0, 1.0, 0.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (-3.0, 4.0, 12.0),
            (0.25, -0.5, -0.75), (0.1, 0.2, 0.97)]
    # 2^500..2^600 and 2^-600..2^-500 pin Vec2Pix's rescaling bounds
    # (0x1p500 / 0x1p-500): squares overflow from about 2^512 and
    # underflow (lose bits, or flush to 0) below about 2^-511.
    exps = [0, 600, 1000, 1019, -600, -1000, -1070,
            501, 511, 513, 520, 550, 589, 590, 599,
            -501, -511, -513, -520, -550, -589, -590, -599]
    cases = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for nside in [4, 8192, 2**29]:
            for b in base:
                ur = int(hp.vec2pix(nside, *b))
                un = int(hp.vec2pix(nside, *b, nest=True))
                for k in exps:
                    v = [math.ldexp(c, k) for c in b]
                    if any(math.ldexp(x, -k) != c for x, c in zip(v, b)):
                        continue  # not exact (subnormal rounding): skip
                    cases.append({"nside": nside, "v": v, "ring": ur, "nest": un,
                                  "healpy_ring": int(hp.vec2pix(nside, *v)),
                                  "healpy_nest": int(hp.vec2pix(nside, *v, nest=True))})
    nwrong = sum(c["healpy_ring"] != c["ring"] or c["healpy_nest"] != c["nest"] for c in cases)
    dump("vec2pix_extreme.json", {"cases": cases})
    print(f"  vec2pix_extreme: {len(cases)} cases (healpy wrong on {nwrong})")


# Ang2Pix switches from the cos(theta) to the sin(theta) formula at
# theta < 0.01 and theta > 3.14159 - 0.01 (the C double expression, which
# equals the Go constant), Vec2Pix at |z| > 0.99. Each band straddles one
# switch-over point and the nearby values a wrong threshold would use
# (pi - 0.01, 0.0099, 3.14159 - 0.011, 0.9901).
SWITCH_ANG_BANDS = [("north_0.01", 0.0098, 0.0102), ("south_3.13159", 3.1304, 3.1318)]
SWITCH_ANG_POINTS = [0.01, 3.14159 - 0.01, math.pi - 0.01]
SWITCH_VEC_BANDS = [("z_0.99", math.acos(0.9902), math.acos(0.9898))]
SWITCH_KEEP = 40


def _cap_tmp(nside, theta):
    """healpix_cxx polar-cap tmp = nside*sqrt(3*(1-|z|)), exact (mpmath)."""
    th = mpmath.mpf(theta)
    if th > mpmath.pi / 2:
        th = mpmath.pi - th
    return nside * mpmath.sqrt(6) * mpmath.sin(th / 2)


def _healpix_tp(phi):
    tt = math.fmod(phi * (2.0 / math.pi), 4.0)
    return tt - int(tt)


def _ulps(x, d):
    for _ in range(abs(d)):
        x = float(np.nextafter(x, np.inf if d > 0 else -np.inf))
    return x


def _switch_ang_candidates(rng, nside, lo, hi, n):
    """Directions with theta in [lo, hi] (either pole) on, and 2 ULP either
    side of, a polar-cap pixel boundary tp*tmp = k or (1-tp)*tmp = k."""
    out = []
    south = lo > 1
    with mpmath.workdps(50):
        s6 = mpmath.sqrt(6)
        for _ in range(n):
            phi = (int(rng.integers(0, 4)) + rng.uniform(0.02, 0.98)) * (math.pi / 2)
            tp = _healpix_tp(phi)
            frac = mpmath.mpf(tp) if rng.uniform() < 0.5 else 1 - mpmath.mpf(tp)
            k = int(mpmath.floor(frac * _cap_tmp(nside, rng.uniform(lo, hi))))
            if k < 1:
                continue
            arg = (mpmath.mpf(k) / frac) / (nside * s6)
            if arg >= 1:
                continue
            t = 2 * mpmath.asin(arg)
            t = mpmath.pi - t if south else t
            if not (lo <= t <= hi):
                continue
            t0 = float(t)
            out += [(_ulps(t0, d), phi) for d in (-2, -1, 0, 1, 2)]
    return out


def _switch_point_candidates(rng, nside, theta, n, dts=range(-3, 4), dps=range(-2, 3)):
    """theta and (by default) 3 ULP either side of it, with phi on, and 2
    ULP either side of, a pixel boundary at that theta."""
    out = []
    with mpmath.workdps(50):
        tmp = _cap_tmp(nside, theta)
        for _ in range(n):
            k = int(rng.integers(1, int(tmp) + 1))
            f = mpmath.mpf(k) / tmp
            if rng.uniform() < 0.5:
                f = 1 - f
            phi0 = float((int(rng.integers(0, 4)) + f) * (mpmath.pi / 2))
            for dt in dts:
                for dp in dps:
                    out.append((_ulps(theta, dt), _ulps(phi0, dp)))
    return out


def write_branch_switch():
    """Ang2Pix / Vec2Pix at the cos -> sin switch-over points.

    Inputs: directions within each SWITCH_*_BANDS band on (and 2 ULP
    around) a polar-cap pixel boundary, and the switch-over thetas
    themselves (+-3 ULP) with phi on a boundary; vectors with z == +-0.99
    exactly. Only inputs where the cos and the sin formula give different
    pixels are kept (hpref C/S/Q/P), SWITCH_KEEP per band and Nside, so a
    test with a threshold moved inside a band gets a wrong pixel.
    Expected: the real healpix_cxx (hpref A/V) fed correctly rounded
    cos/sin/atan2. Inputs where that differs from healpy's own answer
    (glibc not correctly rounded) are dropped, so every expected value is
    also healpy's answer; the generator asserts it."""
    rng = np.random.default_rng(SEED + 123)
    ang, vec = [], []
    dropped = 0
    for nside in POLAR_NSIDES:
        nests = (0, 1) if hp.isnsideok(nside, nest=True) else (0,)
        groups = []
        for name, lo, hi in SWITCH_ANG_BANDS:
            groups.append(("ang", name, _switch_ang_candidates(rng, nside, lo, hi, 600)))
        for t in SWITCH_ANG_POINTS:
            groups.append(("ang", f"theta_{t!r}", _switch_point_candidates(rng, nside, t, 40)))
        for name, lo, hi in SWITCH_VEC_BANDS:
            c = _switch_ang_candidates(rng, nside, lo, hi, 600)
            c += _switch_ang_candidates(rng, nside, math.pi - hi, math.pi - lo, 600)
            groups.append(("vec", name, c))
        exact = []
        for z in (0.99, -0.99):
            for th, ph in _switch_point_candidates(rng, nside, math.acos(z), 1200, (0,), range(-8, 9)):
                exact.append((z, ph))
        groups.append(("vecz", "z_exact_0.99", exact))
        for kind, name, cands in groups:
            seen, items = set(), []
            for a, b in cands:
                if kind == "ang":
                    if not (0 <= a <= math.pi) or (a, b) in seen:
                        continue
                    seen.add((a, b))
                    items.append({"theta": a, "phi": b})
                    continue
                if kind == "vec":
                    v = [cr(lambda t, p: mpmath.sin(t) * mpmath.cos(p), a, b),
                         cr(lambda t, p: mpmath.sin(t) * mpmath.sin(p), a, b), cr(mpmath.cos, a)]
                else:  # z exactly +-0.99, x, y correctly rounded
                    with mpmath.workprec(256):
                        st = mpmath.sqrt(1 - mpmath.mpf(a) ** 2)
                        v = [float(st * mpmath.cos(b)), float(st * mpmath.sin(b)), a]
                    xl = 1.0 / math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
                    if abs(v[2] * xl) != 0.99:
                        continue  # healpix's nz is not exactly 0.99
                if tuple(v) in seen:
                    continue
                seen.add(tuple(v))
                items.append({"v": v})
            q = []
            for it in items:
                if "theta" in it:
                    t, p = it["theta"], it["phi"]
                    c, s_ = cr(mpmath.cos, t).hex(), cr(mpmath.sin, t).hex()
                    q.append(f"C {nside} 0 {t.hex()} {p.hex()} {c} {s_}\n")
                    q.append(f"S {nside} 0 {t.hex()} {p.hex()} {c} {s_}\n")
                else:
                    x, y, z = it["v"]
                    f = cr(mpmath.atan2, y, x).hex()
                    q.append(f"Q {nside} 0 {x.hex()} {y.hex()} {z.hex()} {f}\n")
                    q.append(f"P {nside} 0 {x.hex()} {y.hex()} {z.hex()} {f}\n")
            r = hpref_run(q) if q else []
            items = [it for it, a, b in zip(items, r[0::2], r[1::2]) if a != b]
            kept = []
            for it in items:
                rec = {"nside": nside, "band": name}
                rec.update(it)
                # The correctly rounded libm values healpix_cxx was fed: the
                # Go test runs its post-libm code on them and requires the
                # expected pixel exactly.
                if "theta" in it:
                    rec["cos_cr"] = cr(mpmath.cos, it["theta"])
                    rec["sin_cr"] = cr(mpmath.sin, it["theta"])
                else:
                    rec["atan2_cr"] = cr(mpmath.atan2, it["v"][1], it["v"][0])
                ok = True
                for nest in nests:
                    key = "nest" if nest else "ring"
                    if "theta" in it:
                        t, p = it["theta"], it["phi"]
                        line = f"A {nside} {nest} {t.hex()} {p.hex()} {cr(mpmath.cos, t).hex()} {cr(mpmath.sin, t).hex()}\n"
                        h = int(hp.ang2pix(nside, t, p, nest=bool(nest)))
                        glibc = f"a {nside} {nest} {t.hex()} {p.hex()}\n"
                    else:
                        x, y, z = it["v"]
                        line = f"V {nside} {nest} {x.hex()} {y.hex()} {z.hex()} {cr(mpmath.atan2, y, x).hex()}\n"
                        h = int(hp.vec2pix(nside, x, y, z, nest=bool(nest)))
                        glibc = f"v {nside} {nest} {x.hex()} {y.hex()} {z.hex()}\n"
                    e, g = hpref_run([line, glibc])
                    assert g == h, f"hpref (glibc) != healpy at nside {nside} {it}"
                    if e != h:
                        ok = False
                    rec[key] = e
                if not ok:
                    dropped += 1
                    continue
                kept.append(rec)
                if len(kept) == SWITCH_KEEP:
                    break
            assert len(kept) >= 10, f"branch_switch {name} nside {nside}: only {len(kept)} inputs where the branch matters"
            if "theta" in (kept[0] if kept else {}):
                ang += kept
            else:
                vec += kept
    dump("branch_switch.json", {"ang": ang, "vec": vec})
    print(f"  branch_switch: {len(ang)} Ang2Pix and {len(vec)} Vec2Pix inputs where the cos/sin branch "
          f"changes the pixel; {dropped} dropped where healpy (glibc) differs from correctly rounded healpix_cxx")


# ---------------------------------------------------------------------------
# Added by the round-2 second pass.
# ---------------------------------------------------------------------------
TWOTHIRD_NSIDES = [2**k for k in range(30)]


def write_ang2pix_twothird():
    """Ang2Pix at z = cos(theta) exactly 2/3 (the double nearest 2/3; the
    only such theta is 0.8410686705679303),
    where loc2pix's 'za <= twothird' sends NESTED and RING to the
    equatorial formula and 'za < twothird' would send them to the polar
    one, plus the doubles 1 and 2 ULP either side, for z = +2/3 and -2/3
    and phi = k pi/32, at every power-of-2 Nside 2^0..2^29. Only thetas
    whose glibc cos (what healpy calls) is the correctly rounded cos are
    kept, so healpy's answer is the one for the correctly rounded cos (the
    Go package's math.Cos can be 1 ULP off, which moves z across 2/3: the
    test accepts that only as a boundary-ambiguous difference). Explicit
    lists, not digests, so the test sees every expected value."""
    th = []
    for target in (2.0 / 3.0, -2.0 / 3.0):
        t0 = math.acos(target)
        cands = [t0]
        lo = hi = t0
        for _ in range(8):
            lo, hi = float(np.nextafter(lo, -np.inf)), float(np.nextafter(hi, np.inf))
            cands += [lo, hi]
        # No double theta near acos(-2/3) has a cos that rounds to -2/3
        # (cos steps by about 3 ULP of z per ULP of theta there); its
        # nearest neighbours are used instead.
        exact = [t for t in cands if cr(mpmath.cos, t) == target]
        assert exact or target < 0, f"no double theta with correctly rounded cos == {target}"
        t = exact[0] if exact else t0
        pts = [t]
        a = b = t
        for _ in range(2):
            a, b = float(np.nextafter(a, -np.inf)), float(np.nextafter(b, np.inf))
            pts += [a, b]
        th += sorted(pts)
    th = [t for t in th if math.cos(t) == cr(mpmath.cos, t)]
    phis = [k * math.pi / 32 for k in range(64)]
    theta = [t for t in th for _ in phis]
    phi = [p for _ in th for p in phis]
    results = []
    for nside in TWOTHIRD_NSIDES:
        results.append({"nside": nside,
                         "ring": [int(v) for v in hp.ang2pix(nside, np.array(theta), np.array(phi))],
                         "nest": [int(v) for v in hp.ang2pix(nside, np.array(theta), np.array(phi), nest=True)]})
    dump("ang2pix_twothird.json", {"theta": f64(theta), "phi": f64(phi),
                                   "cos_cr": f64([cr(mpmath.cos, t) for t in theta]), "results": results})
    print(f"  ang2pix_twothird: {len(th)} thetas x {len(phis)} phi x {len(TWOTHIRD_NSIDES)} Nside")


# Every exponent from -490 to -480 (healpy first moves a vector at 2^-482),
# and -495, -500.
SCALED_HEALPY_EXPS = list(range(-490, -479)) + [-495, -500]


def write_vec2pix_scaled_healpy():
    """healpy's own vec2pix on the ang2pix_polar.json vectors scaled by
    2^k (exact): for each Nside, scheme and k, the indices where healpy's
    answer for the scaled vector differs from its answer for the unscaled
    one, with the scaled answer. This is the truth for the Vec2Pix doc's
    statement about healpy (near-pole vectors scaled by 2^-490 move)."""
    import warnings
    th_all, ph_all = polar_directions()
    vx = [cr(lambda a, b: mpmath.sin(a) * mpmath.cos(b), a, b) for a, b in zip(th_all, ph_all)]
    vy = [cr(lambda a, b: mpmath.sin(a) * mpmath.sin(b), a, b) for a, b in zip(th_all, ph_all)]
    vz = [cr(mpmath.cos, a) for a in th_all]
    ax, ay, az = np.array(vx), np.array(vy), np.array(vz)
    cases = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for nside, nest in ((2**29 - 1, False), (2**29, False), (2**29, True)):
            base = hp.vec2pix(nside, ax, ay, az, nest=nest)
            for k in SCALED_HEALPY_EXPS:
                sx, sy, sz = np.ldexp(ax, k), np.ldexp(ay, k), np.ldexp(az, k)
                ok = (np.ldexp(sx, -k) == ax) & (np.ldexp(sy, -k) == ay) & (np.ldexp(sz, -k) == az)
                got = hp.vec2pix(nside, sx, sy, sz, nest=nest)
                diffs = [[int(i), int(got[i])] for i in range(len(vx)) if ok[i] and got[i] != base[i]]
                cases.append({"nside": nside, "nest": nest, "k": k, "exact": int(ok.sum()), "healpy_differs": diffs})
        # The Vec2Pix doc's example: (1e-200, 0, 1e-200), whose squares
        # underflow to 0 in healpy, against (1, 0, 1) (the same direction).
        tiny = []
        for nside in (1, 2, 16, 2**20, 2**29):
            for nest in (False, True):
                tiny.append({"nside": nside, "nest": nest,
                             "unit": int(hp.vec2pix(nside, 1.0, 0.0, 1.0, nest=nest)),
                             "tiny": int(hp.vec2pix(nside, 1e-200, 0.0, 1e-200, nest=nest))})
    dump("vec2pix_scaled_healpy.json", {"cases": cases, "tiny": tiny})
    print("  vec2pix_scaled_healpy: healpy moves " +
          ", ".join(f"{len(c['healpy_differs'])}@{c['nside']}{'n' if c['nest'] else 'r'}2^{c['k']}" for c in cases))


if __name__ == "__main__":
    write_ang2pix_random()
    write_ang2pix_polar()
    write_branch_switch()
    write_vec2pix_extreme()
    write_ang2pix_phiwrap()
    write_pixfuncs_ext()
    write_pix2vec()
    write_vec2pix()
    write_neighbours_edges()
    write_utilities_ext()
    write_ang2pix_twothird()
    write_vec2pix_scaled_healpy()
    nfiles = len([f for f in os.listdir(OUTDIR) if f.endswith(".json")])
    print(f"Done: {nfiles} golden files in {OUTDIR}.")
