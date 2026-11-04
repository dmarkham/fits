#!/usr/bin/env python3
"""Correctly rounded sin, cos, asin, acos, atan2 reference values (mpmath).

Writes libm/trig.json + .bin and libm/atan2.json + .bin (or under
<outroot>/libm; format in
f64pack.py: metadata JSON and the exact doubles, zlib): for a fixed set of
double arguments, the double nearest to the exact sin/cos/asin/acos,
computed by mpmath at 200 bits and rounded once. TestLibmAccuracy
(package transform) bounds Go's math.Sin/Cos and the port's
libAsin/libAcos against these doubles in ULP, TestAtan2Accuracy
math.Atan2 and math.Atan.

The arguments are the ones that matter for the projections: radian values
as wcslib forms them from degrees (fl(deg * fl(pi/180)), e.g. 35, 55,
89.99 deg), arguments within 1e-16..1e-1 of 0, +-pi/2 and +-pi, uniform
ones over [-7, 7], and asin/acos arguments uniform over [-1, 1] and
within 1e-16..1e-1 of +-1. Every random draw comes from a seeded
random.Random, so a re-run reproduces the file byte for byte (checked by
TestLibmGeneratorReproduces).

Usage: gen_libm.py [outroot]     (default: this directory)
Needs mpmath (pinned in PROVENANCE.md; the shared venv has 1.4.1).
"""
import math
import os
import random
import sys

import mpmath

import f64pack

mpmath.mp.prec = 200


def nearest(v):
    # float(mpf) rounds to nearest, ties to even.
    return float(v)


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    rng = random.Random(20261101)
    d2r = math.pi / 180.0  # wcsmath.h D2R = PI/180.0 in double arithmetic

    args = []
    # Degree angles as wcslib passes them to sin/cos.
    for deg in [1, 10, 15, 30, 35, 45, 55, 60, 75, 80, 89, 89.9, 89.99, 89.999,
                89.9999, 89.99999, 89.999999, 89.9999999, 135, 179.9, 179.99999,
                180.0000001, 270.0000001, -89.99, -89.9999999, 0.001, 1e-7]:
        for s in (1, -1):
            args.append(s * deg * d2r)
    for _ in range(400):
        args.append(rng.uniform(-360, 360) * d2r)
    # Near 0, +-pi/2, +-pi.
    for c in (0.0, math.pi / 2, -math.pi / 2, math.pi, -math.pi):
        for _ in range(240):
            off = 10 ** rng.uniform(-16, -1) * rng.choice((1, -1))
            args.append(c + off)
    for _ in range(400):
        args.append(rng.uniform(-7, 7))

    inv = []
    for _ in range(500):
        inv.append(rng.uniform(-1, 1))
    for _ in range(400):
        off = 10 ** rng.uniform(-16, -1)
        inv.append(rng.choice((1, -1)) * (1 - off))
    inv += [1.0, -1.0, 0.0, 0.5, -0.5, 1 - 2 ** -53, -(1 - 2 ** -53), 1e-300]

    trig = []
    for a in args:
        m = mpmath.mpf(a)
        trig.append([a, nearest(mpmath.sin(m)), nearest(mpmath.cos(m))])
    arc = []
    for x in inv:
        m = mpmath.mpf(x)
        arc.append([x, nearest(mpmath.asin(m)), nearest(mpmath.acos(m))])

    meta = {
        "generator": "gen_libm.py",
        "mpmath": mpmath.__version__,
        "precision_bits": 200,
        "layout": "sincos rows [x, sin x, cos x], then asinacos rows [x, asin x, acos x]",
        "sincos": len(trig),
        "asinacos": len(arc),
    }
    values = [v for r in trig + arc for v in r]
    path = f64pack.write(os.path.join(root, "libm"), "trig", meta, values)
    print("wrote %s: %d sin/cos, %d asin/acos" % (path, len(trig), len(arc)))

    # atan2(y, x) (atand and atan2d are built on it). The
    # arguments wcslib passes: (sin phi, cos phi) of degree angles within
    # 1e-12..1e-1 deg of 0, +-90 and +-180 (qscs2x's small-angle zeco takes
    # atan2 of them and subtracts pi/2 or pi) and over the whole circle;
    # atan2d(1, s) and atand(v) = atan2(v, 1) for s, v over 1e-6..1e6;
    # random (y, x) of every sign with magnitudes 1e-3..1e3.
    pairs = []
    for c in (0.0, 90.0, -90.0, 180.0, -180.0):
        for _ in range(240):
            deg = c + 10 ** rng.uniform(-12, -1) * rng.choice((1, -1))
            a = mpmath.mpf(deg * d2r)
            pairs.append((nearest(mpmath.sin(a)), nearest(mpmath.cos(a))))
    for _ in range(300):
        a = mpmath.mpf(rng.uniform(-180, 180) * d2r)
        pairs.append((nearest(mpmath.sin(a)), nearest(mpmath.cos(a))))
    for _ in range(300):
        v = 10 ** rng.uniform(-6, 6) * rng.choice((1, -1))
        pairs.append((1.0, v))
        pairs.append((v, 1.0))
    for _ in range(400):
        pairs.append((10 ** rng.uniform(-3, 3) * rng.choice((1, -1)),
                      10 ** rng.uniform(-3, 3) * rng.choice((1, -1))))
    pairs += [(1.0, 1.0), (-1.0, 1.0), (1.0, -1.0), (-1.0, -1.0), (1e-300, 1.0), (1.0, 1e-300)]
    at = []
    for y, x in pairs:
        at.append([y, x, nearest(mpmath.atan2(mpmath.mpf(y), mpmath.mpf(x)))])
    meta = {
        "generator": "gen_libm.py",
        "mpmath": mpmath.__version__,
        "precision_bits": 200,
        "layout": "rows [y, x, atan2(y, x)]",
        "atan2": len(at),
    }
    path = f64pack.write(os.path.join(root, "libm"), "atan2", meta, [v for r in at for v in r])
    print("wrote %s: %d atan2" % (path, len(at)))


if __name__ == "__main__":
    main()
