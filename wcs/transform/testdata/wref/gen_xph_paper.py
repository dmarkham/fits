#!/usr/bin/env python3
"""XPH truth from the paper's equations, independent of wcslib.

Writes paper/xph_paper.json + .bin (or under <outroot>/paper; format in
f64pack.py: metadata JSON and the exact doubles, zlib): the XPH
("butterfly") projection evaluated straight from the equations of

  Calabretta M. R. & Lowe S. R. 2013, PASA 30, e050,
  "Representing the 'butterfly' projection in FITS - projection code XPH",
  arXiv:1308.6647v1 (PDF sha256
  fd5ff0cd630a713aaedd3029c08b67b6439c2ce28e0c1a62569e7d80c301a9e5),

projection equations (2)-(9) and deprojection equations (10)-(17), in
mpmath at 60 significant digits and rounded once to double. The small-angle
forms (18)/(19) are numerical devices for double precision; at 60 digits
the exact (9) and (16) are used everywhere. No wcslib code, and no code of
the Go package, is involved: TestXPHPaperEquations holds the Go XPH port to
these values.

forward rows [phi, theta, x, y] on a 7.5 deg grid (phi in [-180, 180)), plus
theta = +-(90 - m 10^-k), m = 1, 3, k = 1..7, theta = +-89.95 .. +-89.996
(either side of sigma = 1e-4, where the double-precision code switches to
(19)), and theta = +-(asin(2/3) +- 1e-7 deg) (either side of the zone
boundary).
inverse rows [x, y, phi, theta] at the double-rounded image of every forward
point except the south pole (so every input lies on the projection),
with phi normalised to [-180, 180), plus direct inverse probes at
sigma = m 10^-k (k = 1..10) next to both poles.

Usage: gen_xph_paper.py [outroot]    (needs mpmath; venv 1.4.1)
"""
import os
import sys

import mpmath as mp

import f64pack

mp.mp.dps = 60
D = mp.mpf
KAPPA = mp.sqrt(2) / 2
DEG = mp.pi / 180


def sind(a):
    return mp.sin(a * DEG)


def asind(v):
    return mp.asin(v) / DEG


def forward(phi, theta):
    """Equations (2)-(9); phi in [-180, 180)."""
    phi, theta = D(phi), D(theta)
    psi = mp.fmod(phi + 180, 90)
    thx = asind(D(2) / 3)
    if abs(theta) <= thx:
        xi = psi
        eta = D(135) / 2 * sind(theta)
    else:
        sigma = mp.sqrt(3 * (1 - abs(sind(theta))))
        xi = 45 + (psi - 45) * sigma
        eta = mp.sign(theta) * (90 - 45 * sigma)
    xp, ep = xi - 45, eta - 90
    if phi < -90:
        return KAPPA * (-xp + ep), KAPPA * (-xp - ep)
    if phi < 0:
        return KAPPA * (xp + ep), KAPPA * (-xp + ep)
    if phi < 90:
        return KAPPA * (xp - ep), KAPPA * (xp + ep)
    return KAPPA * (-xp - ep), KAPPA * (xp - ep)


def inverse(x, y):
    """Equations (10)-(17)."""
    x, y = D(x), D(y)
    if x <= 0 and y > 0:
        xp, ep, off = KAPPA * (-x - y), KAPPA * (x - y), -180
    elif x < 0 and y <= 0:
        xp, ep, off = KAPPA * (x - y), KAPPA * (x + y), -90
    elif x >= 0 and y < 0:
        xp, ep, off = KAPPA * (x + y), KAPPA * (-x + y), 0
    else:
        xp, ep, off = KAPPA * (-x + y), KAPPA * (-x - y), 90
    xi, eta = xp + 45, ep + 90
    if abs(eta) <= 45:
        psi = xi
        theta = asind(2 * eta / 135)
    else:
        sigma = (90 - abs(eta)) / 45
        psi = 45 + (xi - 45) / sigma if sigma != 0 else D(45)
        theta = mp.sign(eta) * asind(1 - sigma * sigma / 3)
    phi = psi + off
    if phi >= 180:
        phi -= 360
    if phi < -180:
        phi += 360
    return phi, theta


def fmt(v):
    return repr(float(v))


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    pts = []
    for i in range(48):
        phi = -180 + 7.5 * i
        for j in range(25):
            pts.append((phi, -90 + 7.5 * j))
    thx = float(asind(D(2) / 3))
    for phi in (0.0, 30.0, 97.5, -150.0, -45.0, 179.9):
        for k in range(1, 8):
            for m in (1.0, 3.0):
                for s in (1, -1):
                    pts.append((phi, s * (90 - m * 10.0 ** -k)))
        # Either side of the port's (and wcslib's) switch to the small-angle
        # form (19) at sigma = 1e-4, theta = 89.99532 deg: a moved switch
        # shows here as the truncation error of (19).
        for t in (89.95, 89.96, 89.97, 89.98, 89.994, 89.996):
            for s in (1, -1):
                pts.append((phi, s * t))
        for s in (1, -1):
            pts.append((phi, s * (thx + 1e-7)))
            pts.append((phi, s * (thx - 1e-7)))
    fwd, inv = [], []
    for phi, theta in pts:
        x, y = forward(phi, theta)
        fwd.append([phi, theta, x, y])
        if theta == -90:
            # The south pole is the four outer corners of the butterfly;
            # rounded to double they fall on either side of the map edge,
            # so the inverse there is not a meaningful test.
            continue
        xd, yd = float(x), float(y)
        p, t = inverse(xd, yd)
        inv.append([xd, yd, p, t])
    # Direct inverse probes next to both poles: sigma = m 10^-k (k = 1..10,
    # m = 1, 3) in the quadrant x >= 0, y < 0 (eqs. 10-13 with
    # xi' = x' - 45 = kappa (x + y), eta' = eta - 90 = kappa (-x + y)), on
    # the facet centre line (xi = 45) and halfway to its edge. Near the
    # north pole these reach the small-angle form (18); near the south
    # pole they reach the double-precision limit of (16), which peaks where
    # sigma^2 / 3 approaches 2^-53 and which the images of the forward
    # points do not reach.
    for k in range(1, 11):
        for m in (1.0, 3.0):
            sig = D(m * 10.0 ** -k)
            for sgn in (1, -1):
                eta = sgn * (90 - 45 * sig)
                for xi in (D(45), 45 + D(45) / 2 * sig):
                    xpr, epr = xi - 45, eta - 90
                    xd = float((xpr - epr) / (2 * KAPPA))
                    yd = float((xpr + epr) / (2 * KAPPA))
                    p, t = inverse(xd, yd)
                    inv.append([xd, yd, p, t])
    meta = {
        "source": "Calabretta & Lowe 2013 PASA 30 e050 (arXiv:1308.6647v1) eqs (2)-(17), mpmath %s at 60 digits" % mp.__version__,
        "layout": "forward rows [phi, theta, x, y], then inverse rows [x, y, phi, theta]",
        "forward": len(fwd),
        "inverse": len(inv),
    }
    values = [v for r in fwd + inv for v in r]
    path = f64pack.write(os.path.join(root, "paper"), "xph_paper", meta, values)
    print("wrote %s: %d forward, %d inverse" % (path, len(fwd), len(inv)))


if __name__ == "__main__":
    main()
