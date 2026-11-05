#!/usr/bin/env python3
"""HPX truth from the paper's equations, independent of wcslib.

Writes paper/hpx_paper.json + .bin (or under <outroot>/paper; format in
f64pack.py: metadata JSON and the exact doubles, zlib): the HPX
(HEALPix) projection with H facets in longitude and K in latitude,
evaluated straight from the equations of

  Calabretta M. R. & Roukema B. F. 2007, MNRAS 381, 865,
  "Mapping on the HEALPix grid", arXiv:0704.0626, section 6 (the HPX
  projection for general H and K):

  equatorial zone, |sin theta| <= (K - 1) / K:
      x = phi,   y = (90 K / H) sin theta
  polar zones:
      sigma = sqrt(K (1 - |sin theta|))
      phi_c = -180 + (2 floor((phi + 180) H / 360) + 1) 180 / H
      x = phi_c + (phi - phi_c) sigma
      y = sign(theta) (180 / H) ((K + 1) / 2 - sigma)
  and the inverses
      equatorial, |y| <= 90 (K - 1) / H:
          phi = x,   theta = asin(y H / (90 K))
      polar:
          sigma = (K + 1) / 2 - |y| H / 180
          x_c = -180 + (2 floor((x + 180) H / 360) + 1) 180 / H
          phi = x_c + (x - x_c) / sigma
          theta = sign(y) asin(1 - sigma^2 / K)

in mpmath at 60 significant digits, rounded once to double. Only odd K is
generated: for even K the southern polar facets are shifted by 180 / H, a
layout detail the wcslib goldens cover (hpx_h3_k2); the near-pole
precision, which is what this file is for, is the same formula for every K.
No wcslib code and no Go code is involved: TestHPXPaperEquations holds the
Go HPX port to these values.

Cases: (H, K) = (4, 3) (the HEALPix default) on a 7.5 deg grid, and
(5, 5) (odd H, so facet centres are not at multiples of 45 deg) on a
15 deg grid. Both also get theta = +-(90 - m 10^-k), m = 1, 3, k = 1..8,
at six phi,
and theta = +-(asin((K - 1) / K) +- 1e-7 deg), either side of the zone
boundary.
forward rows [phi, theta, x, y]; phi in [-180, 180).
inverse rows [x, y, phi, theta] at the double-rounded image of every forward
point that is not a pole (the pole is a facet apex, where phi is
undefined) and not on a polar facet edge (phi_c +- 180 / H), plus direct
inverse probes at sigma = m 10^-k (k = 1..10) next to the polar apexes.

Usage: gen_hpx_paper.py [outroot]    (needs mpmath; venv 1.4.1)
"""
import os
import sys

import mpmath as mp

import f64pack

mp.mp.dps = 60
D = mp.mpf
DEG = mp.pi / 180


def sind(a):
    return mp.sin(a * DEG)


def asind(v):
    return mp.asin(v) / DEG


def centre(v, H):
    return -180 + (2 * mp.floor((v + 180) * H / 360) + 1) * D(180) / H


def forward(phi, theta, H, K):
    phi, theta = D(phi), D(theta)
    s = sind(theta)
    if abs(s) <= D(K - 1) / K:
        return phi, D(90) * K / H * s
    sigma = mp.sqrt(K * (1 - abs(s)))
    pc = centre(phi, H)
    return pc + (phi - pc) * sigma, mp.sign(theta) * D(180) / H * (D(K + 1) / 2 - sigma)


def inverse(x, y, H, K):
    x, y = D(x), D(y)
    if abs(y) <= D(90) * (K - 1) / H:
        return x, asind(y * H / (D(90) * K))
    sigma = D(K + 1) / 2 - abs(y) * H / 180
    xc = centre(x, H)
    return xc + (x - xc) / sigma, mp.sign(y) * asind(1 - sigma * sigma / K)


def fmt(v):
    return repr(float(v))


def rows_for(H, K, step):
    pts = []
    n = int(360 / step)
    for i in range(n):
        phi = -180 + step * i
        for j in range(int(180 / step) + 1):
            pts.append((phi, -90 + step * j))
    thx = float(asind(D(K - 1) / K))
    for phi in (0.0, 30.0, 97.5, -150.0, -45.0, 179.9):
        for k in range(1, 9):
            for m in (1.0, 3.0):
                for s in (1, -1):
                    pts.append((phi, s * (90 - m * 10.0 ** -k)))
        if K > 1:
            for s in (1, -1):
                pts.append((phi, s * (thx + 1e-7)))
                pts.append((phi, s * (thx - 1e-7)))
    fwd, inv = [], []
    for phi, theta in pts:
        x, y = forward(phi, theta, H, K)
        fwd.append([phi, theta, x, y])
        if abs(theta) == 90:
            continue
        if abs(sind(theta)) > D(K - 1) / K and ((phi + 180) * H / 360) % 1 == 0:
            # phi on a polar facet edge maps to the edge of the facet; its
            # double-rounded image lies on either side of it (wcslib and
            # the port reject some), a knife edge, not a test.
            continue
        xd, yd = float(x), float(y)
        p, t = inverse(xd, yd, H, K)
        inv.append([xd, yd, p, t])
    # Direct inverse probes next to every polar apex: sigma = m 10^-k
    # (k = 1..10, m = 1, 3), at the facet centre line and halfway to the
    # facet edge. The double-rounding limit of the inverse peaks where
    # sigma^2 / K approaches 2^-53 (sigma ~ 1e-8), which the images of
    # the forward points above do not reach.
    for k in range(1, 11):
        for m in (1.0, 3.0):
            sig = m * 10.0 ** -k
            for s in (1, -1):
                y = float(s * D(180) / H * (D(K + 1) / 2 - D(sig)))
                for xc in (float(centre(D(0), H)), float(centre(D(-170), H))):
                    for frac in (0.0, 0.5):
                        x = float(D(xc) + D(frac) * D(180) / H * D(sig))
                        p, t = inverse(x, y, H, K)
                        inv.append([x, y, p, t])
    return fwd, inv


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    cases, values = [], []
    for H, K, step in ((4, 3, 7.5), (5, 5, 15.0)):
        fwd, inv = rows_for(H, K, step)
        cases.append({"h": H, "k": K, "forward": len(fwd), "inverse": len(inv)})
        for r in fwd + inv:
            values.extend(r)
        print("H=%d K=%d: %d forward, %d inverse" % (H, K, len(fwd), len(inv)))
    meta = {
        "source": "Calabretta & Roukema 2007 MNRAS 381 865 (arXiv:0704.0626) sec. 6, mpmath %s at 60 digits" % mp.__version__,
        "layout": "per case, in order: forward rows [phi, theta, x, y], then inverse rows [x, y, phi, theta]",
        "cases": cases,
    }
    path = f64pack.write(os.path.join(root, "paper"), "hpx_paper", meta, values)
    print("wrote %s" % path)


if __name__ == "__main__":
    main()
