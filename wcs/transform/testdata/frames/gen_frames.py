#!/usr/bin/env python3
"""Generate frames_golden.json: sky-frame conversion truth for frames.go.

Truth sources (no value here comes from the Go code under test):

  * astropy.coordinates (SkyCoord.transform_to) for ICRS/FK5/FK4 <->
    Galactic, Supergalactic, BarycentricMeanEcliptic and FK4 <-> FK5;
  * ERFA (pyerfa, the C SOFA port) for the IAU 1976 FK5 precession:
    erfa.prec76 gives (zeta, z, theta) and the rotation is built and applied
    with erfa.ir / erfa.rz / erfa.ry / erfa.rxp / erfa.s2c / erfa.c2s, i.e.
    exactly SOFA's iauPmat76 recipe generalised to two arbitrary epochs.

Run with the repo venv (versions pinned in PROVENANCE.md):

    venv/bin/python wcs/transform/testdata/frames/gen_frames.py > \
        wcs/transform/testdata/frames/frames_golden.json

Output is deterministic (fixed point list, repr() floats), so re-running
reproduces the committed bytes exactly with the same astropy/pyerfa/numpy.
"""
import json
import sys

import astropy
import astropy.units as u
import erfa
import numpy as np
from astropy.coordinates import (FK4, FK5, ICRS, FK4NoETerms, BarycentricMeanEcliptic,
                                 Galactic, SkyCoord, Supergalactic)
from astropy.time import Time


def points():
    pts = []
    for lat in range(-90, 91, 15):
        lons = [0.0] if abs(lat) == 90 else [float(l) for l in range(0, 360, 30)]
        for lon in lons:
            pts.append((lon, float(lat)))
    # Off-grid points, near-pole points, and the classic reference objects.
    pts += [
        (266.4168, -28.9360),      # near the Galactic centre
        (192.85948, 27.12825),     # Galactic north pole (FK5 J2000, rounded)
        (10.684708, 41.268750),    # M31
        (83.633083, 22.014500),    # M1
        (359.999999, 89.9999),
        (123.456, -89.9999),
        (0.000001, 0.0),
        (201.365063, -43.019112),  # Cen A
        (41.054063, 49.227750),    # theta Persei (Meeus ex. 21.b)
        (299.868, 40.734),         # Cyg A
    ]
    return pts


PTS = points()
LON = np.array([p[0] for p in PTS]) * u.deg
LAT = np.array([p[1] for p in PTS]) * u.deg


def lonlat(sc):
    sph = sc.represent_as("unitspherical")
    return [[float(a), float(b)] for a, b in
            zip(sph.lon.wrap_at(360 * u.deg).deg, sph.lat.deg)]


def case(name, out, **params):
    # The input directions are the same for every case; they are written
    # once, as the top-level "points", and each case holds only its output.
    assert len(out) == len(PTS)
    return {"name": name, "params": params, "out": out}


def jd_tt(jd):
    return Time(jd, format="jd", scale="tt")


def main():
    cases = []
    icrs = SkyCoord(ICRS(ra=LON, dec=LAT))
    fk5 = SkyCoord(FK5(ra=LON, dec=LAT, equinox="J2000"))
    gal = SkyCoord(Galactic(l=LON, b=LAT))
    sgal = SkyCoord(Supergalactic(sgl=LON, sgb=LAT))

    cases.append(case("icrs_to_galactic", lonlat(icrs.transform_to(Galactic()))))
    cases.append(case("galactic_to_icrs", lonlat(gal.transform_to(ICRS()))))
    cases.append(case("fk5_to_galactic", lonlat(fk5.transform_to(Galactic()))))
    cases.append(case("galactic_to_fk5",
                      lonlat(gal.transform_to(FK5(equinox="J2000")))))
    cases.append(case("galactic_to_supergalactic",
                      lonlat(gal.transform_to(Supergalactic()))))
    cases.append(case("supergalactic_to_galactic",
                      lonlat(sgal.transform_to(Galactic()))))
    cases.append(case("icrs_to_supergalactic",
                      lonlat(icrs.transform_to(Supergalactic()))))
    cases.append(case("supergalactic_to_icrs", lonlat(sgal.transform_to(ICRS()))))

    # Mean ecliptic and equinox of date (IAU 2006, erfa.ecm06 inside astropy).
    for jd in (2451545.0, 2433282.42345905, 2460676.5, 2488069.5, 2415020.31352):
        ecl = BarycentricMeanEcliptic(equinox=jd_tt(jd))
        cases.append(case("icrs_to_ecliptic", lonlat(icrs.transform_to(ecl)), jd=jd))
        e = SkyCoord(BarycentricMeanEcliptic(lon=LON, lat=LAT, equinox=jd_tt(jd)))
        cases.append(case("ecliptic_to_icrs", lonlat(e.transform_to(ICRS())), jd=jd))

    # IAU 1976 precession between Julian epochs (ERFA).
    for ep1, ep2 in ((1950.0, 2000.0), (2000.0, 1950.0), (2000.0, 2050.0),
                     (1900.0, 2100.0), (2025.5, 1987.25)):
        d1 = (ep1 - 2000.0) * 365.25
        d2 = (ep2 - 2000.0) * 365.25
        zeta, z, theta = erfa.prec76(2451545.0, d1, 2451545.0, d2)
        r = erfa.ir()
        r = erfa.rz(-zeta, r)
        r = erfa.ry(theta, r)
        r = erfa.rz(-z, r)
        out = []
        for lon, lat in PTS:
            v = erfa.s2c(np.radians(lon), np.radians(lat))
            a, b = erfa.c2s(erfa.rxp(r, v))
            out.append([float(np.degrees(erfa.anp(a))), float(np.degrees(b))])
        cases.append(case("precess_fk5", out, equinox1=ep1, equinox2=ep2))

    # FK4 (equinox B1950, with E-terms) <-> FK5 J2000.
    #
    # The hop through FK4NoETerms(equinox=B1950) is spelled out. A direct
    # FK4 -> FK5 transform_to builds its intermediate FK4NoETerms frame with
    # the FK5 frame's equinox (J2000), so astropy precesses B1950 -> J2000
    # and back with Newcomb's formulas, which are not exact inverses: that
    # detour moves every point by a constant ~2.5e-5 arcsec. The explicit
    # path is astropy's own definition (fk4_to_fk4_no_e, fk4_no_e_to_fk5)
    # without the round trip.
    for obs in (1950.0, 1900.0, 1984.0, 2000.0):
        ot = Time(obs, format="byear", scale="tt")
        noe = FK4NoETerms(equinox="B1950", obstime=ot)
        fk4 = SkyCoord(FK4(ra=LON, dec=LAT, equinox="B1950", obstime=ot))
        cases.append(case("fk4_to_fk5",
                          lonlat(fk4.transform_to(noe).transform_to(FK5(equinox="J2000"))),
                          obstime_byear=obs))
        cases.append(case("fk5_to_fk4",
                          lonlat(fk5.transform_to(noe).transform_to(
                              FK4(equinox="B1950", obstime=ot))),
                          obstime_byear=obs))

    doc = {
        "generator": "wcs/transform/testdata/frames/gen_frames.py",
        "astropy": astropy.__version__,
        "pyerfa": erfa.__version__,
        "erfa": erfa.version.erfa_version,
        "numpy": np.__version__,
        "points": [[p[0], p[1]] for p in PTS],
        "cases": cases,
    }
    # One case per line keeps `git diff` of a regeneration readable; no
    # blanks after separators (compact, the values are unchanged).
    sep = (",", ":")
    head = {k: v for k, v in doc.items() if k != "cases"}
    out = json.dumps(head, separators=sep)[:-1] + ',"cases":[\n'
    out += ",\n".join(json.dumps(c, separators=sep) for c in cases)
    out += "\n]}\n"
    sys.stdout.write(out)


if __name__ == "__main__":
    main()
