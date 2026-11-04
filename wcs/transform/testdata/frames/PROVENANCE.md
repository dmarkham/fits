# Sky-frame truth: provenance

`frames_golden.json` is written by `gen_frames.py` in this directory. No
value in it comes from the Go code under test.

| Tool | Version | Used for |
|---|---|---|
| astropy (`astropy.coordinates`) | 8.0.2 | ICRS / FK5 J2000 / FK4 B1950 / Galactic / Supergalactic / BarycentricMeanEcliptic transforms |
| pyerfa (ERFA C library) | 2.0.1.5 (ERFA 2.0.1) | IAU 1976 precession: `erfa.prec76` angles, matrix built with `erfa.ir/rz/ry`, applied with `erfa.s2c/rxp/c2s` |
| numpy | 2.4.6 | |
| Python | 3.11.11 (pyenv), venv at `fits-proof-audit/venv` | |

Regenerate (from the repo root):

    venv/bin/python wcs/transform/testdata/frames/gen_frames.py \
        > wcs/transform/testdata/frames/frames_golden.json

The output is deterministic; `TestFramesGoldenReproduce` re-runs the
generator and requires byte-identical output (tier 2, `go test -tags
interop`: it fails without `FITS_ASTROPY_PYTHON` pointing at the venv
interpreter). `TestFramesAstropyGolden`
needs only the committed JSON and fails if it is missing or was produced by
other astropy / pyerfa versions.

## What each case is

The same 144 input directions are used by every case; the golden stores
them once (top-level `points`) and each case holds only its `out` array,
in the same order (compact JSON, one case per line): a 30 deg x 15 deg grid including both poles, plus
near-pole points (|lat| = 89.9999), lon just above 0, and named objects
(Galactic centre, FK5 galactic pole, M31, M1, Cen A, Cyg A, theta Persei).

| Case | astropy / ERFA call | Go function |
|---|---|---|
| `icrs_to_galactic`, `galactic_to_icrs` | ICRS <-> Galactic (frame bias + FK5 J2000 galactic rotation) | `EquatorialToGalactic`, `GalacticToEquatorial` |
| `fk5_to_galactic`, `galactic_to_fk5` | FK5(equinox=J2000) <-> Galactic | `FK5ToGalactic`, `GalacticToFK5` |
| `galactic_to_supergalactic`, `supergalactic_to_galactic` | Galactic <-> Supergalactic | `GalacticToSupergalactic`, `SupergalacticToGalactic` |
| `icrs_to_supergalactic`, `supergalactic_to_icrs` | ICRS <-> Supergalactic | `EquatorialToSupergalactic`, `SupergalacticToEquatorial` |
| `icrs_to_ecliptic`, `ecliptic_to_icrs` (jd = J2000, B1950, 2025-01-01, 2100-01-01, B1900; TT) | ICRS <-> BarycentricMeanEcliptic(equinox=jd) (= ERFA `ecm06`) | `EquatorialToEcliptic`, `EclipticToEquatorial` |
| `precess_fk5` (1950->2000, 2000->1950, 2000->2050, 1900->2100, 2025.5->1987.25, Julian epochs) | ERFA `prec76` + Rz(-zeta) Ry(theta) Rz(-z) | `PrecessFK5` |
| `fk4_to_fk5`, `fk5_to_fk4` (obstime B1950, B1900, B1984, B2000) | FK4(B1950, obstime) <-> FK4NoETerms(B1950, obstime) <-> FK5(J2000) | `FK4ToFK5`, `FK5ToFK4` |

Note on `fk4_to_fk5`: the generator goes through FK4NoETerms(equinox=B1950)
explicitly. A one-step `SkyCoord(FK4(...)).transform_to(FK5())` makes
astropy build the intermediate FK4NoETerms frame with the FK5 frame's
equinox (J2000), i.e. it precesses B1950 -> J2000 -> B1950 with Newcomb's
formulas, which are not exact inverses; that detour shifts every point by a
constant 2.5e-5 arcsec. The explicit path is astropy's own definition of the
transform (`fk4_to_fk4_no_e` then `fk4_no_e_to_fk5`) without the round trip.

Measured agreement (Go vs this golden): worst 4.1e-10 arcsec over all
31 cases; budget 1e-8 arcsec.

Published references checked separately in `frames_test.go` (no library):
Meeus, *Astronomical Algorithms* 2nd ed., example 21.b (FK5 precession,
to half a unit of the printed digits); IAU 2006 obliquity at J2000 =
84381.406"; galactic pole / NCP from the FK5 J2000 galactic constants;
supergalactic pole (47.37, +6.32) and origin (137.37, 0).
