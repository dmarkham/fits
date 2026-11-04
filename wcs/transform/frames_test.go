package transform

import (
	"math"
	"testing"
)

// Published-reference and definitional checks for frames.go. The
// point-by-point proof against astropy / ERFA is TestFramesAstropyGolden
// (frames_golden_test.go); these pin the defining points of each frame
// independently of any library.

// lonDiff returns the smallest absolute difference between two longitudes,
// accounting for the 360° wrap-around.
func lonDiff(a, b float64) float64 {
	d := math.Abs(a - b)
	if d > 180 {
		d = 360 - d
	}
	return d
}

// TestRoundTrips: every conversion pair is self-inverse to float64
// rounding (1e-9 arcsec) everywhere, including the poles.
func TestRoundTrips(t *testing.T) {
	type pair struct {
		name     string
		fwd, inv func(a, b float64) (float64, float64)
	}
	const jd1950 = 2433282.42345905
	pairs := []pair{
		{"galactic", EquatorialToGalactic, GalacticToEquatorial},
		{"fk5-galactic", FK5ToGalactic, GalacticToFK5},
		{"supergalactic", GalacticToSupergalactic, SupergalacticToGalactic},
		{"icrs-supergalactic", EquatorialToSupergalactic, SupergalacticToEquatorial},
		{"ecliptic J2000", func(a, d float64) (float64, float64) { return EquatorialToEcliptic(a, d, 2451545) },
			func(a, d float64) (float64, float64) { return EclipticToEquatorial(a, d, 2451545) }},
		{"ecliptic B1950", func(a, d float64) (float64, float64) { return EquatorialToEcliptic(a, d, jd1950) },
			func(a, d float64) (float64, float64) { return EclipticToEquatorial(a, d, jd1950) }},
		{"precess", func(a, d float64) (float64, float64) { return PrecessFK5(a, d, 1950, 2050) },
			func(a, d float64) (float64, float64) { return PrecessFK5(a, d, 2050, 1950) }},
		{"fk4", func(a, d float64) (float64, float64) { return FK4ToFK5(a, d, 1950) },
			func(a, d float64) (float64, float64) { return FK5ToFK4(a, d, 1950) }},
	}
	for _, p := range pairs {
		for _, c := range [][2]float64{{0, 0}, {180, 0}, {90, 45}, {270, -30}, {45, 60}, {359, -89}, {12, 90}, {200, -90}} {
			x, y := p.fwd(c[0], c[1])
			a, d := p.inv(x, y)
			// Precession is not exactly self-inverse (the 1976 angles for
			// 1950->2050 and 2050->1950 are separate series); allow its
			// documented series truncation, 1e-6 arcsec.
			tol := 1e-9
			if p.name == "precess" {
				tol = 1e-6
			}
			if sep := angSepArcsec(a, d, c[0], c[1]); sep > tol {
				t.Errorf("%s round trip %v -> (%v,%v) -> (%v,%v): %.3g arcsec", p.name, c, x, y, a, d, sep)
			}
		}
	}
}

// TestGalacticDefiningPoints: the FK5 J2000 galactic north pole maps to
// b = +90 and the north celestial pole to l = lon0 (astropy
// Galactic._ngp_J2000 / _lon0_J2000, the IAU 1958 definition in FK5).
func TestGalacticDefiningPoints(t *testing.T) {
	if _, b := FK5ToGalactic(galNGPRaJ2000, galNGPDecJ2000); math.Abs(b-90) > 1e-12 {
		t.Errorf("galactic latitude of the NGP = %.15f, want 90", b)
	}
	l, b := FK5ToGalactic(0, 90)
	if math.Abs(l-galLon0J2000) > 1e-11 || math.Abs(b-galNGPDecJ2000) > 1e-11 {
		t.Errorf("NCP -> (%.13f, %.13f), want (%.13f, %.13f)", l, b, galLon0J2000, galNGPDecJ2000)
	}
	// The ICRS version differs only by the frame bias (<= 0.03 arcsec).
	la, ba := EquatorialToGalactic(266.4168, -28.9360)
	lf, bf := FK5ToGalactic(266.4168, -28.9360)
	if sep := angSepArcsec(la, ba, lf, bf); sep < 0.005 || sep > 0.035 {
		t.Errorf("ICRS vs FK5 galactic at the GC differ by %.4f arcsec, want the 0.005-0.035 arcsec frame bias", sep)
	}
}

// TestSupergalacticDefiningPoints: the supergalactic north pole is at
// galactic (47.37, +6.32) and SGL = SGB = 0 at galactic (137.37, 0)
// (de Vaucouleurs et al. 1976; astropy Supergalactic).
func TestSupergalacticDefiningPoints(t *testing.T) {
	if _, sgb := GalacticToSupergalactic(47.37, 6.32); math.Abs(sgb-90) > 1e-9 {
		t.Errorf("SGB of the supergalactic pole = %.12f, want 90", sgb)
	}
	sgl, sgb := GalacticToSupergalactic(137.37, 0)
	if lonDiff(sgl, 0) > 1e-9 || math.Abs(sgb) > 1e-9 {
		t.Errorf("galactic (137.37, 0) -> (%.12f, %.12f), want (0, 0)", sgl, sgb)
	}
	l, b := SupergalacticToGalactic(0, 0)
	if lonDiff(l, 137.37) > 1e-9 || math.Abs(b) > 1e-9 {
		t.Errorf("SG origin -> galactic (%.12f, %.12f), want (137.37, 0)", l, b)
	}
}

// TestObliquityJ2000: the IAU 2006 mean obliquity at J2000.0 is exactly
// 84381.406 arcsec (Capitaine et al. 2003).
func TestObliquityJ2000(t *testing.T) {
	eps := obliquityIAU2006(2451545.0) * 180 * 3600 / math.Pi
	if math.Abs(eps-84381.406) > 1e-9 {
		t.Errorf("obliquity at J2000: %.9f arcsec, want 84381.406", eps)
	}
}

// TestPrecessMeeus21b reproduces Meeus, Astronomical Algorithms (2nd ed.),
// example 21.b: theta Persei, FK5 J2000.0 (2h44m11.986s, +49°13'42.48"),
// proper motion (+0.03425 s/yr, -0.0895"/yr), precessed to 2028 Nov 13.19
// TD (JD 2462088.69) gives (2h46m11.331s, +49°20'54.54"). The published
// values are rounded to 0.001 s and 0.01", so the budget is half a unit
// in the last printed place (0.0005 s, 0.005").
func TestPrecessMeeus21b(t *testing.T) {
	years := (2462088.69 - 2451545.0) / 365.25
	a0 := (2+44.0/60+11.986/3600)*15 + 0.03425*15/3600*years
	d0 := 49 + 13.0/60 + 42.48/3600 - 0.0895/3600*years
	a, d := PrecessFK5(a0, d0, 2000, 2000+years)
	wantA := (2 + 46.0/60 + 11.331/3600) * 15
	wantD := 49 + 20.0/60 + 54.54/3600
	if ds := math.Abs(a-wantA) * 3600 / 15; ds > 0.0005 {
		t.Errorf("RA %.9f deg off Meeus by %.5f s", a, ds)
	}
	if dd := math.Abs(d-wantD) * 3600; dd > 0.005 {
		t.Errorf("Dec %.9f deg off Meeus by %.4f arcsec", d, dd)
	}
}

// TestPrecessIdentity: precessing from an epoch to itself is the identity.
func TestPrecessIdentity(t *testing.T) {
	a, d := PrecessFK5(180, 30, 2000, 2000)
	if math.Abs(a-180) > 1e-12 || math.Abs(d-30) > 1e-12 {
		t.Errorf("self-precession: (%v, %v)", a, d)
	}
}
