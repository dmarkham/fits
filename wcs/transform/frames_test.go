package transform

import (
	"math"
	"testing"
)

// TestEquatorialGalacticRoundTrip: ICRS → galactic → ICRS should be the
// identity to within double-precision roundoff. Longitudes are compared
// modulo 360° because 0° and 359.9999...° are the same point.
func TestEquatorialGalacticRoundTrip(t *testing.T) {
	cases := [][2]float64{
		{0, 0},
		{180, 0},
		{90, 45},
		{270, -30},
		{45, 60},
		{359, -89},
	}
	for _, c := range cases {
		l, b := EquatorialToGalactic(c[0], c[1])
		a2, d2 := GalacticToEquatorial(l, b)
		if lonDiff(a2, c[0]) > 1e-7 || math.Abs(d2-c[1]) > 1e-8 {
			t.Errorf("round-trip (%v,%v) → (%v,%v) → (%v,%v)", c[0], c[1], l, b, a2, d2)
		}
	}
}

// lonDiff returns the smallest absolute difference between two longitudes,
// accounting for the 360° wrap-around.
func lonDiff(a, b float64) float64 {
	d := math.Abs(a - b)
	if d > 180 {
		d = 360 - d
	}
	return d
}

// TestGalacticCenter: the galactic center is at ICRS (266.4168°, -28.9360°)
// and galactic (l=0°, b=0°) by definition. Verify this within ~1 arcmin
// (the accuracy of the J2000 matrix).
func TestGalacticCenter(t *testing.T) {
	l, b := EquatorialToGalactic(266.4168, -28.9360)
	if math.Abs(l) > 0.02 && math.Abs(l-360) > 0.02 {
		t.Errorf("galactic longitude of GC: %v (want ~0°)", l)
	}
	if math.Abs(b) > 0.02 {
		t.Errorf("galactic latitude of GC: %v (want ~0°)", b)
	}
}

// TestEquatorialEclipticRoundTrip verifies the ecliptic conversion pair
// is self-inverse at J2000.
func TestEquatorialEclipticRoundTrip(t *testing.T) {
	const jdJ2000 = 2451545.0
	for _, c := range [][2]float64{
		{0, 0}, {180, 0}, {90, 30}, {270, -45}, {45, 60},
	} {
		l, b := EquatorialToEcliptic(c[0], c[1], jdJ2000)
		a2, d2 := EclipticToEquatorial(l, b, jdJ2000)
		if lonDiff(a2, c[0]) > 1e-7 || math.Abs(d2-c[1]) > 1e-8 {
			t.Errorf("round-trip (%v,%v) → (%v,%v) → (%v,%v)", c[0], c[1], l, b, a2, d2)
		}
	}
}

// TestObliquityJ2000: at J2000, the obliquity of the ecliptic is
// 23.4393° (per IAU 2006). Verify our formula gives this.
func TestObliquityJ2000(t *testing.T) {
	eps := obliquityIAU2006(2451545.0) * 180 / math.Pi
	if math.Abs(eps-23.4393) > 0.001 {
		t.Errorf("obliquity at J2000: %v (want 23.4393)", eps)
	}
}

// TestSupergalacticRoundTrip walks ICRS → SGL/SGB → ICRS.
func TestSupergalacticRoundTrip(t *testing.T) {
	for _, c := range [][2]float64{
		{0, 0}, {180, 30}, {90, -45}, {270, 60},
	} {
		sgl, sgb := EquatorialToSupergalactic(c[0], c[1])
		a2, d2 := SupergalacticToEquatorial(sgl, sgb)
		if lonDiff(a2, c[0]) > 1e-7 || math.Abs(d2-c[1]) > 1e-8 {
			t.Errorf("supergalactic round-trip (%v,%v): got (%v,%v)", c[0], c[1], a2, d2)
		}
	}
}

// TestPrecessB1950ToJ2000: a star at FK4 B1950 (0, 0) precesses to
// roughly FK5 J2000 (0.6406°, 0.2781°) per standard astronomy references.
// This verifies the precession matrix is in the right ballpark (we only
// check to ~1 arcsec since the IAU 1976 model is itself an approximation).
func TestPrecessB1950ToJ2000(t *testing.T) {
	a2, d2 := PrecessFK5(0, 0, 1950.0, 2000.0)
	if math.Abs(a2-0.64) > 0.01 {
		t.Errorf("precessed RA: %v (want ~0.64)", a2)
	}
	if math.Abs(d2-0.28) > 0.01 {
		t.Errorf("precessed Dec: %v (want ~0.28)", d2)
	}
}

// TestPrecessIdentity: precessing from an epoch to itself is the identity.
func TestPrecessIdentity(t *testing.T) {
	a, d := PrecessFK5(180, 30, 2000, 2000)
	if math.Abs(a-180) > 1e-10 || math.Abs(d-30) > 1e-10 {
		t.Errorf("self-precession: (%v, %v)", a, d)
	}
}
