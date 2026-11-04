package transform

import (
	"math"
	"testing"
)

// TestCSCUnderflowGuardPremises checks the numeric premises of the proof
// in cscXfYf that its float32 underflow guards (chi2, psi2, |chi psi| >
// 1e-16) decide nothing, so the mutants that widen or remove them are
// equivalent (the proof is in the comment there). With a = the guarded
// variable's square (or chi2 psi2) at most T = 1e-8 and b the other one,
// over every float32 b in [0, 1] (stepped by 2^-20, plus the end points):
//
//  1. a dropped term is at most |c20|, |c02|, |c11| times T^2, below
//     1.8e-17, and 2^25 times that is below 6.1e-10;
//  2. c00 + c10 a + c20 a^2 >= 0.0439 for a in [0, 1] (the partial sum in
//     front of a 6th-position guarded term);
//  3. wherever c00 + c10 a + c01 b (a <= 1) is within 6.1e-10 of zero,
//     b is in [0.50, 0.80] and c02 b^2 >= 2^-6 (half a float32 ULP
//     there is 2^-30 or more, i.e. >= 9.3e-10 > 6.1e-10).
func TestCSCUnderflowGuardPremises(t *testing.T) {
	const T = 1e-8
	maxc := math.Max(math.Max(math.Abs(float64(cscC20)), math.Abs(float64(cscC02))), math.Abs(float64(cscC11)))
	drop := maxc * T * T
	if !(drop < 1.8e-17 && math.Ldexp(drop, 25) < 6.1e-10) {
		t.Errorf("premise 1: dropped term bound %.3g", drop)
	}
	minF := math.Inf(1)
	for i := 0; i <= 1<<20; i++ {
		a := float32(math.Ldexp(float64(i), -20))
		f := cscC00 + float32(cscC10*a) + float32(float32(cscC20*a)*a)
		minF = math.Min(minF, float64(f))
	}
	if !(minF >= 0.0439) {
		t.Errorf("premise 2: min of c00 + c10 a + c20 a^2 over [0, 1] = %.6g, the proof needs >= 0.0439", minF)
	}
	// Premise 3: the b that make c00 + c10 a + c01 b vanish lie in
	// [-(c00)/c01, -(c00 + c10)/c01]; check those ends and c02 b^2.
	lo := -float64(cscC00) / float64(cscC01)
	hi := -(float64(cscC00) + float64(cscC10)) / float64(cscC01)
	// widen by the 6.1e-10 slack and float32 rounding of the products
	slack := 6.1e-10/math.Abs(float64(cscC01)) + 1e-6
	lo, hi = lo-slack, hi+slack
	if !(lo >= 0.50 && hi <= 0.80) {
		t.Errorf("premise 3: cancelling b in [%.6g, %.6g], the proof needs [0.50, 0.80]", lo, hi)
	}
	for _, b := range []float64{lo, hi} {
		c := float64(cscC02) * b * b
		if !(c >= 1.0/64) {
			t.Errorf("premise 3: c02 b^2 = %.6g at b = %.6g, the proof needs >= 2^-6", c, b)
		}
	}
	if h := math.Ldexp(1, -30); !(h > 6.1e-10) {
		t.Errorf("premise 3: half a float32 ULP at 2^-6 is %.3g", h)
	}
}
