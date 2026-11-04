package transform

import (
	"math"
	"testing"
)

// TestPoleSmallAngleBand pins the small-angle switches next to the native
// pole of the zenithal projections and XPH against wcslib 8.10, to a few
// ULP of the latitude-like coordinate:
//
//   - SIN forward 90 - |theta| < 1e-5 rad and, with a slant (PV2_1, PV2_2
//     not both 0), the inverse x^2 + y^2 < 1e-10 rad^2 (zenithal.go);
//   - SZP inverse x^2 + y^2 < 1e-10 rad^2 (zenithal_extra.go);
//   - AIR forward xi < 1e-4 rad and inverse r < w2 1e-4 (newAIR's tol);
//   - XPH forward theta < 90 - 1e-4 sqrt(2/3) rad and inverse sigma < 1e-4
//     (healpix.go xphW3, xphW5).
//
// On either side of such a switch wcslib's two formulas agree to 1e-10 ..
// 5e-9 deg (the exact branch is ill-conditioned there: theta = asin(s)
// with s ~ 1 loses eps / cos(theta)), so a port that moves the switch
// differs from wcslib only inside a thin annulus around the pole and by
// little more than TestWrefProjections' 1e-11 deg budget: before these
// probes existed the verifier moved the SIN/SZP switch to 1e-11 and
// 0.5e-10, AIR's tol to 1.2e-4 and 2e-4 and XPH's to 1.2e-4 .. 2e-4, and
// every test passed. wref.c pole_fwd now puts forward points at
// theta = +-(90 - d) for d from 1e-5 to 0.1 deg (26% apart) and 0.05%,
// 0.1% and 1% either side of every threshold (theta = w[5] exactly for
// XPH; the 0.05% points were added in round 2, when SIN's forward switch
// moved exactly 0.1% up still survived); their
// wcslib images are the inverse round-trip rows. wref.c tie_inv adds
// inverse inputs ON the SIN/SZP and AIR inverse thresholds (the computed
// r^2 or r equal to it exactly), which pin the comparison itself (< as
// wcslib writes it, not <=); the forward ones cannot be tied
// (TestSmallAngleSwitchTies).
//
// On these rows the port runs exactly wcslib's operations (same branch,
// same arithmetic), so the latitude-like output (forward: x and y, both
// small here; inverse: theta) is held to 4 ULP of wcslib's value plus its
// own libm sensitivity (libm_sens_test.go: how far it moves when every
// sin/cos/atan2/log result moves by 4 ULP, which is what separates Go's
// math functions from glibc's), and the longitude to
// TestWrefProjections' absolute budget. Where the formula is
// well-conditioned that is a few ULP; next to the pole the exact branch
// magnifies the last bit of libm by up to 1e8 (AIR's log(cos xi), STG's
// 1 + sin(theta), SIN/SZP/XPH) and the allowance grows with it, counted
// in the log. A moved switch puts the other formula where wcslib uses
// one, and that still fails here: with the allowance in place this test
// fails for every one of SIN's forward switch moved by +-0.1% and +20%,
// SIN's and SZP's inverse r^2 < 1e-10 changed to 1e-11, 0.5e-10 or <=,
// AIR's tol at 0.8e-4, 1.2e-4 and 2e-4, and XPH's w3 at 0.8e-4, 1.2e-4
// and 2e-4 (checked by hand when the allowance was introduced; the other
// formula is off by more than the sensitivity, which only covers the
// last bits of libm).
func TestPoleSmallAngleBand(t *testing.T) {
	poleP := map[float64]bool{-150.3: true, -41.7: true, 22.9: true, 113.1: true}
	var total int
	for _, name := range wrefCaseManifest {
		g := loadWref(t, name)
		if g.Prjset.Theta0 == nil || *g.Prjset.Theta0 != 90 {
			continue
		}
		proj, err := wrefSelect(g)
		if err != nil {
			t.Fatal(err)
		}
		ulp := func(v float64) float64 {
			a := math.Abs(v)
			return math.Nextafter(a, math.Inf(1)) - a
		}
		band := map[[2]uint64]bool{}
		var n int
		var worst float64
		var cond condTally
		ulpDist := func(inverse bool) func(a0, a1, b0, b1 float64) float64 {
			return func(a0, a1, b0, b1 float64) float64 {
				if inverse {
					return math.Abs(b1-a1) / ulp(a1)
				}
				return math.Max(math.Abs(b0-a0)/ulp(a0), math.Abs(b1-a1)/ulp(a1))
			}
		}
		check := func(inverse bool, r [5]*float64) {
			in0, in1 := *r[0], *r[1]
			if *r[4] != 0 || r[2] == nil || r[3] == nil {
				if o0, o1, ok := wrefEvalDeg(proj, inverse, in0, in1); ok {
					t.Errorf("%s inverse=%v (%.17g, %.17g): wcslib rejects, we return (%.17g, %.17g)", name, inverse, in0, in1, o0, o1)
				}
				return
			}
			w0, w1 := *r[2], *r[3]
			o0, o1, ok := wrefEvalDeg(proj, inverse, in0, in1)
			if !ok {
				t.Errorf("%s inverse=%v (%.17g, %.17g): rejected, wcslib gives (%.17g, %.17g)", name, inverse, in0, in1, w0, w1)
				return
			}
			n++
			var u float64
			if inverse {
				u = math.Abs(o1-w1) / ulp(w1)
				if d := wrefDiff(true, w0, w1, o0, o1); d > wrefTolDeg {
					t.Errorf("%s inverse (%.17g, %.17g): got (%.17g, %.17g), wcslib (%.17g, %.17g), %.3g deg", name, in0, in1, o0, o1, w0, w1, d)
				}
			} else {
				u = math.Max(math.Abs(o0-w0)/ulp(w0), math.Abs(o1-w1)/ulp(w1))
			}
			worst = math.Max(worst, u)
			cond.rows++
			sens := func() float64 {
				return libmSensitivity(func() (float64, float64, bool) { return wrefEvalDeg(proj, inverse, in0, in1) }, ulpDist(inverse))
			}
			if !cond.valueOK(u, 4, sens) {
				t.Errorf("%s inverse=%v (%.17g, %.17g): got (%.17g, %.17g), wcslib (%.17g, %.17g): %.3g ULP, budget 4 + libm sensitivity %.3g",
					name, inverse, in0, in1, o0, o1, w0, w1, u, sens())
			}
		}
		for _, r := range g.Forward {
			if !poleP[*r[0]] || 90-math.Abs(*r[1]) > 0.1001 {
				continue
			}
			check(false, r)
			if *r[4] == 0 && r[2] != nil && r[3] != nil {
				band[[2]uint64{math.Float64bits(*r[2]), math.Float64bits(*r[3])}] = true
			}
		}
		// The inverse rows: the round trips of the band rows, and every
		// other inverse input within 0.2 deg of the pole's image (0, 0):
		// the tie probes (wref.c tie_inv: inputs on SIN's/SZP's
		// r^2 = 1e-10 and AIR's r = w2 1e-4 exactly, and the doubles
		// either side) and the boundary probes around the pole image.
		for _, r := range g.Inverse {
			if band[[2]uint64{math.Float64bits(*r[0]), math.Float64bits(*r[1])}] || math.Hypot(*r[0], *r[1]) <= 0.2 {
				check(true, r)
			}
		}
		if n < 300 {
			t.Errorf("%s: %d accepted pole band rows, want at least 300 (wref.c pole_fwd)", name, n)
		}
		total += n
		t.Logf("%-28s %4d pole band rows, worst %.3g ULP; %d needed their libm sensitivity (worst used %.2g of it)", name, n, worst, cond.conditioned, cond.worstRatio)
	}
	if total < 25*300 {
		t.Errorf("%d pole band rows in all, want at least %d", total, 25*300)
	}
}
