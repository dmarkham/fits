package transform

import (
	"math"
	"testing"
)

// TestPCOSmallAngleBand pins the two small-angle switches of PCO
// (polyconic.go: the inverse uses theta = y / (1 + w3 x^2) for |y| < 1e-4
// deg, the forward uses x = phi cos(theta), y = (1 + w3 x^2) theta for
// |theta| < 1e-4 deg, both as wcslib 8.10 pcox2s / pcos2x do).
//
// Next to the equator the approximation and the exact branch differ by
// only ~1e-13 deg (measured with wcslib: 9.9e-14 at |y| = 1.1e-4, 1.1e-13
// at 9.9e-4, 5.7e-13 at 3e-3, 1.7e-11 at 9.9e-3 deg), below the 1e-11 deg
// absolute budget of TestWrefProjections, so that test alone cannot tell
// a port that moves a switch from a faithful one: the verifier widened
// the inverse switch to 1e-3 deg and every test still passed.
//
// The latitude-like coordinate (inverse theta, forward y) of every
// accepted band probe (wref.c band_inv / band_fwd: |y| or |theta| between
// 5e-5 and 9e-2 deg, off the central meridian) is therefore held to 4 ULP
// of wcslib's value, a relative budget of 9e-16. That is reachable
// because on these rows the port runs exactly wcslib's operations (the
// same branch, the same regula-falsi steps) and its result is measured
// bit-identical to wcslib's at every such row; any moved switch puts
// the approximation where wcslib iterates, or the reverse, and shows up
// at 1e-10 relative or more. With probes at 9.95e-5, 9.99e-5, 1e-4,
// 1.005e-4 and 1.01e-4 deg each switch is pinned to 0.1% either side and
// its comparison (< against <=) at the threshold itself; the mutation
// sweep (docs/VALIDATION.md) kills every move of 0.1% or more. The longitude-like coordinate keeps the
// absolute budget (it is ~170 deg, where 1 ULP is 2.8e-14).
func TestPCOSmallAngleBand(t *testing.T) {
	g := loadWref(t, "pco_default")
	proj, err := wrefSelect(g)
	if err != nil {
		t.Fatal(err)
	}
	var n int
	for _, inverse := range []bool{false, true} {
		rows := g.Forward
		if inverse {
			rows = g.Inverse
		}
		for _, r := range rows {
			in0, in1 := *r[0], *r[1]
			if math.Abs(in1) < 5e-5 || math.Abs(in1) > 0.09 || in0 == 0 {
				continue
			}
			if *r[4] != 0 || r[2] == nil {
				t.Errorf("inverse=%v (%g, %g): wcslib rejected a band probe", inverse, in0, in1)
				continue
			}
			n++
			w0, w1 := *r[2], *r[3]
			o0, o1, ok := wrefEvalDeg(proj, inverse, in0, in1)
			if !ok {
				t.Errorf("inverse=%v (%g, %g): rejected, wcslib gives (%.17g, %.17g)", inverse, in0, in1, w0, w1)
				continue
			}
			ulp := math.Nextafter(math.Abs(w1), math.Inf(1)) - math.Abs(w1)
			if d := math.Abs(o1 - w1); d > 4*ulp {
				t.Errorf("inverse=%v (%g, %g): latitude-like output %.17g, wcslib %.17g (%.3g deg, %.3g ULP; budget 4 ULP)",
					inverse, in0, in1, o1, w1, d, d/ulp)
			}
			if d := math.Abs(o0 - w0); d > wrefTolDeg {
				t.Errorf("inverse=%v (%g, %g): longitude-like output %.17g, wcslib %.17g", inverse, in0, in1, o0, w0)
			}
		}
	}
	// 13 values x 2 signs x 4 per direction (104 + 104), plus the inverse
	// round trips of the forward band probes that land inside the band.
	if n < 208 {
		t.Errorf("%d band probes in pco_default, want at least 208 (wref.c band_inv/band_fwd)", n)
	}
	// theta = 0 exactly: pcos2x returns y = -0.0 (its own theta == 0
	// case), where the small-angle formula would give +0. The value is
	// the same, the sign bit is not; it is pinned bit for bit.
	var nz int
	for _, r := range g.Forward {
		if *r[1] != 0 || *r[4] != 0 || r[3] == nil {
			continue
		}
		nz++
		_, o1, ok := wrefEvalDeg(proj, false, *r[0], *r[1])
		if !ok || math.Float64bits(o1) != math.Float64bits(*r[3]) {
			t.Errorf("forward (%g, 0): y = %v (ok %v), wcslib %v (bits %#x vs %#x)", *r[0], o1, ok, *r[3], math.Float64bits(o1), math.Float64bits(*r[3]))
		}
	}
	if nz < 27 {
		t.Errorf("%d forward rows at theta = 0, want at least 27 (the equator of the regular grid)", nz)
	}
	t.Logf("%d band rows, latitude-like output within 4 ULP of wcslib", n)
}
