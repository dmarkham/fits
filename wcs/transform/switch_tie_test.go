package transform

import (
	"math"
	"testing"
)

// TestSmallAngleSwitchTies proves that three of the small-angle switches
// can never be tied, so writing their comparison as <= instead of <
// (as wcslib 8.10 writes them) would change no result: those mutants
// survive every golden for this reason, not for want of a probe.
// TestPoleSmallAngleBand pins where each switch is (to 0.1%), and the
// inverse SIN/SZP/AIR switches, which CAN be tied, are probed at the tie
// itself (wref.c tie_inv).
//
//   - SIN forward (zenithal.go): t = (90 - |theta|) d2r < 1e-5;
//   - AIR forward (zenithal_extra.go): xi = d2r (90 - theta) / 2 < 1e-4;
//   - XPH inverse (healpix.go): sigma = (90 - |eta|) / 45 < xphW3 = 1e-4.
//
// In each, the angle next to the switch is a double in [64, 128), a
// multiple of 2^-46, and 90 minus it is exact (Sterbenz), so the value
// compared is fl(k 2^-46 c) (c = d2r, d2r / 2, 1/45) for an integer k.
// It can equal the threshold T only if k 2^-46 c is within half an ULP of
// T, a window 1e5 times narrower than the spacing of k 2^-46 c, so only
// the k nearest T / (2^-46 c) can tie. The test evaluates the port's own
// expressions over that k and 1000 either side, and requires that none
// ties and that the value crosses T (the switch is reached from both
// sides).
func TestSmallAngleSwitchTies(t *testing.T) {
	ulp := math.Ldexp(1, -46)
	air, err := newAIR(nil, 2)
	if err != nil {
		t.Fatal(err)
	}
	for _, sw := range []struct {
		name   string
		thresh float64
		scale  float64 // value / (90 - angle), for the starting k
		value  func(angle float64) float64
	}{
		{"SIN forward t", 1.0e-5, d2r, func(theta float64) float64 { return float64((90.0 - math.Abs(theta)) * d2r) }},
		{"AIR forward xi", air.w4, d2r / 2, func(theta float64) float64 { return d2r * (90.0 - theta) / 2.0 }},
		{"XPH inverse sigma", xphW3, 1.0 / 45, func(abseta float64) float64 { return (90.0 - abseta) / 45.0 }},
	} {
		k0 := math.Round(sw.thresh / sw.scale / ulp)
		var below, above bool
		for dk := -1000.0; dk <= 1000; dk++ {
			angle := 90.0 - (k0+dk)*ulp
			if angle < 64 || angle >= 128 || angle-math.Nextafter(angle, 0) != ulp {
				t.Fatalf("%s: angle %v is not in the 2^-46 binade", sw.name, angle)
			}
			v := sw.value(angle)
			if v == sw.thresh {
				t.Errorf("%s: angle %.17g gives exactly the threshold %g: the < / <= choice matters there and needs a wcslib probe", sw.name, angle, sw.thresh)
			}
			below = below || v < sw.thresh
			above = above || v > sw.thresh
		}
		if !below || !above {
			t.Errorf("%s: the value does not cross %g over k0 +- 1000 (below %v, above %v)", sw.name, sw.thresh, below, above)
		}
	}
}
