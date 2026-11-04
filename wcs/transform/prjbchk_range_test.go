package transform

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// TestPrjbchkRangeByConstruction proves that the final prjbchk of TSC,
// CSC, QSC and AIR (x2s) can never clamp or reject, so its tolerance
// decides nothing there and the mutants that widen it are equivalent
// (they survive every golden for this reason, not for want of a probe):
//
//   - phi is atan2d(m, l) (or 0): |math.Atan2| <= fl(pi) (for x < 0 Go
//     returns atan(y/x) +- fl(pi) with the atan term of the opposite
//     sign, and rounding is monotonic), so |atan2d| <= fl(fl(pi) r2d) =
//     180 exactly;
//   - theta is asind(n) (TSC, CSC, QSC): +-90 exactly for |n| >= 1, else
//     |fl(libAsin(n) r2d)| <= fl(fl(pi/2) r2d) = 90 exactly;
//   - AIR: theta = 90 - fl(2 xi) with xi = acosd(cos xi), cos xi in (0, 1]
//     (the regula falsi keeps it between two positive bracket ends), so
//     0 <= xi <= fl(fl(pi/2) r2d) = 90 and theta >= -90 exactly; or
//     xi = r w6 with 0 < r < w5 = w2 1e-4, w6 = r2d / w2, so
//     0 < xi < 0.0058 and theta in (89.98, 90).
//
// So the clamped values these four ever pass prjbchk are inside
// [-180, 180] x [-90, 90] already. (AZP's and XPH's prjbchk, which CAN
// see values a few 1e-14 .. 1e-7 deg out, are pinned by the wref cases
// azp_mu1_gamma45 and xph_default's seam probes, wref.c xph_seam_inv.)
func TestPrjbchkRangeByConstruction(t *testing.T) {
	if v := float64(math.Pi * r2d); v != 180 {
		t.Errorf("fl(fl(pi) r2d) = %.17g, the proof needs exactly 180", v)
	}
	if v := float64((math.Pi / 2) * r2d); v != 90 {
		t.Errorf("fl(fl(pi/2) r2d) = %.17g, the proof needs exactly 90", v)
	}
	// The extremes of the functions themselves.
	for _, c := range []struct {
		name string
		got  float64
		max  float64
	}{
		{"|atan2d(+tiny, -1)|", math.Abs(atan2d(0x1p-1000, -1)), 180},
		{"|atan2d(-tiny, -1)|", math.Abs(atan2d(-0x1p-1000, -1)), 180},
		{"|atan2d(-0, -1)|", math.Abs(atan2d(math.Copysign(0, -1), -1)), 180},
		{"|math.Atan2(+tiny, -1)|", math.Abs(math.Atan2(0x1p-1000, -1)), math.Pi},
		{"|math.Atan2(-tiny, -1)|", math.Abs(math.Atan2(-0x1p-1000, -1)), math.Pi},
		{"asind(1 - ulp)", asind(math.Nextafter(1, 0)), 90},
		{"-asind(-1 + ulp)", -asind(math.Nextafter(-1, 0)), 90},
		{"libAsin(1)", libAsin(1), math.Pi / 2},
		{"acosd(tiny)", acosd(0x1p-1000), 90},
		{"acosd(ulp)", acosd(0x1p-1074), 90},
	} {
		if !(c.got <= c.max) {
			t.Errorf("%s = %.17g > %.17g", c.name, c.got, c.max)
		}
	}
	// AIR's small-radius branch: xi = r w6 < w5 w6.
	for _, tb := range []float64{90, 45, 0, -45, -89} {
		p, err := newAIR(map[wcs.PVKey]float64{{Axis: 2, Index: 1}: tb}, 2)
		if err != nil {
			t.Fatal(err)
		}
		if xi := float64(p.w5 * p.w6); !(xi > 0 && xi < 0.0058) {
			t.Errorf("AIR theta_b = %g: w5 w6 = %.17g, the proof needs (0, 0.0058)", tb, xi)
		}
	}
}
