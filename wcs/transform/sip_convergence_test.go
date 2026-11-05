package transform

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// The two tests below are adapted from Daryl Manning's PR #3 review-fix
// commit da8e66d.

// TestSIPInverseConvergesFarFromCRPIX: the Newton inverse converges for a
// cubic SIP solution at offsets of ~1e4 px from CRPIX, so the convergence
// test is not tighter than the arithmetic allows there.
func TestSIPInverseConvergesFarFromCRPIX(t *testing.T) {
	sip := &SIP{
		A: &SIPPoly{Order: 3, Coeffs: [][]float64{{0, 0, 2e-7, 0}, {0, -1e-7, 0, 0}, {-5e-8, 0, 0, 0}, {3e-11, 0, 0, 0}}},
		B: &SIPPoly{Order: 3, Coeffs: [][]float64{{0, 0, -7e-8, 1e-11}, {0, -2e-7, 0, 0}, {-1.5e-7, 0, 0, 0}, {0, 0, 0, 0}}},
	}
	for _, p := range [][2]float64{{8000, -8000}, {-9500, 9100}, {12000, 300}} {
		up, vp := sip.Forward(p[0], p[1])
		u, v, ok := sip.Inverse(up, vp)
		if !ok {
			t.Fatalf("Inverse did not converge at %v", p)
		}
		if math.Abs(u-p[0]) > 1e-6 || math.Abs(v-p[1]) > 1e-6 {
			t.Errorf("round trip %v -> (%v, %v)", p, u, v)
		}
	}
}

// TestSIPSkyToPixelErrorsWithoutInverse: Copilot's PR #3 example, end to
// end. With A_2_0 = 1e-4 the undistorted offset u' = -10000 px has no SIP
// preimage (u + 1e-4 u^2 >= -2500 for every u), so SkyToPixel must return
// an error rather than the last Newton iterate. The sky point is taken
// from the same header without any SIP keywords (plain TAN).
func TestSIPSkyToPixelErrorsWithoutInverse(t *testing.T) {
	base := []any{"NAXIS", int64(2),
		"CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN",
		"CRVAL1", 150.0, "CRVAL2", 2.0, "CRPIX1", 1.0, "CRPIX2", 1.0,
		"CDELT1", -1e-5, "CDELT2", 1e-5}
	build := func(pairs []any) *Transform {
		w, err := wcs.Parse(newHeader(pairs...))
		if err != nil {
			t.Fatal(err)
		}
		tr, err := New(w)
		if err != nil {
			t.Fatal(err)
		}
		return tr
	}
	// PixelToSky takes 0-based pixels, so CRPIX1 = 1 is pixel 0 and pixel
	// -10000 is offset -10000.
	a, d, err := build(base).PixelToSky(-10000, 0)
	if err != nil {
		t.Fatal(err)
	}
	sipHdr := append([]any{}, base...)
	sipHdr[3], sipHdr[5] = "RA---TAN-SIP", "DEC--TAN-SIP"
	sipHdr = append(sipHdr, "A_ORDER", int64(2), "A_2_0", 1e-4, "B_ORDER", int64(2))
	if x, y, err := build(sipHdr).SkyToPixel(a, d); err == nil {
		t.Errorf("SkyToPixel(%v, %v) = (%v, %v), want an error (no SIP inverse exists)", a, d, x, y)
	}
}
