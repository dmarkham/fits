package transform

import (
	"bytes"
	"fmt"
	"math"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// Policy tests: behaviour that is our choice (documented on New) rather
// than something wcslib computes. Where wcslib has an opinion it is cited.

// TestSIPInverseNonConvergenceNotOK (PR #3 review): when Newton iteration
// cannot converge, Inverse must say so. F(u) = u - 1e-3 u^2 peaks at
// F(500) = 250, so u' = 300 has no real inverse; the old code returned
// ok = true after its last iteration with whatever u it had.
func TestSIPInverseNonConvergenceNotOK(t *testing.T) {
	sip := &SIP{
		A: &SIPPoly{Order: 2, Coeffs: [][]float64{{0, 0, 0}, {0, 0, 0}, {-1e-3, 0, 0}}},
		B: &SIPPoly{Order: 2, Coeffs: [][]float64{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}},
	}
	if u, v, ok := sip.Inverse(300, 10); ok {
		fu, fv := sip.Forward(u, v)
		t.Fatalf("Inverse(300, 10) reported ok with (%g, %g) -> Forward (%g, %g)", u, v, fu, fv)
	}
	// Targets just past the fold maximum (250 + 1e-4 .. 250 + 0.2 px) have
	// no inverse either. Newton creeps towards the fold there with steps
	// that stop shrinking well above round-off; accepting such a plateau
	// (a stall tolerance of 1e-2 instead of 1e-11) reported ok on 182 of
	// these 2000 targets.
	for k := 0; k < 2000; k++ {
		tu := 250 + 1e-4*math.Pow(2000, float64(k)/1999) // 250.0001 .. 250.2, log-spaced
		if u, v, ok := sip.Inverse(tu, 10); ok {
			fu, fv := sip.Forward(u, v)
			t.Fatalf("Inverse(%.10g, 10) beyond the fold reported ok with (%g, %g) -> Forward (%g, %g)", tu, u, v, fu, fv)
		}
	}
	// Same polynomial, a point that does have an inverse, converges.
	u, v, ok := sip.Inverse(200, 10)
	if !ok {
		t.Fatal("Inverse(200, 10) failed")
	}
	if fu, fv := sip.Forward(u, v); math.Abs(fu-200) > 1e-12 || math.Abs(fv-10) > 1e-12 {
		t.Fatalf("Forward(Inverse(200, 10)) = (%g, %g)", fu, fv)
	}

	// End to end: SkyToPixel must turn the non-convergence into an error,
	// not a pixel. Same fold (u' = u - 1e-3 u^2, max 250 px) behind a TAN
	// header at CRVAL (0, 0) with 1e-3 deg/px and RA increasing with x.
	w, err := wcs.Parse(newHeader("NAXIS", int64(2),
		"CTYPE1", "RA---TAN-SIP", "CTYPE2", "DEC--TAN-SIP", "CRVAL1", 0.0, "CRVAL2", 0.0,
		"CRPIX1", 1.0, "CRPIX2", 1.0, "CD1_1", 1e-3, "CD2_2", 1e-3,
		"A_ORDER", int64(2), "A_2_0", -1e-3, "B_ORDER", int64(0)))
	if err != nil {
		t.Fatal(err)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	// Plane x = 0.3 deg needs u' = 300 px > 250: no inverse.
	ra := radToDeg(math.Atan(degToRad(0.3)))
	if x, y, err := tr.SkyToPixel(ra, 0); err == nil {
		t.Fatalf("SkyToPixel through a non-invertible SIP returned (%g, %g)", x, y)
	}
	// Plane x = 0.2 deg (u' = 200) is invertible and must round-trip.
	ra = radToDeg(math.Atan(degToRad(0.2)))
	x, y, err := tr.SkyToPixel(ra, 0)
	if err != nil {
		t.Fatalf("SkyToPixel at u' = 200: %v", err)
	}
	if a, d, err := tr.PixelToSky(x, y); err != nil || math.Abs(a-ra) > 1e-12 || math.Abs(d) > 1e-12 {
		t.Fatalf("PixelToSky(SkyToPixel) = (%v, %v, %v), want (%v, 0)", a, d, err, ra)
	}
}

// TestSIPInverseConvergesAtLargeOffsets (PR #3 review note): an absolute
// 1e-13 px test can never pass at a few thousand px from CRPIX (one ULP of
// 2640 is 4.5e-13), so a correct relative test is needed for real
// detectors to report ok. Every point of a 5280 x 5280 UH88-like detector
// must invert to round-off.
func TestSIPInverseConvergesAtLargeOffsets(t *testing.T) {
	sip := &SIP{
		A: &SIPPoly{Order: 2, Coeffs: [][]float64{{0, 0, -4.51e-8}, {0, -1.02e-7, 0}, {-5.44e-8, 0, 0}}},
		B: &SIPPoly{Order: 2, Coeffs: [][]float64{{0, 0, -7.14e-8}, {0, -2.08e-7, 0}, {-1.52e-7, 0, 0}}},
	}
	for _, u := range []float64{-2640, -1000, 0, 1500, 2640} {
		for _, v := range []float64{-2640, -7, 2000, 2640} {
			fu, fv := sip.Forward(u, v)
			gu, gv, ok := sip.Inverse(fu, fv)
			if !ok {
				t.Fatalf("Inverse at (%g, %g) did not converge", u, v)
			}
			if d := math.Max(math.Abs(gu-u), math.Abs(gv-v)); d > 1e-10 {
				t.Errorf("(%g, %g): round trip off by %.3g px", u, v, d)
			}
		}
	}
}

// TestTPVInverseNonConvergenceNotOK (PR #3 review): xi' = xi - xi^2 peaks
// at 0.25, so xi' = 0.3 has no inverse; Inverse must not report ok, and
// SkyToPixel must return an error instead of a pixel.
func TestTPVInverseNonConvergenceNotOK(t *testing.T) {
	tpv := &TPV{Lon: 1, Lat: 2}
	tpv.Axis1[1], tpv.Axis1[4] = 1, -1
	tpv.Axis2[1] = 1
	if xi, eta, ok := tpv.Inverse(0.3, 0.1); ok {
		t.Fatalf("Inverse(0.3, 0.1) reported ok with (%g, %g)", xi, eta)
	}
	// Just past the fold maximum (0.25 + 1e-7 .. 0.25 + 2e-4 deg): no
	// inverse, so no ok from a slow creep towards the fold.
	for k := 0; k < 2000; k++ {
		tx := 0.25 + 1e-7*math.Pow(2000, float64(k)/1999)
		if xi, eta, ok := tpv.Inverse(tx, 0.1); ok {
			t.Fatalf("Inverse(%.10g, 0.1) beyond the fold reported ok with (%g, %g)", tx, xi, eta)
		}
	}
	xi, eta, ok := tpv.Inverse(0.2, 0.1)
	if !ok {
		t.Fatal("Inverse(0.2, 0.1) failed")
	}
	if fx, fy := tpv.Forward(xi, eta); math.Abs(fx-0.2) > 1e-15 || math.Abs(fy-0.1) > 1e-15 {
		t.Fatalf("Forward(Inverse) = (%g, %g)", fx, fy)
	}

	w, err := wcs.Parse(newHeader("NAXIS", int64(2),
		"CTYPE1", "RA---TPV", "CTYPE2", "DEC--TPV", "CRVAL1", 0.0, "CRVAL2", 0.0,
		"CD1_1", 1e-3, "CD2_2", 1e-3, "PV1_1", 1.0, "PV1_4", -1.0, "PV2_1", 1.0))
	if err != nil {
		t.Fatal(err)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	// TAN position with xi' = 0.3 deg, eta' = 0.
	if x, y, err := tr.SkyToPixel(0.3, 0); err == nil {
		t.Fatalf("SkyToPixel through a non-invertible TPV returned (%g, %g)", x, y)
	}
}

// TestDistortionCodeMustMatchOnBothAxes (PR #3 review): the distortion
// code is checked on BOTH celestial axes. wcslib 8.10 accepts a header
// with "-SIP" on one axis only (and applies the SIP keywords); we refuse
// it as ambiguous, as decided for this repository.
func TestDistortionCodeMustMatchOnBothAxes(t *testing.T) {
	for _, ct := range [][2]string{
		{"RA---TAN-SIP", "DEC--TAN"},
		{"RA---TAN", "DEC--TAN-SIP"},
		{"RA---TAN-SIP", "DEC--TAN-TPV"},
	} {
		w, err := wcs.Parse(newHeader("NAXIS", int64(2), "CTYPE1", ct[0], "CTYPE2", ct[1],
			"CDELT1", -1e-3, "CDELT2", 1e-3,
			"A_ORDER", int64(2), "A_2_0", 1e-5, "B_ORDER", int64(2), "B_0_2", 1e-5))
		if err != nil {
			t.Fatal(err)
		}
		if _, err := New(w); err == nil || !strings.Contains(err.Error(), "inconsistent distortion") {
			t.Errorf("%v: New error = %v, want inconsistent distortion codes", ct, err)
		}
	}
}

// TestNewRefusesUnimplementedDistortions: the unmodified WFC3/UVIS FLT
// header (real_hst_wfc3_uvis_sip4 is the same header with its CPDIS /
// D2IM lookup-table records removed) carries lookup-table distortions we
// do not implement; applying only its SIP part would be silently ~0.1 px
// off, so New refuses. Likewise SIP referenced from an alternate
// description (the ACS FLT header's 'O' WCS): SIP keywords have no alt
// form; wcslib silently applies no distortion there, we refuse.
func TestNewRefusesUnimplementedDistortions(t *testing.T) {
	g := loadHdrGolden(t, "real_hst_wfc3_uvis_sip4")
	cards := append([]string{}, g.Cards...)
	cards = append(cards,
		"CPDIS1  = 'LOOKUP  '           / Prior distortion function type",
		"CPDIS2  = 'LOOKUP  '           / Prior distortion function type",
		"D2IMDIS1= 'LOOKUP  '           / Detector to image correction type",
		"D2IMDIS2= 'LOOKUP  '           / Detector to image correction type")
	w, err := hdrParse(cards, "")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := New(w); err == nil || !strings.Contains(err.Error(), "CPDIS1") {
		t.Errorf("WFC3 FLT with CPDIS/D2IM: New error = %v", err)
	}

	acs := loadHdrGolden(t, "real_hst_acs_wfc_sip4")
	w, err = hdrParse(acs.Cards, "O")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := New(w); err == nil || !strings.Contains(err.Error(), "SIP") {
		t.Errorf("ACS alt WCS 'O' with -SIP: New error = %v", err)
	}
}

// TestNewRefusesCoupledMatrix: a PC matrix that feeds a non-celestial
// pixel axis into RA/DEC cannot be evaluated from (p1, p2) alone.
func TestNewRefusesCoupledMatrix(t *testing.T) {
	w, err := wcs.Parse(newHeader("NAXIS", int64(3),
		"CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN", "CTYPE3", "FREQ",
		"CDELT1", -1e-3, "CDELT2", 1e-3, "CDELT3", 1e6, "PC1_3", 0.01))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := New(w); err == nil || !strings.Contains(err.Error(), "couples pixel axis 3") {
		t.Errorf("New error = %v, want coupling error", err)
	}
	// Coupling only INTO the spectral axis is harmless.
	w, err = wcs.Parse(newHeader("NAXIS", int64(3),
		"CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN", "CTYPE3", "FREQ",
		"CDELT1", -1e-3, "CDELT2", 1e-3, "CDELT3", 1e6, "PC3_1", 0.01))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := New(w); err != nil {
		t.Errorf("PC3_1 only: %v", err)
	}
}

// TestNewRefusesSIPOffAxes12: SIP is defined on pixel axes 1 and 2; a
// "-SIP" pair on axes 2,3 of a cube is refused, not applied to the wrong
// axes.
func TestNewRefusesSIPOffAxes12(t *testing.T) {
	w, err := wcs.Parse(newHeader("NAXIS", int64(3),
		"CTYPE1", "FREQ", "CTYPE2", "RA---TAN-SIP", "CTYPE3", "DEC--TAN-SIP",
		"CDELT2", -1e-3, "CDELT3", 1e-3,
		"A_ORDER", int64(2), "A_2_0", 1e-5, "B_ORDER", int64(2), "B_0_2", 1e-5))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := New(w); err == nil || !strings.Contains(err.Error(), "SIP on celestial pixel axes 2,3") {
		t.Errorf("New error = %v", err)
	}
}

// TestNewRefusesTPVWithNonDegreeUnits: TPV coefficients are in degrees
// (registry); a TPV header with CUNIT arcsec is ambiguous and refused.
func TestNewRefusesTPVWithNonDegreeUnits(t *testing.T) {
	w, err := wcs.Parse(newHeader("NAXIS", int64(2),
		"CTYPE1", "RA---TPV", "CTYPE2", "DEC--TPV", "CUNIT1", "arcsec", "CUNIT2", "arcsec",
		"CD1_1", -1.0, "CD2_2", 1.0, "PV1_1", 1.0, "PV2_1", 1.0))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := New(w); err == nil || !strings.Contains(err.Error(), "TPV") {
		t.Errorf("New error = %v", err)
	}
}

// TestNewtonStallAtRoundOffAccepted pins the other acceptance branch of
// newton2: when the evaluation itself is only good to a quantum q (here
// F(x) = x rounded to q = 1e-12, an exact model of round-off-limited
// evaluation), the Newton step bounces between 0.4q and 0.6q and never
// meets the relative test (1e-14 * scale); once the step stops shrinking
// at that level (below newtonStallTol * scale) the result is accepted, and
// it is within q of the true root. Without the stall branch this returns
// ok = false after newtonMaxIter iterations.
func TestNewtonStallAtRoundOffAccepted(t *testing.T) {
	const q = 1e-12
	f := func(x, y float64) (fx, fy, j11, j12, j21, j22 float64) {
		return math.Round(x/q) * q, y, 1, 0, 0, 1
	}
	tx := 1 + 0.4*q
	x, y, ok := newton2(f, tx, 0, 0, 0, 1)
	if !ok {
		t.Fatalf("newton2 on a round-off plateau: ok = false (x = %.17g)", x)
	}
	if math.Abs(x-tx) > q || y != 0 {
		t.Fatalf("newton2 accepted (%.17g, %g), %g from the root (quantum %g)", x, y, x-tx, q)
	}
	// A plateau far above round-off (q = 1e-6: steps ~ 5e-7 * scale) is
	// not accepted.
	g := func(x, y float64) (fx, fy, j11, j12, j21, j22 float64) {
		return math.Round(x/1e-6) * 1e-6, y, 1, 0, 0, 1
	}
	if x, _, ok := newton2(g, 1+0.4e-6, 0, 0, 0, 1); ok {
		t.Fatalf("newton2 accepted a 1e-6 plateau at %.17g", x)
	}
}

// TestNewtonRefusesNonFiniteJacobian: newton2's two loud-failure guards.
// Each input is built so that, without the guard, a point that is not a
// root would come back with ok = true (checked by exact arithmetic, not by
// our own transform):
//
//   - a determinant that OVERFLOWS although every Jacobian entry is finite
//     (F = (1e200 x, 1e200 y): det = 1e400 = +Inf, target (1, 1)): each
//     Newton step is (finite)/Inf = 0, so without the non-finite-
//     determinant guard the zero step passes the convergence test and the
//     start (0, 0), where F = (0, 0) != (1, 1), would be returned as the
//     root. (A det of exactly 0, NaN, or an infinite ENTRY always makes the
//     step non-finite and is also caught by the iterate guard.)
//   - a determinant so small that the step overflows (F = 1e-300 x,
//     target 1e10: the root 1e310 is not a float64): the iterate becomes
//     +Inf, and without the non-finite-iterate guard Inf <= 1e-14 * Inf
//     passes the relative test and +Inf would come back as the root.
func TestNewtonRefusesNonFiniteJacobian(t *testing.T) {
	big := func(x, y float64) (fx, fy, j11, j12, j21, j22 float64) {
		return 1e200 * x, 1e200 * y, 1e200, 0, 0, 1e200
	}
	if x, y, ok := newton2(big, 1, 1, 0, 0, 1); ok {
		if fx, fy, _, _, _, _ := big(x, y); math.Abs(fx-1) > 1e-12 || math.Abs(fy-1) > 1e-12 {
			t.Errorf("overflowing determinant: newton2 returned ok at (%g, %g), where F = (%g, %g) != (1, 1)", x, y, fx, fy)
		}
	}
	tiny := func(x, y float64) (fx, fy, j11, j12, j21, j22 float64) {
		return 1e-300 * x, y, 1e-300, 0, 0, 1
	}
	if x, y, ok := newton2(tiny, 1e10, 0, 0, 0, 1); ok {
		t.Errorf("overflowing step: newton2 returned ok at (%g, %g); the root 1e310 is not representable", x, y)
	}
}

// TestNewTNXRefusesSIP: a TNX header that also carries SIP keywords (with
// or without -SIP) has two distortion solutions. wcslib 8.10 applies SIP
// with TAN and ignores the WAT correction (hdrcases_stricter/tnx_sip_*,
// accepted there by TestHdrStricterThanWcslib); NewTNX used to apply the
// TNX correction and silently ignore the SIP set it had parsed. Both
// NewTNX and New must refuse the header.
func TestNewTNXRefusesSIP(t *testing.T) {
	for _, name := range []string{"tnx_sip_nosuffix", "tnx_sip_suffix"} {
		cards, _, _ := readCaseCards(t, filepath.Join(hdrStricterDir, name+".case"))
		var buf bytes.Buffer
		for _, c := range cards {
			buf.WriteString(fmt.Sprintf("%-80s", c))
		}
		buf.WriteString(fmt.Sprintf("%-80s", "END"))
		for buf.Len()%2880 != 0 {
			buf.WriteByte(' ')
		}
		parsed, _, err := header.ParseCards(buf.Bytes())
		if err != nil {
			t.Fatal(err)
		}
		h := header.FromCards(parsed)
		if tr, err := NewTNX(h); err == nil {
			a, d, _ := tr.PixelToSky(0, 0)
			t.Errorf("%s: NewTNX built a transform that ignores the SIP keywords (pixel (0,0) -> %.10f, %.10f)", name, a, d)
		} else if !strings.Contains(err.Error(), "TNX header carries SIP keywords") {
			t.Errorf("%s: NewTNX refused, but with %v", name, err)
		}
		w, err := wcs.Parse(h)
		if err == nil {
			_, err = New(w)
		}
		if err == nil || !strings.Contains(err.Error(), "TNX header carries SIP keywords") {
			t.Errorf("%s: Parse/New: %v, want the TNX + SIP refusal", name, err)
		}
	}
}
