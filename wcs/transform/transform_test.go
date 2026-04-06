package transform

import (
	"errors"
	"math"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// newHeader is a tiny test helper: builds a *header.Header from key-value
// pairs so that tests can call wcs.Parse.
func newHeader(pairs ...any) *header.Header {
	h := header.New()
	for i := 0; i < len(pairs); i += 2 {
		k := pairs[i].(string)
		v := pairs[i+1]
		if err := h.Set(k, v, ""); err != nil {
			panic(err)
		}
	}
	return h
}

// makeTANHeader builds a parsed *wcs.Header for a TAN projection with
// straightforward defaults used by several tests.
func makeTANHeader(t *testing.T, crval1, crval2, crpix1, crpix2, cd11, cd12, cd21, cd22 float64) *wcs.Header {
	t.Helper()
	h := newHeader(
		"NAXIS", int64(2),
		"NAXIS1", int64(512),
		"NAXIS2", int64(512),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CUNIT1", "deg",
		"CUNIT2", "deg",
		"CRVAL1", crval1,
		"CRVAL2", crval2,
		"CRPIX1", crpix1,
		"CRPIX2", crpix2,
		"CD1_1", cd11,
		"CD1_2", cd12,
		"CD2_1", cd21,
		"CD2_2", cd22,
	)
	w, err := wcs.Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	return w
}

// TestTANReferencePointIdentity: the reference pixel must map exactly to
// the reference celestial coordinate. This is the "does the plumbing
// work at all" test.
func TestTANReferencePointIdentity(t *testing.T) {
	// CRPIX in FITS 1-based, Go 0-based = CRPIX - 1.
	w := makeTANHeader(t, 180.0, 30.0, 256.5, 256.5, -0.0002778, 0, 0, 0.0002778)
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	// Pixel (255.5, 255.5) in 0-based == CRPIX - 1.
	a, d, err := tr.PixelToSky(255.5, 255.5)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(a-180.0) > 1e-12 {
		t.Fatalf("reference RA: got %.15g, want 180", a)
	}
	if math.Abs(d-30.0) > 1e-12 {
		t.Fatalf("reference Dec: got %.15g, want 30", d)
	}
}

// TestTANSkyToPixelReference: SkyToPixel at (CRVAL1, CRVAL2) must return
// (CRPIX1 - 1, CRPIX2 - 1).
func TestTANSkyToPixelReference(t *testing.T) {
	w := makeTANHeader(t, 180.0, 30.0, 256.5, 256.5, -0.0002778, 0, 0, 0.0002778)
	tr, _ := New(w)
	p1, p2, err := tr.SkyToPixel(180.0, 30.0)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(p1-255.5) > 1e-9 || math.Abs(p2-255.5) > 1e-9 {
		t.Fatalf("got (%.15g, %.15g) want (255.5, 255.5)", p1, p2)
	}
}

// TestTANRoundTrip: pixel → sky → pixel must be an identity to floating-
// point precision across a grid of sample points inside the image.
func TestTANRoundTrip(t *testing.T) {
	w := makeTANHeader(t, 180.0, 30.0, 256.5, 256.5, -0.0002778, 0, 0, 0.0002778)
	tr, _ := New(w)
	for _, p := range []struct{ x, y float64 }{
		{0, 0}, {100, 100}, {255.5, 255.5}, {400, 100}, {511, 511},
		{256, 256}, {1.5, 512.25},
	} {
		a, d, err := tr.PixelToSky(p.x, p.y)
		if err != nil {
			t.Fatalf("PixelToSky(%v,%v): %v", p.x, p.y, err)
		}
		p1, p2, err := tr.SkyToPixel(a, d)
		if err != nil {
			t.Fatalf("SkyToPixel(%v,%v): %v", a, d, err)
		}
		if math.Abs(p1-p.x) > 1e-9 || math.Abs(p2-p.y) > 1e-9 {
			t.Errorf("roundtrip %v → (%.15g,%.15g) → (%.15g,%.15g)", p, a, d, p1, p2)
		}
	}
}

// TestTANAtEquator: on a TAN projection centered at (RA=0, Dec=0), a 1-
// pixel step in the +x direction (CDELT1 < 0) should decrease RA by
// |CDELT1| arcsec to first order. This is a sanity check that the sign
// convention is right.
func TestTANAtEquator(t *testing.T) {
	const scale = 0.001 // 3.6 arcsec per pixel, exact
	w := makeTANHeader(t, 0.0, 0.0, 100.0, 100.0, -scale, 0, 0, scale)
	tr, _ := New(w)

	// Go pixel (100, 99) is one pixel to the right of the reference in the
	// 0-based frame (reference is at 100-1=99).
	a, d, err := tr.PixelToSky(100, 99)
	if err != nil {
		t.Fatal(err)
	}
	// Expect RA ≈ -0.001° → wrapped to 359.999°, Dec ≈ 0.
	if math.Abs(a-359.999) > 1e-6 {
		t.Errorf("RA one pixel right: got %.15g, want 359.999", a)
	}
	if math.Abs(d) > 1e-6 {
		t.Errorf("Dec one pixel right: got %.15g, want 0", d)
	}
}

// TestTANCosDec: at a high declination, the RA scale shrinks by cos(dec).
// Moving 1 pixel in +x at Dec=60° should change RA by CDELT1/cos(60°) =
// 2×|CDELT1|.
func TestTANCosDec(t *testing.T) {
	const scale = 0.001
	w := makeTANHeader(t, 180.0, 60.0, 100.0, 100.0, -scale, 0, 0, scale)
	tr, _ := New(w)
	a, _, err := tr.PixelToSky(100, 99) // one pixel right of reference
	if err != nil {
		t.Fatal(err)
	}
	// Expected RA shift = -CDELT1 / cos(60°) ≈ 0.001 / 0.5 = 0.002.
	// So RA ≈ 180 - 0.002 = 179.998. (Sign: CDELT1 < 0, +x pixel → RA decreases.)
	want := 179.998
	if math.Abs(a-want) > 1e-6 {
		t.Errorf("cos(dec) scaling: got %.15g, want %.15g", a, want)
	}
}

// TestSINRoundTrip: same round-trip discipline for SIN.
func TestSINRoundTrip(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---SIN",
		"CTYPE2", "DEC--SIN",
		"CRVAL1", 45.0,
		"CRVAL2", -20.0,
		"CRPIX1", 128.0,
		"CRPIX2", 128.0,
		"CDELT1", -0.005,
		"CDELT2", 0.005,
	)
	w, _ := wcs.Parse(h)
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	for _, p := range []struct{ x, y float64 }{
		{127, 127}, {50, 200}, {200, 50}, {128, 128},
	} {
		a, d, err := tr.PixelToSky(p.x, p.y)
		if err != nil {
			t.Fatal(err)
		}
		p1, p2, err := tr.SkyToPixel(a, d)
		if err != nil {
			t.Fatal(err)
		}
		if math.Abs(p1-p.x) > 1e-9 || math.Abs(p2-p.y) > 1e-9 {
			t.Errorf("SIN roundtrip %v → (%v,%v)", p, p1, p2)
		}
	}
}

// TestSTGRoundTrip: same for stereographic.
func TestSTGRoundTrip(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---STG",
		"CTYPE2", "DEC--STG",
		"CRVAL1", 90.0,
		"CRVAL2", 45.0,
		"CRPIX1", 512.0,
		"CRPIX2", 512.0,
		"CDELT1", -0.01,
		"CDELT2", 0.01,
	)
	w, _ := wcs.Parse(h)
	tr, _ := New(w)
	a, d, err := tr.PixelToSky(256, 256)
	if err != nil {
		t.Fatal(err)
	}
	p1, p2, err := tr.SkyToPixel(a, d)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(p1-256) > 1e-9 || math.Abs(p2-256) > 1e-9 {
		t.Errorf("STG roundtrip: (%v,%v)", p1, p2)
	}
}

// TestARCRoundTrip and TestZEARoundTrip verify the remaining zenithal members.
func TestARCRoundTrip(t *testing.T) {
	assertRoundTrip(t, "ARC", 200, 10, 256, 256, 0.01)
}
func TestZEARoundTrip(t *testing.T) {
	assertRoundTrip(t, "ZEA", 60, -30, 256, 256, 0.01)
}

// TestCARRoundTrip exercises the CAR cylindrical projection.
func TestCARRoundTrip(t *testing.T) {
	assertRoundTrip(t, "CAR", 0, 0, 360, 180, 1.0)
}

// assertRoundTrip builds a simple header with the given projection code
// and verifies pixel → sky → pixel is an identity at a grid of test
// points.
func assertRoundTrip(t *testing.T, code string, crval1, crval2, crpix1, crpix2, scale float64) {
	t.Helper()
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---"+code,
		"CTYPE2", "DEC--"+code,
		"CRVAL1", crval1,
		"CRVAL2", crval2,
		"CRPIX1", crpix1,
		"CRPIX2", crpix2,
		"CDELT1", -scale,
		"CDELT2", scale,
	)
	// Fix up CTYPE length for short codes (CAR is 3 chars already, so
	// "RA---CAR" works; ARC/ZEA likewise).
	_ = h
	w, err := wcs.Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	points := []struct{ x, y float64 }{
		{crpix1 - 1, crpix2 - 1},
		{crpix1 - 1 + 10, crpix2 - 1 + 10},
		{crpix1 - 1 - 50, crpix2 - 1 + 5},
		{crpix1 - 1 + 100, crpix2 - 1 - 20},
	}
	for _, p := range points {
		a, d, err := tr.PixelToSky(p.x, p.y)
		if err != nil {
			t.Fatalf("%s PixelToSky(%v,%v): %v", code, p.x, p.y, err)
		}
		p1, p2, err := tr.SkyToPixel(a, d)
		if err != nil {
			t.Fatalf("%s SkyToPixel: %v", code, err)
		}
		if math.Abs(p1-p.x) > 1e-9 || math.Abs(p2-p.y) > 1e-9 {
			t.Errorf("%s roundtrip: (%v,%v) → (%v,%v) → (%v,%v)", code, p.x, p.y, a, d, p1, p2)
		}
	}
}

// TestAllPaperIICodesImplemented verifies that every 3-letter code in
// allProjectionCodes constructs successfully via Select (with sensible
// default PV parameters). All 27 Paper II projections plus HPX/XPH
// are implemented; none may return ErrUnsupportedProjection.
func TestAllPaperIICodesImplemented(t *testing.T) {
	// Conic projections require a non-zero PV2_1 (theta_a); supply one.
	conicPV := map[wcs.PVKey]float64{
		{Axis: 2, Index: 1}: 45.0, // theta_a = 45°
	}
	// ZPN requires at least one non-zero coefficient; supply a linear term.
	zpnPV := map[wcs.PVKey]float64{
		{Axis: 2, Index: 1}: 1.0,
	}
	for code := range allProjectionCodes {
		var pv map[wcs.PVKey]float64
		switch code {
		case "COP", "COE", "COD", "COO":
			pv = conicPV
		case "ZPN":
			pv = zpnPV
		}
		if _, err := Select(code, pv, 2); err != nil {
			t.Errorf("Select(%q): %v (every standard code must be implemented)", code, err)
		}
	}
}

// TestUnknownProjectionReturnsError rejects codes that aren't in the FITS
// projection list at all.
func TestUnknownProjectionReturnsError(t *testing.T) {
	for _, code := range []string{"XYZ", "FOO", "", "TAN1"} {
		_, err := Select(code, nil, 2)
		if err == nil {
			t.Errorf("Select(%q) returned nil for unknown code", code)
			continue
		}
		if !errors.Is(err, ErrUnknownProjection) {
			t.Errorf("Select(%q): %v (not wrapping ErrUnknownProjection)", code, err)
		}
	}
}

// TestNewAcceptsMERHeader is a sanity check that New() now builds a
// transform for MER (previously unsupported, now implemented in Phase 2).
func TestNewAcceptsMERHeader(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---MER",
		"CTYPE2", "DEC--MER",
		"CRVAL1", 0.0,
		"CRVAL2", 0.0,
		"CRPIX1", 100.0,
		"CRPIX2", 100.0,
		"CDELT1", -0.1,
		"CDELT2", 0.1,
	)
	w, _ := wcs.Parse(h)
	if _, err := New(w); err != nil {
		t.Fatalf("MER should be supported now: %v", err)
	}
}

// TestNewRejectsNonCelestial: a WCS without lon/lat axes should not
// construct a Transform.
func TestNewRejectsNonCelestial(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "FREQ",
		"CTYPE2", "WAVE",
	)
	w, _ := wcs.Parse(h)
	if _, err := New(w); err == nil {
		t.Fatal("expected error for non-celestial header")
	}
}

// TestPCMatrixPath verifies the PC matrix + CDELT code path (distinct from
// the CD path exercised by the TAN tests).
func TestPCMatrixPath(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CRVAL1", 10.0,
		"CRVAL2", 10.0,
		"CRPIX1", 100.0,
		"CRPIX2", 100.0,
		"CDELT1", -0.001,
		"CDELT2", 0.001,
		"PC1_1", 1.0,
		"PC1_2", 0.0,
		"PC2_1", 0.0,
		"PC2_2", 1.0,
	)
	w, _ := wcs.Parse(h)
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	// Reference pixel → reference sky.
	a, d, err := tr.PixelToSky(99, 99)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(a-10.0) > 1e-12 || math.Abs(d-10.0) > 1e-12 {
		t.Fatalf("PC path reference: (%v,%v)", a, d)
	}
	// Round-trip check.
	p1, p2, _ := tr.SkyToPixel(a, d)
	if math.Abs(p1-99) > 1e-9 || math.Abs(p2-99) > 1e-9 {
		t.Errorf("PC roundtrip: (%v,%v)", p1, p2)
	}
}

// TestSingularLinearMatrixRejected: a singular CD matrix must fail New.
func TestSingularLinearMatrixRejected(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CD1_1", 1.0,
		"CD1_2", 2.0,
		"CD2_1", 2.0,
		"CD2_2", 4.0, // row 2 = 2 × row 1 → singular
	)
	w, _ := wcs.Parse(h)
	if _, err := New(w); err == nil {
		t.Fatal("expected singular matrix error")
	}
}
