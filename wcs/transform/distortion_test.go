package transform

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// TestTPVRoundTrip builds a minimal TPV distortion (identity + small
// quadratic perturbation) and verifies Forward then Inverse recovers the
// original intermediate coordinates.
func TestTPVRoundTrip(t *testing.T) {
	tpv := &TPV{Lon: 1, Lat: 2}
	// Identity baseline.
	tpv.Axis1[1] = 1.0 // coefficient of xi in the x output
	tpv.Axis2[1] = 1.0 // coefficient of eta in the y output (PV2_1 is eta)
	// Small quadratic distortion.
	tpv.Axis1[4] = 1e-4 // xi² term
	tpv.Axis2[4] = 2e-4 // eta² term (PV2_4 is eta²)

	xi, eta := 0.1, 0.15 // degrees
	xiP, etaP := tpv.Forward(xi, eta)
	xi2, eta2, ok := tpv.Inverse(xiP, etaP)
	if !ok {
		t.Fatal("Inverse !ok")
	}
	if math.Abs(xi2-xi) > 1e-9 || math.Abs(eta2-eta) > 1e-9 {
		t.Fatalf("TPV round-trip: (%v,%v) → (%v,%v) → (%v,%v)", xi, eta, xiP, etaP, xi2, eta2)
	}
}

// TestTPVRadialTerms exercises the radial coefficients (slots 3, 11, 23,
// 39 — r, r³, r⁵, r⁷), which a previous bug in rPow initialization had
// silently zeroed. SCAMP-generated TPV headers for wide-field optical
// surveys routinely use r³ for the dominant radial distortion.
func TestTPVRadialTerms(t *testing.T) {
	tpv := &TPV{Lon: 1, Lat: 2}
	// Identity linear baseline + small r and r³ radial distortion.
	tpv.Axis1[1] = 1.0  // xi
	tpv.Axis2[1] = 1.0  // eta (PV2_1)
	tpv.Axis1[3] = 1e-3 // r in the x output
	tpv.Axis2[3] = 1e-3 // r in the y output
	tpv.Axis1[11] = 1e-5 // r³ in the x output
	tpv.Axis2[11] = 1e-5 // r³ in the y output

	// Pick a mid-field point where r ~ 0.3°.
	xi, eta := 0.2, 0.225
	r := math.Hypot(xi, eta)

	xiP, etaP := tpv.Forward(xi, eta)
	// xiP should = xi + 1e-3·r + 1e-5·r³. Verify the radial
	// contribution is nonzero (catches the old rPow bug):
	expectedXiP := xi + 1e-3*r + 1e-5*r*r*r
	if math.Abs(xiP-expectedXiP) > 1e-12 {
		t.Fatalf("TPV radial forward: xiP=%.15g want %.15g (diff %g) — r, r³ probably zeroed",
			xiP, expectedXiP, xiP-expectedXiP)
	}

	// Round-trip.
	xi2, eta2, ok := tpv.Inverse(xiP, etaP)
	if !ok {
		t.Fatal("Inverse !ok")
	}
	if math.Abs(xi2-xi) > 1e-11 || math.Abs(eta2-eta) > 1e-11 {
		t.Fatalf("TPV radial round-trip: (%v,%v) → (%v,%v) → (%v,%v)",
			xi, eta, xiP, etaP, xi2, eta2)
	}
}

// TestTPVAxis2TermOrder pins the axis-2 term layout: the same table as
// axis 1 with xi and eta exchanged (wcslib tpv.c, SCAMP), so PV2_1 is
// eta, PV2_2 is xi, PV2_4 is eta², PV2_6 is xi², PV2_8 is eta²·xi.
func TestTPVAxis2TermOrder(t *testing.T) {
	tpv := &TPV{Lon: 1, Lat: 2}
	tpv.Axis1[1] = 1.0
	tpv.Axis1[2] = 0.5  // eta
	tpv.Axis1[4] = 0.25 // xi²
	tpv.Axis1[8] = 0.125 // xi²·eta
	tpv.Axis2[1] = 1.0
	tpv.Axis2[2] = 0.5  // xi
	tpv.Axis2[4] = 0.25 // eta²
	tpv.Axis2[8] = 0.125 // eta²·xi

	xi, eta := 0.3, -0.2
	xiP, etaP := tpv.Forward(xi, eta)
	wantXi := xi + 0.5*eta + 0.25*xi*xi + 0.125*xi*xi*eta
	wantEta := eta + 0.5*xi + 0.25*eta*eta + 0.125*eta*eta*xi
	if math.Abs(xiP-wantXi) > 1e-15 {
		t.Errorf("axis 1: got %.15g, want %.15g", xiP, wantXi)
	}
	if math.Abs(etaP-wantEta) > 1e-15 {
		t.Errorf("axis 2: got %.15g, want %.15g (xi and eta not exchanged?)", etaP, wantEta)
	}
}

// TestTPVAxis2CrossTermRoundTrip: an axis-2 polynomial dominated by
// cross terms only converges if Inverse uses axis 2's partials in the
// exchanged order.
func TestTPVAxis2CrossTermRoundTrip(t *testing.T) {
	tpv := &TPV{Lon: 1, Lat: 2}
	tpv.Axis1[1] = 1.0
	tpv.Axis2[1] = 1.0
	tpv.Axis2[2] = 0.4  // xi
	tpv.Axis2[6] = 0.3  // xi²
	tpv.Axis2[10] = 0.2 // xi³

	xi, eta := 0.25, 0.1
	xiP, etaP := tpv.Forward(xi, eta)
	xi2, eta2, ok := tpv.Inverse(xiP, etaP)
	if !ok {
		t.Fatal("Inverse !ok")
	}
	if math.Abs(xi2-xi) > 1e-12 || math.Abs(eta2-eta) > 1e-12 {
		t.Fatalf("round-trip: (%v,%v) → (%v,%v) → (%v,%v)", xi, eta, xiP, etaP, xi2, eta2)
	}
}

// TestParseTPVDefaultIsIdentity: ParseTPV defaults PV1_1 and PV2_1 to 1,
// which must leave eta unchanged when only axis-1 terms are given.
func TestParseTPVDefaultIsIdentity(t *testing.T) {
	tpv := ParseTPV(&wcs.Header{
		LonAxis: 1,
		LatAxis: 2,
		PV:      map[wcs.PVKey]float64{{Axis: 1, Index: 4}: 1e-4},
	})
	if tpv == nil {
		t.Fatal("expected TPV, got nil")
	}
	xi, eta := 0.1, 0.2
	xiP, etaP := tpv.Forward(xi, eta)
	if math.Abs(xiP-(xi+1e-4*xi*xi)) > 1e-15 || etaP != eta {
		t.Fatalf("Forward(%v,%v) = (%v,%v), want (%v,%v)", xi, eta, xiP, etaP, xi+1e-4*xi*xi, eta)
	}
}

// TestParseTPV reads TPV coefficients from a *wcs.Header and verifies the
// coefficient map was populated correctly.
func TestParseTPV(t *testing.T) {
	w := &wcs.Header{
		LonAxis: 1,
		LatAxis: 2,
		PV: map[wcs.PVKey]float64{
			{Axis: 1, Index: 1}: 1.0,
			{Axis: 1, Index: 4}: 1e-4,
			{Axis: 2, Index: 2}: 1.0,
			{Axis: 2, Index: 6}: 2e-4,
		},
	}
	tpv := ParseTPV(w)
	if tpv == nil {
		t.Fatal("expected TPV, got nil")
	}
	if tpv.Axis1[1] != 1.0 || tpv.Axis1[4] != 1e-4 {
		t.Fatalf("Axis1: %+v", tpv.Axis1[:5])
	}
	if tpv.Axis2[2] != 1.0 || tpv.Axis2[6] != 2e-4 {
		t.Fatalf("Axis2: %+v", tpv.Axis2[:8])
	}
}

// TestParseTPVAbsent returns nil when no PV coefficients are present.
func TestParseTPVAbsent(t *testing.T) {
	w := &wcs.Header{LonAxis: 1, LatAxis: 2, PV: map[wcs.PVKey]float64{}}
	if ParseTPV(w) != nil {
		t.Fatal("expected nil for empty PV")
	}
}

// TestTNXParse constructs a synthetic TNX WAT header and parses it. Uses
// a simple polynomial (surface type 3) with identity-like coefficients
// and verifies the parser round-trips the key fields.
func TestTNXParse(t *testing.T) {
	h := header.New()
	// wtype=tnx axtype=ra lngcor = "3. 3. 3. 0. -1 1 -1 1  0 0 0 0 0 0 0 0 0"
	// That is: simple polynomial (3), xorder=3, yorder=3, cross=0 (none),
	// domain [-1,1] × [-1,1], then 9 zero coefficients.
	watValue := `wtype=tnx axtype=ra projection=tan lngcor = "3. 3. 3. 0. -1 1 -1 1  0 0 0 0 0 0 0 0 0"`
	h.Set("WAT1_001", watValue, "")

	tnx, err := ParseTNX(h)
	if err != nil {
		t.Fatal(err)
	}
	if tnx == nil {
		t.Fatal("expected TNX, got nil")
	}
	if tnx.Lon.Type != 3 {
		t.Fatalf("surface type: %d", tnx.Lon.Type)
	}
	if tnx.Lon.XOrder != 3 || tnx.Lon.YOrder != 3 {
		t.Fatalf("order: %d %d", tnx.Lon.XOrder, tnx.Lon.YOrder)
	}
	if tnx.Lon.XMin != -1 || tnx.Lon.XMax != 1 {
		t.Fatalf("x domain: [%v, %v]", tnx.Lon.XMin, tnx.Lon.XMax)
	}
}

// TestTNXRoundTrip builds a TNX surface with non-trivial quadratic
// coefficients across all three surface types (simple / Chebyshev /
// Legendre) and verifies Forward then Inverse recovers the original
// intermediate coordinates via the analytic Jacobian path.
func TestTNXRoundTrip(t *testing.T) {
	for _, surfType := range []int{1, 2, 3} {
		name := map[int]string{1: "chebyshev", 2: "legendre", 3: "simple"}[surfType]
		t.Run(name, func(t *testing.T) {
			tnx := &TNX{
				Lon: tnxSurface{
					Type: surfType, XOrder: 3, YOrder: 3, Cross: 1,
					XMin: -0.5, XMax: 0.5, YMin: -0.5, YMax: 0.5,
					// Row-major coefficient list: (i, j) pairs.
					// We use small non-trivial values so the Jacobian
					// is non-identity.
					Coeffs:  []float64{0.0, 0.001, -0.0005, 0.001, 0.0002, 0, -0.0005, 0, 0},
					present: true,
				},
				Lat: tnxSurface{
					Type: surfType, XOrder: 3, YOrder: 3, Cross: 1,
					XMin: -0.5, XMax: 0.5, YMin: -0.5, YMax: 0.5,
					Coeffs:  []float64{0.0, -0.0003, 0.0008, 0.0007, -0.0001, 0, 0.0004, 0, 0},
					present: true,
				},
			}
			xi, eta := 0.1, 0.12
			xiP, etaP := tnx.Forward(xi, eta)
			xi2, eta2, ok := tnx.Inverse(xiP, etaP)
			if !ok {
				t.Fatal("Inverse !ok")
			}
			if math.Abs(xi2-xi) > 1e-11 || math.Abs(eta2-eta) > 1e-11 {
				t.Errorf("round-trip: (%v,%v) → (%v,%v) → (%v,%v)",
					xi, eta, xiP, etaP, xi2, eta2)
			}
		})
	}
}

// TestTNXAbsent returns nil when the header has no WAT keywords.
func TestTNXAbsent(t *testing.T) {
	h := header.New()
	h.Set("NAXIS", int64(2), "")
	tnx, err := ParseTNX(h)
	if err != nil {
		t.Fatal(err)
	}
	if tnx != nil {
		t.Fatal("expected nil for header without WAT")
	}
}
