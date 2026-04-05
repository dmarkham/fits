package transform

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestSIPForwardReverseIdentity constructs a minimal SIP polynomial,
// applies forward and then numerical inverse, and verifies round-trip.
func TestSIPForwardReverseIdentity(t *testing.T) {
	// Build a small SIP polynomial: a quadratic x-distortion and a linear
	// y-distortion.
	//   A(u, v) = 1e-5 * u² + 2e-6 * v²
	//   B(u, v) = -3e-6 * u + 4e-6 * v²
	sip := &SIP{
		A: &SIPPoly{Order: 2, Coeffs: [][]float64{
			{0, 0, 2e-6}, // A_0_0, A_0_1, A_0_2
			{0, 0, 0},    // A_1_0, A_1_1
			{1e-5, 0, 0}, // A_2_0
		}},
		B: &SIPPoly{Order: 2, Coeffs: [][]float64{
			{0, 0, 4e-6},
			{-3e-6, 0, 0},
			{0, 0, 0},
		}},
	}
	// Forward at (100, -50).
	u, v := 100.0, -50.0
	u2, v2 := sip.Forward(u, v)
	// Verify non-identity (distortion is applied).
	if math.Abs(u2-u) < 1e-9 || math.Abs(v2-v) < 1e-9 {
		t.Fatalf("forward had no effect: (%v,%v) → (%v,%v)", u, v, u2, v2)
	}
	// Inverse via Newton.
	u3, v3, ok := sip.Inverse(u2, v2)
	if !ok {
		t.Fatal("Inverse returned !ok")
	}
	if math.Abs(u3-u) > 1e-9 || math.Abs(v3-v) > 1e-9 {
		t.Fatalf("round trip: (%v,%v) → (%v,%v) → (%v,%v)", u, v, u2, v2, u3, v3)
	}
}

// TestParseSIPFromHeader reads SIP coefficients from a *header.Header and
// verifies the polynomial has the expected structure.
func TestParseSIPFromHeader(t *testing.T) {
	h := header.New()
	h.Set("A_ORDER", int64(2), "")
	h.Set("B_ORDER", int64(2), "")
	h.Set("A_0_2", 2e-6, "")
	h.Set("A_2_0", 1e-5, "")
	h.Set("B_1_0", -3e-6, "")
	h.Set("B_0_2", 4e-6, "")

	sip, err := ParseSIP(h)
	if err != nil {
		t.Fatal(err)
	}
	if sip == nil {
		t.Fatal("expected SIP, got nil")
	}
	if sip.A.Order != 2 || sip.B.Order != 2 {
		t.Fatalf("orders: A=%d B=%d", sip.A.Order, sip.B.Order)
	}
	if sip.A.Coeffs[2][0] != 1e-5 {
		t.Fatalf("A_2_0 = %v", sip.A.Coeffs[2][0])
	}
	if sip.B.Coeffs[0][2] != 4e-6 {
		t.Fatalf("B_0_2 = %v", sip.B.Coeffs[0][2])
	}
}

// TestParseSIPAbsent returns nil when A_ORDER is missing.
func TestParseSIPAbsent(t *testing.T) {
	h := header.New()
	h.Set("NAXIS", int64(2), "")
	sip, err := ParseSIP(h)
	if err != nil {
		t.Fatal(err)
	}
	if sip != nil {
		t.Fatal("expected nil SIP for header without A_ORDER")
	}
}

// TestSIPAPInversePolynomial: when AP/BP are present, the inverse uses
// the supplied polynomial instead of Newton iteration.
func TestSIPAPInversePolynomial(t *testing.T) {
	h := header.New()
	h.Set("A_ORDER", int64(1), "")
	h.Set("B_ORDER", int64(1), "")
	h.Set("A_1_0", 1e-3, "")
	h.Set("B_0_1", 2e-3, "")
	h.Set("AP_ORDER", int64(1), "")
	h.Set("BP_ORDER", int64(1), "")
	h.Set("AP_1_0", -1e-3, "") // approximate inverse
	h.Set("BP_0_1", -2e-3, "")

	sip, _ := ParseSIP(h)
	if sip.AP == nil || sip.BP == nil {
		t.Fatal("AP/BP not parsed")
	}
	// Inverse at (100, 50) should use AP/BP directly.
	u, v, ok := sip.Inverse(100, 50)
	if !ok {
		t.Fatal("inverse failed")
	}
	// Expect u = 100 + AP_1_0*100 = 100 - 0.1 = 99.9
	if math.Abs(u-99.9) > 1e-9 {
		t.Fatalf("AP inverse: u=%v", u)
	}
	if math.Abs(v-49.9) > 1e-9 {
		t.Fatalf("BP inverse: v=%v", v)
	}
}
