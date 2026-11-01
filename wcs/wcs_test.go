package wcs

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/header"
)

// newHeader builds a *header.Header from a flat list of (key, value) pairs
// for test convenience.
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

// TestParseMinimal: NAXIS=2 with defaults only. Everything should fall to
// sane defaults per Paper I.
func TestParseMinimal(t *testing.T) {
	h := newHeader("NAXIS", int64(2), "NAXIS1", int64(100), "NAXIS2", int64(100))
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.NAxis != 2 {
		t.Fatalf("NAxis=%d", w.NAxis)
	}
	if w.CDelt[0] != 1 || w.CDelt[1] != 1 {
		t.Fatalf("CDELT default not 1: %v", w.CDelt)
	}
	if w.CRPix[0] != 0 || w.CRPix[1] != 0 {
		t.Fatalf("CRPIX default not 0: %v", w.CRPix)
	}
	// PC should default to identity.
	if len(w.PC) != 2 || w.PC[0][0] != 1 || w.PC[1][1] != 1 || w.PC[0][1] != 0 {
		t.Fatalf("PC default not identity: %v", w.PC)
	}
	if w.CD != nil {
		t.Fatalf("CD should be nil when not present")
	}
	if w.IsCelestial() {
		t.Fatalf("minimal header should not be celestial")
	}
}

// TestParseTANHeader parses a typical HST-style TAN header with all
// standard keywords.
func TestParseTANHeader(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"NAXIS1", int64(512),
		"NAXIS2", int64(512),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CUNIT1", "deg",
		"CUNIT2", "deg",
		"CRVAL1", 180.0,
		"CRVAL2", 30.0,
		"CRPIX1", 256.5,
		"CRPIX2", 256.5,
		"CDELT1", -0.0002777778,
		"CDELT2", 0.0002777778,
		"RADESYS", "ICRS",
		"EQUINOX", 2000.0,
		"DATE-OBS", "2025-04-05T12:34:56",
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() {
		t.Fatalf("not celestial")
	}
	if w.LonAxis != 1 || w.LatAxis != 2 {
		t.Fatalf("axes: lon=%d lat=%d", w.LonAxis, w.LatAxis)
	}
	if w.AxisType[0] != "RA" || w.AxisType[1] != "DEC" {
		t.Fatalf("axis types: %v", w.AxisType)
	}
	if w.ProjCode[0] != "TAN" || w.ProjCode[1] != "TAN" {
		t.Fatalf("proj codes: %v", w.ProjCode)
	}
	if w.CelestialProjCode() != "TAN" {
		t.Fatalf("celestial code")
	}
	if w.CRVal[0] != 180.0 || w.CRVal[1] != 30.0 {
		t.Fatalf("CRVAL: %v", w.CRVal)
	}
	if math.Abs(w.CDelt[0]+0.0002777778) > 1e-15 {
		t.Fatalf("CDELT1: %v", w.CDelt[0])
	}
	if w.RadeSys != "ICRS" {
		t.Fatalf("RADESYS: %q", w.RadeSys)
	}
	if w.Equinox != 2000.0 {
		t.Fatalf("EQUINOX: %v", w.Equinox)
	}
	if w.DateObs != "2025-04-05T12:34:56" {
		t.Fatalf("DATE-OBS: %q", w.DateObs)
	}
}

// TestParseCDMatrix: CD matrix form instead of PC+CDELT.
func TestParseCDMatrix(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CRVAL1", 10.0,
		"CRVAL2", 20.0,
		"CRPIX1", 100.0,
		"CRPIX2", 100.0,
		"CD1_1", 1e-5,
		"CD1_2", 2e-6,
		"CD2_1", 3e-6,
		"CD2_2", -1e-5,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.CD == nil {
		t.Fatalf("CD matrix not populated")
	}
	if w.CD[0][0] != 1e-5 || w.CD[0][1] != 2e-6 || w.CD[1][0] != 3e-6 || w.CD[1][1] != -1e-5 {
		t.Fatalf("CD matrix: %v", w.CD)
	}
}

// TestParseConflictingMatrices: PC and CD both present → error.
func TestParseConflictingMatrices(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"PC1_1", 1.0,
		"CD1_1", 1e-5,
	)
	_, err := Parse(h)
	if err != ErrConflictingMatrix {
		t.Fatalf("expected ErrConflictingMatrix, got %v", err)
	}
}

// TestParsePCMatrix: explicit PC matrix (non-identity).
func TestParsePCMatrix(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CDELT1", 0.01,
		"CDELT2", 0.01,
		"PC1_1", 0.707,
		"PC1_2", -0.707,
		"PC2_1", 0.707,
		"PC2_2", 0.707,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(w.PC[0][0]-0.707) > 1e-15 {
		t.Fatalf("PC[0][0]: %v", w.PC[0][0])
	}
	if math.Abs(w.PC[0][1]+0.707) > 1e-15 {
		t.Fatalf("PC[0][1]: %v", w.PC[0][1])
	}
}

// TestParsePV: PVi_m projection parameters.
func TestParsePV(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---SZP",
		"CTYPE2", "DEC--SZP",
		"PV2_1", 0.5,
		"PV2_2", 1.5,
		"PV2_3", 2.5,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.PV[PVKey{Axis: 2, Index: 1}] != 0.5 {
		t.Fatalf("PV2_1: %v", w.PV[PVKey{Axis: 2, Index: 1}])
	}
	if w.PV[PVKey{Axis: 2, Index: 3}] != 2.5 {
		t.Fatalf("PV2_3: %v", w.PV[PVKey{Axis: 2, Index: 3}])
	}
}

// TestSplitCTypeVariants covers the quirky CTYPE formats seen in the wild.
func TestSplitCTypeVariants(t *testing.T) {
	cases := []struct {
		in       string
		axisType string
		projCode string
		distCode string
	}{
		{"RA---TAN", "RA", "TAN", ""},
		{"DEC--TAN", "DEC", "TAN", ""},
		{"GLON-CAR", "GLON", "CAR", ""},
		{"GLAT-CAR", "GLAT", "CAR", ""},
		{"ELON-AIT", "ELON", "AIT", ""},
		{"RA---SIN", "RA", "SIN", ""},
		{"DEC--SIN", "DEC", "SIN", ""},
		{"RA---TAN-SIP", "RA", "TAN", "SIP"},
		{"DEC--TAN-SIP", "DEC", "TAN", "SIP"},
		{"ra---tan-sip", "RA", "TAN", "SIP"},
		{"RA---TPV", "RA", "TPV", ""},
		{"FREQ", "FREQ", "", ""},
		{"WAVE", "WAVE", "", ""},
		{"", "", "", ""},
	}
	for _, c := range cases {
		at, pc, dc := splitCType(c.in)
		if at != c.axisType || pc != c.projCode || dc != c.distCode {
			t.Errorf("splitCType(%q) = (%q,%q,%q), want (%q,%q,%q)",
				c.in, at, pc, dc, c.axisType, c.projCode, c.distCode)
		}
	}
}

// sipCards is a TAN WCS with a quadratic SIP distortion and its inverse.
var sipCards = []any{
	"NAXIS", int64(2),
	"CRVAL1", 101.5, "CRVAL2", 28.9,
	"CRPIX1", 2640.5, "CRPIX2", 2640.5,
	"CD1_1", 4.5e-05, "CD1_2", 0.0, "CD2_1", 0.0, "CD2_2", -4.5e-05,
	"A_ORDER", int64(2), "A_2_0", -5.4e-08, "A_1_1", -1.0e-07,
	"B_ORDER", int64(2), "B_0_2", -7.1e-08,
	"AP_ORDER", int64(2), "AP_2_0", 5.4e-08,
	"BP_ORDER", int64(2), "BP_0_2", 7.1e-08,
}

// TestParseSIPSuffix: "RA---TAN-SIP" is a celestial TAN axis with a SIP
// distortion, and the SIP coefficients are read into Header.SIP.
func TestParseSIPSuffix(t *testing.T) {
	h := newHeader(append([]any{"CTYPE1", "RA---TAN-SIP", "CTYPE2", "DEC--TAN-SIP"}, sipCards...)...)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() {
		t.Fatalf("RA---TAN-SIP not celestial: AxisType=%v", w.AxisType)
	}
	if got := w.CelestialProjCode(); got != "TAN" {
		t.Errorf("CelestialProjCode() = %q, want TAN", got)
	}
	if w.DistCode[0] != "SIP" || w.DistCode[1] != "SIP" {
		t.Errorf("DistCode = %v, want [SIP SIP]", w.DistCode)
	}
	if w.SIP == nil {
		t.Fatal("SIP coefficients not parsed")
	}
	if w.SIP.A.Order != 2 || w.SIP.A.Coeffs[2][0] != -5.4e-08 || w.SIP.A.Coeffs[1][1] != -1.0e-07 {
		t.Errorf("A = %+v", w.SIP.A)
	}
	if w.SIP.B.Coeffs[0][2] != -7.1e-08 {
		t.Errorf("B_0_2 = %v", w.SIP.B.Coeffs[0][2])
	}
	if w.SIP.AP == nil || w.SIP.AP.Coeffs[2][0] != 5.4e-08 || w.SIP.BP == nil || w.SIP.BP.Coeffs[0][2] != 7.1e-08 {
		t.Errorf("AP/BP = %+v / %+v", w.SIP.AP, w.SIP.BP)
	}
}

// TestParseSIPKeysWithoutSuffix: SIP keywords are only read when CTYPE
// declares "-SIP". ATLAS difference images, for example, declare TPV and
// carry a redundant SIP set that must not be applied on top of it.
func TestParseSIPKeysWithoutSuffix(t *testing.T) {
	for _, ct := range [][2]string{{"RA---TAN", "DEC--TAN"}, {"RA---TPV", "DEC--TPV"}} {
		h := newHeader(append([]any{"CTYPE1", ct[0], "CTYPE2", ct[1]}, sipCards...)...)
		w, err := Parse(h)
		if err != nil {
			t.Fatal(err)
		}
		if w.SIP != nil {
			t.Errorf("%s: SIP parsed without a -SIP suffix", ct[0])
		}
		if w.DistCode[0] != "" {
			t.Errorf("%s: DistCode = %q, want empty", ct[0], w.DistCode[0])
		}
	}
}

// TestParseSIPSuffixWithoutCoefficients: a "-SIP" CTYPE with no A_ORDER /
// B_ORDER parses as celestial with no SIP polynomial (as astropy.wcs).
func TestParseSIPSuffixWithoutCoefficients(t *testing.T) {
	h := newHeader("NAXIS", int64(2), "CTYPE1", "RA---TAN-SIP", "CTYPE2", "DEC--TAN-SIP",
		"CRVAL1", 10.0, "CRVAL2", 20.0, "CDELT1", -1e-3, "CDELT2", 1e-3)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() || w.CelestialProjCode() != "TAN" {
		t.Fatalf("not a celestial TAN header: %v", w)
	}
	if w.DistCode[0] != "SIP" {
		t.Errorf("DistCode = %v", w.DistCode)
	}
	if w.SIP != nil {
		t.Errorf("SIP = %+v, want nil", w.SIP)
	}
}

// TestParseRADECSYSFallback: the pre-1998 RADECSYS spelling should still be
// recognized.
func TestParseRADECSYSFallback(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"RADECSYS", "FK5",
	)
	w, _ := Parse(h)
	if w.RadeSys != "FK5" {
		t.Fatalf("RADECSYS fallback: %q", w.RadeSys)
	}
}

// TestParseEpochFallback: pre-1988 EPOCH is a synonym for EQUINOX.
func TestParseEpochFallback(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"EPOCH", 1950.0,
	)
	w, _ := Parse(h)
	if w.Equinox != 1950.0 {
		t.Fatalf("EPOCH fallback: %v", w.Equinox)
	}
}

// TestParseRejectsInvalidSIPOrder: a "-SIP" header with an order outside
// 0..MaxSIPOrder is an error, never a panic. Without the suffix the SIP
// keywords are not read, so the same orders are harmless.
func TestParseRejectsInvalidSIPOrder(t *testing.T) {
	for _, order := range []int64{-2, MaxSIPOrder + 1, 1 << 40} {
		for _, suffix := range []string{"-SIP", ""} {
			h := newHeader("NAXIS", int64(2),
				"CTYPE1", "RA---TAN"+suffix, "CTYPE2", "DEC--TAN"+suffix,
				"CDELT1", -1e-3, "CDELT2", 1e-3,
				"A_ORDER", order, "B_ORDER", int64(2))
			func() {
				defer func() {
					if r := recover(); r != nil {
						t.Errorf("A_ORDER=%d %q: Parse panicked: %v", order, suffix, r)
					}
				}()
				_, err := Parse(h)
				if suffix != "" && err == nil {
					t.Errorf("A_ORDER=%d: want an error", order)
				}
				if suffix == "" && err != nil {
					t.Errorf("A_ORDER=%d without -SIP: unexpected error %v", order, err)
				}
			}()
		}
	}
}
