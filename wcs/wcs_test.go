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
	}{
		{"RA---TAN", "RA", "TAN"},
		{"DEC--TAN", "DEC", "TAN"},
		{"GLON-CAR", "GLON", "CAR"},
		{"GLAT-CAR", "GLAT", "CAR"},
		{"ELON-AIT", "ELON", "AIT"},
		{"RA---SIN", "RA", "SIN"},
		{"DEC--SIN", "DEC", "SIN"},
		{"FREQ", "FREQ", ""},
		{"WAVE", "WAVE", ""},
		{"", "", ""},
	}
	for _, c := range cases {
		at, pc := splitCType(c.in)
		if at != c.axisType || pc != c.projCode {
			t.Errorf("splitCType(%q) = (%q,%q), want (%q,%q)", c.in, at, pc, c.axisType, c.projCode)
		}
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
