package wcs

import (
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestParseAltDescription: a header with both a primary WCS (no suffix)
// and an alternate 'A' description must let the caller select which one
// to parse via ParseAlt.
func TestParseAltDescription(t *testing.T) {
	h := header.New()
	h.Set("NAXIS", int64(2), "")
	h.Set("CTYPE1", "RA---TAN", "")
	h.Set("CTYPE2", "DEC--TAN", "")
	h.Set("CRVAL1", 180.0, "")
	h.Set("CRVAL2", 30.0, "")
	h.Set("CRPIX1", 256.0, "")
	h.Set("CRPIX2", 256.0, "")
	// Alternate A: pixel-pixel identity (the cfitsio iter_image.fit pattern).
	h.Set("CTYPE1A", "X", "")
	h.Set("CTYPE2A", "Y", "")
	h.Set("CRVAL1A", 1.0, "")
	h.Set("CRVAL2A", 1.0, "")
	h.Set("CD1_1A", 1.0, "")
	h.Set("CD2_2A", 1.0, "")

	// Primary parse.
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.CType[0] != "RA---TAN" {
		t.Fatalf("primary CTYPE1: %q", w.CType[0])
	}
	if w.CRVal[0] != 180.0 || w.CRVal[1] != 30.0 {
		t.Fatalf("primary CRVAL: %v", w.CRVal)
	}
	if !w.IsCelestial() {
		t.Fatal("primary should be celestial")
	}

	// Alternate A parse.
	wA, err := ParseAlt(h, "A")
	if err != nil {
		t.Fatal(err)
	}
	if wA.CType[0] != "X" {
		t.Fatalf("alt A CTYPE1: %q", wA.CType[0])
	}
	if wA.CRVal[0] != 1.0 || wA.CRVal[1] != 1.0 {
		t.Fatalf("alt A CRVAL: %v", wA.CRVal)
	}
	if wA.CD == nil || wA.CD[0][0] != 1.0 || wA.CD[1][1] != 1.0 {
		t.Fatalf("alt A CD: %v", wA.CD)
	}
	// Alt description is not celestial (CTYPE is "X"/"Y" with no projection).
	if wA.IsCelestial() {
		t.Fatal("alt A should not be celestial (CTYPE=X/Y)")
	}
}

// TestParseAltInvalidIdentifier rejects malformed alt identifiers.
func TestParseAltInvalidIdentifier(t *testing.T) {
	h := header.New()
	h.Set("NAXIS", int64(2), "")
	if _, err := ParseAlt(h, "AB"); err == nil {
		t.Error("expected error for 2-char alt")
	}
	if _, err := ParseAlt(h, "1"); err == nil {
		t.Error("expected error for digit alt")
	}
	if _, err := ParseAlt(h, "a"); err == nil {
		t.Error("expected error for lowercase alt")
	}
}
