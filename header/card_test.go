package header

import (
	"math"
	"testing"
)

// mkCard builds an 80-byte card from s, right-padded with spaces.
func mkCard(s string) []byte {
	if len(s) > CardWidth {
		panic("card too long: " + s)
	}
	b := make([]byte, CardWidth)
	copy(b, s)
	for i := len(s); i < CardWidth; i++ {
		b[i] = ' '
	}
	return b
}

func TestDecodeEnd(t *testing.T) {
	c, err := DecodeCard(mkCard("END"))
	if err != nil {
		t.Fatal(err)
	}
	if !c.IsEnd() || c.Key != "END" {
		t.Fatalf("not END: %+v", c)
	}
}

func TestDecodeString(t *testing.T) {
	c, err := DecodeCard(mkCard("OBJECT  = 'NGC 5194'           / galaxy name"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Key != "OBJECT" || c.Type != TypeString {
		t.Fatalf("%+v", c)
	}
	if c.Value.(string) != "NGC 5194" {
		t.Fatalf("value %q", c.Value)
	}
	if c.Comment != "galaxy name" {
		t.Fatalf("comment %q", c.Comment)
	}
}

func TestDecodeStringWithEscapedQuote(t *testing.T) {
	// FITS escapes single quote as '' inside a string value.
	c, err := DecodeCard(mkCard("AUTHOR  = 'O''Malley' / name"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Value.(string) != "O'Malley" {
		t.Fatalf("escape failed: %q", c.Value)
	}
}

func TestDecodeInteger(t *testing.T) {
	c, err := DecodeCard(mkCard("NAXIS   =                    2 / number of axes"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Type != TypeInt || c.Value.(int64) != 2 {
		t.Fatalf("%+v", c)
	}
}

func TestDecodeIntegerNegative(t *testing.T) {
	c, err := DecodeCard(mkCard("BITPIX  =                  -32 / bits per pixel"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Type != TypeInt || c.Value.(int64) != -32 {
		t.Fatalf("%+v", c)
	}
}

func TestDecodeFloat(t *testing.T) {
	c, err := DecodeCard(mkCard("BSCALE  =              1.5E-03 / scale factor"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Type != TypeFloat {
		t.Fatalf("not float: %+v", c)
	}
	if math.Abs(c.Value.(float64)-0.0015) > 1e-12 {
		t.Fatalf("value %v", c.Value)
	}
}

func TestDecodeFloatDExponent(t *testing.T) {
	c, err := DecodeCard(mkCard("BZERO   =              3.14D+02"))
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(c.Value.(float64)-314.0) > 1e-12 {
		t.Fatalf("D-exponent: %v", c.Value)
	}
}

func TestDecodeLogical(t *testing.T) {
	c, err := DecodeCard(mkCard("SIMPLE  =                    T / standard FITS"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Type != TypeLogical || c.Value.(bool) != true {
		t.Fatalf("%+v", c)
	}
	c2, _ := DecodeCard(mkCard("EXTEND  =                    F"))
	if c2.Value.(bool) != false {
		t.Fatalf("false: %+v", c2)
	}
}

func TestDecodeComment(t *testing.T) {
	c, err := DecodeCard(mkCard("COMMENT   This is a user comment."))
	if err != nil {
		t.Fatal(err)
	}
	if c.Key != "COMMENT" || c.Type != TypeEmpty {
		t.Fatalf("%+v", c)
	}
	if c.Comment != "  This is a user comment." {
		// The space after "COMMENT" is part of the user text (col 9 onward).
		// Our decoder keeps leading spaces; just verify trim-right-only.
		t.Fatalf("comment %q", c.Comment)
	}
}

func TestDecodeHistory(t *testing.T) {
	c, _ := DecodeCard(mkCard("HISTORY Calibrated 2026-04-05 by fits"))
	if c.Key != "HISTORY" {
		t.Fatalf("%+v", c)
	}
}

func TestDecodeBlank(t *testing.T) {
	c, _ := DecodeCard(mkCard("          just some text"))
	if c.Key != "" || c.Type != TypeEmpty {
		t.Fatalf("%+v", c)
	}
}

func TestDecodeHierarch(t *testing.T) {
	c, err := DecodeCard(mkCard("HIERARCH ESO DET NAME = 'CCD42-80' / detector"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Key != "HIERARCH ESO DET NAME" {
		t.Fatalf("key %q", c.Key)
	}
	if c.Value.(string) != "CCD42-80" {
		t.Fatalf("value %q", c.Value)
	}
}

func TestDecodeComplex(t *testing.T) {
	c, err := DecodeCard(mkCard("PHASE   = (3.14, -2.72)"))
	if err != nil {
		t.Fatal(err)
	}
	if c.Type != TypeComplexFloat {
		t.Fatalf("%+v", c)
	}
	v := c.Value.(Complex)
	if math.Abs(v.Re-3.14) > 1e-12 || math.Abs(v.Im+2.72) > 1e-12 {
		t.Fatalf("%+v", v)
	}
}

func TestDecodeBadLength(t *testing.T) {
	if _, err := DecodeCard([]byte("too short")); err == nil {
		t.Fatalf("expected error")
	}
}
