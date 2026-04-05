package tform

import (
	"testing"
)

func TestParseBinarySimple(t *testing.T) {
	f, err := ParseBinary("10E")
	if err != nil {
		t.Fatal(err)
	}
	if f.Repeat != 10 || f.Type != BinFloat32 {
		t.Fatalf("%+v", f)
	}
}

func TestParseBinaryDefaultRepeat(t *testing.T) {
	f, err := ParseBinary("J")
	if err != nil {
		t.Fatal(err)
	}
	if f.Repeat != 1 || f.Type != BinInt32 {
		t.Fatalf("%+v", f)
	}
}

func TestParseBinaryVLA(t *testing.T) {
	f, err := ParseBinary("1PE(256)")
	if err != nil {
		t.Fatal(err)
	}
	if f.Type != BinPArrayDesc32 || f.VarType != BinFloat32 || f.VarMaxLen != 256 {
		t.Fatalf("%+v", f)
	}
}

func TestParseBinaryQVLA(t *testing.T) {
	f, err := ParseBinary("1QJ(1024)")
	if err != nil {
		t.Fatal(err)
	}
	if f.Type != BinQArrayDesc64 || f.VarType != BinInt32 || f.VarMaxLen != 1024 {
		t.Fatalf("%+v", f)
	}
}

func TestParseBinaryBadType(t *testing.T) {
	if _, err := ParseBinary("10Z"); err == nil {
		t.Fatal("expected error")
	}
}

func TestParseASCIIInt(t *testing.T) {
	f, err := ParseASCII("I10")
	if err != nil {
		t.Fatal(err)
	}
	if f.Type != AsciiInt || f.Width != 10 || f.Fraction != 0 {
		t.Fatalf("%+v", f)
	}
}

func TestParseASCIIFloat(t *testing.T) {
	f, err := ParseASCII("F10.3")
	if err != nil {
		t.Fatal(err)
	}
	if f.Type != AsciiFloatF || f.Width != 10 || f.Fraction != 3 {
		t.Fatalf("%+v", f)
	}
}

func TestParseASCIIChar(t *testing.T) {
	f, err := ParseASCII("A8")
	if err != nil {
		t.Fatal(err)
	}
	if f.Type != AsciiChar || f.Width != 8 {
		t.Fatalf("%+v", f)
	}
}

func TestParseDim(t *testing.T) {
	d, err := ParseDim("(3,4,5)")
	if err != nil {
		t.Fatal(err)
	}
	if len(d) != 3 || d[0] != 3 || d[1] != 4 || d[2] != 5 {
		t.Fatalf("%v", d)
	}
}

func TestBinElementSize(t *testing.T) {
	cases := []struct {
		t BinaryType
		s int
	}{
		{BinLogical, 1},
		{BinUint8, 1},
		{BinInt16, 2},
		{BinInt32, 4},
		{BinInt64, 8},
		{BinFloat32, 4},
		{BinFloat64, 8},
		{BinComplex64, 8},
		{BinComplex128, 16},
		{BinPArrayDesc32, 8},
		{BinQArrayDesc64, 16},
	}
	for _, c := range cases {
		if got := c.t.ElementSize(); got != c.s {
			t.Errorf("%c: got %d want %d", c.t, got, c.s)
		}
	}
}
