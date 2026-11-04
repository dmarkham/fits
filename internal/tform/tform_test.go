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

// The grammar cases below follow the FITS 4.0 standard text: binary TFORM
// "rTa" (§7.3.2, Table 18; "rAw" substrings per the registered Substring
// Array convention; P/Q descriptors §7.3.5 with repeat 0 or 1), ASCII TFORM
// "Tw[.d]" (§7.2.2, Table 15) and TDIM "(n1,...,nk)" (§7.3.2). cfitsio's
// own fits_binary_tform is not used as the oracle here: it accepts "10E5"
// as repeat 10^6 (sscanf %lf), "1PE(4)x" and "2PE(4)", and leaves the
// repeat of "1E" uninitialised, so it is too lax to judge syntax. The
// values our parser returns for real columns are checked against cfitsio
// on every committed sample in TestParseMatchesCfitsioColumns.

func TestParseBinaryGrammar(t *testing.T) {
	ok := []struct {
		in      string
		repeat  int64
		typ     BinaryType
		width   int64
		varType BinaryType
		maxLen  int64
	}{
		{"1L", 1, BinLogical, 0, 0, 0},
		{"L", 1, BinLogical, 0, 0, 0},
		{"13X", 13, BinBit, 0, 0, 0},
		{"4B", 4, BinUint8, 0, 0, 0},
		{"2I", 2, BinInt16, 0, 0, 0},
		{"J", 1, BinInt32, 0, 0, 0},
		{"1K", 1, BinInt64, 0, 0, 0},
		{"12A", 12, BinChar, 0, 0, 0},
		{"20A5", 20, BinChar, 5, 0, 0},
		{"1E", 1, BinFloat32, 0, 0, 0},
		{"0D", 0, BinFloat64, 0, 0, 0},
		{"2C", 2, BinComplex64, 0, 0, 0},
		{"1M", 1, BinComplex128, 0, 0, 0},
		{"1PE(256)", 1, BinPArrayDesc32, 0, BinFloat32, 256},
		{"PE(256)", 1, BinPArrayDesc32, 0, BinFloat32, 256},
		{"0PJ(3)", 0, BinPArrayDesc32, 0, BinInt32, 3},
		{"1PJ", 1, BinPArrayDesc32, 0, BinInt32, 0},
		{"1PL(20)", 1, BinPArrayDesc32, 0, BinLogical, 20},
		{"1PX(17)", 1, BinPArrayDesc32, 0, BinBit, 17},
		{"1PA(20)", 1, BinPArrayDesc32, 0, BinChar, 20},
		{"1PC(0)", 1, BinPArrayDesc32, 0, BinComplex64, 0},
		{"1QM(4)", 1, BinQArrayDesc64, 0, BinComplex128, 4},
		{"1QJ(1024)", 1, BinQArrayDesc64, 0, BinInt32, 1024},
		{" 1J  ", 1, BinInt32, 0, 0, 0},
	}
	for _, c := range ok {
		f, err := ParseBinary(c.in)
		if err != nil {
			t.Errorf("ParseBinary(%q): %v", c.in, err)
			continue
		}
		if f.Repeat != c.repeat || f.Type != c.typ || f.Width != c.width || f.VarType != c.varType || f.VarMaxLen != c.maxLen {
			t.Errorf("ParseBinary(%q) = %+v", c.in, f)
		}
	}
	bad := []string{"", "  ", "10E5", "1PE(4)x", "2PE(4)", "1PP(4)", "1PQ(4)", "1Z", "1PZ",
		"1PE(4", "1PE(x)", "1PE(-3)", "J1", "1W", "1U", "1V", "1S", "12A0", "1E ,", "-1J"}
	for _, in := range bad {
		if f, err := ParseBinary(in); err == nil {
			t.Errorf("ParseBinary(%q) = %+v, want error", in, f)
		}
	}
}

func TestParseASCIIGrammar(t *testing.T) {
	ok := []struct {
		in    string
		typ   ASCIIType
		width int
		frac  int
	}{
		{"A8", AsciiChar, 8, 0},
		{"I6", AsciiInt, 6, 0},
		{"F8.2", AsciiFloatF, 8, 2},
		{"E11.3", AsciiFloatE, 11, 3},
		{"D21.14", AsciiFloatD, 21, 14},
		{"F8", AsciiFloatF, 8, 0},
		{" I6 ", AsciiInt, 6, 0},
	}
	for _, c := range ok {
		f, err := ParseASCII(c.in)
		if err != nil {
			t.Errorf("ParseASCII(%q): %v", c.in, err)
			continue
		}
		if f.Type != c.typ || f.Width != c.width || f.Fraction != c.frac {
			t.Errorf("ParseASCII(%q) = %+v", c.in, f)
		}
	}
	bad := []string{"", "I10.5", "A8.2", "F.2", "F0.1", "X5", "I", "F8.x", "E12.5x", "I-3", "A0"}
	for _, in := range bad {
		if f, err := ParseASCII(in); err == nil {
			t.Errorf("ParseASCII(%q) = %+v, want error", in, f)
		}
	}
}

func TestParseDimGrammar(t *testing.T) {
	ok := map[string][]int64{"(3)": {3}, "(3,2)": {3, 2}, " (1,2,8) ": {1, 2, 8}, "(0)": {0}, "( 2 , 2 )": {2, 2}}
	for in, want := range ok {
		got, err := ParseDim(in)
		if err != nil || len(got) != len(want) {
			t.Errorf("ParseDim(%q) = %v, %v", in, got, err)
			continue
		}
		for i := range want {
			if got[i] != want[i] {
				t.Errorf("ParseDim(%q) = %v", in, got)
			}
		}
	}
	for _, in := range []string{"()", "(3,)", "(a)", "3,2", "(-1)", "(3,2", "(+3)"} {
		if got, err := ParseDim(in); err == nil {
			t.Errorf("ParseDim(%q) = %v, want error", in, got)
		}
	}
}
