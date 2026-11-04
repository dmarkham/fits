package compress

import (
	"math"
	"testing"
)

// The FitsRandoms table is checked against libcfitsio itself in
// randoms_cref_test.go (TestFitsRandomsMatchCfitsio).

// TestDequantizeNoDither: no random offset, plain linear reconstruction.
func TestDequantizeNoDither(t *testing.T) {
	input := []int32{0, 10, 20, 30, 40}
	out := make([]float64, len(input))
	Dequantize(input, out, 0.5, 100.0, NoDither, 1, 1)
	for i, v := range input {
		want := float64(v)*0.5 + 100.0
		if out[i] != want {
			t.Errorf("pixel %d: got %v want %v", i, out[i], want)
		}
	}
}

// TestDequantizeSubtractive1: verify the dither offset is applied and
// the per-pixel sequence advances through the lookup table.
func TestDequantizeSubtractive1(t *testing.T) {
	randoms := FitsRandoms()
	input := []int32{100, 100, 100, 100}
	out := make([]float64, len(input))
	Dequantize(input, out, 1.0, 0.0, SubtractiveDither1, 1, 1)
	// For tileRow=1, zdither0=1: iseed = (1+1-2) % 10000 = 0.
	// nextrand starts at int(randoms[0] * 500).
	iseed := 0
	nextrand := int(randoms[iseed] * 500)
	for i := range input {
		want := (float64(100) - float64(randoms[nextrand]) + 0.5) * 1.0
		if math.Abs(out[i]-want) > 1e-12 {
			t.Errorf("pixel %d: got %v want %v", i, out[i], want)
		}
		nextrand++
	}
}

// TestDequantizeSubtractive2PreservesZero: a ZeroValue sentinel must map
// to exactly 0.0, regardless of the dither offset.
func TestDequantizeSubtractive2PreservesZero(t *testing.T) {
	input := []int32{50, ZeroValue, 100, ZeroValue, 200}
	out := make([]float64, len(input))
	Dequantize(input, out, 0.25, 10.0, SubtractiveDither2, 1, 1)
	if out[1] != 0.0 {
		t.Errorf("ZeroValue at idx 1 → %v, expected 0.0", out[1])
	}
	if out[3] != 0.0 {
		t.Errorf("ZeroValue at idx 3 → %v, expected 0.0", out[3])
	}
	// Non-zero pixels should still be dithered.
	if out[0] == 50*0.25+10.0 {
		t.Error("idx 0 should have been dithered, not plain-scaled")
	}
}

// TestParseDitherMethod covers the ZQUANTIZ string mapping.
func TestParseDitherMethod(t *testing.T) {
	cases := []struct {
		s string
		m DitherMethod
	}{
		{"", NoDither},
		{"NO_DITHER", NoDither},
		{"SUBTRACTIVE_DITHER_1", SubtractiveDither1},
		{"SUBTRACTIVE_DITHER_2", SubtractiveDither2},
		{"something else", NoDither},
	}
	for _, c := range cases {
		if got := ParseDitherMethod(c.s); got != c.m {
			t.Errorf("ParseDitherMethod(%q) = %v, want %v", c.s, got, c.m)
		}
	}
}
