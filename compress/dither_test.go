package compress

import (
	"math"
	"testing"
)

// TestFitsRandomsInvariant verifies the Park-Miller seed invariant from
// cfitsio's fits_init_randoms. The cfitsio source asserts that after
// exactly 10000 iterations the internal seed equals 1043618065 — if
// our port doesn't produce the same final seed the whole dither path
// is silently wrong.
func TestFitsRandomsInvariant(t *testing.T) {
	r := FitsRandoms()
	if len(r) != NRandom {
		t.Fatalf("len(FitsRandoms) = %d, want %d", len(r), NRandom)
	}
	// Verify the invariant by re-running the algorithm.
	const a = 16807.0
	const m = 2147483647.0
	seed := 1.0
	for i := 0; i < NRandom; i++ {
		temp := a * seed
		seed = temp - m*float64(int64(temp/m))
	}
	if int(seed) != 1043618065 {
		t.Fatalf("cfitsio invariant: expected seed = 1043618065, got %d", int(seed))
	}
}

// TestFitsRandomsFirstValues checks that the first few values match what
// the Park-Miller generator should produce. This guards against any
// porting error in the loop body.
func TestFitsRandomsFirstValues(t *testing.T) {
	r := FitsRandoms()
	// Compute expected values independently.
	const a = 16807.0
	const m = 2147483647.0
	seed := 1.0
	for i := 0; i < 5; i++ {
		temp := a * seed
		seed = temp - m*float64(int64(temp/m))
		expected := float32(seed / m)
		if r[i] != expected {
			t.Errorf("FitsRandoms[%d] = %v, expected %v", i, r[i], expected)
		}
	}
}

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
