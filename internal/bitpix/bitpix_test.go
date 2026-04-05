package bitpix

import "testing"

func TestSize(t *testing.T) {
	cases := []struct {
		b    BITPIX
		want int
	}{
		{Int8, 1},
		{Int16, 2},
		{Int32, 4},
		{Int64, 8},
		{Float32, 4},
		{Float64, 8},
	}
	for _, c := range cases {
		if got := c.b.Size(); got != c.want {
			t.Errorf("%s Size = %d, want %d", c.b, got, c.want)
		}
	}
}

func TestValid(t *testing.T) {
	for _, b := range []BITPIX{Int8, Int16, Int32, Int64, Float32, Float64} {
		if !b.Valid() {
			t.Errorf("%s should be valid", b)
		}
	}
	for _, b := range []BITPIX{0, 1, 7, 9, 17, -1, -33} {
		if b.Valid() {
			t.Errorf("%d should be invalid", int(b))
		}
	}
}

func TestIsFloatInt(t *testing.T) {
	if !Float32.IsFloat() || !Float64.IsFloat() {
		t.Errorf("floats")
	}
	if Int8.IsFloat() || Int16.IsFloat() || Int32.IsFloat() || Int64.IsFloat() {
		t.Errorf("ints marked float")
	}
	if !Int8.IsInt() || !Int16.IsInt() || !Int32.IsInt() || !Int64.IsInt() {
		t.Errorf("ints")
	}
}
