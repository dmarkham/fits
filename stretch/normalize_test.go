package stretch

import (
	"fmt"
	"math"
	"testing"
)

// NormalizeChannels: expected values are worked out by hand from the
// documented rule (shared min/max over the finite pixels of all channels,
// NaN -> 0 when rescaling).
func TestNormalizeChannels(t *testing.T) {
	nan := float32(math.NaN())
	inf := float32(math.Inf(1))
	for _, tc := range []struct {
		name string
		in   [][]float32
		want [][]float32
	}{
		{"shared range", [][]float32{{0, 2}, {4, 8}}, [][]float32{{0, 0.25}, {0.5, 1}}},
		{"negative min", [][]float32{{-2, 0, 2}}, [][]float32{{0, 0.5, 1}}},
		{"NaN first", [][]float32{{nan, 2, 4}, {0, 1, 3}}, [][]float32{{0, 0.5, 1}, {0, 0.25, 0.75}}},
		{"NaN later", [][]float32{{-0.5, 0.5, nan}}, [][]float32{{0, 1, 0}}},
		{"empty first channel", [][]float32{{}, {1, 3}}, [][]float32{{}, {0, 1}}},
		{"all empty", [][]float32{{}, {}}, [][]float32{{}, {}}},
		{"already in [0,1]", [][]float32{{0.2, 0.5}, {nan, 0.9}}, [][]float32{{0.2, 0.5}, {nan, 0.9}}},
		{"constant", [][]float32{{5, 5}, {5}}, [][]float32{{5, 5}, {5}}},
		{"all NaN", [][]float32{{nan}, {nan, nan}}, [][]float32{{nan}, {nan, nan}}},
		{"one above 1", [][]float32{{0.5, 1, 2.5}}, [][]float32{{0, 0.25, 1}}},
		// The [0, 1] skip boundary, both sides (float32 results checked with
		// numpy float32: 1.5*(1/1.5f) = 1, 0.5*(1/1.0000001f) = 0.49999994).
		{"max in (1,2]", [][]float32{{0, 1.5}, {0.75}}, [][]float32{{0, 1}, {0.5}}},
		{"exactly [0,1]", [][]float32{{0, 1}, {0.5}}, [][]float32{{0, 1}, {0.5}}},
		{"one ULP above 1", [][]float32{{0, 0.5, math.Nextafter32(1, 2)}}, [][]float32{{0, 0.49999994, 1}}},
		{"min -0", [][]float32{{float32(math.Copysign(0, -1)), 1}}, [][]float32{{float32(math.Copysign(0, -1)), 1}}},
		{"min just below 0", [][]float32{{-math.SmallestNonzeroFloat32, 1}}, [][]float32{{0, 1}}},
		// ±Inf are not part of the range; when the data is rescaled +Inf
		// becomes 1 and -Inf 0. Before: [0 0 0 NaN] and [NaN NaN NaN NaN]
		// (rng = +Inf, inv = 0).
		{"+Inf", [][]float32{{0, 10, 20, inf}}, [][]float32{{0, 0.5, 1, 1}}},
		{"-Inf", [][]float32{{0, 10, 20, -inf}}, [][]float32{{0, 0.5, 1, 0}}},
		{"±Inf and NaN across channels", [][]float32{{inf, 2}, {-inf, 4, nan}}, [][]float32{{1, 0}, {0, 1, 0}}},
		{"Inf, finite range in [0,1]", [][]float32{{0.25, inf, 0.5}, {-inf}}, [][]float32{{0.25, inf, 0.5}, {-inf}}},
		{"Inf, constant finite", [][]float32{{3, inf, 3}}, [][]float32{{3, inf, 3}}},
		{"only Inf", [][]float32{{inf, -inf}}, [][]float32{{inf, -inf}}},
		// A finite range that overflows float32 (max - min > MaxFloat32) is
		// scaled in float64: float32((v - min) / (max - min)). Before:
		// [0 NaN 0] (rng = +Inf). Values from numpy:
		// np.float32((np.float64(v) - mn) / (np.float64(mx) - mn)).
		{"overflowing range", [][]float32{{-3e38, 3e38, 0}}, [][]float32{{0, 1, 0.5}}},
		{"overflowing range, 4 values", [][]float32{{-3.4e38, 3.4e38, 1e38, -1e38}}, [][]float32{{0, 1, 0.64705884, 0.3529412}}},
		// float32 multiply by the float32 reciprocal, not a division:
		// 1.25 * (1/3f) = 0.41666669, 1.25 / 3 = 0.41666666 (numpy float32).
		{"reciprocal, not division", [][]float32{{0, 1.25, 3}}, [][]float32{{0, 0.41666669, 1}}},
		{"overflowing range and +Inf", [][]float32{{3e38, -3e38}, {inf}}, [][]float32{{1, 0}, {1}}},
	} {
		t.Run(tc.name, func(t *testing.T) {
			NormalizeChannels(tc.in...)
			if fmt.Sprint(tc.in) != fmt.Sprint(tc.want) {
				t.Errorf("got %v, want %v", tc.in, tc.want)
			}
		})
	}
	NormalizeChannels() // no channels: no panic
}
