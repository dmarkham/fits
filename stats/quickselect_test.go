package stats

import (
	"math/rand"
	"slices"
	"testing"
)

// TestQuickSelect_AgainstSort cross-checks against sort.Slice for
// every k across a variety of inputs. The reference impl is "sort and
// index" — slow but obviously correct.
func TestQuickSelect_AgainstSort(t *testing.T) {
	cases := map[string][]float32{
		"distinct_small":  {3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5},
		"already_sorted": {1, 2, 3, 4, 5, 6, 7, 8, 9, 10},
		"reverse_sorted": {10, 9, 8, 7, 6, 5, 4, 3, 2, 1},
		"all_same":       {7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
		"two_values":     {1, 2, 1, 2, 1, 2, 1, 2, 1, 2},
		"single":         {42},
		"pair":           {2, 1},
		"negatives":      {-3, -1, -4, 1, -5, 9, -2, 6, -5, 3, 5},
		"random":         genRandomFloat32(100, 42),
		"large_random":   genRandomFloat32(10000, 7),
	}

	for name, base := range cases {
		t.Run(name, func(t *testing.T) {
			// Sort a copy as the reference.
			ref := slices.Clone(base)
			slices.Sort(ref)

			// QuickSelect for every k.
			for k := 0; k < len(base); k++ {
				work := slices.Clone(base)
				got := QuickSelect(work, k)
				want := ref[k]
				if got != want {
					t.Errorf("k=%d: got %v want %v (input %v)",
						k, got, want, base)
				}
			}
		})
	}
}

// TestQuickSelect_PartitionInvariant verifies the post-condition: after
// QuickSelect(data, k), data[k] is the k-th smallest, and every
// element before k is <= and every element after k is >=. This is
// what makes QuickSelect useful as a partial sort.
func TestQuickSelect_PartitionInvariant(t *testing.T) {
	rng := rand.New(rand.NewSource(42))
	for trial := 0; trial < 50; trial++ {
		n := 50 + rng.Intn(200)
		data := make([]float32, n)
		for i := range data {
			data[i] = rng.Float32()
		}
		k := rng.Intn(n)

		QuickSelect(data, k)

		pivot := data[k]
		for i := 0; i < k; i++ {
			if data[i] > pivot {
				t.Errorf("trial %d k=%d: data[%d]=%v > pivot=%v",
					trial, k, i, data[i], pivot)
			}
		}
		for i := k + 1; i < n; i++ {
			if data[i] < pivot {
				t.Errorf("trial %d k=%d: data[%d]=%v < pivot=%v",
					trial, k, i, data[i], pivot)
			}
		}
	}
}

// TestQuickSelect_EdgeCases covers k=0 (min), k=n-1 (max), and
// out-of-bounds panics.
func TestQuickSelect_EdgeCases(t *testing.T) {
	t.Run("k0_is_min", func(t *testing.T) {
		data := []float32{5, 3, 8, 1, 9, 2, 7, 4, 6}
		got := QuickSelect(data, 0)
		if got != 1 {
			t.Errorf("k=0: got %v want 1", got)
		}
	})

	t.Run("k_n-1_is_max", func(t *testing.T) {
		data := []float32{5, 3, 8, 1, 9, 2, 7, 4, 6}
		got := QuickSelect(data, len(data)-1)
		if got != 9 {
			t.Errorf("k=n-1: got %v want 9", got)
		}
	})

	t.Run("empty_panics", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Fatal("expected panic on empty input")
			}
		}()
		var data []float32
		_ = QuickSelect(data, 0)
	})

	t.Run("k_negative_panics", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Fatal("expected panic on k<0")
			}
		}()
		_ = QuickSelect([]float32{1, 2, 3}, -1)
	})

	t.Run("k_too_large_panics", func(t *testing.T) {
		defer func() {
			if r := recover(); r == nil {
				t.Fatal("expected panic on k>=n")
			}
		}()
		_ = QuickSelect([]float32{1, 2, 3}, 3)
	})
}

// TestQuickSelect_IntegerTypes confirms the generic constraint works
// for the integer Numeric types — important for the planned use case
// where QuickSelect on int16 image pixels gives an exact-element
// median that the histogram method can't.
func TestQuickSelect_IntegerTypes(t *testing.T) {
	t.Run("int16", func(t *testing.T) {
		data := []int16{500, 100, 300, 200, 400}
		got := QuickSelect(data, 2)
		if got != 300 {
			t.Errorf("int16 median: got %v want 300", got)
		}
	})
	t.Run("uint8", func(t *testing.T) {
		data := []uint8{50, 10, 30, 20, 40}
		got := QuickSelect(data, 2)
		if got != 30 {
			t.Errorf("uint8 median: got %v want 30", got)
		}
	})
}

// TestQuickSelect_ZeroAlloc — in-place algorithm must not allocate.
func TestQuickSelect_ZeroAlloc(t *testing.T) {
	original := genRandomFloat32(1000, 99)
	work := make([]float32, len(original))

	// Warm up.
	copy(work, original)
	_ = QuickSelect(work, 500)

	allocs := testing.AllocsPerRun(100, func() {
		copy(work, original)
		_ = QuickSelect(work, 500)
	})
	// AllocsPerRun reports the avg over multiple iterations. Expect
	// 0 — the copy is on a pre-allocated buffer, and QuickSelect is
	// fully in place.
	if allocs > 0 {
		t.Errorf("QuickSelect allocates %v times per call; expected 0", allocs)
	}
}

// genRandomFloat32 returns a deterministic random slice. Used so the
// quickselect tests are reproducible.
func genRandomFloat32(n int, seed int64) []float32 {
	rng := rand.New(rand.NewSource(seed))
	out := make([]float32, n)
	for i := range out {
		out[i] = rng.Float32()
	}
	return out
}
