package stats

import (
	"math"
	"testing"
)

// TestPercentileBuf_MatchesPercentile is the equivalence gate: for
// every fixture in the existing golden suite, PercentileBuf with a
// caller-supplied histogram must return bit-for-bit the same value
// as Percentile. If they ever diverge, the refactor has drifted.
func TestPercentileBuf_MatchesPercentile(t *testing.T) {
	cases := loadGoldens(t)
	histo := make([]uint32, MaxHistoSize)
	for name, c := range cases {
		if c.Len == 0 {
			continue
		}
		t.Run(name, func(t *testing.T) {
			data := toF32Slice(c.Input)
			for _, p := range []float64{0.0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0} {
				want := Percentile(data, p)
				got := PercentileBuf(data, p, histo)
				if got != want {
					t.Errorf("p=%g: PercentileBuf=%g != Percentile=%g", p, got, want)
				}
			}
		})
	}
}

// TestPercentileBuf_DirtyBuffer verifies the function zeroes the
// caller's histogram before use, so dirty (non-zero) input doesn't
// corrupt the result. Without this guarantee, callers would have to
// reset the buffer manually between calls — fragile and easy to forget.
func TestPercentileBuf_DirtyBuffer(t *testing.T) {
	data := []float32{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
	histo := make([]uint32, MaxHistoSize)

	// Pollute the buffer with garbage that would corrupt a percentile
	// walk if not zeroed (large counts at low bins shift the CDF).
	for i := range histo {
		histo[i] = 999999
	}

	want := Percentile(data, 0.5)
	got := PercentileBuf(data, 0.5, histo)
	if got != want {
		t.Errorf("dirty buffer: got %g want %g (buffer was not zeroed)", got, want)
	}
}

// TestPercentileBuf_PanicsOnUndersizedBuffer verifies the contract
// is enforced loudly. Silently truncating a too-small buffer would
// produce wrong answers; we panic instead so callers find the bug.
func TestPercentileBuf_PanicsOnUndersizedBuffer(t *testing.T) {
	data := make([]float32, 100)
	for i := range data {
		data[i] = float32(i)
	}
	// Need at least min(100, MaxHistoSize)=100; pass 50.
	histo := make([]uint32, 50)

	defer func() {
		r := recover()
		if r == nil {
			t.Fatalf("expected panic for undersized buffer; got none")
		}
		// Sanity check: error message should mention sizes.
		msg, ok := r.(string)
		if !ok {
			if e, ok := r.(error); ok {
				msg = e.Error()
			}
		}
		_ = msg // not asserting exact text — just confirming it panicked
	}()

	_ = PercentileBuf(data, 0.5, histo)
}

// TestPercentileBuf_ZeroAlloc verifies the whole point of the API:
// no allocations in steady state when reusing the histo buffer.
func TestPercentileBuf_ZeroAlloc(t *testing.T) {
	data := make([]float32, 1000)
	for i := range data {
		data[i] = float32(i) * 0.001
	}
	histo := make([]uint32, MaxHistoSize)

	// Warm up — first call may have generic-instantiation cost.
	_ = PercentileBuf(data, 0.5, histo)

	allocs := testing.AllocsPerRun(100, func() {
		_ = PercentileBuf(data, 0.5, histo)
	})
	if allocs > 0 {
		t.Errorf("PercentileBuf allocates %v times per call; expected 0", allocs)
	}
}

// TestPercentileBuf_MinSizedBuffer confirms the boundary: a buffer
// of exactly min(n, MaxHistoSize) is accepted (off-by-one safety).
func TestPercentileBuf_MinSizedBuffer(t *testing.T) {
	data := make([]float32, 50)
	for i := range data {
		data[i] = float32(i) * 0.02
	}
	// histoSize will be 50 (n < MaxHistoSize). A buffer of exactly 50
	// should work — the >= check must not be a > check.
	histo := make([]uint32, 50)
	got := PercentileBuf(data, 0.5, histo)
	want := Percentile(data, 0.5)
	if got != want {
		t.Errorf("got %g want %g", got, want)
	}
}

// TestPercentileBuf_Float64Path exercises the generic path (non-float32
// input) to make sure the buf refactor for percentileGenericBuf is
// also bit-equivalent.
func TestPercentileBuf_Float64Path(t *testing.T) {
	data := []float64{0.1, 0.5, 0.3, 0.9, 0.2, 0.7, 0.4, 0.6, 0.8, 1.0}
	histo := make([]uint32, MaxHistoSize)
	for _, p := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
		want := Percentile(data, p)
		got := PercentileBuf(data, p, histo)
		if got != want {
			t.Errorf("p=%g: got %g want %g", p, got, want)
		}
	}
}

// TestPercentileBuf_AllNaN confirms the n=0-after-NaN-skip edge case
// returns 0 without touching histo (so a nil histo would be safe in
// principle — but we still take the slice for API consistency).
func TestPercentileBuf_AllNaN(t *testing.T) {
	data := []float32{float32(math.NaN()), float32(math.NaN())}
	histo := make([]uint32, MaxHistoSize)
	got := PercentileBuf(data, 0.5, histo)
	if got != 0 {
		t.Errorf("all-NaN: got %g want 0", got)
	}
}

// TestMedianBuf_MatchesMedian: MedianBuf must equal Median bit-for-bit
// across the existing fixture suite.
func TestMedianBuf_MatchesMedian(t *testing.T) {
	cases := loadGoldens(t)
	histo := make([]uint32, MaxHistoSize)
	for name, c := range cases {
		if c.Len == 0 {
			continue
		}
		t.Run(name, func(t *testing.T) {
			data := toF32Slice(c.Input)
			want := Median(data)
			got := MedianBuf(data, histo)
			if got != want {
				t.Errorf("MedianBuf=%g != Median=%g", got, want)
			}
		})
	}
}

// TestMedianBuf_ZeroAlloc — the whole point.
func TestMedianBuf_ZeroAlloc(t *testing.T) {
	data := make([]float32, 1000)
	for i := range data {
		data[i] = float32(i) * 0.001
	}
	histo := make([]uint32, MaxHistoSize)

	_ = MedianBuf(data, histo) // warm up

	allocs := testing.AllocsPerRun(100, func() {
		_ = MedianBuf(data, histo)
	})
	if allocs > 0 {
		t.Errorf("MedianBuf allocates %v times per call; expected 0", allocs)
	}
}

// TestMADWithMedianBuf_MatchesMADWithMedian: bit-for-bit equivalence
// gate against the existing MADWithMedian across the fixture suite.
func TestMADWithMedianBuf_MatchesMADWithMedian(t *testing.T) {
	cases := loadGoldens(t)
	histo := make([]uint32, MaxHistoSize)
	absDev := make([]float64, 0, 100000) // big enough for all fixtures
	for name, c := range cases {
		if c.Len == 0 {
			continue
		}
		t.Run(name, func(t *testing.T) {
			data := toF32Slice(c.Input)
			wantMed, wantMad := MADWithMedian(data)
			gotMed, gotMad := MADWithMedianBuf(data, histo, absDev)
			if gotMed != wantMed || gotMad != wantMad {
				t.Errorf("got (%g, %g) want (%g, %g)",
					gotMed, gotMad, wantMed, wantMad)
			}
		})
	}
}

// TestMADWithMedianBuf_ZeroAlloc verifies steady-state zero allocation
// with reused scratch buffers. This is the load-bearing test — if MAD
// allocates per call, the whole exercise was pointless.
func TestMADWithMedianBuf_ZeroAlloc(t *testing.T) {
	data := make([]float32, 1000)
	for i := range data {
		data[i] = float32(i) * 0.001
	}
	histo := make([]uint32, MaxHistoSize)
	absDev := make([]float64, 0, 1000)

	_, _ = MADWithMedianBuf(data, histo, absDev) // warm up

	allocs := testing.AllocsPerRun(100, func() {
		_, _ = MADWithMedianBuf(data, histo, absDev)
	})
	if allocs > 0 {
		t.Errorf("MADWithMedianBuf allocates %v times per call; expected 0", allocs)
	}
}

// TestMADWithMedianBuf_PanicsOnUndersizedAbsDev verifies the absDev
// contract is enforced.
func TestMADWithMedianBuf_PanicsOnUndersizedAbsDev(t *testing.T) {
	data := make([]float32, 100)
	for i := range data {
		data[i] = float32(i)
	}
	histo := make([]uint32, MaxHistoSize)
	absDev := make([]float64, 0, 50) // cap=50, need 100

	defer func() {
		if r := recover(); r == nil {
			t.Fatalf("expected panic for undersized absDev; got none")
		}
	}()
	_, _ = MADWithMedianBuf(data, histo, absDev)
}

// TestMADWithMedianBuf_DoesNotMutateInputAbsDev confirms the function
// resets absDev to len 0 internally — caller's len is irrelevant, only
// cap matters.
func TestMADWithMedianBuf_DoesNotMutateInputAbsDev(t *testing.T) {
	data := []float32{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}
	histo := make([]uint32, MaxHistoSize)
	// Pass a buffer with stale length and garbage values.
	absDev := make([]float64, 5, 100)
	for i := range absDev {
		absDev[i] = -999
	}

	wantMed, wantMad := MADWithMedian(data)
	gotMed, gotMad := MADWithMedianBuf(data, histo, absDev)
	if gotMed != wantMed || gotMad != wantMad {
		t.Errorf("stale absDev: got (%g, %g) want (%g, %g)",
			gotMed, gotMad, wantMed, wantMad)
	}
}
