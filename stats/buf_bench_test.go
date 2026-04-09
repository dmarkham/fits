package stats

import (
	"math/rand"
	"testing"
)

// Benchmarks comparing the allocating Median/MAD/Percentile/MeanStdev
// API to the zero-alloc *Buf variants and to QuickSelect/MeanStdevSiril.
//
// Two sizes:
//
//   - small (n=50): the rejection-step use case from
//     fits-processing/stack/reject.go — millions of calls per run on
//     short slices. Per-call alloc thrashes GC even though each alloc
//     is tiny.
//
//   - large (n=500_000): the overlap_norm use case — fewer calls (~10K
//     per run) but each allocation is hundreds of KB (histogram of
//     65536 uint32s + absDev of 500K float64s = ~4.25 MB per call).
//
// Run all of them:
//   go test -bench=. -benchmem ./stats/

const (
	benchSmallN = 50
	benchLargeN = 500_000
)

func benchData(n int, seed int64) []float32 {
	rng := rand.New(rand.NewSource(seed))
	out := make([]float32, n)
	for i := range out {
		out[i] = rng.Float32()
	}
	return out
}

// ----------------------------------------------------------------------
// Median: alloc vs Buf
// ----------------------------------------------------------------------

func BenchmarkMedian_Small(b *testing.B) {
	data := benchData(benchSmallN, 1)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = Median(data)
	}
}

func BenchmarkMedianBuf_Small(b *testing.B) {
	data := benchData(benchSmallN, 1)
	histo := make([]uint32, MaxHistoSize)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = MedianBuf(data, histo)
	}
}

func BenchmarkMedian_Large(b *testing.B) {
	data := benchData(benchLargeN, 1)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = Median(data)
	}
}

func BenchmarkMedianBuf_Large(b *testing.B) {
	data := benchData(benchLargeN, 1)
	histo := make([]uint32, MaxHistoSize)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = MedianBuf(data, histo)
	}
}

// ----------------------------------------------------------------------
// MAD: alloc vs Buf
// ----------------------------------------------------------------------

func BenchmarkMAD_Small(b *testing.B) {
	data := benchData(benchSmallN, 2)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = MAD(data)
	}
}

func BenchmarkMADWithMedianBuf_Small(b *testing.B) {
	data := benchData(benchSmallN, 2)
	histo := make([]uint32, MaxHistoSize)
	absDev := make([]float64, 0, benchSmallN)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_, _ = MADWithMedianBuf(data, histo, absDev)
	}
}

func BenchmarkMAD_Large(b *testing.B) {
	data := benchData(benchLargeN, 2)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = MAD(data)
	}
}

func BenchmarkMADWithMedianBuf_Large(b *testing.B) {
	data := benchData(benchLargeN, 2)
	histo := make([]uint32, MaxHistoSize)
	absDev := make([]float64, 0, benchLargeN)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_, _ = MADWithMedianBuf(data, histo, absDev)
	}
}

// ----------------------------------------------------------------------
// Percentile: alloc vs Buf (one mid-range size only — covered above
// for the median case, but worth a separate point for non-0.5 p)
// ----------------------------------------------------------------------

func BenchmarkPercentile_p99_Small(b *testing.B) {
	data := benchData(benchSmallN, 3)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = Percentile(data, 0.99)
	}
}

func BenchmarkPercentileBuf_p99_Small(b *testing.B) {
	data := benchData(benchSmallN, 3)
	histo := make([]uint32, MaxHistoSize)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = PercentileBuf(data, 0.99, histo)
	}
}

// ----------------------------------------------------------------------
// MeanStdev: existing generic vs Siril variant
// ----------------------------------------------------------------------

func BenchmarkMeanStdev_Small(b *testing.B) {
	data := benchData(benchSmallN, 4)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_, _ = MeanStdev(data)
	}
}

func BenchmarkMeanStdevSiril_Small(b *testing.B) {
	data := benchData(benchSmallN, 4)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_, _ = MeanStdevSiril(data)
	}
}

func BenchmarkMeanStdev_Large(b *testing.B) {
	data := benchData(benchLargeN, 4)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_, _ = MeanStdev(data)
	}
}

func BenchmarkMeanStdevSiril_Large(b *testing.B) {
	data := benchData(benchLargeN, 4)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_, _ = MeanStdevSiril(data)
	}
}

// ----------------------------------------------------------------------
// QuickSelect: exact-element median for small n (the reject.go use case)
// ----------------------------------------------------------------------

func BenchmarkQuickSelect_Small(b *testing.B) {
	original := benchData(benchSmallN, 5)
	work := make([]float32, len(original))
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		copy(work, original)
		_ = QuickSelect(work, len(work)/2)
	}
}

func BenchmarkMedianBuf_SmallVsQuickSelect(b *testing.B) {
	// Apples-to-apples comparison at the rejection-step size: how
	// does the histogram-based MedianBuf compare to QuickSelect for
	// n=50? For small n the histogram has very few bins and
	// quickselect is usually faster — this benchmark is the evidence.
	data := benchData(benchSmallN, 5)
	histo := make([]uint32, MaxHistoSize)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		_ = MedianBuf(data, histo)
	}
}
