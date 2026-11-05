package stats

import (
	"math"
	"testing"
)

// sirilSDNoVarRound is siril_stats_float_sd with one change: the variance
// is not rounded to float32 before the square root (only the result is).
// It exists only to pin the MeanStdevSiril doc's claim that this rounding
// is load-bearing; it produces no expected value.
func sirilSDNoVarRound(data []float32) float32 {
	n := len(data)
	var acc float64
	for _, v := range data {
		acc += float64(v)
	}
	mean := float32(acc / float64(n))
	acc = 0
	for _, v := range data {
		d := v - mean
		acc += float64(d * d)
	}
	return float32(math.Sqrt(acc / float64(n-1)))
}

// sirilSDAllFloat64 is the textbook two-pass sample stdev entirely in
// float64 (mean, deviations, variance), rounded to float32 only at the end.
func sirilSDAllFloat64(data []float32) float32 {
	n := len(data)
	var acc float64
	for _, v := range data {
		acc += float64(v)
	}
	mean := acc / float64(n)
	acc = 0
	for _, v := range data {
		d := float64(v) - mean
		acc += float64(d * d)
	}
	return float32(math.Sqrt(acc / float64(n-1)))
}

// TestMeanStdevSiril_DocCounts pins the numbers in the MeanStdevSiril doc:
// over the 424 Siril golden inputs, skipping the float32 rounding of the
// variance differs from Siril's stdev on 43 of them, and doing everything
// in float64 on 50. Siril's value is the golden (the verbatim C code).
func TestMeanStdevSiril_DocCounts(t *testing.T) {
	g := loadSirilGolden(t)
	noVarRound, allF64 := 0, 0
	for _, c := range g.F32 {
		raw := c.data(t)
		want := hex32(t, c.SD.Stdev)
		if !same32(sirilSDNoVarRound(raw), want) {
			noVarRound++
		}
		if !same32(sirilSDAllFloat64(raw), want) {
			allF64++
		}
	}
	t.Logf("differ from Siril: variance not rounded %d, all float64 %d (of %d)", noVarRound, allF64, len(g.F32))
	if noVarRound != 43 || allF64 != 50 || len(g.F32) != 424 {
		t.Errorf("variance not rounded differs on %d, all-float64 on %d, of %d cases; the MeanStdevSiril doc says 43 and 50 of 424",
			noVarRound, allF64, len(g.F32))
	}
}
