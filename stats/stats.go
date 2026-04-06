// Package stats provides generic, NaN-aware statistics functions for
// astronomical image data. Every function operates on slices of any
// Numeric type and uses float64 accumulators internally for precision.
//
// The core algorithms are extracted from the library's existing
// cross-validated implementations:
//
//   - Percentile/Median/MAD: ported from Siril/RawTherapee's
//     findMinMaxPercentile (histogram-based interpolated percentile)
//   - MeanStdev: ported from cfitsio's FnMeanSigma
//
// NaN handling: for float32/float64 inputs, NaN values are silently
// skipped in all computations. Integer types cannot be NaN, so the
// generic code fast-paths them with no NaN checks.
package stats

import "math"

// Numeric is the type constraint for all stats functions, matching the
// fits package's pixel type set.
type Numeric interface {
	~uint8 | ~int8 | ~int16 | ~uint16 | ~int32 | ~uint32 | ~int64 | ~uint64 | ~float32 | ~float64
}

// isNaN reports whether v is NaN. For integer types this is always false.
func isNaN[T Numeric](v T) bool {
	// The compiler should optimize this away for integer types.
	return v != v
}

// toFloat64 converts any Numeric value to float64.
func toFloat64[T Numeric](v T) float64 {
	return float64(v)
}

// MinMax returns the minimum and maximum of data, skipping NaN.
// For empty slices, returns (0, 0).
func MinMax[T Numeric](data []T) (min, max T) {
	started := false
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		if !started {
			min, max = v, v
			started = true
			continue
		}
		if v < min {
			min = v
		}
		if v > max {
			max = v
		}
	}
	return min, max
}

// Mean returns the arithmetic mean of data, skipping NaN.
// For empty slices (or all-NaN), returns 0.
func Mean[T Numeric](data []T) float64 {
	var sum float64
	var n int64
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		sum += toFloat64(v)
		n++
	}
	if n == 0 {
		return 0
	}
	return sum / float64(n)
}

// MeanStdev returns the arithmetic mean and population standard
// deviation of data, skipping NaN. Uses float64 accumulators.
// For empty slices, returns (0, 0).
func MeanStdev[T Numeric](data []T) (mean, stdev float64) {
	var sum, sum2 float64
	var n int64
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		x := toFloat64(v)
		sum += x
		sum2 += x * x
		n++
	}
	if n == 0 {
		return 0, 0
	}
	mean = sum / float64(n)
	// One-pass variance (cfitsio convention). Can lose precision due to
	// catastrophic cancellation when values cluster far from zero;
	// the variance < 0 guard handles the worst case.
	variance := sum2/float64(n) - mean*mean
	if variance < 0 {
		variance = 0
	}
	stdev = math.Sqrt(variance)
	return mean, stdev
}
