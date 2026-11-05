// Package stats provides NaN-aware statistics functions for astronomical
// image data. Most functions are generic over the Numeric element types;
// MeanStdevSiril takes []float32 only.
//
// Each function documents what it computes and what it is tested against.
// In summary (truth sources live under testdata/, with their generators):
//
//   - Percentile, Median, MAD: RawTherapee/Siril's histogram-based
//     findMinMaxPercentile and siril_stats_float_mad. Bit-identical to
//     Siril's C code for []float32 (verbatim Siril source compiled in
//     testdata/sref); other types run the same algorithm in float64 and
//     are bit-identical to that source compiled at double precision.
//   - MeanStdevSiril: bit-identical to Siril's siril_stats_float_sd.
//   - MeanStdev: one-pass float64 sample mean/stdev (N-1). Agrees with
//     numpy std(ddof=1) within the rounding-error bound of the one-pass
//     formula, which is large when the data sit far from zero.
//   - MeanStdevTwoPass: compensated two-pass sample mean/stdev; on every
//     tested input the mean is the correctly rounded exact mean and the
//     stdev is within 2 ULP of the exact (rational-arithmetic) value.
//   - SigmaClip: astropy.stats.sigma_clip semantics (asymmetric bounds,
//     maxiters=None), with Median as the median centre; the same survivors
//     and iteration count as astropy given the exact mean and std.
//   - BuildHistogram: identical bins and edges to numpy.histogram.
//   - QuickSelect: exact order statistics (checked against a sort and
//     Siril's quickmedian_float).
//
// NaN handling: unless a function says otherwise, NaN values are skipped
// (integer types cannot be NaN). ±Inf handling is documented per
// function.
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
// For empty (or all-NaN) slices, returns (0, 0). ±Inf values take part
// like any other value.
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

// Mean returns the arithmetic mean of data, skipping NaN, as a plain
// float64 running sum divided by the count. Its rounding error is at
// most about n·ε·mean(|x|) (ε = 2^-53), the standard bound for recursive
// summation. ±Inf propagate by IEEE rules. For empty slices (or all-NaN),
// returns 0.
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

// MeanStdev returns the arithmetic mean and the sample standard deviation
// (Bessel-corrected, N-1 denominator) of data, skipping NaN, using the
// one-pass formula
//
//	variance = (Σx² - n·mean²) / (n-1)
//
// with float64 accumulators (a negative variance from rounding is clamped
// to 0). It agrees with numpy's std(x, ddof=1) to within the rounding
// error of that formula: the absolute error of the variance is up to
// about n·ε·mean(x²) (ε = 2^-53). That is negligible when the spread is
// comparable to the values, but the subtraction cancels catastrophically
// when the data sit far from zero relative to their spread: for 1000
// values alternating 1e9 and 1e9+1 the true sample stdev is 0.50025 and
// MeanStdev returns 95.834358265834 on amd64, or 96.153395548339716 on
// arm64. Use MeanStdevTwoPass for such data, and MeanStdevSiril to
// reproduce Siril.
//
// The low bits of the result depend on the build: the code has no
// explicit float64() conversions, so the Go compiler fuses its
// multiply-adds where the target has FMA instructions (nothing on amd64
// below GOAMD64=v3; the Σx² accumulation on GOAMD64=v3; the accumulation
// and the final subtraction on arm64). Each case is pinned bit for bit
// against a reproduction of that arithmetic (stats/testdata/gen_golden.py,
// TestMeanStdev_OnePassFormula).
//
// For empty (or all-NaN) input it returns (0, 0); for a single value,
// (value, 0). ±Inf propagate by IEEE rules (an Inf value gives a NaN
// stdev).
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
	if n < 2 {
		return mean, 0
	}
	// One-pass sample variance: (sum2 - n*mean^2) / (n-1). Loses
	// precision via catastrophic cancellation when values cluster far
	// from zero (see the doc comment); the variance < 0 guard only
	// stops a negative result.
	variance := (sum2 - float64(n)*mean*mean) / float64(n-1)
	if variance < 0 {
		variance = 0
	}
	stdev = math.Sqrt(variance)
	return mean, stdev
}

// MeanStdevTwoPass returns the arithmetic mean and the sample standard
// deviation (N-1 denominator) of data, skipping NaN, without the
// cancellation problem of MeanStdev. Values are converted to float64
// first, so int64 and uint64 values beyond 2^53 are rounded on the way in.
//
// It is the corrected two-pass algorithm (Björck; Chan, Golub & LeVeque):
// pass one computes the mean from a compensated (Neumaier) sum, divided
// by n with an exact-remainder correction; pass two accumulates
// d = x - mean as both Σd² and Σd with compensated sums and returns
// (Σd² - (Σd)²/n) / (n-1), where the (Σd)²/n term removes the error left
// by the rounded mean.
//
// The mean is the correctly rounded quotient of the compensated sum by n,
// except when that quotient lies within about 2^-50 ULP of a rounding
// boundary (the correction term is itself rounded once). The compensated
// sum is the exact sum unless the compensation rounds; its error is at
// most about 2ε|S| + n·ε²·Σ|x| (ε = 2^-53), so only sums that cancel to
// far below Σ|x|/(n·ε) can be inexact.
//
// Accuracy, tested against exact arithmetic (testdata/gen_golden.py): the
// mean is bit-identical to the correctly rounded exact mean on every
// golden input (including data offset by 1e9 and float32 ADU data offset
// by 3e4) and on 20000 sets of integers of magnitude 2^40 to 2^62, as
// []int64 and []float64; the stdev is within 2 ULP of the exact value.
//
// Conventions match MeanStdev: (0, 0) for empty or all-NaN input and
// (value, 0) for a single value. Non-finite cases follow IEEE arithmetic
// and give what numpy's mean and std(ddof=1) give (tested): ±Inf in the
// data makes the stdev NaN; a sum or Σd² that overflows float64 gives
// ±Inf instead of the exact value.
func MeanStdevTwoPass[T Numeric](data []T) (mean, stdev float64) {
	n, mean, m2 := twoPassMoments(data)
	if n == 0 {
		return 0, 0
	}
	if n < 2 {
		return mean, 0
	}
	variance := m2 / float64(n-1)
	if variance < 0 {
		variance = 0
	}
	return mean, math.Sqrt(variance)
}

// twoPassMoments returns the count of non-NaN values, their mean and the
// corrected sum of squared deviations Σd² - (Σd)²/n (d = x - mean), with
// the compensated sums of MeanStdevTwoPass. For n == 0 the mean is NaN.
func twoPassMoments[T Numeric](data []T) (n int64, mean, m2 float64) {
	var s neumaier
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		s.add(toFloat64(v))
		n++
	}
	if n == 0 {
		return 0, math.NaN(), math.NaN()
	}
	mean = s.mean(float64(n))
	var sd, sd2 neumaier
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		d := toFloat64(v) - mean
		sd.add(d)
		sd2.add(float64(d * d)) // explicit rounding: no fused multiply-add
	}
	m2 = sd2.sum()
	if corr := sd.sum(); !math.IsInf(corr, 0) && !math.IsNaN(corr) {
		m2 -= corr * corr / float64(n)
	}
	return n, mean, m2
}

// neumaier is a compensated (Kahan-Babuška-Neumaier) float64 sum. Once
// the running sum is no longer finite the compensation is meaningless
// (Inf-Inf), so sum then reports the plain IEEE sum.
type neumaier struct{ s, c float64 }

func (k *neumaier) add(x float64) {
	t := k.s + x
	if math.Abs(k.s) >= math.Abs(x) {
		k.c += (k.s - t) + x
	} else {
		k.c += (x - t) + k.s
	}
	k.s = t
}

// mean returns sum/n rounded once (to within a tiny fraction of an ULP
// of the correctly rounded quotient): the running sum and its
// compensation are renormalised into a double-double t + e, and the
// remainder of t/n, computed exactly with a fused multiply-add, corrects
// the quotient. A rounding of sum() followed by a division could be one
// ULP off, which matters when the data are quantised at the mean's ULP
// (integers beyond 2^53): a sigma-clip bound then moves by a whole value.
func (k *neumaier) mean(n float64) float64 {
	if math.IsInf(k.s, 0) || math.IsNaN(k.s) {
		return k.s / n
	}
	t := k.s + k.c
	if math.IsInf(t, 0) {
		return t / n
	}
	// Two-sum: t + e == k.s + k.c exactly.
	bb := t - k.s
	e := (k.s - (t - bb)) + (k.c - bb)
	q := t / n
	r := math.FMA(-q, n, t) + e
	return q + r/n
}

func (k *neumaier) sum() float64 {
	if math.IsInf(k.s, 0) || math.IsNaN(k.s) {
		return k.s
	}
	return k.s + k.c
}
