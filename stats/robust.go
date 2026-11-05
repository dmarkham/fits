package stats

import (
	"fmt"
	"math"
)

// MaxHistoSize is the upper bound on histogram bins used by the
// percentile/median family. Callers using the *Buf zero-alloc
// variants should pre-allocate a histogram slice of this length
// (or len(data), whichever is smaller).
const MaxHistoSize = 65536

// Percentile returns the p-th percentile (p in [0, 1], clamped) of data
// using RawTherapee/Siril's histogram-based findMinMaxPercentile
// (Siril src/rt/rt_algo.cc). The algorithm:
//
//  1. Find min/max of the non-NaN data
//  2. Build a histogram with min(65536, n) bins spanning [min, max]
//  3. Walk the cumulative count to p*n
//  4. Interpolate between the two bins straddling it
//  5. Convert back to the original value range and clamp to [min, max]
//
// It is not numpy's percentile. With n non-NaN values the result is
// within one histogram bin, (max-min)/(min(n, 65536)-1), of the order
// statistic x_(j), j = max(1, ceil(p*n)) (p*n rounded as the code rounds
// it), plus a few ULP of rounding in the final conversion
// (TestPercentile_WithinOneBin, against numpy order statistics). p = 0
// gives the minimum and p = 1 the maximum exactly.
//
// For []float32 input the computation is done in float32 exactly as
// Siril's C code does it, and the result is bit-identical to
// findMinMaxPercentile run on the same values with the NaNs removed
// (tested against the verbatim Siril source in testdata/sref, including
// inputs of 2^24+3 and 2*10^7 elements). For every other element type
// (including named types whose underlying type is float32) the values
// are converted to float64 and the same algorithm runs in float64; the
// result is bit-identical to that Siril source compiled with float
// replaced by double (tested for every element type, including int8,
// uint32, uint64 above 2^63 and a named float32 type, and on float64 and
// uint16 inputs of 2*10^7 elements).
//
// Edge cases: empty or all-NaN input returns 0; p is clamped to [0, 1]
// and a NaN p returns NaN; if the non-NaN data contains +Inf or -Inf, or
// max-min overflows the working precision (float32 for []float32, e.g.
// {-3e38, 3e38}; float64 otherwise), the histogram has no finite width
// and the result is NaN (Siril's code has undefined behaviour there: its
// bin index is a cast of NaN). TestPercentile_NonFinite and
// TestPercentile_RangeOverflow pin these.
//
// Allocates a histogram slice on every call. For hot paths, use
// PercentileBuf with a caller-supplied scratch buffer.
func Percentile[T Numeric](data []T, p float64) float64 {
	histoSize := len(data)
	if histoSize > MaxHistoSize {
		histoSize = MaxHistoSize
	}
	histo := make([]uint32, histoSize)
	return PercentileBuf(data, p, histo)
}

// PercentileBuf is the zero-alloc variant of Percentile. The histo
// slice must have len >= min(len(data), MaxHistoSize); panics
// otherwise. histo is zeroed and overwritten on each call, so callers
// can reuse the same buffer across many calls without resetting it.
//
// Typical usage: pre-allocate one MaxHistoSize histogram per goroutine
// and reuse it for every Percentile/Median/MAD call:
//
//	var histo = make([]uint32, stats.MaxHistoSize)
//	for _, frame := range frames {
//		med := stats.PercentileBuf(frame, 0.5, histo)
//		// ...
//	}
func PercentileBuf[T Numeric](data []T, p float64, histo []uint32) float64 {
	if p != p {
		return math.NaN()
	}
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	if f32, ok := any(data).([]float32); ok {
		return percentileFloat32Buf(f32, float32(p), histo)
	}
	return percentileGenericBuf(data, p, histo)
}

// percentileFloat32Buf is findMinMaxPercentile in float32, operation for
// operation. S is float64 only for the MAD path, whose deviations are
// float32 values stored in the caller's []float64 scratch buffer; every
// value converts to float32 exactly.
//
// Siril details that matter for bit-exactness:
//   - count is an integer (size_t) and is converted to float only where
//     the C code compares or subtracts it (a float32 accumulator drifts
//     once n exceeds 2^24);
//   - the bin index is static_cast<uint16_t> of a float product;
//   - every float32 product is rounded before it is added (explicit
//     float32 conversions forbid the compiler from fusing them).
func percentileFloat32Buf[S float32 | float64](data []S, p float32, histo []uint32) float64 {
	n := 0
	var minVal, maxVal float32
	for _, sv := range data {
		v := float32(sv)
		if v != v { // skip NaN
			continue
		}
		if n == 0 {
			minVal, maxVal = v, v
		} else {
			if v < minVal {
				minVal = v
			}
			if v > maxVal {
				maxVal = v
			}
		}
		n++
	}
	if n == 0 {
		return 0
	}
	if math.IsInf(float64(minVal), 0) || math.IsInf(float64(maxVal), 0) {
		return math.NaN()
	}
	if math.IsInf(float64(maxVal-minVal), 0) {
		// The float32 range overflows: the scale is 0 and Siril's bin
		// index would be a cast of NaN (undefined).
		return math.NaN()
	}
	if maxVal-minVal == 0 {
		return float64(minVal)
	}
	histoSize := n
	if histoSize > MaxHistoSize {
		histoSize = MaxHistoSize
	}
	if len(histo) < histoSize {
		panic(fmt.Sprintf("stats: histo buffer too small: got len %d, need %d", len(histo), histoSize))
	}
	histo = histo[:histoSize]
	for i := range histo {
		histo[i] = 0
	}
	scale := float32(histoSize-1) / (maxVal - minVal)
	for _, sv := range data {
		v := float32(sv)
		if v != v {
			continue
		}
		histo[uint16(float32(scale*(v-minVal)))]++
	}
	thresh := float32(p * float32(n))
	k := 0
	var count uint64
	for float32(count) < thresh && k < len(histo) {
		count += uint64(histo[k])
		k++
	}
	var result float32
	if k > 0 {
		countPrev := count - uint64(histo[k-1])
		c0 := float32(count) - thresh
		c1 := thresh - float32(countPrev)
		result = (float32(c1*float32(k)) + float32(c0*float32(k-1))) / (c0 + c1)
	} else {
		result = float32(k)
	}
	result /= scale
	result += minVal
	if result < minVal {
		result = minVal
	}
	if result > maxVal {
		result = maxVal
	}
	return float64(result)
}

// percentileGenericBuf is findMinMaxPercentile in float64.
func percentileGenericBuf[T Numeric](data []T, p float64, histo []uint32) float64 {
	n := 0
	var minVal, maxVal float64
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		fv := toFloat64(v)
		if n == 0 {
			minVal, maxVal = fv, fv
		} else {
			if fv < minVal {
				minVal = fv
			}
			if fv > maxVal {
				maxVal = fv
			}
		}
		n++
	}
	if n == 0 {
		return 0
	}
	if math.IsInf(minVal, 0) || math.IsInf(maxVal, 0) || math.IsInf(maxVal-minVal, 0) {
		return math.NaN()
	}
	if maxVal-minVal == 0 {
		return minVal
	}

	histoSize := n
	if histoSize > MaxHistoSize {
		histoSize = MaxHistoSize
	}
	if len(histo) < histoSize {
		panic(fmt.Sprintf("stats: histo buffer too small: got len %d, need %d", len(histo), histoSize))
	}
	histo = histo[:histoSize]
	for i := range histo {
		histo[i] = 0
	}
	scale := float64(histoSize-1) / (maxVal - minVal)

	for _, v := range data {
		if isNaN(v) {
			continue
		}
		bin := int(float64(scale * (toFloat64(v) - minVal)))
		if bin < 0 {
			bin = 0
		}
		if bin >= histoSize {
			bin = histoSize - 1
		}
		histo[bin]++
	}

	thresh := float64(p * float64(n))
	k := 0
	var count uint64
	for float64(count) < thresh && k < len(histo) {
		count += uint64(histo[k])
		k++
	}

	var result float64
	if k > 0 {
		countPrev := count - uint64(histo[k-1])
		c0 := float64(count) - thresh
		c1 := thresh - float64(countPrev)
		result = (float64(c1*float64(k)) + float64(c0*float64(k-1))) / (c0 + c1)
	} else {
		result = float64(k)
	}
	result /= scale
	result += minVal
	if result < minVal {
		result = minVal
	}
	if result > maxVal {
		result = maxVal
	}
	return result
}

// Median returns Percentile(data, 0.5): Siril's histogram_median_float
// for []float32 input (bit-identical), the float64 equivalent otherwise.
// It is a histogram approximation, not the exact middle element; use
// QuickSelect for an exact order statistic.
//
// Allocates a histogram on every call. For hot paths, use MedianBuf
// with a caller-supplied scratch buffer.
func Median[T Numeric](data []T) float64 {
	return Percentile(data, 0.5)
}

// MedianBuf is the zero-alloc variant of Median. The histo slice must
// have len >= min(len(data), MaxHistoSize); panics otherwise. histo is
// zeroed and overwritten on each call.
func MedianBuf[T Numeric](data []T, histo []uint32) float64 {
	return PercentileBuf(data, 0.5, histo)
}

// MAD returns the raw Median Absolute Deviation, median(|xi - median(x)|),
// with both medians computed by Median. Multiply by 1.4826 for a
// Gaussian-equivalent sigma.
//
// For []float32 input this is Siril's siril_stats_float_mad applied to
// the NaN-free data with the median from histogram_median_float: the
// deviations are fabsf(x - (float)median) in float32 and the result is
// bit-identical to Siril's (tested against the verbatim Siril source).
// Other element types use float64 deviations and the float64 median.
// NaN values are skipped; if the median is NaN (the data contains ±Inf,
// or max-min overflows; see Percentile) the MAD is NaN. Empty or all-NaN
// input returns 0.
//
// Allocates two slices on every call (histogram + absolute deviations).
// For hot paths, use MADWithMedianBuf with caller-supplied scratch.
func MAD[T Numeric](data []T) float64 {
	_, mad := MADWithMedian(data)
	return mad
}

// MADWithMedian returns the median and the MAD, avoiding double-
// computing the median when the caller needs both.
//
// Allocates two slices on every call. For hot paths, use
// MADWithMedianBuf.
func MADWithMedian[T Numeric](data []T) (median, mad float64) {
	histoSize := len(data)
	if histoSize > MaxHistoSize {
		histoSize = MaxHistoSize
	}
	histo := make([]uint32, histoSize)
	absDev := make([]float64, 0, len(data))
	return MADWithMedianBuf(data, histo, absDev)
}

// MADWithMedianBuf is the zero-alloc variant of MADWithMedian. Caller
// supplies two scratch slices:
//
//   - histo: must have len >= min(len(data), MaxHistoSize); used for
//     both the median pass and the MAD-of-deviations pass. Zeroed and
//     overwritten internally.
//   - absDev: must have cap >= len(data); used as scratch for the
//     |xi - median| array. Length is reset to 0 on entry, then grown
//     by append up to the non-NaN count of data.
//
// Panics if either buffer is too small.
//
// absDev is []float64 (not generic []T) so a single buffer works for
// any T. For []float32 input it holds the float32 deviations (exactly
// representable in float64), so the result is still Siril's float32
// computation.
func MADWithMedianBuf[T Numeric](data []T, histo []uint32, absDev []float64) (median, mad float64) {
	if cap(absDev) < len(data) {
		panic(fmt.Sprintf("stats: absDev buffer too small: got cap %d, need %d", cap(absDev), len(data)))
	}

	median = MedianBuf(data, histo)
	if median != median {
		return median, median
	}

	absDev = absDev[:0]
	if f32, ok := any(data).([]float32); ok {
		// siril_stats_float_mad: tmp[i] = fabsf(data[i] - median) with
		// median converted to float.
		m := float32(median)
		for _, v := range f32 {
			if v != v {
				continue
			}
			absDev = append(absDev, math.Abs(float64(v-m)))
		}
		if len(absDev) == 0 {
			return median, 0
		}
		return median, percentileFloat32Buf(absDev, 0.5, histo)
	}

	for _, v := range data {
		if isNaN(v) {
			continue
		}
		absDev = append(absDev, math.Abs(toFloat64(v)-median))
	}
	if len(absDev) == 0 {
		return median, 0
	}
	// Reuses the same histo scratch: PercentileBuf zeroes it on entry.
	mad = PercentileBuf(absDev, 0.5, histo)
	return median, mad
}
