package stats

import "fmt"

// MaxHistoSize is the upper bound on histogram bins used by the
// percentile/median family. Callers using the *Buf zero-alloc
// variants should pre-allocate a histogram slice of this length
// (or len(data), whichever is smaller).
const MaxHistoSize = 65536

// Percentile returns the p-th percentile (p in [0, 1], clamped) of data using
// a histogram-based interpolated algorithm. This is a line-for-line
// port of RawTherapee/Siril's findMinMaxPercentile (rt_algo.cc:167),
// made generic over any Numeric type.
//
// The algorithm:
//  1. Find min/max of the non-NaN data
//  2. Build a histogram with min(65536, n) bins spanning [min, max]
//  3. Walk the CDF to the percentile
//  4. Interpolate between the two bins straddling the percentile
//  5. Convert back to the original value range
//
// For empty slices (or all-NaN), returns 0.
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
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	// For float32 input, use float32 arithmetic for the histogram
	// to match Siril/RawTherapee's findMinMaxPercentile exactly
	// (it operates on float data with float scale/bin computation).
	if f32, ok := any(data).([]float32); ok {
		return percentileFloat32Buf(f32, float32(p), histo)
	}
	// Generic path for all other types: float64 arithmetic.
	return percentileGenericBuf(data, p, histo)
}

// percentileFloat32Buf is the Siril-matching float32 specialization,
// taking a caller-supplied histogram scratch buffer.
//
// Reimplementation of the findMinMaxPercentile algorithm from
// RawTherapee/Siril using float32 arithmetic for scale, bin
// assignment, and CDF walk — matching the C code's behavior.
// The bin index uses uint16 truncation without clamping, matching
// the C code's static_cast<uint16_t> behavior.
func percentileFloat32Buf(data []float32, p float32, histo []uint32) float64 {
	n := 0
	for _, v := range data {
		if v == v { // skip NaN
			n++
		}
	}
	if n == 0 {
		return 0
	}
	var minVal, maxVal float32
	first := true
	for _, v := range data {
		if v != v {
			continue
		}
		if first {
			minVal, maxVal = v, v
			first = false
		} else {
			if v < minVal {
				minVal = v
			}
			if v > maxVal {
				maxVal = v
			}
		}
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
	for _, v := range data {
		if v != v {
			continue
		}
		bin := uint16(scale * (v - minVal))
		histo[bin]++
	}
	thresh := p * float32(n)
	k := 0
	var count float32
	for count < thresh && k < len(histo) {
		count += float32(histo[k])
		k++
	}
	var result float32
	if k > 0 {
		countPrev := count - float32(histo[k-1])
		c0 := count - thresh
		c1 := thresh - countPrev
		result = (c1*float32(k) + c0*float32(k-1)) / (c0 + c1)
	} else {
		result = float32(k)
	}
	result = result/scale + minVal
	if result < minVal {
		result = minVal
	}
	if result > maxVal {
		result = maxVal
	}
	return float64(result)
}

func percentileGenericBuf[T Numeric](data []T, p float64, histo []uint32) float64 {
	// Collect non-NaN values and find min/max in one pass.
	n := 0
	var minVal, maxVal float64
	started := false
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		fv := toFloat64(v)
		if !started {
			minVal, maxVal = fv, fv
			started = true
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
	if maxVal-minVal == 0 {
		return minVal
	}

	// Histogram with min(65536, n) bins.
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
		bin := int(scale * (toFloat64(v) - minVal))
		if bin < 0 {
			bin = 0
		}
		if bin >= histoSize {
			bin = histoSize - 1
		}
		histo[bin]++
	}

	// Walk CDF.
	thresh := p * float64(n)
	k := 0
	var count float64
	for count < thresh {
		count += float64(histo[k])
		k++
	}

	// Interpolate.
	var result float64
	if k > 0 {
		countPrev := count - float64(histo[k-1])
		c0 := count - thresh
		c1 := thresh - countPrev
		result = (c1*float64(k) + c0*float64(k-1)) / (c0 + c1)
	} else {
		result = float64(k)
	}

	// Back to original range + clamp.
	result = result/scale + minVal
	if result < minVal {
		result = minVal
	}
	if result > maxVal {
		result = maxVal
	}
	return result
}

// Median returns the median of data (the 50th percentile).
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

// MAD returns the Median Absolute Deviation: median(|xi - median(x)|).
// This is the raw MAD — caller multiplies by 1.4826 for Gaussian-
// equivalent sigma if needed.
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
// Note: absDev is []float64 (not generic []T) so a single buffer
// works for any T. The deviations are computed as float64
// regardless of input type, matching the existing MADWithMedian.
func MADWithMedianBuf[T Numeric](data []T, histo []uint32, absDev []float64) (median, mad float64) {
	if cap(absDev) < len(data) {
		panic(fmt.Sprintf("stats: absDev buffer too small: got cap %d, need %d", cap(absDev), len(data)))
	}

	median = MedianBuf(data, histo)

	// Reset absDev length and append |xi - median|.
	absDev = absDev[:0]
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		d := toFloat64(v) - median
		if d < 0 {
			d = -d
		}
		absDev = append(absDev, d)
	}
	if len(absDev) == 0 {
		return median, 0
	}

	// MAD = median of the absolute deviations. Reuses the same histo
	// scratch — PercentileBuf zeroes it on entry, so the previous
	// median pass's contents don't corrupt this one.
	mad = PercentileBuf(absDev, 0.5, histo)
	return median, mad
}

