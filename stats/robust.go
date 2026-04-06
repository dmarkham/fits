package stats


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
func Percentile[T Numeric](data []T, p float64) float64 {
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
		return percentileFloat32(f32, float32(p))
	}
	// Generic path for all other types: float64 arithmetic.
	return percentileGeneric(data, p)
}

// percentileFloat32 is the Siril-matching float32 specialization.
// Reimplementation of the findMinMaxPercentile algorithm from
// RawTherapee/Siril using float32 arithmetic for scale, bin
// assignment, and CDF walk — matching the C code's behavior.
// The bin index uses uint16 truncation without clamping, matching
// the C code's static_cast<uint16_t> behavior.
func percentileFloat32(data []float32, p float32) float64 {
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
	if histoSize > 65536 {
		histoSize = 65536
	}
	scale := float32(histoSize-1) / (maxVal - minVal)
	histo := make([]uint32, histoSize)
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

func percentileGeneric[T Numeric](data []T, p float64) float64 {
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
	if histoSize > 65536 {
		histoSize = 65536
	}
	scale := float64(histoSize-1) / (maxVal - minVal)

	histo := make([]uint32, histoSize)
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
func Median[T Numeric](data []T) float64 {
	return Percentile(data, 0.5)
}

// MAD returns the Median Absolute Deviation: median(|xi - median(x)|).
// This is the raw MAD — caller multiplies by 1.4826 for Gaussian-
// equivalent sigma if needed.
func MAD[T Numeric](data []T) float64 {
	_, mad := MADWithMedian(data)
	return mad
}

// MADWithMedian returns the median and the MAD, avoiding double-
// computing the median when the caller needs both.
func MADWithMedian[T Numeric](data []T) (median, mad float64) {
	median = Median(data)

	// Compute |xi - median| as float64.
	absDev := make([]float64, 0, len(data))
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

	// MAD = median of the absolute deviations. Use the same
	// histogram-based percentile on the float64 deviations.
	mad = Percentile(absDev, 0.5)
	return median, mad
}

