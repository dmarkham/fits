package stretch

import "math"

// MTF applies the Midtone Transfer Function to a single channel.
// This is the core of Siril's autostretch command.
//
//	MTF(x, m) = (m - 1) * x / ((2*m - 1) * x - m)
//
// where x is the input pixel (rescaled to [0, 1] after shadow/highlight
// clipping) and m is the midtone balance (0 < m < 1; m=0.5 is identity).
//
// shadows is the black clip point (pixels <= shadows become 0).
// midtones is the midtone balance parameter.
// highlights is the white clip point (pixels >= highlights become 1).
//
// Port of Siril's MTF function (mtf.c:95).
func MTF(pixels []float32, shadows, midtones, highlights float32) {
	if midtones == 0.5 && shadows == 0 && highlights == 1 {
		return // identity
	}
	invRange := float32(1) / (highlights - shadows)
	for i, v := range pixels {
		if v <= shadows {
			pixels[i] = 0
			continue
		}
		if v >= highlights {
			pixels[i] = 1
			continue
		}
		xp := (v - shadows) * invRange
		// MTF formula: (m-1)*xp / ((2m-1)*xp - m)
		pixels[i] = (midtones - 1) * xp / ((2*midtones-1)*xp - midtones)
	}
}

// AutostretchParams computes the MTF parameters (shadows, midtones,
// highlights) from image statistics and the given shadow clipping
// and target background values.
//
// shadowsClip is negative (e.g. -2.8) — the number of MAD-normalized
// sigma below the median to place the shadow clip point.
// targetBG is the target midtone value (e.g. 0.25).
//
// For multi-channel linked mode, pass the concatenation of all channels
// or call AutostretchLinkedParams.
//
// Port of Siril's find_unlinked_midtones_balance (mtf.c:374).
func AutostretchParams(pixels []float32, shadowsClip, targetBG float32) (shadows, midtones, highlights float32) {
	med := medianFloat32(pixels)
	mad := madFloat32(pixels, med)
	const madNorm = 1.4826
	madN := mad * madNorm

	if madN == 0 {
		madN = 0.001
	}

	c0 := med + shadowsClip*madN
	if c0 < 0 {
		c0 = 0
	}

	m2 := med - c0
	shadows = c0
	highlights = 1.0
	// midtones = MTF(m2, targetBG, 0, 1)
	if m2 <= 0 {
		midtones = 0.5
	} else if m2 >= 1 {
		midtones = 0.5
	} else {
		midtones = (targetBG - 1) * m2 / ((2*targetBG-1)*m2 - targetBG)
	}
	return shadows, midtones, highlights
}

// AutostretchLinkedParams computes linked MTF parameters by averaging
// statistics across all channels. Each channel is a separate []float32.
//
// Port of Siril's find_linked_midtones_balance (mtf.c:259).
func AutostretchLinkedParams(channels [][]float32, shadowsClip, targetBG float32) (shadows, midtones, highlights float32) {
	n := float32(len(channels))
	var c0sum, medSum float32
	const madNorm = 1.4826

	for _, ch := range channels {
		med := medianFloat32(ch)
		mad := madFloat32(ch, med)
		madN := mad * madNorm
		if madN == 0 {
			madN = 0.001
		}
		c0sum += med + shadowsClip*madN
		medSum += med
	}

	c0 := c0sum / n
	if c0 < 0 {
		c0 = 0
	}
	m2 := medSum/n - c0

	shadows = c0
	highlights = 1.0
	if m2 <= 0 || m2 >= 1 {
		midtones = 0.5
	} else {
		midtones = (targetBG - 1) * m2 / ((2*targetBG-1)*m2 - targetBG)
	}
	return shadows, midtones, highlights
}

// histogramPercentile computes the p-th percentile (p in [0,1]) of a
// float32 array using the histogram-based interpolated algorithm from
// Siril/RawTherapee (rt_algo.cc:findMinMaxPercentile). This matches
// Siril's statistics output bit-for-bit, unlike a sort-based median
// which gives a slightly different value due to binning interpolation.
//
// The algorithm:
//  1. Find min/max of the data
//  2. Build a histogram with min(65536, n) bins in [min, max]
//  3. Walk the CDF to find the percentile
//  4. Interpolate between the two bins straddling the percentile
func histogramPercentile(data []float32, p float32) float32 {
	n := len(data)
	if n == 0 {
		return 0
	}
	// Find min, max.
	minVal, maxVal := data[0], data[0]
	for _, v := range data[1:] {
		if v < minVal {
			minVal = v
		}
		if v > maxVal {
			maxVal = v
		}
	}
	if maxVal == minVal {
		return minVal
	}
	// Histogram with min(65536, n) bins.
	histoSize := n
	if histoSize > 65536 {
		histoSize = 65536
	}
	scale := float32(histoSize-1) / (maxVal - minVal)
	histo := make([]int, histoSize)
	for _, v := range data {
		bin := int(scale * (v - minVal))
		if bin < 0 {
			bin = 0
		}
		if bin >= histoSize {
			bin = histoSize - 1
		}
		histo[bin]++
	}
	// Walk CDF to find percentile.
	thresh := p * float32(n)
	k := 0
	count := 0
	for count < int(thresh) {
		count += histo[k]
		k++
	}
	// Interpolate (matching Siril's rt_algo.cc:270-277).
	var result float32
	if k > 0 {
		countPrev := count - histo[k-1]
		c0 := float32(count) - thresh
		c1 := thresh - float32(countPrev)
		result = (c1*float32(k) + c0*float32(k-1)) / (c0 + c1)
	} else {
		result = float32(k)
	}
	// Convert back to original range.
	result = result/scale + minVal
	if result < minVal {
		result = minVal
	}
	if result > maxVal {
		result = maxVal
	}
	return result
}

// medianFloat32 computes the median using Siril's histogram-based
// interpolated percentile algorithm. This matches Siril's
// histogram_median_float exactly.
func medianFloat32(data []float32) float32 {
	return histogramPercentile(data, 0.5)
}

// madFloat32 computes the Median Absolute Deviation: median(|xi - median|).
// The raw MAD (without the 1.4826 normalization factor). Uses the same
// histogram-based median algorithm as medianFloat32.
func madFloat32(data []float32, median float32) float32 {
	if len(data) == 0 {
		return 0
	}
	tmp := make([]float32, len(data))
	for i, v := range data {
		d := v - median
		if d < 0 {
			d = -d
		}
		tmp[i] = d
	}
	return histogramPercentile(tmp, 0.5)
}

// Autostretch applies the full autostretch pipeline to a single channel:
// compute statistics, derive MTF parameters, apply the transfer function.
func Autostretch(pixels []float32, shadowsClip, targetBG float32) {
	shadows, midtones, highlights := AutostretchParams(pixels, shadowsClip, targetBG)
	MTF(pixels, shadows, midtones, highlights)
}

// AutostretchLinked applies the linked autostretch pipeline to an RGB
// image: compute averaged statistics, derive one set of MTF parameters,
// apply the same transfer function to all channels.
func AutostretchLinked(channels [][]float32, shadowsClip, targetBG float32) {
	shadows, midtones, highlights := AutostretchLinkedParams(channels, shadowsClip, targetBG)
	for _, ch := range channels {
		MTF(ch, shadows, midtones, highlights)
	}
}

// MTFValue evaluates the MTF formula at a single point.
// Exported for testing and direct use.
func MTFValue(x, m float32) float32 {
	if x <= 0 {
		return 0
	}
	if x >= 1 {
		return 1
	}
	return (m - 1) * x / ((2*m-1)*x - m)
}

// ensure math is used
var _ = math.Abs
