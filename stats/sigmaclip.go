package stats

import "math"

// CenterFunc selects the center estimator for sigma clipping.
type CenterFunc int

const (
	CenterMean   CenterFunc = iota // Use arithmetic mean as center
	CenterMedian                   // Use median as center (more robust)
)

// SigmaClipResult holds the output of iterative sigma clipping.
type SigmaClipResult struct {
	Mean   float64
	Stdev  float64
	Median float64
	NGood  int // count of values that survived clipping
}

// SigmaClip performs iterative sigma clipping on data.
//
// At each iteration, values outside [center - kLow*sigma,
// center + kHigh*sigma] are rejected. The center is recomputed from
// the surviving values using the chosen CenterFunc. Iteration stops
// when no more values are rejected or maxIter is reached (0 = iterate
// until stable, capped at 100).
//
// Data is copied internally — the caller's slice is never modified.
// NaN values are excluded before clipping begins.
func SigmaClip[T Numeric](data []T, kLow, kHigh float64, maxIter int, center CenterFunc) SigmaClipResult {
	// For float32 input, stay in float32 throughout the clipping loop
	// so that Percentile uses the float32 specialization (matching
	// Siril's findMinMaxPercentile which works in float arithmetic).
	if f32, ok := any(data).([]float32); ok {
		return sigmaClipFloat32(f32, kLow, kHigh, maxIter, center)
	}
	return sigmaClipGeneric(data, kLow, kHigh, maxIter, center)
}

func sigmaClipFloat32(data []float32, kLow, kHigh float64, maxIter int, center CenterFunc) SigmaClipResult {
	work := make([]float32, 0, len(data))
	for _, v := range data {
		if v == v { // skip NaN
			work = append(work, v)
		}
	}
	if len(work) == 0 {
		return SigmaClipResult{}
	}
	if maxIter <= 0 {
		maxIter = 100
	}
	for iter := 0; iter < maxIter; iter++ {
		var c, sigma float64
		switch center {
		case CenterMedian:
			c = Percentile(work, 0.5)
		default:
			var sum float64
			for _, v := range work {
				sum += float64(v)
			}
			c = sum / float64(len(work))
		}
		var sum2 float64
		for _, v := range work {
			d := float64(v) - c
			sum2 += d * d
		}
		sigma = math.Sqrt(sum2 / float64(len(work)))
		if sigma == 0 {
			break
		}
		lo := c - kLow*sigma
		hi := c + kHigh*sigma
		kept := work[:0]
		for _, v := range work {
			fv := float64(v)
			if fv >= lo && fv <= hi {
				kept = append(kept, v)
			}
		}
		if len(kept) == len(work) {
			break
		}
		work = kept
	}
	if len(work) == 0 {
		return SigmaClipResult{}
	}
	var sum float64
	for _, v := range work {
		sum += float64(v)
	}
	m := sum / float64(len(work))
	var sum2 float64
	for _, v := range work {
		d := float64(v) - m
		sum2 += d * d
	}
	return SigmaClipResult{
		Mean:   m,
		Stdev:  math.Sqrt(sum2 / float64(len(work))),
		Median: Percentile(work, 0.5),
		NGood:  len(work),
	}
}

func sigmaClipGeneric[T Numeric](data []T, kLow, kHigh float64, maxIter int, center CenterFunc) SigmaClipResult {
	work := make([]float64, 0, len(data))
	for _, v := range data {
		if !isNaN(v) {
			work = append(work, toFloat64(v))
		}
	}
	if len(work) == 0 {
		return SigmaClipResult{}
	}
	if maxIter <= 0 {
		maxIter = 100
	}
	for iter := 0; iter < maxIter; iter++ {
		var c, sigma float64
		switch center {
		case CenterMedian:
			c = Percentile(work, 0.5)
		default:
			c = meanF64(work)
		}
		sigma = stdevF64(work, c)
		if sigma == 0 {
			break
		}
		lo := c - kLow*sigma
		hi := c + kHigh*sigma
		kept := work[:0]
		for _, v := range work {
			if v >= lo && v <= hi {
				kept = append(kept, v)
			}
		}
		if len(kept) == len(work) {
			break
		}
		work = kept
	}
	if len(work) == 0 {
		return SigmaClipResult{}
	}
	m := meanF64(work)
	s := stdevF64(work, m)
	med := Percentile(work, 0.5)
	return SigmaClipResult{
		Mean:   m,
		Stdev:  s,
		Median: med,
		NGood:  len(work),
	}
}

func meanF64(data []float64) float64 {
	if len(data) == 0 {
		return 0
	}
	var s float64
	for _, v := range data {
		s += v
	}
	return s / float64(len(data))
}

func stdevF64(data []float64, mean float64) float64 {
	if len(data) == 0 {
		return 0
	}
	var s float64
	for _, v := range data {
		d := v - mean
		s += d * d
	}
	return math.Sqrt(s / float64(len(data)))
}
