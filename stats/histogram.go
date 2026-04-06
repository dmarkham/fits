package stats

// Histogram holds the result of a histogram computation.
type Histogram struct {
	Counts []int64    // bin counts, length = NBins
	Edges  []float64  // bin edges, length = NBins + 1
	NBins  int
	Total  int64      // sum of all counts
}

// BuildHistogram bins data into nbins equal-width bins spanning
// [min, max] of the non-NaN data. Returns a Histogram struct.
func BuildHistogram[T Numeric](data []T, nbins int) Histogram {
	if nbins <= 0 {
		nbins = 256
	}

	// Find min/max, count non-NaN.
	var total int64
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
		total++
	}

	h := Histogram{
		Counts: make([]int64, nbins),
		Edges:  make([]float64, nbins+1),
		NBins:  nbins,
		Total:  total,
	}

	if total == 0 {
		return h
	}

	// Build edges.
	binWidth := (maxVal - minVal) / float64(nbins)
	for i := 0; i <= nbins; i++ {
		h.Edges[i] = minVal + float64(i)*binWidth
	}

	// Bin data.
	if binWidth == 0 {
		// All values the same — put everything in bin 0.
		h.Counts[0] = total
		return h
	}

	scale := float64(nbins) / (maxVal - minVal)
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		bin := int(scale * (toFloat64(v) - minVal))
		if bin < 0 {
			bin = 0
		}
		if bin >= nbins {
			bin = nbins - 1
		}
		h.Counts[bin]++
	}
	return h
}

// CDF returns the cumulative distribution function as a float64 slice
// of length NBins. cdf[i] is the fraction of data at or below the
// upper edge of bin i.
func (h *Histogram) CDF() []float64 {
	cdf := make([]float64, h.NBins)
	if h.Total == 0 {
		return cdf
	}
	var cum int64
	for i, c := range h.Counts {
		cum += c
		cdf[i] = float64(cum) / float64(h.Total)
	}
	return cdf
}

// Percentile returns the p-th percentile (p in [0, 1]) by
// interpolating the CDF of this histogram.
func (h *Histogram) Percentile(p float64) float64 {
	if h.Total == 0 || h.NBins == 0 {
		return 0
	}
	target := p * float64(h.Total)
	var cum int64
	for i, c := range h.Counts {
		cum += c
		if float64(cum) >= target {
			// Linear interpolation within this bin.
			if c == 0 {
				return h.Edges[i]
			}
			prevCum := cum - c
			frac := (target - float64(prevCum)) / float64(c)
			return h.Edges[i] + frac*(h.Edges[i+1]-h.Edges[i])
		}
	}
	return h.Edges[h.NBins]
}
