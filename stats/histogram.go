package stats

import (
	"errors"
	"fmt"
	"math"
)

// Histogram holds the result of a histogram computation.
type Histogram struct {
	Counts []int64   // bin counts, length = NBins
	Edges  []float64 // bin edges, length = NBins + 1
	NBins  int
	Total  int64 // sum of all counts (number of finite values binned)
}

// ErrHistogramRange is returned (wrapped) by BuildHistogram where
// numpy.histogram raises "Too many bins for data range": the value range
// cannot be split into nbins finite, strictly increasing bins.
var ErrHistogramRange = errors.New("stats: too many bins for data range")

// BuildHistogram bins data into nbins equal-width bins (256 if nbins <= 0)
// and returns the counts and edges that numpy.histogram(x, bins=nbins)
// returns for x = the finite values of data converted to float64; the
// tests compare both bit-for-bit. Specifically:
//
//   - NaN and ±Inf values are excluded (numpy would raise on ±Inf); Total
//     is the number of values binned, so len(data) - Total of them were
//     excluded.
//   - The range is [min, max] of those values; if min == max the range
//     is [min-0.5, max+0.5], as numpy does.
//   - Edges[i] = min + i*(max-min)/nbins, with Edges[nbins] = max exactly
//     (numpy.linspace).
//   - Bin i holds Edges[i] <= v < Edges[i+1]; the last bin also holds
//     v == Edges[nbins].
//
// With no finite values, all counts are 0, Total is 0 and the edges span
// [0, 1] (numpy's empty-input range).
//
// Where numpy.histogram raises instead of binning, BuildHistogram returns
// a zero Histogram and an error wrapping ErrHistogramRange (its only
// error). That is wherever the edges are not finite and strictly
// increasing: a single-valued range whose ±0.5 padding is lost (|v| >=
// 2^53, or 2^52 with more than one bin), a max-min that overflows float64
// (numpy raises IndexError for one bin, ValueError otherwise), and a range
// too narrow for nbins distinct float64 edges (numpy's "Too many bins for
// data range", e.g. 100 bins over integers clustered at 2^60, where
// float64 values are 256 apart, or 2 bins over [0, 5e-324]).
func BuildHistogram[T Numeric](data []T, nbins int) (Histogram, error) {
	if nbins <= 0 {
		nbins = 256
	}

	var total int64
	var minVal, maxVal float64
	for _, v := range data {
		fv := toFloat64(v)
		if math.IsNaN(fv) || math.IsInf(fv, 0) {
			continue
		}
		if total == 0 {
			minVal, maxVal = fv, fv
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

	first, last := minVal, maxVal
	if total == 0 {
		first, last = 0, 1
	}
	if first == last {
		first -= 0.5
		last += 0.5
	}
	rangeErr := func() (Histogram, error) {
		return Histogram{}, fmt.Errorf("%w: cannot create %d finite-sized bins over [%g, %g]", ErrHistogramRange, nbins, minVal, maxVal)
	}

	delta := last - first
	if delta == 0 || math.IsInf(delta, 0) {
		// |values| >= 2^52 so the ±0.5 padding is lost, or max-min
		// overflows: numpy's linspace edges are not finite and distinct.
		return rangeErr()
	}

	h := Histogram{
		Counts: make([]int64, nbins),
		Edges:  make([]float64, nbins+1),
		NBins:  nbins,
		Total:  total,
	}

	// numpy.linspace(first, last, nbins+1): i*step + first, last forced.
	step := delta / float64(nbins)
	for i := 0; i <= nbins; i++ {
		if step == 0 {
			h.Edges[i] = float64(float64(float64(i)/float64(nbins))*delta) + first
		} else {
			h.Edges[i] = float64(float64(i)*step) + first
		}
	}
	h.Edges[nbins] = last
	for i := 0; i < nbins; i++ {
		if h.Edges[i] >= h.Edges[i+1] {
			// numpy: "Too many bins for data range. Cannot create
			// nbins finite-sized bins."
			return rangeErr()
		}
	}

	if total == 0 {
		return h, nil
	}

	// numpy's uniform-bin index computation and its edge corrections.
	for _, v := range data {
		x := toFloat64(v)
		if math.IsNaN(x) || math.IsInf(x, 0) {
			continue
		}
		f := float64(float64((x-first)/delta) * float64(nbins))
		i := int(f)
		if i == nbins {
			i--
		}
		if x < h.Edges[i] {
			i--
		}
		if x >= h.Edges[i+1] && i != nbins-1 {
			i++
		}
		h.Counts[i]++
	}
	return h, nil
}

// CDF returns the cumulative distribution function as a float64 slice
// of length NBins: cdf[i] = (Counts[0] + ... + Counts[i]) / Total, the
// fraction of the binned values in bins 0..i. All zeros if Total is 0.
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

// Percentile returns the p-th percentile (p in [0, 1], clamped) of the
// histogram by inverting its piecewise-linear CDF: with t = p*Total, it
// finds the first non-empty bin i whose cumulative count reaches t and
// returns Edges[i] + (t - countBelow)/Counts[i] * (Edges[i+1] - Edges[i]).
// For p = 0 that is the lower edge of the first non-empty bin, for p = 1
// the upper edge of the last one. Where the CDF is flat at p (t equals
// the count below an empty bin), every point of the flat part is a
// p-quantile and this is its left end. The tests compare it with that
// definition evaluated in exact rational arithmetic (within 4 ULP of the
// edge magnitudes), and, for 0 < p < 1 off the flat parts, with
// scipy.stats.rv_histogram's ppf of the same counts and edges. Returns 0
// for an empty histogram and NaN for a NaN p.
func (h *Histogram) Percentile(p float64) float64 {
	if p != p {
		return math.NaN()
	}
	if h.Total == 0 || h.NBins == 0 {
		return 0
	}
	if p < 0 {
		p = 0
	}
	if p > 1 {
		p = 1
	}
	target := float64(p * float64(h.Total))
	var cum int64
	for i, c := range h.Counts {
		cum += c
		if float64(cum) >= target && c > 0 {
			prevCum := cum - c
			frac := (target - float64(prevCum)) / float64(c)
			return h.Edges[i] + float64(frac*(h.Edges[i+1]-h.Edges[i]))
		}
	}
	return h.Edges[h.NBins]
}
