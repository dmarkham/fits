package stats

import "math"

// CenterFunc selects the center estimator for sigma clipping.
type CenterFunc int

const (
	CenterMean   CenterFunc = iota // Use arithmetic mean as center
	CenterMedian                   // Use Median (histogram median) as center
)

// SigmaClipResult holds the output of iterative sigma clipping.
type SigmaClipResult struct {
	Mean   float64 // mean of the surviving values
	Stdev  float64 // population standard deviation (ddof=0) of the survivors
	Median float64 // Median of the survivors
	NGood  int     // count of values that survived clipping
	// NIter is the number of clipping iterations run, counted as
	// astropy's SigmaClip counts them: the final iteration that rejects
	// nothing is included.
	NIter int
}

// SigmaClip performs iterative sigma clipping on data with the semantics
// of astropy.stats.sigma_clip(data, sigma_lower=kLow, sigma_upper=kHigh,
// maxiters=maxIter, cenfunc=..., stdfunc="std"):
//
//   - NaN and ±Inf values are removed before clipping begins (astropy
//     does the same and warns); NGood counts finite survivors only.
//   - Each iteration computes the center (the mean, or Median for
//     CenterMedian) and sigma, the population standard deviation (ddof=0)
//     about the mean of the current values, and keeps the values v with
//     center - sigma*kLow <= v <= center + sigma*kHigh. kLow and kHigh
//     are independent (astropy's sigma_lower and sigma_upper).
//   - Iteration stops when an iteration rejects nothing or after maxIter
//     iterations. maxIter <= 0 is astropy's maxiters=None: iterate until
//     an iteration rejects nothing, with no cap (every other iteration
//     removes at least one value, so this terminates).
//
// The mean and sigma are computed accurately (the compensated two-pass
// method of MeanStdevTwoPass: the mean correctly rounded, sigma within 2
// ULP of the exact value on every tested input), and
// the tests run astropy with the exactly rounded mean and population std
// as cenfunc and stdfunc: the surviving set and the iteration count are
// identical to that. numpy's np.mean/np.std (astropy's default "mean" and
// "std") carry their own rounding, which is large for data with a big
// offset, so on such data plain astropy can keep a slightly different set.
//
// The result holds the mean, population stdev and Median of the
// survivors. Median is the histogram median (Siril's for []float32), not
// numpy's exact median, so CenterMedian matches astropy run with that
// median as cenfunc; this is what the tests check (testdata/gen_golden.py
// runs astropy with Siril's median compiled from source).
//
// Values are compared in float64. For []float32 input the working set
// stays float32 so that Median uses the Siril float32 path (a named type
// whose underlying type is float32 uses the float64 path). Data is copied
// internally; the caller's slice is never modified. With no finite values
// the result is the zero SigmaClipResult; if an iteration rejects every
// value, the result is zero except NIter.
func SigmaClip[T Numeric](data []T, kLow, kHigh float64, maxIter int, center CenterFunc) SigmaClipResult {
	if f32, ok := any(data).([]float32); ok {
		work := make([]float32, 0, len(f32))
		for _, v := range f32 {
			if !math.IsNaN(float64(v)) && !math.IsInf(float64(v), 0) {
				work = append(work, v)
			}
		}
		return sigmaClipLoop(work, kLow, kHigh, maxIter, center)
	}
	work := make([]float64, 0, len(data))
	for _, v := range data {
		x := toFloat64(v)
		if !math.IsNaN(x) && !math.IsInf(x, 0) {
			work = append(work, x)
		}
	}
	return sigmaClipLoop(work, kLow, kHigh, maxIter, center)
}

func sigmaClipLoop[S float32 | float64](work []S, kLow, kHigh float64, maxIter int, center CenterFunc) SigmaClipResult {
	if len(work) == 0 {
		return SigmaClipResult{}
	}
	niter := 0
	// astropy's loop: while nchanged != 0 and iteration < maxiters. An
	// iteration on an empty set (everything was rejected) changes nothing
	// and ends it, as in astropy.
	for maxIter <= 0 || niter < maxIter {
		niter++
		m, sigma := meanPopStdev(work)
		c := m
		if center == CenterMedian {
			c = Percentile(work, 0.5)
		}
		lo := c - float64(sigma*kLow)
		hi := c + float64(sigma*kHigh)
		kept := work[:0]
		for _, v := range work {
			fv := float64(v)
			if fv >= lo && fv <= hi {
				kept = append(kept, v)
			}
		}
		changed := len(kept) != len(work)
		work = kept
		if !changed {
			break
		}
	}
	if len(work) == 0 {
		return SigmaClipResult{NIter: niter}
	}
	m, s := meanPopStdev(work)
	return SigmaClipResult{
		Mean:   m,
		Stdev:  s,
		Median: Percentile(work, 0.5),
		NGood:  len(work),
		NIter:  niter,
	}
}

// meanPopStdev returns the mean and the population (ddof=0) standard
// deviation about it, with the compensated two-pass sums of
// MeanStdevTwoPass (NaN for an empty set).
func meanPopStdev[S float32 | float64](data []S) (mean, stdev float64) {
	n, mean, m2 := twoPassMoments(data)
	if n == 0 {
		return mean, math.NaN()
	}
	variance := m2 / float64(n)
	if variance < 0 {
		variance = 0
	}
	return mean, math.Sqrt(variance)
}
