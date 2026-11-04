package stats

import "math"

// MeanStdevSiril computes the mean and sample standard deviation (N-1)
// of data with the two-pass algorithm of Siril's siril_stats_float_sd
// (src/algos/statistics_float.c), and returns exactly what that C
// function returns: both outputs are bit-identical to Siril's, which is
// tested against the verbatim Siril source (testdata/sref) on 424
// inputs from n=0 to n=2*10^7.
//
// Algorithm (operation for operation):
//
//  1. First pass: accumulate sum in float64.
//  2. Round the mean to float32: mean = (float)(sum / n).
//  3. Second pass: subtract in float32 (data[i] - mean), square in
//     float32, accumulate in float64.
//  4. Return sqrtf((float)(acc / (n-1))): the variance is rounded to
//     float32 before the square root.
//
// The float32 mean, the float32 deviations and the float32 variance are
// all load-bearing: doing any of them in float64 changes the result by
// an ULP on a sizeable fraction of inputs, and rejection pipelines
// compare pixels against boundaries derived from this value. On the 424
// golden inputs, keeping the float32 mean and deviations but taking the
// square root of the unrounded float64 variance differs from Siril on 43
// of them; computing everything in float64 and rounding only the result
// differs on 50 (TestMeanStdevSiril_DocCounts).
//
// NaN handling: none, as in Siril. NaN inputs propagate by IEEE rules
// and both outputs are NaN. Callers must pre-filter NaN if needed (cf.
// Siril's reassign_to_non_null_data_float, FilterNonZero).
//
// Edge cases (Siril's behaviour, also tested):
//   - n=0: returns (NaN, -0). Mean is 0/0 = NaN; the deviation
//     accumulator is 0 and 0/(n-1) = 0/-1 = -0.
//   - n=1: returns (data[0], NaN). The single deviation is 0 and
//     0/(n-1) = 0/0 = NaN.
//
// For a friendly default (e.g. (0, 0) for empty input) use MeanStdev or
// MeanStdevTwoPass; both return float64 results that differ from this
// function at the ULP level even on well-formed inputs.
func MeanStdevSiril(data []float32) (mean, stdev float32) {
	n := len(data)
	var acc float64
	for _, v := range data {
		acc += float64(v)
	}
	mean = float32(acc / float64(n))

	acc = 0
	for _, v := range data {
		// (data[i] - mean) * (data[i] - mean): float subtraction and
		// float product, promoted to double for the +=.
		d := v - mean
		acc += float64(d * d)
	}
	return mean, float32(math.Sqrt(float64(float32(acc / float64(n-1)))))
}
