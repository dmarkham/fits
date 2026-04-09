package stats

import "math"

// MeanStdevSiril computes the mean and sample standard deviation (N-1)
// of data using the two-pass algorithm from Siril's siril_stats_float_sd
// (statistics_float.c). Returns float32 to preserve rejection-boundary
// precision — silently promoting to float64 produces ULP-level
// differences at sigma-clip boundaries that affect pixel rejection
// decisions on real astronomical data.
//
// Algorithm (bit-for-bit port of the C source):
//
//  1. First pass: accumulate sum in float64.
//  2. Truncate the mean to float32: mean = float32(sum / n).
//  3. Second pass: subtract in float32 (data[i] - mean), square in
//     float32, accumulate in float64.
//  4. Return float32(sqrt(acc / (n-1))).
//
// The float32 mean truncation and float32 deviation arithmetic are the
// load-bearing parts: the same loops in float64 produce different ULP-
// level results that change pixel rejection decisions. The fits-
// processing rejection pipeline traced a visible artifact (red hotspot
// in stacked output) to this precision gap and validated the two-pass-
// with-truncation algorithm against Siril's actual rejection output on
// real Rosette Nebula data.
//
// NaN handling: NONE — matches Siril. NaN inputs propagate by IEEE
// rules and the result is NaN in both outputs. Callers must pre-filter
// NaN if needed (cf. Siril's reassign_to_non_null_data_float).
//
// Edge cases (also Siril behavior):
//   - n=0: returns (NaN, -0). Mean is 0/0 = NaN; the deviation
//     accumulator is 0 and the (n-1)=-1 division yields -0.
//   - n=1: returns (data[0], NaN). Mean is data[0]; the single
//     deviation is 0 and the (n-1)=0 division yields NaN.
//
// If you want a friendly default (e.g. (0, 0) for empty input), use
// the generic MeanStdev — but be aware its outputs differ from this
// function at the ULP level even for well-formed inputs.
//
// Source: Siril src/algos/statistics_float.c siril_stats_float_sd.
func MeanStdevSiril(data []float32) (mean, stdev float32) {
	n := len(data)
	var acc float64
	for _, v := range data {
		acc += float64(v)
	}
	mean = float32(acc / float64(n))

	acc = 0
	for _, v := range data {
		// float32 subtraction and float32 squaring, matching Siril's
		// (data[i] - mean) * (data[i] - mean) where both operands are
		// float and the result is implicitly promoted to double for
		// the +=.
		d := v - mean
		acc += float64(d * d)
	}
	return mean, float32(math.Sqrt(acc / float64(n-1)))
}
