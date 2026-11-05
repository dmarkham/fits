package stats

// FilterNonZero returns a new slice containing only the non-zero,
// non-NaN elements of data, in order. The predicate (v != 0 && !isnan(v),
// so -0 is excluded too) is the one Siril's
// reassign_to_non_null_data_float uses to exclude blank pixels before
// computing the median and MAD; for []float32 the output is checked
// element-for-element against that Siril function.
//
// For integer types, excludes exact zeros (NaN is not possible).
//
// If no element survives (all zero/NaN, or empty input), FilterNonZero
// returns data itself, not a copy — unlike Siril, which computes no
// statistics in that case. A caller that must tell the two apart can
// check whether the first element of the result is zero or NaN.
func FilterNonZero[T Numeric](data []T) []T {
	out := make([]T, 0, len(data))
	for _, v := range data {
		if isNaN(v) {
			continue
		}
		if v == 0 {
			continue
		}
		out = append(out, v)
	}
	if len(out) == 0 {
		return data
	}
	return out
}
