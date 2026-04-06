package stats

// FilterNonZero returns a new slice containing only the non-zero,
// non-NaN elements of data. This matches the astronomical convention
// for "blank pixel" exclusion — Siril's reassign_to_non_null_data_float
// uses the same logic before computing median and MAD.
//
// For integer types, excludes exact zeros (NaN is not possible).
// For float types, excludes both exact 0.0 and NaN.
//
// If all elements are zero/NaN, returns the original slice as a
// fallback (avoids returning empty data to downstream statistics).
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
