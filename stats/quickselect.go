package stats

// QuickSelect rearranges data so that data[k] contains the k-th
// smallest element (0-indexed) and returns it. It uses Hoare's
// two-pointer partition with the middle element as pivot: O(n) average,
// O(n²) worst case. After the call every element before k is <= data[k]
// and every element after it is >= data[k].
//
// The returned value is the exact order statistic, the same value any
// correct selection returns (tested against slices.Sort and against
// Siril's quickmedian_float, src/algos/sorting.c, compiled from source).
// The partition itself is not Siril's (Siril uses a Lomuto partition and
// a sorting network for n < 9), so the order the other elements are left
// in differs; only data[k] and the partition property are guaranteed.
//
// Unlike Median, this computes the exact order statistic — no
// histogram quantization or sub-bin interpolation. For small n this is
// both faster than Median and more accurate, since the histogram
// method's bin width is comparable to the data spacing in that regime.
//
// Relation to Siril's median for []float32 (both checked bit-for-bit
// against quickmedian_float):
//   - odd n: float64(QuickSelect(data, n/2)) equals quickmedian_float.
//   - even n >= 9: quickmedian_float returns
//     (float64(lo) + float64(hi)) / 2, with lo = QuickSelect(data, n/2-1)
//     and hi = QuickSelect(data, n/2).
//   - even n < 9: Siril's sorting network adds in float32 first:
//     float64(lo + hi) / 2.
//
// QuickSelect(data, len(data)/2) on its own is therefore not the median
// for even n, and not Median (the histogram approximation) for any n.
//
// NaN handling: NONE. NaN comparisons always return false, which
// breaks partitioning — NaN values can end up anywhere in the output
// and corrupt the answer. Callers MUST pre-filter NaN before calling.
//
// Mutation: data is rearranged in place. If callers need to preserve
// the original order, copy the slice first.
//
// Panics if data is empty, or if k < 0 or k >= len(data).
func QuickSelect[T Numeric](data []T, k int) T {
	n := len(data)
	if n == 0 {
		panic("stats.QuickSelect: empty data")
	}
	if k < 0 || k >= n {
		panic("stats.QuickSelect: k out of range")
	}

	lo, hi := 0, n-1
	for lo < hi {
		// Middle-element pivot. No median-of-three.
		// For random astronomical data this is fine; the worst case
		// (sorted/reverse-sorted/adversarial) is rare in practice.
		pivot := data[lo+(hi-lo)/2]
		i, j := lo, hi
		for i <= j {
			for data[i] < pivot {
				i++
			}
			for data[j] > pivot {
				j--
			}
			if i <= j {
				data[i], data[j] = data[j], data[i]
				i++
				j--
			}
		}
		// After the partition: data[lo..j] are <= pivot,
		// data[i..hi] are >= pivot. The gap (j+1..i-1) holds elements
		// equal to the pivot — if k falls in there, we're done.
		if k <= j {
			hi = j
		} else if k >= i {
			lo = i
		} else {
			break
		}
	}
	return data[k]
}
