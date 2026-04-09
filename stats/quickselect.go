package stats

// QuickSelect rearranges data so that data[k] contains the k-th
// smallest element (0-indexed). Returns data[k]. Uses the Hoare-style
// two-pointer partition with the middle element as pivot. O(n)
// average, O(n²) worst case.
//
// This is a generic line-for-line port of Siril's quickmedian_float
// (sorting.c) — same partition, same pivot strategy, same loop
// structure. The fits-processing rejection pipeline already has a
// float32-specific copy (stack/reject.go quickselectF32); this is the
// library version, generic over any Numeric type. Behavioral
// equivalence with the float32 version is intentional, so callers
// can drop the local copy and use this without changing rejection
// decisions.
//
// Unlike Median, this computes the exact order statistic — no
// histogram quantization or sub-bin interpolation. For small n
// (under ~50 elements) this is both faster than Median and more
// accurate, since the histogram method's bin width is comparable to
// the data spacing in that regime. For large n, prefer MedianBuf
// (zero-alloc, O(n)) which is faster despite the histogram overhead.
//
// IMPORTANT differences from Median for the same input:
//   - QuickSelect returns an exact element of data, not an
//     interpolated value between two elements
//   - For even n, "the median" is conventionally the mean of the two
//     middle elements; QuickSelect with k=n/2 returns just the upper
//     of those two
//   - Therefore QuickSelect(data, len(data)/2) is NOT equivalent to
//     Median(data); callers who need a drop-in replacement should
//     compute (QuickSelect(data, n/2-1) + QuickSelect(data, n/2)) / 2
//     for even n
//
// NaN handling: NONE. NaN comparisons always return false, which
// breaks partitioning — NaN values can end up anywhere in the output
// and corrupt the answer. Callers MUST pre-filter NaN before calling.
//
// Mutation: data is rearranged in place. If callers need to preserve
// the original order, copy the slice first.
//
// Panics if data is empty, or if k < 0 or k >= len(data).
//
// Source: Siril src/algos/sorting.c quickmedian_float, mirrored at
// fits-processing/stack/reject.go quickselectF32.
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
		// Middle-element pivot — matches Siril. No median-of-three.
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
