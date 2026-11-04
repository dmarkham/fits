package fits

import (
	"fmt"
)

// ReadSubset reads an n-dimensional hyperslab from h.
//
// lower, upper, and stride each must have length equal to h.NAXIS(). The
// slab selected is the half-open interval [lower[i], upper[i]) along each
// axis, stepping by stride[i] (stride==1 means contiguous). Coordinates are
// 0-based and match the order of h.Shape() (axis 0 varies fastest on disk).
//
// The output slice has length Π ceil((upper[i]-lower[i]) / stride[i]) and
// is laid out row-major matching the same axis order as Shape() — i.e. the
// first axis is the fastest-varying dimension, as on disk.
//
// Unlike ReadPixels, ReadSubset does not issue one big read. It computes a
// per-row byte range along axis 0 and reads each row independently, seeking
// between them. For stride[0]==1 and contiguous axis-0 runs this is
// near-optimal; for large strided axes it issues more seeks than an
// all-at-once read would.
//
// BSCALE/BZERO are applied, and target types accepted or refused, exactly
// as in ReadPixels (integer targets exact, float targets equal to
// cfitsio's fits_read_subset, with ReadPixels' one difference: a value
// beyond the float32 range reads as +-Inf into float32, where cfitsio
// returns +-FLT_MAX). BLANK is not applied, as in ReadPixels. Any
// NAXIS >= 1 is supported.
func ReadSubset[T Numeric](h *ImageHDU, lower, upper, stride []int64) ([]T, error) {
	naxis := h.rec.naxis
	if len(lower) != naxis || len(upper) != naxis || len(stride) != naxis {
		return nil, fmt.Errorf("fits: ReadSubset: lower/upper/stride must have length NAXIS=%d", naxis)
	}
	shape := h.rec.shape
	for i := range naxis {
		if lower[i] < 0 || upper[i] > shape[i] || lower[i] >= upper[i] {
			return nil, fmt.Errorf("fits: ReadSubset: axis %d range [%d,%d) out of bounds (shape %d)", i, lower[i], upper[i], shape[i])
		}
		if stride[i] < 1 {
			return nil, fmt.Errorf("fits: ReadSubset: axis %d stride %d < 1", i, stride[i])
		}
	}
	bp := h.bitpixType()
	conv, err := pixelConverter[T](h)
	if err != nil {
		return nil, err
	}

	// Count output elements along each axis.
	outLen := int64(1)
	perAxisCount := make([]int64, naxis)
	for i := range naxis {
		perAxisCount[i] = (upper[i] - lower[i] + stride[i] - 1) / stride[i]
		outLen *= perAxisCount[i]
	}
	out := make([]T, outLen)
	elemSize := int64(bp.Size())

	// We iterate over every combination of axis 1..(naxis-1) indices (call
	// these "rows") and for each row read an axis-0 slab. If stride[0]==1
	// the slab is contiguous on disk; otherwise we read the full axis-0
	// extent and stride in memory.
	rowReadLen := upper[0] - lower[0]
	rowBytes := rowReadLen * elemSize

	// Strides in file coordinates: fileStride[i] = elemSize * Π_{j<i} shape[j].
	fileStride := make([]int64, naxis)
	var prod int64 = 1
	for i := range naxis {
		fileStride[i] = prod * elemSize
		prod *= shape[i]
	}

	// Output strides: outStride[i] = Π_{j<i} perAxisCount[j].
	outStride := make([]int64, naxis)
	prod = 1
	for i := range naxis {
		outStride[i] = prod
		prod *= perAxisCount[i]
	}

	rowBuf := make([]byte, rowBytes)
	axisIdx := make([]int64, naxis) // current in-subset index
	// Iterate axes 1..naxis-1 via nested counters.
	for {
		// Compute file offset for this row's axis-0 start.
		var fileOff int64 = 0
		for i := 1; i < naxis; i++ {
			fileOff += (lower[i] + axisIdx[i]*stride[i]) * fileStride[i]
		}
		fileOff += lower[0] * fileStride[0]
		if err := h.rec.file.br.ReadRange(h.rec.dataStart+fileOff, rowBuf); err != nil {
			return nil, err
		}
		// Compute output base for this row.
		var outBase int64 = 0
		for i := 1; i < naxis; i++ {
			outBase += axisIdx[i] * outStride[i]
		}
		// Decode axis-0 slab (with stride).
		count := perAxisCount[0]
		for k := range count {
			srcPix := k * stride[0]
			out[outBase+k] = conv.one(rowBuf[srcPix*elemSize : srcPix*elemSize+elemSize])
		}

		// Advance the axis 1.. counters.
		if naxis == 1 {
			break
		}
		i := 1
		for i < naxis {
			axisIdx[i]++
			if axisIdx[i] < perAxisCount[i] {
				break
			}
			axisIdx[i] = 0
			i++
		}
		if i == naxis {
			break
		}
	}
	return out, nil
}
