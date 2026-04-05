package fits

import (
	"fmt"
	"iter"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/tform"
)

// VarColumn stores all rows' variable-length data in a single contiguous
// slab plus an offsets table. This layout mirrors the FITS on-disk heap and
// avoids one allocation per row.
//
// See plan decision 5.5.
type VarColumn[T Numeric] struct {
	values  []T
	offsets []int64 // len = nrows+1; values[offsets[r]:offsets[r+1]] is row r
}

// Len returns the number of rows (not elements).
func (v *VarColumn[T]) Len() int { return len(v.offsets) - 1 }

// At returns a zero-copy view into the backing Values slab for row r.
// The returned slice aliases the backing buffer; if the caller needs an
// independent copy, they must clone it.
func (v *VarColumn[T]) At(row int) []T {
	return v.values[v.offsets[row]:v.offsets[row+1]]
}

// Rows returns an iter.Seq2[int, []T] that yields (rowIndex, zeroCopySlice)
// for each row in insertion order.
func (v *VarColumn[T]) Rows() iter.Seq2[int, []T] {
	return func(yield func(int, []T) bool) {
		for r := 0; r < v.Len(); r++ {
			if !yield(r, v.At(r)) {
				return
			}
		}
	}
}

// Raw returns the backing slab and offsets table directly. The caller may
// pass them to GPU, SIMD, or bulk-transfer code without an intermediate
// copy. Mutating either slice invalidates the VarColumn.
func (v *VarColumn[T]) Raw() (values []T, offsets []int64) {
	return v.values, v.offsets
}

// ReadVarColumn reads a variable-length (P/Q) binary-table column into a
// VarColumn[T]. The declared element type in the TFORM (e.g. "1PE(256)")
// must match T after scale/zero application.
func ReadVarColumn[T Numeric](t *BinaryTableHDU, col int) (*VarColumn[T], error) {
	cols, err := t.Columns()
	if err != nil {
		return nil, err
	}
	if col < 1 || col > len(cols) {
		return nil, fmt.Errorf("fits: column %d out of range [1,%d]", col, len(cols))
	}
	c := cols[col-1]
	if c.bin.Type != tform.BinPArrayDesc32 && c.bin.Type != tform.BinQArrayDesc64 {
		return nil, fmt.Errorf("fits: column %d is not variable-length (TFORM=%s)", col, c.TForm)
	}
	nrows := t.NumRows()
	// Read the fixed row area to pull out (nelem, offset) pairs.
	rowBytes := t.rowBytes
	rowsBuf := make([]byte, int64(rowBytes)*nrows)
	if err := t.rec.file.br.ReadRange(t.rec.dataStart, rowsBuf); err != nil {
		return nil, err
	}

	// Parse descriptors row by row.
	type desc struct{ n, off int64 }
	descs := make([]desc, nrows)
	isQ := c.bin.Type == tform.BinQArrayDesc64
	elemSize := c.bin.VarType.ElementSize()
	var totalElems int64
	for r := int64(0); r < nrows; r++ {
		base := int(r)*rowBytes + c.byteOffset
		if !isQ {
			n := int64(bigendian.Int32(rowsBuf[base:]))
			off := int64(bigendian.Int32(rowsBuf[base+4:]))
			descs[r] = desc{n: n, off: off}
		} else {
			n := bigendian.Int64(rowsBuf[base:])
			off := bigendian.Int64(rowsBuf[base+8:])
			descs[r] = desc{n: n, off: off}
		}
		totalElems += descs[r].n
	}

	// Read the heap slab (the entire heap area; we do not try to read per-row
	// because rows may be out-of-order or share regions).
	heapStart := t.heapStart()
	heapSize := t.rec.dataEnd - heapStart
	if heapSize < 0 {
		heapSize = 0
	}
	heap := make([]byte, heapSize)
	if heapSize > 0 {
		if err := t.rec.file.br.ReadRange(heapStart, heap); err != nil {
			return nil, err
		}
	}

	out := &VarColumn[T]{
		values:  make([]T, totalElems),
		offsets: make([]int64, nrows+1),
	}
	var pos int64
	for r := int64(0); r < nrows; r++ {
		out.offsets[r] = pos
		d := descs[r]
		if d.n == 0 {
			continue
		}
		if d.off < 0 || d.off+d.n*int64(elemSize) > int64(len(heap)) {
			return nil, fmt.Errorf("fits: VLA row %d heap range out of bounds", r)
		}
		slab := heap[d.off : d.off+d.n*int64(elemSize)]
		for i := int64(0); i < d.n; i++ {
			src := slab[i*int64(elemSize) : (i+1)*int64(elemSize)]
			var f float64
			switch c.bin.VarType {
			case tform.BinUint8:
				f = float64(src[0])
			case tform.BinInt16:
				f = float64(bigendian.Int16(src))
			case tform.BinInt32:
				f = float64(bigendian.Int32(src))
			case tform.BinInt64:
				f = float64(bigendian.Int64(src))
			case tform.BinFloat32:
				f = float64(bigendian.Float32(src))
			case tform.BinFloat64:
				f = bigendian.Float64(src)
			}
			if c.Scale != 0 && c.Scale != 1 {
				f *= c.Scale
			}
			f += c.Zero
			out.values[pos+i] = numericFromFloat[T](f)
		}
		pos += d.n
	}
	out.offsets[nrows] = pos
	return out, nil
}
