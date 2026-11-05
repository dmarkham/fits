package fits

import (
	"fmt"
	"iter"
	"math"
	"reflect"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/tform"
)

// VarColumn stores all rows' variable-length data in a single contiguous
// slab plus an offsets table. This layout mirrors the FITS on-disk heap and
// avoids one allocation per row.
//
// See plan decision 5.5.
type VarColumn[T any] struct {
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

// vlaDescriptor is one row's (element count, heap byte offset) pair.
type vlaDescriptor struct{ n, off int64 }

// varCells reads the descriptors of a P/Q column, validates every one
// against the heap (THEAP .. NAXIS1*NAXIS2+PCOUNT, §7.3.5) and calls fn
// with each row's element count and payload bytes. elemBytes gives the
// payload size of n elements.
func (t *BinaryTableHDU) varCells(c Column, elemBytes func(n int64) int64, fn func(r int64, n int64, payload []byte) error) error {
	if !isVarType(c.bin.Type) {
		return wrongReader(c, "not a variable-length column; "+readerFor(c))
	}
	nrows := t.NumRows()
	descs := make([]vlaDescriptor, nrows)
	if c.Repeat == 0 {
		// "0P..." columns occupy no bytes: every row is empty.
		for r := int64(0); r < nrows; r++ {
			if err := fn(r, 0, nil); err != nil {
				return err
			}
		}
		return nil
	}
	heapStart, heapSize, err := t.heapRange()
	if err != nil {
		return err
	}
	isQ := c.bin.Type == tform.BinQArrayDesc64
	lo, hi := int64(-1), int64(0)
	err = t.cells(c, func(r int64, cell []byte) error {
		var d vlaDescriptor
		if isQ {
			d = vlaDescriptor{bigendian.Int64(cell), bigendian.Int64(cell[8:])}
		} else {
			d = vlaDescriptor{int64(bigendian.Int32(cell)), int64(bigendian.Int32(cell[4:]))}
		}
		if d.n < 0 || d.off < 0 {
			return fmt.Errorf("fits: column %d row %d: invalid VLA descriptor (%d, %d)", c.Index, r, d.n, d.off)
		}
		descs[r] = d
		if d.n > 0 {
			// Every element is at most 16 bytes, so after this check
			// elemBytes cannot overflow, and the range test below is
			// written so that off+bytes is never computed past the heap
			// (heapSize and off are >= 0, so heapSize-off cannot
			// overflow; an offset past the heap makes it negative, and
			// bytes >= 1 for n > 0).
			var bytes int64
			if d.n <= math.MaxInt64/16 {
				bytes = elemBytes(d.n)
			}
			if d.n > math.MaxInt64/16 || bytes > heapSize-d.off {
				return fmt.Errorf("fits: column %d row %d: VLA descriptor (%d elements, offset %d) lies outside the heap (%d bytes: PCOUNT minus the THEAP gap)", c.Index, r, d.n, d.off, heapSize)
			}
			if lo < 0 || d.off < lo {
				lo = d.off
			}
			if e := d.off + bytes; e > hi {
				hi = e
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	var heap []byte
	if lo >= 0 {
		heap = make([]byte, hi-lo)
		if err := t.rec.file.br.ReadRange(heapStart+lo, heap); err != nil {
			return fmt.Errorf("fits: column %d: read VLA heap: %w", c.Index, err)
		}
	}
	for r, d := range descs {
		var payload []byte
		if d.n > 0 {
			payload = heap[d.off-lo : d.off-lo+elemBytes(d.n)]
		}
		if err := fn(int64(r), d.n, payload); err != nil {
			return err
		}
	}
	return nil
}

func fixedElemBytes(size int) func(int64) int64 {
	return func(n int64) int64 { return n * int64(size) }
}

// ReadVarColumn reads a variable-length (P/Q) numeric column (element type
// B, I, J, K, E or D) into a VarColumn[T]. TSCALn/TZEROn/TNULLn are applied
// and checked per element exactly as for fixed columns (see ReadColumn):
// integer T must hold every physical value exactly, and a null element is
// an ErrNullValue error for integer T (NaN for float T); use
// ReadVarColumnMasked to get the null flags. Logical, bit, character and
// complex VLAs have their own readers (ReadVarLogicalColumn,
// ReadVarBitColumn, ReadVarStringColumn, ReadVarComplexColumn).
//
// Descriptors may point anywhere in the heap, in any order, and may share
// payloads; the heap starts at THEAP (§7.3.5).
func ReadVarColumn[T Numeric](t *BinaryTableHDU, col int) (*VarColumn[T], error) {
	v, _, err := readVarNumeric[T](t, col, false)
	return v, err
}

// ReadVarColumnMasked is ReadVarColumn plus a VarColumn of validity flags
// with the same shape (false = null element).
func ReadVarColumnMasked[T Numeric](t *BinaryTableHDU, col int) (*VarColumn[T], *VarColumn[bool], error) {
	return readVarNumeric[T](t, col, true)
}

func readVarNumeric[T Numeric](t *BinaryTableHDU, col int, masked bool) (*VarColumn[T], *VarColumn[bool], error) {
	c, err := t.column(col)
	if err != nil {
		return nil, nil, err
	}
	if !isVarType(c.bin.Type) {
		return nil, nil, fmt.Errorf("fits: column %d is not variable-length (TFORM=%s)", col, c.TForm)
	}
	if !isNumericBinType(c.bin.VarType) {
		return nil, nil, wrongReader(c, "not a numeric VLA; "+readerFor(c))
	}
	kind := reflect.TypeFor[T]().Kind()
	es := c.bin.VarType.ElementSize()
	sc := columnScaling(c)
	fd := newFastDec(&sc, c.bin.VarType, kind)
	nrows := t.NumRows()
	out := &VarColumn[T]{offsets: make([]int64, nrows+1)}
	valid := &VarColumn[bool]{offsets: out.offsets}
	err = t.varCells(c, fixedElemBytes(es), func(r, n int64, payload []byte) error {
		out.offsets[r] = int64(len(out.values))
		for i := int64(0); i < n; i++ {
			if fd.ok {
				var x T
				if cellFast(&fd, c.bin.VarType, payload[i*int64(es):], &x) {
					out.values = append(out.values, x)
					valid.values = append(valid.values, true)
					continue
				}
			}
			v := sc.decodeElem(c.bin.VarType, payload[i*int64(es):])
			x, err := numTo[T](v, masked, kind)
			if err != nil {
				return fmt.Errorf("fits: column %d row %d element %d: %w", col, r, i, err)
			}
			out.values = append(out.values, x)
			valid.values = append(valid.values, !v.null)
		}
		return nil
	})
	if err != nil {
		return nil, nil, err
	}
	out.offsets[nrows] = int64(len(out.values))
	if out.values == nil {
		out.values = []T{}
		valid.values = []bool{}
	}
	return out, valid, nil
}

// ReadVarLogicalColumn reads a P/Q logical (L) VLA: values and validity
// flags (false = undefined, byte 0x00) with the same shape. A logical byte
// other than 'T', 'F' or 0 is an error.
func ReadVarLogicalColumn(t *BinaryTableHDU, col int) (values, valid *VarColumn[bool], err error) {
	c, err := t.column(col)
	if err != nil {
		return nil, nil, err
	}
	if !isVarType(c.bin.Type) || c.bin.VarType != tform.BinLogical {
		return nil, nil, wrongReader(c, "not a logical VLA; "+readerFor(c))
	}
	nrows := t.NumRows()
	values = &VarColumn[bool]{offsets: make([]int64, nrows+1), values: []bool{}}
	valid = &VarColumn[bool]{offsets: values.offsets, values: []bool{}}
	err = t.varCells(c, fixedElemBytes(1), func(r, n int64, payload []byte) error {
		values.offsets[r] = int64(len(values.values))
		for i := int64(0); i < n; i++ {
			v, ok, err := decodeLogical(payload[i])
			if err != nil {
				return fmt.Errorf("fits: column %d row %d element %d: %w", col, r, i, err)
			}
			values.values = append(values.values, v)
			valid.values = append(valid.values, ok)
		}
		return nil
	})
	if err != nil {
		return nil, nil, err
	}
	values.offsets[nrows] = int64(len(values.values))
	return values, valid, nil
}

// ReadVarBitColumn reads a P/Q bit (X) VLA. The descriptor count is the
// number of bits; the heap holds ceil(n/8) bytes per row, most significant
// bit first.
func ReadVarBitColumn(t *BinaryTableHDU, col int) (*VarColumn[bool], error) {
	c, err := t.column(col)
	if err != nil {
		return nil, err
	}
	if !isVarType(c.bin.Type) || c.bin.VarType != tform.BinBit {
		return nil, wrongReader(c, "not a bit VLA; "+readerFor(c))
	}
	nrows := t.NumRows()
	out := &VarColumn[bool]{offsets: make([]int64, nrows+1), values: []bool{}}
	err = t.varCells(c, func(n int64) int64 { return (n + 7) / 8 }, func(r, n int64, payload []byte) error {
		out.offsets[r] = int64(len(out.values))
		bitsRow := make([]bool, n)
		unpackBits(bitsRow, payload)
		out.values = append(out.values, bitsRow...)
		return nil
	})
	if err != nil {
		return nil, err
	}
	out.offsets[nrows] = int64(len(out.values))
	return out, nil
}

// ReadVarStringColumn reads a P/Q character (A) VLA: one string per row
// made of the row's n characters, cut at the first NUL, trailing blanks
// removed (as for fixed A columns).
func ReadVarStringColumn(t *BinaryTableHDU, col int) ([]string, error) {
	c, err := t.column(col)
	if err != nil {
		return nil, err
	}
	if !isVarType(c.bin.Type) || c.bin.VarType != tform.BinChar {
		return nil, wrongReader(c, "not a character VLA; "+readerFor(c))
	}
	if k := dimTail(c.Dim); k != 1 {
		// TDIM (w, d2, ...) would make each row an array of strings;
		// returning the row as one string cut at a NUL could drop data.
		return nil, wrongReader(c, fmt.Sprintf("TDIM %v makes each row an array of strings; string arrays in VLA columns are not supported", c.Dim))
	}
	out := make([]string, t.NumRows())
	err = t.varCells(c, fixedElemBytes(1), func(r, n int64, payload []byte) error {
		out[r] = cellString(payload)
		return nil
	})
	if err != nil {
		return nil, err
	}
	return out, nil
}

// ReadVarComplexColumn reads a P/Q complex (C or M) VLA. TSCALn/TZEROn on
// complex columns are not supported (error).
func ReadVarComplexColumn[T Complex](t *BinaryTableHDU, col int) (*VarColumn[T], error) {
	c, err := t.column(col)
	if err != nil {
		return nil, err
	}
	if !isVarType(c.bin.Type) || (c.bin.VarType != tform.BinComplex64 && c.bin.VarType != tform.BinComplex128) {
		return nil, wrongReader(c, "not a complex VLA; "+readerFor(c))
	}
	if !columnScaling(c).identity { // any TSCAL/TZERO literal other than exactly 1 / 0
		return nil, wrongReader(c, "TSCAL/TZERO on complex columns is not supported")
	}
	nrows := t.NumRows()
	out := &VarColumn[T]{offsets: make([]int64, nrows+1), values: []T{}}
	es := c.bin.VarType.ElementSize()
	err = t.varCells(c, fixedElemBytes(es), func(r, n int64, payload []byte) error {
		out.offsets[r] = int64(len(out.values))
		row := make([]T, n)
		decodeComplex(row, payload, c.bin.VarType)
		out.values = append(out.values, row...)
		return nil
	})
	if err != nil {
		return nil, err
	}
	out.offsets[nrows] = int64(len(out.values))
	return out, nil
}

// dimTail is the product of the TDIM axes after the first (1 without a
// TDIM): the number of strings in a character cell or row.
func dimTail(d []int64) int64 {
	k := int64(1)
	if len(d) > 1 {
		for _, x := range d[1:] {
			k *= x
		}
	}
	return k
}
