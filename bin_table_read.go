package fits

import (
	"fmt"
	"math"
	"unsafe"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/tform"
)

// ColumnValue is the constraint for the typed ReadColumn API. It accepts
// every Numeric type plus string and bool for Char/Logical columns.
type ColumnValue interface {
	~uint8 | ~int8 | ~int16 | ~uint16 | ~int32 | ~uint32 |
		~int64 | ~uint64 | ~float32 | ~float64 |
		~string | ~bool
}

// ReadColumn reads a scalar column (Repeat==1) across all rows and returns
// one value per row as a []T. For Char columns (TFORM "A") use T=string.
// For Logical columns use T=bool. For numeric columns use any Numeric T
// that can hold the column type (same rules as CanConvert on images).
//
// TSCAL/TZERO are applied to numeric columns (BSCALE/BZERO analogue).
func ReadColumn[T ColumnValue](t *BinaryTableHDU, col int) ([]T, error) {
	cols, err := t.Columns()
	if err != nil {
		return nil, err
	}
	if col < 1 || col > len(cols) {
		return nil, fmt.Errorf("fits: column %d out of range [1,%d]", col, len(cols))
	}
	c := cols[col-1]
	if c.bin.Type == tform.BinPArrayDesc32 || c.bin.Type == tform.BinQArrayDesc64 {
		return nil, fmt.Errorf("fits: column %d is variable-length; use ReadVarColumn", col)
	}
	if c.Repeat != 1 && c.bin.Type != tform.BinChar {
		return nil, fmt.Errorf("fits: column %d has repeat %d; use ReadVectorColumn", col, c.Repeat)
	}
	nrows := t.NumRows()
	out := make([]T, nrows)

	// Fast path: read the whole column in one pass by slurping all rows.
	// For very wide tables this is wasteful; we bias toward simplicity.
	rowBytes := t.rowBytes
	total := int64(rowBytes) * nrows
	buf := make([]byte, total)
	if total > 0 {
		if err := t.rec.file.br.ReadRange(t.rec.dataStart, buf); err != nil {
			return nil, err
		}
	}

	for r := int64(0); r < nrows; r++ {
		cellBase := int(r)*rowBytes + c.byteOffset
		cell := buf[cellBase : cellBase+c.byteSize]
		v, err := decodeScalarCell[T](cell, c)
		if err != nil {
			return nil, fmt.Errorf("fits: column %d row %d: %w", col, r, err)
		}
		out[r] = v
	}
	return out, nil
}

// decodeScalarCell decodes a single-cell binary table entry. For Char
// columns the entire c.Repeat bytes form the string (trailing spaces
// trimmed). For numeric columns TSCAL/TZERO are applied.
func decodeScalarCell[T ColumnValue](cell []byte, c Column) (T, error) {
	var zero T
	switch c.bin.Type {
	case tform.BinChar:
		// Trim trailing NUL and spaces per FITS §7.3.3.1.
		end := len(cell)
		for end > 0 && (cell[end-1] == ' ' || cell[end-1] == 0) {
			end--
		}
		s := string(cell[:end])
		var out any = s
		if _, ok := out.(T); !ok {
			return zero, fmt.Errorf("TFORM 'A' requires T=string")
		}
		return any(s).(T), nil
	case tform.BinLogical:
		b := cell[0] == 'T'
		if _, ok := any(b).(T); !ok {
			return zero, fmt.Errorf("TFORM 'L' requires T=bool")
		}
		return any(b).(T), nil
	case tform.BinUint8:
		f := float64(cell[0])
		return applyScaleAndCast[T](f, c), nil
	case tform.BinInt16:
		f := float64(bigendian.Int16(cell))
		return applyScaleAndCast[T](f, c), nil
	case tform.BinInt32:
		f := float64(bigendian.Int32(cell))
		return applyScaleAndCast[T](f, c), nil
	case tform.BinInt64:
		f := float64(bigendian.Int64(cell))
		return applyScaleAndCast[T](f, c), nil
	case tform.BinFloat32:
		f := float64(bigendian.Float32(cell))
		return applyScaleAndCast[T](f, c), nil
	case tform.BinFloat64:
		f := bigendian.Float64(cell)
		return applyScaleAndCast[T](f, c), nil
	}
	return zero, fmt.Errorf("fits: unsupported column type %v for scalar read", c.Type)
}

// applyScaleAndCast applies TSCAL/TZERO and casts to T.
func applyScaleAndCast[T ColumnValue](v float64, c Column) T {
	if c.Scale != 0 && c.Scale != 1 {
		v *= c.Scale
	} else if c.Scale == 0 {
		// Missing TSCAL defaults to 1 — decoded as 0 indicates bug; fall
		// through treating as 1.
	}
	v += c.Zero
	return columnCast[T](v)
}

// columnCast casts a float64 to T. Mirrors numericFromFloat but for
// ColumnValue (which additionally includes string and bool — unused here).
func columnCast[T ColumnValue](v float64) T {
	var zero T
	switch any(zero).(type) {
	case float64:
		return *(*T)(unsafe.Pointer(&v))
	case float32:
		f := float32(v)
		return *(*T)(unsafe.Pointer(&f))
	case int8:
		i := int8(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint8:
		i := uint8(v)
		return *(*T)(unsafe.Pointer(&i))
	case int16:
		i := int16(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint16:
		i := uint16(v)
		return *(*T)(unsafe.Pointer(&i))
	case int32:
		i := int32(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint32:
		i := uint32(v)
		return *(*T)(unsafe.Pointer(&i))
	case int64:
		i := int64(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint64:
		i := uint64(v)
		return *(*T)(unsafe.Pointer(&i))
	}
	return zero
}

// ReadVectorColumn reads a fixed-width vector column (Repeat>1) and returns
// a flat slice of length nrows*repeat in row-major order.
func ReadVectorColumn[T Numeric](t *BinaryTableHDU, col int) ([]T, error) {
	cols, err := t.Columns()
	if err != nil {
		return nil, err
	}
	if col < 1 || col > len(cols) {
		return nil, fmt.Errorf("fits: column %d out of range [1,%d]", col, len(cols))
	}
	c := cols[col-1]
	if c.bin.Type == tform.BinPArrayDesc32 || c.bin.Type == tform.BinQArrayDesc64 {
		return nil, fmt.Errorf("fits: column %d is variable-length; use ReadVarColumn", col)
	}
	if c.bin.Type == tform.BinChar || c.bin.Type == tform.BinLogical {
		return nil, fmt.Errorf("fits: column %d has non-numeric type; use ReadColumn/RowBytes", col)
	}
	nrows := t.NumRows()
	repeat := c.Repeat
	total := nrows * repeat
	out := make([]T, total)
	buf := make([]byte, int64(t.rowBytes)*nrows)
	if err := t.rec.file.br.ReadRange(t.rec.dataStart, buf); err != nil {
		return nil, err
	}
	elemSize := int(c.bin.Type.ElementSize())
	for r := int64(0); r < nrows; r++ {
		cellBase := int(r)*t.rowBytes + c.byteOffset
		for k := int64(0); k < repeat; k++ {
			src := buf[cellBase+int(k)*elemSize : cellBase+int(k)*elemSize+elemSize]
			var f float64
			switch c.bin.Type {
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
			out[r*repeat+k] = numericFromFloat[T](f)
		}
	}
	return out, nil
}

// Silence unused-import guard in case future code drops one of these.
var _ = math.Pi
