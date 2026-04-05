package fits

import (
	"fmt"
	"strconv"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// Columns returns the ordered column metadata for the binary table.
// The result is parsed lazily on first call and cached.
func (h *BinaryTableHDU) Columns() ([]Column, error) {
	if h.columnsOnce {
		return h.columns, h.columnsErr
	}
	h.columnsOnce = true
	hdr := h.Header()
	nf, err := hdr.Int(header.KeyTfields)
	if err != nil {
		h.columnsErr = fmt.Errorf("fits: BinaryTable missing TFIELDS: %w", err)
		return nil, h.columnsErr
	}
	cols := make([]Column, nf)
	var rowBytes int64
	for i := int64(1); i <= nf; i++ {
		idx := i
		col := Column{Index: int(idx)}
		// Required: TFORMi
		tf, err := hdr.String(fmt.Sprintf("TFORM%d", idx))
		if err != nil {
			h.columnsErr = fmt.Errorf("fits: missing TFORM%d: %w", idx, err)
			return nil, h.columnsErr
		}
		col.TForm = tf
		bin, err := tform.ParseBinary(tf)
		if err != nil {
			h.columnsErr = fmt.Errorf("fits: TFORM%d: %w", idx, err)
			return nil, h.columnsErr
		}
		col.bin = bin
		col.Repeat = bin.Repeat
		col.Type = binaryColumnType(bin)

		// Optional keywords.
		if s, err := hdr.String(fmt.Sprintf("TTYPE%d", idx)); err == nil {
			col.Name = s
		}
		if s, err := hdr.String(fmt.Sprintf("TUNIT%d", idx)); err == nil {
			col.Unit = s
		}
		if s, err := hdr.String(fmt.Sprintf("TDISP%d", idx)); err == nil {
			col.Display = s
		}
		if f, err := hdr.Float(fmt.Sprintf("TSCAL%d", idx)); err == nil {
			col.Scale = f
		} else {
			col.Scale = 1.0
		}
		if f, err := hdr.Float(fmt.Sprintf("TZERO%d", idx)); err == nil {
			col.Zero = f
		}
		if v, err := hdr.Int(fmt.Sprintf("TNULL%d", idx)); err == nil {
			col.Null = v
		}
		if s, err := hdr.String(fmt.Sprintf("TDIM%d", idx)); err == nil {
			d, err := tform.ParseDim(s)
			if err == nil {
				col.Dim = d
			}
		}

		// Byte layout: per row, each column consumes repeat * element_size
		// bytes. Bit columns (X) use ceil(repeat/8). P/Q descriptors use
		// fixed 8/16 bytes regardless of repeat.
		var colBytes int64
		switch bin.Type {
		case tform.BinBit:
			colBytes = (bin.Repeat + 7) / 8
		case tform.BinPArrayDesc32:
			colBytes = 8
		case tform.BinQArrayDesc64:
			colBytes = 16
		case tform.BinChar:
			colBytes = bin.Repeat
		default:
			colBytes = bin.Repeat * int64(bin.Type.ElementSize())
		}
		col.byteOffset = int(rowBytes)
		col.byteSize = int(colBytes)
		cols[idx-1] = col
		rowBytes += colBytes
	}
	// Cross-check against NAXIS1 (row byte width).
	if h.rec.naxis >= 1 && h.rec.shape[0] != rowBytes {
		h.columnsErr = fmt.Errorf("fits: column total %d != NAXIS1 %d", rowBytes, h.rec.shape[0])
		return nil, h.columnsErr
	}
	h.rowBytes = int(rowBytes)
	h.columns = cols
	return cols, nil
}

// ColumnByName returns the column whose TTYPE matches name, or ok=false.
// The comparison is case-sensitive.
func (h *BinaryTableHDU) ColumnByName(name string) (Column, bool) {
	cols, err := h.Columns()
	if err != nil {
		return Column{}, false
	}
	for _, c := range cols {
		if c.Name == name {
			return c, true
		}
	}
	return Column{}, false
}

// ColumnIndex returns the 1-based column index for name, or 0 if absent.
func (h *BinaryTableHDU) ColumnIndex(name string) int {
	c, ok := h.ColumnByName(name)
	if !ok {
		return 0
	}
	return c.Index
}

// binaryColumnType maps a tform.BinaryType to the exported ColumnType enum.
func binaryColumnType(b tform.BinaryForm) ColumnType {
	switch b.Type {
	case tform.BinLogical:
		return ColLogical
	case tform.BinBit:
		return ColBit
	case tform.BinUint8:
		return ColByte
	case tform.BinInt16:
		return ColInt16
	case tform.BinInt32:
		return ColInt32
	case tform.BinInt64:
		return ColInt64
	case tform.BinChar:
		return ColString
	case tform.BinFloat32:
		return ColFloat32
	case tform.BinFloat64:
		return ColFloat64
	case tform.BinComplex64:
		return ColComplex64
	case tform.BinComplex128:
		return ColComplex128
	case tform.BinPArrayDesc32, tform.BinQArrayDesc64:
		return ColVarArray
	}
	return ColByte
}

// RowBytes returns the raw bytes of one row from the table data area.
// The returned slice is a copy safe to retain.
func (h *BinaryTableHDU) RowBytes(row int64) ([]byte, error) {
	if _, err := h.Columns(); err != nil {
		return nil, err
	}
	nrows := h.NumRows()
	if row < 0 || row >= nrows {
		return nil, fmt.Errorf("fits: row %d out of range [0,%d)", row, nrows)
	}
	buf := make([]byte, h.rowBytes)
	off := h.rec.dataStart + row*int64(h.rowBytes)
	if err := h.rec.file.br.ReadRange(off, buf); err != nil {
		return nil, err
	}
	return buf, nil
}

// heapStart returns the byte offset of the VLA heap area for this table.
// The heap begins at THEAP bytes from the start of the data section; if
// THEAP is absent it defaults to NAXIS1*NAXIS2 (i.e. immediately after the
// fixed rows).
func (h *BinaryTableHDU) heapStart() int64 {
	hdr := h.Header()
	if v, err := hdr.Int("THEAP"); err == nil {
		return h.rec.dataStart + v
	}
	if len(h.rec.shape) >= 2 {
		return h.rec.dataStart + h.rec.shape[0]*h.rec.shape[1]
	}
	return h.rec.dataStart
}

// naxisnKey returns "NAXISn" — small helper local to this file.
func naxisnKey(n int) string { return "NAXIS" + strconv.Itoa(n) }

var _ = naxisnKey // reserved for future use
