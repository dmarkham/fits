package fits

import (
	"fmt"
	"math"
	"math/big"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// Columns returns the ordered column metadata for the binary table.
// The result is parsed lazily on first call and cached.
//
// Columns fails (and every read of the table with it) when a column's
// keywords are not usable: a missing or malformed TFORMn, a malformed
// TDIMn or one whose element count exceeds the repeat count, a repeat
// count or TDIM whose byte or element count overflows int64, TSCALn = 0, a
// TNULLn on an integer column that is not an integer, or column widths
// that do not add up to NAXIS1. TSCALn/TZEROn on L, X and A columns and
// TNULLn on non-integer columns have no meaning in the standard and are
// ignored, as cfitsio does.
func (h *BinaryTableHDU) Columns() ([]Column, error) {
	if h.columnsOnce {
		return h.columns, h.columnsErr
	}
	h.columnsOnce = true
	cols, rowBytes, err := parseBinaryColumns(h.Header())
	if err == nil && h.rec.naxis >= 1 && h.rec.shape[0] != rowBytes {
		err = fmt.Errorf("fits: column total %d != NAXIS1 %d", rowBytes, h.rec.shape[0])
	}
	if err != nil {
		h.columnsErr = err
		return nil, err
	}
	h.rowBytes = int(rowBytes)
	h.columns = cols
	return cols, nil
}

func parseBinaryColumns(hdr *header.Header) ([]Column, int64, error) {
	nf, err := hdr.Int(header.KeyTfields)
	if err != nil {
		return nil, 0, fmt.Errorf("fits: BinaryTable missing TFIELDS: %w", err)
	}
	if nf < 0 || nf > 999 {
		return nil, 0, fmt.Errorf("fits: TFIELDS = %d out of range [0,999]", nf)
	}
	cols := make([]Column, nf)
	var rowBytes int64
	for idx := int64(1); idx <= nf; idx++ {
		col := Column{Index: int(idx)}
		tf, err := hdr.String(fmt.Sprintf("TFORM%d", idx))
		if err != nil {
			return nil, 0, fmt.Errorf("fits: missing TFORM%d: %w", idx, err)
		}
		col.TForm = tf
		bin, err := tform.ParseBinary(tf)
		if err != nil {
			return nil, 0, fmt.Errorf("fits: TFORM%d: %w", idx, err)
		}
		col.bin = bin
		col.Repeat = bin.Repeat
		col.Width = int(bin.Width)
		col.Type = binaryColumnType(bin)

		if s, err := hdr.String(fmt.Sprintf("TTYPE%d", idx)); err == nil {
			col.Name = s
		}
		if s, err := hdr.String(fmt.Sprintf("TUNIT%d", idx)); err == nil {
			col.Unit = s
		}
		if s, err := hdr.String(fmt.Sprintf("TDISP%d", idx)); err == nil {
			col.Display = s
		}
		col.Scale = 1
		if err := readScaleZero(hdr, idx, &col); err != nil {
			return nil, 0, err
		}
		elem := bin.Type
		if isVarType(bin.Type) {
			elem = bin.VarType
		}
		if isIntType(elem) {
			key := fmt.Sprintf("TNULL%d", idx)
			if c, ok := hdr.Get(key); ok {
				v, isInt := c.Value.(int64)
				if !isInt || c.Type != header.TypeInt {
					return nil, 0, fmt.Errorf("fits: %s = %v is not an integer", key, c.Value)
				}
				col.Null = v
				col.hasNull = true
				col.nullInt = v
			}
		}
		if c, ok := hdr.Get(fmt.Sprintf("TDIM%d", idx)); ok {
			s, isStr := c.Value.(string)
			if !isStr {
				return nil, 0, fmt.Errorf("fits: TDIM%d is not a string", idx)
			}
			d, err := tform.ParseDim(s)
			if err != nil {
				return nil, 0, fmt.Errorf("fits: TDIM%d: %w", idx, err)
			}
			n := int64(1)
			for _, ax := range d {
				if ax < 0 || (ax > 0 && n > math.MaxInt64/ax) {
					return nil, 0, fmt.Errorf("fits: TDIM%d %s: the element count overflows int64", idx, s)
				}
				n *= ax
			}
			if !isVarType(bin.Type) {
				if n > bin.Repeat {
					return nil, 0, fmt.Errorf("fits: TDIM%d %s has %d elements, more than the repeat count %d of TFORM%d %q", idx, s, n, bin.Repeat, idx, tf)
				}
			}
			col.Dim = d
		}
		if bin.Type == tform.BinChar && bin.Width > 0 && bin.Width > bin.Repeat {
			return nil, 0, fmt.Errorf("fits: TFORM%d %q: substring width exceeds repeat", idx, tf)
		}

		// Byte layout: per row, each column consumes repeat * element_size
		// bytes. Bit columns (X) use ceil(repeat/8). P/Q descriptors use
		// fixed 8/16 bytes (times a repeat of 0 or 1).
		var colBytes int64
		switch bin.Type {
		case tform.BinBit:
			colBytes = (bin.Repeat + 7) / 8
		case tform.BinChar:
			colBytes = bin.Repeat
		default:
			es := int64(bin.Type.ElementSize())
			if bin.Repeat > math.MaxInt64/es {
				return nil, 0, fmt.Errorf("fits: TFORM%d %q: %d elements of %d bytes overflow int64", idx, tf, bin.Repeat, es)
			}
			colBytes = bin.Repeat * es
		}
		if colBytes > math.MaxInt64-rowBytes {
			return nil, 0, fmt.Errorf("fits: TFORM%d %q: the column widths overflow int64", idx, tf)
		}
		col.byteOffset = int(rowBytes)
		col.byteSize = int(colBytes)
		cols[idx-1] = col
		rowBytes += colBytes
	}
	return cols, rowBytes, nil
}

// readScaleZero fills col.Scale/col.Zero from TSCALn/TZEROn, and keeps
// both literals exactly (col.scaleLit, col.zeroLit) for the exact paths.
func readScaleZero(hdr *header.Header, idx int64, col *Column) error {
	if c, ok := hdr.Get(fmt.Sprintf("TSCAL%d", idx)); ok {
		v, err := cardNumber(c)
		if err != nil {
			return fmt.Errorf("fits: TSCAL%d: %w", idx, err)
		}
		if v == 0 {
			return fmt.Errorf("fits: TSCAL%d = 0", idx)
		}
		col.Scale = v
		col.scaleLit = exactCardValue(c)
	}
	if c, ok := hdr.Get(fmt.Sprintf("TZERO%d", idx)); ok {
		v, err := cardNumber(c)
		if err != nil {
			return fmt.Errorf("fits: TZERO%d: %w", idx, err)
		}
		col.Zero = v
		col.zeroLit = exactCardValue(c)
		// A float literal that is not 2^63 but whose nearest double is
		// (9.223372036854776E+18): cfitsio and astropy read it as the
		// uint64 convention, the literal itself says stored +
		// 9223372036854776000. Integer reads refuse to guess.
		col.zeroNear2p63 = c.Type == header.TypeFloat && v == 9223372036854775808.0 && col.zeroLit.Cmp(rat2p63) != 0
	}
	return nil
}

var rat2p63 = new(big.Rat).SetInt(new(big.Int).Lsh(big.NewInt(1), 63))

// exactCardValue returns a numeric card's value exactly: an integer
// literal as is (of any size), a float literal ("4611686018427387905.0",
// "1.0E3", "0.1") as the exact decimal of its card text. A card built in
// memory (no card text) is the float64 value itself. Infinite or NaN
// values (not FITS numbers) give nil.
func exactCardValue(c header.Card) *big.Rat {
	switch v := c.Value.(type) {
	case int64:
		return new(big.Rat).SetInt64(v)
	case *big.Int:
		return new(big.Rat).SetInt(v)
	case float64:
		if lit, ok := numericLiteral(c); ok {
			if r, ok := new(big.Rat).SetString(lit); ok {
				return r
			}
		}
		if math.IsInf(v, 0) || math.IsNaN(v) {
			return nil
		}
		return new(big.Rat).SetFloat64(v)
	}
	return nil
}

// numericLiteral extracts the value text of a numeric card from its raw
// 80 bytes ("KEYWORD = value / comment"), with a D exponent written as E.
func numericLiteral(c header.Card) (string, bool) {
	raw := string(c.Raw[:])
	if len(raw) < 10 || raw[8:10] != "= " {
		return "", false
	}
	v := raw[10:]
	if i := strings.IndexByte(v, '/'); i >= 0 {
		v = v[:i]
	}
	v = strings.TrimSpace(v)
	if v == "" {
		return "", false
	}
	return strings.Map(func(r rune) rune {
		if r == 'D' || r == 'd' {
			return 'E'
		}
		return r
	}, v), true
}

// cardNumber returns a numeric card value as float64. An integer literal
// outside int64 (header.TypeBigInt, e.g. cfitsio's uint64 TZERO
// 9223372036854775808) is rounded to the nearest double, as cfitsio reads
// it (ffgkyd on a TDOUBLE).
func cardNumber(c header.Card) (float64, error) {
	switch v := c.Value.(type) {
	case int64:
		return float64(v), nil
	case float64:
		return v, nil
	case *big.Int:
		f, _ := new(big.Float).SetInt(v).Float64()
		return f, nil
	}
	return 0, fmt.Errorf("value %v is not a number", c.Value)
}

func isVarType(b tform.BinaryType) bool {
	return b == tform.BinPArrayDesc32 || b == tform.BinQArrayDesc64
}

func isIntType(b tform.BinaryType) bool {
	switch b {
	case tform.BinUint8, tform.BinInt16, tform.BinInt32, tform.BinInt64:
		return true
	}
	return false
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

// heapRange returns the absolute file offset of the VLA heap and its size
// in bytes. The heap begins THEAP bytes after the start of the data unit
// (default NAXIS1*NAXIS2, i.e. right after the fixed rows) and ends at
// NAXIS1*NAXIS2 + PCOUNT (§7.3.5). A THEAP outside that range is an error.
func (h *BinaryTableHDU) heapRange() (start, size int64, err error) {
	var fixed int64
	if len(h.rec.shape) >= 2 {
		fixed = h.rec.shape[0] * h.rec.shape[1]
	}
	end := fixed + h.rec.pcount
	theap := fixed
	if c, ok := h.Header().Get("THEAP"); ok {
		v, isInt := c.Value.(int64)
		if !isInt {
			return 0, 0, fmt.Errorf("fits: THEAP = %v is not an integer", c.Value)
		}
		theap = v
	}
	if theap < fixed || theap > end {
		return 0, 0, fmt.Errorf("fits: THEAP = %d outside the data unit's heap area [%d,%d] (NAXIS1*NAXIS2, NAXIS1*NAXIS2+PCOUNT)", theap, fixed, end)
	}
	return h.rec.dataStart + theap, end - theap, nil
}
