package fits

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// Columns returns the column metadata for the ASCII table.
func (h *ASCIITableHDU) Columns() ([]Column, error) {
	if h.columnsOnce {
		return h.columns, h.columnsErr
	}
	h.columnsOnce = true
	hdr := h.Header()
	nf, err := hdr.Int(header.KeyTfields)
	if err != nil {
		h.columnsErr = fmt.Errorf("fits: ASCIITable missing TFIELDS: %w", err)
		return nil, h.columnsErr
	}
	cols := make([]Column, nf)
	for i := int64(1); i <= nf; i++ {
		idx := i
		col := Column{Index: int(idx)}
		tfStr, err := hdr.String(fmt.Sprintf("TFORM%d", idx))
		if err != nil {
			h.columnsErr = fmt.Errorf("fits: missing TFORM%d: %w", idx, err)
			return nil, h.columnsErr
		}
		col.TForm = tfStr
		af, err := tform.ParseASCII(tfStr)
		if err != nil {
			h.columnsErr = fmt.Errorf("fits: TFORM%d: %w", idx, err)
			return nil, h.columnsErr
		}
		// TBCOLi — 1-based column offset of the field in the row.
		bcol, err := hdr.Int(fmt.Sprintf("TBCOL%d", idx))
		if err != nil {
			h.columnsErr = fmt.Errorf("fits: missing TBCOL%d: %w", idx, err)
			return nil, h.columnsErr
		}
		col.byteOffset = int(bcol - 1)
		col.byteSize = af.Width
		col.Repeat = 1
		switch af.Type {
		case tform.AsciiChar:
			col.Type = ColString
		case tform.AsciiInt:
			col.Type = ColInt64
		case tform.AsciiFloatF, tform.AsciiFloatE, tform.AsciiFloatD:
			col.Type = ColFloat64
		}
		if s, err := hdr.String(fmt.Sprintf("TTYPE%d", idx)); err == nil {
			col.Name = s
		}
		if s, err := hdr.String(fmt.Sprintf("TUNIT%d", idx)); err == nil {
			col.Unit = s
		}
		if f, err := hdr.Float(fmt.Sprintf("TSCAL%d", idx)); err == nil {
			col.Scale = f
		} else {
			col.Scale = 1.0
		}
		if f, err := hdr.Float(fmt.Sprintf("TZERO%d", idx)); err == nil {
			col.Zero = f
		}
		cols[idx-1] = col
	}
	h.columns = cols
	return cols, nil
}

// ColumnByName returns the column whose TTYPE matches name.
func (h *ASCIITableHDU) ColumnByName(name string) (Column, bool) {
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

// ReadASCIIColumn reads a scalar column from an ASCII table. For string
// columns (TFORM 'A') T must be ~string; for integer and float columns T
// must be a Numeric type.
func ReadASCIIColumn[T ColumnValue](t *ASCIITableHDU, col int) ([]T, error) {
	cols, err := t.Columns()
	if err != nil {
		return nil, err
	}
	if col < 1 || col > len(cols) {
		return nil, fmt.Errorf("fits: column %d out of range [1,%d]", col, len(cols))
	}
	c := cols[col-1]
	if len(t.rec.shape) < 2 {
		return nil, fmt.Errorf("fits: ASCII table missing NAXIS2")
	}
	nrows := t.rec.shape[1]
	rowWidth := int(t.rec.shape[0])
	buf := make([]byte, int64(rowWidth)*nrows)
	if err := t.rec.file.br.ReadRange(t.rec.dataStart, buf); err != nil {
		return nil, err
	}
	out := make([]T, nrows)
	for r := int64(0); r < nrows; r++ {
		field := string(buf[int(r)*rowWidth+c.byteOffset : int(r)*rowWidth+c.byteOffset+c.byteSize])
		v, err := decodeASCIIField[T](field, c)
		if err != nil {
			return nil, fmt.Errorf("fits: col %d row %d: %w", col, r, err)
		}
		out[r] = v
	}
	return out, nil
}

// decodeASCIIField parses one ASCII-table cell into T. Leading/trailing
// whitespace is stripped for numeric fields; strings preserve trailing
// whitespace trimmed of pad.
func decodeASCIIField[T ColumnValue](field string, c Column) (T, error) {
	var zero T
	trimmed := strings.TrimSpace(field)
	switch c.Type {
	case ColString:
		s := strings.TrimRight(field, " ")
		if _, ok := any(s).(T); !ok {
			return zero, fmt.Errorf("string column requires T=string")
		}
		return any(s).(T), nil
	case ColInt64:
		n, err := strconv.ParseInt(trimmed, 10, 64)
		if err != nil {
			return zero, err
		}
		f := float64(n)*safeScale(c.Scale) + c.Zero
		return columnCast[T](f), nil
	case ColFloat64:
		// Accept both E and D exponent forms.
		s := trimmed
		for i := 0; i < len(s); i++ {
			if s[i] == 'd' || s[i] == 'D' {
				s = s[:i] + "e" + s[i+1:]
				break
			}
		}
		n, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return zero, err
		}
		n = n*safeScale(c.Scale) + c.Zero
		return columnCast[T](n), nil
	}
	return zero, fmt.Errorf("unsupported ASCII column type %v", c.Type)
}

func safeScale(s float64) float64 {
	if s == 0 {
		return 1
	}
	return s
}
