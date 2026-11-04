package fits

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// ASCIIColumnData describes one column of a new ASCII table (§7.2).
// Exactly one Data* field must be set:
//
//	DataString  -> Aw   (w = TForm's width, else the longest string)
//	DataInt64   -> Iw   (default I20, wide enough for every int64)
//	DataFloat64 -> Fw.d, Ew.d or Dw.d (default E24.16: 17 significant
//	               digits, so every float64 reads back exactly)
//
// TForm may set any ASCII format of the matching kind (e.g. "I6",
// "F10.3", "D25.17", "A8"); a value that does not fit the field width is
// an error, never a truncated or rounded-away field beyond what the
// format's d digits state. Floats are written with an explicit decimal
// point and an E (or D) exponent per the TFORM letter; NaN and ±Inf
// cannot be written (error) unless the element is marked null.
//
// Null sets TNULLn; Valid marks null rows (false = write the TNULL string,
// left-justified). A null row without Null is an error, and so is a
// non-null value whose field text equals the TNULL string.
type ASCIIColumnData struct {
	Name  string
	Unit  string
	TForm string
	Null  string // TNULLn; "" = no TNULL

	DataString  []string
	DataInt64   []int64
	DataFloat64 []float64

	Valid []bool
}

// AppendASCIITable appends a new ASCII TABLE extension. Fields are laid
// out left to right with one blank between neighbouring fields (TBCOLn
// accordingly); NAXIS1 is the total row width. hdr may carry extra cards;
// the structural and column keywords (XTENSION, BITPIX, NAXIS*, PCOUNT,
// GCOUNT, TFIELDS, TBCOLn, TFORMn, TTYPEn, TUNITn, TNULLn, TSCALn, TZEROn,
// TDIMn, THEAP) are written by the writer: a TTYPEn/TUNITn/TNULLn in hdr
// is used when the column leaves the field empty and is an error when it
// conflicts; any other column keyword there is an error. TSCALn/TZEROn
// are not supported by this writer.
func AppendASCIITable(f *File, hdr *header.Header, cols []ASCIIColumnData) (*ASCIITableHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	if len(cols) == 0 {
		return nil, fmt.Errorf("fits: AppendASCIITable: no columns")
	}
	if len(cols) > 999 {
		return nil, fmt.Errorf("fits: AppendASCIITable: %d columns, FITS allows 999", len(cols))
	}
	cols = append([]ASCIIColumnData(nil), cols...)
	extra, err := mergeASCIIKeys(hdr, cols)
	if err != nil {
		return nil, err
	}
	type spec struct {
		form  tform.ASCIIForm
		tbcol int
		n     int
	}
	specs := make([]spec, len(cols))
	nrows := -1
	pos := 1
	for i := range cols {
		c := &cols[i]
		form, n, err := asciiForm(c)
		if err != nil {
			return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q): %w", i+1, c.Name, err)
		}
		if nrows >= 0 && n != nrows {
			return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q): %d rows, earlier columns have %d", i+1, c.Name, n, nrows)
		}
		nrows = n
		if c.Valid != nil && len(c.Valid) != n {
			return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q): Valid has %d flags, data %d", i+1, c.Name, len(c.Valid), n)
		}
		if len(c.Null) > form.Width {
			return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q): TNULL %q is wider than the field (%d)", i+1, c.Name, c.Null, form.Width)
		}
		if err := printableASCII(c.Null); err != nil {
			return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q): TNULL: %w", i+1, c.Name, err)
		}
		if i > 0 {
			pos++ // one blank between fields
		}
		specs[i] = spec{form: form, tbcol: pos, n: n}
		pos += form.Width
	}
	rowWidth := pos - 1
	data := make([]byte, rowWidth*nrows)
	for i := range data {
		data[i] = ' '
	}
	for i := range cols {
		c := &cols[i]
		sp := specs[i]
		for r := 0; r < nrows; r++ {
			var field string
			if c.Valid != nil && !c.Valid[r] {
				if c.Null == "" {
					return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q) row %d: null row but no Null (TNULL) set", i+1, c.Name, r)
				}
				field = c.Null + strings.Repeat(" ", sp.form.Width-len(c.Null))
			} else {
				field, err = formatASCIIField(c, sp.form, r)
				if err != nil {
					return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q) row %d: %w", i+1, c.Name, r, err)
				}
				if c.Null != "" && isASCIINull(field, c.Null) {
					return nil, fmt.Errorf("fits: AppendASCIITable: column %d (%q) row %d: field %q equals TNULL and would read back as null", i+1, c.Name, r, field)
				}
			}
			copy(data[r*rowWidth+sp.tbcol-1:], field)
		}
	}

	hb := newHeaderBuilder()
	// Fixed-format string values must be at least 8 characters (§4.2.1.1),
	// which fitsverify enforces for XTENSION.
	hb.set(header.KeyXtension, "TABLE   ", "ASCII table extension")
	hb.set(header.KeyBitpix, int64(8), "")
	hb.set(header.KeyNaxis, int64(2), "")
	hb.set("NAXIS1", int64(rowWidth), "width of table row in characters")
	hb.set("NAXIS2", int64(nrows), "number of rows")
	hb.set(header.KeyPcount, int64(0), "")
	hb.set(header.KeyGcount, int64(1), "")
	hb.set(header.KeyTfields, int64(len(cols)), "number of columns")
	for i := range cols {
		k := strconv.Itoa(i + 1)
		c := &cols[i]
		hb.set("TBCOL"+k, int64(specs[i].tbcol), "")
		hb.set("TFORM"+k, specs[i].form.Raw, "")
		if c.Name != "" {
			hb.set("TTYPE"+k, c.Name, "")
		}
		if c.Unit != "" {
			hb.set("TUNIT"+k, c.Unit, "")
		}
		if c.Null != "" {
			hb.set("TNULL"+k, c.Null, "")
		}
	}
	for _, c := range extra {
		hb.add(c)
	}
	if hb.err != nil {
		return nil, fmt.Errorf("fits: AppendASCIITable: %w", hb.err)
	}
	h := hb.h
	rec, err := appendTable(f, h, kindASCIITable, "TABLE", int64(rowWidth), int64(nrows), int64(len(cols)), data)
	if err != nil {
		return nil, err
	}
	return &ASCIITableHDU{rec: rec}, nil
}

// asciiForm resolves the TFORM of a column and its row count.
func asciiForm(c *ASCIIColumnData) (tform.ASCIIForm, int, error) {
	set := 0
	var n int
	var def string
	var kinds string
	if c.DataString != nil {
		set++
		n = len(c.DataString)
		w := 1
		for r, s := range c.DataString {
			if err := printableASCII(s); err != nil {
				return tform.ASCIIForm{}, 0, fmt.Errorf("row %d: %w", r, err)
			}
			w = max(w, len(s))
		}
		def = "A" + strconv.Itoa(w)
		kinds = "A"
	}
	if c.DataInt64 != nil {
		set++
		n = len(c.DataInt64)
		def = "I20"
		kinds = "I"
	}
	if c.DataFloat64 != nil {
		set++
		n = len(c.DataFloat64)
		def = "E24.16"
		kinds = "FED"
	}
	if set != 1 {
		return tform.ASCIIForm{}, 0, fmt.Errorf("%d Data fields set, want exactly one", set)
	}
	tf := c.TForm
	if tf == "" {
		tf = def
	}
	form, err := tform.ParseASCII(tf)
	if err != nil {
		return form, 0, err
	}
	if !strings.ContainsRune(kinds, rune(form.Type)) {
		return form, 0, fmt.Errorf("TForm %q does not match the data (want one of %s)", tf, kinds)
	}
	form.Raw = strings.TrimSpace(tf)
	return form, n, nil
}

// formatASCIIField renders row r right-justified (numbers) or
// left-justified (strings) in exactly form.Width characters.
func formatASCIIField(c *ASCIIColumnData, form tform.ASCIIForm, r int) (string, error) {
	var s string
	switch form.Type {
	case tform.AsciiChar:
		v := c.DataString[r]
		if len(v) > form.Width {
			return "", fmt.Errorf("string of %d bytes does not fit %s", len(v), form.Raw)
		}
		return v + strings.Repeat(" ", form.Width-len(v)), nil
	case tform.AsciiInt:
		s = strconv.FormatInt(c.DataInt64[r], 10)
	default:
		v := c.DataFloat64[r]
		if math.IsNaN(v) || math.IsInf(v, 0) {
			return "", fmt.Errorf("%v cannot be written in an ASCII table (mark it null)", v)
		}
		switch form.Type {
		case tform.AsciiFloatF:
			s = strconv.FormatFloat(v, 'f', form.Fraction, 64)
			if form.Fraction == 0 {
				s += "." // keep the decimal point explicit
			}
		default:
			s = strconv.FormatFloat(v, 'E', form.Fraction, 64)
			if form.Fraction == 0 {
				// "1E+05" -> "1.E+05": explicit decimal point.
				s = strings.Replace(s, "E", ".E", 1)
			}
			if form.Type == tform.AsciiFloatD {
				s = strings.Replace(s, "E", "D", 1)
			}
		}
	}
	if len(s) > form.Width {
		return "", fmt.Errorf("value %s does not fit %s", s, form.Raw)
	}
	return strings.Repeat(" ", form.Width-len(s)) + s, nil
}

func mergeASCIIKeys(hdr *header.Header, cols []ASCIIColumnData) ([]header.Card, error) {
	if hdr == nil {
		return nil, nil
	}
	var extra []header.Card
	for _, c := range hdr.Cards() {
		if isStructuralTableKey(c.Key) {
			continue
		}
		root, n, ok := indexedKey(c.Key)
		switch root {
		case "TTYPE", "TUNIT", "TNULL":
		case "TFORM", "TSCAL", "TZERO", "TDIM", "TDISP":
			if ok {
				return nil, fmt.Errorf("fits: AppendASCIITable: hdr %s is not supported (set ASCIIColumnData.TForm; TSCAL/TZERO/TDIM/TDISP are not written)", c.Key)
			}
		default:
			ok = false
		}
		if !ok {
			extra = append(extra, c)
			continue
		}
		if n < 1 || n > len(cols) {
			return nil, fmt.Errorf("fits: AppendASCIITable: hdr has %s but the table has %d columns", c.Key, len(cols))
		}
		v, isStr := c.Value.(string)
		if !isStr {
			return nil, fmt.Errorf("fits: AppendASCIITable: hdr %s is not a string", c.Key)
		}
		col := &cols[n-1]
		dst := map[string]*string{"TTYPE": &col.Name, "TUNIT": &col.Unit, "TNULL": &col.Null}[root]
		if *dst != "" && *dst != v {
			return nil, fmt.Errorf("fits: AppendASCIITable: hdr %s = %q conflicts with column %d's %q", c.Key, v, n, *dst)
		}
		*dst = v
	}
	return extra, nil
}
