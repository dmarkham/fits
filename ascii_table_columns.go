package fits

import (
	"fmt"
	"math"
	"math/big"
	"reflect"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// Columns returns the column metadata for the ASCII table.
//
// It fails when a TFORMn or TBCOLn is missing or malformed, a field does
// not fit inside NAXIS1, or TSCALn = 0. Column.Null holds TNULLn as a
// string; Column.Width is the field width.
func (h *ASCIITableHDU) Columns() ([]Column, error) {
	if h.columnsOnce {
		return h.columns, h.columnsErr
	}
	h.columnsOnce = true
	cols, err := h.parseColumns()
	if err != nil {
		h.columnsErr = err
		return nil, err
	}
	h.columns = cols
	return cols, nil
}

func (h *ASCIITableHDU) parseColumns() ([]Column, error) {
	hdr := h.Header()
	nf, err := hdr.Int(header.KeyTfields)
	if err != nil {
		return nil, fmt.Errorf("fits: ASCIITable missing TFIELDS: %w", err)
	}
	if nf < 0 || nf > 999 {
		return nil, fmt.Errorf("fits: TFIELDS = %d out of range [0,999]", nf)
	}
	var rowWidth int64
	if len(h.rec.shape) >= 1 {
		rowWidth = h.rec.shape[0]
	}
	cols := make([]Column, nf)
	for idx := int64(1); idx <= nf; idx++ {
		col := Column{Index: int(idx)}
		tfStr, err := hdr.String(fmt.Sprintf("TFORM%d", idx))
		if err != nil {
			return nil, fmt.Errorf("fits: missing TFORM%d: %w", idx, err)
		}
		col.TForm = tfStr
		af, err := tform.ParseASCII(tfStr)
		if err != nil {
			return nil, fmt.Errorf("fits: TFORM%d: %w", idx, err)
		}
		col.ascii = af
		bcol, err := hdr.Int(fmt.Sprintf("TBCOL%d", idx))
		if err != nil {
			return nil, fmt.Errorf("fits: missing TBCOL%d: %w", idx, err)
		}
		// Written so that nothing wraps around int64: 0 <= bcol-1 < rowWidth.
		if bcol < 1 || bcol > rowWidth || int64(af.Width) > rowWidth-(bcol-1) {
			return nil, fmt.Errorf("fits: TBCOL%d = %d with width %d does not fit NAXIS1 = %d", idx, bcol, af.Width, rowWidth)
		}
		col.byteOffset = int(bcol - 1)
		col.byteSize = af.Width
		col.Width = af.Width
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
		if s, err := hdr.String(fmt.Sprintf("TDISP%d", idx)); err == nil {
			col.Display = s
		}
		col.Scale = 1
		if af.Type != tform.AsciiChar {
			if err := readScaleZero(hdr, idx, &col); err != nil {
				return nil, err
			}
		}
		if c, ok := hdr.Get(fmt.Sprintf("TNULL%d", idx)); ok {
			s, isStr := c.Value.(string)
			if !isStr {
				return nil, fmt.Errorf("fits: TNULL%d of an ASCII table must be a string, got %v", idx, c.Value)
			}
			col.Null = s
			col.hasNullStr = true
			col.nullStr = s
		}
		cols[idx-1] = col
	}
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

// ReadASCIIColumn reads a column of an ASCII table. For character columns
// (TFORM Aw) T must be ~string; for I, F, E and D columns T must be a
// numeric type.
//
// Fields are decoded the way cfitsio decodes them, except for TNULL
// matching, exact decimal arithmetic and -0 (kept in an unscaled column;
// cfitsio returns +0; under TSCAL/TZERO it is the value 0, so the result
// is exactly TZERO), all below:
//
//   - Leading blanks, blanks between the sign and the digits and blanks
//     embedded in the digits are ignored; E or D introduces the exponent.
//   - A field without a decimal point gets the implied decimal of its
//     TFORM: "12345" in F8.2 is 123.45, "12500E-2" in E11.3 is 0.125.
//   - An all-blank numeric field is 0 (not null), as in cfitsio.
//   - A field with a sign, a decimal point or an exponent but no digits
//     ("+", ".", "E5") is 0, as in cfitsio.
//   - TSCALn/TZEROn are applied as for binary tables, with the literals
//     read exactly. When TSCAL is 1 and the TZEROn literal is an integer
//     (of any size) the field is read as an exact decimal: an integer
//     value of any size is exact (int64, or uint64 above 2^63-1; beyond
//     those only float targets accept it), and a value with a fraction is
//     the correctly rounded float64 (and, for float32 targets, float32) of
//     decimal + TZERO and never fits an integer target, even when that
//     float64 is integral
//     ("9223372036854775807.5" is not 2^63). Otherwise float targets get
//     the field's float64 scaled in float64, and integer targets the
//     exact decimal*TSCAL + TZERO ("4611686018427387905" with TSCAL 2 is
//     9223372036854775810). cfitsio reads every field through a double,
//     so it differs for integers beyond 2^53 and is not always correctly
//     rounded for long fields.
//   - A field with a fraction is never read into an integer target, at
//     any length or exponent ("1.5E-900" is not 0). Beyond 10^400 in
//     magnitude it is +-Inf; below 10^-400 it is +-0, or, under an integer
//     TZERO, TZERO correctly rounded with the field's sign deciding a tie
//     ("1.5E-900" + (2^53+1) is 2^53+2).
//   - A field equal to TNULLn padded with blanks to the field width is
//     null (FITS 4.0 §7.2.2). cfitsio instead treats every field that
//     merely STARTS with the TNULL string as null (so with TNULL '-' it
//     reads "-1.234E-100" as null); astropy follows the standard.
//   - Character fields lose their trailing blanks.
//   - Anything else in a numeric field is an error (exponents of more than
//     6 significant digits included); an exponent within that limit that
//     overflows float64 gives +-Inf (or 0 on underflow) for float targets.
//
// A null field is NaN for float T and an ErrNullValue error otherwise
// (use ReadASCIIColumnMasked). Integer T must hold the value exactly (an
// ErrTypeMismatch error otherwise), as for binary tables.
func ReadASCIIColumn[T ColumnValue](t *ASCIITableHDU, col int) ([]T, error) {
	out, _, err := readASCII[T](t, col, false)
	return out, err
}

// ReadASCIIColumnMasked is ReadASCIIColumn plus a validity mask (false =
// the field matched TNULLn). Null rows hold the zero value of T.
func ReadASCIIColumnMasked[T ColumnValue](t *ASCIITableHDU, col int) ([]T, []bool, error) {
	return readASCII[T](t, col, true)
}

func readASCII[T ColumnValue](t *ASCIITableHDU, col int, masked bool) ([]T, []bool, error) {
	cols, err := t.Columns()
	if err != nil {
		return nil, nil, err
	}
	if col < 1 || col > len(cols) {
		return nil, nil, fmt.Errorf("fits: column %d out of range [1,%d]", col, len(cols))
	}
	c := cols[col-1]
	if len(t.rec.shape) < 2 {
		return nil, nil, fmt.Errorf("fits: ASCII table missing NAXIS2")
	}
	kind := reflect.TypeFor[T]().Kind()
	isChar := c.ascii.Type == tform.AsciiChar
	switch {
	case isChar && kind != reflect.String:
		return nil, nil, wrongReader(c, "a character column needs T=string")
	case !isChar && (kind == reflect.String || kind == reflect.Bool):
		return nil, nil, wrongReader(c, "a numeric ASCII column needs a numeric T")
	}
	nrows := t.rec.shape[1]
	rowWidth := int(t.rec.shape[0])
	buf := make([]byte, int64(rowWidth)*nrows)
	if len(buf) > 0 {
		if err := t.rec.file.br.ReadRange(t.rec.dataStart, buf); err != nil {
			return nil, nil, err
		}
	}
	sc := columnScaling(c)
	out := make([]T, nrows)
	valid := make([]bool, nrows)
	for r := int64(0); r < nrows; r++ {
		base := int(r)*rowWidth + c.byteOffset
		field := string(buf[base : base+c.byteSize])
		null := c.hasNullStr && isASCIINull(field, c.nullStr)
		if isChar {
			if null {
				if !masked {
					return nil, nil, fmt.Errorf("fits: column %d row %d: %w (field equals TNULL; use ReadASCIIColumnMasked)", col, r, ErrNullValue)
				}
				continue
			}
			out[r] = fromString[T](strings.TrimRight(field, " "))
			valid[r] = true
			continue
		}
		var v tval
		if null {
			v = tval{null: true}
		} else {
			v, err = parseASCIINumber(field, c.ascii, &sc)
			if err != nil {
				return nil, nil, fmt.Errorf("fits: column %d row %d: %w", col, r, err)
			}
		}
		x, err := numTo[T](v, masked, kind)
		if err != nil {
			return nil, nil, fmt.Errorf("fits: column %d row %d: %w", col, r, err)
		}
		out[r] = x
		valid[r] = !v.null
	}
	return out, valid, nil
}

// isASCIINull reports whether field equals tnull blank-padded to the
// field width (§7.2.2: "The string is implicitly space filled to the width
// of the field").
func isASCIINull(field, tnull string) bool {
	if len(tnull) > len(field) {
		return false
	}
	if field[:len(tnull)] != tnull {
		return false
	}
	for i := len(tnull); i < len(field); i++ {
		if field[i] != ' ' {
			return false
		}
	}
	return true
}

// parseASCIINumber decodes one numeric ASCII-table field (see
// ReadASCIIColumn) and applies the column scaling.
func parseASCIINumber(field string, af tform.ASCIIForm, sc *colScaling) (tval, error) {
	i, n := 0, len(field)
	skip := func() {
		for i < n && field[i] == ' ' {
			i++
		}
	}
	digits := func() string {
		var b strings.Builder
		for i < n && field[i] >= '0' && field[i] <= '9' {
			b.WriteByte(field[i])
			i++
			skip()
		}
		return b.String()
	}
	skip()
	if i == n {
		return sc.decodeASCIIInt(0), nil // all blank: 0, as cfitsio
	}
	neg := false
	if field[i] == '+' || field[i] == '-' {
		neg = field[i] == '-'
		i++
		skip()
	}
	intDigits := digits()
	var fracDigits string
	hasDot := false
	if i < n && field[i] == '.' {
		hasDot = true
		i++
		skip()
		fracDigits = digits()
	}
	exp := 0
	if i < n && (field[i] == 'E' || field[i] == 'D') {
		i++
		skip()
		eneg := false
		if i < n && (field[i] == '+' || field[i] == '-') {
			eneg = field[i] == '-'
			i++
			skip()
		}
		// Leading zeros are not significant ("E+0000002" is 10^2); more
		// than 6 significant exponent digits is out of every float range
		// and would overflow the integer parse.
		ed := strings.TrimLeft(digits(), "0")
		if len(ed) > 6 {
			return tval{}, fmt.Errorf("fits: ASCII field %q: exponent out of range", field)
		}
		if ed != "" {
			e, _ := strconv.Atoi(ed)
			exp = e
		}
		if eneg {
			exp = -exp
		}
	}
	if i != n {
		return tval{}, fmt.Errorf("fits: cannot read a number from ASCII field %q (TFORM %s)", field, af.Raw)
	}
	mant := intDigits + fracDigits
	exp10 := exp - len(fracDigits)
	if !hasDot {
		exp10 -= af.Fraction // implied decimal point
	}
	mant = strings.TrimLeft(mant, "0")
	if mant == "" {
		if neg && sc.identity {
			return tval{kind: tvFloat, f: math.Copysign(0, -1)}, nil
		}
		// Under scaling -0 is the value 0: 0*TSCAL + TZERO exactly (the
		// float64 arithmetic gives the same double for -0 and +0).
		return sc.decodeASCIIInt(0), nil
	}
	// An integer value (no fraction, or a fraction of zeros): exact, of
	// any size (bounded by the 6-digit exponent check above and the field
	// width; beyond 400 digits it is outside every float range anyway).
	var intDigits10 string
	isInt := false
	switch {
	case exp10 >= 0 && len(mant)+exp10 <= 400:
		intDigits10, isInt = mant+strings.Repeat("0", exp10), true
	case exp10 < 0 && len(mant) > -exp10 && strings.TrimRight(mant[len(mant)+exp10:], "0") == "":
		intDigits10, isInt = mant[:len(mant)+exp10], true
	}
	if isInt {
		if len(intDigits10) <= 18 {
			iv, _ := strconv.ParseInt(intDigits10, 10, 64)
			if neg {
				iv = -iv
			}
			return sc.decodeASCIIInt(iv), nil
		}
		b, _ := new(big.Int).SetString(intDigits10, 10)
		if neg {
			b.Neg(b)
		}
		return sc.decodeASCIIBig(b), nil
	}
	// A value with a fraction. mag bounds its magnitude: 10^(mag-1) <=
	// |decimal| < 10^mag. Within 10^-400 .. 10^400 it is evaluated exactly
	// (the integers involved are bounded by the field width and 400
	// digits); outside, the 6-digit exponent could make them huge, and the
	// result does not need them: beyond 10^400 it is +-Inf after any TZERO
	// a card can hold, and below 10^-400 it changes no rounding of an
	// integer TZERO except the direction of a tie, which its sign decides.
	mag := len(mant) + exp10
	exact := sc.identity || sc.intZero || sc.bigZero != nil
	var r *big.Rat
	switch {
	case mag > 400:
	case mag < -400:
		if exact {
			// A stand-in of the same sign and below every rounding step
			// of decimal + TZERO (|TZERO| < 10^71 is an integer, so the
			// doubles near it are at least 2^-54 apart and only TZERO
			// itself can be a tie).
			r = new(big.Rat).SetFrac(big.NewInt(1), new(big.Int).Exp(big.NewInt(10), big.NewInt(401), nil))
			if neg {
				r.Neg(r)
			}
		}
	default:
		r, _ = new(big.Rat).SetString(mant)
		if exp10 > 0 {
			r.Mul(r, new(big.Rat).SetInt(new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(exp10)), nil)))
		} else {
			r.Quo(r, new(big.Rat).SetInt(new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(-exp10)), nil)))
		}
		if neg {
			r.Neg(r)
		}
	}
	if exact {
		// A decimal with a fraction under exact arithmetic: the correctly
		// rounded value of decimal + TZERO, and never an integer target
		// (the double may be integral: 9223372036854775807.5 -> 2^63).
		if r == nil { // beyond 10^400
			f := math.Inf(1)
			if neg {
				f = math.Inf(-1)
			}
			return tval{kind: tvFloat, f: f, noInt: true}, nil
		}
		return ratTval(r.Add(r, sc.zeroR)), nil
	}
	if r != nil {
		// Float64 path: float targets get the field's float64 scaled in
		// float64, integer targets the exact decimal*TSCAL + TZERO.
		f, _ := r.Float64()
		v := sc.decodeASCIIFloat(f)
		v.src, v.sr, v.sc = srcRat, r, sc
		return v, nil
	}
	// Beyond 10^+-400 on the float64 path: the field's float64 (+-Inf or
	// +-0) scaled in float64, never an integer target.
	f := math.Inf(1)
	if mag < -400 {
		f = 0
	}
	if neg {
		f = -f
	}
	v := sc.decodeASCIIFloat(f)
	v.noInt = true
	return v, nil
}

// decodeASCIIInt scales an exact integer field value.
func (s *colScaling) decodeASCIIInt(v int64) tval {
	saved := s.hasNull
	s.hasNull = false // ASCII nulls are string matches, handled earlier
	t := s.decodeInt(v, false)
	s.hasNull = saved
	return t
}

// decodeASCIIBig scales an exact integer field value outside int64.
func (s *colScaling) decodeASCIIBig(b *big.Int) tval {
	switch {
	case s.identity:
		return bigTval(b)
	case s.intZero:
		return bigTval(b.Add(b, big.NewInt(s.zeroI)))
	case s.bigZero != nil:
		return bigTval(b.Add(b, s.bigZero))
	}
	f, _ := new(big.Float).SetInt(b).Float64()
	v := s.decodeASCIIFloat(f)
	v.src, v.sr, v.sc = srcRat, new(big.Rat).SetInt(b), s
	return v
}

// decodeASCIIFloat scales a float field value (NaN is not null here: an
// ASCII field cannot spell NaN).
func (s *colScaling) decodeASCIIFloat(v float64) tval {
	if !s.identity {
		v = float64(v*s.scale) + s.zero
	}
	return tval{kind: tvFloat, f: v}
}
