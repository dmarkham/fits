package fits

import (
	"fmt"
	"math"
	"reflect"
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

// Complex is the constraint for complex (C, M) column reads.
type Complex interface{ ~complex64 | ~complex128 }

// How binary-table cells are decoded (all readers below):
//
//   - B, I, J, K, E, D: TSCALn/TZEROn are applied as
//     phys = stored*TSCAL + TZERO, where the TSCALn/TZEROn values are the
//     literals on the cards read as exact decimals (TZERO =
//     9223372036854775807 or 4611686018427387905.0 is that integer, not
//     its nearest double; TSCAL = 0.1 is one tenth).
//   - When TSCAL is 1 and TZERO is an integer (of any size) the
//     arithmetic is exact: B/I/J/K get exact integer arithmetic, so K
//     (int64) columns keep every bit and the unsigned conventions of
//     §7.3.2 Table 19 (B with TZERO -128 = int8, I with 32768 = uint16,
//     J with 2147483648 = uint32, K with 9223372036854775808 = uint64)
//     are exact. A result
//     above 2^63-1 is exact as a uint64; one outside both ranges only
//     float targets accept (the correctly rounded value). cfitsio holds
//     TZERO as a double instead, so it is not exact there: it reads
//     TZERO = 9223372036854775809 or 9223372036854775807 as the uint64
//     convention (2^63); we add the literal. E/D values get stored + TZERO
//     correctly rounded: 1.0 + 9007199254740993 is 2^53+2, where adding
//     TZERO's double (cfitsio, astropy) gives 2^53.
//   - Otherwise (TSCAL != 1 or a fractional TZERO) float targets get the
//     float64 arithmetic float64(stored)*TSCAL + TZERO, like cfitsio, and
//     integer targets get the exact value stored*TSCAL + TZERO: TSCAL 2
//     on stored 2^61+1 is 2^62+2 (float64 says 2^62); TZERO 0.5 on 2^60 is
//     not an integer (float64 says 2^60).
//   - A TZERO float literal that is not exactly 2^63 but whose nearest
//     double is (9.223372036854776E+18, which cfitsio and astropy take as
//     the uint64 convention) is ambiguous: integer targets are refused
//     (ErrTypeMismatch), float targets get the float64 arithmetic.
//   - A stored integer equal to TNULLn is null (the comparison is on the
//     stored value, before scaling). E and D: NaN is null.
//   - Integer Go targets must hold the exact physical value: a fraction,
//     a value outside the target's range or a non-finite value is an
//     ErrTypeMismatch error for the whole read, never a rounded, truncated
//     or wrapped number. Float targets get the physical value: on the
//     exact path (TSCAL 1, integer TZERO, and unscaled columns) the
//     float64 or float32 nearest the exact value (a float32 target does
//     not round twice through the float64: 1K 2^60+2^36+1 is 0x5d800001,
//     as cfitsio's TFLOAT read gives), otherwise the float64 arithmetic
//     result, or its float32 conversion (what cfitsio computes).
//   - Non-masked readers return NaN for a null cell when T is a float and
//     ErrNullValue otherwise. The *Masked readers return the zero value
//     and valid[i] = false.
//   - L: 'T' is true, 'F' false, 0x00 undefined (null); any other byte is
//     an error.
//   - A: the cell is cut at its first NUL and trailing blanks are removed.
//     A cell holding several strings (TDIM (w, d2, ...) string arrays, rAw
//     substring arrays) is read with ReadStringArrayColumn; ReadColumn
//     refuses it rather than return the cell as one string.
//   - TSCAL/TZERO on C/M columns are not supported (an error); RowBytes
//     gives the raw row bytes.

// ReadColumn reads a scalar column (Repeat==1) across all rows and returns
// one value per row as a []T. For Char columns (TFORM "rA") use T=string:
// the whole r-character cell is one string (a cell that holds an array of
// strings is an error; use ReadStringArrayColumn). For Logical columns use
// T=bool. For numeric columns (B, I, J, K, E, D) use any Numeric T that
// holds the physical values exactly (see the rules above).
//
// Bit (X) columns are read with ReadBitColumn, complex (C, M) columns with
// ReadComplexColumn, vector columns with ReadVectorColumn and VLAs with
// ReadVarColumn and friends; ReadColumn returns an error for them.
func ReadColumn[T ColumnValue](t *BinaryTableHDU, col int) ([]T, error) {
	out, _, err := readScalarColumn[T](t, col, false)
	return out, err
}

// ReadColumnMasked is ReadColumn plus a validity mask: valid[i] is false
// when row i is null (TNULLn match, NaN float, undefined logical). Null
// rows hold the zero value of T. String columns are always valid.
func ReadColumnMasked[T ColumnValue](t *BinaryTableHDU, col int) ([]T, []bool, error) {
	return readScalarColumn[T](t, col, true)
}

func (t *BinaryTableHDU) column(col int) (Column, error) {
	cols, err := t.Columns()
	if err != nil {
		return Column{}, err
	}
	if col < 1 || col > len(cols) {
		return Column{}, fmt.Errorf("fits: column %d out of range [1,%d]", col, len(cols))
	}
	return cols[col-1], nil
}

// readRows reads the fixed-width row area (NAXIS1*NAXIS2 bytes).
func (t *BinaryTableHDU) readRows() ([]byte, error) {
	total := int64(t.rowBytes) * t.NumRows()
	buf := make([]byte, total)
	if total > 0 {
		if err := t.rec.file.br.ReadRange(t.rec.dataStart, buf); err != nil {
			return nil, err
		}
	}
	return buf, nil
}

// cells calls fn with the bytes of column c in every row.
func (t *BinaryTableHDU) cells(c Column, fn func(r int64, cell []byte) error) error {
	buf, err := t.readRows()
	if err != nil {
		return err
	}
	nrows := t.NumRows()
	for r := int64(0); r < nrows; r++ {
		base := int(r)*t.rowBytes + c.byteOffset
		if err := fn(r, buf[base:base+c.byteSize]); err != nil {
			return err
		}
	}
	return nil
}

func wrongReader(c Column, what string) error {
	return fmt.Errorf("fits: column %d (%q, TFORM %s): %s", c.Index, c.Name, c.TForm, what)
}

// readerFor names the reader to use for a column, for error messages.
func readerFor(c Column) string {
	switch c.bin.Type {
	case tform.BinPArrayDesc32, tform.BinQArrayDesc64:
		switch c.bin.VarType {
		case tform.BinLogical:
			return "use ReadVarLogicalColumn"
		case tform.BinBit:
			return "use ReadVarBitColumn"
		case tform.BinChar:
			return "use ReadVarStringColumn"
		case tform.BinComplex64, tform.BinComplex128:
			return "use ReadVarComplexColumn"
		}
		return "use ReadVarColumn"
	case tform.BinBit:
		return "use ReadBitColumn"
	case tform.BinComplex64, tform.BinComplex128:
		return "use ReadComplexColumn"
	case tform.BinLogical:
		return "use ReadLogicalColumn"
	case tform.BinChar:
		return "use ReadColumn[string]"
	}
	return "use ReadVectorColumn"
}

func readScalarColumn[T ColumnValue](t *BinaryTableHDU, col int, masked bool) ([]T, []bool, error) {
	c, err := t.column(col)
	if err != nil {
		return nil, nil, err
	}
	kind := reflect.TypeFor[T]().Kind()
	nrows := t.NumRows()
	out := make([]T, nrows)
	valid := make([]bool, nrows)
	switch {
	case c.bin.Type == tform.BinChar:
		if kind != reflect.String {
			return nil, nil, wrongReader(c, "a character column needs T=string")
		}
		_, n, lerr := stringLayout(c)
		if lerr != nil {
			return nil, nil, lerr
		}
		if n != 1 {
			return nil, nil, wrongReader(c, fmt.Sprintf("the cell is an array of %d strings (rAw substring array or TDIM string array); use ReadStringArrayColumn", n))
		}
		err = t.cells(c, func(r int64, cell []byte) error {
			out[r] = fromString[T](cellString(cell))
			valid[r] = true
			return nil
		})
	case c.bin.Type == tform.BinLogical && c.Repeat == 1:
		if kind != reflect.Bool {
			return nil, nil, wrongReader(c, "a logical column needs T=bool")
		}
		err = t.cells(c, func(r int64, cell []byte) error {
			v, ok, err := decodeLogical(cell[0])
			if err != nil {
				return fmt.Errorf("fits: column %d row %d: %w", col, r, err)
			}
			if !ok && !masked {
				return fmt.Errorf("fits: column %d row %d: %w (undefined logical; use ReadColumnMasked)", col, r, ErrNullValue)
			}
			out[r] = fromBool[T](v)
			valid[r] = ok
			return nil
		})
	case isNumericBinType(c.bin.Type) && c.Repeat == 1:
		if kind == reflect.String || kind == reflect.Bool {
			return nil, nil, wrongReader(c, "a numeric column needs a numeric T")
		}
		sc := columnScaling(c)
		fd := newFastDec(&sc, c.bin.Type, kind)
		err = t.cells(c, func(r int64, cell []byte) error {
			if fd.ok && cellFast(&fd, c.bin.Type, cell, &out[r]) {
				valid[r] = true
				return nil
			}
			v := sc.decodeElem(c.bin.Type, cell)
			x, err := numTo[T](v, masked, kind)
			if err != nil {
				return fmt.Errorf("fits: column %d row %d: %w", col, r, err)
			}
			out[r] = x
			valid[r] = !v.null
			return nil
		})
	case c.Repeat != 1 && (isNumericBinType(c.bin.Type) || c.bin.Type == tform.BinLogical):
		return nil, nil, wrongReader(c, fmt.Sprintf("repeat %d; use ReadVectorColumn / ReadLogicalColumn", c.Repeat))
	default:
		return nil, nil, wrongReader(c, readerFor(c))
	}
	if err != nil {
		return nil, nil, err
	}
	return out, valid, nil
}

// stringLayout returns the width of one string and the number of strings
// in a cell of character column c (FITS 4.0 7.3.2 and 7.3.3.1): a TDIMn
// (w, d2, d3, ...) makes the cell an array of d2*d3*... strings of w
// characters, and it wins over an rAw width (8A2 with TDIM (8) is one
// string of 8, with TDIM (4,2) two strings of 4: cfitsio 4.7.0 and astropy
// 8.0.2 both read it so, cf_bintable_strarray.fits W8/W42); otherwise TFORM
// rAw makes it r/w strings of w characters; otherwise the cell is one
// string. A TDIM that covers fewer characters than the repeat count and an
// rAw whose r is not a multiple of w are errors: cfitsio splits such cells
// into repeat/width strings while the standard leaves the rest undefined,
// and neither can be proven to be the writer's meaning.
func stringLayout(c Column) (width, n int64, err error) {
	repeat := c.bin.Repeat
	switch {
	case len(c.Dim) > 0:
		w := c.Dim[0]
		k := int64(1)
		for _, d := range c.Dim[1:] {
			k *= d
		}
		if w <= 0 || w*k != repeat {
			return 0, 0, wrongReader(c, fmt.Sprintf("TDIM %v covers %d of the %d characters of the cell; string arrays whose TDIM does not cover the whole cell are not supported", c.Dim, w*k, repeat))
		}
		return w, k, nil
	case c.bin.Width > 0 && c.bin.Width < repeat:
		if repeat%c.bin.Width != 0 {
			return 0, 0, wrongReader(c, fmt.Sprintf("rAw with r = %d not a multiple of w = %d is not supported", repeat, c.bin.Width))
		}
		return c.bin.Width, repeat / c.bin.Width, nil
	}
	return repeat, 1, nil
}

// ReadStringArrayColumn reads a character column (TFORM rA) whose cells
// are arrays of strings, as one []string per row: a TDIMn string array
// (TDIM (w, d2, ...), the first axis being the string length; astropy
// writes 2-D string columns this way) or an rAw substring array (r/w
// strings of w characters). Each string is cut at its first NUL and loses
// its trailing blanks (§7.3.3.1), as ReadColumn[string] does for a whole
// cell. A plain rA column gives one string per row. The layouts that
// stringLayout refuses are an error.
func ReadStringArrayColumn(t *BinaryTableHDU, col int) ([][]string, error) {
	c, err := t.column(col)
	if err != nil {
		return nil, err
	}
	if c.bin.Type != tform.BinChar {
		return nil, wrongReader(c, "not a character column; "+readerFor(c))
	}
	w, n, err := stringLayout(c)
	if err != nil {
		return nil, err
	}
	out := make([][]string, t.NumRows())
	err = t.cells(c, func(r int64, cell []byte) error {
		row := make([]string, n)
		for k := range row {
			row[k] = cellString(cell[int64(k)*w : int64(k+1)*w])
		}
		out[r] = row
		return nil
	})
	if err != nil {
		return nil, err
	}
	return out, nil
}

// cellString decodes a binary-table character cell: it ends at the first
// NUL, and trailing blanks are not significant (§7.3.3.1).
func cellString(cell []byte) string {
	end := len(cell)
	for i, b := range cell {
		if b == 0 {
			end = i
			break
		}
	}
	for end > 0 && cell[end-1] == ' ' {
		end--
	}
	return string(cell[:end])
}

// decodeLogical maps a logical byte: 'T', 'F', or 0 (undefined).
func decodeLogical(b byte) (v, ok bool, err error) {
	switch b {
	case 'T':
		return true, true, nil
	case 'F':
		return false, true, nil
	case 0:
		return false, false, nil
	}
	return false, false, fmt.Errorf("fits: invalid logical byte 0x%02x (want 'T', 'F' or 0)", b)
}

func fromString[T any](s string) T {
	var out T
	*(*string)(unsafe.Pointer(&out)) = s
	return out
}

func fromBool[T any](b bool) T {
	var out T
	*(*bool)(unsafe.Pointer(&out)) = b
	return out
}

// numTo converts a decoded element to T, whose reflect kind is kind (a
// numeric kind), following the rules above ReadColumn.
func numTo[T any](v tval, masked bool, kind reflect.Kind) (T, error) {
	var out T
	p := unsafe.Pointer(&out)
	isFloat := kind == reflect.Float32 || kind == reflect.Float64
	if v.null {
		if isFloat && !masked {
			if kind == reflect.Float64 {
				*(*float64)(p) = math.NaN()
			} else {
				*(*float32)(p) = float32(math.NaN())
			}
			return out, nil
		}
		if masked {
			return out, nil
		}
		return out, fmt.Errorf("%w (use the Masked reader, or a float type to get NaN)", ErrNullValue)
	}
	ok := true
	switch kind {
	case reflect.Float64:
		*(*float64)(p) = v.asFloat()
	case reflect.Float32:
		*(*float32)(p) = v.asFloat32()
	case reflect.Int8:
		var x int64
		x, ok = v.asInt(math.MinInt8, math.MaxInt8)
		*(*int8)(p) = int8(x)
	case reflect.Int16:
		var x int64
		x, ok = v.asInt(math.MinInt16, math.MaxInt16)
		*(*int16)(p) = int16(x)
	case reflect.Int32:
		var x int64
		x, ok = v.asInt(math.MinInt32, math.MaxInt32)
		*(*int32)(p) = int32(x)
	case reflect.Int64:
		var x int64
		x, ok = v.asInt(math.MinInt64, math.MaxInt64)
		*(*int64)(p) = x
	case reflect.Uint8:
		var x uint64
		x, ok = v.asUint(math.MaxUint8)
		*(*uint8)(p) = uint8(x)
	case reflect.Uint16:
		var x uint64
		x, ok = v.asUint(math.MaxUint16)
		*(*uint16)(p) = uint16(x)
	case reflect.Uint32:
		var x uint64
		x, ok = v.asUint(math.MaxUint32)
		*(*uint32)(p) = uint32(x)
	case reflect.Uint64:
		var x uint64
		x, ok = v.asUint(math.MaxUint64)
		*(*uint64)(p) = x
	default:
		return out, fmt.Errorf("fits: unsupported target kind %v", kind)
	}
	if !ok {
		var zero T
		if v.sc != nil && v.sc.ambiguous {
			return zero, fmt.Errorf("%w: TZERO is a float literal near 2^63 but not 2^63 (the uint64 convention needs exactly 9223372036854775808); integer targets refused, read as float", ErrTypeMismatch)
		}
		return zero, fmt.Errorf("%w: table value %s does not fit %v exactly", ErrTypeMismatch, v.String(), reflect.TypeFor[T]())
	}
	return out, nil
}

// ReadVectorColumn reads a numeric (B, I, J, K, E, D) column of any repeat
// count and returns a flat slice of length nrows*repeat in row-major order
// (element k of row r at index r*repeat+k). Decoding follows the rules
// described above ReadColumn; a null element is an ErrNullValue error for
// integer T (use ReadVectorColumnMasked) and NaN for float T.
func ReadVectorColumn[T Numeric](t *BinaryTableHDU, col int) ([]T, error) {
	out, _, err := readVector[T](t, col, false)
	return out, err
}

// ReadVectorColumnMasked is ReadVectorColumn plus a validity mask of the
// same length (false = null element).
func ReadVectorColumnMasked[T Numeric](t *BinaryTableHDU, col int) ([]T, []bool, error) {
	return readVector[T](t, col, true)
}

func readVector[T Numeric](t *BinaryTableHDU, col int, masked bool) ([]T, []bool, error) {
	c, err := t.column(col)
	if err != nil {
		return nil, nil, err
	}
	if !isNumericBinType(c.bin.Type) {
		return nil, nil, wrongReader(c, "not a numeric column; "+readerFor(c))
	}
	kind := reflect.TypeFor[T]().Kind()
	repeat := c.Repeat
	n := t.NumRows() * repeat
	out := make([]T, n)
	valid := make([]bool, n)
	sc := columnScaling(c)
	fd := newFastDec(&sc, c.bin.Type, kind)
	es := c.bin.Type.ElementSize()
	err = t.cells(c, func(r int64, cell []byte) error {
		for k := int64(0); k < repeat; k++ {
			if fd.ok && cellFast(&fd, c.bin.Type, cell[int(k)*es:], &out[r*repeat+k]) {
				valid[r*repeat+k] = true
				continue
			}
			v := sc.decodeElem(c.bin.Type, cell[int(k)*es:])
			x, err := numTo[T](v, masked, kind)
			if err != nil {
				return fmt.Errorf("fits: column %d row %d element %d: %w", col, r, k, err)
			}
			out[r*repeat+k] = x
			valid[r*repeat+k] = !v.null
		}
		return nil
	})
	if err != nil {
		return nil, nil, err
	}
	return out, valid, nil
}

// ReadLogicalColumn reads a logical (L) column of any repeat count as a
// flat nrows*repeat slice plus a validity mask (false = undefined, byte
// 0x00). A byte other than 'T', 'F' or 0 is an error.
func ReadLogicalColumn(t *BinaryTableHDU, col int) (values, valid []bool, err error) {
	c, err := t.column(col)
	if err != nil {
		return nil, nil, err
	}
	if c.bin.Type != tform.BinLogical {
		return nil, nil, wrongReader(c, "not a logical column; "+readerFor(c))
	}
	repeat := c.Repeat
	n := t.NumRows() * repeat
	values = make([]bool, n)
	valid = make([]bool, n)
	err = t.cells(c, func(r int64, cell []byte) error {
		for k := int64(0); k < repeat; k++ {
			v, ok, err := decodeLogical(cell[k])
			if err != nil {
				return fmt.Errorf("fits: column %d row %d element %d: %w", col, r, k, err)
			}
			values[r*repeat+k], valid[r*repeat+k] = v, ok
		}
		return nil
	})
	if err != nil {
		return nil, nil, err
	}
	return values, valid, nil
}

// ReadBitColumn reads a bit (X) column as a flat nrows*repeat slice of
// bits, row-major; within a row bit 0 is the most significant bit of the
// first byte (§7.3.3.2).
func ReadBitColumn(t *BinaryTableHDU, col int) ([]bool, error) {
	c, err := t.column(col)
	if err != nil {
		return nil, err
	}
	if c.bin.Type != tform.BinBit {
		return nil, wrongReader(c, "not a bit column; "+readerFor(c))
	}
	repeat := c.Repeat
	out := make([]bool, t.NumRows()*repeat)
	err = t.cells(c, func(r int64, cell []byte) error {
		unpackBits(out[r*repeat:(r+1)*repeat], cell)
		return nil
	})
	if err != nil {
		return nil, err
	}
	return out, nil
}

func unpackBits(dst []bool, src []byte) {
	for k := range dst {
		dst[k] = src[k/8]&(0x80>>(k%8)) != 0
	}
}

// ReadComplexColumn reads a complex (C or M) column of any repeat count as
// a flat nrows*repeat slice. Reading M into complex64 rounds each part to
// float32. TSCALn/TZEROn on a complex column are not supported (error).
func ReadComplexColumn[T Complex](t *BinaryTableHDU, col int) ([]T, error) {
	c, err := t.column(col)
	if err != nil {
		return nil, err
	}
	if c.bin.Type != tform.BinComplex64 && c.bin.Type != tform.BinComplex128 {
		return nil, wrongReader(c, "not a complex column; "+readerFor(c))
	}
	if !columnScaling(c).identity { // any TSCAL/TZERO literal other than exactly 1 / 0
		return nil, wrongReader(c, "TSCAL/TZERO on complex columns is not supported")
	}
	repeat := c.Repeat
	out := make([]T, t.NumRows()*repeat)
	err = t.cells(c, func(r int64, cell []byte) error {
		decodeComplex(out[r*repeat:(r+1)*repeat], cell, c.bin.Type)
		return nil
	})
	if err != nil {
		return nil, err
	}
	return out, nil
}

func decodeComplex[T Complex](dst []T, src []byte, typ tform.BinaryType) {
	for k := range dst {
		if typ == tform.BinComplex64 {
			re := bigendian.Float32(src[k*8:])
			im := bigendian.Float32(src[k*8+4:])
			dst[k] = T(complex(float64(re), float64(im)))
		} else {
			re := bigendian.Float64(src[k*16:])
			im := bigendian.Float64(src[k*16+8:])
			dst[k] = T(complex(re, im))
		}
	}
}
