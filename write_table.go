package fits

import (
	"fmt"
	"io"
	"math"
	"math/big"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/tform"
)

// ColumnData describes a single column to write into a new binary table.
// Exactly one Data* field must be set; its Go type selects the TFORM type
// code:
//
//	DataUint8 B   DataInt16 I   DataInt32 J   DataInt64 K
//	DataFloat32 E DataFloat64 D DataString A  DataBool L  DataBits X
//	DataComplex64 C               DataComplex128 M
//	DataInt8 B + TZERO -128       DataUint16 I + TZERO 32768
//	DataUint32 J + TZERO 2^31     DataUint64 K + TZERO 2^63
//
// The last four hold PHYSICAL values: the writer stores value-TZERO and
// writes the TZEROn card (the unsigned-integer convention of FITS 4.0
// §7.3.2, Table 19), so readers that apply TZERO (cfitsio, astropy, this
// package) get the values back exactly.
//
// Fixed-width columns take a flat slice of nrows*repeat elements, row
// after row. The repeat count comes from TForm (e.g. "3E"), else from the
// product of Dim, else 1. DataBits is one bool per bit. DataString holds
// one string per row of printable ASCII; the column width is TForm's
// repeat (a longer string is an error) or else the longest string.
//
// Variable-length array (VLA, §7.3.5) columns take one inner slice per
// row (DataVar*), or one string per row (DataVarString). The writer emits
// a 1P<type>(max) descriptor column, or 1Q when TForm says so (e.g.
// "1QE"), and appends the payloads to the heap right after the rows.
//
// TScal/TZero/TNull write TSCALn/TZEROn/TNULLn. With them, the Data*
// values of a column are the STORED values (no inverse scaling is
// applied). TScal/TZero are only allowed on B, I, J, K, E, D columns (and
// their VLAs) and not together with the four TZERO types above. An
// integral value is written as its integer literal, whose exact decimal
// is the float64 (TZERO 2^63 as 9223372036854775808); one of more than 70
// characters (|v| >= 1e70 or so) does not fit a card and is an error. Any
// other value is written as the shortest decimal that reads back as the
// same float64 (0.1 as "0.1", 1/3 as "0.3333333333333333"); a reader that
// takes the literal exactly, as this library does, scales by that decimal
// (stored 10 with TScal 0.1 is exactly 1), not by the float64's binary
// value. TNull only
// on integer columns. Valid (fixed columns, one flag per element) and
// VarValid (VLA columns, one flag per element of each row) mark null
// elements: an integer element is written as TNULLn (which must then be
// set), a float as NaN, a logical as 0x00; A, X, C and M have no null and
// a false flag there is an error. A non-null integer equal to TNULLn is an
// error too (it would read back as null).
type ColumnData struct {
	Name    string
	Unit    string
	Display string
	TForm   string // optional; must agree with the Data* type
	Dim     []int64

	TScal *float64
	TZero *float64
	TNull *int64

	DataUint8      []uint8
	DataInt8       []int8
	DataInt16      []int16
	DataUint16     []uint16
	DataInt32      []int32
	DataUint32     []uint32
	DataInt64      []int64
	DataUint64     []uint64
	DataFloat32    []float32
	DataFloat64    []float64
	DataComplex64  []complex64
	DataComplex128 []complex128
	DataString     []string
	DataBool       []bool
	DataBits       []bool

	// Valid marks null elements of a fixed-width column (nil = all valid).
	Valid []bool

	DataVarUint8      [][]uint8
	DataVarInt16      [][]int16
	DataVarInt32      [][]int32
	DataVarInt64      [][]int64
	DataVarFloat32    [][]float32
	DataVarFloat64    [][]float64
	DataVarComplex64  [][]complex64
	DataVarComplex128 [][]complex128
	DataVarBool       [][]bool
	DataVarBits       [][]bool
	DataVarString     []string

	// VarValid marks null elements of a VLA column (nil = all valid).
	VarValid [][]bool
}

// AppendBinaryTable appends a new BINTABLE HDU to the file. All columns must
// have the same row count. The TFORM code is inferred from the DataXxx slice
// that is populated (and checked against an explicit .TForm if set).
//
// hdr may carry extra cards for the new HDU. Of the column keywords,
// TTYPEn, TUNITn, TDISPn, TDIMn, TSCALn, TZEROn, TNULLn and TFORMn are
// honoured as if they were set in cols[n-1] (a conflicting value is an
// error, and so is an n beyond the number of columns). The structural
// keywords XTENSION, BITPIX, NAXIS, NAXISn, PCOUNT, GCOUNT, TFIELDS,
// TBCOLn, THEAP, SIMPLE, EXTEND and END are always computed by the writer;
// copies in hdr are not written.
//
// Caller cards are copied with header.AddCard. A card read from a file is
// written as its bytes. An edited Cards() copy (Key/Type/Value/Comment
// changed) is written from its fields under the upper-cased keyword, and
// is left out like any caller card for a keyword the writer computes
// itself. An unedited lowercase-keyword (non-HIERARCH) card makes the
// writer return an error wrapping header.ErrInvalidKeyword: cfitsio and
// astropy read such a card differently.
func AppendBinaryTable(f *File, hdr *header.Header, cols []ColumnData) (*BinaryTableHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	if len(cols) == 0 {
		return nil, fmt.Errorf("fits: AppendBinaryTable: no columns")
	}
	if len(cols) > 999 {
		return nil, fmt.Errorf("fits: AppendBinaryTable: %d columns, FITS allows 999", len(cols))
	}
	cols = append([]ColumnData(nil), cols...) // hdr keys are merged into a copy
	extra, err := mergeColumnKeys(hdr, cols)
	if err != nil {
		return nil, err
	}
	specs := make([]*colSpec, len(cols))
	nrows := int64(-1)
	for i := range cols {
		s, err := describeColumn(&cols[i])
		if err != nil {
			return nil, fmt.Errorf("fits: AppendBinaryTable: column %d (%q): %w", i+1, cols[i].Name, err)
		}
		if s.nrows >= 0 {
			if nrows >= 0 && s.nrows != nrows {
				return nil, fmt.Errorf("fits: AppendBinaryTable: column %d (%q): %d rows, earlier columns have %d", i+1, cols[i].Name, s.nrows, nrows)
			}
			nrows = s.nrows
		}
		specs[i] = s
	}
	if nrows < 0 {
		nrows = 0
	}
	var rowBytes int64
	for _, s := range specs {
		s.offset = rowBytes
		rowBytes += s.cellBytes
	}
	data := make([]byte, rowBytes*nrows)
	var heap []byte
	for i, s := range specs {
		for r := int64(0); r < nrows; r++ {
			cell := data[r*rowBytes+s.offset : r*rowBytes+s.offset+s.cellBytes]
			var err error
			if s.vla {
				heap, err = s.encodeVar(cell, r, heap)
			} else {
				err = s.encodeFixed(cell, r)
			}
			if err != nil {
				return nil, fmt.Errorf("fits: AppendBinaryTable: column %d (%q) row %d: %w", i+1, cols[i].Name, r, err)
			}
		}
	}

	// Every card goes through a headerBuilder: a Set/AddCard error (a user
	// card the encoder cannot represent) fails the call instead of leaving
	// the card out, and a card read from a file keeps its original bytes.
	hb := newHeaderBuilder()
	hb.set(header.KeyXtension, "BINTABLE", "binary table extension")
	hb.set(header.KeyBitpix, int64(8), "")
	hb.set(header.KeyNaxis, int64(2), "")
	hb.set("NAXIS1", rowBytes, "width of table row in bytes")
	hb.set("NAXIS2", nrows, "number of rows")
	hb.set(header.KeyPcount, int64(len(heap)), "size of heap area in bytes")
	hb.set(header.KeyGcount, int64(1), "")
	hb.set(header.KeyTfields, int64(len(cols)), "number of columns")
	for i, s := range specs {
		s.writeKeys(hb, i+1)
	}
	for _, c := range extra {
		hb.add(c)
	}
	if hb.err != nil {
		return nil, fmt.Errorf("fits: AppendBinaryTable: %w", hb.err)
	}
	h := hb.h
	rec, err := appendTable(f, h, kindBinTable, "BINTABLE", rowBytes, nrows, int64(len(cols)), append(data, heap...))
	if err != nil {
		return nil, err
	}
	return &BinaryTableHDU{rec: rec}, nil
}

// appendTable writes a table HDU (header h, data unit body) at the end of
// the file and registers it.
func appendTable(f *File, h *header.Header, kind hduKind, xt string, rowBytes, nrows, tfields int64, body []byte) (*hduRecord, error) {
	if f.mode == ModeCreate && len(f.hdus) == 0 {
		// A FITS file starts with a primary HDU; an extension cannot be
		// the first HDU. Write an empty one (BITPIX 16, NAXIS 0), as
		// cfitsio's ffcrtb does when a table is created in an empty file
		// (putkey.c: ffcrim(fptr, 16, 0, ...)).
		if _, err := f.writeImage(pickBitpix[int16](), nil, nil, anyData{bp: pickBitpix[int16](), emit: func(io.Writer) error { return nil }}); err != nil {
			return nil, err
		}
	}
	if f.mode == ModeEdit && len(f.hdus) > 0 {
		last := f.hdus[len(f.hdus)-1]
		if _, err := f.bw.Seek(last.paddedEnd, io.SeekStart); err != nil {
			return nil, err
		}
	} else if f.mode == ModeCreate {
		if _, err := f.bw.Seek(0, io.SeekEnd); err != nil {
			return nil, err
		}
	}
	// In ModeEdit the append runs under the undo journal (journal.go).
	return f.journaledAppend(f.bw.Pos(), func() (*hduRecord, error) {
		return emitTableHDU(f, h, kind, xt, rowBytes, nrows, tfields, body)
	})
}

// emitTableHDU writes the table HDU at the current writer position and
// registers it.
func emitTableHDU(f *File, h *header.Header, kind hduKind, xt string, rowBytes, nrows, tfields int64, body []byte) (*hduRecord, error) {
	start := f.bw.Pos()
	hdrBytes, err := header.Encode(h)
	if err != nil {
		return nil, err
	}
	if err := f.bw.WriteRange(hdrBytes); err != nil {
		return nil, err
	}
	dataStart := f.bw.Pos()
	if len(body) > 0 {
		if _, err := f.bw.Write(body); err != nil {
			return nil, err
		}
	}
	dataEnd := f.bw.Pos()
	pad := byte(0)
	if kind == kindASCIITable {
		pad = ' ' // §7.2.3: ASCII table data is padded with blanks
	}
	if err := f.bw.PadToBlock(pad); err != nil {
		return nil, err
	}
	rec := &hduRecord{
		index:       len(f.hdus),
		kind:        kind,
		headerStart: start,
		dataStart:   dataStart,
		dataEnd:     dataEnd,
		paddedEnd:   f.bw.Pos(),
		xtension:    xt,
		bitpix:      8,
		naxis:       2,
		shape:       []int64{rowBytes, nrows},
		pcount:      int64(len(body)) - rowBytes*nrows,
		gcount:      1,
		tfields:     tfields,
		rawHeader:   hdrBytes,
		file:        f,
	}
	// The header is parsed back from the bytes written (not kept as h),
	// so reads through the returned HDU see the card text a reopen sees
	// (the exact TSCAL/TZERO literals among them).
	f.hdus = append(f.hdus, rec)
	return rec, nil
}

// colSpec is the resolved layout of one column to write.
type colSpec struct {
	col       *ColumnData
	code      tform.BinaryType // element type
	vla       bool
	q         bool  // 1Q descriptors
	repeat    int64 // fixed: elements per row (characters for A, bits for X)
	cellBytes int64
	offset    int64
	nrows     int64 // -1 if this column does not determine it
	maxLen    int
	tform     string
	autoZero  string // TZERO of the four unsigned-convention Data types
}

// dataKind returns the element type code, whether the data is a VLA, the
// number of elements (fixed) or rows (VLA), and the automatic TZERO.
func dataKind(c *ColumnData) (code tform.BinaryType, vla bool, n int, autoZero string, err error) {
	set := 0
	pick := func(nonNil bool, t tform.BinaryType, isVar bool, length int, z string) {
		if nonNil {
			set++
			code, vla, n, autoZero = t, isVar, length, z
		}
	}
	pick(c.DataUint8 != nil, tform.BinUint8, false, len(c.DataUint8), "")
	pick(c.DataInt8 != nil, tform.BinUint8, false, len(c.DataInt8), "-128")
	pick(c.DataInt16 != nil, tform.BinInt16, false, len(c.DataInt16), "")
	pick(c.DataUint16 != nil, tform.BinInt16, false, len(c.DataUint16), "32768")
	pick(c.DataInt32 != nil, tform.BinInt32, false, len(c.DataInt32), "")
	pick(c.DataUint32 != nil, tform.BinInt32, false, len(c.DataUint32), "2147483648")
	pick(c.DataInt64 != nil, tform.BinInt64, false, len(c.DataInt64), "")
	pick(c.DataUint64 != nil, tform.BinInt64, false, len(c.DataUint64), "9223372036854775808")
	pick(c.DataFloat32 != nil, tform.BinFloat32, false, len(c.DataFloat32), "")
	pick(c.DataFloat64 != nil, tform.BinFloat64, false, len(c.DataFloat64), "")
	pick(c.DataComplex64 != nil, tform.BinComplex64, false, len(c.DataComplex64), "")
	pick(c.DataComplex128 != nil, tform.BinComplex128, false, len(c.DataComplex128), "")
	pick(c.DataString != nil, tform.BinChar, false, len(c.DataString), "")
	pick(c.DataBool != nil, tform.BinLogical, false, len(c.DataBool), "")
	pick(c.DataBits != nil, tform.BinBit, false, len(c.DataBits), "")
	pick(c.DataVarUint8 != nil, tform.BinUint8, true, len(c.DataVarUint8), "")
	pick(c.DataVarInt16 != nil, tform.BinInt16, true, len(c.DataVarInt16), "")
	pick(c.DataVarInt32 != nil, tform.BinInt32, true, len(c.DataVarInt32), "")
	pick(c.DataVarInt64 != nil, tform.BinInt64, true, len(c.DataVarInt64), "")
	pick(c.DataVarFloat32 != nil, tform.BinFloat32, true, len(c.DataVarFloat32), "")
	pick(c.DataVarFloat64 != nil, tform.BinFloat64, true, len(c.DataVarFloat64), "")
	pick(c.DataVarComplex64 != nil, tform.BinComplex64, true, len(c.DataVarComplex64), "")
	pick(c.DataVarComplex128 != nil, tform.BinComplex128, true, len(c.DataVarComplex128), "")
	pick(c.DataVarBool != nil, tform.BinLogical, true, len(c.DataVarBool), "")
	pick(c.DataVarBits != nil, tform.BinBit, true, len(c.DataVarBits), "")
	pick(c.DataVarString != nil, tform.BinChar, true, len(c.DataVarString), "")
	switch set {
	case 0:
		err = fmt.Errorf("no data set")
	case 1:
	default:
		err = fmt.Errorf("%d Data fields set, want exactly one", set)
	}
	return
}

func describeColumn(c *ColumnData) (*colSpec, error) {
	code, vla, n, autoZero, err := dataKind(c)
	if err != nil {
		return nil, err
	}
	s := &colSpec{col: c, code: code, vla: vla, autoZero: autoZero}
	if autoZero != "" && (c.TScal != nil || c.TZero != nil) {
		return nil, fmt.Errorf("TScal/TZero cannot be combined with the signed-byte/unsigned Data types (they set TZERO = %s themselves)", autoZero)
	}
	if c.TScal != nil || c.TZero != nil {
		if !isIntType(code) && code != tform.BinFloat32 && code != tform.BinFloat64 {
			return nil, fmt.Errorf("TSCAL/TZERO are not allowed on %c columns", code)
		}
		if c.TScal != nil && (*c.TScal == 0 || math.IsNaN(*c.TScal) || math.IsInf(*c.TScal, 0)) {
			return nil, fmt.Errorf("TSCAL = %v", *c.TScal)
		}
		if c.TZero != nil && (math.IsNaN(*c.TZero) || math.IsInf(*c.TZero, 0)) {
			return nil, fmt.Errorf("TZERO = %v", *c.TZero)
		}
		for _, k := range []struct {
			key string
			v   *float64
		}{{"TSCAL", c.TScal}, {"TZERO", c.TZero}} {
			if k.v == nil {
				continue
			}
			if b, ok := numberValue(*k.v).(*big.Int); ok && len(b.String()) > maxValueChars {
				return nil, fmt.Errorf("%s = %g: no exact decimal literal of this float64 fits a header card (it is an integer of %d characters, a card holds %d), and a shorter literal such as %s is a different number to a reader that takes the literal exactly", k.key, *k.v, len(b.String()), maxValueChars, strconv.FormatFloat(*k.v, 'G', -1, 64))
			}
		}
	}
	if c.TNull != nil {
		if !isIntType(code) {
			return nil, fmt.Errorf("TNULL is only allowed on integer (B, I, J, K) columns, not %c", code)
		}
		if !storedFits(code, *c.TNull) {
			return nil, fmt.Errorf("TNULL %d does not fit the stored %c type", *c.TNull, code)
		}
	}
	var form tform.BinaryForm
	if c.TForm != "" {
		form, err = tform.ParseBinary(c.TForm)
		if err != nil {
			return nil, err
		}
		gotVar := isVarType(form.Type)
		elem := form.Type
		if gotVar {
			elem = form.VarType
		}
		if gotVar != vla || elem != code {
			what := "fixed-width"
			if vla {
				what = "variable-length"
			}
			return nil, fmt.Errorf("TForm %q does not match the %s %c data", c.TForm, what, code)
		}
		if form.Width > 0 && form.Width != form.Repeat {
			return nil, fmt.Errorf("TForm %q: writing rAw substring arrays is not supported", c.TForm)
		}
	}
	if vla {
		if c.TForm != "" && form.Repeat != 1 {
			return nil, fmt.Errorf("TForm %q: the writer needs a VLA repeat of 1", c.TForm)
		}
		s.q = c.TForm != "" && form.Type == tform.BinQArrayDesc64
		s.nrows = int64(n)
		s.cellBytes = 8
		if s.q {
			s.cellBytes = 16
		}
		if c.Valid != nil {
			return nil, fmt.Errorf("use VarValid, not Valid, for a VLA column")
		}
		if c.VarValid != nil {
			if len(c.VarValid) != n {
				return nil, fmt.Errorf("VarValid has %d rows, data %d", len(c.VarValid), n)
			}
			for r := 0; r < n; r++ {
				if len(c.VarValid[r]) != varRowLen(c, r) {
					return nil, fmt.Errorf("VarValid row %d has %d flags, data %d", r, len(c.VarValid[r]), varRowLen(c, r))
				}
			}
		}
		for r := 0; r < n; r++ {
			if l := varRowLen(c, r); l > s.maxLen {
				s.maxLen = l
			}
		}
		if c.TForm != "" && form.VarMaxLen != 0 && int(form.VarMaxLen) < s.maxLen {
			return nil, fmt.Errorf("TForm %q: max length %d < longest row %d", c.TForm, form.VarMaxLen, s.maxLen)
		}
		d := byte('P')
		if s.q {
			d = 'Q'
		}
		s.tform = fmt.Sprintf("1%c%c(%d)", d, code, s.maxLen)
		if len(c.Dim) > 0 {
			return nil, fmt.Errorf("TDIM on a VLA column is not supported by the writer")
		}
		return s, nil
	}
	if c.VarValid != nil {
		return nil, fmt.Errorf("use Valid, not VarValid, for a fixed-width column")
	}
	repeat := int64(1)
	switch {
	case c.TForm != "":
		repeat = form.Repeat
	case len(c.Dim) > 0:
		for _, d := range c.Dim {
			repeat *= d
		}
	case code == tform.BinChar:
		repeat = 0
		for _, str := range c.DataString {
			repeat = max(repeat, int64(len(str)))
		}
		repeat = max(repeat, 1)
	}
	if len(c.Dim) > 0 {
		p := int64(1)
		for _, d := range c.Dim {
			if d < 0 {
				return nil, fmt.Errorf("negative TDIM axis %d", d)
			}
			p *= d
		}
		if p > repeat {
			return nil, fmt.Errorf("Dim %v has %d elements, more than the repeat count %d", c.Dim, p, repeat)
		}
	}
	s.repeat = repeat
	if c.Valid != nil && len(c.Valid) != n {
		return nil, fmt.Errorf("Valid has %d flags, data %d elements", len(c.Valid), n)
	}
	if code == tform.BinChar {
		if len(c.Dim) > 0 && c.Dim[0] != repeat {
			// TDIM (w, d2, ...) with w < repeat declares an array of
			// strings (cfitsio and astropy split the cell); DataString
			// holds one string per row.
			return nil, fmt.Errorf("Dim %v on a character column declares an array of %d-character strings; writing string arrays is not supported (DataString holds one string per row)", c.Dim, c.Dim[0])
		}
		s.nrows = int64(n)
		for r, str := range c.DataString {
			if int64(len(str)) > repeat {
				return nil, fmt.Errorf("row %d: string of %d bytes does not fit TFORM %dA", r, len(str), repeat)
			}
			if err := printableASCII(str); err != nil {
				return nil, fmt.Errorf("row %d: %w", r, err)
			}
			if c.Valid != nil && !c.Valid[r] {
				return nil, fmt.Errorf("row %d: character columns have no null value", r)
			}
		}
		s.tform = strconv.FormatInt(repeat, 10) + "A"
		s.cellBytes = repeat
		return s, nil
	}
	switch {
	case repeat > 0:
		if int64(n)%repeat != 0 {
			return nil, fmt.Errorf("%d elements is not a multiple of the repeat count %d", n, repeat)
		}
		s.nrows = int64(n) / repeat
	case n != 0:
		return nil, fmt.Errorf("repeat count 0 but %d elements given", n)
	default:
		s.nrows = -1 // a zero-width column holds nothing in any row
	}
	s.tform = strconv.FormatInt(repeat, 10) + string(rune(code))
	if code == tform.BinBit {
		s.cellBytes = (repeat + 7) / 8
	} else {
		s.cellBytes = repeat * int64(code.ElementSize())
	}
	return s, nil
}

func printableASCII(s string) error {
	for i := 0; i < len(s); i++ {
		if s[i] < 0x20 || s[i] > 0x7e {
			return fmt.Errorf("byte 0x%02x is not printable ASCII (§7.3.3.1)", s[i])
		}
	}
	return nil
}

func varRowLen(c *ColumnData, r int) int {
	switch {
	case c.DataVarUint8 != nil:
		return len(c.DataVarUint8[r])
	case c.DataVarInt16 != nil:
		return len(c.DataVarInt16[r])
	case c.DataVarInt32 != nil:
		return len(c.DataVarInt32[r])
	case c.DataVarInt64 != nil:
		return len(c.DataVarInt64[r])
	case c.DataVarFloat32 != nil:
		return len(c.DataVarFloat32[r])
	case c.DataVarFloat64 != nil:
		return len(c.DataVarFloat64[r])
	case c.DataVarComplex64 != nil:
		return len(c.DataVarComplex64[r])
	case c.DataVarComplex128 != nil:
		return len(c.DataVarComplex128[r])
	case c.DataVarBool != nil:
		return len(c.DataVarBool[r])
	case c.DataVarBits != nil:
		return len(c.DataVarBits[r])
	case c.DataVarString != nil:
		return len(c.DataVarString[r])
	}
	return 0
}

func storedFits(code tform.BinaryType, v int64) bool {
	switch code {
	case tform.BinUint8:
		return v >= 0 && v <= math.MaxUint8
	case tform.BinInt16:
		return v >= math.MinInt16 && v <= math.MaxInt16
	case tform.BinInt32:
		return v >= math.MinInt32 && v <= math.MaxInt32
	}
	return true
}

// putInt stores one integer element (already a stored value).
func (s *colSpec) putInt(dst []byte, v int64, valid bool) error {
	tnull := s.col.TNull
	if !valid {
		if tnull == nil {
			return fmt.Errorf("null element but the column has no TNULL")
		}
		v = *tnull
	} else if tnull != nil && v == *tnull {
		return fmt.Errorf("stored value %d equals TNULL and would read back as null", v)
	}
	switch s.code {
	case tform.BinUint8:
		dst[0] = byte(v)
	case tform.BinInt16:
		bigendian.PutInt16(dst, int16(v))
	case tform.BinInt32:
		bigendian.PutInt32(dst, int32(v))
	case tform.BinInt64:
		bigendian.PutInt64(dst, v)
	}
	return nil
}

func validAt(valid []bool, i int) bool { return valid == nil || valid[i] }

// encodeFixed writes row r of a fixed-width column into cell.
func (s *colSpec) encodeFixed(cell []byte, r int64) error {
	c := s.col
	if s.code == tform.BinChar {
		str := c.DataString[r]
		copy(cell, str)
		for i := len(str); i < len(cell); i++ {
			cell[i] = ' '
		}
		return nil
	}
	lo := int(r * s.repeat)
	return s.encodeElems(cell, lo, int(s.repeat), c.Valid, func(i int) any { return elemAt(c, i, -1) })
}

// encodeElems writes n elements, the first being element lo of the data
// (as returned by at), into dst.
func (s *colSpec) encodeElems(dst []byte, lo, n int, valid []bool, at func(i int) any) error {
	if s.code == tform.BinBit {
		for k := 0; k < n; k++ {
			if !validAt(valid, lo+k) {
				return fmt.Errorf("bit columns have no null value")
			}
			if at(lo + k).(bool) {
				dst[k/8] |= 0x80 >> (k % 8)
			}
		}
		return nil
	}
	es := s.code.ElementSize()
	for k := 0; k < n; k++ {
		ok := validAt(valid, lo+k)
		out := dst[k*es:]
		var err error
		switch x := at(lo + k).(type) {
		case uint8:
			err = s.putInt(out, int64(x), ok)
		case int8:
			err = s.putInt(out, int64(x)+128, ok)
		case int16:
			err = s.putInt(out, int64(x), ok)
		case uint16:
			err = s.putInt(out, int64(x)-32768, ok)
		case int32:
			err = s.putInt(out, int64(x), ok)
		case uint32:
			err = s.putInt(out, int64(x)-2147483648, ok)
		case int64:
			err = s.putInt(out, x, ok)
		case uint64:
			err = s.putInt(out, int64(x^(1<<63)), ok)
		case float32:
			if !ok {
				x = float32(math.NaN())
			}
			bigendian.PutFloat32(out, x)
		case float64:
			if !ok {
				x = math.NaN()
			}
			bigendian.PutFloat64(out, x)
		case complex64:
			if !ok {
				return fmt.Errorf("complex columns have no null value")
			}
			bigendian.PutFloat32(out, real(x))
			bigendian.PutFloat32(out[4:], imag(x))
		case complex128:
			if !ok {
				return fmt.Errorf("complex columns have no null value")
			}
			bigendian.PutFloat64(out, real(x))
			bigendian.PutFloat64(out[8:], imag(x))
		case bool:
			switch {
			case !ok:
				out[0] = 0
			case x:
				out[0] = 'T'
			default:
				out[0] = 'F'
			}
		}
		if err != nil {
			return fmt.Errorf("element %d: %w", lo+k, err)
		}
	}
	return nil
}

// elemAt returns element i of the column's flat fixed data (row < 0) or
// element i of VLA row `row`.
func elemAt(c *ColumnData, i, row int) any {
	if row < 0 {
		switch {
		case c.DataUint8 != nil:
			return c.DataUint8[i]
		case c.DataInt8 != nil:
			return c.DataInt8[i]
		case c.DataInt16 != nil:
			return c.DataInt16[i]
		case c.DataUint16 != nil:
			return c.DataUint16[i]
		case c.DataInt32 != nil:
			return c.DataInt32[i]
		case c.DataUint32 != nil:
			return c.DataUint32[i]
		case c.DataInt64 != nil:
			return c.DataInt64[i]
		case c.DataUint64 != nil:
			return c.DataUint64[i]
		case c.DataFloat32 != nil:
			return c.DataFloat32[i]
		case c.DataFloat64 != nil:
			return c.DataFloat64[i]
		case c.DataComplex64 != nil:
			return c.DataComplex64[i]
		case c.DataComplex128 != nil:
			return c.DataComplex128[i]
		case c.DataBool != nil:
			return c.DataBool[i]
		case c.DataBits != nil:
			return c.DataBits[i]
		}
		return nil
	}
	switch {
	case c.DataVarUint8 != nil:
		return c.DataVarUint8[row][i]
	case c.DataVarInt16 != nil:
		return c.DataVarInt16[row][i]
	case c.DataVarInt32 != nil:
		return c.DataVarInt32[row][i]
	case c.DataVarInt64 != nil:
		return c.DataVarInt64[row][i]
	case c.DataVarFloat32 != nil:
		return c.DataVarFloat32[row][i]
	case c.DataVarFloat64 != nil:
		return c.DataVarFloat64[row][i]
	case c.DataVarComplex64 != nil:
		return c.DataVarComplex64[row][i]
	case c.DataVarComplex128 != nil:
		return c.DataVarComplex128[row][i]
	case c.DataVarBool != nil:
		return c.DataVarBool[row][i]
	case c.DataVarBits != nil:
		return c.DataVarBits[row][i]
	}
	return nil
}

// encodeVar appends row r's payload to heap and writes its descriptor
// into cell.
func (s *colSpec) encodeVar(cell []byte, r int64, heap []byte) ([]byte, error) {
	c := s.col
	n := varRowLen(c, int(r))
	off := int64(len(heap))
	var payload []byte
	var valid []bool
	if c.VarValid != nil {
		valid = c.VarValid[r]
	}
	switch s.code {
	case tform.BinChar:
		str := c.DataVarString[r]
		if err := printableASCII(str); err != nil {
			return nil, err
		}
		for _, ok := range valid {
			if !ok {
				return nil, fmt.Errorf("character columns have no null value")
			}
		}
		payload = []byte(str)
	case tform.BinBit:
		payload = make([]byte, (n+7)/8)
	default:
		payload = make([]byte, n*s.code.ElementSize())
	}
	if s.code != tform.BinChar {
		row := int(r)
		if err := s.encodeElems(payload, 0, n, valid, func(i int) any { return elemAt(c, i, row) }); err != nil {
			return nil, err
		}
	}
	if n == 0 {
		off = 0
	}
	if err := putDescriptor(cell, s.q, s.code, int64(n), off, int64(len(heap))+int64(len(payload))); err != nil {
		return nil, err
	}
	return append(heap, payload...), nil
}

// putDescriptor writes a P (8-byte) or Q (16-byte) array descriptor for n
// elements at heap offset off; end is the heap size after this row. A 1P
// descriptor holds 32-bit signed fields, so a row that would need a length
// or an offset beyond 2^31-1 is an error, never a wrapped descriptor.
func putDescriptor(cell []byte, q bool, code tform.BinaryType, n, off, end int64) error {
	if q {
		bigendian.PutInt64(cell, n)
		bigendian.PutInt64(cell[8:], off)
		return nil
	}
	if end > math.MaxInt32 || n > math.MaxInt32 {
		return fmt.Errorf("heap of %d bytes exceeds the 1P descriptor range; use TForm \"1Q%c\"", end, code)
	}
	bigendian.PutInt32(cell, int32(n))
	bigendian.PutInt32(cell[4:], int32(off))
	return nil
}

// writeKeys emits column n's TFORMn and descriptive keywords.
func (s *colSpec) writeKeys(h *headerBuilder, n int) {
	c := s.col
	k := strconv.Itoa(n)
	h.set("TFORM"+k, s.tform, "")
	if c.Name != "" {
		h.set("TTYPE"+k, c.Name, "")
	}
	if c.Unit != "" {
		h.set("TUNIT"+k, c.Unit, "")
	}
	if c.Display != "" {
		h.set("TDISP"+k, c.Display, "")
	}
	if len(c.Dim) > 0 {
		dims := make([]string, len(c.Dim))
		for j, d := range c.Dim {
			dims[j] = strconv.FormatInt(d, 10)
		}
		h.set("TDIM"+k, "("+strings.Join(dims, ",")+")", "")
	}
	if c.TNull != nil {
		h.set("TNULL"+k, *c.TNull, "")
	}
	if c.TScal != nil {
		h.set("TSCAL"+k, numberValue(*c.TScal), "")
	}
	switch {
	case s.autoZero == "9223372036854775808":
		// 2^63 does not fit an int64 card value: it is written as the
		// exact integer literal (a TypeBigInt card), as cfitsio and
		// astropy write it. A float literal such as
		// 9.223372036854776E+18 would be a different decimal.
		h.set("TZERO"+k, uint64(1)<<63, "offset for unsigned integers")
	case s.autoZero != "":
		z, _ := strconv.ParseInt(s.autoZero, 10, 64)
		h.set("TZERO"+k, z, "offset for signed/unsigned integers")
	case c.TZero != nil:
		h.set("TZERO"+k, numberValue(*c.TZero), "")
	}
}

// maxValueChars is the room for a card value with no comment (bytes 11-80).
const maxValueChars = 70

// numberValue returns the card value written for a TSCAL/TZERO: an
// integral value as an exact integer literal (int64, or *big.Int beyond
// it), so that the literal's decimal value is the float64 itself, and any
// other value as the float64 (written as its shortest decimal).
func numberValue(v float64) any {
	if v != math.Trunc(v) || math.IsInf(v, 0) {
		return v
	}
	if math.Abs(v) < 1<<63 {
		return int64(v)
	}
	b, _ := new(big.Float).SetFloat64(v).Int(nil)
	return b
}

// writtenLiteral is the exact decimal value of the card numberValue(v)
// produces.
func writtenLiteral(v float64) *big.Rat {
	switch x := numberValue(v).(type) {
	case int64:
		return new(big.Rat).SetInt64(x)
	case *big.Int:
		return new(big.Rat).SetInt(x)
	}
	r, _ := new(big.Rat).SetString(strconv.FormatFloat(v, 'G', -1, 64))
	return r
}

// mergeColumnKeys applies hdr's per-column keywords to cols (in place)
// and returns the remaining cards to copy into the new header.
func mergeColumnKeys(hdr *header.Header, cols []ColumnData) ([]header.Card, error) {
	if hdr == nil {
		return nil, nil
	}
	var extra []header.Card
	for _, c := range hdr.Cards() {
		key := c.Key
		if isStructuralTableKey(key) {
			continue
		}
		root, n, ok := indexedKey(key)
		switch root {
		case "TTYPE", "TUNIT", "TDISP", "TDIM", "TSCAL", "TZERO", "TNULL", "TFORM":
		default:
			ok = false
		}
		if !ok {
			extra = append(extra, c)
			continue
		}
		if n < 1 || n > len(cols) {
			return nil, fmt.Errorf("fits: AppendBinaryTable: hdr has %s but the table has %d columns", key, len(cols))
		}
		col := &cols[n-1]
		conflict := func(field string, have any) error {
			return fmt.Errorf("fits: AppendBinaryTable: hdr %s = %v conflicts with column %d's %s %v", key, c.Value, n, field, have)
		}
		switch root {
		case "TTYPE", "TUNIT", "TDISP", "TFORM", "TDIM":
			v, isStr := c.Value.(string)
			if !isStr {
				return nil, fmt.Errorf("fits: AppendBinaryTable: hdr %s is not a string", key)
			}
			switch root {
			case "TTYPE", "TUNIT", "TDISP":
				dst := map[string]*string{"TTYPE": &col.Name, "TUNIT": &col.Unit, "TDISP": &col.Display}[root]
				if *dst != "" && *dst != v {
					return nil, conflict(root, *dst)
				}
				*dst = v
			case "TFORM":
				v = strings.TrimSpace(v)
				if col.TForm != "" && strings.TrimSpace(col.TForm) != v {
					return nil, conflict("TForm", col.TForm)
				}
				col.TForm = v
			case "TDIM":
				d, err := tform.ParseDim(v)
				if err != nil {
					return nil, fmt.Errorf("fits: AppendBinaryTable: hdr %s: %w", key, err)
				}
				if len(col.Dim) > 0 && !equalDims(col.Dim, d) {
					return nil, conflict("Dim", col.Dim)
				}
				col.Dim = d
			}
		case "TSCAL", "TZERO":
			v, err := cardNumber(c)
			if err != nil {
				return nil, fmt.Errorf("fits: AppendBinaryTable: hdr %s: %w", key, err)
			}
			if lit := exactCardValue(c); lit == nil || lit.Cmp(writtenLiteral(v)) != 0 {
				// e.g. TZERO = 9223372036854775809, whose nearest double is
				// 2^63: writing the double back would change the column.
				return nil, fmt.Errorf("fits: AppendBinaryTable: hdr %s = %v cannot be written back exactly (the writer holds TSCAL/TZERO as float64)", key, c.Value)
			}
			dst := &col.TScal
			if root == "TZERO" {
				dst = &col.TZero
			}
			if *dst != nil && **dst != v {
				return nil, conflict(root, **dst)
			}
			*dst = &v
		case "TNULL":
			v, isInt := c.Value.(int64)
			if !isInt {
				return nil, fmt.Errorf("fits: AppendBinaryTable: hdr %s is not an integer", key)
			}
			if col.TNull != nil && *col.TNull != v {
				return nil, conflict("TNull", *col.TNull)
			}
			col.TNull = &v
		}
	}
	return extra, nil
}

func equalDims(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// isStructuralTableKey reports whether key is computed by the table
// writers (and so never copied from a user header).
func isStructuralTableKey(key string) bool {
	switch key {
	case header.KeyXtension, header.KeyBitpix, header.KeyNaxis, header.KeyPcount,
		header.KeyGcount, header.KeyTfields, header.KeyEnd, "THEAP", "SIMPLE", "EXTEND":
		return true
	}
	root, _, ok := indexedKey(key)
	return ok && (root == "NAXIS" || root == "TBCOL")
}

// indexedKey splits "TFORM12" into ("TFORM", 12).
func indexedKey(key string) (root string, n int, ok bool) {
	i := len(key)
	for i > 0 && key[i-1] >= '0' && key[i-1] <= '9' {
		i--
	}
	if i == len(key) || i == 0 || len(key)-i > 3 {
		return "", 0, false
	}
	n, err := strconv.Atoi(key[i:])
	if err != nil {
		return "", 0, false
	}
	return key[:i], n, true
}
