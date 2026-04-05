// Package tform parses the FITS TFORMn and TDIMn column-format keyword
// values for ASCII and binary tables.
//
// Binary-table TFORM (§7.3.2, Table 18):
//
//	 rTa[(dims)]
//
// where r is an optional repeat count, T is a single-letter type code, and
// "a" is only present for variable-length array descriptors (P/Q), e.g.:
//
//	"10E"        — 10 single-precision floats per row
//	"1PE(256)"   — variable-length float array, max 256 elements
//	"1QJ(1024)"  — variable-length int32 array with 64-bit descriptor
//
// ASCII-table TFORM (§7.2.3, Table 15):
//
//	 Tw[.d]
//
// where T is one of A, I, F, E, D, w is the field width, and d is the number
// of fraction digits (for F, E, D).
package tform

import (
	"errors"
	"fmt"
	"strconv"
	"strings"
)

// BinaryType is a binary-table column type code.
type BinaryType byte

// Binary TFORM type codes per Table 18.
const (
	BinLogical      BinaryType = 'L' // 1-byte boolean
	BinBit          BinaryType = 'X' // bit string
	BinUint8        BinaryType = 'B' // unsigned byte
	BinInt16        BinaryType = 'I' // 16-bit int
	BinInt32        BinaryType = 'J' // 32-bit int
	BinInt64        BinaryType = 'K' // 64-bit int
	BinChar         BinaryType = 'A' // character string
	BinFloat32      BinaryType = 'E' // single-precision float
	BinFloat64      BinaryType = 'D' // double-precision float
	BinComplex64    BinaryType = 'C' // complex single
	BinComplex128   BinaryType = 'M' // complex double
	BinPArrayDesc32 BinaryType = 'P' // 32-bit VLA descriptor
	BinQArrayDesc64 BinaryType = 'Q' // 64-bit VLA descriptor
)

// BinaryForm is the parsed result of a binary-table TFORMn value.
type BinaryForm struct {
	Raw    string     // original TFORM string
	Repeat int64      // leading repeat count; defaults to 1 if absent
	Type   BinaryType // single-letter type code
	// For P/Q descriptors: inner element type and optional max-length hint.
	VarType   BinaryType // element type inside the VLA (only for P/Q)
	VarMaxLen int64      // max length hint from "(n)" suffix, 0 if absent
}

// ElementSize returns the on-disk byte size of one element of the type.
// For P descriptors it returns 8 (two 32-bit ints); for Q descriptors 16
// (two 64-bit ints).
func (b BinaryType) ElementSize() int {
	switch b {
	case BinLogical, BinUint8, BinChar:
		return 1
	case BinInt16:
		return 2
	case BinInt32, BinFloat32:
		return 4
	case BinInt64, BinFloat64, BinComplex64:
		return 8
	case BinComplex128:
		return 16
	case BinBit:
		// Bit columns are measured in bits; callers compute byte cost as
		// ceil(r/8). ElementSize is undefined here.
		return 0
	case BinPArrayDesc32:
		return 8
	case BinQArrayDesc64:
		return 16
	}
	return 0
}

// Valid reports whether b is a known binary TFORM type.
func (b BinaryType) Valid() bool {
	switch b {
	case BinLogical, BinBit, BinUint8, BinInt16, BinInt32, BinInt64,
		BinChar, BinFloat32, BinFloat64, BinComplex64, BinComplex128,
		BinPArrayDesc32, BinQArrayDesc64:
		return true
	}
	return false
}

// ErrBadTForm is returned when a TFORM value cannot be parsed.
var ErrBadTForm = errors.New("fits/internal/tform: malformed TFORM")

// ParseBinary parses a binary-table TFORMn value. Leading whitespace is
// tolerated; trailing whitespace is stripped.
func ParseBinary(s string) (BinaryForm, error) {
	form := BinaryForm{Raw: s}
	trimmed := strings.TrimSpace(s)
	if trimmed == "" {
		return form, fmt.Errorf("%w: empty", ErrBadTForm)
	}
	// Leading repeat count.
	i := 0
	for i < len(trimmed) && trimmed[i] >= '0' && trimmed[i] <= '9' {
		i++
	}
	if i > 0 {
		n, err := strconv.ParseInt(trimmed[:i], 10, 64)
		if err != nil {
			return form, fmt.Errorf("%w: repeat count: %v", ErrBadTForm, err)
		}
		form.Repeat = n
	} else {
		form.Repeat = 1
	}
	if i >= len(trimmed) {
		return form, fmt.Errorf("%w: missing type code", ErrBadTForm)
	}
	form.Type = BinaryType(trimmed[i])
	if !form.Type.Valid() {
		return form, fmt.Errorf("%w: unknown type code %q", ErrBadTForm, trimmed[i])
	}
	i++
	// P/Q descriptors carry an inner type and an optional "(maxlen)" hint.
	if form.Type == BinPArrayDesc32 || form.Type == BinQArrayDesc64 {
		if i >= len(trimmed) {
			return form, fmt.Errorf("%w: P/Q missing element type", ErrBadTForm)
		}
		form.VarType = BinaryType(trimmed[i])
		if !form.VarType.Valid() || form.VarType == BinPArrayDesc32 || form.VarType == BinQArrayDesc64 {
			return form, fmt.Errorf("%w: invalid VLA element type %q", ErrBadTForm, trimmed[i])
		}
		i++
		// Optional "(n)" max length.
		if i < len(trimmed) && trimmed[i] == '(' {
			j := strings.IndexByte(trimmed[i:], ')')
			if j < 0 {
				return form, fmt.Errorf("%w: unterminated VLA max-length", ErrBadTForm)
			}
			inner := trimmed[i+1 : i+j]
			n, err := strconv.ParseInt(strings.TrimSpace(inner), 10, 64)
			if err != nil {
				return form, fmt.Errorf("%w: bad VLA max-length %q", ErrBadTForm, inner)
			}
			form.VarMaxLen = n
			i += j + 1
		}
	}
	return form, nil
}

// ASCIIType is an ASCII-table column type code.
type ASCIIType byte

// ASCII TFORM type codes per Table 15.
const (
	AsciiChar   ASCIIType = 'A'
	AsciiInt    ASCIIType = 'I'
	AsciiFloatF ASCIIType = 'F'
	AsciiFloatE ASCIIType = 'E'
	AsciiFloatD ASCIIType = 'D'
)

// ASCIIForm is the parsed result of an ASCII-table TFORMn value.
type ASCIIForm struct {
	Raw      string
	Type     ASCIIType
	Width    int // field width in characters
	Fraction int // fraction digits; 0 for A, I
}

// ParseASCII parses an ASCII-table TFORMn value.
func ParseASCII(s string) (ASCIIForm, error) {
	form := ASCIIForm{Raw: s}
	trimmed := strings.TrimSpace(s)
	if trimmed == "" {
		return form, fmt.Errorf("%w: empty", ErrBadTForm)
	}
	form.Type = ASCIIType(trimmed[0])
	switch form.Type {
	case AsciiChar, AsciiInt, AsciiFloatF, AsciiFloatE, AsciiFloatD:
	default:
		return form, fmt.Errorf("%w: unknown ASCII type %q", ErrBadTForm, trimmed[0])
	}
	rest := trimmed[1:]
	// Split on '.' for F/E/D. For A/I there is no fraction.
	wStr, fStr, _ := strings.Cut(rest, ".")
	if wStr == "" {
		return form, fmt.Errorf("%w: missing width", ErrBadTForm)
	}
	w, err := strconv.Atoi(strings.TrimSpace(wStr))
	if err != nil {
		return form, fmt.Errorf("%w: bad width %q: %v", ErrBadTForm, wStr, err)
	}
	form.Width = w
	if fStr != "" {
		f, err := strconv.Atoi(strings.TrimSpace(fStr))
		if err != nil {
			return form, fmt.Errorf("%w: bad fraction %q: %v", ErrBadTForm, fStr, err)
		}
		form.Fraction = f
	}
	return form, nil
}

// ParseDim parses a TDIMn value "(n1,n2,...,nk)" into a []int64 of axes.
// Leading and trailing whitespace is tolerated.
func ParseDim(s string) ([]int64, error) {
	t := strings.TrimSpace(s)
	if !strings.HasPrefix(t, "(") || !strings.HasSuffix(t, ")") {
		return nil, fmt.Errorf("fits/internal/tform: bad TDIM %q", s)
	}
	inner := t[1 : len(t)-1]
	parts := strings.Split(inner, ",")
	out := make([]int64, 0, len(parts))
	for _, p := range parts {
		n, err := strconv.ParseInt(strings.TrimSpace(p), 10, 64)
		if err != nil {
			return nil, fmt.Errorf("fits/internal/tform: bad TDIM axis %q", p)
		}
		out = append(out, n)
	}
	return out, nil
}
