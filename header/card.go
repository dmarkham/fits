// Package header implements FITS header-card parsing and serialization.
//
// A FITS header is a sequence of fixed-width 80-byte "cards" grouped into
// 2880-byte blocks. Every card has a keyword in columns 1–8 (or the ESO
// HIERARCH form, "HIERARCH " followed by a longer keyword path), optionally
// followed by a value indicator and a value field that may carry a string,
// logical, integer, floating-point, or complex value, followed by an
// optional "/ comment" field.
//
// Reading follows cfitsio: a card is split into value and comment exactly
// as cfitsio 4.7.0's ffpsvc does, the value of a long string is joined
// across CONTINUE cards as fits_read_key_longstr joins it (its joined
// comment too, except that it is kept whole where cfitsio cuts it at 72
// characters; see PROVENANCE.md), and a card whose value is not valid
// FITS does not fail the header: it is kept (TypeUnparsed, with
// the raw value text) next to every other card. The test suite checks the
// decoded value, type and comment of every card of every sample against
// cfitsio reading the same bytes (header/testdata, testdata/cfref).
//
// Writing: a card read from a file and not replaced is written as the
// bytes it was read from; a copy of such a card whose fields were edited
// is written from its fields (see EncodeCardWithContinuation). A card
// written from its fields has its keyword upper-cased, as cfitsio's
// fits_write_key writes it (FromCards and AddCard store it so). A
// decoded card with a lower-case keyword is not copied into another
// header as its bytes: AddCard refuses it (readers disagree on such a
// card), and a keyword that would only upper-case to END is never
// written as the END card. A card built in memory (Set, Add, AddCard) is
// written so that cfitsio and
// astropy read back the same value and comment (checked by the tests
// through both, plus fitsverify), or Set / Encode return an error: the
// encoder never truncates a value or comment, or rewrites one to make it
// fit. Trailing blanks are not significant in FITS and readers (cfitsio,
// astropy, this package) drop them, so the encoder drops them too: from
// every string value, from the comment of every value card, and the
// commentary text of COMMENT / HISTORY / blank cards reads back without
// them. A RawValue is written verbatim and is accepted only when it reads
// back as the same unparsed text (checked through cfitsio and this
// package; astropy does not read invalid values). Encode also reads its
// own output back and refuses a header whose cards would join or split on
// read.
//
// Reference: FITS Standard v4.0 §4.1, §4.2, Appendix A.
package header

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"
)

// CardWidth is the fixed width of a header card, in bytes.
const CardWidth = 80

// ValueType enumerates the distinct value types a parsed card can carry.
type ValueType int

const (
	// TypeEmpty is a card with no value: commentary (COMMENT, HISTORY,
	// blank keyword, a CONTINUE card that does not continue a string),
	// a card without a value indicator, or a keyword whose value field is
	// blank ("undefined value").
	TypeEmpty ValueType = iota
	// TypeString — single-quoted string value (§4.2.1). Value is a string.
	TypeString
	// TypeLogical — "T" or "F" (§4.2.2). Value is a bool.
	TypeLogical
	// TypeInt — integer value within int64 (§4.2.3). Value is an int64.
	TypeInt
	// TypeFloat — floating-point value including the D-exponent form
	// (§4.2.4). Value is a float64.
	TypeFloat
	// TypeComplexInt — integer complex "(re,im)" (§4.2.5). Value is a Complex.
	TypeComplexInt
	// TypeComplexFloat — floating complex "(re,im)" (§4.2.6). Value is a Complex.
	TypeComplexFloat
	// TypeBigInt — an integer literal outside the int64 range, e.g. the
	// BZERO = 9223372036854775808 cfitsio writes for uint64 images. Value
	// is a non-nil *big.Int holding it exactly; a card written from its
	// fields with a nil Value or one inside int64 (which is TypeInt) is
	// an error.
	TypeBigInt
	// TypeUnparsed — the card has a value field that is not a valid FITS
	// value (an unquoted string such as DATE-OBS= 2020-01-01, "TRUE", a
	// number that overflows float64, NaN, ...). Value is a RawValue with
	// the value text exactly as cfitsio's value parser returns it. The card
	// is kept, re-encoded byte for byte from Raw, and every typed getter
	// refuses it with ErrWrongType.
	TypeUnparsed
)

func (t ValueType) String() string {
	switch t {
	case TypeEmpty:
		return "empty"
	case TypeString:
		return "string"
	case TypeLogical:
		return "logical"
	case TypeInt:
		return "int"
	case TypeFloat:
		return "float"
	case TypeComplexInt:
		return "complex-int"
	case TypeComplexFloat:
		return "complex-float"
	case TypeBigInt:
		return "big-int"
	case TypeUnparsed:
		return "unparsed"
	}
	return "ValueType(" + strconv.Itoa(int(t)) + ")"
}

// RawValue is the Value of a TypeUnparsed card.
type RawValue struct {
	// Text is the value field as cfitsio's value parser (ffpsvc) splits it
	// from the card: the first blank- or '/'-terminated token, or a quoted
	// or parenthesised group.
	Text string
	// Reason says why Text is not a valid FITS value.
	Reason string
}

// Card is a single parsed FITS header card.
//
// Key holds the keyword name as cfitsio's ffgknm reads it: columns 1 up
// to the first blank or '=' (columns 1–8 without trailing blanks for a
// standard card; a malformed card such as "B9  =    5" gives "B9", one
// with no blank or '=' in its first 74 columns gives those 74), or
// "HIERARCH " plus the keyword path for HIERARCH cards. A HIERARCH path is
// stored with single blanks ("HIERARCH ESO DET NAME") even when the card
// has runs of blanks between its words; cfitsio keeps the runs, so here a
// lookup works with the single-blank form only.
//
// Value holds the typed value (see ValueType for the Go type per Type). For
// commentary cards (COMMENT, HISTORY, blank) Value is nil and Comment
// carries the text of columns 9–80.
//
// Raw preserves the original 80 bytes read from disk; for a long string
// continued over CONTINUE cards it is the first card and RawCards returns
// all of them. Encode re-emits the original bytes of a card read from disk
// unless it was replaced through Set, or it is a copy whose Key, Type,
// Value or Comment was changed (a change of letter case in Key, of the
// sign of a float zero, or of one part of a complex counts): such a copy
// is written from its fields (see EncodeCardWithContinuation), with its
// keyword upper-cased; switching Type between TypeComplexInt and
// TypeComplexFloat counts too, and a TypeComplexInt card must hold
// integer parts. Raw is zero for cards built in memory.
//
// An integer card written "-0" decodes to Value int64(0), and Float reads
// it as -0, as cfitsio's fits_read_key_dbl does. A copy of it keeps that
// -0 while its Value is int64(0) (one whose comment alone was edited is
// written "-0"); a copy whose Value was changed holds the new value. To
// write an integer +0, use Set.
type Card struct {
	Key     string
	Value   any
	Comment string
	Type    ValueType
	Raw     [CardWidth]byte

	// cont holds the CONTINUE cards absorbed into this card's string.
	cont [][CardWidth]byte
	// negZero records an integer card written "-0": its int64 value is 0,
	// but cfitsio's fits_read_key_dbl (strtod) reads it as -0, and so does
	// Header.Float.
	negZero bool
}

// RawCards returns every physical card this card was decoded from, in
// file order: Raw followed by the CONTINUE cards of a long string. It
// returns nil for a card built in memory.
func (c *Card) RawCards() [][CardWidth]byte {
	if c.Raw == ([CardWidth]byte{}) {
		return nil
	}
	out := make([][CardWidth]byte, 0, 1+len(c.cont))
	out = append(out, c.Raw)
	return append(out, c.cont...)
}

// IsCommentary reports whether c is a commentary card (COMMENT, HISTORY, or a
// blank keyword). Commentary cards never carry a typed value.
func (c *Card) IsCommentary() bool {
	return c.Key == "COMMENT" || c.Key == "HISTORY" || c.Key == ""
}

// IsEnd reports whether c is the END card that terminates the header.
func (c *Card) IsEnd() bool { return c.Key == "END" }

// IsHierarch reports whether c uses the ESO HIERARCH long-keyword convention.
func (c *Card) IsHierarch() bool {
	return strings.HasPrefix(c.Key, "HIERARCH ")
}

// Complex is the in-memory representation of a FITS complex value. The real
// and imaginary parts are stored as float64 regardless of whether the card
// wrote them as integer or floating literals — callers that need the
// original textual form can inspect Card.Raw.
type Complex struct {
	Re, Im float64
}

// String formats v for diagnostic output; the canonical serializer lives in
// encoder.go.
func (v Complex) String() string {
	return fmt.Sprintf("(%g,%g)", v.Re, v.Im)
}

// ------------------------- card decode -------------------------

// ErrBadCard is returned when a card's bytes cannot be treated as a card at
// all (wrong length). A card whose value is malformed is not an error: it
// decodes as TypeUnparsed.
var ErrBadCard = errors.New("fits/header: malformed card")

// cfitsio buffer limits that shape what its parser returns (fitsio.h):
// FLEN_VALUE-1 value characters, FLEN_COMMENT-1 comment characters.
const (
	maxValueLen   = 70
	maxCommentLen = 72
)

// DecodeCard parses one 80-byte FITS header card and returns the typed Card.
// It returns an error only when len(raw) != CardWidth.
//
// The value and comment are split from the card exactly as cfitsio's
// ffpsvc splits them:
//
//   - COMMENT, HISTORY, blank keyword, END and CONTINUE cards have no value;
//     Comment is columns 9–80 without trailing blanks.
//   - "HIERARCH " cards: the value follows the first '='; without '=' the
//     card has Key "HIERARCH", no value, and Comment = columns 9–80.
//   - otherwise the value follows "= " in columns 9–10, or failing that the
//     first '=' anywhere on the card; a card without '=' has no value and
//     Comment is columns 9–80.
//   - the value is a quoted string (an unterminated quote runs to the end of
//     the card), a parenthesised complex, or the first token ending at a
//     blank or '/'. The comment is the rest of the card after the value,
//     with one leading '/' and one blank after it removed.
//
// The value is then typed: a quoted string is TypeString (doubled quotes
// unescaped, trailing blanks removed); exactly "T" or "F" is TypeLogical;
// an integer literal is TypeInt, or TypeBigInt outside int64; a number
// with '.', 'E', 'e' or an upper-case 'D' exponent is TypeFloat; "(re, im)"
// is a complex type; anything else, or a number that is not finite in
// float64 (overflow, NaN, Inf, underflow to zero or a subnormal, which
// cfitsio reports as overflow), is TypeUnparsed.
//
// Where cfitsio is laxer than FITS, this decoder keeps the text instead of
// guessing: cfitsio reads "TRUE" or "Tx" as logical true, a hexadecimal
// token as a number, and a complex value with an empty part ("(1, )",
// "(,1)") as having a zero there; here they are TypeUnparsed.
func DecodeCard(raw []byte) (Card, error) {
	if len(raw) != CardWidth {
		return Card{}, fmt.Errorf("%w: card length %d != %d", ErrBadCard, len(raw), CardWidth)
	}
	return decodeCard((*[CardWidth]byte)(raw)), nil
}

// decodeCard is DecodeCard for a card image of the right length.
func decodeCard(raw *[CardWidth]byte) Card {
	var c Card
	c.Raw = *raw
	// cfitsio works on the card with trailing blanks removed.
	t := strings.TrimRight(string(raw[:]), " ")
	key8 := string(raw[:8])

	if key8 == "END     " {
		c.Key = "END"
		return c
	}
	if strings.HasPrefix(t, "HIERARCH ") {
		eq := strings.IndexByte(t, '=')
		if eq < 0 {
			// No value indicator: cfitsio names the card "HIERARCH" and
			// reports columns 9-80 as its comment.
			c.Key = "HIERARCH"
			c.Comment = commentaryText(raw[:])
			return c
		}
		c.Key = "HIERARCH " + strings.Join(strings.Fields(t[9:eq]), " ")
		decodeValueField(t, eq+1, &c)
		return c
	}
	c.Key = keywordName(t)
	switch key8 {
	case "COMMENT ", "HISTORY ", "        ", "CONTINUE":
		c.Comment = commentaryText(raw[:])
		return c
	}
	valpos := 10
	if len(t) < 10 || t[8:10] != "= " {
		eq := strings.IndexByte(t, '=')
		if eq < 0 {
			c.Comment = commentaryText(raw[:])
			return c
		}
		valpos = eq + 1
	}
	decodeValueField(t, valpos, &c)
	return c
}

// keywordName is cfitsio's ffgknm for a card that is not HIERARCH: the
// name runs from column 1 to the first blank or '=', for at most
// FLEN_KEYWORD-1 (74) characters. For a standard card that is columns 1-8
// without trailing blanks; for a malformed one it can be shorter
// ("B9  =    5" names B9) or longer than 8 characters.
func keywordName(t string) string {
	const maxKeyLen = 74
	n := min(len(t), maxKeyLen)
	if i := strings.IndexAny(t[:n], " ="); i >= 0 {
		return t[:i]
	}
	return t[:n]
}

func commentaryText(raw []byte) string {
	s := strings.TrimRight(string(raw[8:]), " ")
	if len(s) > maxCommentLen {
		s = s[:maxCommentLen]
	}
	return s
}

// splitValue is cfitsio's ffpsvc from the value position on: it returns the
// value string (quotes kept for strings) and the comment. t is the card
// without trailing blanks. ok is false where ffpsvc reports an error (a
// complex value without a closing parenthesis).
func splitValue(t string, valpos int) (value, comment string, ok bool) {
	n := len(t)
	ii := valpos
	for ii < n && t[ii] == ' ' {
		ii++
	}
	if ii >= n {
		return "", "", true // undefined value
	}
	switch t[ii] {
	case '/':
		ii++
	case '\'':
		var v [maxValueLen + 1]byte
		v[0] = '\''
		jj := 1
		for ii++; ii < n && jj < maxValueLen; ii, jj = ii+1, jj+1 {
			if t[ii] == '\'' {
				if ii+1 < n && t[ii+1] == '\'' {
					v[jj] = t[ii]
					ii++
					jj++
				} else {
					v[jj] = t[ii]
					break
				}
			}
			if jj < maxValueLen {
				v[jj] = t[ii]
			}
		}
		if ii == n || jj >= maxValueLen {
			// No closing quote: cfitsio closes the string at the end of the
			// card (since 2008 this is not an error).
			jj = min(jj, maxValueLen-1)
			v[jj] = '\''
			value = string(v[:jj+1])
		} else {
			value = string(v[:jj+1])
			ii++
		}
	case '(':
		k := strings.IndexByte(t[ii:], ')')
		if k < 0 || k >= maxValueLen {
			return "", "", false
		}
		value = t[ii : ii+k+1]
		ii += k + 1
	default:
		k := strings.IndexAny(t[ii:], " /")
		if k < 0 {
			k = n - ii
		}
		if k > maxValueLen {
			k = maxValueLen
		}
		value = t[ii : ii+k]
		ii += k
	}
	for ii < n && t[ii] == ' ' {
		ii++
	}
	if ii < n {
		if t[ii] == '/' {
			ii++
			if ii < n && t[ii] == ' ' {
				ii++
			}
		}
		comment = strings.TrimRight(t[ii:], " ")
		if len(comment) > maxCommentLen {
			comment = comment[:maxCommentLen]
		}
	}
	return value, comment, true
}

// decodeValueField splits and types the value of a keyed card.
func decodeValueField(t string, valpos int, c *Card) {
	value, comment, ok := splitValue(t, valpos)
	if !ok {
		c.Type = TypeUnparsed
		c.Value = RawValue{Text: strings.TrimSpace(t[valpos:]), Reason: "complex value without a closing ')'"}
		return
	}
	c.Comment = comment
	typeValue(value, c)
}

// typeValue sets c.Type and c.Value from cfitsio's value string.
func typeValue(v string, c *Card) {
	unparsed := func(reason string) {
		c.Type = TypeUnparsed
		c.Value = RawValue{Text: v, Reason: reason}
	}
	switch {
	case v == "":
		c.Type = TypeEmpty
		c.Value = nil
	case v[0] == '\'':
		c.Type = TypeString
		c.Value = unquote(v)
	case v == "T" || v == "F":
		c.Type = TypeLogical
		c.Value = v == "T"
	case v[0] == '(':
		cv, ct, err := parseComplex(v)
		if err != nil {
			unparsed(err.Error())
			return
		}
		c.Type, c.Value = ct, cv
	case isIntegerLiteral(v):
		if i, err := strconv.ParseInt(v, 10, 64); err == nil {
			c.Type, c.Value = TypeInt, i
			c.negZero = i == 0 && v[0] == '-'
			return
		}
		b, _ := new(big.Int).SetString(strings.TrimPrefix(v, "+"), 10)
		c.Type, c.Value = TypeBigInt, b
	default:
		f, err := parseFITSFloat(v)
		if err != nil {
			unparsed(err.Error())
			return
		}
		c.Type, c.Value = TypeFloat, f
	}
}

// unquote is cfitsio's ffc2s for a value that starts with a quote (splitValue
// guarantees a closing quote): strip the quotes, turn each doubled quote into one, and drop
// trailing blanks.
func unquote(v string) string {
	var sb strings.Builder
	for i := 1; i < len(v); i++ {
		if v[i] == '\'' {
			if i+1 < len(v) && v[i+1] == '\'' {
				i++
			} else {
				break
			}
		}
		sb.WriteByte(v[i])
	}
	return strings.TrimRight(sb.String(), " ")
}

// parseComplex parses "(re, im)". Integer parts give TypeComplexInt.
func parseComplex(v string) (Complex, ValueType, error) {
	inner := v[1 : len(v)-1]
	parts := strings.Split(inner, ",")
	if len(parts) != 2 {
		return Complex{}, 0, errors.New("complex value needs exactly one ','")
	}
	reTok, imTok := strings.TrimSpace(parts[0]), strings.TrimSpace(parts[1])
	if isIntegerLiteral(reTok) && isIntegerLiteral(imTok) {
		re, err1 := strconv.ParseInt(reTok, 10, 64)
		im, err2 := strconv.ParseInt(imTok, 10, 64)
		if err1 == nil && err2 == nil {
			return Complex{Re: float64(re), Im: float64(im)}, TypeComplexInt, nil
		}
	}
	re, err1 := parseFITSFloat(reTok)
	im, err2 := parseFITSFloat(imTok)
	if err1 != nil || err2 != nil {
		return Complex{}, 0, fmt.Errorf("complex parts %q, %q are not numbers", reTok, imTok)
	}
	return Complex{Re: re, Im: im}, TypeComplexFloat, nil
}

// isIntegerLiteral reports whether tok is a pure integer literal (optional
// leading + or -, then one or more digits).
func isIntegerLiteral(tok string) bool {
	if tok == "" {
		return false
	}
	i := 0
	if tok[0] == '+' || tok[0] == '-' {
		i = 1
	}
	if i == len(tok) {
		return false
	}
	for ; i < len(tok); i++ {
		if tok[i] < '0' || tok[i] > '9' {
			return false
		}
	}
	return true
}

// parseFITSFloat parses a FITS floating-point literal: decimal digits with
// an optional sign, point and E/e/D exponent (cfitsio replaces only an
// upper-case 'D'). The result must be a finite, normal (or zero from an
// all-zero mantissa) float64; cfitsio reports overflow, NaN, Inf and
// underflow as a conversion error and so does this function.
func parseFITSFloat(tok string) (float64, error) {
	for i := 0; i < len(tok); i++ {
		switch b := tok[i]; {
		case b >= '0' && b <= '9', b == '+', b == '-', b == '.', b == 'E', b == 'e', b == 'D':
		default:
			return 0, fmt.Errorf("%q is not a FITS number", tok)
		}
	}
	t := strings.Replace(tok, "D", "E", 1)
	f, err := strconv.ParseFloat(t, 64)
	if err != nil {
		return 0, fmt.Errorf("%q is not a FITS number", tok)
	}
	if math.IsInf(f, 0) || math.IsNaN(f) {
		return 0, fmt.Errorf("%q overflows float64", tok)
	}
	if a := math.Abs(f); a != 0 && a < 0x1p-1022 {
		return 0, fmt.Errorf("%q underflows float64", tok)
	}
	if f == 0 && hasNonzeroMantissa(t) {
		return 0, fmt.Errorf("%q underflows float64", tok)
	}
	return f, nil
}

func hasNonzeroMantissa(t string) bool {
	for i := 0; i < len(t); i++ {
		switch t[i] {
		case 'E', 'e':
			return false
		case '1', '2', '3', '4', '5', '6', '7', '8', '9':
			return true
		}
	}
	return false
}
