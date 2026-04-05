// Package header implements FITS header-card parsing and serialization.
//
// A FITS header is a sequence of fixed-width 80-byte "cards" grouped into
// 2880-byte blocks. Every card has a keyword in columns 1–8 (or the extended
// HIERARCH form in columns 1–8 + columns 9..), optionally followed by the
// value indicator "= " in columns 9–10 and a value field that may carry a
// string, logical, integer, floating-point, or complex value, followed by
// an optional "/ comment" field.
//
// Reference: FITS Standard v4.0 §4.1, §4.2, Appendix A.
package header

import (
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"
)

// CardWidth is the fixed width of a header card, in bytes.
const CardWidth = 80

// ValueType enumerates the distinct value types a parsed card can carry.
type ValueType int

const (
	// TypeEmpty is a card with no value (COMMENT, HISTORY, blank, END, or a
	// commentary card with no "= " indicator).
	TypeEmpty ValueType = iota
	// TypeString — single-quoted string value (§4.2.1).
	TypeString
	// TypeLogical — "T" or "F" (§4.2.2).
	TypeLogical
	// TypeInt — integer value (§4.2.3).
	TypeInt
	// TypeFloat — floating-point value including D-exponent form (§4.2.4).
	TypeFloat
	// TypeComplexInt — integer complex "(re,im)" (§4.2.5).
	TypeComplexInt
	// TypeComplexFloat — floating complex "(re,im)" (§4.2.6).
	TypeComplexFloat
)

// Card is a single parsed FITS header card.
//
// Key holds the keyword in uppercase for normal cards, or the full HIERARCH
// keyword path ("HIERARCH ESO DET NAME") for HIERARCH cards.
//
// Value holds the typed value. For commentary cards (COMMENT, HISTORY, blank)
// Value is nil and Comment carries the text. For keyword cards without a
// value (i.e. badly-formed user cards that still contain "= "), Value is nil
// and Type is TypeEmpty.
//
// Raw preserves the original 80 bytes read from disk so that unmodified cards
// can be round-tripped bit-for-bit through Encode() without re-serialization.
// Raw is zero for cards built in memory.
type Card struct {
	Key     string
	Value   any
	Comment string
	Type    ValueType
	Raw     [CardWidth]byte
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

// ErrBadCard is returned when a card's 80 bytes cannot be parsed. The error
// text includes the keyword and a short description of the problem.
var ErrBadCard = errors.New("fits/header: malformed card")

// DecodeCard parses one 80-byte FITS header card and returns the typed Card.
// The caller guarantees len(raw) == CardWidth.
//
// Layout (§4.1, §4.2):
//
//	columns 1..8  : keyword (left-justified, space-padded)
//	columns 9..10 : value indicator "= " (present iff a value follows)
//	columns 11..80: value + comment field, or free-text for commentary cards
//
// HIERARCH cards use "HIERARCH" in columns 1..8 followed by a dotted/space
// keyword path and "= " elsewhere on the card (§4.1.2.2).
func DecodeCard(raw []byte) (Card, error) {
	if len(raw) != CardWidth {
		return Card{}, fmt.Errorf("%w: card length %d != %d", ErrBadCard, len(raw), CardWidth)
	}
	var c Card
	copy(c.Raw[:], raw)

	// Examine columns 1..8 as keyword.
	keyField := raw[0:8]
	key := strings.TrimRight(string(keyField), " ")

	// END is a special case: columns 1..3 = "END" and the rest is blanks. We
	// accept it if the keyword matches exactly "END".
	if key == "END" {
		c.Key = "END"
		c.Type = TypeEmpty
		return c, nil
	}

	// HIERARCH: key field reads "HIERARCH" and the real keyword follows.
	if key == "HIERARCH" {
		return decodeHierarch(raw, c)
	}

	// CONTINUE long-string continuation: "CONTINUE" in cols 1..8, then the
	// string value starts at column 11 with no "= " indicator (§4.2.1.2).
	if key == "CONTINUE" {
		c.Key = "CONTINUE"
		if err := decodeValue(raw[8:], &c); err != nil {
			return c, fmt.Errorf("%w: CONTINUE card: %v", ErrBadCard, err)
		}
		return c, nil
	}

	// Commentary cards: COMMENT, HISTORY, or entirely-blank keyword with free
	// text in columns 9..80.
	if key == "COMMENT" || key == "HISTORY" || key == "" {
		c.Key = key
		c.Type = TypeEmpty
		c.Comment = strings.TrimRight(string(raw[8:]), " ")
		// Leading spaces are part of the user text in principle; many writers
		// put one space after the keyword. We keep the text exactly as seen,
		// minus trailing blanks.
		return c, nil
	}

	// Normal keyed card. Columns 9–10 must contain "= " for the card to carry
	// a value; otherwise it is treated as commentary under §4.1.2.3.
	if raw[8] != '=' || raw[9] != ' ' {
		c.Key = key
		c.Type = TypeEmpty
		c.Comment = strings.TrimRight(string(raw[8:]), " ")
		return c, nil
	}

	c.Key = key
	if err := decodeValue(raw[10:], &c); err != nil {
		return c, fmt.Errorf("%w: key %q: %v", ErrBadCard, key, err)
	}
	return c, nil
}

// decodeHierarch decodes a HIERARCH card. The keyword path is the run of
// non-space tokens starting at column 9 and ending before the "=" sign. The
// value field starts at the first byte after "= ".
func decodeHierarch(raw []byte, c Card) (Card, error) {
	// Scan from column 8 (0-based: index 8) for the key path.
	body := string(raw[8:])
	before, after, ok := strings.Cut(body, "=")
	if !ok {
		// No value indicator: treat as commentary-ish HIERARCH card.
		c.Key = "HIERARCH " + strings.TrimRight(body, " ")
		c.Type = TypeEmpty
		return c, nil
	}
	path := strings.TrimSpace(before)
	// Collapse runs of spaces in the path to single spaces so that comparison
	// is canonical. FITS permits variable whitespace in HIERARCH paths.
	fields := strings.Fields(path)
	c.Key = "HIERARCH " + strings.Join(fields, " ")

	// valueField may lead with zero or more spaces.
	valueField := strings.TrimLeft(after, " ")
	// Re-wrap valueField into the byte slice form decodeValue expects by
	// prefixing enough spaces so the "comment boundary" math downstream has a
	// stable offset — simpler is to just pass a byte slice.
	if err := decodeValue([]byte(valueField), &c); err != nil {
		return c, fmt.Errorf("%w: HIERARCH %s: %v", ErrBadCard, c.Key, err)
	}
	return c, nil
}

// decodeValue parses the value+comment half of a card. The input is the
// bytes from column 11 onward (normal card) or from the first non-space byte
// after "=" (HIERARCH card).
func decodeValue(buf []byte, c *Card) error {
	s := string(buf)
	i := 0
	// Skip leading spaces (free-format values may be anywhere in the field).
	for i < len(s) && s[i] == ' ' {
		i++
	}
	if i >= len(s) {
		c.Type = TypeEmpty
		return nil
	}

	switch s[i] {
	case '\'':
		return decodeStringValue(s, i, c)
	case 'T', 'F':
		// Logical is a single character with optional trailing whitespace and
		// comment. Fixed format requires it at column 30 but free format allows
		// anywhere — we only check that the next char is a value terminator.
		if i+1 >= len(s) || s[i+1] == ' ' || s[i+1] == '/' {
			c.Type = TypeLogical
			c.Value = s[i] == 'T'
			c.Comment = parseInlineComment(s, i+1)
			return nil
		}
		// Otherwise fall through to numeric parse; rare but legal ("TRUE"
		// isn't a FITS logical, but some malformed cards put text here — we
		// reject below by treating as a numeric parse that fails).
		fallthrough
	default:
		return decodeNumericValue(s, i, c)
	}
}

// decodeStringValue parses a single-quoted FITS string value. Embedded
// single quotes are escaped by doubling ("''"). Trailing blanks before the
// closing quote are significant per §4.2.1.1 only if the string is meant to
// be length-preserving; we keep them as-is and let callers trim if desired.
func decodeStringValue(s string, i int, c *Card) error {
	if s[i] != '\'' {
		return errors.New("expected opening quote")
	}
	j := i + 1
	var sb strings.Builder
	for j < len(s) {
		if s[j] == '\'' {
			if j+1 < len(s) && s[j+1] == '\'' {
				sb.WriteByte('\'')
				j += 2
				continue
			}
			// Closing quote found.
			c.Type = TypeString
			c.Value = strings.TrimRight(sb.String(), " ")
			c.Comment = parseInlineComment(s, j+1)
			return nil
		}
		sb.WriteByte(s[j])
		j++
	}
	return errors.New("unterminated string value")
}

// decodeNumericValue parses an integer, float (including D-exponent), or
// complex value at s[i:]. It dispatches on the shape of the token rather
// than attempting every parse in sequence.
func decodeNumericValue(s string, i int, c *Card) error {
	// Complex literal: "(re,im)" with optional whitespace.
	if s[i] == '(' {
		return decodeComplexValue(s, i, c)
	}

	// Scan the token up to a comment slash or end-of-card.
	j := i
	for j < len(s) && s[j] != '/' {
		j++
	}
	tok := strings.TrimSpace(s[i:j])
	if tok == "" {
		c.Type = TypeEmpty
		c.Comment = parseInlineComment(s, j)
		return nil
	}

	// Integer: only digits (with optional leading sign).
	if isIntegerLiteral(tok) {
		v, err := strconv.ParseInt(tok, 10, 64)
		if err != nil {
			return fmt.Errorf("parse int %q: %w", tok, err)
		}
		c.Type = TypeInt
		c.Value = v
		c.Comment = parseInlineComment(s, j)
		return nil
	}

	// Float: allow "1.5", ".5", "1.", "1e5", "1E5", "1d5", "1D5".
	f, err := parseFITSFloat(tok)
	if err != nil {
		return fmt.Errorf("parse float %q: %w", tok, err)
	}
	c.Type = TypeFloat
	c.Value = f
	c.Comment = parseInlineComment(s, j)
	return nil
}

// decodeComplexValue parses "(re, im)" at s[i:]. Either form (integer or
// float components) is accepted; the resulting Type reflects which.
func decodeComplexValue(s string, i int, c *Card) error {
	// Find the closing ')'.
	j := strings.IndexByte(s[i:], ')')
	if j < 0 {
		return errors.New("unterminated complex value")
	}
	inner := s[i+1 : i+j]
	parts := strings.SplitN(inner, ",", 2)
	if len(parts) != 2 {
		return errors.New("complex value missing comma")
	}
	reTok := strings.TrimSpace(parts[0])
	imTok := strings.TrimSpace(parts[1])

	intForm := isIntegerLiteral(reTok) && isIntegerLiteral(imTok)
	if intForm {
		re, err1 := strconv.ParseInt(reTok, 10, 64)
		im, err2 := strconv.ParseInt(imTok, 10, 64)
		if err1 != nil || err2 != nil {
			return fmt.Errorf("parse integer complex (%s,%s)", reTok, imTok)
		}
		c.Type = TypeComplexInt
		c.Value = Complex{Re: float64(re), Im: float64(im)}
	} else {
		re, err1 := parseFITSFloat(reTok)
		im, err2 := parseFITSFloat(imTok)
		if err1 != nil || err2 != nil {
			return fmt.Errorf("parse float complex (%s,%s)", reTok, imTok)
		}
		c.Type = TypeComplexFloat
		c.Value = Complex{Re: re, Im: im}
	}
	c.Comment = parseInlineComment(s, i+j+1)
	return nil
}

// parseInlineComment extracts the "/ comment" tail starting at s[i:]. Leading
// whitespace is stripped; a single leading space after the slash is the
// common convention but not required.
func parseInlineComment(s string, i int) string {
	for i < len(s) && s[i] == ' ' {
		i++
	}
	if i >= len(s) || s[i] != '/' {
		return ""
	}
	rest := s[i+1:]
	return strings.TrimRight(strings.TrimLeft(rest, " "), " ")
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

// parseFITSFloat parses a FITS floating-point literal, accepting both the
// standard Go 'e'/'E' exponent form and the Fortran-style 'd'/'D' exponent
// form permitted by §4.2.4.
func parseFITSFloat(tok string) (float64, error) {
	t := tok
	// Replace a single d/D exponent with e.
	for i := 0; i < len(t); i++ {
		if t[i] == 'd' || t[i] == 'D' {
			t = t[:i] + "e" + t[i+1:]
			break
		}
	}
	f, err := strconv.ParseFloat(t, 64)
	if err != nil {
		return 0, err
	}
	if math.IsInf(f, 0) {
		return 0, fmt.Errorf("overflow parsing %q", tok)
	}
	return f, nil
}
