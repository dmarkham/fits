package header

import (
	"errors"
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"
)

// commentaryWidth is the text width of a COMMENT/HISTORY/blank card
// (columns 9–80). Longer text is split over several cards, as cfitsio's
// fits_write_comment / fits_write_history do.
const commentaryWidth = CardWidth - 8

// contHead is the start of a continuation card.
const contHead = "CONTINUE  "

// EncodeCard serializes a card that occupies exactly one 80-byte card.
//
// A card read from disk (Raw non-zero) whose fields are still the ones its
// bytes decode to is returned as its original bytes; an edited copy is
// encoded from its fields (see EncodeCardWithContinuation).
// A card that needs more than one physical card (a string continued with
// CONTINUE, commentary text longer than 72 characters) is an error here;
// use EncodeCardWithContinuation.
func EncodeCard(c Card) ([CardWidth]byte, error) {
	out, err := EncodeCardWithContinuation(c)
	if err != nil {
		return [CardWidth]byte{}, err
	}
	if len(out) != 1 {
		return [CardWidth]byte{}, fmt.Errorf("fits/header: key %q needs %d cards; use EncodeCardWithContinuation", c.Key, len(out))
	}
	return out[0], nil
}

// EncodeCardWithContinuation serializes c into one or more 80-byte cards.
//
// A card read from disk (Raw non-zero) is returned as its original bytes,
// including the CONTINUE cards of a long string (Card.RawCards), as long
// as those bytes decode to the card's Key, Type, Value and Comment. A copy
// of a decoded card whose fields were changed (Raw left as it was) is
// encoded from its fields like a card built in memory, so what is written
// is what the Header getters return; if the edited card could not be
// built by Set (a keyword such as "DATE.OBS") that is an error. Fields
// are compared exactly: key letter case, float bits (so the sign of a
// zero), both parts of a complex. The one
// difference accepted between the fields and the bytes is a TypeFloat card
// over an integer literal of exactly the same value with the same sign
// (the BZERO = 9223372036854775808 the image writer builds): its Float is
// the same either way.
//
// Otherwise:
//   - commentary text (COMMENT, HISTORY, blank keyword) longer than 72
//     characters is split into consecutive 72-character cards, as cfitsio
//     does; readers return each card's text without trailing blanks;
//   - a string that does not fit on one card with its comment uses the
//     CONTINUE convention (§4.2.1.2): every card but the last ends its
//     piece with '&', an embedded quote is never split from its escape, and
//     a comment that does not fit after the last piece goes on following
//     "CONTINUE  '&' / ..." cards, split only at single blanks (a comment
//     that cannot be split that way is an error); such a string must not
//     end in '&' (an error: cfitsio keeps that '&', astropy drops it);
//   - any other value must fit, with its comment, on one card, or an error
//     is returned.
//
// Trailing blanks of a string value are not significant in FITS
// (§4.2.1.1) and readers (cfitsio, astropy, this package) drop them, so
// they do not survive a write; leading and embedded blanks do. The same
// holds for the comment of a value card (written without its trailing
// blanks) and for commentary text (readers return each card's text
// without them). A string
// that needs CONTINUE cards is written without its trailing blanks, so the
// chain never ends in a blank-only piece that readers take for an empty
// CONTINUE.
//
// Every keyword, value and comment is validated first (see Set); invalid
// input is an error, never a truncated or altered card. As Set does, and
// as cfitsio's fits_write_key does, the keyword is written upper-cased.
// A card whose keyword upper-cases to END is written only as the bare END
// card, from Card{Key: "END"}: one with a comment, or with a lower-case
// key such as a decoded "end" card (which cfitsio reads as an ordinary
// card, not as the end of the header), is an error.
func EncodeCardWithContinuation(c Card) ([][CardWidth]byte, error) {
	out, _, err := encodeCard(c)
	return out, err
}

// encodeCard is EncodeCardWithContinuation; fromRaw reports whether the
// card's original bytes were used.
func encodeCard(c Card) (out [][CardWidth]byte, fromRaw bool, err error) {
	if c.Raw != ([CardWidth]byte{}) {
		if rawMatches(&c) {
			return c.RawCards(), true, nil
		}
	}
	out, err = encodeFields(c)
	return out, false, err
}

// encodeFields encodes a card from its Key, Type, Value and Comment.
func encodeFields(c Card) ([][CardWidth]byte, error) {
	// Keyword names are written upper-cased, as Set and AddCard store them
	// and as cfitsio's fits_write_key writes them.
	key := canonKey(c.Key)
	if key == KeyEnd && (c.Key != KeyEnd || c.Comment != "") {
		// Upper-casing would turn an ordinary card ("end", which cfitsio
		// reads as a card like any other) into the END card, and the END
		// card has no room for a comment: either way something is lost.
		return nil, fmt.Errorf("%w: %q with comment %q would be written as a bare END card, which ends the header", ErrInvalidKeyword, c.Key, c.Comment)
	}
	c.Key = key
	if err := validateCard(&c); err != nil {
		return nil, err
	}
	if c.Key == KeyEnd {
		return [][CardWidth]byte{newCard("END")}, nil
	}
	if c.IsCommentary() {
		return encodeCommentary(c.Key, c.Comment), nil
	}
	comment := strings.TrimRight(c.Comment, " ")
	head := cardHead(c.Key)
	if c.Type == TypeString {
		return encodeString(c.Key, head, c.Value.(string), comment)
	}
	val, err := formatValue(c)
	if err != nil {
		return nil, err
	}
	if c.Type == TypeEmpty && (strings.HasPrefix(comment, " ") || strings.HasPrefix(comment, "/")) {
		// After the '/' of an undefined value readers skip every blank and
		// a second '/', so such a comment would not read back.
		return nil, fmt.Errorf("fits/header: key %q: the comment of an undefined value cannot start with a blank or '/'", c.Key)
	}
	line := head
	if !c.IsHierarch() && len(head)+len(val) < 30 {
		// Fixed format: numbers and logicals right-justified to column 30.
		line += strings.Repeat(" ", 30-len(head)-len(val))
	}
	line += val
	if comment != "" {
		line += " / " + comment
	}
	if len(line) > CardWidth {
		return nil, fmt.Errorf("fits/header: key %q: value and comment need %d columns, a card has %d (shorten the comment)", c.Key, len(line), CardWidth)
	}
	return [][CardWidth]byte{newCard(line)}, nil
}

// cardHead is the keyword and value indicator: "KEYWORD = " or
// "HIERARCH path = ".
func cardHead(key string) string {
	if strings.HasPrefix(key, "HIERARCH ") {
		return key + " = "
	}
	return padKey(key) + "= "
}

func newCard(s string) [CardWidth]byte {
	var out [CardWidth]byte
	n := copy(out[:], s)
	for i := n; i < CardWidth; i++ {
		out[i] = ' '
	}
	return out
}

func encodeCommentary(key, text string) [][CardWidth]byte {
	k := padKey(key)
	var out [][CardWidth]byte
	for {
		n := min(len(text), commentaryWidth)
		out = append(out, newCard(k+text[:n]))
		text = text[n:]
		if text == "" {
			return out
		}
	}
}

// encodeString writes a string value on one card if it fits with its
// comment, else as a CONTINUE chain.
//
// A string that needs CONTINUE cards loses its trailing blanks before it
// is split: they are not significant (FITS §4.2.1.1), and a last piece of
// only blanks with no comment would be written as CONTINUE '   ', which
// cfitsio and this package read as an empty CONTINUE that ends the chain
// and leaves the '&' marker in the value. A string that needs CONTINUE
// cards and ends in '&' is refused: cfitsio keeps that '&', astropy takes
// it for a continuation marker and drops it.
func encodeString(key, head, s, comment string) ([][CardWidth]byte, error) {
	esc := fitsEscape(s)
	single := head + "'" + esc + "'"
	if comment != "" {
		single += " / " + comment
	}
	if len(single) <= CardWidth {
		return [][CardWidth]byte{newCard(single)}, nil
	}
	if t := strings.TrimRight(s, " "); t != s {
		return encodeString(key, head, t, comment)
	}
	if strings.HasSuffix(s, "&") {
		return nil, fmt.Errorf("fits/header: key %q: a string that needs CONTINUE cards cannot end in '&' (readers disagree on whether it is a continuation marker)", key)
	}
	firstCap := CardWidth - len(head) - 3 // quotes and '&'
	if firstCap < 0 {
		return nil, fmt.Errorf("fits/header: HIERARCH key %q leaves no room for a string value", key)
	}
	const contCap = CardWidth - len(contHead) - 3
	pieces := splitEscaped(s, firstCap, contCap)

	var out [][CardWidth]byte
	for i, p := range pieces {
		h := contHead
		if i == 0 {
			h = head
		}
		last := i == len(pieces)-1
		body := h + "'" + p
		if !last {
			body += "&'"
			out = append(out, newCard(body))
			continue
		}
		// Last piece: the comment goes here if it fits whole.
		if comment == "" || len(body)+len("' / ")+len(comment) <= CardWidth {
			body += "'"
			if comment != "" {
				body += " / " + comment
			}
			out = append(out, newCard(body))
			return out, nil
		}
		out = append(out, newCard(body+"&'"))
	}
	// Comment cards. Each one is kept below 80 columns so cfitsio rejoins
	// the pieces with one blank, as astropy does.
	const cmtHead = contHead + "'&' / "
	const cmtCap = CardWidth - len(cmtHead) - 1
	words, err := splitComment(comment, cmtCap)
	if err != nil {
		return nil, fmt.Errorf("fits/header: key %q: %w", key, err)
	}
	for i, w := range words {
		if i == len(words)-1 {
			out = append(out, newCard(contHead+"'' / "+w))
		} else {
			out = append(out, newCard(cmtHead+w))
		}
	}
	return out, nil
}

// splitEscaped cuts s into pieces whose escaped length is at most firstCap
// for the first piece and contCap for the rest, never splitting a quote
// from its escape. The last piece may be empty only when s is empty.
func splitEscaped(s string, firstCap, contCap int) []string {
	var out []string
	capacity := firstCap
	var cur strings.Builder
	n := 0
	for i := 0; i < len(s); i++ {
		w := 1
		if s[i] == '\'' {
			w = 2
		}
		if n+w > capacity {
			out = append(out, cur.String())
			cur.Reset()
			n = 0
			capacity = contCap
		}
		if s[i] == '\'' {
			cur.WriteString("''")
		} else {
			cur.WriteByte(s[i])
		}
		n += w
	}
	return append(out, cur.String())
}

// splitComment splits a comment into pieces of at most capacity bytes,
// breaking only at a single blank (which is dropped and restored by
// readers as the joining blank).
func splitComment(c string, capacity int) ([]string, error) {
	var out []string
	for len(c) > capacity {
		cut := -1
		for k := min(capacity, len(c)-2); k > 0; k-- {
			if c[k] == ' ' && c[k-1] != ' ' && c[k+1] != ' ' {
				cut = k
				break
			}
		}
		if cut < 0 {
			return nil, fmt.Errorf("comment %q cannot be continued: no single blank to break it at within %d characters", c, capacity)
		}
		out = append(out, c[:cut])
		c = c[cut+1:]
	}
	return append(out, c), nil
}

// padKey returns name padded with spaces to 8 characters. Callers have
// validated that name is at most 8 bytes.
func padKey(name string) string {
	return name + strings.Repeat(" ", 8-len(name))
}

// fitsEscape returns s with embedded single quotes doubled, as required for
// a FITS string literal.
func fitsEscape(s string) string {
	if !strings.ContainsRune(s, '\'') {
		return s
	}
	return strings.ReplaceAll(s, "'", "''")
}

// formatFloat writes f so that it reads back bit-exactly and is typed as a
// float by FITS readers: the shortest round-trip form, with ".0" added when
// it would otherwise have neither '.' nor an exponent.
func formatFloat(f float64) (string, error) {
	if math.IsNaN(f) || math.IsInf(f, 0) {
		return "", fmt.Errorf("fits/header: %v cannot be written as a FITS number", f)
	}
	s := strconv.FormatFloat(f, 'G', -1, 64)
	if !strings.ContainsAny(s, ".E") {
		s += ".0"
	}
	return s, nil
}

// formatValue formats a non-commentary, non-string card's value — the text
// between the value indicator and the optional "/ comment".
func formatValue(c Card) (string, error) {
	switch c.Type {
	case TypeString:
		s, _ := c.Value.(string)
		return "'" + fitsEscape(s) + "'", nil
	case TypeLogical:
		if c.Value.(bool) {
			return "T", nil
		}
		return "F", nil
	case TypeInt:
		v := c.Value.(int64)
		if v == 0 && c.negZero {
			// A decoded "-0" card (Float reads -0, as cfitsio does)
			// whose comment was edited still holds -0.
			return "-0", nil
		}
		return strconv.FormatInt(v, 10), nil
	case TypeBigInt:
		return c.Value.(*big.Int).String(), nil
	case TypeFloat:
		return formatFloat(c.Value.(float64))
	case TypeComplexInt:
		v := c.Value.(Complex)
		if !integerComplex(v) {
			return "", fmt.Errorf("fits/header: integer complex %v has non-integer parts", v)
		}
		return fmt.Sprintf("(%d, %d)", int64(v.Re), int64(v.Im)), nil
	case TypeComplexFloat:
		v := c.Value.(Complex)
		re, err := formatFloat(v.Re)
		if err != nil {
			return "", err
		}
		im, err := formatFloat(v.Im)
		if err != nil {
			return "", err
		}
		return "(" + re + ", " + im + ")", nil
	case TypeUnparsed:
		return c.Value.(RawValue).Text, nil
	case TypeEmpty:
		return "", nil
	}
	return "", fmt.Errorf("fits/header: unknown value type %v", c.Type)
}

// integerComplex reports whether both parts of v are integers that a
// TypeComplexInt card can hold exactly.
func integerComplex(v Complex) bool {
	return v.Re == math.Trunc(v.Re) && v.Im == math.Trunc(v.Im) && math.Abs(v.Re) <= 1<<53 && math.Abs(v.Im) <= 1<<53
}

// ErrInvalidKeyword is returned for a keyword name that cannot be written.
var ErrInvalidKeyword = errors.New("fits/header: invalid keyword")

// ErrInvalidText is returned for a value or comment that is not printable
// ASCII (FITS 4.0 §4.1.2: header text is restricted to 0x20–0x7E).
var ErrInvalidText = errors.New("fits/header: text is not printable ASCII")

// validateKey checks a keyword for writing. Legal names: "" (blank
// commentary); up to 8 characters from A–Z, 0–9, '-' and '_' (the set
// cfitsio's fits_test_keyword accepts); or "HIERARCH " followed by a
// non-empty path of printable ASCII without '=', its words separated by
// single blanks.
func validateKey(key string) error {
	if p, ok := strings.CutPrefix(key, "HIERARCH "); ok {
		if strings.TrimSpace(p) == "" {
			return fmt.Errorf("%w: %q: empty HIERARCH keyword path", ErrInvalidKeyword, key)
		}
		if norm := strings.Join(strings.Fields(p), " "); norm != p {
			// The decoder (like the Card doc) keeps single blanks between
			// path words, so any other spacing would change on a round trip.
			return fmt.Errorf("%w: %q: HIERARCH path words must be separated by single blanks, with none leading or trailing (write %q)", ErrInvalidKeyword, key, "HIERARCH "+norm)
		}
		for i := 0; i < len(p); i++ {
			if p[i] < 0x20 || p[i] > 0x7E || p[i] == '=' {
				return fmt.Errorf("%w: %q: HIERARCH path byte %q is not allowed", ErrInvalidKeyword, key, p[i])
			}
		}
		if len(key)+len(" = ") >= CardWidth {
			return fmt.Errorf("%w: %q: HIERARCH keyword leaves no room for a value", ErrInvalidKeyword, key)
		}
		return nil
	}
	if len(key) > 8 {
		return fmt.Errorf("%w: %q has %d characters; FITS keywords have at most 8 (use \"HIERARCH %s\" for a long name)", ErrInvalidKeyword, key, len(key), key)
	}
	for i := 0; i < len(key); i++ {
		b := key[i]
		if !(b >= 'A' && b <= 'Z' || b >= '0' && b <= '9' || b == '-' || b == '_') {
			return fmt.Errorf("%w: %q: character %q is not allowed (A-Z, 0-9, '-', '_')", ErrInvalidKeyword, key, b)
		}
	}
	return nil
}

func validateText(what, s string) error {
	for i := 0; i < len(s); i++ {
		if s[i] < 0x20 || s[i] > 0x7E {
			return fmt.Errorf("%w: %s byte %d is 0x%02X", ErrInvalidText, what, i, s[i])
		}
	}
	return nil
}

// validateRawText checks that a RawValue text, written verbatim after
// "KEY = ", reads back as the same unparsed text: one token as cfitsio's
// ffpsvc splits it (so no blank or '/' outside parentheses, no leading
// quote, at most 70 characters) that does not type as a FITS value.
func validateRawText(key, text string) error {
	bad := func(why string) error {
		return fmt.Errorf("fits/header: key %q: RawValue %q %s, so it would not read back as written", key, text, why)
	}
	if text == "" {
		return bad("is empty (an undefined value)")
	}
	v, _, ok := splitValue(text+" / c", 0)
	if !ok {
		return bad("is a '(' group without ')'")
	}
	if v != text {
		return bad(fmt.Sprintf("is read as the value %q followed by more text", v))
	}
	var d Card
	typeValue(v, &d)
	if d.Type != TypeUnparsed {
		return bad("is a valid FITS " + d.Type.String() + " value (set it as that type)")
	}
	return nil
}

// validateCard checks an in-memory card before encoding.
func validateCard(c *Card) error {
	if err := validateKey(c.Key); err != nil {
		return err
	}
	if err := validateText("comment of "+c.Key, c.Comment); err != nil {
		return err
	}
	if c.Key == KeyEnd && c.Type != TypeEmpty {
		return fmt.Errorf("%w: END cannot carry a value", ErrInvalidKeyword)
	}
	if c.Key == "CONTINUE" {
		// A hand-made CONTINUE card would be joined on read into a
		// preceding string that ends in '&'.
		return fmt.Errorf("%w: CONTINUE cards are written by the encoder for long strings; set the string on its own keyword", ErrInvalidKeyword)
	}
	if (c.IsCommentary() || c.Key == KeyEnd) && c.Type != TypeEmpty {
		return fmt.Errorf("fits/header: commentary key %q cannot carry a typed value", c.Key)
	}
	switch c.Type {
	case TypeString:
		s, ok := c.Value.(string)
		if !ok {
			return fmt.Errorf("fits/header: key %q: TypeString with value %T", c.Key, c.Value)
		}
		return validateText("string value of "+c.Key, s)
	case TypeUnparsed:
		rv, ok := c.Value.(RawValue)
		if !ok {
			return fmt.Errorf("fits/header: key %q: TypeUnparsed with value %T", c.Key, c.Value)
		}
		if err := validateText("value of "+c.Key, rv.Text); err != nil {
			return err
		}
		return validateRawText(c.Key, rv.Text)
	case TypeBigInt:
		b, ok := c.Value.(*big.Int)
		if !ok {
			return fmt.Errorf("fits/header: key %q: TypeBigInt with value %T", c.Key, c.Value)
		}
		if b == nil {
			return fmt.Errorf("fits/header: key %q: nil *big.Int", c.Key)
		}
		if b.IsInt64() {
			// Int would report ErrIntOverflow for a value every reader
			// takes for an ordinary integer.
			return fmt.Errorf("fits/header: key %q: TypeBigInt value %s is inside int64; use TypeInt", c.Key, b)
		}
	case TypeInt:
		if _, ok := c.Value.(int64); !ok {
			return fmt.Errorf("fits/header: key %q: TypeInt with value %T", c.Key, c.Value)
		}
	case TypeFloat:
		f, ok := c.Value.(float64)
		if !ok {
			return fmt.Errorf("fits/header: key %q: TypeFloat with value %T", c.Key, c.Value)
		}
		if math.IsNaN(f) || math.IsInf(f, 0) {
			return fmt.Errorf("fits/header: key %q: %v cannot be written as a FITS number", c.Key, f)
		}
	case TypeLogical:
		if _, ok := c.Value.(bool); !ok {
			return fmt.Errorf("fits/header: key %q: TypeLogical with value %T", c.Key, c.Value)
		}
	case TypeComplexInt, TypeComplexFloat:
		v, ok := c.Value.(Complex)
		if !ok {
			return fmt.Errorf("fits/header: key %q: complex type with value %T", c.Key, c.Value)
		}
		for _, f := range []float64{v.Re, v.Im} {
			if math.IsNaN(f) || math.IsInf(f, 0) {
				return fmt.Errorf("fits/header: key %q: complex part %v cannot be written as a FITS number", c.Key, f)
			}
		}
		if c.Type == TypeComplexInt && !integerComplex(v) {
			return fmt.Errorf("fits/header: key %q: integer complex %v has non-integer parts (use TypeComplexFloat)", c.Key, v)
		}
	case TypeEmpty:
		if c.Value != nil {
			return fmt.Errorf("fits/header: key %q: TypeEmpty with value %T", c.Key, c.Value)
		}
	default:
		return fmt.Errorf("fits/header: key %q: unknown value type %v", c.Key, c.Type)
	}
	return nil
}

// Encode serializes a full Header into one or more 2880-byte blocks,
// terminated by an END card and padded with blank cards to the next block
// boundary. Cards read from disk and not replaced are emitted byte for
// byte, so an unmodified parsed header encodes to the bytes it was parsed
// from (when, as FITS requires, the space after END is blank).
func Encode(h *Header) ([]byte, error) {
	var cards [][CardWidth]byte
	// owner[j] is the index in h.cards of the card that read-back card j
	// comes from: a long commentary text built in memory reads back as
	// one card per 72 characters, every other card as one card.
	var owner []int
	for i := range h.cards {
		enc, fromRaw, err := encodeCard(h.cards[i])
		if err != nil {
			return nil, fmt.Errorf("encoding card %d (%s): %w", i, h.cards[i].Key, err)
		}
		cards = append(cards, enc...)
		// Commentary text written from memory is n cards of 72 characters,
		// read back as n cards. Bytes from a file read back as one card,
		// even under a blank key: a malformed card such as
		// " 0000= '...&'" has Key "" and can carry CONTINUE cards.
		n := 1
		if !fromRaw && h.cards[i].IsCommentary() {
			n = len(enc)
		}
		for range n {
			owner = append(owner, i)
		}
	}
	cards = append(cards, newCard("END"))

	const cardsPerBlock = 2880 / CardWidth
	for len(cards)%cardsPerBlock != 0 {
		cards = append(cards, newCard(""))
	}
	if err := checkCardBoundaries(h, cards, owner); err != nil {
		return nil, err
	}
	out := make([]byte, 0, len(cards)*CardWidth)
	for i := range cards {
		out = append(out, cards[i][:]...)
	}
	return out, nil
}

// checkCardBoundaries reads the encoded cards back as ParseCards would
// (decoding each one, joining CONTINUE cards, stopping at the first END)
// and requires them to be the cards they were built from, one for one, in
// order and under the same keys. It catches cards that are fine one by
// one but not side by side: a one-card string ending in '&' followed by a
// CONTINUE card (taken from another header) would be read back as one
// joined string. The error names the first card that would not read back
// as itself: for such a join, the CONTINUE card that is swallowed.
func checkCardBoundaries(h *Header, cards [][CardWidth]byte, owner []int) error {
	var j cardJoiner
	for i := range cards {
		d := decodeCard(&cards[i])
		if d.IsEnd() {
			break
		}
		j.add(d, cards[i][:])
	}
	back := j.cards
	for k := range max(len(back), len(owner)) {
		if k >= len(back) || k >= len(owner) || canonKey(back[k].Key) != canonKey(h.cards[owner[k]].Key) {
			i := owner[min(k, len(owner)-1)]
			got := "nothing"
			if k < len(back) {
				got = fmt.Sprintf("%q", back[k].Key)
			}
			return fmt.Errorf("fits/header: card %d (%s) would read back as %s (the encoded header reads back as %d cards, not %d): the cards around it join or split on read", i, h.cards[i].Key, got, len(back), len(owner))
		}
	}
	return nil
}

// rawMatches reports whether c's original bytes (Raw and the CONTINUE
// cards of a long string) decode, on their own, to exactly one card with
// c's Key, Comment and value, so that writing them writes c.
func rawMatches(c *Card) bool {
	var j cardJoiner
	for _, r := range c.RawCards() {
		j.add(decodeCard(&r), r[:])
	}
	if len(j.cards) != 1 {
		return false
	}
	d := &j.cards[0]
	return d.Key == c.Key && d.Comment == c.Comment && sameDecodedValue(c, d)
}

// sameDecodedValue reports whether c holds the value d was decoded with.
// Values of the same type must be identical (floats bit for bit). The one
// cross-type match is a TypeFloat c over an integer d of exactly the same
// value, whose Float (as Header.Float reads d, -0 included) is the same
// float64.
func sameDecodedValue(c, d *Card) bool {
	if c.Type != d.Type {
		f, ok := c.Value.(float64)
		if c.Type != TypeFloat || !ok {
			return false
		}
		var exact *big.Int
		switch d.Type {
		case TypeInt:
			exact = big.NewInt(d.Value.(int64))
		case TypeBigInt:
			exact = d.Value.(*big.Int)
		default:
			return false
		}
		df, acc := new(big.Float).SetInt(exact).Float64()
		if d.negZero {
			df = math.Copysign(0, -1)
		}
		return acc == big.Exact && math.Float64bits(df) == math.Float64bits(f)
	}
	switch c.Type {
	case TypeEmpty:
		return c.Value == nil
	case TypeString:
		v, ok := c.Value.(string)
		return ok && v == d.Value.(string)
	case TypeLogical:
		v, ok := c.Value.(bool)
		return ok && v == d.Value.(bool)
	case TypeInt:
		v, ok := c.Value.(int64)
		return ok && v == d.Value.(int64) && c.negZero == d.negZero
	case TypeBigInt:
		v, ok := c.Value.(*big.Int)
		return ok && v != nil && v.Cmp(d.Value.(*big.Int)) == 0
	case TypeFloat:
		v, ok := c.Value.(float64)
		return ok && math.Float64bits(v) == math.Float64bits(d.Value.(float64))
	case TypeComplexInt, TypeComplexFloat:
		v, ok := c.Value.(Complex)
		w := d.Value.(Complex)
		return ok && math.Float64bits(v.Re) == math.Float64bits(w.Re) && math.Float64bits(v.Im) == math.Float64bits(w.Im)
	case TypeUnparsed:
		v, ok := c.Value.(RawValue)
		return ok && v.Text == d.Value.(RawValue).Text
	}
	return false
}
