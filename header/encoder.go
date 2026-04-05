package header

import (
	"fmt"
	"strconv"
	"strings"
)

// EncodeCard serializes a single Card into its 80-byte on-disk form.
//
// If the card's Raw field is non-zero (i.e. it was parsed from disk and
// never mutated through Set/Add), the raw bytes are returned as-is to
// guarantee byte-for-byte round-trip fidelity for unchanged cards.
//
// For strings that exceed the 68-character fixed-format value field, this
// function returns an error — callers that need CONTINUE splitting must use
// EncodeCardWithContinuation, which returns one or more 80-byte cards.
func EncodeCard(c Card) ([CardWidth]byte, error) {
	if c.Raw != ([CardWidth]byte{}) {
		return c.Raw, nil
	}
	buf, err := encodeCardOnce(c)
	if err != nil {
		return [CardWidth]byte{}, err
	}
	return buf, nil
}

// EncodeCardWithContinuation serializes c into one or more 80-byte cards,
// applying the CONTINUE convention (§4.2.1.2) if the string value is too
// long for a single card.
//
// The caller's comment is placed on the LAST emitted card per convention
// (cfitsio ffmkky); leading cards carry an empty comment.
func EncodeCardWithContinuation(c Card) ([][CardWidth]byte, error) {
	if c.Raw != ([CardWidth]byte{}) {
		return [][CardWidth]byte{c.Raw}, nil
	}
	// Strings that fit in a single card go through the normal path.
	if c.Type != TypeString {
		b, err := encodeCardOnce(c)
		if err != nil {
			return nil, err
		}
		return [][CardWidth]byte{b}, nil
	}
	s := c.Value.(string)
	// Measure whether the value fits on one card. The fixed-format string
	// value field occupies columns 11..80 (70 chars), of which two are the
	// enclosing quotes. Each embedded quote character doubles, so we account
	// for escaped length.
	escLen := escapedLen(s)
	// Available payload on a primary card: column 11 to (up to) column 80,
	// leaving room for the comment. The string payload (inside quotes) can
	// use up to 68 characters in the fully-filled case (columns 12..79,
	// reserving cols 80 for space). Comment, if any, needs at least 3 more
	// chars (" / "). If both fit we stay single-card.
	const maxSinglePayload = 68
	const maxSingleWithComment = 60 // leaves "= '...' / comment " room
	fitsOne := escLen <= maxSinglePayload && (c.Comment == "" || escLen <= maxSingleWithComment)
	if fitsOne {
		b, err := encodeCardOnce(c)
		if err != nil {
			return nil, err
		}
		return [][CardWidth]byte{b}, nil
	}
	return encodeCONTINUE(c.Key, s, c.Comment)
}

// encodeCardOnce writes a single 80-byte card for values that fit in one
// card. It does NOT apply CONTINUE; callers must pre-check.
func encodeCardOnce(c Card) ([CardWidth]byte, error) {
	var out [CardWidth]byte
	for i := range out {
		out[i] = ' '
	}

	if c.Key == KeyEnd {
		copy(out[0:], "END")
		return out, nil
	}

	// Commentary card (COMMENT, HISTORY, blank): free text starts at col 9.
	if c.IsCommentary() {
		if c.Key != "" {
			copy(out[0:8], padKey(c.Key))
		}
		txt := c.Comment
		if len(txt) > CardWidth-8 {
			txt = txt[:CardWidth-8]
		}
		copy(out[8:], txt)
		return out, nil
	}

	// HIERARCH: "HIERARCH path = value / comment".
	if c.IsHierarch() {
		return encodeHierarchCard(c)
	}

	// Normal keyed card: "KEYWORD = value / comment".
	copy(out[0:8], padKey(c.Key))
	out[8] = '='
	out[9] = ' '

	valStr, err := formatValue(c)
	if err != nil {
		return out, err
	}
	// Fixed-format placement puts numeric values right-justified in cols
	// 11..30 and strings left-justified starting at column 11 with the
	// opening quote at column 11. We follow cfitsio here: strings begin at
	// col 11, numerics right-justify to column 30 when they fit.
	pos := 10 // column 11 (0-based index 10)
	if c.Type == TypeInt || c.Type == TypeFloat || c.Type == TypeLogical {
		// Right-justify to column 30 (index 29).
		right := 30
		if right-pos >= len(valStr) {
			pos = right - len(valStr)
		}
	}
	if pos+len(valStr) > CardWidth {
		return out, fmt.Errorf("fits/header: value field too long for key %q", c.Key)
	}
	copy(out[pos:], valStr)
	end := pos + len(valStr)

	// Inline comment.
	if c.Comment != "" {
		// " / comment" starting at end+1 (leave one space).
		if end+3 >= CardWidth {
			// No room for comment — silently drop. Round-trip via Raw.
			return out, nil
		}
		out[end+1] = '/'
		out[end+2] = ' '
		cmt := c.Comment
		avail := CardWidth - (end + 3)
		if len(cmt) > avail {
			cmt = cmt[:avail]
		}
		copy(out[end+3:], cmt)
	}
	return out, nil
}

// encodeHierarchCard serializes a HIERARCH card: "HIERARCH path = value / comment".
// The value is placed free-form (no column 30 right-justification).
func encodeHierarchCard(c Card) ([CardWidth]byte, error) {
	var out [CardWidth]byte
	for i := range out {
		out[i] = ' '
	}
	// Key is "HIERARCH path" literally.
	if !strings.HasPrefix(c.Key, "HIERARCH ") {
		return out, fmt.Errorf("fits/header: HIERARCH card must begin with \"HIERARCH \"")
	}
	head := c.Key + " = "
	if len(head) > CardWidth-2 {
		return out, fmt.Errorf("fits/header: HIERARCH key too long: %q", c.Key)
	}
	copy(out[0:], head)
	pos := len(head)
	valStr, err := formatValue(c)
	if err != nil {
		return out, err
	}
	if pos+len(valStr) > CardWidth {
		return out, fmt.Errorf("fits/header: HIERARCH value too long: %q", c.Key)
	}
	copy(out[pos:], valStr)
	end := pos + len(valStr)
	if c.Comment != "" && end+3 < CardWidth {
		out[end+1] = '/'
		out[end+2] = ' '
		avail := CardWidth - (end + 3)
		cmt := c.Comment
		if len(cmt) > avail {
			cmt = cmt[:avail]
		}
		copy(out[end+3:], cmt)
	}
	return out, nil
}

// padKey returns name padded with spaces to exactly 8 characters. name is
// assumed to be ≤ 8 bytes; longer names must go through the HIERARCH path.
func padKey(name string) string {
	if len(name) > 8 {
		return name[:8]
	}
	return name + strings.Repeat(" ", 8-len(name))
}

// escapedLen returns the length of s after FITS string-escaping (embedded
// single quotes doubled).
func escapedLen(s string) int {
	n := len(s)
	for i := 0; i < len(s); i++ {
		if s[i] == '\'' {
			n++
		}
	}
	return n
}

// fitsEscape returns s with embedded single quotes doubled, as required for
// a FITS string literal.
func fitsEscape(s string) string {
	if !strings.ContainsRune(s, '\'') {
		return s
	}
	return strings.ReplaceAll(s, "'", "''")
}

// formatValue formats a non-commentary card's typed value to its textual
// form — the string between the "= " and the optional "/ comment".
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
		return strconv.FormatInt(c.Value.(int64), 10), nil
	case TypeFloat:
		f := c.Value.(float64)
		// Use a representation that round-trips via ParseFloat.
		return strconv.FormatFloat(f, 'G', -1, 64), nil
	case TypeComplexInt:
		v := c.Value.(Complex)
		return fmt.Sprintf("(%d, %d)", int64(v.Re), int64(v.Im)), nil
	case TypeComplexFloat:
		v := c.Value.(Complex)
		return fmt.Sprintf("(%s, %s)",
			strconv.FormatFloat(v.Re, 'G', -1, 64),
			strconv.FormatFloat(v.Im, 'G', -1, 64)), nil
	case TypeEmpty:
		return "", nil
	}
	return "", fmt.Errorf("fits/header: unknown value type %v", c.Type)
}

// encodeCONTINUE implements §4.2.1.2 CONTINUE long-string splitting. The
// first card has the real keyword; every following card has keyword
// "CONTINUE" and is a pure commentary-style extension.
//
// Chunk sizing: the first card must fit "KEYWORD = 'head&'" — columns 11..80
// hold the quoted payload, 70 chars, minus 2 quotes, minus 1 for '&' = 67
// chars of user payload on the first card (no room for a leading/trailing
// comment on the first card). Continuation cards use "CONTINUE  'piece&'"
// where "CONTINUE" occupies cols 1..8, then two spaces, then the string —
// same 67-char capacity. The last card omits the '&' and may include a
// trailing "/ comment" if room allows.
func encodeCONTINUE(key, value, comment string) ([][CardWidth]byte, error) {
	const chunkCap = 67 // max chars inside the quotes per card (not counting &)
	var chunks []string
	rest := value
	for len(rest) > chunkCap {
		chunks = append(chunks, rest[:chunkCap])
		rest = rest[chunkCap:]
	}
	// Last chunk: the remainder (may be empty if the value landed exactly
	// on a boundary, which is a legitimate edge case handled here).
	chunks = append(chunks, rest)

	out := make([][CardWidth]byte, 0, len(chunks))
	for i, ch := range chunks {
		last := i == len(chunks)-1
		var card [CardWidth]byte
		for j := range card {
			card[j] = ' '
		}
		if i == 0 {
			copy(card[0:8], padKey(key))
			card[8] = '='
			card[9] = ' '
		} else {
			copy(card[0:8], "CONTINUE")
			// Two spaces at columns 9..10 per the convention — CONTINUE
			// cards have NO "= " indicator, only the quoted string.
			card[8] = ' '
			card[9] = ' '
		}
		text := fitsEscape(ch)
		if !last {
			text += "&"
		}
		// Quote the text.
		card[10] = '\''
		write := 11
		if write+len(text)+1 > CardWidth {
			return nil, fmt.Errorf("fits/header: CONTINUE chunk overflow for key %q", key)
		}
		copy(card[write:], text)
		write += len(text)
		card[write] = '\''
		write++
		// Comment only on the last card, space permitting.
		if last && comment != "" && write+3 < CardWidth {
			card[write+1] = '/'
			card[write+2] = ' '
			avail := CardWidth - (write + 3)
			cmt := comment
			if len(cmt) > avail {
				cmt = cmt[:avail]
			}
			copy(card[write+3:], cmt)
		}
		out = append(out, card)
	}
	return out, nil
}

// Encode serializes a full Header into one or more 2880-byte blocks,
// terminated by an END card and padded with spaces to the next block
// boundary.
func Encode(h *Header) ([]byte, error) {
	// Pre-compute: walk cards, emit one or more physical 80-byte cards per
	// logical card (CONTINUE for long strings), then append END and pad.
	var cards [][CardWidth]byte
	for i := range h.cards {
		enc, err := EncodeCardWithContinuation(h.cards[i])
		if err != nil {
			return nil, fmt.Errorf("encoding card %d (%s): %w", i, h.cards[i].Key, err)
		}
		cards = append(cards, enc...)
	}
	// END card.
	var endCard [CardWidth]byte
	for i := range endCard {
		endCard[i] = ' '
	}
	copy(endCard[0:], "END")
	cards = append(cards, endCard)

	// Pad up to block boundary with all-space cards.
	const cardsPerBlock = 2880 / CardWidth
	rem := len(cards) % cardsPerBlock
	if rem != 0 {
		pad := cardsPerBlock - rem
		blank := [CardWidth]byte{}
		for i := range blank {
			blank[i] = ' '
		}
		for range pad {
			cards = append(cards, blank)
		}
	}

	out := make([]byte, 0, len(cards)*CardWidth)
	for i := range cards {
		out = append(out, cards[i][:]...)
	}
	return out, nil
}
