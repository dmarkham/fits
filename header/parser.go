package header

import (
	"errors"
	"fmt"
	"strings"
)

// ParseError describes a failure while parsing a header byte stream.
type ParseError struct {
	Offset int    // byte offset inside the input where the problem was detected
	Msg    string // human-readable description
}

func (e *ParseError) Error() string {
	return fmt.Sprintf("fits/header: parse error at byte %d: %s", e.Offset, e.Msg)
}

// ErrNoEnd is returned when a header byte slice does not contain an END
// keyword.
var ErrNoEnd = errors.New("fits/header: no END card found")

// ParseCards parses a header byte slice (one or more 2880-byte blocks) into
// an ordered slice of Cards. The returned slice includes every keyword,
// commentary, and COMMENT/HISTORY card but NOT the terminating END (the END
// position is returned separately as endIndex, counted in 80-byte units from
// the start of the input). If no END is found ErrNoEnd is returned. A card
// with a malformed value does not fail the parse (see DecodeCard).
//
// Long strings (the CONTINUE convention, FITS 4.0 §4.2.1.2) are joined the
// way cfitsio's fits_read_key_longstr joins them: while the string decoded
// so far ends in '&' and the next card starts with "CONTINUE" followed by
// two blanks, the '&' is removed and that card's string (or its unquoted
// value text) is appended; a CONTINUE card whose string is empty and that
// has no comment ends the chain and the '&' stays. The joined comment is the
// comments of the chain joined by one blank, except that no blank is added
// after a card whose text fills all 80 columns (cfitsio truncates the joined
// comment to 72 characters; this parser keeps all of it).
//
// The absorbed CONTINUE cards do not appear as separate Cards; Card.RawCards
// returns them, and Encode re-emits them byte for byte, so an unmodified
// header re-encodes to exactly the input. A CONTINUE card that does not
// continue a string (no '&' before it, after a non-string card, or with a
// value indicator) is kept as its own Card with Key "CONTINUE", TypeEmpty
// and Comment = columns 9–80, as cfitsio reports it.
func ParseCards(buf []byte) (cards []Card, endIndex int, err error) {
	if len(buf)%CardWidth != 0 {
		return nil, 0, &ParseError{Offset: 0, Msg: fmt.Sprintf("byte length %d not a multiple of %d", len(buf), CardWidth)}
	}
	nCards := len(buf) / CardWidth
	j := cardJoiner{cards: make([]Card, 0, nCards)}
	for i := range nCards {
		off := i * CardWidth
		raw := buf[off : off+CardWidth]
		c, decErr := DecodeCard(raw)
		if decErr != nil {
			return nil, 0, &ParseError{Offset: off, Msg: decErr.Error()}
		}
		if c.IsEnd() {
			return j.cards, i, nil
		}
		j.add(c, raw)
	}
	return nil, 0, ErrNoEnd
}

// cardJoiner collects decoded cards, joining CONTINUE cards into the long
// string before them as fits_read_key_longstr does (see ParseCards).
type cardJoiner struct {
	cards []Card
	// addDelim mirrors cfitsio's addCommDelim: whether the last physical
	// card of the current chain was shorter than 80 columns.
	addDelim bool
}

// add appends c, decoded from raw, or joins it into the previous card.
func (j *cardJoiner) add(c Card, raw []byte) {
	short := len(strings.TrimRight(string(raw), " ")) < CardWidth
	if n := len(j.cards); n > 0 && string(raw[:10]) == "CONTINUE  " {
		prev := &j.cards[n-1]
		if s, ok := prev.Value.(string); ok && prev.Type == TypeString && strings.HasSuffix(s, "&") {
			v, cm := continuation(raw)
			if v != "" || cm != "" {
				prev.Value = s[:len(s)-1] + v
				if cm != "" {
					if prev.Comment != "" && j.addDelim {
						prev.Comment += " "
					}
					prev.Comment += cm
				}
				prev.cont = append(prev.cont, c.Raw)
				j.addDelim = short
				return
			}
		}
	}
	j.cards = append(j.cards, c)
	j.addDelim = short
}

// continuation is cfitsio's ffgcnt for a card that starts "CONTINUE  ": the
// card is read as if it were "D2345678= " + columns 11–80 and its value is
// unquoted (ffc2s: a value without a leading quote is returned as is).
func continuation(raw []byte) (value, comment string) {
	t := strings.TrimRight("D2345678= "+string(raw[10:]), " ")
	v, cm, ok := splitValue(t, 10)
	if !ok {
		return "", ""
	}
	if v != "" && v[0] == '\'' {
		v = unquote(v)
	}
	return v, cm
}

// BlockCountFor returns the number of 2880-byte blocks required to hold the
// given number of cards plus a terminating END card.
func BlockCountFor(nCards int) int {
	const cardsPerBlock = 2880 / CardWidth // 36
	total := nCards + 1                    // +1 for END
	return (total + cardsPerBlock - 1) / cardsPerBlock
}
