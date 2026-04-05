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
// the start of the input). If no END is found ErrNoEnd is returned.
//
// ParseCards applies the CONTINUE long-string convention (§4.2.1.2): when a
// string value ends in "&" and is immediately followed by a CONTINUE card,
// the string is joined and the pair emitted as one logical Card. The joined
// card retains the original bytes of the first card in Card.Raw; subsequent
// CONTINUE cards are consumed and do not appear separately in the output.
// Callers that need to round-trip continuation intact must re-serialize via
// Encode, which re-emits CONTINUE cards as needed.
func ParseCards(buf []byte) (cards []Card, endIndex int, err error) {
	if len(buf)%CardWidth != 0 {
		return nil, 0, &ParseError{Offset: 0, Msg: fmt.Sprintf("byte length %d not a multiple of %d", len(buf), CardWidth)}
	}
	nCards := len(buf) / CardWidth
	cards = make([]Card, 0, nCards)
	for i := range nCards {
		off := i * CardWidth
		raw := buf[off : off+CardWidth]
		c, decErr := DecodeCard(raw)
		if decErr != nil {
			return nil, 0, &ParseError{Offset: off, Msg: decErr.Error()}
		}
		if c.IsEnd() {
			return cards, i, nil
		}
		// CONTINUE continuation of the previous string value.
		if c.Key == "CONTINUE" {
			if len(cards) == 0 {
				return nil, 0, &ParseError{Offset: off, Msg: "CONTINUE with no preceding card"}
			}
			prev := &cards[len(cards)-1]
			if prev.Type != TypeString {
				return nil, 0, &ParseError{Offset: off, Msg: "CONTINUE following non-string card"}
			}
			prevStr, _ := prev.Value.(string)
			// Only join if the previous string ended with '&'.
			if !strings.HasSuffix(prevStr, "&") {
				return nil, 0, &ParseError{Offset: off, Msg: "CONTINUE without '&' marker on previous card"}
			}
			contStr, _ := c.Value.(string)
			prev.Value = strings.TrimSuffix(prevStr, "&") + contStr
			if c.Comment != "" {
				if prev.Comment == "" {
					prev.Comment = c.Comment
				} else {
					prev.Comment = prev.Comment + " " + c.Comment
				}
			}
			continue
		}
		cards = append(cards, c)
	}
	return nil, 0, ErrNoEnd
}

// BlockCountFor returns the number of 2880-byte blocks required to hold the
// given number of cards plus a terminating END card.
func BlockCountFor(nCards int) int {
	const cardsPerBlock = 2880 / CardWidth // 36
	total := nCards + 1                    // +1 for END
	return (total + cardsPerBlock - 1) / cardsPerBlock
}

