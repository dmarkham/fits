package fits

import (
	"bytes"
	"errors"
)

// ErrMalformedEND is returned by Open when a header holds a card that
// readers take for the END card, but that is not one: "END" followed by
// something other than 77 blanks ("END  = 1", "END=1", "END     junk",
// "END" padded with NUL bytes, "HIERARCH END = 1"). cfitsio names a card by
// the text up to the first blank or '=' (ffgknm), and a HIERARCH card by
// the text between "HIERARCH " and '=' with blanks trimmed, so for it such
// a card ends the header (and
// one where a required keyword must be is an error), while the END search
// of the header scan and the FITS Standard (4.0, 4.4.1.1: columns 9-80
// blank) need the exact card. With such a card the header and the data
// unit have no single reading, so the file is refused.
var ErrMalformedEND = errors.New("fits: a card named END that is not the END card (END followed by 77 blanks)")

// hierarchPrefix starts a HIERARCH card (cfitsio compares it
// case-sensitively, FSTRNCMP(card, "HIERARCH ", 9)).
var hierarchPrefix = []byte("HIERARCH ")

// endCard is the one valid END card.
var endCard = append([]byte("END"), bytes.Repeat([]byte(" "), 77)...)

// malformedEND reports whether the 80-byte card raw is named END the way
// cfitsio names cards (columns up to the first blank, '=' or NUL; for a
// HIERARCH card the text before '=') but is not the END card.
func malformedEND(raw []byte) bool {
	if bytes.HasPrefix(raw, hierarchPrefix) {
		// ffgknm's HIERARCH branch: the name is the text between
		// "HIERARCH " and the first '=' (strchr, so a NUL before it
		// hides it), leading and trailing blanks removed; ffgphd and
		// ffbinit compare it with "END" case-sensitively.
		s := raw
		if i := bytes.IndexByte(s, 0); i >= 0 {
			s = s[:i]
		}
		eq := bytes.IndexByte(s, '=')
		if eq < 0 {
			return false
		}
		return string(bytes.Trim(s[len(hierarchPrefix):eq], " ")) == "END"
	}
	if len(raw) < 4 || raw[0] != 'E' || raw[1] != 'N' || raw[2] != 'D' {
		return false
	}
	switch raw[3] {
	case ' ', '=', 0:
		return !bytes.Equal(raw, endCard)
	}
	return false
}
