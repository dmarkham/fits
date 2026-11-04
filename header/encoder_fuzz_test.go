package header

import (
	"strings"
	"testing"
)

// FuzzEncodeString: any string value and comment that Set and Encode
// accept must parse back through this package as the same string without
// its trailing blanks (FITS §4.2.1.1) and the same comment without its
// trailing blanks. What the written cards mean to cfitsio and astropy is
// pinned by TestEncodeReadByCfitsio / TestEncodeReadByAstropy on the
// encodeCases; this fuzzer searches for inputs the encoder writes that do
// not read back at all (for example a chain ending in a blank-only
// CONTINUE piece, which leaves the '&' in the value).
func FuzzEncodeString(f *testing.F) {
	for _, c := range encodeCases() {
		if v, ok := c.value.(string); ok {
			f.Add(v, c.comment)
		}
	}
	f.Add(rep("x", 67)+"     ", "")
	f.Add(rep("x", 134)+"  ", "")
	f.Add(rep(" ", 100), "")
	f.Add(rep("a", 70)+"&", "")
	f.Add("short&", rep("comment ", 96)+"end")
	f.Fuzz(func(t *testing.T, v, comment string) {
		h := New()
		if err := h.Set("FUZZSTR", v, comment); err != nil {
			return
		}
		b, err := Encode(h)
		if err != nil {
			return
		}
		cards, _, err := ParseCards(b)
		if err != nil {
			t.Fatalf("ParseCards of encoded header: %v\n%s", err, cardsText(b))
		}
		got := FromCards(cards)
		s, err := got.String("FUZZSTR")
		if want := strings.TrimRight(v, " "); err != nil || s != want {
			t.Fatalf("value read back %q, %v; want %q\n%s", s, err, want, cardsText(b))
		}
		c, _ := got.Get("FUZZSTR")
		if want := strings.TrimRight(comment, " "); c.Comment != want {
			t.Fatalf("comment read back %q; want %q\n%s", c.Comment, want, cardsText(b))
		}
	})
}
