package header

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// rawGood are RawValue texts that a card can hold: each is one value token
// that cfitsio's ffpsvc splits off whole and that does not parse as a FITS
// value, so it reads back as the same TypeUnparsed text.
var rawGood = []struct{ text, comment string }{
	{"1.5E3D", "malformed exponent"},
	{"TRUE", ""},
	{"0x1F", "hex"},
	{"1_000", ""},
	{"++5", "two signs"},
	{"(1, )", "empty imaginary part"},
	{"(1,2,3)", ""},
	{"1E", ""},
	{"abc'", "quote inside"},
	{"x(", ""},
	{"NaN", ""},
	{"-Inf", "c"},
	{strings.Repeat("9", 30) + ".x", ""},
	{strings.Repeat("u", 70), ""},
}

// rawBad are RawValue texts that would be written verbatim and read back as
// something else: another value, another type, or a different comment.
var rawBad = []string{
	"",                        // reads back as an undefined value
	" x",                      // leading blank: still "x", but not the text given
	"x ",                      // trailing blank
	"a b",                     // value "a", comment "b"
	"5 / injected",            // TypeInt 5 with comment "injected"
	"x/y",                     // value "x", comment "y"
	"'str'",                   // TypeString
	"'unterminated",           // TypeString
	"(1,2)",                   // TypeComplexInt
	"(1 2",                    // no ')': ffpsvc fails
	"(1, 2) x",                // value "(1, 2)"
	"T",                       // TypeLogical
	"F",                       // TypeLogical
	"5",                       // TypeInt
	"-0",                      // TypeInt
	"1.5",                     // TypeFloat
	"1.0D3",                   // TypeFloat
	"99999999999999999999999", // TypeBigInt
	strings.Repeat("u", 71),   // cut at 70 by cfitsio
}

// TestRawValueRefusedUnlessItReadsBack: Set refuses a RawValue that would
// not read back as the same unparsed text, and Encode refuses such a card
// built in memory. Before the fix Set and Encode accepted every printable
// ASCII text and wrote it verbatim (for example RawValue{"5 / injected"}
// read back as TypeInt 5 with comment "injected").
func TestRawValueRefusedUnlessItReadsBack(t *testing.T) {
	for _, txt := range rawBad {
		h := New()
		if err := h.Set("RAW", RawValue{Text: txt}, "c"); err == nil {
			enc, _ := Encode(h)
			cards, _, _ := ParseCards(enc)
			t.Errorf("Set(RawValue{%q}) accepted it; it reads back as %+v", txt, cards[0])
		}
		c := Card{Key: "RAW", Type: TypeUnparsed, Value: RawValue{Text: txt}}
		if _, err := Encode(FromCards([]Card{c})); err == nil {
			t.Errorf("Encode accepted in-memory RawValue{%q}", txt)
		}
	}
}

// TestRawValueReadBackByGoAndCfitsio: every rawGood text written through
// Set reads back as the same TypeUnparsed text and comment through this
// package, and cfitsio (cfref dump) reports the same value text and
// comment.
func TestRawValueReadBackByGoAndCfitsio(t *testing.T) {
	h := New()
	mustSet(t, h, "SIMPLE", true, "")
	mustSet(t, h, "BITPIX", int64(8), "")
	mustSet(t, h, "NAXIS", int64(0), "")
	key := func(i int) string { return fmt.Sprintf("RAW%02d", i) }
	for i, g := range rawGood {
		mustSet(t, h, key(i), RawValue{Text: g.text}, g.comment)
	}
	b, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	cards, _, err := ParseCards(b)
	if err != nil {
		t.Fatal(err)
	}
	got := FromCards(cards)
	for i, g := range rawGood {
		c, ok := got.Get(key(i))
		rv, isRaw := c.Value.(RawValue)
		if !ok || c.Type != TypeUnparsed || !isRaw || rv.Text != g.text || c.Comment != g.comment {
			t.Errorf("%s: RawValue{%q} / %q read back as %+v", key(i), g.text, g.comment, c)
		}
	}
	p := filepath.Join(t.TempDir(), "raw.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	cf := testutil.CfrefDump(t, p)
	if len(cf.HDUs) != 1 {
		t.Fatalf("cfitsio: %d HDUs", len(cf.HDUs))
	}
	for i, g := range rawGood {
		w := cf.HDUs[0].Card(key(i))
		if w == nil || string(w.Value) != g.text || string(w.Comment) != g.comment {
			t.Errorf("%s: cfitsio reads %+v, want value %q comment %q", key(i), w, g.text, g.comment)
		}
	}
}

// FuzzEncodeRawValue: any RawValue and comment that Set and Encode accept
// read back as the same TypeUnparsed text and comment (without trailing
// blanks).
func FuzzEncodeRawValue(f *testing.F) {
	for _, g := range rawGood {
		f.Add(g.text, g.comment)
	}
	for _, b := range rawBad {
		f.Add(b, "c")
	}
	f.Fuzz(func(t *testing.T, text, comment string) {
		h := New()
		if err := h.Set("FUZZRAW", RawValue{Text: text}, comment); err != nil {
			return
		}
		b, err := Encode(h)
		if err != nil {
			return
		}
		cards, _, err := ParseCards(b)
		if err != nil {
			t.Fatalf("ParseCards: %v\n%s", err, cardsText(b))
		}
		if len(cards) != 1 {
			t.Fatalf("%d cards read back, want 1\n%s", len(cards), cardsText(b))
		}
		c := cards[0]
		rv, ok := c.Value.(RawValue)
		if c.Key != "FUZZRAW" || c.Type != TypeUnparsed || !ok || rv.Text != text || c.Comment != strings.TrimRight(comment, " ") {
			t.Fatalf("RawValue{%q} / %q read back as %+v\n%s", text, comment, c, cardsText(b))
		}
	})
}
