package header

import (
	"fmt"
	"os"
	"strings"
	"testing"
)

// decodeLines decodes card images as ParseCards reads them from a file.
func decodeLines(t *testing.T, lines ...string) []Card {
	t.Helper()
	var buf []byte
	for _, l := range lines {
		buf = append(buf, fmt.Sprintf("%-80s", l)...)
	}
	buf = append(buf, fmt.Sprintf("%-80s", "END")...)
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	cards, _, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	return cards
}

// TestSetRefusesWhatWouldNotReadBack: keys and values that Set used to
// write but that read back as something else are refused.
func TestSetRefusesWhatWouldNotReadBack(t *testing.T) {
	for _, b := range []struct {
		name, key string
		value     any
		comment   string
	}{
		// The decoder collapses blank runs in a HIERARCH path (as the Card
		// doc says), so these keys would change on a round trip.
		{"HIERARCH blank run", "HIERARCH ESO  DBL", 5, ""},
		{"HIERARCH leading blank", "HIERARCH  ESO X", 5, ""},
		{"HIERARCH trailing blank", "HIERARCH ESO X ", 5, ""},
		// A CONTINUE card built by hand joins a preceding string ending in
		// '&' on read; CONTINUE cards are written by the encoder only.
		{"CONTINUE commentary", "CONTINUE", nil, "  'def'"},
		{"CONTINUE commentary string", "CONTINUE", "'def'", ""},
		// A commentary card has one text; a value and a comment would
		// silently lose the value.
		{"COMMENT value and comment", "COMMENT", "x", "y"},
		{"HISTORY value and comment", "HISTORY", "x", "y"},
	} {
		h := New()
		if err := h.Set(b.key, b.value, b.comment); err == nil {
			enc, _ := Encode(h)
			t.Errorf("%s: Set(%q, %#v, %q) accepted it; written as:\n%s", b.name, b.key, b.value, b.comment, cardsText(enc))
		}
		if err := New().Add(b.key, b.value, b.comment); err == nil {
			t.Errorf("%s: Add(%q, %#v, %q) accepted it", b.name, b.key, b.value, b.comment)
		}
	}
	// The single-blank HIERARCH form is accepted and reads back unchanged.
	h := New()
	mustSet(t, h, "HIERARCH ESO DBL", 5, "")
	enc, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	cards, _, _ := ParseCards(enc)
	if len(cards) != 1 || cards[0].Key != "HIERARCH ESO DBL" {
		t.Errorf("HIERARCH ESO DBL read back as %+v", cards)
	}
}

// TestEncodeRefusesJoiningCards: Encode refuses a header whose cards would
// read back as a different set of cards. A one-card string ending in '&'
// followed by a lone CONTINUE card (here decoded from a file) is written as
// two cards but would read back as one joined string. Before the fix
// A = "abc&" + CONTINUE 'def' was written and read back as A = "abcdef".
func TestEncodeRefusesJoiningCards(t *testing.T) {
	lone := decodeLines(t, "CONTINUE  'def'")
	if len(lone) != 1 || lone[0].Key != "CONTINUE" {
		t.Fatalf("lone CONTINUE decoded as %+v", lone)
	}
	h := New()
	mustSet(t, h, "A", "abc&", "")
	if err := h.AddCard(lone[0]); err != nil {
		t.Fatal(err)
	}
	if enc, err := Encode(h); err == nil {
		cards, _, _ := ParseCards(enc)
		t.Errorf("Encode accepted it; written:\n%s\nread back as %d cards: %+v", cardsText(enc), len(cards), cards)
	}
	// The same cards as read from a file (both decoded, the CONTINUE empty
	// so it was not joined) are re-emitted as they were.
	file := decodeLines(t, "A       = 'abc&'", "CONTINUE")
	if len(file) != 2 {
		t.Fatalf("decoded %d cards", len(file))
	}
	if _, err := Encode(FromCards(file)); err != nil {
		t.Errorf("Encode of an unmodified decoded header: %v", err)
	}
}

// TestAddCard: a card decoded from a file is appended with its original
// bytes, even when it is not a card Set could build; a card built in
// memory is validated as Set validates it.
func TestAddCard(t *testing.T) {
	src := decodeLines(t,
		"KEYLONG10= 5 / cfitsio names this KEYLONG10",
		"DATE.OBS= '2020-01-01'",
		"LONGSTR = 'first part of a long string value that goes over the card&'",
		"CONTINUE  'and continues here' / joined comment",
		"HISTORY a history line",
	)
	h := New()
	for _, c := range src {
		if err := h.AddCard(c); err != nil {
			t.Fatalf("AddCard(%s): %v", c.Key, err)
		}
	}
	enc, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	var want []byte
	for i := range src {
		for _, r := range src[i].RawCards() {
			want = append(want, r[:]...)
		}
	}
	if string(enc[:len(want)]) != string(want) {
		t.Errorf("AddCard did not keep the source bytes:\n got %s\nwant %s", cardsText(enc[:len(want)]), cardsText(want))
	}
	for _, k := range []string{"KEYLONG10", "DATE.OBS", "LONGSTR"} {
		if !h.Has(k) {
			t.Errorf("Has(%q) false after AddCard", k)
		}
	}
	if got := h.History(); len(got) != 1 || got[0] != "a history line" {
		t.Errorf("History %q", got)
	}
	for _, c := range []Card{
		{Key: "DATE.OBS", Type: TypeString, Value: "x"},
		{Key: "KEYLONG10", Type: TypeInt, Value: int64(1)},
		{Key: "OK", Type: TypeString, Value: "x", Comment: "caf\xe9"},
		{Key: "END"},
	} {
		if err := New().AddCard(c); err == nil {
			t.Errorf("AddCard accepted in-memory card %+v", c)
		}
	}
	if err := h.AddCard(Card{Key: "lower", Type: TypeInt, Value: int64(3)}); err != nil {
		t.Errorf("AddCard(lower): %v", err)
	} else if v, err := h.Int("LOWER"); err != nil || v != 3 {
		t.Errorf("Int(LOWER) = %d, %v", v, err)
	}
}

// TestJoinedContinueCommentIsWholeComment: past 72 characters a joined
// CONTINUE comment cannot be checked against cfitsio, which truncates it
// (fits_read_key_longstr). The truth is the comment hdrref.c gave cfitsio
// to write in cfwritten.fits: c120 = rep("comment text ", 120), the
// pattern repeated to 120 characters. Go must read back exactly that.
func TestJoinedContinueCommentIsWholeComment(t *testing.T) {
	p := headerTestdata(t, "cfwritten.fits")
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	cards, _, err := ParseCards(b)
	if err != nil {
		t.Fatal(err)
	}
	c, ok := FromCards(cards).Get("LONGSTR")
	want := rep("comment text ", 120)
	if !ok || c.Comment != want {
		t.Errorf("LONGSTR comment\n got %q\nwant %q (hdrref.c c120)", c.Comment, want)
	}
	if !strings.HasPrefix(want, "comment text") || len(want) != 120 {
		t.Fatal("bad want")
	}
}
