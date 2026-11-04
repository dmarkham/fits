package header

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Write rules of the encoder and lookup rules of Header that each had a
// mutant no test killed. Every test states its truth: cfitsio 4.7.0 writing
// or reading the same keywords (hdrref, cfref), or a refusal whose reason
// is shown by cfitsio reading the bytes that would otherwise be written.

func fileCards(gh *testutil.CfrefHDU) map[string]string {
	out := map[string]string{}
	for _, c := range gh.Cards {
		out[string(c.Key)] = string(c.Raw)
	}
	return out
}

// TestFixedFormatMatchesCfitsioBytes: integer and logical cards are
// written exactly as cfitsio's fits_write_key writes them, value
// right-justified to column 30 for every printed width from 1 to 20
// characters, with and without a comment. Truth: "hdrref writekeys"
// (cfitsio writes the same keywords); the card bytes must be equal.
func TestFixedFormatMatchesCfitsioBytes(t *testing.T) {
	cfPath := filepath.Join(t.TempDir(), "cfkeys.fits")
	runHdrrefWriting(t, []string{cfPath}, "writekeys", cfPath)
	cf := fileCards(&testutil.CfrefDump(t, cfPath).HDUs[0])

	h := New()
	cmt := func(odd bool) string {
		if odd {
			return "c"
		}
		return ""
	}
	var keys []string
	set := func(k string, v any, c string) {
		mustSet(t, h, k, v, c)
		keys = append(keys, k)
	}
	p := int64(0)
	for k := 1; k <= 18; k++ {
		p = p*10 + 9
		set(fmt.Sprintf("P%02d", k), p, cmt(k%2 == 1))
		set(fmt.Sprintf("N%02d", k), -p, cmt(k%2 == 1))
	}
	set("IMAX", int64(math.MaxInt64), "c")
	set("IMIN", int64(math.MinInt64), "")
	set("NBIG", int64(-1000000000000000000), "")
	set("UMAX", uint64(math.MaxUint64), "c")
	set("LT", true, "c")
	set("LF", false, "")
	for _, k := range keys {
		c, _ := h.Get(k)
		enc, err := EncodeCard(c)
		if err != nil {
			t.Fatal(err)
		}
		w, ok := cf[k]
		if !ok {
			t.Fatalf("cfitsio did not write %s", k)
		}
		if got := trimCard(enc); got != w {
			t.Errorf("%s:\n got %q\nwant %q (cfitsio fits_write_key)", k, got, w)
		}
	}
}

// TestInMemoryContinueRefusedOnEveryPath: a card named CONTINUE built in
// memory, under that name in either case, is refused by AddCard,
// EncodeCard(WithContinuation) and Encode
// (FromCards does not validate), whatever its type. The reason, from
// cfitsio: the bytes such a card would be written as ("CONTINUE= 'abc'")
// carry no value for cfitsio, which reads them as a CONTINUE card with the
// comment "= 'abc'", so the value would be lost; after a string ending in
// '&' they would be joined into it.
func TestInMemoryContinueRefusedOnEveryPath(t *testing.T) {
	lines := append(append([]string(nil), primaryLines...),
		"CONTINUE= 'abc'",
		"CONTINUE=                    5 / five",
	)
	var buf []byte
	for _, l := range append(lines, "END") {
		buf = append(buf, padCard(l)...)
	}
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	gh := testutil.CfrefDump(t, writeBytes(t, buf, "continue_eq.fits")).HDUs[0]
	for _, n := range []int{4, 5} {
		c := gh.Cards[n-1]
		if string(c.Key) != "CONTINUE" || c.Type != "" || c.Value != "" {
			t.Fatalf("cfitsio reads card %d (%q) as key %q value %q type %q; the refusal assumes it has no value", n, c.Raw, c.Key, c.Value, c.Type)
		}
	}
	if w := string(gh.Cards[3].Comment); w != "= 'abc'" {
		t.Fatalf("cfitsio comment of CONTINUE= 'abc' is %q", w)
	}

	for _, key := range []string{"CONTINUE", "continue"} {
		for _, c := range []Card{
			{Key: key, Type: TypeString, Value: "abc"},
			{Key: key, Type: TypeInt, Value: int64(5), Comment: "five"},
			{Key: key, Type: TypeFloat, Value: 1.5},
			{Key: key, Type: TypeLogical, Value: true},
			{Key: key, Type: TypeEmpty, Comment: "  'abc'"},
		} {
			if err := New().AddCard(c); !errors.Is(err, ErrInvalidKeyword) {
				t.Errorf("AddCard(%+v): err %v, want ErrInvalidKeyword", c, err)
			}
			if _, err := EncodeCardWithContinuation(c); !errors.Is(err, ErrInvalidKeyword) {
				t.Errorf("EncodeCardWithContinuation(%+v): err %v, want ErrInvalidKeyword", c, err)
			}
			if b, err := Encode(FromCards([]Card{c})); !errors.Is(err, ErrInvalidKeyword) {
				t.Errorf("Encode(FromCards(%+v)): err %v, want ErrInvalidKeyword; written:\n%s", c, err, cardsText(b))
			}
		}
	}
}

// TestEncodeJoinErrorNamesTheSwallowedCard: a one-card string ending in
// '&', a lone CONTINUE card decoded from a file and another card after
// it. cfitsio reads those bytes as ONE string, "abcdef" (shown below on
// the bytes written by hand), so Encode refuses them, and its error names
// the CONTINUE card that would be swallowed, not a later one.
func TestEncodeJoinErrorNamesTheSwallowedCard(t *testing.T) {
	lone := decodeLines(t, "CONTINUE  'def'")
	h := FromCards(decodeLines(t, primaryLines...))
	mustSet(t, h, "A", "abc&", "")
	if err := h.AddCard(lone[0]); err != nil {
		t.Fatal(err)
	}
	mustSet(t, h, "B", int64(1), "")

	// What would be written without the check.
	var buf []byte
	for _, c := range h.Cards() {
		enc, err := EncodeCard(c)
		if err != nil {
			t.Fatal(err)
		}
		buf = append(buf, enc[:]...)
	}
	buf = append(buf, padCard("END")...)
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	path := writeBytes(t, buf, "joined.fits")
	gh := testutil.CfrefDump(t, path).HDUs[0]
	a := cfCard(t, &gh, "A")
	if a.Str == nil || string(*a.Str) != "abcdef" {
		t.Fatalf("cfitsio reads A as %v; the refusal assumes the join %q", a.Str, "abcdef")
	}

	_, err := Encode(h)
	if err == nil {
		t.Fatal("Encode accepted cards that cfitsio reads back joined")
	}
	if want := "card 4 (CONTINUE)"; !strings.Contains(err.Error(), want) {
		t.Errorf("Encode error %q does not name %s, the card that is swallowed", err, want)
	}
}

// TestLongCommentWithDoubleBlanksRefused: the comment of a CONTINUE'd
// string can only be broken at a single blank, which readers restore as
// the joining blank. A comment with only double blanks is refused. The
// reason, from cfitsio: a card's trailing blanks are not part of its
// comment ("X = 1 / ab  " reads back as "ab"), so a break after the first
// blank of a pair would lose it.
func TestLongCommentWithDoubleBlanksRefused(t *testing.T) {
	var buf []byte
	for _, l := range append(append([]string(nil), primaryLines...), "X       =                    1 / ab  ", "END") {
		buf = append(buf, padCard(l)...)
	}
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	gh := testutil.CfrefDump(t, writeBytes(t, buf, "trail.fits")).HDUs[0]
	if c := cfCard(t, &gh, "X"); string(c.Comment) != "ab" {
		t.Fatalf("cfitsio comment %q, want \"ab\"", c.Comment)
	}

	h := New()
	mustSet(t, h, "LONGSTR", rep("v", 100), rep("ab  ", 120))
	b, err := Encode(h)
	if err == nil {
		cards, _, _ := ParseCards(b)
		got := ""
		if len(cards) > 0 {
			got = cards[0].Comment
		}
		t.Fatalf("Encode accepted the comment; it reads back as %q, want %q", got, rep("ab  ", 120))
	}
	if !strings.Contains(err.Error(), "cannot be continued") {
		t.Errorf("error %q, want the comment-continuation refusal", err)
	}
}

// dupLines is a primary header with duplicated keywords: three DUPA, two
// DUPS, and commentary cards around them.
var dupLines = append(append([]string(nil), primaryLines...),
	"HISTORY h1",
	"DUPA    =                    1 / first",
	"DUPS    = 'x'",
	"HISTORY h2",
	"DUPA    =                    2 / second",
	"COMMENT c1",
	"DUPS    = 'y'",
	"DUPA    =                    3 / third",
	"HISTORY h3",
)

func linesFile(t *testing.T, lines []string, name string) (string, []byte) {
	t.Helper()
	var buf []byte
	for _, l := range append(append([]string(nil), lines...), "END") {
		buf = append(buf, padCard(l)...)
	}
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	return writeBytes(t, buf, name), buf
}

// TestDuplicateKeywordFirstCardMatchesCfitsio: for a keyword that appears
// more than once in a header READ FROM A FILE (FromCards), or copied card
// by card with AddCard, every getter returns the first card, as cfitsio's
// fits_read_keyword does from the top of the header; after Delete, which
// removes the first card as fits_delete_key does, lookups find the next
// one. Truth: hdrref readkey / delkey on the same bytes, deleting twice.
func TestDuplicateKeywordFirstCardMatchesCfitsio(t *testing.T) {
	path, buf := linesFile(t, dupLines, "dup.fits")
	cards, _, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	fromFile := FromCards(cards)
	added := New()
	for _, c := range cards {
		if err := added.AddCard(c); err != nil {
			t.Fatal(err)
		}
	}
	checkFirst := func(label string, h *Header, path string) {
		t.Helper()
		cf := hdrrefReadkey(t, path, "dupa", "DUPS")
		for _, k := range []string{"dupa", "DUPS"} {
			r := cf[k]
			c, ok := h.Get(k)
			if r.Status != 0 {
				if ok {
					t.Errorf("%s: Get(%q) finds %q; cfitsio status %d", label, k, trimCard(c.Raw), r.Status)
				}
				continue
			}
			if !ok {
				t.Errorf("%s: Get(%q) finds nothing; cfitsio %q", label, k, *r.Value)
				continue
			}
			enc, err := EncodeCard(c)
			if err != nil {
				t.Fatal(err)
			}
			vc, cm, _ := valueAndComment(enc)
			if vc != *r.Value || cm != *r.Comment {
				t.Errorf("%s: Get(%q) gives %q / %q; cfitsio's first card %q / %q", label, k, vc, cm, *r.Value, *r.Comment)
			}
		}
	}
	checkFirst("FromCards", fromFile, path)
	checkFirst("AddCard", added, path)
	if v, err := fromFile.Int("DUPA"); err != nil || v != 1 {
		t.Errorf("FromCards: Int(DUPA) = %d, %v", v, err)
	}
	if v, err := added.String("dups"); err != nil || v != "x" {
		t.Errorf("AddCard: String(dups) = %q, %v", v, err)
	}

	// Delete twice, each time against cfitsio deleting from the same file.
	for _, h := range []*Header{fromFile, added} {
		cur := path
		for round := 1; round <= 2; round++ {
			if err := h.Delete("DupA"); err != nil {
				t.Fatal(err)
			}
			next := filepath.Join(t.TempDir(), fmt.Sprintf("del%d.fits", round))
			runHdrrefWriting(t, []string{next}, "delkey", cur, next, "DUPA")
			want, err := os.ReadFile(next)
			if err != nil {
				t.Fatal(err)
			}
			ours, err := Encode(h)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(ours, want) {
				t.Fatalf("after %d Delete(DupA): header differs from cfitsio fits_delete_key: %s", round, firstCardDiff(ours, want))
			}
			checkFirst(fmt.Sprintf("after %d deletes", round), h, next)
			cur = next
		}
		if v, err := h.Int("dupa"); err != nil || v != 3 {
			t.Errorf("after two deletes Int(dupa) = %d, %v; want the third card's 3", v, err)
		}
	}
}

// TestDeleteCommentaryRemovesFirstMatchesCfitsio: Delete("HISTORY")
// removes the first HISTORY card, as cfitsio's fits_delete_key does, also
// when it is the very first card of the header (an extension header built
// from user cards). Truth: hdrref delkey on the encoded file; the header
// without the primary cards must lose the same card.
func TestDeleteCommentaryRemovesFirstMatchesCfitsio(t *testing.T) {
	path, buf := linesFile(t, dupLines, "hist.fits")
	cards, _, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	full := FromCards(cards)
	for round := 1; round <= 2; round++ {
		if err := full.Delete("history"); err != nil {
			t.Fatal(err)
		}
		next := filepath.Join(t.TempDir(), fmt.Sprintf("hist%d.fits", round))
		runHdrrefWriting(t, []string{next}, "delkey", path, next, "HISTORY")
		want, err := os.ReadFile(next)
		if err != nil {
			t.Fatal(err)
		}
		ours, err := Encode(full)
		if err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(ours, want) {
			t.Fatalf("after %d Delete(history): differs from cfitsio fits_delete_key: %s", round, firstCardDiff(ours, want))
		}
		path = next
	}
	if got := full.History(); len(got) != 1 || got[0] != "h3" {
		t.Errorf("History after two deletes %q, want [h3]", got)
	}

	// The same cards without the primary keywords: HISTORY h1 is card 0.
	bare := FromCards(decodeLines(t, dupLines[len(primaryLines):]...))
	for range 2 {
		if err := bare.Delete("HISTORY"); err != nil {
			t.Fatal(err)
		}
	}
	got, want := bare.Cards(), full.Cards()[len(primaryLines):]
	if len(got) != len(want) {
		t.Fatalf("bare header has %d cards after deletes, the full one %d", len(got), len(want))
	}
	for i := range got {
		if got[i].Raw != want[i].Raw {
			t.Errorf("card %d: bare %q, cfitsio-checked %q", i, trimCard(got[i].Raw), trimCard(want[i].Raw))
		}
	}
}

// TestRawValueRefusalNamesItsReason: the two RawValue texts with their own
// refusal reason are refused with that reason (an empty text reads back
// as an undefined value, "(1 2" is a complex group cfitsio's ffpsvc cannot
// close; cfitsio's reading of both is in TestRawValueReadBackByGoAndCfitsio
// and rawBad).
func TestRawValueRefusalNamesItsReason(t *testing.T) {
	for text, reason := range map[string]string{
		"":     "is empty",
		"(1 2": "without ')'",
	} {
		err := New().Set("RAW", RawValue{Text: text}, "")
		if err == nil || !strings.Contains(err.Error(), reason) {
			t.Errorf("Set(RawValue %q): %v, want the reason %q", text, err, reason)
		}
	}
}
