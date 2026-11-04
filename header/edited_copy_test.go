package header

import (
	"bytes"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// A Card taken from a decoded header keeps its Raw bytes when it is copied.
// The tests here edit such copies (Cards() of a header read from a file)
// and check, through cfitsio reading the encoded bytes, that what is
// written is the edited card, never its old bytes: the Header getters
// return the edited fields, so writing the old bytes would make the file
// and the in-memory header disagree without an error.

// primaryLines start every sample here so cfitsio opens the file.
var primaryLines = []string{
	"SIMPLE  =                    T",
	"BITPIX  =                    8",
	"NAXIS   =                    0",
}

func writeBytes(t *testing.T, b []byte, name string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// cfCard returns cfitsio's reading of the only card named key in the
// primary header, failing when there is none or more than one.
func cfCard(t *testing.T, gh *testutil.CfrefHDU, key string) *testutil.CfrefCard {
	t.Helper()
	var found *testutil.CfrefCard
	for i := range gh.Cards {
		if string(gh.Cards[i].Key) == key {
			if found != nil {
				t.Fatalf("cfitsio sees %s twice (cards %d and %d)", key, found.N, gh.Cards[i].N)
			}
			found = &gh.Cards[i]
		}
	}
	if found == nil {
		t.Fatalf("cfitsio does not see %s", key)
	}
	return found
}

// TestEditedDecodedCardIsWrittenAsEdited: a copy of a decoded card whose
// Value, Comment, Type or Key was changed (Raw left as it was) is written
// from its fields, through both FromCards and AddCard; an untouched card
// keeps its bytes. Before the fix Encode wrote the old bytes for all of
// them: EXTIME read back as 5 with comment "s" while h.Float said 7.5.
func TestEditedDecodedCardIsWrittenAsEdited(t *testing.T) {
	src := decodeLines(t, append(primaryLines,
		"EXTIME  =                  5.0 / s",
		"FLAG    =                    T / flag",
		"OLDKEY  =                    5 / renamed",
		"LONGSTR = 'first part of a long string value that goes over the card&'",
		"CONTINUE  'and continues here' / joined comment",
		"HISTORY original history",
		"KEEP    =                   42 / untouched",
		"FLAGB   =                    T / same comment",
		"RATIO   =                 0.25 / same comment",
		"BIG     = 99999999999999999999 / same comment",
		"CPX     =           (1.5, 2.5) / same comment",
		"BADNUM  = 1.5E3D / same comment",
	)...)
	orig := FromCards(src)
	cards := append([]Card(nil), orig.Cards()...)
	for i := range cards {
		c := &cards[i]
		switch c.Key {
		case "EXTIME":
			c.Value, c.Comment = 7.5, "edited"
		case "FLAG":
			c.Type, c.Value = TypeString, "F"
		case "OLDKEY":
			c.Key = "NEWKEY"
		case "LONGSTR":
			c.Value = "short now"
		case "HISTORY":
			c.Comment = "edited history"
		// Value only, same type and comment.
		case "FLAGB":
			c.Value = false
		case "RATIO":
			c.Value = 0.5
		case "BIG":
			c.Value, _ = new(big.Int).SetString("99999999999999999998", 10)
		case "CPX":
			c.Value = Complex{Re: 3.5, Im: -4.5}
		case "BADNUM":
			c.Value = RawValue{Text: "2.5E3D"}
		}
	}
	viaFrom := FromCards(cards)
	viaAdd := New()
	for _, c := range cards {
		if err := viaAdd.AddCard(c); err != nil {
			t.Fatalf("AddCard(%s): %v", c.Key, err)
		}
	}
	encFrom, err := Encode(viaFrom)
	if err != nil {
		t.Fatalf("Encode(FromCards): %v", err)
	}
	encAdd, err := Encode(viaAdd)
	if err != nil {
		t.Fatalf("Encode(AddCard): %v", err)
	}
	if !bytes.Equal(encFrom, encAdd) {
		t.Fatalf("FromCards and AddCard encode differently: %s", firstCardDiff(encFrom, encAdd))
	}
	// The header the copy came from is unchanged.
	if v, err := orig.Float("EXTIME"); err != nil || v != 5 {
		t.Errorf("source header EXTIME = %v, %v; editing a copy changed it", v, err)
	}

	path := writeBytes(t, encFrom, "edited.fits")
	cf := testutil.CfrefDump(t, path)
	gh := &cf.HDUs[0]

	ex := cfCard(t, gh, "EXTIME")
	if ex.Type != "F" || ex.Float == nil || float64(*ex.Float) != 7.5 || string(ex.Comment) != "edited" {
		t.Errorf("cfitsio reads EXTIME as %q (type %s) comment %q; the header says 7.5 / edited", ex.Value, ex.Type, ex.Comment)
	}
	if v, err := viaFrom.Float("EXTIME"); err != nil || v != 7.5 {
		t.Errorf("Float(EXTIME) = %v, %v", v, err)
	}
	fl := cfCard(t, gh, "FLAG")
	if fl.Type != "C" || fl.Str == nil || string(*fl.Str) != "F" || string(fl.Comment) != "flag" {
		t.Errorf("cfitsio reads FLAG as %q (type %s); the header says the string 'F'", fl.Value, fl.Type)
	}
	nk := cfCard(t, gh, "NEWKEY")
	if nk.Int == nil || *nk.Int != 5 || string(nk.Comment) != "renamed" {
		t.Errorf("cfitsio reads NEWKEY as %q / %q", nk.Value, nk.Comment)
	}
	if gh.Card("OLDKEY") != nil {
		t.Errorf("cfitsio still sees OLDKEY: the renamed card was written under its old keyword")
	}
	ls := cfCard(t, gh, "LONGSTR")
	if ls.Str == nil || string(*ls.Str) != "short now" {
		t.Errorf("cfitsio reads LONGSTR as %v; the header says %q", ls.Str, "short now")
	}
	lsl := lsHDU(t, path).card(ls.N)
	if lsl == nil || lsl.Comment == nil || string(*lsl.Comment) != "joined comment" {
		t.Errorf("cfitsio longstr comment of LONGSTR: %+v, want %q", lsl, "joined comment")
	}
	if gh.Card("CONTINUE") != nil {
		t.Error("the old CONTINUE card of LONGSTR was written after the edited value")
	}
	hi := cfCard(t, gh, "HISTORY")
	if string(hi.Comment) != "edited history" {
		t.Errorf("cfitsio reads HISTORY %q, want %q", hi.Comment, "edited history")
	}
	if c := cfCard(t, gh, "FLAGB"); c.Logical == nil || *c.Logical {
		t.Errorf("cfitsio reads FLAGB as %q; the header says F", c.Value)
	}
	if c := cfCard(t, gh, "RATIO"); c.Float == nil || float64(*c.Float) != 0.5 {
		t.Errorf("cfitsio reads RATIO as %q; the header says 0.5", c.Value)
	}
	if c := cfCard(t, gh, "BIG"); string(c.Value) != "99999999999999999998" {
		t.Errorf("cfitsio reads BIG as %q; the header says 99999999999999999998", c.Value)
	}
	if c := lsHDU(t, path).card(cfCard(t, gh, "CPX").N); c == nil || c.Cplx == nil || *c.Cplx != [2]float64{3.5, -4.5} {
		t.Errorf("cfitsio reads CPX as %+v; the header says (3.5, -4.5)", c)
	}
	if c := cfCard(t, gh, "BADNUM"); string(c.Value) != "2.5E3D" || c.Float != nil {
		t.Errorf("cfitsio reads BADNUM as %q; the header says the unparsed text 2.5E3D", c.Value)
	}
	for _, k := range []string{"FLAGB", "RATIO", "BIG", "CPX", "BADNUM"} {
		if c := cfCard(t, gh, k); string(c.Comment) != "same comment" {
			t.Errorf("%s: cfitsio comment %q", k, c.Comment)
		}
	}
	// Untouched cards keep their bytes.
	for _, k := range []string{"KEEP", "SIMPLE"} {
		c, _ := orig.Get(k)
		if w := cfCard(t, gh, k); string(w.Raw) != trimCard(c.Raw) {
			t.Errorf("%s rewritten: %q, source %q", k, w.Raw, trimCard(c.Raw))
		}
	}
	// fitsverify's only complaint is the deliberately malformed BADNUM.
	if r := testutil.Fitsverify(t, path); r.Warnings != 0 || r.Errors != 1 || !strings.Contains(r.Report, `BADNUM: Bad numerical value "2.5E3D"`) {
		t.Errorf("fitsverify:\n%s", r.Report)
	}
}

func lsHDU(t *testing.T, path string) *longstrHDU {
	t.Helper()
	return &hdrrefLongstr(t, path).HDUs[0]
}

// TestReplacedRawIsNotWrittenWithOldContinueCards: a copy of a decoded
// long string whose Raw was replaced by the bytes of a short string card
// still carries the CONTINUE cards of the old chain. Those bytes would
// read back as two cards, so they are not used: the card is written from
// its fields, and cfitsio reads one LONGSTR = 'short' and no CONTINUE.
func TestReplacedRawIsNotWrittenWithOldContinueCards(t *testing.T) {
	src := decodeLines(t, append(primaryLines,
		"LONGSTR = 'first part of a long string value that goes over the card&'",
		"CONTINUE  'and continues here' / joined comment",
	)...)
	c := src[len(src)-1]
	copy(c.Raw[:], padCard("LONGSTR = 'short'"))
	c.Value, c.Comment = "short", ""
	src[len(src)-1] = c
	b, err := Encode(FromCards(src))
	if err != nil {
		t.Fatal(err)
	}
	gh := &testutil.CfrefDump(t, writeBytes(t, b, "replaced.fits")).HDUs[0]
	if w := cfCard(t, gh, "LONGSTR"); w.Str == nil || string(*w.Str) != "short" {
		t.Errorf("cfitsio reads LONGSTR %+v", w)
	}
	if gh.Card("CONTINUE") != nil {
		t.Error("the old CONTINUE card was written")
	}
}

// TestEditedCardThatCannotBeWrittenIsAnError: an edited copy of a decoded
// card that Set could not build (here a keyword cfitsio reads as
// "DATE.OBS") is refused by Encode, not written as its old bytes; the
// unedited card is still written as it was read.
func TestEditedCardThatCannotBeWrittenIsAnError(t *testing.T) {
	src := decodeLines(t, "DATE.OBS= '2020-01-01' / odd keyword")
	if _, err := Encode(FromCards(src)); err != nil {
		t.Fatalf("unedited card: %v", err)
	}
	c := src[0]
	c.Value = "2021-02-02"
	if b, err := Encode(FromCards([]Card{c})); err == nil {
		t.Fatalf("Encode wrote an edited DATE.OBS:\n%s", cardsText(b))
	}
	if err := New().AddCard(c); err == nil {
		t.Fatal("AddCard accepted an edited DATE.OBS that cannot be written")
	}
	// An undefined value given a Value without a type: the value cannot be
	// written as TypeEmpty, and must not be dropped silently.
	u := decodeLines(t, "UNDEF   =")[0]
	if u.Type != TypeEmpty {
		t.Fatalf("UNDEF decoded as %v", u.Type)
	}
	u.Value = int64(5)
	if b, err := Encode(FromCards([]Card{u})); err == nil {
		t.Fatalf("Encode dropped the value of an edited undefined card:\n%s", cardsText(b))
	}
}

// TestRawKeptOnlyWhenItReadsBackAsTheCard: a card with Raw bytes whose
// fields differ from what the bytes decode to is written from its fields.
// The one accepted difference is a float card over an integer literal of
// exactly the same value (the writers' BZERO = 9223372036854775808 and
// BSCALE = 1 cards): Float gives the same float64 from both, BigInt
// gives more on disk, nothing reads back different. The truth is
// cfitsio's reading of the encoded header, compared with the Go getters.
func TestRawKeptOnlyWhenItReadsBackAsTheCard(t *testing.T) {
	rawCard := func(key string, typ ValueType, v any, line string) Card {
		var r [CardWidth]byte
		copy(r[:], padCard(line))
		return Card{Key: key, Type: typ, Value: v, Raw: r}
	}
	two63 := math.Ldexp(1, 63)
	cases := []struct {
		card    Card
		keepRaw bool
	}{
		// What write_image.go builds for a uint64 / unsigned image.
		{rawCard("BZERO", TypeFloat, two63, "BZERO   =  9223372036854775808"), true},
		{rawCard("BSCALE", TypeFloat, 1.0, "BSCALE  =                    1"), true},
		{rawCard("ZEROI", TypeFloat, 32768.0, "ZEROI   =                32768"), true},
		// 2^53+1 is not the float64 2^53: the bytes would read back as
		// another integer than the header holds.
		{rawCard("NEAR53", TypeFloat, 9007199254740992.0, "NEAR53  =     9007199254740993"), false},
		// An integer card over float bytes: Int would fail on disk.
		{rawCard("INTF", TypeInt, int64(1), "INTF    =                  1.0"), false},
		// +0 over "-0": Float of the bytes is -0.
		{rawCard("ZSIGN", TypeFloat, 0.0, "ZSIGN   =                   -0"), false},
		// An integer 0 over "-0": Float of the bytes is -0, of the card +0.
		{rawCard("ZINT", TypeInt, int64(0), "ZINT    =                   -0"), false},
		// Same value, other comment.
		{rawCard("CMT", TypeInt, int64(3), "CMT     =                    3 / old"), false},
	}
	h := FromCards(decodeLines(t, primaryLines...))
	for _, c := range cases {
		if err := h.AddCard(c.card); err != nil {
			t.Fatalf("AddCard(%s): %v", c.card.Key, err)
		}
	}
	enc, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	cf := testutil.CfrefDump(t, writeBytes(t, enc, "rawkeep.fits"))
	gh := &cf.HDUs[0]
	for _, c := range cases {
		k := c.card.Key
		w := cfCard(t, gh, k)
		if kept := string(w.Raw) == trimCard(c.card.Raw); kept != c.keepRaw {
			t.Errorf("%s: written as %q; keep the source bytes = %v, want %v", k, w.Raw, kept, c.keepRaw)
		}
		if w.Float == nil {
			t.Errorf("%s: cfitsio cannot read %q as a number: %s", k, w.Value, w.ReadError)
			continue
		}
		f, err := h.Float(k)
		if err != nil || math.Float64bits(f) != math.Float64bits(float64(*w.Float)) {
			t.Errorf("%s: Float = %v, %v; cfitsio fits_read_key_dbl %v", k, f, err, float64(*w.Float))
		}
		if string(w.Comment) != c.card.Comment {
			t.Errorf("%s: cfitsio comment %q, the header says %q", k, w.Comment, c.card.Comment)
		}
		// The disk value is exactly the in-memory value.
		mem := new(big.Float).SetFloat64(f)
		if iv, err := h.BigInt(k); err == nil {
			mem = new(big.Float).SetInt(iv)
		}
		disk, ok := new(big.Float).SetPrec(256).SetString(strings.Replace(string(w.Value), "D", "E", 1))
		if !ok || disk.Cmp(mem) != 0 {
			t.Errorf("%s: written as %q, the header holds %s exactly", k, w.Value, mem.Text('g', 25))
		}
		if c.card.Type == TypeInt && (w.Type != "I" || w.Int == nil) {
			t.Errorf("%s: in memory an integer, cfitsio reads type %s", k, w.Type)
		}
	}
}

// TestBlankKeyLongStringReencodes: a malformed card whose keyword cfitsio
// names "" (a blank in column 1) but which has a value indicator and a
// string ending in '&' takes the CONTINUE card after it into its value.
// Unmodified, the header must re-encode to its own bytes. (Found by
// FuzzParseHeader: Encode counted the two physical cards of the
// blank-keyed card as two read-back cards and refused the header.)
func TestBlankKeyLongStringReencodes(t *testing.T) {
	lines := []string{
		" 0000000000000000000000000= '00000000000000&'00000000000000000000000000000000000",
		"CONTINUE  0000000000000000000000000000000000000000000000000000000000000000000000",
	}
	cards := decodeLines(t, lines...)
	if len(cards) != 1 || cards[0].Key != "" || len(cards[0].RawCards()) != 2 {
		t.Fatalf("decoded %d cards, first key %q with %d physical cards; want one blank-keyed card over two", len(cards), cards[0].Key, len(cards[0].RawCards()))
	}
	b, err := Encode(FromCards(cards))
	if err != nil {
		t.Fatal(err)
	}
	want := padCard(lines[0]) + padCard(lines[1]) + padCard("END")
	if string(b[:len(want)]) != want {
		t.Errorf("re-encoded as\n%s", cardsText(b))
	}
}

func padCard(s string) string { return s + strings.Repeat(" ", CardWidth-len(s)) }
