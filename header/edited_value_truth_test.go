package header

import (
	"bytes"
	"math"
	"math/big"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Edited copies of decoded cards, one test per value rule the encoder's
// "do the bytes still decode to the card" check depends on. In every test
// the truth is cfitsio 4.7.0 reading the encoded header (cfref dump,
// hdrref longstr) compared with the Header getters of the header that was
// encoded, bit for bit: a card written as its old bytes, or a getter that
// says something else than the bytes, makes them disagree.

// editedHeaders decodes primaryLines+lines, applies edit to a copy of the
// cards and returns the edited header built both ways (FromCards and
// AddCard), after checking that both encode to the same bytes, and the
// path of the encoded file.
func editedHeaders(t *testing.T, edit func(c *Card), lines ...string) (fromCards, added *Header, path string) {
	t.Helper()
	src := FromCards(decodeLines(t, append(append([]string(nil), primaryLines...), lines...)...))
	cards := append([]Card(nil), src.Cards()...)
	for i := range cards {
		edit(&cards[i])
	}
	fromCards = FromCards(cards)
	added = New()
	for _, c := range cards {
		if err := added.AddCard(c); err != nil {
			t.Fatalf("AddCard(%s): %v", c.Key, err)
		}
	}
	encFrom, err := Encode(fromCards)
	if err != nil {
		t.Fatalf("Encode(FromCards): %v", err)
	}
	encAdd, err := Encode(added)
	if err != nil {
		t.Fatalf("Encode(AddCard): %v", err)
	}
	if !bytes.Equal(encFrom, encAdd) {
		t.Fatalf("FromCards and AddCard encode differently: %s", firstCardDiff(encFrom, encAdd))
	}
	return fromCards, added, writeBytes(t, encFrom, "edited.fits")
}

// checkNumberAsCfitsio: Float (and Int, where the card is an integer) of
// key in every header equal cfitsio's fits_read_key_dbl / _lnglng of the
// encoded card, bit for bit.
func checkNumberAsCfitsio(t *testing.T, gh *testutil.CfrefHDU, key string, hs ...*Header) {
	t.Helper()
	w := cfCard(t, gh, key)
	if w.Float == nil {
		t.Errorf("%s: cfitsio cannot read %q as a number: %s", key, w.Raw, w.ReadError)
		return
	}
	for i, h := range hs {
		f, err := h.Float(key)
		if err != nil || math.Float64bits(f) != math.Float64bits(float64(*w.Float)) {
			t.Errorf("%s (header %d): Float = %v (sign bit %v), %v; cfitsio reads the written card %q as %v (sign bit %v)",
				key, i, f, math.Signbit(f), err, w.Raw, float64(*w.Float), math.Signbit(float64(*w.Float)))
		}
		iv, err := h.Int(key)
		switch {
		case err == nil && (w.Type != "I" || w.Int == nil || *w.Int != iv):
			t.Errorf("%s (header %d): Int = %d; cfitsio reads %q as type %s int %v", key, i, iv, w.Raw, w.Type, w.Int)
		case err != nil && w.Type == "I":
			t.Errorf("%s (header %d): Int: %v; cfitsio reads the integer %q", key, i, err, w.Raw)
		}
	}
}

// TestEditedNegZeroIntegerReadsBackAsHeld: a decoded integer card written
// "-0" (cfitsio's fits_read_key_dbl reads it as -0, and so does Float).
// A copy whose Value was edited to 5000 must read 5000 everywhere; a copy
// whose comment alone was edited still holds -0 and must be written so
// cfitsio reads -0; an untouched copy keeps its bytes. Before the fix the
// edited copy read Float = -0 in memory while 5000 was written, so the
// image writer (which scales with Float(BZERO)) wrote every pixel 5000 off;
// the comment-only copy was written "0", which cfitsio reads as +0.
func TestEditedNegZeroIntegerReadsBackAsHeld(t *testing.T) {
	lines := []string{
		"NZEDIT  =                   -0 / offset",
		"NZCMT   =                   -0 / offset",
		"NZKEEP  =                   -0 / offset",
	}
	from, added, path := editedHeaders(t, func(c *Card) {
		switch c.Key {
		case "NZEDIT":
			c.Value = int64(5000)
		case "NZCMT":
			c.Comment = "edited"
		}
	}, lines...)
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	for _, k := range []string{"NZEDIT", "NZCMT", "NZKEEP"} {
		checkNumberAsCfitsio(t, gh, k, from, added)
	}
	if v, err := from.Float("NZEDIT"); err != nil || v != 5000 {
		t.Errorf("Float(NZEDIT) = %v, %v; the card holds 5000", v, err)
	}
	if w := cfCard(t, gh, "NZKEEP"); string(w.Raw) != strings.TrimRight(lines[2], " ") {
		t.Errorf("untouched NZKEEP written as %q", w.Raw)
	}
	if w := cfCard(t, gh, "NZCMT"); string(w.Comment) != "edited" {
		t.Errorf("NZCMT comment %q", w.Comment)
	}
	if r := testutil.Fitsverify(t, path); !r.OK() {
		t.Errorf("fitsverify:\n%s", r.Report)
	}
}

// TestEditedFloatZeroSignReadsBackAsHeld: the encoder compares float
// values by their bits, so a decoded "-0.0" card edited to +0 (and a "0.0"
// card edited to -0) is written from its fields: -0 == +0 as numbers, and
// a numeric compare would keep the old bytes, which cfitsio reads with the
// other sign.
func TestEditedFloatZeroSignReadsBackAsHeld(t *testing.T) {
	from, added, path := editedHeaders(t, func(c *Card) {
		switch c.Key {
		case "FNZ":
			c.Value = 0.0
		case "FPZ":
			c.Value = math.Copysign(0, -1)
		}
	},
		"FNZ     =                 -0.0 / was -0",
		"FPZ     =                  0.0 / was +0",
	)
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	for _, k := range []string{"FNZ", "FPZ"} {
		checkNumberAsCfitsio(t, gh, k, from, added)
	}
	if f, _ := from.Float("FNZ"); math.Signbit(f) {
		t.Error("FNZ: the edit to +0 was lost in memory")
	}
	if f, _ := from.Float("FPZ"); !math.Signbit(f) {
		t.Error("FPZ: the edit to -0 was lost in memory")
	}
}

// TestEditedComplexPartReadsBackAsHeld: a decoded complex card whose
// imaginary part alone (or real part alone) was edited is written with
// the edit. Truth: cfitsio fits_read_key(TDBLCOMPLEX) of the written card.
func TestEditedComplexPartReadsBackAsHeld(t *testing.T) {
	want := map[string]Complex{"CPXIM": {Re: 1.5, Im: 9.5}, "CPXRE": {Re: 7.5, Im: 2.5}, "CPXINT": {Re: 3, Im: -8}}
	from, _, path := editedHeaders(t, func(c *Card) {
		if v, ok := want[c.Key]; ok {
			c.Value = v
		}
	},
		"CPXIM   =           (1.5, 2.5) / im edited",
		"CPXRE   =           (1.5, 2.5) / re edited",
		"CPXINT  =              (3, 4) / im edited",
	)
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	ls := lsHDU(t, path)
	for k, v := range want {
		c := ls.card(cfCard(t, gh, k).N)
		if c == nil || c.Cplx == nil || math.Float64bits(c.Cplx[0]) != math.Float64bits(v.Re) || math.Float64bits(c.Cplx[1]) != math.Float64bits(v.Im) {
			t.Errorf("%s: cfitsio reads %+v (card %q); the header holds %v", k, c, cfCard(t, gh, k).Raw, v)
		}
		if got, _ := from.Get(k); got.Value != v {
			t.Errorf("%s: header holds %v", k, got.Value)
		}
	}
}

// TestEditedKeyCaseIsWritten: a decoded lower-case keyword ("object", as
// cfitsio's ffgknm names it) whose copy had its Key upper-cased is written
// under the new key; keeping the old bytes would write the lower-case
// name, which fitsverify rejects. Truth: cfitsio's keyword name of the
// written card, and fitsverify on the source bytes (1 error) and on the
// written file (none).
func TestEditedKeyCaseIsWritten(t *testing.T) {
	srcPath, _ := linesFile(t, append(append([]string(nil), primaryLines...), "object  = 'M31' / lower"), "lower.fits")
	if r := testutil.Fitsverify(t, srcPath); r.Errors != 1 || !strings.Contains(r.Report, `"object"`) {
		t.Fatalf("fitsverify on the source bytes (the refusal assumes it rejects the lower-case name):\n%s", r.Report)
	}
	_, _, path := editedHeaders(t, func(c *Card) {
		if c.Key == "object" {
			c.Key = "OBJECT"
		}
	}, "object  = 'M31' / lower")
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	w := cfCard(t, gh, "OBJECT")
	if w.Str == nil || string(*w.Str) != "M31" || string(w.Comment) != "lower" {
		t.Errorf("cfitsio reads OBJECT %+v", w)
	}
	if gh.Card("object") != nil {
		t.Error("the copy renamed to OBJECT was written under its old lower-case name")
	}
	if r := testutil.Fitsverify(t, path); !r.OK() {
		t.Errorf("fitsverify:\n%s", r.Report)
	}
}

// TestLowerCaseKeyWrittenUpperCasedLikeCfitsio: a card built in memory or
// an edited copy of a decoded lower-case card is written with its keyword
// upper-cased by every path (EncodeCard, Encode of FromCards, AddCard), as
// Set upper-cases a name and as cfitsio's fits_write_key does: hdrref
// writekeys passes the name "lowkey" and cfitsio writes
// "LOWKEY  =                    5 / c", which must equal our card byte for
// byte. Before the fix Encode(FromCards) refused the edited copy ("invalid
// keyword ... 'o' is not allowed") while AddCard wrote it upper-cased.
func TestLowerCaseKeyWrittenUpperCasedLikeCfitsio(t *testing.T) {
	cfPath := filepath.Join(t.TempDir(), "cfkeys.fits")
	runHdrrefWriting(t, []string{cfPath}, "writekeys", cfPath)
	cf := fileCards(&testutil.CfrefDump(t, cfPath).HDUs[0])
	wantCard, ok := cf["LOWKEY"]
	if !ok {
		t.Fatalf("cfitsio did not write LOWKEY: %v", cf)
	}
	mem := Card{Key: "lowkey", Type: TypeInt, Value: int64(5), Comment: "c"}
	if enc, err := EncodeCard(mem); err != nil {
		t.Errorf("EncodeCard(lowkey): %v", err)
	} else if got := trimCard(enc); got != wantCard {
		t.Errorf("EncodeCard(lowkey):\n got %q\nwant %q (cfitsio fits_write_key)", got, wantCard)
	}
	for label, h := range map[string]*Header{"FromCards": FromCards([]Card{mem}), "AddCard": func() *Header {
		h := New()
		if err := h.AddCard(mem); err != nil {
			t.Fatal(err)
		}
		return h
	}()} {
		b, err := Encode(h)
		if err != nil {
			t.Errorf("%s: Encode: %v", label, err)
			continue
		}
		if got := strings.TrimRight(string(b[:CardWidth]), " "); got != wantCard {
			t.Errorf("%s: Encode wrote %q, cfitsio %q", label, got, wantCard)
		}
	}

	// An edited copy of a decoded lower-case card, both ways.
	from, _, path := editedHeaders(t, func(c *Card) {
		if c.Key == "target" {
			c.Value = "M33"
		}
	}, "target  = 'M31' / lower")
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	if w := cfCard(t, gh, "TARGET"); w.Str == nil || string(*w.Str) != "M33" || string(w.Comment) != "lower" {
		t.Errorf("cfitsio reads TARGET %+v", w)
	}
	if v, err := from.String("TARGET"); err != nil || v != "M33" {
		t.Errorf("String(TARGET) = %q, %v", v, err)
	}
	if r := testutil.Fitsverify(t, path); !r.OK() {
		t.Errorf("fitsverify:\n%s", r.Report)
	}
}

// TestEditedCommentaryCopyPastOneCard: a decoded HISTORY card whose copy
// was given 100 characters of text is written as two HISTORY cards of 72
// and 28 characters (the split cfitsio's fits_write_history makes, see
// TestEncodeReadByCfitsio), by Encode of FromCards and through AddCard.
// Truth: cfitsio's reading of the written cards.
func TestEditedCommentaryCopyPastOneCard(t *testing.T) {
	text := strings.Repeat("h", 100)
	_, _, path := editedHeaders(t, func(c *Card) {
		if c.Key == "HISTORY" {
			c.Comment = text
		}
	}, "HISTORY original", "KEEP    =                    1")
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	var got []string
	for _, c := range gh.Cards {
		if string(c.Key) == "HISTORY" {
			got = append(got, string(c.Comment))
		}
	}
	if len(got) != 2 || got[0] != text[:72] || got[1] != text[72:] {
		t.Errorf("cfitsio reads HISTORY cards %q; want the 100 characters as 72 + 28", got)
	}
	if gh.Card("KEEP") == nil {
		t.Error("the card after the edited HISTORY is missing")
	}
}

// TestEditedCopyWithMismatchedValueRefused: a copy of a decoded integer
// card whose Type was changed while its Value became a float64 that
// equals the integer bytes is an invalid card (a TypeString or
// TypeLogical card with a float64 value): it is refused by Encode,
// EncodeCardWithContinuation and AddCard, not written as its old bytes.
// The old bytes would not be the card: cfitsio reads them as the integer
// 5 (type I), while String / Bool of the header fail.
func TestEditedCopyWithMismatchedValueRefused(t *testing.T) {
	path, _ := linesFile(t, append(append([]string(nil), primaryLines...), "NUM     =                    5"), "num.fits")
	if w := cfCard(t, &testutil.CfrefDump(t, path).HDUs[0], "NUM"); w.Type != "I" || w.Int == nil || *w.Int != 5 {
		t.Fatalf("cfitsio reads NUM as %+v; the refusal assumes the integer 5", w)
	}
	for _, typ := range []ValueType{TypeString, TypeLogical, TypeInt, TypeComplexFloat, TypeEmpty} {
		c := decodeLines(t, "NUM     =                    5")[0]
		c.Type, c.Value = typ, 5.0
		if _, err := EncodeCardWithContinuation(c); err == nil {
			t.Errorf("%v with a float64 value: EncodeCardWithContinuation accepted it", typ)
		}
		if b, err := Encode(FromCards([]Card{c})); err == nil {
			t.Errorf("%v with a float64 value: Encode wrote\n%s", typ, cardsText(b))
		}
		if err := New().AddCard(c); err == nil {
			t.Errorf("%v with a float64 value: AddCard accepted it", typ)
		}
	}
}

// TestBigIntCardValueRefusedUnlessOutsideInt64: a TypeBigInt card must
// hold a non-nil *big.Int outside the int64 range, as the decoder makes
// it. A nil one (an edited copy of a decoded BIG card) would be written
// "<nil>", which cfitsio cannot read as a number (shown below on those
// bytes); one inside int64 would be written as an integer cfitsio reads
// with fits_read_key_lnglng while Int of the header in memory reports
// ErrIntOverflow. Both are refused by EncodeCardWithContinuation, Encode
// and AddCard, as Set refuses a nil *big.Int (and turns an in-range one
// into TypeInt).
func TestBigIntCardValueRefusedUnlessOutsideInt64(t *testing.T) {
	path, _ := linesFile(t, append(append([]string(nil), primaryLines...),
		"BIG     =                <nil> / b",
		"SMALL   =                    5 / b",
	), "nilbig.fits")
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	if w := cfCard(t, gh, "BIG"); w.Float != nil || w.Int != nil {
		t.Fatalf("cfitsio reads %q as a number (%+v); the refusal assumes it cannot", w.Raw, w)
	}
	if w := cfCard(t, gh, "SMALL"); w.Type != "I" || w.Int == nil || *w.Int != 5 {
		t.Fatalf("cfitsio reads SMALL as %+v", w)
	}
	if err := New().Set("BIG", (*big.Int)(nil), "b"); err == nil {
		t.Fatal("Set accepted a nil *big.Int")
	}
	if c, err := buildCard("SMALL", big.NewInt(5), "b"); err != nil || c.Type != TypeInt {
		t.Fatalf("Set of an in-range *big.Int: type %v, %v; want TypeInt", c.Type, err)
	}

	decoded := decodeLines(t, "BIG     = 99999999999999999999 / b")[0]
	nilCopy := decoded
	nilCopy.Value = (*big.Int)(nil)
	smallCopy := decoded
	smallCopy.Value = big.NewInt(5)
	for label, c := range map[string]Card{
		"edited copy, nil":      nilCopy,
		"edited copy, in range": smallCopy,
		"in memory, nil":        {Key: "BIG", Type: TypeBigInt, Value: (*big.Int)(nil)},
		"in memory, in range":   {Key: "BIG", Type: TypeBigInt, Value: big.NewInt(5)},
	} {
		if out, err := EncodeCardWithContinuation(c); err == nil {
			t.Errorf("%s: EncodeCardWithContinuation wrote %q", label, trimCard(out[0]))
		}
		if b, err := Encode(FromCards([]Card{c})); err == nil {
			t.Errorf("%s: Encode wrote\n%s", label, cardsText(b))
		}
		if err := New().AddCard(c); err == nil {
			t.Errorf("%s: AddCard accepted it", label)
		}
	}
	// The decoded card itself and an in-range edit made as TypeInt are fine.
	if _, err := Encode(FromCards([]Card{decoded})); err != nil {
		t.Errorf("unedited BIG: %v", err)
	}
	asInt := decoded
	asInt.Type, asInt.Value = TypeInt, int64(5)
	if _, err := Encode(FromCards([]Card{asInt})); err != nil {
		t.Errorf("BIG edited to TypeInt 5: %v", err)
	}
}
