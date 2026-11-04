package header

import (
	"bytes"
	"errors"
	"math"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Round-2 rules for keyword letter case, the END keyword and complex type
// switches. The truth in every test is an external reader of the bytes:
// cfitsio 4.7.0 (cfref dump, hdrref longstr), astropy and fitsverify.

// lowerSourceLines is a 3-pixel BITPIX 16 primary header with one extra
// card; the data block (zeros) follows in lowerSourceFile.
var lowerSourceLines = []string{
	"SIMPLE  =                    T",
	"BITPIX  =                   16",
	"NAXIS   =                    1",
	"NAXIS1  =                    3",
}

// lowerSourceFile writes lowerSourceLines + extra and a data unit whose
// stored int16 pixels are 1, 2, 3.
func lowerSourceFile(t *testing.T, extra string) string {
	t.Helper()
	_, hdr := linesFile(t, append(append([]string(nil), lowerSourceLines...), extra), "src.fits")
	data := make([]byte, 2880)
	for i, v := range []byte{1, 2, 3} {
		data[2*i+1] = v
	}
	return writeBytes(t, append(hdr, data...), "src.fits")
}

// astropyPixels returns astropy's scaled reading of the primary data.
func astropyPixels(t *testing.T, path string) string {
	t.Helper()
	out := testutil.AstropyCode(t, `
import sys, warnings
warnings.simplefilter("ignore")
from astropy.io import fits
with fits.open(sys.argv[1]) as f:
    print(" ".join(repr(float(v)) for v in f[0].data))
`, path)
	return strings.TrimSpace(string(out))
}

// TestAddCardRefusesUneditedLowerCaseKeyword: a decoded card whose
// (non-HIERARCH) keyword is not upper case is illegal FITS, and readers
// disagree on what it means, so it has no copy that every reader reads
// the same: on the source bytes fitsverify reports an error, cfitsio's
// HDU parser (ffgphd, case-sensitive) does not apply a "bzero" card to
// the pixels (fits_read_key still finds it by any case) while astropy
// applies it, and astropy takes "naxis1" for NAXIS1 while
// cfitsio's parser takes the NAXIS1 card. AddCard (which every writer
// uses to copy a caller's cards into the new HDU) therefore refuses it,
// with an error that says how to copy it (upper-case the copy's Key, so
// it is written from its fields, or drop it). Before the fix AddCard
// kept the card and the writers copied its bytes: a "bzero = 5000" card
// made WriteImage scale by 5000 while cfitsio read the file unscaled
// (pixels [-4999 -4998 -4997] for [1 2 3]), and a "naxis1" card was
// written as a second NAXIS1.
//
// A lower-case HIERARCH path is not refused: fitsverify accepts it and
// astropy and cfitsio both find it by any case.
func TestAddCardRefusesUneditedLowerCaseKeyword(t *testing.T) {
	// The disagreement the refusal rests on.
	bz := lowerSourceFile(t, "bzero   =                 5000 / lower")
	if r := testutil.Fitsverify(t, bz); r.Errors == 0 {
		t.Errorf("fitsverify accepts a lower-case bzero card:\n%s", r.Report)
	}
	cf := testutil.CfrefDump(t, bz).HDUs[0]
	if cf.Image == nil {
		t.Fatal("cfitsio sees no image")
	}
	// fits_read_key finds "bzero" by any case (cf.Image.Bzero is 5000),
	// but the pixel reader does not apply it: TDOUBLE reads 1, 2, 3.
	cfPix := []testutil.NullFloat{{V: 1, Valid: true}, {V: 2, Valid: true}, {V: 3, Valid: true}}
	if !cf.Image.Phys.Check(t, "cfitsio TDOUBLE of the source (bzero not applied)", cfPix) {
		t.FailNow()
	}
	if got, want := astropyPixels(t, bz), "5001.0 5002.0 5003.0"; got != want {
		t.Fatalf("astropy reads %q, want %q (it applies the lower-case bzero)", got, want)
	}

	for _, line := range []string{
		"bzero   =                 5000 / lower",
		"naxis1  =                    7 / lower",
		"bscale  =                    2",
		"blank   =                    2",
		"object  = 'M31'",
		"end      note here",
		"date.obs= '2020-01-01'",
	} {
		c := decodeLines(t, line)[0]
		h := New()
		err := h.AddCard(c)
		if !errors.Is(err, ErrInvalidKeyword) || !strings.Contains(err.Error(), "upper case") {
			t.Errorf("AddCard(%q) = %v, want an ErrInvalidKeyword about upper case", line, err)
		}
		if h.Len() != 0 {
			t.Errorf("AddCard(%q) failed but left %d cards", line, h.Len())
		}
		// The way out the error names: upper-case the copy's Key. It is
		// then written from its fields (a key Set cannot build stays an
		// error, at Encode).
		c.Key = strings.ToUpper(c.Key)
		if c.Key == KeyEnd || c.Key == "DATE.OBS" {
			continue
		}
		if err := h.AddCard(c); err != nil {
			t.Errorf("AddCard(%q upper-cased): %v", line, err)
		}
	}

	// Lower-case HIERARCH path: accepted, and both readers find it.
	line := "HIERARCH ESO det NAME = 'x' / c"
	path, _ := linesFile(t, append(append([]string(nil), primaryLines...), line), "hier.fits")
	if r := testutil.Fitsverify(t, path); !r.OK() {
		t.Fatalf("fitsverify on a lower-case HIERARCH path:\n%s", r.Report)
	}
	h := New()
	c := decodeLines(t, line)[0]
	if err := h.AddCard(c); err != nil {
		t.Fatalf("AddCard(%q): %v", line, err)
	}
	b, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(b[:CardWidth], c.Raw[:]) {
		t.Errorf("the HIERARCH card was not copied as its bytes: %q", b[:CardWidth])
	}
}

// TestFromCardsUpperCasesEditedKeys: FromCards stores the keyword of a
// card it will write from its fields (built in memory, or an edited copy
// of a decoded card) upper-cased, as AddCard stores it and as the card
// is written, so Cards(), Get and every caller that reads Key (the
// writers' structural-key filters) see the keyword the file will carry.
// A decoded card that is still its bytes keeps the decoded name. Before
// the fix FromCards kept an edited "naxis1" copy under Key "naxis1"
// while Encode wrote it "NAXIS1": WriteImage's filter missed it and
// wrote a second NAXIS1 (astropy then read 2 of 4 pixels). Truth: cfitsio
// reading the encoded header names the card as Cards() does.
func TestFromCardsUpperCasesEditedKeys(t *testing.T) {
	src := decodeLines(t, append(append([]string(nil), primaryLines...),
		"naxis1  =                    2 / lower",
		"object  = 'M31' / lower",
	)...)
	cards := append([]Card(nil), src...)
	for i := range cards {
		if cards[i].Key == "naxis1" {
			cards[i].Comment = "edited"
		}
	}
	cards = append(cards, Card{Key: "lowmem", Type: TypeInt, Value: int64(1)})
	h := FromCards(cards)
	var keys []string
	for _, c := range h.Cards() {
		keys = append(keys, c.Key)
	}
	want := []string{"SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "object", "LOWMEM"}
	if strings.Join(keys, ",") != strings.Join(want, ",") {
		t.Fatalf("FromCards keys %v, want %v", keys, want)
	}
	added := New()
	for _, c := range cards {
		if c.Key == "object" {
			continue // an unedited lower-case card: AddCard refuses it
		}
		if err := added.AddCard(c); err != nil {
			t.Fatal(err)
		}
	}
	if c, _ := added.Get("NAXIS1"); c.Key != "NAXIS1" {
		t.Errorf("AddCard stored the edited copy as %q", c.Key)
	}
	b, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	gh := &testutil.CfrefDump(t, writeBytes(t, b, "fc.fits")).HDUs[0]
	for i, c := range h.Cards() {
		if got := string(gh.Cards[i].Key); got != c.Key {
			t.Errorf("card %d: cfitsio names it %q, Cards() %q", i, got, c.Key)
		}
	}
	if w := cfCard(t, gh, "NAXIS1"); string(w.Comment) != "edited" {
		t.Errorf("cfitsio reads NAXIS1 %+v", w)
	}
}

// TestEndKeywordNeverWrittenFromFields: a card that would be written from
// its fields as an END card with anything else on it, or whose keyword
// only upper-cases to END, is refused by EncodeCard,
// EncodeCardWithContinuation, Encode and AddCard. cfitsio reads a
// lower-case "end" card as an ordinary card (the header goes on after it,
// shown on the source bytes) while the END card ends the header, and an
// END card has no comment field, so writing either from its fields would
// drop data. A bare Card{Key: "END"} still encodes as the END card. Before
// the fix EncodeCardWithContinuation wrote Card{Key: "end", Comment:
// "note"}, Card{Key: "END", Comment: "x"} and an edited copy of a decoded
// "end      note here" card as a bare END card, without an error.
func TestEndKeywordNeverWrittenFromFields(t *testing.T) {
	path, _ := linesFile(t, append(append([]string(nil), primaryLines...),
		"end      note here", "KEEP    =                    1"), "end.fits")
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	if gh.Card("KEEP") == nil || cfCard(t, gh, "end") == nil {
		t.Fatalf("cfitsio does not read the card after a lower-case end: %+v", gh.Cards)
	}

	decoded := decodeLines(t, "end      note here")[0]
	if decoded.Key != "end" {
		t.Fatalf("decoded key %q", decoded.Key)
	}
	edited := decoded
	edited.Comment = "edited"
	for name, c := range map[string]Card{
		"lower-case end":           {Key: "end", Comment: "note"},
		"lower-case end, no text":  {Key: "end"},
		"END with a comment":       {Key: "END", Comment: "x"},
		"edited decoded end":       edited,
		"decoded end, key renamed": func() Card { c := decoded; c.Key = "END"; return c }(),
	} {
		if out, err := EncodeCardWithContinuation(c); err == nil {
			t.Errorf("%s: EncodeCardWithContinuation wrote %q", name, out)
		}
		if out, err := EncodeCard(c); err == nil {
			t.Errorf("%s: EncodeCard wrote %q", name, out)
		}
		if b, err := Encode(FromCards([]Card{c})); err == nil {
			t.Errorf("%s: Encode wrote\n%s", name, cardsText(b))
		}
		if err := New().AddCard(c); err == nil {
			t.Errorf("%s: AddCard accepted it", name)
		}
	}
	end, err := EncodeCard(Card{Key: "END"})
	if err != nil || trimCard(end) != "END" {
		t.Errorf("EncodeCard(END) = %q, %v", end, err)
	}
	// The decoded card itself, unedited, is still its bytes for Encode.
	b, err := Encode(FromCards(decodeLines(t, append(append([]string(nil), primaryLines...), "end      note here")...)))
	if err != nil || !strings.Contains(cardsText(b), "end      note here") {
		t.Errorf("Encode of the unedited header: %v\n%s", err, cardsText(b))
	}
}

// TestEditedComplexTypeSwitchReadsBackAsHeld: a decoded complex card whose
// copy had its Type switched between TypeComplexInt and TypeComplexFloat
// is not its bytes any more. "CF = (1.5, 2.0)" switched to ComplexInt
// holds non-integer parts and is refused (cfitsio would read the old
// bytes as the float parts 1.5, 2.0 the Type says cannot be there);
// "CI = (1, 2)" switched to ComplexFloat is written as a float complex
// that reads back, through cfitsio (both parts, bit for bit) and this
// package, as the held value and type. A regression that let either copy
// keep its old bytes (mutant V17) wrote "(1.5, 2.0)" for the ComplexInt
// card and "(1, 2)", read back as ComplexInt, for the ComplexFloat one.
func TestEditedComplexTypeSwitchReadsBackAsHeld(t *testing.T) {
	cf := decodeLines(t, "CF      =           (1.5, 2.0) / c")[0]
	if cf.Type != TypeComplexFloat {
		t.Fatalf("CF decodes as %v", cf.Type)
	}
	cf.Type = TypeComplexInt
	if out, err := EncodeCardWithContinuation(cf); err == nil {
		t.Errorf("ComplexInt (1.5,2): EncodeCardWithContinuation wrote %q", out)
	}
	if b, err := Encode(FromCards([]Card{cf})); err == nil {
		t.Errorf("ComplexInt (1.5,2): Encode wrote\n%s", cardsText(b))
	}
	if err := New().AddCard(cf); err == nil {
		t.Error("ComplexInt (1.5,2): AddCard accepted it")
	}

	from, _, path := editedHeaders(t, func(c *Card) {
		if c.Key == "CI" {
			c.Type = TypeComplexFloat
		}
	}, "CI      =               (1, 2) / c")
	held, _ := from.Get("CI")
	if held.Type != TypeComplexFloat {
		t.Fatalf("held type %v", held.Type)
	}
	gh := &testutil.CfrefDump(t, path).HDUs[0]
	w := cfCard(t, gh, "CI")
	ls := lsHDU(t, path).card(w.N)
	if ls == nil || ls.Cplx == nil || math.Float64bits(ls.Cplx[0]) != math.Float64bits(1) || math.Float64bits(ls.Cplx[1]) != math.Float64bits(2) {
		t.Errorf("cfitsio reads CI %+v (card %q)", ls, w.Raw)
	}
	var raw [CardWidth]byte
	copy(raw[:], padCard(string(w.Raw)))
	back := decodeCard(&raw)
	if back.Type != TypeComplexFloat || back.Value != held.Value {
		t.Errorf("the written card %q reads back as %v %v; the header holds %v %v", w.Raw, back.Type, back.Value, held.Type, held.Value)
	}
}
