package header

import (
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// encodeCase is one keyword the Go encoder writes and cfitsio / astropy
// must read back with exactly this value and comment.
type encodeCase struct {
	key     string // as passed to Set ("HIERARCH ..." for HIERARCH)
	value   any
	comment string
}

func rep(p string, n int) string {
	var b strings.Builder
	for b.Len() < n {
		b.WriteString(p)
	}
	return b.String()[:n]
}

// encodeCases covers long strings around every CONTINUE boundary, quote
// escaping across card boundaries, comments too long for one card,
// HIERARCH long strings and every number form.
func encodeCases() []encodeCase {
	cs := []encodeCase{
		{"LONGSTR", rep("0123456789abcdef'ghij", 200), rep("comment text ", 119) + "end"},
		{"QEDGE", rep("a", 66) + "'" + rep("b", 80), "quote at the first chunk boundary"},
		{"QEDGE2", rep("a", 65) + "''" + rep("c", 70), ""},
		{"ALLQUOTE", rep("'", 100), "only quotes"},
		{"LEADSP", "   leading blanks", "kept"},
		{"MIDSP", rep("x", 66) + "      " + rep("y", 10), "blanks at a chunk end"},
		{"EMPTY", "", "empty string"},
		{"EMPTYLC", "", rep("long comment on an empty string ", 4) + "end"},
		{"SHORTLC", "abc", rep("long comment on a short string ", 5) + "end"},
		{"AMPEND", "ends with &", "a literal ampersand"},
		{"HIERARCH ESO OBS LONGNAME", rep("hierarch-value ", 150) + "end", "hierarch comment"},
		{"HIERARCH ESO DET NAME", "short", "hierarch short"},
		{"HIERARCH ESO VERY LONG KEYWORD PATH NAME FOR TESTING", rep("v", 120), "c"},
		{"HIERARCH ESO TEL ALT", 45.5, "hierarch float"},
		{"FLT120", 120.0, "integral double"},
		{"FINT3", 3.0, ""},
		{"FSMALL", 1e-7, "small"},
		{"FBIG", 1.5e300, "big"},
		{"FTENTH", 0.1, "0.1"},
		{"FNEGZERO", math.Copysign(0, -1), "negative zero"},
		{"FMAX", math.MaxFloat64, "largest float64"},
		{"FMINNORM", 2.2250738585072014e-308, "smallest normal float64"},
		{"F1E21", 1e21, ""},
		{"F123456", 123456.0, ""},
		{"F32", float32(0.1), "float32 widened exactly"},
		{"I64MIN", int64(math.MinInt64), "int64 min"},
		{"I64MAX", int64(math.MaxInt64), "int64 max"},
		{"U64MAX", uint64(math.MaxUint64), "2**64-1"},
		{"U63", uint64(1) << 63, "2**63 (uint64 BZERO)"},
		{"BIGNEG", new(big.Int).Sub(big.NewInt(math.MinInt64), big.NewInt(1)), "int64 min - 1"},
		{"LOGT", true, "true"},
		{"LOGF", false, ""},
		{"CPX", Complex{Re: 1.5, Im: -2.25}, "complex"},
		{"UNDEF", nil, "undefined value"},
		// Trailing blanks on a string that needs CONTINUE cards: they are
		// not significant, but a last piece of only blanks with no comment
		// would be an empty CONTINUE, which ends the chain and leaves the
		// '&' in the value.
		{"TRAILBL", rep("x", 67) + "     ", ""},
		{"TRAILBL2", rep("x", 134) + "  ", ""},
		{"TRAILBLC", rep("x", 67) + "     ", "c"},
		{"BLANKS", rep(" ", 100), ""},
	}
	// Every string length across the first two CONTINUE boundaries, with
	// and without a comment.
	for n := 60; n <= 140; n += 1 {
		cs = append(cs, encodeCase{fmt.Sprintf("S%03d", n), rep("z", n), ""})
		cs = append(cs, encodeCase{fmt.Sprintf("C%03d", n), rep("w", n), "c"})
		cs = append(cs, encodeCase{fmt.Sprintf("T%03d", n), rep("t", n) + "   ", ""})
	}
	return cs
}

var (
	longHistory = rep("history words ", 199) + "end"
	longComment = rep("comment words ", 149) + "end"
	longBlank   = rep("blank keyword text ", 6) + "end"
)

// writeEncodedHeader builds the encodeCases header with the Go encoder
// and writes it as a NAXIS=0 primary HDU.
func writeEncodedHeader(t *testing.T) (string, *Header) {
	t.Helper()
	h := New()
	mustSet(t, h, "SIMPLE", true, "conforms to FITS standard")
	mustSet(t, h, "BITPIX", int64(8), "")
	mustSet(t, h, "NAXIS", int64(0), "")
	// The OGIP long-string convention recommends announcing itself.
	mustSet(t, h, "LONGSTRN", "OGIP 1.0", "The HEASARC Long String Convention may be used.")
	for _, c := range encodeCases() {
		mustSet(t, h, c.key, c.value, c.comment)
	}
	if err := h.Add("HISTORY", longHistory, ""); err != nil {
		t.Fatal(err)
	}
	if err := h.Add("COMMENT", longComment, ""); err != nil {
		t.Fatal(err)
	}
	if err := h.Add("", longBlank, ""); err != nil {
		t.Fatal(err)
	}
	b, err := Encode(h)
	if err != nil {
		t.Fatalf("Encode: %v", err)
	}
	p := filepath.Join(t.TempDir(), "goenc.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return p, h
}

func mustSet(t *testing.T, h *Header, k string, v any, c string) {
	t.Helper()
	if err := h.Set(k, v, c); err != nil {
		t.Fatalf("Set(%s): %v", k, err)
	}
}

// chunks72 is how cfitsio (fits_write_history / fits_write_comment) and
// this encoder split long commentary text; readers drop each card's
// trailing blanks.
func chunks72(s string) []string {
	var out []string
	for len(s) > 72 {
		out = append(out, strings.TrimRight(s[:72], " "))
		s = s[72:]
	}
	return append(out, strings.TrimRight(s, " "))
}

func cfitsioKey(k string) string {
	if strings.HasPrefix(k, "HIERARCH ") {
		return strings.TrimPrefix(k, "HIERARCH ")
	}
	return k
}

// TestEncodeReadByCfitsio writes every encodeCase with the Go encoder and
// requires cfitsio to read back exactly the value and comment given, and
// fitsverify to find no warning or error.
func TestEncodeReadByCfitsio(t *testing.T) {
	path, _ := writeEncodedHeader(t)
	cf := testutil.CfrefDump(t, path)
	ls := hdrrefLongstr(t, path)
	if len(cf.HDUs) != 1 || len(cf.HDUs[0].Errors) != 0 {
		t.Fatalf("cfitsio: %d HDUs, errors %v", len(cf.HDUs), cf.HDUs[0].Errors)
	}
	gh, lh := &cf.HDUs[0], &ls.HDUs[0]
	for _, c := range encodeCases() {
		k := cfitsioKey(c.key)
		w := gh.Card(k)
		if w == nil {
			t.Errorf("%s: cfitsio does not see the keyword", c.key)
			continue
		}
		l := lh.card(w.N)
		switch v := c.value.(type) {
		case string:
			if w.Str == nil || string(*w.Str) != strings.TrimRight(v, " ") {
				t.Errorf("%s: cfitsio string %v (%s), want %q", c.key, w.Str, w.ReadError, v)
			}
			wantC := c.comment
			if len(wantC) > 72 {
				wantC = wantC[:72] // fits_read_key_longstr truncates the joined comment
			}
			if l == nil || l.Comment == nil || string(*l.Comment) != wantC {
				t.Errorf("%s: cfitsio longstr comment %v, want %q", c.key, l, wantC)
			}
			continue
		case float64, float32:
			f := toF64(v)
			if w.Type != "F" || w.Float == nil || math.Float64bits(float64(*w.Float)) != math.Float64bits(f) {
				t.Errorf("%s: cfitsio type %q float %v (%s), want F %v (raw %q)", c.key, w.Type, w.Float, w.ReadError, f, w.Raw)
			}
		case int64:
			if w.Type != "I" || w.Int == nil || *w.Int != v {
				t.Errorf("%s: cfitsio type %q int %v (%s), want %d", c.key, w.Type, w.Int, w.ReadError, v)
			}
		case uint64:
			if l == nil || l.U64 == nil || *l.U64 != v {
				t.Errorf("%s: cfitsio TULONGLONG %+v, want %d (raw %q)", c.key, l, v, w.Raw)
			}
		case *big.Int:
			if string(w.Value) != v.String() || w.Float == nil {
				t.Errorf("%s: cfitsio value %q float %v, want %s", c.key, w.Value, w.Float, v)
			}
		case bool:
			if w.Logical == nil || *w.Logical != v {
				t.Errorf("%s: cfitsio logical %v, want %v", c.key, w.Logical, v)
			}
		case Complex:
			if l == nil || l.Cplx == nil || l.Cplx[0] != v.Re || l.Cplx[1] != v.Im {
				t.Errorf("%s: cfitsio complex %+v, want %v", c.key, l, v)
			}
		case nil:
			if w.Value != "" || w.Type != "" {
				t.Errorf("%s: cfitsio value %q type %q, want undefined", c.key, w.Value, w.Type)
			}
		default:
			t.Fatalf("%s: no check for %T", c.key, v)
		}
		if string(w.Comment) != c.comment {
			t.Errorf("%s: cfitsio comment %q, want %q", c.key, w.Comment, c.comment)
		}
	}
	// Long commentary: split into 72-character cards, in order.
	for _, cc := range []struct{ key, text string }{{"HISTORY", longHistory}, {"COMMENT", longComment}, {"", longBlank}} {
		var got []string
		for _, card := range gh.Cards {
			if string(card.Key) == cc.key && card.Value == "" && strings.HasPrefix(string(card.Raw), padKey(cc.key)) {
				got = append(got, string(card.Comment))
			}
		}
		want := chunks72(cc.text)
		if strings.Join(got, "\n") != strings.Join(want, "\n") {
			t.Errorf("%q cards as cfitsio reads them:\n got %q\nwant %q", cc.key, got, want)
		}
	}
	// The Go decoder reads its own output exactly as cfitsio does.
	b, _ := os.ReadFile(path)
	cards, _, err := ParseCards(b)
	if err != nil {
		t.Fatal(err)
	}
	compareWithCfitsio(t, "goenc.fits", cards, gh, lh)

	// fitsverify: no errors; the one warning is the deliberate undefined
	// value of UNDEF.
	r := testutil.Fitsverify(t, path)
	if r.Errors != 0 || r.Warnings != 1 || !strings.Contains(r.Report, "UNDEF has a null value") {
		t.Errorf("fitsverify: %d warnings, %d errors (want only the UNDEF null-value warning)\n%s", r.Warnings, r.Errors, r.Report)
	}
}

// TestEncodeStringsReadBackByGo: every encodeCase string, written by the
// encoder and parsed again by this package, reads back as the string given
// without its trailing blanks (FITS §4.2.1.1), with the full comment. It
// needs no external tool; TestEncodeReadByCfitsio and
// TestEncodeReadByAstropy check the same file with cfitsio and astropy.
func TestEncodeStringsReadBackByGo(t *testing.T) {
	path, _ := writeEncodedHeader(t)
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	cards, _, err := ParseCards(b)
	if err != nil {
		t.Fatal(err)
	}
	h := FromCards(cards)
	for _, c := range encodeCases() {
		v, ok := c.value.(string)
		if !ok {
			continue
		}
		got, err := h.String(c.key)
		if err != nil || got != strings.TrimRight(v, " ") {
			t.Errorf("%s: read back %q, %v; want %q", c.key, got, err, strings.TrimRight(v, " "))
		}
		if card, ok := h.Get(c.key); !ok || card.Comment != c.comment {
			t.Errorf("%s: comment %q (found %v), want %q", c.key, card.Comment, ok, c.comment)
		}
	}
}

func toF64(v any) float64 {
	switch x := v.(type) {
	case float32:
		return float64(x)
	case float64:
		return x
	}
	panic("not a float")
}

const astropyDumpScript = `
import json, sys, warnings
from astropy.io import fits
warnings.simplefilter("error")
h = fits.getheader(sys.argv[1])
out = {"cards": {}, "HISTORY": [str(x) for x in h["HISTORY"]], "COMMENT": [str(x) for x in h["COMMENT"]], "BLANK": [str(x) for x in h[""]]}
for k in h.keys():
    if k in ("HISTORY", "COMMENT", ""):
        continue
    v = h[k]
    if isinstance(v, bool):
        e = {"t": "bool", "v": v}
    elif isinstance(v, int):
        e = {"t": "int", "v": str(v)}
    elif isinstance(v, float):
        e = {"t": "float", "v": v.hex()}
    elif isinstance(v, complex):
        e = {"t": "complex", "v": [v.real.hex(), v.imag.hex()]}
    elif isinstance(v, str):
        e = {"t": "str", "v": v}
    elif v is None or type(v).__name__ == "Undefined":
        e = {"t": "undef"}
    else:
        e = {"t": type(v).__name__, "v": repr(v)}
    e["c"] = h.comments[k]
    out["cards"][k] = e
json.dump(out, sys.stdout)
`

type astropyCard struct {
	T string          `json:"t"`
	V json.RawMessage `json:"v"`
	C string          `json:"c"`
}

// TestEncodeReadByAstropy is the second, independent reader of the Go
// encoder's output: astropy (with warnings turned into errors) must read
// every value and the FULL comment back exactly.
func TestEncodeReadByAstropy(t *testing.T) {
	path, _ := writeEncodedHeader(t)
	out := testutil.AstropyCode(t, astropyDumpScript, path)
	var got struct {
		Cards   map[string]astropyCard `json:"cards"`
		History []string               `json:"HISTORY"`
		Comment []string               `json:"COMMENT"`
		Blank   []string               `json:"BLANK"`
	}
	if err := json.Unmarshal(out, &got); err != nil {
		t.Fatalf("astropy output: %v\n%s", err, out)
	}
	for _, c := range encodeCases() {
		k := cfitsioKey(c.key)
		a, ok := got.Cards[k]
		if !ok {
			t.Errorf("%s: astropy does not see the keyword", c.key)
			continue
		}
		if a.C != c.comment {
			t.Errorf("%s: astropy comment %q, want %q", c.key, a.C, c.comment)
		}
		var sv string
		_ = json.Unmarshal(a.V, &sv)
		switch v := c.value.(type) {
		case string:
			if a.T != "str" || sv != strings.TrimRight(v, " ") {
				t.Errorf("%s: astropy %s %q, want %q", c.key, a.T, sv, v)
			}
		case float64, float32:
			f, err := strconv.ParseFloat(sv, 64)
			if a.T != "float" || err != nil || math.Float64bits(f) != math.Float64bits(toF64(v)) {
				t.Errorf("%s: astropy %s %s, want %v", c.key, a.T, a.V, v)
			}
		case int64, uint64, *big.Int:
			if a.T != "int" || sv != fmt.Sprint(v) {
				t.Errorf("%s: astropy %s %s, want %v", c.key, a.T, a.V, v)
			}
		case bool:
			var b bool
			_ = json.Unmarshal(a.V, &b)
			if a.T != "bool" || b != v {
				t.Errorf("%s: astropy %s %s, want %v", c.key, a.T, a.V, v)
			}
		case Complex:
			var p [2]string
			_ = json.Unmarshal(a.V, &p)
			re, _ := strconv.ParseFloat(p[0], 64)
			im, _ := strconv.ParseFloat(p[1], 64)
			if a.T != "complex" || re != v.Re || im != v.Im {
				t.Errorf("%s: astropy %s %s, want %v", c.key, a.T, a.V, v)
			}
		case nil:
			if a.T != "undef" {
				t.Errorf("%s: astropy %s %s, want undefined", c.key, a.T, a.V)
			}
		}
	}
	for _, cc := range []struct {
		name      string
		got, want []string
	}{{"HISTORY", got.History, chunks72(longHistory)}, {"COMMENT", got.Comment, chunks72(longComment)}, {"blank", got.Blank, chunks72(longBlank)}} {
		if strings.Join(cc.got, "\n") != strings.Join(cc.want, "\n") {
			t.Errorf("astropy %s:\n got %q\nwant %q", cc.name, cc.got, cc.want)
		}
	}
}

// TestEncodeFloatAlwaysFloatSyntax: every float64 the encoder writes is
// typed 'F' by cfitsio (it carries '.' or 'E'), so a float that happens to
// be integral is not re-read as an integer. Checked on the encoded text and
// against cfitsio's typing in TestEncodeReadByCfitsio.
func TestEncodeFloatAlwaysFloatSyntax(t *testing.T) {
	for _, f := range []float64{0, 1, -1, 120, 3, 1e6, 1e15, 1e16, 1e21, 123456, 2.5, math.MaxFloat64, -math.MaxFloat64, 2.2250738585072014e-308, 1e-7} {
		s, err := formatValue(Card{Type: TypeFloat, Value: f})
		if err != nil {
			t.Fatalf("%v: %v", f, err)
		}
		if !strings.ContainsAny(s, ".E") {
			t.Errorf("float %v encoded as %q: no '.' or 'E', cfitsio would type it as an integer", f, s)
		}
		back, err := parseFITSFloat(s)
		if err != nil || math.Float64bits(back) != math.Float64bits(f) {
			t.Errorf("float %v encoded as %q does not round-trip: %v %v", f, s, back, err)
		}
	}
}

// TestEncodeRejectsWhatFITSCannotHold: non-finite floats, keywords that
// are not legal FITS keywords, and text outside printable ASCII return an
// error instead of writing an invalid or truncated card.
func TestEncodeRejectsWhatFITSCannotHold(t *testing.T) {
	bad := []struct {
		name    string
		key     string
		value   any
		comment string
	}{
		{"NaN", "FNAN", math.NaN(), ""},
		{"+Inf", "FINF", math.Inf(1), ""},
		{"-Inf", "FNINF", math.Inf(-1), ""},
		{"float32 NaN", "F32NAN", float32(math.NaN()), ""},
		{"complex NaN", "CNAN", Complex{Re: math.NaN()}, ""},
		{"complex Inf", "CINF", Complex{Im: math.Inf(1)}, ""},
		{"9-char key", "TOOLONGKY", 1, ""},
		{"long key", "THISKEYWORDISLONG", "x", ""},
		{"non-ASCII string", "UTF8", "caf\u00e9", ""},
		{"control char in string", "CTRL", "a\tb", ""},
		{"non-ASCII comment", "CMT8", 1, "\u00e9"},
		{"newline in HISTORY", "HISTORY", "line1\nline2", ""},
		{"comment too long for a number card", "NUMCMT", 1, rep("c", 60)},
		{"END keyword", "END", nil, ""},
		{"HIERARCH with '='", "HIERARCH A=B", 1, ""},
		{"empty HIERARCH path", "HIERARCH ", 1, ""},
		{"HIERARCH key too long for any value", "HIERARCH " + rep("K", 72), 1, ""},
		{"unsupported Go type", "MAP", map[string]int{}, ""},
		// A string that needs CONTINUE cards and ends in '&': cfitsio keeps
		// the '&', astropy drops it (seen with the AMPTRAIL case).
		{"CONTINUE'd string ending in '&'", "AMPLONG", rep("a", 70) + "&", ""},
		{"CONTINUE'd string ending in '&' and blanks", "AMPTRAIL", rep("a", 70) + "&  ", ""},
		{"string ending in '&' with a CONTINUE'd comment", "AMPCMT", "short&", rep("comment ", 96) + "end"},
	}
	// Set validates the keyword, the value and the text itself; only the
	// layout (does it fit, can it be continued) is left to Encode.
	encodeOnly := map[string]bool{"NUMCMT": true, "AMPLONG": true, "AMPTRAIL": true, "AMPCMT": true}
	for _, b := range bad {
		h := New()
		err := h.Set(b.key, b.value, b.comment)
		if !encodeOnly[b.key] {
			if err == nil {
				t.Errorf("%s: Set(%q, %#v) accepted it", b.name, b.key, b.value)
			}
			continue
		}
		if err != nil {
			t.Errorf("%s: Set(%q) refused a value it should leave to Encode: %v", b.name, b.key, err)
			continue
		}
		if enc, eerr := Encode(h); eerr == nil {
			t.Errorf("%s: Encode accepted Set(%q, %#v); header:\n%s", b.name, b.key, b.value, cardsText(enc))
		}
	}
	// Each remaining guard of the encoder and card validation, with the
	// error it must give (each could be removed before without a failing
	// test; with the ComplexInt guard off, (1.5, 2) was written as (1, 2)).
	for _, tc := range []struct {
		name string
		c    Card
		want string
	}{
		{"integer complex with a fractional part", Card{Key: "CINT", Type: TypeComplexInt, Value: Complex{Re: 1.5, Im: 2}}, "non-integer parts"},
		{"integer complex beyond 2^53", Card{Key: "CINT", Type: TypeComplexInt, Value: Complex{Re: 1 << 54, Im: 0}}, "non-integer parts"},
		{"HIERARCH key with no room for a string", Card{Key: "HIERARCH " + rep("K", 66), Type: TypeString, Value: "abc"}, "leaves no room for a string value"},
		{"END with a value", Card{Key: "END", Type: TypeInt, Value: int64(1)}, "END cannot carry a value"},
		{"COMMENT with a typed value", Card{Key: "COMMENT", Type: TypeInt, Value: int64(1)}, "cannot carry a typed value"},
		{"HISTORY with a typed value", Card{Key: "HISTORY", Type: TypeLogical, Value: true}, "cannot carry a typed value"},
		{"TypeBigInt holding int64", Card{Key: "BIG", Type: TypeBigInt, Value: int64(1)}, "TypeBigInt with value int64"},
		{"TypeInt holding float64", Card{Key: "INT", Type: TypeInt, Value: 1.5}, "TypeInt with value float64"},
		{"TypeFloat holding string", Card{Key: "FLT", Type: TypeFloat, Value: "1.5"}, "TypeFloat with value string"},
		{"TypeLogical holding int64", Card{Key: "LOG", Type: TypeLogical, Value: int64(1)}, "TypeLogical with value int64"},
		{"TypeEmpty holding a value", Card{Key: "EMP", Type: TypeEmpty, Value: int64(1)}, "TypeEmpty with value int64"},
	} {
		_, err := EncodeCardWithContinuation(tc.c)
		if err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%s: EncodeCardWithContinuation err = %v, want one containing %q", tc.name, err, tc.want)
		}
		if _, err := Encode(FromCards([]Card{tc.c})); err == nil {
			t.Errorf("%s: Encode accepted %+v", tc.name, tc.c)
		}
	}
	// EncodeCard is for one physical card: a string that needs CONTINUE
	// cards is an error there (EncodeCardWithContinuation writes it).
	long := Card{Key: "LONGSTR", Type: TypeString, Value: rep("x", 100)}
	if _, err := EncodeCard(long); err == nil || !strings.Contains(err.Error(), "needs 2 cards; use EncodeCardWithContinuation") {
		t.Errorf("EncodeCard on a 100-character string: err = %v", err)
	}
	if out, err := EncodeCardWithContinuation(long); err != nil || len(out) != 2 {
		t.Errorf("EncodeCardWithContinuation on a 100-character string: %d cards, %v", len(out), err)
	}
	// Bool on a card that is not logical is ErrWrongType, not a panic or false.
	hb := New()
	if err := hb.Set("NOTBOOL", int64(1), ""); err != nil {
		t.Fatal(err)
	}
	if _, err := hb.Bool("NOTBOOL"); !errors.Is(err, ErrWrongType) {
		t.Errorf("Bool on an integer card: err = %v, want ErrWrongType", err)
	}
	// Cards built directly (not through Set) are validated by Encode. (A
	// lower-case name is written upper-cased, as cfitsio writes it: see
	// TestLowerCaseKeyWrittenUpperCasedLikeCfitsio.)
	for _, c := range []Card{
		{Key: "FNAN", Type: TypeFloat, Value: math.NaN()},
		{Key: "LONGKEYWORD", Type: TypeInt, Value: int64(1)},
		{Key: "BAD KEY", Type: TypeInt, Value: int64(1)},
		{Key: "HISTORY", Comment: "\x00"},
	} {
		if _, err := Encode(FromCards([]Card{c})); err == nil {
			t.Errorf("Encode accepted in-memory card %+v", c)
		}
	}
}

func cardsText(b []byte) string {
	var s []string
	for i := 0; i+CardWidth <= len(b); i += CardWidth {
		if l := strings.TrimRight(string(b[i:i+CardWidth]), " "); l != "" {
			s = append(s, l)
		}
	}
	return strings.Join(s, "\n")
}

// TestKeywordValidationMatchesCfitsio: for every name in
// golden/testkey.json (at most 8 characters), Set accepts it exactly when
// cfitsio's fits_test_keyword does; names longer than 8 characters are
// rejected unless written with the explicit "HIERARCH " prefix (cfitsio's
// writer would silently switch them to HIERARCH; this package does not).
func TestKeywordValidationMatchesCfitsio(t *testing.T) {
	p := headerTestdata(t, "golden", "testkey.json")
	testutil.RequireFile(t, p)
	b, _ := os.ReadFile(p)
	var keys []struct {
		Key    testutil.Latin1 `json:"key"`
		Status int             `json:"status"`
	}
	if err := json.Unmarshal(b, &keys); err != nil {
		t.Fatal(err)
	}
	if len(keys) < 10 {
		t.Fatalf("testkey golden has %d keys", len(keys))
	}
	for _, k := range keys {
		name := string(k.Key)
		if len(name) > 8 {
			t.Fatalf("testkey golden key %q longer than 8", name)
		}
		err := New().Set(name, int64(1), "")
		if (err == nil) != (k.Status == 0) {
			t.Errorf("Set(%q) err = %v; cfitsio fits_test_keyword status %d", name, err, k.Status)
		}
		// The documented sentinel (README, ErrInvalidKeyword doc).
		if k.Status != 0 && err != nil && !errors.Is(err, ErrInvalidKeyword) {
			t.Errorf("Set(%q) err = %v, want errors.Is ErrInvalidKeyword", name, err)
		}
	}
	for _, name := range []string{"NAXIS1234", "TOOLONGKEY", "DATE-OBS-X"} {
		if err := New().Set(name, int64(1), ""); !errors.Is(err, ErrInvalidKeyword) {
			t.Errorf("Set(%q): err %v, want ErrInvalidKeyword for a %d-character keyword without HIERARCH", name, err, len(name))
		}
		if err := New().Set("HIERARCH "+name, int64(1), ""); err != nil {
			t.Errorf("Set(HIERARCH %s): %v", name, err)
		}
	}
}

// TestNonASCIIIsInvalidFITS backs the printable-ASCII rule with fitsverify:
// a card holding a byte >= 0x80 in a string value is reported as an error.
func TestNonASCIIIsInvalidFITS(t *testing.T) {
	h := New()
	mustSet(t, h, "SIMPLE", true, "")
	mustSet(t, h, "BITPIX", int64(8), "")
	mustSet(t, h, "NAXIS", int64(0), "")
	mustSet(t, h, "OBJECT", "caf_", "")
	b, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	i := strings.Index(string(b), "caf_")
	b[i+3] = 0xE9 // Latin-1 e-acute, written behind the encoder's back
	p := filepath.Join(t.TempDir(), "nonascii.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	if r := testutil.Fitsverify(t, p); r.Errors == 0 {
		t.Fatalf("fitsverify accepts a non-ASCII string value:\n%s", r.Report)
	}
	if err := New().Set("OBJECT", "caf\xe9", ""); !errors.Is(err, ErrInvalidText) {
		t.Errorf("Set of the same non-ASCII string: err %v, want ErrInvalidText", err)
	}
	if err := New().Set("OBJECT", "cafe", "caf\xe9"); !errors.Is(err, ErrInvalidText) {
		t.Errorf("Set of a non-ASCII comment: err %v, want ErrInvalidText", err)
	}
}

// TestCommentWordTooLongToContinue: a comment that must be continued but
// has no single blank to break at within a card is refused rather than
// split mid-word (readers would rejoin it with an extra blank).
func TestCommentWordTooLongToContinue(t *testing.T) {
	h := New()
	if err := h.Set("LONGSTR", rep("v", 100), rep("w", 90)); err == nil {
		if _, err := Encode(h); err == nil {
			t.Fatal("encoded a 90-character comment word across CONTINUE cards")
		}
	}
	if err := h.Set("LONGSTR", rep("v", 100), rep("w", 60)+" "+rep("u", 60)); err != nil {
		t.Fatal(err)
	}
	if _, err := Encode(h); err != nil {
		t.Fatalf("a comment with a blank to break at: %v", err)
	}
}

// TestIntOverflowGetters: an in-memory uint64 above MaxInt64 is kept
// exactly (TypeBigInt), Int refuses it, Float and BigInt give it.
func TestIntOverflowGetters(t *testing.T) {
	h := New()
	mustSet(t, h, "U64MAX", uint64(math.MaxUint64), "")
	if _, err := h.Int("U64MAX"); !errors.Is(err, ErrIntOverflow) {
		t.Errorf("Int: %v, want ErrIntOverflow", err)
	}
	if b, err := h.BigInt("U64MAX"); err != nil || !b.IsUint64() || b.Uint64() != math.MaxUint64 {
		t.Errorf("BigInt: %v %v", b, err)
	}
	if f, err := h.Float("U64MAX"); err != nil || f != 18446744073709551615.0 {
		t.Errorf("Float: %v %v", f, err)
	}
}
