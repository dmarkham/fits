package header

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// This file holds the helpers that compare the header package against
// cfitsio's reading of the same bytes. Two harnesses provide the truth:
//
//   - testdata/cfref (shared): "cfref dump" lists every physical card with
//     cfitsio's keyword name, value string, comment and typed reads.
//   - header/testdata/hdrref: "hdrref longstr" adds what cfref does not
//     record: the comment fits_read_key_longstr returns for a CONTINUE'd
//     string, uint64 reads and complex reads.
//
// Committed goldens are read with LoadCfrefGolden / loadLongstrGolden;
// files written by the Go encoder at test time go through testutil.Run:
// the recorded reading of those exact bytes in tier 1, the built
// binaries live in tier 2 (-tags interop).

// longstrFile is one "hdrref longstr" result.
type longstrFile struct {
	File    string       `json:"file"`
	Cfitsio string       `json:"cfitsio"`
	HDUs    []longstrHDU `json:"hdus"`
}

type longstrHDU struct {
	Index int           `json:"index"`
	Keys  []longstrCard `json:"keys"`
}

type longstrCard struct {
	N       int              `json:"n"` // 1-based physical card number
	Key     testutil.Latin1  `json:"key"`
	Str     *testutil.Latin1 `json:"str"`
	Comment *testutil.Latin1 `json:"comment"`
	U64     *uint64          `json:"u64"`
	Cplx    *[2]float64      `json:"cplx"`
	Skipped string           `json:"skipped"`
	Error   string           `json:"error"`
}

func (h *longstrHDU) card(n int) *longstrCard {
	for i := range h.Keys {
		if h.Keys[i].N == n {
			return &h.Keys[i]
		}
	}
	return nil
}

func decodeLongstr(b []byte) (*longstrFile, error) {
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	var f longstrFile
	if err := dec.Decode(&f); err != nil {
		return nil, err
	}
	if f.Cfitsio != "4.7.0" {
		return nil, fmt.Errorf("longstr golden made with cfitsio %q, want 4.7.0", f.Cfitsio)
	}
	return &f, nil
}

// headerTestdata is <repo>/header/testdata.
func headerTestdata(t testing.TB, elem ...string) string {
	t.Helper()
	return filepath.Join(append([]string{testutil.RepoRoot(t), "header", "testdata"}, elem...)...)
}

func loadLongstrGolden(t testing.TB, path string) *longstrFile {
	t.Helper()
	testutil.RequireFile(t, path)
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	f, err := decodeLongstr(b)
	if err != nil {
		t.Fatalf("decode %s: %v", path, err)
	}
	return f
}

// toolHdrref is the header harness header/testdata/hdrref (built by
// `make -C header/testdata hdrref/hdrref` for tier 2).
var toolHdrref = testutil.Tool{Name: "hdrref", Harness: "header/testdata/hdrref/hdrref", Sources: []string{"header/testdata/hdrref/hdrref.c"}}

// runHdrref is hdrref's output for args; every argument that names an
// existing file is an input of the recorded reading.
func runHdrref(t testing.TB, args ...string) []byte {
	t.Helper()
	return runHdrrefWriting(t, nil, args...)
}

// runHdrrefWriting is runHdrref for a call that writes the files outputs.
func runHdrrefWriting(t testing.TB, outputs []string, args ...string) []byte {
	t.Helper()
	var inputs []string
	for _, a := range args {
		if st, err := os.Stat(a); err == nil && st.Mode().IsRegular() {
			inputs = append(inputs, a)
		}
	}
	return testutil.RunOK(t, toolHdrref, testutil.Call{Args: args, Inputs: inputs, Outputs: outputs}).Stdout
}

func hdrrefLongstr(t testing.TB, path string) *longstrFile {
	t.Helper()
	f, err := decodeLongstr(runHdrref(t, "longstr", path))
	if err != nil {
		t.Fatalf("decode hdrref longstr %s: %v", path, err)
	}
	return f
}

func trimCard(b [CardWidth]byte) string { return strings.TrimRight(string(b[:]), " ") }

// wantGoKey is the Key the Go decoder must report for a card cfitsio
// names k (fits_get_keyname). Non-HIERARCH names must be equal byte for
// byte. For a HIERARCH card cfitsio returns the path alone, with any runs
// of blanks the card has; the Go Key is "HIERARCH " plus the path with
// single blanks (documented on Card), and nothing else is accepted.
func wantGoKey(c *Card, k string) string {
	if c.IsHierarch() {
		return "HIERARCH " + strings.Join(strings.Fields(k), " ")
	}
	return k
}

// compareWithCfitsio asserts that the logical cards Go decoded from one
// header equal cfitsio's reading of the same bytes. Each Go card stands for
// 1 + (number of CONTINUE cards it absorbed) physical cards; together they
// must cover every physical card cfitsio lists, in order.
//
// Per card (cfitsio typing from fits_get_keytype, values from its typed
// reads):
//
//	no value          -> TypeEmpty, nil value
//	parse error       -> TypeUnparsed
//	'C' string        -> TypeString == fits_read_key_longstr
//	'L' "T"/"F"       -> TypeLogical == fits_read_key_log
//	'L' other text    -> TypeUnparsed, Text == cfitsio's value text
//	                     (cfitsio reads "TRUE"/"Tx" as true; Go keeps the
//	                     text and does not guess)
//	'I' int ok        -> TypeInt == fits_read_key_lnglng
//	'I' overflow      -> TypeBigInt; float64 == fits_read_key_dbl, exact
//	                     value == fits_read_key(TULONGLONG) where in range
//	'I'/'F' float ok  -> TypeFloat == fits_read_key_dbl (bit-exact)
//	'I'/'F' no read   -> TypeUnparsed, Text == cfitsio's value text
//	'X' complex       -> TypeComplex*, == fits_read_key(TDBLCOMPLEX)
//
// Comments equal cfitsio's per-card comment; for a CONTINUE'd string they
// equal fits_read_key_longstr's joined comment, which cfitsio truncates to
// 72 characters (Go keeps the whole comment, so cfitsio's is its prefix).
func compareWithCfitsio(t *testing.T, label string, cards []Card, gh *testutil.CfrefHDU, lh *longstrHDU) {
	t.Helper()
	g := 0
	for ci := range cards {
		c := &cards[ci]
		if g == len(gh.Cards) && trailingBlanks(cards, ci) {
			break
		}
		phys := c.RawCards()
		if g+len(phys) > len(gh.Cards) {
			t.Errorf("%s: Go card %d (%s) needs %d physical cards, cfitsio has only %d left", label, ci, c.Key, len(phys), len(gh.Cards)-g)
			return
		}
		want := &gh.Cards[g]
		where := fmt.Sprintf("%s card %d (%s)", label, want.N, strings.TrimSpace(string(want.Raw)))
		for j, raw := range phys {
			if got, w := trimCard(raw), string(gh.Cards[g+j].Raw); got != w {
				t.Errorf("%s: physical card %d raw\n got %q\nwant %q", where, j, got, w)
			}
		}
		if got, w := c.Key, wantGoKey(c, string(want.Key)); got != w {
			t.Errorf("%s: key %q, cfitsio %q", where, got, w)
		}
		ls := lh.card(want.N)
		// Comment.
		if len(phys) == 1 {
			if c.Comment != string(want.Comment) {
				t.Errorf("%s: comment %q, cfitsio %q", where, c.Comment, want.Comment)
			}
		} else {
			if ls == nil || ls.Comment == nil {
				t.Errorf("%s: continued string but no longstr golden comment", where)
			} else {
				w := string(*ls.Comment)
				got := c.Comment
				if len(got) > 72 {
					got = got[:72]
				}
				if got != w {
					t.Errorf("%s: joined comment %q (first 72 bytes %q), cfitsio %q", where, c.Comment, got, w)
				}
			}
		}
		checkValue(t, where, c, want, ls)
		g += len(phys)
	}
	if g != len(gh.Cards) {
		t.Errorf("%s: Go cards cover %d physical cards, cfitsio has %d", label, g, len(gh.Cards))
	}
}

// trailingBlanks reports whether every card from i on is an all-blank
// card. cfitsio does not count blank cards just before END as keywords
// (it treats them as free header space); Go keeps them as blank
// commentary cards so the header re-encodes byte for byte.
func trailingBlanks(cards []Card, i int) bool {
	for ; i < len(cards); i++ {
		if cards[i].Raw != newCard("") {
			return false
		}
	}
	return true
}

func checkValue(t *testing.T, where string, c *Card, want *testutil.CfrefCard, ls *longstrCard) {
	t.Helper()
	text := string(want.Value)
	unparsed := func(why string) {
		t.Helper()
		rv, ok := c.Value.(RawValue)
		if c.Type != TypeUnparsed || !ok {
			t.Errorf("%s: %s: want TypeUnparsed, got type %v value %#v", where, why, c.Type, c.Value)
			return
		}
		if want.ParseError == "" && rv.Text != text {
			t.Errorf("%s: unparsed text %q, cfitsio value %q", where, rv.Text, text)
		}
	}
	if want.ParseError != "" {
		unparsed("cfitsio parse error " + want.ParseError)
		return
	}
	if text == "" {
		if c.Type != TypeEmpty || c.Value != nil {
			t.Errorf("%s: cfitsio sees no value; got type %v value %#v", where, c.Type, c.Value)
		}
		return
	}
	switch want.Type {
	case "C":
		if want.Str == nil {
			t.Errorf("%s: cfitsio could not read the string: %s", where, want.ReadError)
			return
		}
		if s, ok := c.Value.(string); c.Type != TypeString || !ok || s != string(*want.Str) {
			t.Errorf("%s: string %#v (type %v), cfitsio %q", where, c.Value, c.Type, *want.Str)
		}
	case "L":
		if text != "T" && text != "F" {
			unparsed("cfitsio coerces " + text + " to a logical")
			return
		}
		if b, ok := c.Value.(bool); c.Type != TypeLogical || !ok || want.Logical == nil || b != *want.Logical {
			t.Errorf("%s: logical %#v (type %v), cfitsio %v", where, c.Value, c.Type, want.Logical)
		}
	case "I", "F":
		switch {
		case strings.ContainsAny(text, "xX"):
			// C strtod reads hexadecimal; FITS has no hexadecimal numbers.
			unparsed("hexadecimal token (cfitsio's strtod accepts it)")
		case want.Int != nil:
			if v, ok := c.Value.(int64); c.Type != TypeInt || !ok || v != *want.Int {
				t.Errorf("%s: int %#v (type %v), cfitsio %d", where, c.Value, c.Type, *want.Int)
			}
		case want.Type == "I" && want.Float != nil && isIntegerLiteral(text):
			// Integer literal outside int64: cfitsio's integer read
			// overflows, its double read succeeds.
			b, ok := c.Value.(*big.Int)
			if c.Type != TypeBigInt || !ok {
				t.Errorf("%s: want TypeBigInt, got type %v value %#v", where, c.Type, c.Value)
				return
			}
			f, _ := new(big.Float).SetInt(b).Float64()
			if math.Float64bits(f) != math.Float64bits(float64(*want.Float)) {
				t.Errorf("%s: big int as float64 %v, cfitsio %v", where, f, float64(*want.Float))
			}
			if ls != nil && ls.U64 != nil && (!b.IsUint64() || b.Uint64() != *ls.U64) {
				t.Errorf("%s: big int %s, cfitsio TULONGLONG %d", where, b, *ls.U64)
			}
			if exact, _ := new(big.Int).SetString(strings.TrimPrefix(text, "+"), 10); exact == nil || exact.Cmp(b) != 0 {
				t.Errorf("%s: big int %s does not equal the card text %q", where, b, text)
			}
		case want.Float != nil:
			v, ok := c.Value.(float64)
			if c.Type != TypeFloat || !ok || math.Float64bits(v) != math.Float64bits(float64(*want.Float)) {
				t.Errorf("%s: float %#v (type %v), cfitsio %v", where, c.Value, c.Type, float64(*want.Float))
			}
		default:
			unparsed("cfitsio cannot read it as a number (" + want.ReadError + ")")
		}
	case "X":
		if ls == nil {
			t.Errorf("%s: no cfitsio complex read in longstr golden", where)
			return
		}
		if ls.Cplx == nil {
			unparsed("cfitsio cannot read it as a complex (" + ls.Error + ")")
			return
		}
		if inner := strings.Split(text[1:len(text)-1], ","); len(inner) == 2 &&
			(strings.TrimSpace(inner[0]) == "" || strings.TrimSpace(inner[1]) == "") {
			// cfitsio reads an empty part as 0; FITS has no such value.
			unparsed("complex with an empty part (cfitsio reads it as 0)")
			return
		}
		v, ok := c.Value.(Complex)
		if (c.Type != TypeComplexInt && c.Type != TypeComplexFloat) || !ok ||
			math.Float64bits(v.Re) != math.Float64bits(ls.Cplx[0]) || math.Float64bits(v.Im) != math.Float64bits(ls.Cplx[1]) {
			t.Errorf("%s: complex %#v (type %v), cfitsio %v", where, c.Value, c.Type, *ls.Cplx)
		}
	default:
		t.Errorf("%s: unexpected cfitsio type %q", where, want.Type)
	}
}
