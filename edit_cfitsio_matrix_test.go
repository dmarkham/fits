package fits

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// editMatrixCases are header edits pushed through Flush and read back by
// cfitsio (the cfref harness: fits_read_key_longstr for strings,
// fits_get_keytype plus the typed read for numbers and logicals). The
// guard in encodeVerified compares the encoder with this package's own
// parser; this matrix checks the result against cfitsio instead, so an
// edit that both of ours agree on but cfitsio reads differently fails
// here.
//
// A case with refused set is one the encoder cannot write exactly: Flush
// must refuse it and leave the file byte-identical, and cfitsio must read
// the card the encoder would have written as something other than the
// value that was set (so the refusal is not over-cautious).
var editMatrixCases = []struct {
	name    string
	key     string // as passed to Set
	cfKey   string // fits_get_keyname (HIERARCH prefix removed)
	value   any
	comment string
	refused bool
}{
	{name: "hierarch-string", key: "HIERARCH ESO DET NAME", cfKey: "ESO DET NAME", value: "ccd-42", comment: "detector"},
	{name: "hierarch-float", key: "HIERARCH ESO TEL ALT", cfKey: "ESO TEL ALT", value: 61.25, comment: "deg"},
	{name: "ampersand-no-continue", key: "AMPSTR", value: "abc&", comment: "not a long string"},
	{name: "slash-in-value", key: "SLASHV", value: "a/b // c", comment: "x"},
	{name: "slash-in-comment", key: "SLASHC", value: int64(3), comment: "one / two // three"},
	{name: "quotes", key: "QUOTES", value: "it's ''q''", comment: "a 'quoted' comment"},
	{name: "int-2^53+1", key: "BIGINT", value: int64(1<<53 + 1), comment: ""},
	{name: "int-min", key: "MININT", value: int64(math.MinInt64), comment: ""},
	{name: "int-max", key: "MAXINT", value: int64(math.MaxInt64), comment: ""},
	{name: "float32-widened", key: "F32WIDE", value: float64(float32(0.1)), comment: ""},
	{name: "float-huge", key: "FHUGE", value: 1.7976931348623157e308, comment: ""},
	{name: "float-min-normal", key: "FMINNORM", value: 2.2250738585072014e-308, comment: ""},
	{name: "float-negative-tiny", key: "FNEG", value: -1.5e-300, comment: ""},
	{name: "string-80-chars", key: "LONG80", value: strings.Repeat("0123456789", 8), comment: "becomes CONTINUE"},
	{name: "string-equals-x", key: "EQX", value: "=x", comment: ""},
	{name: "string-one-blank", key: "BLANK1", value: " ", comment: ""},
	{name: "string-empty", key: "EMPTY", value: "", comment: ""},
	{name: "logical-false", key: "LOGF", value: false, comment: "no"},
	{name: "logical-true", key: "LOGT", value: true, comment: ""},

	// Subnormal floats underflow strtod (cfitsio: NUM_OVERFLOW); our parser
	// rejects them too, so the card would not read back.
	{name: "float-subnormal", key: "TINY", value: 5e-324, refused: true},
	{name: "float-subnormal-2", key: "TINY2", value: -1e-310, refused: true},
	// -0.0 keeps its sign (written "-0.0"; cfitsio's strtod reads -0).
	{name: "float-negative-zero", key: "NEG0", value: math.Copysign(0, -1)},
}

// cfValueMatches reports whether cfitsio's reading of card c is the value
// v that was set (strings: trailing blanks not significant; floats: the
// same bits, so the sign of zero counts).
func cfValueMatches(c *testutil.CfrefCard, v any) bool {
	switch v := v.(type) {
	case string:
		return c.Type == "C" && c.Str != nil && strings.TrimRight(string(*c.Str), " ") == strings.TrimRight(v, " ")
	case int64:
		return c.Type == "I" && c.Int != nil && *c.Int == v
	case float64:
		return c.Type == "F" && c.Float != nil && math.Float64bits(float64(*c.Float)) == math.Float64bits(v)
	case bool:
		return c.Type == "L" && c.Logical != nil && *c.Logical == v
	}
	return false
}

func TestFlushEditMatrixMatchesCfitsio(t *testing.T) {
	for _, tc := range editMatrixCases {
		t.Run(tc.name, func(t *testing.T) {
			cfKey := tc.cfKey
			if cfKey == "" {
				cfKey = tc.key
			}
			path, orig := copySample(t, "testprog.std")
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			h := hduHeader(t, f, 2)
			if err := h.Set(tc.key, tc.value, tc.comment); err != nil {
				t.Fatalf("Set: %v", err)
			}
			encoded, encErr := header.Encode(h)
			ferr := f.Flush()
			f.Close()
			if tc.refused {
				if ferr == nil {
					t.Fatalf("Flush wrote %s = %v; it cannot be written exactly", tc.key, tc.value)
				}
				if !strings.Contains(ferr.Error(), "cannot represent") {
					t.Errorf("Flush error %q does not name the card the encoder cannot represent", ferr)
				}
				if !bytes.Equal(mustRead(t, path), orig) {
					t.Fatal("refused Flush changed the file")
				}
				if encErr != nil {
					return // the encoder itself refuses: nothing would be written
				}
				// What cfitsio would have read from the encoder's card.
				c := cfitsioReadsEncodedCard(t, encoded, tc.key, cfKey)
				if c != nil && cfValueMatches(c, tc.value) {
					t.Fatalf("cfitsio reads the encoded card %q as the value that was set; the refusal is not needed", c.Raw)
				}
				return
			}
			if ferr != nil {
				t.Fatalf("Flush: %v", ferr)
			}
			cf := testutil.CfrefDump(t, path)
			c := cf.HDUs[2].Card(cfKey)
			if c == nil {
				t.Fatalf("cfitsio finds no %s card", cfKey)
			}
			if !cfValueMatches(c, tc.value) {
				t.Fatalf("cfitsio reads %s as type %q (str %v int %v float %v logical %v, raw %q), set %v (%T)",
					cfKey, c.Type, c.Str, c.Int, c.Float, c.Logical, c.Raw, tc.value, tc.value)
			}
			// A long string: fits_read_key_longstr takes the comment from the
			// CONTINUE cards that follow, which cfref lists as their own cards.
			var last *testutil.CfrefCard
			for i := range cf.HDUs[2].Cards {
				if cc := &cf.HDUs[2].Cards[i]; cc.Key == "CONTINUE" && cc.N == c.N+1+countContinue(last, c) {
					last = cc
				}
			}
			if _, isStr := tc.value.(string); isStr && len(tc.value.(string)) > 68 {
				if last == nil || !strings.HasSuffix(strings.TrimRight(string(last.Raw), " "), "/ "+tc.comment) {
					t.Errorf("last CONTINUE card %+v does not carry the comment %q", last, tc.comment)
				}
			} else if last != nil {
				t.Errorf("a short value is followed by CONTINUE card %q", last.Raw)
			} else if got := strings.TrimSpace(string(c.Comment)); got != tc.comment {
				t.Errorf("cfitsio comment %q, set %q", got, tc.comment)
			}
			if v := testutil.Fitsverify(t, path); v.Errors != 0 {
				t.Errorf("fitsverify: %d errors\n%s", v.Errors, v.Report)
			}
		})
	}
}

// cfitsioReadsEncodedCard puts the cards for key out of an encoded header
// into a minimal primary HDU and returns cfitsio's reading of the card.
func cfitsioReadsEncodedCard(t *testing.T, encoded []byte, key, cfKey string) *testutil.CfrefCard {
	t.Helper()
	var cards []string
	for i := 0; i+80 <= len(encoded); i += 80 {
		c := string(encoded[i : i+80])
		k := strings.TrimRight(c[:8], " ")
		if strings.HasPrefix(c, key) && (len(key) >= 8 || k == key) {
			cards = append(cards, c)
		}
	}
	if len(cards) == 0 {
		t.Fatalf("encoded header has no %s card", key)
	}
	var b bytes.Buffer
	for _, c := range append([]string{"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0"}, cards...) {
		fmt.Fprintf(&b, "%-80s", c)
	}
	fmt.Fprintf(&b, "%-80s", "END")
	for b.Len()%2880 != 0 {
		b.WriteByte(' ')
	}
	p := filepath.Join(t.TempDir(), "card.fits")
	if err := os.WriteFile(p, b.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	cf, err := testutil.CfrefDumpErr(t, p)
	if err != nil {
		t.Logf("cfitsio cannot read the encoded card %q: %v", cards[0], err)
		return nil
	}
	c := cf.HDUs[0].Card(cfKey)
	if c != nil {
		t.Logf("cfitsio reads the encoded card %q as type %q float %v (error %q %q)", c.Raw, c.Type, c.Float, c.Error, c.ReadError)
	}
	return c
}

// countContinue is how many CONTINUE cards after c have been taken so far.
func countContinue(last, c *testutil.CfrefCard) int {
	if last == nil {
		return 0
	}
	return last.N - c.N
}
