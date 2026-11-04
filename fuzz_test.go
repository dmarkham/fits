package fits

import (
	"bytes"
	"fmt"
	"math"
	"math/big"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
	"github.com/dmarkham/fits/internal/tform"
)

// The header fuzzers check values, not counts. Every decoded card must
//
//  1. re-encode, unmodified, to exactly the bytes it was decoded from
//     (Raw, including the CONTINUE cards of a long string), and
//  2. when rebuilt in memory from its decoded Key/Type/Value/Comment alone
//     and accepted by the encoder, decode again to the same key, type,
//     value (floats bit for bit, big integers exactly) and comment.
//
// What a card means is pinned to cfitsio by the golden tests in the
// header package (TestDecodeMatchesCfitsio*); the fuzzers extend the
// encode/decode agreement to arbitrary bytes. Seeds are every header card
// of the committed samples (CONTINUE chains, malformed values, HIERARCH).

// FuzzParseHeader exercises ParseCards on arbitrary 2880-byte multiples.
func FuzzParseHeader(f *testing.F) {
	f.Add(mkSeedBlock([]string{
		"SIMPLE  =                    T",
		"BITPIX  =                    8",
		"NAXIS   =                    0",
		"END",
	}))
	for _, hdr := range seedHeaders(f) {
		f.Add(hdr)
	}
	f.Fuzz(func(t *testing.T, in []byte) {
		// Any input is cut to whole cards and blank-padded to whole blocks,
		// so (nearly) every mutation reaches the parser.
		n := len(in) / header.CardWidth * header.CardWidth
		if n == 0 || n > 8*2880 {
			return
		}
		buf := append([]byte(nil), in[:n]...)
		for len(buf)%2880 != 0 {
			buf = append(buf, ' ')
		}
		cards, end, err := header.ParseCards(buf)
		if err != nil {
			return // only ErrNoEnd can happen for card-aligned input
		}
		// 1. Unmodified re-encode reproduces every card before END.
		enc, err := header.Encode(header.FromCards(cards))
		if err != nil {
			t.Fatalf("Encode of an unmodified parsed header: %v", err)
		}
		n = end * header.CardWidth
		if len(enc) < n+header.CardWidth || !bytes.Equal(enc[:n], buf[:n]) || string(enc[n:n+8]) != "END     " {
			t.Fatalf("unmodified re-encode differs from the input before END")
		}
		// 2. Value-level round trip of the in-memory rebuild. FromCards and
		// the encoder store and write keywords upper-cased (as cfitsio's
		// fits_write_key does), so the expected key is the upper-cased one.
		rebuilt := make([]header.Card, 0, len(cards))
		for _, c := range cards {
			rebuilt = append(rebuilt, header.Card{Key: strings.ToUpper(c.Key), Type: c.Type, Value: c.Value, Comment: c.Comment})
		}
		enc2, err := header.Encode(header.FromCards(rebuilt))
		if err != nil {
			return // the encoder refuses what it cannot write faithfully
		}
		cards2, _, err := header.ParseCards(enc2)
		if err != nil {
			t.Fatalf("re-parse of encoder output: %v", err)
		}
		cards2 = joinCommentary(cards2, rebuilt)
		if len(cards2) != len(rebuilt) {
			t.Fatalf("rebuilt %d cards, re-parsed %d", len(rebuilt), len(cards2))
		}
		for i := range rebuilt {
			if d := cardDiff(rebuilt[i], cards2[i]); d != "" {
				t.Fatalf("card %d: %s\nencoded:\n%s", i, d, enc2)
			}
		}
	})
}

// FuzzCardDecode exercises DecodeCard on single 80-byte cards.
func FuzzCardDecode(f *testing.F) {
	for _, hdr := range seedHeaders(f) {
		for i := 0; i+header.CardWidth <= len(hdr); i += header.CardWidth {
			f.Add(hdr[i : i+header.CardWidth])
		}
	}
	f.Fuzz(func(t *testing.T, raw []byte) {
		c, err := header.DecodeCard(raw)
		if len(raw) != header.CardWidth {
			if err == nil {
				t.Fatalf("DecodeCard accepted %d bytes", len(raw))
			}
			return
		}
		if err != nil {
			t.Fatalf("DecodeCard on an 80-byte card: %v", err)
		}
		if !bytes.Equal(c.Raw[:], raw) {
			t.Fatal("Raw is not the input")
		}
		if err := typeMatchesValue(c); err != nil {
			t.Fatal(err)
		}
		if rc := c.RawCards(); len(rc) != 1 || rc[0] != c.Raw {
			t.Fatalf("RawCards of a single decoded card: %d cards", len(rc))
		}
		out, err := header.EncodeCard(c)
		if err != nil || !bytes.Equal(out[:], raw) {
			t.Fatalf("EncodeCard of a decoded card does not reproduce it: %v", err)
		}
		if c.IsEnd() {
			return
		}
		mem := header.Card{Key: c.Key, Type: c.Type, Value: c.Value, Comment: c.Comment}
		enc, err := header.EncodeCardWithContinuation(mem)
		if err != nil {
			return // refused: not writable faithfully
		}
		var buf []byte
		for _, e := range enc {
			buf = append(buf, e[:]...)
		}
		buf = append(buf, mkSeedBlock([]string{"END"})...)
		buf = buf[:len(buf)/2880*2880]
		cards, _, err := header.ParseCards(buf)
		if err != nil {
			t.Fatalf("re-parse: %v", err)
		}
		cards = joinCommentary(cards, []header.Card{mem})
		if len(cards) != 1 {
			t.Fatalf("one card encoded as %d logical cards:\n%s", len(cards), buf)
		}
		// The encoder writes the keyword upper-cased (as cfitsio's
		// fits_write_key does, TestLowerCaseKeyWrittenUpperCasedLikeCfitsio).
		want := mem
		want.Key = strings.ToUpper(mem.Key)
		if d := cardDiff(want, cards[0]); d != "" {
			t.Fatalf("%s\ninput %q\nencoded %q", d, raw, buf[:len(enc)*80])
		}
	})
}

// joinCommentary merges the 72-character pieces the encoder splits long
// commentary text into, so they line up with the in-memory cards.
func joinCommentary(parsed, want []header.Card) []header.Card {
	var out []header.Card
	j := 0
	for i := range want {
		if j >= len(parsed) {
			break
		}
		c := parsed[j]
		j++
		if want[i].IsCommentary() || want[i].Key == "CONTINUE" {
			for len(c.Comment) < len(want[i].Comment) && j < len(parsed) && parsed[j].Key == c.Key && len(strings.TrimRight(string(parsed[j-1].Raw[:]), " ")) == header.CardWidth {
				c.Comment += parsed[j].Comment
				j++
			}
		}
		out = append(out, c)
	}
	return append(out, parsed[j:]...)
}

func cardDiff(a, b header.Card) string {
	if a.Key != b.Key || a.Type != b.Type {
		return fmt.Sprintf("key/type %q/%v -> %q/%v", a.Key, a.Type, b.Key, b.Type)
	}
	if strings.TrimRight(a.Comment, " ") != b.Comment {
		return fmt.Sprintf("%s: comment %q -> %q", a.Key, a.Comment, b.Comment)
	}
	switch av := a.Value.(type) {
	case float64:
		if bv, ok := b.Value.(float64); !ok || math.Float64bits(av) != math.Float64bits(bv) {
			return fmt.Sprintf("%s: float %v -> %v", a.Key, av, b.Value)
		}
	case header.Complex:
		bv, ok := b.Value.(header.Complex)
		if !ok || math.Float64bits(av.Re) != math.Float64bits(bv.Re) || math.Float64bits(av.Im) != math.Float64bits(bv.Im) {
			return fmt.Sprintf("%s: complex %v -> %v", a.Key, av, b.Value)
		}
	case string:
		// Trailing blanks of a FITS string are not significant (§4.2.1.1);
		// a joined CONTINUE string can end in one, the encoder drops it.
		if bv, ok := b.Value.(string); !ok || strings.TrimRight(av, " ") != strings.TrimRight(bv, " ") {
			return fmt.Sprintf("%s: string %q -> %q", a.Key, av, b.Value)
		}
	case *big.Int:
		if bv, ok := b.Value.(*big.Int); !ok || av.Cmp(bv) != 0 {
			return fmt.Sprintf("%s: big int %v -> %v", a.Key, av, b.Value)
		}
	default:
		if a.Value != b.Value {
			return fmt.Sprintf("%s: value %#v -> %#v", a.Key, a.Value, b.Value)
		}
	}
	return ""
}

// typeMatchesValue checks the documented Go type of Value for each Type.
func typeMatchesValue(c header.Card) error {
	ok := false
	switch c.Type {
	case header.TypeEmpty:
		ok = c.Value == nil
	case header.TypeString:
		_, ok = c.Value.(string)
	case header.TypeLogical:
		_, ok = c.Value.(bool)
	case header.TypeInt:
		_, ok = c.Value.(int64)
	case header.TypeBigInt:
		var b *big.Int
		b, ok = c.Value.(*big.Int)
		ok = ok && b != nil && !b.IsInt64()
	case header.TypeFloat:
		var f float64
		f, ok = c.Value.(float64)
		ok = ok && !math.IsNaN(f) && !math.IsInf(f, 0)
	case header.TypeComplexInt, header.TypeComplexFloat:
		_, ok = c.Value.(header.Complex)
	case header.TypeUnparsed:
		_, ok = c.Value.(header.RawValue)
	}
	if !ok {
		return fmt.Errorf("type %v with value %#v", c.Type, c.Value)
	}
	return nil
}

// seedFiles are the committed samples whose headers seed the header
// fuzzers: the header-package samples and three cfitsio test files. Each
// one is required by name, so a missing sample fails the fuzzers instead
// of silently dropping its seeds.
var seedFiles = []string{
	"header/testdata/adversarial.fits",
	"header/testdata/cfwritten.fits",
	"header/testdata/continue.fits",
	"header/testdata/malformed.fits",
	"testdata/testprog.std",
	"testdata/testf77.std",
	"testdata/iter_c.fit",
}

// seedHeaders returns the header bytes of every HDU of seedFiles.
func seedHeaders(f *testing.F) [][]byte {
	paths := seedFiles
	for _, p := range paths {
		testutil.RequireFile(f, p)
	}
	if m, _ := filepath.Glob("header/testdata/*.fits"); len(m) != 4 {
		f.Fatalf("fuzz seeds: header/testdata has %d .fits samples, seedFiles lists 4; add the new one", len(m))
	}
	var out [][]byte
	for _, p := range paths {
		b, err := os.ReadFile(p)
		if err != nil {
			f.Fatal(err)
		}
		fl, err := OpenReader(bytes.NewReader(b))
		if err != nil {
			f.Fatalf("seed %s: %v", p, err)
		}
		for _, r := range fl.hdus {
			out = append(out, b[r.headerStart:r.dataStart])
		}
	}
	return out
}

// FuzzTForm exercises the binary-table TFORM parser. Must never panic.
func FuzzTForm(f *testing.F) {
	f.Add("1PE(256)")
	f.Add("10E")
	f.Add("1J")
	f.Add("128A")
	f.Add("1QJ(1024)")
	f.Fuzz(func(t *testing.T, s string) {
		_, _ = tform.ParseBinary(s) // no-panic contract only
	})
}

// FuzzBlockFraming exercises Open on arbitrary byte sequences. It must not
// panic; any well-formed input either parses or returns a typed error.
func FuzzBlockFraming(f *testing.F) {
	seed, _ := makeMinimalFITS()
	f.Add(seed)
	f.Fuzz(func(t *testing.T, buf []byte) {
		_, _ = OpenReader(bytes.NewReader(buf))
	})
}

// mkSeedBlock builds a 2880-byte block padded with spaces.
func mkSeedBlock(cards []string) []byte {
	const size = 2880
	const cw = 80
	out := make([]byte, 0, size)
	for _, s := range cards {
		c := make([]byte, cw)
		for i := range c {
			c[i] = ' '
		}
		copy(c, s)
		out = append(out, c...)
	}
	for len(out) < size {
		out = append(out, ' ')
	}
	return out
}

// makeMinimalFITS returns a minimal valid FITS file (NAXIS=0 primary + END).
func makeMinimalFITS() ([]byte, error) {
	h := header.New()
	h.Set("SIMPLE", true, "")
	h.Set("BITPIX", int64(8), "")
	h.Set("NAXIS", int64(0), "")
	return header.Encode(h)
}
