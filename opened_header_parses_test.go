package fits

import (
	"bytes"
	"errors"
	"fmt"
	"math/rand/v2"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// Every header of a file Open accepted parses, and parses to the same
// cards the scan used: HeaderOf never fails, so the Header() methods (and
// Columns, CompressionType, which call them) cannot panic on it, and the
// parsed header ends at the END card the scan found, so Header() holds
// every card behind Shape() and the data offset. The argument, which this
// test checks on mutated real files:
//
//   - Open keeps a header only when extractStructural found the END card
//     ("END" and 77 blanks; its loop stops there; with none it fails), and
//     the header bytes are a whole number of 80-byte cards. A card before it
//     that is named END but is not the END card ("END  = 1", "END=1",
//     "END     junk", "END" with NULs: the name runs to the first blank, '='
//     or NUL, as in cfitsio's ffgknm) makes Open fail with ErrMalformedEND;
//   - header.ParseCards fails only on a length that is not a multiple of
//     80, on a DecodeCard error (which happens only for a card that is not
//     80 bytes), or when no card decodes as END; DecodeCard decodes a card
//     as END only when its first 8 bytes are "END     " or its name (up to
//     the first blank or '=') is END, and Open has refused every such card
//     before the END card.
//
// So cmd/fitsdump has no undecodable header to report.

// openedHeadersParse opens b and, if Open accepts it, checks that every
// HDU's header parses. It reports whether Open accepted the file.
func openedHeadersParse(t testing.TB, dir string, b []byte) bool {
	p := filepath.Join(dir, "f.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := Open(p)
	if err != nil {
		return false
	}
	defer f.Close()
	for i, rec := range f.hdus {
		if _, err := rec.loadHeader(); err != nil {
			t.Fatalf("HDU %d of a file Open accepted: header does not parse: %v", i, err)
		}
		// The parse must stop at the END card the scan stopped at.
		_, endIdx, err := header.ParseCards(rec.rawHeader)
		if err != nil {
			t.Fatal(err)
		}
		scanEnd := -1
		for k := 0; k+80 <= len(rec.rawHeader); k += 80 {
			if bytes.Equal(rec.rawHeader[k:k+80], endCard) {
				scanEnd = k / 80
				break
			}
		}
		if endIdx != scanEnd {
			t.Fatalf("HDU %d of a file Open accepted: the parsed header ends at card %d, the scan's END card is card %d", i, endIdx, scanEnd)
		}
		h, err := f.HDU(i)
		if err != nil {
			continue // random groups, unsupported extension: reported by HDU
		}
		_ = h.Header()
		_ = h.CompressionType()
		if bt, ok := h.(*BinaryTableHDU); ok {
			_, _ = bt.Columns()
		}
	}
	return true
}

// TestOpenedHeadersAlwaysParse: crafted END and CONTINUE variants, then
// 3000 seeded mutations of the headers of committed samples (a random card
// replaced by END variants, NULs, lone CONTINUE, non-ASCII or random
// bytes).
func TestOpenedHeadersAlwaysParse(t *testing.T) {
	dir := t.TempDir()
	card := func(s string) []byte {
		b := bytes.Repeat([]byte(" "), 80)
		copy(b, s)
		return b
	}
	endNUL := make([]byte, 80)
	copy(endNUL, "END")
	variants := [][]byte{
		card("END"), card("END     junk"), card("END     = 1"), endNUL, card("ENDX"),
		card("CONTINUE  'def'"), card("CONTINUE  'abc&'"), card("CONTINUE"),
		card("LONGSTR = 'abc&'"), card("HIERARCH"), card("HIERARCH A B = "),
		card("KEY     = '\xff\x80'"), card("KEY     = 'unterminated"), card("=       = 1"),
		bytes.Repeat([]byte{0}, 80), bytes.Repeat([]byte{0xff}, 80),
	}
	prim := func(extra ...[]byte) []byte {
		var b []byte
		for _, c := range []string{"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0"} {
			b = append(b, card(c)...)
		}
		for _, c := range extra {
			b = append(b, c...)
		}
		b = append(b, card("END")...)
		for len(b)%2880 != 0 {
			b = append(b, ' ')
		}
		return b
	}
	accepted := 0
	for _, v := range variants {
		for _, b := range [][]byte{prim(v), prim(v, v), prim(card("LONGSTR = 'abc&'"), v)} {
			if openedHeadersParse(t, dir, b) {
				accepted++
			}
		}
	}

	rng := rand.New(rand.NewPCG(1, 2))
	samples := []string{"testprog.std", "iter_a.fit", "comp_rice_i16.fits", "vari.fits"}
	for _, s := range samples {
		src := testutil.TestdataPath(t, s)
		testutil.RequireFile(t, src)
		orig, err := os.ReadFile(src)
		if err != nil {
			t.Fatal(err)
		}
		f, err := Open(src)
		if err != nil {
			t.Fatal(err)
		}
		type span struct{ start, end int64 }
		var hdrs []span
		for _, r := range f.hdus {
			hdrs = append(hdrs, span{r.headerStart, r.dataStart})
		}
		f.Close()
		for range 750 {
			b := append([]byte(nil), orig...)
			for range 1 + rng.IntN(3) {
				h := hdrs[rng.IntN(len(hdrs))]
				off := h.start + 80*rng.Int64N((h.end-h.start)/80)
				var c []byte
				if rng.IntN(4) == 0 {
					c = make([]byte, 80)
					for i := range c {
						c[i] = byte(rng.IntN(256))
					}
				} else {
					c = variants[rng.IntN(len(variants))]
				}
				copy(b[off:off+80], c)
			}
			if openedHeadersParse(t, dir, b) {
				accepted++
			}
		}
	}
	if accepted < 1000 {
		t.Fatalf("only %d mutated files were accepted by Open; the property is barely exercised", accepted)
	}
	t.Logf("%d accepted files, every header parsed", accepted)
}

// FuzzOpenedHeadersParse is the same property over arbitrary bytes.
func FuzzOpenedHeadersParse(f *testing.F) {
	for _, s := range []string{"iter_a.fit", "testprog.std"} {
		b, err := os.ReadFile(testutil.TestdataPath(f, s))
		if err != nil {
			f.Fatal(err)
		}
		f.Add(b)
	}
	dir := f.TempDir()
	f.Fuzz(func(t *testing.T, b []byte) {
		openedHeadersParse(t, dir, b)
	})
}

// A card named END that is not the END card. Where a required keyword must
// be, cfitsio refuses the file ("Expected positive integer keyword NAXIS1,
// but instead found keyword END"); before this check Open accepted it and
// read pixels with NAXISn cards taken from after the fake END, while
// Header() stopped at it. Elsewhere cfitsio ends the header at that card,
// which puts the data unit at a different offset than the END card does
// when they are in different blocks. Open refuses every such card
// (ErrMalformedEND). Truth: cfitsio (testdata/cfref dump) for the refusal
// at a required position; the FITS Standard (4.0, 4.4.1.1) for the rest.
func TestMalformedENDRefused(t *testing.T) {
	card := func(s string) string { return fmt.Sprintf("%-80s", s) }
	build := func(cards ...string) []byte {
		var b []byte
		for _, c := range append(cards, "END") {
			b = append(b, card(c)...)
		}
		for len(b)%2880 != 0 {
			b = append(b, ' ')
		}
		b = append(b, 1, 2, 3, 4, 5, 6, 7, 8)
		for len(b)%2880 != 0 {
			b = append(b, 0)
		}
		return b
	}
	head := []string{"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    2"}
	axes := []string{"NAXIS1  =                    4", "NAXIS2  =                    2"}
	nul := "END" + strings.Repeat("\x00", 77)
	cases := []struct {
		name string
		fake string
		// required: the fake END is where NAXIS1 must be (cfitsio refuses).
		required bool
		// later: 40 COMMENT cards between the fake END and the END card, so
		// the END card is in the next block.
		later bool
	}{
		{"END  = 1/required", "END  = 1", true, false},
		{"END=1/required", "END=1", true, false},
		{"END     junk/required", "END     junk", true, false},
		{"END  = 1/after-axes", "END  = 1", false, false},
		{"END     junk/after-axes", "END     junk", false, false},
		{"END     = 1/after-axes", "END     = 1", false, false},
		{"END=1/after-axes-next-block", "END=1", false, true},
		{"END NUL-padded/after-axes", nul, false, false},
		// cfitsio's ffgknm names a HIERARCH card by the text between
		// "HIERARCH " and '=' (blanks trimmed), so these are END for it.
		{"HIERARCH END = 1/required", "HIERARCH END = 1", true, false},
		{"HIERARCH END = 1/after-axes", "HIERARCH END = 1", false, false},
		{"HIERARCH END = 1/after-axes-next-block", "HIERARCH END = 1", false, true},
		{"HIERARCH  END= 'x'/after-axes-next-block", "HIERARCH  END= 'x'", false, true},
		{"HIERARCH END    =/after-axes", "HIERARCH END    =", false, false},
	}
	dir := t.TempDir()
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			cards := append([]string(nil), head...)
			if c.required {
				cards = append(cards, c.fake)
				cards = append(cards, axes...)
			} else {
				cards = append(cards, axes...)
				cards = append(cards, c.fake)
			}
			if c.later {
				for k := range 40 {
					cards = append(cards, fmt.Sprintf("COMMENT filler %d", k))
				}
			}
			b := build(cards...)
			p := filepath.Join(dir, strings.NewReplacer(" ", "_", "=", "_", "/", "_").Replace(c.name)+".fits")
			if err := os.WriteFile(p, b, 0o644); err != nil {
				t.Fatal(err)
			}
			f, err := Open(p)
			if err == nil {
				img, _ := f.Primary()
				var shape []int64
				if img != nil {
					shape = img.Shape()
				}
				f.Close()
				t.Fatalf("Open accepted a header with %q (shape %v)", c.fake, shape)
			}
			if !errors.Is(err, ErrMalformedEND) {
				t.Fatalf("Open: %v, want ErrMalformedEND", err)
			}
			if c.later {
				// cfitsio ends the header at the fake END: its data unit
				// starts in the block the END card is in.
				cf, err := testutil.CfrefDumpErr(t, p)
				if err != nil {
					t.Fatalf("cfitsio: %v", err)
				}
				if cf.HDUs[0].DataStart != 2880 {
					t.Fatalf("cfitsio puts the data at %d, not at the block after the fake END (2880); the test's premise is wrong", cf.HDUs[0].DataStart)
				}
			}
			if c.required {
				if _, err := testutil.CfrefDumpErr(t, p); err == nil {
					t.Fatalf("cfitsio accepted the file; the test's premise is wrong")
				} else {
					t.Logf("cfitsio refuses too: %v", err)
				}
			}
		})
	}
	// The END card itself, and names that only start with END, are not
	// affected.
	// HIERARCH paths that are not exactly END (cfitsio compares the name
	// case-sensitively) are ordinary cards for cfitsio too.
	for _, c := range []string{"ENDX    = 1", "ENDING  = 'x'", "HIERARCH ENDX = 1", "HIERARCH ESO END = 1", "HIERARCH end = 1", "HIERARCH END"} {
		b := build(append(append(append([]string(nil), head...), axes...), c)...)
		p := filepath.Join(dir, "ok.fits")
		if err := os.WriteFile(p, b, 0o644); err != nil {
			t.Fatal(err)
		}
		f, err := Open(p)
		if err != nil {
			t.Fatalf("Open with %q: %v", c, err)
		}
		f.Close()
		cf, err := testutil.CfrefDumpErr(t, p)
		if err != nil {
			t.Fatalf("cfitsio with %q: %v", c, err)
		}
		if cf.HDUs[0].NKeys != 6 {
			t.Fatalf("cfitsio reads %d keywords with %q, want 6 (it is not END for cfitsio)", cf.HDUs[0].NKeys, c)
		}
	}
}
