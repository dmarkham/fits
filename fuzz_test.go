package fits

import (
	"bytes"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// FuzzParseHeader exercises the header parser against random 2880-byte
// blocks. The parser must never panic and any successful parse must
// round-trip through Encode+ParseCards again.
func FuzzParseHeader(f *testing.F) {
	// Seed with a minimal valid header.
	seed := mkSeedBlock([]string{
		"SIMPLE  =                    T",
		"BITPIX  =                    8",
		"NAXIS   =                    0",
		"END",
	})
	f.Add(seed)

	f.Fuzz(func(t *testing.T, buf []byte) {
		if len(buf)%2880 != 0 {
			return
		}
		cards, _, err := header.ParseCards(buf)
		if err != nil {
			// Errors are fine; we only require no panic.
			return
		}
		// Round-trip: encode what we parsed, re-parse, compare keys.
		h := header.FromCards(cards)
		encoded, err := header.Encode(h)
		if err != nil {
			t.Fatalf("Encode failed on parsed input: %v", err)
		}
		cards2, _, err := header.ParseCards(encoded)
		if err != nil {
			t.Fatalf("re-parse failed: %v", err)
		}
		if len(cards) != len(cards2) {
			t.Fatalf("card count changed: %d → %d", len(cards), len(cards2))
		}
	})
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
