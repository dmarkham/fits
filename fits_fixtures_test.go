package fits

import (
	"bytes"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// allSamples returns every committed FITS sample in testdata/ (the same
// set testdata/cfref has goldens for).
func allSamples(t *testing.T) []string {
	t.Helper()
	var out []string
	for _, pat := range []string{"*.fits", "*.fit", "*.std"} {
		m, err := filepath.Glob(testutil.TestdataPath(t, pat))
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, m...)
	}
	sort.Strings(out)
	if len(out) < 43 {
		t.Fatalf("found %d samples in testdata, want at least 43", len(out))
	}
	return out
}

// TestCopyToByteIdentical: for every committed sample, (*File).CopyTo
// reproduces the input byte for byte.
//
// This proves only the framing pass-through (HDU boundaries found and every
// header and data byte copied in order). CopyTo does not decode any value,
// so it says nothing about reading pixels, table cells or keyword values;
// those are proven against the cfitsio goldens by the image, table and
// header tests. (This test used to be called TestRoundTripFixtures and was
// cited as validating the library against cfitsio's test files.)
func TestCopyToByteIdentical(t *testing.T) {
	for _, path := range allSamples(t) {
		name := filepath.Base(path)
		for _, op := range openPaths {
			t.Run(name+"/"+op.name, func(t *testing.T) {
				orig, err := os.ReadFile(path)
				if err != nil {
					t.Fatalf("committed sample unreadable: %v", err)
				}
				f, err := op.open(t, path)
				if err != nil {
					t.Fatalf("Open: %v", err)
				}
				defer f.Close()
				var buf bytes.Buffer
				if err := f.CopyTo(&buf); err != nil {
					t.Fatalf("CopyTo: %v", err)
				}
				if out := buf.Bytes(); !bytes.Equal(out, orig) {
					t.Fatalf("byte mismatch: lens %d (out) vs %d (orig), first diff at offset %d", len(out), len(orig), firstDiff(out, orig))
				}
			})
		}
	}
}

// CopyTo writes the HDUs, not the file: bytes after the last HDU
// (TrailingBytes) are not copied, and a last data unit that lacks its final
// padding is written with it (zeros, FITS 4.0 section 3.3.2), so the output
// is the input cut at the end of the last HDU and padded to whole blocks.
// Truth: testprog.std (the copy of testprog.std with junk after it must be
// testprog.std again) and cfitsio's own f32_plain.fits (the copy of it cut
// inside its last block must be cfitsio's padded file).
func TestCopyToDropsTrailingBytesAndPads(t *testing.T) {
	read := func(elem ...string) []byte {
		p := testutil.TestdataPath(t, elem...)
		testutil.RequireFile(t, p)
		b, err := os.ReadFile(p)
		if err != nil {
			t.Fatal(err)
		}
		return b
	}
	prog := read("testprog.std")
	f32 := read("images", "cfwrite", "f32_plain.fits")
	cases := []struct {
		name     string
		in, want []byte
		trailing int64
	}{
		{"junk-after-last-HDU", append(append([]byte(nil), prog...), []byte("thirteen junk")...), prog, 13},
		{"zero-block-after-last-HDU", append(append([]byte(nil), prog...), make([]byte, 2880)...), prog, 2880},
		{"last-block-cut", f32[:2924], f32, 0},
		{"last-block-cut-to-data-end", f32[:2880+256], f32, 0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			p := filepath.Join(t.TempDir(), "in.fits")
			if err := os.WriteFile(p, c.in, 0o644); err != nil {
				t.Fatal(err)
			}
			f, err := Open(p)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			if f.TrailingBytes() != c.trailing {
				t.Fatalf("TrailingBytes = %d, want %d", f.TrailingBytes(), c.trailing)
			}
			var buf bytes.Buffer
			if err := f.CopyTo(&buf); err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(buf.Bytes(), c.want) {
				t.Fatalf("CopyTo output: len %d, want %d, first diff at %d", buf.Len(), len(c.want), firstDiff(buf.Bytes(), c.want))
			}
		})
	}
}

// TestHDULayoutMatchesCfitsio: for every committed sample, the scanner
// finds the same HDUs at the same byte offsets as cfitsio
// (fits_get_hduaddrll, recorded in testdata/cfref/golden), classifies them
// as the same kind, and splits each header into the same cards (cfitsio's
// fits_read_record text, trailing blanks stripped; cfitsio stops at the
// last non-blank card, so reserved blank cards before END are not counted).
func TestHDULayoutMatchesCfitsio(t *testing.T) {
	for _, path := range allSamples(t) {
		name := filepath.Base(path)
		t.Run(name, func(t *testing.T) {
			want := testutil.LoadCfrefGolden(t, name)
			f, err := Open(path)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer f.Close()
			if f.NumHDU() != want.NHDU {
				t.Fatalf("NumHDU = %d, cfitsio = %d", f.NumHDU(), want.NHDU)
			}
			for i, w := range want.HDUs {
				r := f.hdus[i]
				if r.headerStart != w.HeaderStart || r.dataStart != w.DataStart || r.paddedEnd != w.DataEnd {
					t.Errorf("HDU %d: header/data/end = %d/%d/%d, cfitsio %d/%d/%d", i, r.headerStart, r.dataStart, r.paddedEnd, w.HeaderStart, w.DataStart, w.DataEnd)
				}
				h, err := f.HDU(i)
				if err != nil {
					t.Fatalf("HDU(%d): %v", i, err)
				}
				var kind string
				switch h.Type() {
				case TypeImage:
					kind = "IMAGE"
				case TypeBinaryTable:
					kind = "BINARY_TBL"
				case TypeASCIITable:
					kind = "ASCII_TBL"
				}
				// cfitsio reports a tile-compressed image as an image.
				if h.Compressed() {
					kind = "IMAGE"
				}
				if kind != w.Type {
					t.Errorf("HDU %d: type %v (%s), cfitsio %s", i, h.Type(), kind, w.Type)
				}
				cards, err := splitCards(r.rawHeader)
				if err != nil {
					t.Fatalf("HDU %d: %v", i, err)
				}
				for len(cards) > 0 && strings.TrimSpace(cards[len(cards)-1]) == "" {
					cards = cards[:len(cards)-1]
				}
				if len(cards) != len(w.Cards) {
					t.Errorf("HDU %d: %d cards, cfitsio %d", i, len(cards), len(w.Cards))
					continue
				}
				for j, c := range cards {
					if got := strings.TrimRight(c, " "); got != string(w.Cards[j].Raw) {
						t.Errorf("HDU %d card %d = %q, cfitsio %q", i, j+1, got, w.Cards[j].Raw)
					}
				}
			}
		})
	}
}

// scaledHDU is an HDU with BSCALE / BZERO accessors (ImageHDU,
// CompressedImageHDU).
type scaledHDU interface {
	BSCALE() (float64, error)
	BZERO() (float64, error)
}

// bscaleOf is h.BSCALE() for a test that expects the cards to be valid.
func bscaleOf(t testing.TB, h scaledHDU) float64 {
	t.Helper()
	v, err := h.BSCALE()
	if err != nil {
		t.Fatalf("BSCALE: %v", err)
	}
	return v
}

// bzeroOf is h.BZERO() for a test that expects the cards to be valid.
func bzeroOf(t testing.TB, h scaledHDU) float64 {
	t.Helper()
	v, err := h.BZERO()
	if err != nil {
		t.Fatalf("BZERO: %v", err)
	}
	return v
}
