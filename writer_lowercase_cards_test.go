package fits

import (
	"errors"
	"fmt"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// Lower-case keywords in a caller's header, through every writer. A
// lower-case keyword is illegal FITS (§4.1.2.1) and readers disagree on
// it (cfitsio's HDU parser matches BSCALE/BZERO/BLANK and the structural
// keywords case-sensitively, astropy upper-cases every keyword), so:
//
//   - an edited copy is written from its fields under the upper-cased
//     keyword, and FromCards / AddCard store it under that keyword, so a
//     writer that computes the keyword itself (NAXIS1, BITPIX, TFIELDS,
//     ...) leaves the copy out as it leaves out an upper-case one;
//   - an unedited decoded card is refused (header.AddCard, which every
//     writer copies the caller's cards with): no copy of its bytes reads
//     the same in every reader.
//
// The helpers here are this file's own, so the tests do not depend on
// other areas' test files.

// lcDecoded decodes card lines as a header read from a file.
func lcDecoded(t *testing.T, lines ...string) []header.Card {
	t.Helper()
	var buf []byte
	for _, l := range append(append([]string(nil), lines...), "END") {
		buf = append(buf, fmt.Sprintf("%-80s", l)...)
	}
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	cards, _, err := header.ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	return cards
}

// lcWriter is one writer that copies caller header cards.
type lcWriter struct {
	name  string
	write func(f *File, h *header.Header) error
	hdu   int // the HDU the caller's cards go to
}

var lcWriters = []lcWriter{
	{"WriteImage", func(f *File, h *header.Header) error {
		_, err := WriteImage[int16](f, h, []int64{4}, []int16{1, 2, 3, 4})
		return err
	}, 0},
	{"AppendBinaryTable", func(f *File, h *header.Header) error {
		if _, err := WriteImage[int16](f, nil, []int64{1}, []int16{1}); err != nil {
			return fmt.Errorf("primary: %w", err)
		}
		_, err := AppendBinaryTable(f, h, []ColumnData{{Name: "A", DataInt32: []int32{1, 2}}})
		return err
	}, 1},
	{"AppendCompressedImage", func(f *File, h *header.Header) error {
		if _, err := WriteImage[int16](f, nil, []int64{1}, []int16{1}); err != nil {
			return fmt.Errorf("primary: %w", err)
		}
		_, err := AppendCompressedImage[int16](f, h, []int64{4, 4}, make([]int16, 16), CompressOptions{})
		return err
	}, 1},
}

// lcWrite runs w into a new file and returns its path and w's error
// (or Close's).
func lcWrite(t *testing.T, w lcWriter, h *header.Header) (string, error) {
	t.Helper()
	p := filepath.Join(t.TempDir(), "out.fits")
	f, err := Create(p)
	if err != nil {
		t.Fatal(err)
	}
	werr := w.write(f, h)
	if cerr := f.Close(); werr == nil {
		werr = cerr
	}
	return p, werr
}

// astropyImage returns astropy's reading of HDU 0's data, scaled.
func astropyImage(t *testing.T, path string) string {
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

// TestWritersLeaveOutEditedLowerCaseStructuralCopies: an edited copy
// (comment changed) of a decoded lower-case structural card, through
// FromCards and through AddCard, is left out by the writer that computes
// the keyword, as an upper-case one is; the file is the one the writer
// would write without it. Truth: cfitsio (each structural keyword once,
// the pixels / rows), astropy (the image data), fitsverify (no error or
// warning) and this package's Open. Before the fix FromCards kept the
// copy under its lower-case Key, the writers' case-sensitive filter let
// it through and it was written as a second, upper-case card: for
// 'naxis1 = 2' WriteImage returned nil and wrote NAXIS1 = 4 and then
// NAXIS1 = 2, astropy read 2 of the 4 pixels, a 'bitpix = 8' copy was
// read as uint8 [0 1 0 2], a 'naxis = 2' copy failed in astropy, and
// fitsverify failed on all of them.
func TestWritersLeaveOutEditedLowerCaseStructuralCopies(t *testing.T) {
	cases := []struct {
		w    lcWriter
		line string
		key  string
	}{
		{lcWriters[0], "naxis1  =                    2 / lower", "NAXIS1"},
		{lcWriters[0], "bitpix  =                    8 / lower", "BITPIX"},
		{lcWriters[0], "naxis   =                    2 / lower", "NAXIS"},
		{lcWriters[1], "naxis2  =                    9 / lower", "NAXIS2"},
		{lcWriters[1], "tfields =                    1 / lower", "TFIELDS"},
		{lcWriters[1], "naxis1  =                    7 / lower", "NAXIS1"},
	}
	for _, tc := range cases {
		for _, via := range []string{"FromCards", "AddCard"} {
			t.Run(tc.w.name+"/"+tc.key+"/"+via, func(t *testing.T) {
				edited := lcDecoded(t, tc.line, "OBJECT  = 'M31'")
				edited[0].Comment = "edited"
				var h *header.Header
				if via == "FromCards" {
					h = header.FromCards(edited)
				} else {
					h = header.New()
					for _, c := range edited {
						if err := h.AddCard(c); err != nil {
							t.Fatal(err)
						}
					}
				}
				if c := h.Cards()[0]; c.Key != tc.key {
					t.Fatalf("%s stores the edited copy as %q", via, c.Key)
				}
				p, err := lcWrite(t, tc.w, h)
				if err != nil {
					t.Fatalf("write: %v", err)
				}
				if r := testutil.Fitsverify(t, p); !r.OK() {
					t.Errorf("fitsverify:\n%s", r.Report)
				}
				cf := testutil.CfrefDump(t, p)
				gh := &cf.HDUs[tc.w.hdu]
				n := 0
				for _, c := range gh.Cards {
					if strings.EqualFold(string(c.Key), tc.key) {
						n++
						if strings.Contains(string(c.Comment), "edited") {
							t.Errorf("the edited copy was written: %q", c.Raw)
						}
					}
				}
				if n != 1 {
					t.Errorf("cfitsio sees %d %s cards", n, tc.key)
				}
				if gh.Card("OBJECT") == nil {
					t.Error("the OBJECT card after it is missing")
				}
				f, err := Open(p)
				if err != nil {
					t.Fatalf("Open: %v", err)
				}
				defer f.Close()
				switch tc.w.hdu {
				case 0:
					want := []testutil.NullFloat{{V: 1, Valid: true}, {V: 2, Valid: true}, {V: 3, Valid: true}, {V: 4, Valid: true}}
					gh.Image.Phys.Check(t, "cfitsio TDOUBLE", want)
					if got := astropyImage(t, p); got != "1.0 2.0 3.0 4.0" {
						t.Errorf("astropy reads %q", got)
					}
					img, err := f.Primary()
					if err != nil {
						t.Fatal(err)
					}
					px, err := ReadPixels[int16](img)
					if err != nil || fmt.Sprint(px) != "[1 2 3 4]" {
						t.Errorf("ReadPixels = %v, %v", px, err)
					}
				case 1:
					if gh.Table == nil || gh.Table.NRows != 2 || gh.Table.NCols != 1 {
						t.Errorf("cfitsio reads the table as %+v", gh.Table)
					}
				}
			})
		}
	}
}

// TestWritersRefuseUneditedLowerCaseCards: a caller header holding a
// decoded lower-case card that was not edited (a header read from a
// file, or FromCards of its cards) makes every writer return an error
// wrapping header.ErrInvalidKeyword. Truth for why no copy is right
// (header.TestAddCardRefusesUneditedLowerCaseKeyword shows it on the
// source bytes): for 'bzero = 5000', cfitsio does not apply it to the
// pixels and astropy does. Before the fix the writers returned nil and
// copied the bytes: WriteImage with 'bzero = 5000' scaled by 5000 and
// cfitsio (and ReadPixels) read [1 2 3] back as [-4999 -4998 -4997];
// 'bscale = 2' read [1 1 2]; 'naxis1 = 7' made astropy read 7 pixels
// [1 2 3 0 0 0 0]; 'blank = 2' made astropy turn pixel 2 into NaN while
// cfitsio kept it; fitsverify failed on every file.
func TestWritersRefuseUneditedLowerCaseCards(t *testing.T) {
	for _, line := range []string{
		"bzero   =                 5000 / lower",
		"bscale  =                    2 / lower",
		"blank   =                    2 / lower",
		"naxis1  =                    7 / lower",
		"bitpix  =                    8 / lower",
		"naxis2  =                    9 / lower",
		"tform1  = '1J      '           / lower",
		"object  = 'M31'                / lower",
	} {
		for _, w := range lcWriters {
			t.Run(w.name+"/"+strings.Fields(line)[0], func(t *testing.T) {
				src := lcDecoded(t, line)
				for name, h := range map[string]*header.Header{"FromCards": header.FromCards(src)} {
					_, err := lcWrite(t, w, h)
					// The compressed writer refuses any caller BSCALE /
					// BZERO before it copies cards; that error is enough.
					early := w.name == "AppendCompressedImage" && (strings.HasPrefix(line, "bzero") || strings.HasPrefix(line, "bscale"))
					if !errors.Is(err, header.ErrInvalidKeyword) && !(early && err != nil) {
						t.Errorf("%s: write returned %v, want header.ErrInvalidKeyword", name, err)
					}
				}
			})
		}
	}
}
