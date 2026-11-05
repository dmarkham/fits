package fits

import (
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// TestWritersWriteEditedCardCopies: a caller reads a header, edits copies
// of its cards (Cards() then FromCards) and hands the result to a writer.
// The written file must carry the edited cards, so cfitsio reads what the
// in-memory header says, and the image writer's scaling (which uses the
// header's getters) must agree with the BZERO it writes. Before the fix
// every writer wrote the old bytes of the edited cards: WriteImage wrote
// BZERO = 100 and OBJECT = 'M31' while scaling the data with BZERO 5000,
// and cfitsio read the pixels [1 2 3] back as [-4899 -4898 -4897].
func TestWritersWriteEditedCardCopies(t *testing.T) {
	src := lcDecoded(t,
		"BZERO   =                  100 / offset",
		"OBJECT  = 'M31'",
		"EXPTIME =                  5.0 / s",
	)
	edited := append([]header.Card(nil), header.FromCards(src).Cards()...)
	for i := range edited {
		switch edited[i].Key {
		case "BZERO":
			edited[i].Value = int64(5000)
		case "OBJECT":
			edited[i].Value = "M33"
		case "EXPTIME":
			edited[i].Value, edited[i].Comment = 7.5, "edited"
		}
	}

	// The image writer: cards and pixels through cfitsio.
	p := filepath.Join(t.TempDir(), "img.fits")
	f, err := Create(p)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage[int16](f, header.FromCards(edited), []int64{3}, []int16{1, 2, 3}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	cf := testutil.CfrefDump(t, p)
	gh := &cf.HDUs[0]
	checkEditedCards(t, "WriteImage", gh)
	if gh.Image == nil {
		t.Fatal("cfitsio sees no image")
	}
	if float64(gh.Image.Bzero) != 5000 {
		t.Errorf("cfitsio BZERO %v, want the edited 5000", float64(gh.Image.Bzero))
	}
	want := []testutil.NullFloat{{V: 1, Valid: true}, {V: 2, Valid: true}, {V: 3, Valid: true}}
	gh.Image.Phys.Check(t, "WriteImage pixels (cfitsio TDOUBLE)", want)

	// The table and compressed writers: cards through cfitsio.
	for _, w := range lcWriters[1:] {
		t.Run(w.name, func(t *testing.T) {
			hdr := header.FromCards(edited)
			if err := hdr.Delete("BZERO"); err != nil { // BZERO is an image keyword
				t.Fatal(err)
			}
			p, werr := lcWrite(t, w, hdr)
			if werr != nil {
				t.Fatal(werr)
			}
			cf := testutil.CfrefDump(t, p)
			if len(cf.HDUs) != 2 {
				t.Fatalf("cfitsio sees %d HDUs", len(cf.HDUs))
			}
			checkEditedCards(t, w.name, &cf.HDUs[1])
		})
	}
}

func checkEditedCards(t *testing.T, where string, gh *testutil.CfrefHDU) {
	t.Helper()
	if c := gh.Card("OBJECT"); c == nil || c.Str == nil || string(*c.Str) != "M33" {
		t.Errorf("%s: cfitsio reads OBJECT %+v, want the edited 'M33'", where, c)
	}
	if c := gh.Card("EXPTIME"); c == nil || c.Float == nil || float64(*c.Float) != 7.5 || string(c.Comment) != "edited" {
		t.Errorf("%s: cfitsio reads EXPTIME %+v, want 7.5 / edited", where, c)
	}
}

// TestWriteImageEditedNegZeroBZERO: the source header has BZERO = -0 (an
// integer card cfitsio's fits_read_key_dbl reads as -0). A copy edited to
// 5000 must be written as 5000 AND used as 5000 for the scaling, so
// cfitsio reads the pixels [1 2 3] back; a copy whose comment alone was
// edited still holds -0. Before the fix the header's Float(BZERO) kept
// returning -0 for the edited copy: WriteImage wrote BZERO = 5000 but
// scaled with 0, and cfitsio read [5001 5002 5003] without any error.
func TestWriteImageEditedNegZeroBZERO(t *testing.T) {
	for _, tc := range []struct {
		name  string
		edit  func(c *header.Card)
		bzero float64
	}{
		{"value edited", func(c *header.Card) { c.Value = int64(5000) }, 5000},
		{"comment edited", func(c *header.Card) { c.Comment = "edited" }, 0},
	} {
		t.Run(tc.name, func(t *testing.T) {
			edited := append([]header.Card(nil), lcDecoded(t, "BZERO   =                   -0 / offset")...)
			tc.edit(&edited[0])
			p := filepath.Join(t.TempDir(), "img.fits")
			f, err := Create(p)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := WriteImage[int16](f, header.FromCards(edited), []int64{3}, []int16{1, 2, 3}); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}
			gh := &testutil.CfrefDump(t, p).HDUs[0]
			if gh.Image == nil {
				t.Fatal("cfitsio sees no image")
			}
			if float64(gh.Image.Bzero) != tc.bzero {
				t.Errorf("cfitsio BZERO %v, want %v", float64(gh.Image.Bzero), tc.bzero)
			}
			want := []testutil.NullFloat{{V: 1, Valid: true}, {V: 2, Valid: true}, {V: 3, Valid: true}}
			gh.Image.Phys.Check(t, "pixels (cfitsio TDOUBLE)", want)
		})
	}
}
