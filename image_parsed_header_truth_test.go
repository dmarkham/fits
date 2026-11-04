package fits

// Reads and overwrites that take BSCALE/BZERO/BLANK from the fully parsed
// (and possibly edited in memory) header instead of the raw on-disk cards
// (hduRecord.imageKeywords, r.parsed branch). Truth: cfitsio's goldens and
// cfitsio-written data units, as for the raw-card path.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// TestReadsAfterHeaderParseMatchCfitsio re-runs the corpus read checks
// (float64, float32, integer targets, masked) with the full header of
// every HDU parsed first, so BSCALE/BZERO/BLANK are taken from the parsed
// header. Every result must still equal cfitsio's.
func TestReadsAfterHeaderParseMatchCfitsio(t *testing.T) {
	var parsed, withBscale, withBzero, withBlank int
	parseHeaderFirst = func(t *testing.T, s imageSample, i int, img *ImageHDU) {
		if _, err := HeaderOf(img); err != nil {
			t.Errorf("%s HDU %d: HeaderOf: %v", s.name, i, err)
			return
		}
		if img.rec.parsed == nil {
			t.Fatalf("%s HDU %d: HeaderOf succeeded but the parsed header is not cached", s.name, i)
		}
		parsed++
		im := s.golden.HDUs[i].Image
		if float64(im.Bscale) != 1 {
			withBscale++
		}
		if float64(im.Bzero) != 0 {
			withBzero++
		}
		if im.Blank != nil {
			withBlank++
		}
	}
	defer func() { parseHeaderFirst = nil }()
	t.Run("float64", TestReadPixelsFloat64MatchesCfitsio)
	t.Run("float32", TestReadPixelsFloat32MatchesCfitsio)
	t.Run("integer", TestReadPixelsIntegerTargets)
	t.Run("masked", TestReadPixelsMaskedMatchesCfitsio)
	// The parsed path must really have been exercised on scaled and
	// BLANK data (4 sub-tests per HDU).
	if parsed == 0 || withBscale == 0 || withBzero == 0 || withBlank == 0 {
		t.Errorf("parsed-header reads: %d HDUs, %d with BSCALE != 1, %d with BZERO != 0, %d with BLANK; want all > 0",
			parsed, withBscale, withBzero, withBlank)
	}
	t.Logf("parsed-header reads: %d HDU reads, %d with BSCALE != 1, %d with BZERO != 0, %d with BLANK", parsed, withBscale, withBzero, withBlank)
}

// TestOverwritePixelsAfterInMemoryHeaderEdit: in ModeEdit, set BSCALE and
// BZERO on the parsed header (not yet flushed), then OverwritePixels; the
// values must be stored with the edited scaling, so the data unit equals
// the one cfitsio wrote for the same values and cards, and after
// CloseAndFlush the file's BSCALE/BZERO equal cfitsio's.
func TestOverwritePixelsAfterInMemoryHeaderEdit(t *testing.T) {
	dir := t.TempDir()
	ran := 0
	for _, c := range loadCases(t) {
		if c.Bscale == nil || c.Bzero == nil || c.Blank != nil {
			continue
		}
		t.Run(c.Name, func(t *testing.T) {
			path := filepath.Join(dir, c.Name+".fits")
			seed, _ := caseData(t, c, true)
			f, err := Create(path)
			if err != nil {
				t.Fatal(err)
			}
			// Seed without the case's cards: the raw header has no BSCALE
			// and, for the unsigned types, only the convention BZERO.
			if _, err := seed(f, nil, c.Naxes); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}
			e, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			img, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			h, err := HeaderOf(img)
			if err != nil {
				t.Fatal(err)
			}
			if err := h.Set("BSCALE", *c.Bscale, ""); err != nil {
				t.Fatal(err)
			}
			if err := h.Set("BZERO", *c.Bzero, ""); err != nil {
				t.Fatal(err)
			}
			if bscaleOf(t, img) != *c.Bscale || bzeroOf(t, img) != *c.Bzero {
				t.Fatalf("after the in-memory edit BSCALE/BZERO = %v/%v, want %v/%v", bscaleOf(t, img), bzeroOf(t, img), *c.Bscale, *c.Bzero)
			}
			_, over := caseData(t, c, false)
			if err := over(img); err != nil {
				t.Fatalf("OverwritePixels: %v", err)
			}
			// Close does not flush header edits; CloseAndFlush does.
			if err := e.CloseAndFlush(); err != nil {
				t.Fatal(err)
			}
			_, gd := dataUnit(t, path)
			wd := expectedCaseData(t, c)
			if !bytes.Equal(gd, wd) {
				t.Errorf("data stored with the edited BSCALE/BZERO differs from cfitsio's at byte %d\n go:      % x\n cfitsio: % x",
					firstDiff(gd, wd), head(gd, 48), head(wd, 48))
			}
			g := openImage(t, path)
			w := openImage(t, imagesPath(t, "cfwrite", c.Name+".fits"))
			if g.BITPIX() != w.BITPIX() || bscaleOf(t, g) != bscaleOf(t, w) || bzeroOf(t, g) != bzeroOf(t, w) {
				t.Errorf("flushed BITPIX/BSCALE/BZERO = %d/%v/%v, cfitsio %d/%v/%v",
					g.BITPIX(), bscaleOf(t, g), bzeroOf(t, g), w.BITPIX(), bscaleOf(t, w), bzeroOf(t, w))
			}
			ran++
		})
	}
	if ran < 10 {
		t.Errorf("only %d cases exercised the in-memory edit path", ran)
	}
}

// TestUint64ImagesReadWithoutHeaderParse: uint64-convention images (integer
// card BZERO = 9223372036854775808) read exactly as cfitsio's TULONGLONG
// reading, before and after the full header is parsed (the header package
// reads that card as TypeBigInt; a parse failure is a test failure).
func TestUint64ImagesReadWithoutHeaderParse(t *testing.T) {
	samples := []string{"cfwrite/u64_conv.fits", "cfwrite/u64_bzero_user.fits", "astropy_written/u64.fits", "astropy_written/u64_cube.fits"}
	for _, s := range samples {
		t.Run(s, func(t *testing.T) {
			p := imagesPath(t, filepath.FromSlash(s))
			if _, err := os.Stat(p); err != nil {
				t.Fatal(err)
			}
			wp := loadPixGolden(t, s).HDUs[0]
			if wp.U64Error != "" || wp.U64.Len() == 0 {
				t.Fatalf("golden has no TULONGLONG reading: %s", wp.U64Error)
			}
			want := wp.U64.Values(t)
			img := openImage(t, p)
			check := func(when string) {
				got, err := ReadPixels[uint64](img)
				if err != nil {
					t.Fatalf("%s: ReadPixels[uint64]: %v", when, err)
				}
				if len(got) != len(want) {
					t.Fatalf("%s: %d pixels, cfitsio %d", when, len(got), len(want))
				}
				for k, w := range want {
					if got[k] != w.V {
						t.Fatalf("%s: pixel %d = %d, cfitsio TULONGLONG %d", when, k, got[k], w.V)
					}
				}
			}
			check("before HeaderOf")
			if _, err := HeaderOf(img); err != nil {
				t.Fatalf("HeaderOf: %v", err)
			}
			check("after HeaderOf")
		})
	}
}
