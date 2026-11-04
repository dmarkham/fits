package fits

import (
	"errors"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCompressedNullOutsideZBITPIX: testdata/comp_bad/i16_gzip1_w4_nullout
// .fits holds 4-byte GZIP tiles on a 16-bit image whose pixel 0 stores
// 37232, with BLANK = 37232. cfitsio's fits_read_imgnull flags it null
// (the golden) while astropy 8.0.2 reads the BLANK against the int16 data
// and flags nothing (pixel 0 = -28304): a BLANK outside the ZBITPIX range
// is refused by the scaling-card rule (scalingCards), so every pixel read
// returns ErrInvalidScaling, with the HDU's structure still readable.
func TestCompressedNullOutsideZBITPIX(t *testing.T) {
	const name = "i16_gzip1_w4_nullout.fits"
	g := loadGoldenFile(t, "testdata/comp_bad/golden/"+name+".json")
	if g.NHDU != 2 || g.HDUs[1].Image == nil || len(g.HDUs[1].Errors) != 0 {
		t.Fatalf("cfitsio golden: %d HDUs, errors %v", g.NHDU, g.HDUs[1].Errors)
	}
	gh := &g.HDUs[1]
	if b, _ := goldenCardValue(gh, "BLANK"); b != "37232" {
		t.Fatalf("golden BLANK %q, want 37232", b)
	}
	if w := inflatedTileWidth(t, filepath.Join("comp_bad", name)); w != 4 {
		t.Fatalf("tile 1 inflates to %d bytes per pixel, want 4", w)
	}
	if want := gh.Image.Phys.Values(t); want[0].Valid || gh.Image.Phys.NullCount() != 1 {
		t.Fatalf("cfitsio golden: pixel 0 valid %v, %d nulls; want pixel 0 the only null", want[0].Valid, gh.Image.Phys.NullCount())
	}
	f, h := openCompressed(t, filepath.Join("testdata", "comp_bad", name))
	defer f.Close()
	if h.BITPIX() != 16 || len(h.Shape()) != 2 {
		t.Errorf("BITPIX %d Shape %v: the structure must stay readable", h.BITPIX(), h.Shape())
	}
	refused := func(what string, err error) {
		t.Helper()
		if !errors.Is(err, ErrInvalidScaling) || !strings.Contains(err.Error(), "BLANK = 37232") {
			t.Errorf("%s: %v, want ErrInvalidScaling naming BLANK = 37232", what, err)
		}
	}
	_, _, err := ReadPixelsCompressedMasked[float64](h)
	refused("ReadPixelsCompressedMasked[float64]", err)
	_, _, err = ReadPixelsCompressedMasked[int16](h)
	refused("ReadPixelsCompressedMasked[int16]", err)
	_, err = ReadPixelsCompressed[float64](h)
	refused("ReadPixelsCompressed[float64]", err)
	_, err = ReadFloat32(h)
	refused("ReadFloat32", err)

	// Two such null pixels in two tiles (comp_bad/i16_gzip1_w4_nullout2:
	// per-tile ZBLANK column 70000, 70001 and pixel 0 of each tile holding
	// it): cfitsio flags both (golden), funpack fails (412); the masked
	// read equals the golden and the unmasked read names the first one.
	const name2 = "i16_gzip1_w4_nullout2.fits"
	g2 := loadGoldenFile(t, "testdata/comp_bad/golden/"+name2+".json")
	if g2.NHDU != 2 || g2.HDUs[1].Image == nil || len(g2.HDUs[1].Errors) != 0 {
		t.Fatalf("%s: cfitsio golden: %d HDUs", name2, g2.NHDU)
	}
	want2 := g2.HDUs[1].Image.Phys.Values(t)
	if want2[0].Valid || want2[40].Valid || g2.HDUs[1].Image.Phys.NullCount() != 2 {
		t.Fatalf("%s: cfitsio golden: want pixels 0 and 40 the only nulls", name2)
	}
	f2, h2 := openCompressed(t, filepath.Join("testdata", "comp_bad", name2))
	defer f2.Close()
	v2, m2, err := ReadPixelsCompressedMasked[int16](h2)
	if err != nil {
		t.Fatalf("%s: ReadPixelsCompressedMasked[int16]: %v", name2, err)
	}
	for k, w := range want2 {
		if m2[k] != w.Valid || (w.Valid && float64(v2[k]) != w.V) {
			t.Fatalf("%s pixel %d: %d valid %v; cfitsio %v valid %v", name2, k, v2[k], m2[k], w.V, w.Valid)
		}
	}
	const msg2 = "tile 0: null tile value 70000 does not fit ZBITPIX 16"
	if _, err := ReadPixelsCompressed[float64](h2); err == nil || !strings.Contains(err.Error(), msg2) {
		t.Errorf("%s: ReadPixelsCompressed[float64]: err %v, want one mentioning %q", name2, err, msg2)
	}
}

// TestCompressedBlankRangeAstropy: astropy 8.0.2 refuses to decode the
// comp_bad files whose BLANK does not fit a 32-bit int (OverflowError
// "BLANK value ... is too large"), the second truth behind the library's
// refusal (TestCompressedBadSamplesRefused; cfitsio ignores such a BLANK).
func TestCompressedBlankRangeAstropy(t *testing.T) {
	const script = `import sys
import astropy.io.fits as fits
try:
    fits.open(sys.argv[1])[1].data
except OverflowError as e:
    print("OverflowError:", e)
else:
    print("decoded")
`
	for file, value := range map[string]string{
		"i32_rice_blank2p32.fits":  "4294967296",
		"i64_gzip1_blank2p40.fits": "1099511627783",
	} {
		p := filepath.Join("testdata", "comp_bad", file)
		r := testutil.RunOK(t, testutil.ToolAstropy, testutil.Call{Args: []string{"-c", script, p}, Inputs: []string{p}})
		out := append(r.Stdout, r.Stderr...)
		if want := "OverflowError: BLANK value " + value + " is too large"; strings.TrimSpace(string(out)) != want {
			t.Errorf("%s: astropy says %q, want %q", file, strings.TrimSpace(string(out)), want)
		}
	}
}

// TestCompressedUnmaskedNullsAreStoredValues: the unmasked reads return a
// null pixel of an integer image as its stored value, scaled (the doc of
// ReadPixelsCompressed). The truth is the funpack twin, which writes the
// stored values to a plain image: its value where cfitsio's plain-image
// read of the twin has one, else (the twin's BLANK flags the pixel too)
// the twin's BLANK scaled by its BSCALE/BZERO, which is the stored value
// by definition. Every integer corpus case with null pixels is checked,
// except the fallback rows where cfitsio's own fits_read_img differs
// (rawRows).
func TestCompressedUnmaskedNullsAreStoredValues(t *testing.T) {
	checked := 0
	for _, c := range fzCorpus {
		if c.zbitpix < 0 || c.i64Truth {
			continue
		}
		_, null := caseTruth(t, c, false)
		n := 0
		for _, x := range null {
			if x {
				n++
			}
		}
		if n == 0 {
			continue
		}
		t.Run(c.name, func(t *testing.T) {
			if checkCorpusScalingRefused(t, c) {
				return
			}
			th, tv, tn := cfrefPixels(t, twinFile(t, c.name), false)
			gh, _, _ := cfrefPixels(t, "comp_fz_"+c.name+".fits", false)
			raw := map[int]bool{}
			nx := int(gh.Image.Naxes[0])
			for _, r := range c.rawRows {
				for i := r * nx; i < (r+1)*nx; i++ {
					raw[i] = true
				}
			}
			f, h := openCompressed(t, filepathTestdata(t, "comp_fz_"+c.name+".fits"))
			defer f.Close()
			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatalf("ReadPixelsCompressed[float64]: %v", err)
			}
			for k, isNull := range null {
				if !isNull || raw[k] {
					continue
				}
				want := tv[k]
				if tn[k] {
					b, ok := goldenCardValue(th, "BLANK")
					bv, err := strconv.ParseInt(b, 10, 64)
					if !ok || err != nil {
						t.Fatalf("pixel %d: twin null but twin BLANK %q", k, b)
					}
					want = float64(th.Image.Bzero) + float64(th.Image.Bscale)*float64(bv)
				}
				if got[k] != want {
					t.Fatalf("null pixel %d: unmasked %v, stored value per the funpack twin %v", k, got[k], want)
				}
			}
		})
		checked++
	}
	if checked < 10 {
		t.Fatalf("only %d integer corpus cases with nulls", checked)
	}
}
