package fits

import (
	"math"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCompressedFixtures checks the astropy-written fixtures made by
// testdata/generate_compressed.py (astropy 8.0.2; see
// testdata/COMP_FZ_PROVENANCE.md). Each compressed fixture is decoded and
// compared BIT-EXACT against what libcfitsio decodes from the same file
// (testdata/cfref/golden/comp_<name>.fits.json; float32 images compare to
// float32 of cfitsio's value). Lossless fixtures must in addition equal
// their comp_raw_ twin exactly.
//
// Quantized float fixtures are NOT compared to their raw twin for
// correctness: that holds the original floats, and a tolerance against it
// cannot tell a correct dequantizer from one that ignores the dither (the
// earlier version of this test passed with dithering disabled).
//
// NOCOMPRESS has no cfitsio truth: cfitsio 4.7.0 (and funpack) fail on it
// with "unknown compression algorithm" (recorded in its golden), so only
// the raw-twin comparison applies: for integers the raw twin is the input
// image, for floats (astropy quantizes NOCOMPRESS floats: ZSCALE/ZZERO and
// GZIP_COMPRESSED_DATA columns, SUBTRACTIVE_DITHER_1) it is astropy
// 8.0.2's own decode of the file.
//
// Raw twins that would be byte-identical are written once (raw names the
// fixture whose comp_raw_ file holds the image; generate_compressed.py
// raw_of).
func TestCompressedFixtures(t *testing.T) {
	cases := []struct {
		name      string // comp_<name>.fits / comp_raw_<raw or name>.fits
		raw       string
		algo      string
		bitpix    int
		quantiz   string // ZQUANTIZ the fixture must carry ("" for integer)
		lossless  bool   // decoded pixels must equal the raw twin
		noCfitsio bool   // cfitsio cannot decode it (NOCOMPRESS)
	}{
		{"rice_i16", "", "RICE_1", 16, "", true, false},
		{"rice_i32", "", "RICE_1", 32, "", true, false},
		{"gzip1_i16", "rice_i16", "GZIP_1", 16, "", true, false},
		{"gzip1_i32", "rice_i32", "GZIP_1", 32, "", true, false},
		{"gzip2_i16", "rice_i16", "GZIP_2", 16, "", true, false},
		{"gzip2_i32", "rice_i32", "GZIP_2", 32, "", true, false},
		{"nocompress_i16", "rice_i16", "NOCOMPRESS", 16, "", true, true},
		{"plio_mask", "", "PLIO_1", 16, "", true, false},
		{"hcompress_i16", "rice_i16", "HCOMPRESS_1", 16, "", true, false},
		// astropy writes no ZQUANTIZ for NO_DITHER; absent means NO_DITHER
		// (cfitsio imcomp_get_compressed_image_par).
		{"hcompress_f32", "rice_f32_nodither", "HCOMPRESS_1", -32, "", false, false},
		{"rice_f32_nodither", "", "RICE_1", -32, "", false, false},
		{"rice_f32_dither1", "rice_f32_nodither", "RICE_1", -32, "SUBTRACTIVE_DITHER_1", false, false},
		{"rice_f32_dither2", "", "RICE_1", -32, "SUBTRACTIVE_DITHER_2", false, false},
		{"gzip1_f32", "rice_f32_nodither", "GZIP_1", -32, "", false, false},
		{"gzip2_f32", "rice_f32_nodither", "GZIP_2", -32, "SUBTRACTIVE_DITHER_1", false, false},
		{"rice_multitile", "rice_i16", "RICE_1", 16, "", true, false},
		{"rice_onetile", "rice_i16", "RICE_1", 16, "", true, false},
		{"rice_odd", "", "RICE_1", 16, "", true, false},
		{"nocompress_u8", "", "NOCOMPRESS", 8, "", true, true},
		{"nocompress_i32", "", "NOCOMPRESS", 32, "", true, true},
		{"nocompress_i64", "", "NOCOMPRESS", 64, "", true, true},
		{"nocompress_f32", "", "NOCOMPRESS", -32, "SUBTRACTIVE_DITHER_1", true, true},
		{"nocompress_f64", "", "NOCOMPRESS", -64, "SUBTRACTIVE_DITHER_1", true, true},
	}
	raws := map[string]bool{}
	for _, c := range cases {
		if c.raw != "" {
			raws[c.raw] = true
		} else {
			raws[c.name] = true
		}
		t.Run(c.name, func(t *testing.T) {
			compName := "comp_" + c.name + ".fits"
			comp, err := Open(testutil.TestdataPath(t, compName))
			if err != nil {
				t.Fatalf("open compressed: %v", err)
			}
			defer comp.Close()
			cHDU, err := comp.HDU(1)
			if err != nil {
				t.Fatalf("HDU(1): %v", err)
			}
			cimg, ok := cHDU.(*CompressedImageHDU)
			if !ok {
				t.Fatalf("HDU(1) is %T, expected *CompressedImageHDU", cHDU)
			}
			if got := cimg.CompressionType(); got != c.algo {
				t.Errorf("CompressionType = %q, want %q", got, c.algo)
			}
			if cimg.BITPIX() != c.bitpix {
				t.Fatalf("BITPIX = %d, want %d", cimg.BITPIX(), c.bitpix)
			}
			if q, _ := cimg.Header().String("ZQUANTIZ"); q != c.quantiz {
				t.Fatalf("fixture ZQUANTIZ = %q, want %q (generator selects the method explicitly)", q, c.quantiz)
			}
			got, err := ReadPixelsCompressed[float64](cimg)
			if err != nil {
				t.Fatalf("decompress: %v", err)
			}
			if c.bitpix == -32 {
				g32, err := ReadPixelsCompressed[float32](cimg)
				if err != nil {
					t.Fatalf("decompress float32: %v", err)
				}
				for i := range g32 {
					if !sameFloat(float64(g32[i]), got[i]) {
						t.Fatalf("pixel %d: float32 read %v, float64 read %v", i, g32[i], got[i])
					}
				}
			}

			if !c.noCfitsio {
				_, want, _ := cfrefPixels(t, compName, c.bitpix == -32)
				if d := diffPix(got, want, nil); d.n > 0 {
					t.Errorf("vs cfitsio decode: %s", d)
				}
			} else {
				g := testutil.LoadCfrefGolden(t, compName)
				if len(g.HDUs) < 2 || len(g.HDUs[1].Errors) == 0 {
					t.Errorf("golden no longer records the cfitsio failure on %s; compare against it", compName)
				}
			}

			rawName := c.raw
			if rawName == "" {
				rawName = c.name
			}
			raw, err := Open(testutil.TestdataPath(t, "comp_raw_"+rawName+".fits"))
			if err != nil {
				t.Fatalf("open raw: %v", err)
			}
			defer raw.Close()
			rHDU, err := raw.HDU(1)
			if err != nil {
				t.Fatalf("HDU(1) raw: %v", err)
			}
			rimg := rHDU.(*ImageHDU)
			cs, rs := cimg.Shape(), rimg.Shape()
			if len(cs) != len(rs) {
				t.Fatalf("NAXIS mismatch: comp=%d raw=%d", len(cs), len(rs))
			}
			for i := range cs {
				if cs[i] != rs[i] {
					t.Fatalf("NAXIS%d mismatch: comp=%d raw=%d", i+1, cs[i], rs[i])
				}
			}
			want, err := ReadPixels[float64](rimg)
			if err != nil {
				t.Fatal(err)
			}
			if c.lossless {
				if d := diffPix(got, want, nil); d.n > 0 {
					t.Errorf("lossless fixture vs raw twin: %s", d)
				}
				if c.bitpix == 64 {
					// past 2^53 float64 cannot tell values apart: compare the
					// integers with cfitsio's read of the raw twin
					gi, err := ReadPixelsCompressed[int64](cimg)
					if err != nil {
						t.Fatal(err)
					}
					wi := testutil.LoadCfrefGolden(t, "comp_raw_"+rawName+".fits").HDUs[1].Image.I64.Values(t)
					if len(wi) != len(gi) {
						t.Fatalf("%d pixels, raw twin golden %d", len(gi), len(wi))
					}
					for k := range gi {
						if !wi[k].Valid || gi[k] != wi[k].V {
							t.Fatalf("pixel %d: %d, cfitsio read of the raw twin %+v", k, gi[k], wi[k])
						}
					}
				}
				return
			}
			// Sanity only (the exact check is the cfitsio one above): a
			// quantized pixel stays within one quantization step (the
			// largest tile ZSCALE) of the original.
			zs, err := ReadColumn[float64](cimg.tbl, findCompColumn(mustCols(t, cimg), "ZSCALE"))
			if err != nil {
				t.Fatal(err)
			}
			step := 0.0
			for _, z := range zs {
				step = math.Max(step, z)
			}
			for i := range want {
				if d := math.Abs(got[i] - want[i]); d > step {
					t.Fatalf("pixel %d: |decoded-original| %g > ZSCALE %g", i, d, step)
				}
			}
		})
	}
	files, _ := filepath.Glob(filepath.Join("testdata", "comp_raw_*.fits"))
	for _, f := range files {
		if !raws[strings.TrimSuffix(strings.TrimPrefix(filepath.Base(f), "comp_raw_"), ".fits")] {
			t.Errorf("%s is no fixture's raw twin", filepath.Base(f))
		}
	}
	if len(files) != len(raws) {
		t.Errorf("%d comp_raw_ files, the cases name %d", len(files), len(raws))
	}
}

// TestCompressedFixtureSetComplete: every comp_raw_ twin committed is one
// a TestCompressedFixtures case reads, and every astropy fixture has a case.
func TestCompressedFixtureSetComplete(t *testing.T) {
	raws, _ := filepath.Glob(filepath.Join("testdata", "comp_raw_*.fits"))
	comps, _ := filepath.Glob(filepath.Join("testdata", "comp_*.fits"))
	n := 0
	for _, c := range comps {
		b := filepath.Base(c)
		if !strings.HasPrefix(b, "comp_fz") && !strings.HasPrefix(b, "comp_raw_") {
			n++
		}
	}
	if len(raws) != 11 || n != 23 {
		t.Fatalf("%d comp_raw_ twins and %d astropy fixtures, want 11 and 23 (the TestCompressedFixtures cases)", len(raws), n)
	}
}

func mustCols(t *testing.T, h *CompressedImageHDU) []Column {
	t.Helper()
	cols, err := h.tbl.Columns()
	if err != nil {
		t.Fatal(err)
	}
	return cols
}
