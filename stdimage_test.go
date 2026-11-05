package fits

import (
	"bytes"
	"encoding/json"
	"errors"
	"image"
	"image/color"
	"math"
	"path/filepath"
	"sort"
	"testing"
)

// Tests for AsImage, FloatGray and the image.Decode registration against
// testdata/imaging/stdimage_expected.json (pixel values read by astropy;
// see testdata/imaging/PROVENANCE.md).

type stdimageExpected struct {
	Bitpix int      `json:"bitpix"`
	Naxis  int      `json:"naxis"`
	Kind   string   `json:"kind"`
	W      int      `json:"w"`
	H      int      `json:"h"`
	Values pixArray `json:"values"` // a digest past 1024 values (NaN bits as 0x7fc00000)
}

func loadStdimageExpected(t *testing.T) map[string]stdimageExpected {
	t.Helper()
	var g struct {
		Files map[string]stdimageExpected `json:"files"`
	}
	if err := json.Unmarshal(mustReadFile(t, imagingPath("stdimage_expected.json")), &g); err != nil {
		t.Fatal(err)
	}
	if len(g.Files) < 35 {
		t.Fatalf("stdimage_expected.json lists %d files, want >= 35", len(g.Files))
	}
	return g.Files
}

// checkImage asserts img against the expected kind and values. Pixel
// (x, y) of the image.Image is FITS pixel (x, y) (0-based, x along NAXIS1):
// values[y*w + x].
func checkImage(t *testing.T, img image.Image, e stdimageExpected) {
	t.Helper()
	b := img.Bounds()
	if b != image.Rect(0, 0, e.W, e.H) {
		t.Fatalf("bounds %v, want %dx%d", b, e.W, e.H)
	}
	float := e.Kind == "float" || e.Kind == "floatphys"
	got := make([]int64, 0, e.W*e.H)
	for y := 0; y < e.H; y++ {
		for x := 0; x < e.W; x++ {
			switch e.Kind {
			case "gray":
				g, ok := img.(*image.Gray)
				if !ok {
					t.Fatalf("got %T, want *image.Gray", img)
				}
				got = append(got, int64(g.GrayAt(x, y).Y))
			case "gray16":
				g, ok := img.(*image.Gray16)
				if !ok {
					t.Fatalf("got %T, want *image.Gray16", img)
				}
				got = append(got, int64(g.Gray16At(x, y).Y))
			case "float", "floatphys":
				fg, ok := img.(*FloatGray)
				if !ok {
					t.Fatalf("got %T, want *FloatGray", img)
				}
				got = append(got, f32Bits([]float32{fg.At(x, y).(FloatColor).V}, true)[0])
			default:
				t.Fatalf("kind %q", e.Kind)
			}
		}
	}
	if !e.Values.Digested() && float {
		// Explicit float values keep astropy's NaN bits; any NaN matches.
		want := e.Values.Explicit(t)
		if len(want) != len(got) {
			t.Fatalf("%d expected values for %d pixels", len(want), len(got))
		}
		for i, w := range want {
			if w != nil && math.IsNaN(float64(math.Float32frombits(uint32(*w)))) && got[i] == 0x7fc00000 {
				got[i] = *w
			}
		}
	}
	e.Values.Check(t, "AsImage pixels (x fastest)", got)
}

// TestAsImage_Matrix: every uncompressed sample through AsImage — the
// lossless cases return the exact stored/physical values, everything that
// cannot be represented without loss returns ErrNotImageCompatible.
func TestAsImage_Matrix(t *testing.T) {
	exp := loadStdimageExpected(t)
	names := make([]string, 0, len(exp))
	for n := range exp {
		names = append(names, n)
	}
	sort.Strings(names)
	kinds := map[string]int{}
	for _, name := range names {
		e := exp[name]
		kinds[e.Kind]++
		t.Run(name, func(t *testing.T) {
			f, err := Open(imagingPath(name))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			p, err := f.Primary()
			if err != nil {
				t.Fatal(err)
			}
			img, err := p.AsImage()
			if e.Kind == "reject" {
				if !errors.Is(err, ErrNotImageCompatible) {
					t.Fatalf("AsImage = %T, %v; want ErrNotImageCompatible (BITPIX=%d NAXIS=%d)", img, err, e.Bitpix, e.Naxis)
				}
				return
			}
			if err != nil {
				t.Fatalf("AsImage: %v", err)
			}
			checkImage(t, img, e)
		})
	}
	for _, k := range []string{"gray", "gray16", "float", "floatphys", "reject"} {
		if kinds[k] == 0 {
			t.Errorf("no sample of kind %q", k)
		}
	}
}

// TestImageDecode_Registered: the blank import of this package registers
// "fits" with the image package; image.Decode / DecodeConfig on the bytes
// of a FITS file return the same image as AsImage and matching config, and
// the rejected samples fail through image.Decode as well.
func TestImageDecode_Registered(t *testing.T) {
	exp := loadStdimageExpected(t)
	for _, name := range []string{"u8_all.fits", "u16_all.fits", "asi533_u16_crop.fits", "f32_2d.fits", "f64_exact_2d.fits", "f64_2d.fits", "i16_all.fits", "u16_blank_2d.fits", "u16_rgb.fits", "i32.fits"} {
		e := exp[name]
		t.Run(name, func(t *testing.T) {
			data := mustReadFile(t, imagingPath(name))
			img, format, err := image.Decode(bytes.NewReader(data))
			cfg, cfgFormat, cfgErr := image.DecodeConfig(bytes.NewReader(data))
			if e.Kind == "reject" {
				if err == nil {
					t.Fatalf("image.Decode accepted it as %T", img)
				}
				if format != "fits" && !errors.Is(err, ErrNotImageCompatible) {
					t.Fatalf("image.Decode: format %q err %v", format, err)
				}
				if name != "u16_blank_2d.fits" && name != "f64_2d.fits" && cfgErr == nil {
					// DecodeConfig reads the header only, so it cannot see
					// BLANK pixels or inexact float values; everything else
					// must agree with Decode.
					t.Fatalf("DecodeConfig accepted a file AsImage rejects: %+v", cfg)
				}
				return
			}
			if err != nil || format != "fits" {
				t.Fatalf("image.Decode: format %q err %v", format, err)
			}
			checkImage(t, img, e)
			if cfgErr != nil || cfgFormat != "fits" {
				t.Fatalf("image.DecodeConfig: format %q err %v", cfgFormat, cfgErr)
			}
			if cfg.Width != e.W || cfg.Height != e.H || cfg.ColorModel != img.ColorModel() {
				t.Errorf("DecodeConfig %dx%d model %T, image %dx%d model %T", cfg.Width, cfg.Height, cfg.ColorModel, e.W, e.H, img.ColorModel())
			}
		})
	}
	// A tile-compressed file has an empty (NAXIS=0) primary HDU; the image
	// is in HDU 2 (cfitsio: testdata/imaging/cfnull/golden image_hdu = 2).
	// Decode uses the primary only, so it must reject the file, through
	// image.Decode too, rather than decode the compressed table.
	for _, name := range []string{"comp_asi533_u16_crop.fits.fz", "comp_u8_blank.fits.fz"} {
		data := mustReadFile(t, imagingPath(name))
		img, format, err := image.Decode(bytes.NewReader(data))
		if err == nil || !errors.Is(err, ErrNotImageCompatible) || format != "fits" {
			t.Errorf("%s: image.Decode = %T, %q, %v; want format fits and ErrNotImageCompatible", name, img, format, err)
		}
		if _, _, err := image.DecodeConfig(bytes.NewReader(data)); !errors.Is(err, ErrNotImageCompatible) {
			t.Errorf("%s: image.DecodeConfig err %v; want ErrNotImageCompatible", name, err)
		}
		var cf struct {
			ImageHDU int `json:"image_hdu"`
		}
		if err := json.Unmarshal(mustReadFile(t, imagingPath(filepath.Join("cfnull", "golden", name+".json"))), &cf); err != nil || cf.ImageHDU != 2 {
			t.Errorf("%s: cfitsio golden image_hdu = %d (%v), want 2", name, cf.ImageHDU, err)
		}
	}
	if _, format, err := image.Decode(bytes.NewReader([]byte("\x89PNG\r\n\x1a\nnot a fits file"))); format == "fits" || err == nil {
		t.Errorf("non-FITS input decoded as %q (%v)", format, err)
	}
}

// TestAsImage_RowOrientation pins the documented orientation: image row
// y is FITS row y+1 (the first stored row is y=0, the top of the
// image.Image), so an encoder shows the frame flipped vertically relative
// to the FITS convention (y up).
func TestAsImage_RowOrientation(t *testing.T) {
	e := loadStdimageExpected(t)["u16_orient.fits"]
	f, err := Open(imagingPath("u16_orient.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	p, _ := f.Primary()
	img, err := p.AsImage()
	if err != nil {
		t.Fatal(err)
	}
	g := img.(*image.Gray16)
	for y := 0; y < e.H; y++ {
		for x := 0; x < e.W; x++ {
			if want := uint16(1000 * (4*y + x)); g.Gray16At(x, y).Y != want {
				t.Errorf("(%d,%d) = %d, want %d", x, y, g.Gray16At(x, y).Y, want)
			}
		}
	}
}

// TestFloatColor_RGBA pins the documented placeholder mapping, including
// NaN and out-of-range inputs.
func TestFloatColor_RGBA(t *testing.T) {
	cases := []struct {
		v    float32
		want uint32
	}{
		{0, 0}, {1, 0xffff}, {0.5, 0x7fff}, {-3, 0}, {7, 0xffff},
		// just above 1, and inside (1, 2]: clamped to 1, never above 0xffff
		{math.Nextafter32(1, 2), 0xffff}, {1.5, 0xffff}, {2, 0xffff},
		{float32(math.Inf(1)), 0xffff}, {float32(math.Inf(-1)), 0}, {float32(math.NaN()), 0},
	}
	for _, c := range cases {
		r, g, b, a := FloatColor{V: c.v}.RGBA()
		if r != c.want || g != c.want || b != c.want || a != 0xffff {
			t.Errorf("FloatColor{%v}.RGBA() = %#x,%#x,%#x,%#x; want %#x x3, 0xffff", c.v, r, g, b, a, c.want)
		}
	}
	if got := FloatGrayModel.Convert(color.Gray16{Y: 0xffff}).(FloatColor).V; got != 1 {
		t.Errorf("Convert(Gray16 max) = %v, want 1", got)
	}
	fc := FloatColor{V: 123.5}
	if FloatGrayModel.Convert(fc) != fc {
		t.Errorf("Convert must pass a FloatColor through unchanged")
	}
	fg := &FloatGray{Rect: image.Rect(0, 0, 2, 1), Stride: 2, Pix: []float32{1.5, -2}}
	if v := fg.At(1, 0).(FloatColor).V; v != -2 {
		t.Errorf("At(1,0) = %v", v)
	}
	if v := fg.At(2, 0).(FloatColor).V; v != 0 {
		t.Errorf("out-of-bounds At = %v, want 0", v)
	}
}
