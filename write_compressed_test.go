package fits

import (
	"errors"
	"math"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/header"
)

// TestWriteCompressedRICE: create a compressed file via AppendCompressedImage
// with RICE_1, reopen it via Open, decompress via ReadPixelsCompressed,
// verify pixel-for-pixel equality.
func TestWriteCompressedRICE(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "write_rice.fits")

	// Write.
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	// Primary placeholder.
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	// 64x48 int16 image with a predictable pattern.
	shape := []int64{64, 48}
	pixels := make([]int16, 64*48)
	for i := range pixels {
		pixels[i] = int16(i % 1000)
	}
	if _, err := AppendCompressedImage(f, nil, shape, pixels, CompressOptions{
		Algorithm: compress.RICE1,
	}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	// Read back.
	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()

	if fr.NumHDU() != 2 {
		t.Fatalf("NumHDU = %d, want 2", fr.NumHDU())
	}
	h, _ := fr.HDU(1)
	cimg, ok := h.(*CompressedImageHDU)
	if !ok {
		t.Fatalf("HDU(1) is %T, want *CompressedImageHDU", h)
	}
	if cimg.CompressionType() != "RICE_1" {
		t.Fatalf("algorithm %q", cimg.CompressionType())
	}
	if cimg.BITPIX() != 16 {
		t.Fatalf("BITPIX %d", cimg.BITPIX())
	}
	sh := cimg.Shape()
	if len(sh) != 2 || sh[0] != 64 || sh[1] != 48 {
		t.Fatalf("shape %v", sh)
	}

	got, err := ReadPixelsCompressed[int16](cimg)
	if err != nil {
		t.Fatalf("decompress: %v", err)
	}
	if len(got) != len(pixels) {
		t.Fatalf("len %d vs %d", len(got), len(pixels))
	}
	for i := range pixels {
		if got[i] != pixels[i] {
			t.Fatalf("pixel %d: got %d want %d", i, got[i], pixels[i])
		}
	}
}

// TestWriteCompressedGZIP1: same but with GZIP_1.
func TestWriteCompressedGZIP1(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "write_gzip1.fits")
	writeReadCompressedInt16Test(t, path, compress.GZIP1, []int64{32, 32})
}

// TestWriteCompressedGZIP2: GZIP_2 with byte shuffle.
func TestWriteCompressedGZIP2(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "write_gzip2.fits")
	writeReadCompressedInt16Test(t, path, compress.GZIP2, []int64{40, 40})
}

// TestWriteCompressedRefusals: inputs cfitsio cannot write or read are
// refused with an error instead of producing a file other tools reject:
//
//   - NOCOMPRESS: funpack/cfitsio 4.7.0 fail on ZCMPTYPE='NOCOMPRESS'
//     ("FITSIO status = 414: error uncompressing image / unknown
//     compression algorithm"), so the writer no longer emits it;
//   - 64-bit integers: fpack refuses BITPIX 64 ("Bad image datatype")
//     and funpack cannot read a compressed one ("unsupported datatype
//     for uncompressing image");
//   - HCOMPRESS on 1-D images or tiles thinner than 4 pixels (cfitsio's
//     imcomp_init_table limits);
//   - PLIO pixels outside 0..2^24-1 (cfitsio: "data out of range for
//     PLIO compression");
//   - inconsistent float options and caller scaling keywords.
//
// refusalReason is the text each TestWriteCompressedRefusals case must be
// refused with, so a case cannot pass on an unrelated error.
var refusalReason = map[string]string{
	"NOCOMPRESS int16":                   "NOCOMPRESS is not written",
	"NOCOMPRESS float32":                 "NOCOMPRESS is not written",
	"int64":                              "int64 pixels: cfitsio cannot write or read tile-compressed 64-bit integer images",
	"uint64":                             "uint64 pixels: cfitsio cannot write or read tile-compressed 64-bit integer images",
	"int64 GZIP_1":                       "int64 pixels: cfitsio cannot write or read tile-compressed 64-bit integer images",
	"uint64 GZIP_2":                      "uint64 pixels: cfitsio cannot write or read tile-compressed 64-bit integer images",
	"int64 PLIO in range":                "int64 pixels: cfitsio cannot write or read tile-compressed 64-bit integer images",
	"PLIO uint16 in range":               "PLIO_1 with uint16 pixels is not written",
	"float32 through the integer writer": "integer type required, got float32",
	"HCOMPRESS 1-D":                      "HCOMPRESS_1 cannot be used with 1-dimensional images",
	"HCOMPRESS row tiles":                "HCOMPRESS_1 minimum tile dimension is 4 pixels (ZTILE2 = 1)",
	"HCOMPRESS 3-D tile":                 "HCOMPRESS_1 tiles must be 2-D over the first two axes (ZTILE3 = 2)",
	"HCOMPRESS thin edge tile":           "HCOMPRESS_1 last tile along axis 1 would have 2 < 4 pixels",
	"HCOMPRESS negative scale":           "HCompressScale -1 < 0",
	"HCOMPRESS tile of 3":                "HCOMPRESS_1 minimum tile dimension is 4 pixels (ZTILE1 = 3)",
	"HCOMPRESS edge tile of 3":           "HCOMPRESS_1 last tile along axis 1 would have 3 < 4 pixels",
	"HCOMPRESS edge tile of 3 on axis 2": "HCOMPRESS_1 last tile along axis 2 would have 3 < 4 pixels",
	"unknown DitherMethod":               "unknown DitherMethod 3",
	"PLIO negative":                      "data out of range for PLIO compression (0 - 2**24): pixel 3 = -1",
	"PLIO above 2^24":                    "data out of range for PLIO compression (0 - 2**24): pixel 3 = 16777216",
	"float HCOMPRESS":                    "algorithm HCOMPRESS_1 is not supported for float input",
	"float64 PLIO":                       "algorithm PLIO_1 is not supported for float input",
	"float64 NOCOMPRESS":                 "NOCOMPRESS",
	"caller BLANK on float64":            "BLANK is not allowed in a floating-point image header",
	"NoDither with DitherMethod":         "NoDither and DitherMethod 2 are mutually exclusive",
	"ZDither0 > 10000":                   "ZDither0 10001 out of range 1..10000",
	"ZDither0 < 0":                       "ZDither0 -5 out of range 1..10000",
	"caller BZERO on uint16":             "BZERO in the header must be 32768 for this pixel type",
	"caller BSCALE on int16":             "BSCALE in the header is not supported for compressed writes",
	"float64 through the integer writer": "float64 pixels: use AppendCompressedFloat64Image",
	"caller BLANK on float32":            "BLANK is not allowed in a floating-point image header",
	"caller SIMPLE not logical":          "SIMPLE in the header is not a logical value",
}

func TestWriteCompressedRefusals(t *testing.T) {
	px16 := make([]int16, 64)
	cases := []struct {
		name string
		run  func(f *File) error
	}{
		{"NOCOMPRESS int16", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, px16, CompressOptions{Algorithm: compress.NoCompress})
			return err
		}},
		{"NOCOMPRESS float32", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{8, 8}, make([]float32, 64), CompressFloatOptions{Algorithm: compress.NoCompress})
			return err
		}},
		{"int64", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]int64, 64), CompressOptions{})
			return err
		}},
		{"uint64", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]uint64, 64), CompressOptions{})
			return err
		}},
		// The 64-bit refusal must come from the writer itself, not from
		// a codec that happens to reject BYTEPIX 8 (GZIP would encode it).
		{"int64 GZIP_1", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]int64, 64), CompressOptions{Algorithm: compress.GZIP1})
			return err
		}},
		{"uint64 GZIP_2", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]uint64, 64), CompressOptions{Algorithm: compress.GZIP2})
			return err
		}},
		{"int64 PLIO in range", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]int64, 64), CompressOptions{Algorithm: compress.PLIO1})
			return err
		}},
		// uint16 PLIO with every value >= 32768: the stored v-32768 is
		// inside PLIO's range, so only the writer's own refusal stops it.
		{"PLIO uint16 in range", func(f *File) error {
			v := make([]uint16, 64)
			for i := range v {
				v[i] = 32768 + uint16(i%3)
			}
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, v, CompressOptions{Algorithm: compress.PLIO1})
			return err
		}},
		{"float32 through the integer writer", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]float32, 64), CompressOptions{})
			return err
		}},
		{"HCOMPRESS 1-D", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{64}, px16, CompressOptions{Algorithm: compress.HCOMPRESS1})
			return err
		}},
		{"HCOMPRESS row tiles", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, px16, CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{8, 1}})
			return err
		}},
		{"HCOMPRESS 3-D tile", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{4, 4, 4}, px16, CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{4, 4, 2}})
			return err
		}},
		{"HCOMPRESS thin edge tile", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{10, 8}, make([]int16, 80), CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{8, 8}})
			return err
		}},
		// The limits are cfitsio's: imcomp_init_table refuses a tile
		// dimension < 4 (golden hcomp_i16_t3x8 records cfitsio's refusal
		// "Hcompress minimum tile dimension is 4 pixels") and treats an
		// edge tile of 1..3 pixels specially (golden hcomp_i16_11x8_t8x8,
		// where cfitsio enlarges the tile; we refuse). 4-pixel tiles and
		// 4-pixel edges are accepted (goldens hcomp_i16_12x24_t4x12,
		// hcomp_i16_12x9_t8x5).
		{"HCOMPRESS tile of 3", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{40, 25}, make([]int16, 1000), CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{3, 8}})
			return err
		}},
		{"HCOMPRESS edge tile of 3", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{11, 8}, make([]int16, 88), CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{8, 8}})
			return err
		}},
		{"HCOMPRESS edge tile of 3 on axis 2", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 11}, make([]int16, 88), CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{8, 8}})
			return err
		}},
		// cfitsio's fits_set_quantize_method refuses methods outside
		// NO_DITHER..SUBTRACTIVE_DITHER_2 ("illegal dithering value")
		{"unknown DitherMethod", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{8, 8}, make([]float32, 64), CompressFloatOptions{DitherMethod: 3})
			return err
		}},
		{"HCOMPRESS negative scale", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, px16, CompressOptions{Algorithm: compress.HCOMPRESS1, HCompressScale: -1})
			return err
		}},
		{"PLIO negative", func(f *File) error {
			v := make([]int16, 64)
			v[3] = -1
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, v, CompressOptions{Algorithm: compress.PLIO1})
			return err
		}},
		{"PLIO above 2^24", func(f *File) error {
			v := make([]int32, 64)
			v[3] = 1 << 24
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, v, CompressOptions{Algorithm: compress.PLIO1})
			return err
		}},
		{"float HCOMPRESS", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{8, 8}, make([]float32, 64), CompressFloatOptions{Algorithm: compress.HCOMPRESS1})
			return err
		}},
		{"float64 PLIO", func(f *File) error {
			_, err := AppendCompressedFloat64Image(f, nil, []int64{8, 8}, make([]float64, 64), CompressFloatOptions{Algorithm: compress.PLIO1})
			return err
		}},
		{"float64 NOCOMPRESS", func(f *File) error {
			_, err := AppendCompressedFloat64Image(f, nil, []int64{8, 8}, make([]float64, 64), CompressFloatOptions{Algorithm: compress.NoCompress})
			return err
		}},
		{"caller BLANK on float64", func(f *File) error {
			h := header.New()
			h.Set("BLANK", int64(5), "")
			_, err := AppendCompressedFloat64Image(f, h, []int64{8, 8}, make([]float64, 64), CompressFloatOptions{})
			return err
		}},
		{"NoDither with DitherMethod", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{8, 8}, make([]float32, 64), CompressFloatOptions{NoDither: true, DitherMethod: compress.SubtractiveDither2})
			return err
		}},
		{"ZDither0 > 10000", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{8, 8}, make([]float32, 64), CompressFloatOptions{ZDither0: 10001})
			return err
		}},
		{"ZDither0 < 0", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{8, 8}, make([]float32, 64), CompressFloatOptions{ZDither0: -5})
			return err
		}},
		{"caller BZERO on uint16", func(f *File) error {
			h := header.New()
			h.Set("BZERO", 1000.0, "")
			_, err := AppendCompressedImage(f, h, []int64{8, 8}, make([]uint16, 64), CompressOptions{})
			return err
		}},
		{"caller BSCALE on int16", func(f *File) error {
			h := header.New()
			h.Set("BSCALE", 2.0, "")
			_, err := AppendCompressedImage(f, h, []int64{8, 8}, px16, CompressOptions{})
			return err
		}},
		// cfitsio can write ZBITPIX -64 (fits_quantize_double); this
		// writer has no float64 path and must say so, not misreport it.
		{"float64 through the integer writer", func(f *File) error {
			_, err := AppendCompressedImage(f, nil, []int64{8, 8}, make([]float64, 64), CompressOptions{})
			return err
		}},
		// BLANK with BITPIX < 0 is illegal (FITS 4.0 sect. 4.4.2.5); copied
		// into a ZBITPIX -32 header without ZBLANK, every reader takes it
		// as the null of the quantized integers and nulls real pixels
		// (golden rice_f32_blank_hdr, TestFloatBLANKNullsRealPixels).
		{"caller BLANK on float32", func(f *File) error {
			h := header.New()
			h.Set("BLANK", int64(5), "")
			v := make([]float32, 64)
			for i := range v {
				v[i] = float32(100 + i%7)
			}
			_, err := AppendCompressedFloat32Image(f, h, []int64{8, 8}, v, CompressFloatOptions{})
			return err
		}},
		// SIMPLE must be the logical T or F (FITS 4.0 sect. 4.4.1.1); a
		// caller value the writer would have to replace is refused.
		{"caller SIMPLE not logical", func(f *File) error {
			h := header.New()
			h.Set("SIMPLE", "X", "")
			_, err := AppendCompressedImage(f, h, []int64{8, 8}, px16, CompressOptions{})
			return err
		}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "x.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
				t.Fatal(err)
			}
			err = c.run(f)
			if err == nil {
				t.Fatalf("%s: written without error; cfitsio cannot write/read it", c.name)
			}
			want, ok := refusalReason[c.name]
			if !ok {
				t.Fatalf("no expected refusal reason for %q", c.name)
			}
			if !strings.Contains(err.Error(), want) {
				t.Fatalf("refused for the wrong reason:\n  got  %v\n  want an error containing %q", err, want)
			}
			t.Logf("refused: %v", err)
		})
	}
	// The matching unsigned convention in the caller header is accepted.
	f, err := Create(filepath.Join(t.TempDir(), "ok.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	WriteImage(f, nil, []int64{}, []uint8{})
	h := header.New()
	h.Set("BZERO", 32768.0, "")
	h.Set("BSCALE", 1.0, "")
	if _, err := AppendCompressedImage(f, h, []int64{8, 8}, make([]uint16, 64), CompressOptions{}); err != nil {
		t.Fatalf("uint16 with BZERO=32768/BSCALE=1 refused: %v", err)
	}
}

// TestWriteCompressedPLIO: mask data through PLIO_1 encode → decode.
func TestWriteCompressedPLIO(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "write_plio.fits")

	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{})
	// Sparse mask.
	shape := []int64{50, 50}
	pixels := make([]int16, 50*50)
	for row := 10; row < 20; row++ {
		for col := 15; col < 35; col++ {
			pixels[row*50+col] = 1
		}
	}
	for row := 30; row < 40; row++ {
		for col := 5; col < 45; col++ {
			pixels[row*50+col] = 2
		}
	}
	if _, err := AppendCompressedImage(f, nil, shape, pixels, CompressOptions{
		Algorithm: compress.PLIO1,
	}); err != nil {
		t.Fatal(err)
	}
	f.Close()

	fr, _ := Open(path)
	defer fr.Close()
	h, _ := fr.HDU(1)
	cimg := h.(*CompressedImageHDU)
	if cimg.CompressionType() != "PLIO_1" {
		t.Fatalf("algorithm %q", cimg.CompressionType())
	}
	got, err := ReadPixelsCompressed[int16](cimg)
	if err != nil {
		t.Fatal(err)
	}
	for i := range pixels {
		if got[i] != pixels[i] {
			t.Fatalf("pixel %d: got %d want %d", i, got[i], pixels[i])
		}
	}
}

// TestWriteCompressedHCOMPRESS: 2D image through the forward H-transform
// write path and the inverse read path.
func TestWriteCompressedHCOMPRESS(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "write_hcompress.fits")

	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{})
	shape := []int64{32, 32}
	pixels := make([]int16, 32*32)
	for i := range pixels {
		pixels[i] = int16((i * 3) % 500)
	}
	// HCOMPRESS needs a 2D tile; use the whole image as one tile.
	if _, err := AppendCompressedImage(f, nil, shape, pixels, CompressOptions{
		Algorithm: compress.HCOMPRESS1,
		TileShape: []int64{32, 32},
	}); err != nil {
		t.Fatal(err)
	}
	f.Close()

	fr, _ := Open(path)
	defer fr.Close()
	h, _ := fr.HDU(1)
	cimg := h.(*CompressedImageHDU)
	if cimg.CompressionType() != "HCOMPRESS_1" {
		t.Fatalf("algorithm %q", cimg.CompressionType())
	}
	got, err := ReadPixelsCompressed[int16](cimg)
	if err != nil {
		t.Fatalf("decompress: %v", err)
	}
	for i := range pixels {
		if got[i] != pixels[i] {
			t.Fatalf("pixel %d: got %d want %d", i, got[i], pixels[i])
			break
		}
	}
}

// writeReadCompressedInt16Test is a shared helper for the algorithm
// variants that don't need special per-algorithm setup.
func writeReadCompressedInt16Test(t *testing.T, path string, algo compress.Algorithm, shape []int64) {
	t.Helper()
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	n := int(shape[0] * shape[1])
	pixels := make([]int16, n)
	for i := range pixels {
		pixels[i] = int16(i)
	}
	if _, err := AppendCompressedImage(f, nil, shape, pixels, CompressOptions{
		Algorithm: algo,
	}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()
	h, err := fr.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	cimg, ok := h.(*CompressedImageHDU)
	if !ok {
		t.Fatalf("HDU(1) is %T, want *CompressedImageHDU", h)
	}
	if cimg.CompressionType() != algo.String() {
		t.Errorf("CompressionType: got %q want %q", cimg.CompressionType(), algo.String())
	}
	got, err := ReadPixelsCompressed[int16](cimg)
	if err != nil {
		t.Fatalf("decompress: %v", err)
	}
	for i := range pixels {
		if got[i] != pixels[i] {
			t.Fatalf("pixel %d: got %d want %d", i, got[i], pixels[i])
		}
	}
}

// TestCompressedWriterShapeRefusals: both compressed writers refuse a
// shape they cannot write before anything is written: data shorter or
// longer than the shape product (ErrShapeMismatch, as WriteImage), an
// empty shape, a zero or negative axis, more than 6 axes (cfitsio's
// MAX_COMPRESS_DIM; fits_img_compress refuses 7 axes with status 212,
// golden img/rice_i16_7d, and funpack/fitsverify reject such a file), a
// tile shape of the wrong rank and a zero or negative tile axis. Each of
// these guards could be removed before without a failing test (two of
// them then panicked, the others wrote a file).
func TestCompressedWriterShapeRefusals(t *testing.T) {
	type tc struct {
		name          string
		shape         []int64
		n             int
		tile          []int64
		want          string
		shapeMismatch bool
	}
	cases := []tc{
		{"short data", []int64{4, 4}, 8, nil, "data length 8 != product of shape 16", true},
		{"long data", []int64{4, 4}, 20, nil, "data length 20 != product of shape 16", true},
		{"empty shape", []int64{}, 0, nil, "empty shape", false},
		{"zero axis", []int64{4, 0}, 0, nil, "non-positive axis 0", false},
		{"negative axis", []int64{4, -2}, 8, nil, "non-positive axis -2", false},
		{"7 axes", []int64{2, 2, 2, 2, 2, 2, 2}, 128, nil, "7 axes", false},
		{"tile rank", []int64{4, 4}, 16, []int64{4}, "tile shape has 1 axes, image has 2", false},
		{"zero tile axis", []int64{4, 4}, 16, []int64{4, 0}, "non-positive tile axis 0", false},
		{"negative tile axis", []int64{4, 4}, 16, []int64{-1, 4}, "non-positive tile axis -1", false},
	}
	writers := []struct {
		name string
		run  func(f *File, c tc) error
	}{
		{"AppendCompressedImage", func(f *File, c tc) error {
			_, err := AppendCompressedImage(f, nil, c.shape, make([]int16, c.n), CompressOptions{TileShape: c.tile})
			return err
		}},
		{"AppendCompressedFloat32Image", func(f *File, c tc) error {
			_, err := AppendCompressedFloat32Image(f, nil, c.shape, make([]float32, c.n), CompressFloatOptions{TileShape: c.tile})
			return err
		}},
	}
	for _, w := range writers {
		for _, c := range cases {
			t.Run(w.name+"/"+c.name, func(t *testing.T) {
				path := filepath.Join(t.TempDir(), "x.fits")
				f, err := Create(path)
				if err != nil {
					t.Fatal(err)
				}
				if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
					t.Fatal(err)
				}
				func() {
					defer func() {
						if r := recover(); r != nil {
							t.Fatalf("panicked: %v", r)
						}
					}()
					err = w.run(f, c)
				}()
				if err == nil || !strings.Contains(err.Error(), c.want) {
					t.Fatalf("err = %v, want one containing %q", err, c.want)
				}
				if errors.Is(err, ErrShapeMismatch) != c.shapeMismatch {
					t.Errorf("errors.Is(err, ErrShapeMismatch) = %v, want %v", !c.shapeMismatch, c.shapeMismatch)
				}
				if f.NumHDU() != 1 {
					t.Errorf("%d HDUs after the refusal, want only the primary", f.NumHDU())
				}
				f.Close()
			})
		}
	}
}

// TestCompressedFloatQLevelNaNInfRefused: a QLevel of NaN or -Inf has no
// quantization step (NaN gives delta NaN; -Inf gives delta +Inf and a
// zero point 0*Inf = NaN). cfitsio does not check either: its range
// tests are false for NaN, so it writes garbage integers and NaN/Inf
// ZSCALE/ZZERO, and every pixel reads back NaN. The writer refuses both
// instead of losing the data. QLevel +Inf (delta 0) and a tiny negative
// QLevel stay accepted: the tile's range (as here) or its zero point
// (constant tiles; cfitsio overflows there and writes wrong pixels,
// TestQuantizerOverflowStoredLosslessly) does not fit the integers, and
// the tile is stored losslessly.
func TestCompressedFloatQLevelNaNInfRefused(t *testing.T) {
	data := make([]float32, 100)
	for i := range data {
		data[i] = float32(100 + i%10)
	}
	for _, q := range []float32{float32(math.NaN()), float32(math.Inf(-1))} {
		f, err := Create(filepath.Join(t.TempDir(), "q.fits"))
		if err != nil {
			t.Fatal(err)
		}
		if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
			t.Fatal(err)
		}
		_, werr := AppendCompressedFloat32Image(f, nil, []int64{10, 10}, data, CompressFloatOptions{QLevel: q})
		f.Close()
		if werr == nil {
			t.Errorf("QLevel %v: write accepted (every pixel would read back NaN)", q)
		} else {
			t.Logf("QLevel %v refused: %v", q, werr)
		}
	}
	for _, q := range []float32{float32(math.Inf(1)), -1e-40} {
		path := filepath.Join(t.TempDir(), "q.fits")
		f, err := Create(path)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
			t.Fatal(err)
		}
		if _, err := AppendCompressedFloat32Image(f, nil, []int64{10, 10}, data, CompressFloatOptions{QLevel: q}); err != nil {
			t.Fatalf("QLevel %v: %v", q, err)
		}
		if err := f.Close(); err != nil {
			t.Fatal(err)
		}
		ff, h := openCompressed(t, path)
		got, err := ReadPixelsCompressed[float32](h)
		ff.Close()
		if err != nil {
			t.Fatal(err)
		}
		for i := range data {
			if got[i] != data[i] {
				t.Fatalf("QLevel %v: pixel %d = %v, want %v (lossless fallback)", q, i, got[i], data[i])
			}
		}
	}
}
