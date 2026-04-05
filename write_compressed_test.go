package fits

import (
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/compress"
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

// TestWriteCompressedNOCOMPRESS: the identity path, good for sanity.
func TestWriteCompressedNOCOMPRESS(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "write_nocompress.fits")
	writeReadCompressedInt16Test(t, path, compress.NoCompress, []int64{16, 16})
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
