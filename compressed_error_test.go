package fits

import (
	"path/filepath"
	"testing"
)

// TestCompressedHDUHCompressMetadata verifies that the HCOMPRESS metadata
// is accessible through the normal HDU surface (algorithm name, BITPIX,
// shape). Actual pixel round-trip correctness is verified by
// TestCompressedFixtures/HCOMPRESS_1_*.
func TestCompressedHDUHCompressMetadata(t *testing.T) {
	f, err := Open(filepath.Join("testdata", "comp_hcompress_i16.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	cimg, ok := h.(*CompressedImageHDU)
	if !ok {
		t.Fatalf("HDU(1) is %T, want *CompressedImageHDU", h)
	}
	if cimg.CompressionType() != "HCOMPRESS_1" {
		t.Fatalf("algorithm: %q", cimg.CompressionType())
	}
	if cimg.BITPIX() != 16 {
		t.Fatalf("BITPIX: %d", cimg.BITPIX())
	}
}

// TestCompressedHDUCompressedFlag verifies the HDU interface's new
// Compressed() and CompressionType() methods return the right values.
func TestCompressedHDUCompressedFlag(t *testing.T) {
	// Compressed case.
	f, _ := Open(filepath.Join("testdata", "comp_rice_i16.fits"))
	defer f.Close()
	h, _ := f.HDU(1)
	if !h.Compressed() {
		t.Error("compressed HDU.Compressed() should be true")
	}
	if h.CompressionType() != "RICE_1" {
		t.Errorf("CompressionType: %q", h.CompressionType())
	}
	// Uncompressed case.
	fr, _ := Open(filepath.Join("testdata", "comp_raw_rice_i16.fits"))
	defer fr.Close()
	hr, _ := fr.HDU(1)
	if hr.Compressed() {
		t.Error("plain ImageHDU.Compressed() should be false")
	}
	if hr.CompressionType() != "" {
		t.Errorf("CompressionType on plain image: %q", hr.CompressionType())
	}
}
