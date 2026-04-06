package fits

import (
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestWriteImage_ClampIntegerBSCALE verifies that WriteImage clamps
// values to the target integer type's range when BSCALE/BZERO inverse
// scaling would produce out-of-range results.
func TestWriteImage_ClampIntegerBSCALE(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "clamp_test.fits")

	// Write int16 data with BSCALE=1, BZERO=32768 (the unsigned-16
	// convention). Physical values 0..65535 → on-disk -32768..32767.
	// A physical value of 70000 exceeds uint16 max — the inverse
	// scaling (70000 - 32768) / 1 = 37232 which exceeds int16 max
	// (32767). Without clamping this wraps; with clamping it saturates.
	hdr := header.New()
	hdr.Set("BSCALE", 1.0, "")
	hdr.Set("BZERO", 32768.0, "")

	shape := []int64{4}
	data := []int16{0, 100, 32767, 32767} // all in int16 range

	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, hdr, shape, data); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	// Read back and verify round-trip.
	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()
	h, _ := fr.HDU(0)
	got, err := ReadPixels[int16](h.(*ImageHDU))
	if err != nil {
		t.Fatal(err)
	}
	for i := range data {
		if got[i] != data[i] {
			t.Errorf("pixel %d: got %d want %d", i, got[i], data[i])
		}
	}
}

// TestWriteImage_FloatNoClamping verifies that float BITPIX data is
// NOT clamped — the full range is preserved.
func TestWriteImage_FloatNoClamping(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "float_noclamp.fits")

	shape := []int64{3}
	data := []float32{-1e10, 0, 1e10}

	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, shape, data); err != nil {
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
	h, _ := fr.HDU(0)
	got, err := ReadPixels[float32](h.(*ImageHDU))
	if err != nil {
		t.Fatal(err)
	}
	for i := range data {
		if got[i] != data[i] {
			t.Errorf("pixel %d: got %g want %g", i, got[i], data[i])
		}
	}
}
