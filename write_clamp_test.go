package fits

import (
	"errors"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestWriteImage_IntegerBZERORoundTrip: int16 data written with a caller
// BZERO = 32768 (raw = v - 32768) round-trips when read into a type that
// holds the physical range [0, 65535] (int16 itself is refused: it cannot
// hold 65535). A value whose raw form leaves the stored range is refused
// with ErrNumOverflow, never clamped (TestWriteImageOverflowMatchesCfitsio;
// here: -1, whose raw form is -32769).
func TestWriteImage_IntegerBZERORoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "clamp_test.fits")

	hdr := header.New()
	hdr.Set("BSCALE", 1.0, "")
	hdr.Set("BZERO", 32768.0, "")

	shape := []int64{4}
	data := []int16{0, 100, 32767, 32767}

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
	if _, err := ReadPixels[int16](h.(*ImageHDU)); !errors.Is(err, ErrTypeMismatch) {
		t.Errorf("ReadPixels[int16] on physical range [0, 65535]: err %v, want ErrTypeMismatch", err)
	}
	got, err := ReadPixels[int32](h.(*ImageHDU))
	if err != nil {
		t.Fatal(err)
	}
	for i := range data {
		if got[i] != int32(data[i]) {
			t.Errorf("pixel %d: got %d want %d", i, got[i], data[i])
		}
	}

	f2, err := Create(filepath.Join(dir, "overflow.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	if _, err := WriteImage(f2, hdr, shape, []int16{0, 100, -1, 32767}); !errors.Is(err, ErrNumOverflow) || f2.NumHDU() != 0 {
		t.Errorf("WriteImage of -1 with BZERO 32768: err %v, %d HDUs; want ErrNumOverflow and none", err, f2.NumHDU())
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
