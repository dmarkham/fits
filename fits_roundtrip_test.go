package fits

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestRoundTripFloat32 writes a 4x3 float32 image to a temp file, reopens
// it, reads back the pixels, and checks they match exactly.
func TestRoundTripFloat32(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "float32.fits")

	// Write
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	h := header.New()
	h.Set("OBJECT", "unit test", "synthetic image")
	shape := []int64{4, 3} // NAXIS1=4, NAXIS2=3 → 12 pixels
	data := []float32{
		1, 2, 3, 4,
		5, 6, 7, 8,
		9, 10, 11, 12,
	}
	if _, err := WriteImage(f, h, shape, data); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	// Read
	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	if f2.NumHDU() != 1 {
		t.Fatalf("NumHDU=%d", f2.NumHDU())
	}
	img, err := f2.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if img.BITPIX() != -32 {
		t.Fatalf("BITPIX %d", img.BITPIX())
	}
	if img.NAXIS() != 2 {
		t.Fatalf("NAXIS %d", img.NAXIS())
	}
	sh := img.Shape()
	if len(sh) != 2 || sh[0] != 4 || sh[1] != 3 {
		t.Fatalf("shape %v", sh)
	}
	obj, err := img.Header().String("OBJECT")
	if err != nil || obj != "unit test" {
		t.Fatalf("OBJECT %q %v", obj, err)
	}

	got, err := ReadPixels[float32](img)
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != len(data) {
		t.Fatalf("len mismatch %d vs %d", len(got), len(data))
	}
	for i := range got {
		if got[i] != data[i] {
			t.Fatalf("pixel %d: got %v want %v", i, got[i], data[i])
		}
	}
}

func TestRoundTripInt16(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "int16.fits")

	f, _ := Create(path)
	shape := []int64{8, 2}
	data := []int16{-100, -50, 0, 50, 100, 200, 300, 400,
		1, 2, 3, 4, 5, 6, 7, 8}
	if _, err := WriteImage(f, nil, shape, data); err != nil {
		t.Fatal(err)
	}
	f.Close()

	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	img, _ := f2.Primary()
	if img.BITPIX() != 16 {
		t.Fatalf("BITPIX %d", img.BITPIX())
	}
	got, err := ReadPixels[int16](img)
	if err != nil {
		t.Fatal(err)
	}
	for i := range got {
		if got[i] != data[i] {
			t.Fatalf("pixel %d: %v vs %v", i, got[i], data[i])
		}
	}
}

func TestReadNonexistent(t *testing.T) {
	if _, err := Open(filepath.Join(t.TempDir(), "no.fits")); err == nil {
		t.Fatal("expected open error")
	}
}

func TestEmptyFile(t *testing.T) {
	path := filepath.Join(t.TempDir(), "empty.fits")
	if err := os.WriteFile(path, nil, 0644); err != nil {
		t.Fatal(err)
	}
	_, err := Open(path)
	if err != ErrEmptyFile {
		t.Fatalf("want ErrEmptyFile, got %v", err)
	}
}

func TestNotFITS(t *testing.T) {
	path := filepath.Join(t.TempDir(), "bad.fits")
	// Exactly one block of garbage.
	buf := make([]byte, 2880)
	for i := range buf {
		buf[i] = 'X'
	}
	if err := os.WriteFile(path, buf, 0644); err != nil {
		t.Fatal(err)
	}
	_, err := Open(path)
	if err == nil {
		t.Fatal("expected error")
	}
}
