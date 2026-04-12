package fits

import (
	"math"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestReadFloat32_Int16_Signed verifies normalization for signed 16-bit
// data (BSCALE=1 BZERO=0, the default). Raw [-32768, 32767] is mapped
// to [0, 1].
func TestReadFloat32_Int16_Signed(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "int16_signed.fits")

	hdr := header.New()
	data := []int16{-32768, 0, 32767}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	WriteImage(f, hdr, []int64{3}, data)
	f.Close()

	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	prim, _ := f2.Primary()
	got, err := ReadFloat32(prim)
	if err != nil {
		t.Fatal(err)
	}

	// physMin=-32768, physMax=32767, span=65535.
	// [-32768 → 0.0, 0 → 32768/65535≈0.50001, 32767 → 1.0]
	const tol = 1e-5
	if got[0] != 0.0 {
		t.Errorf("pixel 0 (min): got %v want 0.0", got[0])
	}
	if got[2] != 1.0 {
		t.Errorf("pixel 2 (max): got %v want 1.0", got[2])
	}
	if math.Abs(float64(got[1]-0.5)) > tol {
		t.Errorf("pixel 1 (mid): got %v want ~0.5", got[1])
	}
}

// TestReadFloat32_Uint8 verifies 8-bit normalization: [0, 255] → [0, 1].
func TestReadFloat32_Uint8(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "uint8.fits")

	hdr := header.New()
	data := []uint8{0, 128, 255}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	WriteImage(f, hdr, []int64{3}, data)
	f.Close()

	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	prim, _ := f2.Primary()
	got, err := ReadFloat32(prim)
	if err != nil {
		t.Fatal(err)
	}

	const tol = 1e-5
	if got[0] != 0.0 {
		t.Errorf("pixel 0: got %v want 0.0", got[0])
	}
	if math.Abs(float64(got[2]-1.0)) > tol {
		t.Errorf("pixel 2: got %v want 1.0", got[2])
	}
	if math.Abs(float64(got[1]-128.0/255.0)) > tol {
		t.Errorf("pixel 1: got %v want %v", got[1], 128.0/255.0)
	}
}

// TestReadFloat32_Float32_Passthrough verifies float BITPIX returns
// values as-is with no normalization applied (assumed already [0,1]).
func TestReadFloat32_Float32_Passthrough(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "float32.fits")

	hdr := header.New()
	data := []float32{0.0, 0.25, 0.5, 0.75, 1.0}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	WriteImage(f, hdr, []int64{5}, data)
	f.Close()

	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	prim, _ := f2.Primary()
	got, err := ReadFloat32(prim)
	if err != nil {
		t.Fatal(err)
	}
	for i, w := range data {
		if got[i] != w {
			t.Errorf("pixel %d: got %v want %v", i, got[i], w)
		}
	}
}

// TestReadFloat32_CompressedImage verifies ReadFloat32 dispatches
// correctly for *CompressedImageHDU via the HDU interface, and that
// the normalized output matches the normalized output from the raw
// uncompressed counterpart.
func TestReadFloat32_CompressedImage(t *testing.T) {
	compPath := filepath.Join("testdata", "comp_rice_i16.fits")
	rawPath := filepath.Join("testdata", "comp_raw_rice_i16.fits")

	// The actual image is in HDU 1 for both fixtures.
	rawF, err := Open(rawPath)
	if err != nil {
		t.Skipf("raw fixture not found: %v", err)
	}
	defer rawF.Close()
	rawHDU, err := rawF.HDU(1)
	if err != nil {
		t.Fatalf("raw HDU(1): %v", err)
	}
	rawPix, err := ReadFloat32(rawHDU)
	if err != nil {
		t.Fatalf("raw ReadFloat32: %v", err)
	}
	if len(rawPix) == 0 {
		t.Fatal("raw image has 0 pixels")
	}

	compF, err := Open(compPath)
	if err != nil {
		t.Skipf("compressed fixture not found: %v", err)
	}
	defer compF.Close()
	compHDU, err := compF.HDU(1)
	if err != nil {
		t.Fatalf("comp HDU(1): %v", err)
	}
	compPix, err := ReadFloat32(compHDU)
	if err != nil {
		t.Fatalf("compressed ReadFloat32: %v", err)
	}

	if len(compPix) != len(rawPix) {
		t.Fatalf("pixel count mismatch: compressed %d vs raw %d", len(compPix), len(rawPix))
	}

	diffs := 0
	for i := range rawPix {
		if rawPix[i] != compPix[i] {
			diffs++
			if diffs <= 3 {
				t.Errorf("pixel %d: raw=%v compressed=%v", i, rawPix[i], compPix[i])
			}
		}
	}
	if diffs > 3 {
		t.Errorf("... and %d more differing pixels", diffs-3)
	}
}

// TestReadFloat32_RangeInvariant verifies every pixel from an integer-
// BITPIX image is in [0, 1] after normalization. Writes a gradient
// spanning the full int16 range.
func TestReadFloat32_RangeInvariant(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "range.fits")

	hdr := header.New()
	data := make([]int16, 1000)
	for i := range data {
		data[i] = int16(-32768 + int(i)*65535/999)
	}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	WriteImage(f, hdr, []int64{1000}, data)
	f.Close()

	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	prim, _ := f2.Primary()
	got, err := ReadFloat32(prim)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range got {
		if v < 0 || v > 1 {
			t.Errorf("pixel %d out of [0,1]: %v", i, v)
		}
	}
}

// TestReadFloat32_UnsupportedHDU verifies that passing a non-image HDU
// returns a clear error.
func TestReadFloat32_UnsupportedHDU(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "table.fits")

	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	hdr := header.New()
	WriteImage(f, hdr, []int64{1}, []float32{1.0})
	cols := []ColumnData{{Name: "X", DataInt32: []int32{1, 2, 3}}}
	AppendBinaryTable(f, header.New(), cols)
	f.Close()

	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	tbl, _ := f2.HDU(1)
	_, err = ReadFloat32(tbl)
	if err == nil {
		t.Fatal("expected error for non-image HDU")
	}
}
