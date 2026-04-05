package fits

import (
	"path/filepath"
	"testing"
)

func TestEditHeaderInPlace(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "edit.fits")

	// Create a file.
	f, _ := Create(path)
	shape := []int64{4, 4}
	data := make([]float32, 16)
	for i := range data {
		data[i] = float32(i)
	}
	WriteImage(f, nil, shape, data)
	f.Close()

	// Open for edit and mutate a header keyword.
	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	p, _ := fe.Primary()
	h := p.Header()
	if err := h.Set("OBJECT", "test object", "added by test"); err != nil {
		t.Fatal(err)
	}
	if err := fe.Flush(); err != nil {
		t.Fatal(err)
	}
	if err := fe.Close(); err != nil {
		t.Fatal(err)
	}

	// Reopen and verify the mutation persisted.
	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()
	p2, _ := fr.Primary()
	obj, err := p2.Header().String("OBJECT")
	if err != nil {
		t.Fatal(err)
	}
	if obj != "test object" {
		t.Fatalf("object = %q", obj)
	}
	// Pixels unchanged.
	pix, _ := ReadPixels[float32](p2)
	for i := range data {
		if pix[i] != data[i] {
			t.Fatalf("pixel %d changed after edit", i)
		}
	}
}

func TestOverwritePixels(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ow.fits")
	f, _ := Create(path)
	shape := []int64{8, 2}
	initial := make([]int16, 16)
	WriteImage(f, nil, shape, initial)
	f.Close()

	fe, _ := OpenForEdit(path)
	p, _ := fe.Primary()
	newData := []int16{100, 101, 102, 103, 104, 105, 106, 107,
		200, 201, 202, 203, 204, 205, 206, 207}
	if err := OverwritePixels(p, newData); err != nil {
		t.Fatal(err)
	}
	fe.Close()

	fr, _ := Open(path)
	defer fr.Close()
	p2, _ := fr.Primary()
	got, _ := ReadPixels[int16](p2)
	for i := range newData {
		if got[i] != newData[i] {
			t.Fatalf("pixel %d: %d vs %d", i, got[i], newData[i])
		}
	}
}

func TestOverwritePixelsShapeMismatch(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "mis.fits")
	f, _ := Create(path)
	WriteImage(f, nil, []int64{4, 4}, make([]int32, 16))
	f.Close()

	fe, _ := OpenForEdit(path)
	p, _ := fe.Primary()
	err := OverwritePixels(p, []int32{1, 2, 3}) // wrong length
	if err == nil {
		t.Fatal("expected shape mismatch")
	}
	fe.Close()
}

func TestAppendImageInEdit(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "append.fits")

	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{}) // placeholder primary
	f.Close()

	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	data := []float32{1, 2, 3, 4}
	if _, err := AppendImage(fe, nil, []int64{2, 2}, data); err != nil {
		t.Fatal(err)
	}
	fe.Close()

	fr, _ := Open(path)
	defer fr.Close()
	if fr.NumHDU() != 2 {
		t.Fatalf("NumHDU=%d", fr.NumHDU())
	}
	h, _ := fr.HDU(1)
	img := h.(*ImageHDU)
	if img.NAXIS() != 2 {
		t.Fatalf("appended NAXIS %d", img.NAXIS())
	}
	got, _ := ReadPixels[float32](img)
	for i := range data {
		if got[i] != data[i] {
			t.Fatalf("pixel %d: %v vs %v", i, got[i], data[i])
		}
	}
}
