package fits

import (
	"path/filepath"
	"testing"
)

func TestReadSubset2D(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "sub.fits")

	// 5×4 int32 image: pixel[y][x] = y*10 + x
	// Row-major on disk: axis 0 = x (NAXIS1), axis 1 = y (NAXIS2)
	shape := []int64{5, 4}
	data := make([]int32, 20)
	for y := range int64(4) {
		for x := range int64(5) {
			data[y*5+x] = int32(y*10 + x)
		}
	}
	f, _ := Create(path)
	WriteImage(f, nil, shape, data)
	f.Close()

	f2, _ := Open(path)
	defer f2.Close()
	img, _ := f2.Primary()

	// Read a 3×2 subset: x=[1,4), y=[1,3), stride [1,1]
	got, err := ReadSubset[int32](img, []int64{1, 1}, []int64{4, 3}, []int64{1, 1})
	if err != nil {
		t.Fatal(err)
	}
	// Expected: y=1 [11 12 13], y=2 [21 22 23]
	want := []int32{11, 12, 13, 21, 22, 23}
	if len(got) != len(want) {
		t.Fatalf("len %d want %d", len(got), len(want))
	}
	for i := range got {
		if got[i] != want[i] {
			t.Fatalf("idx %d: got %d want %d (full %v)", i, got[i], want[i], got)
		}
	}
}

func TestReadSubsetStride(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "stride.fits")
	shape := []int64{8, 1}
	data := []int32{10, 11, 12, 13, 14, 15, 16, 17}
	f, _ := Create(path)
	WriteImage(f, nil, shape, data)
	f.Close()

	f2, _ := Open(path)
	defer f2.Close()
	img, _ := f2.Primary()
	got, err := ReadSubset[int32](img, []int64{0, 0}, []int64{8, 1}, []int64{2, 1})
	if err != nil {
		t.Fatal(err)
	}
	want := []int32{10, 12, 14, 16}
	if len(got) != len(want) {
		t.Fatalf("len %d", len(got))
	}
	for i := range got {
		if got[i] != want[i] {
			t.Fatalf("idx %d: %d vs %d", i, got[i], want[i])
		}
	}
}
