package fits

import (
	"path/filepath"
	"testing"
)

func TestBinaryTableRoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "tbl.fits")

	// Primary placeholder + one BINTABLE.
	f, _ := Create(path)
	// Primary: 0-axis image (required by FITS — every file starts with a primary).
	_, err := WriteImage(f, nil, []int64{}, []float32{})
	if err != nil {
		t.Fatal(err)
	}

	cols := []ColumnData{
		{Name: "ID", DataInt32: []int32{1, 2, 3, 4}},
		{Name: "FLUX", Unit: "erg/s", DataFloat64: []float64{1.5, 2.5, 3.5, 4.5}},
		{Name: "LABEL", DataString: []string{"A", "B", "C", "D"}},
	}
	if _, err := AppendBinaryTable(f, nil, cols); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	// Read back.
	f2, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f2.Close()
	if f2.NumHDU() != 2 {
		t.Fatalf("NumHDU=%d", f2.NumHDU())
	}
	h, err := f2.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	tbl, ok := h.(*BinaryTableHDU)
	if !ok {
		t.Fatalf("HDU 1 is %T", h)
	}
	if tbl.NumRows() != 4 {
		t.Fatalf("rows=%d", tbl.NumRows())
	}
	cs, err := tbl.Columns()
	if err != nil {
		t.Fatal(err)
	}
	if len(cs) != 3 {
		t.Fatalf("cols=%d", len(cs))
	}
	if cs[0].Name != "ID" || cs[1].Name != "FLUX" || cs[2].Name != "LABEL" {
		t.Fatalf("names: %v %v %v", cs[0].Name, cs[1].Name, cs[2].Name)
	}

	ids, err := ReadColumn[int32](tbl, 1)
	if err != nil {
		t.Fatal(err)
	}
	want := []int32{1, 2, 3, 4}
	for i, v := range ids {
		if v != want[i] {
			t.Fatalf("id %d: %d vs %d", i, v, want[i])
		}
	}

	flux, err := ReadColumn[float64](tbl, 2)
	if err != nil {
		t.Fatal(err)
	}
	wantF := []float64{1.5, 2.5, 3.5, 4.5}
	for i, v := range flux {
		if v != wantF[i] {
			t.Fatalf("flux %d: %v vs %v", i, v, wantF[i])
		}
	}

	labels, err := ReadColumn[string](tbl, 3)
	if err != nil {
		t.Fatal(err)
	}
	wantS := []string{"A", "B", "C", "D"}
	for i, v := range labels {
		if v != wantS[i] {
			t.Fatalf("label %d: %q vs %q", i, v, wantS[i])
		}
	}
}
