package fits

import (
	"path/filepath"
	"testing"
)

// TestVarColumnRoundTripFloat32 writes a binary table with a single
// variable-length float32 column whose row lengths differ, reopens the
// file, reads the column back through ReadVarColumn, and asserts every
// row's contents match byte-for-byte.
func TestVarColumnRoundTripFloat32(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "vla_f32.fits")

	rows := [][]float32{
		{1.0, 2.0, 3.0},
		{},
		{10.5, 20.5},
		{100.25, 200.25, 300.25, 400.25, 500.25},
		{-1, -2, -3, -4},
	}

	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	cols := []ColumnData{
		{Name: "SPEC", Unit: "counts", DataVarFloat32: rows},
	}
	if _, err := AppendBinaryTable(f, nil, cols); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	fr, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer fr.Close()

	h, err := fr.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	tbl, ok := h.(*BinaryTableHDU)
	if !ok {
		t.Fatalf("HDU 1 is %T", h)
	}
	// Header must declare PCOUNT > 0 (heap is non-empty).
	pcount, err := tbl.Header().Int("PCOUNT")
	if err != nil {
		t.Fatal(err)
	}
	if pcount == 0 {
		t.Fatalf("PCOUNT is zero for non-empty VLA heap")
	}

	// Columns metadata.
	cs, err := tbl.Columns()
	if err != nil {
		t.Fatal(err)
	}
	if len(cs) != 1 || cs[0].Type != ColVarArray {
		t.Fatalf("unexpected column metadata: %+v", cs)
	}

	// Read back.
	vc, err := ReadVarColumn[float32](tbl, 1)
	if err != nil {
		t.Fatal(err)
	}
	if vc.Len() != len(rows) {
		t.Fatalf("row count %d != %d", vc.Len(), len(rows))
	}
	for r := range rows {
		got := vc.At(r)
		want := rows[r]
		if len(got) != len(want) {
			t.Fatalf("row %d: len %d != %d", r, len(got), len(want))
		}
		for i := range want {
			if got[i] != want[i] {
				t.Fatalf("row %d elem %d: %v vs %v", r, i, got[i], want[i])
			}
		}
	}
}

// TestVarColumnRoundTripInt32 exercises the int32 VLA path.
func TestVarColumnRoundTripInt32(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "vla_i32.fits")

	rows := [][]int32{
		{1, 2, 3},
		{42},
		{},
		{-1000, -2000, -3000, -4000},
	}

	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{})
	AppendBinaryTable(f, nil, []ColumnData{
		{Name: "IDS", DataVarInt32: rows},
	})
	f.Close()

	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()

	tbl := mustBinTable(t, fr, 1)
	vc, err := ReadVarColumn[int32](tbl, 1)
	if err != nil {
		t.Fatal(err)
	}
	if vc.Len() != len(rows) {
		t.Fatalf("row count")
	}
	for r := range rows {
		got := vc.At(r)
		if len(got) != len(rows[r]) {
			t.Fatalf("row %d len %d != %d", r, len(got), len(rows[r]))
		}
		for i := range rows[r] {
			if got[i] != rows[r][i] {
				t.Fatalf("row %d elem %d: %d vs %d", r, i, got[i], rows[r][i])
			}
		}
	}
}

// TestMixedFixedAndVarColumns writes a table that combines a fixed-width
// column with a variable-length one and verifies both read back correctly —
// this catches any off-by-one error in row byte layout when VLA descriptors
// sit alongside scalar cells.
func TestMixedFixedAndVarColumns(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "mixed.fits")

	ids := []int32{10, 20, 30}
	spectra := [][]float32{
		{1.1, 2.2, 3.3, 4.4},
		{5.5},
		{6.6, 7.7, 8.8},
	}

	f, _ := Create(path)
	WriteImage(f, nil, []int64{}, []float32{})
	AppendBinaryTable(f, nil, []ColumnData{
		{Name: "ID", DataInt32: ids},
		{Name: "SPEC", DataVarFloat32: spectra},
	})
	f.Close()

	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()

	tbl := mustBinTable(t, fr, 1)

	gotIDs, err := ReadColumn[int32](tbl, 1)
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range ids {
		if gotIDs[i] != v {
			t.Fatalf("id %d: %d vs %d", i, gotIDs[i], v)
		}
	}

	gotSpec, err := ReadVarColumn[float32](tbl, 2)
	if err != nil {
		t.Fatal(err)
	}
	for r := range spectra {
		got := gotSpec.At(r)
		if len(got) != len(spectra[r]) {
			t.Fatalf("spec %d len", r)
		}
		for i := range spectra[r] {
			if got[i] != spectra[r][i] {
				t.Fatalf("spec %d[%d]: %v vs %v", r, i, got[i], spectra[r][i])
			}
		}
	}
}

// mustBinTable is a tiny helper that fetches an HDU and asserts it is a
// binary table; keeps the tests free of the same boilerplate.
func mustBinTable(t *testing.T, f *File, idx int) *BinaryTableHDU {
	t.Helper()
	h, err := f.HDU(idx)
	if err != nil {
		t.Fatalf("HDU(%d): %v", idx, err)
	}
	tbl, ok := h.(*BinaryTableHDU)
	if !ok {
		t.Fatalf("HDU(%d) is %T, want *BinaryTableHDU", idx, h)
	}
	return tbl
}
