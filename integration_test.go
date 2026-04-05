package fits

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestIntegrationWriteThenReadSemantic is the end-to-end sanity test that
// limitation #2 of v1 was about: we build a file entirely through the
// library's writer API, exercising every HDU type and every column type,
// close it, reopen it through the library's reader API, and assert that
// every header keyword and every data element matches the original.
//
// This proves the encoder/decoder are self-consistent. It does NOT prove
// third-party compatibility (that's limitation #3, which needs fitsverify
// and astropy). It DOES prove that any file our writer produces can be
// round-tripped through our reader with byte-exact fidelity on the content.
func TestIntegrationWriteThenReadSemantic(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "integration.fits")

	// -------- Inputs we'll cross-check after read-back. --------
	primaryData := make([]float32, 4*3)
	for i := range primaryData {
		primaryData[i] = float32(i) * 0.5
	}
	imgExtData := []int16{-100, -50, 0, 50, 100, 200, 300, 400}
	tblIDs := []int32{1, 2, 3, 4, 5}
	tblFlux := []float64{1.5, 2.5, math.Pi, math.E, -0.125}
	tblLabels := []string{"alpha", "beta", "gamma", "delta", "omega"}
	tblFlags := []bool{true, false, true, true, false}
	tblSpectra := [][]float32{
		{1, 2, 3},
		{4, 5, 6, 7, 8},
		{},
		{9},
		{10, 11, 12, 13, 14, 15},
	}

	// -------- Write phase. --------
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}

	// HDU 0: primary image (float32, 2D) with a user keyword.
	primHdr := header.New()
	if err := primHdr.Set("OBJECT", "integration", "test object"); err != nil {
		t.Fatal(err)
	}
	if err := primHdr.Set("OBSERVER", "unit test", ""); err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, primHdr, []int64{4, 3}, primaryData); err != nil {
		t.Fatal(err)
	}

	// HDU 1: image extension (int16, 1D) with EXTNAME.
	extHdr := header.New()
	extHdr.Set("EXTNAME", "SIGNAL", "")
	if _, err := AppendImage(f, extHdr, []int64{8}, imgExtData); err != nil {
		t.Fatal(err)
	}

	// HDU 2: binary table mixing fixed-width and variable-length columns.
	tblHdr := header.New()
	tblHdr.Set("EXTNAME", "CATALOG", "")
	tblHdr.Set("AUTHOR", "test", "")
	cols := []ColumnData{
		{Name: "ID", DataInt32: tblIDs},
		{Name: "FLUX", Unit: "erg/s/cm2", DataFloat64: tblFlux},
		{Name: "LABEL", DataString: tblLabels},
		{Name: "VALID", DataBool: tblFlags},
		{Name: "SPEC", Unit: "counts", DataVarFloat32: tblSpectra},
	}
	if _, err := AppendBinaryTable(f, tblHdr, cols); err != nil {
		t.Fatal(err)
	}

	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	// -------- Read phase. --------
	fr, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer fr.Close()

	if fr.NumHDU() != 3 {
		t.Fatalf("NumHDU = %d, want 3", fr.NumHDU())
	}

	// -------- HDU 0: primary image. --------
	prim, err := fr.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if prim.BITPIX() != -32 {
		t.Fatalf("primary BITPIX = %d, want -32", prim.BITPIX())
	}
	if prim.NAXIS() != 2 {
		t.Fatalf("primary NAXIS = %d", prim.NAXIS())
	}
	sh := prim.Shape()
	if len(sh) != 2 || sh[0] != 4 || sh[1] != 3 {
		t.Fatalf("primary shape = %v", sh)
	}
	if obj, err := prim.Header().String("OBJECT"); err != nil || obj != "integration" {
		t.Fatalf("OBJECT = %q, err %v", obj, err)
	}
	if obs, err := prim.Header().String("OBSERVER"); err != nil || obs != "unit test" {
		t.Fatalf("OBSERVER = %q, err %v", obs, err)
	}
	gotPrim, err := ReadPixels[float32](prim)
	if err != nil {
		t.Fatal(err)
	}
	if len(gotPrim) != len(primaryData) {
		t.Fatalf("primary pixel count = %d, want %d", len(gotPrim), len(primaryData))
	}
	for i := range primaryData {
		if gotPrim[i] != primaryData[i] {
			t.Fatalf("primary pixel %d: %v vs %v", i, gotPrim[i], primaryData[i])
		}
	}

	// -------- HDU 1: image extension by EXTNAME. --------
	ext, err := fr.HDUByName("SIGNAL")
	if err != nil {
		t.Fatal(err)
	}
	extImg, ok := ext.(*ImageHDU)
	if !ok {
		t.Fatalf("SIGNAL HDU is %T", ext)
	}
	if extImg.BITPIX() != 16 {
		t.Fatalf("SIGNAL BITPIX = %d", extImg.BITPIX())
	}
	gotExt, err := ReadPixels[int16](extImg)
	if err != nil {
		t.Fatal(err)
	}
	if len(gotExt) != len(imgExtData) {
		t.Fatalf("SIGNAL len")
	}
	for i := range imgExtData {
		if gotExt[i] != imgExtData[i] {
			t.Fatalf("SIGNAL pixel %d: %d vs %d", i, gotExt[i], imgExtData[i])
		}
	}

	// -------- HDU 2: binary table with all column types. --------
	catAny, err := fr.HDUByName("CATALOG")
	if err != nil {
		t.Fatal(err)
	}
	cat, ok := catAny.(*BinaryTableHDU)
	if !ok {
		t.Fatalf("CATALOG HDU is %T", catAny)
	}
	if cat.NumRows() != int64(len(tblIDs)) {
		t.Fatalf("CATALOG rows = %d", cat.NumRows())
	}
	if auth, err := cat.Header().String("AUTHOR"); err != nil || auth != "test" {
		t.Fatalf("AUTHOR = %q %v", auth, err)
	}

	// ID
	gotIDs, err := ReadColumn[int32](cat, cat.ColumnIndex("ID"))
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range tblIDs {
		if gotIDs[i] != v {
			t.Fatalf("ID[%d]: %d vs %d", i, gotIDs[i], v)
		}
	}

	// FLUX
	gotFlux, err := ReadColumn[float64](cat, cat.ColumnIndex("FLUX"))
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range tblFlux {
		if gotFlux[i] != v {
			t.Fatalf("FLUX[%d]: %v vs %v", i, gotFlux[i], v)
		}
	}

	// Verify TUNIT survived.
	if unit, err := cat.Header().String("TUNIT2"); err != nil || unit != "erg/s/cm2" {
		t.Fatalf("TUNIT2 = %q %v", unit, err)
	}

	// LABEL (strings, padded in storage)
	gotLabels, err := ReadColumn[string](cat, cat.ColumnIndex("LABEL"))
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range tblLabels {
		if gotLabels[i] != v {
			t.Fatalf("LABEL[%d]: %q vs %q", i, gotLabels[i], v)
		}
	}

	// VALID (logical)
	gotFlags, err := ReadColumn[bool](cat, cat.ColumnIndex("VALID"))
	if err != nil {
		t.Fatal(err)
	}
	for i, v := range tblFlags {
		if gotFlags[i] != v {
			t.Fatalf("VALID[%d]: %v vs %v", i, gotFlags[i], v)
		}
	}

	// SPEC (VLA float32)
	gotSpec, err := ReadVarColumn[float32](cat, cat.ColumnIndex("SPEC"))
	if err != nil {
		t.Fatal(err)
	}
	if gotSpec.Len() != len(tblSpectra) {
		t.Fatalf("SPEC rows")
	}
	for r := range tblSpectra {
		got := gotSpec.At(r)
		if len(got) != len(tblSpectra[r]) {
			t.Fatalf("SPEC row %d len: %d vs %d", r, len(got), len(tblSpectra[r]))
		}
		for i := range tblSpectra[r] {
			if got[i] != tblSpectra[r][i] {
				t.Fatalf("SPEC row %d elem %d: %v vs %v", r, i, got[i], tblSpectra[r][i])
			}
		}
	}
}

// TestIntegrationWriterIsDeterministic writes the same logical content to
// two separate files and asserts the resulting bytes are identical. This
// is the v1 encoder's determinism guarantee — without it, the fixture
// round-trip regression tests would eventually flake on keyword ordering,
// float formatting, or map-iteration order.
func TestIntegrationWriterIsDeterministic(t *testing.T) {
	dir := t.TempDir()
	pathA := filepath.Join(dir, "a.fits")
	pathB := filepath.Join(dir, "b.fits")

	write := func(path string) {
		f, err := Create(path)
		if err != nil {
			t.Fatal(err)
		}
		hdr := header.New()
		hdr.Set("OBJECT", "det test", "")
		hdr.Set("EXPTIME", 120.5, "seconds")
		hdr.Set("IMAGETYP", "LIGHT", "")
		data := make([]float32, 64*48)
		for i := range data {
			data[i] = float32(i) * 0.25
		}
		if _, err := WriteImage(f, hdr, []int64{64, 48}, data); err != nil {
			t.Fatal(err)
		}
		tblHdr := header.New()
		tblHdr.Set("EXTNAME", "META", "")
		if _, err := AppendBinaryTable(f, tblHdr, []ColumnData{
			{Name: "IDX", DataInt32: []int32{1, 2, 3}},
			{Name: "TAG", DataString: []string{"aa", "bb", "cc"}},
			{Name: "VLA", DataVarInt32: [][]int32{{1, 2}, {3}, {4, 5, 6}}},
		}); err != nil {
			t.Fatal(err)
		}
		if err := f.Close(); err != nil {
			t.Fatal(err)
		}
	}
	write(pathA)
	write(pathB)

	bytesA, err := os.ReadFile(pathA)
	if err != nil {
		t.Fatal(err)
	}
	bytesB, err := os.ReadFile(pathB)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(bytesA, bytesB) {
		if len(bytesA) != len(bytesB) {
			t.Fatalf("writer not deterministic: file sizes %d vs %d", len(bytesA), len(bytesB))
		}
		for i := range bytesA {
			if bytesA[i] != bytesB[i] {
				t.Fatalf("writer not deterministic: first diff at byte %d: %#x vs %#x", i, bytesA[i], bytesB[i])
			}
		}
	}
}

// TestIntegrationWriteReadCopyByteExact writes a file via the library,
// then opens it and streams it through (*File).CopyTo to a fresh buffer,
// and asserts the bytes match. This proves the reader's CopyTo path is
// byte-exact for files produced by our own writer (complementing the
// cfitsio fixtures test which uses externally-produced inputs).
func TestIntegrationWriteReadCopyByteExact(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "wrc.fits")

	f, _ := Create(path)
	WriteImage(f, nil, []int64{16, 16}, makeFloat32Ramp(16*16))
	AppendBinaryTable(f, nil, []ColumnData{
		{Name: "A", DataInt32: []int32{10, 20, 30}},
		{Name: "B", DataVarFloat64: [][]float64{{1.1, 2.2}, {}, {3.3}}},
	})
	f.Close()

	orig, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}

	fr, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()

	var buf bytes.Buffer
	if err := fr.CopyTo(&buf); err != nil {
		t.Fatal(err)
	}
	out := buf.Bytes()
	if !bytes.Equal(orig, out) {
		n := len(orig)
		if len(out) < n {
			n = len(out)
		}
		for i := 0; i < n; i++ {
			if orig[i] != out[i] {
				t.Fatalf("CopyTo diverges at byte %d: %#x vs %#x (lens %d, %d)", i, orig[i], out[i], len(orig), len(out))
			}
		}
		t.Fatalf("CopyTo lengths differ: %d vs %d", len(orig), len(out))
	}
}

// makeFloat32Ramp returns n floats from 0..n-1.
func makeFloat32Ramp(n int) []float32 {
	out := make([]float32, n)
	for i := range out {
		out[i] = float32(i)
	}
	return out
}
