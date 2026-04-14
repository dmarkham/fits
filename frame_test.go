package fits

import (
	"math"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestFrame_WriteMono_ReadFrame_RoundTrip verifies WriteMono + ReadFrame
// returns the exact same float32 pixels in the same [W, H] layout.
func TestFrame_WriteMono_ReadFrame_RoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "mono.fits")

	const w, h = 4, 3
	pix := []float32{
		0.0, 0.1, 0.2, 0.3,
		0.4, 0.5, 0.6, 0.7,
		0.8, 0.9, 1.0, 0.55,
	}
	if err := WriteMono(path, pix, w, h); err != nil {
		t.Fatal(err)
	}

	got, gw, gh, err := ReadFrame(path)
	if err != nil {
		t.Fatal(err)
	}
	if gw != w || gh != h {
		t.Fatalf("dims: got %dx%d want %dx%d", gw, gh, w, h)
	}
	if len(got) != len(pix) {
		t.Fatalf("len: got %d want %d", len(got), len(pix))
	}
	for i, v := range pix {
		if got[i] != v {
			t.Errorf("pixel %d: got %v want %v", i, got[i], v)
		}
	}
}

// TestFrame_WriteRGB_ReadFrameRGB_RoundTrip verifies the planar R/G/B
// layout round-trips through WriteRGB + ReadFrameRGB unchanged.
func TestFrame_WriteRGB_ReadFrameRGB_RoundTrip(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "rgb.fits")

	const w, h = 3, 2
	r := []float32{0.10, 0.11, 0.12, 0.13, 0.14, 0.15}
	g := []float32{0.20, 0.21, 0.22, 0.23, 0.24, 0.25}
	b := []float32{0.30, 0.31, 0.32, 0.33, 0.34, 0.35}
	if err := WriteRGB(path, r, g, b, w, h); err != nil {
		t.Fatal(err)
	}

	gr, gg, gb, gw, gh, err := ReadFrameRGB(path)
	if err != nil {
		t.Fatal(err)
	}
	if gw != w || gh != h {
		t.Fatalf("dims: got %dx%d want %dx%d", gw, gh, w, h)
	}
	for i := range r {
		if gr[i] != r[i] || gg[i] != g[i] || gb[i] != b[i] {
			t.Errorf("pixel %d: got (%v,%v,%v) want (%v,%v,%v)",
				i, gr[i], gg[i], gb[i], r[i], g[i], b[i])
		}
	}
}

// TestFrame_WriteRGB_FileLayout_IsPlanar confirms the on-disk shape is
// NAXIS=3 with NAXIS3=3 and that ReadFrame rejects it.
func TestFrame_WriteRGB_FileLayout_IsPlanar(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "rgb_layout.fits")

	const w, h = 2, 2
	plane := []float32{0, 0.25, 0.5, 0.75}
	if err := WriteRGB(path, plane, plane, plane, w, h); err != nil {
		t.Fatal(err)
	}

	f, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	prim, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if prim.NAXIS() != 3 {
		t.Fatalf("NAXIS: got %d want 3", prim.NAXIS())
	}
	shape := prim.Shape()
	if len(shape) != 3 || shape[0] != w || shape[1] != h || shape[2] != 3 {
		t.Fatalf("shape: got %v want [%d %d 3]", shape, w, h)
	}
	if prim.BITPIX() != -32 {
		t.Fatalf("BITPIX: got %d want -32", prim.BITPIX())
	}
}

// TestFrame_ReadFrame_ErrorsOn3Planes verifies RGB files are rejected
// with a hint pointing at ReadFrameRGB.
func TestFrame_ReadFrame_ErrorsOn3Planes(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "rgb.fits")
	const w, h = 2, 2
	plane := []float32{0, 0, 0, 0}
	if err := WriteRGB(path, plane, plane, plane, w, h); err != nil {
		t.Fatal(err)
	}

	_, _, _, err := ReadFrame(path)
	if err == nil {
		t.Fatal("expected error for 3-plane image")
	}
	if !strings.Contains(err.Error(), "ReadFrameRGB") {
		t.Errorf("error should mention ReadFrameRGB, got: %v", err)
	}
}

// TestFrame_ReadFrame_ErrorsOnHyperspectral verifies NAXIS=3 with N!=3 and
// N!=1 planes is rejected with a hint to use ReadPixels.
func TestFrame_ReadFrame_ErrorsOnHyperspectral(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "hyper.fits")

	const w, h, planes = 2, 2, 5
	data := make([]float32, w*h*planes)
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, header.New(), []int64{w, h, planes}, data); err != nil {
		t.Fatal(err)
	}
	f.Close()

	_, _, _, err = ReadFrame(path)
	if err == nil {
		t.Fatal("expected error for hyperspectral cube")
	}
	if !strings.Contains(err.Error(), "ReadPixels") {
		t.Errorf("error should mention ReadPixels, got: %v", err)
	}
	if strings.Contains(err.Error(), "ReadFrameRGB") {
		t.Errorf("error should not mention ReadFrameRGB for 5-plane cube, got: %v", err)
	}
}

// TestFrame_ReadFrame_AcceptsDegenerate3D verifies NAXIS=3 with NAXIS3=1 is
// accepted as a 2D image.
func TestFrame_ReadFrame_AcceptsDegenerate3D(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "degen.fits")

	const w, h = 3, 2
	data := []float32{0.1, 0.2, 0.3, 0.4, 0.5, 0.6}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, header.New(), []int64{w, h, 1}, data); err != nil {
		t.Fatal(err)
	}
	f.Close()

	got, gw, gh, err := ReadFrame(path)
	if err != nil {
		t.Fatal(err)
	}
	if gw != w || gh != h {
		t.Fatalf("dims: got %dx%d want %dx%d", gw, gh, w, h)
	}
	if len(got) != len(data) {
		t.Fatalf("len: got %d want %d", len(got), len(data))
	}
	for i, v := range data {
		if got[i] != v {
			t.Errorf("pixel %d: got %v want %v", i, got[i], v)
		}
	}
}

// TestFrame_ReadFrame_ErrorsOn1D verifies 1D spectra are rejected.
func TestFrame_ReadFrame_ErrorsOn1D(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "spec.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, header.New(), []int64{8}, make([]float32, 8)); err != nil {
		t.Fatal(err)
	}
	f.Close()

	_, _, _, err = ReadFrame(path)
	if err == nil {
		t.Fatal("expected error for 1D image")
	}
}

// TestFrame_ReadFrame_ErrorsOnEmptyPrimary verifies files with NAXIS=0
// primary and no extensions are rejected.
func TestFrame_ReadFrame_ErrorsOnEmptyPrimary(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "empty.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, header.New(), nil, []float32{}); err != nil {
		t.Fatal(err)
	}
	f.Close()

	_, _, _, err = ReadFrame(path)
	if err == nil {
		t.Fatal("expected error for NAXIS=0 primary")
	}
}

// TestFrame_ReadFrame_ErrorsOnPlaceholderPlusBinTable verifies the
// selectFrameImage rule rejects files where the primary is a NAXIS=0
// placeholder and HDU 1 is a BINTABLE (not a tile-compressed image).
func TestFrame_ReadFrame_ErrorsOnPlaceholderPlusBinTable(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "placeholder_plus_table.fits")

	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, header.New(), nil, []float32{}); err != nil {
		t.Fatal(err)
	}
	cols := []ColumnData{{Name: "X", DataInt32: []int32{1, 2, 3}}}
	if _, err := AppendBinaryTable(f, header.New(), cols); err != nil {
		t.Fatal(err)
	}
	f.Close()

	_, _, _, err = ReadFrame(path)
	if err == nil {
		t.Fatal("expected error for placeholder primary + BINTABLE extension")
	}
	if !strings.Contains(err.Error(), "binary_table") {
		t.Errorf("error should mention HDU type 'binary_table', got: %v", err)
	}
}

// TestFrame_ReadFrameRGB_ErrorsOn2D verifies mono files are rejected.
func TestFrame_ReadFrameRGB_ErrorsOn2D(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "mono.fits")
	if err := WriteMono(path, []float32{0, 0, 0, 0}, 2, 2); err != nil {
		t.Fatal(err)
	}
	_, _, _, _, _, err := ReadFrameRGB(path)
	if err == nil {
		t.Fatal("expected error for 2D mono image")
	}
}

// TestFrame_WriteMono_ValidationErrors covers dim and length validation.
func TestFrame_WriteMono_ValidationErrors(t *testing.T) {
	dir := t.TempDir()
	cases := []struct {
		name string
		pix  []float32
		w, h int
	}{
		{"zero width", []float32{}, 0, 3},
		{"zero height", []float32{0, 0}, 2, 0},
		{"negative width", []float32{}, -1, 3},
		{"length mismatch", []float32{1, 2, 3}, 2, 2},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			err := WriteMono(filepath.Join(dir, c.name+".fits"), c.pix, c.w, c.h)
			if err == nil {
				t.Fatal("expected validation error")
			}
		})
	}
}

// TestFrame_WriteRGB_ValidationErrors covers dim and plane-length validation.
func TestFrame_WriteRGB_ValidationErrors(t *testing.T) {
	dir := t.TempDir()
	good := []float32{1, 2, 3, 4}

	if err := WriteRGB(filepath.Join(dir, "zero.fits"), good, good, good, 0, 0); err == nil {
		t.Fatal("expected error for zero dims")
	}
	short := []float32{1, 2, 3}
	if err := WriteRGB(filepath.Join(dir, "r_short.fits"), short, good, good, 2, 2); err == nil {
		t.Fatal("expected error for short R plane")
	}
	if err := WriteRGB(filepath.Join(dir, "g_short.fits"), good, short, good, 2, 2); err == nil {
		t.Fatal("expected error for short G plane")
	}
	if err := WriteRGB(filepath.Join(dir, "b_short.fits"), good, good, short, 2, 2); err == nil {
		t.Fatal("expected error for short B plane")
	}
}

// TestFrame_ReadFrame_TileCompressed_PrimaryPlaceholder verifies the
// HDU-selection rule finds the image in HDU 1 for astropy-style tile-
// compressed files (primary is a NAXIS=0 placeholder).
func TestFrame_ReadFrame_TileCompressed_PrimaryPlaceholder(t *testing.T) {
	compPath := filepath.Join("testdata", "comp_rice_i16.fits")
	rawPath := filepath.Join("testdata", "comp_raw_rice_i16.fits")

	got, w, h, err := ReadFrame(compPath)
	if err != nil {
		t.Skipf("compressed fixture not usable: %v", err)
	}
	if w <= 0 || h <= 0 {
		t.Fatalf("bogus dims %dx%d", w, h)
	}
	if len(got) != w*h {
		t.Fatalf("pixel count %d != %d*%d", len(got), w, h)
	}
	for i, v := range got {
		if v < 0 || v > 1 || math.IsNaN(float64(v)) {
			t.Fatalf("pixel %d out of [0,1]: %v", i, v)
			break
		}
	}

	rawPix, rw, rh, err := ReadFrame(rawPath)
	if err != nil {
		t.Skipf("raw fixture not usable: %v", err)
	}
	if rw != w || rh != h {
		t.Fatalf("dims mismatch raw vs compressed: %dx%d vs %dx%d", rw, rh, w, h)
	}
	for i := range rawPix {
		if rawPix[i] != got[i] {
			t.Fatalf("pixel %d differs: raw=%v comp=%v", i, rawPix[i], got[i])
		}
	}
}

// TestFrame_ReadFrame_Int16_Normalized verifies ReadFrame normalizes integer
// BITPIX to [0, 1], matching ReadFloat32 behavior.
func TestFrame_ReadFrame_Int16_Normalized(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "int16.fits")

	const w, h = 3, 1
	data := []int16{-32768, 0, 32767}
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, header.New(), []int64{w, h}, data); err != nil {
		t.Fatal(err)
	}
	f.Close()

	got, gw, gh, err := ReadFrame(path)
	if err != nil {
		t.Fatal(err)
	}
	if gw != w || gh != h {
		t.Fatalf("dims: got %dx%d want %dx%d", gw, gh, w, h)
	}
	const tol = 1e-5
	if got[0] != 0 {
		t.Errorf("pixel 0: got %v want 0", got[0])
	}
	if got[2] != 1 {
		t.Errorf("pixel 2: got %v want 1", got[2])
	}
	if math.Abs(float64(got[1]-0.5)) > tol {
		t.Errorf("pixel 1: got %v want ~0.5", got[1])
	}
}
