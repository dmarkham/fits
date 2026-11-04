package fits

import (
	"bytes"
	"encoding/json"
	"fmt"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestWriteImageMatchesCfitsioBytes is the core of the write-interop
// matrix: for every case in testdata/images/cases.json (every BITPIX, the
// int8/uint16/uint32/uint64 conventions, explicit BSCALE/BZERO with
// rounding up to the edges of the stored range, BLANK, 1-D to 4-D), the
// data unit the Go library writes must be byte-identical to the one
// cfitsio wrote for the same values (committed cfwrite/<case>.fits), and
// BITPIX/NAXISn/BZERO/BSCALE must agree with cfitsio's header. (Values
// outside the stored range are refused: TestWriteImageOverflowMatchesCfitsio.)
func TestWriteImageMatchesCfitsioBytes(t *testing.T) {
	dir := t.TempDir()
	cases := loadCases(t)
	// Caller headers with only one of BSCALE/BZERO must stay in the
	// matrix (the other defaults to 1 / 0; using the convention instead
	// would store the data unscaled).
	need := map[string]bool{"i16_bzero_only": true, "i16_bscale_only": true, "f64_bzero_only": true}
	for _, c := range cases {
		delete(need, c.Name)
	}
	for n := range need {
		t.Errorf("cases.json has no case %s", n)
	}
	for _, c := range cases {
		t.Run(c.Name, func(t *testing.T) {
			ref := imagesPath(t, "cfwrite", c.Name+".fits")
			testutil.RequireFile(t, ref)
			got := writeCaseFile(t, dir, c)
			_, gotData := dataUnit(t, got)
			wantData := expectedCaseData(t, c)
			if !bytes.Equal(gotData, wantData) {
				t.Errorf("data unit differs from cfitsio's (%s): first diff at byte %d\n go:      % x\n cfitsio: % x",
					c.Note, firstDiff(gotData, wantData), head(gotData, 48), head(wantData, 48))
			}
			// Header semantics through the Go reader of both files.
			g := openImage(t, got)
			w := openImage(t, ref)
			if g.BITPIX() != w.BITPIX() || bscaleOf(t, g) != bscaleOf(t, w) || bzeroOf(t, g) != bzeroOf(t, w) ||
				!equalShape(g.Shape(), w.Shape()) {
				t.Errorf("BITPIX/BSCALE/BZERO/shape = %d/%v/%v/%v, cfitsio %d/%v/%v/%v",
					g.BITPIX(), bscaleOf(t, g), bzeroOf(t, g), g.Shape(), w.BITPIX(), bscaleOf(t, w), bzeroOf(t, w), w.Shape())
			}
		})
	}
}

// TestWriteImageMatrixInRange: cfitsio wrote every matrix case with
// status 0 (no NUM_OVERFLOW), so the byte comparisons above compare values
// both writers store; the cases that reach the edges of the stored range
// are present by name (BITPIX 8 0..255 from both sides, lo-0.5+1 rounding
// to lo on BITPIX 8/16/32, int64 edges under BSCALE 0.5, negative and
// non-integer BSCALE on the unsigned and int8 conventions). The values
// beyond those edges are refused (TestWriteImageOverflowMatchesCfitsio).
func TestWriteImageMatrixInRange(t *testing.T) {
	m := loadManifest(t)
	cases := loadCases(t)
	have := map[string]bool{}
	for _, c := range cases {
		have[c.Name] = true
		st, ok := m.CfwriteStatus[c.Name]
		if !ok || st.Overflow || st.Status != 0 {
			t.Errorf("%s: cfitsio write status %+v, want 0", c.Name, st)
		}
	}
	for _, name := range []string{"u8_bzero_neg100", "u8_bzero_10", "i16_bzero1000", "i32_bscale_milli",
		"u8_bzero_half", "i16_bzero_half", "i32_bzero_half", "u16_bscale_1p1", "u32_bscale_1p1", "u64_bscale_1p1",
		"u64_bscale_neg1p5", "i8_bscale_neg1", "i64_bscale_half", "i64_bzero1000", "i16_bscale_neg1"} {
		if !have[name] {
			t.Errorf("cases.json has no case %s", name)
		}
	}
}

// TestWriteImageCfitsioReadsBack runs the cfitsio harness on every file
// the Go library wrote and requires cfitsio to see exactly what it sees in
// its own reference file: BITPIX, equivalent type, axes, BSCALE/BZERO,
// BLANK and every pixel (scaled double, int64, raw).
func TestWriteImageCfitsioReadsBack(t *testing.T) {
	dir := t.TempDir()
	for _, c := range loadCases(t) {
		t.Run(c.Name, func(t *testing.T) {
			got := testutil.CfrefDump(t, writeCaseFile(t, dir, c))
			want := loadImagesCfref(t, "cfwrite/"+c.Name+".fits")
			if got.NHDU != 1 || want.NHDU != 1 {
				t.Fatalf("nhdu = %d, cfitsio-written %d", got.NHDU, want.NHDU)
			}
			if len(got.HDUs[0].Errors) != len(want.HDUs[0].Errors) {
				t.Errorf("cfitsio errors on Go file %v, on its own file %v", got.HDUs[0].Errors, want.HDUs[0].Errors)
			}
			for _, d := range compareCfrefImages(got.HDUs[0].Image, want.HDUs[0].Image) {
				t.Error(d)
			}
			// uint64 values: cfref reads BITPIX 64 as TLONGLONG only, so
			// compare cfitsio's TULONGLONG reading (imgref pix) too.
			gp := imgrefPix(t, filepath.Join(dir, c.Name+".fits"))
			wp := loadPixGolden(t, "cfwrite/"+c.Name+".fits")
			if d := pixGoldenDiff(gp, wp); d != "" {
				t.Errorf("cfitsio TFLOAT/TULONGLONG reading of the Go file differs from its own file: %s", d)
			}
		})
	}
}

// TestWriteImageFitsverify requires HEASARC fitsverify to report no
// warnings and no errors for every Go-written matrix file.
func TestWriteImageFitsverify(t *testing.T) {
	dir := t.TempDir()
	for _, c := range loadCases(t) {
		p := writeCaseFile(t, dir, c)
		r := testutil.Fitsverify(t, p)
		if !r.OK() {
			t.Errorf("%s: fitsverify %d warnings, %d errors:\n%s", c.Name, r.Warnings, r.Errors, r.Report)
		}
	}
}

// TestWriteImageAstropyReadsBack: astropy (second independent reader)
// must decode every Go-written file to the same dtype and values as the
// cfitsio-written reference, with the same BITPIX/BZERO/BSCALE.
func TestWriteImageAstropyReadsBack(t *testing.T) {
	dir := t.TempDir()
	type pair struct {
		Name string `json:"name"`
		Go   string `json:"go"`
		Ref  string `json:"ref"`
	}
	var pairs []pair
	script := imagesPath(t, "astropy_compare.py")
	testutil.RequireFile(t, script)
	inputs := []string{script}
	for _, c := range loadCases(t) {
		pairs = append(pairs, pair{c.Name, writeCaseFile(t, dir, c), expectedCaseFile(t, dir, c)})
		inputs = append(inputs, pairs[len(pairs)-1].Go, pairs[len(pairs)-1].Ref)
	}
	in, _ := json.Marshal(pairs)
	out := testutil.RunOK(t, testutil.ToolAstropy, testutil.Call{Args: []string{script}, Inputs: inputs, Stdin: in}).Stdout
	var res map[string]struct {
		Equal  bool   `json:"equal"`
		Detail string `json:"detail"`
	}
	if err := json.Unmarshal(out, &res); err != nil {
		t.Fatalf("decode astropy output: %v\n%s", err, out)
	}
	if len(res) != len(pairs) {
		t.Fatalf("astropy compared %d files, want %d", len(res), len(pairs))
	}
	for _, p := range pairs {
		r, ok := res[p.Name]
		if !ok || !r.Equal {
			t.Errorf("%s: astropy reads the Go file differently from cfitsio's: %s", p.Name, r.Detail)
		}
	}
}

// TestOverwritePixelsMatchesCfitsio: OverwritePixels must apply the HDU's
// BSCALE/BZERO inverse (and the unsigned conventions) exactly like
// WriteImage. Each case is first written with its values reversed, then
// overwritten in place with the real values; the data unit must equal the
// cfitsio reference.
func TestOverwritePixelsMatchesCfitsio(t *testing.T) {
	dir := t.TempDir()
	for _, c := range loadCases(t) {
		t.Run(c.Name, func(t *testing.T) {
			path := filepath.Join(dir, c.Name+".fits")
			seed, _ := caseData(t, c, true)
			f, err := Create(path)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := seed(f, caseHeader(t, c), c.Naxes); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}
			e, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			img, err := e.Primary()
			if err != nil {
				t.Fatal(err)
			}
			_, over := caseData(t, c, false)
			if err := over(img); err != nil {
				t.Fatalf("OverwritePixels: %v", err)
			}
			if err := e.Close(); err != nil {
				t.Fatal(err)
			}
			_, gd := dataUnit(t, path)
			wd := expectedCaseData(t, c)
			if !bytes.Equal(gd, wd) {
				t.Errorf("overwritten data differs from cfitsio's at byte %d\n go:      % x\n cfitsio: % x",
					firstDiff(gd, wd), head(gd, 48), head(wd, 48))
			}
		})
	}
}

// --- helpers ---

func imgrefPix(t testing.TB, path string) *pixGolden {
	t.Helper()
	r := testutil.RunOK(t, toolImgref, testutil.Call{Args: []string{"pix", path}, Inputs: []string{path}})
	var g pixGolden
	dec := json.NewDecoder(bytes.NewReader(r.Stdout))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&g); err != nil {
		t.Fatalf("decode imgref pix: %v", err)
	}
	return &g
}

// pixGoldenDiff compares two imgref pix readings (either array form) and
// describes the first difference, "" when they agree.
func pixGoldenDiff(a, b *pixGolden) string {
	if a.Cfitsio != b.Cfitsio || len(a.HDUs) != len(b.HDUs) {
		return fmt.Sprintf("cfitsio %s / %d HDUs vs %s / %d HDUs", a.Cfitsio, len(a.HDUs), b.Cfitsio, len(b.HDUs))
	}
	for i := range a.HDUs {
		x, y := a.HDUs[i], b.HDUs[i]
		if x.Index != y.Index || x.F32Error != y.F32Error || x.F32Ovf != y.F32Ovf || x.U64Error != y.U64Error || x.F64NNError != y.F64NNError {
			return fmt.Sprintf("HDU %d: index/errors/overflow %+v vs %+v", i,
				[]any{x.Index, x.F32Error, x.F32Ovf, x.U64Error, x.F64NNError}, []any{y.Index, y.F32Error, y.F32Ovf, y.U64Error, y.F64NNError})
		}
		for _, d := range []struct{ what, diff string }{
			{"f32", x.F32.Same(y.F32)}, {"u64", x.U64.Same(y.U64)}, {"f64nn", x.F64NN.Same(y.F64NN)},
		} {
			if d.diff != "" {
				return fmt.Sprintf("HDU %d %s: %s", i, d.what, d.diff)
			}
		}
	}
	return ""
}

func openImage(t testing.TB, path string) *ImageHDU {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		t.Fatalf("Open %s: %v", path, err)
	}
	t.Cleanup(func() { f.Close() })
	img, err := f.Primary()
	if err != nil {
		t.Fatalf("Primary %s: %v", path, err)
	}
	return img
}

func equalShape(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func firstDiff(a, b []byte) int {
	n := min(len(a), len(b))
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return i
		}
	}
	if len(a) != len(b) {
		return n
	}
	return -1
}

func head(b []byte, n int) []byte { return b[:min(n, len(b))] }
