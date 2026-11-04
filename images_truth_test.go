package fits

// Shared loaders for the image / checksum truth corpus in testdata/images
// (samples written by cfitsio, by astropy, copied from astropy's test data,
// or byte-assembled; goldens produced by the cfitsio harnesses imgref and
// cfref). See testdata/images/PROVENANCE.md. Nothing in the expected values
// below comes from this library.

import (
	"bytes"
	"encoding/binary"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

func imagesPath(t testing.TB, elem ...string) string {
	t.Helper()
	return testutil.TestdataPath(t, append([]string{"images"}, elem...)...)
}

// imgCase is one entry of testdata/images/cases.json.
type imgCase struct {
	Name       string   `json:"name"`
	GoType     string   `json:"gotype"`
	CfImgtype  int      `json:"cf_imgtype"`
	CfDatatype int      `json:"cf_datatype"`
	Bscale     *float64 `json:"bscale"`
	Bzero      *float64 `json:"bzero"`
	Blank      *int64   `json:"blank"`
	Naxes      []int64  `json:"naxes"`
	Values     []string `json:"values"`
	Note       string   `json:"note"`
}

func loadCases(t testing.TB) []imgCase {
	t.Helper()
	p := imagesPath(t, "cases.json")
	testutil.RequireFile(t, p)
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	var doc struct {
		Comment string    `json:"comment"`
		Cases   []imgCase `json:"cases"`
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&doc); err != nil {
		t.Fatalf("decode %s: %v", p, err)
	}
	if len(doc.Cases) == 0 {
		t.Fatal("cases.json has no cases")
	}
	return doc.Cases
}

// imagesManifest is testdata/images/golden/manifest.json.
type imagesManifest struct {
	Astropy       string            `json:"astropy"`
	Numpy         string            `json:"numpy"`
	AstropySHA256 map[string]string `json:"astropy_data_sha256"`
	CfwriteStatus map[string]struct {
		Status   int  `json:"status"`
		Overflow bool `json:"overflow"`
	} `json:"cfwrite_status"`
	Samples []string `json:"samples"`
}

func loadManifest(t testing.TB) imagesManifest {
	t.Helper()
	var m imagesManifest
	loadJSON(t, imagesPath(t, "golden", "manifest.json"), &m)
	if len(m.Samples) == 0 {
		t.Fatal("manifest lists no samples")
	}
	return m
}

func loadJSON(t testing.TB, p string, v any) {
	t.Helper()
	b := testutil.MustReadGolden(t, p)
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	if err := dec.Decode(v); err != nil {
		t.Fatalf("decode %s: %v", p, err)
	}
}

// goldenName maps "cfwrite/u16_conv.fits" to its golden base name.
func goldenName(sample string) string { return strings.ReplaceAll(sample, "/", "__") }

func loadImagesCfref(t testing.TB, sample string) *testutil.CfrefFile {
	t.Helper()
	p := imagesPath(t, "golden", goldenName(sample)+".json")
	b := testutil.MustReadGolden(t, p)
	cf, err := testutil.DecodeCfref(b)
	if err != nil {
		t.Fatalf("decode %s: %v", p, err)
	}
	if cf.Cfitsio != "4.7.0" {
		t.Fatalf("%s made with cfitsio %s, want 4.7.0", p, cf.Cfitsio)
	}
	return cf
}

// chksumGolden is "imgref chksum" output.
type chksumGolden struct {
	File    string `json:"file"`
	Cfitsio string `json:"cfitsio"`
	HDUs    []struct {
		Index         int    `json:"index"`
		Datasum       uint32 `json:"datasum"`
		HDUSum        uint32 `json:"hdusum"`
		VerifyDatasum int    `json:"verify_datasum"`
		VerifyHDUSum  int    `json:"verify_hdusum"`
		VerifyError   string `json:"verify_error"` // fits_verify_chksum failed (e.g. 204)
		Error         string `json:"error"`
	} `json:"hdus"`
}

func decodeChksum(t testing.TB, b []byte, what string) chksumGolden {
	t.Helper()
	var g chksumGolden
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&g); err != nil {
		t.Fatalf("decode %s: %v", what, err)
	}
	return g
}

func loadChksumGolden(t testing.TB, sample string) chksumGolden {
	t.Helper()
	p := imagesPath(t, "golden", goldenName(sample)+".chksum.json")
	return decodeChksum(t, testutil.MustReadGolden(t, p), p)
}

// pixGolden is "imgref pix" output: cfitsio's TFLOAT reading of every
// image (each float32 as its exact double), its TULONGLONG reading of
// BITPIX 64 images, and, for integer images with null pixels, its TDOUBLE
// reading with null checking off (the scaled stored value of every pixel,
// BLANK ones included). Arrays over 1024 values are exact digests
// (testutil.Values; resolve them with a candidate reading).
type pixGolden struct {
	Cfitsio string   `json:"cfitsio"`
	HDUs    []pixHDU `json:"hdus"`
}

// pixHDU is one HDU of a pixGolden.
type pixHDU struct {
	Index      int                                  `json:"index"`
	F32        testutil.Values[testutil.NullFloat]  `json:"f32"`
	F32Error   string                               `json:"f32_error"`
	F32Ovf     bool                                 `json:"f32_overflow"`
	U64        testutil.Values[testutil.NullUint64] `json:"u64"`
	U64Error   string                               `json:"u64_error"`
	F64NN      testutil.Values[testutil.NullFloat]  `json:"f64nn"`
	F64NNError string                               `json:"f64nn_error"`
}

func loadPixGolden(t testing.TB, sample string) *pixGolden {
	t.Helper()
	var g pixGolden
	loadJSON(t, imagesPath(t, "golden", goldenName(sample)+".pix.json"), &g)
	if g.Cfitsio != "4.7.0" {
		t.Fatalf("%s pix golden made with cfitsio %s, want 4.7.0", sample, g.Cfitsio)
	}
	return &g
}

// cfitsioF32View is our masked float32 reading mapped into cfitsio's
// TFLOAT terms (the Resolve candidate of a digested "f32"): BLANK and NaN
// pixels null; for float BITPIX a stored +-Inf is null and a stored +-0 or
// subnormal reads as BZERO (fnan() underflow); when cfitsio reported
// NUM_OVERFLOW a +-Inf that came from a finite stored value is +-FLT_MAX.
// The digest comparison checks this mapping for every pixel.
func cfitsioF32View(t testing.TB, img *ImageHDU, overflow bool) []testutil.NullFloat {
	t.Helper()
	vals, mask, err := ReadPixelsMasked[float32](img)
	if err != nil {
		return nil
	}
	out := make([]testutil.NullFloat, len(vals))
	for k := range vals {
		if mask[k] {
			out[k] = testutil.NullFloat{V: float64(vals[k]), Valid: true}
		}
	}
	bp := img.BITPIX()
	var buf []byte
	size := 0
	if bp < 0 {
		size = -bp / 8
		buf = make([]byte, len(vals)*size)
		if _, err := img.rec.file.rs.Seek(img.rec.dataStart, 0); err != nil {
			t.Fatal(err)
		}
		if _, err := io.ReadFull(img.rec.file.rs, buf); err != nil {
			t.Fatal(err)
		}
	}
	for k := range out {
		if size > 0 {
			var exp, max uint64
			if size == 4 {
				exp, max = uint64(binary.BigEndian.Uint32(buf[4*k:])>>23)&0xff, 0xff
			} else {
				exp, max = (binary.BigEndian.Uint64(buf[8*k:])>>52)&0x7ff, 0x7ff
			}
			switch exp {
			case max:
				out[k] = testutil.NullFloat{}
				continue
			case 0:
				out[k] = testutil.NullFloat{V: float64(float32(bzeroOf(t, img))), Valid: true}
				continue
			}
		}
		if overflow && out[k].Valid && math.IsInf(out[k].V, 0) {
			out[k].V = math.Copysign(math.MaxFloat32, out[k].V)
		}
	}
	return out
}

// toolImgref is the cfitsio harness testdata/images/imgref (recorded
// readings in tier 1, testutil.Run).
var toolImgref = testutil.Tool{Name: "imgref", Harness: "testdata/images/imgref", Sources: []string{"testdata/images/imgref.c"}}

// imgrefChksum is "imgref chksum" on a file the Go code wrote.
func imgrefChksum(t testing.TB, path string) chksumGolden {
	t.Helper()
	r := testutil.RunOK(t, toolImgref, testutil.Call{Args: []string{"chksum", path}, Inputs: []string{path}})
	return decodeChksum(t, r.Stdout, path)
}

// dataUnit returns the bytes of the first HDU's data unit (header located by
// scanning for the END card on 80-byte boundaries, independent of the
// library's scanner) through the end of its padding.
func dataUnit(t testing.TB, path string) (hdr, data []byte) {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for off := 0; off+80 <= len(b); off += 80 {
		if string(b[off:off+80]) == "END"+strings.Repeat(" ", 77) {
			end := (off/2880 + 1) * 2880
			return b[:end], b[end:]
		}
	}
	t.Fatalf("%s: no END card", path)
	return nil, nil
}

// caseHeader builds the user header for a case: BSCALE/BZERO/BLANK only.
func caseHeader(t testing.TB, c imgCase) *header.Header {
	t.Helper()
	if c.Bscale == nil && c.Bzero == nil && c.Blank == nil {
		return nil
	}
	h := header.New()
	if c.Bscale != nil {
		if err := h.Set("BSCALE", *c.Bscale, ""); err != nil {
			t.Fatal(err)
		}
	}
	if c.Bzero != nil {
		if err := h.Set("BZERO", *c.Bzero, ""); err != nil {
			t.Fatal(err)
		}
	}
	if c.Blank != nil {
		if err := h.Set("BLANK", *c.Blank, ""); err != nil {
			t.Fatal(err)
		}
	}
	return h
}

func parseIntVals[T int8 | int16 | int32 | int64](t testing.TB, vals []string, bits int) []T {
	t.Helper()
	out := make([]T, len(vals))
	for i, s := range vals {
		v, err := strconv.ParseInt(s, 10, bits)
		if err != nil {
			t.Fatalf("value %q: %v", s, err)
		}
		out[i] = T(v)
	}
	return out
}

func parseUintVals[T uint8 | uint16 | uint32 | uint64](t testing.TB, vals []string, bits int) []T {
	t.Helper()
	out := make([]T, len(vals))
	for i, s := range vals {
		v, err := strconv.ParseUint(s, 10, bits)
		if err != nil {
			t.Fatalf("value %q: %v", s, err)
		}
		out[i] = T(v)
	}
	return out
}

func parseFloatVals[T float32 | float64](t testing.TB, vals []string, bits int) []T {
	t.Helper()
	out := make([]T, len(vals))
	for i, s := range vals {
		if h, ok := strings.CutPrefix(s, "bits:"); ok {
			// The exact IEEE bit pattern (signaling NaNs, payloads).
			u, err := strconv.ParseUint(strings.TrimPrefix(h, "0x"), 16, bits)
			if err != nil {
				t.Fatalf("value %q: %v", s, err)
			}
			if bits == 32 {
				out[i] = T(math.Float32frombits(uint32(u)))
			} else {
				out[i] = T(math.Float64frombits(u))
			}
			continue
		}
		if s == "nan" {
			// The quiet NaN C's strtod/strtof return (Go's math.NaN has a
			// different payload); the matrix compares bytes.
			if bits == 32 {
				out[i] = T(math.Float32frombits(0x7FC00000))
			} else {
				out[i] = T(math.Float64frombits(0x7FF8000000000000))
			}
			continue
		}
		v, err := strconv.ParseFloat(s, bits)
		if err != nil {
			t.Fatalf("value %q: %v", s, err)
		}
		out[i] = T(v)
	}
	return out
}

// caseWriter writes the case with the Go library through write(T-typed data).
type caseWriter func(dst WriteTarget, hdr *header.Header, shape []int64) (*ImageHDU, error)

// caseData returns a writer bound to the case's parsed values. When
// reversed is true the values are written in reverse order (used to seed a
// file that OverwritePixels then fixes).
func caseData(t testing.TB, c imgCase, reversed bool) (caseWriter, func(h *ImageHDU) error) {
	t.Helper()
	vals := append([]string(nil), c.Values...)
	if reversed {
		for i, j := 0, len(vals)-1; i < j; i, j = i+1, j-1 {
			vals[i], vals[j] = vals[j], vals[i]
		}
	}
	switch c.GoType {
	case "uint8":
		return bindCase(parseUintVals[uint8](t, vals, 8))
	case "int8":
		return bindCase(parseIntVals[int8](t, vals, 8))
	case "int16":
		return bindCase(parseIntVals[int16](t, vals, 16))
	case "uint16":
		return bindCase(parseUintVals[uint16](t, vals, 16))
	case "int32":
		return bindCase(parseIntVals[int32](t, vals, 32))
	case "uint32":
		return bindCase(parseUintVals[uint32](t, vals, 32))
	case "int64":
		return bindCase(parseIntVals[int64](t, vals, 64))
	case "uint64":
		return bindCase(parseUintVals[uint64](t, vals, 64))
	case "float32":
		return bindCase(parseFloatVals[float32](t, vals, 32))
	case "float64":
		return bindCase(parseFloatVals[float64](t, vals, 64))
	}
	t.Fatalf("case %s: unknown gotype %s", c.Name, c.GoType)
	return nil, nil
}

func bindCase[T Numeric](data []T) (caseWriter, func(h *ImageHDU) error) {
	return func(dst WriteTarget, hdr *header.Header, shape []int64) (*ImageHDU, error) {
			return WriteImage(dst, hdr, shape, data)
		}, func(h *ImageHDU) error {
			return OverwritePixels(h, data)
		}
}

// writeCaseFile writes case c to a new file in dir with the Go library.
func writeCaseFile(t testing.TB, dir string, c imgCase) string {
	t.Helper()
	path := filepath.Join(dir, c.Name+".fits")
	w, _ := caseData(t, c, false)
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := w(f, caseHeader(t, c), c.Naxes); err != nil {
		f.Close()
		t.Fatalf("%s: WriteImage: %v", c.Name, err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

// sameFloat compares physical values: equal, or both NaN.
func sameFloat(a, b float64) bool {
	return a == b && math.Signbit(a) == math.Signbit(b) || math.IsNaN(a) && math.IsNaN(b)
}

// compareCfrefImages reports the differences between two cfitsio views of
// an image (bitpix, equivalent type, axes, scaling, every pixel).
func compareCfrefImages(got, want *testutil.CfrefImage) []string {
	var d []string
	if got == nil || want == nil {
		return []string{fmt.Sprintf("image missing: got %v want %v", got != nil, want != nil)}
	}
	if got.Bitpix != want.Bitpix || got.Equivtype != want.Equivtype || got.Naxis != want.Naxis ||
		fmt.Sprint(got.Naxes) != fmt.Sprint(want.Naxes) {
		d = append(d, fmt.Sprintf("bitpix/equivtype/naxes = %d/%d/%v, cfitsio-written %d/%d/%v",
			got.Bitpix, got.Equivtype, got.Naxes, want.Bitpix, want.Equivtype, want.Naxes))
	}
	if got.Bscale != want.Bscale || got.Bzero != want.Bzero {
		d = append(d, fmt.Sprintf("bscale/bzero = %v/%v, cfitsio-written %v/%v", got.Bscale, got.Bzero, want.Bscale, want.Bzero))
	}
	if (got.Blank == nil) != (want.Blank == nil) || got.Blank != nil && *got.Blank != *want.Blank {
		d = append(d, "BLANK differs")
	}
	// Every pixel, exactly (explicit values or the sha256 digest of all of
	// them; floats bit for bit, every NaN alike).
	for _, c := range []struct{ name, diff string }{
		{"phys", got.Phys.Same(want.Phys)},
		{"i64", got.I64.Same(want.I64)},
		{"u64", got.U64.Same(want.U64)},
		{"raw", got.Raw.Same(want.Raw)},
	} {
		if c.diff != "" {
			d = append(d, fmt.Sprintf("%s: %s (cfitsio-written)", c.name, c.diff))
		}
	}
	return d
}

// imageSamples lists every committed testdata/images sample plus the
// top-level testdata samples, as (display name, path, cfitsio golden).
type imageSample struct {
	name   string // e.g. "cfwrite/u16_conv.fits" or "testprog.std"
	path   string
	golden *testutil.CfrefFile
	pix    *pixGolden
}

func allImageSamples(t testing.TB) []imageSample {
	t.Helper()
	var out []imageSample
	for _, s := range loadManifest(t).Samples {
		p := imagesPath(t, filepath.FromSlash(s))
		testutil.RequireFile(t, p)
		out = append(out, imageSample{name: s, path: p, golden: loadImagesCfref(t, s), pix: loadPixGolden(t, s)})
	}
	entries, err := os.ReadDir(testutil.TestdataPath(t))
	if err != nil {
		t.Fatal(err)
	}
	var top []string
	for _, e := range entries {
		n := e.Name()
		if !e.IsDir() && (strings.HasSuffix(n, ".fits") || strings.HasSuffix(n, ".fit") || strings.HasSuffix(n, ".std")) {
			top = append(top, n)
		}
	}
	sort.Strings(top)
	if len(top) == 0 {
		t.Fatal("no top-level testdata samples")
	}
	for _, n := range top {
		s := imageSample{name: n, path: testutil.TestdataPath(t, n), golden: testutil.LoadCfrefGolden(t, n)}
		// gen.py writes a pix golden only for a sample with a plain image
		// that has pixels (the HDUs forEachImage reads); the pixels of
		// tile-compressed HDUs are proven by the compressed-read corpus.
		pp := imagesPath(t, "golden", goldenName("testdata/"+n)+".pix.json")
		if hasPlainImagePixels(s.golden) {
			s.pix = loadPixGolden(t, "testdata/"+n)
		} else if testutil.GoldenExists(pp) {
			t.Errorf("%s: pix golden exists but the sample has no plain image with pixels (gen.py out of sync)", pp)
		}
		out = append(out, s)
	}
	return out
}

// hasPlainImagePixels reports whether a cfref golden has an image HDU that
// is not tile-compressed and has pixels (the HDUs forEachImage visits).
func hasPlainImagePixels(g *testutil.CfrefFile) bool {
	for _, h := range g.HDUs {
		if h.Image != nil && !h.Image.Compressed && h.Image.Npix > 0 {
			return true
		}
	}
	return false
}

// rawFloatPixel decodes float pixel k of an image straight from the file
// bytes (big-endian IEEE, independent of the library's decoder) and applies
// BSCALE/BZERO per the FITS definition. It is used where cfitsio's
// null-checked reads do not report the stored value: cfitsio flags +-Inf as
// null like NaN and returns +0 for -0 and subnormals ("underflow", fnan()
// in fitsio2.h, used by fffr4r8/fffr8r8).
func rawFloatPixel(t testing.TB, img *ImageHDU, k int) (raw, phys float64) {
	t.Helper()
	size := 4
	if img.BITPIX() == -64 {
		size = 8
	}
	buf := make([]byte, size)
	if _, err := img.rec.file.rs.Seek(img.rec.dataStart+int64(k*size), 0); err != nil {
		t.Fatal(err)
	}
	if _, err := io.ReadFull(img.rec.file.rs, buf); err != nil {
		t.Fatal(err)
	}
	if size == 4 {
		raw = float64(math.Float32frombits(binary.BigEndian.Uint32(buf)))
	} else {
		raw = math.Float64frombits(binary.BigEndian.Uint64(buf))
	}
	if bscaleOf(t, img) == 1 && bzeroOf(t, img) == 0 {
		return raw, raw // no scaling defined: the stored value (keeps -0)
	}
	return raw, raw*bscaleOf(t, img) + bzeroOf(t, img)
}

// cfitsioPhysView is our masked float64 reading of a plain image mapped
// into cfitsio's TDOUBLE terms (testutil.Values.Resolve candidate for a
// digested golden): BLANK and NaN pixels null; for float BITPIX, cfitsio
// also flags +-Inf null and returns BZERO for a stored +-0 or subnormal
// ("underflow", fnan() in fitsio2.h, used by fffr4r8/fffr8r8). The
// resolution compares every pixel with the golden's sha256, so this
// mapping is checked, not assumed.
func cfitsioPhysView(t testing.TB, img *ImageHDU) []testutil.NullFloat {
	t.Helper()
	vals, mask, err := ReadPixelsMasked[float64](img)
	if err != nil {
		return nil
	}
	out := testutil.NullFloats(vals, mask)
	bp := img.BITPIX()
	if bp > 0 {
		return out
	}
	size := -bp / 8
	buf := make([]byte, len(vals)*size)
	if _, err := img.rec.file.rs.Seek(img.rec.dataStart, 0); err != nil {
		t.Fatal(err)
	}
	if _, err := io.ReadFull(img.rec.file.rs, buf); err != nil {
		t.Fatal(err)
	}
	for k := range out {
		var exp, max uint64
		if size == 4 {
			exp, max = uint64(binary.BigEndian.Uint32(buf[4*k:])>>23)&0xff, 0xff
		} else {
			exp, max = (binary.BigEndian.Uint64(buf[8*k:])>>52)&0x7ff, 0x7ff
		}
		switch exp {
		case max:
			out[k] = testutil.NullFloat{}
		case 0:
			out[k] = testutil.NullFloat{V: bzeroOf(t, img), Valid: true}
		}
	}
	return out
}

// int64View is our masked int64 reading (cfitsio TLONGLONG terms), nil
// when the library refuses int64 for the image.
func int64View(img *ImageHDU) []testutil.NullInt64 {
	v, ok, err := ReadPixelsMasked[int64](img)
	if err != nil {
		return nil
	}
	return testutil.NullInt64s(v, ok)
}

// cfitsioFloatQuirk reports whether a float pixel that cfitsio returned as
// null or as +0 is one of cfitsio's documented substitutions (Inf -> null,
// -0/subnormal -> 0), returning the stored physical value to compare with.
func cfitsioFloatQuirk(t testing.TB, img *ImageHDU, k int, want testutil.NullFloat) (phys float64, quirk bool) {
	t.Helper()
	if img.BITPIX() > 0 {
		return 0, false
	}
	raw, phys := rawFloatPixel(t, img, k)
	switch {
	case !want.Valid && math.IsInf(raw, 0):
		return phys, true
	case want.Valid && want.V == 0 && !math.Signbit(want.V) && raw != 0 || want.Valid && want.V == 0 && math.Signbit(raw):
		// -0 or a subnormal stored; cfitsio reports +0
		if raw == 0 || math.Abs(raw) < 0x1p-126 {
			return phys, true
		}
	}
	return 0, false
}

// expectedCaseData returns the data unit the library must write for case
// c: cfitsio's committed data unit (cfwrite/<case>.fits). Every matrix
// value is storable (TestWriteImageMatrixInRange), so there is no element
// where the library and cfitsio differ.
func expectedCaseData(t testing.TB, c imgCase) []byte {
	t.Helper()
	ref := expectedCaseFile(t, "", c)
	_, data := dataUnit(t, ref)
	return data
}

// expectedCaseFile is the committed cfitsio reference of case c, which
// the read-back comparisons read live.
func expectedCaseFile(t testing.TB, _ string, c imgCase) string {
	t.Helper()
	ref := imagesPath(t, "cfwrite", c.Name+".fits")
	testutil.RequireFile(t, ref)
	return ref
}
