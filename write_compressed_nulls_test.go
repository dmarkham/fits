package fits

import (
	"bytes"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/internal/crefgen"
	"github.com/dmarkham/fits/internal/testutil"
)

// Null handling and header rules of the compressed writers that the
// golden parity test alone does not explain.

// TestFloatBLANKNullsRealPixels is why AppendCompressedFloat32Image
// refuses a caller BLANK (TestWriteCompressedRefusals "caller BLANK on
// float32"). Golden rice_f32_blank_hdr is cfitsio's own compression of a
// float32 source image carrying BLANK = 5 and no NaN: cfitsio copied
// BLANK into the ZBITPIX -32 header (no ZBLANK, as no tile has a null),
// and its decode of its own file nulls every real pixel whose quantized
// integer is 5. astropy reads the same file with the same nulls. A
// writer that passed BLANK through would lose those pixels silently.
func TestFloatBLANKNullsRealPixels(t *testing.T) {
	g := imgGoldenByName(t, "rice_f32_blank_hdr")
	for i, v := range g.floatInput(t) {
		if math.IsNaN(float64(v)) || math.IsInf(float64(v), 0) {
			t.Fatalf("input pixel %d = %v: the case needs an input without nulls", i, v)
		}
	}
	nulls := g.Dec.NullCount()
	if nulls == 0 {
		t.Fatalf("cfitsio's decode of its BLANK file has no nulls; the hazard this refusal guards against is gone")
	}
	f, h := openCompressed(t, filepath.Join(imgGoldenDir, g.Name+".fits"))
	defer f.Close()
	if c, ok := h.Header().Get("BLANK"); !ok || h.Header().Has("ZBLANK") {
		t.Fatalf("cfitsio's file: BLANK %v (present %v), ZBLANK present %v; want BLANK copied and no ZBLANK", c.Value, ok, h.Header().Has("ZBLANK"))
	}
	t.Logf("cfitsio nulls %d/%d real pixels of its own BLANK = 5 float file", nulls, g.npix())

	const script = `
import sys, json, warnings
import numpy as np
from astropy.io import fits
warnings.simplefilter("ignore")
with fits.open(sys.argv[1]) as h:
    d = np.asarray(h[1].data, dtype=np.float64).ravel(order="C")
json.dump([int(i) for i in np.flatnonzero(np.isnan(d))], sys.stdout)
`
	out := testutil.AstropyCode(t, script, goldenFITS(t, filepath.Join(testutil.RepoRoot(t), imgGoldenDir, g.Name+".fits")))
	var idx []int
	if err := json.Unmarshal(out, &idx); err != nil {
		t.Fatalf("parse astropy output: %v", err)
	}
	if len(idx) != nulls {
		t.Fatalf("astropy nulls %d pixels, cfitsio %d", len(idx), nulls)
	}
	for _, k := range idx {
		if !g.Dec.IsNull(k) {
			t.Fatalf("astropy nulls pixel %d, cfitsio does not", k)
		}
	}
}

// floatNullCase is a 64-pixel float32 row with a NaN at 3, cfitsio's
// FLOATNULLVALUE (-9.11912e-36) as a real pixel at 10 and real -1e30
// pixels at 20 and 40: one tile that cannot be quantized (the range is
// 1e30), so it is stored losslessly in GZIP_COMPRESSED_DATA.
func floatNullCase() []float32 {
	v := make([]float32, 64)
	for i := range v {
		v[i] = 100 + float32(i%7)
	}
	v[3] = float32(math.NaN())
	v[10] = crefgen.FloatNullValue
	v[20], v[40] = -1e30, -1e30
	return v
}

// TestFloatNullValuePixelStaysValue pins the documented difference from
// cfitsio: a real pixel equal to FLOATNULLVALUE in a tile with nulls is
// kept by this writer (cfitsio's own writer nulls it: golden
// rice_f32_floatnullvalue_t40x5, knownFloatNullValue). The pixel sits in
// a GZIP_COMPRESSED_DATA (lossless) tile here, so our reader, libcfitsio
// (fits_read_img with a NaN null value: cref "decode") and astropy all
// return it exactly; in quantized tiles it is a value too, within the
// quantization step (TestCompressedWriterFloatNullValueQuantizedMatchesCfitsio), and
// the real -1e30 pixels (the writer's default NullValue, Nullcheck off)
// too. funpack is the exception: it writes its output with
// FLOATNULLVALUE as the null value, so the pixel comes back NaN there;
// if funpack changes, this test says so.
func TestFloatNullValuePixelStaysValue(t *testing.T) {
	in := floatNullCase()
	path := filepath.Join(t.TempDir(), "floatnull.fits")
	if err := writeOne(path, func(f *File) error {
		_, err := AppendCompressedFloat32Image(f, nil, []int64{64}, in, CompressFloatOptions{})
		return err
	}); err != nil {
		t.Fatal(err)
	}
	check := func(who string, px []*float64) {
		t.Helper()
		if len(px) != len(in) {
			t.Fatalf("%s: %d pixels, want %d", who, len(px), len(in))
		}
		for i, v := range in {
			p := px[i]
			if i == 3 {
				if p != nil && !math.IsNaN(*p) {
					t.Errorf("%s: pixel 3 (NaN input) = %v, want null", who, *p)
				}
				continue
			}
			if p == nil || float32(*p) != v {
				t.Errorf("%s: pixel %d = %v, input %v (lossless tile: must be exact)", who, i, fmtp(p), v)
			}
		}
	}
	f, h := openCompressed(t, path)
	got, err := ReadPixelsCompressed[float32](h)
	f.Close()
	if err != nil {
		t.Fatal(err)
	}
	ours := make([]*float64, len(got))
	for i, v := range got {
		x := float64(v)
		ours[i] = &x
	}
	check("our reader", ours)

	check("libcfitsio", cfitsioPixels(t, path))
	a := astropyDecode(t, []interopCase{{path: path}})[path]
	if a.Error != "" {
		t.Fatalf("astropy: %s", a.Error)
	}
	check("astropy", a.Pixels)

	out := path + ".funpacked"
	funpackTo(t, out, path)
	fp := cfitsioPixels(t, out)
	if fp[10] != nil && !math.IsNaN(*fp[10]) {
		t.Fatalf("funpack now keeps the FLOATNULLVALUE pixel (%v); update the funpack note in CompressFloatOptions and docs/VALIDATION.md", *fp[10])
	}
	for _, i := range []int{20, 40} {
		if fp[i] == nil || float32(*fp[i]) != in[i] {
			t.Errorf("funpack: pixel %d = %v, input %v", i, fmtp(fp[i]), in[i])
		}
	}
}

// TestCfitsioImgnullStaleFallbackFlags documents why the interop tests
// read cfitsio's decode with fits_read_img and a null value, not with
// fits_read_imgnull: on cfitsio's own file golden
// rice_f32_nan_realnull_rows (whole-row tiles, most stored in
// GZIP_COMPRESSED_DATA with a NaN) libcfitsio 4.7.0's fits_read_imgnull
// keeps the null flags of earlier fallback tiles, so real pixels are
// flagged null, while fits_read_img (the golden's "dec" digest), our
// reader (TestReaderDecodesCfitsioWrittenGoldens) and astropy agree that
// only the NaN inputs are null. If a future cfitsio fixes this, the test
// says so.
func TestCfitsioImgnullStaleFallbackFlags(t *testing.T) {
	g := imgGoldenByName(t, "rice_f32_nan_realnull_rows")
	in := g.floatInput(t)
	for i, v := range in {
		if isNaN := math.IsNaN(float64(v)); isNaN != g.Dec.IsNull(i) {
			t.Fatalf("pixel %d: input %v, fits_read_img null %v: the golden no longer shows nulls exactly at the NaNs", i, v, g.Dec.IsNull(i))
		}
	}
	src := goldenFITS(t, filepath.Join(testutil.RepoRoot(t), imgGoldenDir, g.Name+".fits"))
	flagged := cfitsioDecode(t, "decodenull", src)
	stale := 0
	for i, p := range flagged {
		if p == nil && !math.IsNaN(float64(in[i])) {
			stale++
		}
	}
	if stale == 0 {
		t.Fatalf("fits_read_imgnull now flags only the NaN pixels of %s; cref decode may use it again", g.Name)
	}
	t.Logf("fits_read_imgnull flags %d real pixels of cfitsio's own file as null (fits_read_img: %d nulls, the NaN inputs)", stale, g.Dec.NullCount())
}

// TestCompressedWriterBlockSize16MatchesFitsRcomp: CompressOptions.BlockSize
// 16 is really used and declared. cfitsio's fits_img_compress always
// writes 32-pixel blocks, so the truth is its codec: a one-tile 1-D image
// of each pixel width must store exactly the stream fits_rcomp* produced
// with blocksize 16 (golden rice/rice_{b,s,i}_bs16), declared as ZNAME1 =
// 'BLOCKSIZE', ZVAL1 = 16 (the keyword every reader decodes the stream
// with), and our reader (cfitsio-proven) must give the input back.
func TestCompressedWriterBlockSize16MatchesFitsRcomp(t *testing.T) {
	for _, tag := range []string{"b", "s", "i"} {
		name := "rice_" + tag + "_bs16"
		t.Run(name, func(t *testing.T) {
			dir := filepath.Join("compress", "testdata", "cref", "golden", "rice")
			b, err := os.ReadFile(filepath.Join(dir, name+".json"))
			if err != nil {
				t.Fatalf("committed golden missing (failure, not a skip): %v", err)
			}
			var meta struct {
				BytePix   int    `json:"bytepix"`
				N         int    `json:"n"`
				BlockSize int    `json:"blocksize"`
				Kind      int    `json:"kind"`
				GenSeed   uint64 `json:"gen_seed"`
				InSHA256  string `json:"in_sha256"`
			}
			if err := json.Unmarshal(b, &meta); err != nil {
				t.Fatal(err)
			}
			want, err := os.ReadFile(filepath.Join(dir, name+".rice"))
			if err != nil {
				t.Fatalf("committed golden missing (failure, not a skip): %v", err)
			}
			if meta.BlockSize != 16 {
				t.Fatalf("golden blocksize %d", meta.BlockSize)
			}
			vals := crefgen.IntKind(int64(meta.N), meta.BytePix, meta.Kind, meta.GenSeed)
			if err := crefgen.Check(name, crefgen.LEInts(vals, meta.BytePix), meta.InSHA256); err != nil {
				t.Fatal(err)
			}
			path := filepath.Join(t.TempDir(), name+".fits")
			opts := CompressOptions{Algorithm: compress.RICE1, BlockSize: 16}
			shape := []int64{int64(meta.N)}
			if err := writeOne(path, func(f *File) error {
				switch meta.BytePix {
				case 1:
					return appendConv(f, shape, vals, opts, func(x int64) uint8 { return uint8(x) })
				case 2:
					return appendConv(f, shape, vals, opts, func(x int64) int16 { return int16(x) })
				}
				return appendConv(f, shape, vals, opts, func(x int64) int32 { return int32(x) })
			}); err != nil {
				t.Fatal(err)
			}
			f, h := openCompressed(t, path)
			defer f.Close()
			hd := h.Header()
			if n1, _ := hd.String("ZNAME1"); strings.TrimSpace(n1) != "BLOCKSIZE" {
				t.Errorf("ZNAME1 = %q, want BLOCKSIZE", n1)
			}
			if v, err := hd.Int("ZVAL1"); err != nil || v != 16 {
				t.Errorf("ZVAL1 = %v (%v), want 16: readers would decode the 16-pixel blocks as 32", v, err)
			}
			if v, err := hd.Int("ZVAL2"); err != nil || v != int64(meta.BytePix) {
				t.Errorf("ZVAL2 = %v (%v), want %d", v, err, meta.BytePix)
			}
			cells := readTableCols(t, h).data["COMPRESSED_DATA"]
			if len(cells) != 1 || !bytes.Equal(cells[0], want) {
				got := []byte(nil)
				if len(cells) > 0 {
					got = cells[0]
				}
				t.Errorf("COMPRESSED_DATA differs from fits_rcomp blocksize 16: %s", diffBytes(got, want))
			}
			px, err := ReadPixelsCompressed[int64](h)
			if err != nil {
				t.Fatalf("read back: %v", err)
			}
			for i := range vals {
				if px[i] != vals[i] {
					t.Fatalf("pixel %d read back %d, input %d", i, px[i], vals[i])
				}
			}
		})
	}
}

// TestCompressedWriterNoDuplicateKeys: a keyword appears once in a header
// (fitsverify reports "Keyword ... is duplicated" as an error; checked on
// the same header by TestCompressedWriterInterop i16_rice_caller_zkeys).
// A caller header with EXTEND and ZEXTEND, CHECKSUM and ZHECKSUM, BLOCKED
// and ZBLOCKED, and stale compression keywords must give one card per
// keyword, ours winning for the compression keywords.
func TestCompressedWriterNoDuplicateKeys(t *testing.T) {
	hdr := primaryStyleHeader(16, []int64{8, 8})
	hdr.Set("ZEXTEND", true, "already renamed")
	hdr.Set("ZHECKSUM", "BBBBBBBBBBBBBBBB", "already renamed")
	hdr.Set("ZBLOCKED", true, "already renamed")
	hdr.Set("ZBITPIX", int64(99), "stale")
	hdr.Set("ZCMPTYPE", "GZIP_1", "stale")
	hdr.Set("ZTILE1", int64(7), "stale")
	hdr.Set("ZBLANK", int64(-1), "stale")
	path := filepath.Join(t.TempDir(), "dup.fits")
	px := make([]int16, 64)
	if err := writeOne(path, func(f *File) error {
		_, err := AppendCompressedImage(f, hdr, []int64{8, 8}, px, CompressOptions{})
		return err
	}); err != nil {
		t.Fatal(err)
	}
	f, h := openCompressed(t, path)
	defer f.Close()
	seen := map[string]int{}
	for _, c := range h.Header().Cards() {
		if c.Key == "COMMENT" || c.Key == "HISTORY" || c.Key == "" {
			continue
		}
		seen[c.Key]++
	}
	for k, n := range seen {
		if n > 1 {
			t.Errorf("keyword %s appears %d times", k, n)
		}
	}
	hd := h.Header()
	if v, _ := hd.Int("ZBITPIX"); v != 16 {
		t.Errorf("ZBITPIX = %d, want 16 (ours, not the caller's stale 99)", v)
	}
	if v, _ := hd.Int("ZTILE1"); v != 8 {
		t.Errorf("ZTILE1 = %d, want 8 (ours, not the caller's stale 7)", v)
	}
	if v, _ := hd.String("ZCMPTYPE"); strings.TrimSpace(v) != "RICE_1" {
		t.Errorf("ZCMPTYPE = %q, want RICE_1 (ours, not the caller's stale GZIP_1)", v)
	}
	if hd.Has("ZBLANK") {
		t.Errorf("stale caller ZBLANK copied: readers would null pixel value -1")
	}
}
