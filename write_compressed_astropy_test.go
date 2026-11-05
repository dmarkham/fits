package fits

import (
	"encoding/binary"
	"encoding/json"
	"fmt"
	"math"
	"path/filepath"
	"regexp"
	"testing"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// TestCompressedWriterInterop is the end-to-end interop matrix for the
// compressed-image writers: every algorithm x pixel type x dither method
// x tiling the writers accept is written, then judged only by external
// tools:
//
//  1. fitsverify (HEASARC, cfitsio 4.7.0): no warnings, no errors;
//  2. libcfitsio, through the testdata/cfref harness: the file opens and
//     decodes;
//  3. funpack: decompresses it; the funpacked image (dumped again by
//     cfitsio) has exactly the pixels cfitsio decoded directly;
//  4. astropy.io.fits (its own tile-compression implementation): decodes
//     exactly the same pixels.
//
// Integer images that are lossless (all but HCOMPRESS_1 with scale > 1)
// must decode to the input exactly. Float32 images must decode to the
// same float32 bits in all three readers, NaN input pixels must come back
// as nulls/NaN, tiles stored through the GZIP_COMPRESSED_DATA fallback
// must be exact, and every other pixel must be within half a
// quantization step (0.5*ZSCALE of its tile, plus float32 rounding) of
// the input -- the bound fits_quantize_float guarantees for both
// NO_DITHER and subtractive dithering. A writer that dithered with the
// wrong seed, or claimed SUBTRACTIVE_DITHER_1 while not dithering,
// breaks that bound. Under SUBTRACTIVE_DITHER_2 input zeros must decode
// to exactly 0. With Nullcheck, pixels equal to NullValue must come back
// null. The matrix also covers 1-D images, HCOMPRESS_1/PLIO_1 cubes and
// a caller header copied from a primary HDU (SIMPLE/EXTEND/BLOCKED/
// CHECKSUM/DATASUM must not reach the extension as such: fitsverify).
//
// Bit-for-bit parity of the same writer with cfitsio's own writer is
// TestCompressedWriterMatchesCfitsio; this test proves the files are
// read correctly by everyone.
//
// Tier 1 replays the tools' recorded readings of these exact bytes
// (testutil.Run); tier 2 runs fitsverify and funpack (PATH), the built
// cref harness and FITS_ASTROPY_PYTHON live.
func TestCompressedWriterInterop(t *testing.T) {

	dir := t.TempDir()
	cases := interopCases()
	for i := range cases {
		cases[i].path = filepath.Join(dir, fmt.Sprintf("%02d_%s.fits", i, cases[i].name))
		if err := cases[i].write(cases[i].path); err != nil {
			t.Fatalf("%s: write: %v", cases[i].name, err)
		}
	}
	// Candidates for the recorded readings (cfitsioPixelsAs): this
	// package's own reading, which the readers must equal exactly.
	ours := map[string][]*float64{}
	for i := range cases {
		ours[cases[i].path] = ourPixels(t, cases[i].path, cases[i].float)
	}
	astro := astropyDecodeAs(t, cases, ours)

	for i := range cases {
		c := &cases[i]
		t.Run(c.name, func(t *testing.T) {
			fv := testutil.Fitsverify(t, c.path)
			if !fv.OK() {
				t.Errorf("fitsverify: %d warning(s), %d error(s):\n%s", fv.Warnings, fv.Errors, fv.Report)
			}
			direct := cfitsioPixelsAs(t, c.path, ours[c.path], c.float)

			out := c.path + ".funpacked"
			funpackTo(t, out, c.path)
			unpacked := cfitsioPixelsAs(t, out, direct, c.float)

			a := astro[c.path]
			if a.Error != "" {
				t.Fatalf("astropy failed: %s", a.Error)
			}
			if len(direct) != c.npix || len(unpacked) != c.npix || len(a.Pixels) != c.npix {
				t.Fatalf("pixel counts: cfitsio %d, funpack %d, astropy %d, want %d", len(direct), len(unpacked), len(a.Pixels), c.npix)
			}
			bad := 0
			report := func(format string, args ...any) {
				if bad < 5 {
					t.Errorf(format, args...)
				}
				bad++
			}
			// funpack 1.7.0 (CFITSIO 4.7.0) mis-decodes float images that have null pixels
			// and tiles narrower than NAXIS1 -- also cfitsio's own files
			// (TestFunpackMisdecodesNullTiledFloats). For those its pixels
			// are not compared; it must still succeed.
			funpackTrusted := !((c.float || c.double) && c.hasNaN() && c.tile[0] < c.shape[0])
			if !funpackTrusted {
				t.Logf("funpack pixels not compared: known funpack 1.7.0 (CFITSIO 4.7.0) bug for null + narrow-tile float images")
			}
			for k := 0; k < c.npix; k++ {
				d, u, ap := direct[k], unpacked[k], a.Pixels[k]
				if !funpackTrusted {
					u = d
				}
				if c.double {
					if !sameDouble(d, u) || !sameDouble(d, ap) {
						report("pixel %d: cfitsio %v, funpack %v, astropy %v disagree", k, fmtp(d), fmtp(u), fmtp(ap))
						continue
					}
					c.checkDouble(k, d, a.ZScale, report)
					continue
				}
				if !samePixel(d, u, c.float) || !samePixel(d, ap, c.float) {
					report("pixel %d: cfitsio %v, funpack %v, astropy %v disagree", k, fmtp(d), fmtp(u), fmtp(ap))
					continue
				}
				if !c.float {
					if c.lossless && (d == nil || *d != c.ints[k]) {
						report("pixel %d: decoded %v, input %v (lossless)", k, fmtp(d), c.ints[k])
					}
					continue
				}
				in, null := cfitsioFloatInput(c.floats[k])
				if null || (c.nullcheck && in == c.nullValue) {
					if d != nil && !math.IsNaN(*d) {
						report("pixel %d: input %v (a null) decoded as %v", k, c.floats[k], *d)
					}
					continue
				}
				if d == nil || math.IsNaN(*d) {
					report("pixel %d: input %v decoded as null", k, in)
					continue
				}
				got := float32(*d)
				if c.wantQuantiz == "SUBTRACTIVE_DITHER_2" && in == 0 && got != 0 {
					// SUBTRACTIVE_DITHER_2 stores exact zeros as ZERO_VALUE;
					// a writer quantizing them like other pixels breaks this.
					report("pixel %d: input 0 decoded as %v (SUBTRACTIVE_DITHER_2 keeps zeros exact)", k, got)
					continue
				}
				zs := a.ZScale[c.tileOf(k)]
				if zs < 0 {
					report("pixel %d: tile %d has a non-finite ZSCALE", k, c.tileOf(k)+1)
					continue
				}
				if zs == 0 { // GZIP_COMPRESSED_DATA fallback tile: lossless
					if got != in {
						report("pixel %d: fallback tile decoded %v, input %v (must be exact)", k, got, in)
					}
					continue
				}
				bound := 0.5*zs*(1+1e-12) + float64(ulp32(in))
				if diff := math.Abs(float64(got) - float64(in)); diff > bound {
					report("pixel %d: |decoded %v - input %v| = %g > half a quantization step %g", k, got, in, diff, bound)
				}
			}
			if bad > 0 {
				t.Errorf("%d/%d pixels wrong", bad, c.npix)
			}
			if (c.float || c.double) && a.ZQuantiz != c.wantQuantiz {
				t.Errorf("ZQUANTIZ %q, want %q", a.ZQuantiz, c.wantQuantiz)
			}
		})
	}
}

// cfitsioFloatInput is what cfitsio's writer sees for a float pixel
// (fffr4r4's fnan test when fits_img_compress reads a tile): an all-ones
// exponent (NaN, +-Inf) is a null, an all-zero exponent (+-0, denormals)
// is 0. Golden-proven in TestCompressedWriterMatchesCfitsio
// (rice_f32_inf_*, gzip2_f32_inf_fallback_t40x5).
func cfitsioFloatInput(v float32) (float32, bool) {
	switch math.Float32bits(v) & 0x7F800000 {
	case 0x7F800000:
		return 0, true
	case 0:
		return 0, false
	}
	return v, false
}

func fmtp(p *float64) string {
	if p == nil {
		return "null"
	}
	return fmt.Sprint(*p)
}

func ulp32(v float32) float32 {
	a := float32(math.Abs(float64(v)))
	return math.Nextafter32(a, float32(math.Inf(1))) - a
}

// sameValue compares two decoded pixels; nil is a null. Float pixels are
// compared as float32 bits (NaN equals NaN).
func samePixel(a, b *float64, isFloat bool) bool {
	if a == nil || b == nil {
		return a == nil && b == nil
	}
	if isFloat {
		x, y := float32(*a), float32(*b)
		if x != x || y != y {
			return x != x && y != y
		}
		return math.Float32bits(x) == math.Float32bits(y)
	}
	return *a == *b
}

// toolCref is the compress/testdata/cref harness (`make -C
// compress/testdata/cref` builds it for tier 2; tier 1 replays its
// recorded readings, testutil.Run).
var toolCref = testutil.Tool{Name: "cref", Harness: "compress/testdata/cref/cref", Sources: []string{"compress/testdata/cref/cref.c"}}

// funpackTo is funpack -O out src. The decompressed file is only read
// back by the cfitsio harness, so only its digest is recorded.
func funpackTo(t *testing.T, out, src string) {
	t.Helper()
	testutil.RunOK(t, testutil.ToolFunpack, testutil.Call{Args: []string{"-O", out, src}, Inputs: []string{src}, Outputs: []string{out},
		OutputsByDigest: true, FixOutput: unstampChecksums})
}

// checksumStamp is the wall-clock time cfitsio writes into the comment of
// the CHECKSUM and DATASUM cards ("HDU checksum updated <time>").
var checksumStamp = regexp.MustCompile(`checksum updated \d{4}-\d\d-\d\dT\d\d:\d\d:\d\d`)

// unstampChecksums makes a funpacked file the same on every run: the
// time in the CHECKSUM/DATASUM comments becomes 2000-01-01T00:00:00 and
// the CHECKSUM value, which covers that comment, sixteen zeros. Only the
// pixels of the file are read back (cfitsio's decode checks no checksum).
func unstampChecksums(b []byte) []byte {
	out := append([]byte(nil), b...)
	for off := 0; off+80 <= len(out); off += 80 {
		card := out[off : off+80]
		if loc := checksumStamp.FindIndex(card); loc != nil {
			copy(card[loc[0]:loc[1]], "checksum updated 2000-01-01T00:00:00")
			if string(card[:10]) == "CHECKSUM= " && card[10] == '\'' && card[27] == '\'' {
				copy(card[11:27], "0000000000000000")
			}
		}
	}
	return out
}

// cfitsioPixels decodes HDU 1 of path (compressed or not) with libcfitsio
// (`cref decode`: fits_read_img TDOUBLE with a NaN null value, scaling
// applied), nil = null. fits_read_imgnull is not used: libcfitsio 4.7.0
// flags real pixels null in float images with several null-holding
// GZIP_COMPRESSED_DATA tiles (TestCfitsioImgnullStaleFallbackFlags).
func cfitsioPixels(t *testing.T, path string) []*float64 {
	t.Helper()
	return cfitsioDecode(t, "decode", path)
}

// cfitsioDecode runs `cref <mode> path` (decode or decodenull).
func cfitsioDecode(t *testing.T, mode, path string) []*float64 {
	t.Helper()
	return cfitsioDecodeAs(t, mode, path, nil, false)
}

// cfitsioPixelsAs is cfitsioPixels for a reading expected to equal cand
// (nil: no expectation). When cfitsio's live reading is exactly cand, only
// its digest is recorded (testutil.FloatsDigest; the full reading of the
// interop matrix would not fit the repository's size budget), and on
// replay cand must match that digest exactly; otherwise the reading is
// recorded value by value. Either way the result is cfitsio's reading.
//
// f32: the image is float32 and the caller judges each pixel as
// float32(value) only (samePixel, the quantization bound), so the digest
// is over the values rounded to float32, and so is the result.
func cfitsioPixelsAs(t *testing.T, path string, cand []*float64, f32 bool) []*float64 {
	t.Helper()
	return cfitsioDecodeAs(t, "decode", path, cand, f32)
}

// roundPixF32 is v with every value rounded to float32 (nil kept).
func roundPixF32(v []*float64) []*float64 {
	out := make([]*float64, len(v))
	for i, p := range v {
		if p != nil {
			x := float64(float32(*p))
			out[i] = &x
		}
	}
	return out
}

// crefDecoded is the output of `cref decode`, or its digested form.
type crefDecoded struct {
	Cfitsio string                 `json:"cfitsio"`
	Naxes   []int64                `json:"naxes"`
	Pixels  []*float64             `json:"pixels,omitempty"`
	Digest  *testutil.FloatsDigest `json:"digest,omitempty"`
}

func parseCrefDecoded(b []byte) (crefDecoded, error) {
	var r crefDecoded
	if err := json.Unmarshal(b, &r); err != nil {
		return r, err
	}
	for _, p := range r.Pixels {
		if p != nil && math.Abs(*p) == math.MaxFloat64 { // +-Inf, see cref.c run_decode
			*p = math.Inf(int(math.Copysign(1, *p)))
		}
	}
	return r, nil
}

func cfitsioDecodeAs(t *testing.T, mode, path string, cand []*float64, f32 bool) []*float64 {
	t.Helper()
	call := testutil.Call{Args: []string{mode, path}, Inputs: []string{path}}
	if cand != nil {
		call.Compact = func(b []byte) []byte {
			r, err := parseCrefDecoded(b)
			if err != nil || len(r.Pixels) != len(cand) {
				return b
			}
			if f32 {
				r.Pixels = roundPixF32(r.Pixels)
			}
			d := testutil.DigestFloats(r.Pixels)
			if ok, _ := d.Matches(cand); !ok {
				return b
			}
			out, _ := json.Marshal(crefDecoded{Cfitsio: r.Cfitsio, Naxes: r.Naxes, Digest: &d})
			return out
		}
	}
	res := testutil.Run(t, toolCref, call)
	if res.Exit != 0 {
		t.Fatalf("cfitsio cannot decode %s: exit status %d: %s", path, res.Exit, res.Stderr)
	}
	r, err := parseCrefDecoded(res.Stdout)
	if err != nil {
		t.Fatalf("parse cref decode output: %v", err)
	}
	if r.Cfitsio != "4.070" {
		t.Fatalf("cref linked against cfitsio %s, goldens and expectations are for 4.070 (4.7.0)", r.Cfitsio)
	}
	if r.Digest != nil {
		if ok, why := r.Digest.Matches(cand); !ok {
			t.Fatalf("cfitsio's recorded reading of %s is not the expected pixels: %s", path, why)
		}
		return cand
	}
	if f32 {
		return roundPixF32(r.Pixels)
	}
	return r.Pixels
}

// ourPixels is this package's reading of HDU 1 of a compressed file, in
// the form of cfitsioPixels (nil = null; a float32 image read as float32),
// or nil when it cannot read it. It is only a candidate for
// cfitsioPixelsAs: the truth stays cfitsio's.
func ourPixels(t *testing.T, path string, f32 bool) []*float64 {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		return nil
	}
	defer f.Close()
	hdu, err := f.HDU(1)
	if err != nil {
		return nil
	}
	h, ok := hdu.(*CompressedImageHDU)
	if !ok {
		return nil
	}
	var v []float64
	var mask []bool
	if f32 {
		v32, m, err := ReadPixelsCompressedMasked[float32](h)
		if err != nil {
			return nil
		}
		for _, x := range v32 {
			v = append(v, float64(x))
		}
		mask = m
	} else if v, mask, err = ReadPixelsCompressedMasked[float64](h); err != nil {
		return nil
	}
	out := make([]*float64, len(v))
	for i := range v {
		// cref decode (fits_read_img with a NaN null value) and the
		// astropy script report every NaN as null.
		if mask[i] && !math.IsNaN(v[i]) {
			x := v[i]
			out[i] = &x
		}
	}
	return out
}

type astropyResult struct {
	Pixels   []*float64             `json:"pixels,omitempty"`
	Digest   *testutil.FloatsDigest `json:"digest,omitempty"`
	ZScale   []float64              `json:"zscale"`
	ZQuantiz string                 `json:"zquantiz"`
	Error    string                 `json:"error"`
}

// astropyDecode reads every case file with astropy in one interpreter.
// Any astropy warning about a file is turned into an error.
func astropyDecode(t *testing.T, cases []interopCase) map[string]astropyResult {
	t.Helper()
	return astropyDecodeAs(t, cases, nil)
}

// astropyDecodeAs is astropyDecode where astropy's pixels for a file are
// expected to equal cand[path]: then only their digest is recorded, as in
// cfitsioPixelsAs.
//
// For a float32 case (interopCase.float) the caller judges astropy's
// pixels as float32 only, so they are digested and returned rounded to
// float32, as in cfitsioPixelsAs.
func astropyDecodeAs(t *testing.T, cases []interopCase, cand map[string][]*float64) map[string]astropyResult {
	t.Helper()
	const script = `
import sys, json, math, warnings
import numpy as np
from astropy.io import fits
warnings.simplefilter("error")
out = {}
for p in sys.argv[1:]:
    try:
        with fits.open(p) as h:
            d = np.asarray(h[1].data, dtype=np.float64).ravel(order="C")
            # JSON has no Inf: +-Inf travel as +-DBL_MAX (turned back into
            # +-Inf on the Go side) so a reader returning Inf is visible
            vals = [None if math.isnan(v) else (math.copysign(sys.float_info.max, v) if math.isinf(v) else v) for v in d.tolist()]
        with fits.open(p, disable_image_compression=True) as h:
            tb = h[1].data
            # a non-finite ZSCALE is reported as -1 (JSON has no NaN)
            zs = [z if math.isfinite(z) else -1.0 for z in tb["ZSCALE"].tolist()] if "ZSCALE" in tb.names else []
            zq = h[1].header.get("ZQUANTIZ", "")
        out[p] = {"pixels": vals, "zscale": zs, "zquantiz": zq, "error": ""}
    except Exception as e:
        out[p] = {"pixels": [], "zscale": [], "zquantiz": "", "error": repr(e)}
json.dump(out, sys.stdout, allow_nan=False)
`
	var paths []string
	f32 := map[string]bool{}
	for _, c := range cases {
		paths = append(paths, c.path)
		f32[c.path] = c.float && cand != nil
	}
	parse := func(b []byte) (map[string]astropyResult, error) {
		res := map[string]astropyResult{}
		if err := json.Unmarshal(b, &res); err != nil {
			return nil, err
		}
		for p, r := range res {
			for _, v := range r.Pixels {
				if v != nil && math.Abs(*v) == math.MaxFloat64 {
					*v = math.Inf(int(math.Copysign(1, *v)))
				}
			}
			if f32[p] {
				r.Pixels = roundPixF32(r.Pixels)
				res[p] = r
			}
		}
		return res, nil
	}
	call := testutil.Call{Args: append([]string{"-c", script}, paths...), Inputs: paths}
	if cand != nil {
		call.Compact = func(b []byte) []byte {
			res, err := parse(b)
			if err != nil {
				return b
			}
			for p, r := range res {
				if c, ok := cand[p]; ok && c != nil && r.Error == "" {
					d := testutil.DigestFloats(r.Pixels)
					if ok, _ := d.Matches(c); ok {
						r.Pixels, r.Digest = nil, &d
						res[p] = r
					}
				}
			}
			out, _ := json.Marshal(res)
			return out
		}
	}
	b := testutil.RunOK(t, testutil.ToolAstropy, call).Stdout
	res, err := parse(b)
	if err != nil {
		t.Fatalf("parse astropy output: %v", err)
	}
	for p, r := range res {
		if r.Digest != nil {
			if ok, why := r.Digest.Matches(cand[p]); !ok {
				t.Fatalf("astropy's recorded reading of %s is not the expected pixels: %s", p, why)
			}
			r.Pixels, r.Digest = cand[p], nil
			res[p] = r
		}
	}
	return res
}

// interopCase is one written file of the matrix.
type interopCase struct {
	name        string
	path        string
	shape       []int64
	tile        []int64 // resolved tile shape (maps pixels to ZSCALE rows)
	npix        int
	float       bool
	lossless    bool
	ints        []float64 // integer input as physical values
	floats      []float32
	wantQuantiz string
	nullcheck   bool // float: pixels equal to nullValue must decode as null
	nullValue   float32
	double      bool // float64 image (AppendCompressedFloat64Image)
	doubles     []float64
	nullValue64 float64
	write       func(path string) error
}

func (c *interopCase) hasNaN() bool {
	for _, v := range c.floats {
		if x, null := cfitsioFloatInput(v); null || (c.nullcheck && x == c.nullValue) {
			return true
		}
	}
	for _, v := range c.doubles {
		if x, null := cfitsioDoubleInput(v); null || (c.nullcheck && x == c.nullValue64) {
			return true
		}
	}
	return false
}

// cfitsioDoubleInput is cfitsioFloatInput for float64 (dnan): an
// all-ones exponent is a null, an all-zero exponent is 0 (golden
// img64/rice_f64_inf_sd1_z1, gzip2_f64_inf_fallback_t40x5).
func cfitsioDoubleInput(v float64) (float64, bool) {
	switch math.Float64bits(v) & 0x7FF0000000000000 {
	case 0x7FF0000000000000:
		return 0, true
	case 0:
		return 0, false
	}
	return v, false
}

// sameDouble compares two decoded float64 pixels bit for bit (nil is a
// null; NaN equals NaN).
func sameDouble(a, b *float64) bool {
	if a == nil || b == nil {
		return a == nil && b == nil
	}
	if *a != *a || *b != *b {
		return *a != *a && *b != *b
	}
	return math.Float64bits(*a) == math.Float64bits(*b)
}

// checkDouble judges decoded pixel k of a float64 case against its input:
// nulls null, fallback tiles exact, SUBTRACTIVE_DITHER_2 zeros exact, the
// rest within half a quantization step (plus one float64 ulp).
func (c *interopCase) checkDouble(k int, d *float64, zscale []float64, report func(string, ...any)) {
	in, null := cfitsioDoubleInput(c.doubles[k])
	if null || (c.nullcheck && in == c.nullValue64) {
		if d != nil && !math.IsNaN(*d) {
			report("pixel %d: input %v (a null) decoded as %v", k, c.doubles[k], *d)
		}
		return
	}
	if d == nil || math.IsNaN(*d) {
		report("pixel %d: input %v decoded as null", k, in)
		return
	}
	got := *d
	if c.wantQuantiz == "SUBTRACTIVE_DITHER_2" && in == 0 && got != 0 {
		report("pixel %d: input 0 decoded as %v (SUBTRACTIVE_DITHER_2 keeps zeros exact)", k, got)
		return
	}
	zs := zscale[c.tileOf(k)]
	if zs < 0 {
		report("pixel %d: tile %d has a non-finite ZSCALE", k, c.tileOf(k)+1)
		return
	}
	if zs == 0 {
		if got != in {
			report("pixel %d: fallback tile decoded %v, input %v (must be exact)", k, got, in)
		}
		return
	}
	ulp := math.Nextafter(math.Abs(in), math.Inf(1)) - math.Abs(in)
	bound := 0.5*zs*(1+1e-12) + ulp
	if diff := math.Abs(got - in); diff > bound {
		report("pixel %d: |decoded %v - input %v| = %g > half a quantization step %g", k, got, in, diff, bound)
	}
}

// tileOf returns the 0-based tile (table row) holding pixel k.
func (c *interopCase) tileOf(k int) int {
	idx, mul, rem := 0, 1, k
	for d := range c.shape {
		n, tl := int(c.shape[d]), int(c.tile[d])
		coord := rem % n
		rem /= n
		idx += (coord / tl) * mul
		mul *= (n + tl - 1) / tl
	}
	return idx
}

func writeOne(path string, fn func(f *File) error) error {
	f, err := Create(path)
	if err != nil {
		return err
	}
	if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
		f.Close()
		return err
	}
	if err := fn(f); err != nil {
		f.Close()
		return err
	}
	return f.Close()
}

// lcg is a small deterministic generator for test inputs (the inputs are
// not truth; the external decoders are).
type lcg uint64

func (r *lcg) next() uint64 {
	*r = *r*6364136223846793005 + 1442695040888963407
	return uint64(*r) >> 11
}
func (r *lcg) unif() float64 { return float64(r.next()) / float64(1<<53) }
func (r *lcg) gauss() float64 {
	u1, u2 := r.unif(), r.unif()
	if u1 < 1e-12 {
		u1 = 1e-12
	}
	return math.Sqrt(-2*math.Log(u1)) * math.Cos(2*math.Pi*u2)
}

// intPattern returns adversarial integer pixels in [lo, hi]: bands of
// smooth ramp, full-range random, the range extremes cycled (every
// difference overflows the pixel width) and constant runs. mask gives
// PLIO-style runs in [0, min(hi, 2^24-1)] instead.
func intPattern(n, nx int, lo, hi int64, mask bool, seed uint64) []int64 {
	r := lcg(seed)
	out := make([]int64, n)
	mid := lo + (hi-lo)/2
	for i := range out {
		row, col := i/nx, i%nx
		if mask {
			top := hi
			if top > 16777215 {
				top = 16777215
			}
			vals := []int64{0, 1, 2, top, top/2 + 1, 0}
			if (row/3)%2 == 0 && col > 3 {
				out[i] = vals[(col/4)%len(vals)]
			}
			continue
		}
		switch (row / 4) % 4 {
		case 0:
			out[i] = min(lo+(hi-lo)/3+int64(col*3+row), hi)
		case 1:
			out[i] = lo + int64(r.next()%uint64(hi-lo+1))
		case 2:
			out[i] = []int64{lo, hi, mid, lo + 1}[i%4]
		default:
			out[i] = mid
		}
	}
	return out
}

// noisePattern is a smooth noisy sky well inside [lo, hi].
func noisePattern(n, nx int, lo, hi int64, seed uint64) []int64 {
	r := lcg(seed)
	span := float64(hi - lo)
	mid := float64(lo) + span/2
	sig := span / 200
	out := make([]int64, n)
	for i := range out {
		v := mid + sig*r.gauss() + float64(i%nx)*sig/10
		out[i] = int64(math.Round(math.Max(float64(lo), math.Min(float64(hi), v))))
	}
	return out
}

func appendConv[T Numeric](f *File, shape []int64, v []int64, o CompressOptions, conv func(int64) T) error {
	data := make([]T, len(v))
	for i, x := range v {
		data[i] = conv(x)
	}
	_, err := AppendCompressedImage(f, nil, shape, data, o)
	return err
}

func interopCases() []interopCase {
	var cases []interopCase
	type intType struct {
		name   string
		lo, hi int64
		write  func(f *File, shape []int64, v []int64, o CompressOptions) error
	}
	types := []intType{
		{"u8", 0, 255, func(f *File, s []int64, v []int64, o CompressOptions) error {
			return appendConv(f, s, v, o, func(x int64) uint8 { return uint8(x) })
		}},
		{"s8", -128, 127, func(f *File, s []int64, v []int64, o CompressOptions) error {
			return appendConv(f, s, v, o, func(x int64) int8 { return int8(x) })
		}},
		{"i16", -32768, 32767, func(f *File, s []int64, v []int64, o CompressOptions) error {
			return appendConv(f, s, v, o, func(x int64) int16 { return int16(x) })
		}},
		{"u16", 0, 65535, func(f *File, s []int64, v []int64, o CompressOptions) error {
			return appendConv(f, s, v, o, func(x int64) uint16 { return uint16(x) })
		}},
		{"i32", math.MinInt32, math.MaxInt32, func(f *File, s []int64, v []int64, o CompressOptions) error {
			return appendConv(f, s, v, o, func(x int64) int32 { return int32(x) })
		}},
		{"u32", 0, math.MaxUint32, func(f *File, s []int64, v []int64, o CompressOptions) error {
			return appendConv(f, s, v, o, func(x int64) uint32 { return uint32(x) })
		}},
	}
	type intVariant struct {
		name  string
		shape []int64
		opts  CompressOptions
		tile  []int64
		mask  bool
		lossy bool
	}
	variants := []intVariant{
		{"rice", []int64{36, 22}, CompressOptions{Algorithm: compress.RICE1}, []int64{36, 1}, false, false},
		{"rice_bs16_t16x8", []int64{36, 22}, CompressOptions{Algorithm: compress.RICE1, BlockSize: 16, TileShape: []int64{16, 8}}, []int64{16, 8}, false, false},
		{"gzip1", []int64{36, 22}, CompressOptions{Algorithm: compress.GZIP1}, []int64{36, 1}, false, false},
		{"gzip2_t16x8", []int64{36, 22}, CompressOptions{Algorithm: compress.GZIP2, TileShape: []int64{16, 8}}, []int64{16, 8}, false, false},
		{"hcomp", []int64{36, 22}, CompressOptions{Algorithm: compress.HCOMPRESS1}, []int64{36, 22}, false, false},
		{"hcomp_t12x11", []int64{36, 22}, CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{12, 11}}, []int64{12, 11}, false, false},
		{"hcomp_s3", []int64{36, 22}, CompressOptions{Algorithm: compress.HCOMPRESS1, HCompressScale: 3}, []int64{36, 22}, false, true},
		{"hcomp_40rows", []int64{20, 40}, CompressOptions{Algorithm: compress.HCOMPRESS1}, []int64{20, 16}, false, false},
		{"plio", []int64{36, 22}, CompressOptions{Algorithm: compress.PLIO1}, []int64{36, 1}, true, false},
		{"rice_cube_t5x4x2", []int64{9, 7, 3}, CompressOptions{Algorithm: compress.RICE1, TileShape: []int64{5, 4, 2}}, []int64{5, 4, 2}, false, false},
		{"gzip2_cube", []int64{9, 7, 3}, CompressOptions{Algorithm: compress.GZIP2}, []int64{9, 1, 1}, false, false},
		{"rice_1d", []int64{300}, CompressOptions{Algorithm: compress.RICE1}, []int64{300}, false, false},
		{"gzip1_1d", []int64{300}, CompressOptions{Algorithm: compress.GZIP1}, []int64{300}, false, false},
		{"hcomp_cube", []int64{12, 9, 3}, CompressOptions{Algorithm: compress.HCOMPRESS1}, []int64{12, 9, 1}, false, false},
		{"plio_cube", []int64{12, 9, 3}, CompressOptions{Algorithm: compress.PLIO1}, []int64{12, 1, 1}, true, false},
		// tile longer than the image axis: ZTILE1 = 100 > NAXIS1 = 36,
		// as cfitsio writes it (golden rice_i16_11x8_t100x3)
		{"rice_t100x8", []int64{36, 22}, CompressOptions{Algorithm: compress.RICE1, TileShape: []int64{100, 8}}, []int64{100, 8}, false, false},
		// 4-pixel HCOMPRESS tiles and 4-pixel edge tiles (the smallest
		// cfitsio accepts; golden hcomp_i16_12x24_t4x12, hcomp_i16_12x9_t8x5)
		{"hcomp_t4x12", []int64{12, 24}, CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{4, 12}}, []int64{4, 12}, false, false},
		{"hcomp_edge4_t8x5", []int64{12, 9}, CompressOptions{Algorithm: compress.HCOMPRESS1, TileShape: []int64{8, 5}}, []int64{8, 5}, false, false},
	}
	for ti, ty := range types {
		for vi, v := range variants {
			if v.mask && (ty.name == "u32" || ty.name == "u16") {
				// Refused by the writer: uint32 is stored as v-2^31, out of
				// PLIO range (as in cfitsio, golden plio_u32); uint16 PLIO is
				// read differently by cfitsio and astropy
				// (TestPLIOUint16AstropyMisreadsCfitsio).
				continue
			}
			n := 1
			for _, d := range v.shape {
				n *= int(d)
			}
			var vals []int64
			if v.lossy {
				// Lossy HCOMPRESS of full-range data overshoots the type's
				// range on decode (funpack: "datatype conversion overflow",
				// cfitsio-written files alike); use a noisy sky instead.
				vals = noisePattern(n, int(v.shape[0]), ty.lo, ty.hi, uint64(100*ti+vi+1))
			} else {
				vals = intPattern(n, int(v.shape[0]), ty.lo, ty.hi, v.mask, uint64(100*ti+vi+1))
			}
			ints := make([]float64, n)
			for i, x := range vals {
				ints[i] = float64(x)
			}
			ty, v := ty, v
			cases = append(cases, interopCase{
				name: ty.name + "_" + v.name, shape: v.shape, tile: v.tile, npix: n,
				lossless: !v.lossy, ints: ints,
				write: func(path string) error {
					return writeOne(path, func(f *File) error { return ty.write(f, v.shape, vals, v.opts) })
				},
			})
		}
	}

	{
		// A caller header copied from an uncompressed primary HDU (SIMPLE,
		// EXTEND, BLOCKED, CHECKSUM, DATASUM, structural keys): renamed or
		// dropped like cfitsio does, so fitsverify accepts the extension.
		shape := []int64{36, 22}
		vals := intPattern(36*22, 36, -32768, 32767, false, 999)
		ints := make([]float64, len(vals))
		for i, x := range vals {
			ints[i] = float64(x)
		}
		hdr := primaryStyleHeader(16, shape)
		cases = append(cases, interopCase{
			name: "i16_rice_caller_header", shape: shape, tile: []int64{36, 1}, npix: len(vals),
			lossless: true, ints: ints,
			write: func(path string) error {
				return writeOne(path, func(f *File) error {
					data := make([]int16, len(vals))
					for i, x := range vals {
						data[i] = int16(x)
					}
					_, err := AppendCompressedImage(f, hdr, shape, data, CompressOptions{})
					return err
				})
			},
		})
	}

	{
		// A caller header that already carries compression keywords
		// (ZBITPIX, ZCMPTYPE, ZTILE1, ZNAXIS1, ZQUANTIZ, ZBLANK) and both EXTEND and
		// ZEXTEND / CHECKSUM and ZHECKSUM: ours must win and no keyword may
		// appear twice (fitsverify: "Keyword ... is duplicated").
		shape := []int64{36, 22}
		vals := intPattern(36*22, 36, -32768, 32767, false, 998)
		ints := make([]float64, len(vals))
		for i, x := range vals {
			ints[i] = float64(x)
		}
		hdr := primaryStyleHeader(16, shape)
		hdr.Set("ZBITPIX", int64(99), "stale")
		hdr.Set("ZCMPTYPE", "GZIP_1", "stale")
		hdr.Set("ZTILE1", int64(7), "stale")
		hdr.Set("ZNAXIS1", int64(3), "stale")
		hdr.Set("ZQUANTIZ", "NO_DITHER", "stale")
		// a stale ZBLANK would make every reader turn pixel value -1 (the
		// pattern's mid value, present in the data) into a null
		hdr.Set("ZBLANK", int64(-1), "stale")
		hdr.Set("ZEXTEND", true, "already renamed")
		hdr.Set("ZHECKSUM", "BBBBBBBBBBBBBBBB", "already renamed")
		cases = append(cases, interopCase{
			name: "i16_rice_caller_zkeys", shape: shape, tile: []int64{36, 1}, npix: len(vals),
			lossless: true, ints: ints,
			write: func(path string) error {
				return writeOne(path, func(f *File) error {
					data := make([]int16, len(vals))
					for i, x := range vals {
						data[i] = int16(x)
					}
					_, err := AppendCompressedImage(f, hdr, shape, data, CompressOptions{})
					return err
				})
			},
		})
	}

	type fVariant struct {
		name    string
		shape   []int64
		tile    []int64
		opts    CompressFloatOptions
		quantiz string
		content string // "noise", "const" (constant rows -> fallback tiles), "nullval", "nulldef", "realnull", "specials", "huge"
		hdr     bool   // write with a primary-style caller header
	}
	var fvars []fVariant
	dithers := []struct {
		name    string
		opts    CompressFloatOptions
		quantiz string
	}{
		{"nodither", CompressFloatOptions{NoDither: true}, "NO_DITHER"},
		{"sd1_z1", CompressFloatOptions{DitherMethod: compress.SubtractiveDither1, ZDither0: 1}, "SUBTRACTIVE_DITHER_1"},
		{"sd1_z77", CompressFloatOptions{DitherMethod: compress.SubtractiveDither1, ZDither0: 77}, "SUBTRACTIVE_DITHER_1"},
		{"sd2_z9999", CompressFloatOptions{DitherMethod: compress.SubtractiveDither2, ZDither0: 9999}, "SUBTRACTIVE_DITHER_2"},
	}
	for _, a := range []compress.Algorithm{compress.RICE1, compress.GZIP1, compress.GZIP2} {
		for _, d := range dithers {
			o := d.opts
			o.Algorithm = a
			fvars = append(fvars, fVariant{a.String() + "_" + d.name, []int64{36, 22}, []int64{36, 1}, o, d.quantiz, "noise", false})
			o.TileShape = []int64{12, 11}
			fvars = append(fvars, fVariant{a.String() + "_" + d.name + "_t12x11", []int64{36, 22}, []int64{12, 11}, o, d.quantiz, "noise", false})
		}
	}
	fvars = append(fvars,
		fVariant{"RICE_1_default", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"RICE_1_abs_q", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{QLevel: -0.25, ZDither0: 5}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"RICE_1_fallback", []int64{36, 22}, []int64{36, 4}, CompressFloatOptions{TileShape: []int64{36, 4}}, "SUBTRACTIVE_DITHER_1", "const", false},
		fVariant{"GZIP_2_cube", []int64{9, 7, 3}, []int64{9, 7, 1}, CompressFloatOptions{Algorithm: compress.GZIP2, TileShape: []int64{9, 7, 1}, ZDither0: 3}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"RICE_1_wrap_100x100", []int64{100, 100}, []int64{100, 100}, CompressFloatOptions{TileShape: []int64{100, 100}, ZDither0: 9000}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"RICE_1_wrapin_z10000", []int64{100, 100}, []int64{100, 100}, CompressFloatOptions{TileShape: []int64{100, 100}, ZDither0: 10000}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"RICE_1_1d", []int64{300}, []int64{300}, CompressFloatOptions{ZDither0: 2}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"GZIP_1_1d_sd2", []int64{300}, []int64{300}, CompressFloatOptions{Algorithm: compress.GZIP1, DitherMethod: compress.SubtractiveDither2, ZDither0: 4}, "SUBTRACTIVE_DITHER_2", "noise", false},
		fVariant{"RICE_1_nullcheck", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{Nullcheck: true, NullValue: -999}, "SUBTRACTIVE_DITHER_1", "nullval", false},
		fVariant{"RICE_1_nullcheck_t12x11", []int64{36, 22}, []int64{12, 11}, CompressFloatOptions{Nullcheck: true, NullValue: -999, TileShape: []int64{12, 11}}, "SUBTRACTIVE_DITHER_1", "nullval", false},
		fVariant{"RICE_1_caller_header", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{}, "SUBTRACTIVE_DITHER_1", "noise", true},
		// BlockSize is honoured: the stream must really be coded in
		// 16-pixel blocks, as ZVAL1 = 16 tells the readers
		fVariant{"RICE_1_bs16", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{BlockSize: 16, ZDither0: 6}, "SUBTRACTIVE_DITHER_1", "noise", false},
		fVariant{"RICE_1_bs16_t12x11_sd2", []int64{36, 22}, []int64{12, 11}, CompressFloatOptions{BlockSize: 16, TileShape: []int64{12, 11}, DitherMethod: compress.SubtractiveDither2, ZDither0: 8}, "SUBTRACTIVE_DITHER_2", "noise", false},
		// Nullcheck with NullValue left 0: the documented default -1e30
		// is the null value (and 0 is not)
		fVariant{"RICE_1_nullcheck_default", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{Nullcheck: true}, "SUBTRACTIVE_DITHER_1", "nulldef", false},
		// IEEE specials: +-Inf are nulls like NaN; +-0 and denormals are 0
		fVariant{"RICE_1_specials", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{ZDither0: 2}, "SUBTRACTIVE_DITHER_1", "specials", false},
		fVariant{"GZIP_2_specials_sd2_t12x11", []int64{36, 22}, []int64{12, 11}, CompressFloatOptions{Algorithm: compress.GZIP2, TileShape: []int64{12, 11}, DitherMethod: compress.SubtractiveDither2, ZDither0: 3}, "SUBTRACTIVE_DITHER_2", "specials", false},
		fVariant{"RICE_1_specials_nodither", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{NoDither: true}, "NO_DITHER", "specials", false},
		fVariant{"RICE_1_specials_fallback", []int64{36, 22}, []int64{36, 4}, CompressFloatOptions{TileShape: []int64{36, 4}}, "SUBTRACTIVE_DITHER_1", "huge", false},
		// a real -1e30 (the default NullValue) in tiles that also hold a
		// NaN, Nullcheck off: it must stay a value (the null sentinel the
		// quantizer is given may equal no real pixel; golden
		// rice_f32_nan_realnull_t40x5)
		fVariant{"RICE_1_realnull", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{ZDither0: 4}, "SUBTRACTIVE_DITHER_1", "realnull", false},
		fVariant{"GZIP_2_realnull_t12x11", []int64{36, 22}, []int64{12, 11}, CompressFloatOptions{Algorithm: compress.GZIP2, TileShape: []int64{12, 11}}, "SUBTRACTIVE_DITHER_1", "realnull", false},
		// tile longer than the image axis (ZTILE1 = 100 > NAXIS1)
		fVariant{"GZIP_1_t100x5", []int64{36, 22}, []int64{100, 5}, CompressFloatOptions{Algorithm: compress.GZIP1, TileShape: []int64{100, 5}}, "SUBTRACTIVE_DITHER_1", "noise", false},
	)
	for vi, v := range fvars {
		n := 1
		for _, d := range v.shape {
			n *= int(d)
		}
		r := lcg(uint64(7000 + vi))
		vals := make([]float32, n)
		for i := range vals {
			row, col := i/int(v.shape[0]), i%int(v.shape[0])
			x := float32(100 + float64(0.3*float64(col)) + float64(0.2*float64(row)) + float64(5*r.gauss()))
			if v.content == "const" && row < 8 {
				x = 42.5
			}
			if v.content == "nullval" && i%41 == 7 {
				x = -999
			}
			if (v.content == "nulldef" || v.content == "realnull") && i%41 == 7 {
				x = -1e30
			}
			if v.content == "huge" && row < 4 && i%2 == 1 {
				x += 1e20 // tile 1 cannot be quantized: GZIP_COMPRESSED_DATA
			}
			switch {
			case i%37 == 5:
				x = float32(math.NaN())
			case i%29 == 3:
				x = 0
			}
			if v.content == "specials" || v.content == "huge" {
				switch {
				case i%31 == 0:
					x = float32(math.Inf(1 - 2*((i/31)%2)))
				case i%17 == 0:
					x = []float32{1e-40, -3e-42}[(i/17)%2]
				case i%13 == 0:
					x = float32(math.Copysign(0, -1))
				case i%19 == 0:
					x = []float32{math.SmallestNonzeroFloat32 * (1 << 23), -math.SmallestNonzeroFloat32 * (1 << 23)}[(i/19)%2]
				}
			}
			vals[i] = x
		}
		v := v
		nullValue := v.opts.NullValue
		if nullValue == 0 {
			nullValue = -1e30 // documented default of CompressFloatOptions.NullValue
		}
		var hdr *header.Header
		if v.hdr {
			hdr = primaryStyleHeader(-32, v.shape)
		}
		cases = append(cases, interopCase{
			name: "f32_" + v.name, shape: v.shape, tile: v.tile, npix: n, float: true,
			floats: vals, wantQuantiz: v.quantiz,
			nullcheck: v.opts.Nullcheck, nullValue: nullValue,
			write: func(path string) error {
				return writeOne(path, func(f *File) error {
					_, err := AppendCompressedFloat32Image(f, hdr, v.shape, vals, v.opts)
					return err
				})
			},
		})
	}
	// float64 images (AppendCompressedFloat64Image, ZBITPIX -64)
	dvars := []fVariant{
		{"RICE_1_default", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{}, "SUBTRACTIVE_DITHER_1", "noise", false},
		{"GZIP_2_sd2_t12x11", []int64{36, 22}, []int64{12, 11}, CompressFloatOptions{Algorithm: compress.GZIP2, TileShape: []int64{12, 11}, DitherMethod: compress.SubtractiveDither2, ZDither0: 3}, "SUBTRACTIVE_DITHER_2", "noise", false},
		{"GZIP_1_nodither", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{Algorithm: compress.GZIP1, NoDither: true}, "NO_DITHER", "noise", false},
		{"RICE_1_abs_q", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{QLevel: -0.25, ZDither0: 5}, "SUBTRACTIVE_DITHER_1", "noise", false},
		{"RICE_1_fallback", []int64{36, 22}, []int64{36, 4}, CompressFloatOptions{TileShape: []int64{36, 4}}, "SUBTRACTIVE_DITHER_1", "const", false},
		{"RICE_1_specials", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{ZDither0: 2}, "SUBTRACTIVE_DITHER_1", "specials", false},
		{"RICE_1_specials_fallback", []int64{36, 22}, []int64{36, 4}, CompressFloatOptions{TileShape: []int64{36, 4}}, "SUBTRACTIVE_DITHER_1", "huge", false},
		{"RICE_1_nullcheck", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{Nullcheck: true, NullValue64: -999.25}, "SUBTRACTIVE_DITHER_1", "nullval", false},
		{"RICE_1_realnull", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{ZDither0: 4}, "SUBTRACTIVE_DITHER_1", "realnull", false},
		{"RICE_1_precise_1e8", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{}, "SUBTRACTIVE_DITHER_1", "precise", false},
		{"RICE_1_caller_header", []int64{36, 22}, []int64{36, 1}, CompressFloatOptions{}, "SUBTRACTIVE_DITHER_1", "noise", true},
	}
	for vi, v := range dvars {
		n := 1
		for _, d := range v.shape {
			n *= int(d)
		}
		r := lcg(uint64(9000 + vi))
		vals := make([]float64, n)
		for i := range vals {
			row, col := i/int(v.shape[0]), i%int(v.shape[0])
			// Each product is rounded on its own (float64(...)): without
			// the conversions Go fuses them into FMAs on arm64 and the
			// input, the written bytes and their recorded readings would
			// differ by platform.
			x := 100 + float64(0.3*float64(col)) + float64(0.2*float64(row)) + float64(5*r.gauss())
			switch v.content {
			case "const":
				if row < 8 {
					x = 42.5
				}
			case "nullval":
				if i%41 == 7 {
					x = -999.25
				}
			case "realnull":
				if i%41 == 7 {
					x = -1e30
				}
			case "huge":
				if row < 4 && i%2 == 1 {
					x += 1e20
				}
			case "precise":
				x = 1e8 + float64(1e-4*float64(col)) + float64(0.01*r.gauss())
			}
			switch {
			case i%37 == 5:
				x = math.NaN()
			case i%29 == 3:
				x = 0
			}
			if v.content == "specials" || v.content == "huge" {
				switch {
				case i%31 == 0:
					x = math.Inf(1 - 2*((i/31)%2))
				case i%17 == 0:
					x = []float64{1e-310, -4e-320}[(i/17)%2]
				case i%13 == 0:
					x = math.Copysign(0, -1)
				case i%19 == 0:
					x = []float64{0x1p-1022, -0x1p-1022}[(i/19)%2]
				}
			}
			vals[i] = x
		}
		v := v
		nullValue := v.opts.NullValue64
		if nullValue == 0 {
			nullValue = -1e30
		}
		var hdr *header.Header
		if v.hdr {
			hdr = primaryStyleHeader(-64, v.shape)
		}
		cases = append(cases, interopCase{
			name: "f64_" + v.name, shape: v.shape, tile: v.tile, npix: n, double: true,
			doubles: vals, wantQuantiz: v.quantiz,
			nullcheck: v.opts.Nullcheck, nullValue64: nullValue,
			write: func(path string) error {
				return writeOne(path, func(f *File) error {
					_, err := AppendCompressedFloat64Image(f, hdr, v.shape, vals, v.opts)
					return err
				})
			},
		})
	}
	// inputs on which cfitsio's quantizer overflows: stored losslessly
	// (TestQuantizerOverflowStoredLosslessly), so every reader must give
	// them back exactly
	cases = append(cases, quantOverflowInteropCases()...)
	return cases
}

// primaryStyleHeader is the header a caller would copy from an
// uncompressed primary image HDU.
func primaryStyleHeader(bp int64, shape []int64) *header.Header {
	h := header.New()
	h.Set("SIMPLE", true, "file does conform to FITS standard")
	h.Set("BITPIX", bp, "")
	h.Set("NAXIS", int64(len(shape)), "")
	for i, n := range shape {
		h.Set(fmt.Sprintf("NAXIS%d", i+1), n, "")
	}
	h.Set("EXTEND", true, "FITS dataset may contain extensions")
	h.Set("BLOCKED", true, "deprecated")
	h.Set("OBJECT", "NGC 1234", "target")
	h.Set("CHECKSUM", "AAAAAAAAAAAAAAAA", "of the uncompressed HDU")
	h.Set("DATASUM", "1234567890", "of the uncompressed HDU")
	h.Add("HISTORY", nil, "copied from the uncompressed image")
	return h
}

// TestFunpackMisdecodesNullTiledFloats documents the funpack defect the
// interop matrix works around, on a file cfitsio itself wrote
// (golden rice_f32_nan_z1_t16x8: float32, NaNs, 16x8 tiles on a 40x25
// image): funpack 1.7.0 (CFITSIO 4.7.0)'s output differs from cfitsio's fits_read_img and
// from astropy, which agree with each other. If a future funpack fixes
// this, this test fails and the exclusion in TestCompressedWriterInterop
// should be removed.
func TestFunpackMisdecodesNullTiledFloats(t *testing.T) {
	src := goldenFITS(t, filepath.Join(testutil.RepoRoot(t), imgGoldenDir, "rice_f32_nan_z1_t16x8.fits"))
	out := filepath.Join(t.TempDir(), "unpacked.fits")
	funpackTo(t, out, src)
	direct := cfitsioPixels(t, src)
	unpacked := cfitsioPixels(t, out)
	differ := 0
	for i := range direct {
		if !samePixel(direct[i], unpacked[i], true) {
			differ++
		}
	}
	if differ == 0 {
		t.Fatalf("funpack now decodes null + narrow-tile float images like fits_read_img; drop the funpack exclusion in TestCompressedWriterInterop")
	}
	t.Logf("funpack differs from fits_read_img on %d/%d pixels of cfitsio's own file", differ, len(direct))
}

// TestPLIOUint16AstropyMisreadsCfitsio is why AppendCompressedImage
// refuses PLIO_1 for uint16: on cfitsio's own PLIO_1 uint16 file (golden
// plio_u16) cfitsio decodes the input exactly while astropy returns every
// pixel 32768 too high.
func TestPLIOUint16AstropyMisreadsCfitsio(t *testing.T) {
	root := testutil.RepoRoot(t)
	src := goldenFITS(t, filepath.Join(root, imgGoldenDir, "plio_u16.fits"))
	in := imgGoldenByName(t, "plio_u16").readInput(t)
	direct := cfitsioPixels(t, src)
	for i := range direct {
		want := float64(binary.LittleEndian.Uint16(in[2*i:]))
		if direct[i] == nil || *direct[i] != want {
			t.Fatalf("cfitsio pixel %d = %v, input %v", i, fmtp(direct[i]), want)
		}
	}
	a := astropyDecode(t, []interopCase{{path: src}})[src]
	if a.Error != "" {
		t.Fatalf("astropy: %s", a.Error)
	}
	same, low, lowOff := 0, 0, 0
	for i, p := range a.Pixels {
		if p != nil && *p == *direct[i] {
			same++
		}
		if *direct[i] < 32768 {
			low++
			if p != nil && *p == *direct[i]+32768 {
				lowOff++
			}
		}
	}
	if same != 0 || low == 0 || lowOff != low {
		t.Fatalf("astropy: %d/%d pixels equal cfitsio, %d/%d values < 32768 read +32768; if astropy now agrees with cfitsio, PLIO uint16 writes can be enabled", same, len(direct), lowOff, low)
	}
	t.Logf("astropy agrees with cfitsio on 0/%d pixels; every value < 32768 (%d) comes back +32768", len(direct), low)
}
