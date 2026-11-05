package fits

import (
	"math"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/internal/crefgen"
)

// Limits of the float quantizer that cfitsio 4.7.0 does not guard, the
// caller null value, and the caller SIMPLE card.

// quantOverflowCase is an input on which fits_quantize_float/_double's
// arithmetic overflows: the iqfactor (int64(zeropt/delta + 0.5)) when
// |minval| / delta >= 2^63, or the int32 pixel when the centered zero point
// (minval + maxval)/2, summed in float, lands so far from the true center
// that the maximum quantizes past 2^31-1, or a null-branch tile whose zero
// point rounds so far that the minimum lands on the null sentinels. cfitsio
// writes such tiles with status 0 and they decode to wrong values. Our
// writer must store them so that every pixel reads back exactly.
type quantOverflowCase struct {
	name   string
	shape  []int64
	f32    []float32
	f64    []float64
	opts   CompressFloatOptions
	golden string // the cfitsio golden of the same input, if any
	// exact is the number of leading pixels (the overflowing tiles) that
	// must read back exactly; 0 means all. The rest are quantized
	// normally and must be within half a quantization step (-QLevel).
	exact int
}

func quantOverflowCases() []quantOverflowCase {
	fill32 := func(n int, v float32) []float32 {
		o := make([]float32, n)
		for i := range o {
			o[i] = v
		}
		return o
	}
	ramp := crefgen.ImgFloat(64, 64, 20, 0)
	big1e17nan := fill32(100, 1e17)
	big1e17nan[17] = float32(math.NaN())
	f64ramp := crefgen.ImgDouble(64, 64, 19, 0)
	f64nan := append([]float64(nil), f64ramp...)
	f64nan[5] = math.NaN()
	return []quantOverflowCase{
		// golden img2/rice_f32_zpt_overflow_abs: cfitsio decodes 100 as -0.0922
		{"f32_const100_q-1e-20", []int64{64}, crefgen.ImgFloat(64, 64, 19, 0), nil,
			CompressFloatOptions{QLevel: -1e-20, NoDither: true}, "img2/rice_f32_zpt_overflow_abs", 0},
		// golden img2/rice_f32_center_overflow: cfitsio decodes the maximum
		// 14294966272 as 9999998976
		{"f32_center_ramp_q-1", []int64{64}, ramp, nil,
			CompressFloatOptions{QLevel: -1, NoDither: true}, "img2/rice_f32_center_overflow", 0},
		{"f32_center_ramp_q-1_sd1", []int64{64}, ramp, nil,
			CompressFloatOptions{QLevel: -1, ZDither0: 3}, "", 0},
		// golden img2/rice_f32_null_sentinel_hit: row 1's minimum quantizes
		// to NULL_VALUE and cfitsio decodes the real pixel as a null
		{"f32_null_sentinel_hit", []int64{64, 2}, crefgen.ImgFloat(128, 64, 21, 0), nil,
			CompressFloatOptions{QLevel: -0x1.528c8cp+0, NoDither: true}, "img2/rice_f32_null_sentinel_hit", 64},
		// golden img64/rice_f64_zpt_overflow_abs: cfitsio decodes 1e10 as -9223372.002
		{"f64_1e10_q-1e-12", []int64{64}, nil, f64ramp,
			CompressFloatOptions{QLevel: -1e-12, NoDither: true}, "img64/rice_f64_zpt_overflow_abs", 0},
		{"f64_1e10_q-1e-12_nan_sd2", []int64{64}, nil, f64nan,
			CompressFloatOptions{QLevel: -1e-12, DitherMethod: compress.SubtractiveDither2, ZDither0: 9}, "", 0},
		// the round-2 verifier's probes
		{"f32_const100_q-1e-40", []int64{10, 10}, fill32(100, 100), nil, CompressFloatOptions{QLevel: -1e-40}, "", 0},
		{"f32_const1e17_q-0.001", []int64{10, 10}, fill32(100, 1e17), nil, CompressFloatOptions{QLevel: -0.001}, "", 0},
		// null branch: zeropt = 1e17 + 2147483.637 rounds to a multiple of
		// 16 in float64, thousands of steps off
		{"f32_const1e17_nan_q-0.001", []int64{10, 10}, big1e17nan, nil, CompressFloatOptions{QLevel: -0.001}, "", 0},
	}
}

// writeQuantOverflow writes c with our writer and returns the path.
func (c quantOverflowCase) write(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "q.fits")
	err := writeOne(path, func(f *File) error {
		var err error
		if c.f64 != nil {
			_, err = AppendCompressedFloat64Image(f, nil, c.shape, c.f64, c.opts)
		} else {
			_, err = AppendCompressedFloat32Image(f, nil, c.shape, c.f32, c.opts)
		}
		return err
	})
	if err != nil {
		t.Fatalf("write: %v", err)
	}
	return path
}

// input returns the case's pixels as cfitsio's writer sees them, widened
// (NaN for a null).
func (c quantOverflowCase) input() []float64 {
	if c.f64 != nil {
		return append([]float64(nil), c.f64...)
	}
	out := make([]float64, len(c.f32))
	for i, v := range c.f32 {
		out[i] = float64(v)
	}
	return out
}

// TestQuantizerOverflowStoredLosslessly: where cfitsio's quantizer
// overflows (see quantOverflowCase) our writer must not write wrong
// pixels. The truth is the input itself: the tile is stored in
// GZIP_COMPRESSED_DATA (lossless; the format is proven bit-exact against
// cfitsio by the fallback goldens) and reads back bit for bit. For the
// cases with a cfitsio golden, the golden (cfitsio's decode of its own
// file) is checked to really be wrong, so the deliberate difference stays
// pinned. TestCompressedWriterInterop judges the same inputs with
// cfitsio, funpack and astropy (the "quantovf_" cases).
func TestQuantizerOverflowStoredLosslessly(t *testing.T) {
	goldens := map[string]imgGolden{}
	for _, g := range append(img2GoldenCases(t), img64GoldenCases(t)...) {
		goldens[filepath.Base(g.dir)+"/"+g.Name] = g
	}
	for _, c := range quantOverflowCases() {
		t.Run(c.name, func(t *testing.T) {
			in := c.input()
			if c.golden != "" {
				g, ok := goldens[c.golden]
				if !ok {
					t.Fatalf("committed golden %s missing (failure, not a skip)", c.golden)
				}
				if _, ok := knownCfitsioOverflow[g.Name]; !ok {
					t.Fatalf("%s not listed in knownCfitsioOverflow", g.Name)
				}
				var gin []float64
				if g.Type == "f64" {
					gin = g.doubleInput(t)
				} else {
					for _, v := range g.floatInput(t) {
						gin = append(gin, float64(v))
					}
				}
				for i := range in {
					if math.Float64bits(gin[i]) != math.Float64bits(in[i]) {
						t.Fatalf("golden input pixel %d = %v, case input %v", i, gin[i], in[i])
					}
				}
				dec := g.decTruth(t, "cfitsio's decode", g.decCandidate(t))
				wrong := 0
				for i := range in {
					if math.IsNaN(dec[i]) != math.IsNaN(in[i]) || math.Abs(dec[i]-in[i]) > 1e3*math.Abs(g.QLevel) {
						wrong++
					}
				}
				if wrong == 0 {
					t.Errorf("cfitsio's decode of its file is within the quantization step everywhere: the defect (%s) is gone, revisit knownCfitsioOverflow", knownCfitsioOverflow[g.Name])
				} else {
					t.Logf("cfitsio 4.7.0: %d/%d pixels wrong in its own file (%s)", wrong, len(in), knownCfitsioOverflow[g.Name])
				}
			}

			path := c.write(t)
			f, h := openCompressed(t, path)
			defer f.Close()
			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatal(err)
			}
			bad := 0
			for i := range in {
				w := in[i]
				if c.f64 == nil {
					w = float64(float32(w))
				}
				tol := 0.0
				if c.exact > 0 && i >= c.exact {
					tol = 0.5 * math.Abs(float64(c.opts.QLevel)) * (1 + 1e-6)
				}
				if math.IsNaN(w) != math.IsNaN(got[i]) || (!math.IsNaN(w) && math.Abs(got[i]-w) > tol) {
					if bad < 3 {
						t.Errorf("pixel %d = %v, input %v (tolerance %g; 0: must be stored exactly)", i, got[i], w, tol)
					}
					bad++
				}
			}
			if bad > 0 {
				t.Errorf("%d/%d pixels wrong", bad, len(in))
			}
		})
	}
}

// decCandidate is a placeholder reading for decTruth when the golden holds
// explicit values (decTruth returns them after checking cand against
// them); for a digest it would have to be a real reading. The overflow
// goldens are small (explicit), so the cfitsio-written file is read with
// our reader -- whose agreement with cfitsio on that file is
// TestReaderDecodesCfitsioImg2Goldens.
func (g imgGolden) decCandidate(t *testing.T) []float64 {
	t.Helper()
	f, h := openCompressed(t, filepath.Join(g.dir, g.Name+".fits"))
	defer f.Close()
	v, err := ReadPixelsCompressed[float64](h)
	if err != nil {
		t.Fatal(err)
	}
	return v
}

// quantOverflowInteropCases are the quantOverflowCases as interop cases.
func quantOverflowInteropCases() []interopCase {
	var out []interopCase
	for _, c := range quantOverflowCases() {
		c := c
		n := len(c.f32) + len(c.f64)
		tile := []int64{c.shape[0], 1}
		if len(c.shape) == 1 {
			tile = []int64{c.shape[0]}
		}
		quantiz := "SUBTRACTIVE_DITHER_1"
		switch {
		case c.opts.NoDither:
			quantiz = "NO_DITHER"
		case c.opts.DitherMethod == compress.SubtractiveDither2:
			quantiz = "SUBTRACTIVE_DITHER_2"
		}
		ic := interopCase{
			name: "quantovf_" + c.name, shape: c.shape, tile: tile, npix: n,
			wantQuantiz: quantiz, nullValue: -1e30, nullValue64: -1e30,
			write: func(path string) error {
				return writeOne(path, func(f *File) error {
					var err error
					if c.f64 != nil {
						_, err = AppendCompressedFloat64Image(f, nil, c.shape, c.f64, c.opts)
					} else {
						_, err = AppendCompressedFloat32Image(f, nil, c.shape, c.f32, c.opts)
					}
					return err
				})
			},
		}
		if c.f64 != nil {
			ic.double, ic.doubles = true, c.f64
		} else {
			ic.float, ic.floats = true, c.f32
		}
		out = append(out, ic)
	}
	return out
}

// TestCompressedWriterNullValueMatchesCfitsio: a caller null value
// (Nullcheck with NullValue / NullValue64) against cfitsio's
// fits_write_imgnull with the same null value (golden/imgw nullval cases:
// -999 for float32 and float64, and the denormal 1.04488e-40, which
// cfitsio compares raw before any denormal is read as 0): the same tiles,
// ZSCALE and ZZERO as cfitsio's file, and cfitsio's decode of its file has
// a null exactly where the input is NaN or equals the null value.
func TestCompressedWriterNullValueMatchesCfitsio(t *testing.T) {
	ran := 0
	for _, g := range imgwGoldenCases(t) {
		if g.NullValue == 0 {
			continue
		}
		ran++
		t.Run(g.Name, func(t *testing.T) {
			matchCfitsioGolden(t, g)
			var in []float64
			if g.Type == "f64" {
				in = g.doubleInput(t)
			} else {
				for _, v := range g.floatInput(t) {
					in = append(in, float64(v))
				}
			}
			path, err := g.writeOurs(t, t.TempDir())
			if err != nil {
				t.Fatal(err)
			}
			ff, h := openCompressed(t, path)
			defer ff.Close()
			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatal(err)
			}
			truth := g.decTruth(t, "our file read as float64", got)
			nnull := 0
			for i, v := range in {
				want := math.IsNaN(v) || v == g.NullValue
				if want {
					nnull++
				}
				if math.IsNaN(truth[i]) != want {
					t.Errorf("pixel %d (input %v, null value %v): cfitsio null = %v, want %v", i, v, g.NullValue, math.IsNaN(truth[i]), want)
				}
			}
			if nnull == 0 {
				t.Fatal("input holds no null")
			}
		})
	}
	if ran != 3 {
		t.Fatalf("%d imgw null-value cases, want 3", ran)
	}
}

// TestCompressFloatNullValueOfOtherWidthRefused: the float32 writer takes
// NullValue and the float64 writer NullValue64. A non-zero value in the
// other field would be ignored, and the pixels the caller asked to null
// kept as values; both writers refuse it instead.
func TestCompressFloatNullValueOfOtherWidthRefused(t *testing.T) {
	data32 := make([]float32, 100)
	data64 := make([]float64, 100)
	for i := range data32 {
		data32[i] = float32(100 + i%10)
		data64[i] = float64(100 + i%10)
	}
	data32[5], data64[5] = -999, -999
	cases := []struct {
		name string
		run  func(f *File) error
		want string
	}{
		{"float64 writer with NullValue", func(f *File) error {
			_, err := AppendCompressedFloat64Image(f, nil, []int64{10, 10}, data64, CompressFloatOptions{Nullcheck: true, NullValue: -999})
			return err
		}, "NullValue64"},
		{"float32 writer with NullValue64", func(f *File) error {
			_, err := AppendCompressedFloat32Image(f, nil, []int64{10, 10}, data32, CompressFloatOptions{Nullcheck: true, NullValue64: -999})
			return err
		}, "NullValue"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			err := writeOne(filepath.Join(t.TempDir(), "n.fits"), c.run)
			if err == nil || !strings.Contains(err.Error(), c.want) {
				t.Fatalf("err = %v, want a refusal naming %s", err, c.want)
			}
		})
	}
}

// TestCompressedWriterSimpleFMatchesCfitsio: a caller header with SIMPLE =
// F is written as ZSIMPLE = F, as cfitsio's fits_img_compress writes it
// from a SIMPLE = F source image (golden img2/rice_i16_simple_f, its
// source image rice_i16_simple_f.src.fits; the full keyword comparison is
// TestCompressedWriterMatchesCfitsio/img2/rice_i16_simple_f).
func TestCompressedWriterSimpleFMatchesCfitsio(t *testing.T) {
	var g imgGolden
	for _, c := range img2GoldenCases(t) {
		if c.Name == "rice_i16_simple_f" {
			g = c
		}
	}
	if g.Name == "" {
		t.Fatal("committed golden img2/rice_i16_simple_f missing (failure, not a skip)")
	}
	gf, gc := openCompressed(t, filepath.Join(g.dir, g.Name+".fits"))
	defer gf.Close()
	if v, err := gc.Header().Bool("ZSIMPLE"); err != nil || v {
		t.Fatalf("cfitsio's file: ZSIMPLE = %v (%v), want F", v, err)
	}
	path, err := g.writeOurs(t, t.TempDir())
	if err != nil {
		t.Fatal(err)
	}
	ff, fc := openCompressed(t, path)
	defer ff.Close()
	if v, err := fc.Header().Bool("ZSIMPLE"); err != nil || v {
		t.Fatalf("our file: ZSIMPLE = %v (%v), cfitsio writes F", v, err)
	}
}
