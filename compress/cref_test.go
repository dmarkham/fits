package compress

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/crefgen"
	"github.com/dmarkham/fits/internal/testutil"
)

// C reference validation. These tests load the golden fixtures produced
// by testdata/cref/cref.c (which links against libcfitsio) and diff the
// Go port output against the real cfitsio implementation.
//
// The fixtures are committed so the tests need neither a C compiler nor
// libcfitsio. The inputs are not committed: cref.c generates them with a
// fixed-seed LCG and records their sha256; internal/crefgen regenerates
// them and must reproduce that hash. cfitsio's quantized integers are an
// exact digest (testutil.Values: sha256 over every value plus a sample).
// To regenerate (and check reproducibility) locally:
//
//	make -C compress/testdata/cref goldens   # rewrite in place
//	make -C compress/testdata/cref verify    # regenerate + byte-compare
//
// A missing or extra golden is a test FAILURE, never a skip: a partial
// deletion of the fixtures must not silently shrink the coverage.

const cgoldenDir = "testdata/cref/golden"

// crefStats mirrors the fields cref.c writes in each .out.json file.
type crefCase struct {
	Name         string  `json:"name"`
	NX           int64   `json:"nx"`
	NY           int64   `json:"ny"`
	Nullcheck    int     `json:"nullcheck"`
	NullValue    float64 `json:"nullvalue"`
	QLevel       float64 `json:"qlevel"`
	DitherMethod int     `json:"dither_method"`
	Row          int64   `json:"row"`
	Fill         string  `json:"fill"`      // cref.c generator (internal/crefgen.Fill)
	InSHA256     string  `json:"in_sha256"` // sha256 of the generated float32 LE input
	Stats        struct {
		NGood  int64   `json:"ngood"`
		Min    float64 `json:"min"`
		Max    float64 `json:"max"`
		Mean   float64 `json:"mean"`
		Sigma  float64 `json:"sigma"`
		Noise1 float64 `json:"noise1"`
		Noise2 float64 `json:"noise2"`
		Noise3 float64 `json:"noise3"`
		Noise5 float64 `json:"noise5"`
	} `json:"stats"`
	Quantize struct {
		Return  int     `json:"return"`
		BScale  float64 `json:"bscale"`
		BZero   float64 `json:"bzero"`
		IMinVal int32   `json:"iminval"`
		IMaxVal int32   `json:"imaxval"`
	} `json:"quantize"`
	// Out is fits_quantize_float's idata (absent/null: not quantized).
	Out testutil.Values[testutil.NullInt64] `json:"out"`
}

// loadCrefCase reads the golden JSON of a case and regenerates its input
// float32 array, which must hash to the input cfitsio was given. A
// missing fixture or a drifted generator fails the test.
func loadCrefCase(t *testing.T, name string) ([]float32, *crefCase) {
	t.Helper()
	jsonPath := filepath.Join(cgoldenDir, name+".out.json")
	jsonData, err := os.ReadFile(jsonPath)
	if err != nil {
		t.Fatalf("committed golden missing (failure, not a skip): %v", err)
	}
	var c crefCase
	if err := json.Unmarshal(jsonData, &c); err != nil {
		t.Fatalf("parse %s: %v", jsonPath, err)
	}
	arr, err := crefgen.Fill(c.Fill, c.NX*c.NY)
	if err != nil {
		t.Fatal(err)
	}
	if err := crefgen.Check(name, crefgen.LEFloat32(arr), c.InSHA256); err != nil {
		t.Fatal(err)
	}
	return arr, &c
}

// crefAllCases is the canonical list maintained in sync with cref.c
// (`cref list` prints the same names).
var crefAllCases = []string{
	"gauss_64x64_q4_nodither",
	"gauss_64x64_q4_dither1",
	"gauss_64x64_q4_dither1_row7",
	"gauss_64x64_q4_dither2_zeros",
	"gauss_64x64_q4_nulls",
	"gauss_64x64_q4_real_nans",
	"gauss_64x64_q16_dither1",
	"gauss_64x64_abs_step",
	"gradient_64x64_q4_dither1",
	"constant_16x16",
	"tiny_6x1",
	"runs3_64x64_q4_dither1",
	"runs3_8x40_q4_nulls_dither2",
	"runs3_stale_200x2_q4_dither1",
	"gauss_64x64_q0_dither1",
	"tiny_2x8",
	"tiny_4x5_abs",
	"tiny_2x2_abs",
	"range_fsub_256x1",
	"wide_center_64x1",
	"center_fsum_64x1",
	"margin_64x1",
	"halfmargin_64x1",
	"n5_one_9x3",
	"n5_two_10x3",
	"n5_zero_nulls_9x3",
	"n1_short_nulls_5x3",
	"one_good_3x3",
	"clip_tie_27x1",
	"abs_nulls_10x4",
	"abs_nulls_10x4_dither1",
	"all_null_4x4_q4_dither1",
	"all_null_12x3_q4_nodither",
	"all_null_4x4_abs_dither2",
	"all_nan_10x10_q16_dither2",
	"margin_far_64x1",
}

// TestCrefGoldensComplete fails when the committed quantizer goldens and
// crefAllCases disagree in either direction (a golden deleted, or a case
// added to cref.c without its golden being committed and listed here).
func TestCrefGoldensComplete(t *testing.T) {
	want := map[string]bool{}
	for _, n := range crefAllCases {
		want[n] = true
		if _, err := os.Stat(filepath.Join(cgoldenDir, n+".out.json")); err != nil {
			t.Errorf("golden %s.out.json missing: %v", n, err)
		}
	}
	got, err := filepath.Glob(filepath.Join(cgoldenDir, "*.out.json"))
	if err != nil {
		t.Fatal(err)
	}
	if len(got) != len(crefAllCases) {
		t.Errorf("%d *.out.json goldens committed, crefAllCases lists %d", len(got), len(crefAllCases))
	}
	for _, p := range got {
		n := strings.TrimSuffix(filepath.Base(p), ".out.json")
		if !want[n] {
			t.Errorf("golden %s is not in crefAllCases (untested)", n)
		}
	}
}

// TestCrefRealNaNsInputHasNaNs pins the fixture itself: the
// gauss_64x64_q4_real_nans input must contain the raw NaNs (every 41st
// pixel) so the NaN -> sentinel -> NULL_VALUE path is really exercised.
// cref.c used to hash/write the input after the substitution, leaving
// zero NaNs.
func TestCrefRealNaNsInputHasNaNs(t *testing.T) {
	arr, c := loadCrefCase(t, "gauss_64x64_q4_real_nans")
	nan := 0
	for i, v := range arr {
		isNaN := math.IsNaN(float64(v))
		if isNaN {
			nan++
		}
		if isNaN != (i%41 == 0) {
			t.Fatalf("pixel %d: NaN=%v, want NaN exactly at i%%41==0", i, isNaN)
		}
	}
	if want := int((c.NX*c.NY + 40) / 41); nan != want {
		t.Fatalf("%d NaNs in input, want %d", nan, want)
	}
}

// TestReplaceSentinelWithNaNBits verifies that the lossless-path NaN
// substitution produces a bit pattern that is (a) NaN when reinterpreted
// as float32, and (b) specifically the all-ones pattern cfitsio writes
// (imcomp_float2nan uses outdata[ii] = -1). Non-sentinel pixels are
// untouched.
func TestReplaceSentinelWithNaNBits(t *testing.T) {
	buf := []float32{1.0, -1e30, 2.5, -1e30, 3.5}
	n := ReplaceSentinelWithNaNBits(buf, -1e30)
	if n != 2 {
		t.Fatalf("substituted %d, want 2", n)
	}
	// Sentinel positions must now be NaN.
	if !math.IsNaN(float64(buf[1])) {
		t.Errorf("buf[1] not NaN: %v", buf[1])
	}
	if !math.IsNaN(float64(buf[3])) {
		t.Errorf("buf[3] not NaN: %v", buf[3])
	}
	// And specifically the all-ones bit pattern cfitsio uses.
	if bits := math.Float32bits(buf[1]); bits != 0xFFFFFFFF {
		t.Errorf("buf[1] bits = %#x, want 0xFFFFFFFF (cfitsio all-ones NaN)", bits)
	}
	// Non-sentinel values must be byte-identical.
	if buf[0] != 1.0 || buf[2] != 2.5 || buf[4] != 3.5 {
		t.Errorf("non-null pixels mutated: %v", buf)
	}
}

// TestSubstituteNaN32 verifies raw-NaN → sentinel substitution.
func TestSubstituteNaN32(t *testing.T) {
	buf := []float32{1.0, float32(math.NaN()), 2.5, float32(math.NaN()), 3.5}
	n := SubstituteNaN32(buf, -1e30)
	if n != 2 {
		t.Fatalf("substituted %d, want 2", n)
	}
	if buf[1] != -1e30 || buf[3] != -1e30 {
		t.Errorf("NaN positions not replaced: %v", buf)
	}
	if buf[0] != 1.0 || buf[2] != 2.5 || buf[4] != 3.5 {
		t.Errorf("non-NaN pixels mutated: %v", buf)
	}
}

// TestImgStatsFloat32_CRef validates fnMeanSigma / fnNoise1 / fnNoise5
// against live cfitsio output byte-for-byte. The noise functions do
// a lot of float32 arithmetic, so bit-exact agreement is the bar.
func TestImgStatsFloat32_CRef(t *testing.T) {
	for _, name := range crefAllCases {
		t.Run(name, func(t *testing.T) {
			arr, c := loadCrefCase(t, name)

			// nullcheck==0 -> no nulls
			// nullcheck==1 -> sentinel already in input, pass through
			// nullcheck==2 -> input has raw NaNs; substitute before calling
			nullcheck := c.Nullcheck != 0
			nullvalue := float32(c.NullValue)
			if c.Nullcheck == 2 {
				SubstituteNaN32(arr, nullvalue)
			}

			got := ImgStatsFloat32(arr, c.NX, c.NY, nullcheck, nullvalue)

			if got.NGood != c.Stats.NGood {
				t.Errorf("ngood: got %d want %d", got.NGood, c.Stats.NGood)
			}
			if float64(got.Min) != c.Stats.Min {
				t.Errorf("min: got %v want %v", got.Min, c.Stats.Min)
			}
			if float64(got.Max) != c.Stats.Max {
				t.Errorf("max: got %v want %v", got.Max, c.Stats.Max)
			}
			if got.Mean != c.Stats.Mean {
				t.Errorf("mean: got %.17g want %.17g (diff %g)", got.Mean, c.Stats.Mean, got.Mean-c.Stats.Mean)
			}
			if got.Sigma != c.Stats.Sigma {
				t.Errorf("sigma: got %.17g want %.17g (diff %g)", got.Sigma, c.Stats.Sigma, got.Sigma-c.Stats.Sigma)
			}
			if got.Noise1 != c.Stats.Noise1 {
				t.Errorf("noise1: got %.17g want %.17g (diff %g)", got.Noise1, c.Stats.Noise1, got.Noise1-c.Stats.Noise1)
			}
			if got.Noise2 != c.Stats.Noise2 {
				t.Errorf("noise2: got %.17g want %.17g (diff %g)", got.Noise2, c.Stats.Noise2, got.Noise2-c.Stats.Noise2)
			}
			if got.Noise3 != c.Stats.Noise3 {
				t.Errorf("noise3: got %.17g want %.17g (diff %g)", got.Noise3, c.Stats.Noise3, got.Noise3-c.Stats.Noise3)
			}
			if got.Noise5 != c.Stats.Noise5 {
				t.Errorf("noise5: got %.17g want %.17g (diff %g)", got.Noise5, c.Stats.Noise5, got.Noise5-c.Stats.Noise5)
			}
		})
	}
}

// TestQuantizeFloat32_CRef validates fits_quantize_float bit-exactly
// against cfitsio. Every non-constant case should produce identical
// int32 output, identical bscale/bzero, and identical iminval/imaxval.
// Constant cases and ranges past the int32 margin must return
// quantized=false, matching the C quantize_return=0 ("out": null). The
// range cases pin cfitsio's float arithmetic: maxval - minval and
// minval + maxval are float sums (range_fsub_256x1, center_fsum_64x1),
// and the N_RESERVED_VALUES margins of both range tests (margin_64x1,
// halfmargin_64x1, reached with a positive qlevel).
func TestQuantizeFloat32_CRef(t *testing.T) {
	for _, name := range crefAllCases {
		t.Run(name, func(t *testing.T) {
			arr, c := loadCrefCase(t, name)

			nullcheck := c.Nullcheck != 0
			var nullvalue float32
			if nullcheck {
				nullvalue = float32(c.NullValue)
			}
			// nullcheck==2: the input holds raw NaNs; the writer pipeline
			// substitutes the sentinel before quantizing (as cref.c does).
			if c.Nullcheck == 2 {
				SubstituteNaN32(arr, nullvalue)
			}

			n := c.NX * c.NY
			idata := make([]int32, n)
			res, ok := QuantizeFloat32(c.Row, arr, c.NX, c.NY, nullcheck,
				nullvalue, float32(c.QLevel), DitherMethod(c.DitherMethod), idata)

			wantOK := c.Quantize.Return == 1
			if ok != wantOK {
				t.Fatalf("quantized=%v want %v", ok, wantOK)
			}
			if !ok {
				if c.Out.Present() {
					t.Fatalf("golden has quantized output but quantize returned 0")
				}
				return
			}

			if res.BScale != c.Quantize.BScale {
				t.Errorf("bscale: got %.17g want %.17g (diff %g)",
					res.BScale, c.Quantize.BScale, res.BScale-c.Quantize.BScale)
			}
			if res.BZero != c.Quantize.BZero {
				t.Errorf("bzero: got %.17g want %.17g (diff %g)",
					res.BZero, c.Quantize.BZero, res.BZero-c.Quantize.BZero)
			}
			if res.IMinVal != c.Quantize.IMinVal {
				t.Errorf("iminval: got %d want %d", res.IMinVal, c.Quantize.IMinVal)
			}
			if res.IMaxVal != c.Quantize.IMaxVal {
				t.Errorf("imaxval: got %d want %d", res.IMaxVal, c.Quantize.IMaxVal)
			}

			// Every quantized integer against cfitsio's (exact digest).
			if !c.Out.Present() {
				t.Fatalf("golden has no quantized output but quantize returned 1")
			}
			got := make([]testutil.NullInt64, n)
			for i := range got {
				got[i] = testutil.NullInt64{V: int64(idata[i]), Valid: true}
			}
			c.Out.Check(t, "idata vs fits_quantize_float", got)
		})
	}
}
