package compress

import (
	"encoding/binary"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"
)

// C reference validation. These tests load the golden fixtures produced
// by testdata/cref/cref.c (which links against libcfitsio) and diff the
// Go port output against the real cfitsio implementation.
//
// The fixtures are checked into the repo so CI doesn't need a C
// compiler or libcfitsio. To regenerate them locally:
//
//	cd compress/testdata/cref && make clean && make goldens
//
// If the golden directory is missing these tests are skipped rather
// than failing, so a fresh clone without the fixtures is still green.

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
}

// loadCrefCase reads the input float32 array and the golden JSON metadata
// for a test case. Returns nil slice + nil error + true to skip if the
// fixture is missing.
func loadCrefCase(t *testing.T, name string) ([]float32, *crefCase, bool) {
	t.Helper()
	jsonPath := filepath.Join(cgoldenDir, name+".out.json")
	binPath := filepath.Join(cgoldenDir, name+".in.bin")

	jsonData, err := os.ReadFile(jsonPath)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, nil, true
		}
		t.Fatalf("read %s: %v", jsonPath, err)
	}
	var c crefCase
	if err := json.Unmarshal(jsonData, &c); err != nil {
		t.Fatalf("parse %s: %v", jsonPath, err)
	}

	binData, err := os.ReadFile(binPath)
	if err != nil {
		t.Fatalf("read %s: %v", binPath, err)
	}
	n := c.NX * c.NY
	if int64(len(binData)) != n*4 {
		t.Fatalf("%s: expected %d float32 bytes, got %d", binPath, n*4, len(binData))
	}
	arr := make([]float32, n)
	for i := int64(0); i < n; i++ {
		arr[i] = math.Float32frombits(binary.LittleEndian.Uint32(binData[i*4:]))
	}
	return arr, &c, false
}

// loadCrefInt32Out reads the quantized int32 array for a case. Empty
// file = cfitsio returned 0 (not quantized).
func loadCrefInt32Out(t *testing.T, name string, n int64) []int32 {
	t.Helper()
	path := filepath.Join(cgoldenDir, name+".out.bin")
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	if len(data) == 0 {
		return nil
	}
	if int64(len(data)) != n*4 {
		t.Fatalf("%s: expected %d int32 bytes, got %d", path, n*4, len(data))
	}
	out := make([]int32, n)
	for i := int64(0); i < n; i++ {
		out[i] = int32(binary.LittleEndian.Uint32(data[i*4:]))
	}
	return out
}

// crefAllCases is the canonical list maintained in sync with cref.c.
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
	anyLoaded := false
	for _, name := range crefAllCases {
		t.Run(name, func(t *testing.T) {
			arr, c, skip := loadCrefCase(t, name)
			if skip {
				t.Skipf("fixture not found — run `make -C compress/testdata/cref goldens`")
				return
			}
			anyLoaded = true

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
	if !anyLoaded {
		t.Skip("no cref fixtures available — regenerate via `make -C compress/testdata/cref goldens`")
	}
}

// TestQuantizeFloat32_CRef validates fits_quantize_float bit-exactly
// against cfitsio. Every non-constant case should produce identical
// int32 output, identical bscale/bzero, and identical iminval/imaxval.
// Constant / tiny cases should return quantized=false and leave idata
// unused, matching the C quantize_return=0 marker (empty .out.bin).
func TestQuantizeFloat32_CRef(t *testing.T) {
	anyLoaded := false
	for _, name := range crefAllCases {
		t.Run(name, func(t *testing.T) {
			arr, c, skip := loadCrefCase(t, name)
			if skip {
				t.Skipf("fixture not found — run `make -C compress/testdata/cref goldens`")
				return
			}
			anyLoaded = true

			nullcheck := c.Nullcheck != 0
			var nullvalue float32
			if nullcheck {
				nullvalue = -1e30
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
				// Nothing more to diff for the not-quantized cases.
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

			// Byte-for-byte compare against the golden int32 stream.
			want := loadCrefInt32Out(t, name, n)
			if want == nil {
				t.Fatalf("empty golden .out.bin but quantize returned 1")
			}
			mismatches := 0
			firstMis := -1
			for i := int64(0); i < n; i++ {
				if idata[i] != want[i] {
					if mismatches == 0 {
						firstMis = int(i)
					}
					mismatches++
				}
			}
			if mismatches > 0 {
				t.Errorf("%d/%d pixels differ (first at %d: got %d want %d)",
					mismatches, n, firstMis, idata[firstMis], want[firstMis])
			}
		})
	}
	if !anyLoaded {
		t.Skip("no cref fixtures available — regenerate via `make -C compress/testdata/cref goldens`")
	}
}
