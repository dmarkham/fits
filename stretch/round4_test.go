package stretch

import (
	"errors"
	"fmt"
	"math"
	"testing"
)

// Round-4 cases (gen_siril_goldens.py round4_cases), against the verbatim
// Siril 1.4.0 oracle; each case named here is also run through a real
// siril-cli 1.4.0 when it has a command (TestSirilCLIReproduce).
// Tolerance: none (canonical float32 bits), as in oracle_edge_test.go.

// TestCLAHEClipConversion_SirilOracle: OpenCV converts the clip limit with
// static_cast<int>(clipLimit*tileSizeTotal/histSize), which on x86-64 is
// cvttsd2si: 2^31 and above (and +Inf) give INT_MIN, which max(., 1) turns
// into the strongest limit, 1. On a 32x32 image with tile 1 the boundary
// is clip = 2^37; 2^37 - 64 is 2^31 - 1 (no bin reaches it, so no limit);
// a NaN clip fails OpenCV's clipLimit > 0 test (no limit). The generator
// checked each output equals the oracle at clip 1 (limit 1) or clip 0 (no
// limit) as stated.
func TestCLAHEClipConversion_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"clahe_clipconv_below", "clahe_clipconv_at", "clahe_clipconv_1e12", "clahe_clipconv_1e30",
		"clahe_clipconv_inf", "clahe_clipconv_nan", "clahe_clipconv_rgb_at",
	} {
		runEdge(t, name)
	}
}

// TestMTFEdgeParams_SirilOracle: MTF with the finite parameters Siril's
// autostretch can hand it but siril-cli's mtf command refuses (midtones 0
// and 1, shadows >= 1, highlights <= 0, shadows below 0 with highlights
// above 1), and two that siril-cli accepts (shadows above or equal to
// highlights).
func TestMTFEdgeParams_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"mtf_m0", "mtf_m1", "mtf_lo_gt_hi", "mtf_lo_eq_hi", "mtf_lo_ge_1", "mtf_hi_le_0", "mtf_wide",
	} {
		runEdge(t, name)
	}
}

// TestAsinhEdgeParams_SirilOracle: asinh at the edge of the accepted
// range: offset 0.99 (just below the refused 1), beta = 3e38 (beta*x
// overflows to +Inf for input above 1) and beta = 0.3 with -human (the
// GUI's range starts at 0; siril-cli refuses beta < 1).
func TestAsinhEdgeParams_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"asinh_off099_mono", "asinh_off099_rgb", "asinh_bigbeta_mono", "asinh_bigbeta_rgb", "asinh_b03_rgb_human",
	} {
		runEdge(t, name)
	}
}

// TestGHTD0_SirilOracle: D = 0, which Siril's GUI accepts (siril-cli
// refuses it): the output is the clamped input, for each B family, with
// NaN, -0, negative and above-1 pixels.
func TestGHTD0_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_d0_b0", "ght_d0_bpos", "ght_d0_bneg1", "ght_d0_bneg"} {
		runEdge(t, name)
	}
}

// cliArgWidth is the width of golden/siril_cli_args.json's input image
// (gen_golden.py arg_input: 8x4).
const cliArgWidth = 8

func asinhMono(px []float32, beta, offset float32) error { return Asinh(px, beta, offset) }

func mtfApply(px []float32, s, m, h float32) error { return MTF(px, s, m, h) }

// claheMono runs CLAHE; its argument panic (tile < 1) counts as a refusal.
func claheMono(px []float32, ny, tile int, clip float64) (err error) {
	defer func() {
		if r := recover(); r != nil {
			err = fmt.Errorf("CLAHE panicked: %v", r)
		}
	}()
	CLAHE(px, cliArgWidth, ny, tile, clip)
	return nil
}

// TestAsinhMTF_ParamErrors: the parameters Asinh, AsinhRGB, AsinhRGBHuman
// and MTF refuse give an error wrapping ErrAsinhParams / ErrMTFParams and
// leave every channel unchanged, and the edges of the accepted ranges run.
// TestSirilCLIArgs compares the mono decisions with siril-cli; this covers
// the RGB functions and the exact range edges.
func TestAsinhMTF_ParamErrors(t *testing.T) {
	nan, inf := float32(math.NaN()), float32(math.Inf(1))
	below1 := math.Nextafter32(1, 0)
	in := []float32{0, 0.1, 0.5, 0.9, 1}
	asinhFns := map[string]func(r, g, b []float32, beta, offset float32) error{
		"Asinh":         func(r, _, _ []float32, beta, offset float32) error { return Asinh(r, beta, offset) },
		"AsinhRGB":      AsinhRGB,
		"AsinhRGBHuman": AsinhRGBHuman,
	}
	for name, f := range asinhFns {
		for _, a := range [][2]float32{
			{nan, 0}, {inf, 0}, {-inf, 0}, {-1, 0}, {-math.SmallestNonzeroFloat32, 0},
			{10, nan}, {10, inf}, {10, -inf}, {10, 1}, {10, 1.5},
		} {
			r, g, b := append([]float32(nil), in...), append([]float32(nil), in...), append([]float32(nil), in...)
			err := f(r, g, b, a[0], a[1])
			if !errors.Is(err, ErrAsinhParams) {
				t.Errorf("%s(beta %v, offset %v): err = %v, want ErrAsinhParams", name, a[0], a[1], err)
			}
			if sha32(r, g, b) != sha32(in, in, in) {
				t.Errorf("%s(beta %v, offset %v): refused but changed the pixels", name, a[0], a[1])
			}
		}
		for _, a := range [][2]float32{{0, 0}, {0.5, 0}, {math.MaxFloat32, 0}, {10, below1}, {10, -math.MaxFloat32}} {
			r, g, b := append([]float32(nil), in...), append([]float32(nil), in...), append([]float32(nil), in...)
			if err := f(r, g, b, a[0], a[1]); err != nil {
				t.Errorf("%s(beta %v, offset %v): %v, want accepted", name, a[0], a[1], err)
			}
		}
	}
	for _, a := range [][3]float32{
		{nan, 0.5, 1}, {0, nan, 1}, {0, 0.5, nan}, {inf, 0.5, 1}, {-inf, 0.5, 1}, {0, 0.5, inf}, {0, 0.5, -inf},
		{0, -math.SmallestNonzeroFloat32, 1}, {0, math.Nextafter32(1, 2), 1}, {0, inf, 1},
	} {
		px := append([]float32(nil), in...)
		if err := MTF(px, a[0], a[1], a[2]); !errors.Is(err, ErrMTFParams) {
			t.Errorf("MTF(%v): err = %v, want ErrMTFParams", a, err)
		}
		if sha32(px) != sha32(in) {
			t.Errorf("MTF(%v): refused but changed the pixels", a)
		}
	}
	for _, a := range [][3]float32{{0, 0, 1}, {0, 1, 1}, {-math.MaxFloat32, 0.5, math.MaxFloat32}, {2, 0.5, -1}} {
		if err := MTF(append([]float32(nil), in...), a[0], a[1], a[2]); err != nil {
			t.Errorf("MTF(%v): %v, want accepted", a, err)
		}
	}
}
