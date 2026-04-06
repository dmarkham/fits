package stretch

import (
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"testing"
)

const goldenDir = "testdata/golden"

type goldenInput struct {
	Shape  []int     `json:"shape"`
	Pixels []float64 `json:"pixels"`
}

type goldenStretch struct {
	Name    string    `json:"name"`
	Command string    `json:"command"`
	Shape   []int     `json:"shape"`
	Pixels  []float64 `json:"pixels"`
}

func loadInput(t *testing.T) (r, g, b []float32, nx, ny int) {
	t.Helper()
	data, err := os.ReadFile(filepath.Join(goldenDir, "input.json"))
	if err != nil {
		t.Skipf("golden input not found — run gen_golden.py")
	}
	var inp goldenInput
	if err := json.Unmarshal(data, &inp); err != nil {
		t.Fatal(err)
	}
	nch, ny2, nx2 := inp.Shape[0], inp.Shape[1], inp.Shape[2]
	if nch != 3 {
		t.Fatalf("expected 3 channels, got %d", nch)
	}
	chSize := ny2 * nx2
	r = make([]float32, chSize)
	g = make([]float32, chSize)
	b = make([]float32, chSize)
	for i := 0; i < chSize; i++ {
		r[i] = float32(inp.Pixels[i])
		g[i] = float32(inp.Pixels[chSize+i])
		b[i] = float32(inp.Pixels[2*chSize+i])
	}
	return r, g, b, nx2, ny2
}

func loadGolden(t *testing.T, name string) (r, g, b []float32) {
	t.Helper()
	path := filepath.Join(goldenDir, fmt.Sprintf("stretch_%s.json", name))
	data, err := os.ReadFile(path)
	if err != nil {
		t.Skipf("golden %s not found", name)
	}
	var gs goldenStretch
	if err := json.Unmarshal(data, &gs); err != nil {
		t.Fatal(err)
	}
	chSize := gs.Shape[1] * gs.Shape[2]
	r = make([]float32, chSize)
	g = make([]float32, chSize)
	b = make([]float32, chSize)
	for i := 0; i < chSize; i++ {
		r[i] = float32(gs.Pixels[i])
		g[i] = float32(gs.Pixels[chSize+i])
		b[i] = float32(gs.Pixels[2*chSize+i])
	}
	return r, g, b
}

// diffChannels compares two channel arrays and returns max absolute
// difference plus the count of pixels exceeding the given tolerance.
func diffChannels(got, want []float32, tol float32) (maxDiff float32, nbad int) {
	for i := range got {
		d := got[i] - want[i]
		if d < 0 {
			d = -d
		}
		if d > maxDiff {
			maxDiff = d
		}
		if d > tol {
			nbad++
		}
	}
	return maxDiff, nbad
}

// copyChannel returns a copy of a float32 slice.
func copyChannel(src []float32) []float32 {
	dst := make([]float32, len(src))
	copy(dst, src)
	return dst
}

// ---------- Linear stretch test ----------

func TestLinear_Siril(t *testing.T) {
	r, g, b, _, _ := loadInput(t)
	wantR, wantG, wantB := loadGolden(t, "linear")

	Linear(r, 0.05)
	Linear(g, 0.05)
	Linear(b, 0.05)

	const tol = 1e-5
	for _, tc := range []struct {
		name string
		got  []float32
		want []float32
	}{
		{"R", r, wantR},
		{"G", g, wantG},
		{"B", b, wantB},
	} {
		maxDiff, nbad := diffChannels(tc.got, tc.want, tol)
		if nbad > 0 {
			t.Errorf("%s: %d/%d pixels exceed tolerance %g (max diff %g)",
				tc.name, nbad, len(tc.got), tol, maxDiff)
		}
	}
}

// ---------- Asinh stretch tests ----------

func TestAsinh_Siril(t *testing.T) {
	cases := []struct {
		name   string
		beta   float32
		offset float32
	}{
		{"asinh_light", 1, 0},
		{"asinh_medium", 50, 0},
		{"asinh_strong", 500, 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r, g, b, _, _ := loadInput(t)
			wantR, wantG, wantB := loadGolden(t, tc.name)

			AsinhRGB(r, g, b, tc.beta, tc.offset)

			const tol = 1e-4
			for _, ch := range []struct {
				name string
				got  []float32
				want []float32
			}{
				{"R", r, wantR},
				{"G", g, wantG},
				{"B", b, wantB},
			} {
				maxDiff, nbad := diffChannels(ch.got, ch.want, tol)
				if nbad > 0 {
					t.Errorf("%s: %d/%d exceed tol %g (max %g)",
						ch.name, nbad, len(ch.got), tol, maxDiff)
				}
			}
		})
	}
}

// ---------- Autostretch (MTF) tests ----------

func TestAutostretch_Siril(t *testing.T) {
	cases := []struct {
		name       string
		shadowClip float32
		targetBG   float32
	}{
		{"kstars1", -2.8, 0.25},
		{"kstars2", -1.5, 0.25},
		{"kstars3", -1.0, 0.125},
		{"kstars4", -0.5, 0.125},
		{"kstars5", -2.0, 0.125},
		{"kstars6", -4.0, 0.125},
		{"kstars7", -5.5, 0.25},
		{"custom_dark", -6.5, 0.03},
		{"custom_extreme", -7.5, 0.02},
		{"custom_ultra", -9.0, 0.01},
		{"custom_mask", -12.0, 0.005},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r, g, b, _, _ := loadInput(t)
			wantR, wantG, wantB := loadGolden(t, tc.name)

			Autostretch(r, tc.shadowClip, tc.targetBG)
			Autostretch(g, tc.shadowClip, tc.targetBG)
			Autostretch(b, tc.shadowClip, tc.targetBG)

			// Tolerance: mono autostretch matches Siril within 3e-4. RGB
			// has a known channel-statistics discrepancy (up to 0.025
			// on the R channel, which has many zeros from gradient
			// clipping). The core MTF formula is correct; the
			// remaining diff is from per-channel median/MAD computation
			// on 3D FITS data. Tracked for future investigation.
			const tol float32 = 5e-4
			for _, ch := range []struct {
				name string
				got  []float32
				want []float32
			}{
				{"R", r, wantR},
				{"G", g, wantG},
				{"B", b, wantB},
			} {
				maxDiff, nbad := diffChannels(ch.got, ch.want, tol)
				if nbad > 0 {
					t.Errorf("%s: %d/%d exceed tol %g (max %g)",
						ch.name, nbad, len(ch.got), tol, maxDiff)
				}
			}
		})
	}
}

func TestAutostretchLinked_Siril(t *testing.T) {
	r, g, b, _, _ := loadInput(t)
	wantR, wantG, wantB := loadGolden(t, "auto_linked")

	AutostretchLinked([][]float32{r, g, b}, -2.8, 0.25)

	const tol float32 = 5e-4
	for _, ch := range []struct {
		name string
		got  []float32
		want []float32
	}{
		{"R", r, wantR},
		{"G", g, wantG},
		{"B", b, wantB},
	} {
		maxDiff, nbad := diffChannels(ch.got, ch.want, tol)
		if nbad > 0 {
			t.Errorf("%s: %d/%d exceed tol %g (max %g)",
				ch.name, nbad, len(ch.got), tol, maxDiff)
		}
	}
}

// ---------- GHT tests ----------

func TestGHT_Siril(t *testing.T) {
	cases := []struct {
		name string
		D    float64
		B    float64
		SP   float64
	}{
		{"ght_mild", 2.0, 0.0, 0.5},
		{"ght_strong", 5.0, 2.0, 0.3},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r, g, b, _, _ := loadInput(t)
			wantR, wantG, wantB := loadGolden(t, tc.name)

			params := GHTParams{D: tc.D, B: tc.B, SP: tc.SP, HP: 1.0}
			GHT(r, params)
			GHT(g, params)
			GHT(b, params)

			const tol = 1e-4
			for _, ch := range []struct {
				name string
				got  []float32
				want []float32
			}{
				{"R", r, wantR},
				{"G", g, wantG},
				{"B", b, wantB},
			} {
				maxDiff, nbad := diffChannels(ch.got, ch.want, tol)
				if nbad > 0 {
					t.Errorf("%s: %d/%d exceed tol %g (max %g)",
						ch.name, nbad, len(ch.got), tol, maxDiff)
				}
			}
		})
	}
}

// ---------- MTF unit test ----------

func TestMTFValue(t *testing.T) {
	// MTF(0.5, 0.5) = identity → 0.5
	if got := MTFValue(0.5, 0.5); math.Abs(float64(got-0.5)) > 1e-10 {
		t.Errorf("MTF(0.5, 0.5) = %g, want 0.5", got)
	}
	// MTF(0, m) = 0 for any m
	if got := MTFValue(0, 0.3); got != 0 {
		t.Errorf("MTF(0, 0.3) = %g, want 0", got)
	}
	// MTF(1, m) = 1 for any m
	if got := MTFValue(1, 0.3); got != 1 {
		t.Errorf("MTF(1, 0.3) = %g, want 1", got)
	}
}
