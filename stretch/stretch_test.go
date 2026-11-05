package stretch

import (
	"fmt"
	"math"
	"testing"
)

// ---------- siril-cli goldens (testdata/golden/stretch_*.json) ----------
//
// The 18 DefaultPresets commands on inputs/synth.f32. Each golden's sha256
// equals the digest of the original siril-cli output (TestSirilCLIDigests);
// the files are written by the verbatim Siril 1.4.0 oracle, which
// reproduces that output bit for bit. Tolerance: none.

func TestLinear_Siril(t *testing.T) {
	c := loadCLI(t, "linear")
	if c.Command != "linstretch -BP=0.05" {
		t.Fatalf("golden command %q", c.Command)
	}
	ch, _, _ := c.input(t)
	for _, p := range ch {
		Linear(p, 0.05)
	}
	c.check(t, "linear", ch)
}

func TestAsinh_Siril(t *testing.T) {
	for _, tc := range []struct {
		name string
		beta float32
	}{{"asinh_light", 1}, {"asinh_medium", 50}, {"asinh_strong", 500}} {
		t.Run(tc.name, func(t *testing.T) {
			c := loadCLI(t, tc.name)
			ch, _, _ := c.input(t)
			if err := AsinhRGB(ch[0], ch[1], ch[2], tc.beta, 0); err != nil {
				t.Fatal(err)
			}
			c.check(t, tc.name, ch)
		})
	}
}

var autostretchCLI = []struct {
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

// TestAutostretch_Siril: siril-cli `autostretch` on a 3-channel image is
// Siril's unlinked autostretch (AutostretchUnlinked). On this dark image no
// channel is bright, so per-channel mono Autostretch gives the same pixels;
// that is checked too.
func TestAutostretch_Siril(t *testing.T) {
	for _, tc := range autostretchCLI {
		t.Run(tc.name, func(t *testing.T) {
			c := loadCLI(t, tc.name)
			ch, _, _ := c.input(t)
			mono := copyPlanes(ch)
			AutostretchUnlinked(ch, tc.shadowClip, tc.targetBG)
			c.check(t, tc.name, ch)
			for _, p := range mono {
				Autostretch(p, tc.shadowClip, tc.targetBG)
			}
			c.check(t, tc.name+" per-channel", mono)
		})
	}
}

func TestAutostretchLinked_Siril(t *testing.T) {
	c := loadCLI(t, "auto_linked")
	if c.Command != "autostretch -linked" {
		t.Fatalf("golden command %q", c.Command)
	}
	ch, _, _ := c.input(t)
	AutostretchLinked(ch, -2.8, 0.25)
	c.check(t, "auto_linked", ch)
}

// ---------- MTF unit checks (values follow from the formula) ----------

func TestMTFValue(t *testing.T) {
	// MTF(x, 0.5) = (-0.5x)/(-0.5) = x exactly
	for _, x := range []float32{0.1, 0.25, 0.5, 0.9} {
		if got := MTFValue(x, 0.5); got != x {
			t.Errorf("MTFValue(%g, 0.5) = %g, want %g", x, got, x)
		}
	}
	// Siril MTF(): x <= lo -> 0, x >= hi -> 1, for any m.
	for _, m := range []float32{0, 0.3, 1} {
		if got := MTFValue(0, m); got != 0 {
			t.Errorf("MTFValue(0, %g) = %g, want 0", m, got)
		}
		if got := MTFValue(-1, m); got != 0 {
			t.Errorf("MTFValue(-1, %g) = %g, want 0", m, got)
		}
		if got := MTFValue(1, m); got != 1 {
			t.Errorf("MTFValue(1, %g) = %g, want 1", m, got)
		}
	}
	// m = 0 sends every x in (0,1) to 1, m = 1 to 0 (limits of the formula).
	if got := MTFValue(0.3, 0); got != 1 {
		t.Errorf("MTFValue(0.3, 0) = %g, want 1", got)
	}
	if got := MTFValue(0.3, 1); got != 0 {
		t.Errorf("MTFValue(0.3, 1) = %g, want 0", got)
	}
	if got := MTFValue(float32(math.NaN()), 0.3); got == got {
		t.Errorf("MTFValue(NaN) = %g, want NaN (Siril propagates NaN)", got)
	}
}

// TestLinear_Basics: properties stated in Linear's doc (clamp to [0, 1]
// first, 0 and NaN map to 0, bp = 0 is the clamp itself).
func TestLinear_Basics(t *testing.T) {
	nan := float32(math.NaN())
	px := []float32{-1, 0, nan, 0.5, 1, 3}
	Linear(px, 0)
	if want := []float32{0, 0, 0, 0.5, 1, 1}; fmt.Sprint(px) != fmt.Sprint(want) {
		t.Errorf("bp=0: %v, want %v", px, want)
	}
	px = []float32{0.25, 0.5, 1, 2}
	Linear(px, 0.5)
	if want := []float32{0, 0, 1, 1}; fmt.Sprint(px) != fmt.Sprint(want) {
		t.Errorf("bp=0.5: %v, want %v", px, want)
	}
}
