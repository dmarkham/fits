package stretch

import (
	"testing"
)

// Tests against the verbatim Siril 1.4.0 oracle (testdata/sref,
// goldens in testdata/golden/sref, generator gen_siril_goldens.py). The
// oracle is the Siril release source compiled unmodified; it reproduces the
// 18 siril-cli goldens bit for bit. Tolerance: none (identical float32
// bits, sha256 over every pixel).

func runSrefCase(t *testing.T, name string, apply func(c srefCase, ch [][]float32, nx, ny int)) {
	t.Helper()
	t.Run(name, func(t *testing.T) {
		c := loadSref(t, name)
		ch, nx, ny := c.input(t)
		apply(c, ch, nx, ny)
		c.check(t, name+" ("+c.Note+")", ch)
	})
}

func f32(v float64) float32 { return float32(v) }

func TestAsinh_SirilOracle(t *testing.T) {
	mono := func(c srefCase, ch [][]float32, nx, ny int) { applySrefOp(t, c, ch, nx, ny) }
	rgb := mono
	runSrefCase(t, "asinh_mono_b50", mono)
	runSrefCase(t, "asinh_mono_b50_off", mono)
	runSrefCase(t, "asinh_rgb_b20_off", rgb)
	runSrefCase(t, "asinh_rgb_b20_human", rgb)
	runSrefCase(t, "asinh_rgb_b200_human_off", rgb)
	runSrefCase(t, "asinh_bright_b500", rgb)
}

// autostretchCase applies `autostretch [-linked] clip target` the way the
// siril command does: linked -> AutostretchLinked, otherwise unlinked over
// all channels of the image (AutostretchUnlinked; for a mono image that is
// Autostretch).
func autostretchCase(c srefCase, ch [][]float32, _, _ int) {
	sc, tb := f32(c.Args[1]), f32(c.Args[2])
	switch {
	case c.Args[0] != 0:
		AutostretchLinked(ch, sc, tb)
	case len(ch) == 1:
		Autostretch(ch[0], sc, tb)
	default:
		AutostretchUnlinked(ch, sc, tb)
	}
}

func TestAutostretch_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"autostretch_bright", "autostretch_bright_linked", "autostretch_bright_mono",
		"autostretch_bright_mono_strong", "autostretch_bright_c1_clip", "autostretch_mixed",
		"autostretch_mixed_linked", "autostretch_m2_le_0", "autostretch_const", "autostretch_allzero",
		"autostretch_zero_g", "autostretch_zero_g_linked", "autostretch_bright_zero_b",
		"autostretch_bright_zero_b_linked",
	} {
		runSrefCase(t, name, autostretchCase)
	}
}

// TestAutostretchParams_SirilOracle pins the (shadows, midtones,
// highlights) Siril computes, including the inverted branch (median >
// 0.5), the MTF edge values for m2 <= 0, the MAD == 0 guard and the
// no-good-pixel fallback (0, 0.2, 1): for a mono all-zero image, and for
// RGB images (unlinked and linked, normal and bright) where one channel is
// all zero, in which case Siril applies the fallback to every channel.
// StatsOK is checked against Siril's return value for the unlinked
// functions; AutostretchLinkedParams has no such output.
func TestAutostretchParams_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"autostretch_bright", "autostretch_bright_linked", "autostretch_bright_mono",
		"autostretch_bright_mono_strong", "autostretch_bright_c1_clip", "autostretch_mixed",
		"autostretch_mixed_linked", "autostretch_m2_le_0", "autostretch_const", "autostretch_allzero",
		"autostretch_zero_g", "autostretch_zero_g_linked", "autostretch_bright_zero_b",
		"autostretch_bright_zero_b_linked", "nan_channel_g", "nan_channel_g_linked", "nan_allnan_mono",
	} {
		t.Run(name, func(t *testing.T) {
			c := loadSref(t, name)
			ch, _, _ := c.input(t)
			sc, tb := f32(c.Args[1]), f32(c.Args[2])
			var got []MTFParams
			switch {
			case c.Args[0] != 0:
				s, m, h := AutostretchLinkedParams(ch, sc, tb)
				got = []MTFParams{{Shadows: s, Midtones: m, Highlights: h}}
			case len(ch) == 1:
				s, m, h := AutostretchParams(ch[0], sc, tb)
				got = AutostretchUnlinkedParams(ch, sc, tb)
				if !same32(got[0].Shadows, s) || !same32(got[0].Midtones, m) || !same32(got[0].Highlights, h) {
					t.Errorf("AutostretchParams %v %v %v != AutostretchUnlinkedParams %+v", s, m, h, got[0])
				}
			default:
				got = AutostretchUnlinkedParams(ch, sc, tb)
			}
			if len(got) != len(c.Params) {
				t.Fatalf("%d parameter sets, want %d", len(got), len(c.Params))
			}
			for i, w := range c.Params {
				g := got[i]
				if !sameParams(g, w) {
					t.Errorf("set %d: got (s=%v m=%v h=%v), Siril (s=%v m=%v h=%v)", i,
						g.Shadows, g.Midtones, g.Highlights, w.Shadows, w.Midtones, w.Highlights)
				}
				if c.Args[0] == 0 && (w.Retval != 0) != !g.StatsOK {
					t.Errorf("set %d: StatsOK=%v, Siril retval %d", i, g.StatsOK, w.Retval)
				}
			}
		})
	}
}

// sameParams: bit-identical shadows, midtones and highlights.
func sameParams(g MTFParams, w asParam) bool {
	return same32(g.Shadows, w.Shadows) && same32(g.Midtones, w.Midtones) && same32(g.Highlights, w.Highlights)
}

func TestMTF_SirilOracle(t *testing.T) {
	apply := func(c srefCase, ch [][]float32, nx, ny int) { applySrefOp(t, c, ch, nx, ny) }
	runSrefCase(t, "mtf_hl_08", apply)
	runSrefCase(t, "mtf_hl_03", apply)
}

func TestLinear_SirilOracle(t *testing.T) {
	runSrefCase(t, "linstretch_bright", func(c srefCase, ch [][]float32, _, _ int) {
		for _, p := range ch {
			Linear(p, f32(c.Args[0]))
		}
	})
}
