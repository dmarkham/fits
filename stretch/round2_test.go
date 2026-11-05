package stretch

import (
	"fmt"
	"math"
	"strconv"
	"testing"
)

// Inputs no earlier golden had (gen_siril_goldens.py round2_cases and the
// sweeps): negative pixels and -0, +-Inf, a median of exactly 0.5,
// B = +-0.001 exactly, random GHT parameter sets and random autostretch
// images. Truth: the verbatim Siril 1.4.0 oracle. Tolerance: none, except
// the listed libm-sensitive GHT pixels (ghtLibm).

// asParamsOf runs the Go parameter function matching the case's
// autostretch arguments (linked, clip, target) on ch.
func asParamsOf(c srefCase, ch [][]float32) []MTFParams {
	sc, tb := f32(c.Args[1]), f32(c.Args[2])
	if c.Args[0] != 0 {
		s, m, h := AutostretchLinkedParams(ch, sc, tb)
		return []MTFParams{{Shadows: s, Midtones: m, Highlights: h, StatsOK: true}}
	}
	return AutostretchUnlinkedParams(ch, sc, tb)
}

// checkASParams compares Go's autostretch parameters with Siril's, bit for
// bit; for unlinked cases StatsOK must match Siril's return value.
func checkASParams(t *testing.T, label string, linked bool, got []MTFParams, want []asParam) {
	t.Helper()
	if len(got) != len(want) {
		t.Fatalf("%s: %d parameter sets, want %d", label, len(got), len(want))
	}
	for i, w := range want {
		g := got[i]
		if !sameParams(g, w) {
			t.Errorf("%s set %d: got (s=%v m=%v h=%v), Siril (s=%v m=%v h=%v)", label, i,
				g.Shadows, g.Midtones, g.Highlights, w.Shadows, w.Midtones, w.Highlights)
		}
		if !linked && (w.Retval != 0) != !g.StatsOK {
			t.Errorf("%s set %d: StatsOK=%v, Siril retval %d", label, i, g.StatsOK, w.Retval)
		}
	}
}

// runCase runs an oracle case (pixels, and parameters for autostretch).
func runCase(t *testing.T, name string) {
	t.Helper()
	t.Run(name, func(t *testing.T) {
		c := loadSref(t, name)
		if c.Op == "autostretch" {
			ch, _, _ := c.input(t)
			checkASParams(t, name, c.Args[0] != 0, asParamsOf(c, ch), c.Params)
		}
		if c.Op == "ght" {
			runGHT(t, name+" ("+c.Note+")", loadGHT(t, "sref/"+name+".json"))
			return
		}
		ch, nx, ny := c.input(t)
		applySrefOp(t, c, ch, nx, ny)
		c.check(t, name+" ("+c.Note+")", ch)
	})
}

// TestNegative_SirilOracle: about a quarter of the pixels negative (and
// some -0) for every operation. Autostretch statistics exclude only 0 and
// NaN, so negative pixels count (mono, unlinked, linked, inverted branch);
// asinh, MTF, linstretch, GHT and CLAHE follow Siril on them, and the sign
// of every output zero matches (linstretch -BP=0 maps -0 to +0).
func TestNegative_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"neg_autostretch", "neg_autostretch_linked", "neg_autostretch_mono", "neg_autostretch_bright",
		"neg_autostretch_bright_linked", "neg_asinh_rgb", "neg_asinh_rgb_off", "neg_asinh_human_off",
		"neg_asinh_mono", "neg_linstretch", "neg_linstretch_bp0", "neg_mtf", "neg_ght", "neg_clahe_rgb",
		"neg_clahe_mono",
	} {
		runCase(t, name)
	}
	// The input really has negatives (the generator requires >= 20%).
	c := loadSref(t, "neg_autostretch")
	ch, _, _ := c.input(t)
	n, negz := 0, 0
	for _, p := range ch {
		for _, v := range p {
			if v < 0 {
				n++
			}
			if math.Float32bits(v) == 0x80000000 {
				negz++
			}
		}
	}
	if n*5 < len(ch)*len(ch[0]) || negz == 0 {
		t.Fatalf("neg input: %d negative, %d -0 of %d", n, negz, len(ch)*len(ch[0]))
	}
}

// TestLinear_NegativeZero states the sign-of-zero rule directly: Siril's
// linstretch with BP = 0 returns +0 for a -0 input (x == 0 -> 0.0f), and so
// does Linear.
func TestLinear_NegativeZero(t *testing.T) {
	c := loadSref(t, "neg_linstretch_bp0")
	ch, _, _ := c.input(t)
	for _, p := range ch {
		Linear(p, 0)
	}
	c.check(t, "neg_linstretch_bp0", ch)
	in, _, _ := c.input(t)
	seen := 0
	for k := range in {
		for i, v := range in[k] {
			if math.Float32bits(v) == 0x80000000 {
				seen++
				if math.Float32bits(ch[k][i]) != 0 {
					t.Errorf("ch%d pixel %d: Linear(-0) = %v (%08x), Siril +0", k, i, ch[k][i], math.Float32bits(ch[k][i]))
				}
			}
		}
	}
	if seen == 0 {
		t.Fatal("input has no -0 pixel")
	}
}

// TestAutostretchMedianHalf_SirilOracle: Siril takes the inverted branch
// only for median > 0.5. A constant 0.5 image (mono, RGB, linked) has a
// median of exactly 0.5 and takes the normal branch; one float32 step
// above takes the inverted one (the generator checked the branch Siril
// took).
func TestAutostretchMedianHalf_SirilOracle(t *testing.T) {
	for _, name := range []string{"autostretch_const_05", "autostretch_const_05_rgb", "autostretch_const_05_linked",
		"autostretch_const_05_up", "autostretch_const_05_up_linked"} {
		runCase(t, name)
	}
}

// TestGHTBThreshold_SirilOracle: B = 0.001 and -0.001 exactly. siril-cli
// compares fabs(B) < 1.e-3f with B in double (command.c:2328), so the user
// value 0.001 becomes 0 and the exponential branch is used.
func TestGHTBThreshold_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_b_0001", "ght_b_neg0001"} {
		runCase(t, name)
	}
}

// infCase is an inf_autostretch* golden: Params and Out are the truth for
// the Go rule (Siril on the image with the Inf pixels set to 0, Siril's
// MTF on the Inf image); SirilInfParams is what Siril itself computes on
// the Inf image.
type infCase struct {
	srefCase
	SirilInfParams []asParam `json:"siril_inf_params"`
}

// TestInf_SirilOracle: +-Inf pixels. Siril's autostretch statistics on a
// channel with an infinite pixel depend on undefined behaviour
// (findMinMaxPercentile casts an infinite or NaN histogram position to
// uint16) and come out NaN or meaningless (recorded in the golden). The
// Go autostretch excludes +-Inf from the statistics like 0 and NaN: its
// parameters must equal Siril's for the same image with the Inf pixels
// set to 0, and its pixels Siril's MTF with those parameters applied to
// the image with the Inf pixels (+Inf -> 1, -Inf -> 0). Every other
// operation follows Siril on Inf input directly.
func TestInf_SirilOracle(t *testing.T) {
	for _, name := range []string{"inf_autostretch", "inf_autostretch_linked", "inf_autostretch_mono",
		"inf_autostretch_mono_neg", "inf_autostretch_bright"} {
		t.Run(name, func(t *testing.T) {
			var c infCase
			readGolden(t, "sref/"+name+".json", &c)
			ch, nx, ny := c.input(t)
			got := asParamsOf(c.srefCase, ch)
			checkASParams(t, name, c.Args[0] != 0, got, c.Params)
			for i, p := range got {
				if p.Shadows != p.Shadows || p.Midtones != p.Midtones || p.Highlights != p.Highlights {
					t.Errorf("set %d: NaN parameter %+v", i, p)
				}
			}
			differs := false
			for i, w := range c.SirilInfParams {
				differs = differs || !sameParams(got[i], w)
			}
			if !differs {
				t.Errorf("Siril's own parameters on the Inf image %+v equal the Inf-excluded ones; the documented divergence is gone", c.SirilInfParams)
			}
			applySrefOp(t, c.srefCase, ch, nx, ny)
			c.check(t, name, ch)
		})
	}
	for _, name := range []string{"inf_asinh_rgb", "inf_asinh_mono", "inf_asinh_rgb_off", "inf_linstretch",
		"inf_mtf", "inf_ght", "inf_clahe_rgb", "inf_clahe_mono"} {
		runCase(t, name)
	}
}

// ghtSweep is one golden/sweep/ght_<family>.json.
type ghtSweep struct {
	Name   string   `json:"name"`
	Fields []string `json:"coef_fields"`
	Cases  []struct {
		Args  []float64  `json:"args"`
		Coef  string     `json:"coef"`
		Input *goldInput `json:"input"`
		Out   *goldArray `json:"out"`
		Libm  *libmList  `json:"libm"`
	} `json:"cases"`
}

// TestGHTSweep_SirilOracle: 256 random parameter sets per B family (B =
// -1, B < 0, B = 0 incl. |B| < 0.001, B > 0) inside the accepted range,
// including D = 10, B = -5, B = 15, LP = 0, LP = SP, SP = HP and HP = 1.
// The GHTsetup coefficients must equal Siril's bit for bit (every float /
// double promotion of filters/ght.c), and on 64 sets per family the
// pixels 0, -0, negative, NaN, exactly LP, SP, HP and one float32 step
// either side, 1, 1.5 and random values must equal Siril's GHT (branch
// comparisons `<`, the "0 stays 0" rule).
func TestGHTSweep_SirilOracle(t *testing.T) {
	want := []string{"b1", "a2", "b2", "c2", "d2", "e2", "a3", "b3", "c3", "d3", "e3", "a4", "b4"}
	total, withPixels := 0, 0
	for _, fam := range []string{"bneg1", "bneg", "b0", "bpos"} {
		var g ghtSweep
		readGolden(t, "sweep/ght_"+fam+".json", &g)
		if fmt.Sprint(g.Fields) != fmt.Sprint(want) || len(g.Cases) != 256 {
			t.Fatalf("%s: fields %v, %d cases", fam, g.Fields, len(g.Cases))
		}
		for i, c := range g.Cases {
			label := fmt.Sprintf("%s set %d %v", fam, i, c.Args)
			p := GHTParams{D: c.Args[0], B: c.Args[1], LP: c.Args[2], SP: c.Args[3], HP: c.Args[4]}
			co := prepareGHT(p).c
			got := []float32{co.b1, co.a2, co.b2, co.c2, co.d2, co.e2, co.a3, co.b3, co.c3, co.d3, co.e3, co.a4, co.b4}
			if len(c.Coef) != 8*len(got) {
				t.Fatalf("%s: coef has %d hex digits", label, len(c.Coef))
			}
			for k, v := range got {
				u, _ := strconv.ParseUint(c.Coef[8*k:8*k+8], 16, 32)
				if w := math.Float32frombits(uint32(u)); !same32(v, w) {
					t.Errorf("%s: %s = %v (%08x), Siril GHTsetup %v (%08x)", label, want[k], v, math.Float32bits(v), w, u)
				}
			}
			total++
			if c.Input == nil {
				continue
			}
			withPixels++
			in := c.Input.flat(t)
			out := append([]float32(nil), in...)
			mustGHT(t, out, p)
			ghtLibm(t, label, in, out, *c.Libm, p)
			checkF32(t, label, *c.Out, out)
		}
	}
	if total != 1024 || withPixels != 256 {
		t.Fatalf("%d parameter sets, %d with pixels", total, withPixels)
	}
}

// TestAutostretchSweep_SirilOracle: Siril's unlinked and linked
// autostretch parameters on 96 scaled and offset copies of the three base
// images with random shadows clip and target, both branches, with the
// linked shadows above 0 or highlights below 1 in many of them so the
// float32 averaging c0/n and c1/n is visible.
func TestAutostretchSweep_SirilOracle(t *testing.T) {
	var g struct {
		Cases []struct {
			Input    goldInput `json:"input"`
			Clip     float64   `json:"clip"`
			Target   float64   `json:"target"`
			Unlinked []asParam `json:"unlinked"`
			Linked   []asParam `json:"linked"`
		} `json:"cases"`
	}
	readGolden(t, "sweep/autostretch.json", &g)
	if len(g.Cases) != 96 {
		t.Fatalf("%d cases, want 96", len(g.Cases))
	}
	for i, c := range g.Cases {
		ch, _, _ := c.Input.f32(t)
		sc, tb := f32(c.Clip), f32(c.Target)
		label := fmt.Sprintf("case %d (%s clip %v target %v)", i, opsString(c.Input.Ops), sc, tb)
		checkASParams(t, label+" unlinked", false, AutostretchUnlinkedParams(ch, sc, tb), c.Unlinked)
		s, m, h := AutostretchLinkedParams(ch, sc, tb)
		checkASParams(t, label+" linked", true, []MTFParams{{Shadows: s, Midtones: m, Highlights: h}}, c.Linked)
	}
}
