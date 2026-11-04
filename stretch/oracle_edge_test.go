package stretch

import (
	"fmt"
	"math"
	"testing"
)

// Edge cases against the verbatim Siril 1.4.0 oracle (gen_siril_goldens.py
// edge_cases). Tolerance: none (canonical float32 bits; NaN must be NaN
// where Siril gives NaN, and the sign of zero must match).

// applySrefOp runs the Go function that corresponds to an oracle op, the
// way the siril command applies it to an image with these channels.
func applySrefOp(t *testing.T, c srefCase, ch [][]float32, nx, ny int) {
	t.Helper()
	var err error
	defer func() {
		if err != nil {
			t.Fatalf("%s: %v", c.Name, err)
		}
	}()
	switch c.Op {
	case "asinh":
		switch {
		case len(ch) == 1:
			err = Asinh(ch[0], f32(c.Args[0]), f32(c.Args[1]))
		case c.Args[2] != 0:
			err = AsinhRGBHuman(ch[0], ch[1], ch[2], f32(c.Args[0]), f32(c.Args[1]))
		default:
			err = AsinhRGB(ch[0], ch[1], ch[2], f32(c.Args[0]), f32(c.Args[1]))
		}
	case "autostretch":
		autostretchCase(c, ch, nx, ny)
	case "mtf":
		for _, p := range ch {
			if err = MTF(p, f32(c.Args[0]), f32(c.Args[1]), f32(c.Args[2])); err != nil {
				return
			}
		}
	case "linstretch":
		for _, p := range ch {
			if err = Linear(p, f32(c.Args[0])); err != nil {
				return
			}
		}
	case "clahe":
		if c.OpenCV != "5.0.0" {
			t.Fatalf("golden from OpenCV %q, the port follows 5.0.0", c.OpenCV)
		}
		if len(ch) == 1 {
			CLAHE(ch[0], nx, ny, int(c.Args[1]), c.Args[0])
			return
		}
		CLAHERGB(ch[0], ch[1], ch[2], nx, ny, int(c.Args[1]), c.Args[0])
	default:
		t.Fatalf("op %q has no Go counterpart here", c.Op)
	}
}

func ghtParamsOf(c srefCase) GHTParams {
	// The oracle takes HP literally, so HP = 0 there is GHTParams.HPZero.
	return GHTParams{D: c.Args[0], B: c.Args[1], LP: c.Args[2], SP: c.Args[3], HP: c.Args[4], HPZero: c.Args[4] == 0}
}

func runEdge(t *testing.T, name string) {
	t.Run(name, func(t *testing.T) {
		c := loadSref(t, name)
		if c.Op == "ght" {
			runGHT(t, name+" ("+c.Note+")", loadGHT(t, "sref/"+name+".json"))
			return
		}
		ch, nx, ny := c.input(t)
		applySrefOp(t, c, ch, nx, ny)
		c.check(t, name+" ("+c.Note+")", ch)
	})
}

// TestNaN_SirilOracle: NaN pixels at different positions in each channel
// plus one pixel NaN in all channels. Pins: autostretch statistics exclude
// NaN (mono, unlinked, linked, inverted branch); asinh propagates NaN
// through Siril's max() macro (a NaN in one channel makes the pixel's
// luminance, and so every channel of that pixel, NaN); linstretch and GHT
// turn NaN into 0; MTF propagates NaN; CLAHE follows OpenCV's SSE min/max
// in the Lab conversions and saturate_cast. A channel that is all NaN has
// no good pixel, so Siril falls back to (0, 0.2, 1) for every channel
// (parameters checked in TestAutostretchParams_SirilOracle).
func TestNaN_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"nan_autostretch", "nan_autostretch_linked", "nan_autostretch_mono", "nan_autostretch_bright",
		"nan_asinh_rgb", "nan_asinh_rgb_human", "nan_asinh_mono", "nan_linstretch", "nan_mtf",
		"nan_ght", "nan_clahe_rgb", "nan_clahe_mono",
		"nan_channel_g", "nan_channel_g_linked", "nan_allnan_mono",
	} {
		runEdge(t, name)
	}
}

// TestNaN_AutostretchParams_SirilOracle: the parameters themselves, so a
// NaN that leaked into the median or MAD is reported as such.
func TestNaN_AutostretchParams_SirilOracle(t *testing.T) {
	for _, name := range []string{"nan_autostretch", "nan_autostretch_linked", "nan_autostretch_mono",
		"nan_autostretch_bright", "autostretch_median_049", "autostretch_median_0495",
		"autostretch_median_0499", "autostretch_median_04999", "autostretch_median_05"} {
		t.Run(name, func(t *testing.T) {
			c := loadSref(t, name)
			ch, _, _ := c.input(t)
			sc, tb := f32(c.Args[1]), f32(c.Args[2])
			var got []MTFParams
			if c.Args[0] != 0 {
				s, m, h := AutostretchLinkedParams(ch, sc, tb)
				got = []MTFParams{{Shadows: s, Midtones: m, Highlights: h, StatsOK: true}}
			} else {
				got = AutostretchUnlinkedParams(ch, sc, tb)
			}
			if len(got) != len(c.Params) {
				t.Fatalf("%d parameter sets, want %d", len(got), len(c.Params))
			}
			for i, w := range c.Params {
				g := got[i]
				if w.Retval != 0 || !g.StatsOK {
					t.Errorf("set %d: Siril retval %d, StatsOK %v; both should succeed", i, w.Retval, g.StatsOK)
				}
				if !sameParams(g, w) {
					t.Errorf("set %d: got (s=%v m=%v h=%v), Siril (s=%v m=%v h=%v)", i,
						g.Shadows, g.Midtones, g.Highlights, w.Shadows, w.Midtones, w.Highlights)
				}
			}
		})
	}
}

// TestAutostretchInvertedBoundary_SirilOracle: mono images whose median
// sits just below and just above Siril's inverted-branch edge (median >
// 0.5 on Siril's histogram median). The generator checked which branch
// the oracle took for each.
func TestAutostretchInvertedBoundary_SirilOracle(t *testing.T) {
	for _, name := range []string{"autostretch_median_049", "autostretch_median_0495",
		"autostretch_median_0499", "autostretch_median_04999", "autostretch_median_05"} {
		runEdge(t, name)
	}
}

// TestGHTEdges_SirilOracle: output above 1 (Siril does not clamp GHT's
// output; the generator required at least one value > 1), small D (where
// float32(expm1(D)) differs from float32(exp(D))-1), B on both sides of
// the |B| < 0.001 threshold, B = -0.35 (where (B-1)/B and 1-1/B round
// differently in float32) and B = -1 with many pixels below LP.
func TestGHTEdges_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_out_above_1", "ght_out_above_1_b0", "ght_small_d_005",
		"ght_small_d_001", "ght_b_tiny_pos", "ght_b_tiny_neg", "ght_b_above_thresh",
		"ght_bneg035_lp_hp", "ght_bneg1_big_lp"} {
		runEdge(t, name)
	}
}

// TestGHTOutputAbove1_SirilOracle states the claim directly: Siril's
// output has values above 1, and GHT returns those same values.
func TestGHTOutputAbove1_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_out_above_1", "ght_out_above_1_b0"} {
		c := loadGHT(t, "sref/"+name+".json")
		in, _, _ := c.input(t)
		ch := copyPlanes(in)
		for _, p := range ch {
			mustGHT(t, p, ghtParamsOf(c.srefCase))
		}
		checkGHT(t, name, c, in, ch, ghtParamsOf(c.srefCase))
		// The sampled pixels include the image maximum; Siril's is > 1.
		above := 0
		words := c.Out.words(t, 8)
		for i, k := range c.Out.Idx {
			if w := math.Float32frombits(words[i]); w > 1 {
				above++
				if g := ch[k/len(ch[0])][k%len(ch[0])]; g != w {
					t.Errorf("%s pixel %d: GHT %v, Siril %v (> 1)", name, k, g, w)
				}
			}
		}
		if above == 0 {
			t.Fatalf("%s: golden sample has no value above 1", name)
		}
	}
}

// TestAsinhEdges_SirilOracle: beta = 0 (k = 1; siril-cli rejects beta < 1,
// the oracle calls command_asinh directly) for mono and RGB, and mono
// input above 1, where min(1, .) clamps.
func TestAsinhEdges_SirilOracle(t *testing.T) {
	for _, name := range []string{"asinh_beta0_mono", "asinh_beta0_rgb", "asinh_mono_above_1"} {
		runEdge(t, name)
	}
	c := loadSref(t, "asinh_mono_above_1")
	in, _, _ := c.input(t)
	n := 0
	for _, v := range in[0] {
		if v > 1 {
			n++
		}
	}
	if n == 0 {
		t.Fatal("asinh_mono_above_1: no input above 1")
	}
	t.Logf("asinh_mono_above_1: %d inputs above 1", n)
}

// TestLabNaN_cv2: OpenCV's float BGR2Lab and Lab2BGR with NaN in every
// combination of components, in rows of width 1, 13 and 64 (scalar loop,
// SIMD body and tail). Exact, NaN must match NaN.
func TestLabNaN_cv2(t *testing.T) {
	var g struct {
		Forward []struct {
			Width int       `json:"width"`
			Input goldInput `json:"input"`
			Out   goldArray `json:"out"`
		} `json:"forward"`
		Inverse []struct {
			Width int       `json:"width"`
			Input goldInput `json:"input"`
			Out   goldArray `json:"out"`
		} `json:"inverse"`
	}
	readGolden(t, "nan_cv2.json", &g)
	if len(g.Forward) == 0 || len(g.Inverse) == 0 {
		t.Fatal("nan_cv2.json has no vectors")
	}
	for _, v := range g.Forward {
		checkF32(t, fmt.Sprintf("BGR2Lab NaN width %d", v.Width), v.Out, convertRows(v.Input.flat(t), v.Width, bgr2labRow))
	}
	for _, v := range g.Inverse {
		checkF32(t, fmt.Sprintf("Lab2BGR NaN width %d", v.Width), v.Out, convertRows(v.Input.flat(t), v.Width, lab2bgrRow))
	}
}

// TestCvRound_cv2: cvRound (saturate_cast<int>(float)) against OpenCV:
// half to even, INT_MIN for NaN, +-inf and values outside int32.
func TestCvRound_cv2(t *testing.T) {
	var g struct {
		RoundIn  goldInput `json:"round_in"`
		RoundOut []int64   `json:"round_out"`
	}
	readGolden(t, "nan_cv2.json", &g)
	in := g.RoundIn.flat(t)
	if len(in) == 0 || len(in) != len(g.RoundOut) {
		t.Fatalf("%d inputs, %d outputs", len(in), len(g.RoundOut))
	}
	for i, v := range in {
		if got := int64(cvRound(v)); got != g.RoundOut[i] {
			t.Errorf("cvRound(%v) = %d, OpenCV %d", v, got, g.RoundOut[i])
		}
	}
}
