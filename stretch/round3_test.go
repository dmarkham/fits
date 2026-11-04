package stretch

import (
	"errors"
	"math"
	"strconv"
	"strings"
	"testing"
)

// TestGHTBThresholdF32_SirilOracle: B exactly at siril-cli's threshold,
// (double)1.e-3f = 0.0010000000474974513 (a float32 0.001 widened to
// float64), and its negative. process_ght_args zeroes B only when
// fabs(B) < 1.e-3f, so this B is kept and the power family is used; a
// threshold test with <= would turn it into the exponential branch. The
// positive case is also a siril-cli command (gen_golden.py checks it).
func TestGHTBThresholdF32_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_b_f32thr", "ght_b_negf32thr"} {
		runCase(t, name)
	}
}

// TestGHTHPZero_SirilOracle: HP = 0, which siril-cli accepts with
// LP = SP = 0 (`ght -D=1 -SP=0 -HP=0`). The oracle gives the clamp alone
// (every pixel >= HP) for B = 0, B > 0, B < 0 and B = -1.
func TestGHTHPZero_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_hp0", "ght_hp0_bpos", "ght_hp0_bneg", "ght_hp0_bm1"} {
		runCase(t, name)
	}
}

// TestLinearBP1_SirilOracle: linstretch -BP=1, the upper edge siril-cli
// accepts: (x - 1) / 0 is -Inf, NaN or +Inf, and fmaxf(0, .) makes every
// pixel +0, including x = 1 (0/0 = NaN) and -0.
func TestLinearBP1_SirilOracle(t *testing.T) {
	runCase(t, "linstretch_bp1")
}

// sirilCLIArgs is golden/siril_cli_args.json, written by gen_golden.py
// --write-args from a real siril-cli 1.4.0 (build_siril_cli.sh).
type sirilCLIArgs struct {
	SirilCLI string      `json:"siril_cli"`
	Input    goldArray   `json:"input"`
	Cases    []cliArgRow `json:"cases"`
}

type cliArgRow struct {
	Command   string     `json:"command"`
	Accepted  bool       `json:"accepted"`
	Error     string     `json:"error"`
	Out       *goldArray `json:"out"`
	OutVaries bool       `json:"out_varies"`
}

// cliDivergence lists the commands where the Go function deliberately
// does not follow siril-cli's accept/refuse decision, with the reason.
// Every entry must be a row of the table, so the list cannot go stale.
var cliDivergence = map[string]string{
	// Go refuses what siril-cli runs: NaN parameters pass Siril's range
	// checks (every comparison with NaN is false) and turn the image into
	// NaN, or, for B, read GHTsetup coefficients it never set (the table
	// records out_varies: the pixels change from run to run).
	"autostretch -2.8 nan":         "NaN targetBG: Siril's output is NaN",
	"autostretch nan 0.25":         "NaN shadowsClip: Siril's output is NaN",
	"autostretch -linked nan 0.25": "NaN shadowsClip: Siril's output is NaN",
	"ght -D=2 -B=nan -SP=0.5":      "NaN B: Siril reads uninitialised coefficients",
	"ght -D=2 -SP=nan":             "NaN SP: Siril's output is NaN",
	"ght -D=2 -LP=nan -SP=0.2":     "NaN LP: Siril's output is NaN",
	"ght -D=2 -SP=0.5 -HP=nan":     "NaN HP: Siril's output is NaN",
	// Go runs what siril-cli refuses: Siril's GUI accepts these and runs
	// the same GHT code (TestGHT_Basics, TestGHTSweep_SirilOracle).
	"ght -D=0 -SP=0.5":         "D = 0: the GUI's identity after the clamp",
	"ght -D=2 -B=-0.5 -SP=0.5": "B < 0: the GUI's range is [-5, 15]",

	// asinh. Go refuses a NaN or infinite beta or offset (the output is
	// NaN, or 0 for every pixel) and an offset >= 1 (1 - offset <= 0:
	// division by zero, NaN at x = 1, or an inverted image clamped to 1).
	"asinh nan":     "NaN beta: Siril's output is NaN",
	"asinh inf":     "infinite beta: asinh(inf)/inf is NaN",
	"asinh 10 nan":  "NaN offset: Siril's output is NaN",
	"asinh 10 inf":  "infinite offset: Siril's output is NaN",
	"asinh 10 -inf": "infinite offset: Siril's output is NaN",
	"asinh 10 1":    "offset 1: (x - 1)/0, NaN at x = 1",
	"asinh 10 1.5":  "offset > 1: negative denominator, inverted image",
	// Go runs what siril-cli refuses: Siril's GUI accepts beta in [0, 1)
	// (asinh_dialog.ui adjustment 0..1000) and runs the same code
	// (asinh_beta0_*, asinh_tiny_*_b05, asinh_b03_rgb_human).
	"asinh 0.999": "beta in [0, 1): the GUI's range starts at 0",
	"asinh 0":     "beta = 0: the GUI's range starts at 0",

	// mtf. Go refuses NaN parameters (Siril's output is NaN).
	"mtf nan 0.5 1": "NaN shadows: Siril's output is NaN",
	"mtf 0 nan 1":   "NaN midtones: Siril's output is NaN",
	"mtf 0 0.5 nan": "NaN highlights: Siril's output is NaN",
	// Go runs every finite shadows/highlights and 0 <= midtones <= 1,
	// which Siril's autostretch hands MTF (find_unlinked_midtones_balance:
	// midtones 0 or 1, shadows above 1, highlights below 0) and which
	// mtf_m0, mtf_m1, mtf_lo_ge_1, mtf_hi_le_0, mtf_wide pin.
	"mtf 0 0 1":      "midtones 0: autostretch's midtones when c0 >= median",
	"mtf 0 1 1":      "midtones 1: autostretch's inverted-branch midtones when c1 <= median",
	"mtf -0.1 0.5 1": "shadows < 0: finite, Siril's MTF() takes it",
	"mtf 1 0.5 1":    "shadows >= 1: autostretch's c0 can exceed 1",
	"mtf 0 0.5 0":    "highlights <= 0: autostretch's c1 can be below 0",
	"mtf 0 0.5 1.1":  "highlights > 1: finite, Siril's MTF() takes it",

	// clahe. CLAHE has no error result; clip <= 0 is OpenCV's "no limit"
	// (TestCLAHECore16bit_cv2), as documented on CLAHE.
	"clahe 0 1":  "clip 0: OpenCV's no limit",
	"clahe -1 1": "clip < 0: OpenCV's no limit",
}

// applyCLI runs the Go counterpart of a siril command on a mono image.
func applyCLI(t *testing.T, command string, px []float32) error {
	t.Helper()
	w := strings.Fields(command)
	num := func(s string) float64 {
		v, err := strconv.ParseFloat(s, 64)
		if err != nil {
			t.Fatalf("%q: %v", command, err)
		}
		return v
	}
	switch w[0] {
	case "linstretch":
		return Linear(px, float32(num(strings.TrimPrefix(w[1], "-BP="))))
	case "autostretch":
		linked := len(w) > 1 && w[1] == "-linked"
		if linked {
			w = w[1:]
		}
		sc, tb := float32(-2.8), float32(0.25) // process_autostretch defaults
		if len(w) > 1 {
			sc = float32(num(w[1]))
		}
		if len(w) > 2 {
			tb = float32(num(w[2]))
		}
		if linked {
			return AutostretchLinked([][]float32{px}, sc, tb)
		}
		return AutostretchUnlinked([][]float32{px}, sc, tb)
	case "ght":
		p := GHTParams{HP: 1} // process_ght_args defaults: B = LP = SP = 0, HP = 1
		for _, kv := range w[1:] {
			k, v, _ := strings.Cut(kv, "=")
			x := num(v)
			switch k {
			case "-D":
				p.D = x
			case "-B":
				p.B = x
			case "-LP":
				p.LP = x
			case "-SP":
				p.SP = x
			case "-HP":
				p.HP, p.HPZero = x, x == 0
			default:
				t.Fatalf("%q: unknown option %s", command, k)
			}
		}
		return GHT(px, p)
	case "asinh":
		if len(w) > 1 && w[1] == "-human" { // mono: Siril ignores -human
			w = w[1:]
		}
		beta, offset := float32(num(w[1])), float32(0)
		if len(w) > 2 {
			offset = float32(num(w[2]))
		}
		return asinhMono(px, beta, offset)
	case "mtf":
		return mtfApply(px, float32(num(w[1])), float32(num(w[2])), float32(num(w[3])))
	case "clahe":
		tile, err := strconv.Atoi(w[2])
		if err != nil {
			t.Fatalf("%q: %v", command, err)
		}
		return claheMono(px, len(px)/cliArgWidth, tile, num(w[1]))
	}
	t.Fatalf("no Go counterpart for %q", command)
	return nil
}

// TestSirilCLIArgs: for every command of golden/siril_cli_args.json the Go
// function returns an error exactly when siril-cli refused the command
// (cliDivergence lists the deliberate exceptions), leaves the pixels
// unchanged when it does, and otherwise produces siril-cli's pixels bit
// for bit. The table covers the edges of linstretch's BP, autostretch's
// targetBG and shadowsClip and every GHT range check, NaN and ±Inf
// included.
func TestSirilCLIArgs(t *testing.T) {
	var tab sirilCLIArgs
	readGolden(t, "siril_cli_args.json", &tab)
	if tab.SirilCLI != "siril 1.4.0" || len(tab.Cases) < 50 {
		t.Fatalf("table from %q with %d cases", tab.SirilCLI, len(tab.Cases))
	}
	in := tab.Input.explicit32(t)
	seen := map[string]bool{}
	for _, c := range tab.Cases {
		px := append([]float32(nil), in...)
		err := applyCLI(t, c.Command, px)
		if why, ok := cliDivergence[c.Command]; ok {
			seen[c.Command] = true
			if (err == nil) == c.Accepted {
				t.Errorf("%q: listed as a divergence (%s) but Go agrees with siril-cli (accepted=%v, err=%v)", c.Command, why, c.Accepted, err)
			}
			if err != nil && sha32(px) != sha32(in) {
				t.Errorf("%q: refused but changed the pixels", c.Command)
			}
			continue
		}
		switch {
		case c.Accepted && err != nil:
			t.Errorf("%q: siril-cli runs it, Go: %v", c.Command, err)
		case !c.Accepted && err == nil:
			t.Errorf("%q: siril-cli refuses it (%s), Go ran it: %v", c.Command, c.Error, px)
		case !c.Accepted:
			if sha32(px) != sha32(in) {
				t.Errorf("%q: refused but changed the pixels", c.Command)
			}
		case c.Out == nil:
			t.Errorf("%q: accepted without recorded pixels (out_varies=%v) and not a listed divergence", c.Command, c.OutVaries)
		default:
			checkF32(t, c.Command, *c.Out, px)
		}
	}
	for cmd := range cliDivergence {
		if !seen[cmd] {
			t.Errorf("cliDivergence entry %q is not in the table", cmd)
		}
	}
}

// TestAutostretchParams_PanicOnBadArgs: the *Params functions keep their
// error-free signatures, so they panic (with an error wrapping
// ErrAutostretchParams) on the arguments the Autostretch functions refuse
// (TestSirilCLIArgs), instead of returning parameters for them.
func TestAutostretchParams_PanicOnBadArgs(t *testing.T) {
	px := []float32{0.1, 0.2, 0.3}
	nan := float32(math.NaN())
	for _, a := range [][2]float32{{-2.8, 1.5}, {-2.8, -0.01}, {-2.8, nan}, {nan, 0.25}} {
		for name, f := range map[string]func(){
			"AutostretchParams":         func() { AutostretchParams(px, a[0], a[1]) },
			"AutostretchLinkedParams":   func() { AutostretchLinkedParams([][]float32{px}, a[0], a[1]) },
			"AutostretchUnlinkedParams": func() { AutostretchUnlinkedParams([][]float32{px}, a[0], a[1]) },
		} {
			func() {
				defer func() {
					r := recover()
					if err, ok := r.(error); !ok || !errors.Is(err, ErrAutostretchParams) {
						t.Errorf("%s(sc=%v, tb=%v): recovered %v, want a panic with ErrAutostretchParams", name, a[0], a[1], r)
					}
				}()
				f()
			}()
		}
		if err := Autostretch(append([]float32(nil), px...), a[0], a[1]); !errors.Is(err, ErrAutostretchParams) {
			t.Errorf("Autostretch(sc=%v, tb=%v) = %v, want ErrAutostretchParams", a[0], a[1], err)
		}
	}
	if err := Linear([]float32{0.5}, 1.5); !errors.Is(err, ErrLinearParams) {
		t.Errorf("Linear(bp=1.5) = %v, want ErrLinearParams", err)
	}
}

// TestAsinhTiny_SirilOracle: tiny and subnormal pixels (1e-30, 1e-38,
// 1.4e-45, -1e-30, ...), mono and RGB. asinhf of a tiny argument is the
// argument itself, so k = asinh(beta*x)/(x*asinh(beta)) stays finite; with
// beta = 0.5 the smallest subnormal makes 0/0 and that pixel is NaN in
// Siril too.
func TestAsinhTiny_SirilOracle(t *testing.T) {
	for _, name := range []string{"asinh_tiny_mono", "asinh_tiny_mono_b05", "asinh_tiny_rgb", "asinh_tiny_rgb_b05"} {
		runCase(t, name)
	}
}

// TestCLAHETiny_SirilOracle: images smaller than the tile grid (5x3 with
// an 8x8 grid, 3x2 with 4x4, 1x6 with 4x4). OpenCV pads them with
// BORDER_REFLECT_101 by more than their own size, so the reflection runs
// past index 0, and a 1-pixel-wide image reflects onto its only column.
func TestCLAHETiny_SirilOracle(t *testing.T) {
	for _, name := range []string{"clahe_tiny_mono", "clahe_tiny_rgb", "clahe_1col_mono"} {
		runCase(t, name)
	}
}
