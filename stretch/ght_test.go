package stretch

import (
	"errors"
	"fmt"
	"math"
	"strconv"
	"strings"
	"testing"
)

// GHT tolerance. Siril evaluates its GHT pieces with glibc's logf, powf
// and expf, which are not correctly rounded: on this x86-64/glibc 2.44
// machine they differ from the correctly rounded float32 on ~0.06%
// (powf, expf) and ~0.9% (logf) of arguments (testdata/sref/libm_probe.c;
// glibc picks FMA code paths at run time, so Siril's own output varies
// with the CPU). Go has no float32 libm and rounds the float64 result.
//
// Each GHT golden therefore lists ("libm") the pixels where the verbatim
// oracle and the same Siril code built with double-rounded logf/powf/expf
// (sref_libm, which only locates them) disagree, with Siril's value. Only
// at a listed pixel may the Go value differ, and only if re-running the Go
// GHT on that pixel with the transcendental's float32 result moved by
// exactly one ULP (up or down) reproduces Siril's value bit for bit; such
// a pixel is then replaced by Siril's value and the sha256 over all pixels
// must match. Any other difference fails.

// libmList is a GHT golden's list of libm-sensitive pixels.
type libmList struct {
	Idx  []int  `json:"idx"`
	Bits string `json:"bits"`
}

// ghtLibm accepts, at the listed pixels only, a Go value that a 1-ULP
// different logf/powf/expf result explains, by replacing it in got (flat,
// channel-major) with Siril's value. It returns how many it replaced.
func ghtLibm(t *testing.T, label string, input, got []float32, list libmList, p GHTParams) int {
	t.Helper()
	if len(list.Bits) != 8*len(list.Idx) {
		t.Fatalf("%s: libm list has %d indices and %d hex digits", label, len(list.Idx), len(list.Bits))
	}
	g := prepareGHT(p)
	nudge := func(dir float32) *float32Libm {
		n := func(v float32) float32 { return float32(math.Nextafter32(v, dir)) }
		return &float32Libm{
			logf: func(x float32) float32 { return n(logf(x)) },
			powf: func(x, y float32) float32 { return n(powf(x, y)) },
			expf: func(x float32) float32 { return n(expf(x)) },
		}
	}
	up, down := nudge(float32(math.Inf(1))), nudge(float32(math.Inf(-1)))
	replaced := 0
	for j, k := range list.Idx {
		u, err := strconv.ParseUint(list.Bits[8*j:8*j+8], 16, 32)
		if err != nil {
			t.Fatal(err)
		}
		want := math.Float32frombits(uint32(u))
		if same32(got[k], want) {
			continue
		}
		x := fminf(1, fmaxf(0, input[k]))
		if x != 0 && (g.eval(x, up) == want || g.eval(x, down) == want) {
			got[k] = want
			replaced++
			continue
		}
		t.Errorf("%s: pixel %d (in %v): got %v, Siril %v (%d ULP), not explained by a 1-ULP libm difference",
			label, k, input[k], got[k], want, ulpDiff(got[k], want))
	}
	if replaced > 0 {
		t.Logf("%s: %d/%d pixels match Siril only with glibc's 1-ULP-different logf/powf/expf result (%d listed)",
			label, replaced, len(got), len(list.Idx))
	}
	return replaced
}

// mustGHT applies GHT and fails the test on an error.
func mustGHT(t *testing.T, px []float32, p GHTParams) {
	t.Helper()
	if err := GHT(px, p); err != nil {
		t.Fatalf("GHT(%+v): %v", p, err)
	}
}

// ghtCase is a GHT golden: an oracle case plus its libm list.
type ghtCase struct {
	srefCase
	Libm *libmList `json:"libm"`
}

func loadGHT(t *testing.T, rel string) ghtCase {
	t.Helper()
	var c ghtCase
	readGolden(t, rel, &c)
	if c.Op != "ght" || c.Libm == nil {
		t.Fatalf("%s: not a GHT golden with a libm list", rel)
	}
	return c
}

// checkGHT compares GHT output got (applied to in) with the golden.
func checkGHT(t *testing.T, label string, c ghtCase, in, got [][]float32, p GHTParams) {
	t.Helper()
	var fin, fgot []float32
	for i := range in {
		fin = append(fin, in[i]...)
		fgot = append(fgot, got[i]...)
	}
	ghtLibm(t, label, fin, fgot, *c.Libm, p)
	checkF32(t, label, c.Out, fgot)
}

// runGHT applies GHT with the case's parameters to every channel of its
// input and checks the result.
func runGHT(t *testing.T, label string, c ghtCase) {
	t.Helper()
	in, _, _ := c.input(t)
	got := copyPlanes(in)
	p := ghtParamsOf(c.srefCase)
	for _, ch := range got {
		mustGHT(t, ch, p)
	}
	checkGHT(t, label, c, in, got, p)
}

func TestGHT_Siril(t *testing.T) {
	for _, tc := range []struct {
		name     string
		D, B, SP float64
		cmd      string
	}{
		{"ght_mild", 2.0, 0.0, 0.5, "ght -D=2.0 -B=0.0 -SP=0.5"},
		{"ght_strong", 5.0, 2.0, 0.3, "ght -D=5.0 -B=2.0 -SP=0.3"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			c := loadGHT(t, "stretch_"+tc.name+".json")
			if c.Command != tc.cmd {
				t.Fatalf("golden command %q, want %q", c.Command, tc.cmd)
			}
			in, _, _ := c.input(t)
			got := copyPlanes(in)
			p := GHTParams{D: tc.D, B: tc.B, SP: tc.SP, HP: 1.0}
			for _, ch := range got {
				mustGHT(t, ch, p)
			}
			checkGHT(t, tc.name, c, in, got, p)
		})
	}
}

func TestGHT_SirilOracle(t *testing.T) {
	for _, name := range []string{"ght_b_neg1", "ght_b_neg05", "ght_lp_hp", "ght_b0_lp_hp", "ght_bneg1_lp_hp", "ght_bright_hp"} {
		t.Run(name, func(t *testing.T) {
			c := loadGHT(t, "sref/"+name+".json")
			runGHT(t, name+" ("+c.Note+")", c)
		})
	}
}

// TestGHT_Basics: properties that follow from Siril's code, independent of
// any golden: 0 stays 0, inputs are clamped to [0, 1] first, D = 0 is the
// clamp alone, and HP = 0 means HP = 1.
func TestGHT_Basics(t *testing.T) {
	px := []float32{-0.5, 0, 0.25, 1, 1.5, float32(math.NaN())}
	p := GHTParams{D: 2, B: 1, SP: 0.2, HP: 1}
	got := append([]float32(nil), px...)
	mustGHT(t, got, p)
	ref := []float32{0, 0, 0.25, 1, 1, 0}
	mustGHT(t, ref, p)
	for i := range got {
		if got[i] != ref[i] {
			t.Errorf("pixel %d (%v): %v, want same as clamped input %v -> %v", i, px[i], got[i], ref[i], ref[i])
		}
	}
	if got[0] != 0 || got[1] != 0 || got[5] != 0 {
		t.Errorf("0, negative and NaN must map to 0: %v", got)
	}
	id := append([]float32(nil), px...)
	mustGHT(t, id, GHTParams{D: 0, SP: 0.5})
	if want := []float32{0, 0, 0.25, 1, 1, 0}; fmt.Sprint(id) != fmt.Sprint(want) {
		t.Errorf("D=0: %v, want clamp only %v", id, want)
	}
	a := []float32{0.1, 0.3, 0.7}
	bb := append([]float32(nil), a...)
	mustGHT(t, a, GHTParams{D: 3, B: 0.5, SP: 0.2})
	mustGHT(t, bb, GHTParams{D: 3, B: 0.5, SP: 0.2, HP: 1})
	if fmt.Sprint(a) != fmt.Sprint(bb) {
		t.Errorf("HP=0 should equal HP=1: %v vs %v", a, bb)
	}
}

// ghtOutOfRange are parameter sets siril-cli's process_ght_args refuses
// (core/command.c:2363-2387: D must be in (0, 10], B in [0, 15], SP in
// [0, 1], LP in [0, SP], HP in [SP, 1]; the GUI widens B to [-5, 15],
// gui/uifiles/histogram_dialog.ui:77, and allows D = 0), plus NaN.
var ghtOutOfRange = []GHTParams{
	{D: -0.5, B: 1, SP: 0.2, HP: 1},
	{D: 10.0001, B: 1, SP: 0.2, HP: 1},
	{D: 2, B: -5.0001, SP: 0.2, HP: 1},
	{D: 2, B: 15.0001, SP: 0.2, HP: 1},
	{D: 2, B: 1, SP: -0.1, HP: 1},
	{D: 2, B: 1, SP: 1.1, HP: 1.2},
	{D: 2, B: 1, LP: -0.01, SP: 0.2, HP: 1},
	{D: 2, B: 1, LP: 0.3, SP: 0.2, HP: 1},
	{D: 2, B: 1, SP: 0.5, HP: 0.4},
	{D: 2, B: 1, SP: 0.5, HP: 1.01},
	{D: math.NaN(), B: 1, SP: 0.2, HP: 1},
	{D: 2, B: math.NaN(), SP: 0.2, HP: 1},
	{D: 2, B: 1, LP: math.NaN(), SP: 0.2, HP: 1},
	{D: 2, B: 1, SP: math.NaN(), HP: 1},
	{D: 2, B: 1, SP: 0.2, HP: math.NaN()},
	{D: math.Inf(1), B: 1, SP: 0.2, HP: 1},
	{D: 1, SP: 0.1, HPZero: true},        // -SP=0.1 -HP=0: HP < SP
	{D: 1, SP: 0, HP: 0.5, HPZero: true}, // HPZero with a non-zero HP
}

// TestGHT_RejectsOutOfRange: parameters outside the range Siril accepts are
// an error and leave the pixels untouched (they used to be run through the
// formulas, which siril-cli never does). The edges themselves are accepted
// and follow Siril: D = 10, B = -5 and 15, LP = 0, LP = SP, SP = HP and
// HP = 1 are among TestGHTSweep_SirilOracle's parameter sets, and D = 0
// is the clamp alone (TestGHT_Basics).
func TestGHT_RejectsOutOfRange(t *testing.T) {
	for _, p := range ghtOutOfRange {
		px := []float32{0.1, 0.3, 0.7}
		err := GHT(px, p)
		if !errors.Is(err, ErrGHTParams) {
			t.Errorf("GHT(%+v) error = %v, want ErrGHTParams", p, err)
		}
		if fmt.Sprint(px) != fmt.Sprint([]float32{0.1, 0.3, 0.7}) {
			t.Errorf("GHT(%+v) changed the pixels to %v", p, px)
		}
	}
	// Each refusal names the parameter Siril's own check names (the SP
	// check is otherwise implied by the LP and HP checks).
	for _, c := range []struct {
		p    GHTParams
		name string
	}{
		{GHTParams{D: -1, SP: 0.5}, ": D = "}, {GHTParams{D: 1, B: 16, SP: 0.5}, ": B = "},
		{GHTParams{D: 1, SP: -0.1}, ": SP = "}, {GHTParams{D: 1, SP: 1.1, HP: 1.2}, ": SP = "},
		{GHTParams{D: 1, SP: math.NaN()}, ": SP = "}, {GHTParams{D: 1, LP: 0.3, SP: 0.2}, ": LP = "},
		{GHTParams{D: 1, SP: 0.5, HP: 0.4}, ": HP = "},
	} {
		if err := GHT([]float32{0.5}, c.p); err == nil || !strings.Contains(err.Error(), c.name) {
			t.Errorf("GHT(%+v) error %v, want one naming %q", c.p, err, c.name)
		}
	}
	// the accepted edges
	for _, p := range []GHTParams{{D: 0, SP: 0.5}, {D: 10, B: 15, SP: 1}, {D: 10, B: -5, LP: 0.3, SP: 0.3, HP: 0.3},
		{D: 1, B: 1e-4, SP: 0, HP: 0}, {D: 3, B: -1, LP: 0.2, SP: 0.6, HP: 1}, {D: 1, SP: 0, HPZero: true}} {
		if err := GHT([]float32{0.1, 0.5, 0.9}, p); err != nil {
			t.Errorf("GHT(%+v): %v", p, err)
		}
	}
}
