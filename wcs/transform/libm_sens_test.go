package transform

import (
	"math"
	"testing"
)

// Libm conditioning and boundary ambiguity.
//
// The port calls Go's math functions where wcslib calls glibc's; both are
// within a couple of ULP of the exact value (TestLibmAccuracy,
// TestAtan2Accuracy) but they are not the same functions, so an output
// can differ from wcslib's wherever the projection magnifies the last bit
// of a sin, cos, atan2, ... (1 - sin(theta) next to a pole, log(cos xi)
// for small xi, asin of a value next to 1, a cancellation next to a cube
// face centre). Two allowances, both measured on our own code and both
// counted by the tests that use them:
//
//   - Values: libmSensitivity runs the computation again with every libm
//     result moved by libmSensULP ULP (all up, all down, and alternating
//     in both phases, through the lm* hook in wcstrig.go) and returns the
//     largest change of the output. A row whose difference from wcslib
//     exceeds the budget is accepted when the difference is within the
//     budget plus that sensitivity: then a few ULP of libm explain it.
//     libmSensULP = 4 covers our up-to-2-ULP error plus glibc's 1 ULP.
//   - Decisions (accept/reject): a mismatch is accepted only at a
//     boundary-ambiguous point, where our own decision flips under a nudge
//     of each input by up to nudgeULP ULP or under the libm perturbation.
//
// A well-conditioned row has a sensitivity of a few ULP, far below every
// budget, so the allowance changes nothing there.
const (
	libmSensULP = 4
	nudgeULP    = 4
)

// withLibmPerturbed runs f with the libm results moved by k ULP; alt
// alternates the sign from call to call, starting with -k when phase is
// 0 and +k when it is 1.
func withLibmPerturbed(k int64, alt bool, phase int64, f func()) {
	libmPerturbULP, libmPerturbAlt, libmPerturbN = k, alt, phase
	defer func() { libmPerturbULP, libmPerturbAlt, libmPerturbN = 0, false, 0 }()
	f()
}

// libmRuns evaluates eval under the four perturbation patterns.
func libmRuns(eval func() (float64, float64, bool)) (outs [4][2]float64, oks [4]bool) {
	pats := [4]struct {
		k     int64
		alt   bool
		phase int64
	}{{libmSensULP, false, 0}, {-libmSensULP, false, 0}, {libmSensULP, true, 0}, {libmSensULP, true, 1}}
	for i, p := range pats {
		withLibmPerturbed(p.k, p.alt, p.phase, func() {
			outs[i][0], outs[i][1], oks[i] = eval()
		})
	}
	return outs, oks
}

// libmSensitivity returns the largest dist between eval's unperturbed
// output and its output under each perturbation pattern (+Inf if a
// perturbed run is rejected while the unperturbed one is not).
func libmSensitivity(eval func() (float64, float64, bool), dist func(a0, a1, b0, b1 float64) float64) float64 {
	o0, o1, ok := eval()
	if !ok {
		return math.Inf(1)
	}
	outs, oks := libmRuns(eval)
	s := 0.0
	for i := range outs {
		if !oks[i] {
			return math.Inf(1)
		}
		s = math.Max(s, dist(o0, o1, outs[i][0], outs[i][1]))
	}
	return s
}

// nudgedF moves x by k ULP (k < 0 toward -Inf).
func nudgedF(x float64, k int) float64 {
	for ; k > 0; k-- {
		x = math.Nextafter(x, math.Inf(1))
	}
	for ; k < 0; k++ {
		x = math.Nextafter(x, math.Inf(-1))
	}
	return x
}

// decisionAmbiguous reports whether eval's accept/reject decision at
// (in0, in1) changes when either input is nudged by up to nudgeULP ULP or
// the libm results are perturbed.
func decisionAmbiguous(eval func(in0, in1 float64) (float64, float64, bool), in0, in1 float64) bool {
	_, _, ok := eval(in0, in1)
	for d0 := -nudgeULP; d0 <= nudgeULP; d0++ {
		for d1 := -nudgeULP; d1 <= nudgeULP; d1++ {
			if _, _, k := eval(nudgedF(in0, d0), nudgedF(in1, d1)); k != ok {
				return true
			}
		}
	}
	_, oks := libmRuns(func() (float64, float64, bool) { return eval(in0, in1) })
	for _, k := range oks {
		if k != ok {
			return true
		}
	}
	return false
}

// condTally counts rows accepted through either allowance.
type condTally struct {
	rows, conditioned, ambiguous int
	worstRatio                   float64 // largest (difference - budget) / sensitivity accepted
}

func (c *condTally) log(t testing.TB, where string) {
	t.Helper()
	if c.conditioned+c.ambiguous > 0 {
		t.Logf("%s: %d rows; %d within budget only with their libm sensitivity (worst used %.2g of it), %d boundary-ambiguous decisions",
			where, c.rows, c.conditioned, c.worstRatio, c.ambiguous)
	}
}

// valueOK reports whether a difference d from the reference is
// acceptable: d <= budget, or d <= budget + sensitivity() (counted).
func (c *condTally) valueOK(d, budget float64, sensitivity func() float64) bool {
	if d <= budget {
		return true
	}
	s := sensitivity()
	if d <= budget+s {
		c.conditioned++
		c.worstRatio = math.Max(c.worstRatio, (d-budget)/s)
		return true
	}
	return false
}

// TestLibmPerturbHook: the hook moves results by the requested ULP and
// leaves production untouched when unset; the sensitivity of a
// well-conditioned projection point is a few ULP, of an ill-conditioned
// one (STG next to its divergent anti-pole) far more.
func TestLibmPerturbHook(t *testing.T) {
	x := 0.7
	if lmSin(x) != math.Sin(x) {
		t.Fatal("lmSin differs from math.Sin with the hook unset")
	}
	withLibmPerturbed(3, false, 0, func() {
		if got, want := lmSin(x), nudgedF(math.Sin(x), 3); got != want {
			t.Errorf("perturbed lmSin = %v, want %v", got, want)
		}
		if got, want := lmCos(-2.5), nudgedF(math.Cos(-2.5), -3); got != want {
			t.Errorf("perturbed lmCos (negative) = %v, want %v (|v| up 3 ULP)", got, want)
		}
		if lmSin(0) != 0 {
			t.Error("exact zero perturbed")
		}
	})
	if lmSin(x) != math.Sin(x) || libmPerturbULP != 0 {
		t.Fatal("hook left set")
	}
	stg, err := Select("STG", nil, 2)
	if err != nil {
		t.Fatal(err)
	}
	ulpDist := func(a0, a1, b0, b1 float64) float64 {
		u := func(a, b float64) float64 {
			v := math.Abs(a)
			return math.Abs(a-b) / (math.Nextafter(v, math.Inf(1)) - v)
		}
		return math.Max(u(a0, b0), u(a1, b1))
	}
	well := libmSensitivity(func() (float64, float64, bool) { return wrefEvalDeg(stg, false, 30, 40) }, ulpDist)
	ill := libmSensitivity(func() (float64, float64, bool) { return wrefEvalDeg(stg, false, 30, -89.99) }, ulpDist)
	t.Logf("STG forward sensitivity to %d ULP of libm: %.3g ULP at theta = 40, %.3g ULP at theta = -89.99", libmSensULP, well, ill)
	if !(well <= 20) || !(ill > 1e6) {
		t.Errorf("STG sensitivity %.3g ULP (theta 40, want <= 20) and %.3g ULP (theta -89.99, want > 1e6)", well, ill)
	}
}
