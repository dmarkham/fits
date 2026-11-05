package transform

import (
	"fmt"
	"math"
)

// solveNewton finds x such that f(x) = target to within tolerance, using
// Newton's method with bracketed bisection fallback. fprime is the
// derivative of f.
//
// NOT used by any projection or distortion in this package: the
// projections without closed-form directions (MOL forward; PCO, AIR and
// ZPN degree >= 3 inverse) use wcslib's own bisection / weighted
// regula-falsi loops, ported verbatim so their stopping rules match
// wcslib's (see pseudo_cylindrical.go, polyconic.go, zenithal_extra.go),
// and TPV/TNX/SIP invert with their own iterations (tpv.go, tnx.go,
// sip.go). This general solver is kept, with its own tests
// (solver_test.go), for future use.
//
// x0 is the initial guess; [lo, hi] brackets the root (both must evaluate
// to values bracketing target — i.e. (f(lo)-target) and (f(hi)-target)
// have opposite signs). If Newton produces a step outside [lo, hi] or
// fails to converge, the function falls back to bisection on the bracket.
//
// Returns an error if convergence fails within maxIter iterations (which
// should never happen for well-behaved problems but is detected rather
// than silently returning a bad value).
func solveNewton(f, fprime func(float64) float64, target, x0, lo, hi, tol float64, maxIter int) (float64, error) {
	x := x0
	// Ensure the bracket is valid.
	fLo := f(lo) - target
	fHi := f(hi) - target
	if float64(fLo*fHi) > 0 {
		return 0, fmt.Errorf("wcs/transform: solver bracket [%g, %g] does not straddle target %g (f values %g, %g)", lo, hi, target, fLo+target, fHi+target)
	}
	for i := 0; i < maxIter; i++ {
		fx := f(x) - target
		if math.Abs(fx) < tol {
			return x, nil
		}
		dfx := fprime(x)
		// Newton step, clipped to stay inside the bracket.
		if dfx != 0 {
			next := x - fx/dfx
			if next > lo && next < hi {
				x = next
				continue
			}
		}
		// Bisection fallback: tighten the bracket using the current x.
		if float64(fx*fLo) < 0 {
			hi = x
			fHi = fx
		} else {
			lo = x
			fLo = fx
		}
		x = (lo + hi) / 2
	}
	return 0, fmt.Errorf("wcs/transform: solver failed to converge after %d iterations (last x=%g, residual=%g)", maxIter, x, f(x)-target)
}

// solveBisection finds x such that f(x) = target on the bracket [lo, hi].
// Slower than solveNewton but needs no derivative. Returns an error if the
// bracket does not straddle target, or if maxIter halvings leave neither
// |f(x) - target| nor the bracket width below tol (the caller asked for
// more precision than maxIter can deliver; the midpoint is NOT returned
// as if it were a root).
func solveBisection(f func(float64) float64, target, lo, hi, tol float64, maxIter int) (float64, error) {
	fLo := f(lo) - target
	fHi := f(hi) - target
	if float64(fLo*fHi) > 0 {
		return 0, fmt.Errorf("wcs/transform: bisection bracket [%g, %g] does not straddle target %g", lo, hi, target)
	}
	for i := 0; i < maxIter; i++ {
		mid := (lo + hi) / 2
		fMid := f(mid) - target
		if math.Abs(fMid) < tol || (hi-lo) < tol {
			return mid, nil
		}
		if float64(fMid*fLo) < 0 {
			hi = mid
			fHi = fMid
		} else {
			lo = mid
			fLo = fMid
		}
	}
	return 0, fmt.Errorf("wcs/transform: bisection did not reach tolerance %g in %d iterations (bracket [%g, %g])", tol, maxIter, lo, hi)
}
