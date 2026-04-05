package transform

import (
	"fmt"
	"math"
)

// solveNewton finds x such that f(x) = target to within tolerance, using
// Newton's method with bracketed bisection fallback. fprime is the
// derivative of f.
//
// Used by projections without closed-form inverses: MOL (Mollweide) needs
// to solve 2*psi + sin(2*psi) = pi*sin(theta) for psi, PCO (polyconic)
// needs a similar 1D root find, and distortion inversions (TPV/TNX) fall
// back on iterative solve when the inverse polynomial is not supplied.
//
// x0 is the initial guess; [lo, hi] brackets the root (both must evaluate
// to values bracketing target — i.e. (f(lo)-target) and (f(hi)-target)
// have opposite signs). If Newton produces a step outside [lo, hi] or
// fails to converge, the function falls back to bisection on the bracket
// to guarantee convergence.
//
// Returns an error if convergence fails within maxIter iterations (which
// should never happen for well-behaved problems but is detected rather
// than silently returning a bad value).
func solveNewton(f, fprime func(float64) float64, target, x0, lo, hi, tol float64, maxIter int) (float64, error) {
	x := x0
	// Ensure the bracket is valid.
	fLo := f(lo) - target
	fHi := f(hi) - target
	if fLo*fHi > 0 {
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
		if fx*fLo < 0 {
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
// Always converges if the bracket is valid; slower than solveNewton but
// does not need a derivative.
func solveBisection(f func(float64) float64, target, lo, hi, tol float64, maxIter int) (float64, error) {
	fLo := f(lo) - target
	fHi := f(hi) - target
	if fLo*fHi > 0 {
		return 0, fmt.Errorf("wcs/transform: bisection bracket [%g, %g] does not straddle target %g", lo, hi, target)
	}
	for i := 0; i < maxIter; i++ {
		mid := (lo + hi) / 2
		fMid := f(mid) - target
		if math.Abs(fMid) < tol || (hi-lo) < tol {
			return mid, nil
		}
		if fMid*fLo < 0 {
			hi = mid
			fHi = fMid
		} else {
			lo = mid
			fLo = fMid
		}
	}
	return (lo + hi) / 2, nil
}
