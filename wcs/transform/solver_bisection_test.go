package transform

import (
	"math"
	"testing"
)

// TestSolveBisectionNonConvergenceIsAnError pins the solveBisection
// contract: when maxIter is exhausted before the bracket shrinks below tol
// it must return an error, not the bracket midpoint dressed up as a root.
// (Before the fix it returned (midpoint, nil) after maxIter.)
func TestSolveBisectionNonConvergenceIsAnError(t *testing.T) {
	f := func(x float64) float64 { return x * x }
	// Root of x^2 = 2 on [0, 2]; 5 halvings leave a bracket of 2/32.
	x, err := solveBisection(f, 2.0, 0, 2, 1e-13, 5)
	if err == nil {
		t.Fatalf("solveBisection with 5 iterations returned x=%.17g (|x-sqrt2|=%.3g), nil error; want a non-convergence error", x, math.Abs(x-math.Sqrt2))
	}
	// And with enough iterations it converges with no error.
	x, err = solveBisection(f, 2.0, 0, 2, 1e-13, 200)
	if err != nil || math.Abs(x-math.Sqrt2) > 1e-12 {
		t.Fatalf("solveBisection converged run: x=%.17g err=%v", x, err)
	}
}
