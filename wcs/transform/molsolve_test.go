package transform

import (
	"math"
	"math/rand"
	"testing"
)

// molBisectPlain is wcslib 8.10 mols2x's bisection exactly as
// prj.c writes it, evaluating the residual at every step.
func molBisectPlain(u, tol float64) float64 {
	v0, v1, v := -wpi, wpi, u
	for k := 0; k < 100; k++ {
		resid := (v - u) + lmSin(v)
		if resid < 0.0 {
			if resid > -tol {
				break
			}
			v0 = v
		} else {
			if resid < tol {
				break
			}
			v1 = v
		}
		v = (v0 + v1) / 2.0
	}
	return v
}

// TestMolSolveMatchesBisection: molSolve, which skips the sin evaluations
// far from the root, stops at exactly the iterate of the plain bisection
// (bit for bit) for latitudes over the whole sphere, dense next to the
// poles (90 - |theta| from 1e-12 to 1 deg) and the equator, and under the
// libm perturbation the tests use (so the skip bounds stay valid when
// sin moves by a few ULP).
func TestMolSolveMatchesBisection(t *testing.T) {
	const tol = 1.0e-13
	rng := rand.New(rand.NewSource(20261101))
	var thetas []float64
	for i := 0; i < 20000; i++ {
		thetas = append(thetas, rng.Float64()*180-90)
	}
	for e := -12.0; e <= 0; e += 0.01 {
		d := math.Pow(10, e)
		thetas = append(thetas, 90-d, d-90, d, -d)
	}
	n, bad := 0, 0
	check := func(label string) {
		for _, th := range thetas {
			u := float64(wpi * sind(th))
			got, want := molSolve(u, tol), molBisectPlain(u, tol)
			n++
			if math.Float64bits(got) != math.Float64bits(want) {
				if bad++; bad <= 10 {
					t.Errorf("%s theta %.17g (u %.17g): molSolve %.17g, plain bisection %.17g", label, th, u, got, want)
				}
			}
		}
	}
	check("libm")
	for _, k := range []int64{1, -1, 4, -4} {
		withLibmPerturbed(k, false, 0, func() { check("perturbed") })
	}
	if n < 100000 {
		t.Fatalf("only %d solves compared", n)
	}
}
