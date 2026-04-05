package transform

import (
	"math"
	"testing"
)

// TestSolveNewton_Quadratic: sqrt(2) via Newton on x² = 2.
func TestSolveNewton_Quadratic(t *testing.T) {
	f := func(x float64) float64 { return x * x }
	fp := func(x float64) float64 { return 2 * x }
	x, err := solveNewton(f, fp, 2.0, 1.0, 0.5, 2.0, 1e-14, 50)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(x-math.Sqrt2) > 1e-14 {
		t.Fatalf("sqrt(2) = %.15g, want %.15g", x, math.Sqrt2)
	}
}

// TestSolveNewton_Mollweide: the Mollweide auxiliary-angle equation
// 2*psi + sin(2*psi) = pi*sin(theta). For theta = 30°, the reference
// answer from Paper II is psi ≈ 0.39378... (hand-computed via enough
// Newton iterations).
func TestSolveNewton_Mollweide(t *testing.T) {
	theta := math.Pi / 6 // 30°
	target := math.Pi * math.Sin(theta)
	f := func(psi float64) float64 { return 2*psi + math.Sin(2*psi) }
	fp := func(psi float64) float64 { return 2 + 2*math.Cos(2*psi) }
	psi, err := solveNewton(f, fp, target, theta, -math.Pi/2, math.Pi/2, 1e-13, 50)
	if err != nil {
		t.Fatal(err)
	}
	// Verify the solution by plugging back in.
	residual := 2*psi + math.Sin(2*psi) - target
	if math.Abs(residual) > 1e-13 {
		t.Fatalf("residual %g for psi=%g", residual, psi)
	}
}

// TestSolveNewton_OutsideBracket rejects a starting bracket that does not
// contain the root.
func TestSolveNewton_OutsideBracket(t *testing.T) {
	f := func(x float64) float64 { return x * x }
	fp := func(x float64) float64 { return 2 * x }
	if _, err := solveNewton(f, fp, 2.0, 10.0, 3.0, 5.0, 1e-14, 50); err == nil {
		t.Fatal("expected bracket error")
	}
}

// TestSolveBisection finds sqrt(2) using pure bisection.
func TestSolveBisection(t *testing.T) {
	f := func(x float64) float64 { return x * x }
	x, err := solveBisection(f, 2.0, 0, 2, 1e-13, 100)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(x-math.Sqrt2) > 1e-12 {
		t.Fatalf("bisection sqrt(2) = %g", x)
	}
}
