package transform

import "math"

// This file implements the pseudo-cylindrical and conventional projections
// from Paper II §5.3.
//
// All have theta_0 = 0 (reference point on the equator). The characteristic
// of pseudo-cylindrical is that parallels are horizontal lines (as in
// cylindrical) but meridians are curves.

// SFL — Sanson-Flamsteed sinusoidal (§5.3.1).
//
//	x = phi * cos(theta)
//	y = theta
//
// Equal-area. No projection parameters. Poles map to single points.
type sflProjection struct{}

func (sflProjection) Code() string    { return "SFL" }
func (sflProjection) Theta0() float64 { return 0 }

func (sflProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	return phi * math.Cos(theta), theta, true
}

func (sflProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	theta = y
	c := math.Cos(theta)
	if c == 0 {
		// At the pole, phi is indeterminate; pick 0.
		return 0, theta, true
	}
	phi = x / c
	return phi, theta, true
}

// PAR — parabolic (§5.3.2).
//
//	x = phi * (2*cos(2*theta/3) - 1)
//	y = pi * sin(theta/3)
//
// Equal-area. No parameters.
type parProjection struct{}

func (parProjection) Code() string    { return "PAR" }
func (parProjection) Theta0() float64 { return 0 }

func (parProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	x = phi * (2*math.Cos(2*theta/3) - 1)
	y = math.Pi * math.Sin(theta/3)
	return x, y, true
}

func (parProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	arg := y / math.Pi
	if arg < -1 || arg > 1 {
		return 0, 0, false
	}
	theta = 3 * math.Asin(arg)
	d := 2*math.Cos(2*theta/3) - 1
	if d == 0 {
		return 0, theta, true
	}
	phi = x / d
	return phi, theta, true
}

// MOL — Mollweide (§5.3.3).
//
// Auxiliary angle psi satisfies 2*psi + sin(2*psi) = pi*sin(theta).
// Then:
//
//	x = 2*sqrt(2)/pi * phi * cos(psi)
//	y = sqrt(2) * sin(psi)
//
// Equal-area. Uses iterative solver for psi in the forward direction;
// inverse is closed-form.
type molProjection struct{}

func (molProjection) Code() string    { return "MOL" }
func (molProjection) Theta0() float64 { return 0 }

func (molProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	// Solve 2*psi + sin(2*psi) = pi*sin(theta).
	target := math.Pi * math.Sin(theta)
	f := func(psi float64) float64 { return 2*psi + math.Sin(2*psi) }
	fp := func(psi float64) float64 { return 2 + 2*math.Cos(2*psi) }
	psi, err := solveNewton(f, fp, target, theta, -math.Pi/2, math.Pi/2, 1e-13, 50)
	if err != nil {
		return 0, 0, false
	}
	x = 2 * math.Sqrt2 / math.Pi * phi * math.Cos(psi)
	y = math.Sqrt2 * math.Sin(psi)
	return x, y, true
}

func (molProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// y = sqrt(2)*sin(psi) → sin(psi) = y/sqrt(2).
	sinPsi := y / math.Sqrt2
	if sinPsi < -1 || sinPsi > 1 {
		return 0, 0, false
	}
	psi := math.Asin(sinPsi)
	cosPsi := math.Cos(psi)
	// sin(theta) = (2*psi + sin(2*psi)) / pi
	sinTheta := (2*psi + math.Sin(2*psi)) / math.Pi
	if sinTheta < -1 || sinTheta > 1 {
		if sinTheta > 1 {
			sinTheta = 1
		} else {
			sinTheta = -1
		}
	}
	theta = math.Asin(sinTheta)
	if cosPsi == 0 {
		return 0, theta, true
	}
	phi = math.Pi * x / (2 * math.Sqrt2 * cosPsi)
	return phi, theta, true
}

// AIT — Hammer-Aitoff (§5.3.4).
//
//	gamma = sqrt(2 / (1 + cos(theta)*cos(phi/2)))
//	x = 2 * gamma * cos(theta) * sin(phi/2)
//	y = gamma * sin(theta)
//
// Equal-area, 2:1 aspect. The canonical all-sky astronomy projection.
type aitProjection struct{}

func (aitProjection) Code() string    { return "AIT" }
func (aitProjection) Theta0() float64 { return 0 }

func (aitProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	cosT := math.Cos(theta)
	denom := 1 + cosT*math.Cos(phi/2)
	if denom <= 0 {
		return 0, 0, false
	}
	gamma := math.Sqrt(2 / denom)
	x = 2 * gamma * cosT * math.Sin(phi/2)
	y = gamma * math.Sin(theta)
	return x, y, true
}

func (aitProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// Closed-form inverse: let z = sqrt(1 - (x/4)² - (y/2)²).
	// Then theta = asin(y*z), phi = 2*atan2(z*x/2, 2*z² - 1).
	a := x / 4
	b := y / 2
	u := 1 - a*a - b*b
	if u < 0 {
		return 0, 0, false
	}
	z := math.Sqrt(u)
	theta = math.Asin(y * z)
	phi = 2 * math.Atan2(z*x/2, 2*z*z-1)
	return phi, theta, true
}
