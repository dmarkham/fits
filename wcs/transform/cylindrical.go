package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file implements the cylindrical family of projections from Paper II
// §5.2. All have theta_0 = 0 (reference point on the native equator).

// CAR — plate carrée (§5.2.3).
//
//	x = phi
//	y = theta
//
// The simplest projection: longitude is x, latitude is y, both in radians.
// Valid for the entire sphere. Used for all-sky atlases.
type carProjection struct{}

func (carProjection) Code() string    { return "CAR" }
func (carProjection) Theta0() float64 { return 0 }

func (carProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	return phi, theta, true
}

func (carProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	return x, y, true
}

// CYP — cylindrical perspective (Paper II §5.2.1).
//
// Parameters (latitude axis):
//
//	PV2_1 = mu     (distance from axis to projection point, default 1)
//	PV2_2 = lambda (cylinder radius / projection scale, default 1)
//
// Forward: x = lambda * phi, y = (mu + 1) * sin(theta) / (mu + cos(theta)).
// Special cases include Gall's stereographic (mu=1, lambda=sqrt(2)/2),
// central cylindrical (mu=0, lambda=1), and Miller-style variants.
type cypProjection struct {
	mu     float64
	lambda float64
}

func newCYP(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return cypProjection{
		mu:     pvFloat(pv, latAxis, 1, 1),
		lambda: pvFloat(pv, latAxis, 2, 1),
	}
}

func (cypProjection) Code() string    { return "CYP" }
func (cypProjection) Theta0() float64 { return 0 }

func (p cypProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	sinT, cosT := math.Sincos(theta)
	denom := p.mu + cosT
	if denom == 0 {
		return 0, 0, false
	}
	x = p.lambda * phi
	y = (p.mu + 1) * sinT / denom
	return x, y, true
}

func (p cypProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	if p.lambda == 0 {
		return 0, 0, false
	}
	phi = x / p.lambda
	// Solve y = (mu+1)*sin(theta) / (mu + cos(theta)) for theta.
	// Let eta = y / (mu+1). Then sin(theta) = eta*(mu + cos(theta)).
	// Square: sin²(theta) = eta²*(mu + cos(theta))². With sin²=1-cos².
	// This is a quadratic in cos(theta). The stable form uses atan2:
	// tan(theta/2) is more direct via half-angle identity.
	// Paper II gives the closed form:
	//   theta = atan2(eta, 1) + asin(eta*mu / sqrt(eta² + 1))
	eta := y / (p.mu + 1)
	denom := math.Sqrt(eta*eta + 1)
	arg := eta * p.mu / denom
	if arg < -1 || arg > 1 {
		return 0, 0, false
	}
	theta = math.Atan2(eta, 1) + math.Asin(arg)
	return phi, theta, true
}

// CEA — cylindrical equal-area (Paper II §5.2.2).
//
// Parameter (latitude axis):
//
//	PV2_1 = lambda (squash factor, default 1 for Lambert's projection)
//
// Forward: x = phi, y = sin(theta) / lambda.
// Inverse: theta = asin(y * lambda).
type ceaProjection struct {
	lambda float64
}

func newCEA(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return ceaProjection{
		lambda: pvFloat(pv, latAxis, 1, 1),
	}
}

func (ceaProjection) Code() string    { return "CEA" }
func (ceaProjection) Theta0() float64 { return 0 }

func (p ceaProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if p.lambda == 0 {
		return 0, 0, false
	}
	return phi, math.Sin(theta) / p.lambda, true
}

func (p ceaProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	arg := y * p.lambda
	if arg < -1 || arg > 1 {
		return 0, 0, false
	}
	return x, math.Asin(arg), true
}

// MER — Mercator (Paper II §5.2.4).
//
//	x = phi
//	y = ln(tan(pi/4 + theta/2))
//
// Diverges at theta = ±pi/2. Used for navigation charts and historical
// sky atlases.
type merProjection struct{}

func (merProjection) Code() string    { return "MER" }
func (merProjection) Theta0() float64 { return 0 }

func (merProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if theta >= math.Pi/2 || theta <= -math.Pi/2 {
		return 0, 0, false
	}
	return phi, math.Log(math.Tan(math.Pi/4 + theta/2)), true
}

func (merProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	return x, 2*math.Atan(math.Exp(y)) - math.Pi/2, true
}
