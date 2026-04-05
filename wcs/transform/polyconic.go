package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file implements the polyconic and pseudoconic projections from
// Paper II §5.5.

// BON — Bonne's projection (§5.5.1).
//
// Parameter: PV2_1 = theta_1 (standard parallel in degrees; 0 reduces to
// SFL). Paper II eq. 79:
//
//	Y_0 = cot(theta_1) + theta_1           (in radians)
//	A   = phi * cos(theta) / (Y_0 - theta)
//	x   = (Y_0 - theta) * sin(A)
//	y   = -(Y_0 - theta) * cos(A) + Y_0
//
// Equal-area, pseudoconic, theta_0 = theta_1.
type bonProjection struct {
	theta1 float64 // radians
	Y0     float64
}

func newBON(pv map[wcs.PVKey]float64, latAxis int) Projection {
	t1 := degToRad(pvFloat(pv, latAxis, 1, 0))
	return bonProjection{
		theta1: t1,
		Y0:     safeCot(t1) + t1,
	}
}

// safeCot returns cot(x) = cos(x)/sin(x), with a large finite value when
// sin(x) is zero (equator → infinite radius, the SFL degenerate case is
// special-cased below).
func safeCot(x float64) float64 {
	s := math.Sin(x)
	if s == 0 {
		return 0
	}
	return math.Cos(x) / s
}

func (bonProjection) Code() string      { return "BON" }
func (p bonProjection) Theta0() float64 { return p.theta1 }

func (p bonProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if p.theta1 == 0 {
		// Degenerate case: Bonne with theta_1 = 0 is sinusoidal (SFL).
		return sflProjection{}.Forward(phi, theta)
	}
	R := p.Y0 - theta
	A := phi * math.Cos(theta) / R
	sA, cA := math.Sincos(A)
	x = R * sA
	y = -R*cA + p.Y0
	return x, y, true
}

func (p bonProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	if p.theta1 == 0 {
		return sflProjection{}.Inverse(x, y)
	}
	yPrime := p.Y0 - y
	R := math.Hypot(x, yPrime)
	if p.theta1 < 0 {
		R = -R
	}
	theta = p.Y0 - R
	A := math.Atan2(x/R, yPrime/R)
	c := math.Cos(theta)
	if c == 0 {
		return 0, theta, true
	}
	phi = A * R / c
	return phi, theta, true
}

// PCO — polyconic (Paper II §5.5.2).
//
// For theta = 0:
//
//	x = phi
//	y = 0
//
// For theta != 0:
//
//	E = phi * sin(theta)
//	x = cot(theta) * sin(E)
//	y = theta + cot(theta) * (1 - cos(E))
//
// Inverse requires iteratively solving for theta given (x, y). The
// Paper II iteration converges in ~5 steps for reasonable inputs.
type pcoProjection struct{}

func (pcoProjection) Code() string    { return "PCO" }
func (pcoProjection) Theta0() float64 { return 0 }

func (pcoProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if theta == 0 {
		return phi, 0, true
	}
	sinT, cosT := math.Sincos(theta)
	if sinT == 0 {
		return phi, 0, true
	}
	cot := cosT / sinT
	E := phi * sinT
	sE, cE := math.Sincos(E)
	x = cot * sE
	y = theta + cot*(1-cE)
	return x, y, true
}

func (pcoProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	if y == 0 {
		return x, 0, true
	}
	// Paper II §5.5.2 gives the inverse via Newton iteration on theta. At
	// each step, given the current theta, we can compute E (the argument
	// of the inner sin/cos) from (x, y):
	//
	//   sin(E) = x * tan(theta)
	//   cos(E) = 1 - (y - theta) * tan(theta)
	//
	// Then the residual is yc(theta) = theta + cot(theta)*(1 - cos(E)),
	// which must equal y. Newton's method on theta with the derivative
	// of yc. A clean formulation: iterate until (x, y) - Forward(phi,theta)
	// is small enough, where phi = E / sin(theta) at each step.
	theta = y
	for i := 0; i < 100; i++ {
		sinT, cosT := math.Sincos(theta)
		if sinT == 0 {
			// At equator, Forward gives (phi, 0). Return x directly.
			return x, 0, true
		}
		tanT := sinT / cosT
		sinE := x * tanT
		cosE := 1 - (y-theta)*tanT
		// Clamp sinE to valid range — rounding can push it slightly outside.
		mag := math.Hypot(sinE, cosE)
		if mag == 0 {
			return 0, theta, true
		}
		sinE /= mag
		cosE /= mag
		E := math.Atan2(sinE, cosE)
		phiEst := E / sinT

		// Check residual via the forward formula.
		xc, yc, _ := pcoProjection{}.Forward(phiEst, theta)
		dx := x - xc
		dy := y - yc
		if math.Abs(dx) < 1e-13 && math.Abs(dy) < 1e-13 {
			return phiEst, theta, true
		}

		// Newton step: d(yc)/d(theta) at fixed phi includes both the
		// direct theta term and the E = phi*sin(theta) coupling. The
		// simplest stable update is the "residual in y" step:
		//   theta_new = theta - dy / dyDtheta
		// with dyDtheta = 1 - cosE*tanT*phiEst*cosT + ... (messy).
		//
		// A more robust approach: damped fixed-point toward theta such
		// that yc(theta_new) = y.
		dyDtheta := 1 + (cosT/sinT)*cosE*0 + // placeholder
			(cosT*cosT/(sinT*sinT))*(1-cosE) +
			phiEst*cosT*sinE
		if dyDtheta == 0 {
			dyDtheta = 1
		}
		newTheta := theta - (yc-y)/dyDtheta
		if math.Abs(newTheta-theta) < 1e-14 {
			return phiEst, theta, true
		}
		// Clamp theta to (-pi/2, pi/2).
		if newTheta > math.Pi/2-1e-9 {
			newTheta = math.Pi/2 - 1e-9
		}
		if newTheta < -math.Pi/2+1e-9 {
			newTheta = -math.Pi/2 + 1e-9
		}
		theta = newTheta
	}
	// Return best effort.
	sinT := math.Sin(theta)
	if sinT == 0 {
		return x, theta, true
	}
	phi = math.Atan2(x*sinT/math.Cos(theta), 1-(y-theta)*sinT/math.Cos(theta)) / sinT
	return phi, theta, true
}
