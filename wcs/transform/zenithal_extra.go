package transform

import (
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file contains the zenithal perspective family (AZP, SZP) and the
// zenithal polynomial (ZPN), separated from the simpler closed-form
// zenithal projections in zenithal.go.

// AZP — zenithal perspective projection (Paper II §5.1.1).
//
// Parameters (on the latitude axis, conventionally axis 2):
//
//	PV2_1 = mu    (dimensionless distance from sphere center to
//	              projection source, in sphere radii; default 0)
//	PV2_2 = gamma (tilt angle in degrees, default 0)
//
// Forward: Paper II eq. 20.
//
//	R_theta = (mu + 1) * cos(theta) / (mu + sin(theta))
//
// (The gamma tilt adds a cos(gamma) denominator term when rotating the
// projection plane; see eq. 19.) For gamma = 0 this reduces to the
// un-tilted zenithal perspective.
//
// Inverse: quadratic in sin(theta). We derive it below from the forward.
type azpProjection struct {
	mu    float64
	gamma float64 // tilt, radians
}

func newAZP(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return azpProjection{
		mu:    pvFloat(pv, latAxis, 1, 0),
		gamma: degToRad(pvFloat(pv, latAxis, 2, 0)),
	}
}

func (azpProjection) Code() string    { return "AZP" }
func (azpProjection) Theta0() float64 { return math.Pi / 2 }

func (p azpProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	sinPhi, cosPhi := math.Sincos(phi)
	sinTheta, cosTheta := math.Sincos(theta)
	sinGamma, cosGamma := math.Sincos(p.gamma)
	denom := p.mu + sinTheta + cosTheta*cosPhi*sinGamma
	if denom == 0 {
		return 0, 0, false
	}
	r := (p.mu + 1) * cosTheta / denom
	x = r * sinPhi
	y = -r * cosPhi / cosGamma
	return x, y, true
}

func (p azpProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// For gamma = 0 (the common case with no tilt) we have a closed-form
	// quadratic inverse. Setting u = sin(theta) and substituting into
	//
	//	R = (mu+1) * cos(theta) / (mu + sin(theta))
	//
	// gives the quadratic
	//
	//	u² * (R² + (mu+1)²) + 2*u*R²*mu + (R²*mu² - (mu+1)²) = 0
	//
	// whose larger root is the near-hemisphere solution. For gamma != 0
	// the formula still has a closed form but involves cos(phi) coupling
	// into the denominator; we handle gamma != 0 by scaling y first then
	// applying the gamma=0 derivation (Paper II eq. 23 with the tilt
	// correction absorbed into y).
	cosGamma := math.Cos(p.gamma)
	if cosGamma == 0 {
		return 0, 0, false
	}
	yTilt := y * cosGamma
	r := math.Hypot(x, yTilt)
	if r == 0 {
		return 0, math.Pi / 2, true
	}
	phi = math.Atan2(x, -yTilt)

	// Quadratic coefficients. Without tilt, this is exact; with gamma != 0
	// it is an approximation that is exact along the principal meridians.
	muPlus1 := p.mu + 1
	a := r*r + muPlus1*muPlus1
	b := 2 * r * r * p.mu
	c := r*r*p.mu*p.mu - muPlus1*muPlus1
	disc := b*b - 4*a*c
	if disc < 0 {
		return 0, 0, false
	}
	sqrtDisc := math.Sqrt(disc)
	// Larger root → near hemisphere.
	u := (-b + sqrtDisc) / (2 * a)
	if u < -1 || u > 1 {
		u = (-b - sqrtDisc) / (2 * a)
		if u < -1 || u > 1 {
			return 0, 0, false
		}
	}
	sinTheta := u
	cosTheta := math.Sqrt((1 - u) * (1 + u))
	theta = math.Atan2(sinTheta, cosTheta)
	return phi, theta, true
}

// SZP — slant zenithal perspective (Paper II §5.1.2).
//
// Parameters (latitude axis):
//
//	PV2_1 = mu      (source distance in sphere radii; default 0)
//	PV2_2 = phi_c   (native longitude of the projection source, degrees)
//	PV2_3 = theta_c (native latitude of the projection source, degrees,
//	                default 90)
//
// This is the most general zenithal perspective: the projection source
// can be anywhere along the axis of a cone whose apex is at theta_c.
// Used by some planetary stereo-pair instruments. When mu = 0 and
// theta_c = 90 it reduces to TAN; when theta_c = 90 and mu != 0 it
// reduces to AZP.
type szpProjection struct {
	mu     float64
	phiC   float64 // radians
	thetaC float64 // radians
}

func newSZP(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return szpProjection{
		mu:     pvFloat(pv, latAxis, 1, 0),
		phiC:   degToRad(pvFloat(pv, latAxis, 2, 0)),
		thetaC: degToRad(pvFloat(pv, latAxis, 3, 90)),
	}
}

func (szpProjection) Code() string    { return "SZP" }
func (szpProjection) Theta0() float64 { return math.Pi / 2 }

func (p szpProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	// Paper II eq. 24:
	//   X_p = mu * cos(theta_c) * sin(phi_c)
	//   Y_p = -mu * cos(theta_c) * cos(phi_c)
	//   Z_p = mu * sin(theta_c) + 1
	// And the projection (Paper II eq. 25):
	//   X = (Z_p * cos(theta)*sin(phi) - X_p*(1 - sin(theta))) / (Z_p - (1 - sin(theta)))
	//   Y = -(Z_p * cos(theta)*cos(phi) + Y_p*(1 - sin(theta))) / (Z_p - (1 - sin(theta)))
	sinPhiC, cosPhiC := math.Sincos(p.phiC)
	sinThC, cosThC := math.Sincos(p.thetaC)
	Xp := p.mu * cosThC * sinPhiC
	Yp := -p.mu * cosThC * cosPhiC
	Zp := p.mu*sinThC + 1

	sinTheta, cosTheta := math.Sincos(theta)
	sinPhi, cosPhi := math.Sincos(phi)
	s := 1 - sinTheta
	denom := Zp - s
	if denom == 0 {
		return 0, 0, false
	}
	x = (Zp*cosTheta*sinPhi - Xp*s) / denom
	y = -(Zp*cosTheta*cosPhi + Yp*s) / denom
	return x, y, true
}

func (p szpProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// The closed-form inverse of SZP requires solving a quadratic in
	// sin(theta). We implement it via Newton iteration against Forward,
	// which is robust and reuses the forward math verbatim. Initial
	// guess: assume the source is at infinity (like TAN) for a first
	// approximation of theta and phi.
	r := math.Hypot(x, y)
	phi0 := math.Atan2(x, -y)
	theta0 := math.Atan2(1, r)

	// Use two-variable Newton by reducing to 1D via fixed-point: given
	// the current (phi, theta), project via Forward to get (x', y'), then
	// adjust theta so that |(x, y) - (x', y')| shrinks. A simpler robust
	// strategy: solve a 2x2 Newton system. For SZP the system is well-
	// conditioned away from the source point.
	phi, theta = phi0, theta0
	for i := 0; i < 50; i++ {
		xc, yc, fok := p.Forward(phi, theta)
		if !fok {
			return 0, 0, false
		}
		dx := x - xc
		dy := y - yc
		if math.Hypot(dx, dy) < 1e-13 {
			return phi, theta, true
		}
		// Numerical Jacobian.
		const h = 1e-8
		x1, y1, _ := p.Forward(phi+h, theta)
		x2, y2, _ := p.Forward(phi, theta+h)
		j11 := (x1 - xc) / h
		j12 := (x2 - xc) / h
		j21 := (y1 - yc) / h
		j22 := (y2 - yc) / h
		det := j11*j22 - j12*j21
		if det == 0 {
			return 0, 0, false
		}
		dphi := (j22*dx - j12*dy) / det
		dtheta := (-j21*dx + j11*dy) / det
		phi += dphi
		theta += dtheta
		// Clamp theta to valid range.
		if theta > math.Pi/2 {
			theta = math.Pi / 2
		}
		if theta < -math.Pi/2 {
			theta = -math.Pi / 2
		}
	}
	return 0, 0, false
}

// ZPN — zenithal polynomial projection (Paper II §5.1.8).
//
// R_theta is a user-defined polynomial in (pi/2 - theta):
//
//	R_theta = PV2_0 + PV2_1 * zeta + PV2_2 * zeta² + ... + PV2_N * zeta^N
//
// where zeta = pi/2 - theta (in radians; Paper II uses degrees but we
// convert). Up to N = 20 coefficients.
//
// Forward is trivial polynomial evaluation. Inverse requires finding the
// real root of the polynomial in [0, pi] that matches R. We use
// solveNewton with the polynomial derivative.
type zpnProjection struct {
	coeffs [21]float64 // PV2_0 .. PV2_20
	degree int         // highest non-zero index
}

func newZPN(pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	p := zpnProjection{}
	for i := 0; i <= 20; i++ {
		p.coeffs[i] = pvFloat(pv, latAxis, i, 0)
		if p.coeffs[i] != 0 {
			p.degree = i
		}
	}
	if p.degree == 0 && p.coeffs[0] == 0 {
		return nil, fmt.Errorf("wcs/transform: ZPN requires at least one non-zero PV%d_i coefficient", latAxis)
	}
	return p, nil
}

func (zpnProjection) Code() string    { return "ZPN" }
func (zpnProjection) Theta0() float64 { return math.Pi / 2 }

// polyEval evaluates the ZPN polynomial at zeta.
func (p zpnProjection) polyEval(zeta float64) float64 {
	// Horner form.
	v := p.coeffs[p.degree]
	for i := p.degree - 1; i >= 0; i-- {
		v = v*zeta + p.coeffs[i]
	}
	return v
}

func (p zpnProjection) polyDeriv(zeta float64) float64 {
	if p.degree == 0 {
		return 0
	}
	v := float64(p.degree) * p.coeffs[p.degree]
	for i := p.degree - 1; i >= 1; i-- {
		v = v*zeta + float64(i)*p.coeffs[i]
	}
	return v
}

func (p zpnProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	zeta := math.Pi/2 - theta
	r := p.polyEval(zeta)
	return r * math.Sin(phi), -r * math.Cos(phi), true
}

func (p zpnProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := math.Hypot(x, y)
	if r == 0 {
		// At the origin, find zeta such that polyEval(zeta) = 0.
		// Usually zeta = 0 (coeffs[0] = 0) but not guaranteed.
		if p.coeffs[0] == 0 {
			return 0, math.Pi / 2, true
		}
	}
	phi = math.Atan2(x, -y)
	// Initial guess: assume linear leading term, i.e. zeta ≈ r / coeffs[1]
	// if coeffs[1] != 0, else r / coeffs[p.degree].
	var guess float64
	if p.coeffs[1] != 0 {
		guess = r / p.coeffs[1]
	} else {
		guess = math.Pow(r/p.coeffs[p.degree], 1/float64(p.degree))
	}
	if guess < 0 {
		guess = 0
	}
	if guess > math.Pi {
		guess = math.Pi
	}
	zeta, err := solveNewton(p.polyEval, p.polyDeriv, r, guess, 0, math.Pi, 1e-13, 50)
	if err != nil {
		return 0, 0, false
	}
	theta = math.Pi/2 - zeta
	return phi, theta, true
}
