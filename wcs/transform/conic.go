package transform

import (
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file implements the conic projection family from Paper II §5.4.
// All four conic projections share the same parameter set:
//
//	PV2_1 = theta_a (midpoint of the two standard parallels, degrees)
//	PV2_2 = eta    (half the spread between the standard parallels,
//	                 degrees; default 0 for a single-parallel tangent case)
//
// The standard parallels are theta_1 = theta_a - eta and
// theta_2 = theta_a + eta.
//
// All conic projections have theta_0 = theta_a (not 0). The formulas
// share the following constants:
//
//	C = sin(theta_a)             (cone constant)
//	Y_0 = (180/pi) * cot(theta_a)  — in radians: cot(theta_a)
//
// Then:
//
//	x = R_theta * sin(C * phi)
//	y = -R_theta * cos(C * phi) + Y_0
//
// with R_theta differing per projection.

// conicBase holds the common parameters and provides helpers for the four
// concrete conic projections.
type conicBase struct {
	thetaA float64 // reference parallel, radians
	eta    float64 // spread, radians
	C      float64 // cone constant = sin(thetaA)
	Y0     float64 // cot(thetaA), radians
}

func newConicBase(pv map[wcs.PVKey]float64, latAxis int) (conicBase, error) {
	ta := degToRad(pvFloat(pv, latAxis, 1, 0))
	if ta == 0 {
		return conicBase{}, fmt.Errorf("wcs/transform: conic projection requires PV%d_1 (theta_a) != 0", latAxis)
	}
	eta := degToRad(pvFloat(pv, latAxis, 2, 0))
	C := math.Sin(ta)
	return conicBase{
		thetaA: ta,
		eta:    eta,
		C:      C,
		Y0:     math.Cos(ta) / math.Sin(ta),
	}, nil
}

// projectFromR takes an R_theta value and returns (x, y) using the shared
// conic (x, y) formula.
func (c conicBase) projectFromR(phi, rTheta float64) (x, y float64) {
	cPhi := c.C * phi
	s, co := math.Sincos(cPhi)
	return rTheta * s, -rTheta*co + c.Y0
}

// inversePhi extracts phi from (x, y) and an already-known R_theta.
func (c conicBase) inversePhi(x, y, rTheta float64) float64 {
	if rTheta == 0 || c.C == 0 {
		return 0
	}
	// x = R*sin(C*phi), y' = -R*cos(C*phi)  where y' = y - Y0.
	yPrime := y - c.Y0
	cPhi := math.Atan2(x/rTheta, -yPrime/rTheta)
	return cPhi / c.C
}

// inverseR returns the R value given (x, y): the radial distance from the
// cone's apex, preserving sign convention (R always positive for northern
// conics, always negative for southern).
func (c conicBase) inverseR(x, y float64) float64 {
	yPrime := y - c.Y0
	r := math.Hypot(x, yPrime)
	// Sign: matches sign of -yPrime because the base formula y = Y0 - R*cos(C*phi).
	// For small phi, cos(C*phi) ≈ 1, so y - Y0 ≈ -R, hence R ≈ -(y-Y0) = -yPrime.
	// But we also need negative R for southern conics. cfitsio uses the
	// convention R > 0 when C > 0 and handles southern separately. We mirror
	// that: use hypot for magnitude, and treat the sign via C.
	if c.C < 0 {
		return -r
	}
	return r
}

// COP — conic perspective (Paper II §5.4.1).
//
//	R_theta = cos(eta) * [cot(theta_a) - tan(theta - theta_a)]
type copProjection struct{ conicBase }

func newCOP(pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	b, err := newConicBase(pv, latAxis)
	if err != nil {
		return nil, err
	}
	return copProjection{b}, nil
}

func (copProjection) Code() string      { return "COP" }
func (p copProjection) Theta0() float64 { return p.thetaA }

func (p copProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	diff := theta - p.thetaA
	if math.Abs(diff) >= math.Pi/2 {
		return 0, 0, false
	}
	r := math.Cos(p.eta) * (p.Y0 - math.Tan(diff))
	x, y = p.projectFromR(phi, r)
	return x, y, true
}

func (p copProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := p.inverseR(x, y)
	phi = p.inversePhi(x, y, r)
	diff := math.Atan(p.Y0 - r/math.Cos(p.eta))
	theta = diff + p.thetaA
	return phi, theta, true
}

// COE — conic equal-area / Albers (Paper II §5.4.2).
//
//	R_theta = 2/(C_coe) * sqrt(1 + sin(theta_1)*sin(theta_2) - (sin(theta_1)+sin(theta_2))*sin(theta))
//
// where C_coe = (sin(theta_1) + sin(theta_2)) / 2. The cone constant
// differs from sin(theta_a) for COE; Paper II eq. 69.
type coeProjection struct {
	conicBase
	c2    float64 // 2*C_coe = sin(theta_1) + sin(theta_2)
	gamma float64 // sin(theta_1) + sin(theta_2) factor, cached
}

func newCOE(pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	b, err := newConicBase(pv, latAxis)
	if err != nil {
		return nil, err
	}
	s1 := math.Sin(b.thetaA - b.eta)
	s2 := math.Sin(b.thetaA + b.eta)
	// COE uses C = (sin(theta_1) + sin(theta_2)) / 2 as the cone constant.
	c := (s1 + s2) / 2
	b.C = c
	p := coeProjection{conicBase: b, c2: s1 + s2, gamma: 1 + s1*s2}
	// Y0 for COE is different: Paper II gives the radius at theta_a.
	p.Y0 = math.Sqrt(p.gamma-p.c2*math.Sin(b.thetaA)) / c
	return p, nil
}

func (coeProjection) Code() string      { return "COE" }
func (p coeProjection) Theta0() float64 { return p.thetaA }

func (p coeProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	u := p.gamma - p.c2*math.Sin(theta)
	if u < 0 {
		return 0, 0, false
	}
	r := math.Sqrt(u) / p.C
	x, y = p.projectFromR(phi, r)
	return x, y, true
}

func (p coeProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := p.inverseR(x, y)
	phi = p.inversePhi(x, y, r)
	arg := (p.gamma - (r*p.C)*(r*p.C)) / p.c2
	if arg < -1 || arg > 1 {
		return 0, 0, false
	}
	theta = math.Asin(arg)
	return phi, theta, true
}

// COD — conic equidistant (Paper II §5.4.3).
//
//	R_theta = theta_a - theta + cot(theta_a)
//
// (The Y_0 offset in projectFromR is cot(theta_a); the extra theta_a
// shift here combines into the full formula from Paper II eq. 73.)
type codProjection struct{ conicBase }

func newCOD(pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	b, err := newConicBase(pv, latAxis)
	if err != nil {
		return nil, err
	}
	// COD's cone constant is C = sin(theta_a) * sin(eta)/eta (from Paper II
	// eq. 74), which simplifies to sin(theta_a) when eta = 0.
	if b.eta != 0 {
		b.C = math.Sin(b.thetaA) * math.Sin(b.eta) / b.eta
	}
	return codProjection{b}, nil
}

func (codProjection) Code() string      { return "COD" }
func (p codProjection) Theta0() float64 { return p.thetaA }

func (p codProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	r := p.thetaA - theta + p.Y0
	x, y = p.projectFromR(phi, r)
	return x, y, true
}

func (p codProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := p.inverseR(x, y)
	phi = p.inversePhi(x, y, r)
	theta = p.thetaA + p.Y0 - r
	return phi, theta, true
}

// COO — conic orthomorphic (Lambert conformal conic, Paper II §5.4.4).
//
//	R_theta = psi * tan^C(pi/4 - theta/2)
//
// where psi is a constant chosen so that the two standard parallels are
// true scale. Conformal (angle-preserving). Used by aeronautical charts.
type cooProjection struct {
	conicBase
	psi float64
}

func newCOO(pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	b, err := newConicBase(pv, latAxis)
	if err != nil {
		return nil, err
	}
	// Paper II gives:
	//   C = ln(cos(theta_2)/cos(theta_1)) / ln(tan((pi/2 - theta_1)/2) / tan((pi/2 - theta_2)/2))
	//   psi = cos(theta_1) / (C * tan^C((pi/2 - theta_1)/2))
	// For eta = 0 (both parallels coincident at theta_a), the degenerate
	// form is C = sin(theta_a), psi = cos(theta_a) * cot^sin(theta_a)(pi/4 - theta_a/2).
	t1 := b.thetaA - b.eta
	t2 := b.thetaA + b.eta
	var c float64
	if b.eta == 0 {
		c = math.Sin(b.thetaA)
	} else {
		c = math.Log(math.Cos(t2)/math.Cos(t1)) / math.Log(math.Tan((math.Pi/2-t2)/2)/math.Tan((math.Pi/2-t1)/2))
	}
	if c == 0 {
		return nil, fmt.Errorf("wcs/transform: COO cone constant C = 0 (theta_a = 0?)")
	}
	b.C = c
	psi := math.Cos(t1) / (c * math.Pow(math.Tan((math.Pi/2-t1)/2), c))
	// Override Y0 to the radius at theta_a.
	b.Y0 = psi * math.Pow(math.Tan((math.Pi/2-b.thetaA)/2), c)
	return cooProjection{conicBase: b, psi: psi}, nil
}

func (cooProjection) Code() string      { return "COO" }
func (p cooProjection) Theta0() float64 { return p.thetaA }

func (p cooProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	t := math.Tan((math.Pi/2 - theta) / 2)
	if t <= 0 {
		return 0, 0, false
	}
	r := p.psi * math.Pow(t, p.C)
	x, y = p.projectFromR(phi, r)
	return x, y, true
}

func (p cooProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := p.inverseR(x, y)
	phi = p.inversePhi(x, y, r)
	if r <= 0 {
		theta = math.Pi / 2
		return phi, theta, true
	}
	t := math.Pow(r/p.psi, 1/p.C)
	theta = math.Pi/2 - 2*math.Atan(t)
	return phi, theta, true
}
