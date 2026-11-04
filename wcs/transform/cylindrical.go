package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Cylindrical projections, Paper II §5.2: CYP, CEA, CAR, MER. All have
// theta_0 = 0. Ports of the wcslib 8.10 prj.c functions of the same names,
// in degrees (see projection.go). The inverse rejects any plane point whose
// native longitude would fall outside [-180, 180] or latitude outside
// [-90, 90] (wcslib's prjbchk), i.e. every point off the map.

// CAR — plate carrée: x = phi, y = theta.
type carProjection struct{}

func (carProjection) Code() string    { return "CAR" }
func (carProjection) theta0() float64 { return 0 }

func (carProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	return phi, theta, true
}

func (carProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	phi, theta = x, y
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// CYP — cylindrical perspective.
//
//	PV2_1 = mu     (default 1)
//	PV2_2 = lambda (default 1)
//
//	x = lambda phi
//	y = r0 (mu + lambda) sin(theta) / (mu + cos(theta))
//
// (the previous Go code used (mu + 1) for (mu + lambda)). Gall's
// stereographic is mu = 1, lambda = sqrt(2)/2; the central cylindrical
// mu = 0, lambda = 1.
type cypProjection struct {
	mu             float64
	w0, w1, w2, w3 float64
}

func newCYP(pv map[wcs.PVKey]float64, latAxis int) (cypProjection, error) {
	p := cypProjection{mu: pvFloat(pv, latAxis, 1, 1)}
	lambda := pvFloat(pv, latAxis, 2, 1)
	p.w0 = lambda
	if p.w0 == 0.0 {
		return p, errBadParam("CYP", "lambda = PV%d_2 = 0", latAxis)
	}
	p.w1 = 1.0 / p.w0
	p.w2 = float64(r2d * (p.mu + lambda))
	if p.w2 == 0.0 {
		return p, errBadParam("CYP", "mu + lambda = 0")
	}
	p.w3 = 1.0 / p.w2
	return p, nil
}

func (cypProjection) Code() string    { return "CYP" }
func (cypProjection) theta0() float64 { return 0 }

func (p cypProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	eta := p.mu + cosd(theta)
	if eta == 0.0 {
		return 0, 0, false
	}
	return float64(p.w0 * phi), p.w2 * sind(theta) / eta, true
}

func (p cypProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	phi = float64(p.w1 * x)
	eta := float64(p.w3 * y)
	theta = atan2d(eta, 1.0) + asind(eta*p.mu/math.Sqrt(float64(eta*eta)+1.0))
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// CEA — cylindrical equal area. PV2_1 = lambda in (0, 1]
// (default 1, Lambert's projection): x = phi, y = r0 sin(theta) / lambda.
type ceaProjection struct {
	w2, w3 float64
}

func newCEA(pv map[wcs.PVKey]float64, latAxis int) (ceaProjection, error) {
	lambda := pvFloat(pv, latAxis, 1, 1)
	if lambda <= 0.0 || lambda > 1.0 {
		return ceaProjection{}, errBadParam("CEA", "lambda = PV%d_1 = %g not in (0, 1]", latAxis, lambda)
	}
	return ceaProjection{w2: r2d / lambda, w3: lambda / r2d}, nil
}

func (ceaProjection) Code() string    { return "CEA" }
func (ceaProjection) theta0() float64 { return 0 }

func (p ceaProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	return phi, float64(p.w2 * sind(theta)), true
}

func (p ceaProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	phi = x
	s := float64(p.w3 * y)
	if math.Abs(s) > 1.0 {
		if math.Abs(s) > 1.0+tol {
			return 0, 0, false
		}
		theta = math.Copysign(90.0, s)
	} else {
		theta = asind(s)
	}
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// MER — Mercator: x = phi, y = r0 ln(tan(45 + theta/2)). The
// poles are at infinity and are rejected.
type merProjection struct{}

func (merProjection) Code() string    { return "MER" }
func (merProjection) theta0() float64 { return 0 }

func (merProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	if theta <= -90.0 || theta >= 90.0 {
		return 0, 0, false
	}
	return phi, float64(r2d * lmLog(tand((theta+90.0)/2.0))), true
}

func (merProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	phi = x
	theta = float64(2.0*atand(lmExp(y/r2d))) - 90.0
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}
