package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file implements the zenithal (azimuthal) projection family from
// Calabretta & Greisen 2002 Paper II §5.1.
//
// Every zenithal projection has its reference (fiducial) point at the
// native north pole (theta_0 = 90°). The forward mapping from native
// (phi, theta) to intermediate (x, y) is
//
//	x = R_theta * sin(phi)
//	y = -R_theta * cos(phi)
//
// where R_theta is a projection-specific function of theta that describes
// how colatitude (90° - theta) is scaled into the plane. The inverse is
//
//	R_theta = sqrt(x^2 + y^2)
//	phi     = atan2(x, -y)
//	theta   = inverse-R(R_theta)
//
// All calculations in this file are in radians. When Paper II writes its
// formulas with the "(180/pi)" prefactor it is converting degrees to
// radians in a unit-naive context; our formulas drop that prefactor
// because every angle is already in radians throughout.

// TAN — gnomonic projection (§5.1.3).
//
//	R_theta = cot(theta) = cos(theta) / sin(theta)
//
// Valid for theta > 0 (the near hemisphere). Points on or south of the
// equator are outside the projection's image.
type tanProjection struct{}

func (tanProjection) Code() string    { return "TAN" }
func (tanProjection) Theta0() float64 { return math.Pi / 2 }

func (tanProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	s := math.Sin(theta)
	if s <= 0 {
		// theta ≤ 0: the point is on the far hemisphere, infinitely far on TAN.
		return 0, 0, false
	}
	r := math.Cos(theta) / s
	return r * math.Sin(phi), -r * math.Cos(phi), true
}

func (tanProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := math.Hypot(x, y)
	if r == 0 {
		return 0, math.Pi / 2, true
	}
	phi = math.Atan2(x, -y)
	theta = math.Atan2(1, r) // theta = atan(1/R) since R = cot(theta)
	return phi, theta, true
}

// SIN — orthographic / slant orthographic projection (§5.1.5).
//
// Plain SIN (xi = eta = 0):
//
//	x = cos(theta) * sin(phi)
//	y = -cos(theta) * cos(phi)
//	R_theta = cos(theta), valid for theta >= 0
//
// Slant SIN (the generalized form, Paper II eq. 43 with PV2_1 = xi and
// PV2_2 = eta):
//
//	x = cos(theta)*sin(phi) + xi*(1 - sin(theta))
//	y = -cos(theta)*cos(phi) + eta*(1 - sin(theta))
//
// The slant form is also known as NCP ("north celestial pole") when
// xi = 0 and eta = cot(delta_0); this is the historical radio-
// interferometry projection convention still used by many VLA/ALMA data
// products.
//
// The inverse of the slant form requires solving a quadratic; we pick the
// root that lies in the near hemisphere (theta >= 0).
type sinProjection struct {
	xi  float64
	eta float64
}

func newSIN(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return sinProjection{
		xi:  pvFloat(pv, latAxis, 1, 0),
		eta: pvFloat(pv, latAxis, 2, 0),
	}
}

func (sinProjection) Code() string    { return "SIN" }
func (sinProjection) Theta0() float64 { return math.Pi / 2 }

func (p sinProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if theta < 0 {
		return 0, 0, false
	}
	sinTheta, cosTheta := math.Sincos(theta)
	sinPhi, cosPhi := math.Sincos(phi)
	x = cosTheta*sinPhi + p.xi*(1-sinTheta)
	y = -cosTheta*cosPhi + p.eta*(1-sinTheta)
	return x, y, true
}

func (p sinProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	if p.xi == 0 && p.eta == 0 {
		// Plain SIN fast path: closed-form inverse.
		r := math.Hypot(x, y)
		if r > 1 {
			return 0, 0, false
		}
		if r == 0 {
			return 0, math.Pi / 2, true
		}
		phi = math.Atan2(x, -y)
		theta = math.Acos(r)
		return phi, theta, true
	}
	// Slant SIN: let s = 1 - sin(theta). Then
	//   x - xi*s = cos(theta)*sin(phi)
	//   y - eta*s = -cos(theta)*cos(phi)
	// Squaring and adding:
	//   (x - xi*s)² + (y - eta*s)² = cos²(theta) = 1 - sin²(theta) = (1-sin(theta))*(1+sin(theta)) = s*(2 - s)
	// Expand:
	//   x² - 2*xi*x*s + xi²*s² + y² - 2*eta*y*s + eta²*s² = 2*s - s²
	//   (xi² + eta² + 1)*s² - 2*(xi*x + eta*y + 1)*s + (x² + y²) = 0
	a := p.xi*p.xi + p.eta*p.eta + 1
	b := -2 * (p.xi*x + p.eta*y + 1)
	c := x*x + y*y
	disc := b*b - 4*a*c
	if disc < 0 {
		return 0, 0, false
	}
	sqrtDisc := math.Sqrt(disc)
	// Choose the smaller root (corresponds to larger sin(theta), i.e. near
	// hemisphere).
	s := (-b - sqrtDisc) / (2 * a)
	if s < 0 || s > 2 {
		s = (-b + sqrtDisc) / (2 * a)
	}
	if s < 0 || s > 2 {
		return 0, 0, false
	}
	sinTheta := 1 - s
	cosTheta := math.Sqrt((1 - sinTheta) * (1 + sinTheta))
	// Back-substitute to get phi.
	sinPhi := (x - p.xi*s)
	cosPhi := -(y - p.eta*s)
	phi = math.Atan2(sinPhi, cosPhi)
	theta = math.Atan2(sinTheta, cosTheta)
	return phi, theta, true
}

// STG — stereographic projection (§5.1.4).
//
//	R_theta = 2 * tan((90° - theta) / 2) = 2 * (1 - sin(theta)) / cos(theta)
//
// Alternative form using the half-angle identity:
//
//	R_theta = 2 * cos(theta) / (1 + sin(theta))
//
// Valid for all theta > -90° (the entire sphere minus the native south
// pole, which maps to infinity). STG is the unique zenithal projection
// that is conformal (angle-preserving).
type stgProjection struct{}

func (stgProjection) Code() string    { return "STG" }
func (stgProjection) Theta0() float64 { return math.Pi / 2 }

func (stgProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	s := math.Sin(theta)
	if 1+s == 0 {
		return 0, 0, false
	}
	r := 2 * math.Cos(theta) / (1 + s)
	return r * math.Sin(phi), -r * math.Cos(phi), true
}

func (stgProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := math.Hypot(x, y)
	if r == 0 {
		return 0, math.Pi / 2, true
	}
	phi = math.Atan2(x, -y)
	// From R = 2 cos(theta) / (1 + sin(theta)) we get
	//   theta = pi/2 - 2 * atan(R/2)
	theta = math.Pi/2 - 2*math.Atan(r/2)
	return phi, theta, true
}

// ARC — zenithal equidistant projection (§5.1.6).
//
//	R_theta = (90° - theta) measured in the same unit (radians here)
//
// Distance along meridians is preserved. Valid for the entire sphere.
type arcProjection struct{}

func (arcProjection) Code() string    { return "ARC" }
func (arcProjection) Theta0() float64 { return math.Pi / 2 }

func (arcProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	r := math.Pi/2 - theta
	return r * math.Sin(phi), -r * math.Cos(phi), true
}

func (arcProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := math.Hypot(x, y)
	if r > math.Pi {
		// Outside the 180° disk that covers the full sphere.
		return 0, 0, false
	}
	if r == 0 {
		return 0, math.Pi / 2, true
	}
	phi = math.Atan2(x, -y)
	theta = math.Pi/2 - r
	return phi, theta, true
}

// ZEA — zenithal equal-area projection (§5.1.7).
//
//	R_theta = 2 * sin((90° - theta) / 2)
//
// Areas are preserved. Valid for the entire sphere; the native south pole
// maps to R_theta = 2.
type zeaProjection struct{}

func (zeaProjection) Code() string    { return "ZEA" }
func (zeaProjection) Theta0() float64 { return math.Pi / 2 }

func (zeaProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	// R = 2 sin((pi/2 - theta)/2) = 2 sin(pi/4 - theta/2)
	r := 2 * math.Sin(math.Pi/4-theta/2)
	return r * math.Sin(phi), -r * math.Cos(phi), true
}

func (zeaProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := math.Hypot(x, y)
	if r > 2 {
		return 0, 0, false
	}
	if r == 0 {
		return 0, math.Pi / 2, true
	}
	phi = math.Atan2(x, -y)
	// theta = pi/2 - 2*asin(R/2)
	theta = math.Pi/2 - 2*math.Asin(r/2)
	return phi, theta, true
}
