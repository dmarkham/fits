package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Polyconic and pseudoconic projections, Paper II §5.5: BON, PCO. Both have
// theta_0 = 0. Ports of the wcslib 8.10 prj.c functions of the same names,
// in degrees (see projection.go).

// BON — Bonne's equal area. PV2_1 = theta_1, the standard
// parallel, degrees (REQUIRED). theta_1 = 0 is Sanson-Flamsteed (SFL).
//
//	Y_0 = r0 cot(theta_1) + theta_1,   R_theta = Y_0 - theta
//	A   = r0 phi cos(theta) / R_theta
//	x   = R_theta sin(A),   y = -R_theta cos(A) + Y_0
//
// The fiducial point is (phi_0, theta_0) = (0, 0) whatever theta_1 is
// (wcslib bonset; checked against wcslib's prjset theta0 in
// TestWrefProjections), not (0, theta_1).
type bonProjection struct {
	theta1 float64
	w2     float64 // Y_0
}

func newBON(pv map[wcs.PVKey]float64, latAxis int) (bonProjection, error) {
	if !pvSet(pv, latAxis, 1) {
		return bonProjection{}, errBadParam("BON", "PV%d_1 (theta_1) is required", latAxis)
	}
	t1 := pvFloat(pv, latAxis, 1, 0)
	p := bonProjection{theta1: t1}
	if t1 != 0.0 {
		p.w2 = r2d*cosd(t1)/sind(t1) + t1
	}
	return p, nil
}

func (bonProjection) Code() string    { return "BON" }
func (bonProjection) theta0() float64 { return 0 }

func (p bonProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	if p.theta1 == 0.0 {
		return sflProjection{}.x2s(x, y)
	}
	alpha, r := conicAlpha(x, p.w2-y, p.theta1 < 0.0)
	t := p.w2 - r
	costhe := cosd(t)
	if costhe != 0.0 {
		phi = alpha * (r / r2d) / costhe
	}
	theta = t
	return phi, theta, prjbchk(1.0e-11, &phi, &theta)
}

func (p bonProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	if p.theta1 == 0.0 {
		return sflProjection{}.s2x(phi, theta)
	}
	r := p.w2 - theta
	s := cosd(theta) / r
	sinalpha, cosalpha := sincosd(float64(s * (r2d * phi)))
	return float64(r * sinalpha), float64(-r*cosalpha) + p.w2, true
}

// PCO — polyconic:
//
//	theta = 0: x = phi, y = 0
//	otherwise: E = phi sin(theta), x = r0 cot(theta) sin(E),
//	           y = r0 [cot(theta)(1 - cos(E)) + theta]
//
// The inverse has no closed form; it is solved by weighted regula falsi as
// in wcslib's pcox2s. Plane points off the map are rejected (or, where
// wcslib's iteration produces NaN, refused by degAdapter).
type pcoProjection struct{}

var (
	pcoW2 = 360.0 / wpi
	pcoW3 = d2r / pcoW2
)

func (pcoProjection) Code() string    { return "PCO" }
func (pcoProjection) theta0() float64 { return 0 }

func (pcoProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	w := math.Abs(y)
	switch {
	case w < tol:
		phi, theta = x, 0.0
	case math.Abs(w-90.0) < tol:
		phi, theta = 0.0, math.Copysign(90.0, y)
	default:
		var the, ymthe, tanthe float64
		if w < 1.0e-4 {
			the = y / (1.0 + float64(pcoW3*x*x))
			ymthe = y - the
			tanthe = tand(the)
		} else {
			thepos := y
			theneg := 0.0
			xx := float64(x * x)
			fpos := xx
			fneg := -xx
			for k := 0; k < 64; k++ {
				lambda := fpos / (fpos - fneg)
				if lambda < 0.1 {
					lambda = 0.1
				} else if lambda > 0.9 {
					lambda = 0.9
				}
				the = thepos - float64(lambda*(thepos-theneg))
				ymthe = y - the
				tanthe = tand(the)
				f := xx + float64(ymthe*(ymthe-pcoW2/tanthe))
				if math.Abs(f) < tol {
					break
				}
				if math.Abs(thepos-theneg) < tol {
					break
				}
				if f > 0.0 {
					thepos, fpos = the, f
				} else {
					theneg, fneg = the, f
				}
			}
		}
		x1 := r2d - float64(ymthe*tanthe)
		y1 := float64(x * tanthe)
		if x1 != 0.0 || y1 != 0.0 {
			phi = atan2d(y1, x1) / sind(the)
		}
		theta = the
	}
	return phi, theta, prjbchk(1.0e-12, &phi, &theta)
}

func (pcoProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	switch {
	case theta == 0.0:
		// wcslib returns y = -0.0 here. In Go the constant -0.0 is +0
		// (constants have no negative zero), so the sign is set explicitly;
		// TestPCOSmallAngleBand pins the bits.
		return phi, math.Copysign(0, -1), true
	case math.Abs(theta) < 1.0e-4:
		x = float64(phi * cosd(theta))
		y = float64((1.0 + float64(pcoW3*x*x)) * theta)
		return x, y, true
	}
	therad := float64(theta * d2r)
	sinthe, costhe := sincosd(theta)
	sinpsi, cospsi := sincosd(float64(phi * sinthe))
	cotthe := costhe / sinthe
	return float64(r2d * cotthe * sinpsi), float64(r2d * (float64(cotthe*(1.0-cospsi)) + therad)), true
}
