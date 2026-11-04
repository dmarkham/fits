package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Conic projections, Paper II §5.4: COP, COE, COD, COO. All four take
//
//	PV2_1 = theta_a (REQUIRED; Paper II: no default) — the mean of the two
//	        standard parallels, degrees
//	PV2_2 = eta     (default 0) — half their separation, degrees
//
// so the standard parallels are theta_1,2 = theta_a -+ eta, and the
// fiducial point is (phi_0, theta_0) = (0, theta_a). With C the cone
// constant and Y_0 the plane offset that puts the fiducial point at the
// origin,
//
//	x = R_theta sin(C phi),  y = -R_theta cos(C phi) + Y_0
//
// Ports of the wcslib 8.10 prj.c functions of the same names, in degrees
// (see projection.go). Paper II's formulas, which these implement:
//
//	COP: C = sin(theta_a), R = r0 cos(eta) [cot(theta_a) - tan(theta - theta_a)],
//	     Y_0 = r0 cos(eta) cot(theta_a) (the previous Go code dropped
//	     cos(eta) from Y_0)
//	COE: C = (sin theta_1 + sin theta_2)/2,
//	     R = (r0/C) sqrt(1 + sin theta_1 sin theta_2 - 2C sin theta)
//	COD: C = sin(theta_a) sin(eta)/eta, eta in radians (C = sin theta_a
//	     for eta = 0), R = theta_a - theta + eta cot(eta) cot(theta_a)
//	     in degrees (the previous Go code dropped eta cot(eta))
//	COO: C = ln(cos theta_2/cos theta_1) / ln(tan((90-theta_2)/2)/tan((90-theta_1)/2)),
//	     R = psi tan^C((90 - theta)/2), psi = r0 cos(theta_1)/(C tan^C((90-theta_1)/2))

func conicThetaA(code string, pv map[wcs.PVKey]float64, latAxis int) (thetaA, eta float64, err error) {
	if !pvSet(pv, latAxis, 1) {
		return 0, 0, errBadParam(code, "PV%d_1 (theta_a) is required", latAxis)
	}
	return pvFloat(pv, latAxis, 1, 0), pvFloat(pv, latAxis, 2, 0), nil
}

// conicAlpha returns the cone angle alpha (C phi in the forward formula)
// of plane point (x, dy) where dy = Y_0 - y, and the signed radius.
func conicAlpha(x, dy float64, southern bool) (alpha, r float64) {
	r = math.Sqrt(float64(x*x) + float64(dy*dy))
	if southern {
		r = -r
	}
	if r != 0.0 {
		alpha = atan2d(x/r, dy/r)
	}
	return alpha, r
}

// COP — conic perspective.
type copProjection struct {
	thetaA                 float64
	w0, w1, w2, w3, w4, w5 float64
}

func newCOP(pv map[wcs.PVKey]float64, latAxis int) (copProjection, error) {
	ta, eta, err := conicThetaA("COP", pv, latAxis)
	if err != nil {
		return copProjection{}, err
	}
	p := copProjection{thetaA: ta}
	p.w0 = sind(ta)
	if p.w0 == 0.0 {
		return p, errBadParam("COP", "theta_a = 0")
	}
	p.w1 = 1.0 / p.w0
	p.w3 = float64(r2d * cosd(eta))
	if p.w3 == 0.0 {
		return p, errBadParam("COP", "cos(eta) = 0")
	}
	p.w4 = 1.0 / p.w3
	p.w5 = 1.0 / tand(ta)
	p.w2 = float64(p.w3 * p.w5)
	return p, nil
}

func (copProjection) Code() string      { return "COP" }
func (p copProjection) theta0() float64 { return p.thetaA }

func (p copProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	alpha, r := conicAlpha(x, p.w2-y, p.thetaA < 0.0)
	phi = float64(alpha * p.w1)
	theta = p.thetaA + atand(p.w5-float64(r*p.w4))
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p copProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinalpha, cosalpha := sincosd(float64(p.w0 * phi))
	t := theta - p.thetaA
	s := cosd(t)
	var r float64
	if s == 0.0 {
		return 0, 0, false
	} else if math.Abs(theta) == 90.0 {
		r = 0.0
		if (theta < 0.0) != (p.thetaA < 0.0) {
			return 0, 0, false
		}
	} else {
		r = p.w2 - p.w3*sind(t)/s
		if float64(r*p.w0) < 0.0 {
			return 0, 0, false
		}
	}
	return float64(r * sinalpha), float64(-r*cosalpha) + p.w2, true
}

// COE — conic equal area.
type coeProjection struct {
	thetaA                         float64
	w0, w1, w2, w3, w4, w5, w6, w7 float64
	w8                             float64
}

func newCOE(pv map[wcs.PVKey]float64, latAxis int) (coeProjection, error) {
	ta, eta, err := conicThetaA("COE", pv, latAxis)
	if err != nil {
		return coeProjection{}, err
	}
	p := coeProjection{thetaA: ta}
	theta1 := ta - eta
	theta2 := ta + eta
	p.w0 = (sind(theta1) + sind(theta2)) / 2.0
	if p.w0 == 0.0 {
		return p, errBadParam("COE", "sin(theta_1) + sin(theta_2) = 0")
	}
	p.w1 = 1.0 / p.w0
	p.w3 = r2d / p.w0
	p.w4 = 1.0 + float64(sind(theta1)*sind(theta2))
	p.w5 = float64(2.0 * p.w0)
	p.w6 = float64(p.w3 * p.w3 * p.w4)
	p.w7 = 1.0 / (2.0 * r2d * p.w3)
	p.w8 = float64(p.w3 * math.Sqrt(p.w4+p.w5))
	p.w2 = float64(p.w3 * math.Sqrt(p.w4-float64(p.w5*sind(ta))))
	return p, nil
}

func (coeProjection) Code() string      { return "COE" }
func (p coeProjection) theta0() float64 { return p.thetaA }

func (p coeProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	alpha, r := conicAlpha(x, p.w2-y, p.thetaA < 0.0)
	if math.Abs(r-p.w8) < tol {
		theta = -90.0
	} else {
		w := float64((p.w6 - float64(r*r)) * p.w7)
		if math.Abs(w) > 1.0 {
			if math.Abs(w-1.0) < tol {
				theta = 90.0
			} else if math.Abs(w+1.0) < tol {
				theta = -90.0
			} else {
				return 0, 0, false
			}
		} else {
			theta = asind(w)
		}
	}
	phi = float64(alpha * p.w1)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p coeProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinalpha, cosalpha := sincosd(float64(p.w0 * phi))
	var r float64
	if theta == -90.0 {
		r = p.w8
	} else {
		r = float64(p.w3 * math.Sqrt(p.w4-float64(p.w5*sind(theta))))
	}
	return float64(r * sinalpha), float64(-r*cosalpha) + p.w2, true
}

// COD — conic equidistant.
type codProjection struct {
	thetaA         float64
	w0, w1, w2, w3 float64
}

func newCOD(pv map[wcs.PVKey]float64, latAxis int) (codProjection, error) {
	ta, eta, err := conicThetaA("COD", pv, latAxis)
	if err != nil {
		return codProjection{}, err
	}
	p := codProjection{thetaA: ta}
	if eta == 0.0 {
		p.w0 = float64(r2d * sind(ta) * d2r)
	} else {
		p.w0 = r2d * sind(ta) * sind(eta) / eta
	}
	if p.w0 == 0.0 {
		return p, errBadParam("COD", "theta_a = 0")
	}
	p.w1 = 1.0 / p.w0
	p.w2 = r2d * cosd(eta) * cosd(ta) / p.w0
	p.w3 = p.w2 + ta
	return p, nil
}

func (codProjection) Code() string      { return "COD" }
func (p codProjection) theta0() float64 { return p.thetaA }

func (p codProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	alpha, r := conicAlpha(x, p.w2-y, p.thetaA < 0.0)
	phi = float64(alpha * p.w1)
	theta = p.w3 - r
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p codProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinalpha, cosalpha := sincosd(float64(p.w0 * phi))
	r := p.w3 - theta
	return float64(r * sinalpha), float64(-r*cosalpha) + p.w2, true
}

// COO — conic orthomorphic (Lambert conformal conic). The apex of
// the cone (theta = -90 for a northern cone, +90 for a southern one where
// wcslib's own forward gives infinity) is not on the projection.
type cooProjection struct {
	thetaA             float64
	w0, w1, w2, w3, w4 float64
}

func newCOO(pv map[wcs.PVKey]float64, latAxis int) (cooProjection, error) {
	ta, eta, err := conicThetaA("COO", pv, latAxis)
	if err != nil {
		return cooProjection{}, err
	}
	p := cooProjection{thetaA: ta}
	theta1 := ta - eta
	theta2 := ta + eta
	tan1 := tand((90.0 - theta1) / 2.0)
	cos1 := cosd(theta1)
	if theta1 == theta2 {
		p.w0 = sind(theta1)
	} else {
		tan2 := tand((90.0 - theta2) / 2.0)
		cos2 := cosd(theta2)
		p.w0 = lmLog(cos2/cos1) / lmLog(tan2/tan1)
	}
	if p.w0 == 0.0 {
		return p, errBadParam("COO", "cone constant C = 0")
	}
	p.w1 = 1.0 / p.w0
	p.w3 = r2d * (cos1 / p.w0) / lmPow(tan1, p.w0)
	if p.w3 == 0.0 {
		return p, errBadParam("COO", "psi = 0")
	}
	p.w2 = float64(p.w3 * lmPow(tand((90.0-ta)/2.0), p.w0))
	p.w4 = 1.0 / p.w3
	return p, nil
}

func (cooProjection) Code() string      { return "COO" }
func (p cooProjection) theta0() float64 { return p.thetaA }

func (p cooProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	alpha, r := conicAlpha(x, p.w2-y, p.thetaA < 0.0)
	if r == 0.0 {
		if p.w0 < 0.0 {
			theta = -90.0
		} else {
			return 0, 0, false
		}
	} else {
		theta = 90.0 - float64(2.0*atand(lmPow(float64(r*p.w4), p.w1)))
	}
	phi = float64(alpha * p.w1)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p cooProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinalpha, cosalpha := sincosd(float64(p.w0 * phi))
	var r float64
	if theta == -90.0 {
		r = 0.0
		if p.w0 >= 0.0 {
			return 0, 0, false
		}
	} else {
		r = float64(p.w3 * lmPow(tand((90.0-theta)/2.0), p.w0))
	}
	return float64(r * sinalpha), float64(-r*cosalpha) + p.w2, true
}
