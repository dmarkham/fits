package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// The zenithal perspective projections (AZP, SZP), the zenithal polynomial
// (ZPN) and Airy's projection (AIR), Paper II §5.1. Ports of the wcslib
// 8.10 prj.c functions of the same names, in degrees (see projection.go).

// AZP — zenithal/azimuthal perspective.
//
//	PV2_1 = mu    distance of the point of projection from the sphere
//	              centre, in sphere radii (default 0)
//	PV2_2 = gamma tilt of the projection plane, degrees (default 0)
//
// Both directions include the tilt exactly (the previous Go inverse was
// an approximation that was exact only along the principal meridians,
// up to 1950 arcsec off). The forward rejects points on the far side of
// the horizon (theta < asin(-1/mu) for |mu| > 1) and, for a tilted plane,
// points beyond the plane's own horizon.
type azpProjection struct {
	mu                         float64
	w0, w1, w2, w3, w4, w5, w7 float64
}

func newAZP(pv map[wcs.PVKey]float64, latAxis int) (azpProjection, error) {
	mu := pvFloat(pv, latAxis, 1, 0)
	gamma := pvFloat(pv, latAxis, 2, 0)
	p := azpProjection{mu: mu}
	p.w0 = float64(r2d * (mu + 1.0))
	if p.w0 == 0.0 {
		return p, errBadParam("AZP", "mu = PV%d_1 = -1", latAxis)
	}
	p.w3 = cosd(gamma)
	if p.w3 == 0.0 {
		return p, errBadParam("AZP", "gamma = PV%d_2 = %g gives cos(gamma) = 0", latAxis, gamma)
	}
	p.w2 = 1.0 / p.w3
	p.w4 = sind(gamma)
	p.w1 = p.w4 / p.w3
	if math.Abs(mu) > 1.0 {
		p.w5 = asind(-1.0 / mu)
	} else {
		p.w5 = -90.0
	}
	w6 := float64(mu * p.w3)
	if math.Abs(w6) < 1.0 {
		p.w7 = 1.0
	}
	return p, nil
}

func (azpProjection) Code() string    { return "AZP" }
func (azpProjection) theta0() float64 { return 90 }

func (p azpProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	yc := float64(y * p.w3)
	q := p.w0 + float64(y*p.w4)
	r := math.Sqrt(float64(x*x) + float64(yc*yc))
	if r == 0.0 {
		phi, theta = 0.0, 90.0
	} else {
		phi = atan2d(x, -yc)
		s := r / q
		t := s * p.mu / math.Sqrt(float64(s*s)+1.0)
		s = atan2d(1.0, s)
		if math.Abs(t) > 1.0 {
			if math.Abs(t) > 1.0+tol {
				return 0, 0, false
			}
			t = math.Copysign(90.0, t)
		} else {
			t = asind(t)
		}
		a := s - t
		b := s + t + 180.0
		if a > 90.0 {
			a -= 360.0
		}
		if b > 90.0 {
			b -= 360.0
		}
		theta = math.Max(a, b)
	}
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p azpProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	sinthe, costhe := sincosd(theta)
	s := float64(p.w1 * cosphi)
	t := (p.mu + sinthe) + float64(costhe*s)
	if t == 0.0 {
		return 0, 0, false
	}
	r := p.w0 * costhe / t
	if theta < p.w5 {
		return 0, 0, false
	} else if p.w7 > 0.0 {
		t = p.mu / math.Sqrt(1.0+float64(s*s))
		if math.Abs(t) <= 1.0 {
			s = atand(-s)
			t = asind(t)
			a := s - t
			b := s + t + 180.0
			if a > 90.0 {
				a -= 360.0
			}
			if b > 90.0 {
				b -= 360.0
			}
			if theta < math.Max(a, b) {
				return 0, 0, false
			}
		}
	}
	return float64(r * sinphi), float64(-r * cosphi * p.w2), true
}

// SZP — slant zenithal perspective.
//
//	PV2_1 = mu      distance of the point of projection from the sphere
//	                centre, in sphere radii (default 0)
//	PV2_2 = phi_c   native longitude of the point of projection (default 0)
//	PV2_3 = theta_c native latitude of the point of projection
//	                (default 90)
//
// theta_c = 90 reduces to AZP with gamma = 0; mu = 0 to TAN.
type szpProjection struct {
	mu                             float64
	w0, w1, w2, w3, w4, w5, w6, w7 float64
	w8                             float64
}

func newSZP(pv map[wcs.PVKey]float64, latAxis int) (szpProjection, error) {
	mu := pvFloat(pv, latAxis, 1, 0)
	phic := pvFloat(pv, latAxis, 2, 0)
	thetac := pvFloat(pv, latAxis, 3, 90)
	p := szpProjection{mu: mu}
	p.w0 = 1.0 / r2d
	p.w3 = float64(mu*sind(thetac)) + 1.0
	if p.w3 == 0.0 {
		return p, errBadParam("SZP", "mu sin(theta_c) = -1")
	}
	p.w1 = float64(-mu * cosd(thetac) * sind(phic))
	p.w2 = float64(mu * cosd(thetac) * cosd(phic))
	p.w4 = float64(r2d * p.w1)
	p.w5 = float64(r2d * p.w2)
	p.w6 = float64(r2d * p.w3)
	p.w7 = float64((p.w3-1.0)*p.w3) - 1.0
	if math.Abs(p.w3-1.0) < 1.0 {
		p.w8 = asind(1.0 - p.w3)
	} else {
		p.w8 = -90.0
	}
	return p, nil
}

func (szpProjection) Code() string    { return "SZP" }
func (szpProjection) theta0() float64 { return 90 }

func (p szpProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	xr := float64(x * p.w0)
	yr := float64(y * p.w0)
	r2 := float64(xr*xr) + float64(yr*yr)
	x1 := (xr - p.w1) / p.w3
	y1 := (yr - p.w2) / p.w3
	xy := float64(xr*x1) + float64(yr*y1)
	var z float64
	if r2 < 1.0e-10 {
		z = r2 / 2.0
		theta = 90.0 - float64(r2d*math.Sqrt(r2/(1.0+xy)))
	} else {
		t := float64(x1*x1) + float64(y1*y1)
		a := t + 1.0
		b := xy - t
		c := r2 - xy - xy + t - 1.0
		d := float64(b*b) - float64(a*c)
		if d < 0.0 {
			return 0, 0, false
		}
		d = math.Sqrt(d)
		sinth1 := (-b + d) / a
		sinth2 := (-b - d) / a
		sinthe := math.Max(sinth1, sinth2)
		if sinthe > 1.0 {
			if sinthe-1.0 < tol {
				sinthe = 1.0
			} else {
				sinthe = math.Min(sinth1, sinth2)
			}
		}
		if sinthe < -1.0 && sinthe+1.0 > -tol {
			sinthe = -1.0
		}
		if sinthe > 1.0 || sinthe < -1.0 {
			return 0, 0, false
		}
		theta = asind(sinthe)
		z = 1.0 - sinthe
	}
	phi = atan2d(xr-float64(x1*z), -(yr - float64(y1*z)))
	return phi, theta, prjbchk(2.0e-13, &phi, &theta)
}

func (p szpProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	s := 1.0 - sind(theta)
	t := p.w3 - s
	if t == 0.0 {
		return 0, 0, false
	}
	r := p.w6 * cosd(theta) / t
	u := p.w4 * s / t
	v := p.w5 * s / t
	if theta < p.w8 {
		return 0, 0, false
	} else if math.Abs(p.mu) > 1.0 {
		s := float64(p.w1*sinphi) - float64(p.w2*cosphi)
		t := 1.0 / math.Sqrt(p.w7+float64(s*s))
		if math.Abs(t) <= 1.0 {
			s = atan2d(s, p.w3-1.0)
			t = asind(t)
			a := s - t
			b := s + t + 180.0
			if a > 90.0 {
				a -= 360.0
			}
			if b > 90.0 {
				b -= 360.0
			}
			if theta < math.Max(a, b) {
				return 0, 0, false
			}
		}
	}
	return float64(r*sinphi) - u, float64(-r*cosphi) - v, true
}

// ZPN — zenithal polynomial.
//
//	R_theta = r0 sum_{m=0}^{29} PV2_m zeta^m,  zeta = (90 - theta) in radians
//
// (all PV2_m default 0; at least one must be non-zero). For degree >= 2
// the projection is defined only out to the first zero of dR/dzeta (found
// as in wcslib's zpnset); beyond it both directions reject. The inverse
// is closed-form for degree 1 and 2 and weighted regula falsi otherwise.
type zpnProjection struct {
	pv     [30]float64
	n      int     // degree: highest non-zero index
	w0, w1 float64 // zeta and R/r0 at the first stationary point
}

func newZPN(pv map[wcs.PVKey]float64, latAxis int) (zpnProjection, error) {
	var p zpnProjection
	for m := 0; m < 30; m++ {
		p.pv[m] = pvFloat(pv, latAxis, m, 0)
	}
	const tol = 1.0e-13
	k := 29
	for k >= 0 && p.pv[k] == 0.0 {
		k--
	}
	if k < 0 {
		return p, errBadParam("ZPN", "every coefficient PV%d_m is zero", latAxis)
	}
	p.n = k
	if k < 2 {
		p.w0 = wpi
		return p, nil
	}
	d1 := p.pv[1]
	if d1 <= 0.0 {
		return p, errBadParam("ZPN", "PV%d_1 = %g must be > 0 for degree %d", latAxis, d1, k)
	}
	var zd, zd1, zd2, d, d2 float64
	j := 0
	for ; j < 180; j++ {
		zd2 = float64(float64(j) * d2r)
		d2 = 0.0
		for m := k; m > 0; m-- {
			d2 = float64(d2*zd2) + float64(float64(m)*p.pv[m])
		}
		if d2 <= 0.0 {
			break
		}
		zd1 = zd2
		d1 = d2
	}
	if j == 180 {
		zd = wpi
	} else {
		for j = 1; j <= 10; j++ {
			zd = zd1 - d1*(zd2-zd1)/(d2-d1)
			d = 0.0
			for m := k; m > 0; m-- {
				d = float64(d*zd) + float64(float64(m)*p.pv[m])
			}
			if math.Abs(d) < tol {
				break
			}
			if d < 0.0 {
				zd2, d2 = zd, d
			} else {
				zd1, d1 = zd, d
			}
		}
	}
	r := 0.0
	for m := k; m >= 0; m-- {
		r = float64(r*zd) + p.pv[m]
	}
	p.w0 = zd
	p.w1 = r
	return p, nil
}

func (zpnProjection) Code() string    { return "ZPN" }
func (zpnProjection) theta0() float64 { return 90 }

func (p zpnProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	k := p.n
	r := math.Sqrt(float64(x*x)+float64(y*y)) / r2d
	if r != 0.0 {
		phi = atan2d(x, -y)
	}
	var zd float64
	switch {
	case k < 1:
		return 0, 0, false
	case k == 1:
		zd = (r - p.pv[0]) / p.pv[1]
	case k == 2:
		a := p.pv[2]
		b := p.pv[1]
		c := p.pv[0] - r
		d := float64(b*b) - float64(4.0*a*c)
		if d < 0.0 {
			return 0, 0, false
		}
		d = math.Sqrt(d)
		zd1 := (-b + d) / (2.0 * a)
		zd2 := (-b - d) / (2.0 * a)
		zd = math.Min(zd1, zd2)
		if zd < -tol {
			zd = math.Max(zd1, zd2)
		}
		if zd < 0.0 {
			if zd < -tol {
				return 0, 0, false
			}
			zd = 0.0
		} else if zd > wpi {
			if zd > wpi+tol {
				return 0, 0, false
			}
			zd = wpi
		}
	default:
		zd1 := 0.0
		r1 := p.pv[0]
		zd2 := p.w0
		r2 := p.w1
		if r < r1 {
			if r < r1-tol {
				return 0, 0, false
			}
			zd = zd1
		} else if r > r2 {
			if r > r2+tol {
				return 0, 0, false
			}
			zd = zd2
		} else {
			for j := 0; j < 100; j++ {
				lambda := (r2 - r) / (r2 - r1)
				if lambda < 0.1 {
					lambda = 0.1
				} else if lambda > 0.9 {
					lambda = 0.9
				}
				zd = zd2 - float64(lambda*(zd2-zd1))
				rt := 0.0
				for m := k; m >= 0; m-- {
					rt = (float64(rt * zd)) + p.pv[m]
				}
				if rt < r {
					if r-rt < tol {
						break
					}
					r1, zd1 = rt, zd
				} else {
					if rt-r < tol {
						break
					}
					r2, zd2 = rt, zd
				}
				if math.Abs(zd2-zd1) < tol {
					break
				}
			}
		}
	}
	theta = 90.0 - float64(zd*r2d)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p zpnProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	s := float64((90.0 - theta) * d2r)
	r := 0.0
	for m := p.n; m >= 0; m-- {
		r = float64(r*s) + p.pv[m]
	}
	r *= r2d
	if s > p.w0 {
		return 0, 0, false
	}
	return float64(r * sinphi), float64(-r * cosphi), true
}

// AIR — Airy's zenithal. PV2_1 = theta_b, the latitude at which
// the error is minimized, degrees (default 90). The inverse has no closed
// form; it is solved by bracketing and weighted regula falsi as in wcslib.
type airProjection struct {
	w0, w1, w3, w4, w5, w6 float64
}

func newAIR(pv map[wcs.PVKey]float64, latAxis int) (airProjection, error) {
	const tol = 1.0e-4
	thetab := pvFloat(pv, latAxis, 1, 90)
	var p airProjection
	p.w0 = float64(2.0 * r2d)
	var w2 float64
	if thetab == 90.0 {
		p.w1 = -0.5
		w2 = 1.0
	} else if thetab > -90.0 {
		cosxi := cosd((90.0 - thetab) / 2.0)
		p.w1 = lmLog(cosxi) * (cosxi * cosxi) / (1.0 - float64(cosxi*cosxi))
		w2 = 0.5 - p.w1
	} else {
		return p, errBadParam("AIR", "theta_b = PV%d_1 = %g must be > -90", latAxis, thetab)
	}
	p.w3 = float64(p.w0 * w2)
	p.w4 = tol
	p.w5 = float64(w2 * tol)
	p.w6 = r2d / w2
	return p, nil
}

func (airProjection) Code() string    { return "AIR" }
func (airProjection) theta0() float64 { return 90 }

func (p airProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	r := math.Sqrt(float64(x*x)+float64(y*y)) / p.w0
	if r != 0.0 {
		phi = atan2d(x, -y)
	}
	var xi float64
	if r == 0.0 {
		xi = 0.0
	} else if r < p.w5 {
		xi = float64(r * p.w6)
	} else {
		x1, x2 := 1.0, 1.0
		r1, r2 := 0.0, 0.0
		k := 0
		for ; k < 30; k++ {
			x2 = x1 / 2.0
			tanxi := math.Sqrt(1.0-float64(x2*x2)) / x2
			r2 = -(lmLog(x2)/tanxi + float64(p.w1*tanxi))
			if r2 >= r {
				break
			}
			x1 = x2
			r1 = r2
		}
		// Unreachable as a distinct decision (see "Guards that no input
		// reaches" in projection.go); kept because wcslib has it.
		if k == 30 {
			return 0, 0, false
		}
		var cosxi float64
		for k = 0; k < 100; k++ {
			lambda := (r2 - r) / (r2 - r1)
			if lambda < 0.1 {
				lambda = 0.1
			} else if lambda > 0.9 {
				lambda = 0.9
			}
			cosxi = x2 - float64(lambda*(x2-x1))
			tanxi := math.Sqrt(1.0-float64(cosxi*cosxi)) / cosxi
			rt := -(lmLog(cosxi)/tanxi + float64(p.w1*tanxi))
			if rt < r {
				if r-rt < tol {
					break
				}
				r1 = rt
				x1 = cosxi
			} else {
				if rt-r < tol {
					break
				}
				r2 = rt
				x2 = cosxi
			}
		}
		if k == 100 {
			return 0, 0, false
		}
		xi = acosd(cosxi)
	}
	theta = 90.0 - float64(2.0*xi)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

func (p airProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	var r float64
	if theta == 90.0 {
		r = 0.0
	} else if theta > -90.0 {
		xi := d2r * (90.0 - theta) / 2.0
		if xi < p.w4 {
			r = float64(xi * p.w3)
		} else {
			cosxi := cosd((90.0 - theta) / 2.0)
			tanxi := math.Sqrt(1.0-float64(cosxi*cosxi)) / cosxi
			r = float64(-p.w0 * (lmLog(cosxi)/tanxi + float64(p.w1*tanxi)))
		}
	} else {
		return 0, 0, false
	}
	return float64(r * sinphi), float64(-r * cosphi), true
}
