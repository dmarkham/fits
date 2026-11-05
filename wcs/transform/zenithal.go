package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Zenithal (azimuthal) projections, Calabretta & Greisen 2002 (Paper II)
// §5.1: TAN, STG, SIN, ARC, ZEA here; AZP, SZP, ZPN and AIR in
// zenithal_extra.go.
//
// Every zenithal projection has its fiducial point at the native pole
// (theta_0 = 90). With R_theta a projection-specific function of theta,
//
//	x =  R_theta sin(phi)
//	y = -R_theta cos(phi)
//
// Each s2x / x2s below is a port of the wcslib 8.10 prj.c function of the
// same name, in wcslib's degree units (r0 = 180/pi); see projection.go.

// TAN — gnomonic: R_theta = r0 cot(theta). Only the hemisphere
// theta > 0 is on the projection.
type tanProjection struct{}

func (tanProjection) Code() string    { return "TAN" }
func (tanProjection) theta0() float64 { return 90 }

func (tanProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	s := sind(theta)
	if s <= 0.0 {
		return 0, 0, false
	}
	r := r2d * cosd(theta) / s
	return float64(r * sinphi), float64(-r * cosphi), true
}

func (tanProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	r := math.Sqrt(float64(x*x) + float64(y*y))
	if r != 0.0 {
		phi = atan2d(x, -y)
	}
	theta = atan2d(r2d, r)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// STG — stereographic: R_theta = 2 r0 tan((90 - theta)/2). The
// whole sphere except the native south pole.
type stgProjection struct{}

func (stgProjection) Code() string    { return "STG" }
func (stgProjection) theta0() float64 { return 90 }

func (stgProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	s := 1.0 + sind(theta)
	if s == 0.0 {
		return 0, 0, false
	}
	r := (360.0 / wpi) * cosd(theta) / s
	return float64(r * sinphi), float64(-r * cosphi), true
}

func (stgProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	r := math.Sqrt(float64(x*x) + float64(y*y))
	if r != 0.0 {
		phi = atan2d(x, -y)
	}
	theta = 90.0 - float64(2.0*atand(float64(r*(wpi/360.0))))
	return phi, theta, true
}

// SIN — orthographic / synthesis, with Paper II's slant
// generalization (PV2_1 = xi, PV2_2 = eta, default 0):
//
//	x =  r0 cos(theta) sin(phi) + r0 xi  (1 - sin(theta))
//	y = -r0 cos(theta) cos(phi) + r0 eta (1 - sin(theta))
//
// NCP is the special case xi = 0, eta = cot(delta_0). A point is on the
// projection when it is on the visible side of the slanted plane
// (theta >= -atan(xi sin(phi) - eta cos(phi))); for plain SIN that is
// theta >= 0.
type sinProjection struct {
	xi, eta float64
	w1      float64 // xi^2 + eta^2
}

func newSIN(pv map[wcs.PVKey]float64, latAxis int) sinProjection {
	xi := pvFloat(pv, latAxis, 1, 0)
	eta := pvFloat(pv, latAxis, 2, 0)
	return sinProjection{xi: xi, eta: eta, w1: float64(xi*xi) + float64(eta*eta)}
}

func (sinProjection) Code() string    { return "SIN" }
func (sinProjection) theta0() float64 { return 90 }

func (p sinProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	var costhe, z float64
	t := float64((90.0 - math.Abs(theta)) * d2r)
	if t < 1.0e-5 {
		if theta > 0.0 {
			z = t * t / 2.0
		} else {
			z = 2.0 - t*t/2.0
		}
		costhe = t
	} else {
		z = 1.0 - sind(theta)
		costhe = cosd(theta)
	}
	r := float64(r2d * costhe)
	if p.w1 == 0.0 {
		if theta < 0.0 {
			return 0, 0, false
		}
		return float64(r * sinphi), float64(-r * cosphi), true
	}
	z *= r2d
	if theta < -atand(float64(p.xi*sinphi)-float64(p.eta*cosphi)) {
		return 0, 0, false
	}
	return float64(r*sinphi) + float64(p.xi*z), float64(-r*cosphi) + float64(p.eta*z), true
}

func (p sinProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	x0 := float64(x * d2r)
	y0 := float64(y * d2r)
	r2 := float64(x0*x0) + float64(y0*y0)
	if p.w1 == 0.0 {
		if r2 != 0.0 {
			phi = atan2d(x0, -y0)
		}
		switch {
		case r2 < 0.5:
			theta = acosd(math.Sqrt(r2))
		case r2 <= 1.0:
			theta = asind(math.Sqrt(1.0 - r2))
		default:
			return 0, 0, false
		}
		return phi, theta, prjbchk(1.0e-13, &phi, &theta)
	}
	xy := float64(x0*p.xi) + float64(y0*p.eta)
	var z float64
	if r2 < 1.0e-10 {
		z = r2 / 2.0
		theta = 90.0 - float64(r2d*math.Sqrt(r2/(1.0+xy)))
	} else {
		a := p.w1 + 1.0
		b := xy - p.w1
		c := r2 - xy - xy + (p.w1 - 1.0)
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
	x1 := -y0 + float64(p.eta*z)
	y1 := x0 - float64(p.xi*z)
	if x1 != 0.0 || y1 != 0.0 {
		phi = atan2d(y1, x1)
	}
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// ARC — zenithal equidistant: R_theta = 90 - theta (degrees).
// The whole sphere; R > 180 is off the projection.
type arcProjection struct{}

func (arcProjection) Code() string    { return "ARC" }
func (arcProjection) theta0() float64 { return 90 }

func (arcProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	r := 90.0 - theta
	return float64(r * sinphi), float64(-r * cosphi), true
}

func (arcProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	r := math.Sqrt(float64(x*x) + float64(y*y))
	if r == 0.0 {
		theta = 90.0
	} else {
		phi = atan2d(x, -y)
		theta = 90.0 - r
	}
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// ZEA — zenithal equal-area: R_theta = 2 r0 sin((90 - theta)/2).
// The whole sphere; R > 2 r0 is off the projection.
type zeaProjection struct{}

func (zeaProjection) Code() string    { return "ZEA" }
func (zeaProjection) theta0() float64 { return 90 }

func (zeaProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	r := float64((360.0 / wpi) * sind((90.0-theta)/2.0))
	return float64(r * sinphi), float64(-r * cosphi), true
}

func (zeaProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	r := math.Sqrt(float64(x*x) + float64(y*y))
	if r != 0.0 {
		phi = atan2d(x, -y)
	}
	s := float64(r * (wpi / 360.0))
	if math.Abs(s) > 1.0 {
		if math.Abs(r-360.0/wpi) < tol {
			theta = -90.0
		} else {
			return 0, 0, false
		}
	} else {
		theta = 90.0 - float64(2.0*asind(s))
	}
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}
