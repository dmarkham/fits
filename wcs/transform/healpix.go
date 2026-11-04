package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// HEALPix projections: HPX (Calabretta & Roukema 2007, MNRAS 381, 865)
// and XPH, the polar "butterfly" HEALPix (Calabretta & Lowe 2013, PASA 30,
// e050). Both are equal-area by construction. Ports of the wcslib 8.10 prj.c
// functions of the same names, in degrees (see projection.go).

// HPX — HEALPix (Calabretta & Roukema 2007).
//
//	PV2_1 = H  number of facets in longitude (default 4)
//	PV2_2 = K  number of facets in latitude  (default 3)
//
// Equatorial zone, |sin(theta)| <= (K-1)/K:
//
//	x = phi,  y = (90 K / H) sin(theta)
//
// Polar caps, with sigma = sqrt(K (1 - |sin theta|)) and phi_c the
// longitude of the centre of the facet containing phi:
//
//	x = phi_c + (phi - phi_c) sigma
//	y = +-(180/H) ((K + 1)/2 - sigma)
//
// For odd K the polar facets of the two caps line up; for even K the
// southern cap is shifted by half a facet (wcslib's "offset"). Points in
// the gaps between the polar triangles are rejected.
type hpxProjection struct {
	K                                  float64
	m, n                               int
	w2, w3, w4, w5, w6, w7, w8, w9, w0 float64
}

func newHPX(pv map[wcs.PVKey]float64, latAxis int) (hpxProjection, error) {
	H := pvFloat(pv, latAxis, 1, 4)
	K := pvFloat(pv, latAxis, 2, 3)
	if H <= 0.0 || K <= 0.0 {
		return hpxProjection{}, errBadParam("HPX", "H = PV%d_1 = %g and K = PV%d_2 = %g must be > 0", latAxis, H, latAxis, K)
	}
	p := hpxProjection{K: K}
	p.m = int(H+0.5) % 2
	p.n = int(K+0.5) % 2
	p.w0 = 1.0
	p.w2 = (K - 1.0) / K
	p.w3 = 90.0 * K / H
	p.w4 = (K + 1.0) / 2.0
	p.w5 = 90.0 * (K - 1.0) / H
	p.w6 = 180.0 / H
	p.w7 = H / 360.0
	p.w8 = float64(p.w3 * p.w0)
	p.w9 = float64(p.w6 * p.w0)
	return p, nil
}

func (hpxProjection) Code() string    { return "HPX" }
func (hpxProjection) theta0() float64 { return 0 }

func (p hpxProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	s := x
	t := -180.0 + float64((float64(2.0*math.Floor(float64((x+180.0)*p.w7)))+1.0)*p.w6)
	t = x - t
	slim := p.w6 + 1e-12
	ylim := float64(p.w9 * p.w4)
	absy := math.Abs(y)
	switch {
	case absy <= p.w5:
		phi = s
		theta = asind(y / p.w3)
	case absy <= ylim:
		offset := !(p.n != 0 || y > 0.0)
		sigma := p.w4 - absy/p.w6
		var s2, tt float64
		if sigma == 0.0 {
			s2 = 1e9
			tt = 90.0
		} else {
			tt = 1.0 - sigma*sigma/p.K
			if tt < -1.0 {
				return 0, 0, false
			}
			s2 = 1.0 / sigma
			tt = asind(tt)
		}
		if y < 0.0 {
			tt = -tt
		}
		if offset {
			h := int(math.Floor(s/p.w6)) + p.m
			if h%2 != 0 {
				t -= p.w6
			} else {
				t += p.w6
			}
		}
		r := float64(s2 * t)
		if slim <= math.Abs(r) {
			return 0, 0, false
		}
		if r != 0.0 {
			r -= t
		}
		phi = s + r
		theta = tt
	default:
		return 0, 0, false
	}
	return phi, theta, prjbchk(1.0e-12, &phi, &theta)
}

func (p hpxProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	xi := float64(p.w0 * phi)
	t := -180.0 + float64((float64(2.0*math.Floor(float64((phi+180.0)*p.w7)))+1.0)*p.w6)
	t = float64(p.w0 * (phi - t))
	sinthe := sind(theta)
	abssin := math.Abs(sinthe)
	if abssin <= p.w2 {
		return xi, float64(p.w8 * sinthe), true
	}
	offset := !(p.n != 0 || theta > 0.0)
	sigma := math.Sqrt(float64(p.K * (1.0 - abssin)))
	xi2 := sigma - 1.0
	eta := float64(p.w9 * (p.w4 - sigma))
	if theta < 0 {
		eta = -eta
	}
	if offset {
		h := int(math.Floor(xi/p.w9)) + p.m
		if h%2 != 0 {
			t -= p.w9
		} else {
			t += p.w9
		}
	}
	x = xi + float64(t*xi2)
	if 180.0 < x {
		x = 360.0 - x
	}
	return x, eta, true
}

// XPH — HEALPix polar, "butterfly" (Calabretta & Lowe 2013). The sphere is
// cut along the four meridians phi = -180, -90, 0, 90 and laid out as four
// HPX quadrants rotated by 45 deg about the native pole, which is the
// fiducial point (theta_0 = 90). Plane points outside the butterfly are
// rejected.
type xphProjection struct{}

var (
	xphW0 = 1.0 / math.Sqrt(2.0)
	xphW2 = 2.0 / 3.0
	xphW3 = 1e-4
	xphW4 = float64(math.Sqrt(xphW2) * r2d)
	xphW5 = 90.0 - float64(xphW3*xphW4)
	xphW6 = float64(math.Sqrt(1.5) * d2r)
)

func (xphProjection) Code() string    { return "XPH" }
func (xphProjection) theta0() float64 { return 90 }

func (xphProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	xr := float64(x * xphW0)
	yr := float64(y * xphW0)
	var xi1, eta1 float64
	switch {
	case xr <= 0.0 && 0.0 < yr:
		xi1 = -xr - yr
		eta1 = xr - yr
		phi = -180.0
	case xr < 0.0 && yr <= 0.0:
		xi1 = xr - yr
		eta1 = xr + yr
		phi = -90.0
	case 0.0 <= xr && yr < 0.0:
		xi1 = xr + yr
		eta1 = -xr + yr
		phi = 0.0
	default:
		xi1 = -xr + yr
		eta1 = -xr - yr
		phi = 90.0
	}
	xi := xi1 + 45.0
	eta := eta1 + 90.0
	abseta := math.Abs(eta)
	if abseta > 90.0 {
		return 0, 0, false
	}
	if abseta <= 45.0 {
		phi += xi
		theta = asind(eta / 67.5)
		if 45.0+tol < math.Abs(xi1) {
			return 0, 0, false
		}
	} else {
		sigma := (90.0 - abseta) / 45.0
		switch {
		case xr == 0.0:
			if yr <= 0.0 {
				phi = 0.0
			} else {
				phi = 180.0
			}
		case yr == 0.0:
			if xr < 0.0 {
				phi = -90.0
			} else {
				phi = 90.0
			}
		default:
			phi += 45.0 + xi1/sigma
		}
		if sigma < xphW3 {
			theta = 90.0 - float64(sigma*xphW4)
		} else {
			theta = asind(1.0 - sigma*sigma/3.0)
		}
		if eta < 0.0 {
			theta = -theta
		}
		if eta < -45.0 && eta+90.0+tol < math.Abs(xi1) {
			return 0, 0, false
		}
	}
	return phi, theta, prjbchk(1.0e-12, &phi, &theta)
}

func (xphProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	chi := phi
	if 180.0 <= math.Abs(chi) {
		chi = math.Mod(chi, 360.0)
		if chi < -180.0 {
			chi += 360.0
		} else if 180.0 <= chi {
			chi -= 360.0
		}
	}
	chi += 180.0
	psi := math.Mod(chi, 90.0)
	yq := chi - 180.0 // selects the quadrant

	sinthe := sind(theta)
	abssin := math.Abs(sinthe)
	var xi, eta float64
	if abssin <= xphW2 {
		xi = psi
		eta = float64(67.5 * sinthe)
	} else {
		var sigma float64
		if theta < xphW5 {
			sigma = math.Sqrt(float64(3.0 * (1.0 - abssin)))
		} else {
			sigma = float64((90.0 - theta) * xphW6)
		}
		xi = 45.0 + float64((psi-45.0)*sigma)
		eta = float64(45.0 * (2.0 - sigma))
		if theta < 0.0 {
			eta = -eta
		}
	}
	xi -= 45.0
	eta -= 90.0
	switch {
	case yq < -90.0:
		return float64(xphW0 * (-xi + eta)), float64(xphW0 * (-xi - eta)), true
	case yq < 0.0:
		return float64(xphW0 * (+xi + eta)), float64(xphW0 * (-xi + eta)), true
	case yq < 90.0:
		return float64(xphW0 * (xi - eta)), float64(xphW0 * (xi + eta)), true
	}
	return float64(xphW0 * (-xi - eta)), float64(xphW0 * (xi - eta)), true
}
