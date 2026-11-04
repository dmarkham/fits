package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// dssPoly is the DSS plate polynomial as wcslib 8.10 applies it (dis.c
// dssset): a sequent distortion on the intermediate pixel coordinates
// (m1, m2) = PC (p - CRPIX), between the matrix and CDELT, translated
// into wcslib's TPD form with auxiliary variables. For axis j:
//
//	u = aux[j][0] + aux[j][1] m1 + aux[j][2] m2
//	v = aux[j][3] + aux[j][4] m1 + aux[j][5] m2
//	m'_j = sum_k c[j][k] T_k(u, v)
//
// where T_k is the TPD (= TPV) term layout without the radial terms, of
// degree 3, or 5 when AMD13 is non-zero. The auxiliary variables undo the
// plate rotation and scale so that (u, v) are the plate coordinates in mm
// (x, y on axis 1; y, x on axis 2), and c are the AMD coefficients
// divided by -S (axis 1) or S (axis 2), so that CDELT = -S/3600, S/3600
// turns m' into degrees.
type dssPoly struct {
	aux    [2][6]float64
	c      [2][40]float64 // TPD slots 0..23 used (tpvTerms layout)
	degree int            // 3 or 5
}

// newDSSPoly ports dssset coefficient for coefficient.
func newDSSPoly(d *wcs.DSS) *dssPoly {
	a1, a2, a3 := d.AMDX[0], d.AMDX[1], d.AMDX[2]
	b1, b2, b3 := d.AMDY[0], d.AMDY[1], d.AMDY[2]
	det := float64(a1*b1) - float64(a2*b2)
	x0 := (float64(a2*b3) - float64(a3*b1)) / det
	y0 := (float64(a3*b2) - float64(a1*b3)) / det
	s := math.Sqrt(math.Abs(det))
	p := &dssPoly{degree: 3}
	p.aux[0] = [6]float64{x0, -b1 / s, -a2 / s, y0, b2 / s, a1 / s}
	p.aux[1] = [6]float64{y0, b2 / s, a1 / s, x0, -b1 / s, -a2 / s}
	for j, amd := range [2]*[13]float64{&d.AMDX, &d.AMDY} {
		sj := s
		if j == 0 {
			sj = -s // the A coefficients are scaled by -S
		}
		c := &p.c[j]
		for m := 1; m <= 13; m++ {
			coeff := amd[m-1]
			if coeff == 0 {
				continue
			}
			coeff /= sj
			switch m {
			case 1:
				c[1] = coeff
			case 2:
				c[2] = coeff
			case 3:
				c[0] = coeff
			case 4:
				c[4] += coeff
			case 5:
				c[5] = coeff
			case 6:
				c[6] += coeff
			case 7:
				c[4] += coeff
				c[6] += coeff
			case 8:
				c[7] += coeff
			case 9:
				c[8] = coeff
			case 10:
				c[9] += coeff
			case 11:
				c[10] = coeff
			case 12:
				c[7] += coeff
				c[9] += coeff
			case 13:
				c[17] = coeff
				c[19] = float64(coeff * 2.0)
				c[21] = coeff
				p.degree = 5
			}
		}
	}
	return p
}

// uv is the auxiliary-variable step of wcslib tpd3/tpd5.
func (p *dssPoly) uv(j int, m1, m2 float64) (u, v float64) {
	a := &p.aux[j]
	u = a[0] + float64(a[1]*m1) + float64(a[2]*m2)
	v = a[3] + float64(a[4]*m1) + float64(a[5]*m2)
	return u, v
}

// eval is wcslib tpd3 / tpd5 (no radial terms), Horner form for form.
func (p *dssPoly) eval(j int, m1, m2 float64) float64 {
	u, v := p.uv(j, m1, m2)
	c := &p.c[j]
	if p.degree == 3 {
		d := c[0] + float64(u*(c[1]+float64(u*(c[4]+float64(u*c[7])))))
		d += float64(v*(c[2]+float64(v*(c[6]+float64(v*c[10]))))) +
			float64(float64(u*(c[5]+float64(v*c[9])+float64(u*c[8])))*v)
		return d
	}
	d := c[0] + float64(u*(c[1]+float64(u*(c[4]+float64(u*(c[7]+float64(u*(c[12]+float64(u*c[17])))))))))
	d += float64(v*(c[2]+float64(v*(c[6]+float64(v*(c[10]+float64(v*(c[16]+float64(v*c[22]))))))))) +
		float64(float64(u*(c[5]+float64(v*(c[9]+float64(v*(c[15]+float64(v*c[21])))))+
			float64(u*(c[8]+float64(v*(c[14]+float64(v*c[20])))+
				float64(u*(c[13]+float64(v*c[19])+
					float64(u*c[18])))))))*v)
	return d
}

// Forward maps intermediate pixel coordinates to the distorted ones.
func (p *dssPoly) Forward(m1, m2 float64) (float64, float64) {
	return p.eval(0, m1, m2), p.eval(1, m1, m2)
}

// jac returns axis j's value and its partials with respect to (m1, m2),
// from the monomial form of the same polynomial (tpvTerms layout).
func (p *dssPoly) jac(j int, m1, m2 float64) (val, d1, d2 float64) {
	u, v := p.uv(j, m1, m2)
	_, du, dv := evalTPVWithJacobian(&p.c[j], u, v)
	a := &p.aux[j]
	return p.eval(j, m1, m2), float64(du*a[1]) + float64(dv*a[4]), float64(du*a[2]) + float64(dv*a[5])
}

// Inverse inverts the plate polynomial by Newton iteration (wcslib
// disx2p iterates on the forward polynomial too; both converge to the
// same root). ok is false when the iteration does not converge. The start
// point (the distorted coordinates themselves) and the Jacobian only
// steer the iteration: the result is the root of the forward polynomial,
// which the goldens pin; the Jacobian is pinned by TestDSSJacobian.
func (p *dssPoly) Inverse(t1, t2 float64) (m1, m2 float64, ok bool) {
	return newton2(func(x, y float64) (fx, fy, j11, j12, j21, j22 float64) {
		fx, j11, j12 = p.jac(0, x, y)
		fy, j21, j22 = p.jac(1, x, y)
		return
	}, t1, t2, t1, t2, dssNewtonFloor)
}

// dssNewtonFloor: intermediate pixel coordinates are pixels; below 1 px
// convergence is tested absolutely (1e-14 px).
const dssNewtonFloor = 1
