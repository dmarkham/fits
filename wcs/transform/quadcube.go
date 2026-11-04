package transform

import "math"

// Quad-cube (spherical cube) projections from Paper II §5.6:
//
//   TSC — tangential spherical cube
//   CSC — COBE quadrilateralized cube
//   QSC — quadrilateralized cube
//
// Ported from wcslib prj.c in wcslib's degree units, operation for
// operation (CSC also in its single precision; see projection.go). Until
// round 2 of the proof audit TSC and QSC were radian-native ports, and QSC
// was up to 4e-11 deg off wcslib next to the face centres. All three agree
// with wcslib 8.10, with identical accept/reject decisions, on the wref
// grids that run past every face boundary (TestWrefProjections,
// testdata/wref/), and to 4 ULP next to the face centres
// (TestQuadcubeFaceCentreBand) except where QSC's own cancellation next to
// the centres of faces 2, 3, 4 and the pole face magnifies the last bit of
// Go's atan2 against glibc's (there within the point's libm sensitivity,
// libm_sens_test.go). With the correctly rounded libm the port used
// during the proof audit, a dense scan (per code 3.8e5 forward points,
// uniform and within 1e-6..0.1 deg of every face centre, and 2e5 inverse)
// agreed to 5.7e-14 deg forward (one ULP at 270 deg) and 2.8e-14 deg
// inverse. The forward
// face-coordinate overshoot rejections of tscs2x, qscs2x and cscs2x are
// ported too; no (phi, theta) reaches them in wcslib or here (see
// quadcubeFaceBound and the CSC notes), so TestQuadcubeFaceBound pins the
// TSC/QSC threshold directly.
//
// # Face layout (wcslib convention)
//
// The six cube faces are unfolded into a cross pattern in the output
// plane. In "face-cell" units where each face is a 2×2 square:
//
//	              +------+
//	              | f=0  |  top    (+z direction)
//	+------+------+------+------+
//	| f=4  | f=1  | f=2  | f=3  |  equatorial band (left→right)
//	+------+------+------+------+
//	              | f=5  |  bottom (-z direction)
//	              +------+
//
// With face-cell offsets:
//
//	f=0: (0,+2)   f=1: (0,0)    f=2: (+2,0)
//	f=3: (+4,0)   f=4: (+6,0)   f=5: (0,-2)
//
// Face 4 (left, -y direction) sits at x = +6 in the cross; points in
// the image plane with xf < -1 are mapped there via xf += 8.
//

// Each face is 90 x 90 deg in the plane: wcslib's w[0] = 45 deg per
// half-face (r0 = 180/pi), face coordinates xf, yf in [-1, 1].

// quadcubeFaceXY is the face selection of wcslib 8.10 tscx2s / qscx2s: the
// plane point (x, y) in degrees to the face and its face coordinates, or
// ok = false off the unfolded cube.
func quadcubeFaceXY(x, y float64) (face int, xf, yf float64, ok bool) {
	xf = float64(x * (1.0 / 45.0))
	yf = float64(y * (1.0 / 45.0))
	if math.Abs(xf) <= 1.0 {
		if math.Abs(yf) > 3.0 {
			return 0, 0, 0, false
		}
	} else if math.Abs(xf) > 7.0 || math.Abs(yf) > 1.0 {
		return 0, 0, 0, false
	}
	if xf < -1.0 {
		xf += 8.0
	}
	switch {
	case xf > 5.0:
		face, xf = 4, xf-6.0
	case xf > 3.0:
		face, xf = 3, xf-4.0
	case xf > 1.0:
		face, xf = 2, xf-2.0
	case yf > 1.0:
		face, yf = 0, yf-2.0
	case yf < -1.0:
		face, yf = 5, yf+2.0
	default:
		face = 1
	}
	return face, xf, yf, true
}

// quadcubeFace returns the face index (0..5) for a point on the unit
// sphere with direction cosines (l, m, n) = (cosθ·cosφ, cosθ·sinφ,
// sinθ). The face is whichever signed direction cosine is the
// maximum — ties broken in a deterministic order (0→1→2→3→4→5)
// matching wcslib's sequential if-else in tscs2x/cscs2x/qscs2x.
func quadcubeFace(l, m, n float64) (face int, zeta float64) {
	face = 0
	zeta = n
	if l > zeta {
		face = 1
		zeta = l
	}
	if m > zeta {
		face = 2
		zeta = m
	}
	if -l > zeta {
		face = 3
		zeta = -l
	}
	if -m > zeta {
		face = 4
		zeta = -m
	}
	if -n > zeta {
		face = 5
		zeta = -n
	}
	return face, zeta
}

// --------------------------------------------------------------------
// TSC — tangential spherical cube (Paper II, wcslib prj.c)
// --------------------------------------------------------------------

type tscProjection struct{}

func (tscProjection) Code() string    { return "TSC" }
func (tscProjection) theta0() float64 { return 0 }

// s2x is wcslib tscs2x.
func (tscProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	sinthe, costhe := sincosd(theta)
	l := float64(costhe * cosphi)
	m := float64(costhe * sinphi)
	n := sinthe

	face, zeta := quadcubeFace(l, m, n)
	var xf, yf, x0, y0 float64
	switch face {
	case 1:
		xf, yf, x0, y0 = m/zeta, n/zeta, 0, 0
	case 2:
		xf, yf, x0, y0 = -l/zeta, n/zeta, 2, 0
	case 3:
		xf, yf, x0, y0 = -m/zeta, n/zeta, 4, 0
	case 4:
		xf, yf, x0, y0 = l/zeta, n/zeta, 6, 0
	case 5:
		xf, yf, x0, y0 = m/zeta, l/zeta, 0, -2
	default:
		xf, yf, x0, y0 = m/zeta, -l/zeta, 0, 2
	}
	xf, yf, ok = quadcubeFaceBound(xf, yf)
	if !ok {
		return 0, 0, false
	}
	return float64(45.0 * (xf + x0)), float64(45.0 * (yf + y0)), true
}

// quadcubeFaceBound is the face-coordinate check of wcslib 8.10 tscs2x and
// qscs2x: a face coordinate past +-1 by at most quadcubeFaceTol is clamped
// to +-1, one further out rejects the point (PRJERR_BAD_WORLD_SET).
//
// No (phi, theta) reaches the rejection, in wcslib or here: for TSC the
// face is the one whose direction cosine zeta is largest in magnitude, so
// |xf| = |m| / zeta <= 1 exactly in IEEE arithmetic (division is
// monotonic); QSC's area-equalised xf, yf can exceed 1 only by rounding
// (a few ULP, which the clamp absorbs). A wcslib scan of 4.62e8 forward
// points per code (1e-6 deg spacing around all 16 face corners and along
// the face edges, the verifier's qscan.c) found no rejection for TSC, QSC
// or CSC. The check is ported anyway so that the decision rule is
// wcslib's, and TestQuadcubeFaceBound pins its threshold directly.
func quadcubeFaceBound(xf, yf float64) (float64, float64, bool) {
	ok := true
	if math.Abs(xf) > 1.0 {
		if math.Abs(xf) > 1.0+quadcubeFaceTol {
			ok = false
		}
		xf = math.Copysign(1.0, xf)
	}
	if math.Abs(yf) > 1.0 {
		if math.Abs(yf) > 1.0+quadcubeFaceTol {
			ok = false
		}
		yf = math.Copysign(1.0, yf)
	}
	return xf, yf, ok
}

// quadcubeFaceTol is the tol of wcslib 8.10 tscs2x / qscs2x (1.0e-12).
const quadcubeFaceTol = 1.0e-12

// x2s is wcslib tscx2s.
func (tscProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	face, xf, yf, ok := quadcubeFaceXY(x, y)
	if !ok {
		return 0, 0, false
	}
	r := 1.0 / math.Sqrt(1.0+float64(xf*xf)+float64(yf*yf))
	var l, m, n float64
	switch face {
	case 4:
		m = -r
		l = float64(-m * xf)
		n = float64(-m * yf)
	case 3:
		l = -r
		m = float64(l * xf)
		n = float64(-l * yf)
	case 2:
		m = r
		l = float64(-m * xf)
		n = float64(m * yf)
	case 0:
		n = r
		l = float64(-n * yf)
		m = float64(n * xf)
	case 5:
		n = -r
		l = float64(-n * yf)
		m = float64(-n * xf)
	default:
		l = r
		m = float64(l * xf)
		n = float64(l * yf)
	}
	if l != 0.0 || m != 0.0 {
		phi = atan2d(m, l)
	}
	theta = asind(n)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// --------------------------------------------------------------------
// CSC — COBE quadrilateralized spherical cube (Paper II)
// --------------------------------------------------------------------
//
// CSC uses two INDEPENDENT polynomial fits (Chan & O'Neill 1975; O'Neill
// & Laubscher 1976), one for each direction; they are not exact inverses
// of each other: a forward then inverse round trip through wcslib's own
// CSC drifts by up to 0.012 deg on the wref grid, and ours drifts
// identically. The port is in wcslib's degree units and single-precision
// arithmetic (see projection.go). wcslib's forward rejects a point whose
// face coordinate overshoots +-1 by more than 1e-7, and the port keeps that
// check, but no point reaches it: wcslib's own cscs2x rejects none of the
// 6.5e8 points of a 0.01 deg (phi, theta) grid, so the check is untested
// (a mutant that removes it survives every test).

// cscChiPsi evaluates the CSC inverse polynomial P(xf, yf) → chi, psi.
// Ported from wcslib cscx2s. The polynomial is a
// degree-7 bivariate in xx = xf² and yy = yf² with 28 coefficients,
// evaluated as a sum of row polynomials z0..z6 then composed.
//
// wcslib evaluates these polynomials in single precision (float) and so
// do we, with the same float/double conversion points as cscx2s/cscs2x.
func cscChiPsi(xf, yf float32) (chi, psi float32) {
	const (
		p00 = float32(-0.27292696)
		p10 = float32(-0.07629969)
		p20 = float32(-0.22797056)
		p30 = float32(0.54852384)
		p40 = float32(-0.62930065)
		p50 = float32(0.25795794)
		p60 = float32(0.02584375)
		p01 = float32(-0.02819452)
		p11 = float32(-0.01471565)
		p21 = float32(0.48051509)
		p31 = float32(-1.74114454)
		p41 = float32(1.71547508)
		p51 = float32(-0.53022337)
		p02 = float32(0.27058160)
		p12 = float32(-0.56800938)
		p22 = float32(0.30803317)
		p32 = float32(0.98938102)
		p42 = float32(-0.83180469)
		p03 = float32(-0.60441560)
		p13 = float32(1.50880086)
		p23 = float32(-0.93678576)
		p33 = float32(0.08693841)
		p04 = float32(0.93412077)
		p14 = float32(-1.41601920)
		p24 = float32(0.33887446)
		p05 = float32(-0.63915306)
		p15 = float32(0.52032238)
		p06 = float32(0.14381585)
	)
	xx := float32(xf * xf)
	yy := float32(yf * yf)

	// chi polynomial: bivariate in (xx, yy).
	z0 := p00 + float32(xx*(p10+float32(xx*(p20+float32(xx*(p30+float32(xx*(p40+float32(xx*(p50+float32(xx*p60)))))))))))
	z1 := p01 + float32(xx*(p11+float32(xx*(p21+float32(xx*(p31+float32(xx*(p41+float32(xx*p51)))))))))
	z2 := p02 + float32(xx*(p12+float32(xx*(p22+float32(xx*(p32+float32(xx*p42)))))))
	z3 := p03 + float32(xx*(p13+float32(xx*(p23+float32(xx*p33)))))
	z4 := p04 + float32(xx*(p14+float32(xx*p24)))
	z5 := p05 + float32(xx*p15)
	z6 := p06
	chi32 := z0 + float32(yy*(z1+float32(yy*(z2+float32(yy*(z3+float32(yy*(z4+float32(yy*(z5+float32(yy*z6)))))))))))
	chi32 = xf + float32(xf*(1-xx)*chi32)

	// psi polynomial: same structure with xx and yy swapped.
	z0 = p00 + float32(yy*(p10+float32(yy*(p20+float32(yy*(p30+float32(yy*(p40+float32(yy*(p50+float32(yy*p60)))))))))))
	z1 = p01 + float32(yy*(p11+float32(yy*(p21+float32(yy*(p31+float32(yy*(p41+float32(yy*p51)))))))))
	z2 = p02 + float32(yy*(p12+float32(yy*(p22+float32(yy*(p32+float32(yy*p42)))))))
	z3 = p03 + float32(yy*(p13+float32(yy*(p23+float32(yy*p33)))))
	z4 = p04 + float32(yy*(p14+float32(yy*p24)))
	z5 = p05 + float32(yy*p15)
	z6 = p06
	psi32 := z0 + float32(xx*(z1+float32(xx*(z2+float32(xx*(z3+float32(xx*(z4+float32(xx*(z5+float32(xx*z6)))))))))))
	psi32 = yf + float32(yf*(1-yy)*psi32)

	return chi32, psi32
}

// cscXfYf's coefficients (wcslib cscs2x, single precision).
const (
	cscGstar  = float32(1.37484847732)
	cscMm     = float32(0.004869491981)
	cscGamma  = float32(-0.13161671474)
	cscOmega1 = float32(-0.159596235474)
	cscD0     = float32(0.0759196200467)
	cscD1     = float32(-0.0217762490699)
	cscC00    = float32(0.141189631152)
	cscC10    = float32(0.0809701286525)
	cscC01    = float32(-0.281528535557)
	cscC11    = float32(0.15384112876)
	cscC20    = float32(-0.178251207466)
	cscC02    = float32(0.106959469314)
)

// cscXfYf evaluates the CSC forward polynomial (chi, psi) → (xf, yf).
// Different coefficient set from the inverse — this is a separate
// least-squares fit. Ported from wcslib cscs2x.
func cscXfYf(chi, psi float32) (xf, yf float32) {
	chi2 := float32(chi * chi)
	psi2 := float32(psi * psi)
	chi2co := 1 - chi2
	psi2co := 1 - psi2

	// Avoid FP underflow (wcslib's guards, kept as written). They decide
	// nothing: dropping chi4 (chi2 <= T), psi4 (psi2 <= T) or chi2psi2
	// (|chi psi| <= T) for any T <= 1e-8 never changes xf or yf, so the
	// 1e-16 here could be 0 or 1e-8 (TestCSCUnderflowGuardPremises checks
	// the bounds this uses). Each dropped term t is at most
	// 0.18 T^2 <= 1.8e-17 in magnitude and enters the inner sum
	// P = cscC00 + cscC10 a + cscC01 b + cscC11 ab + cscC20 a^2 + cscC02 b^2 (a, b = chi2,
	// psi2 for xf; psi2, chi2 for yf), evaluated left to right in float32.
	// fl(S + t) = S whenever |S| > 2^25 |t| (6e-10), S the partial sum
	// before t. For the chi4/psi4 term in 6th position (the guard on the
	// variable entering as b = a tiny value) S = cscC00 + cscC10 a + cscC20 a^2 +
	// tiny >= 0.0439 on a in [0, 1]. For a term in 4th or 5th position
	// with S below 6e-10, cscC01 b must cancel cscC00 (+ cscC10 a), so b is in
	// [0.50, 0.80] and the 6th term cscC02 b^2 >= 0.026: S + t and S are both
	// below half an ULP of it (a float32 >= 2^-6 has half an ULP
	// >= 9.3e-10) and the sum rounds to the same float32 either way. So P, and xf and
	// yf, are identical with and without the guards.
	chipsi := float32(math.Abs(float64(chi * psi)))
	var chi4, psi4, chi2psi2 float32
	if chi2 > 1e-16 {
		chi4 = float32(chi2 * chi2)
	}
	if psi2 > 1e-16 {
		psi4 = float32(psi2 * psi2)
	}
	if chipsi > 1e-16 {
		chi2psi2 = float32(chi2 * psi2)
	}

	xf32 := float32(chi * (chi2 + float32(chi2co*(cscGstar+float32(psi2*(float32(cscGamma*chi2co)+float32(cscMm*chi2)+
		float32(psi2co*(cscC00+float32(cscC10*chi2)+float32(cscC01*psi2)+float32(cscC11*chi2psi2)+float32(cscC20*chi4)+float32(cscC02*psi4)))))+
		float32(chi2*(cscOmega1-float32(chi2co*(cscD0+float32(cscD1*chi2)))))))))
	yf32 := float32(psi * (psi2 + float32(psi2co*(cscGstar+float32(chi2*(float32(cscGamma*psi2co)+float32(cscMm*psi2)+
		float32(chi2co*(cscC00+float32(cscC10*psi2)+float32(cscC01*chi2)+float32(cscC11*chi2psi2)+float32(cscC20*psi4)+float32(cscC02*chi4)))))+
		float32(psi2*(cscOmega1-float32(psi2co*(cscD0+float32(cscD1*psi2)))))))))
	return xf32, yf32
}

type cscProjection struct{}

func (cscProjection) Code() string    { return "CSC" }
func (cscProjection) theta0() float64 { return 0 }

// s2x is wcslib cscs2x.
func (cscProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	const tol = 1.0e-7
	sinphi, cosphi := sincosd(phi)
	sinthe, costhe := sincosd(theta)
	l := float64(costhe * cosphi)
	m := float64(costhe * sinphi)
	n := sinthe

	face := 0
	zeta := n
	if l > zeta {
		face, zeta = 1, l
	}
	if m > zeta {
		face, zeta = 2, m
	}
	if -l > zeta {
		face, zeta = 3, -l
	}
	if -m > zeta {
		face, zeta = 4, -m
	}
	if -n > zeta {
		face, zeta = 5, -n
	}

	var xi, eta float64
	var x0, y0 float32
	switch face {
	case 1:
		xi, eta, x0, y0 = m, n, 0, 0
	case 2:
		xi, eta, x0, y0 = -l, n, 2, 0
	case 3:
		xi, eta, x0, y0 = -m, n, 4, 0
	case 4:
		xi, eta, x0, y0 = l, n, 6, 0
	case 5:
		xi, eta, x0, y0 = m, l, 0, -2
	default:
		xi, eta, x0, y0 = m, -l, 0, 2
	}

	xf, yf := cscXfYf(float32(xi/zeta), float32(eta/zeta))
	if math.Abs(float64(xf)) > 1.0 {
		if math.Abs(float64(xf)) > 1.0+tol {
			return 0, 0, false
		}
		xf = float32(math.Copysign(1.0, float64(xf)))
	}
	if math.Abs(float64(yf)) > 1.0 {
		if math.Abs(float64(yf)) > 1.0+tol {
			return 0, 0, false
		}
		yf = float32(math.Copysign(1.0, float64(yf)))
	}
	return float64(45.0 * float64(xf+x0)), float64(45.0 * float64(yf+y0)), true
}

// x2s is wcslib cscx2s.
func (cscProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	xf := float32(x * (1.0 / 45.0))
	yf := float32(y * (1.0 / 45.0))
	if math.Abs(float64(xf)) <= 1.0 {
		if math.Abs(float64(yf)) > 3.0 {
			return 0, 0, false
		}
	} else if math.Abs(float64(xf)) > 7.0 || math.Abs(float64(yf)) > 1.0 {
		return 0, 0, false
	}
	if xf < -1.0 {
		xf += 8.0
	}
	face := 0
	switch {
	case xf > 5.0:
		face, xf = 4, xf-6.0
	case xf > 3.0:
		face, xf = 3, xf-4.0
	case xf > 1.0:
		face, xf = 2, xf-2.0
	case yf > 1.0:
		face, yf = 0, yf-2.0
	case yf < -1.0:
		face, yf = 5, yf+2.0
	default:
		face = 1
	}

	chi32, psi32 := cscChiPsi(xf, yf)
	chi, psi := float64(chi32), float64(psi32)
	t := 1.0 / math.Sqrt(float64(float32(chi32*chi32)+float32(psi32*psi32))+1.0)
	var l, m, n float64
	switch face {
	case 1:
		l = t
		m = float64(chi * l)
		n = float64(psi * l)
	case 2:
		m = t
		l = float64(-chi * m)
		n = float64(psi * m)
	case 3:
		l = -t
		m = float64(chi * l)
		n = float64(-psi * l)
	case 4:
		m = -t
		l = float64(-chi * m)
		n = float64(-psi * m)
	case 5:
		n = -t
		l = float64(-psi * n)
		m = float64(-chi * n)
	default:
		n = t
		l = float64(-psi * n)
		m = float64(chi * n)
	}
	if l != 0.0 || m != 0.0 {
		phi = atan2d(m, l)
	}
	theta = asind(n)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}

// --------------------------------------------------------------------
// QSC — quadrilateralized spherical cube (Paper II, wcslib prj.c)
// --------------------------------------------------------------------
//
// QSC is equal-area by construction (unlike CSC which is an
// approximation). The forward and inverse are closed-form but each
// has four sub-quadrant branches within every face, giving 24 total
// code paths. Ported from wcslib qscs2x and qscx2s in wcslib's degree
// units, operation for operation (the earlier radian port differed from
// wcslib by up to 4e-11 deg next to the face centres, where the
// small-angle formula's atan2 and 90 - theta are computed).

// qscSqrt2Inv is wcsmath.h's SQRT2INV = 1.0/SQRT2, evaluated in double
// arithmetic on the rounded SQRT2 as C does: 0x1.6a09e667f3bccp-1, one ULP
// below the correctly rounded 1/sqrt(2) (0x1.6a09e667f3bcdp-1).
var qscSqrt2Inv = 1.0 / wsqrt2

// wsqrt2 is wcsmath.h's SQRT2 (a variable so the division above is done
// in double arithmetic, see wpi).
var wsqrt2 = 1.4142135623730950488

type qscProjection struct{}

func (qscProjection) Code() string    { return "QSC" }
func (qscProjection) theta0() float64 { return 0 }

// qscZecoSmall is the zeco below which qscs2x replaces 1 - zeta by its
// small-angle form (wcslib 8.10: if (zeco < 1.0e-8)). The comparison can
// never be tied: 1 - zeta is exact (zeta is in [1/2, 1] there, Sterbenz)
// and a multiple of 2^-53, and 1e-8 is not (TestQSCZecoSwitchTie).
const qscZecoSmall = 1.0e-8

// s2x is wcslib qscs2x.
func (qscProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi)
	sinthe, costhe := sincosd(theta)
	if math.Abs(theta) == 90.0 {
		// wcslib: x = -prj->x0 (= -0.0), y = copysign(2 w[0], theta) - y0.
		return math.Copysign(0, -1), math.Copysign(2.0*45.0, theta), true
	}
	l := float64(costhe * cosphi)
	m := float64(costhe * sinphi)
	n := sinthe

	face, zeta := quadcubeFace(l, m, n)
	zeco := 1.0 - zeta
	var xi, eta, x0, y0 float64
	switch face {
	case 1:
		xi, eta = m, n
		if zeco < qscZecoSmall {
			t := float64(theta * d2r)
			p := lmAtan2(sinphi, cosphi)
			zeco = (float64(p*p) + float64(t*t)) / 2.0
		}
		x0, y0 = 0.0, 0.0
	case 2:
		xi, eta = -l, n
		if zeco < qscZecoSmall {
			t := float64(theta * d2r)
			p := lmAtan2(sinphi, cosphi) - wpi/2.0
			zeco = (float64(p*p) + float64(t*t)) / 2.0
		}
		x0, y0 = 2.0, 0.0
	case 3:
		xi, eta = -m, n
		if zeco < qscZecoSmall {
			t := float64(theta * d2r)
			p := lmAtan2(sinphi, cosphi)
			p -= math.Copysign(wpi, p)
			zeco = (float64(p*p) + float64(t*t)) / 2.0
		}
		x0, y0 = 4.0, 0.0
	case 4:
		xi, eta = l, n
		if zeco < qscZecoSmall {
			t := float64(theta * d2r)
			p := lmAtan2(sinphi, cosphi) + wpi/2.0
			zeco = (float64(p*p) + float64(t*t)) / 2.0
		}
		x0, y0 = 6.0, 0.0
	case 5:
		xi, eta = m, l
		if zeco < qscZecoSmall {
			t := float64((theta + 90.0) * d2r)
			zeco = float64(t*t) / 2.0
		}
		x0, y0 = 0.0, -2.0
	default:
		xi, eta = m, -l
		if zeco < qscZecoSmall {
			t := float64((90.0 - theta) * d2r)
			zeco = float64(t*t) / 2.0
		}
		x0, y0 = 0.0, 2.0
	}

	var xf, yf float64
	if xi != 0.0 || eta != 0.0 {
		switch {
		case -xi > math.Abs(eta):
			omega := eta / xi
			tau := 1.0 + float64(omega*omega)
			xf = -math.Sqrt(zeco / (1.0 - 1.0/math.Sqrt(1.0+tau)))
			yf = float64((xf / 15.0) * (atand(omega) - asind(omega/math.Sqrt(tau+tau))))
		case xi > math.Abs(eta):
			omega := eta / xi
			tau := 1.0 + float64(omega*omega)
			xf = math.Sqrt(zeco / (1.0 - 1.0/math.Sqrt(1.0+tau)))
			yf = float64((xf / 15.0) * (atand(omega) - asind(omega/math.Sqrt(tau+tau))))
		case -eta >= math.Abs(xi):
			omega := xi / eta
			tau := 1.0 + float64(omega*omega)
			yf = -math.Sqrt(zeco / (1.0 - 1.0/math.Sqrt(1.0+tau)))
			xf = float64((yf / 15.0) * (atand(omega) - asind(omega/math.Sqrt(tau+tau))))
		case eta >= math.Abs(xi):
			omega := xi / eta
			tau := 1.0 + float64(omega*omega)
			yf = math.Sqrt(zeco / (1.0 - 1.0/math.Sqrt(1.0+tau)))
			xf = float64((yf / 15.0) * (atand(omega) - asind(omega/math.Sqrt(tau+tau))))
		}
	}

	// The face overshoot check (unreachable, see quadcubeFaceBound).
	xf, yf, ok = quadcubeFaceBound(xf, yf)
	if !ok {
		return 0, 0, false
	}
	return float64(45.0 * (xf + x0)), float64(45.0 * (yf + y0)), true
}

// qscInverseZeta is the first half of wcslib qscx2s: from the face
// coordinates (|xf|, |yf| <= 1, as quadcubeFaceXY returns them) the
// branch, omega, tau and zeta = 1 - zeco, the face-normal direction
// cosine.
//
// zeta is always in (0, 1], so qscx2s's "zeta < -1" rejection (and its
// 1e-12 tolerance) is unreachable: w = 15 v / u with |v| <= |u| is in
// [-15, 15] deg, so cos w >= cos 15 = 0.966 > SQRT2INV and
// |omega| = |sin w| / (cos w - SQRT2INV) <= 1.0 (0.259 / 0.259 at
// w = 15); tau = 1 + omega^2 is in [1, 2]; 1 / sqrt(1 + tau) is in
// [1/sqrt(3), 1/sqrt(2)], so 1 - 1/sqrt(1 + tau) is in [0.29, 0.43]; and
// zeco = u^2 (that) with u^2 <= 1 is in [0, 0.43]. Every step is a
// correctly rounded operation on values well inside those intervals, so
// the computed zeco is in [0, 0.5] and zeta = 1 - zeco in [0.5, 1].
// TestQSCInverseZetaRange checks it at the face corners, edges and centre
// and on a dense grid.
func qscInverseZeta(xf, yf float64) (direct bool, omega, tau, zeta, zeco float64) {
	direct = math.Abs(xf) > math.Abs(yf)
	if direct {
		if xf == 0.0 {
			return direct, 0.0, 1.0, 1.0, 0.0
		}
		w := float64(15.0*yf) / xf
		omega = sind(w) / (cosd(w) - qscSqrt2Inv)
		tau = 1.0 + float64(omega*omega)
		zeco = float64(float64(xf*xf) * (1.0 - 1.0/math.Sqrt(1.0+tau)))
		return direct, omega, tau, 1.0 - zeco, zeco
	}
	if yf == 0.0 {
		return direct, 0.0, 1.0, 1.0, 0.0
	}
	w := float64(15.0*xf) / yf
	sinw, cosw := sincosd(w)
	omega = sinw / (cosw - qscSqrt2Inv)
	tau = 1.0 + float64(omega*omega)
	zeco = float64(float64(yf*yf) * (1.0 - 1.0/math.Sqrt(1.0+tau)))
	return direct, omega, tau, 1.0 - zeco, zeco
}

// x2s is wcslib qscx2s.
func (qscProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	face, xf, yf, ok := quadcubeFaceXY(x, y)
	if !ok {
		return 0, 0, false
	}

	direct, omega, tau, zeta, zeco := qscInverseZeta(xf, yf)
	var w float64
	// zeta < -1 cannot happen (qscInverseZeta, TestQSCInverseZetaRange),
	// so this check, and its tolerance, decide nothing; wcslib's check is
	// kept as written.
	if zeta < -1.0 {
		if zeta < -1.0-tol {
			return 0, 0, false
		}
		zeta = -1.0
		w = 0.0
	} else {
		w = math.Sqrt(zeco * (2.0 - zeco) / tau)
	}

	var l, m, n float64
	switch face {
	case 1:
		l = zeta
		if direct {
			m = w
			if xf < 0.0 {
				m = -m
			}
			n = float64(m * omega)
		} else {
			n = w
			if yf < 0.0 {
				n = -n
			}
			m = float64(n * omega)
		}
	case 2:
		m = zeta
		if direct {
			l = w
			if xf > 0.0 {
				l = -l
			}
			n = float64(-l * omega)
		} else {
			n = w
			if yf < 0.0 {
				n = -n
			}
			l = float64(-n * omega)
		}
	case 3:
		l = -zeta
		if direct {
			m = w
			if xf > 0.0 {
				m = -m
			}
			n = float64(-m * omega)
		} else {
			n = w
			if yf < 0.0 {
				n = -n
			}
			m = float64(-n * omega)
		}
	case 4:
		m = -zeta
		if direct {
			l = w
			if xf < 0.0 {
				l = -l
			}
			n = float64(l * omega)
		} else {
			n = w
			if yf < 0.0 {
				n = -n
			}
			l = float64(n * omega)
		}
	case 5:
		n = -zeta
		if direct {
			m = w
			if xf < 0.0 {
				m = -m
			}
			l = float64(m * omega)
		} else {
			l = w
			if yf < 0.0 {
				l = -l
			}
			m = float64(l * omega)
		}
	default:
		n = zeta
		if direct {
			m = w
			if xf < 0.0 {
				m = -m
			}
			l = float64(-m * omega)
		} else {
			l = w
			if yf > 0.0 {
				l = -l
			}
			m = float64(-l * omega)
		}
	}
	if l != 0.0 || m != 0.0 {
		phi = atan2d(m, l)
	}
	theta = asind(n)
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}
