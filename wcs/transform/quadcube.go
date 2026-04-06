package transform

import "math"

// Quad-cube (spherical cube) projections from Paper II §5.6:
//
//   TSC — tangential spherical cube      (§5.6.1)
//   CSC — COBE quadrilateralized cube    (§5.6.2)
//   QSC — quadrilateralized cube         (§5.6.3)
//
// Ported from wcslib (reference: Siril subprojects/wcslib/prj.c:6648+).
// The face-layout helpers in this file mirror wcslib's convention so
// that forward/inverse outputs agree bit-for-bit with wcslib on every
// pixel (validated by wref goldens in testdata/wref/).
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
// The face-cell axis scale is w0 = π/4 radians per half-face, so each
// face is π/2 × π/2 radians. In wcslib's default r0 = 180/π scheme
// w0 = 45°; we work purely in radians so w0 = π/4.

const quadW0 = math.Pi / 4
const quadW1 = 4.0 / math.Pi // = 1/quadW0

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

// quadcubeFaceFromXY maps an output-plane point (x, y) in radians to
// a face index and its local (xf, yf) coordinates in [-1, 1]. Returns
// ok=false for points outside the cross pattern.
func quadcubeFaceFromXY(x, y float64) (face int, xf, yf float64, ok bool) {
	xf = x * quadW1
	yf = y * quadW1

	// Bounds check (wcslib tscx2s:6746).
	if math.Abs(xf) <= 1 {
		if math.Abs(yf) > 3 {
			return 0, 0, 0, false
		}
	} else {
		if math.Abs(xf) > 7 || math.Abs(yf) > 1 {
			return 0, 0, 0, false
		}
	}

	// Map the leftmost face strip (xf < -1) to the wrapped position xf+8.
	if xf < -1 {
		xf += 8
	}

	switch {
	case xf > 5:
		face = 4
		xf -= 6
	case xf > 3:
		face = 3
		xf -= 4
	case xf > 1:
		face = 2
		xf -= 2
	case yf > 1:
		face = 0
		yf -= 2
	case yf < -1:
		face = 5
		yf += 2
	default:
		face = 1
	}
	return face, xf, yf, true
}

// quadcubeFaceOffset returns the (x0, y0) offset in face-cell units
// (i.e. multiples of 2) for the given face.
func quadcubeFaceOffset(face int) (x0, y0 float64) {
	switch face {
	case 0:
		return 0, 2
	case 1:
		return 0, 0
	case 2:
		return 2, 0
	case 3:
		return 4, 0
	case 4:
		return 6, 0
	case 5:
		return 0, -2
	}
	return 0, 0
}

// --------------------------------------------------------------------
// TSC — tangential spherical cube (Paper II §5.6.1, wcslib prj.c:6648)
// --------------------------------------------------------------------

type tscProjection struct{}

func (tscProjection) Code() string    { return "TSC" }
func (tscProjection) Theta0() float64 { return 0 }

func (tscProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	sinTheta, cosTheta := math.Sincos(theta)
	sinPhi, cosPhi := math.Sincos(phi)
	l := cosTheta * cosPhi
	m := cosTheta * sinPhi
	n := sinTheta

	face, zeta := quadcubeFace(l, m, n)
	if zeta == 0 {
		return 0, 0, false
	}

	var xf, yf float64
	switch face {
	case 1:
		xf = m / zeta
		yf = n / zeta
	case 2:
		xf = -l / zeta
		yf = n / zeta
	case 3:
		xf = -m / zeta
		yf = n / zeta
	case 4:
		xf = l / zeta
		yf = n / zeta
	case 5:
		xf = m / zeta
		yf = l / zeta
	default: // face 0
		xf = m / zeta
		yf = -l / zeta
	}

	x0, y0 := quadcubeFaceOffset(face)
	x = quadW0 * (xf + x0)
	y = quadW0 * (yf + y0)
	return x, y, true
}

func (tscProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	face, xf, yf, fOK := quadcubeFaceFromXY(x, y)
	if !fOK {
		return 0, 0, false
	}

	// Reconstruct direction cosines (l, m, n) on the unit sphere.
	// The face-normal component is 1/sqrt(1+xf²+yf²) (gnomonic
	// projection onto the face plane).
	var l, m, n float64
	switch face {
	case 1:
		l = 1 / math.Sqrt(1+xf*xf+yf*yf)
		m = l * xf
		n = l * yf
	case 2:
		m = 1 / math.Sqrt(1+xf*xf+yf*yf)
		l = -m * xf
		n = m * yf
	case 3:
		l = -1 / math.Sqrt(1+xf*xf+yf*yf)
		m = l * xf
		n = -l * yf
	case 4:
		m = -1 / math.Sqrt(1+xf*xf+yf*yf)
		l = -m * xf
		n = -m * yf
	case 5:
		n = -1 / math.Sqrt(1+xf*xf+yf*yf)
		l = -n * yf
		m = -n * xf
	default: // face 0
		n = 1 / math.Sqrt(1+xf*xf+yf*yf)
		l = -n * yf
		m = n * xf
	}

	if l == 0 && m == 0 {
		phi = 0
	} else {
		phi = atan2Safe(m, l)
	}
	theta = math.Asin(n)
	return phi, theta, true
}

// --------------------------------------------------------------------
// CSC — COBE quadrilateralized (Paper II §5.6.2, wcslib prj.c:7004)
// --------------------------------------------------------------------
//
// CSC uses two INDEPENDENT polynomial fits — one for forward, one for
// inverse — both from O'Neill & Laubscher 1976. The two polynomials
// are not exact inverses of each other; the FITS/wcslib spec
// documents a ~0.02 pixel residual between them by design (it's an
// area-equalization approximation, not an analytic inverse).
//
// We replicate wcslib's coefficients verbatim so our output matches
// the reference bit-for-bit. Users doing forward/inverse round trips
// through CSC should expect the same ~0.02 px drift that wcslib has.

// cscChiPsi evaluates the CSC inverse polynomial P(xf, yf) → chi, psi.
// Ported from wcslib cscx2s (prj.c:7173-7198). The polynomial is a
// degree-7 bivariate in xx = xf² and yy = yf² with 28 coefficients,
// evaluated as a sum of row polynomials z0..z6 then composed.
//
// wcslib uses single-precision (float) for the inner polynomial
// evaluation — we keep this to match wcslib bit-exactly, even though
// float64 would be more precise. CSC's ~0.02 px residual swamps the
// float/double difference by several orders of magnitude.
func cscChiPsi(xfIn, yfIn float64) (chi, psi float64) {
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
	xf := float32(xfIn)
	yf := float32(yfIn)
	xx := xf * xf
	yy := yf * yf

	// chi polynomial: bivariate in (xx, yy).
	z0 := p00 + xx*(p10+xx*(p20+xx*(p30+xx*(p40+xx*(p50+xx*p60)))))
	z1 := p01 + xx*(p11+xx*(p21+xx*(p31+xx*(p41+xx*p51))))
	z2 := p02 + xx*(p12+xx*(p22+xx*(p32+xx*p42)))
	z3 := p03 + xx*(p13+xx*(p23+xx*p33))
	z4 := p04 + xx*(p14+xx*p24)
	z5 := p05 + xx*p15
	z6 := p06
	chi32 := z0 + yy*(z1+yy*(z2+yy*(z3+yy*(z4+yy*(z5+yy*z6)))))
	chi32 = xf + xf*(1-xx)*chi32

	// psi polynomial: same structure with xx and yy swapped.
	z0 = p00 + yy*(p10+yy*(p20+yy*(p30+yy*(p40+yy*(p50+yy*p60)))))
	z1 = p01 + yy*(p11+yy*(p21+yy*(p31+yy*(p41+yy*p51))))
	z2 = p02 + yy*(p12+yy*(p22+yy*(p32+yy*p42)))
	z3 = p03 + yy*(p13+yy*(p23+yy*p33))
	z4 = p04 + yy*(p14+yy*p24)
	z5 = p05 + yy*p15
	z6 = p06
	psi32 := z0 + xx*(z1+xx*(z2+xx*(z3+xx*(z4+xx*(z5+xx*z6)))))
	psi32 = yf + yf*(1-yy)*psi32

	return float64(chi32), float64(psi32)
}

// cscXfYf evaluates the CSC forward polynomial (chi, psi) → (xf, yf).
// Different coefficient set from the inverse — this is a separate
// least-squares fit. Ported from wcslib cscs2x (prj.c:7279-7424).
func cscXfYf(chiIn, psiIn float64) (xf, yf float64) {
	const (
		gstar  = float32(1.37484847732)
		mm     = float32(0.004869491981)
		gamma  = float32(-0.13161671474)
		omega1 = float32(-0.159596235474)
		d0     = float32(0.0759196200467)
		d1     = float32(-0.0217762490699)
		c00    = float32(0.141189631152)
		c10    = float32(0.0809701286525)
		c01    = float32(-0.281528535557)
		c11    = float32(0.15384112876)
		c20    = float32(-0.178251207466)
		c02    = float32(0.106959469314)
	)
	chi := float32(chiIn)
	psi := float32(psiIn)
	chi2 := chi * chi
	psi2 := psi * psi
	chi2co := 1 - chi2
	psi2co := 1 - psi2

	// Avoid FP underflow.
	chipsi := float32(math.Abs(float64(chi * psi)))
	var chi4, psi4, chi2psi2 float32
	if chi2 > 1e-16 {
		chi4 = chi2 * chi2
	}
	if psi2 > 1e-16 {
		psi4 = psi2 * psi2
	}
	if chipsi > 1e-16 {
		chi2psi2 = chi2 * psi2
	}

	xf32 := chi * (chi2 + chi2co*(gstar+psi2*(gamma*chi2co+mm*chi2+
		psi2co*(c00+c10*chi2+c01*psi2+c11*chi2psi2+c20*chi4+c02*psi4))+
		chi2*(omega1-chi2co*(d0+d1*chi2))))
	yf32 := psi * (psi2 + psi2co*(gstar+chi2*(gamma*psi2co+mm*psi2+
		chi2co*(c00+c10*psi2+c01*chi2+c11*chi2psi2+c20*psi4+c02*chi4))+
		psi2*(omega1-psi2co*(d0+d1*psi2))))
	return float64(xf32), float64(yf32)
}

type cscProjection struct{}

func (cscProjection) Code() string    { return "CSC" }
func (cscProjection) Theta0() float64 { return 0 }

func (cscProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	sinTheta, cosTheta := math.Sincos(theta)
	sinPhi, cosPhi := math.Sincos(phi)
	l := cosTheta * cosPhi
	m := cosTheta * sinPhi
	n := sinTheta

	face, zeta := quadcubeFace(l, m, n)
	if zeta == 0 {
		return 0, 0, false
	}

	// Per-face xi, eta (from wcslib cscs2x:7365-7403).
	var xi, eta float64
	switch face {
	case 1:
		xi, eta = m, n
	case 2:
		xi, eta = -l, n
	case 3:
		xi, eta = -m, n
	case 4:
		xi, eta = l, n
	case 5:
		xi, eta = m, l
	default: // face 0
		xi, eta = m, -l
	}

	chi := xi / zeta
	psi := eta / zeta
	xf, yf := cscXfYf(chi, psi)

	// Clamp to [-1, 1].
	if xf > 1 {
		xf = 1
	} else if xf < -1 {
		xf = -1
	}
	if yf > 1 {
		yf = 1
	} else if yf < -1 {
		yf = -1
	}

	x0, y0 := quadcubeFaceOffset(face)
	x = quadW0 * (xf + x0)
	y = quadW0 * (yf + y0)
	return x, y, true
}

func (cscProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	face, xf, yf, fOK := quadcubeFaceFromXY(x, y)
	if !fOK {
		return 0, 0, false
	}

	chi, psi := cscChiPsi(xf, yf)

	// Reconstruct direction cosines (from wcslib cscx2s:7200-7232).
	t := 1 / math.Sqrt(chi*chi+psi*psi+1)
	var l, m, n float64
	switch face {
	case 1:
		l = t
		m = chi * l
		n = psi * l
	case 2:
		m = t
		l = -chi * m
		n = psi * m
	case 3:
		l = -t
		m = chi * l
		n = -psi * l
	case 4:
		m = -t
		l = -chi * m
		n = -psi * m
	case 5:
		n = -t
		l = -psi * n
		m = -chi * n
	default: // face 0
		n = t
		l = -psi * n
		m = chi * n
	}

	if l == 0 && m == 0 {
		phi = 0
	} else {
		phi = atan2Safe(m, l)
	}
	theta = math.Asin(n)
	return phi, theta, true
}

// --------------------------------------------------------------------
// QSC — quadrilateralized spherical cube (Paper II §5.6.3,
//       wcslib prj.c:7470)
// --------------------------------------------------------------------
//
// QSC is an exact equal-area projection (unlike CSC which is an
// approximation). The forward and inverse are closed-form but each
// has four sub-quadrant branches within every face, giving 24 total
// code paths. Ported from wcslib qscs2x (prj.c:7755) and qscx2s
// (prj.c:7504).

const sqrt2Inv = 0.70710678118654752440 // 1/sqrt(2)

type qscProjection struct{}

func (qscProjection) Code() string    { return "QSC" }
func (qscProjection) Theta0() float64 { return 0 }

func (qscProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	sinTheta, cosTheta := math.Sincos(theta)
	sinPhi, cosPhi := math.Sincos(phi)

	// Pole fast-path.
	if math.Abs(theta) == math.Pi/2 {
		return 0, math.Copysign(2*quadW0, theta), true
	}

	l := cosTheta * cosPhi
	m := cosTheta * sinPhi
	n := sinTheta

	face, zeta := quadcubeFace(l, m, n)

	// Per-face (xi, eta) assignment and small-angle zeco formula
	// (from wcslib qscs2x:7858-7931).
	var xi, eta, x0, y0 float64
	var zeco float64 // = 1 - zeta, used in the QSC area formula
	switch face {
	case 1:
		xi, eta = m, n
		x0, y0 = 0, 0
	case 2:
		xi, eta = -l, n
		x0, y0 = 2, 0
	case 3:
		xi, eta = -m, n
		x0, y0 = 4, 0
	case 4:
		xi, eta = l, n
		x0, y0 = 6, 0
	case 5:
		xi, eta = m, l
		x0, y0 = 0, -2
	default: // face 0
		xi, eta = m, -l
		x0, y0 = 0, 2
	}
	zeco = 1 - zeta

	// Small-angle formula: when zeta is near 1 (point near the face
	// center) the 1-zeta subtraction loses precision. Use a
	// second-order Taylor based on the face-local polar coordinates.
	if zeco < 1e-8 {
		var p float64
		switch face {
		case 1:
			t := theta
			p = math.Atan2(m, l) // = phi for face 1 since l = cos(θ)·cos(φ)
			zeco = (p*p + t*t) / 2
		case 2:
			t := theta
			p = math.Atan2(m, l) - math.Pi/2
			zeco = (p*p + t*t) / 2
		case 3:
			t := theta
			p = math.Atan2(m, l)
			p -= math.Copysign(math.Pi, p)
			zeco = (p*p + t*t) / 2
		case 4:
			t := theta
			p = math.Atan2(m, l) + math.Pi/2
			zeco = (p*p + t*t) / 2
		case 5:
			t := theta + math.Pi/2
			zeco = t * t / 2
		default: // face 0
			t := math.Pi/2 - theta
			zeco = t * t / 2
		}
	}
	// Apply the QSC area-equalization (wcslib qscs2x:7935-7957).
	// Four branches depending on which of ±xi, ±eta dominates.
	var xf, yf float64
	if xi != 0 || eta != 0 {
		switch {
		case -xi > math.Abs(eta):
			omega := eta / xi
			tau := 1 + omega*omega
			xf = -math.Sqrt(zeco / (1 - 1/math.Sqrt(1+tau)))
			yf = (xf / 15) * (radToDeg(math.Atan(omega)) - radToDeg(math.Asin(omega/math.Sqrt(tau+tau))))
		case xi > math.Abs(eta):
			omega := eta / xi
			tau := 1 + omega*omega
			xf = math.Sqrt(zeco / (1 - 1/math.Sqrt(1+tau)))
			yf = (xf / 15) * (radToDeg(math.Atan(omega)) - radToDeg(math.Asin(omega/math.Sqrt(tau+tau))))
		case -eta >= math.Abs(xi):
			omega := xi / eta
			tau := 1 + omega*omega
			yf = -math.Sqrt(zeco / (1 - 1/math.Sqrt(1+tau)))
			xf = (yf / 15) * (radToDeg(math.Atan(omega)) - radToDeg(math.Asin(omega/math.Sqrt(tau+tau))))
		case eta >= math.Abs(xi):
			omega := xi / eta
			tau := 1 + omega*omega
			yf = math.Sqrt(zeco / (1 - 1/math.Sqrt(1+tau)))
			xf = (yf / 15) * (radToDeg(math.Atan(omega)) - radToDeg(math.Asin(omega/math.Sqrt(tau+tau))))
		}
	}

	// Clamp to [-1, 1] face square.
	if math.Abs(xf) > 1 {
		xf = math.Copysign(1, xf)
	}
	if math.Abs(yf) > 1 {
		yf = math.Copysign(1, yf)
	}

	x = quadW0 * (xf + x0)
	y = quadW0 * (yf + y0)
	return x, y, true
}

func (qscProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	face, xf, yf, fOK := quadcubeFaceFromXY(x, y)
	if !fOK {
		return 0, 0, false
	}

	// Determine which of xf or yf drives the QSC parameter (the
	// "direct" branch uses xf, the indirect uses yf).
	direct := math.Abs(xf) > math.Abs(yf)
	var omega, tau, zeco, zeta, w float64
	const tol = 1e-12

	if direct {
		if xf == 0 {
			omega = 0
			tau = 1
			zeta = 1
			zeco = 0
		} else {
			// w is in wcslib's degree space where the face
			// half-width is 15 (= 45/3). In radians, we pass
			// (15 * yf/xf) directly through sind/cosd — the
			// numeric value of sind(15·yf/xf) in wcslib is the
			// same as sin(15·yf/xf · π/180) in our radians.
			// We use math.Sincos with the argument already
			// converted.
			wDeg := 15 * yf / xf
			sinW, cosW := math.Sincos(degToRad(wDeg))
			omega = sinW / (cosW - sqrt2Inv)
			tau = 1 + omega*omega
			zeco = xf * xf * (1 - 1/math.Sqrt(1+tau))
			zeta = 1 - zeco
		}
	} else {
		if yf == 0 {
			omega = 0
			tau = 1
			zeta = 1
			zeco = 0
		} else {
			wDeg := 15 * xf / yf
			sinW, cosW := math.Sincos(degToRad(wDeg))
			omega = sinW / (cosW - sqrt2Inv)
			tau = 1 + omega*omega
			zeco = yf * yf * (1 - 1/math.Sqrt(1+tau))
			zeta = 1 - zeco
		}
	}

	// Clamp zeta to [-1, 1].
	if zeta < -1 {
		if zeta < -1-tol {
			return 0, 0, false
		}
		zeta = -1
		zeco = 2
		w = 0
	} else {
		w = math.Sqrt(zeco * (2 - zeco) / tau)
	}

	// Reconstruct (l, m, n) per face (wcslib qscx2s:7657-7731).
	var l, m, n float64
	switch face {
	case 1:
		l = zeta
		if direct {
			m = w
			if xf < 0 {
				m = -m
			}
			n = m * omega
		} else {
			n = w
			if yf < 0 {
				n = -n
			}
			m = n * omega
		}
	case 2:
		m = zeta
		if direct {
			l = w
			if xf > 0 {
				l = -l
			}
			n = -l * omega
		} else {
			n = w
			if yf < 0 {
				n = -n
			}
			l = -n * omega
		}
	case 3:
		l = -zeta
		if direct {
			m = w
			if xf > 0 {
				m = -m
			}
			n = -m * omega
		} else {
			n = w
			if yf < 0 {
				n = -n
			}
			m = -n * omega
		}
	case 4:
		m = -zeta
		if direct {
			l = w
			if xf < 0 {
				l = -l
			}
			n = l * omega
		} else {
			n = w
			if yf < 0 {
				n = -n
			}
			l = n * omega
		}
	case 5:
		n = -zeta
		if direct {
			m = w
			if xf < 0 {
				m = -m
			}
			l = m * omega
		} else {
			l = w
			if yf < 0 {
				l = -l
			}
			m = l * omega
		}
	default: // face 0
		n = zeta
		if direct {
			m = w
			if xf < 0 {
				m = -m
			}
			l = -m * omega
		} else {
			l = w
			if yf > 0 {
				l = -l
			}
			m = -l * omega
		}
	}

	if l == 0 && m == 0 {
		phi = 0
	} else {
		phi = atan2Safe(m, l)
	}
	theta = math.Asin(n)
	return phi, theta, true
}
