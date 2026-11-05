package transform

import "math"

// Pseudo-cylindrical and conventional projections, Paper II §5.3: SFL,
// PAR, MOL, AIT. All have theta_0 = 0 and map the sphere onto a bounded
// region of the plane; the inverse rejects every plane point outside that
// region (the corners of an all-sky map, beyond the ellipse of MOL/AIT,
// beyond the meridian bounds of SFL/PAR). Ports of the wcslib 8.10 prj.c
// functions of the same names, in degrees (see projection.go).

// SFL — Sanson-Flamsteed: x = phi cos(theta), y = theta.
type sflProjection struct{}

func (sflProjection) Code() string    { return "SFL" }
func (sflProjection) theta0() float64 { return 0 }

func (sflProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	return float64(cosd(theta) * phi), theta, true
}

func (sflProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	s := lmCos(y / r2d)
	if s == 0.0 {
		return 0, 0, false
	}
	phi = float64(x * (1.0 / s))
	theta = y
	return phi, theta, prjbchk(1.0e-12, &phi, &theta)
}

// PAR — parabolic:
//
//	x = phi (2 cos(2 theta/3) - 1),  y = 180 sin(theta/3)
type parProjection struct{}

func (parProjection) Code() string    { return "PAR" }
func (parProjection) theta0() float64 { return 0 }

func (parProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	s := sind(theta / 3.0)
	return float64((1.0 - float64(4.0*s*s)) * phi), float64(180.0 * s), true
}

func (parProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	r := float64(y * (1.0 / 180.0))
	if r > 1.0 || r < -1.0 {
		return 0, 0, false
	}
	s := 1.0 - float64(4.0*r*r)
	if s == 0.0 {
		// The poles: on the projection only at x = 0.
		if math.Abs(x)-tol >= 0.0 {
			return 0, 0, false
		}
	} else {
		s = 1.0 / s
	}
	phi = float64(x * s)
	theta = float64(3.0 * asind(r))
	return phi, theta, prjbchk(1.0e-12, &phi, &theta)
}

// MOL — Mollweide:
//
//	x = (2 sqrt(2)/pi) phi cos(gamma),  y = sqrt(2) r0 sin(gamma)
//	where 2 gamma + sin(2 gamma) = pi sin(theta)
//
// The forward solves for gamma by bisection to 1e-13, as wcslib does.
type molProjection struct{}

var (
	molW0 = float64(math.Sqrt2 * r2d) // SQRT2*r0
	molW1 = molW0 / 90.0
	molW2 = 1.0 / molW0
	molW3 = 90.0 / r2d
	molW4 = 2.0 / wpi
)

func (molProjection) Code() string    { return "MOL" }
func (molProjection) theta0() float64 { return 0 }

func (molProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	const tol = 1.0e-13
	var xi, eta float64
	switch {
	case math.Abs(theta) == 90.0:
		xi = 0.0
		eta = math.Copysign(molW0, theta)
	case theta == 0.0:
		xi = 1.0
		eta = 0.0
	default:
		v := molSolve(float64(wpi*sind(theta)), tol)
		gamma := v / 2.0
		s, c := lmSincos(gamma)
		xi = c
		eta = float64(molW0 * s)
	}
	return float64(xi * (molW1 * phi)), eta, true
}

// molSolve is wcslib's MOL forward bisection for v = 2 gamma,
// v + sin(v) = u, step for step: start at v = u, then bisect [-pi, pi]
// until |(v - u) + sin(v)| < tol. Where v is outside (lo, hi) the branch
// it takes is known without evaluating the residual (molSkipBounds), so
// only the last few steps call sin; the iterate it stops at is the same
// (TestMolSolveMatchesBisection).
func molSolve(u, tol float64) float64 {
	lo, hi := molSkipBounds(u, tol)
	v0 := -wpi
	v1 := wpi
	v := u
	for k := 0; k < 100; k++ {
		switch {
		case v >= hi:
			v1 = v
		case v <= lo:
			v0 = v
		default:
			resid := (v - u) + lmSin(v)
			if resid < 0.0 {
				if resid > -tol {
					return v
				}
				v0 = v
			} else {
				if resid < tol {
					return v
				}
				v1 = v
			}
		}
		v = (v0 + v1) / 2.0
	}
	return v
}

// molSkipBounds returns lo < hi such that wcslib's MOL bisection on
// resid(v) = (v - u) + sin(v) takes the "v0 = v" branch at every v <= lo
// and the "v1 = v" branch at every v >= hi without stopping, so those
// steps need no sin. The bisection itself (about 46 halvings of
// [-pi, pi] down to |resid| < tol) is kept exactly, so the solution is
// the same iterate wcslib stops at; evaluating sin only next to the root
// is what makes it cheap.
//
// Why it is safe: g(v) = v + sin(v) is non-decreasing (g' = 1 + cos v >= 0),
// and the computed resid is within 4e-15 of g(v) - u for |v| <= pi. The
// bounds are taken around a Newton estimate of the root and then checked
// by evaluating resid there: resid(hi) >= 2 tol implies resid(v) > tol
// for every v >= hi, and resid(lo) <= -2 tol implies resid(v) < -tol for
// every v <= lo, whatever the accuracy of the estimate. If the check
// fails the margin is widened, and as a last resort (-inf, +inf) is
// returned, which evaluates every step.
func molSkipBounds(u, tol float64) (lo, hi float64) {
	resid := func(v float64) float64 { return (v - u) + math.Sin(v) }
	// Newton on g(v) = u. From u/2, which is at or below the root for
	// u >= 0 (sin v <= v), g's concavity on [0, pi] makes the iterates
	// increase monotonically; near the pole (|u| > 3) start from the
	// cubic expansion g(v) ~ pi - (pi - v)^3/6.
	au := math.Abs(u)
	v := float64(0.5 * au)
	if au > 3 {
		v = wpi - math.Cbrt(float64(6*(wpi-au)))
	}
	for i := 0; i < 8; i++ {
		s, c := math.Sincos(v)
		d := 1 + c
		if d < 1e-300 {
			break
		}
		dv := ((v - au) + s) / d
		v -= dv
		if v > wpi {
			v = wpi
		}
		if math.Abs(dv) < 1e-15 {
			break
		}
	}
	if u < 0 {
		v = -v
	}
	m := 4 * tol / (1 + math.Cos(v))
	if cm := math.Cbrt(float64(48 * tol)); !(m < cm) {
		m = cm // next to +-pi, where g - u grows like (v - root)^3/6
	}
	for try := 0; try < 6; try++ {
		lo, hi = v-m, v+m
		okLo := lo <= -wpi || resid(lo) <= float64(-2*tol)
		okHi := hi >= wpi || resid(hi) >= float64(2*tol)
		if okLo && okHi {
			if lo <= -wpi {
				lo = math.Inf(-1)
			}
			if hi >= wpi {
				hi = math.Inf(1)
			}
			return lo, hi
		}
		m *= 4
	}
	return math.Inf(-1), math.Inf(1)
}

func (molProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-12
	sx := float64(molW3 * x)
	tx := math.Abs(x) - tol
	y0 := y / r2d
	r := 2.0 - float64(y0*y0)
	istat := 0
	var s float64
	if r <= tol {
		if r < -tol {
			istat = 1
		} else {
			istat = -1
		}
		r = 0.0
		s = 0.0
	} else {
		r = math.Sqrt(r)
		s = 1.0 / r
	}
	z := float64(y * molW2)
	if math.Abs(z) > 1.0 {
		if math.Abs(z) > 1.0+tol {
			z = 0.0
			istat = 1
		} else {
			z = math.Copysign(1.0, z) + y0*r/wpi
		}
	} else {
		z = float64(libAsin(z)*molW4) + y0*r/wpi
	}
	if math.Abs(z) > 1.0 {
		if math.Abs(z) > 1.0+tol {
			z = 0.0
			istat = 1
		} else {
			z = math.Copysign(1.0, z)
		}
	}
	if istat > 0 || (istat < 0 && tx >= 0.0) {
		return 0, 0, false
	}
	phi = float64(sx * s)
	theta = asind(z)
	return phi, theta, prjbchk(1.0e-11, &phi, &theta)
}

// AIT — Hammer-Aitoff:
//
//	gamma = r0 sqrt(2 / (1 + cos(theta) cos(phi/2)))
//	x = 2 gamma cos(theta) sin(phi/2),  y = gamma sin(theta)
//
// The inverse rejects every point outside the bounding ellipse
// (x/(4 r0))^2 + (y/(2 r0))^2 > 1/2.
type aitProjection struct{}

var (
	aitW0 = float64(2.0 * r2d * r2d)
	aitW1 = 1.0 / (2.0 * aitW0)
	aitW2 = aitW1 / 4.0
	aitW3 = 1.0 / (2.0 * r2d)
)

func (aitProjection) Code() string    { return "AIT" }
func (aitProjection) theta0() float64 { return 0 }

func (aitProjection) s2x(phi, theta float64) (x, y float64, ok bool) {
	sinphi, cosphi := sincosd(phi / 2.0)
	sinthe, costhe := sincosd(theta)
	w := math.Sqrt(aitW0 / (1.0 + float64(costhe*cosphi)))
	return float64(2.0 * w * costhe * sinphi), float64(w * sinthe), true
}

func (aitProjection) x2s(x, y float64) (phi, theta float64, ok bool) {
	const tol = 1.0e-13
	s := 1.0 - float64(x*x*aitW2)
	t := float64(x * aitW3)
	s -= float64(y * y * aitW1)
	if s < 0.5 {
		if s < 0.5-tol {
			return 0, 0, false
		}
		s = 0.5
	}
	z := math.Sqrt(s)
	x0 := float64(2.0*z*z) - 1.0
	y0 := float64(z * t)
	if x0 != 0.0 || y0 != 0.0 {
		phi = float64(2.0 * atan2d(y0, x0))
	}
	t = z * y / r2d
	if math.Abs(t) > 1.0 {
		if math.Abs(t) > 1.0+tol {
			return 0, 0, false
		}
		theta = math.Copysign(90.0, t)
	} else {
		theta = asind(t)
	}
	return phi, theta, prjbchk(1.0e-13, &phi, &theta)
}
