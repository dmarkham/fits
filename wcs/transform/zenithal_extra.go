package transform

import (
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file contains the zenithal perspective family (AZP, SZP) and the
// zenithal polynomial (ZPN), separated from the simpler closed-form
// zenithal projections in zenithal.go.

// AZP — zenithal perspective projection (Paper II §5.1.1).
//
// Parameters (on the latitude axis, conventionally axis 2):
//
//	PV2_1 = mu    (dimensionless distance from sphere center to
//	              projection source, in sphere radii; default 0)
//	PV2_2 = gamma (tilt angle in degrees, default 0)
//
// Forward: Paper II eq. 20.
//
//	R_theta = (mu + 1) * cos(theta) / (mu + sin(theta))
//
// (The gamma tilt adds a cos(gamma) denominator term when rotating the
// projection plane; see eq. 19.) For gamma = 0 this reduces to the
// un-tilted zenithal perspective.
//
// Inverse: quadratic in sin(theta). We derive it below from the forward.
type azpProjection struct {
	mu    float64
	gamma float64 // tilt, radians
}

func newAZP(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return azpProjection{
		mu:    pvFloat(pv, latAxis, 1, 0),
		gamma: degToRad(pvFloat(pv, latAxis, 2, 0)),
	}
}

func (azpProjection) Code() string    { return "AZP" }
func (azpProjection) Theta0() float64 { return math.Pi / 2 }

func (p azpProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	sinPhi, cosPhi := math.Sincos(phi)
	sinTheta, cosTheta := math.Sincos(theta)
	sinGamma, cosGamma := math.Sincos(p.gamma)
	denom := p.mu + sinTheta + cosTheta*cosPhi*sinGamma
	if denom == 0 {
		return 0, 0, false
	}
	r := (p.mu + 1) * cosTheta / denom
	x = r * sinPhi
	y = -r * cosPhi / cosGamma
	return x, y, true
}

func (p azpProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// For gamma = 0 (the common case with no tilt) we have a closed-form
	// quadratic inverse. Setting u = sin(theta) and substituting into
	//
	//	R = (mu+1) * cos(theta) / (mu + sin(theta))
	//
	// gives the quadratic
	//
	//	u² * (R² + (mu+1)²) + 2*u*R²*mu + (R²*mu² - (mu+1)²) = 0
	//
	// whose larger root is the near-hemisphere solution. For gamma != 0
	// the formula still has a closed form but involves cos(phi) coupling
	// into the denominator; we handle gamma != 0 by scaling y first then
	// applying the gamma=0 derivation (Paper II eq. 23 with the tilt
	// correction absorbed into y).
	cosGamma := math.Cos(p.gamma)
	if cosGamma == 0 {
		return 0, 0, false
	}
	yTilt := y * cosGamma
	r := math.Hypot(x, yTilt)
	if r == 0 {
		return 0, math.Pi / 2, true
	}
	phi = math.Atan2(x, -yTilt)

	// Quadratic coefficients. Without tilt, this is exact; with gamma != 0
	// it is an approximation that is exact along the principal meridians.
	muPlus1 := p.mu + 1
	a := r*r + muPlus1*muPlus1
	b := 2 * r * r * p.mu
	c := r*r*p.mu*p.mu - muPlus1*muPlus1
	disc := b*b - 4*a*c
	if disc < 0 {
		return 0, 0, false
	}
	sqrtDisc := math.Sqrt(disc)
	// Larger root → near hemisphere.
	u := (-b + sqrtDisc) / (2 * a)
	if u < -1 || u > 1 {
		u = (-b - sqrtDisc) / (2 * a)
		if u < -1 || u > 1 {
			return 0, 0, false
		}
	}
	sinTheta := u
	cosTheta := math.Sqrt((1 - u) * (1 + u))
	theta = math.Atan2(sinTheta, cosTheta)
	return phi, theta, true
}

// SZP — slant zenithal perspective (Paper II §5.1.2).
//
// Parameters (latitude axis):
//
//	PV2_1 = mu      (source distance in sphere radii; default 0)
//	PV2_2 = phi_c   (native longitude of the projection source, degrees)
//	PV2_3 = theta_c (native latitude of the projection source, degrees,
//	                default 90)
//
// This is the most general zenithal perspective: the projection source
// can be anywhere along the axis of a cone whose apex is at theta_c.
// Used by some planetary stereo-pair instruments. When mu = 0 and
// theta_c = 90 it reduces to TAN; when theta_c = 90 and mu != 0 it
// reduces to AZP.
type szpProjection struct {
	mu     float64
	phiC   float64 // radians
	thetaC float64 // radians
}

func newSZP(pv map[wcs.PVKey]float64, latAxis int) Projection {
	return szpProjection{
		mu:     pvFloat(pv, latAxis, 1, 0),
		phiC:   degToRad(pvFloat(pv, latAxis, 2, 0)),
		thetaC: degToRad(pvFloat(pv, latAxis, 3, 90)),
	}
}

func (szpProjection) Code() string    { return "SZP" }
func (szpProjection) Theta0() float64 { return math.Pi / 2 }

func (p szpProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	// Paper II eq. 24:
	//   X_p = mu * cos(theta_c) * sin(phi_c)
	//   Y_p = -mu * cos(theta_c) * cos(phi_c)
	//   Z_p = mu * sin(theta_c) + 1
	// And the projection (Paper II eq. 25):
	//   X = (Z_p * cos(theta)*sin(phi) - X_p*(1 - sin(theta))) / (Z_p - (1 - sin(theta)))
	//   Y = -(Z_p * cos(theta)*cos(phi) + Y_p*(1 - sin(theta))) / (Z_p - (1 - sin(theta)))
	sinPhiC, cosPhiC := math.Sincos(p.phiC)
	sinThC, cosThC := math.Sincos(p.thetaC)
	Xp := p.mu * cosThC * sinPhiC
	Yp := -p.mu * cosThC * cosPhiC
	Zp := p.mu*sinThC + 1

	sinTheta, cosTheta := math.Sincos(theta)
	sinPhi, cosPhi := math.Sincos(phi)
	s := 1 - sinTheta
	denom := Zp - s
	if denom == 0 {
		return 0, 0, false
	}
	x = (Zp*cosTheta*sinPhi - Xp*s) / denom
	y = -(Zp*cosTheta*cosPhi + Yp*s) / denom
	return x, y, true
}

func (p szpProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// Closed-form quadratic inverse. Port of wcslib's szpx2s
	// (reference/cfitsio not applicable; see Siril subprojects/wcslib/
	// prj.c:955). The SZP inverse solves a quadratic in sin(theta):
	//
	//	a * sin²(θ) + 2b * sin(θ) + c = 0
	//
	// with coefficients derived from the direction cosines of the
	// projection source relative to the reference point. In wcslib's
	// default unit system r0 = 180/π, giving w[0] = π/180 (the
	// deg→rad factor). Our library is already radian-native, so we
	// set r0 = 1 which collapses w[0] = 1 and the degree-scale
	// constants w[4..6] to their unscaled counterparts w[1..3].
	//
	// The three SZP-specific constants (direction cosines of the
	// projection source on the unit sphere) are:
	//
	//	w1 = -μ * cos(θc) * sin(φc)   (= -Xp in Paper II eq. 24)
	//	w2 =  μ * cos(θc) * cos(φc)   (= -Yp)
	//	w3 =  μ * sin(θc) + 1         (=  Zp)
	sinPhiC, cosPhiC := math.Sincos(p.phiC)
	sinThC, cosThC := math.Sincos(p.thetaC)
	w1 := -p.mu * cosThC * sinPhiC
	w2 := p.mu * cosThC * cosPhiC
	w3 := p.mu*sinThC + 1
	if w3 == 0 {
		return 0, 0, false
	}

	xr, yr := x, y
	r2 := xr*xr + yr*yr

	x1 := (xr - w1) / w3
	y1 := (yr - w2) / w3
	xy := xr*x1 + yr*y1

	const tol = 1e-13
	var z, sinthe float64
	if r2 < 1e-10 {
		// Small-angle formula near the pole: avoids the
		// 1-sinthe cancellation by computing z = r²/2 directly
		// (exact leading term of the Taylor expansion). Matches
		// wcslib's "Use small angle formula" branch.
		theta = math.Pi/2 - math.Sqrt(r2/(1+xy))
		sinthe = math.Sin(theta)
		z = r2 / 2
	} else {
		// Full quadratic.
		t := x1*x1 + y1*y1
		a := t + 1
		b := xy - t
		c := r2 - 2*xy + t - 1
		d := b*b - a*c
		if d < 0 {
			return 0, 0, false
		}
		d = math.Sqrt(d)

		// Two candidate solutions. wcslib picks the one closest
		// to the pole (larger sin(θ)).
		sinth1 := (-b + d) / a
		sinth2 := (-b - d) / a
		if sinth1 > sinth2 {
			sinthe = sinth1
		} else {
			sinthe = sinth2
		}
		if sinthe > 1 {
			if sinthe-1 < tol {
				sinthe = 1
			} else {
				// Numerical overshoot above 1: try the other root.
				if sinth1 < sinth2 {
					sinthe = sinth1
				} else {
					sinthe = sinth2
				}
			}
		}
		if sinthe < -1 {
			if sinthe+1 > -tol {
				sinthe = -1
			}
		}
		if sinthe > 1 || sinthe < -1 {
			return 0, 0, false
		}
		theta = math.Asin(sinthe)
		z = 1 - sinthe
	}

	phi = atan2Safe(xr-x1*z, -(yr - y1*z))
	return phi, theta, true
}

// ZPN — zenithal polynomial projection (Paper II §5.1.8).
//
// R_theta is a user-defined polynomial in zeta = π/2 - theta:
//
//	R(ζ) = pv[0] + pv[1]·ζ + pv[2]·ζ² + ... + pv[n]·ζⁿ
//
// with up to 30 coefficients. The polynomial is monotonic in ζ from
// ζ=0 (pole) out to an inflection point where the derivative first
// becomes non-positive; beyond that the inverse is undefined.
//
// Forward is straight polynomial evaluation. Inverse is a port of
// wcslib's zpnx2s / zpnset (reference implementation, prj.c:2241 and
// prj.c:2342): three-case dispatch on degree with closed-form
// solutions for k=1 and k=2, and weighted false-position bracketing
// for k≥3. The bracketing range is established by zpnset, which walks
// the polynomial's derivative outward from ζ=0 to find the first
// inflection point and caches it as the valid upper bound.
type zpnProjection struct {
	coeffs [30]float64 // pv[0..29]
	degree int         // highest non-zero index (n in the math above)

	// Inflection-point constants cached by setup (wcslib w[0], w[1]).
	// zdMax is the ζ at the inflection (upper bracket for the inverse
	// search); rMax is R(zdMax), the upper R bound. For degree < 2
	// there is no inflection and zdMax is π (the whole hemisphere).
	zdMax float64
	rMax  float64
}

func newZPN(pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	p := zpnProjection{}
	for i := 0; i < 30; i++ {
		p.coeffs[i] = pvFloat(pv, latAxis, i, 0)
		if p.coeffs[i] != 0 {
			p.degree = i
		}
	}
	if p.degree == 0 && p.coeffs[0] == 0 {
		return nil, fmt.Errorf("wcs/transform: ZPN requires at least one non-zero PV%d_i coefficient", latAxis)
	}
	if err := p.setup(); err != nil {
		return nil, err
	}
	return p, nil
}

// setup mirrors wcslib's zpnset (prj.c:2241). For degree < 2 there's no
// inflection, so the whole ζ ∈ [0, π] range is valid. For degree ≥ 2
// we walk the derivative outward from ζ=0 in 1° steps until it first
// goes non-positive, then refine the inflection with secant iteration,
// then cache (zdMax, R(zdMax)) as the upper bound for the inverse
// bracket.
func (p *zpnProjection) setup() error {
	k := p.degree
	if k < 2 {
		p.zdMax = math.Pi
		return nil
	}
	if p.coeffs[1] <= 0 {
		return fmt.Errorf("wcs/transform: ZPN pv[1] must be > 0 (got %g)", p.coeffs[1])
	}

	// Walk in 1° steps looking for the first derivative sign change.
	var zd1, d1 float64
	var zd2, d2 float64
	d1 = p.coeffs[1]
	zd1 = 0
	found := false
	for j := 0; j < 180; j++ {
		zd2 = float64(j) * math.Pi / 180
		d2 = 0
		for m := k; m > 0; m-- {
			d2 = d2*zd2 + float64(m)*p.coeffs[m]
		}
		if d2 <= 0 {
			found = true
			break
		}
		zd1 = zd2
		d1 = d2
	}

	var zd, r float64
	if !found {
		// Derivative never went negative within [0, π] — whole
		// hemisphere is valid.
		zd = math.Pi
	} else {
		// Refine the zero-of-derivative by secant iteration.
		const tol = 1e-13
		for j := 1; j <= 10; j++ {
			zd = zd1 - d1*(zd2-zd1)/(d2-d1)
			d := 0.0
			for m := k; m > 0; m-- {
				d = d*zd + float64(m)*p.coeffs[m]
			}
			if math.Abs(d) < tol {
				break
			}
			if d < 0 {
				zd2 = zd
				d2 = d
			} else {
				zd1 = zd
				d1 = d
			}
		}
	}
	// R at the inflection (upper R bound for inverse).
	r = 0
	for m := k; m >= 0; m-- {
		r = r*zd + p.coeffs[m]
	}
	p.zdMax = zd
	p.rMax = r
	return nil
}

func (zpnProjection) Code() string    { return "ZPN" }
func (zpnProjection) Theta0() float64 { return math.Pi / 2 }

// polyEval evaluates the ZPN polynomial R(ζ) via Horner's rule.
func (p zpnProjection) polyEval(zeta float64) float64 {
	v := p.coeffs[p.degree]
	for i := p.degree - 1; i >= 0; i-- {
		v = v*zeta + p.coeffs[i]
	}
	return v
}

func (p zpnProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	zeta := math.Pi/2 - theta
	r := p.polyEval(zeta)
	return r * math.Sin(phi), -r * math.Cos(phi), true
}

func (p zpnProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	r := math.Hypot(x, y)
	if r == 0 {
		phi = 0
	} else {
		phi = atan2Safe(x, -y)
	}

	k := p.degree
	var zd float64
	const tol = 1e-13

	switch {
	case k < 1:
		// Constant polynomial — no solution.
		return 0, 0, false

	case k == 1:
		// Linear: R = pv[0] + pv[1]·ζ
		zd = (r - p.coeffs[0]) / p.coeffs[1]

	case k == 2:
		// Quadratic: pv[2]·ζ² + pv[1]·ζ + (pv[0] - r) = 0
		a := p.coeffs[2]
		b := p.coeffs[1]
		c := p.coeffs[0] - r
		d := b*b - 4*a*c
		if d < 0 {
			return 0, 0, false
		}
		d = math.Sqrt(d)
		// Choose the root closest to the pole (smaller ζ).
		zd1 := (-b + d) / (2 * a)
		zd2 := (-b - d) / (2 * a)
		if zd1 < zd2 {
			zd = zd1
		} else {
			zd = zd2
		}
		if zd < -tol {
			// Try the other root.
			if zd1 > zd2 {
				zd = zd1
			} else {
				zd = zd2
			}
		}
		if zd < 0 {
			if zd < -tol {
				return 0, 0, false
			}
			zd = 0
		} else if zd > math.Pi {
			if zd > math.Pi+tol {
				return 0, 0, false
			}
			zd = math.Pi
		}

	default:
		// Higher order: weighted false-position bracketing on
		// [0, zdMax]. The bracket is valid because the polynomial
		// is monotonic on this range by construction (zdMax is
		// the first inflection point, located by setup()).
		zd1 := 0.0
		r1 := p.coeffs[0]
		zd2 := p.zdMax
		r2 := p.rMax

		// Outside the valid monotonic range → no solution.
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
			// Dissect the interval with false-position,
			// clamping lambda to avoid stalls. wcslib caps
			// the loop at 100 iterations.
			for j := 0; j < 100; j++ {
				lambda := (r2 - r) / (r2 - r1)
				if lambda < 0.1 {
					lambda = 0.1
				} else if lambda > 0.9 {
					lambda = 0.9
				}
				zd = zd2 - lambda*(zd2-zd1)
				rt := p.polyEval(zd)
				if rt < r {
					if r-rt < tol {
						break
					}
					r1 = rt
					zd1 = zd
				} else {
					if rt-r < tol {
						break
					}
					r2 = rt
					zd2 = zd
				}
				if math.Abs(zd2-zd1) < tol {
					break
				}
			}
		}
	}

	theta = math.Pi/2 - zd
	return phi, theta, true
}
