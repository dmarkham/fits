package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file implements the polyconic and pseudoconic projections from
// Paper II §5.5.

// BON — Bonne's projection (§5.5.1).
//
// Parameter: PV2_1 = theta_1 (standard parallel in degrees; 0 reduces to
// SFL). Paper II eq. 79:
//
//	Y_0 = cot(theta_1) + theta_1           (in radians)
//	A   = phi * cos(theta) / (Y_0 - theta)
//	x   = (Y_0 - theta) * sin(A)
//	y   = -(Y_0 - theta) * cos(A) + Y_0
//
// Equal-area, pseudoconic, theta_0 = theta_1.
type bonProjection struct {
	theta1 float64 // radians
	Y0     float64
}

func newBON(pv map[wcs.PVKey]float64, latAxis int) Projection {
	t1 := degToRad(pvFloat(pv, latAxis, 1, 0))
	return bonProjection{
		theta1: t1,
		Y0:     safeCot(t1) + t1,
	}
}

// safeCot returns cot(x) = cos(x)/sin(x), or 0 when sin(x) is zero.
// Returning 0 is safe because BON's Forward/Inverse special-case
// theta1==0 to delegate to SFL before safeCot is called with x=0.
func safeCot(x float64) float64 {
	s := math.Sin(x)
	if s == 0 {
		return 0
	}
	return math.Cos(x) / s
}

func (bonProjection) Code() string      { return "BON" }
func (p bonProjection) Theta0() float64 { return p.theta1 }

func (p bonProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if p.theta1 == 0 {
		// Degenerate case: Bonne with theta_1 = 0 is sinusoidal (SFL).
		return sflProjection{}.Forward(phi, theta)
	}
	R := p.Y0 - theta
	A := phi * math.Cos(theta) / R
	sA, cA := math.Sincos(A)
	x = R * sA
	y = -R*cA + p.Y0
	return x, y, true
}

func (p bonProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	if p.theta1 == 0 {
		return sflProjection{}.Inverse(x, y)
	}
	yPrime := p.Y0 - y
	R := math.Hypot(x, yPrime)
	if p.theta1 < 0 {
		R = -R
	}
	theta = p.Y0 - R
	A := math.Atan2(x/R, yPrime/R)
	c := math.Cos(theta)
	if c == 0 {
		return 0, theta, true
	}
	phi = A * R / c
	return phi, theta, true
}

// PCO — polyconic (Paper II §5.5.2).
//
// Forward map:
//
//	theta == 0:  x = phi,  y = 0
//	otherwise:   E = phi * sin(theta)
//	             x = cot(theta) * sin(E)
//	             y = theta + cot(theta) * (1 - cos(E))
//
// Inverse map: there is no closed form. We solve for theta via weighted
// false-position bracketing (regula falsi with clamped lambda), matching
// wcslib's pcox2s algorithm in prj.c:6395. This approach is preferred over
// Newton iteration because:
//
//   - It is guaranteed to converge (the solution is always bracketed).
//   - It needs no Jacobian — only the sign of the residue function.
//   - A closed-form Taylor branch handles the cot(theta) singularity at
//     theta = 0 without the general solver ever touching it.
//
// The residue we drive to zero is derived from the forward map:
//
//	f(theta) = x^2 + (y-theta) * ((y-theta) - 2*cot(theta))
//
// Equivalent to 0 iff (phi, theta) is the correct inverse. (Derivation:
// square x, add (y-theta)^2, use the identity sin^2 E + (1-cos E)^2 =
// 2*(1-cos E), then substitute cot(theta)*(1-cos E) = y-theta.)
//
// The library's coordinates are in radians throughout, so wcslib's
// scaling constants w[0]..w[3] collapse to w[0]=1, w[2]=2, w[3]=0.5 in
// our unit system.
type pcoProjection struct{}

func (pcoProjection) Code() string    { return "PCO" }
func (pcoProjection) Theta0() float64 { return 0 }

func (pcoProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	if theta == 0 {
		return phi, 0, true
	}
	sinT, cosT := math.Sincos(theta)
	if sinT == 0 {
		return phi, 0, true
	}
	cot := cosT / sinT
	E := phi * sinT
	sE, cE := math.Sincos(E)
	x = cot * sE
	y = theta + cot*(1-cE)
	return x, y, true
}

func (pcoProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	// Exact-equator fast path.
	if y == 0 {
		return x, 0, true
	}
	// Pole: |y| == pi/2 exactly. wcslib does this check in degrees with
	// tol = 1e-12; in radians the equivalent threshold is ~1.7e-14.
	const tol = 1e-12
	const nearPoleTol = 1.7e-14
	w := math.Abs(y)
	if math.Abs(w-math.Pi/2) < nearPoleTol {
		return 0, math.Copysign(math.Pi/2, y), true
	}

	xx := x * x
	var ymthe, tanthe float64

	if w < 1.0e-6 {
		// Closed-form Taylor branch for small |theta|. Avoids the
		// cot(theta) singularity in the general solver. Derivation:
		// for small theta, y ~ theta * (1 + 0.5*x^2), so
		//    theta ~ y / (1 + 0.5*x^2)
		// matches wcslib's w[0] + w[3]*xj*xj in radian-space.
		theta = y / (1 + 0.5*xx)
		ymthe = y - theta
		tanthe = math.Tan(theta)
	} else {
		// Weighted false-position bracketing. Bracket:
		//    thepos = y       (upper bound: polyconic correction >= 0 so theta <= y)
		//    theneg = 0
		// Seed fneg = -fpos = -xx so the first iteration is a pure
		// midpoint split (lambda = 0.5). This avoids ever evaluating
		// the residue at theta = 0, where cot(theta) is undefined.
		thepos := y
		theneg := 0.0
		fpos := xx
		fneg := -xx

		for k := 0; k < 64; k++ {
			// Regula falsi weight, clamped to [0.1, 0.9] to avoid
			// degenerate cases where one endpoint dominates and the
			// iteration stalls near a non-root.
			lambda := fpos / (fpos - fneg)
			if lambda < 0.1 {
				lambda = 0.1
			} else if lambda > 0.9 {
				lambda = 0.9
			}
			theta = thepos - lambda*(thepos-theneg)
			ymthe = y - theta
			tanthe = math.Tan(theta)
			// Residue f(theta) = x^2 + (y-theta)*((y-theta) - 2*cot(theta))
			f := xx + ymthe*(ymthe-2.0/tanthe)
			if math.Abs(f) < tol {
				break
			}
			if math.Abs(thepos-theneg) < tol {
				break
			}
			if f > 0 {
				thepos = theta
				fpos = f
			} else {
				theneg = theta
				fneg = f
			}
		}
	}

	// Back out phi from the forward-map identities:
	//    x1 = r0 - (y-theta)*tan(theta)   [with r0 = 1 in radian units]
	//    y1 = x*tan(theta)
	//    phi = atan2(y1, x1) / sin(theta)
	x1 := 1.0 - ymthe*tanthe
	y1 := x * tanthe
	if x1 == 0 && y1 == 0 {
		phi = 0
	} else {
		phi = math.Atan2(y1, x1) / math.Sin(theta)
	}
	return phi, theta, true
}
