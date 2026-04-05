package transform

import "math"

// This file implements the native-spherical ↔ celestial-spherical rotation
// from Paper II §2.4. All angles are in radians.
//
// The rotation is parameterized by three numbers:
//
//   alphaP, deltaP : celestial coordinates of the native reference point
//   phiP           : native longitude of the celestial north pole (LONPOLE)
//
// # Numerical stability
//
// A naive implementation of the rotation equations uses asin() to extract
// theta (or delta) from its sine. When the target angle is close to ±pi/2
// (the image is near a celestial pole, OR — more commonly — the projection
// places the native pole at the reference point and the query pixel is
// close to the reference), the asin is ill-conditioned: its derivative is
// 1/sqrt(1 - y²), which is unbounded as |y| → 1. One ULP of error in
// sinTheta becomes amplified by ~1/sqrt(eps_mach) ≈ 1e8 in theta.
//
// The fix is to compute cos(theta) directly via hypot(num, den), where num
// and den are the numerator and denominator of the paired atan2 that yields
// phi. The identity num² + den² = cos²(theta) can be verified by
// substitution from Paper II eq. 2 / eq. 5 — the three equations express
// (phi, theta) as a unit vector decomposition, so squaring and summing
// gives the remaining component. With cos(theta) in hand, we use
// theta = atan2(sinTheta, cosTheta), which is stable to ULPs everywhere
// on the sphere.

// nativeToCelestial rotates from native (phi, theta) to celestial
// (alpha, delta). All angles in radians. Stable near both native and
// celestial poles.
func nativeToCelestial(phi, theta, alphaP, deltaP, phiP float64) (alpha, delta float64) {
	dphi := phi - phiP
	sinTheta, cosTheta := math.Sincos(theta)
	sinDp, cosDp := math.Sincos(deltaP)
	sinDphi, cosDphi := math.Sincos(dphi)

	// Paper II eq. 2, rearranged for stability.
	sinDelta := sinTheta*sinDp + cosTheta*cosDp*cosDphi

	num := -cosTheta * sinDphi
	den := sinTheta*cosDp - cosTheta*sinDp*cosDphi

	// cos(delta) = sqrt(num² + den²) — see package comment above for the
	// derivation. math.Hypot is the numerically stable way to compute this.
	cosDelta := math.Hypot(num, den)

	// atan2 here is stable regardless of whether delta is near ±pi/2.
	delta = math.Atan2(sinDelta, cosDelta)

	alpha = alphaP + math.Atan2(num, den)
	return alpha, delta
}

// celestialToNative is the inverse rotation. Uses the same hypot-based
// stability trick so that round-trips are precise to ULPs near the pole.
func celestialToNative(alpha, delta, alphaP, deltaP, phiP float64) (phi, theta float64) {
	dalpha := alpha - alphaP
	sinD, cosD := math.Sincos(delta)
	sinDp, cosDp := math.Sincos(deltaP)
	sinDA, cosDA := math.Sincos(dalpha)

	sinTheta := sinD*sinDp + cosD*cosDp*cosDA

	num := -cosD * sinDA
	den := sinD*cosDp - cosD*sinDp*cosDA

	cosTheta := math.Hypot(num, den)
	theta = math.Atan2(sinTheta, cosTheta)

	phi = phiP + math.Atan2(num, den)
	// Normalize phi into [-pi, pi] so downstream projection forward calls
	// see a bounded longitude. Without this, conic and pseudo-cylindrical
	// projections (which multiply phi by a non-unit factor) produce wildly
	// wrong intermediate coordinates for alpha values that wrap past 2*pi.
	for phi > math.Pi {
		phi -= 2 * math.Pi
	}
	for phi < -math.Pi {
		phi += 2 * math.Pi
	}
	return phi, theta
}
