package transform

import (
	"errors"
	"fmt"
	"math"
)

// This file implements the native-spherical ↔ celestial-spherical rotation
// from Paper II §2.4. All angles are in radians.
//
// The rotation is parameterized by three numbers:
//
//   alphaP, deltaP : celestial coordinates of the native pole (wcslib
//                    euler[0] and 90 - euler[1]); see celPole
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
// celestial poles. sinDp, cosDp are sin and cos of deltaP.
func nativeToCelestial(phi, theta, alphaP, sinDp, cosDp, phiP float64) (alpha, delta float64) {
	dphi := phi - phiP
	sinTheta, cosTheta := lmSincos(theta)
	sinDphi, cosDphi := lmSincos(dphi)

	// Paper II eq. 2, rearranged for stability.
	sinDelta := float64(sinTheta*sinDp) + float64(cosTheta*cosDp*cosDphi)

	num := float64(-cosTheta * sinDphi)
	den := float64(sinTheta*cosDp) - float64(cosTheta*sinDp*cosDphi)

	// cos(delta) = sqrt(num² + den²) — see package comment above for the
	// derivation. |num|, |den| <= 1, so the squares cannot overflow, and a
	// square that underflows moves delta only within 1e-150 rad of +-pi/2;
	// math.Hypot, which guards against both, costs twice as much here.
	cosDelta := math.Sqrt(float64(num*num) + float64(den*den))

	// atan2 here is stable regardless of whether delta is near ±pi/2.
	delta = lmAtan2(sinDelta, cosDelta)

	alpha = alphaP + lmAtan2(num, den)
	return alpha, delta
}

// celestialToNative is the inverse rotation. Uses the same hypot-based
// stability trick so that round-trips are precise to ULPs near the pole.
func celestialToNative(alpha, delta, alphaP, sinDp, cosDp, phiP float64) (phi, theta float64) {
	dalpha := alpha - alphaP
	sinD, cosD := lmSincos(delta)
	sinDA, cosDA := lmSincos(dalpha)

	sinTheta := float64(sinD*sinDp) + float64(cosD*cosDp*cosDA)

	num := float64(-cosD * sinDA)
	den := float64(sinD*cosDp) - float64(cosD*sinDp*cosDA)

	cosTheta := math.Sqrt(float64(num*num) + float64(den*den)) // see nativeToCelestial
	theta = lmAtan2(sinTheta, cosTheta)

	phi = phiP + lmAtan2(num, den)
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

// ---------------------------------------------------------------------------
// Native-pole solution (Paper II §2.4), ported from wcslib 8.10 cel.c
// celset() so that every branch, default and rejection matches the
// reference implementation. Inputs and outputs are in DEGREES:
//
//	lng0, lat0     celestial coordinates of the fiducial point (CRVAL)
//	phi0, theta0   native coordinates of the fiducial point
//	lonpole        LONPOLE (phi_p); ignored unless lonpoleSet
//	latpole        LATPOLE, used to choose between two solutions (default 90)
//
// It returns the celestial coordinates of the native pole (lngp, latp) and
// the native longitude of the celestial pole phip, which feed
// nativeToCelestial / celestialToNative as (alphaP, deltaP, phiP).
// It fails, as wcslib does, when no native-pole latitude puts the fiducial
// point at (phi0, theta0) (e.g. CAR with CRVAL2 = 30 and LONPOLE = 180).

const celTol = 1e-10 // wcslib cel.c tol

func celPole(lng0, lat0, phi0, theta0, lonpole float64, lonpoleSet bool, latpole float64) (lngp, latp, phip float64, err error) {
	phip = lonpole
	if !lonpoleSet || lonpole == 999 {
		if lat0 < theta0 {
			phip = 180
		} else {
			phip = 0
		}
		phip += phi0
		if phip < -180 {
			phip += 360
		} else if phip > 180 {
			phip -= 360
		}
	}
	thetap := latpole

	if theta0 == 90 {
		// Fiducial point at the native pole.
		return lng0, lat0, phip, nil
	}

	slat0, clat0 := sincosd(lat0)
	sthe0, cthe0 := sincosd(theta0)
	var sphip, cphip, u, v float64
	latpreq := 0
	if phip == phi0 {
		sphip, cphip = 0, 1
		u = theta0
		v = 90 - lat0
	} else {
		sphip, cphip = sincosd(phip - phi0)
		x := float64(cthe0 * cphip)
		y := sthe0
		z := math.Sqrt(float64(x*x) + float64(y*y))
		if z == 0 {
			if slat0 != 0 {
				return 0, 0, 0, fmt.Errorf("wcs/transform: invalid celestial parameters: |lat0| must be 0 when theta0 = 0 and |phip - phi0| = 90 (lat0 = %g)", lat0)
			}
			// latp determined solely by LATPOLE in this case.
			latpreq = 2
			switch {
			case thetap > 90:
				latp = 90
			case thetap < -90:
				latp = -90
			default:
				latp = thetap
			}
		} else {
			slz := slat0 / z
			if math.Abs(slz) > 1 {
				if math.Abs(slz)-1 < celTol {
					slz = math.Copysign(1, slz)
				} else {
					return 0, 0, 0, fmt.Errorf("wcs/transform: ill-conditioned celestial parameters: |slz| = %g > 1 (CRVAL lat %g, phip %g, phi0 %g, theta0 %g)", math.Abs(slz), lat0, phip, phi0, theta0)
				}
			}
			u = atan2d(y, x)
			v = acosd(slz)
		}
	}

	if latpreq == 0 {
		latp1 := u + v
		if latp1 > 180 {
			latp1 -= 360
		} else if latp1 < -180 {
			latp1 += 360
		}
		latp2 := u - v
		if latp2 > 180 {
			latp2 -= 360
		} else if latp2 < -180 {
			latp2 += 360
		}
		if math.Abs(thetap-latp1) < math.Abs(thetap-latp2) {
			if math.Abs(latp1) < 90+celTol {
				latp = latp1
			} else {
				latp = latp2
			}
		} else {
			if math.Abs(latp2) < 90+celTol {
				latp = latp2
			} else {
				latp = latp1
			}
		}
		// Account for rounding error.
		if math.Abs(latp) < 90+celTol {
			if latp > 90 {
				latp = 90
			} else if latp < -90 {
				latp = -90
			}
		}
	}

	z := float64(cosd(latp) * clat0)
	if math.Abs(z) < celTol {
		switch {
		case math.Abs(clat0) < celTol:
			// Celestial pole at the fiducial point.
			lngp = lng0
		case latp > 0:
			// Celestial north pole at the native pole.
			lngp = lng0 + phip - phi0 - 180
		default:
			// Celestial south pole at the native pole.
			lngp = lng0 - phip + phi0
		}
	} else {
		x := (sthe0 - float64(sind(latp)*slat0)) / z
		y := sphip * cthe0 / clat0
		if x == 0 && y == 0 {
			// wcslib's guard; no input reaching it is known: y = 0 needs
			// phip - phi0 = 0 or 180, and x = 0 then needs
			// cos(latp) cos(lat0) = 0 (the |z| < tol branch above),
			// up to round-off.
			return 0, 0, 0, errors.New("wcs/transform: invalid celestial parameters (native pole undetermined)")
		}
		lngp = lng0 - atan2d(y, x)
	}

	// Make celestial longitude of the native pole the same sign as at the
	// fiducial point (as wcslib; a 360 deg shift of lngp is not visible
	// through PixelToSky, which normalises alpha, or SkyToPixel).
	if lng0 >= 0 {
		if lngp < 0 {
			lngp += 360
		} else if lngp > 360 {
			lngp -= 360
		}
	} else {
		if lngp > 0 {
			lngp -= 360
		} else if lngp < -360 {
			lngp += 360
		}
	}

	if math.Abs(latp) > 90+celTol {
		return 0, 0, 0, fmt.Errorf("wcs/transform: ill-conditioned celestial parameters: no valid native-pole latitude for LONPOLE %g, phi0 %g, theta0 %g, CRVAL lat %g", phip, phi0, theta0, lat0)
	}
	return lngp, latp, phip, nil
}

// eulerAngles are wcslib's celprm euler[5] (cel.c celset): the celestial
// longitude of the native pole, its colatitude 90 - latp, the native
// longitude of the celestial pole, and the cosine and sine of the second,
// all in degrees.
type eulerAngles [5]float64

func newEulerAngles(lngp, latp, phip float64) eulerAngles {
	e := eulerAngles{lngp, 90 - latp, phip}
	e[4], e[3] = sincosd(e[1])
	return e
}

// sphTol is wcslib sph.c tol.
const sphTol = 1.0e-5

// sphS2X is wcslib 8.10 sphs2x for one point, in degrees: celestial (lng,
// lat) to native (phi, theta), with phi normalised to [-180, 180]. It is
// ported branch for branch (the aligned-pole special cases, the
// rearranged formula for |x| < 1e-5, the change of longitude origin when
// the point is at the native pole, and the acos form of theta near the
// poles) so that sky->pixel matches wcslib also where the native longitude
// is degenerate: at a native pole every phi is the same point, and wcslib
// picks one by convention (a cylindrical projection maps the native pole
// to a whole row of pixels, so the convention decides the pixel).
// Not every branch is observable: the exact dlng mod 180 = 0 latitude
// form and the sphTol threshold of the rearranged formula only change
// round-off (mutants that drop or widen them pass every golden), and the
// eul[1] >= 90 native-pole longitude origin is equivalent (a point at the
// native pole has dlng = 0 or 180 there, so -dlng and dlng give the same
// phi); see docs/VALIDATION.md.
func sphS2X(eul *eulerAngles, lng, lat float64) (phi, theta float64) {
	norm := func(p float64) float64 {
		if p > 180 {
			p -= 360
		} else if p < -180 {
			p += 360
		}
		return p
	}
	if eul[4] == 0 {
		if eul[1] == 0 {
			// Simple change in origin of longitude.
			dphi := math.Mod(eul[2]-180-eul[0], 360)
			return norm(math.Mod(lng+dphi, 360)), lat
		}
		// Pole-flip with change in origin of longitude.
		dphi := math.Mod(eul[2]+eul[0], 360)
		return norm(math.Mod(dphi-lng, 360)), -lat
	}
	dlng := lng - eul[0]
	sinlat, coslat := sincosd(lat)
	coslat3 := float64(coslat * eul[3])
	coslat4 := float64(coslat * eul[4])
	sinlat3 := float64(sinlat * eul[3])
	sinlat4 := float64(sinlat * eul[4])
	sinlng, coslng := sincosd(dlng)

	// Native longitude.
	x := sinlat4 - float64(coslat3*coslng)
	if math.Abs(x) < sphTol {
		// Rearranged formula to reduce roundoff errors.
		x = -cosd(lat+eul[1]) + float64(coslat3*(1-coslng))
	}
	y := float64(-coslat * sinlng)
	var dphi float64
	if x != 0 || y != 0 {
		dphi = atan2d(y, x)
	} else if eul[1] < 90 {
		// Change of origin of longitude.
		dphi = dlng - 180
	} else {
		dphi = -dlng
	}
	phi = norm(math.Mod(eul[2]+dphi, 360))

	// Native latitude.
	if math.Mod(dlng, 180) == 0 {
		theta = lat + float64(coslng*eul[1])
		if theta > 90 {
			theta = 180 - theta
		}
		if theta < -90 {
			theta = -180 - theta
		}
	} else {
		z := sinlat3 + float64(coslat4*coslng)
		if math.Abs(z) > 0.99 {
			// Alternative formula for greater accuracy.
			theta = math.Copysign(acosd(math.Sqrt(float64(x*x)+float64(y*y))), z)
		} else {
			theta = asind(z)
		}
	}
	return phi, theta
}
