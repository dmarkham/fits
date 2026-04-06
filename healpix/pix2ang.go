package healpix

import "math"

// Pix2Ang returns the (theta, phi) direction of the center of pixel pix
// in the given ordering scheme.
//
//   - theta: colatitude in radians, 0 at north pole, π at south pole
//   - phi: longitude in radians, [0, 2π)
//
// Port of astrometry.net hp_to_xyz (healpix.c:1006).
func Pix2Ang(nside int, pix int64, scheme Scheme) (theta, phi float64) {
	vx, vy, vz := Pix2Vec(nside, pix, scheme)
	// For z near ±1 (polar pixels), acos(z) loses precision because
	// its derivative diverges at z=1. Use atan2(sin_theta, z) instead,
	// which is stable everywhere. sin_theta = sqrt(vx²+vy²) = the
	// cylindrical radius, which we already computed in Pix2Vec.
	sintheta := math.Sqrt(vx*vx + vy*vy)
	theta = math.Atan2(sintheta, vz)
	phi = math.Atan2(vy, vx)
	if phi < 0 {
		phi += 2 * math.Pi
	}
	return theta, phi
}

// Pix2Vec returns the (x, y, z) unit vector pointing to the center of
// pixel pix in the given ordering scheme.
func Pix2Vec(nside int, pix int64, scheme Scheme) (vx, vy, vz float64) {
	var xy int64
	switch scheme {
	case Nested:
		xy = nestedToXY(pix, nside)
	default:
		xy = ringToXY(pix, nside)
	}
	return xyToVec(xy, nside)
}

// xyToVec converts an XY-scheme pixel index to a unit vector (pixel
// center). Faithful port of astrometry.net hp_to_xyz (healpix.c:1006),
// called with dx=0.5, dy=0.5 for pixel centers.
func xyToVec(xy int64, nside int) (vx, vy, vz float64) {
	hp := xyDecompose(xy, nside)
	chp := hp.Face
	ns := float64(nside)

	// Pixel-center position: xp + 0.5, yp + 0.5.
	x := float64(hp.X) + 0.5
	y := float64(hp.Y) + 0.5

	equatorial := true
	zfactor := 1.0

	// Check polar/equatorial boundary (healpix.c:1023-1034).
	// North polar faces (0-3): polar region is where x+y > nside.
	// South polar faces (8-11): polar region is where x+y < nside.
	if chp <= 3 {
		if x+y > ns {
			equatorial = false
			zfactor = 1.0
		}
	}
	if chp >= 8 {
		if x+y < ns {
			equatorial = false
			zfactor = -1.0
		}
	}

	var z, phi float64

	if equatorial {
		// Equatorial path (healpix.c:1036-1060).
		var zoff, phioff float64
		x /= ns
		y /= ns

		if chp <= 3 {
			// North polar face, equatorial part.
			phioff = 1.0
			// zoff stays 0, chp unchanged.
		} else if chp <= 7 {
			// Equatorial face.
			zoff = -1.0
			chp -= 4
		} else {
			// South polar face, equatorial part.
			phioff = 1.0
			zoff = -2.0
			chp -= 8
		}

		z = (2.0 / 3.0) * (x + y + zoff)
		phi = (math.Pi / 4.0) * (x - y + phioff + 2*float64(chp))
	} else {
		// Polar path (healpix.c:1062-1104).
		if zfactor == -1 {
			// South polar: swap x,y and mirror.
			x, y = ns-y, ns-x
		}

		var phiT float64
		if y == ns && x == ns {
			phiT = 0
		} else {
			phiT = math.Pi * (ns - y) / (2.0 * ((ns - x) + (ns - y)))
		}

		// Solve for z from either equation depending on phi_t.
		if phiT < math.Pi/4 {
			z = 1.0 - sqr(math.Pi*(ns-x)/((2*phiT-math.Pi)*ns))/3.0
		} else {
			z = 1.0 - sqr(math.Pi*(ns-y)/(2*phiT*ns))/3.0
		}
		z *= zfactor

		// Base face determines the phi offset.
		if hp.Face >= 8 {
			phi = math.Pi/2*float64(hp.Face-8) + phiT
		} else {
			phi = math.Pi/2*float64(hp.Face) + phiT
		}
	}

	// Normalize phi to [0, 2π).
	if phi < 0 {
		phi += 2 * math.Pi
	}

	rad := math.Sqrt(1 - z*z)
	sinPhi, cosPhi := math.Sincos(phi)
	return rad * cosPhi, rad * sinPhi, z
}
