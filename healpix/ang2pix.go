package healpix

import "math"

// Ang2Pix returns the pixel index containing the direction (theta, phi)
// on the unit sphere, in the given ordering scheme.
//
//   - theta: colatitude in radians, 0 at north pole, π at south pole
//   - phi: longitude in radians, [0, 2π)
//   - nside: HEALPix resolution parameter (any positive int for RING,
//     power-of-2 for NESTED)
//
// Port of astrometry.net xyztohp (healpix.c:771), which works via the
// intermediate XY scheme.
func Ang2Pix(nside int, theta, phi float64, scheme Scheme) int64 {
	sinTheta := math.Sin(theta)
	vx := sinTheta * math.Cos(phi)
	vy := sinTheta * math.Sin(phi)
	vz := math.Cos(theta)
	return Vec2Pix(nside, vx, vy, vz, scheme)
}

// Vec2Pix returns the pixel index containing the direction (vx, vy, vz)
// on the unit sphere (the vector need not be unit-length; only the
// direction matters).
func Vec2Pix(nside int, vx, vy, vz float64, scheme Scheme) int64 {
	xy := vec2xy(nside, vx, vy, vz)
	switch scheme {
	case Nested:
		return xyToNested(xy, nside)
	default:
		return xyToRing(xy, nside)
	}
}

// vec2xy converts a direction (vx, vy, vz) to an XY-scheme pixel index.
// Faithful port of astrometry.net xyztohp (healpix.c:771).
func vec2xy(nside int, vx, vy, vz float64) int64 {
	ns := float64(nside)
	twothirds := 2.0 / 3.0
	halfpi := math.Pi / 2

	// Normalize the direction vector.
	r := math.Sqrt(vx*vx + vy*vy + vz*vz)
	if r > 0 {
		vz /= r
	}

	// Compute phi in [0, 2π).
	phi := math.Atan2(vy, vx)
	if phi < 0 {
		phi += 2 * math.Pi
	}
	phiT := math.Mod(phi, halfpi) // phi within sector [0, π/2)

	var basehp, x, y int

	if vz >= twothirds || vz <= -twothirds {
		// ---------- Polar cap ----------
		north := vz >= twothirds
		zfactor := 1.0
		if !north {
			zfactor = -1.0
		}

		// Solve HEALPix paper eqns 19-20.
		root := (1.0 - vz*zfactor) * 3.0 * sqr(ns*(2.0*phiT-math.Pi)/math.Pi)
		kx := 0.0
		if root > 0 {
			kx = math.Sqrt(root)
		}
		root = (1.0 - vz*zfactor) * 3.0 * sqr(ns*2.0*phiT/math.Pi)
		ky := 0.0
		if root > 0 {
			ky = math.Sqrt(root)
		}

		var xx, yy float64
		if north {
			xx = ns - kx
			yy = ns - ky
		} else {
			xx = ky
			yy = kx
		}

		x = clampInt(int(math.Floor(xx)), 0, nside-1)
		y = clampInt(int(math.Floor(yy)), 0, nside-1)

		// Determine the sector column from phi.
		sector := (phi - phiT) / halfpi
		column := int(math.Round(sector))
		column = ((column % 4) + 4) % 4 // positive modulo: result in [0, 3]

		if north {
			basehp = column
		} else {
			basehp = 8 + column
		}
	} else {
		// ---------- Equatorial region ----------
		// Project into the unit square z=[-2/3, 2/3], phi=[0, π/2].
		zunits := (vz + twothirds) / (4.0 / 3.0)
		phiunits := phiT / halfpi
		// Diagonal coordinates in [0, 2].
		u1 := zunits + phiunits
		u2 := zunits - phiunits + 1.0
		// x is NE, y is NW.
		xx := u1 * ns
		yy := u2 * ns

		// Which sector column.
		sector := (phi - phiT) / halfpi
		offset := int(math.Round(sector))
		offset = ((offset % 4) + 4) % 4 // positive modulo: result in [0, 3]

		// Determine which quadrant (north polar / left eq / right eq / south polar).
		if xx >= ns {
			xx -= ns
			if yy >= ns {
				// North polar face.
				yy -= ns
				basehp = offset
			} else {
				// Right equatorial face.
				basehp = (offset+1)%4 + 4
			}
		} else {
			if yy >= ns {
				// Left equatorial face.
				yy -= ns
				basehp = offset + 4
			} else {
				// South polar face.
				basehp = 8 + offset
			}
		}

		x = clampInt(int(math.Floor(xx)), 0, nside-1)
		y = clampInt(int(math.Floor(yy)), 0, nside-1)
	}

	return xyCompose(basehp, x, y, nside)
}

func sqr(x float64) float64 { return x * x }

func clampInt(v, lo, hi int) int {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}
