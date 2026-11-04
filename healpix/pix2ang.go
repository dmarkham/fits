package healpix

import "math"

// Pix2Ang returns the (theta, phi) direction of the center of pixel pix.
//
//   - theta: colatitude in radians, [0, π]
//   - phi: longitude in radians, [0, 2π)
//
// Pixel centres are computed with the same double-precision formulas
// healpy uses; theta and phi agree with healpy.pix2ang to within 1 ULP
// (the libm acos/atan2 rounding difference), proven over every committed
// golden. Panics on an invalid nside/scheme or a pixel outside
// [0, 12·nside²).
func Pix2Ang(nside int, pix int64, scheme Scheme) (theta, phi float64) {
	g := newGeom("Pix2Ang", nside, scheme)
	g.checkPix("Pix2Ang", pix)
	z, phi, sth, haveSth := g.pix2loc(pix, scheme)
	if haveSth {
		return math.Atan2(sth, z), phi
	}
	return acos(z), phi
}

// acos is arccos accurate to about 1 ULP over [-1, 1]. Go's math.Acos is
// computed as π/2 - Asin(z), which loses several bits when the result is
// small (|z| close to 1) — up to 9 ULP at the ring-1 centre of Nside 5
// against healpy's glibc acos. For |z| > 0.5 the angle is instead taken
// from atan2(sqrt((1-z)(1+z)), |z|), where 1-|z| is exact (Sterbenz).
func acos(z float64) float64 {
	az := math.Abs(z)
	if az <= 0.5 {
		return math.Acos(z)
	}
	t := math.Atan2(math.Sqrt(float64((1-az)*(1+az))), az)
	if z < 0 {
		return pi - t
	}
	return t
}

// Pix2Vec returns the unit vector (x, y, z) pointing to the center of
// pixel pix. Components agree with healpy.pix2vec to within 2·2^-52
// absolute. Panics like Pix2Ang.
func Pix2Vec(nside int, pix int64, scheme Scheme) (vx, vy, vz float64) {
	g := newGeom("Pix2Vec", nside, scheme)
	g.checkPix("Pix2Vec", pix)
	z, phi, sth, haveSth := g.pix2loc(pix, scheme)
	if !haveSth {
		sth = math.Sqrt(float64((1 - z) * (1 + z)))
	}
	return float64(sth * math.Cos(phi)), float64(sth * math.Sin(phi)), z
}

// pix2loc returns z = cos(theta) and phi of the centre of pix. Near the
// poles (|z| > 0.99) it also returns sth = sin(theta) computed without
// cancellation, and haveSth = true.
func (g geom) pix2loc(pix int64, scheme Scheme) (z, phi, sth float64, haveSth bool) {
	ns := g.nside
	if scheme == Ring {
		switch {
		case pix < g.ncap: // north polar cap
			iring := (1 + isqrt(1+2*pix)) >> 1
			iphi := (pix + 1) - 2*iring*(iring-1)
			tmp := float64(float64(iring*iring) * g.fact2)
			z = 1.0 - tmp
			if z > 0.99 {
				sth, haveSth = math.Sqrt(float64(tmp*(2.-tmp))), true
			}
			phi = (float64(iphi) - 0.5) * halfpi / float64(iring)
		case pix < g.npix-g.ncap: // equatorial belt
			nl4 := 4 * ns
			ip := pix - g.ncap
			tmp := ip / nl4
			iring := tmp + ns
			iphi := ip - nl4*tmp + 1
			fodd := 0.5
			if (iring+ns)&1 != 0 {
				fodd = 1
			}
			z = float64(float64(2*ns-iring) * g.fact1)
			phi = float64((float64(iphi) - fodd) * pi * 0.75 * g.fact1)
		default: // south polar cap
			ip := g.npix - pix
			iring := (1 + isqrt(2*ip-1)) >> 1
			iphi := 4*iring + 1 - (ip - 2*iring*(iring-1))
			tmp := float64(float64(iring*iring) * g.fact2)
			z = tmp - 1.0
			if z < -0.99 {
				sth, haveSth = math.Sqrt(float64(tmp*(2.-tmp))), true
			}
			phi = (float64(iphi) - 0.5) * halfpi / float64(iring)
		}
		return z, phi, sth, haveSth
	}

	ix, iy, face := g.nest2xyf(pix)
	jr := jrll[face]<<g.order - int64(ix) - int64(iy) - 1
	var nr int64
	switch {
	case jr < ns:
		nr = jr
		tmp := float64(float64(nr*nr) * g.fact2)
		z = 1 - tmp
		if z > 0.99 {
			sth, haveSth = math.Sqrt(float64(tmp*(2.-tmp))), true
		}
	case jr > 3*ns:
		nr = ns*4 - jr
		tmp := float64(float64(nr*nr) * g.fact2)
		z = tmp - 1
		if z < -0.99 {
			sth, haveSth = math.Sqrt(float64(tmp*(2.-tmp))), true
		}
	default:
		nr = ns
		z = float64(float64(2*ns-jr) * g.fact1)
	}
	tmp := jpll[face]*nr + int64(ix) - int64(iy)
	if tmp < 0 {
		tmp += 8 * nr
	}
	if nr == ns {
		phi = float64(0.75 * halfpi * float64(tmp) * g.fact1)
	} else {
		phi = (0.5 * halfpi * float64(tmp)) / float64(nr)
	}
	return z, phi, sth, haveSth
}
