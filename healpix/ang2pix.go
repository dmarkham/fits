package healpix

import (
	"fmt"
	"math"
)

// Ang2Pix returns the pixel index containing the direction (theta, phi)
// on the unit sphere, in the given ordering scheme. It evaluates
// healpy.ang2pix's formulas with Go's math.Cos and math.Sin, so it agrees
// with healpy exactly except at directions within about an ULP of a pixel
// boundary, where the last bit of cos or sin decides the pixel (see the
// package documentation for the test matrix).
//
//   - nside: 1..2^29 for RING (any value), power of 2 for NESTED
//   - theta: colatitude in radians, [0, π]
//   - phi: longitude in radians; any finite value (taken modulo 2π)
//
// Ang2Pix panics on an invalid nside, NESTED with a non-power-of-2 nside,
// theta outside [0, π] (including NaN) or a non-finite phi.
func Ang2Pix(nside int, theta, phi float64, scheme Scheme) int64 {
	g := newGeom("Ang2Pix", nside, scheme)
	if !(theta >= 0 && theta <= pi) {
		panic(fmt.Sprintf("healpix.Ang2Pix: theta %g outside [0, pi]", theta))
	}
	if math.IsNaN(phi) || math.IsInf(phi, 0) {
		panic(fmt.Sprintf("healpix.Ang2Pix: phi %g is not finite", phi))
	}
	sth := 0.0
	if angPolar(theta) {
		sth = math.Sin(theta)
	}
	return g.ang2pixCS(theta, phi, math.Cos(theta), sth, scheme)
}

// angPolar reports whether Ang2Pix takes the sin(theta) formula: close to
// the poles cos(theta) loses the information that sin(theta) still
// carries (the same switch-over points healpix_cxx / healpy use).
func angPolar(theta float64) bool { return theta < 0.01 || theta > 3.14159-0.01 }

// ang2pixCS is Ang2Pix after the libm calls: cth = cos(theta) and, where
// angPolar(theta), sth = sin(theta). Tests feed it correctly rounded
// values to check the arithmetic and the branch choice exactly.
func (g geom) ang2pixCS(theta, phi, cth, sth float64, scheme Scheme) int64 {
	if angPolar(theta) {
		return g.loc2pix(cth, phi, sth, true, scheme)
	}
	return g.loc2pix(cth, phi, 0, false, scheme)
}

// Vec2Pix returns the pixel index containing the direction (vx, vy, vz).
// The vector need not be unit length; only its direction matters:
// Vec2Pix(v * 2^k) == Vec2Pix(v) for every k that scales v exactly. It
// first scales the vector by a power of 2 so its largest component is in
// [0.5, 1) (exact, the direction is unchanged), then evaluates healpy's
// vec2pix formula, with phi = math.Atan2(vy, vx) from the unscaled
// components. For vectors of about unit length this gives
// healpy.vec2pix's answer except within about an ULP of a pixel boundary,
// where the last bit of atan2 decides the pixel; healpy itself, which does not
// rescale, returns wrong pixels once squaring the components overflows or
// loses bits in the subnormal range (e.g. for (1e-200, 0, 1e-200) RING
// gives pixel 0 and NESTED another wrong pixel at every Nside > 1 instead
// of the pixel of (1, 0, 1); near-pole vectors scaled by 2^-482 start to
// move at Nside 2^29 - 1, and by tens of pixels at 2^-490). It panics on the
// zero vector, a NaN or infinite component, or an invalid nside/scheme
// (see Ang2Pix).
func Vec2Pix(nside int, vx, vy, vz float64, scheme Scheme) int64 {
	g := newGeom("Vec2Pix", nside, scheme)
	m := math.Max(math.Abs(vx), math.Max(math.Abs(vy), math.Abs(vz)))
	if math.IsNaN(vx) || math.IsNaN(vy) || math.IsNaN(vz) || math.IsInf(m, 0) {
		panic(fmt.Sprintf("healpix.Vec2Pix: vector (%g, %g, %g) is not finite", vx, vy, vz))
	}
	if m == 0 {
		panic("healpix.Vec2Pix: the zero vector has no direction")
	}
	phi := 0.0
	if vx != 0 || vy != 0 {
		phi = math.Atan2(vy, vx) // scale-free, so from the exact inputs
	}
	return g.vec2pixPhi(vx, vy, vz, m, phi, scheme)
}

// vec2pixPhi is Vec2Pix after its one libm call: phi = atan2(vy, vx) (0
// when vx = vy = 0) and m the largest |component|. Tests feed it a
// correctly rounded phi to check the arithmetic and the branch choice
// exactly.
func (g geom) vec2pixPhi(vx, vy, vz, m, phi float64, scheme Scheme) int64 {
	// Rescaling every vector is simpler than deciding when it is needed.
	// (Skipping it for m in [2^-480, 2^480] would give the same pixels
	// for Nside <= 2^29: no square overflows, and x²+y² loses bits only
	// when sin(theta) < 2^-31.5, deep inside polar ring 1, where the pixel
	// depends on phi alone, which comes from the unscaled components.
	// healpy, which never rescales, moves 0 pixels at 2^-480 and 2^-481,
	// one at 2^-482 and tens at 2^-490: vec2pix_scaled_healpy.json. A
	// Vec2Pix that rescales small vectors only when m < 2^-482 (or any
	// lower threshold) fails TestVec2PixScaledPolar, which checks every
	// exponent -512..-470; a threshold of 2^-481 or higher is equivalent.)
	_, e := math.Frexp(m)
	vx, vy, vz = math.Ldexp(vx, -e), math.Ldexp(vy, -e), math.Ldexp(vz, -e)
	length := math.Sqrt(float64(vx*vx) + float64(vy*vy) + float64(vz*vz))
	xl := 1 / length
	nz := float64(vz * xl)
	if math.Abs(nz) > 0.99 {
		return g.loc2pix(nz, phi, float64(math.Sqrt(float64(vx*vx)+float64(vy*vy))*xl), true, scheme)
	}
	return g.loc2pix(nz, phi, 0, false, scheme)
}

// fmodulo returns v1 mod v2 in [0, v2) for v2 > 0.
func fmodulo(v1, v2 float64) float64 {
	if v1 >= 0 {
		if v1 < v2 {
			return v1
		}
		return math.Mod(v1, v2)
	}
	tmp := math.Mod(v1, v2) + v2
	if tmp == v2 {
		return 0
	}
	return tmp
}

// loc2pix maps (z = cos theta, phi) to a pixel. sth = sin theta is used
// in the polar caps when haveSth is set (|z| close to 1).
func (g geom) loc2pix(z, phi, sth float64, haveSth bool, scheme Scheme) int64 {
	ns := g.nside
	fns := float64(ns)
	za := math.Abs(z)
	tt := fmodulo(float64(phi*invHalfpi), 4.0) // in [0, 4)

	if scheme == Ring {
		if za <= twothird { // equatorial belt
			nl4 := 4 * ns
			temp1 := float64(fns * (0.5 + tt))
			temp2 := float64(fns * z * 0.75)
			jp := int64(temp1 - temp2) // ascending edge line index
			jm := int64(temp1 + temp2) // descending edge line index
			ir := ns + 1 + jp - jm     // ring number counted from z=2/3, in 1..2n+1
			kshift := 1 - (ir & 1)
			t1 := jp + jm - ns + kshift + 1 + nl4 + nl4
			ip := (t1 >> 1) % nl4
			return g.ncap + (ir-1)*nl4 + ip
		}
		// polar caps
		tp := tt - float64(int64(tt))
		var tmp float64
		if za < 0.99 || !haveSth {
			tmp = float64(fns * math.Sqrt(float64(3*(1-za))))
		} else {
			tmp = fns * sth / math.Sqrt((1+za)/3)
		}
		jp := int64(tp * tmp)
		jm := int64((1.0 - tp) * tmp)
		ir := jp + jm + 1 // ring number counted from the closest pole
		ip := int64(tt * float64(ir))
		if z > 0 {
			return 2*ir*(ir-1) + ip
		}
		return g.npix - 2*ir*(ir+1) + ip
	}

	// NESTED
	if za <= twothird {
		temp1 := float64(fns * (0.5 + tt))
		temp2 := float64(fns * (z * 0.75))
		jp := int64(temp1 - temp2)
		jm := int64(temp1 + temp2)
		ifp := jp >> g.order
		ifm := jm >> g.order
		var face int
		switch {
		case ifp == ifm:
			face = int(ifp | 4)
		case ifp < ifm:
			face = int(ifp)
		default:
			face = int(ifm + 8)
		}
		ix := int(jm & (ns - 1))
		iy := int(ns - (jp & (ns - 1)) - 1)
		return g.xyf2nest(ix, iy, face)
	}
	ntt := int(tt)
	if ntt > 3 {
		ntt = 3
	}
	tp := tt - float64(ntt)
	var tmp float64
	if za < 0.99 || !haveSth {
		tmp = float64(fns * math.Sqrt(float64(3*(1-za))))
	} else {
		tmp = fns * sth / math.Sqrt((1+za)/3)
	}
	jp := int64(tp * tmp)
	jm := int64((1.0 - tp) * tmp)
	// healpix_cxx clamps jp/jm to nside-1 "for points too close to the
	// boundary". With IEEE doubles the clamp cannot trigger: in the z
	// branch za > twothird gives 3*(1-za) <= 1-2^-52, so tmp <= nside*(1-2^-53),
	// and tp <= 1-2^-53, so tp*tmp rounds to at most nside*(1-2^-52) <
	// nside; in the sin(theta) branch tmp < 0.18*nside. It is kept so the
	// code stays line-for-line comparable with healpix_cxx.
	if jp > ns-1 {
		jp = ns - 1
	}
	if jm > ns-1 {
		jm = ns - 1
	}
	if z >= 0 {
		return g.xyf2nest(int(ns-jm-1), int(ns-jp-1), ntt)
	}
	return g.xyf2nest(int(jp), int(jm), ntt+8)
}
