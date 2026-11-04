package healpix

import (
	"fmt"
	"math"
	"math/bits"
)

// MaxNside is the largest supported Nside, 2^29 (healpy / healpix_cxx
// order_max for 64-bit pixel indices; 12·4^29 ≈ 3.46e18 pixels).
const MaxNside = 1 << 29

// Floating-point constants used by the pixel formulas. They are package
// variables, not constants, on purpose: Go evaluates constant expressions
// exactly, while the reference (healpix_cxx, used by healpy) evaluates
// products such as 0.75*halfpi in IEEE double. Using variables makes Go
// do the same double-precision operations in the same order.
var (
	pi        = 3.141592653589793
	halfpi    = 1.5707963267948966
	invHalfpi = 0.6366197723675814
	twothird  = 2.0 / 3.0
)

// Face layout tables: jrll[f] is the ring index (in units of nside) of the
// southernmost corner of base face f, jpll[f] its longitude index (in
// units of pi/4). Gorski et al. 2005, section 4.
var (
	jrll = [12]int64{2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4}
	jpll = [12]int64{1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7}
)

// geom holds the per-Nside quantities every pixel formula needs.
type geom struct {
	nside  int64
	order  int // log2(nside), or -1 when nside is not a power of 2
	npface int64
	ncap   int64 // pixels in the north polar cap
	npix   int64
	fact1  float64
	fact2  float64
}

// newGeom validates nside for scheme and returns its geometry. Invalid
// input panics: these functions have no error return, and a silently
// wrong pixel index is worse than a crash.
func newGeom(fn string, nside int, scheme Scheme) geom {
	if nside <= 0 || nside > MaxNside {
		panic(fmt.Sprintf("healpix.%s: nside %d out of range [1, 2^29]", fn, nside))
	}
	order := -1
	if IsNsidePow2(nside) {
		order = bits.TrailingZeros64(uint64(nside))
	}
	if scheme == Nested && order < 0 {
		panic(fmt.Sprintf("healpix.%s: NESTED ordering requires a power-of-2 nside, got %d", fn, nside))
	}
	if scheme != Nested && scheme != Ring {
		panic(fmt.Sprintf("healpix.%s: unknown scheme %d", fn, int(scheme)))
	}
	ns := int64(nside)
	g := geom{nside: ns, order: order, npface: ns * ns}
	g.ncap = (g.npface - ns) << 1
	g.npix = 12 * g.npface
	g.fact2 = 4 / float64(g.npix)
	g.fact1 = float64(float64(ns<<1) * g.fact2)
	return g
}

func (g geom) checkPix(fn string, pix int64) {
	if pix < 0 || pix >= g.npix {
		panic(fmt.Sprintf("healpix.%s: pixel %d out of range [0, %d) for nside %d", fn, pix, g.npix, g.nside))
	}
}

// isqrt returns floor(sqrt(v)) exactly for 0 <= v < 2^63. The float
// square root is exact below 2^50; above that it is corrected by one step.
func isqrt(v int64) int64 {
	r := int64(math.Sqrt(float64(v) + 0.5))
	if v < 1<<50 {
		return r
	}
	if r*r > v {
		r--
	} else if (r+1)*(r+1) <= v {
		r++
	}
	return r
}

// ringInfoSmall returns the first pixel index, the pixel count and the
// shifted flag of 1-based ring ring.
func (g geom) ringInfoSmall(ring int64) (startpix, ringpix int64, shifted bool) {
	switch {
	case ring < g.nside:
		return 2 * ring * (ring - 1), 4 * ring, true
	case ring < 3*g.nside:
		ringpix = 4 * g.nside
		return g.ncap + (ring-g.nside)*ringpix, ringpix, (ring-g.nside)&1 == 0
	default:
		nr := 4*g.nside - ring
		return g.npix - 2*nr*(nr+1), 4 * nr, true
	}
}

// ring2xyf decomposes a RING index into (x, y, face) in O(1): the ring
// number comes from an integer square root in the polar caps and a
// division in the equatorial belt.
func (g geom) ring2xyf(pix int64) (ix, iy, face int) {
	ns := g.nside
	nl2 := 2 * ns
	var iring, iphi, kshift, nr int64
	switch {
	case pix < g.ncap: // north polar cap
		iring = (1 + isqrt(1+2*pix)) >> 1
		iphi = (pix + 1) - 2*iring*(iring-1)
		kshift = 0
		nr = iring
		face = int((iphi - 1) / nr)
	case pix < g.npix-g.ncap: // equatorial belt
		ip := pix - g.ncap
		tmp := ip / (4 * ns)
		iring = tmp + ns
		iphi = ip - tmp*4*ns + 1
		kshift = (iring + ns) & 1
		nr = ns
		ire := tmp + 1
		irm := nl2 + 1 - tmp
		ifm := (iphi - (ire >> 1) + ns - 1) / ns
		ifp := (iphi - (irm >> 1) + ns - 1) / ns
		switch {
		case ifp == ifm:
			face = int(ifp | 4)
		case ifp < ifm:
			face = int(ifp)
		default:
			face = int(ifm + 8)
		}
	default: // south polar cap
		ip := g.npix - pix
		iring = (1 + isqrt(2*ip-1)) >> 1
		iphi = 4*iring + 1 - (ip - 2*iring*(iring-1))
		kshift = 0
		nr = iring
		iring = 2*nl2 - iring
		face = int((iphi-1)/nr) + 8
	}
	irt := iring - (2+int64(face>>2))*ns + 1
	ipt := 2*iphi - jpll[face]*nr - kshift - 1
	if ipt >= nl2 {
		ipt -= 8 * ns
	}
	return int((ipt - irt) >> 1), int((-ipt - irt) >> 1), face
}

// xyf2ring is the inverse of ring2xyf.
func (g geom) xyf2ring(ix, iy, face int) int64 {
	nl4 := 4 * g.nside
	jr := jrll[face]*g.nside - int64(ix) - int64(iy) - 1
	nBefore, nr, shifted := g.ringInfoSmall(jr)
	nr >>= 2
	kshift := int64(1)
	if shifted {
		kshift = 0
	}
	jp := (jpll[face]*nr + int64(ix) - int64(iy) + 1 + kshift) / 2
	if jp < 1 {
		jp += nl4
	}
	return nBefore + jp - 1
}

// spreadBits moves bit k of v to bit 2k.
func spreadBits(v int) int64 {
	x := uint64(uint32(v))
	x = (x | x<<16) & 0x0000ffff0000ffff
	x = (x | x<<8) & 0x00ff00ff00ff00ff
	x = (x | x<<4) & 0x0f0f0f0f0f0f0f0f
	x = (x | x<<2) & 0x3333333333333333
	x = (x | x<<1) & 0x5555555555555555
	return int64(x)
}

// compressBits is the inverse of spreadBits (keeps the even bits).
func compressBits(v int64) int {
	x := uint64(v) & 0x5555555555555555
	x = (x | x>>1) & 0x3333333333333333
	x = (x | x>>2) & 0x0f0f0f0f0f0f0f0f
	x = (x | x>>4) & 0x00ff00ff00ff00ff
	x = (x | x>>8) & 0x0000ffff0000ffff
	x = (x | x>>16) & 0x00000000ffffffff
	return int(x)
}

func (g geom) nest2xyf(pix int64) (ix, iy, face int) {
	face = int(pix >> (2 * g.order))
	pix &= g.npface - 1
	return compressBits(pix), compressBits(pix >> 1), face
}

func (g geom) xyf2nest(ix, iy, face int) int64 {
	return int64(face)<<(2*g.order) + spreadBits(ix) + spreadBits(iy)<<1
}

func (g geom) pix2xyf(pix int64, scheme Scheme) (ix, iy, face int) {
	if scheme == Nested {
		return g.nest2xyf(pix)
	}
	return g.ring2xyf(pix)
}

func (g geom) xyf2pix(ix, iy, face int, scheme Scheme) int64 {
	if scheme == Nested {
		return g.xyf2nest(ix, iy, face)
	}
	return g.xyf2ring(ix, iy, face)
}

// pix2xyf returns the (x, y, face) decomposition of pix, matching
// healpy.pix2xyf.
func pix2xyf(nside int, pix int64, scheme Scheme) (x, y, face int) {
	g := newGeom("pix2xyf", nside, scheme)
	g.checkPix("pix2xyf", pix)
	return g.pix2xyf(pix, scheme)
}

// xyf2pix is the inverse of pix2xyf, matching healpy.xyf2pix.
func xyf2pix(nside, x, y, face int, scheme Scheme) int64 {
	g := newGeom("xyf2pix", nside, scheme)
	if face < 0 || face > 11 || x < 0 || y < 0 || int64(x) >= g.nside || int64(y) >= g.nside {
		panic(fmt.Sprintf("healpix.xyf2pix: (x=%d, y=%d, face=%d) out of range for nside %d", x, y, face, nside))
	}
	return g.xyf2pix(x, y, face, scheme)
}
