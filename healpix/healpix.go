// Package healpix implements HEALPix pixel indexing for the
// Hierarchical Equal-Area isoLatitude Pixelisation of the sphere
// (Gorski et al. 2005).
//
// It provides the pixel <-> direction operations (Ang2Pix, Vec2Pix,
// Pix2Ang, Pix2Vec), RING <-> NESTED conversion, the 8-neighbour lookup and
// the per-Nside utilities. It does not read or write map files itself: a
// HEALPix map in FITS is an ordinary binary table (PIXTYPE='HEALPIX',
// ORDERING, NSIDE in the header), read with the fits package's table API
// and indexed with this package (see TestRealMapHeader_WMAP and
// TestHealpyWrittenMap for both steps checked end to end).
//
// RING ordering works for any Nside in [1, 2^29]; NESTED requires a power
// of 2. Pixel indices are int64 (12·4^29 ≈ 3.46e18 pixels at Nside 2^29).
// Invalid input (bad Nside, NESTED with a non-power-of-2 Nside, a pixel
// index out of range, theta outside [0, π], non-finite angles, the zero
// vector) panics instead of returning a wrong pixel.
//
// # Coordinate convention
//
// Angles follow the physics / healpy convention:
//
//   - theta (colatitude): 0 at the north pole, π at the south pole
//   - phi (longitude): 0 to 2π, increasing eastward
//
// # Validation
//
// Every exported function is checked against healpy 1.20.1 (which wraps
// healpix_cxx) on the committed goldens in healpix/testdata/golden, made by
// testdata/gen_golden.py (provenance in testdata/PROVENANCE.md). The
// per-function goldens (pix2ang, pix2vec, ang2pix, neighbours) hold the
// power-of-2 Nside 1, 2, 4, 8, 16, 64, 256, 1024 and 8192 (not 32, 128,
// 512, 2048 or 4096); the pixel-function goldens add 2^20 and 2^29 in both
// orderings and the non-power-of-2 RING Nside 3, 5, 6, 7, 12, 100, 1000,
// 3^10, 3·2^27 and 2^29-1:
//
//   - Integer results (Ang2Pix, Vec2Pix, Ring2Nest, Nest2Ring,
//     GetNeighbors including slot order, Npix2Nside, IsNsidePow2/NsideOK)
//     agree exactly, with the documented exceptions where healpy is wrong
//     (Npix2Nside; Vec2Pix of vectors whose squared components overflow
//     or lose bits in the subnormal range, which Vec2Pix rescales) and,
//     for Ang2Pix and Vec2Pix, boundary-ambiguous directions (next item).
//   - Ang2Pix and Vec2Pix call Go's math.Cos, math.Sin and math.Atan2,
//     which can differ from glibc's (what healpy calls) in the last bit.
//     For a direction within about one ULP of a pixel boundary that bit
//     decides the pixel, and such a direction has no meaningful correct
//     side. So the pixel equals healpy's everywhere except at
//     boundary-ambiguous directions, where Ang2Pix's own answer changes
//     under a nudge of a few ULP of its input or of the cos/sin/atan2
//     values. The tests on directions placed on and within 2 ULP of
//     polar-cap boundaries, the cos/sin switch-over and z = 2/3 at Nside up
//     to 2^29 (TestAng2PixPolar_Healpy, TestBranchSwitch_Healpy,
//     TestAng2PixTwoThird_Healpy) accept a different pixel only at such a
//     direction and count them (2.7%, 12% and 2.4% of those constructed
//     inputs); the random and standard direction sets agree exactly. The
//     arithmetic after the libm calls, including every branch choice, is
//     checked exactly: fed the correctly rounded cos/sin/atan2 the
//     expected values were made with, it reproduces them bit for bit.
//   - PixelArea and PixelResol agree bit for bit.
//   - Pix2Ang theta and phi agree to within 1 ULP (the measured maximum,
//     enforced by TestAngleULPBudget_Measured); the formulas and their
//     operation order are healpy's, so the difference is libm rounding.
//   - MaxPixRad is within 2 ULP of the exact value of healpix_cxx's
//     definition (60-digit mpmath) and within 2^-52 rad of healpy, whose
//     own formula cancels (see MaxPixRad).
//   - Pix2Vec components agree to within 2·2^-52 absolute.
//
// # Reference
//
// The pixel arithmetic follows Gorski et al. 2005 (ApJ 622, 759) and the
// algorithms of healpix_cxx (healpix_base.cc, as used by healpy): O(1)
// RING decomposition by integer square root, the z/phi formulas of
// pix2loc/loc2pix including the sin(theta) branch near the poles, the
// face-neighbour tables, and max_pixrad.
package healpix

import "math"

// Scheme identifies the pixel ordering convention.
type Scheme int

const (
	Ring   Scheme = 0 // RING ordering: pixels numbered along isolatitude rings
	Nested Scheme = 1 // NESTED ordering: pixels numbered by hierarchical quad-tree
)

// String returns "RING" or "NESTED".
func (s Scheme) String() string {
	if s == Nested {
		return "NESTED"
	}
	return "RING"
}

// Nside2Npix returns the total number of pixels for the given Nside
// parameter: npix = 12 * nside². The result is meaningless for nside <= 0
// or nside > 2^29 (use NsideOK to check).
func Nside2Npix(nside int) int64 {
	return 12 * int64(nside) * int64(nside)
}

// Npix2Nside returns the Nside parameter corresponding to npix total
// pixels. Returns 0 and false if npix is not a valid HEALPix pixel
// count (i.e. not of the form 12·nside² with 1 <= nside <= 2^29). It
// agrees with healpy.isnpixok / npix2nside except where healpy is wrong
// or has no limit: healpy accepts npix 0, accepts 12·2^58+12 (it rounds
// npix through a double), and accepts 12·nside² for nside > 2^29 (e.g.
// 12·(2^29+1)², which healpix_cxx cannot index); all three are rejected
// here.
func Npix2Nside(npix int64) (int, bool) {
	if npix <= 0 || npix%12 != 0 {
		return 0, false
	}
	ns2 := npix / 12
	nside := isqrt(ns2)
	if nside*nside != ns2 || nside > MaxNside {
		return 0, false
	}
	return int(nside), true
}

// IsNsidePow2 reports whether nside is a positive power of 2 (required
// for NESTED ordering). It does not check the 2^29 upper limit; NsideOK
// does.
func IsNsidePow2(nside int) bool {
	return nside > 0 && nside&(nside-1) == 0
}

// NsideOK reports whether nside is usable with scheme: 1 <= nside <= 2^29,
// and a power of 2 for NESTED. Equivalent to healpy.isnsideok(nside,
// nest=(scheme == Nested)).
func NsideOK(nside int, scheme Scheme) bool {
	if nside <= 0 || nside > MaxNside {
		return false
	}
	return scheme != Nested || IsNsidePow2(nside)
}

// PixelArea returns the solid angle of one HEALPix pixel in steradians.
// All pixels have equal area by construction: area = 4π / (12·nside²),
// bit-identical to healpy.nside2pixarea.
func PixelArea(nside int) float64 {
	return 4 * math.Pi / float64(Nside2Npix(nside))
}

// PixelResol returns the pixel "resolution" — the square root of the
// pixel area, in radians. Bit-identical to healpy.nside2resol.
func PixelResol(nside int) float64 {
	return math.Sqrt(PixelArea(nside))
}

// MaxPixRad returns the maximum angular distance (radians) between any
// pixel centre and its corners at the given Nside — healpy.max_pixrad.
//
// The quantity is the healpix_cxx definition: the angle between the point
// A (z=2/3, phi=π/(4·nside)) and the point B (z=1-(1-1/nside)²/3, phi=0).
// healpix_cxx evaluates it as atan2(|A×B|, A·B) from rounded Cartesian
// components, which cancels badly for large Nside (A and B nearly
// coincide). MaxPixRad instead forms the chord |A-B| from the coordinate
// differences directly, so it is within 2 ULP of the exact value
// (checked against a 60-digit mpmath evaluation of the same definition;
// 2 ULP is measured at nside 16384),
// and agrees with healpy to within 2^-52 rad absolute — the accuracy of
// healpy's own cancelling formula. Panics on an invalid nside.
func MaxPixRad(nside int) float64 {
	newGeom("MaxPixRad", nside, Ring)
	n := float64(nside)
	phi := pi / (4 * n)
	q := (n - 1) / n
	t1 := float64(q * q)                          // (1-1/n)^2; 1 - zB = t1/3
	dz := (float64(2*n) - 1) / (3 * n * n)        // zB - zA = (1 - t1)/3, exact form
	sa := math.Sqrt(5) / 3                        // sin(theta_A), zA = 2/3
	sb := math.Sqrt(float64(t1 / 3 * (2 - t1/3))) // sin(theta_B) = sqrt((1-zB)(1+zB))
	sum := (5 - t1) / 3                           // zB + zA
	dsab := dz * sum / (sa + sb)                  // sa - sb = (zB^2 - zA^2)/(sa + sb)
	sh := math.Sin(phi / 2)
	dx := dsab - float64(2*sa*sh*sh) // xA - xB = (sa - sb) - sa(1 - cos phi)
	dy := float64(sa * math.Sin(phi))
	chord := math.Sqrt(float64(dx*dx) + float64(dy*dy) + float64(dz*dz))
	return float64(2 * math.Asin(chord/2))
}
