// Package healpix implements HEALPix pixel indexing for the
// Hierarchical Equal-Area isoLatitude Pixelisation of the sphere
// (Gorski et al. 2005).
//
// This package provides the core pixel↔direction operations needed to
// read, write, and query HEALPix maps stored in FITS binary tables
// (the standard format for Planck, WMAP, Gaia, and other sky surveys).
//
// All functions support both RING and NESTED pixel ordering and work
// for any valid Nside value (NESTED requires power-of-2 Nside).
//
// # Coordinate convention
//
// Angles follow the physics / healpy convention:
//
//   - theta (colatitude): 0 at the north pole, π at the south pole
//   - phi (longitude): 0 to 2π, increasing eastward
//
// Pixel indices are int64 to support Nside up to 2^29 (the HEALPix
// practical maximum, giving ~3.4 trillion pixels).
//
// # Validation
//
// Every function is cross-validated against healpy 1.19.0 on a
// comprehensive golden fixture set (healpix/testdata/golden/). Integer
// pixel indices must agree exactly; floating-point angles must agree
// to machine precision. See healpix_test.go for the full test matrix.
//
// # Reference
//
// The implementation is a clean-room port of the algorithms described
// in:
//
//   - Gorski et al. 2005 (ApJ 622, 759) — HEALPix paper
//   - astrometry.net healpix.c — C reference implementation
//   - healpix_cxx (Siril subprojects) — C++ reference
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
// parameter: npix = 12 * nside². Behavior is undefined for nside <= 0.
func Nside2Npix(nside int) int64 {
	return 12 * int64(nside) * int64(nside)
}

// Npix2Nside returns the Nside parameter corresponding to npix total
// pixels. Returns 0 and false if npix is not a valid HEALPix pixel
// count (i.e. not of the form 12 * nside²).
func Npix2Nside(npix int64) (int, bool) {
	if npix <= 0 || npix%12 != 0 {
		return 0, false
	}
	ns2 := npix / 12
	nside := int(math.Sqrt(float64(ns2)) + 0.5)
	if int64(nside)*int64(nside) != ns2 {
		return 0, false
	}
	return nside, true
}

// IsNsidePow2 reports whether nside is a positive power of 2 (required
// for NESTED ordering).
func IsNsidePow2(nside int) bool {
	return nside > 0 && nside&(nside-1) == 0
}

// PixelArea returns the solid angle of one HEALPix pixel in steradians.
// All pixels have equal area by construction: area = 4π / (12·nside²).
func PixelArea(nside int) float64 {
	return 4 * math.Pi / float64(Nside2Npix(nside))
}

// PixelResol returns the pixel "resolution" — the square root of the
// pixel area, in radians. Matches healpy.nside2resol.
func PixelResol(nside int) float64 {
	return math.Sqrt(PixelArea(nside))
}

// MaxPixRad returns the maximum angular radius (in radians) of any
// HEALPix pixel at the given Nside. This is the radius of the smallest
// circle that encloses the largest pixel. The worst case is always the
// polar cap pixel (ring 1), whose farthest corner from its center is
// the diagonal toward (theta=0, phi=±π/2) or the equatorward boundary.
//
// The computation checks three candidate distances for the ring-1 pixel
// and returns the maximum: center-to-pole, center-to-equatorward-edge,
// and center-to-diagonal-corner. Validated against healpy.max_pixrad.
func MaxPixRad(nside int) float64 {
	ns := float64(nside)
	t1 := 1.0 - 1.0/(3.0*ns*ns)   // cos(theta) of ring-1 center
	t2 := 1.0 - 2.0/(3.0*ns*ns)   // cos(theta) of ring-1/2 boundary
	centerTheta := math.Acos(t1)
	boundaryTheta := math.Acos(t2)
	toPole := centerTheta
	toEquator := boundaryTheta - centerTheta
	// Diagonal corner: ring-1 pixel spans phi=[0, π/2], so the far
	// corner is at (boundaryTheta, π/2). Center is at (centerTheta, π/4).
	cornerDist := math.Acos(
		math.Sin(centerTheta)*math.Sin(boundaryTheta)*math.Cos(math.Pi/4) +
			math.Cos(centerTheta)*math.Cos(boundaryTheta))
	r := toPole
	if toEquator > r {
		r = toEquator
	}
	if cornerDist > r {
		r = cornerDist
	}
	return r
}
