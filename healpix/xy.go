package healpix

// Internal XY scheme for HEALPix pixel indexing.
//
// Every pixel is identified by a base healpix face (0-11) plus an
// (x, y) coordinate within that face, where both x and y are in
// [0, nside-1]. The flat XY index is:
//
//	xyIndex = face * nside² + x * nside + y
//
// This representation is the computational core — both RING and
// NESTED ordering are derived from it, and ang2pix/pix2ang operate
// in XY space internally.
//
// # Face numbering (Gorski et al. 2005)
//
//	0-3:  north polar cap (faces touching the north pole)
//	4-7:  equatorial band
//	8-11: south polar cap (faces touching the south pole)
//
// Within each face, x increases toward the east (NE on the sphere) and
// y increases toward the north (NW on the sphere).
//
// # Reference
//
// Port of astrometry.net healpix.c struct hp_t and the
// healpix_compose_xy / healpix_decompose_xy functions.

// hpXY is the internal (face, x, y) representation of a pixel.
type hpXY struct {
	Face int // 0-11, the base HEALPix face
	X    int // [0, nside-1], increasing NE
	Y    int // [0, nside-1], increasing NW
}

// xyCompose converts (face, x, y) to a flat XY-scheme pixel index.
func xyCompose(face, x, y, nside int) int64 {
	ns := int64(nside)
	return int64(face)*ns*ns + int64(x)*ns + int64(y)
}

// xyDecompose converts a flat XY-scheme pixel index to (face, x, y).
func xyDecompose(hp int64, nside int) hpXY {
	ns := int64(nside)
	ns2 := ns * ns
	face := int(hp / ns2)
	rem := hp % ns2
	x := int(rem / ns)
	y := int(rem % ns)
	return hpXY{Face: face, X: x, Y: y}
}
