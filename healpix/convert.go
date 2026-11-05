package healpix

// Ring2Nest converts a RING-ordered pixel index to NESTED ordering.
// Panics unless nside is a power of 2 in [1, 2^29] and ring is a valid
// pixel index. Agrees exactly with healpy.ring2nest.
func Ring2Nest(nside int, ring int64) int64 {
	g := newGeom("Ring2Nest", nside, Nested)
	g.checkPix("Ring2Nest", ring)
	ix, iy, f := g.ring2xyf(ring)
	return g.xyf2nest(ix, iy, f)
}

// Nest2Ring converts a NESTED pixel index to RING ordering. Panics like
// Ring2Nest. Agrees exactly with healpy.nest2ring.
func Nest2Ring(nside int, nested int64) int64 {
	g := newGeom("Nest2Ring", nside, Nested)
	g.checkPix("Nest2Ring", nested)
	ix, iy, f := g.nest2xyf(nested)
	return g.xyf2ring(ix, iy, f)
}
