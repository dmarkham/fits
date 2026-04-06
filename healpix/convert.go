package healpix

// Ring2Nest converts a RING-ordered pixel index to NESTED ordering.
// Requires nside to be a power of 2; panics otherwise.
func Ring2Nest(nside int, ring int64) int64 {
	xy := ringToXY(ring, nside)
	return xyToNested(xy, nside)
}

// Nest2Ring converts a NESTED pixel index to RING ordering.
// Requires nside to be a power of 2; panics otherwise.
func Nest2Ring(nside int, nested int64) int64 {
	xy := nestedToXY(nested, nside)
	return xyToRing(xy, nside)
}
