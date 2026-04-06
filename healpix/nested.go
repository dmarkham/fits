package healpix

// NESTED ↔ XY conversion via bit-interleaving.
//
// The NESTED scheme numbers pixels within each base face by
// interleaving the bits of the x and y coordinates:
//
//	nested_subindex bit 2k   = x bit k
//	nested_subindex bit 2k+1 = y bit k
//
// Only valid for Nside that is a power of 2 (required for the
// hierarchical quad-tree to have integer depth).
//
// Port of astrometry.net healpix_xy_to_nested / healpix_nested_to_xy.

// xyToNested converts an XY-scheme pixel index to NESTED ordering.
// Panics if nside is not a power of 2.
func xyToNested(xy int64, nside int) int64 {
	if !IsNsidePow2(nside) {
		panic("healpix: xyToNested requires power-of-2 nside")
	}
	hp := xyDecompose(xy, nside)
	return int64(hp.Face)*int64(nside)*int64(nside) + interleave(hp.X, hp.Y)
}

// nestedToXY converts a NESTED pixel index to XY-scheme ordering.
// Panics if nside is not a power of 2.
func nestedToXY(nested int64, nside int) int64 {
	ns2 := int64(nside) * int64(nside)
	face := int(nested / ns2)
	subidx := nested % ns2
	x, y := deinterleave(subidx)
	return xyCompose(face, x, y, nside)
}

// interleave interleaves the bits of x and y into a single index:
// result bit 2k = x bit k, result bit 2k+1 = y bit k.
func interleave(x, y int) int64 {
	var result int64
	for i := 0; x > 0 || y > 0; i++ {
		result |= int64(x&1) << (2 * i)
		result |= int64(y&1) << (2*i + 1)
		x >>= 1
		y >>= 1
	}
	return result
}

// deinterleave extracts x and y from an interleaved index.
func deinterleave(idx int64) (x, y int) {
	for i := 0; idx > 0; i++ {
		x |= int(idx&1) << i
		idx >>= 1
		y |= int(idx&1) << i
		idx >>= 1
	}
	return x, y
}
