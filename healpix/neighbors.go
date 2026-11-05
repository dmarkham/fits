package healpix

// Neighbour offsets in (x, y) within a face, in output order
// SW, W, NW, N, NE, E, SE, S.
var (
	nbXOffset = [8]int{-1, -1, 0, 1, 1, 1, 0, -1}
	nbYOffset = [8]int{0, 1, 1, 1, 0, -1, -1, -1}
)

// nbFaceArray[d][f] is the face reached from face f when stepping off it
// in direction d (index 3*dy + dx + 4 for dx, dy in {-1, 0, 1}); -1 where
// no face exists (the 3-valent corners at the poles and on the equator).
var nbFaceArray = [9][12]int{
	{8, 9, 10, 11, -1, -1, -1, -1, 10, 11, 8, 9}, // S
	{5, 6, 7, 4, 8, 9, 10, 11, 9, 10, 11, 8},     // SE
	{-1, -1, -1, -1, 5, 6, 7, 4, -1, -1, -1, -1}, // E
	{4, 5, 6, 7, 11, 8, 9, 10, 11, 8, 9, 10},     // SW
	{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},       // centre
	{1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4},         // NE
	{-1, -1, -1, -1, 7, 4, 5, 6, -1, -1, -1, -1}, // W
	{3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7},         // NW
	{2, 3, 0, 1, -1, -1, -1, -1, 0, 1, 2, 3},     // N
}

// nbSwapArray[d][f/4]: bit 1 mirrors x, bit 2 mirrors y, bit 4 swaps x
// and y when entering the neighbouring face.
var nbSwapArray = [9][3]int{
	{0, 0, 3}, // S
	{0, 0, 6}, // SE
	{0, 0, 0}, // E
	{0, 0, 5}, // SW
	{0, 0, 0}, // centre
	{5, 0, 0}, // NE
	{0, 0, 0}, // W
	{6, 0, 0}, // NW
	{3, 0, 0}, // N
}

// GetNeighbors returns the 8 neighbours of pixel pix in the order
// SW, W, NW, N, NE, E, SE, S — the order of healpy.get_all_neighbours
// (and healpix_cxx neighbors). A slot is -1 where that neighbour does not
// exist: for nside >= 2 exactly 24 pixels — the three around each of the
// eight 3-valent vertices on the |z| = 2/3 boundary — lack one of W, N, E
// or S; for nside 1 every pixel lacks two. The result equals healpy's slot
// for slot, including the position of -1.
//
// Panics on an invalid nside/scheme or a pixel outside [0, 12·nside²).
func GetNeighbors(nside int, pix int64, scheme Scheme) [8]int64 {
	g := newGeom("GetNeighbors", nside, scheme)
	g.checkPix("GetNeighbors", pix)
	ix, iy, face := g.pix2xyf(pix, scheme)
	ns := int(g.nside)
	var result [8]int64
	for i := 0; i < 8; i++ {
		x, y := ix+nbXOffset[i], iy+nbYOffset[i]
		nbnum := 4
		if x < 0 {
			x += ns
			nbnum--
		} else if x >= ns {
			x -= ns
			nbnum++
		}
		if y < 0 {
			y += ns
			nbnum -= 3
		} else if y >= ns {
			y -= ns
			nbnum += 3
		}
		f := nbFaceArray[nbnum][face]
		if f < 0 {
			result[i] = -1
			continue
		}
		b := nbSwapArray[nbnum][face>>2]
		if b&1 != 0 {
			x = ns - x - 1
		}
		if b&2 != 0 {
			y = ns - y - 1
		}
		if b&4 != 0 {
			x, y = y, x
		}
		result[i] = g.xyf2pix(x, y, f, scheme)
	}
	return result
}
