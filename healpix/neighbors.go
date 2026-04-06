package healpix

// GetNeighbors returns the pixel indices of up to 8 neighbors of
// pixel pix. The returned slice contains 6-8 valid pixel indices;
// some polar pixels have fewer than 8 neighbors because the quad-tree
// structure has 3-valent vertices at the poles.
//
// The neighbor order matches healpy's get_all_neighbours convention
// (set-equal comparison; order within the slice is implementation-
// defined). -1 entries indicate absent neighbors.
//
// Faithful port of astrometry.net get_neighbours + healpix_get_neighbour
// (healpix.c:498-745).
func GetNeighbors(nside int, pix int64, scheme Scheme) [8]int64 {
	var xy int64
	switch scheme {
	case Nested:
		xy = nestedToXY(pix, nside)
	default:
		xy = ringToXY(pix, nside)
	}
	hp := xyDecompose(xy, nside)
	ns := nside
	base := hp.Face
	x := hp.X
	y := hp.Y

	var result [8]int64
	for i := range result {
		result[i] = -1
	}
	nn := 0

	emit := func(nbase, nx, ny int) {
		if nbase < 0 || nn >= 8 {
			return
		}
		nxy := xyCompose(nbase, nx, ny, nside)
		switch scheme {
		case Nested:
			result[nn] = xyToNested(nxy, nside)
		default:
			result[nn] = xyToRing(nxy, nside)
		}
		nn++
	}

	var nbase, nx, ny int

	// (+1, 0) — E direction
	nx = (x + 1) % ns
	ny = y
	if x == ns-1 {
		nbase = baseNeighbor(base, 1, 0)
		if isNorthPolar(base) {
			nx = x
			nx, ny = ny, nx
		}
	} else {
		nbase = base
	}
	emit(nbase, nx, ny)

	// (+1, +1) — NE direction
	nx = (x + 1) % ns
	ny = (y + 1) % ns
	if x == ns-1 && y == ns-1 {
		if isPolar(base) {
			nbase = baseNeighbor(base, 1, 1)
		} else {
			nbase = -1
		}
	} else if x == ns-1 {
		nbase = baseNeighbor(base, 1, 0)
	} else if y == ns-1 {
		nbase = baseNeighbor(base, 0, 1)
	} else {
		nbase = base
	}
	if isNorthPolar(base) {
		if x == ns-1 {
			nx = ns - 1
		}
		if y == ns-1 {
			ny = ns - 1
		}
		if x == ns-1 || y == ns-1 {
			nx, ny = ny, nx
		}
	}
	if nbase >= 0 {
		emit(nbase, nx, ny)
	}

	// (0, +1) — N direction
	nx = x
	ny = (y + 1) % ns
	if y == ns-1 {
		nbase = baseNeighbor(base, 0, 1)
		if isNorthPolar(base) {
			ny = y
			nx, ny = ny, nx
		}
	} else {
		nbase = base
	}
	emit(nbase, nx, ny)

	// (-1, +1) — NW direction
	nx = (x + ns - 1) % ns
	ny = (y + 1) % ns
	if x == 0 && y == ns-1 {
		if isEquatorial(base) {
			nbase = baseNeighbor(base, -1, 1)
		} else {
			nbase = -1
		}
	} else if x == 0 {
		nbase = baseNeighbor(base, -1, 0)
		if isSouthPolar(base) {
			nx = 0
			nx, ny = ny, nx
		}
	} else if y == ns-1 {
		nbase = baseNeighbor(base, 0, 1)
		if isNorthPolar(base) {
			ny = y
			nx, ny = ny, nx
		}
	} else {
		nbase = base
	}
	if nbase >= 0 {
		emit(nbase, nx, ny)
	}

	// (-1, 0) — W direction
	nx = (x + ns - 1) % ns
	ny = y
	if x == 0 {
		nbase = baseNeighbor(base, -1, 0)
		if isSouthPolar(base) {
			nx = 0
			nx, ny = ny, nx
		}
	} else {
		nbase = base
	}
	emit(nbase, nx, ny)

	// (-1, -1) — SW direction
	nx = (x + ns - 1) % ns
	ny = (y + ns - 1) % ns
	if x == 0 && y == 0 {
		if isPolar(base) {
			nbase = baseNeighbor(base, -1, -1)
		} else {
			nbase = -1
		}
	} else if x == 0 {
		nbase = baseNeighbor(base, -1, 0)
	} else if y == 0 {
		nbase = baseNeighbor(base, 0, -1)
	} else {
		nbase = base
	}
	if isSouthPolar(base) {
		if x == 0 {
			nx = 0
		}
		if y == 0 {
			ny = 0
		}
		if x == 0 || y == 0 {
			nx, ny = ny, nx
		}
	}
	if nbase >= 0 {
		emit(nbase, nx, ny)
	}

	// (0, -1) — S direction
	ny = (y + ns - 1) % ns
	nx = x
	if y == 0 {
		nbase = baseNeighbor(base, 0, -1)
		if isSouthPolar(base) {
			ny = y
			nx, ny = ny, nx
		}
	} else {
		nbase = base
	}
	emit(nbase, nx, ny)

	// (+1, -1) — SE direction
	nx = (x + 1) % ns
	ny = (y + ns - 1) % ns
	if x == ns-1 && y == 0 {
		if isEquatorial(base) {
			nbase = baseNeighbor(base, 1, -1)
		} else {
			nbase = -1
		}
	} else if x == ns-1 {
		nbase = baseNeighbor(base, 1, 0)
		if isNorthPolar(base) {
			nx = x
			nx, ny = ny, nx
		}
	} else if y == 0 {
		nbase = baseNeighbor(base, 0, -1)
		if isSouthPolar(base) {
			ny = y
			nx, ny = ny, nx
		}
	} else {
		nbase = base
	}
	if nbase >= 0 {
		emit(nbase, nx, ny)
	}

	return result
}

// baseNeighbor returns the neighboring base healpix face when moving
// from face hp in direction (dx, dy). Returns -1 if no neighbor exists
// in that direction. Port of astrometry.net healpix_get_neighbour
// (healpix.c:498).
func baseNeighbor(hp, dx, dy int) int {
	if isNorthPolar(hp) {
		switch {
		case dx == 1 && dy == 0:
			return (hp + 1) % 4
		case dx == 0 && dy == 1:
			return (hp + 3) % 4
		case dx == 1 && dy == 1:
			return (hp + 2) % 4
		case dx == -1 && dy == 0:
			return hp + 4
		case dx == 0 && dy == -1:
			return 4 + (hp+1)%4
		case dx == -1 && dy == -1:
			return hp + 8
		}
		return -1
	}
	if isSouthPolar(hp) {
		sp := hp - 8
		switch {
		case dx == 1 && dy == 0:
			return 4 + (sp+1)%4
		case dx == 0 && dy == 1:
			return hp - 4
		case dx == -1 && dy == 0:
			return 8 + (sp+3)%4
		case dx == 0 && dy == -1:
			return 8 + (sp+1)%4
		case dx == -1 && dy == -1:
			return 8 + (sp+2)%4
		case dx == 1 && dy == 1:
			return hp - 8
		}
		return -1
	}
	// Equatorial (4-7).
	eq := hp - 4
	switch {
	case dx == 1 && dy == 0:
		return hp - 4
	case dx == 0 && dy == 1:
		return (eq + 3) % 4
	case dx == -1 && dy == 0:
		return 8 + (eq+3)%4
	case dx == 0 && dy == -1:
		return hp + 4
	case dx == 1 && dy == -1:
		return 4 + (eq+1)%4
	case dx == -1 && dy == 1:
		return 4 + ((eq - 1 + 4) % 4)
	}
	return -1
}

func isNorthPolar(hp int) bool  { return hp <= 3 }
func isSouthPolar(hp int) bool  { return hp >= 8 }
func isEquatorial(hp int) bool  { return hp >= 4 && hp <= 7 }
func isPolar(hp int) bool       { return hp <= 3 || hp >= 8 }
