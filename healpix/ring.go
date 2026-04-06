package healpix

// RING ↔ XY conversion.
//
// Faithful port of astrometry.net healpix_xy_to_ring (healpix.c:262)
// and healpix_ring_to_xy (healpix.c:163). The ring scheme numbers
// pixels along isolatitude rings from the north pole (ring 1) to the
// south pole (ring 4·nside-1). Three zones:
//
//	North polar cap:  rings 1..nside            (4·ring pixels per ring)
//	Equatorial belt:  rings nside+1..3·nside-1  (4·nside pixels per ring)
//	South polar cap:  rings 3·nside..4·nside-1  (mirrored north)

// decomposeRing finds the 1-based ring number and the longitude index
// within that ring for a RING-ordered pixel index. Port of
// astrometry.net healpix_decompose_ring (healpix.c:126).
func decomposeRing(ringIdx int64, nside int) (ring, longind int64) {
	ns := int64(nside)
	offset := int64(0)
	ring = 1
	// North polar cap.
	for ring <= ns {
		npixRing := ring * 4
		if offset+npixRing > ringIdx {
			return ring, ringIdx - offset
		}
		offset += npixRing
		ring++
	}
	// Equatorial belt.
	for ring < 3*ns {
		npixRing := ns * 4
		if offset+npixRing > ringIdx {
			return ring, ringIdx - offset
		}
		offset += npixRing
		ring++
	}
	// South polar cap.
	for ring < 4*ns {
		npixRing := (4*ns - ring) * 4
		if offset+npixRing > ringIdx {
			return ring, ringIdx - offset
		}
		offset += npixRing
		ring++
	}
	return 0, 0 // invalid input
}

// xyToRing converts an XY-scheme pixel index to RING ordering.
// Port of astrometry.net healpix_xy_to_ring (healpix.c:262).
func xyToRing(xy int64, nside int) int64 {
	hp := xyDecompose(xy, nside)
	ns := int64(nside)
	face := hp.Face
	x := int64(hp.X)
	y := int64(hp.Y)

	frow := int64(face / 4)
	F1 := frow + 2
	v := x + y
	ring := F1*ns - v - 1 // 1-based ring number

	var index int64
	if ring <= ns {
		// North polar.
		index = (ns - 1 - y) + int64(face%4)*ring
		index += ring * (ring - 1) * 2
	} else if ring >= 3*ns {
		// South polar. Flip so labeling starts in SE corner, then
		// subtract from total.
		ri := 4*ns - ring
		index = (ri - 1) - x
		index += int64(3-face%4) * ri
		index += ri * (ri - 1) * 2
		index = 12*ns*ns - 1 - index
	} else {
		// Equatorial.
		s := (ring - ns) % 2
		F2 := 2*int64(face%4) - (frow % 2) + 1
		h := x - y
		index = (F2*ns + h + s) / 2
		index += ns * (ns - 1) * 2
		index += ns * 4 * (ring - ns)
		// Handle healpix #4 wrap-around.
		if face == 4 && y > x {
			index += 4*ns - 1
		}
	}
	return index
}

// ringToXY converts a RING pixel index to XY-scheme ordering.
// Port of astrometry.net healpix_ring_to_xy (healpix.c:163), which
// calls healpix_decompose_ring internally.
func ringToXY(ringIdx int64, nside int) int64 {
	ns := int64(nside)

	// Step 1: decompose ringIdx into (ring number, longitude index).
	ring, longind := decomposeRing(ringIdx, nside)

	// Step 2: convert (ring, longind) to (face, x, y).
	var face, x, y int

	if ring <= ns {
		// North polar. Port of healpix.c:167-179.
		face = int(longind / ring)
		ind := longind - int64(face)*ring
		y = nside - 1 - int(ind)
		frow := face / 4
		F1 := frow + 2
		v := F1*nside - int(ring) - 1
		x = v - y
	} else if ring >= 3*ns {
		// South polar. Port of healpix.c:244-258.
		ri := 4*ns - ring
		face = 8 + int(longind/ri)
		ind := longind - int64(face%4)*ri
		y = int(ri-1) - int(ind)
		frow := face / 4
		F1 := frow + 2
		v := F1*nside - int(ring) - 1
		x = v - y
	} else {
		// Equatorial. Port of healpix.c:180-243.
		panel := int(longind / ns)
		ind := int(longind % ns)
		bottomleft := ind < (int(ring-ns)+1)/2
		topleft := ind < (3*nside-int(ring)+1)/2

		R := 0
		switch {
		case !bottomleft && topleft:
			// Top row (north polar faces).
			face = panel
		case bottomleft && !topleft:
			// Bottom row (south polar faces).
			face = 8 + panel
		case bottomleft && topleft:
			// Left side (equatorial faces).
			face = 4 + panel
		default: // !bottomleft && !topleft
			// Right side.
			face = 4 + (panel+1)%4
			if face == 4 {
				longind -= 4*ns - 1
				R = 1
			}
		}

		frow := face / 4
		F1 := frow + 2
		F2 := 2*(face%4) - (frow % 2) + 1
		s := (int(ring) - nside) % 2
		v := F1*nside - int(ring) - 1
		h := 2*int(longind) - s - F2*nside
		if R != 0 {
			h--
		}
		x = (v + h) / 2
		y = (v - h) / 2

		// Fixup for integer-division rounding (healpix.c:230-240).
		if v != x+y || h != x-y {
			h++
			x = (v + h) / 2
			y = (v - h) / 2
		}
	}

	return xyCompose(face, x, y, nside)
}
