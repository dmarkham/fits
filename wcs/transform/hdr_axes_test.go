package transform

// Pixel-order adapter between hdr2sky goldens and the Go API.
//
// hdr2sky drives pixel axis (lng+1) with the case's "x" and pixel axis
// (lat+1) with "y" (pixel axis i carries world axis i), so for a header
// whose longitude axis comes AFTER its latitude axis (CTYPE1 = DEC,
// CTYPE2 = RA) golden "x" is pixel axis 2. Transform.PixelToSky(p1, p2)
// takes the two celestial pixel axes in increasing axis order, so the
// pair is swapped when Lng > Lat. For every other case this is the
// identity (minus the 1-based to 0-based shift).

func hdrPixelToSky(g *hdrGolden, tr *Transform, x, y float64) (float64, float64, error) {
	if g.WCS.Lng > g.WCS.Lat {
		x, y = y, x
	}
	return tr.PixelToSky(x-1, y-1)
}

func hdrSkyToPixel(g *hdrGolden, tr *Transform, lng, lat float64) (float64, float64, error) {
	p1, p2, err := tr.SkyToPixel(lng, lat)
	if g.WCS.Lng > g.WCS.Lat {
		p1, p2 = p2, p1
	}
	return p1 + 1, p2 + 1, err
}
