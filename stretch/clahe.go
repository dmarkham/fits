package stretch

// CLAHE — Contrast Limited Adaptive Histogram Equalization.
//
// Implements the Zuiderveld 1994 algorithm. For RGB images, CLAHE is
// applied to the luminance channel only (after converting to a
// luminance + chrominance representation), preserving color.
//
// This is a from-scratch implementation because Siril delegates CLAHE
// entirely to OpenCV's createCLAHE(). We validate against Siril's
// output on our golden fixtures.
//
// Parameters:
//   - tileGridSize: number of tiles along each axis (e.g. 8 means 8×8 grid)
//   - clipLimit: contrast limit per histogram bin (e.g. 2.0)
//
// Reference: Karel Zuiderveld, "Contrast Limited Adaptive Histogram
// Equalization," in Graphics Gems IV, Academic Press, 1994.

// CLAHE applies Contrast Limited Adaptive Histogram Equalization to a
// single-channel image. pixels is a row-major float32 array in [0, 1].
// nx and ny are the image dimensions.
func CLAHE(pixels []float32, nx, ny, tileGridSize int, clipLimit float64) {
	if tileGridSize <= 0 || clipLimit <= 0 {
		return
	}
	nBins := 256

	// Tile dimensions.
	tileW := nx / tileGridSize
	tileH := ny / tileGridSize
	if tileW < 1 {
		tileW = 1
	}
	if tileH < 1 {
		tileH = 1
	}

	// Effective grid (may differ from tileGridSize if image doesn't
	// divide evenly — we use the integer tile size and let the last
	// tile extend to the edge).
	gridX := nx / tileW
	gridY := ny / tileH
	if gridX < 1 {
		gridX = 1
	}
	if gridY < 1 {
		gridY = 1
	}

	// Precompute the CDF (mapping) for each tile.
	// maps[gy][gx] is a [nBins]float32 lookup table.
	maps := make([][][]float32, gridY)
	for gy := range maps {
		maps[gy] = make([][]float32, gridX)
		for gx := range maps[gy] {
			// Tile bounds.
			x0 := gx * tileW
			y0 := gy * tileH
			x1 := x0 + tileW
			y1 := y0 + tileH
			if gx == gridX-1 {
				x1 = nx
			}
			if gy == gridY-1 {
				y1 = ny
			}

			maps[gy][gx] = claheTileMap(pixels, nx, x0, y0, x1, y1, nBins, clipLimit)
		}
	}

	// Apply: for each pixel, bilinear interpolate between the 4
	// surrounding tile CDFs.
	for py := 0; py < ny; py++ {
		// Which tile row, and fractional position within it.
		tyf := (float64(py) + 0.5) / float64(tileH) - 0.5
		gy0 := int(tyf)
		if gy0 < 0 {
			gy0 = 0
		}
		gy1 := gy0 + 1
		if gy1 >= gridY {
			gy1 = gridY - 1
		}
		if gy0 >= gridY {
			gy0 = gridY - 1
		}
		fy := float32(tyf) - float32(gy0)
		if fy < 0 {
			fy = 0
		}
		if fy > 1 {
			fy = 1
		}

		for px := 0; px < nx; px++ {
			txf := (float64(px) + 0.5) / float64(tileW) - 0.5
			gx0 := int(txf)
			if gx0 < 0 {
				gx0 = 0
			}
			gx1 := gx0 + 1
			if gx1 >= gridX {
				gx1 = gridX - 1
			}
			if gx0 >= gridX {
				gx0 = gridX - 1
			}
			fx := float32(txf) - float32(gx0)
			if fx < 0 {
				fx = 0
			}
			if fx > 1 {
				fx = 1
			}

			idx := py*nx + px
			v := pixels[idx]
			bin := int(v * float32(nBins-1))
			if bin < 0 {
				bin = 0
			}
			if bin >= nBins {
				bin = nBins - 1
			}

			// Bilinear interpolation of the 4 surrounding tile CDFs.
			v00 := maps[gy0][gx0][bin]
			v10 := maps[gy0][gx1][bin]
			v01 := maps[gy1][gx0][bin]
			v11 := maps[gy1][gx1][bin]
			top := v00*(1-fx) + v10*fx
			bot := v01*(1-fx) + v11*fx
			pixels[idx] = top*(1-fy) + bot*fy
		}
	}
}

// claheTileMap computes the clipped + redistributed CDF for one tile.
func claheTileMap(pixels []float32, stride, x0, y0, x1, y1, nBins int, clipLimit float64) []float32 {
	tilePixels := (x1 - x0) * (y1 - y0)
	if tilePixels == 0 {
		out := make([]float32, nBins)
		for i := range out {
			out[i] = float32(i) / float32(nBins-1)
		}
		return out
	}

	// Build histogram.
	hist := make([]int, nBins)
	for py := y0; py < y1; py++ {
		for px := x0; px < x1; px++ {
			v := pixels[py*stride+px]
			bin := int(v * float32(nBins-1))
			if bin < 0 {
				bin = 0
			}
			if bin >= nBins {
				bin = nBins - 1
			}
			hist[bin]++
		}
	}

	// Clip histogram and redistribute excess.
	limit := int(clipLimit * float64(tilePixels) / float64(nBins))
	if limit < 1 {
		limit = 1
	}
	excess := 0
	for i := range hist {
		if hist[i] > limit {
			excess += hist[i] - limit
			hist[i] = limit
		}
	}
	// Distribute excess evenly.
	perBin := excess / nBins
	remainder := excess - perBin*nBins
	for i := range hist {
		hist[i] += perBin
	}
	// Distribute the remainder one-per-bin starting from 0.
	for i := 0; i < remainder; i++ {
		hist[i]++
	}

	// Build CDF and normalize to [0, 1].
	cdf := make([]float32, nBins)
	sum := 0
	for i := range hist {
		sum += hist[i]
		cdf[i] = float32(sum) / float32(tilePixels)
	}
	return cdf
}

// CLAHERGB applies CLAHE to a 3-channel RGB image by converting to
// a luminance representation, applying CLAHE to the luminance, and
// scaling the original channels proportionally to preserve color.
//
// This matches Siril's approach (which converts to CIE Lab via OpenCV,
// applies CLAHE to the L channel, and converts back).
func CLAHERGB(r, g, b []float32, nx, ny, tileGridSize int, clipLimit float64) {
	// Compute luminance (equal weights, matching Siril's even-lum convention).
	n := len(r)
	lum := make([]float32, n)
	for i := range lum {
		lum[i] = (r[i] + g[i] + b[i]) / 3
	}

	// Apply CLAHE to luminance.
	lumOrig := make([]float32, n)
	copy(lumOrig, lum)
	CLAHE(lum, nx, ny, tileGridSize, clipLimit)

	// Scale RGB channels by the luminance ratio.
	for i := range r {
		if lumOrig[i] <= 0 {
			continue
		}
		scale := lum[i] / lumOrig[i]
		r[i] = clampf32(r[i]*scale, 0, 1)
		g[i] = clampf32(g[i]*scale, 0, 1)
		b[i] = clampf32(b[i]*scale, 0, 1)
	}
}
