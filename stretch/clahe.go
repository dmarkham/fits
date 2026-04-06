package stretch

import "math"

// CLAHE — Contrast Limited Adaptive Histogram Equalization.
//
// Implements the Zuiderveld 1994 algorithm matching OpenCV's pipeline
// (which Siril delegates to via createCLAHE). For RGB images, CLAHE
// is applied to the L channel of CIE Lab color space, preserving
// chrominance — exactly the same path Siril/OpenCV uses.
//
// The implementation quantizes to uint8 internally (matching OpenCV's
// behavior for float input) and uses 256 histogram bins.

// CLAHE applies Contrast Limited Adaptive Histogram Equalization to a
// single-channel uint8 image and returns the result as uint8.
func claheUint8(pixels []uint8, nx, ny, tileGridSize int, clipLimit float64) {
	nBins := 256
	tileW := nx / tileGridSize
	tileH := ny / tileGridSize
	if tileW < 1 {
		tileW = 1
	}
	if tileH < 1 {
		tileH = 1
	}
	gridX := (nx + tileW - 1) / tileW
	gridY := (ny + tileH - 1) / tileH

	// Precompute CDF for each tile.
	maps := make([][][]uint8, gridY)
	for gy := range maps {
		maps[gy] = make([][]uint8, gridX)
		for gx := range maps[gy] {
			x0, y0 := gx*tileW, gy*tileH
			x1, y1 := x0+tileW, y0+tileH
			if x1 > nx {
				x1 = nx
			}
			if y1 > ny {
				y1 = ny
			}
			maps[gy][gx] = claheTileMapU8(pixels, nx, x0, y0, x1, y1, nBins, clipLimit)
		}
	}

	// Bilinear interpolation between 4 surrounding tile CDFs.
	for py := 0; py < ny; py++ {
		tyf := (float64(py)+0.5)/float64(tileH) - 0.5
		gy0 := int(tyf)
		if gy0 < 0 {
			gy0 = 0
		}
		gy1 := gy0 + 1
		if gy0 >= gridY {
			gy0 = gridY - 1
		}
		if gy1 >= gridY {
			gy1 = gridY - 1
		}
		fy := float32(tyf) - float32(gy0)
		if fy < 0 {
			fy = 0
		}
		if fy > 1 {
			fy = 1
		}
		for px := 0; px < nx; px++ {
			txf := (float64(px)+0.5)/float64(tileW) - 0.5
			gx0 := int(txf)
			if gx0 < 0 {
				gx0 = 0
			}
			gx1 := gx0 + 1
			if gx0 >= gridX {
				gx0 = gridX - 1
			}
			if gx1 >= gridX {
				gx1 = gridX - 1
			}
			fx := float32(txf) - float32(gx0)
			if fx < 0 {
				fx = 0
			}
			if fx > 1 {
				fx = 1
			}
			idx := py*nx + px
			bin := pixels[idx]
			v00 := float32(maps[gy0][gx0][bin])
			v10 := float32(maps[gy0][gx1][bin])
			v01 := float32(maps[gy1][gx0][bin])
			v11 := float32(maps[gy1][gx1][bin])
			top := v00*(1-fx) + v10*fx
			bot := v01*(1-fx) + v11*fx
			val := top*(1-fy) + bot*fy
			pixels[idx] = uint8(val + 0.5)
		}
	}
}

func claheTileMapU8(pixels []uint8, stride, x0, y0, x1, y1, nBins int, clipLimit float64) []uint8 {
	tilePix := (x1 - x0) * (y1 - y0)
	if tilePix == 0 {
		out := make([]uint8, nBins)
		for i := range out {
			out[i] = uint8(i)
		}
		return out
	}
	hist := make([]int, nBins)
	for py := y0; py < y1; py++ {
		for px := x0; px < x1; px++ {
			hist[pixels[py*stride+px]]++
		}
	}
	limit := int(clipLimit * float64(tilePix) / float64(nBins))
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
	perBin := excess / nBins
	remainder := excess - perBin*nBins
	for i := range hist {
		hist[i] += perBin
	}
	for i := 0; i < remainder; i++ {
		hist[i]++
	}
	cdf := make([]uint8, nBins)
	sum := 0
	for i := range hist {
		sum += hist[i]
		cdf[i] = uint8(float64(sum) / float64(tilePix) * 255 + 0.5)
	}
	return cdf
}

// ---------- CIE Lab conversion (sRGB D65) ----------
// Matches OpenCV's BGR2Lab / Lab2BGR for uint8 data.

func srgbToLinear(v float64) float64 {
	if v <= 0.04045 {
		return v / 12.92
	}
	return math.Pow((v+0.055)/1.055, 2.4)
}

func linearToSrgb(v float64) float64 {
	if v <= 0.0031308 {
		return v * 12.92
	}
	return 1.055*math.Pow(v, 1.0/2.4) - 0.055
}

func labF(t float64) float64 {
	if t > 0.008856 {
		return math.Cbrt(t)
	}
	return 7.787*t + 16.0/116.0
}

func labFInv(t float64) float64 {
	if t > 0.206893 {
		return t * t * t
	}
	return (t - 16.0/116.0) / 7.787
}

// D65 reference white.
const (
	xn = 0.950456
	yn = 1.0
	zn = 1.088754
)

func rgbToLab(r, g, b float64) (L, a, bOut float64) {
	// sRGB [0,1] → linear
	rl := srgbToLinear(r)
	gl := srgbToLinear(g)
	bl := srgbToLinear(b)
	// linear RGB → XYZ (sRGB D65 matrix)
	x := 0.4124564*rl + 0.3575761*gl + 0.1804375*bl
	y := 0.2126729*rl + 0.7151522*gl + 0.0721750*bl
	z := 0.0193339*rl + 0.1191920*gl + 0.9503041*bl
	// XYZ → Lab
	fx := labF(x / xn)
	fy := labF(y / yn)
	fz := labF(z / zn)
	L = 116*fy - 16
	a = 500 * (fx - fy)
	bOut = 200 * (fy - fz)
	return
}

func labToRGB(L, a, bIn float64) (r, g, b float64) {
	fy := (L + 16) / 116
	fx := a/500 + fy
	fz := fy - bIn/200
	x := xn * labFInv(fx)
	y := yn * labFInv(fy)
	z := zn * labFInv(fz)
	// XYZ → linear sRGB (inverse matrix)
	rl := 3.2404542*x - 1.5371385*y - 0.4985314*z
	gl := -0.9692660*x + 1.8760108*y + 0.0415560*z
	bl := 0.0556434*x - 0.2040259*y + 1.0572252*z
	r = linearToSrgb(math.Max(0, math.Min(1, rl)))
	g = linearToSrgb(math.Max(0, math.Min(1, gl)))
	b = linearToSrgb(math.Max(0, math.Min(1, bl)))
	return
}

// CLAHERGB applies CLAHE to a 3-channel RGB image using the CIE Lab
// color space, matching the OpenCV/Siril pipeline. The L channel is
// quantized to uint8, CLAHE'd, then the result is converted back to
// float32 RGB.
func CLAHERGB(r, g, b []float32, nx, ny, tileGridSize int, clipLimit float64) {
	n := len(r)
	// Convert RGB → Lab, extract L as uint8.
	lCh := make([]uint8, n)
	aCh := make([]float64, n)
	bCh := make([]float64, n)
	for i := range r {
		L, a, bv := rgbToLab(float64(r[i]), float64(g[i]), float64(b[i]))
		// OpenCV's uint8 Lab: L in [0, 255] mapped from [0, 100]
		lVal := L * 255 / 100
		if lVal < 0 {
			lVal = 0
		}
		if lVal > 255 {
			lVal = 255
		}
		lCh[i] = uint8(lVal + 0.5)
		aCh[i] = a
		bCh[i] = bv
	}

	// Apply CLAHE to L channel.
	claheUint8(lCh, nx, ny, tileGridSize, clipLimit)

	// Convert back: replace L, keep a/b, convert Lab → RGB.
	for i := range r {
		L := float64(lCh[i]) * 100 / 255
		rv, gv, bv := labToRGB(L, aCh[i], bCh[i])
		r[i] = float32(rv)
		g[i] = float32(gv)
		b[i] = float32(bv)
	}
}
