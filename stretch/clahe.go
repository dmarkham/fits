package stretch

import (
	"fmt"
	"math"
	"math/big"
	"sync"
)

// CLAHE — Contrast Limited Adaptive Histogram Equalization, as Siril's
// `clahe cliplimit tilesize` command computes it for a 32-bit float image
// (opencv/opencv.cpp:806 cvClahe_float). Siril delegates the work to
// OpenCV; this file ports the OpenCV 5.0.0 pieces that pipeline uses:
//
//   - imgproc/src/clahe.cpp: CLAHE_Impl::apply (tile LUTs with clipping and
//     redistribution, bilinear interpolation, BORDER_REFLECT_101 padding
//     when the image is not a multiple of the tile grid), 8- and 16-bit;
//   - imgproc/src/color_lab.cpp: float BGR2Lab (RGB2Lab_f, the trilinear
//     16-bit LUT path OpenCV uses for sRGB input) and float Lab2BGR
//     (Lab2RGBfloat, including its separate SIMD-block and scalar-tail
//     arithmetic as built for x86-64 SSE2 without FMA);
//   - core: Mat::convertTo between 32F and 16U (float32 scale,
//     round-half-even, saturation).
//
// Tests compare the result with the verbatim Siril function linked
// against OpenCV 5.0.0 and the building blocks with opencv-python 5.0.0,
// bit for bit.

// CLAHE applies Siril's CLAHE to a mono image (nx*ny float32 pixels,
// nominally in [0, 1]): each pixel is converted to 16 bits
// (round(v*65535), saturated), equalized by OpenCV's CLAHE with 65536
// histogram bins on a tileGridSize×tileGridSize grid, and converted back
// (v/65535). As in OpenCV, the per-tile histogram is clipped at
// max(1, int32(clipLimit*tilePixels/65536)) counts, where tilePixels is the
// tile area after padding the image to a multiple of the grid; so on small
// tiles every clipLimit > 0 acts like the minimum limit of 1. The int32
// conversion is x86-64's: when clipLimit*tilePixels/65536 is 2^31 or more
// (or clipLimit is +Inf) it gives INT_MIN, so the limit is 1, the
// strongest, as in Siril (siril-cli accepts such values). clipLimit <= 0
// or NaN disables contrast limiting (OpenCV's clipLimit > 0 test; Siril's
// clahe command rejects clipLimit <= 0 and runs NaN). A NaN
// pixel converts to 16-bit 0 (OpenCV's cvRound gives INT_MIN for NaN,
// saturated to 0), so it is equalized as 0 and the output has no NaN.
//
// It panics if len(pixels) != nx*ny or tileGridSize < 1.
func CLAHE(pixels []float32, nx, ny, tileGridSize int, clipLimit float64) {
	checkCLAHEArgs(len(pixels), nx, ny, tileGridSize)
	u := make([]uint16, len(pixels))
	scale := float32(65535.0)
	for i, v := range pixels {
		u[i] = saturateU16(float32(v * scale))
	}
	out := claheApply(u, nx, ny, tileGridSize, tileGridSize, clipLimit, 65536)
	back := float32(1.0 / 65535.0)
	for i, v := range out {
		pixels[i] = float32(float32(v) * back)
	}
}

// CLAHERGB applies Siril's CLAHE to an RGB image: the pixels are converted
// to CIE L*a*b* with OpenCV's float BGR2Lab, L (0..100) is scaled to 16 bits
// (round(L*655.35)), equalized by OpenCV's 16-bit CLAHE, scaled back, and
// the image is converted back with OpenCV's float Lab2BGR; a* and b* are
// untouched. Inputs are clamped to [0, 1] by the Lab conversion and the
// output is in [0, 1]; NaN components go through OpenCV's clamps (SSE
// min/max semantics), so the output has no NaN.
//
// It panics if the channel lengths differ from nx*ny or tileGridSize < 1.
func CLAHERGB(r, g, b []float32, nx, ny, tileGridSize int, clipLimit float64) {
	n := nx * ny
	if len(g) != len(r) || len(b) != len(r) {
		panic(fmt.Sprintf("stretch: CLAHERGB channel lengths %d, %d, %d differ", len(r), len(g), len(b)))
	}
	checkCLAHEArgs(len(r), nx, ny, tileGridSize)
	bgr := make([]float32, 3*n)
	for i := 0; i < n; i++ { // fits_to_bgrbgr_float
		bgr[3*i], bgr[3*i+1], bgr[3*i+2] = b[i], g[i], r[i]
	}
	lab := make([]float32, 3*n)
	for y := 0; y < ny; y++ {
		bgr2labRow(bgr[3*y*nx:3*(y+1)*nx], lab[3*y*nx:3*(y+1)*nx])
	}
	L := make([]uint16, n)
	toU16 := float32(65535.0 / 100.0)
	for i := 0; i < n; i++ {
		L[i] = saturateU16(float32(lab[3*i] * toU16))
	}
	L = claheApply(L, nx, ny, tileGridSize, tileGridSize, clipLimit, 65536)
	fromU16 := float32(100.0 / 65535.0)
	for i := 0; i < n; i++ {
		lab[3*i] = float32(float32(L[i]) * fromU16)
	}
	for y := 0; y < ny; y++ {
		lab2bgrRow(lab[3*y*nx:3*(y+1)*nx], bgr[3*y*nx:3*(y+1)*nx])
	}
	for i := 0; i < n; i++ {
		b[i], g[i], r[i] = bgr[3*i], bgr[3*i+1], bgr[3*i+2]
	}
}

func checkCLAHEArgs(n, nx, ny, tile int) {
	if nx < 1 || ny < 1 || n != nx*ny {
		panic(fmt.Sprintf("stretch: CLAHE %d pixels for a %dx%d image", n, nx, ny))
	}
	if tile < 1 {
		panic(fmt.Sprintf("stretch: CLAHE tileGridSize %d < 1", tile))
	}
}

// cvRound is OpenCV's cvRound(float) on x86-64 (cvtss2si): round half to
// even; NaN and values outside int32 give INT_MIN.
func cvRound(v float32) int32 {
	if !(v >= -2147483648 && v < 2147483648) { // also NaN
		return math.MinInt32
	}
	return int32(math.RoundToEven(float64(v)))
}

// saturateU16 is saturate_cast<ushort>(float).
func saturateU16(v float32) uint16 {
	i := cvRound(v)
	if i < 0 {
		return 0
	}
	if i > 65535 {
		return 65535
	}
	return uint16(i)
}

// borderReflect101 is OpenCV's borderInterpolate(p, n, BORDER_REFLECT_101).
func borderReflect101(p, n int) int {
	if n == 1 {
		return 0
	}
	for p < 0 || p >= n {
		if p < 0 {
			p = -p
		} else {
			p = 2*n - 2 - p
		}
	}
	return p
}

// cvttsd2si is static_cast<int>(double) as compiled for x86-64 (the
// SSE2 instruction cvttsd2si): truncation toward zero, and the "integer
// indefinite" value INT_MIN when the truncated value does not fit in an
// int32 (|v| >= 2^31 after truncation, ±Inf, NaN). OpenCV's CLAHE converts
// its clip limit this way, so a clip limit whose clip*tilePixels/histSize
// is 2^31 or more becomes INT_MIN and then max(INT_MIN, 1) = 1. Go's own
// float-to-int conversion of such a value is implementation-specific
// (it differs between amd64 and arm64), so the result is spelled out.
func cvttsd2si(v float64) int32 {
	if v > -2147483649 && v < 2147483648 {
		return int32(v)
	}
	return math.MinInt32
}

// claheApply is OpenCV's CLAHE_Impl::apply with bitShift 0 on a w×h
// image whose values are < histSize (256 for 8-bit, 65536 for 16-bit).
func claheApply(src []uint16, w, h, tilesX, tilesY int, clipLimitIn float64, histSize int) []uint16 {
	maxVal := histSize - 1
	// Extend the image to a multiple of the tile grid when needed.
	lutSrc, ew := src, w
	if w%tilesX != 0 || h%tilesY != 0 {
		ew = w + tilesX - w%tilesX
		eh := h + tilesY - h%tilesY
		lutSrc = make([]uint16, ew*eh)
		for y := 0; y < eh; y++ {
			sy := borderReflect101(y, h)
			for x := 0; x < ew; x++ {
				lutSrc[y*ew+x] = src[sy*w+borderReflect101(x, w)]
			}
		}
	}
	tw := ew / tilesX
	th := (len(lutSrc) / ew) / tilesY
	tileTotal := tw * th
	lutScale := float32(histSize-1) / float32(tileTotal)
	clipLimit := 0
	if clipLimitIn > 0 {
		clipLimit = int(cvttsd2si(clipLimitIn * float64(tileTotal) / float64(histSize)))
		if clipLimit < 1 {
			clipLimit = 1
		}
	}

	// CLAHE_CalcLut_Body
	lut := make([]uint16, tilesX*tilesY*histSize)
	hist := make([]int, histSize)
	for k := 0; k < tilesX*tilesY; k++ {
		ty, tx := k/tilesX, k%tilesX
		for i := range hist {
			hist[i] = 0
		}
		for y := ty * th; y < (ty+1)*th; y++ {
			row := lutSrc[y*ew+tx*tw : y*ew+(tx+1)*tw]
			for _, v := range row {
				hist[v]++
			}
		}
		if clipLimit > 0 {
			clipped := 0
			for i := range hist {
				if hist[i] > clipLimit {
					clipped += hist[i] - clipLimit
					hist[i] = clipLimit
				}
			}
			redistBatch := clipped / histSize
			residual := clipped - redistBatch*histSize
			for i := range hist {
				hist[i] += redistBatch
			}
			if residual != 0 {
				residualStep := histSize / residual
				if residualStep < 1 {
					residualStep = 1
				}
				for i := 0; i < histSize && residual > 0; i, residual = i+residualStep, residual-1 {
					hist[i]++
				}
			}
		}
		sum := 0
		tl := lut[k*histSize : (k+1)*histSize]
		for i, c := range hist {
			sum += c
			tl[i] = saturateTo(cvRound(float32(float32(sum)*lutScale)), maxVal)
		}
	}

	// CLAHE_Interpolation_Body
	dst := make([]uint16, w*h)
	invTW := 1 / float32(tw)
	invTH := 1 / float32(th)
	ind1 := make([]int, w)
	ind2 := make([]int, w)
	xa := make([]float32, w)
	xa1 := make([]float32, w)
	for x := 0; x < w; x++ {
		txf := float32(float32(x)*invTW) - 0.5
		tx1 := int(math.Floor(float64(txf)))
		tx2 := tx1 + 1
		xa[x] = txf - float32(tx1)
		xa1[x] = 1 - xa[x]
		if tx1 < 0 {
			tx1 = 0
		}
		if tx2 > tilesX-1 {
			tx2 = tilesX - 1
		}
		ind1[x] = tx1 * histSize
		ind2[x] = tx2 * histSize
	}
	for y := 0; y < h; y++ {
		tyf := float32(float32(y)*invTH) - 0.5
		ty1 := int(math.Floor(float64(tyf)))
		ty2 := ty1 + 1
		ya := tyf - float32(ty1)
		ya1 := 1 - ya
		if ty1 < 0 {
			ty1 = 0
		}
		if ty2 > tilesY-1 {
			ty2 = tilesY - 1
		}
		p1 := lut[ty1*tilesX*histSize:]
		p2 := lut[ty2*tilesX*histSize:]
		for x := 0; x < w; x++ {
			v := int(src[y*w+x])
			i1, i2 := ind1[x]+v, ind2[x]+v
			top := float32(float32(p1[i1])*xa1[x]) + float32(float32(p1[i2])*xa[x])
			bot := float32(float32(p2[i1])*xa1[x]) + float32(float32(p2[i2])*xa[x])
			res := float32(top*ya1) + float32(bot*ya)
			dst[y*w+x] = saturateTo(cvRound(res), maxVal)
		}
	}
	return dst
}

func saturateTo(v int32, maxVal int) uint16 {
	if v < 0 {
		return 0
	}
	if int(v) > maxVal {
		return uint16(maxVal)
	}
	return uint16(v)
}

// ---------- OpenCV float BGR <-> Lab (imgproc/src/color_lab.cpp) ----------

// OpenCV's constants, bit for bit (color_lab.cpp: sRGB2XYZ_D65,
// XYZ2sRGB_D65, D65 as raw doubles).
var (
	ocvSRGB2XYZ = [9]float64{
		f64bits(0x3fda65a14488c60d), f64bits(0x3fd6e297396d0918), f64bits(0x3fc71819d2391d58),
		f64bits(0x3fcb38cda6e75ff6), f64bits(0x3fe6e297396d0918), f64bits(0x3fb279aae6c8f755),
		f64bits(0x3f93cc4ac6cdaf4b), f64bits(0x3fbe836eb4e98138), f64bits(0x3fee68427418d691),
	}
	ocvXYZ2SRGB = [9]float64{
		f64bits(0x4009ec804102ff8f), f64bits(0xbff8982a9930be0e), f64bits(0xbfdfe7ff583a53b9),
		f64bits(0xbfef042528ae74f3), f64bits(0x3ffe040f23897204), f64bits(0x3fa546d3f9e7b80b),
		f64bits(0x3fac7de5082cf52c), f64bits(0xbfca1e14bdfd2631), f64bits(0x3ff0eabef06b3786),
	}
	ocvD65 = [3]float64{f64bits(0x3fee6a22b3892ee8), 1, f64bits(0x3ff16b8950763a19)}
)

func f64bits(b uint64) float64 { return math.Float64frombits(b) }

const (
	labLUTDim       = 33 // (1 << lab_lut_shift) + 1
	labBase         = 1 << 14
	trilinearShift  = 4 // 8 - lab_lut_shift + 1
	trilinearBase   = 1 << trilinearShift
	gammaTabSize    = 1024
	gammaThreshold  = 809.0 / 20000.0
	gammaInvThresh  = 7827.0 / 2500000.0
	gammaLowScale   = 323.0 / 25.0
	gammaPower      = 12.0 / 5.0
	gammaXshift     = 11.0 / 200.0
	labCvtBlockSize = 8 // SIMD block: 4 float32 lanes (SSE2) x 2
)

var (
	labOnce   sync.Once
	labLUT    []int16   // [p][q][r][3] for the 33^3 grid, p = channel 0 (B)
	invGamma  []float32 // sRGBInvGammaTab: 1024 cubic spline segments x 4
	lab2bgrC  [9]float32
	labFThres float32
)

// fma32 is softfloat mulAdd: a*b+c rounded once to float32.
func fma32(a, b, c float32) float32 {
	x := new(big.Float).SetPrec(256).SetFloat64(float64(a))
	x.Mul(x, new(big.Float).SetPrec(256).SetFloat64(float64(b)))
	x.Add(x, new(big.Float).SetPrec(256).SetFloat64(float64(c)))
	f, _ := x.Float32()
	return f
}

// softCbrt is OpenCV's softfloat cbrt (core/src/softfloat.cpp f32_cbrt,
// Turkowski's method): a rational approximation evaluated in double on the
// mantissa scaled to [0.125, 1), whose result mantissa is TRUNCATED to 23
// bits. It is often 1 ULP below the correctly rounded cube root, and the
// RGB->Lab LUT depends on it.
func softCbrt(x float32) float32 {
	if x != x {
		return float32(math.NaN())
	}
	if math.IsInf(float64(x), 0) {
		return x
	}
	ux := math.Float32bits(x)
	s := ux >> 31
	ex := int(ux>>23&0xff) - 127
	shx := ex % 3
	if shx >= 0 {
		shx -= 3
	}
	ex = (ex-shx)/3 - 1
	fr := math.Float64frombits(uint64(shx+1023)<<52 | uint64(ux&0x7fffff)<<29)
	A := [10]float64{
		f64bits(0x4046a09e6653ba70), f64bits(0x406808f46c6116e0), f64bits(0x405dca97439cae14),
		f64bits(0x402add70d2827500), f64bits(0x3fc4f15f83f55d2d), f64bits(0x402d9e20660edb21),
		f64bits(0x4062ff15c0285815), f64bits(0x406510d06a8112ce), f64bits(0x4040fecbc9e2c375), 1,
	}
	num := float64(float64(float64(float64(float64(float64(float64(float64(A[0]*fr)+A[1])*fr)+A[2])*fr)+A[3])*fr) + A[4])
	den := float64(float64(float64(float64(float64(float64(float64(float64(A[5]*fr)+A[6])*fr)+A[7])*fr)+A[8])*fr) + A[9])
	fr = num / den
	if ux&0x7fffffff == 0 {
		return 0
	}
	// packToF32UI(s, ex+127, frac): the fields are added, so an exponent
	// carry is possible exactly as in the C code.
	frac := uint32(math.Float64bits(fr) & (1<<52 - 1) >> 29)
	return math.Float32frombits(s<<31 + uint32(ex+127)<<23 + frac)
}

// applyGamma / applyInvGamma are color_lab.cpp's sRGB transfer functions
// (computed in double, result rounded to float).
func applyGamma(x float32) float32 {
	xd := float64(x)
	if xd <= gammaThreshold {
		return float32(xd / gammaLowScale)
	}
	return float32(math.Pow((xd+gammaXshift)/(1+gammaXshift), gammaPower))
}

func applyInvGamma(x float32) float32 {
	xd := float64(x)
	if xd <= gammaInvThresh {
		return float32(xd * gammaLowScale)
	}
	return float32(float64(math.Pow(xd, 1/gammaPower)*(1+gammaXshift)) - gammaXshift)
}

// splineBuild is color_lab.cpp splineBuild (float32 arithmetic).
func splineBuild(f []float32, n int) []float32 {
	tab := make([]float32, n*4)
	for i := 1; i < n; i++ {
		t := float32(float32(float32(f[i+1]-float32(f[i]*2))+f[i-1]) * 3)
		l := 1 / (4 - tab[(i-1)*4])
		tab[i*4] = l
		tab[i*4+1] = float32((t - tab[(i-1)*4+1]) * l)
	}
	var cn float32
	for j := 0; j < n; j++ {
		i := n - j - 1
		c := tab[i*4+1] - float32(tab[i*4]*cn)
		b := (f[i+1] - f[i]) - (cn+float32(c*2))/3
		d := (cn - c) / 3
		tab[i*4], tab[i*4+1], tab[i*4+2], tab[i*4+3] = f[i], b, c, d
		cn = c
	}
	return tab
}

// splineInterpolate is color_lab.cpp splineInterpolate (scalar and SIMD
// versions compute the same expression).
func splineInterpolate(x float32, tab []float32, n int) float32 {
	ix := int(x) // truncation, as int(x) / v_trunc
	if ix < 0 {
		ix = 0
	}
	if ix > n-1 {
		ix = n - 1
	}
	x -= float32(ix)
	t := tab[ix*4:]
	return float32(float32(float32(float32(float32(t[3]*x)+t[2])*x)+t[1])*x) + t[0]
}

func initLab() {
	// RGB2Lab LUT (initLUTforLABLUVs16, Lab part). Channel order is the
	// BGR input order: p indexes channel 0 (blue), so its coefficients are
	// the blue column of sRGB2XYZ.
	scaleWhite := [3]float64{1 / ocvD65[0], 1, 1 / ocvD65[2]}
	var S [9]float32
	for i := 0; i < 3; i++ {
		S[i*3+0] = float32(ocvSRGB2XYZ[i*3+2] * scaleWhite[i])
		S[i*3+1] = float32(ocvSRGB2XYZ[i*3+1] * scaleWhite[i])
		S[i*3+2] = float32(ocvSRGB2XYZ[i*3+0] * scaleWhite[i])
	}
	lthresh := float32(216) / float32(24389)
	lscale := float32(841) / float32(108)
	lbias := float32(16) / float32(116)
	f9033 := float32(29*29*29) / float32(27)
	lab := func(X float32) float32 {
		if X > lthresh {
			return softCbrt(X)
		}
		return fma32(X, lscale, lbias)
	}
	labLUT = make([]int16, labLUTDim*labLUTDim*labLUTDim*3)
	for p := 0; p < labLUTDim; p++ {
		for q := 0; q < labLUTDim; q++ {
			for r := 0; r < labLUTDim; r++ {
				c0 := applyGamma(float32(p) / float32(labLUTDim-1))
				c1 := applyGamma(float32(q) / float32(labLUTDim-1))
				c2 := applyGamma(float32(r) / float32(labLUTDim-1))
				X := float32(float32(c0*S[0])+float32(c1*S[1])) + float32(c2*S[2])
				Y := float32(float32(c0*S[3])+float32(c1*S[4])) + float32(c2*S[5])
				Z := float32(float32(c0*S[6])+float32(c1*S[7])) + float32(c2*S[8])
				FX, FY, FZ := lab(X), lab(Y), lab(Z)
				var L float32
				if Y > lthresh {
					L = float32(116*FY) - 16
				} else {
					L = float32(f9033 * Y)
				}
				a := float32(500 * (FX - FY))
				bb := float32(200 * (FY - FZ))
				idx := ((p*labLUTDim+q)*labLUTDim + r) * 3
				labLUT[idx] = int16(cvRound(float32(labBase*L) / 100))
				labLUT[idx+1] = int16(cvRound(float32(labBase*(a+128)) / 256))
				labLUT[idx+2] = int16(cvRound(float32(labBase*(bb+128)) / 256))
			}
		}
	}

	// sRGBInvGammaTab (createLabTabs)
	ig := make([]float32, gammaTabSize+1)
	scale := 1 / float32(gammaTabSize)
	for i := range ig {
		ig[i] = applyInvGamma(float32(scale * float32(i)))
	}
	invGamma = splineBuild(ig, gammaTabSize)

	// Lab2RGBfloat coefficients for Lab2BGR (blueIdx 0): row k of the
	// output is B, G, R.
	for i := 0; i < 3; i++ {
		lab2bgrC[i+6] = float32(ocvXYZ2SRGB[i] * ocvD65[i])
		lab2bgrC[i+3] = float32(ocvXYZ2SRGB[i+3] * ocvD65[i])
		lab2bgrC[i] = float32(ocvXYZ2SRGB[i+6] * ocvD65[i])
	}
	labFThres = float32(6) / float32(29)
}

// trilinear is color_lab.cpp trilinearInterpolate on labLUT.
func trilinear(cx, cy, cz int) (L, a, b int) {
	tx := cx >> (14 - 5)
	ty := cy >> (14 - 5)
	tz := cz >> (14 - 5)
	const mask = trilinearBase - 1
	x := (cx >> (14 - 8 - 1)) & mask
	y := (cy >> (14 - 8 - 1)) & mask
	z := (cz >> (14 - 8 - 1)) & mask
	var acc [3]int
	for i := 0; i < 8; i++ {
		dp, dq, dr := i>>2&1, i>>1&1, i&1
		w := 1
		for _, f := range [3][2]int{{dp, x}, {dq, y}, {dr, z}} {
			if f[0] == 1 {
				w *= f[1]
			} else {
				w *= trilinearBase - f[1]
			}
		}
		p := min(tx+dp, labLUTDim-1)
		q := min(ty+dq, labLUTDim-1)
		r := min(tz+dr, labLUTDim-1)
		idx := ((p*labLUTDim+q)*labLUTDim + r) * 3
		for k := 0; k < 3; k++ {
			acc[k] += int(labLUT[idx+k]) * w
		}
	}
	const sh = trilinearShift * 3
	return (acc[0] + 1<<(sh-1)) >> sh, (acc[1] + 1<<(sh-1)) >> sh, (acc[2] + 1<<(sh-1)) >> sh
}

// bgr2labRow is cvtColor(COLOR_BGR2Lab) on one row of float32 BGR pixels.
func bgr2labRow(src, dst []float32) {
	labOnce.Do(initLab)
	n := len(src) / 3
	for i := 0; i < n; i++ {
		simd := i < n/labCvtBlockSize*labCvtBlockSize
		var ic [3]int
		for k := 0; k < 3; k++ {
			v := src[3*i+k]
			if simd { // v_min(v_max(v, 0), 1) with SSE NaN semantics
				if !(v > 0) {
					v = 0
				}
				if !(v < 1) {
					v = 1
				}
			} else if v < 0 { // clip(): value < 0 ? 0 : value <= 1 ? value : 1
				v = 0
			} else if !(v <= 1) {
				v = 1
			}
			ic[k] = int(cvRound(float32(v * labBase)))
		}
		L, a, b := trilinear(ic[0], ic[1], ic[2])
		dst[3*i] = float32(float32(float32(L)*1/labBase) * 100)
		dst[3*i+1] = float32(float32(float32(a)*1/labBase)*256) - 128
		dst[3*i+2] = float32(float32(float32(b)*1/labBase)*256) - 128
	}
}

// lab2bgrRow is cvtColor(COLOR_Lab2BGR) on one row of float32 Lab pixels:
// blocks of 8 pixels use the SIMD formulas of Lab2RGBfloat::operator(),
// the remaining pixels its scalar formulas (they round differently).
func lab2bgrRow(src, dst []float32) {
	labOnce.Do(initLab)
	n := len(src) / 3
	C := lab2bgrC
	const lThresh = float32(8)
	var (
		inv903   = float32(1) / float32(903.3)
		v7787    = float32(7.787)
		v16_116  = float32(16) / float32(116)
		inv116   = float32(1) / float32(116)
		pinv500  = float32(1) / float32(500)
		ninv200  = float32(-1) / float32(200)
		inv7787  = float32(1) / float32(7.787)
		f9033    = float32(903.3)
		gscale   = float32(gammaTabSize)
		fThresh  = labFThres
		blockEnd = n / labCvtBlockSize * labCvtBlockSize
	)
	for i := 0; i < n; i++ {
		li, ai, bi := src[3*i], src[3*i+1], src[3*i+2]
		var y, fy float32
		var out [3]float32
		if i < blockEnd {
			if li <= lThresh {
				y = float32(li * inv903)
				fy = float32(v7787*y) + v16_116
			} else {
				fy = float32((li + 16) * inv116)
				y = float32(float32(fy*fy) * fy)
			}
			fxz := [2]float32{float32(ai*pinv500) + fy, float32(bi*ninv200) + fy}
			for j := range fxz {
				if f := fxz[j]; f <= fThresh {
					fxz[j] = float32((f - v16_116) * inv7787)
				} else {
					fxz[j] = float32(float32(f*f) * f)
				}
			}
			x, z := fxz[0], fxz[1]
			for k := 0; k < 3; k++ {
				v := float32(C[3*k]*x) + float32(float32(C[3*k+1]*y)+float32(C[3*k+2]*z))
				// v_max(0, v_min(v, 1)), SSE NaN semantics
				if !(v < 1) {
					v = 1
				}
				if !(v > 0) {
					v = 0
				}
				out[k] = splineInterpolate(float32(v*gscale), invGamma, gammaTabSize)
			}
		} else {
			if li <= lThresh {
				y = li / f9033
				fy = float32(v7787*y) + v16_116
			} else {
				fy = (li + 16) / 116
				y = float32(float32(fy*fy) * fy)
			}
			fxz := [2]float32{ai/500 + fy, fy - bi/200}
			for j := range fxz {
				if f := fxz[j]; f <= fThresh {
					fxz[j] = (f - v16_116) / v7787
				} else {
					fxz[j] = float32(float32(f*f) * f)
				}
			}
			x, z := fxz[0], fxz[1]
			for k := 0; k < 3; k++ {
				v := float32(float32(C[3*k]*x)+float32(C[3*k+1]*y)) + float32(C[3*k+2]*z)
				if v < 0 { // clip()
					v = 0
				} else if !(v <= 1) {
					v = 1
				}
				out[k] = splineInterpolate(float32(v*gscale), invGamma, gammaTabSize)
			}
		}
		dst[3*i], dst[3*i+1], dst[3*i+2] = out[0], out[1], out[2]
	}
}
