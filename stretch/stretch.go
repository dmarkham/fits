// Package stretch implements image stretching operations for FITS
// astronomical images. These are the same algorithms used by Siril,
// ported to pure Go and cross-validated against siril-cli output on
// a comprehensive golden fixture set.
//
// The package operates on normalized float32 pixel data in [0, 1],
// where 0 is black and 1 is the maximum signal. Multi-channel images
// are represented as separate float32 slices per channel.
//
// # Operations
//
//   - [Linear] — clip at black point, linear scale to [0, 1]
//   - [Asinh] — inverse hyperbolic sine stretch (preserves color ratios)
//   - [MTF] — Midtone Transfer Function (KStars/PixInsight autostretch)
//   - [GHT] — Generalized Hyperbolic Transformation (Strasser 2022)
//   - [CLAHE] — Contrast Limited Adaptive Histogram Equalization
//
// # Validation
//
// Every operation is cross-validated against siril-cli 1.5.0. The
// golden fixtures in stretch/testdata/golden/ contain per-pixel float32
// output from Siril for a synthetic 3-channel 64×48 test image.
package stretch

import "math"

// Linear applies a linear stretch with black point clipping.
// For each pixel: output = max(0, (input - bp) / (1 - bp)).
//
// Port of Siril's linstretch command (ght.c:291, a special case of
// GHT with D=0, B=0). Applied independently per channel.
func Linear(pixels []float32, bp float32) {
	if bp <= 0 {
		return // identity
	}
	scale := 1.0 / (1.0 - bp)
	for i, v := range pixels {
		v = (v - bp) * scale
		if v < 0 {
			v = 0
		}
		pixels[i] = v
	}
}

// Asinh applies an inverse-hyperbolic-sine stretch to a single channel
// (mono image). For each pixel:
//
//	xprime = max(0, (x - offset) / (1 - offset))
//	output = clamp(asinh(beta * xprime) / asinh(beta), 0, 1)
//
// Port of Siril's asinh command (asinh.c). For multi-channel images,
// use AsinhRGB which preserves color ratios via luminance-based scaling.
func Asinh(pixels []float32, beta, offset float32) {
	if beta == 0 {
		return // identity
	}
	invOneMinusOff := float32(1) / (1 - offset)
	asinhBeta := float32(math.Asinh(float64(beta)))
	for i, v := range pixels {
		xp := (v - offset) * invOneMinusOff
		if xp <= 0 {
			pixels[i] = 0
			continue
		}
		out := float32(math.Asinh(float64(beta*xp))) / asinhBeta
		if out > 1 {
			out = 1
		}
		pixels[i] = out
	}
}

// AsinhRGB applies an asinh stretch to a 3-channel RGB image, preserving
// color ratios by stretching the luminance and scaling all channels
// proportionally. This matches Siril's default asinh behavior for color
// images (RGBBLEND mode with equal luminance weights).
//
// r, g, b must have the same length (one element per pixel).
func AsinhRGB(r, g, b []float32, beta, offset float32) {
	if beta == 0 {
		return
	}
	invOneMinusOff := float32(1) / (1 - offset)
	asinhBeta := float32(math.Asinh(float64(beta)))
	const third = float32(1.0 / 3.0)
	for i := range r {
		rp := maxf32(0, (r[i]-offset)*invOneMinusOff)
		gp := maxf32(0, (g[i]-offset)*invOneMinusOff)
		bp := maxf32(0, (b[i]-offset)*invOneMinusOff)
		lum := third*rp + third*gp + third*bp
		if lum <= 0 {
			r[i], g[i], b[i] = 0, 0, 0
			continue
		}
		k := float32(math.Asinh(float64(beta*lum))) / (lum * asinhBeta)
		r[i] = clampf32(rp*k, 0, 1)
		g[i] = clampf32(gp*k, 0, 1)
		b[i] = clampf32(bp*k, 0, 1)
	}
}

func maxf32(a, b float32) float32 {
	if a > b {
		return a
	}
	return b
}

func clampf32(v, lo, hi float32) float32 {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}
