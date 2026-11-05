// Package stretch implements image stretching operations for FITS
// astronomical images, ported from Siril 1.4.0 to pure Go.
//
// The package operates on normalized float32 pixel data in [0, 1],
// where 0 is black and 1 is the maximum signal. Multi-channel images
// are represented as separate float32 slices per channel; index i of
// every slice is the same pixel, and a 2-D image is stored row by row
// (index y*nx + x) in the order the caller supplies (nothing is flipped).
//
// # Operations
//
//   - [Linear] — clip at a black point, rescale to [0, 1] (Siril linstretch)
//   - [Asinh], [AsinhRGB], [AsinhRGBHuman] — inverse hyperbolic sine stretch
//     (Siril asinh, asinh -human)
//   - [MTF], [Autostretch], [AutostretchUnlinked], [AutostretchLinked] —
//     Midtone Transfer Function and Siril's autostretch [-linked]
//   - [GHT] — Generalized Hyperbolic Stretch (Siril ght, normal mode,
//     independent channels)
//   - [CLAHE], [CLAHERGB] — Contrast Limited Adaptive Histogram
//     Equalization (Siril clahe, an OpenCV 5 pipeline)
//
// The functions port Siril's float32 (32-bit image) code paths operation by
// operation, so they produce the same float32 values Siril does. Siril's
// 16-bit (ushort) image paths quantize through lookup tables and are not
// what these functions compute.
//
// # Validation
//
// Every exported operation is tested pixel by pixel against an external
// truth source (testdata/PROVENANCE.md):
//
//   - 18 goldens of the [DefaultPresets] commands on a synthetic 64×48 RGB
//     image (testdata/golden/stretch_*.json) whose digests are those of
//     the original siril-cli runs. The Siril 1.4.0 source reproduces all of
//     them bit for bit and writes them (gen_siril_goldens.py), and a
//     siril-cli 1.4.0 built from the release tarball
//     (testdata/build_siril_cli.sh) reproduces them too, together with
//     every oracle golden that is a siril command, ±Inf input included
//     (gen_golden.py).
//   - The accept/refuse decision and output of siril-cli 1.4.0 for
//     commands at the edges of every parameter range
//     (testdata/golden/siril_cli_args.json): [Linear], the Autostretch
//     functions, [GHT], [Asinh] and [MTF] return an error exactly where
//     siril-cli refuses a command, except where the tests list a
//     deliberate difference (NaN parameters and what siril-cli runs into
//     NaN or inverted output are refused; GHT accepts the GUI's B < 0 and
//     D = 0, Asinh the GUI's beta in [0, 1), MTF the finite values Siril's
//     autostretch hands it, CLAHE clip <= 0 as OpenCV's no limit).
//   - The Siril 1.4.0 release source compiled verbatim (testdata/sref,
//     gen_siril_goldens.py), linked against OpenCV 5.0.0 for CLAHE. It
//     covers what the presets do not: mono images, asinh offset and -human,
//     bright images (autostretch's inverted branch), degenerate statistics,
//     MTF with highlights < 1, GHT with B < 0 and with LP > 0 / HP < 1,
//     CLAHE in RGB and mono including sizes that are not a multiple of the
//     tile grid, and edge cases: NaN, ±Inf and negative input for every
//     operation, GHT output above 1, small D, the |B| < float32(0.001)
//     threshold on both sides and exactly at it, HP = 0, linstretch
//     BP = 1, asinh beta = 0 and input above 1, and images whose median
//     sits on either side of autostretch's inverted-branch edge.
//   - opencv-python 5.0.0 for CLAHE's building blocks (8- and 16-bit OpenCV
//     CLAHE, float BGR↔Lab conversion), gen_clahe_cv2.py.
//
// [NormalizeChannels] and [ToRGBImage] are checked against values worked
// out by hand (and with numpy) from their documented rules, and [EncodeJPEG] output is
// decoded by libjpeg-turbo (cv2.imdecode).
//
// The tolerance is none (identical float32 bits, a sha256 over every
// pixel) for every operation except GHT: Siril evaluates GHT with glibc's
// logf/powf/expf, which are not correctly rounded, so a GHT pixel may
// differ where glibc's result is one float32 ULP away from the rounded
// double this package uses. Each GHT golden lists those pixels (0 to 58 of
// 9216 in the committed set), and only there is the value accepted, if
// the one-ULP libm change explains it; see ght_test.go. The goldens were produced
// on x86-64 (glibc 2.44, no FMA in the Siril or OpenCV code paths); the
// Go code uses explicit float32 conversions so the compiler may not fuse
// multiply-adds; the default test run re-runs the package with GOAMD64=v3
// (where the compiler does fuse) and built for arm64 under qemu-aarch64
// (TestFMAGuards_GOAMD64v3, TestCrossArch_arm64; docs/VALIDATION.md,
// Platforms).
package stretch

import (
	"errors"
	"fmt"
	"math"
)

// Linear applies Siril's linstretch (the STRETCH_LINEAR case of
// filters/ght.c, applied by apply_linked_ght_to_fbuf_indep) to one
// channel:
//
//	x      = clamp(input, 0, 1)
//	output = 0                              if x == 0
//	output = max(0, (x - bp) / (1 - bp))    otherwise
//
// bp is the black point, 0 <= bp <= 1, the range siril-cli's
// `linstretch -BP=` accepts (process_ght_args, core/command.c:2349-2356);
// bp = 1 makes every pixel 0, as in Siril. Inputs outside [0, 1] are
// clamped first, as Siril does; NaN becomes 0.
//
// A bp outside [0, 1], or NaN, is an error wrapping ErrLinearParams
// (siril-cli refuses it) and the pixels are left unchanged.
func Linear(pixels []float32, bp float32) error {
	if !(bp >= 0 && bp <= 1) {
		return fmt.Errorf("%w: bp = %v, want 0 <= bp <= 1", ErrLinearParams, bp)
	}
	for i, v := range pixels {
		x := fminf(1, fmaxf(0, v))
		if x == 0 {
			pixels[i] = 0
			continue
		}
		pixels[i] = fmaxf(0, (x-bp)/(1-bp))
	}
	return nil
}

// ErrLinearParams is wrapped by the error Linear returns for a black point
// outside [0, 1].
var ErrLinearParams = errors.New("stretch: linstretch black point out of range")

// Asinh applies Siril's asinh stretch to a single channel (mono image),
// with the float32 arithmetic of asinhlut_float (filters/asinh.c):
//
//	x' = max(0, (x - offset) / (1 - offset))
//	k  = asinh(beta*x') / (x' * asinh(beta))   (0 if x' == 0; 1 if beta == 0)
//	out = min(1, max(0, x' * k))
//
// NaN pixels stay NaN (Siril's max() macro returns its NaN operand).
//
// For RGB images use AsinhRGB or AsinhRGBHuman, which stretch a
// luminance and scale the channels by a common factor.
//
// beta must be finite and >= 0 and offset finite and < 1, else the error
// wraps ErrAsinhParams and the pixels are left unchanged. That range
// includes values siril-cli's `asinh` refuses (beta in [0, 1): Siril's
// GUI range is 0..1000) and gives Siril's result there; it leaves out
// what siril-cli runs into NaN or meaningless output: a NaN or infinite
// beta or offset (NaN or zero pixels) and offset >= 1 (1 - offset <= 0
// divides by zero or inverts the image).
func Asinh(pixels []float32, beta, offset float32) error {
	if err := asinhArgs(beta, offset); err != nil {
		return err
	}
	asinhBeta := asinhf(beta)
	for i, x := range pixels {
		xp := sirilMax(0, (x-offset)/(1-offset))
		var k float32
		switch {
		case xp == 0:
			k = 0
		case beta == 0:
			k = 1
		default:
			k = asinhf(float32(beta*xp)) / (xp * asinhBeta)
		}
		pixels[i] = sirilMin(1, sirilMax(0, float32(xp*k)))
	}
	return nil
}

// ErrAsinhParams is wrapped by the error the Asinh functions return for
// parameters they refuse.
var ErrAsinhParams = errors.New("stretch: asinh parameters out of range")

func asinhArgs(beta, offset float32) error {
	if !(beta >= 0 && beta <= math.MaxFloat32) || !(offset < 1 && offset >= -math.MaxFloat32) {
		return fmt.Errorf("%w: beta %v, offset %v; want finite beta >= 0 and finite offset < 1", ErrAsinhParams, beta, offset)
	}
	return nil
}

// AsinhRGB applies Siril's default RGB asinh stretch
// (`asinh beta [offset]`): the luminance
//
//	lum = 0.3333*r' + 0.3333*g' + 0.3333*b'    (x' = max(0, (x-offset)/(1-offset)))
//
// is stretched and every channel is multiplied by the same factor
// k = asinh(beta*lum) / (lum*asinh(beta)), preserving color ratios. The
// weights are Siril's literal 0.3333f, not 1/3. Each channel is clamped to
// [0, 1] and then passed through Siril's default RGBBLEND clip mode
// (core/arithm.h rgbblend with m_CB = 1, ported operation by operation).
//
// A NaN in any channel of a pixel makes that pixel's luminance NaN, so all
// three output channels of the pixel are NaN (Siril's max() macro
// propagates it).
//
// r, g, b must have the same length (one element per pixel); a mismatch
// is an error wrapping ErrAsinhParams. beta and offset are checked as in
// Asinh (error wrapping ErrAsinhParams). On any error the pixels are
// unchanged.
func AsinhRGB(r, g, b []float32, beta, offset float32) error {
	return asinhRGB(r, g, b, beta, offset, 0.3333, 0.3333, 0.3333)
}

// AsinhRGBHuman is AsinhRGB with Siril's `-human` luminance weights
// 0.2126, 0.7152, 0.0722 (`asinh -human beta [offset]`).
func AsinhRGBHuman(r, g, b []float32, beta, offset float32) error {
	return asinhRGB(r, g, b, beta, offset, 0.2126, 0.7152, 0.0722)
}

func asinhRGB(r, g, b []float32, beta, offset, fr, fg, fb float32) error {
	if err := asinhArgs(beta, offset); err != nil {
		return err
	}
	if len(g) != len(r) || len(b) != len(r) {
		return fmt.Errorf("%w: channel lengths r=%d g=%d b=%d differ", ErrAsinhParams, len(r), len(g), len(b))
	}
	asinhBeta := asinhf(beta)
	for i := range r {
		val := [3]float32{r[i], g[i], b[i]}
		var prime, sf, tf [3]float32
		for c := 0; c < 3; c++ {
			prime[c] = sirilMax(0, (val[c]-offset)/(1-offset))
		}
		x := float32(float32(fr*prime[0])+float32(fg*prime[1])) + float32(fb*prime[2])
		var k float32
		switch {
		case x == 0:
			k = 0
		case beta == 0:
			k = 1
		default:
			k = asinhf(float32(beta*x)) / (x * asinhBeta)
		}
		for c := 0; c < 3; c++ {
			sf[c] = sirilMin(1, sirilMax(0, float32(prime[c]*k)))
			switch {
			case prime[c] == 0:
				tf[c] = 0
			case beta == 0:
				tf[c] = 1
			default:
				tf[c] = asinhf(float32(beta*prime[c])) / (prime[c] * asinhBeta)
			}
		}
		out := rgbblend(sf, tf, val, 1)
		r[i], g[i], b[i] = out[0], out[1], out[2]
	}
	return nil
}

// rgbblend is Siril's core/arithm.h rgbblend with every channel enabled,
// written out operation by operation so the float32 rounding matches.
//
// With m_CB = 1 (the only value used) two of its comparisons cannot change
// the result, so no test can pin them: tf_c = asinh(beta*x)/(x*asinh(beta))
// and sf_c = min(1, x*k) <= 1. If tfmax >= 1, A = 1 - tfmax is exact and
// d = sfmax - tfmax rounds to at most A, so tfmax + d <= 1: "> 1" is never
// true by rounding, and where it equals 1, d = A and candidate = A/d = 1.
// If tfmax < 1 every non-zero channel is above 1, so sfmax = 1, d = 1 -
// tfmax rounded the same way as candidate's numerator, and candidate = 1.
// Either way "> 1" versus ">= 1", and min(m_CB, candidate) versus
// candidate, give the same k (a search of 1.8e7 random pixels found no
// counterexample).
func rgbblend(sf, tf, old [3]float32, mCB float32) [3]float32 {
	fastmax := func(a, b float32) float32 {
		if a > b {
			return a
		}
		return b
	}
	sfmax := fastmax(fastmax(sf[0], sf[1]), sf[2])
	tfmax := fastmax(fastmax(tf[0], tf[1]), tf[2])
	d := sfmax - tfmax
	var maskCond, maskDnz float32
	if tfmax+float32(mCB*d) > 1 {
		maskCond = 1
	}
	if d != 0 {
		maskDnz = 1
	}
	full := float32(maskCond * maskDnz)
	safeD := float32(maskDnz*d) + float32((1-maskDnz)*1)
	candidate := (1 - tfmax) / safeD
	limited := fminf(mCB, candidate)
	k := float32(full*limited) + float32((1-full)*mCB)
	omk := 1 - k
	var out [3]float32
	for c := 0; c < 3; c++ {
		// channel enabled (mr = 1): 1*(omk*tf + k*sf) + (1-1)*old
		v := float32(omk*tf[c]) + float32(k*sf[c])
		out[c] = float32(1*v) + float32(0*old[c])
	}
	return out
}

// NormalizeChannels scales all channels to [0, 1] using a single
// shared min/max across every channel. This preserves the relative
// intensity between channels, which matters for stretches that
// operate on color ratios (linked autostretch, asinh RGB).
//
// The range is that of the finite pixels: NaN and ±Inf are ignored when
// finding it. When the data is rescaled, each finite pixel becomes
// float32((v - min) * (1 / (max - min))) in float32 arithmetic, NaN
// becomes 0, +Inf becomes 1 and -Inf 0. If max - min overflows float32
// (a finite range wider than MaxFloat32) the scaling is done in float64
// instead, float32((v - min) / (max - min)), so no finite pixel is lost.
// If the finite data is already within [0, 1], all finite pixels are
// equal, or there are no finite pixels, the channels are left unchanged
// (NaN and ±Inf included). Empty channels are allowed.
func NormalizeChannels(channels ...[]float32) {
	found := false
	var mn, mx float32
	for _, ch := range channels {
		for _, v := range ch {
			if v != v || math.IsInf(float64(v), 0) {
				continue
			}
			if !found {
				mn, mx = v, v
				found = true
				continue
			}
			if v < mn {
				mn = v
			}
			if v > mx {
				mx = v
			}
		}
	}
	if !found {
		return
	}
	rng := mx - mn
	if rng <= 0 {
		return
	}
	// Already in [0, 1] — skip the work.
	if mn >= 0 && mx <= 1 {
		return
	}
	wide := math.IsInf(float64(rng), 1)
	inv := 1.0 / rng
	for _, ch := range channels {
		for i, v := range ch {
			switch {
			case v != v: // NaN
				ch[i] = 0
			case math.IsInf(float64(v), 1):
				ch[i] = 1
			case math.IsInf(float64(v), -1):
				ch[i] = 0
			case wide:
				ch[i] = float32((float64(v) - float64(mn)) / (float64(mx) - float64(mn)))
			default:
				ch[i] = float32((v - mn) * inv)
			}
		}
	}
}

// asinhf is float32 asinh. glibc's asinhf agreed with the rounded double
// asinh on all 2e6 probe arguments (testdata/sref/libm_probe.c) and on
// every golden pixel.
func asinhf(x float32) float32 { return float32(math.Asinh(float64(x))) }

func logf(x float32) float32    { return float32(math.Log(float64(x))) }
func expf(x float32) float32    { return float32(math.Exp(float64(x))) }
func powf(x, y float32) float32 { return float32(math.Pow(float64(x), float64(y))) }

// sirilMax / sirilMin are Siril's max()/min() macros (core/siril.h:47):
// a > b ? a : b. Unlike fmaxf they propagate a NaN second operand.
func sirilMax(a, b float32) float32 {
	if a > b {
		return a
	}
	return b
}

func sirilMin(a, b float32) float32 {
	if a < b {
		return a
	}
	return b
}

// fmaxf / fminf follow C99: a NaN operand is ignored.
func fmaxf(a, b float32) float32 {
	if a != a {
		return b
	}
	if b != b || a > b {
		return a
	}
	return b
}

func fminf(a, b float32) float32 {
	if a != a {
		return b
	}
	if b != b || a < b {
		return a
	}
	return b
}
