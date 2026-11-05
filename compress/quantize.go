package compress

import (
	"math"
)

// Float quantization for FITS tile compression.
//
// Pure-Go port of cfitsio's fits_quantize_float / fits_quantize_double
// from reference/cfitsio/quantize.c:84. Produces an int32 array + a
// per-tile bscale/bzero, which can then be fed to any integer tile
// compressor (RICE_1, GZIP_1, GZIP_2) and reconstructed on the read
// side via Dequantize / DequantizeFloat32.
//
// Reference: reference/cfitsio/quantize.c:84  fits_quantize_float
//            reference/cfitsio/quantize.c:273 fits_quantize_double

// NullValueInt32 is the sentinel stored in the quantized integer stream
// in place of a null input pixel. cfitsio's quantize.c defines this as
// -2147483647.
const NullValueInt32 int32 = -2147483647

// nReservedValues is the number of reserved integer values at the
// bottom of the int32 range (starting from NullValueInt32). Quantized
// pixel values must not collide with these.
const nReservedValues = 10

// QuantizeResult holds the output of a QuantizeFloat32 call when the
// quantizer succeeded.
type QuantizeResult struct {
	BScale  float64 // FITS BSCALE: one unit in idata == BScale in fdata
	BZero   float64 // FITS BZERO: integer 0 in idata == BZero in fdata
	IMinVal int32   // minimum value actually written to idata (excluding sentinels)
	IMaxVal int32   // maximum value actually written to idata (excluding sentinels)
}

// QuantizeFloat32 quantizes a float32 tile to int32, producing a per-tile
// bscale/bzero. The function returns true if the tile was quantized and
// idata is populated; false if quantization was skipped (constant tile,
// delta overflow, or too few pixels).
//
// Parameters mirror fits_quantize_float exactly:
//
//	row           if > 0, the pixels are dithered and row seeds the
//	              position in the random table: iseed = (row-1) % NRandom,
//	              as in cfitsio. 0 = no dithering.
//	fdata         input pixels, row-major with first axis varying fastest
//	nxpix, nypix  image shape
//	nullcheck     if true, fdata[i] == inNullValue is treated as null
//	inNullValue   null sentinel (cfitsio's approach: caller pre-substitutes
//	              NaN with a chosen sentinel before calling)
//	qlevel        positive: delta = stdev / qlevel  (qlevel==0 => stdev/4)
//	              negative: delta = -qlevel         (absolute step)
//	ditherMethod  NoDither / SubtractiveDither1 / SubtractiveDither2
//	idata         caller-allocated output buffer, len(idata) >= nxpix*nypix
//
// The output bscale/bzero are written into the returned QuantizeResult.
//
// Bit-exact against libcfitsio 4.7.0's fits_quantize_float (idata,
// bscale, bzero, iminval, imaxval) on the cases in testdata/cref
// (TestQuantizeFloat32_CRef), with one deliberate difference: where
// cfitsio's arithmetic overflows it returns false (the tile is not
// quantized) instead of wrong integers. That is when zeropt/delta does not
// fit the int64 iqfactor (a tiny absolute qlevel: cfitsio decodes a tile
// of 100 at qlevel -1e-20 as -0.0922), or when a pixel's quantized value
// is outside int32 or lands on the NULL_VALUE (or, under
// SUBTRACTIVE_DITHER_2, ZERO_VALUE) sentinel: the centered branch's float
// sum minval+maxval can round so far that the maximum quantizes to 2^31,
// which wraps, or the minimum to NULL_VALUE. Goldens
// golden/img2/rice_f32_zpt_overflow_abs, rice_f32_center_overflow,
// rice_f32_null_sentinel_hit and golden/img64/rice_f64_zpt_overflow_abs
// pin cfitsio's wrong output;
// TestQuantizerOverflowStoredLosslessly the writer's lossless storage.
// Callers dithering tile number t (1-based)
// of an image with ZDITHER0 = z pass row = t + z - 1, as cfitsio's
// imcomp_convert_tile_tfloat does; row = 0 means NO_DITHER.
func QuantizeFloat32(row int64, fdata []float32, nxpix, nypix int64, nullcheck bool,
	inNullValue float32, qlevel float32, ditherMethod DitherMethod, idata []int32) (QuantizeResult, bool) {
	return quantize(row, fdata, nxpix, nypix, nullcheck, inNullValue, qlevel, ditherMethod, idata)
}

// QuantizeFloat64 is QuantizeFloat32 for float64 pixels: a port of
// cfitsio's fits_quantize_double, which is fits_quantize_float with the
// pixel type (and FnNoise5_double / FnNoise3_double) changed; qlevel is
// a float32 in both. Bit-exact against libcfitsio 4.7.0 through whole
// ZBITPIX -64 images compressed by fits_img_compress
// (compress/testdata/cref/golden/img64, TestCompressedWriterFloat64MatchesCfitsio).
func QuantizeFloat64(row int64, fdata []float64, nxpix, nypix int64, nullcheck bool,
	inNullValue float64, qlevel float32, ditherMethod DitherMethod, idata []int32) (QuantizeResult, bool) {
	return quantize(row, fdata, nxpix, nypix, nullcheck, inNullValue, qlevel, ditherMethod, idata)
}

// quantize is fits_quantize_float (T float32) / fits_quantize_double
// (T float64). Pixel arithmetic (maxval - minval, minval + maxval) is in
// T, as in C.
func quantize[T floatPix](row int64, fdata []T, nxpix, nypix int64, nullcheck bool,
	inNullValue T, qlevel float32, ditherMethod DitherMethod, idata []int32) (QuantizeResult, bool) {

	var res QuantizeResult
	nx := nxpix * nypix
	if nx <= 1 {
		res.BScale = 1
		res.BZero = 0
		return res, false
	}

	var ngood int64
	var minval, maxval T
	var stdev float64
	var delta float64

	if qlevel >= 0 {
		// Positive qlevel: estimate noise, delta = stdev / qlevel.
		g, mn, mx, n2, n3, n5 := fnNoise5(fdata, nxpix, nypix, nullcheck, inNullValue)
		ngood = g
		minval = mn
		maxval = mx

		if nullcheck && ngood == 0 {
			// All-null tile: dummy values so the downstream
			// arithmetic still produces a usable bscale/bzero
			// (though the real decision is to just fill idata
			// with NullValueInt32). cfitsio does this too.
			minval = 0
			maxval = 1
			stdev = 1
		} else {
			// stdev = min(noise2, noise3, noise5) among the
			// non-zero candidates. noise3 is the primary.
			stdev = n3
			if n2 != 0 && n2 < stdev {
				stdev = n2
			}
			if n5 != 0 && n5 < stdev {
				stdev = n5
			}
		}

		if qlevel == 0 {
			delta = stdev / 4.0
		} else {
			delta = stdev / float64(qlevel)
		}
		if delta == 0 {
			return res, false // can't quantize
		}
	} else {
		// Negative qlevel: absolute step.
		delta = -float64(qlevel)
		g, mn, mx := fnNoise3MinMax(fdata, nxpix, nypix, nullcheck, inNullValue)
		ngood = g
		minval = mn
		maxval = mx
	}

	// Range check: cannot quantize if the span is too wide for int32.
	// cfitsio subtracts the float minval/maxval in float (C float
	// arithmetic), then divides in double: span is that float result
	// (goldens range_fsub_256x1, center_fsum_64x1).
	span := float64(maxval - minval)
	if span/delta > 2.0*2147483647.0-float64(nReservedValues) {
		return res, false
	}

	// Initialize PRNG for dither if requested.
	var nextrand int
	var iseed int
	var randoms []float32
	if row > 0 {
		randoms = FitsRandoms()
		iseed = int((row - 1) % NRandom)
		if iseed < 0 {
			iseed += NRandom
		}
		nextrand = int(randoms[iseed] * 500)
	}

	var zeropt float64
	if ngood == nx {
		// No nulls: aim for all-positive integers where possible.
		if ditherMethod == SubtractiveDither2 {
			// Shift so that integer 0 on the output side maps to
			// near the NULL_VALUE sentinel, freeing the ZERO_VALUE
			// sentinel to round-trip exact zeros.
			zeropt = float64(minval) - float64(delta*(float64(NullValueInt32)+float64(nReservedValues)))
		} else if span/delta < 2147483647.0-float64(nReservedValues) {
			zeropt = float64(minval)
			// Fudge zeropt to be an integer multiple of delta.
			// This makes fpack/funpack idempotent. cfitsio converts
			// without a range check; a quotient outside int64 (a tiny
			// absolute step) gives a meaningless zero point there, so the
			// tile is not quantized instead.
			q := float64(zeropt/delta + 0.5)
			if !(q > -9.2233720368547758e18 && q < 9.2233720368547758e18) {
				return res, false
			}
			iqfactor := int64(q)
			zeropt = float64(float64(iqfactor) * delta)
		} else {
			// Range almost fills int32 — center around zero. The sum is a
			// float sum in cfitsio, as above.
			zeropt = float64(minval+maxval) / 2
		}
	} else {
		// Nulls present: shift so null sentinels are well clear of
		// the data range.
		zeropt = float64(minval) - float64(delta*(float64(NullValueInt32)+float64(nReservedValues)))
	}

	// Fill idata. Four distinct branches on (has_nulls, dithered).
	// Every quantized pixel goes through q.nint, which records a value
	// cfitsio would wrap (outside int32) or that a reader would take for a
	// sentinel (NULL_VALUE, ZERO_VALUE); such a tile is not quantized.
	q := quantChecker{sd2: ditherMethod == SubtractiveDither2 && row > 0}
	if ngood == nx {
		if row > 0 {
			for i := int64(0); i < nx; i++ {
				if ditherMethod == SubtractiveDither2 && fdata[i] == 0.0 {
					idata[i] = ZeroValue
				} else {
					idata[i] = q.nint((float64(fdata[i])-zeropt)/delta + float64(randoms[nextrand]) - 0.5)
				}
				nextrand++
				if nextrand == NRandom {
					iseed++
					if iseed == NRandom {
						iseed = 0
					}
					nextrand = int(randoms[iseed] * 500)
				}
			}
		} else {
			for i := int64(0); i < nx; i++ {
				idata[i] = q.nint((float64(fdata[i]) - zeropt) / delta)
			}
		}
	} else {
		if row > 0 {
			for i := int64(0); i < nx; i++ {
				if fdata[i] != inNullValue {
					if ditherMethod == SubtractiveDither2 && fdata[i] == 0.0 {
						idata[i] = ZeroValue
					} else {
						idata[i] = q.nint((float64(fdata[i])-zeropt)/delta + float64(randoms[nextrand]) - 0.5)
					}
				} else {
					idata[i] = NullValueInt32
				}
				nextrand++
				if nextrand == NRandom {
					iseed++
					if iseed == NRandom {
						iseed = 0
					}
					nextrand = int(randoms[iseed] * 500)
				}
			}
		} else {
			for i := int64(0); i < nx; i++ {
				if fdata[i] != inNullValue {
					idata[i] = q.nint((float64(fdata[i]) - zeropt) / delta)
				} else {
					idata[i] = NullValueInt32
				}
			}
		}
	}

	if q.bad {
		return QuantizeResult{}, false
	}

	// Final min/max as reported to the caller (pre-sentinel,
	// matches cfitsio's iminval/imaxval).
	res.IMinVal = nint((float64(minval) - zeropt) / delta)
	res.IMaxVal = nint((float64(maxval) - zeropt) / delta)
	res.BScale = delta
	res.BZero = zeropt
	return res, true
}

// quantChecker is NINT for quantized pixels, recording any result that
// would not survive the round trip: outside int32 (C's conversion is
// undefined there; x86 gives INT_MIN), equal to NULL_VALUE (a null to
// every reader once the image has ZBLANK; golden
// img2/rice_f32_null_sentinel_hit) or, under SUBTRACTIVE_DITHER_2, equal
// to ZERO_VALUE (read back as an exact 0; no input reaching it was found,
// see docs/VALIDATION.md).
type quantChecker struct{ bad, sd2 bool }

func (q *quantChecker) nint(x float64) int32 {
	if !(x > -2147483648.5 && x < 2147483647.5) {
		q.bad = true
		return 0
	}
	v := nint(x)
	if v == NullValueInt32 || (q.sd2 && v == ZeroValue) {
		q.bad = true
	}
	return v
}

// nint is the C macro NINT: nearest integer, halves round away from zero.
func nint(x float64) int32 {
	if x >= 0 {
		return int32(x + 0.5)
	}
	return int32(x - 0.5)
}

// SubstituteNaN32 replaces NaN pixels in place with the given sentinel
// value, returning the number of substitutions performed. Callers
// should invoke this before QuantizeFloat32 on any tile that might
// contain NaN, because the quantizer detects nulls by float equality
// and NaN != NaN in IEEE 754 — a raw NaN would slip through nullcheck
// and produce undefined int32 output.
//
// It substitutes NaN only. cfitsio's writer also treats ±Inf as null
// (fffr4r4's fnan test), so the fits writers no longer use this helper;
// AppendCompressedFloat32Image classifies NaN and ±Inf itself.
func SubstituteNaN32(data []float32, sentinel float32) int {
	n := 0
	for i, v := range data {
		if v != v { // NaN != NaN
			data[i] = sentinel
			n++
		}
	}
	return n
}

// ReplaceSentinelWithNaNBits walks a float32 byte buffer and, wherever
// the 4-byte value equals the sentinel, overwrites those bytes with
// the NaN bit pattern. For the lossless NO_QUANTIZE write path where
// the float tile is gzipped as raw bytes. Ports imcomp_float2nan from
// reference/cfitsio/imcompress.c:7944.
//
// buf is a []float32 in memory; the bytes the fits writer stores in
// GZIP_COMPRESSED_DATA are FITS big-endian (pinned by the
// gzip2_f32_inf_fallback_t40x5 cfitsio golden).
func ReplaceSentinelWithNaNBits(buf []float32, sentinel float32) int {
	n := 0
	for i, v := range buf {
		if v == sentinel {
			// Overwrite with NaN bit pattern (-1 as int32).
			buf[i] = math.Float32frombits(0xFFFFFFFF)
			n++
		}
	}
	return n
}
