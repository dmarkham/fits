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

// QuantizeResult holds the output of a QuantizeFloat32 / QuantizeFloat64
// call when the quantizer succeeded.
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
//	row           if > 0 and ditherMethod > 0, used to seed the
//	              Park-Miller PRNG offset for dither. Matches
//	              `(row - 1) % NRandom` in cfitsio.
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
// This function is intended to be bit-exact against cfitsio. See
// testdata/cref for the validation harness.
func QuantizeFloat32(row int64, fdata []float32, nxpix, nypix int64, nullcheck bool,
	inNullValue float32, qlevel float32, ditherMethod DitherMethod, idata []int32) (QuantizeResult, bool) {

	var res QuantizeResult
	nx := nxpix * nypix
	if nx <= 1 {
		res.BScale = 1
		res.BZero = 0
		return res, false
	}

	var ngood int64
	var minval, maxval float32
	var stdev float64
	var delta float64

	if qlevel >= 0 {
		// Positive qlevel: estimate noise, delta = stdev / qlevel.
		g, mn, mx, n2, n3, n5 := fnNoise5Float32(fdata, nxpix, nypix, nullcheck, inNullValue)
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
		g, mn, mx := fnNoise3Float32MinMax(fdata, nxpix, nypix, nullcheck, inNullValue)
		ngood = g
		minval = mn
		maxval = mx
	}

	// Range check: cannot quantize if the span is too wide for int32.
	if (float64(maxval)-float64(minval))/delta > 2.0*2147483647.0-float64(nReservedValues) {
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
			zeropt = float64(minval) - delta*(float64(NullValueInt32)+float64(nReservedValues))
		} else if (float64(maxval)-float64(minval))/delta < 2147483647.0-float64(nReservedValues) {
			zeropt = float64(minval)
			// Fudge zeropt to be an integer multiple of delta.
			// This makes fpack/funpack idempotent.
			iqfactor := int64(zeropt/delta + 0.5)
			zeropt = float64(iqfactor) * delta
		} else {
			// Range almost fills int32 — center around zero.
			zeropt = (float64(minval) + float64(maxval)) / 2
		}
	} else {
		// Nulls present: shift so null sentinels are well clear of
		// the data range.
		zeropt = float64(minval) - delta*(float64(NullValueInt32)+float64(nReservedValues))
	}

	// Fill idata. Four distinct branches on (has_nulls, dithered).
	if ngood == nx {
		if row > 0 {
			for i := int64(0); i < nx; i++ {
				if ditherMethod == SubtractiveDither2 && fdata[i] == 0.0 {
					idata[i] = ZeroValue
				} else {
					idata[i] = nint((float64(fdata[i])-zeropt)/delta + float64(randoms[nextrand]) - 0.5)
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
				idata[i] = nint((float64(fdata[i]) - zeropt) / delta)
			}
		}
	} else {
		if row > 0 {
			for i := int64(0); i < nx; i++ {
				if fdata[i] != inNullValue {
					if ditherMethod == SubtractiveDither2 && fdata[i] == 0.0 {
						idata[i] = ZeroValue
					} else {
						idata[i] = nint((float64(fdata[i])-zeropt)/delta + float64(randoms[nextrand]) - 0.5)
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
					idata[i] = nint((float64(fdata[i]) - zeropt) / delta)
				} else {
					idata[i] = NullValueInt32
				}
			}
		}
	}

	// Final min/max as reported to the caller (pre-sentinel,
	// matches cfitsio's iminval/imaxval).
	res.IMinVal = nint((float64(minval) - zeropt) / delta)
	res.IMaxVal = nint((float64(maxval) - zeropt) / delta)
	res.BScale = delta
	res.BZero = zeropt
	return res, true
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
// The two-phase "NaN -> sentinel -> NullValueInt32" approach mirrors
// cfitsio's writer (reference/cfitsio/imcompress.c:2729
// imcomp_convert_tile_tfloat, which feeds FLOATNULLVALUE into
// fits_quantize_float).
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

// SubstituteNaN64 is the float64 variant of SubstituteNaN32.
func SubstituteNaN64(data []float64, sentinel float64) int {
	n := 0
	for i, v := range data {
		if v != v {
			data[i] = sentinel
			n++
		}
	}
	return n
}

// Float32NanInt32Bits is the int32 value whose bit pattern equals an
// IEEE 754 float32 NaN (specifically the one cfitsio uses: all-ones,
// i.e. int32(-1)). This is used by the NO_QUANTIZE lossless path,
// where the float tile is gzipped as raw bytes and null positions
// must be marked as NaN in the byte stream.
const Float32NanInt32Bits int32 = -1

// Float64NanInt64Bits is the int64 value whose bit pattern equals an
// IEEE 754 float64 NaN (all-ones, i.e. int64(-1)).
const Float64NanInt64Bits int64 = -1

// ReplaceSentinelWithNaNBits walks a float32 byte buffer and, wherever
// the 4-byte value equals the sentinel, overwrites those bytes with
// the NaN bit pattern. For the lossless NO_QUANTIZE write path where
// the float tile is gzipped as raw bytes. Ports imcomp_float2nan from
// reference/cfitsio/imcompress.c:7944.
//
// buf is interpreted as a host-order float32 array (same layout
// cfitsio uses — it compresses in host byte order and lets the
// reader byte-swap on decode).
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

