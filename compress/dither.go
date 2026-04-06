package compress

import "sync"

// Float-quantization dither support for FITS tile compression.
//
// Reference: cfitsio/imcompress.c fits_init_randoms() and
// unquantize_i*r4/r8. The algorithm generates a deterministic pseudo-
// random table using a Park-Miller multiplicative congruential
// generator, then uses it to drive a per-tile dither offset that is
// subtracted from each quantized integer pixel before scaling back to
// float.
//
// The PRNG is:
//
//	seed_{n+1} = (16807 * seed_n) mod 2147483647
//
// starting from seed = 1, run NRandom = 10000 times. The resulting
// floats (seed / m) form a lookup table that is initialized once and
// reused across every tile.
//
// Sanity invariant: after 10000 iterations seed == 1043618065. The
// FitsRandoms() function verifies this on first call and panics if not,
// catching any porting error.
//
// Per-tile usage matches cfitsio's unquantize path:
//
//	iseed    = (tileRow - 1) % NRandom
//	nextrand = int(randoms[iseed] * 500)
//	for each pixel:
//	    output = (input - randoms[nextrand] + 0.5) * scale + zero
//	    nextrand++
//	    if nextrand == NRandom:
//	        iseed = (iseed + 1) % NRandom
//	        nextrand = int(randoms[iseed] * 500)
//
// where tileRow is the 1-based tile index adjusted by ZDITHER0:
//
//	tileRow = nrow + ZDITHER0 - 1
//
// (nrow is the 1-based tile number; ZDITHER0 is the per-image dither
// seed keyword, typically 1).

// NRandom is the size of the PRNG lookup table. cfitsio's
// fitsio2.h defines this as 10000.
const NRandom = 10000

// ZeroValue is the sentinel integer stored in place of a true 0.0 float
// when SUBTRACTIVE_DITHER_2 is in effect. The dequantizer maps it back
// to exact 0.0.
const ZeroValue int32 = -2147483646

// DitherMethod identifies the quantization / dither convention used
// when a float image was compressed to integer coefficients.
type DitherMethod int

const (
	// NoDither — plain quantize, no random offset applied.
	NoDither DitherMethod = 0
	// SubtractiveDither1 — subtract a random [0,1) offset on encode;
	// add it back on decode. Default for RICE_1 float tiles.
	SubtractiveDither1 DitherMethod = 1
	// SubtractiveDither2 — same as SubtractiveDither1 but preserves
	// exact-zero input pixels (encoded as ZeroValue, decoded as 0.0).
	SubtractiveDither2 DitherMethod = 2
)

// ParseDitherMethod maps the ZQUANTIZ keyword string to a DitherMethod.
// Unknown or empty values default to NoDither.
func ParseDitherMethod(s string) DitherMethod {
	switch s {
	case "SUBTRACTIVE_DITHER_1":
		return SubtractiveDither1
	case "SUBTRACTIVE_DITHER_2":
		return SubtractiveDither2
	}
	return NoDither
}

// fitsRandomsOnce guards the one-time generation of the PRNG table.
// Using sync.Once avoids a data race if two goroutines call
// FitsRandoms concurrently on first use.
var (
	fitsRandoms     []float32
	fitsRandomsOnce sync.Once
)

// FitsRandoms returns the Park-Miller random lookup table used by the
// dither algorithm. On first call it generates the sequence and verifies
// the cfitsio invariant that the internal seed equals 1043618065 after
// exactly NRandom iterations. Thread-safe via sync.Once.
func FitsRandoms() []float32 {
	fitsRandomsOnce.Do(func() {
		const a = 16807.0
		const m = 2147483647.0
		out := make([]float32, NRandom)
		seed := 1.0
		for i := 0; i < NRandom; i++ {
			temp := a * seed
			seed = temp - m*float64(int64(temp/m))
			out[i] = float32(seed / m)
		}
		if int(seed) != 1043618065 {
			panic("compress: fits_init_randoms invariant violated (got seed != 1043618065)")
		}
		fitsRandoms = out
	})
	return fitsRandoms
}

// Dequantize applies the per-pixel dither + scale + zero reconstruction
// to a slice of quantized integers, writing float64 values to out.
//
// tileRow is the 1-based tile number in its binary table (nrow), NOT
// yet adjusted for ZDITHER0 — this function does the adjustment
// internally to match cfitsio's convention.
//
// For NoDither this degenerates to the plain linear reconstruction
// input*scale + zero.
//
// For SubtractiveDither2, integer values equal to ZeroValue map to
// exactly 0.0 on output, preserving the encoder's zero-exact guarantee.
func Dequantize(input []int32, out []float64, scale, zero float64, method DitherMethod, tileRow, zdither0 int) {
	if method == NoDither {
		for i, v := range input {
			out[i] = float64(v)*scale + zero
		}
		return
	}
	randoms := FitsRandoms()
	row := tileRow + zdither0 - 1 // matches cfitsio (row - 1) indexing below
	iseed := (row - 1) % NRandom
	if iseed < 0 {
		iseed += NRandom
	}
	nextrand := int(randoms[iseed] * 500)
	for i, v := range input {
		if method == SubtractiveDither2 && v == ZeroValue {
			out[i] = 0.0
		} else {
			out[i] = (float64(v) - float64(randoms[nextrand]) + 0.5) * scale + zero
		}
		nextrand++
		if nextrand == NRandom {
			iseed = (iseed + 1) % NRandom
			nextrand = int(randoms[iseed] * 500)
		}
	}
}

// DequantizeFloat32 is the float32 variant, writing directly to a
// float32 slice without an intermediate float64 copy. Used by the
// FITS-side tile pipeline where the final pixel type is float32.
func DequantizeFloat32(input []int32, out []float32, scale, zero float64, method DitherMethod, tileRow, zdither0 int) {
	if method == NoDither {
		for i, v := range input {
			out[i] = float32(float64(v)*scale + zero)
		}
		return
	}
	randoms := FitsRandoms()
	row := tileRow + zdither0 - 1
	iseed := (row - 1) % NRandom
	if iseed < 0 {
		iseed += NRandom
	}
	nextrand := int(randoms[iseed] * 500)
	for i, v := range input {
		if method == SubtractiveDither2 && v == ZeroValue {
			out[i] = 0.0
		} else {
			out[i] = float32((float64(v) - float64(randoms[nextrand]) + 0.5) * scale + zero)
		}
		nextrand++
		if nextrand == NRandom {
			iseed = (iseed + 1) % NRandom
			nextrand = int(randoms[iseed] * 500)
		}
	}
}
