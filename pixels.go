package fits

import (
	"context"
	"fmt"
	"math"
	"unsafe"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// Numeric is the constraint for generic pixel and column access. For
// plain images (*ImageHDU: ReadPixels, ReadPixelsMasked, ReadSubset,
// CanConvert), which T a read accepts is decided per HDU from BITPIX,
// BSCALE and BZERO: float targets accept every image (rounded like
// cfitsio), integer targets only images whose every physical value they
// hold exactly; anything else is ErrTypeMismatch, never a silently wrapped
// value. Tile-compressed images (ReadPixelsCompressed,
// ReadPixelsCompressedMasked) follow the same rule with ZBITPIX in place
// of BITPIX.
type Numeric interface {
	~uint8 | ~int8 | ~int16 | ~uint16 | ~int32 | ~uint32 |
		~int64 | ~uint64 | ~float32 | ~float64
}

// ReadPixels reads every pixel of h into a single flat slice of type T,
// applying BSCALE/BZERO (physical = BZERO + BSCALE * stored). The order is
// the FITS on-disk order: the first axis (NAXIS1) varies fastest.
//
// Which T are accepted (the same rule as CanConvert):
//
//   - float32 and float64 accept every BITPIX. Values equal cfitsio's
//     TFLOAT / TDOUBLE reads bit for bit (one rounding to the target, so
//     BITPIX 64 into float32/float64, and huge BSCALE results, are rounded
//     approximations exactly as in cfitsio). One documented difference:
//     a value beyond the float32 range reads as ±Inf into float32, where
//     cfitsio returns ±FLT_MAX with NUM_OVERFLOW.
//   - An integer T is accepted only for integer BITPIX with integral BSCALE
//     and BZERO, and only when T's range holds the physical value of every
//     possible stored value (raw range x BSCALE + BZERO). Accepted reads
//     are exact (computed in 128-bit integer arithmetic, never through
//     float64), so int64 data and the uint64 BZERO = 2^63 convention read
//     exactly. Otherwise ReadPixels returns a *TypeMismatchError (wrapping
//     ErrTypeMismatch) instead of a wrapped or truncated value: for example
//     ReadPixels[int16] on unsigned 16-bit data (BZERO = 32768) is refused;
//     read it as uint16, int32, int64 or a float type.
//
// BSCALE, BZERO and BLANK are taken from the header by the scaling-card
// rule every image read and write of this package applies (scalingCards in
// scaling.go; ImageHDU.BSCALE): the last card of each wins, as cfitsio and
// astropy apply it, and a card the two references read differently (a
// value that is not a FITS number, a lowercase or misplaced '=' card, a
// HIERARCH BSCALE/BZERO/BLANK, a BLANK outside the BITPIX range, ...)
// makes the read fail with a wrapped ErrInvalidScaling, the same before or
// after the header is parsed (TestScalingKeywordsRule).
//
// BLANK is not applied (BLANK pixels come back as their scaled stored
// value, as cfitsio's reads without null checking return them); use
// ReadPixelsMasked to tell them apart.
func ReadPixels[T Numeric](h *ImageHDU) ([]T, error) {
	return ReadPixelsContext[T](context.Background(), h)
}

// ReadPixelsContext is ReadPixels with an explicit context. The context is
// checked once, before the read starts (a cancelled context returns its
// error and reads nothing); a read already in progress is not
// interrupted. TestContextFunctions.
func ReadPixelsContext[T Numeric](ctx context.Context, h *ImageHDU) ([]T, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	conv, err := pixelConverter[T](h)
	if err != nil {
		return nil, err
	}
	n := h.pixelCount()
	if n == 0 {
		return []T{}, nil
	}
	raw, err := readRawBytes(h, n)
	if err != nil {
		return nil, err
	}
	out := make([]T, n)
	conv.all(raw, out)
	return out, nil
}

// pixelConverter checks that T can hold h's physical values and returns
// the converter for them.
func pixelConverter[T Numeric](h *ImageHDU) (*pixConv[T], error) {
	bp := h.bitpixType()
	s, _, _, err := h.rec.imageKeywords()
	if err != nil {
		return nil, fmt.Errorf("fits: HDU %d: %w", h.rec.index, err)
	}
	if !canHold[T](bp, s) {
		var zero T
		return nil, &TypeMismatchError{Requested: fmt.Sprintf("%T", zero), BITPIX: int(bp), Lossy: true,
			BSCALE: s.bscale, BZERO: s.bzero}
	}
	return newPixConv[T](bp, s), nil
}

// ReadPixelsMasked is like ReadPixels but also returns a per-pixel validity
// mask. A pixel is invalid (mask[i] == false) if its stored value equals
// BLANK (integer BITPIX; BLANK is compared with the stored, unscaled value,
// as cfitsio compares it: for BITPIX 32 as a 32-bit value, and a BLANK of
// 1234554321, cfitsio's NULL_UNDEFINED, marks nothing) or is an IEEE NaN
// (float BITPIX). Invalid pixels are left at the zero
// value of T; valid pixels are converted exactly as by ReadPixels. The mask
// equals cfitsio's null flags, except that cfitsio also flags +-Inf as
// null; here an infinity is a valid value.
func ReadPixelsMasked[T Numeric](h *ImageHDU) ([]T, []bool, error) {
	conv, err := pixelConverter[T](h)
	if err != nil {
		return nil, nil, err
	}
	n := h.pixelCount()
	if n == 0 {
		return []T{}, []bool{}, nil
	}
	raw, err := readRawBytes(h, n)
	if err != nil {
		return nil, nil, err
	}
	_, blank, hasBlank, err := h.rec.imageKeywords()
	if err != nil {
		return nil, nil, err // not reached: pixelConverter read the same cards
	}
	out := make([]T, n)
	mask := make([]bool, n)
	for i := range out {
		b := raw[i*conv.size:]
		if conv.valid(b, blank, hasBlank) {
			mask[i] = true
			out[i] = conv.one(b)
		}
	}
	return out, mask, nil
}

// CanConvert reports whether ReadPixels[T] on h would succeed without
// returning ErrTypeMismatch (see ReadPixels for the rule). It reads only
// the header, not the pixels.
func CanConvert[T Numeric](h *ImageHDU) bool {
	s, _, _, err := h.rec.imageKeywords()
	return err == nil && canHold[T](h.bitpixType(), s)
}

// readRawBytes reads npix*|BITPIX|/8 bytes from the HDU data area into a
// fresh byte slice (big-endian, unscaled).
func readRawBytes(h *ImageHDU, npix int64) ([]byte, error) {
	bpb := int64(h.bitpixType().Size())
	nbytes := npix * bpb
	buf := make([]byte, nbytes)
	if err := h.rec.file.br.ReadRange(h.rec.dataStart, buf); err != nil {
		return nil, fmt.Errorf("fits: read image data at HDU %d: %w", h.rec.index, err)
	}
	return buf, nil
}

// decodePixels fills out from the big-endian raw byte slice. It applies
// BSCALE/BZERO inline when applyBScale is true.
//
// The implementation is a straight switch on BITPIX. For performance we
// could bulk-swap and re-interpret, but clarity beats branch elimination at
// this level — hot loops will optimize under the compiler's escape/bounds
// analysis.
func decodePixels[T Numeric](raw []byte, out []T, bp bitpix.BITPIX, applyBScale bool, bscale, bzero float64) error {
	switch bp {
	case bitpix.Int8:
		for i := range out {
			v := float64(raw[i]) // BITPIX=8 is unsigned byte on disk
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int16:
		for i := range out {
			u := bigendian.Int16(raw[i*2:])
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int32:
		for i := range out {
			u := bigendian.Int32(raw[i*4:])
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int64:
		for i := range out {
			u := bigendian.Int64(raw[i*8:])
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Float32:
		for i := range out {
			f := bigendian.Float32(raw[i*4:])
			v := float64(f)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Float64:
		for i := range out {
			f := bigendian.Float64(raw[i*8:])
			if applyBScale {
				f = float64(f*bscale) + bzero
			}
			out[i] = numericFromFloat[T](f)
		}
	default:
		return fmt.Errorf("fits: unsupported BITPIX %d", int(bp))
	}
	return nil
}

// decodePixelsMasked is the masked variant. Invalid pixels (BLANK match or
// NaN) are left as T's zero value and their mask entry is false.
func decodePixelsMasked[T Numeric](raw []byte, out []T, mask []bool, bp bitpix.BITPIX, applyBScale bool, bscale, bzero float64, blankVal int64, haveBlank bool) error {
	switch bp {
	case bitpix.Int8:
		for i := range out {
			u := int64(raw[i])
			if haveBlank && u == blankVal {
				continue
			}
			mask[i] = true
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int16:
		for i := range out {
			u := int64(bigendian.Int16(raw[i*2:]))
			if haveBlank && u == blankVal {
				continue
			}
			mask[i] = true
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int32:
		for i := range out {
			u := int64(bigendian.Int32(raw[i*4:]))
			if haveBlank && u == blankVal {
				continue
			}
			mask[i] = true
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int64:
		for i := range out {
			u := bigendian.Int64(raw[i*8:])
			if haveBlank && u == blankVal {
				continue
			}
			mask[i] = true
			v := float64(u)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Float32:
		for i := range out {
			f := bigendian.Float32(raw[i*4:])
			if f != f { // NaN
				continue
			}
			mask[i] = true
			v := float64(f)
			if applyBScale {
				v = float64(v*bscale) + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Float64:
		for i := range out {
			f := bigendian.Float64(raw[i*8:])
			if math.IsNaN(f) {
				continue
			}
			mask[i] = true
			if applyBScale {
				f = float64(f*bscale) + bzero
			}
			out[i] = numericFromFloat[T](f)
		}
	default:
		return fmt.Errorf("fits: unsupported BITPIX %d", int(bp))
	}
	return nil
}

// numericFromFloat casts a float64 value into the generic target type T
// without rounding (Go's float→int conversion truncates, which is the
// expected behavior for values that already represent exact integers via
// BSCALE/BZERO).
//
// We use an unsafe.Pointer switch to avoid a generic type-switch explosion.
func numericFromFloat[T Numeric](v float64) T {
	var zero T
	switch any(zero).(type) {
	case float64:
		return *(*T)(unsafe.Pointer(&v))
	case float32:
		f := float32(v)
		return *(*T)(unsafe.Pointer(&f))
	case int8:
		i := int8(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint8:
		i := uint8(v)
		return *(*T)(unsafe.Pointer(&i))
	case int16:
		i := int16(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint16:
		i := uint16(v)
		return *(*T)(unsafe.Pointer(&i))
	case int32:
		i := int32(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint32:
		i := uint32(v)
		return *(*T)(unsafe.Pointer(&i))
	case int64:
		i := int64(v)
		return *(*T)(unsafe.Pointer(&i))
	case uint64:
		i := uint64(v)
		return *(*T)(unsafe.Pointer(&i))
	}
	return zero
}
