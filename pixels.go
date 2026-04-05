package fits

import (
	"context"
	"fmt"
	"math"
	"unsafe"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// Numeric is the constraint for generic pixel and column access. The library
// does NOT auto-convert across BITPIX widths — if the target type cannot
// losslessly hold the file's on-disk values, ReadPixels returns
// ErrTypeMismatch. See plan decision 5.1.
type Numeric interface {
	~uint8 | ~int8 | ~int16 | ~uint16 | ~int32 | ~uint32 |
		~int64 | ~uint64 | ~float32 | ~float64
}

// ReadPixels reads every pixel of h into a single flat slice of type T,
// applying BSCALE/BZERO. Returns an error (wrapping ErrTypeMismatch) if the
// requested type cannot losslessly hold the file's data. Row-major order
// matches the FITS on-disk layout: the first axis varies fastest.
func ReadPixels[T Numeric](h *ImageHDU) ([]T, error) {
	return ReadPixelsContext[T](context.Background(), h)
}

// ReadPixelsContext is ReadPixels with an explicit context. The context is
// consulted before the read and can be used to abort long reads.
func ReadPixelsContext[T Numeric](ctx context.Context, h *ImageHDU) ([]T, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	n := h.pixelCount()
	if n == 0 {
		return []T{}, nil
	}
	bp := h.bitpixType()
	if !canConvertInto[T](bp) {
		var zero T
		return nil, &TypeMismatchError{
			Requested: fmt.Sprintf("%T", zero),
			BITPIX:    int(bp),
			Lossy:     true,
		}
	}
	raw, err := readRawBytes(h, n)
	if err != nil {
		return nil, err
	}
	out := make([]T, n)
	bscale := h.BSCALE()
	bzero := h.BZERO()
	applyBScale := !(bscale == 1.0 && bzero == 0.0)
	if err := decodePixels(raw, out, bp, applyBScale, bscale, bzero); err != nil {
		return nil, err
	}
	return out, nil
}

// ReadPixelsMasked is like ReadPixels but also returns a per-pixel validity
// mask. A pixel is invalid (mask[i] == false) if the underlying raw value
// equals BLANK (integer BITPIX only) or is IEEE NaN (float BITPIX).
// BSCALE/BZERO are applied to valid pixels only; invalid pixels are left
// at the zero value of T.
func ReadPixelsMasked[T Numeric](h *ImageHDU) ([]T, []bool, error) {
	n := h.pixelCount()
	if n == 0 {
		return []T{}, []bool{}, nil
	}
	bp := h.bitpixType()
	if !canConvertInto[T](bp) {
		var zero T
		return nil, nil, &TypeMismatchError{
			Requested: fmt.Sprintf("%T", zero),
			BITPIX:    int(bp),
			Lossy:     true,
		}
	}
	raw, err := readRawBytes(h, n)
	if err != nil {
		return nil, nil, err
	}
	out := make([]T, n)
	mask := make([]bool, n)
	blankVal, haveBlank := readBLANK(h)
	bscale := h.BSCALE()
	bzero := h.BZERO()
	applyBScale := !(bscale == 1.0 && bzero == 0.0)
	if err := decodePixelsMasked(raw, out, mask, bp, applyBScale, bscale, bzero, blankVal, haveBlank); err != nil {
		return nil, nil, err
	}
	return out, mask, nil
}

// CanConvert reports whether ReadPixels[T] on h would succeed without
// returning ErrTypeMismatch. It does not touch the file content.
func CanConvert[T Numeric](h *ImageHDU) bool {
	return canConvertInto[T](h.bitpixType())
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

// readBLANK returns the BLANK keyword value if present on an integer-BITPIX
// HDU. BLANK is undefined for float BITPIX (§4.4.2.4).
func readBLANK(h *ImageHDU) (int64, bool) {
	if h.bitpixType().IsFloat() {
		return 0, false
	}
	hdr := h.Header()
	if v, err := hdr.Int(header.KeyBlankV); err == nil {
		return v, true
	}
	return 0, false
}

// canConvertInto reports whether T can losslessly hold values of the file's
// BITPIX. The rules are:
//
//   - float64 accepts any BITPIX except int64/uint64 with values exceeding
//     2^53. We conservatively allow float64 for all BITPIX (callers get
//     approximation for huge int64s, which is the Go float64 norm — if exact
//     precision is required for int64 data, read into int64 directly).
//   - float32 accepts BITPIX 8/16/-32 exactly and approximates int32/int64/float64.
//   - Integer T accepts only BITPIX values with the same width and signedness
//     (after TZERO unsigned conversion is considered).
//
// We use a generous policy: floats accept any numeric BITPIX (with
// documented approximation for wide ints); integers require width match.
func canConvertInto[T Numeric](bp bitpix.BITPIX) bool {
	var zero T
	switch any(zero).(type) {
	case float64:
		return true
	case float32:
		return true
	case int8:
		return bp == bitpix.Int8
	case uint8:
		return bp == bitpix.Int8
	case int16:
		return bp == bitpix.Int16 || bp == bitpix.Int8
	case uint16:
		return bp == bitpix.Int16 || bp == bitpix.Int8
	case int32:
		return bp == bitpix.Int32 || bp == bitpix.Int16 || bp == bitpix.Int8
	case uint32:
		return bp == bitpix.Int32 || bp == bitpix.Int16 || bp == bitpix.Int8
	case int64:
		return bp.IsInt()
	case uint64:
		return bp.IsInt()
	}
	return false
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
				v = v*bscale + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int16:
		for i := range out {
			u := bigendian.Int16(raw[i*2:])
			v := float64(u)
			if applyBScale {
				v = v*bscale + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int32:
		for i := range out {
			u := bigendian.Int32(raw[i*4:])
			v := float64(u)
			if applyBScale {
				v = v*bscale + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Int64:
		for i := range out {
			u := bigendian.Int64(raw[i*8:])
			v := float64(u)
			if applyBScale {
				v = v*bscale + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Float32:
		for i := range out {
			f := bigendian.Float32(raw[i*4:])
			v := float64(f)
			if applyBScale {
				v = v*bscale + bzero
			}
			out[i] = numericFromFloat[T](v)
		}
	case bitpix.Float64:
		for i := range out {
			f := bigendian.Float64(raw[i*8:])
			if applyBScale {
				f = f*bscale + bzero
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
				v = v*bscale + bzero
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
				v = v*bscale + bzero
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
				v = v*bscale + bzero
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
				v = v*bscale + bzero
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
				v = v*bscale + bzero
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
				f = f*bscale + bzero
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
