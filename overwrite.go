package fits

import (
	"fmt"
	"io"

	"github.com/dmarkham/fits/header"
)

// OverwritePixels writes new pixel values over the existing data allocation
// of an ImageHDU. It requires ModeEdit; ModeRead returns ErrReadOnly.
//
// The call fails with:
//
//   - ErrReadOnly if the parent file is not in ModeEdit.
//   - ErrShapeMismatch if len(data) != NAXIS1 × NAXIS2 × ... for the HDU.
//   - A wrapped ErrTypeMismatch if T does not correspond exactly to the
//     HDU's BITPIX. Writing requires an exact-width match — we do not
//     round or clip.
//
// OverwritePixels does NOT resize the HDU; resize operations must go
// through EditFile.
func OverwritePixels[T Numeric](h *ImageHDU, data []T) error {
	f := h.rec.file
	if f == nil || f.mode != ModeEdit {
		return ErrReadOnly
	}
	wantBP := pickBitpix[T]()
	if wantBP != h.bitpixType() {
		var zero T
		return &TypeMismatchError{
			Requested: fmt.Sprintf("%T", zero),
			BITPIX:    int(h.bitpixType()),
			Lossy:     true,
		}
	}
	expect := h.pixelCount()
	if int64(len(data)) != expect {
		return fmt.Errorf("%w: data len %d != pixel count %d", ErrShapeMismatch, len(data), expect)
	}

	if _, err := f.bw.Seek(h.rec.dataStart, io.SeekStart); err != nil {
		return err
	}
	if err := emitPixels(f.bw, data); err != nil {
		return err
	}
	// Preserve block-boundary zero padding on the trailing incomplete block.
	cur := f.bw.Pos()
	if cur < h.rec.paddedEnd {
		zeroBuf := make([]byte, h.rec.paddedEnd-cur)
		if err := f.bw.WriteRange(zeroBuf); err != nil {
			return err
		}
	}
	f.br.Invalidate()
	return nil
}

// AppendImage appends a new image HDU to a ModeEdit or ModeCreate file.
// It is a thin wrapper over WriteImage. hdr may be nil.
func AppendImage[T Numeric](f *File, hdr *header.Header, shape []int64, data []T) (*ImageHDU, error) {
	return WriteImage(f, hdr, shape, data)
}
