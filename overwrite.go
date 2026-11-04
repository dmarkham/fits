package fits

import (
	"fmt"
	"io"
	"os"

	"github.com/dmarkham/fits/header"
)

// OverwritePixels writes new pixel values over the existing data allocation
// of an ImageHDU. It requires ModeEdit; ModeRead returns ErrReadOnly.
//
// Values are stored with the HDU's own BSCALE/BZERO exactly as WriteImage
// stores them (stored = (value - BZERO) / BSCALE, rounded half away from
// zero for integer BITPIX; exact for integer data with BSCALE = 1 and
// integral BZERO, which covers the int8/uint16/uint32/uint64
// conventions), so overwriting an HDU written by WriteImage[T] with the
// same values reproduces its bytes; this is checked against cfitsio for
// every case in testdata/images/cases.json and every value of
// golden/write_overflow.json.
//
// The call fails with:
//
//   - ErrReadOnly if the parent file is not in ModeEdit.
//   - ErrShapeMismatch if len(data) != NAXIS1 × NAXIS2 × ... for the HDU.
//   - A wrapped ErrInvalidScaling if the HDU's scaling cards (including
//     an in-memory edit, as it would be written) are refused by the
//     scaling-card rule every read applies (scalingCards: e.g. a BSCALE
//     that is not a FITS number, a HIERARCH BSCALE), if an in-memory
//     BSCALE/BZERO/BLANK card cannot be encoded (so it would never reach
//     the disk), or if BSCALE is 0: every value would be lost. With
//     duplicate cards the last applies, as for every read.
//   - A wrapped ErrNumOverflow if a value cannot be stored under the HDU's
//     BITPIX and scaling (the WriteImage rules; e.g. every nonzero value
//     under BSCALE 2^-1074): nothing is written, the file is unchanged.
//   - A wrapped ErrTypeMismatch if T is not the type WriteImage would use
//     for the HDU's BITPIX (uint8/int8 for 8, int16/uint16 for 16,
//     int32/uint32 for 32, int64/uint64 for 64, float32 for -32, float64
//     for -64).
//
// OverwritePixels does NOT resize the HDU; resize operations must go
// through EditFile.
//
// The write is crash safe: the old data unit is first saved to the undo
// journal (<name>.journal, see journal.go), so it costs one extra copy of
// the data unit. If the process dies before OverwritePixels returns, the
// next OpenForEdit or EditFile puts the old bytes back, and until then a
// read-only Open refuses with ErrPendingJournal instead of reading a data
// unit that is partly new and partly old. When the data unit on disk lacks
// its final padding the write grows the file, and the rollback also cuts it
// back to its old size (TestInPlaceWriteCrashAtEveryStep cases
// OverwritePixels and OverwritePixelsUnpadded). A real write or fsync error
// is rolled back the same way before it is returned
// (TestInPlaceWriteIOErrorRollsBack).
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

	// A scaling card of an edited header that cannot be encoded (a value
	// changed to NaN in the Cards() slice) never reaches the disk, so no
	// scaling would make the stored data read back as written.
	if p := h.rec.parsed; p != nil {
		for _, c := range p.Cards() {
			if !isScalingKey(c.Key) {
				continue
			}
			if _, err := header.EncodeCardWithContinuation(c); err != nil {
				return fmt.Errorf("%w: OverwritePixels: the %s card cannot be written: %v", ErrInvalidScaling, c.Key, err)
			}
		}
	}
	s, _, _, err := h.rec.imageKeywords()
	if err != nil {
		return fmt.Errorf("fits: OverwritePixels: %w", err)
	}
	if s.bscale == 0 {
		return fmt.Errorf("%w: OverwritePixels: the HDU's BSCALE is 0, so every value would be stored as BZERO", ErrInvalidScaling)
	}
	enc, err := encodeOrOverflow(h.bitpixType(), s, data)
	if err != nil {
		return fmt.Errorf("fits: OverwritePixels: %w", err)
	}
	rec := h.rec
	return f.journaledWrite([]undoRegion{{rec.dataStart, rec.paddedEnd - rec.dataStart}}, true, func(*os.File) error {
		if err := ioFault("write"); err != nil {
			return err
		}
		if _, err := f.bw.Seek(rec.dataStart, io.SeekStart); err != nil {
			return err
		}
		if err := f.bw.WriteRange(enc); err != nil {
			return err
		}
		if err := faultPoint("pixels-written"); err != nil {
			return err
		}
		// The rest of the last block is zero padding (FITS 4.0, 3.3.2),
		// also when the file lacked it or it held other bytes.
		if cur := f.bw.Pos(); cur < rec.paddedEnd {
			if err := f.bw.WriteRange(make([]byte, rec.paddedEnd-cur)); err != nil {
				return err
			}
		}
		return nil
	})
}

// AppendImage appends a new image HDU to a ModeEdit or ModeCreate file.
// It is a thin wrapper over WriteImage. hdr may be nil.
func AppendImage[T Numeric](f *File, hdr *header.Header, shape []int64, data []T) (*ImageHDU, error) {
	return WriteImage(f, hdr, shape, data)
}
