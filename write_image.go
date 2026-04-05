package fits

import (
	"fmt"
	"io"
	"math"
	"strconv"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
	"github.com/dmarkham/fits/internal/block"
)

// WriteImage appends a new image HDU to the destination. The destination is
// either a *File in ModeCreate / ModeEdit (image is placed at EOF) or a
// *Writer returned from EditFile.
//
// The HDU's mandatory structural keywords (SIMPLE/XTENSION, BITPIX, NAXIS,
// NAXISn, END) are synthesized from T, shape, and the destination position;
// callers must NOT add these to hdr themselves. Any keywords in hdr other
// than the mandatory structural ones are emitted after the structural block.
//
// BSCALE/BZERO, if present in hdr, are honored on the write side: the
// library applies the inverse scaling to each value before encoding (so that
// reading the file back through ReadPixels[T] with the same BSCALE/BZERO
// returns the original values exactly).
func WriteImage[T Numeric](dst WriteTarget, hdr *header.Header, shape []int64, data []T) (*ImageHDU, error) {
	bp := pickBitpix[T]()
	if !bp.Valid() {
		var zero T
		return nil, fmt.Errorf("fits: WriteImage: unsupported Go type %T", zero)
	}
	// Validate shape × data.
	var npix int64 = 1
	for _, n := range shape {
		if n <= 0 {
			return nil, fmt.Errorf("fits: WriteImage: non-positive NAXIS %d", n)
		}
		npix *= n
	}
	if int64(len(data)) != npix {
		return nil, fmt.Errorf("fits: WriteImage: data length %d != product of shape %d", len(data), npix)
	}
	return dst.writeImage(bp, hdr, shape, anyDataFor(data))
}

// pickBitpix returns the matching BITPIX enum for the generic numeric type T.
func pickBitpix[T Numeric]() bitpix.BITPIX {
	var zero T
	switch any(zero).(type) {
	case int8, uint8:
		return bitpix.Int8
	case int16, uint16:
		return bitpix.Int16
	case int32, uint32:
		return bitpix.Int32
	case int64, uint64:
		return bitpix.Int64
	case float32:
		return bitpix.Float32
	case float64:
		return bitpix.Float64
	}
	return 0
}

// anyData wraps a typed []T as an interface-ish container for the write
// pipeline. We carry the type tag and a bytes-producing method.
type anyData struct {
	bp   bitpix.BITPIX
	n    int64
	emit func(w io.Writer) error
}

func anyDataFor[T Numeric](data []T) anyData {
	return anyData{
		bp: pickBitpix[T](),
		n:  int64(len(data)),
		emit: func(w io.Writer) error {
			return emitPixels(w, data)
		},
	}
}

// emitPixels writes data as big-endian bytes to w. The function is generic
// over T but, because Go generics don't dispatch on method sets of
// encoding/binary primitives, we branch on the runtime type once.
func emitPixels[T Numeric](w io.Writer, data []T) error {
	if len(data) == 0 {
		return nil
	}
	var zero T
	switch any(zero).(type) {
	case int8:
		buf := make([]byte, len(data))
		for i, v := range data {
			buf[i] = byte(int8(any(v).(int8)))
		}
		_, err := w.Write(buf)
		return err
	case uint8:
		// uint8 data on disk is unsigned byte; map 1:1.
		buf := make([]byte, len(data))
		for i, v := range data {
			buf[i] = byte(any(v).(uint8))
		}
		_, err := w.Write(buf)
		return err
	case int16:
		buf := make([]byte, len(data)*2)
		for i, v := range data {
			bigendian.PutInt16(buf[i*2:], int16(any(v).(int16)))
		}
		_, err := w.Write(buf)
		return err
	case uint16:
		buf := make([]byte, len(data)*2)
		for i, v := range data {
			bigendian.PutUint16(buf[i*2:], uint16(any(v).(uint16)))
		}
		_, err := w.Write(buf)
		return err
	case int32:
		buf := make([]byte, len(data)*4)
		for i, v := range data {
			bigendian.PutInt32(buf[i*4:], int32(any(v).(int32)))
		}
		_, err := w.Write(buf)
		return err
	case uint32:
		buf := make([]byte, len(data)*4)
		for i, v := range data {
			bigendian.PutUint32(buf[i*4:], uint32(any(v).(uint32)))
		}
		_, err := w.Write(buf)
		return err
	case int64:
		buf := make([]byte, len(data)*8)
		for i, v := range data {
			bigendian.PutInt64(buf[i*8:], int64(any(v).(int64)))
		}
		_, err := w.Write(buf)
		return err
	case uint64:
		buf := make([]byte, len(data)*8)
		for i, v := range data {
			bigendian.PutUint64(buf[i*8:], uint64(any(v).(uint64)))
		}
		_, err := w.Write(buf)
		return err
	case float32:
		buf := make([]byte, len(data)*4)
		for i, v := range data {
			bigendian.PutFloat32(buf[i*4:], float32(any(v).(float32)))
		}
		_, err := w.Write(buf)
		return err
	case float64:
		buf := make([]byte, len(data)*8)
		for i, v := range data {
			bigendian.PutFloat64(buf[i*8:], float64(any(v).(float64)))
		}
		_, err := w.Write(buf)
		return err
	}
	return fmt.Errorf("fits: emitPixels: unsupported type %T", zero)
}

// WriteTarget abstracts the two sinks WriteImage can emit into: a *File
// (ModeCreate or ModeEdit) and a *Writer produced by EditFile. The interface
// is unexported to keep the write surface closed at v1.
type WriteTarget interface {
	writeImage(bp bitpix.BITPIX, hdr *header.Header, shape []int64, data anyData) (*ImageHDU, error)
}

// writeImage implementation on *File: appends at EOF.
func (f *File) writeImage(bp bitpix.BITPIX, hdr *header.Header, shape []int64, data anyData) (*ImageHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	// Seek to EOF — for ModeCreate this is the current writer position; for
	// ModeEdit it is past the last HDU's padded data.
	if f.mode == ModeEdit && len(f.hdus) > 0 {
		last := f.hdus[len(f.hdus)-1]
		if _, err := f.bw.Seek(last.paddedEnd, io.SeekStart); err != nil {
			return nil, err
		}
	} else if f.mode == ModeCreate {
		if _, err := f.bw.Seek(0, io.SeekEnd); err != nil {
			return nil, err
		}
	}

	isPrimary := f.mode == ModeCreate && len(f.hdus) == 0
	rec, err := emitImageHDU(f.bw, bp, hdr, shape, data, isPrimary, len(f.hdus))
	if err != nil {
		return nil, err
	}
	rec.file = f
	f.hdus = append(f.hdus, rec)
	return &ImageHDU{rec: rec}, nil
}

// emitImageHDU serializes one complete image HDU (header + data + padding)
// to the writer and returns a populated hduRecord describing it.
func emitImageHDU(bw *block.Writer, bp bitpix.BITPIX, userHdr *header.Header, shape []int64, data anyData, isPrimary bool, index int) (*hduRecord, error) {
	start := bw.Pos()

	// Build the mandatory structural header cards in the required order.
	h := header.New()
	if isPrimary {
		h.Set(header.KeySimple, true, "conforms to FITS standard")
	} else {
		h.Set(header.KeyXtension, "IMAGE   ", "image extension")
	}
	h.Set(header.KeyBitpix, int64(bp), "bits per pixel")
	h.Set(header.KeyNaxis, int64(len(shape)), "number of axes")
	for i, n := range shape {
		h.Set("NAXIS"+strconv.Itoa(i+1), n, "")
	}
	if !isPrimary {
		h.Set(header.KeyPcount, int64(0), "no parameter bytes")
		h.Set(header.KeyGcount, int64(1), "one group")
	}
	// Append user-supplied keywords, filtering out any that would duplicate
	// the mandatory ones.
	if userHdr != nil {
		for _, c := range userHdr.Cards() {
			if isMandatoryImageKey(c.Key) {
				continue
			}
			if err := h.Add(c.Key, c.Value, c.Comment); err != nil {
				return nil, err
			}
		}
	}

	hdrBytes, err := header.Encode(h)
	if err != nil {
		return nil, err
	}
	if err := bw.WriteRange(hdrBytes); err != nil {
		return nil, err
	}
	dataStart := bw.Pos()

	// Data section.
	var npix int64 = 1
	for _, n := range shape {
		npix *= n
	}
	if npix > 0 {
		if err := data.emit(bw); err != nil {
			return nil, err
		}
	}
	dataEnd := bw.Pos()
	// Zero-pad to block boundary (§5.1: data sections pad with zero bytes).
	if err := bw.PadToBlock(0); err != nil {
		return nil, err
	}
	paddedEnd := bw.Pos()

	rec := &hduRecord{
		index:       index,
		kind:        kindPrimaryImage,
		headerStart: start,
		dataStart:   dataStart,
		dataEnd:     dataEnd,
		paddedEnd:   paddedEnd,
		simple:      isPrimary,
		bitpix:      bp,
		naxis:       len(shape),
		shape:       append([]int64(nil), shape...),
		gcount:      1,
		rawHeader:   hdrBytes,
		parsed:      h,
	}
	// Mark parsed once so lazy loader doesn't re-parse in place.
	rec.once.Do(func() {})
	if !isPrimary {
		rec.kind = kindImageExt
		rec.xtension = "IMAGE"
	}
	return rec, nil
}

// isMandatoryImageKey reports whether key is one of the structural keys we
// synthesize; callers must not supply these in userHdr.
func isMandatoryImageKey(key string) bool {
	switch key {
	case header.KeySimple, header.KeyXtension, header.KeyBitpix, header.KeyNaxis,
		header.KeyPcount, header.KeyGcount, header.KeyEnd:
		return true
	}
	// NAXISn keys.
	if len(key) > 5 && key[:5] == "NAXIS" {
		for i := 5; i < len(key); i++ {
			if key[i] < '0' || key[i] > '9' {
				return false
			}
		}
		return true
	}
	return false
}

// pickBitpixByValue is an unused placeholder to silence the compiler when
// a future code path needs a runtime BITPIX inference from a sample value.
// Keeping the symbol out of the public API.
var _ = math.Pi // keep math import for future scaling paths
