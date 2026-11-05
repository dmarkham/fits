package fits

import (
	"fmt"
	"io"
	"strconv"

	"github.com/dmarkham/fits/header"
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
// Caller cards are copied with header.AddCard. A card read from a file is
// written as its bytes. An edited Cards() copy (Key/Type/Value/Comment
// changed) is written from its fields under the upper-cased keyword, and
// is left out like any caller card for a keyword the writer computes
// itself. An unedited lowercase-keyword (non-HIERARCH) card makes the
// writer return an error wrapping header.ErrInvalidKeyword: cfitsio and
// astropy read such a card differently.
//
// BITPIX follows T: uint8 -> 8, int16 -> 16, int32 -> 32, int64 -> 64,
// float32 -> -32, float64 -> -64. Types FITS has no BITPIX for are written
// with the standard offset conventions (FITS 4.0 §5.3, Table 11), exactly
// as cfitsio and astropy write them: int8 -> BITPIX 8 with BZERO = -128,
// uint16 -> 16 with BZERO = 32768, uint32 -> 32 with BZERO = 2147483648,
// uint64 -> 64 with BZERO = 9223372036854775808 (and BSCALE = 1), so that
// every FITS reader sees the original values.
//
// When hdr carries BSCALE and/or BZERO (a missing one defaults to 1 / 0)
// those are used instead of the convention, read from the encoded cards by
// the scaling-card rule every reader applies (scalingCards, scaling.go;
// D exponents are read, a card the rule refuses is refused here; proven
// by TestWriteImageKeywordsMatchCfitsio), and each value is stored as
// (value - BZERO) / BSCALE: for integer BITPIX rounded half away from
// zero, for float BITPIX stored as the float result. This is cfitsio's
// fits_write_img arithmetic, and the bytes written equal cfitsio's for
// the same values and cards, with one deliberate exception: with BSCALE =
// 1 and an integral BZERO, integer data is offset exactly, so for BITPIX
// 64 with an integral BZERO other than 0 and 2^63 a value that no double
// holds is stored as v - BZERO exactly, where cfitsio rounds v - BZERO
// through a double (v = 2^62 + 1, BZERO 1000: v - 1000 here,
// 4611686018427386880 in cfitsio). Reading the file back with
// ReadPixels[T] returns the original values whenever they were
// representable (not rounded).
//
// A value that cannot be stored is never clamped: WriteImage returns an
// error wrapping ErrNumOverflow, naming the first such value, and writes
// nothing. For integer BITPIX that is a stored value outside the stored
// range (BITPIX 8: 0..255, 16: int16, 32: int32, 64: int64) by more than
// cfitsio's 0.49 rounding margin, exactly where cfitsio's fits_write_img
// stores the range limit and reports NUM_OVERFLOW (412), and in two places
// where cfitsio stores a wrong value without reporting it: a BITPIX 64
// scaled value of exactly 2^63 (cfitsio's out-of-range cast stores
// INT64_MIN) and an exact int64 offset beyond the range that cfitsio's
// double arithmetic absorbs (v = -2^63 with BZERO 1000). For float BITPIX
// it is a finite value whose scaled form is infinite in the stored type
// (cfitsio stores the infinity). Proven per value against cfitsio by
// TestWriteImageOverflowMatchesCfitsio.
//
// WriteImage returns an error wrapping ErrInvalidScaling, and writes
// nothing, when hdr's BSCALE/BZERO/BLANK cards cannot store the data so
// that cfitsio reads it back: BSCALE = 0 (cfitsio's fits_set_bscale
// refuses it, ZERO_SCALE), a BSCALE, BZERO or BLANK given twice (cfitsio
// applies the last card), a value cfitsio cannot convert (a quoted string,
// a lowercase 'd' exponent, Inf, NaN, an over- or underflowing number; it
// would substitute the default), a BLANK that is not an integer,
// BLANK = 1234554321 (cfitsio's NULL_UNDEFINED), and a scaling card that
// cannot be encoded (it would never reach the disk).
//
// The written data are proven byte-identical to cfitsio's (outside the
// BITPIX 64 exception above), and read back identically by cfitsio,
// fitsverify and astropy, for every BITPIX and convention
// (testdata/images/cases.json).
func WriteImage[T Numeric](dst WriteTarget, hdr *header.Header, shape []int64, data []T) (*ImageHDU, error) {
	bp := pickBitpix[T]()
	if !bp.Valid() {
		var zero T
		return nil, fmt.Errorf("fits: WriteImage: unsupported Go type %T", zero)
	}
	// Validate shape × data. NAXIS=0 means zero-pixel primary placeholder.
	var npix int64
	if len(shape) == 0 {
		npix = 0
	} else {
		npix = 1
		for _, n := range shape {
			if n <= 0 {
				return nil, fmt.Errorf("fits: WriteImage: non-positive NAXIS %d", n)
			}
			npix *= n
		}
	}
	if int64(len(data)) != npix {
		return nil, fmt.Errorf("fits: WriteImage: %w: data length %d != product of shape %d", ErrShapeMismatch, len(data), npix)
	}
	s, cards, err := writeScaling[T](hdr)
	if err != nil {
		return nil, err
	}
	// Encode before anything is written, so that a value that cannot be
	// stored refuses the whole write.
	enc, err := encodeOrOverflow(bp, s, data)
	if err != nil {
		return nil, fmt.Errorf("fits: WriteImage: %w", err)
	}
	d := anyData{
		bp:    bp,
		n:     int64(len(data)),
		cards: cards,
		emit: func(w io.Writer) error {
			if len(enc) == 0 {
				return nil
			}
			_, err := w.Write(enc)
			return err
		},
	}
	return dst.writeImage(bp, hdr, shape, d)
}

// writeScaling picks the scaling WriteImage stores T with: the caller's
// BSCALE/BZERO when hdr has either, else the int8/unsigned offset
// convention for T (returned with the BZERO/BSCALE cards to add), else none.
//
// The caller's cards are read from their encoded bytes by the scaling-card
// rule every reader applies (scalingCards), so the values are stored with
// the scaling the written file is read with. A card the rule refuses is
// refused here (a wrapped ErrInvalidScaling), and so are a keyword given
// twice (the readers apply the last; a caller who wrote two meant one of
// them) and BSCALE = 0.
func writeScaling[T Numeric](hdr *header.Header) (scaling, []header.Card, error) {
	if hdr != nil {
		k, err := scalingCards(encodedCards(hdr, isMandatoryImageKey), pickBitpix[T](), false)
		if err != nil {
			return scaling{}, nil, err
		}
		if err := checkWriteScaling(k); err != nil {
			return scaling{}, nil, err
		}
		if k.count[keyBscale] > 0 || k.count[keyBzero] > 0 {
			return newScaling(k.bscale, k.bzero), nil, nil
		}
	}
	var z T
	var off float64
	var text, comment string
	switch any(z).(type) {
	case int8:
		off, text, comment = -128, "-128", "offset data range to that of signed byte"
	case uint16:
		off, text, comment = 32768, "32768", "offset data range to that of unsigned short"
	case uint32:
		off, text, comment = 2147483648, "2147483648", "offset data range to that of unsigned long"
	case uint64:
		off, text, comment = 0x1p63, "9223372036854775808", "offset data range to that of unsigned long long"
	default:
		return newScaling(1, 0), nil, nil
	}
	return newScaling(1, off), []header.Card{
		rawNumberCard(header.KeyBzero, text, off, comment),
		rawNumberCard(header.KeyBscale, "1", 1, "default scaling factor"),
	}, nil
}

// checkWriteScaling refuses scaling cards that the rule accepts but that
// cannot store the data as the caller meant (see writeScaling).
func checkWriteScaling(k imageCards) error {
	for i, name := range cardKeys {
		if k.count[i] > 1 {
			return fmt.Errorf("%w: %d %s cards (the readers apply the last one)", ErrInvalidScaling, k.count[i], name)
		}
	}
	if k.count[keyBscale] > 0 && k.bscale == 0 {
		return fmt.Errorf("%w: BSCALE = 0 would store every value as BZERO (cfitsio refuses it: ZERO_SCALE)", ErrInvalidScaling)
	}
	return nil
}

// rawNumberCard builds a numeric card whose on-disk text is exactly text
// (right-justified in columns 11-30, as cfitsio writes it). The integer
// 9223372036854775808 does not fit an int64 card value, so the card keeps
// its text in Raw and its value as the float64 the reader uses.
func rawNumberCard(key, text string, v float64, comment string) header.Card {
	line := fmt.Sprintf("%-8s= %20s / %s", key, text, comment)
	var raw [header.CardWidth]byte
	copy(raw[:], fmt.Sprintf("%-80s", line))
	return header.Card{Key: key, Value: v, Type: header.TypeFloat, Comment: comment, Raw: raw}
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
// pipeline. We carry the type tag, extra header cards (the BZERO/BSCALE of
// an offset convention) and a bytes-producing method.
type anyData struct {
	bp    bitpix.BITPIX
	n     int64
	cards []header.Card
	emit  func(w io.Writer) error
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
	// In ModeEdit the append runs under the undo journal (journal.go).
	rec, err := f.journaledAppend(f.bw.Pos(), func() (*hduRecord, error) {
		return emitImageHDU(f.bw, bp, hdr, shape, data, isPrimary, len(f.hdus))
	})
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
	hb := newHeaderBuilder()
	if isPrimary {
		hb.set(header.KeySimple, true, "conforms to FITS standard")
	} else {
		hb.set(header.KeyXtension, "IMAGE   ", "image extension")
	}
	hb.set(header.KeyBitpix, int64(bp), "bits per pixel")
	hb.set(header.KeyNaxis, int64(len(shape)), "number of axes")
	for i, n := range shape {
		hb.set("NAXIS"+strconv.Itoa(i+1), n, "")
	}
	if !isPrimary {
		hb.set(header.KeyPcount, int64(0), "no parameter bytes")
		hb.set(header.KeyGcount, int64(1), "one group")
	}
	// The offset-convention BZERO/BSCALE cards (raw bytes, kept verbatim).
	for _, c := range data.cards {
		hb.add(c)
	}
	// Append user-supplied keywords, filtering out any that would duplicate
	// the mandatory ones.
	if userHdr != nil {
		for _, c := range userHdr.Cards() {
			if isMandatoryImageKey(c.Key) {
				continue
			}
			hb.add(c)
		}
	}

	if hb.err != nil {
		return nil, fmt.Errorf("fits: image header: %w", hb.err)
	}
	hdrBytes, err := header.Encode(hb.h)
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
		parsed:      hb.h,
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
