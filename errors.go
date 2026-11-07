// Package fits implements a pure-Go reader/writer for the FITS
// (Flexible Image Transport System) file format (FITS Standard v4.0).
// Random-groups HDUs (§6) are located, so the HDUs after them are read,
// but not decoded: (*File).HDU returns ErrRandomGroups for them.
//
// This file defines the public error sentinels and typed errors.
package fits

import (
	"errors"
	"fmt"
)

// Sentinel errors — match via errors.Is.
var (
	// ErrEmptyFile is returned when an input stream is zero bytes.
	ErrEmptyFile = errors.New("fits: empty file")

	// ErrNotFITS is returned when the first block does not begin with
	// "SIMPLE  =".
	ErrNotFITS = errors.New("fits: not a FITS file")

	// ErrKeyNotFound is reserved for keyword lookups at the File level. No
	// function in this package returns it today: the typed getters of
	// package header return header.ErrKeyNotFound.
	ErrKeyNotFound = errors.New("fits: key not found")

	// ErrNotImageHDU is reserved for image operations on non-image HDUs.
	// No function returns it today (image functions take *ImageHDU, so the
	// mismatch is a compile-time error).
	ErrNotImageHDU = errors.New("fits: not an image HDU")

	// ErrNotImageCompatible is returned by AsImage, Decode and DecodeConfig
	// for HDUs that cannot be losslessly represented as an image.Image:
	// NAXIS != 2, signed or scaled integer data (other than BITPIX 16 with
	// BZERO 32768), BLANK pixels, float values float32 cannot hold exactly,
	// BSCALE = 0, and tile-compressed files (the NAXIS = 0 primary).
	ErrNotImageCompatible = errors.New("fits: not image.Image-compatible")

	// ErrTypeMismatch is returned when a generic read is asked for a target
	// type that cannot losslessly represent the on-disk type.
	ErrTypeMismatch = errors.New("fits: type mismatch")

	// ErrReadOnly is returned when a write or edit operation is attempted on
	// a read-only *File.
	ErrReadOnly = errors.New("fits: file is read-only")

	// ErrShapeMismatch is returned (wrapped) when OverwritePixels,
	// WriteImage or AppendImage receive data whose length is not the
	// product of the image shape.
	ErrShapeMismatch = errors.New("fits: shape mismatch")

	// ErrRandomGroups is returned for random-groups (§6) HDUs, which
	// this library does not support by design — random groups are a
	// legacy radio-interferometry format superseded by binary tables.
	ErrRandomGroups = errors.New("fits: random-groups HDUs not supported")

	// ErrTiledTableCompression is returned by (*File).HDU for a binary
	// table with ZTABLE = T (tiled table compression, FITS 4.0 §10.3, as
	// written by fpack -table), which this library does not decode: read
	// as a plain table its rows and columns would be the compressed
	// bytes, not the data.
	ErrTiledTableCompression = errors.New("fits: tiled table compression (ZTABLE = T) not supported")

	// ErrUnsupportedExtension is returned (as an *UnsupportedExtensionError)
	// by (*File).HDU for an extension whose XTENSION value is not one this
	// library decodes (IMAGE, TABLE, BINTABLE and the legacy aliases
	// IUEIMAGE, A3DTABLE, 3DTABLE), e.g. 'FOREIGN' or 'DUMP'. The HDU is
	// still framed (its size follows FITS 4.0 eq. 2, as for any
	// extension), so the HDUs after it are read normally. cfitsio reads
	// such an HDU as an IMAGE_HDU; this library does not treat its bytes
	// as pixels.
	ErrUnsupportedExtension = errors.New("fits: unsupported extension type")
)

// UnsupportedExtensionError reports an extension of a type this library
// does not decode. errors.Is(err, ErrUnsupportedExtension) is true for it.
type UnsupportedExtensionError struct {
	Index    int    // HDU index (0 = primary)
	XTension string // the XTENSION value, trailing blanks removed
}

func (e *UnsupportedExtensionError) Error() string {
	return fmt.Sprintf("fits: HDU %d: unsupported extension type XTENSION = '%s'", e.Index, e.XTension)
}

// Is reports whether target is ErrUnsupportedExtension.
func (e *UnsupportedExtensionError) Is(target error) bool { return target == ErrUnsupportedExtension }

// ErrUnterminatedHeader describes a header that does not contain an END card
// within a generous block limit.
type ErrUnterminatedHeader struct {
	HDU          int
	BytesScanned int64
}

func (e *ErrUnterminatedHeader) Error() string {
	return fmt.Sprintf("fits: HDU %d: unterminated header (scanned %d bytes)", e.HDU, e.BytesScanned)
}

// ErrTruncatedData describes a data section that ends before the declared
// data size.
type ErrTruncatedData struct {
	HDU      int
	Expected int64
	Got      int64
}

func (e *ErrTruncatedData) Error() string {
	return fmt.Sprintf("fits: HDU %d: data truncated: expected %d bytes, got %d", e.HDU, e.Expected, e.Got)
}

// ErrDuplicateKeyword is matched (errors.Is) by *DuplicateKeywordError.
var ErrDuplicateKeyword = errors.New("fits: duplicate structural keyword")

// DuplicateKeywordError reports a header that carries a structural keyword
// (SIMPLE, XTENSION, BITPIX, NAXIS, NAXISn for n <= NAXIS, PCOUNT, GCOUNT,
// TFIELDS) more than once with different values. cfitsio reads BITPIX,
// NAXIS and NAXISn by card position and takes the last PCOUNT/GCOUNT of an
// image, so such a file has no single reading; Open refuses it instead of
// guessing the element size or the HDU length. A keyword repeated with the
// same value reads the same either way and is accepted.
type DuplicateKeywordError struct {
	HDU     int
	Keyword string
}

func (e *DuplicateKeywordError) Error() string {
	return fmt.Sprintf("fits: HDU %d: structural keyword %q appears more than once", e.HDU, e.Keyword)
}

func (e *DuplicateKeywordError) Is(target error) bool { return target == ErrDuplicateKeyword }

// ErrInvalidScaling is returned (wrapped, with the reason) when the
// BSCALE/BZERO/BLANK cards of an image are refused by the scaling-card
// rule (scalingCards in scaling.go: a card cfitsio and astropy read
// differently, e.g. a value that is not a FITS number, a HIERARCH BSCALE,
// a BLANK outside the BITPIX range or equal to 1234554321), by every read
// of the pixels (ReadPixels, ReadPixelsMasked, ReadFloat32, AsImage, the
// tile-compressed readers), by ImageHDU.BSCALE / BZERO / BLANK and by
// OverwritePixels.
//
// WriteImage also returns it for the caller's header when, beyond that
// rule, BSCALE = 0 (cfitsio's fits_set_bscale refuses it with
// ZERO_SCALE), a keyword is given twice, or a card cannot be encoded;
// OverwritePixels when the HDU's BSCALE is 0.
var ErrInvalidScaling = errors.New("fits: invalid image scaling keywords")

// ErrUndefinedChecksum is returned by ChecksumStatus and VerifyChecksum
// for an HDU whose DATASUM or CHECKSUM card has no value at all (cfitsio's
// fits_verify_chksum fails with VALUE_UNDEFINED there).
var ErrUndefinedChecksum = errors.New("fits: DATASUM or CHECKSUM card has no value")

// ErrMissingRequiredKeyword describes a missing mandatory keyword.
type ErrMissingRequiredKeyword struct {
	HDU     int
	Keyword string
}

func (e *ErrMissingRequiredKeyword) Error() string {
	return fmt.Sprintf("fits: HDU %d: missing required keyword %q", e.HDU, e.Keyword)
}

// TypeMismatchError carries detail about an ErrTypeMismatch failure.
// BSCALE and BZERO are the image's scaling (1 and 0 when absent or not
// applicable).
type TypeMismatchError struct {
	Requested string
	BITPIX    int
	Lossy     bool
	BSCALE    float64
	BZERO     float64
}

func (e *TypeMismatchError) Error() string {
	if e.BSCALE != 0 && (e.BSCALE != 1 || e.BZERO != 0) {
		return fmt.Sprintf("fits: cannot represent BITPIX=%d data with BSCALE=%v BZERO=%v as %s (lossy=%v)", e.BITPIX, e.BSCALE, e.BZERO, e.Requested, e.Lossy)
	}
	return fmt.Sprintf("fits: cannot represent BITPIX=%d data as %s (lossy=%v)", e.BITPIX, e.Requested, e.Lossy)
}

func (e *TypeMismatchError) Is(target error) bool {
	return target == ErrTypeMismatch
}

// KeyNotFoundError carries detail about an ErrKeyNotFound failure (not
// returned by any function today; see ErrKeyNotFound).
type KeyNotFoundError struct {
	Key string
	HDU int
}

func (e *KeyNotFoundError) Error() string {
	return fmt.Sprintf("fits: HDU %d: key %q not found", e.HDU, e.Key)
}

func (e *KeyNotFoundError) Is(target error) bool {
	return target == ErrKeyNotFound
}

// ErrNumOverflow is returned by WriteImage and OverwritePixels, which then
// write nothing, when a value cannot be stored under the image's BITPIX
// and BSCALE/BZERO: for integer BITPIX, (value - BZERO) / BSCALE lies
// outside the stored range (BITPIX 8: 0..255, 16/32/64: the signed type)
// by more than cfitsio's 0.49 rounding margin, which is where cfitsio's
// fits_write_img reports NUM_OVERFLOW (412) after storing the range limit;
// for float BITPIX, a finite value whose scaled form is not finite in the
// stored type (cfitsio stores the infinity without an error).
var ErrNumOverflow = errors.New("fits: value out of range for the stored type")
