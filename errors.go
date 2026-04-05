// Package fits implements a pure-Go reader/writer for the FITS
// (Flexible Image Transport System) file format as defined by the
// FITS Standard v4.0.
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

	// ErrKeyNotFound is returned by typed keyword getters when the key is
	// absent from a header.
	ErrKeyNotFound = errors.New("fits: key not found")

	// ErrNotImageHDU is returned when an image operation is attempted on a
	// non-image HDU.
	ErrNotImageHDU = errors.New("fits: not an image HDU")

	// ErrNotImageCompatible is returned by AsImage() for HDUs that cannot be
	// losslessly represented as an image.Image (multi-dim, complex BITPIX,
	// etc.).
	ErrNotImageCompatible = errors.New("fits: not image.Image-compatible")

	// ErrTypeMismatch is returned when a generic read is asked for a target
	// type that cannot losslessly represent the on-disk type.
	ErrTypeMismatch = errors.New("fits: type mismatch")

	// ErrReadOnly is returned when a write or edit operation is attempted on
	// a read-only *File.
	ErrReadOnly = errors.New("fits: file is read-only")

	// ErrShapeMismatch is returned when OverwritePixels or a similar
	// shape-preserving call receives data of the wrong length or dimension.
	ErrShapeMismatch = errors.New("fits: shape mismatch")

	// ErrRandomGroups is returned for random-groups (§6) HDUs, which
	// this library does not support by design — random groups are a
	// legacy radio-interferometry format superseded by binary tables.
	ErrRandomGroups = errors.New("fits: random-groups HDUs not supported")
)

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

// ErrMissingRequiredKeyword describes a missing mandatory keyword.
type ErrMissingRequiredKeyword struct {
	HDU     int
	Keyword string
}

func (e *ErrMissingRequiredKeyword) Error() string {
	return fmt.Sprintf("fits: HDU %d: missing required keyword %q", e.HDU, e.Keyword)
}

// TypeMismatchError carries detail about an ErrTypeMismatch failure.
type TypeMismatchError struct {
	Requested string
	BITPIX    int
	Lossy     bool
}

func (e *TypeMismatchError) Error() string {
	return fmt.Sprintf("fits: cannot represent BITPIX=%d data as %s (lossy=%v)", e.BITPIX, e.Requested, e.Lossy)
}

func (e *TypeMismatchError) Is(target error) bool {
	return target == ErrTypeMismatch
}

// KeyNotFoundError carries detail about an ErrKeyNotFound failure.
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
