// Package bitpix defines the BITPIX enumeration and helpers.
//
// The FITS BITPIX keyword identifies the physical on-disk type of image and
// random-group data and the implicit cell type of binary-table image cells
// (FITS Standard v4.0 §4.4.1.1, Table 8).
package bitpix

import "fmt"

// BITPIX is the FITS data-type encoding stored in the BITPIX header keyword.
// Positive values are integer types (8 is an unsigned byte, 16/32/64 are
// two's-complement signed), negative values are IEEE-754 floating point.
// The constant names give the width only: Int8 is BITPIX 8, which stores
// unsigned bytes 0..255.
type BITPIX int

// Mandatory BITPIX values defined by the FITS standard.
const (
	Int8    BITPIX = 8
	Int16   BITPIX = 16
	Int32   BITPIX = 32
	Int64   BITPIX = 64
	Float32 BITPIX = -32
	Float64 BITPIX = -64
)

// Size returns the on-disk byte size of one element of type b.
// It panics for values that are not a valid BITPIX.
func (b BITPIX) Size() int {
	switch b {
	case Int8:
		return 1
	case Int16:
		return 2
	case Int32, Float32:
		return 4
	case Int64, Float64:
		return 8
	}
	panic(fmt.Sprintf("bitpix: invalid BITPIX value %d", int(b)))
}

// Valid reports whether b is one of the BITPIX values defined by the FITS
// standard.
func (b BITPIX) Valid() bool {
	switch b {
	case Int8, Int16, Int32, Int64, Float32, Float64:
		return true
	}
	return false
}

// IsFloat reports whether b encodes an IEEE-754 floating-point type.
func (b BITPIX) IsFloat() bool { return b < 0 }

// IsInt reports whether b encodes an integer type.
func (b BITPIX) IsInt() bool { return b > 0 }

// Signed reports whether the raw on-disk representation of b is signed.
// BITPIX 16, 32 and 64 are signed; BITPIX 8 is an unsigned byte (FITS 4.0
// Table 8). Float types are always signed.
//
// Note: signed bytes and unsigned 16/32/64-bit integers are expressed with
// the BZERO (TZERO for table columns) offset convention, not a separate
// BITPIX: BITPIX 8 with BZERO = -128 holds int8, BITPIX 16/32/64 with
// BZERO = 2^15/2^31/2^63 hold uint16/uint32/uint64 (FITS 4.0 Table 11).
func (b BITPIX) Signed() bool {
	switch b {
	case Int16, Int32, Int64, Float32, Float64:
		return true
	case Int8:
		// BITPIX=8 on disk is unsigned 8-bit. We return false here so callers
		// that want "signedness of the raw on-disk type" get the right answer.
		return false
	}
	panic(fmt.Sprintf("bitpix: invalid BITPIX value %d", int(b)))
}

// String returns a short human-readable name for b.
func (b BITPIX) String() string {
	switch b {
	case Int8:
		return "uint8"
	case Int16:
		return "int16"
	case Int32:
		return "int32"
	case Int64:
		return "int64"
	case Float32:
		return "float32"
	case Float64:
		return "float64"
	}
	return fmt.Sprintf("BITPIX(%d)", int(b))
}
