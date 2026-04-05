// Package compress implements decoders for the FITS tile-compression
// algorithms defined by the "Compressed Images" FITS convention
// (Pence et al. 2010). It is a pure bytes-in / bytes-out library that
// knows nothing about FITS headers or binary tables — the FITS-side
// wrapping lives in the root fits package as the CompressedImageHDU
// type.
//
// # Supported algorithms
//
//	RICE_1      — Rice coding with per-block k parameter (Pence 1998,
//	              White & Greenfield 1998). The workhorse for JWST and
//	              most HST CCD data.
//	GZIP_1      — Raw pixel bytes through zlib DEFLATE.
//	GZIP_2      — Byte-shuffled GZIP: bytes are grouped by significance
//	              before gzip. Improves ratio on smooth imagery.
//	NOCOMPRESS  — Identity pass-through (used by the tile-level
//	              fallback chain when compression would inflate a tile).
//	PLIO_1      — IRAF Planio run-length encoding for integer masks.
//	HCOMPRESS_1 — H-transform with quadtree coding (White 1991).
//
// # Coverage status
//
// All six algorithms listed above are recognized by Select. RICE_1,
// GZIP_1, GZIP_2, NOCOMPRESS, and PLIO_1 have full decoders. HCOMPRESS_1
// currently returns ErrUnsupportedCompression and will ship in a
// follow-up — it is less common in modern data products than the
// others.
package compress

import (
	"errors"
	"fmt"
	"strings"
)

// Algorithm identifies a tile-compression algorithm.
type Algorithm int

const (
	// Unknown is the zero value; returned by ParseAlgorithm for unrecognized
	// strings.
	Unknown Algorithm = iota
	// RICE1 is the RICE_1 algorithm (Rice coding with per-block k).
	RICE1
	// GZIP1 is plain zlib DEFLATE over the raw pixel bytes.
	GZIP1
	// GZIP2 is zlib DEFLATE over byte-shuffled pixel data.
	GZIP2
	// HCOMPRESS1 is H-transform + quadtree coding (White 1991).
	HCOMPRESS1
	// PLIO1 is IRAF Planio run-length encoding for integer masks.
	PLIO1
	// NoCompress is an identity pass-through. Used by cfitsio when a tile
	// would inflate under the primary compressor.
	NoCompress
)

// String returns the canonical FITS ZCMPTYPE string for this algorithm.
func (a Algorithm) String() string {
	switch a {
	case RICE1:
		return "RICE_1"
	case GZIP1:
		return "GZIP_1"
	case GZIP2:
		return "GZIP_2"
	case HCOMPRESS1:
		return "HCOMPRESS_1"
	case PLIO1:
		return "PLIO_1"
	case NoCompress:
		return "NOCOMPRESS"
	}
	return "UNKNOWN"
}

// ParseAlgorithm maps a ZCMPTYPE value to an Algorithm. Leading and
// trailing whitespace is trimmed; the comparison is case-insensitive.
// Returns Unknown for values not in the standard set.
func ParseAlgorithm(s string) Algorithm {
	switch strings.ToUpper(strings.TrimSpace(s)) {
	case "RICE_1", "RICE_ONE":
		return RICE1
	case "GZIP_1":
		return GZIP1
	case "GZIP_2":
		return GZIP2
	case "HCOMPRESS_1":
		return HCOMPRESS1
	case "PLIO_1":
		return PLIO1
	case "NOCOMPRESS":
		return NoCompress
	}
	return Unknown
}

// Params carries the per-algorithm compression parameters (ZNAMEi=ZVALi
// pairs from the header). Each algorithm extracts what it needs:
//
//	RICE_1:      BLOCKSIZE (default 32), BYTEPIX (1/2/4/8)
//	GZIP_1:      none (byte stream is self-describing)
//	GZIP_2:      BYTEPIX (derived from ZBITPIX if absent)
//	HCOMPRESS_1: SCALE, SMOOTH
//	PLIO_1:      none
//
// Values are int64 to match the FITS integer-card representation.
type Params map[string]int64

// Get returns the parameter value or def if absent.
func (p Params) Get(name string, def int64) int64 {
	if p == nil {
		return def
	}
	if v, ok := p[name]; ok {
		return v
	}
	return def
}

// Decoder decompresses the bytes of a single tile. The same Decoder
// instance may be reused across multiple tiles in a FITS image — it
// is stateless with respect to the FITS layer.
//
// src is the compressed bytes for one tile (the contents of one row
// of the COMPRESSED_DATA variable-length column). dst is a pre-allocated
// buffer that must be large enough for the uncompressed output: for an
// integer tile, len(dst) == nelem * elemSize bytes.
//
// On success, dst is filled with big-endian uncompressed pixel bytes
// in the same byte order FITS uses on disk. nelem is the number of
// pixels in the tile; elemSize is the on-disk byte size of each pixel
// (1, 2, 4, or 8) as given by ZBITPIX / BYTEPIX.
type Decoder interface {
	Decode(src, dst []byte, nelem, elemSize int) error
}

// Sentinel errors.
var (
	// ErrUnsupportedCompression is returned when the algorithm is known
	// but not yet implemented by this package. The v1.1 build implements
	// all algorithms except HCOMPRESS_1.
	ErrUnsupportedCompression = errors.New("fits/compress: unsupported compression algorithm")

	// ErrUnknownCompression is returned when the algorithm string does
	// not match any standard FITS compression type.
	ErrUnknownCompression = errors.New("fits/compress: unknown compression algorithm")

	// ErrCorrupt is returned when a compressed stream is malformed or
	// truncated.
	ErrCorrupt = errors.New("fits/compress: corrupt compressed stream")
)

// Select returns the Decoder for the given algorithm, configured from
// params. Returns ErrUnknownCompression for Unknown, or
// ErrUnsupportedCompression for algorithms this package recognizes but
// does not yet implement.
func Select(algo Algorithm, params Params) (Decoder, error) {
	switch algo {
	case RICE1:
		return newRICE1(params), nil
	case GZIP1:
		return newGZIP1(), nil
	case GZIP2:
		return newGZIP2(params), nil
	case NoCompress:
		return newNoCompress(), nil
	case PLIO1:
		return newPLIO1(), nil
	case HCOMPRESS1:
		return newHCOMPRESS1(params), nil
	case Unknown:
		return nil, ErrUnknownCompression
	}
	return nil, fmt.Errorf("%w: algorithm %d", ErrUnknownCompression, int(algo))
}
