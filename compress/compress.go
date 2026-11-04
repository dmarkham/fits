// Package compress implements the FITS tile-compression codecs of the
// "Compressed Images" convention (Pence et al. 2010, FITS standard 4.0
// section 10). It is a pure bytes-in / bytes-out library that knows
// nothing about FITS headers or binary tables — the FITS-side wrapping
// lives in the root fits package (CompressedImageHDU for reading,
// AppendCompressedImage / AppendCompressedFloat32Image for writing).
//
// # Algorithms
//
//	RICE_1      — Rice coding with per-block k parameter (Pence 1998,
//	              White & Greenfield 1998).
//	GZIP_1      — Raw pixel bytes through DEFLATE.
//	GZIP_2      — Byte-shuffled GZIP: bytes are grouped by significance
//	              before DEFLATE.
//	PLIO_1      — IRAF PLIO run-length encoding for integer masks
//	              (pixel values 0 .. 2^24-1).
//	HCOMPRESS_1 — H-transform with quadtree coding (White 1991).
//	NOCOMPRESS  — Identity tiles. Decoded for completeness (astropy can
//	              write it), but cfitsio 4.7.0 and funpack cannot read
//	              ZCMPTYPE = 'NOCOMPRESS' ("unknown compression
//	              algorithm"), so the fits writers refuse to produce it.
//	              cfitsio itself never falls back to NOCOMPRESS; its
//	              per-tile fallback is the GZIP_COMPRESSED_DATA column
//	              for float tiles that cannot be quantized.
//
// Select returns a decoder for each of the six. SelectEncoder returns
// an encoder for each; the RICE_1, HCOMPRESS_1 and PLIO_1 encoders are
// ports of cfitsio's fits_rcomp/_short/_byte, fits_hcompress/
// fits_hcompress64 and pl_p2li whose output is byte-identical to
// libcfitsio 4.7.0 on adversarial inputs (encode_cref_test.go, goldens
// from testdata/cref). The RICE_1 and HCOMPRESS_1 encoders accept
// BYTEPIX 1, 2 and 4 only, as cfitsio does.
//
// The float quantizer (QuantizeFloat32, QuantizeFloat64, ImgStatsFloat32)
// is a port of cfitsio's fits_quantize_float / fits_quantize_double and
// fits_img_stats_float, bit-exact against libcfitsio (cref_test.go; the
// float64 path through whole BITPIX -64 images, golden img64), and FitsRandoms is cfitsio's
// fits_rand_value dither table (randoms_cref_test.go). Decoder accuracy
// against real compressed files is tested separately from the encoders.
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
	// NoCompress is the identity "algorithm" (ZCMPTYPE = 'NOCOMPRESS').
	// Readable here, but not by cfitsio/funpack, so never written by the
	// fits writers.
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

// ParseAlgorithm maps a ZCMPTYPE value to an Algorithm the way cfitsio
// recognises it (imcompress.c imcomp_get_compressed_image_par): trailing
// blanks are insignificant (FITS strings), the comparison is
// case-sensitive, and 'RICE_ONE' is an alias for RICE_1. So 'rice_1' is
// Unknown, as cfitsio refuses it (TestZCMPTYPEComparedLikeCfitsio in the
// root package). Returns Unknown for values not in the standard set.
func ParseAlgorithm(s string) Algorithm {
	switch strings.TrimRight(s, " ") {
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

// Params carries the per-algorithm compression parameters. On read, the
// fits package fills them the way cfitsio does: RICE_1 BLOCKSIZE and
// BYTEPIX and HCOMPRESS_1 SMOOTH are taken from ZVAL1/ZVAL2 by POSITION,
// not from the ZNAMEi names (compressed.go parseCompressedMetadata,
// proven by TestCompressedCorpusRead). The HCOMPRESS_1 SCALE is not
// passed on read (the stream carries it). Each algorithm extracts what
// it needs:
//
//	RICE_1:      BLOCKSIZE (default 32), BYTEPIX 1/2/4 (8 is refused by
//	             the decoder and the encoder: cfitsio has no 64-bit Rice
//	             coder)
//	GZIP_1:      none (byte stream is self-describing)
//	GZIP_2:      BYTEPIX (the shuffle width; must equal elemSize)
//	HCOMPRESS_1: SCALE (write only), SMOOTH
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
	// ErrUnsupportedCompression is returned when an algorithm is
	// recognized but the requested use is not supported: an encoder
	// asked for a pixel width cfitsio does not support (RICE_1 and
	// HCOMPRESS_1 BYTEPIX 8), the RICE_1 decoder asked for BYTEPIX 8, or
	// PLIO_1 pixels outside 0..2^24-1
	// (ErrPLIORange wraps it).
	ErrUnsupportedCompression = errors.New("fits/compress: unsupported compression algorithm")

	// ErrUnknownCompression is returned when the algorithm string does
	// not match any standard FITS compression type.
	ErrUnknownCompression = errors.New("fits/compress: unknown compression algorithm")

	// ErrCorrupt is returned when a compressed stream is malformed or
	// truncated.
	ErrCorrupt = errors.New("fits/compress: corrupt compressed stream")
)

// Select returns the Decoder for the given algorithm, configured from
// params. Returns ErrUnknownCompression for unknown algorithm strings,
// or ErrUnsupportedCompression if a future FITS standard revision adds
// an algorithm this package does not support.
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
