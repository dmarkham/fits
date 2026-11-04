package compress

import (
	"fmt"
)

// Encoder is the compression counterpart of Decoder. Given a tile's
// uncompressed bytes in FITS big-endian order, it produces the
// compressed bytes that should be stored in the COMPRESSED_DATA VLA
// column for that row.
//
// The caller allocates dst large enough to hold the worst case; the
// encoder returns the actual number of bytes written. If dst is too
// small (including one byte short of the exact stream length), it
// returns ErrBufferTooSmall and the caller retries with a larger buffer.
// There is no NOCOMPRESS fallback: cfitsio cannot read NOCOMPRESS tiles,
// so the fits writers never produce them (see the package doc).
type Encoder interface {
	Encode(src, dst []byte, nelem, elemSize int) (int, error)
}

// ErrBufferTooSmall is returned by Encoder.Encode when the caller's
// destination buffer cannot hold the compressed output. Every encoder
// returns it rather than a truncated stream, including RICE_1 when only
// the final partial byte does not fit (checkExactBuffer in
// encode_cref_test.go). Callers retry with a larger buffer.
var ErrBufferTooSmall = fmt.Errorf("fits/compress: output buffer too small")

// SelectEncoder returns the Encoder for the given algorithm, configured
// from params. Mirrors Select() for decoders.
func SelectEncoder(algo Algorithm, params Params) (Encoder, error) {
	switch algo {
	case RICE1:
		return newRICE1Encoder(params), nil
	case GZIP1:
		return newGZIP1Encoder(), nil
	case GZIP2:
		return newGZIP2Encoder(params), nil
	case NoCompress:
		return nocompressEncoder{}, nil
	case PLIO1:
		return plio1Encoder{}, nil
	case HCOMPRESS1:
		return newHCOMPRESS1Encoder(params), nil
	case Unknown:
		return nil, ErrUnknownCompression
	}
	return nil, fmt.Errorf("%w: algorithm %d", ErrUnknownCompression, int(algo))
}

// ---------- NOCOMPRESS ----------

type nocompressEncoder struct{}

func (nocompressEncoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	if len(src) != nelem*elemSize {
		return 0, fmt.Errorf("%w: NOCOMPRESS src length %d != expected %d", ErrCorrupt, len(src), nelem*elemSize)
	}
	if len(dst) < len(src) {
		return 0, ErrBufferTooSmall
	}
	return copy(dst, src), nil
}
