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
// small, it returns ErrBufferTooSmall so the caller can retry with
// a larger buffer or fall back to NOCOMPRESS for that tile.
type Encoder interface {
	Encode(src, dst []byte, nelem, elemSize int) (int, error)
}

// ErrBufferTooSmall is returned by Encoder.Encode when the caller's
// destination buffer cannot hold the compressed output. Callers that
// hit this should retry with a larger buffer or use the NOCOMPRESS
// fallback for this specific tile.
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
