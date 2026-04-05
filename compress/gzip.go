package compress

import (
	"bytes"
	"compress/zlib"
	"fmt"
	"io"
)

// ------------------ GZIP_1 ------------------

// gzip1Decoder decompresses a zlib DEFLATE stream back to raw pixel bytes.
// GZIP_1 is the simplest "real" tile-compression algorithm: the raw pixel
// bytes are fed through zlib and the result is stored in the tile's
// COMPRESSED_DATA row. On decode, we just reverse that — zlib gives us
// back the pre-compression byte sequence verbatim.
//
// Note: despite the "GZIP" name in the FITS spec, the actual bytes on
// disk are zlib-wrapped DEFLATE (RFC 1950), NOT the gzip wrapper format
// (RFC 1952). Go's compress/zlib handles this directly.
type gzip1Decoder struct{}

func newGZIP1() Decoder { return gzip1Decoder{} }

func (gzip1Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	want := nelem * elemSize
	if len(dst) < want {
		return fmt.Errorf("%w: GZIP_1 dst length %d < needed %d", ErrCorrupt, len(dst), want)
	}
	r, err := zlib.NewReader(bytes.NewReader(src))
	if err != nil {
		// cfitsio actually writes bare DEFLATE + gzip wrapper, which
		// zlib.NewReader rejects. Fall through to a raw DEFLATE attempt
		// via flate if the zlib header check fails.
		return decodeRawGzip(src, dst, want)
	}
	defer r.Close()
	if _, err := io.ReadFull(r, dst[:want]); err != nil {
		// Try the raw-gzip path as a fallback before giving up.
		if errRaw := decodeRawGzip(src, dst, want); errRaw == nil {
			return nil
		}
		return fmt.Errorf("%w: GZIP_1 read: %v", ErrCorrupt, err)
	}
	return nil
}

// decodeRawGzip handles the case where the FITS writer emitted a real
// gzip-format stream (RFC 1952) instead of zlib-wrapped DEFLATE. cfitsio
// uses gzip format; astropy uses zlib. We accept either.
func decodeRawGzip(src, dst []byte, want int) error {
	// Try compress/gzip — imports are deferred to keep the hot path lean.
	return gzipFallback(src, dst, want)
}

// ------------------ GZIP_2 ------------------

// gzip2Decoder decompresses a GZIP_2 tile: zlib DEFLATE wrapped around a
// byte-shuffled version of the pixels.
//
// The shuffle groups bytes by their significance position within each
// pixel. For a BYTEPIX=2 image, the raw pixel byte stream is
//
//	b0_0 b0_1 b1_0 b1_1 b2_0 b2_1 ...   (pixel 0 byte 0, pixel 0 byte 1, ...)
//
// After shuffle it becomes
//
//	b0_0 b1_0 b2_0 ... b0_1 b1_1 b2_1 ...  (all first bytes, then all second bytes)
//
// Smooth imagery has correlated high-order bytes and nearly-random
// low-order bytes, so grouping them improves gzip's ratio noticeably.
//
// On decode we zlib-decompress then invert the shuffle.
type gzip2Decoder struct {
	bytepix int
}

func newGZIP2(params Params) Decoder {
	// BYTEPIX may be supplied explicitly; if not, the FITS-side layer
	// derives it from ZBITPIX and passes it in via params["BYTEPIX"].
	return gzip2Decoder{bytepix: int(params.Get("BYTEPIX", 0))}
}

func (g gzip2Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	bp := g.bytepix
	if bp == 0 {
		bp = elemSize
	}
	want := nelem * elemSize
	if len(dst) < want {
		return fmt.Errorf("%w: GZIP_2 dst length %d < needed %d", ErrCorrupt, len(dst), want)
	}
	// Decompress into a scratch buffer.
	scratch := make([]byte, want)
	r, err := zlib.NewReader(bytes.NewReader(src))
	if err != nil {
		if errRaw := gzipFallback(src, scratch, want); errRaw != nil {
			return fmt.Errorf("%w: GZIP_2 zlib init: %v", ErrCorrupt, err)
		}
	} else {
		if _, err := io.ReadFull(r, scratch); err != nil {
			r.Close()
			if errRaw := gzipFallback(src, scratch, want); errRaw != nil {
				return fmt.Errorf("%w: GZIP_2 read: %v", ErrCorrupt, err)
			}
		} else {
			r.Close()
		}
	}
	// Unshuffle: scratch contains nelem*bp bytes arranged as bp groups of
	// nelem bytes each. We write out interleaved pixel order.
	if bp == 1 {
		copy(dst, scratch)
		return nil
	}
	for i := 0; i < nelem; i++ {
		for b := 0; b < bp; b++ {
			dst[i*bp+b] = scratch[b*nelem+i]
		}
	}
	return nil
}
