package compress

import (
	"bytes"
	"compress/gzip"
	"compress/zlib"
	"fmt"
)

// GZIP_1 encoder — gzip-format (RFC 1952) DEFLATE over the raw pixel
// bytes. astropy.io.fits expects the gzip wrapper with the 0x1f 0x8b
// magic, not the bare zlib wrapper. Our decoder accepts both formats
// via a fallback path, but the encoder emits canonical gzip for maximum
// interop.
type gzip1Encoder struct{}

func newGZIP1Encoder() Encoder { return gzip1Encoder{} }

func (gzip1Encoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	want := nelem * elemSize
	if len(src) < want {
		return 0, fmt.Errorf("%w: GZIP_1 src length %d < expected %d", ErrCorrupt, len(src), want)
	}
	var buf bytes.Buffer
	w, err := gzip.NewWriterLevel(&buf, gzip.DefaultCompression)
	if err != nil {
		return 0, err
	}
	// Fixed metadata for reproducible output (no filename, zero mtime).
	w.Name = ""
	w.ModTime = zeroTime
	if _, err := w.Write(src[:want]); err != nil {
		return 0, err
	}
	if err := w.Close(); err != nil {
		return 0, err
	}
	if buf.Len() > len(dst) {
		return 0, ErrBufferTooSmall
	}
	return copy(dst, buf.Bytes()), nil
}

// Silence unused-import warning in case the zlib import is ever pruned.
var _ = zlib.DefaultCompression

// GZIP_2 encoder — byte-shuffle the pixel bytes by significance
// position (all 1st bytes, then all 2nd bytes, etc.) before zlib
// compression. Symmetric counterpart of the GZIP_2 decoder; preserves
// the exact on-disk bytes written by cfitsio.
type gzip2Encoder struct {
	bytepix int
}

func newGZIP2Encoder(params Params) Encoder {
	return gzip2Encoder{bytepix: int(params.Get("BYTEPIX", 0))}
}

func (g gzip2Encoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	bp := g.bytepix
	if bp == 0 {
		bp = elemSize
	}
	want := nelem * elemSize
	if len(src) < want {
		return 0, fmt.Errorf("%w: GZIP_2 src length %d < expected %d", ErrCorrupt, len(src), want)
	}
	// Shuffle into bp groups of nelem bytes each.
	var shuffled []byte
	if bp == 1 {
		shuffled = src[:want]
	} else {
		shuffled = make([]byte, want)
		for i := 0; i < nelem; i++ {
			for b := 0; b < bp; b++ {
				shuffled[b*nelem+i] = src[i*bp+b]
			}
		}
	}
	var buf bytes.Buffer
	w, err := gzip.NewWriterLevel(&buf, gzip.DefaultCompression)
	if err != nil {
		return 0, err
	}
	w.Name = ""
	w.ModTime = zeroTime
	if _, err := w.Write(shuffled); err != nil {
		return 0, err
	}
	if err := w.Close(); err != nil {
		return 0, err
	}
	if buf.Len() > len(dst) {
		return 0, ErrBufferTooSmall
	}
	return copy(dst, buf.Bytes()), nil
}
