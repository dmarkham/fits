package compress

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
)

// gzipFallback decompresses a real RFC 1952 gzip stream (as opposed to
// zlib-wrapped DEFLATE). cfitsio-written compressed FITS files use gzip
// format; astropy uses zlib. We accept either by trying zlib first and
// falling back here on zlib header mismatch.
func gzipFallback(src, dst []byte, want int) error {
	r, err := gzip.NewReader(bytes.NewReader(src))
	if err != nil {
		return fmt.Errorf("%w: gzip fallback: %v", ErrCorrupt, err)
	}
	defer r.Close()
	n, err := io.ReadFull(r, dst[:want])
	if err != nil && err != io.EOF && err != io.ErrUnexpectedEOF {
		return fmt.Errorf("%w: gzip fallback read: %v", ErrCorrupt, err)
	}
	if n < want {
		// Fewer bytes than expected — still succeed if EOF; caller will
		// zero-fill any trailing pixels.
		return nil
	}
	return nil
}
