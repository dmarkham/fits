package compress

import (
	"bytes"
	"compress/gzip"
	"compress/zlib"
	"errors"
	"io"
)

// inflate decompresses a tile written either as a gzip (RFC 1952) stream,
// as cfitsio writes it, or as zlib-wrapped DEFLATE (RFC 1950), as astropy
// writes it. At most limit bytes are returned: a caller that needs an
// exact length n passes n+1 and so sees a too-long stream as n+1 bytes.
// A truncated or corrupt stream is an error, never a short result.
func inflate(src []byte, limit int) ([]byte, error) {
	var r io.ReadCloser
	if gz, err := gzip.NewReader(bytes.NewReader(src)); err == nil {
		gz.Multistream(false)
		r = gz
	} else if zr, zerr := zlib.NewReader(bytes.NewReader(src)); zerr == nil {
		r = zr
	} else {
		return nil, errors.New("tile is neither a gzip nor a zlib stream")
	}
	defer r.Close()
	out, err := io.ReadAll(io.LimitReader(r, int64(limit)))
	if err != nil {
		return nil, err
	}
	return out, nil
}
