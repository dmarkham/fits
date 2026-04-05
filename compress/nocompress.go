package compress

import "fmt"

// nocompressDecoder is the identity decoder: src and dst are the same
// bytes. Used by cfitsio's per-tile fallback chain when the primary
// compressor would inflate a tile relative to its raw bytes.
type nocompressDecoder struct{}

func newNoCompress() Decoder { return nocompressDecoder{} }

func (nocompressDecoder) Decode(src, dst []byte, nelem, elemSize int) error {
	want := nelem * elemSize
	if len(src) != want {
		return fmt.Errorf("%w: NOCOMPRESS src length %d != expected %d", ErrCorrupt, len(src), want)
	}
	if len(dst) < want {
		return fmt.Errorf("%w: NOCOMPRESS dst length %d < needed %d", ErrCorrupt, len(dst), want)
	}
	copy(dst, src)
	return nil
}
