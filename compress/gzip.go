package compress

import (
	"fmt"
)

// ------------------ GZIP_1 ------------------

// gzip1Decoder inflates a GZIP_1 tile back to raw big-endian pixel bytes.
//
// cfitsio writes the tile as a gzip (RFC 1952) stream; astropy writes
// zlib-wrapped DEFLATE (RFC 1950). Both are accepted (inflate).
//
// The tile must inflate to exactly nelem*elemSize bytes. A tile that
// inflates to any other length, longer or shorter, is refused with
// ErrCorrupt, never truncated or zero-filled: cfitsio picks the tile's
// integer type from the inflated length (1, 2, 4 or 8 bytes per pixel,
// imcompress.c:6372-6417), so a length that is not nelem*elemSize means
// the stored values are not elemSize wide. Use InflateTile to read such a
// tile the way cfitsio does.
type gzip1Decoder struct{}

func newGZIP1() Decoder { return gzip1Decoder{} }

func (gzip1Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	want := nelem * elemSize
	if len(dst) < want {
		return fmt.Errorf("%w: GZIP_1 dst length %d < needed %d", ErrCorrupt, len(dst), want)
	}
	raw, err := inflate(src, want+1)
	if err != nil {
		return fmt.Errorf("%w: GZIP_1: %v", ErrCorrupt, err)
	}
	if len(raw) != want {
		return fmt.Errorf("%w: GZIP_1 tile inflates to %s bytes, want %d (%d pixels x %d bytes)", ErrCorrupt, lenString(raw, want), want, nelem, elemSize)
	}
	copy(dst, raw)
	return nil
}

// lenString formats an inflated length that was read with a limit of
// want+1 bytes.
func lenString(raw []byte, want int) string {
	if len(raw) > want {
		return fmt.Sprintf("more than %d", want)
	}
	return fmt.Sprint(len(raw))
}

// ------------------ GZIP_2 ------------------

// gzip2Decoder decompresses a GZIP_2 tile: DEFLATE wrapped around a
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
// On decode we inflate, then invert the shuffle. The inflated length must
// be exactly nelem*elemSize, as for GZIP_1.
type gzip2Decoder struct {
	bytepix int
}

func newGZIP2(params Params) Decoder {
	// BYTEPIX may be supplied explicitly; if not, the shuffle width is
	// elemSize.
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
	if bp != elemSize {
		return fmt.Errorf("%w: GZIP_2 BYTEPIX %d differs from the %d-byte pixels", ErrCorrupt, bp, elemSize)
	}
	raw, err := inflate(src, want+1)
	if err != nil {
		return fmt.Errorf("%w: GZIP_2: %v", ErrCorrupt, err)
	}
	if len(raw) != want {
		return fmt.Errorf("%w: GZIP_2 tile inflates to %s bytes, want %d (%d pixels x %d bytes)", ErrCorrupt, lenString(raw, want), want, nelem, elemSize)
	}
	unshuffle(raw, dst, nelem, bp)
	return nil
}

// unshuffle inverts the GZIP_2 byte shuffle of nelem values of width bp
// (cfitsio fits_unshuffle_2bytes/_4bytes/_8bytes): src holds bp groups of
// nelem bytes, dst gets the interleaved values.
func unshuffle(src, dst []byte, nelem, bp int) {
	if bp == 1 {
		copy(dst, src[:nelem])
		return
	}
	for i := 0; i < nelem; i++ {
		for b := 0; b < bp; b++ {
			dst[i*bp+b] = src[b*nelem+i]
		}
	}
}

// InflateTile inflates one GZIP_1 or GZIP_2 tile of nelem pixels the way
// cfitsio's imcomp_decompress_tile does (imcompress.c:6372-6417): the
// width of the stored values is taken from the inflated length, which
// must be exactly 1, 2, 4 or 8 times nelem; GZIP_2 is then un-shuffled by
// that width. Any other length (including a tile with leftover bytes) is
// refused with ErrCorrupt, as cfitsio refuses it ("uncompressed tile has
// wrong size").
//
// The result is width-byte big-endian values: unsigned for width 1,
// two's-complement signed for 2, 4 and 8, as cfitsio reads them (TBYTE,
// TSHORT, TINT, TLONGLONG). Old cfitsio versions wrote 16-bit images as
// 4-byte GZIP tiles (astropy's m13_gzip.fits is one); the caller converts
// the values to the image's ZBITPIX type.
func InflateTile(algo Algorithm, src []byte, nelem int) (data []byte, width int, err error) {
	if algo != GZIP1 && algo != GZIP2 {
		return nil, 0, fmt.Errorf("%w: InflateTile: %v is not GZIP_1 or GZIP_2", ErrUnsupportedCompression, algo)
	}
	if nelem <= 0 {
		return nil, 0, fmt.Errorf("%w: InflateTile: %d pixels", ErrCorrupt, nelem)
	}
	raw, err := inflate(src, 8*nelem+1)
	if err != nil {
		return nil, 0, fmt.Errorf("%w: %v: %v", ErrCorrupt, algo, err)
	}
	switch len(raw) {
	case nelem:
		width = 1
	case 2 * nelem:
		width = 2
	case 4 * nelem:
		width = 4
	case 8 * nelem:
		width = 8
	default:
		return nil, 0, fmt.Errorf("%w: %v tile inflates to %s bytes, not 1, 2, 4 or 8 x %d pixels", ErrCorrupt, algo, lenString(raw, 8*nelem), nelem)
	}
	if algo == GZIP1 || width == 1 {
		return raw, width, nil
	}
	out := make([]byte, len(raw))
	unshuffle(raw, out, nelem, width)
	return out, width, nil
}
