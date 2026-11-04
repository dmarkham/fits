package fits

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	"io"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// AsImage returns the 2D image as a native image.Image:
//
//   - BITPIX=8 with BSCALE=1, BZERO=0 → *image.Gray holding the stored
//     bytes (lossless)
//   - BITPIX=16 with BSCALE=1, BZERO=32768 (unsigned 16-bit, the camera
//     format) → *image.Gray16 holding the physical values 0..65535
//     (lossless)
//   - BITPIX=-32 with BSCALE=1, BZERO=0 → *FloatGray holding the stored
//     floats bit for bit (NaN included)
//   - BITPIX=-64, or -32 with BSCALE/BZERO → *FloatGray holding the
//     physical values raw·BSCALE + BZERO (computed in float64) only if every
//     one of them is exactly representable in float32 (NaN allowed);
//     otherwise ErrNotImageCompatible, since FloatGray is float32
//
// Everything else returns ErrNotImageCompatible instead of a lossy
// conversion: a BSCALE or BZERO card whose value is not a finite FITS
// number, quoted numbers such as '10.0' included (rather than guessing
// the default, as cfitsio does), a last BSCALE / BZERO / BLANK card that
// is not a valid value card (the card rule of ReadFloat32), NAXIS != 2 (Slice
// first), signed 16-bit (BZERO=0), signed
// bytes (BZERO=-128) or any other BSCALE/BZERO on integer data, BITPIX=32
// and 64, an integer image with a pixel equal to BLANK (image.Gray has no
// null), and a float image with a value float32 cannot hold exactly. No
// auto-rescaling — see the guiding principle in the README.
//
// Orientation: image pixel (x, y) is FITS pixel (x+1, y+1); the first
// stored row is y=0, the top row of an image.Image. FITS conventionally
// draws y upward, so an image encoder (png.Encode, ...) shows the frame
// flipped vertically relative to DS9/Siril display.
//
// Proven against astropy for every sample in testdata/imaging
// (TestAsImage_Matrix, TestAsImage_RowOrientation).
func (h *ImageHDU) AsImage() (image.Image, error) {
	if h.NAXIS() != 2 {
		return nil, fmt.Errorf("%w: NAXIS=%d (need 2)", ErrNotImageCompatible, h.NAXIS())
	}
	cards, err := h.rec.imageCards()
	if err != nil {
		return nil, fmt.Errorf("%w: %v", ErrNotImageCompatible, err)
	}
	bscale, bzero := cards.bscale, cards.bzero
	kind, err := imageKind(h.bitpixType(), bscale, bzero)
	if err != nil {
		return nil, err
	}
	shape := h.Shape()
	w := int(shape[0])
	ht := int(shape[1])
	rect := image.Rect(0, 0, w, ht)
	n := h.pixelCount()

	switch kind {
	case kindGray, kindGray16:
		raw, err := readRawBytes(h, n)
		if err != nil {
			return nil, err
		}
		if blank := cards.blank; cards.haveBlank {
			for i := int64(0); i < n; i++ {
				var v int64
				if kind == kindGray {
					v = int64(raw[i])
				} else {
					v = int64(int16(uint16(raw[2*i])<<8 | uint16(raw[2*i+1])))
				}
				if v == blank {
					return nil, fmt.Errorf("%w: pixel %d equals BLANK=%d (image.Image has no null value; use ReadPixelsMasked)", ErrNotImageCompatible, i, blank)
				}
			}
		}
		if kind == kindGray {
			return &image.Gray{Pix: raw, Stride: w, Rect: rect}, nil
		}
		// Stored int16 + 32768 is the physical uint16; flipping the sign
		// bit of the big-endian high byte gives its big-endian bytes,
		// which is exactly Gray16's layout.
		for i := int64(0); i < n; i++ {
			raw[2*i] ^= 0x80
		}
		return &image.Gray16{Pix: raw, Stride: 2 * w, Rect: rect}, nil
	default: // float
		// Decoded here from the stored bytes with the BSCALE / BZERO read
		// above (not through ReadPixels, whose header reading differs).
		raw, err := readRawBytes(h, n)
		if err != nil {
			return nil, err
		}
		bp := h.bitpixType()
		pix := make([]float32, n)
		for i := range pix {
			if bp == bitpix.Float32 && bscale == 1 && bzero == 0 {
				pix[i] = bigendian.Float32(raw[4*i:]) // stored bits, NaN payloads included
				continue
			}
			var v float64
			if bp == bitpix.Float32 {
				v = float64(bigendian.Float32(raw[4*i:]))
			} else {
				v = bigendian.Float64(raw[8*i:])
			}
			if bscale != 1 || bzero != 0 {
				v = float64(v*bscale) + bzero
			}
			pix[i] = float32(v)
			if float64(pix[i]) != v && v == v {
				return nil, fmt.Errorf("%w: pixel %d value %v is not exactly representable as float32", ErrNotImageCompatible, i, v)
			}
		}
		return &FloatGray{Rect: rect, Stride: w, Pix: pix}, nil
	}
}

type stdImageKind int

const (
	kindGray stdImageKind = iota
	kindGray16
	kindFloat
)

// imageKind decides from the header alone which image.Image type a 2D
// image maps to (see AsImage).
func imageKind(bp bitpix.BITPIX, bscale, bzero float64) (stdImageKind, error) {
	switch {
	case bp == bitpix.Int8 && bscale == 1 && bzero == 0:
		return kindGray, nil
	case bp == bitpix.Int16 && bscale == 1 && bzero == 32768:
		return kindGray16, nil
	case bp == bitpix.Float32 || bp == bitpix.Float64:
		if bscale == 0 {
			return 0, fmt.Errorf("%w: BSCALE=0", ErrNotImageCompatible)
		}
		return kindFloat, nil
	}
	return 0, fmt.Errorf("%w: BITPIX=%d BSCALE=%g BZERO=%g has no lossless image.Image form", ErrNotImageCompatible, int(bp), bscale, bzero)
}

// FloatGray is an image.Image over float32 pixels. It does NOT rescale
// — rendering through png.Encode without a caller-provided transform
// produces garbage, which is intentional. The purpose of this type is to
// let float FITS data participate in the image.Image ecosystem without
// changing the values (bit for bit for unscaled BITPIX=-32; see AsImage).
// Pixel (x, y) is Pix[(y-Rect.Min.Y)*Stride + (x-Rect.Min.X)]; At outside
// Rect returns FloatColor{0}.
type FloatGray struct {
	Rect   image.Rectangle
	Stride int
	Pix    []float32
}

// ColorModel returns a pass-through float color model.
func (f *FloatGray) ColorModel() color.Model { return FloatGrayModel }

// Bounds returns the rectangle.
func (f *FloatGray) Bounds() image.Rectangle { return f.Rect }

// At returns the raw float32 at (x, y), boxed in a FloatColor.
func (f *FloatGray) At(x, y int) color.Color {
	if !(image.Point{x, y}.In(f.Rect)) {
		return FloatColor{V: 0}
	}
	return FloatColor{V: f.Pix[(y-f.Rect.Min.Y)*f.Stride+(x-f.Rect.Min.X)]}
}

// FloatColor wraps a float32 as a color.Color. RGBA returns the float value
// scaled to the 16-bit RGBA space via a lossy saturating cast. This is the
// only place the library does anything that could be called "scaling", and
// it exists solely to satisfy the color.Color interface contract. Callers
// that need faithful display MUST NOT rely on this mapping — they must
// build their own transform from the raw Pix[] slab.
type FloatColor struct {
	V float32
}

// RGBA implements color.Color. The float value is clamped to [0, 1] (NaN
// and negatives give 0, values above 1 and +Inf give 0xFFFF) and
// truncated to 0..0xFFFF — a cheap placeholder that exists only to satisfy
// the interface. For real rendering, read Pix[] directly and transform.
func (c FloatColor) RGBA() (r, g, b, a uint32) {
	v := c.V
	if !(v > 0) { // zero, negative or NaN: no separate NaN branch, so the
		// result never depends on Go's implementation-defined uint32(NaN)
		return 0, 0, 0, 0xFFFF
	}
	if v > 1 {
		v = 1
	}
	u := uint32(v * 0xFFFF)
	return u, u, u, 0xFFFF
}

// floatGrayModel implements color.Model for FloatGray.
type floatGrayModel struct{}

func (floatGrayModel) Convert(c color.Color) color.Color {
	if fc, ok := c.(FloatColor); ok {
		return fc
	}
	r, _, _, _ := c.RGBA()
	return FloatColor{V: float32(r) / 0xFFFF}
}

// FloatGrayModel is the color model for FloatGray images.
var FloatGrayModel color.Model = floatGrayModel{}

// Decode implements image.Decode for FITS files: it reads the whole
// stream and returns the primary HDU's AsImage (same rules and errors).
// Only the primary HDU is used, so a tile-compressed file (NAXIS=0
// primary) returns ErrNotImageCompatible.
func Decode(r io.Reader) (image.Image, error) {
	// image.Decode passes a bufio.Reader; we need a seeker. We slurp.
	f, err := OpenReadAll(bufio.NewReader(r))
	if err != nil {
		return nil, err
	}
	defer f.Close()
	prim, err := f.Primary()
	if err != nil {
		return nil, err
	}
	return prim.AsImage()
}

// DecodeConfig implements image.DecodeConfig for FITS files. It applies
// AsImage's rules from the header alone, so the cases it cannot see are
// data-dependent: an integer pixel equal to BLANK, or a -64/scaled float
// value float32 cannot hold exactly (Decode then fails).
func DecodeConfig(r io.Reader) (image.Config, error) {
	f, err := OpenReadAll(bufio.NewReader(r))
	if err != nil {
		return image.Config{}, err
	}
	defer f.Close()
	prim, err := f.Primary()
	if err != nil {
		return image.Config{}, err
	}
	if prim.NAXIS() != 2 {
		return image.Config{}, fmt.Errorf("%w: NAXIS=%d (need 2)", ErrNotImageCompatible, prim.NAXIS())
	}
	cards, err := prim.rec.imageCards()
	if err != nil {
		return image.Config{}, fmt.Errorf("%w: %v", ErrNotImageCompatible, err)
	}
	kind, err := imageKind(prim.bitpixType(), cards.bscale, cards.bzero)
	if err != nil {
		return image.Config{}, err
	}
	sh := prim.Shape()
	m := map[stdImageKind]color.Model{kindGray: color.GrayModel, kindGray16: color.Gray16Model, kindFloat: FloatGrayModel}[kind]
	return image.Config{ColorModel: m, Width: int(sh[0]), Height: int(sh[1])}, nil
}

func init() {
	// Register "fits" as a recognizable image format. The magic string is
	// the start of every valid FITS file ("SIMPLE  =").
	image.RegisterFormat("fits", "SIMPLE  =", Decode, DecodeConfig)
}
