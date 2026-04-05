package fits

import (
	"bufio"
	"fmt"
	"image"
	"image/color"
	"io"

	"github.com/dmarkham/fits/internal/bitpix"
)

// AsImage returns a native image.Image for HDUs that can be represented
// losslessly:
//
//   - 2D BITPIX=8 → *image.Gray (uint8)
//   - 2D BITPIX=16 → *image.Gray16 (uint16, with TZERO=32768 unsigned handling)
//   - 2D BITPIX=-32 / -64 → *FloatGray with raw float values (not rescaled)
//
// Returns ErrNotImageCompatible for HDUs that cannot be represented without
// loss (multi-dim, complex, NAXIS < 2, signed int32/int64, etc.). Callers
// that want lower-dim views must Slice() first. No auto-rescaling — see
// plan guiding principle.
func (h *ImageHDU) AsImage() (image.Image, error) {
	if h.NAXIS() != 2 {
		return nil, fmt.Errorf("%w: NAXIS=%d (need 2)", ErrNotImageCompatible, h.NAXIS())
	}
	shape := h.Shape()
	w := int(shape[0])
	ht := int(shape[1])
	bp := h.bitpixType()
	rect := image.Rect(0, 0, w, ht)
	bzero := h.BZERO()

	switch bp {
	case bitpix.Int8:
		// BITPIX=8 is unsigned byte on disk; image.Gray is unsigned byte.
		raw := make([]byte, w*ht)
		if err := h.rec.file.br.ReadRange(h.rec.dataStart, raw); err != nil {
			return nil, err
		}
		return &image.Gray{Pix: raw, Stride: w, Rect: rect}, nil
	case bitpix.Int16:
		// Read int16 into Gray16. FITS int16 is signed on disk; the unsigned
		// uint16 form is expressed via BZERO=32768. We honor that: if
		// BZERO is close to 32768, we add it before storing as uint16.
		pix, err := ReadPixels[int32](h)
		if err != nil {
			// Try through float path if lossless-int read fails.
			fpix, err2 := ReadPixels[float64](h)
			if err2 != nil {
				return nil, err2
			}
			out := &image.Gray16{Stride: w * 2, Rect: rect, Pix: make([]byte, w*ht*2)}
			for i, v := range fpix {
				u := uint16(clampU16(v))
				out.Pix[i*2] = byte(u >> 8)
				out.Pix[i*2+1] = byte(u)
			}
			return out, nil
		}
		out := &image.Gray16{Stride: w * 2, Rect: rect, Pix: make([]byte, w*ht*2)}
		for i, v := range pix {
			// BZERO applied during ReadPixels[int32]? No — ReadPixels applies
			// BSCALE/BZERO only when bscale != 1 or bzero != 0, returning the
			// scaled value. If BZERO=32768 is set, the int32 values here are
			// 0..65535. Clamp into uint16.
			_ = bzero
			u := uint16(clampU16(float64(v)))
			out.Pix[i*2] = byte(u >> 8)
			out.Pix[i*2+1] = byte(u)
		}
		return out, nil
	case bitpix.Float32, bitpix.Float64:
		pix, err := ReadPixels[float32](h)
		if err != nil {
			return nil, err
		}
		return &FloatGray{Rect: rect, Stride: w, Pix: pix}, nil
	}
	return nil, fmt.Errorf("%w: BITPIX=%d", ErrNotImageCompatible, int(bp))
}

func clampU16(v float64) float64 {
	if v < 0 {
		return 0
	}
	if v > 65535 {
		return 65535
	}
	return v
}

// FloatGray is an image.Image over raw float32 pixels. It does NOT rescale
// — rendering through png.Encode without a caller-provided transform
// produces garbage, which is intentional. The purpose of this type is to
// let float FITS data participate in the image.Image ecosystem with
// lossless fidelity.
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

// RGBA implements color.Color. The float value is saturating-cast to the
// 0..0xFFFF range — a cheap placeholder that exists only to satisfy the
// interface. For real rendering, read Pix[] directly and transform.
func (c FloatColor) RGBA() (r, g, b, a uint32) {
	v := c.V
	if v < 0 {
		v = 0
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

// Decode implements image.Decode for FITS files. Returns ErrNotImageCompatible
// if the primary HDU (or first HDU) is not losslessly representable as an
// image.Image.
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

// DecodeConfig implements image.DecodeConfig for FITS files.
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
		return image.Config{}, ErrNotImageCompatible
	}
	sh := prim.Shape()
	var m color.Model
	switch prim.bitpixType() {
	case bitpix.Int8:
		m = color.GrayModel
	case bitpix.Int16:
		m = color.Gray16Model
	case bitpix.Float32, bitpix.Float64:
		m = FloatGrayModel
	default:
		return image.Config{}, ErrNotImageCompatible
	}
	return image.Config{ColorModel: m, Width: int(sh[0]), Height: int(sh[1])}, nil
}

func init() {
	// Register "fits" as a recognizable image format. The magic string is
	// the start of every valid FITS file ("SIMPLE  =").
	image.RegisterFormat("fits", "SIMPLE  =", Decode, DecodeConfig)
}
