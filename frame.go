package fits

import (
	"fmt"

	"github.com/dmarkham/fits/header"
)

// frameImage is the minimal HDU surface the one-shot frame helpers need.
// Satisfied by both *ImageHDU and *CompressedImageHDU.
type frameImage interface {
	HDU
	NAXIS() int
	Shape() []int64
}

// selectFrameImage picks the HDU the one-shot frame helpers read from:
//
//   - primary HDU if it has image data (NAXIS >= 2)
//   - otherwise HDU 1 if it is an image or tile-compressed image
//     (the standard astropy/fpack layout where the primary is a NAXIS=0
//     placeholder and the image lives in the first extension)
//   - otherwise an error
func selectFrameImage(f *File) (frameImage, error) {
	prim, err := f.Primary()
	if err != nil {
		return nil, err
	}
	if prim.NAXIS() >= 2 {
		return prim, nil
	}
	if f.NumHDU() < 2 {
		return nil, fmt.Errorf("fits: no image data: primary NAXIS=%d and file has no extensions", prim.NAXIS())
	}
	h, err := f.HDU(1)
	if err != nil {
		return nil, err
	}
	if img, ok := h.(frameImage); ok {
		return img, nil
	}
	return nil, fmt.Errorf("fits: no image data: primary is NAXIS=%d placeholder and HDU 1 is %s", prim.NAXIS(), h.Type())
}

// ReadFrame reads a 2D mono image from a FITS file path and returns the
// pixels normalized to [0, 1] via ReadFloat32.
//
// Returned dims follow the FITS convention: width == NAXIS1 (columns, the
// fastest-varying axis), height == NAXIS2 (rows). Pixels are row-major:
// index pixel (x, y) as pix[y*width + x].
//
// Accepts:
//   - NAXIS == 2
//   - NAXIS == 3 with NAXIS3 == 1 (degenerate single-plane cube)
//
// The image is taken from the primary HDU if it has image data, otherwise
// from HDU 1 (the standard tile-compressed layout). Tile-compressed images
// are supported transparently. Integer BITPIX is normalized to [0, 1];
// float BITPIX passes through unchanged. NaN and BLANK pixels are not
// masked — callers who need that should use ReadPixelsMasked directly.
//
// Returns an error for NAXIS != 2 (or NAXIS == 3 with NAXIS3 != 1),
// with a hint pointing at ReadFrameRGB for 3-plane color cubes.
func ReadFrame(path string) (pix []float32, width, height int, err error) {
	f, err := Open(path)
	if err != nil {
		return nil, 0, 0, err
	}
	// Surface Close errors only if the read itself succeeded.
	defer func() {
		if cerr := f.Close(); cerr != nil && err == nil {
			err = cerr
		}
	}()

	img, err := selectFrameImage(f)
	if err != nil {
		return nil, 0, 0, err
	}
	shape := img.Shape()
	switch {
	case img.NAXIS() == 2:
		// ok
	case img.NAXIS() == 3 && len(shape) >= 3 && shape[2] == 1:
		// ok — degenerate single-plane cube
	case img.NAXIS() == 3 && len(shape) >= 3 && shape[2] == 3:
		return nil, 0, 0, fmt.Errorf("fits.ReadFrame: image is NAXIS=3 with 3 planes; use ReadFrameRGB")
	case img.NAXIS() == 3 && len(shape) >= 3:
		return nil, 0, 0, fmt.Errorf("fits.ReadFrame: image is NAXIS=3 with %d planes, not a 2D or RGB image; use ReadPixels directly", shape[2])
	default:
		return nil, 0, 0, fmt.Errorf("fits.ReadFrame: expected NAXIS=2, got NAXIS=%d", img.NAXIS())
	}

	pix, err = ReadFloat32(img)
	if err != nil {
		return nil, 0, 0, err
	}
	return pix, int(shape[0]), int(shape[1]), nil
}

// ReadFrameRGB reads a 3D planar RGB image (NAXIS=3, NAXIS3=3) from a FITS
// file path and returns the three color planes normalized to [0, 1] via
// ReadFloat32.
//
// Returned dims follow the FITS convention: width == NAXIS1 (columns, the
// fastest-varying axis), height == NAXIS2 (rows). Each plane is row-major:
// index pixel (x, y) as r[y*width + x] (and likewise for g, b).
//
// Layout on disk is R plane, then G plane, then B plane (NAXIS1 fastest,
// NAXIS3=3 slowest) — the standard FITS/Siril/astropy convention.
//
// HDU selection, normalization, and masking semantics match ReadFrame.
// Returns an error for NAXIS != 3 or NAXIS3 != 3.
func ReadFrameRGB(path string) (r, g, b []float32, width, height int, err error) {
	f, err := Open(path)
	if err != nil {
		return nil, nil, nil, 0, 0, err
	}
	// Surface Close errors only if the read itself succeeded.
	defer func() {
		if cerr := f.Close(); cerr != nil && err == nil {
			err = cerr
		}
	}()

	img, err := selectFrameImage(f)
	if err != nil {
		return nil, nil, nil, 0, 0, err
	}
	shape := img.Shape()
	if img.NAXIS() != 3 || len(shape) < 3 || shape[2] != 3 {
		return nil, nil, nil, 0, 0, fmt.Errorf("fits.ReadFrameRGB: expected NAXIS=3 with 3 planes, got NAXIS=%d shape=%v", img.NAXIS(), shape)
	}

	pix, err := ReadFloat32(img)
	if err != nil {
		return nil, nil, nil, 0, 0, err
	}
	w, h := int(shape[0]), int(shape[1])
	npix := w * h
	if len(pix) != 3*npix {
		return nil, nil, nil, 0, 0, fmt.Errorf("fits.ReadFrameRGB: pixel count %d does not match 3 * %d * %d", len(pix), w, h)
	}
	r = make([]float32, npix)
	g = make([]float32, npix)
	b = make([]float32, npix)
	copy(r, pix[:npix])
	copy(g, pix[npix:2*npix])
	copy(b, pix[2*npix:3*npix])
	return r, g, b, w, h, nil
}

// WriteMono writes a 2D mono float32 image to a FITS file at path.
//
// The output is a BITPIX=-32 primary HDU with shape [width, height]
// (NAXIS1=width, NAXIS2=height). Input pixels are row-major — pixel (x, y)
// at pix[y*width + x] — matching the layout ReadFrame returns. Values are
// written verbatim (no scaling, no clamping). Any existing file at path
// is truncated.
func WriteMono(path string, pix []float32, width, height int) error {
	if width <= 0 || height <= 0 {
		return fmt.Errorf("fits.WriteMono: width=%d height=%d, both must be > 0", width, height)
	}
	if len(pix) != width*height {
		return fmt.Errorf("fits.WriteMono: len(pix)=%d, want width*height=%d", len(pix), width*height)
	}

	f, err := Create(path)
	if err != nil {
		return err
	}
	if _, err := WriteImage(f, header.New(), []int64{int64(width), int64(height)}, pix); err != nil {
		_ = f.Close()
		return err
	}
	return f.Close()
}

// WriteRGB writes a 3D planar RGB float32 image to a FITS file at path.
//
// The output is a BITPIX=-32 primary HDU with shape [width, height, 3]
// (NAXIS1=width, NAXIS2=height, NAXIS3=3). The three planes are
// concatenated in R, G, B order (NAXIS1 fastest, NAXIS3=3 slowest) — the
// standard FITS/Siril/astropy convention. Each input plane is row-major —
// pixel (x, y) at r[y*width + x] (and likewise for g, b) — matching the
// layout ReadFrameRGB returns. Values are written verbatim (no scaling,
// no clamping). Any existing file at path is truncated.
func WriteRGB(path string, r, g, b []float32, width, height int) error {
	if width <= 0 || height <= 0 {
		return fmt.Errorf("fits.WriteRGB: width=%d height=%d, both must be > 0", width, height)
	}
	npix := width * height
	if len(r) != npix || len(g) != npix || len(b) != npix {
		return fmt.Errorf("fits.WriteRGB: len(r)=%d len(g)=%d len(b)=%d, want all == width*height=%d",
			len(r), len(g), len(b), npix)
	}

	data := make([]float32, 3*npix)
	copy(data[:npix], r)
	copy(data[npix:2*npix], g)
	copy(data[2*npix:3*npix], b)

	f, err := Create(path)
	if err != nil {
		return err
	}
	if _, err := WriteImage(f, header.New(), []int64{int64(width), int64(height), 3}, data); err != nil {
		_ = f.Close()
		return err
	}
	return f.Close()
}
