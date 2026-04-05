package fits

import (
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
)

// ImageHDU represents a FITS primary HDU or an IMAGE extension.
type ImageHDU struct {
	rec *hduRecord
}

// Type returns TypeImage.
func (h *ImageHDU) Type() HDUType { return TypeImage }

// Index returns the 0-based HDU index.
func (h *ImageHDU) Index() int { return h.rec.index }

// Header returns the parsed header, loading it lazily on first call.
// Errors during parse panic here; callers that need to handle parse errors
// explicitly must call (*File).HDU followed by direct access via the header
// package — in practice the scan pass already validated structure, so a
// lazy parse failure indicates programmer error or disk corruption between
// open and access.
func (h *ImageHDU) Header() *header.Header {
	hdr, err := h.rec.loadHeader()
	if err != nil {
		// Structural parse errors at this stage indicate the raw header
		// bytes changed underfoot. Surface via panic — callers should not
		// be modifying the underlying stream.
		panic(err)
	}
	return hdr
}

// BITPIX returns the integer BITPIX value from the HDU's structural metadata.
func (h *ImageHDU) BITPIX() int { return int(h.rec.bitpix) }

// bitpixType returns the typed enum form used internally.
func (h *ImageHDU) bitpixType() bitpix.BITPIX { return h.rec.bitpix }

// NAXIS returns the number of axes.
func (h *ImageHDU) NAXIS() int { return h.rec.naxis }

// Shape returns the axis lengths: [NAXIS1, NAXIS2, ...]. The returned slice
// is safe to hold — a copy is made.
func (h *ImageHDU) Shape() []int64 {
	out := make([]int64, len(h.rec.shape))
	copy(out, h.rec.shape)
	return out
}

// BSCALE returns the BSCALE keyword value, or 1.0 if absent (§4.4.2.5).
func (h *ImageHDU) BSCALE() float64 {
	hdr := h.Header()
	if v, err := hdr.Float(header.KeyBscale); err == nil {
		return v
	}
	return 1.0
}

// BZERO returns the BZERO keyword value, or 0.0 if absent (§4.4.2.5).
func (h *ImageHDU) BZERO() float64 {
	hdr := h.Header()
	if v, err := hdr.Float(header.KeyBzero); err == nil {
		return v
	}
	return 0.0
}

// pixelCount returns the total number of pixels = NAXIS1*...*NAXISm, or 0
// for NAXIS=0 HDUs.
func (h *ImageHDU) pixelCount() int64 {
	if h.rec.naxis == 0 {
		return 0
	}
	var n int64 = 1
	for _, s := range h.rec.shape {
		n *= s
	}
	return n
}
