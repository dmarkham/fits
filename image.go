package fits

import (
	"fmt"

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

// Compressed returns false — plain ImageHDUs are raw, not tile-compressed.
func (*ImageHDU) Compressed() bool { return false }

// CompressionType returns the empty string for uncompressed images.
func (*ImageHDU) CompressionType() string { return "" }

// Header returns the parsed header, loading it lazily on first call.
//
// It panics if package header cannot parse a card. Open decodes only the
// structural cards (SIMPLE/XTENSION, BITPIX, NAXISn, PCOUNT, GCOUNT, ...),
// so a file can open, and its pixels can be read (BSCALE/BZERO/BLANK are
// read from the stored cards directly), even though some other card is
// malformed: use HeaderOf to get the parse error instead of a panic.
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

// BSCALE returns the scaling factor every read and write of this package
// applies to the image, by the scaling-card rule (scalingCards in
// scaling.go): 1 when there is no BSCALE card; the LAST of duplicate
// BSCALE cards, as cfitsio and astropy apply it. It returns a wrapped
// ErrInvalidScaling, as every pixel read of the HDU does, when that card
// is not a valid value card with a FITS number that cfitsio converts to
// the same value (a quoted string, logical, complex, Inf/NaN, hex,
// lowercase d exponent, out-of-range or junk value, "BSCALE  =2.0", a
// lowercase "bscale" card), or when a HIERARCH BSCALE / BZERO / BLANK
// card is present (cfitsio applies it, astropy ignores it). It is read
// from the header as parsed (including unflushed in-memory edits, as they
// would be written) or, before the full header is parsed, straight from
// the stored cards, with the same answer either way. Proven against
// cfitsio and astropy by TestScalingKeywordsRule.
func (h *ImageHDU) BSCALE() (float64, error) {
	s, _, _, err := h.rec.imageKeywords()
	if err != nil {
		return 0, fmt.Errorf("fits: HDU %d: %w", h.rec.index, err)
	}
	return s.bscale, nil
}

// BZERO returns the offset applied to the image, read like BSCALE (0 when
// there is no BZERO card). The value is a float64, as in cfitsio; the
// uint64 convention's 9223372036854775808 (2^63) is exact.
func (h *ImageHDU) BZERO() (float64, error) {
	s, _, _, err := h.rec.imageKeywords()
	if err != nil {
		return 0, fmt.Errorf("fits: HDU %d: %w", h.rec.index, err)
	}
	return s.bzero, nil
}

// BLANK returns the null value of an integer image (the stored value that
// marks an undefined pixel) by the same rule, and whether there is one. A
// float image has none (FITS 4.0 §4.4.2.4: NaN marks nulls).
func (h *ImageHDU) BLANK() (int64, bool, error) {
	_, b, ok, err := h.rec.imageKeywords()
	if err != nil {
		return 0, false, fmt.Errorf("fits: HDU %d: %w", h.rec.index, err)
	}
	return b, ok, nil
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
