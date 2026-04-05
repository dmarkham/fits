package fits

import (
	"fmt"
	"strconv"
	"sync"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
)

// CompressedImageHDU represents a tile-compressed image stored as a
// binary table HDU with ZIMAGE=T per the "Compressed Images" FITS
// convention (Pence et al. 2010).
//
// From the user's point of view, a CompressedImageHDU behaves exactly
// like an ImageHDU: the same ReadPixels[T] / ReadPixelsContext[T] /
// ReadPixelsMasked[T] / BITPIX / NAXIS / Shape / Header methods are
// available, and the decompression pipeline runs transparently under
// the covers. Callers who care can ask Compressed() to tell them it's
// backed by tiles rather than raw bytes.
//
// The physical backing is a binary table with:
//   - COMPRESSED_DATA VLA column — one row per tile, carrying the
//     compressed byte stream for that tile.
//   - UNCOMPRESSED_DATA VLA column (optional) — fallback raw bytes for
//     tiles that could not be compressed better than raw.
//   - GZIP_COMPRESSED_DATA VLA column (optional) — secondary fallback
//     when the primary algorithm fails and gzip is used instead.
//   - ZSCALE, ZZERO columns (for float data) — per-tile quantization
//     parameters.
//   - ZBLANK column (optional) — per-tile null sentinel.
type CompressedImageHDU struct {
	rec *hduRecord
	tbl *BinaryTableHDU // the underlying physical binary table

	metaOnce sync.Once
	meta     *compressedMetadata
	metaErr  error
}

// compressedMetadata holds the parsed Z* keywords describing the
// uncompressed image's shape, type, and tile/algorithm parameters.
type compressedMetadata struct {
	Bitpix   bitpix.BITPIX
	Naxis    int
	Shape    []int64 // ZNAXIS1..ZNAXISn
	Tile     []int64 // ZTILE1..ZTILEn
	CmpType  compress.Algorithm
	CmpParams compress.Params
	Quantize string // ZQUANTIZ
	Dither0  int64  // ZDITHER0
	Blank    int64  // ZBLANK keyword
	BlankSet bool

	// BSCALE/BZERO for the uncompressed image live in ZBSCALE/ZBZERO (or
	// BSCALE/BZERO if not prefixed). They apply after decompression.
	BScale float64
	BZero  float64

	// Whole-image ZSCALE/ZZERO fallback when per-tile columns are absent.
	ZScaleFixed float64
	ZZeroFixed  float64
	HasZScale   bool
}

// Type returns TypeImage — the compressed HDU presents an image-shaped
// view to callers regardless of its physical binary-table backing.
func (h *CompressedImageHDU) Type() HDUType { return TypeImage }

// Index returns the 0-based HDU index within the parent file.
func (h *CompressedImageHDU) Index() int { return h.rec.index }

// Header returns the parsed header of the underlying binary table.
// Callers see every keyword on the HDU, including the Z* compression
// keywords. To access the "logical" image header (what the uncompressed
// image would have looked like), reconstruct it from the Z-prefixed
// keywords via the compressedMetadata fields below.
func (h *CompressedImageHDU) Header() *header.Header {
	hdr, err := h.rec.loadHeader()
	if err != nil {
		panic(err)
	}
	return hdr
}

// BITPIX returns the logical BITPIX of the uncompressed image (from
// ZBITPIX), NOT the BITPIX=8 of the wrapping binary table.
func (h *CompressedImageHDU) BITPIX() int {
	m, _ := h.metadata()
	if m == nil {
		return 0
	}
	return int(m.Bitpix)
}

// NAXIS returns the logical NAXIS of the uncompressed image (ZNAXIS).
func (h *CompressedImageHDU) NAXIS() int {
	m, _ := h.metadata()
	if m == nil {
		return 0
	}
	return m.Naxis
}

// Shape returns the logical image dimensions (ZNAXIS1..ZNAXISn).
func (h *CompressedImageHDU) Shape() []int64 {
	m, _ := h.metadata()
	if m == nil {
		return nil
	}
	out := make([]int64, len(m.Shape))
	copy(out, m.Shape)
	return out
}

// BSCALE returns the BSCALE from the compressed image's Z-prefixed
// keyword (ZBSCALE) or falls back to BSCALE. Defaults to 1.0.
func (h *CompressedImageHDU) BSCALE() float64 {
	m, _ := h.metadata()
	if m == nil || m.BScale == 0 {
		return 1.0
	}
	return m.BScale
}

// BZERO returns the BZERO for the uncompressed image. Defaults to 0.
func (h *CompressedImageHDU) BZERO() float64 {
	m, _ := h.metadata()
	if m == nil {
		return 0
	}
	return m.BZero
}

// Compressed always returns true for a CompressedImageHDU.
func (h *CompressedImageHDU) Compressed() bool { return true }

// CompressionType returns the ZCMPTYPE value (e.g. "RICE_1").
func (h *CompressedImageHDU) CompressionType() string {
	m, _ := h.metadata()
	if m == nil {
		return ""
	}
	return m.CmpType.String()
}

// metadata returns the parsed Z* keyword set, loading it lazily on first
// access.
func (h *CompressedImageHDU) metadata() (*compressedMetadata, error) {
	h.metaOnce.Do(func() {
		h.meta, h.metaErr = parseCompressedMetadata(h.Header())
	})
	return h.meta, h.metaErr
}

// parseCompressedMetadata reads the Z* keywords from the given header
// into a compressedMetadata struct.
func parseCompressedMetadata(hdr *header.Header) (*compressedMetadata, error) {
	if v, _ := hdr.Bool("ZIMAGE"); !v {
		return nil, fmt.Errorf("fits: CompressedImageHDU: ZIMAGE not true")
	}
	m := &compressedMetadata{
		BScale:    1.0,
		CmpParams: make(compress.Params),
	}
	if v, err := hdr.Int("ZBITPIX"); err == nil {
		m.Bitpix = bitpix.BITPIX(v)
	} else {
		return nil, fmt.Errorf("fits: CompressedImageHDU: missing ZBITPIX: %w", err)
	}
	if v, err := hdr.Int("ZNAXIS"); err == nil {
		m.Naxis = int(v)
	} else {
		return nil, fmt.Errorf("fits: CompressedImageHDU: missing ZNAXIS: %w", err)
	}
	m.Shape = make([]int64, m.Naxis)
	m.Tile = make([]int64, m.Naxis)
	for i := 1; i <= m.Naxis; i++ {
		if v, err := hdr.Int("ZNAXIS" + strconv.Itoa(i)); err == nil {
			m.Shape[i-1] = v
		} else {
			return nil, fmt.Errorf("fits: CompressedImageHDU: missing ZNAXIS%d", i)
		}
		if v, err := hdr.Int("ZTILE" + strconv.Itoa(i)); err == nil {
			m.Tile[i-1] = v
		} else {
			// ZTILE defaults: the first axis defaults to ZNAXIS1, the
			// rest default to 1 (whole-row tiles are the cfitsio default).
			if i == 1 {
				m.Tile[0] = m.Shape[0]
			} else {
				m.Tile[i-1] = 1
			}
		}
	}
	if s, err := hdr.String("ZCMPTYPE"); err == nil {
		m.CmpType = compress.ParseAlgorithm(s)
	} else {
		return nil, fmt.Errorf("fits: CompressedImageHDU: missing ZCMPTYPE")
	}
	// Parse ZNAMEi / ZVALi parameter pairs (up to i=999 in theory, but in
	// practice up to ~4).
	for i := 1; i < 100; i++ {
		name, err := hdr.String("ZNAME" + strconv.Itoa(i))
		if err != nil {
			break
		}
		val, err := hdr.Int("ZVAL" + strconv.Itoa(i))
		if err != nil {
			continue
		}
		m.CmpParams[name] = val
	}
	// BYTEPIX is sometimes absent for GZIP_2; derive from ZBITPIX if so.
	if _, ok := m.CmpParams["BYTEPIX"]; !ok {
		m.CmpParams["BYTEPIX"] = int64(m.Bitpix.Size())
	}
	if s, err := hdr.String("ZQUANTIZ"); err == nil {
		m.Quantize = s
	}
	if v, err := hdr.Int("ZDITHER0"); err == nil {
		m.Dither0 = v
	}
	if v, err := hdr.Int("ZBLANK"); err == nil {
		m.Blank = v
		m.BlankSet = true
	}
	if v, err := hdr.Float("ZBSCALE"); err == nil {
		m.BScale = v
	}
	if v, err := hdr.Float("ZBZERO"); err == nil {
		m.BZero = v
	}
	if v, err := hdr.Float("ZSCALE"); err == nil {
		m.ZScaleFixed = v
		m.HasZScale = true
	}
	if v, err := hdr.Float("ZZERO"); err == nil {
		m.ZZeroFixed = v
	}
	return m, nil
}
