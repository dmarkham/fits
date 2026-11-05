package fits

import (
	"fmt"
	"strconv"
	"strings"
	"sync"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
)

// CompressedImageHDU represents a tile-compressed image stored as a
// binary table HDU with ZIMAGE=T per the "Compressed Images" FITS
// convention (Pence et al. 2010), as written by fpack/cfitsio.
//
// It presents the logical image: BITPIX, NAXIS and Shape come from ZBITPIX
// and ZNAXISn, and BSCALE/BZERO from the BSCALE/BZERO keywords of the
// compressed header (cfitsio keeps them under their plain names). Pixels
// are read with ReadPixelsCompressed[T] / ReadPixelsCompressedContext[T] /
// ReadPixelsCompressedMasked[T]; the ImageHDU functions ReadPixels[T] and
// ReadPixelsMasked[T] take a *ImageHDU and do not accept this type.
// ReadFloat32 accepts both. Compressed() reports true.
//
// The physical backing is a binary table with:
//   - COMPRESSED_DATA VLA column — one row per tile, carrying the
//     compressed byte stream for that tile.
//   - GZIP_COMPRESSED_DATA VLA column (optional) — gzipped raw floats for
//     float tiles cfitsio could not quantize (used when COMPRESSED_DATA is
//     empty for that row).
//   - UNCOMPRESSED_DATA VLA column (optional, pre-2011 cfitsio) — the raw
//     tile values, used when COMPRESSED_DATA is empty for that row.
//   - ZSCALE, ZZERO columns (quantized float data) — per-tile quantization
//     parameters.
//   - ZBLANK column (optional) — per-tile null value.
//
// Read support and its proof are described on ReadPixelsCompressed.
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
	Bitpix    bitpix.BITPIX
	Naxis     int
	Shape     []int64 // ZNAXIS1..ZNAXISn
	Tile      []int64 // ZTILE1..ZTILEn
	CmpType   compress.Algorithm
	CmpParams compress.Params
	Quantize  string                // ZQUANTIZ as written ("" if absent)
	Dither    compress.DitherMethod // from ZQUANTIZ
	Dither0   int64                 // ZDITHER0; 1 when absent (cfitsio default)
	// Quantized reports whether float tiles hold quantized integers: true
	// unless the table has neither a ZSCALE nor a ZZERO column, or
	// ZQUANTIZ = 'NONE' (cfitsio imcomp_get_compressed_image_par). Without
	// those columns ZQUANTIZ is not checked at all, as in cfitsio. With the
	// columns and ZQUANTIZ = 'NONE' (cfitsio NO_QUANTIZE) the 4-byte tile
	// integers are IEEE floats, which cfitsio still scales by the tile's
	// ZSCALE/ZZERO (comp_fz_f32_rice_zqnone).
	Quantized bool
	// ZScaleCols: a float image with the per-tile ZSCALE/ZZERO columns.
	// Every tile then leaves decompression in physical units (the image
	// BSCALE/BZERO folded into the tile scaling, as cfitsio does).
	ZScaleCols bool

	// Null value: the ZBLANK keyword, else the BLANK keyword (a per-tile
	// ZBLANK column, when present, overrides both per tile). BlankKey
	// names the keyword it came from.
	Blank    int64
	BlankSet bool
	BlankKey string

	// BSCALE/BZERO of the uncompressed image, from the BSCALE/BZERO
	// keywords of the compressed header. They apply after decompression.
	BScale float64
	BZero  float64
	// scalingErr: the scaling-card rule refused the BSCALE / BZERO /
	// BLANK cards (pixel reads and BSCALE / BZERO return it).
	scalingErr error

	// Whole-image ZSCALE/ZZERO keywords of an unquantized float image (no
	// ZSCALE/ZZERO columns). When set they replace BSCALE/BZERO, as in
	// cfitsio.
	FloatScale    float64
	FloatZero     float64
	HasFloatScale bool
}

// Type returns TypeImage — the compressed HDU presents an image-shaped
// view to callers regardless of its physical binary-table backing.
func (h *CompressedImageHDU) Type() HDUType { return TypeImage }

// Index returns the 0-based HDU index within the parent file.
func (h *CompressedImageHDU) Index() int { return h.rec.index }

// Header returns the parsed header of the underlying binary table.
// Callers see every keyword on the HDU, including the Z* compression
// keywords. No logical (uncompressed-image) header is reconstructed; the
// logical image parameters are available from BITPIX, NAXIS, Shape,
// BSCALE and BZERO.
//
// If the header does not parse, Header returns an empty header (it does
// not panic); HeaderOf and Validate return the parse error, and the pixel
// readers fail with it (TestCompressedHeaderParseErrorNoPanic).
func (h *CompressedImageHDU) Header() *header.Header {
	hdr, err := h.rec.loadHeader()
	if err != nil {
		return header.New()
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

// BSCALE returns the BSCALE of the uncompressed image: the BSCALE keyword
// of the compressed header, as cfitsio writes and reads it, by the
// scaling-card rule of every image read (ImageHDU.BSCALE; on a
// tile-compressed header the first of duplicate cards, as cfitsio and
// astropy read it). 1 when there is none. The error is the HDU's metadata
// error (a refused card wraps ErrInvalidScaling).
func (h *CompressedImageHDU) BSCALE() (float64, error) {
	m, err := h.pixelMetadata()
	if err != nil {
		return 0, err
	}
	return m.BScale, nil
}

// BZERO returns the BZERO of the uncompressed image (the BZERO keyword of
// the compressed header, read like BSCALE; e.g. 32768 for an fpack'd
// uint16 image). 0 when there is none.
func (h *CompressedImageHDU) BZERO() (float64, error) {
	m, err := h.pixelMetadata()
	if err != nil {
		return 0, err
	}
	return m.BZero, nil
}

// Compressed always returns true for a CompressedImageHDU.
func (h *CompressedImageHDU) Compressed() bool { return true }

// CompressionType returns the name of the algorithm the ZCMPTYPE card
// selects, as cfitsio recognises it (case-sensitive): "RICE_1" (also for
// cfitsio's alias 'RICE_ONE'), "GZIP_1", "GZIP_2", "HCOMPRESS_1",
// "PLIO_1", "NOCOMPRESS", or "UNKNOWN" for any other value (the pixels of
// such an HDU cannot be read). The card itself is in Header().
func (h *CompressedImageHDU) CompressionType() string {
	m, _ := h.metadata()
	if m == nil {
		return ""
	}
	return m.CmpType.String()
}

// Validate returns the error that makes the HDU's Z* keywords unreadable
// (a header that does not parse, a ZNAXIS outside 1..6, a ZNAXISn below 1,
// a ZBITPIX that is no BITPIX, a missing ZBITPIX or ZCMPTYPE parameter, an
// unknown ZQUANTIZ, ...), or nil. BITPIX, NAXIS and Shape return zero
// values for such an HDU; the pixel readers return this error. BSCALE /
// BZERO / BLANK cards the scaling-card rule refuses (ImageHDU.BSCALE) are
// not a structural error: they leave Validate nil and BITPIX, NAXIS and
// Shape readable, and fail the pixel reads and BSCALE / BZERO.
func (h *CompressedImageHDU) Validate() error {
	_, err := h.metadata()
	return err
}

// pixelMetadata is metadata for the pixel readers: it also fails when the
// scaling-card rule refused the BSCALE / BZERO / BLANK cards.
func (h *CompressedImageHDU) pixelMetadata() (*compressedMetadata, error) {
	m, err := h.metadata()
	if err != nil {
		return nil, err
	}
	if m.scalingErr != nil {
		return nil, m.scalingErr
	}
	return m, nil
}

// metadata returns the parsed Z* keyword set, loading it lazily on first
// access.
func (h *CompressedImageHDU) metadata() (*compressedMetadata, error) {
	h.metaOnce.Do(func() {
		hdr, err := h.rec.loadHeader()
		if err != nil {
			h.metaErr = fmt.Errorf("fits: CompressedImageHDU: header: %w", err)
			return
		}
		h.meta, h.metaErr = parseCompressedMetadata(hdr)
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
		// cfitsio does not check ZBITPIX here and decodes the tiles of a
		// ZBITPIX 12 or 0 image anyway (comp_bad/u16_rice_zbitpix12,
		// _zbitpix0 goldens); with no pixel type there is nothing to prove a
		// reading against: refused.
		if !m.Bitpix.Valid() {
			return nil, fmt.Errorf("fits: CompressedImageHDU: ZBITPIX %d is not a valid BITPIX (8, 16, 32, 64, -32, -64)", v)
		}
	} else {
		return nil, fmt.Errorf("fits: CompressedImageHDU: missing ZBITPIX: %w", err)
	}
	if v, err := hdr.Int("ZNAXIS"); err == nil {
		// cfitsio refuses ZNAXIS < 1 ("Compressed image has no data") and
		// ZNAXIS > MAX_COMPRESS_DIM = 6 with BAD_NAXIS
		// (imcompress.c:5370-5379); so do we (comp_bad/u16_rice_znaxis0,
		// u16_rice_znaxism1, gzip1_7d).
		if v < 1 || v > 6 {
			return nil, fmt.Errorf("fits: CompressedImageHDU: ZNAXIS %d is outside 1..6 (cfitsio: BAD_NAXIS)", v)
		}
		m.Naxis = int(v)
	} else {
		return nil, fmt.Errorf("fits: CompressedImageHDU: missing ZNAXIS: %w", err)
	}
	m.Shape = make([]int64, m.Naxis)
	m.Tile = make([]int64, m.Naxis)
	for i := 1; i <= m.Naxis; i++ {
		if v, err := hdr.Int("ZNAXIS" + strconv.Itoa(i)); err == nil {
			// A negative or zero axis length: cfitsio refuses negative ones
			// (status 414 even when the tile count matches the rows,
			// comp_bad/u16_rice_znaxisneg); a zero one has no tiles.
			if v < 1 {
				return nil, fmt.Errorf("fits: CompressedImageHDU: ZNAXIS%d = %d (the image has no pixels on that axis)", i, v)
			}
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
	// Algorithm parameters. cfitsio (imcomp_get_compressed_image_par)
	// takes them by POSITION, whatever ZNAMEn says: RICE_1 BLOCKSIZE is
	// ZVAL1 (required) and BYTEPIX is ZVAL2 (4 when absent or when ZNAME2
	// is the pre-2009 'NOISEBIT'; the two are swapped back when BLOCKSIZE
	// < 16 and BYTEPIX > 8); HCOMPRESS_1 SCALE is ZVAL1 (required) and
	// SMOOTH is ZVAL2 (0 when absent). Other algorithms take none.
	switch m.CmpType {
	case compress.RICE1:
		bs, err := hdr.Int("ZVAL1")
		if err != nil {
			return nil, fmt.Errorf("fits: CompressedImageHDU: RICE_1 needs ZVAL1 (BLOCKSIZE): %w", err)
		}
		bpx := int64(4)
		if n2, err := hdr.String("ZNAME2"); err != nil || strings.TrimSpace(n2) != "NOISEBIT" {
			if v, err := hdr.Int("ZVAL2"); err == nil {
				bpx = v
			}
		}
		if bs < 16 && bpx > 8 {
			bs, bpx = bpx, bs
		}
		m.CmpParams["BLOCKSIZE"] = bs
		m.CmpParams["BYTEPIX"] = bpx
	case compress.HCOMPRESS1:
		// SCALE itself is not needed to decode (the stream carries it).
		if _, err := hdr.Float("ZVAL1"); err != nil {
			return nil, fmt.Errorf("fits: CompressedImageHDU: HCOMPRESS_1 needs ZVAL1 (SCALE): %w", err)
		}
		if v, err := hdr.Int("ZVAL2"); err == nil {
			m.CmpParams["SMOOTH"] = v
		}
	}
	if s, err := hdr.String("ZQUANTIZ"); err == nil {
		m.Quantize = strings.TrimSpace(s)
	}
	hasZScaleCol, hasZZeroCol := hasColumn(hdr, "ZSCALE"), hasColumn(hdr, "ZZERO")
	// cfitsio: a float image with neither a ZSCALE nor a ZZERO column is
	// not quantized, whatever ZQUANTIZ says (it is not even parsed).
	m.ZScaleCols = m.Bitpix.IsFloat() && (hasZScaleCol || hasZZeroCol)
	m.Quantized = m.ZScaleCols && m.Quantize != "NONE"
	if !m.Bitpix.IsFloat() || hasZScaleCol || hasZZeroCol {
		switch m.Quantize {
		case "", "NO_DITHER", "NONE":
			m.Dither = compress.NoDither
		case "SUBTRACTIVE_DITHER_1", "SUBTRACTIVE_DITHER_2":
			m.Dither = compress.ParseDitherMethod(m.Quantize)
		default:
			// cfitsio refuses unknown quantization types; so do we rather
			// than guess.
			return nil, fmt.Errorf("fits: CompressedImageHDU: unknown ZQUANTIZ %q", m.Quantize)
		}
	}
	m.Dither0 = 1
	if v, err := hdr.Int("ZDITHER0"); err == nil {
		m.Dither0 = v
	}
	// BSCALE, BZERO and (integer images) BLANK of the compressed header by
	// the scaling-card rule of every image read (scalingCards), the FIRST
	// card winning as cfitsio (imcomp_get_compressed_image_par) and
	// astropy read a tile-compressed header.
	// A refused card leaves the structure readable (Shape, BITPIX, ...)
	// and fails the pixel reads and BSCALE / BZERO (pixelMetadata).
	sc, err := scalingCards(encodedCards(hdr, nil), m.Bitpix, true)
	if err != nil {
		m.scalingErr = fmt.Errorf("fits: CompressedImageHDU: %w", err)
	}
	m.BScale, m.BZero = sc.bscale, sc.bzero
	if err != nil {
		m.BScale = 1
	}
	if v, err := hdr.Int("ZBLANK"); err == nil {
		m.Blank, m.BlankSet, m.BlankKey = v, true, "ZBLANK"
	} else if sc.haveBlank {
		m.Blank, m.BlankSet, m.BlankKey = sc.blank, true, "BLANK"
	} else if v, err := hdr.Int("BLANK"); err == nil && m.Bitpix.IsFloat() {
		// A float image's BLANK is not a scaling card (FITS 4.0
		// sect. 4.4.2.4), but cfitsio takes it as the null value of the
		// quantized integers when there is no ZBLANK.
		m.Blank, m.BlankSet, m.BlankKey = v, true, "BLANK"
	}
	_, zsErr := hdr.Float("ZSCALE")
	_, zzErr := hdr.Float("ZZERO")
	hasZScaleKw, hasZZeroKw := zsErr == nil, zzErr == nil
	switch {
	case !m.Bitpix.IsFloat():
		// cfitsio would use a ZSCALE keyword or column in place of
		// BSCALE/BZERO on an integer image. No writer makes one and there
		// is no test file, so refuse it.
		if hasZScaleCol || hasZZeroCol || hasZScaleKw || hasZZeroKw {
			return nil, fmt.Errorf("fits: CompressedImageHDU: ZSCALE/ZZERO on an integer image (ZBITPIX %d) are not supported", int(m.Bitpix))
		}
	case hasZScaleCol != hasZZeroCol:
		// cfitsio's handling of one column without the other is
		// inconsistent (it mixes the column with the keywords, or fails);
		// no writer makes such a table.
		return nil, fmt.Errorf("fits: CompressedImageHDU: quantized float image needs both ZSCALE and ZZERO columns")
	case !hasZScaleCol && hasZScaleKw:
		// Unquantized float with whole-image ZSCALE/ZZERO keywords: cfitsio
		// applies them to the stored floats in place of BSCALE/BZERO
		// (imcomp_decompress_tile, cn_zscale == -1), falling back to
		// BSCALE/BZERO when they are 1 and 0.
		if !hasZZeroKw {
			return nil, fmt.Errorf("fits: CompressedImageHDU: ZSCALE keyword without ZZERO is not supported")
		}
		zs, _ := hdr.Float("ZSCALE")
		zz, _ := hdr.Float("ZZERO")
		if zs != 1 || zz != 0 {
			m.FloatScale, m.FloatZero, m.HasFloatScale = zs, zz, true
		}
	}
	return m, nil
}

// hasColumn reports whether a binary-table header declares a column whose
// TTYPEn is name (case-insensitive, as cfitsio's CASEINSEN lookup).
func hasColumn(hdr *header.Header, name string) bool {
	n, err := hdr.Int("TFIELDS")
	if err != nil {
		return false
	}
	for i := 1; i <= int(n); i++ {
		if s, err := hdr.String("TTYPE" + strconv.Itoa(i)); err == nil && equalFold(strings.TrimSpace(s), name) {
			return true
		}
	}
	return false
}
