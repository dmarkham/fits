package fits

import (
	"context"
	"fmt"
	"math"
	"slices"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// ReadPixelsCompressed reads and decompresses every tile of h, returning
// a flat slice of pixels in row-major order matching the logical image
// shape (Shape()). The image's BSCALE/BZERO (stored under those plain
// names in the compressed header, as cfitsio writes and reads them) are
// applied as for a regular ImageHDU, so an fpack'd uint16 image (BZERO =
// 32768) reads back as its uint16 values.
//
// Which T are accepted is the ReadPixels rule with ZBITPIX for BITPIX:
// float32 and float64 accept every image; an integer T is accepted only
// for an integer ZBITPIX with integral BSCALE and BZERO, and only when T
// holds the physical value of every possible stored value. Otherwise the
// read returns a *TypeMismatchError (wrapping ErrTypeMismatch) instead of
// wrapped or truncated values: ReadPixelsCompressed[int16] on an fpack'd
// uint16 image is refused. Accepted integer reads are exact (no float64
// on the way). Proven over every corpus file and every Go integer type by
// TestReadPixelsCompressedIntegerTargets.
//
// Null pixels: float images return NaN for null pixels (quantized value
// equal to the tile's ZBLANK, or a NaN stored losslessly). Integer images
// return the stored null value, scaled, exactly like ReadPixels on an
// *ImageHDU (TestCompressedUnmaskedNullsAreStoredValues); use
// ReadPixelsCompressedMasked to identify them. A null pixel whose stored
// tile value the ZBITPIX type cannot hold (4-byte tiles of a 16-bit image)
// makes this read fail, as cfitsio's read as the image type and funpack
// do (status 412); ReadPixelsCompressedMasked and ReadFloat32 read such an
// image (comp_bad/i16_gzip1_w4_nullout.fits); the error names the first
// such pixel (i16_gzip1_w4_nullout2). A BLANK, ZBLANK or ZBLANK column value
// outside the int32 range is refused (cfitsio reads it as a C int and
// ignores it or fails, astropy 8.0.2 refuses the file); -2^31 and 2^31 - 1
// themselves are null values (comp_fz_i32_gzip1_blankmin, _blankmax,
// _zblankcolmin), and a ZBLANK column wins over a BLANK keyword cfitsio
// could not read (comp_fz_i16_zblankcol_blank2p32).
//
// HCOMPRESS_1 and PLIO_1 integer tiles decode to C ints, as in cfitsio.
// Lossy HCOMPRESS_1 reconstructs pixels of a saturated image slightly
// outside the ZBITPIX range. For a 16-bit image they are clipped to the
// range, as funpack and cfitsio's TSHORT read do (imcompress.c resets the
// overflow for HCOMPRESS_1 after fffi4i2); cfitsio's TDOUBLE fits_read_img
// returns them unclipped, a difference the corpus pins
// (comp_fz_u16sat_hcomp_s2, i16sat_hcomp_s2), and its TUSHORT read of the
// uint16 image refuses them with status 412 (the fffi4u2 branch keeps the
// overflow). For an 8-bit image they are refused, as funpack refuses
// them (status 412). A reconstruction outside the int32 range is refused at
// every ZBITPIX: cfitsio's (int) cast wraps it, so cfitsio and funpack
// return a saturated 32-bit pixel as about -2^31 (comp_bad/i32sat_hcomp_s2,
// i64_hcomp_2p31, i64_hcomp_m2p31m1). A PLIO_1 value the ZBITPIX type
// cannot hold is refused (cfitsio's read as the image type and funpack
// fail with 412; comp_bad/u16_plio_bzero0, i32_plio_zbitpix8,
// i32_plio_zbitpix16_bzero32768).
//
// A BITPIX -32 image yields float32 pixels (what funpack writes, and what
// cfitsio returns for a TFLOAT read); read as float64 they are those
// float32 values widened, not cfitsio's TDOUBLE intermediate, which keeps
// the dequantization in double precision. The exception is an unquantized
// -32 image with ZSCALE/ZZERO keywords: like an uncompressed -32 image
// with BSCALE, its float64 read scales the stored float32 in double
// precision (equal to cfitsio's TDOUBLE read).
//
// PLIO_1 with BZERO = 32768 (uint16): cfitsio stores value + 32768 in the
// PLIO stream and subtracts it on read (imcompress.c:2343, 6477); the
// library does the same. For ZBITPIX 8 or 32 with BZERO = 32768 cfitsio's
// writer and reader disagree, so that combination is refused.
//
// Fallback tiles: a float tile cfitsio could not quantize is stored as raw
// values in GZIP_COMPRESSED_DATA (or, before 2011, UNCOMPRESSED_DATA, which
// may hold integer tiles too). They get the image's BSCALE/BZERO and BLANK
// like every other tile, as funpack and astropy do. libcfitsio's
// fits_read_img differs here: it returns those tiles unscaled and without
// null checking while scaling the rest of the image. The corpus pins that
// difference (rawRows in compressed_corpus_test.go) with the funpack twin
// as the truth for those tiles.
//
// Tile integers of another width: cfitsio picks a GZIP_1/GZIP_2 tile's
// integer type from its inflated length (1, 2, 4 or 8 bytes per pixel)
// and decodes RICE_1 with the BYTEPIX width (1, 2, else 4: a ZVAL2 of 8
// means 4), then converts to the image type (imcompress.c:6335-6417); the
// lowest and highest value of the ZBITPIX type are read (comp_fz_*_edge),
// one past either is refused. Quantized float GZIP tiles of 1 or 2 bytes
// are widened the same way (unquantize_i1r4 / _i2r4). Old
// cfitsio wrote 16-bit images as 4-byte GZIP tiles and RICE_1 without
// ZVAL2 (BYTEPIX 4); astropy ships real examples (M13). The library does
// the same, and a committed sample proves every width on every integer
// ZBITPIX (GZIP 2/4/8 bytes on 8-bit, 1/4/8 on 16-bit, 1/2/8 on 32-bit;
// RICE_1 BYTEPIX 2/4 on 8-bit, 1/4 on 16-bit, 1/2 on 32-bit, 1/2/8 on
// quantized floats; quantized GZIP 1/2). A null value is compared with the
// tile integers, as cfitsio does, after cfitsio's cast of it to the tile
// type (see tileNullValue).
//
// ZQUANTIZ = 'NONE' on a file with ZSCALE/ZZERO columns: cfitsio
// (NO_QUANTIZE) takes the 4-byte tile integers as IEEE float32 (widened for
// -64) and scales them by the tile's ZSCALE/ZZERO; so does the library
// (astropy 8.0.2 dequantizes such a file instead).
//
// This is the compressed counterpart to ReadPixels[T], which takes a
// *ImageHDU: callers that may get either kind of HDU must type-switch (or
// use ReadFloat32, which accepts both).
//
// What is checked, bit-exact against libcfitsio 4.7.0:
//
//   - TestCompressedCorpusRead over the fpack/libcfitsio corpus
//     testdata/comp_fz_* (every RICE_1, GZIP_1, GZIP_2, HCOMPRESS_1 and
//     PLIO_1 cell of ZBITPIX 8/16/32/-32/-64 x the uint16/uint32/int8
//     BZERO conventions x NO_DITHER / SUBTRACTIVE_DITHER_1 and _2 /
//     ZQUANTIZ 'NONE' that cfitsio can write, enforced by
//     TestCompressedCorpusComplete; HCOMPRESS SMOOTH; ZBLANK keyword over
//     BLANK, BLANK, per-tile ZBLANK column; RICE_1 and HCOMPRESS_1
//     parameters from ZVAL1/ZVAL2 by position; ZSCALE/ZZERO keywords on an
//     unquantized float image; ZQUANTIZ 'NONE' with ZSCALE/ZZERO columns;
//     the GZIP_COMPRESSED_DATA fallback for -32 and -64 and
//     UNCOMPRESSED_DATA as 1PB/1PI/1PJ/1PE/1PD, including NaN inside the
//     float fallback tiles (flagged by ReadPixelsCompressedMasked); 3D, 4D
//     and 6D tiles; the 'RICE_ONE' alias; every tile width above; the real
//     M13 files in GZIP_1 with 4-byte tiles, RICE_1 BYTEPIX 4, HCOMPRESS_1
//     and PLIO_1; the -64 twins of the scaled -32 layouts; the dither
//     random-sequence wrap inside one 10200-pixel -64 tile). Where
//     libcfitsio's fits_read_img is not self-consistent the test pins the
//     difference and takes the independent cfitsio path as the truth: the
//     fallback tiles above and a NaN in a tile of an image without ZBLANK
//     (fits_read_img leaves the first unscaled and does not null-check
//     either; funpack is the truth there), 6-D images (fits_read_compressed_img misreads
//     them), and 8-byte integer tiles (its TDOUBLE read returns memory it
//     never wrote, its TLONG read, funpack and astropy the values).
//   - TestReaderDecodesCfitsioWrittenGoldens over the 83 images written
//     by cfitsio's fits_img_compress in compress/testdata/cref/golden/img
//     (including PLIO uint16 and the -32 dither random-sequence wrap
//     inside 10000-pixel tiles), against cfitsio's own decode.
//
// ZBITPIX 64: fpack and fits_img_compress refuse 64-bit images and funpack
// cannot decompress them, but astropy 8.0.2 writes them with GZIP_1,
// GZIP_2 and NOCOMPRESS, and libcfitsio's TLONG fits_read_img reads the
// GZIP ones (its TDOUBLE read returns unwritten memory for 8-byte tiles, its
// uint64 reads overflow). They are read exactly, int64 and uint64 (BZERO
// 2^63) alike, including GZIP tiles of 1, 2 and 4 bytes and BLANK, against cfitsio's
// stored values and astropy's decode (TestCompressedInt64Read). RICE_1
// (BYTEPIX 4), HCOMPRESS_1 and PLIO_1 at ZBITPIX 64, which cfitsio decodes
// into C ints and widens, are read the same way (comp_fz_i64_rice,
// i64_hcomp, i64_plio, written from libcfitsio's own coders).
//
// Refused with an error: a ZBITPIX that is no BITPIX and a ZNAXISn below 1
// (cfitsio reads the first anyway and refuses the second with 414), a
// ZCMPTYPE cfitsio does not recognise (compared
// case-sensitively, as cfitsio does: 'rice_1' is refused;
// TestZCMPTYPEComparedLikeCfitsio), an unknown ZQUANTIZ on a quantized
// image, ZNAXIS outside 1..6 (cfitsio BAD_NAXIS), RICE_1 tiles of an
// unquantized float image other than 4-byte tiles, HCOMPRESS_1/PLIO_1
// tiles of an unquantized float image, quantized GZIP tiles of 8 bytes, a
// BLANK/ZBLANK outside int32, a GZIP tile
// whose inflated length is not 1, 2, 4 or 8 times its pixel count
// (including leftover bytes; cfitsio: status 414), a tile value the
// ZBITPIX type cannot hold (cfitsio fails a read as the image type with
// status 412; it is refused here for every T), a tile with no payload
// (cfitsio status 415), an integer tile in GZIP_COMPRESSED_DATA (cfitsio
// status 414), a table row count other than the tile count (cfitsio will
// not open the HDU), an UNCOMPRESSED_DATA tile of the wrong length
// (cfitsio reads it without an error, leaving pixels undefined), only one
// of the ZSCALE/ZZERO columns, and ZSCALE/ZZERO on an integer image
// (TestCompressedBadSamplesRefused, TestCompressedReadRejects).
// NOCOMPRESS tiles are decoded as raw values; cfitsio 4.7.0 cannot read
// NOCOMPRESS, so astropy 8.0.2 is the truth (TestCompressedFixtures: u8,
// i16, i32 and the quantized f32/f64 astropy writes).
func ReadPixelsCompressed[T Numeric](h *CompressedImageHDU) ([]T, error) {
	return ReadPixelsCompressedContext[T](context.Background(), h)
}

// ReadPixelsCompressedContext is ReadPixelsCompressed with an explicit
// context, consulted before the read and between tiles.
func ReadPixelsCompressedContext[T Numeric](ctx context.Context, h *CompressedImageHDU) ([]T, error) {
	out, _, err := readPixelsCompressed[T](ctx, h, false)
	return out, err
}

// ReadPixelsCompressedMasked is ReadPixelsCompressed plus a per-pixel
// validity mask, the compressed counterpart of ReadPixelsMasked. A pixel
// is invalid (mask[i] == false) when cfitsio would flag it null: for
// integer images the stored value equals the null value of its tile (the
// per-tile ZBLANK column, else the ZBLANK keyword, else the BLANK keyword);
// for quantized float images the quantized value equals that null value;
// for unquantized floats the value is NaN (also in a tile of an image
// without ZBLANK, which cfitsio's fits_read_img does not null-check; the
// funpack decode and the plain-image rule flag it). For PLIO_1 uint16 images
// (BZERO = 32768) the null value is compared with the PLIO-decoded value
// (stored + 32768), as cfitsio's fits_read_img does (funpack writes the
// stored values to a plain image, where the plain-image rule nulls them). Invalid pixels are left
// at T's zero value; BSCALE/BZERO are applied to valid pixels only.
func ReadPixelsCompressedMasked[T Numeric](h *CompressedImageHDU) ([]T, []bool, error) {
	out, null, err := readPixelsCompressed[T](context.Background(), h, true)
	if err != nil {
		return nil, nil, err
	}
	mask := make([]bool, len(out))
	var zero T
	for i := range mask {
		if null != nil && null[i] {
			out[i] = zero
			continue
		}
		mask[i] = true
	}
	return out, mask, nil
}

func readPixelsCompressed[T Numeric](ctx context.Context, h *CompressedImageHDU, wantNull bool) ([]T, []bool, error) {
	if err := ctx.Err(); err != nil {
		return nil, nil, err
	}
	meta, err := h.pixelMetadata()
	if err != nil {
		return nil, nil, err
	}
	bp := meta.Bitpix
	// The plain-image rule (ReadPixels): an integer T only for an integer
	// ZBITPIX with integral BSCALE/BZERO, and only when T holds the
	// physical value of every possible stored value.
	scale := newScaling(meta.BScale, meta.BZero)
	if !canHold[T](bp, scale) {
		var zero T
		return nil, nil, &TypeMismatchError{
			Requested: fmt.Sprintf("%T", zero),
			BITPIX:    int(bp),
			Lossy:     true,
			BSCALE:    meta.BScale,
			BZERO:     meta.BZero,
		}
	}
	// Decompress every tile into a flat byte buffer with the image layout.
	var nullOverflow error
	raw, null, err := decodeTiles(ctx, h, meta, &nullOverflow)
	if err != nil {
		return nil, nil, err
	}
	if nullOverflow != nil && !wantNull {
		return nil, nil, nullOverflow
	}
	var npix int64 = 1
	for _, n := range meta.Shape {
		npix *= n
	}
	out := make([]T, npix)
	// For quantized float tiles the image BSCALE/BZERO were already folded
	// into the per-tile scale (cfitsio imcomp_decompress_tile); they apply
	// here for integer images and losslessly stored floats.
	bscale, bzero := meta.BScale, meta.BZero
	switch {
	case bp.IsFloat() && meta.ZScaleCols:
		bscale, bzero = 1, 0
	case bp.IsFloat() && meta.HasFloatScale:
		bscale, bzero = meta.FloatScale, meta.FloatZero
	}
	if _, _, isInt, _ := intRange[T](); isInt || !bp.IsFloat() {
		// An integer image (canHold refused integer T on float images):
		// through the same converter as ReadPixels, so the compressed
		// conversion is the plain-image one (exact for integer T, and for
		// float T e.g. uint64 with BZERO 2^63 as cfitsio does it).
		newPixConv[T](bp, scale).all(raw, out)
		if !wantNull {
			null = nil
		}
		return out, null, nil
	}
	applyBScale := !(bscale == 1.0 && bzero == 0.0)
	if err := decodePixels(raw, out, bp, applyBScale, bscale, bzero); err != nil {
		return nil, nil, err
	}
	if !wantNull {
		null = nil
	}
	return out, null, nil
}

// decompressAllTiles walks the tile grid, decompresses each tile via the
// selected algorithm, dequantizes float tiles, and stitches the tiles back
// into a full image byte buffer (big-endian, BITPIX layout, unscaled for
// integer images) so it can be fed into the same decodePixels path used by
// regular ImageHDU reads. It also returns the null flags (nil when the
// image has no null pixels).
//
// The per-tile rules follow cfitsio imcompress.c
// (imcomp_get_compressed_image_par, imcomp_decompress_tile):
//   - payload: COMPRESSED_DATA when non-empty, else GZIP_COMPRESSED_DATA
//     (gzipped raw floats), else UNCOMPRESSED_DATA (raw values);
//   - a float image is quantized unless it has neither a ZSCALE nor a
//     ZZERO column, or ZQUANTIZ = 'NONE'; quantized tiles hold 4-byte
//     integers for both BITPIX -32 and -64;
//   - the null value is the tile's ZBLANK column entry, else the ZBLANK
//     keyword, else the BLANK keyword;
//   - dither row = tile number + ZDITHER0 - 1 (ZDITHER0 defaults to 1).
func decompressAllTiles(ctx context.Context, h *CompressedImageHDU, meta *compressedMetadata) ([]byte, []bool, error) {
	return decodeTiles(ctx, h, meta, nil)
}

// decodeTiles is decompressAllTiles. When nullOverflow is not nil it is
// set to the error an unmasked read must report: a null pixel whose stored
// tile value the ZBITPIX type cannot hold (left nil if there is none). The
// masked reads and ReadFloat32 do not need such a pixel's value; an
// unmasked read would have to return it as the image type, which cfitsio's
// read as the image type and funpack refuse (status 412).
func decodeTiles(ctx context.Context, h *CompressedImageHDU, meta *compressedMetadata, nullOverflow *error) ([]byte, []bool, error) {
	if meta.scalingErr != nil {
		return nil, nil, meta.scalingErr
	}
	bp := meta.Bitpix
	elemSize := bp.Size()
	var npix int64 = 1
	for _, n := range meta.Shape {
		npix *= n
	}
	out := make([]byte, npix*int64(elemSize))
	var null []bool // allocated on the first null pixel

	// Compute the tile grid. tileCounts[i] = ceil(shape[i] / tile[i]).
	tileCounts := make([]int64, meta.Naxis)
	totalTiles := int64(1)
	for i := 0; i < meta.Naxis; i++ {
		if meta.Tile[i] <= 0 {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: ZTILE%d = %d", i+1, meta.Tile[i])
		}
		tc := (meta.Shape[i] + meta.Tile[i] - 1) / meta.Tile[i]
		tileCounts[i] = tc
		totalTiles *= tc
	}

	// Resolve the column indices for the per-row compressed payloads.
	cols, err := h.tbl.Columns()
	if err != nil {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: parse columns: %w", err)
	}
	compIdx := findCompColumn(cols, "COMPRESSED_DATA")
	uncompIdx := findCompColumn(cols, "UNCOMPRESSED_DATA")
	gzipIdx := findCompColumn(cols, "GZIP_COMPRESSED_DATA")
	zscaleIdx := findCompColumn(cols, "ZSCALE")
	zzeroIdx := findCompColumn(cols, "ZZERO")
	zblankIdx := findCompColumn(cols, "ZBLANK")
	if compIdx == 0 {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: missing COMPRESSED_DATA column")
	}

	nrows := h.tbl.NumRows()
	if nrows != totalTiles {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile count mismatch: %d rows vs %d tiles", nrows, totalTiles)
	}

	compPayloads, err := readVarBytes(h.tbl, compIdx, nrows)
	if err != nil {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: read COMPRESSED_DATA: %w", err)
	}
	var uncompPayloads, gzipPayloads [][]byte
	if uncompIdx != 0 {
		if uncompPayloads, err = readVarBytes(h.tbl, uncompIdx, nrows); err != nil {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: read UNCOMPRESSED_DATA: %w", err)
		}
	}
	if gzipIdx != 0 {
		if gzipPayloads, err = readVarBytes(h.tbl, gzipIdx, nrows); err != nil {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: read GZIP_COMPRESSED_DATA: %w", err)
		}
	}
	// Per-tile ZSCALE/ZZERO (quantized float) and ZBLANK columns.
	var zscaleVals, zzeroVals []float64
	var zblankVals []int64
	if zscaleIdx != 0 {
		if zscaleVals, err = ReadColumn[float64](h.tbl, zscaleIdx); err != nil {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: read ZSCALE: %w", err)
		}
	}
	if zzeroIdx != 0 {
		if zzeroVals, err = ReadColumn[float64](h.tbl, zzeroIdx); err != nil {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: read ZZERO: %w", err)
		}
	}
	if zblankIdx != 0 {
		if zblankVals, err = ReadColumn[int64](h.tbl, zblankIdx); err != nil {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: read ZBLANK: %w", err)
		}
	}
	// cfitsio reads the null value as a C int (ZBLANK/BLANK with ffgky
	// TINT, the ZBLANK column with ffgcvk TINT, imcompress.c:5534-5549 and
	// 6234-6245). One outside the int32 range overflows: cfitsio then
	// ignores the keyword (falling back from ZBLANK to BLANK) or fails the
	// read (the column), while astropy 8.0.2 refuses the file. Refused:
	// comp_bad/i32_rice_blank2p32, i64_gzip1_blank2p40, i16_rice_zblank2p33,
	// i16_rice_zblankcol1k.
	if meta.BlankSet && zblankVals == nil && (meta.Blank < math.MinInt32 || meta.Blank > math.MaxInt32) {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: %s %d is outside the int32 range (cfitsio reads it as a C int and ignores it, astropy refuses it); not supported", meta.BlankKey, meta.Blank)
	}
	for t, v := range zblankVals {
		if v < math.MinInt32 || v > math.MaxInt32 {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: ZBLANK column value %d is outside the int32 range (cfitsio: status 412)", t, v)
		}
	}
	quantized := bp.IsFloat() && meta.Quantized
	// zcols: per-tile ZSCALE/ZZERO columns, so the image BSCALE/BZERO are
	// folded into every tile (see compressedMetadata.ZScaleCols).
	zcols := bp.IsFloat() && meta.ZScaleCols

	// riceWidth: the width of RICE_1 tile integers. cfitsio decodes RICE_1
	// with fits_rdecomp_byte for BYTEPIX 1 (unsigned bytes), with
	// fits_rdecomp_short for BYTEPIX 2 and with the 4-byte coder for every
	// other ZVAL2 (imcompress.c:6335-6346), then converts the tile to the
	// image type: integer images by value (with cfitsio's overflow check),
	// quantized floats through unquantize_i1r4/_i2r4. An old file with no
	// ZVAL2 (BYTEPIX 4) on a 16-bit image (astropy's m13_rice.fits) holds
	// 4-byte tile integers. Each width is proven per ZBITPIX by a corpus
	// sample (comp_fz_*_bp1/_bp2/_bp4/_bp8, real_m13_rice).
	riceWidth := 4
	if meta.CmpType == compress.RICE1 {
		switch meta.CmpParams["BYTEPIX"] {
		case 1:
			riceWidth = 1
		case 2:
			riceWidth = 2
		}
		// An unquantized float image in RICE_1 (ZQUANTIZ 'NONE' with
		// ZSCALE/ZZERO columns): cfitsio takes 4-byte tile integers as
		// IEEE float32, for -32 and -64 alike (fffr4r4 / fffr4r8;
		// comp_fz_f32_rice_zqnone, comp_fz_f64_rice_zqnone), and
		// narrower ones as floats read past the end of the tile (no
		// defined value): refused.
		if bp.IsFloat() && !quantized && riceWidth != 4 {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: RICE_1 BYTEPIX %d tiles of an unquantized ZBITPIX %d image are not supported", meta.CmpParams["BYTEPIX"], int(bp))
		}
	}
	// An unquantized float image in HCOMPRESS_1 or PLIO_1 (only ZQUANTIZ
	// 'NONE' edited into a quantized file makes one; fpack refuses
	// '-q 0' with them): no sample, no reference: refused.
	if bp.IsFloat() && !quantized && (meta.CmpType == compress.HCOMPRESS1 || meta.CmpType == compress.PLIO1) {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: %v tiles of an unquantized ZBITPIX %d image are not supported", meta.CmpType, int(bp))
	}
	// GZIP_2 un-shuffles by the width of the stored values (cfitsio uses
	// the tile data type, not a header keyword).
	params := make(compress.Params, len(meta.CmpParams))
	for k, v := range meta.CmpParams {
		params[k] = v
	}
	if meta.CmpType == compress.GZIP2 {
		delete(params, "BYTEPIX")
	}
	decoder, err := compress.Select(meta.CmpType, params)
	if err != nil {
		return nil, nil, fmt.Errorf("fits: CompressedImageHDU: %w", err)
	}
	gzipDecoder, _ := compress.Select(compress.GZIP1, nil)
	rawDecoder, _ := compress.Select(compress.NoCompress, nil)
	dither := meta.Dither
	cnScale, cnZero := meta.BScale, meta.BZero

	// PLIO_1 with BZERO = 32768 (see the tile loop). cfitsio's reader
	// applies the -32768 shift whatever ZBITPIX is, but its writer adds
	// +32768 only to 16-bit images, so for ZBITPIX 8 or 32 the two
	// disagree and there is no file to prove either reading: refuse.
	plioU16 := false
	if meta.CmpType == compress.PLIO1 && meta.BZero == 32768 {
		if bp != bitpix.Int16 {
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: PLIO_1 with BZERO = 32768 on ZBITPIX %d is not supported (cfitsio's writer and reader disagree there)", int(bp))
		}
		plioU16 = true
	}

	tileIdx := make([]int64, meta.Naxis)
	for t := int64(0); t < totalTiles; t++ {
		if err := ctx.Err(); err != nil {
			return nil, nil, err
		}
		box := tileBox(meta.Shape, meta.Tile, tileIdx)
		tileNpix := int64(1)
		for _, n := range box.size {
			tileNpix *= n
		}
		tileBytes := make([]byte, tileNpix*int64(elemSize))
		tileNull := make([]bool, tileNpix)
		anyNull := false

		// nullOv: the first stored value of a null pixel that the ZBITPIX
		// type cannot hold (integer tiles of another width, HCOMPRESS_1,
		// PLIO_1); the first one in the image is what an unmasked read
		// reports (comp_bad/i16_gzip1_w4_nullout2).
		var nullOv *int64
		haveBlank, blank := meta.BlankSet, meta.Blank
		if zblankVals != nil {
			haveBlank, blank = true, zblankVals[t]
		}

		switch {
		case len(compPayloads[t]) > 0 && !bp.IsFloat() && (meta.CmpType == compress.HCOMPRESS1 || meta.CmpType == compress.PLIO1):
			// HCOMPRESS_1 and PLIO_1 tiles decode to C ints (cfitsio
			// imcomp_decompress_tile: tiledatatype TINT for every ZBITPIX),
			// compared with the int null value uncast and then converted to
			// the image type (see hcompPlioTile).
			var err error
			if nullOv, err = hcompPlioTile(meta.CmpType, decoder, compPayloads[t], int(tileNpix), tileBytes, bp, plioU16, haveBlank, blank, tileNull); err != nil {
				return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
			}
			anyNull = slices.Contains(tileNull, true)
		case len(compPayloads[t]) > 0:
			if !quantized && !bp.IsFloat() && (meta.CmpType == compress.GZIP1 || meta.CmpType == compress.GZIP2 || (meta.CmpType == compress.RICE1 && riceWidth != elemSize)) {
				// Integer tile whose stored integers may be wider or
				// narrower than ZBITPIX (see riceWidth and readWideTile).
				vals, w, err := readWideTile(meta.CmpType, decoder, compPayloads[t], int(tileNpix), riceWidth)
				if err != nil {
					return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
				}
				// cfitsio compares the null value with the tile integers
				// (the fffi1i2/fffi2i2/fffi4i2/fffi8i2 tnull check), before
				// conversion, with tnull cast to the tile type.
				if haveBlank {
					anyNull = markBlank(vals, widthBitpix(w), tileNullValue(blank, w), tileNull)
				}
				if nullOv, err = narrowTileInts(vals, w, tileBytes, bp, tileNull); err != nil {
					return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
				}
				break
			}
			if !quantized {
				if bp == bitpix.Float64 && meta.CmpType == compress.RICE1 {
					// 4-byte RICE_1 tile integers of a -64 image: IEEE
					// float32, widened (cfitsio fffr4r8).
					f4 := make([]byte, tileNpix*4)
					if err := decoder.Decode(compPayloads[t], f4, int(tileNpix), 4); err != nil {
						return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
					}
					for i := int64(0); i < tileNpix; i++ {
						bigendian.PutFloat64(tileBytes[i*8:], float64(bigendian.Float32(f4[i*4:])))
					}
				} else if err := decoder.Decode(compPayloads[t], tileBytes, int(tileNpix), elemSize); err != nil {
					return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
				}
				if bp.IsFloat() {
					if zcols {
						// ZQUANTIZ 'NONE' with ZSCALE/ZZERO columns: cfitsio
						// scales the stored floats by the tile's ZSCALE/ZZERO
						// (fffr4r4 with the column values).
						zs, zz := tileScale(zscaleVals[t], zzeroVals[t], cnScale, cnZero)
						applyFloatScale(tileBytes, bp, zs, zz)
					}
					anyNull = markNaN(tileBytes, bp, tileNull)
				} else if haveBlank {
					// Native-width RICE_1 tile: tnull cast to the tile type.
					anyNull = markBlank(tileBytes, bp, tileNullValue(blank, elemSize), tileNull)
				}
				break
			}
			// Quantized float: 4-byte tile integers, the RICE_1 BYTEPIX
			// width, or for GZIP the width of the inflated tile (cfitsio
			// unquantize_i1r4/_i2r4 widen 1- and 2-byte ones, and compare
			// tnull cast to the tile type).
			w := 4
			var narrow []byte
			switch meta.CmpType {
			case compress.GZIP1, compress.GZIP2:
				var err error
				if narrow, w, err = compress.InflateTile(meta.CmpType, compPayloads[t], int(tileNpix)); err != nil {
					return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
				}
				if w == 8 {
					// cfitsio has no 8-byte branch for quantized tiles: it
					// returns memory it never wrote, funpack and astropy
					// 8.0.2 zeros (comp_bad/f32_gzip1_q_w8.fits).
					return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: quantized %v tile inflates to 8 bytes per pixel; not supported (cfitsio returns unwritten memory, funpack and astropy zeros)", t, meta.CmpType)
				}
			default:
				if meta.CmpType == compress.RICE1 {
					w = riceWidth
				}
				narrow = make([]byte, tileNpix*int64(w))
				if err := decoder.Decode(compPayloads[t], narrow, int(tileNpix), w); err != nil {
					return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
				}
			}
			ints := widenToInt32(narrow, w)
			// Both columns exist (parseCompressedMetadata refuses one alone).
			zscale, zzero := tileScale(zscaleVals[t], zzeroVals[t], cnScale, cnZero)
			anyNull = dequantizeTile(ints, tileBytes, bp, int(tileNpix), zscale, zzero, tileNullValue(blank, w), haveBlank, dither, int(t+1), int(meta.Dither0), tileNull)
		case gzipPayloads != nil && len(gzipPayloads[t]) > 0:
			// cfitsio's escape hatch for float tiles it could not quantize:
			// the raw IEEE values, gzipped.
			if !bp.IsFloat() {
				return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: GZIP_COMPRESSED_DATA on an integer image (ZBITPIX %d)", t, int(bp))
			}
			if err := gzipDecoder.Decode(gzipPayloads[t], tileBytes, int(tileNpix), elemSize); err != nil {
				return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d GZIP_COMPRESSED_DATA: %w", t, err)
			}
			if zcols && (cnScale != 1 || cnZero != 0) {
				applyFloatScale(tileBytes, bp, cnScale, cnZero)
			}
			anyNull = markNaN(tileBytes, bp, tileNull)
		case uncompPayloads != nil && len(uncompPayloads[t]) > 0:
			// Pre-2011 cfitsio: the raw values of the tile, column-typed.
			p := uncompPayloads[t]
			if len(p) != len(tileBytes) {
				return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d UNCOMPRESSED_DATA has %d bytes, want %d", t, len(p), len(tileBytes))
			}
			if err := rawDecoder.Decode(p, tileBytes, int(tileNpix), elemSize); err != nil {
				return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d UNCOMPRESSED_DATA: %w", t, err)
			}
			if bp.IsFloat() {
				if zcols && (cnScale != 1 || cnZero != 0) {
					applyFloatScale(tileBytes, bp, cnScale, cnZero)
				}
				anyNull = markNaN(tileBytes, bp, tileNull)
			} else if haveBlank {
				anyNull = markBlank(tileBytes, bp, blank, tileNull)
			}
		default:
			return nil, nil, fmt.Errorf("fits: CompressedImageHDU: tile %d has no payload", t)
		}

		if nullOv != nil && nullOverflow != nil && *nullOverflow == nil {
			*nullOverflow = fmt.Errorf("fits: CompressedImageHDU: tile %d: null tile value %d does not fit ZBITPIX %d (cfitsio: status 412 for a read as the image type); ReadPixelsCompressedMasked reads the image", t, *nullOv, int(bp))
		}
		copyTileIntoImage(out, tileBytes, meta.Shape, box, elemSize)
		if anyNull {
			if null == nil {
				null = make([]bool, npix)
			}
			copyTileNulls(null, tileNull, meta.Shape, box)
		}

		// Advance multi-dim tile counter.
		for d := 0; d < meta.Naxis; d++ {
			tileIdx[d]++
			if tileIdx[d] < tileCounts[d] {
				break
			}
			tileIdx[d] = 0
		}
	}
	return out, null, nil
}

// readWideTile decodes one integer GZIP_1/GZIP_2/RICE_1 tile into its
// stored integers and returns them with their width: for GZIP the width
// cfitsio takes from the inflated length (compress.InflateTile: 1, 2, 4
// or 8 bytes per pixel), for RICE_1 the BYTEPIX width riceWidth. Every
// width on every integer ZBITPIX is proven by a corpus sample: 2/4/8-byte
// tiles on 8-bit (u8_gzip2_w2, u8_gzip1_w4, u8_gzip2_w8, u8_rice_bp2,
// u8_rice_bp4), 1/4/8 on 16-bit (i16small_gzip1_w1, i16small_rice_bp1,
// real_m13_gzip, real_m13_rice, u16_gzip2_w4, i16blank_gzip1_w8) and 1/2/8
// on 32-bit (i32small_gzip2_w1, i32small_gzip2_w2, i32small_rice_bp1,
// i32small_rice_bp2, i32_gzip2_w8).
func readWideTile(algo compress.Algorithm, dec compress.Decoder, src []byte, n, riceWidth int) ([]byte, int, error) {
	if algo == compress.RICE1 {
		vals := make([]byte, n*riceWidth)
		if err := dec.Decode(src, vals, n, riceWidth); err != nil {
			return nil, 0, err
		}
		return vals, riceWidth, nil
	}
	return compress.InflateTile(algo, src, n)
}

// tileNullValue is cfitsio's null value for w-byte tile integers: tnull
// is a C int, passed to the 1- and 2-byte converters cast to unsigned
// char and short (imcompress.c imcomp_decompress_tile, fffi1i2 /
// fffi2i2 / unquantize_i1r4 / unquantize_i2r4 ...), and compared as it is
// with 4- and 8-byte tile integers. So ZBLANK 406 on 1-byte quantized
// tiles flags the 150s (comp_fz_f32_rice_bp1_zblank406: a signed cast,
// -106, or none would flag nothing), BLANK 65541 on 1-byte integer tiles
// the fives (comp_fz_i32small_rice_bp1) and BLANK 65536 on 2-byte tiles
// the zeros (comp_fz_i16small_rice_blankcast), while BLANK 65536 on the
// 4-byte ints of HCOMPRESS_1 / PLIO_1 flags nothing
// (comp_fz_i16small_hcomp_blank65536, _plio_blank65536). The null value
// is in the int32 range (decodeTiles refuses any other).
func tileNullValue(blank int64, w int) int64 {
	switch w {
	case 1:
		return int64(uint8(blank))
	case 2:
		return int64(int16(blank))
	}
	return blank
}

// tileScale is the linear scaling cfitsio applies to a tile with ZSCALE /
// ZZERO columns: the column values combined with a non-default image
// BSCALE/BZERO (imcomp_decompress_tile).
func tileScale(zscale, zzero, cnScale, cnZero float64) (float64, float64) {
	if cnScale != 1 || cnZero != 0 {
		return float64(zscale * cnScale), float64(zzero*cnScale) + cnZero
	}
	return zscale, zzero
}

// widenToInt32 returns w-byte big-endian tile integers (unsigned for w =
// 1, signed for 2 and 4) as big-endian int32.
func widenToInt32(vals []byte, w int) []byte {
	if w == 4 {
		return vals
	}
	n := len(vals) / w
	out := make([]byte, n*4)
	for i := 0; i < n; i++ {
		var v int32
		if w == 1 {
			v = int32(vals[i])
		} else {
			v = int32(bigendian.Int16(vals[i*2:]))
		}
		bigendian.PutInt32(out[i*4:], v)
	}
	return out
}

// widthBitpix is the BITPIX whose pixels are w bytes wide (1 = unsigned
// byte, as cfitsio's TBYTE tile type).
func widthBitpix(w int) bitpix.BITPIX {
	switch w {
	case 1:
		return bitpix.Int8
	case 2:
		return bitpix.Int16
	case 4:
		return bitpix.Int32
	}
	return bitpix.Int64
}

// narrowTileInts converts w-byte tile integers (unsigned for w = 1,
// signed otherwise, big-endian) to the big-endian ZBITPIX integers of an
// integer image. A non-null value outside the ZBITPIX range is an error:
// cfitsio's converters (fffi4i2, fffi8i2, ...) report NUM_OVERFLOW (status
// 412) for it and funpack refuses the file
// (comp_bad/i16_gzip1_w4_overflow.fits and the *_lo / *_hi samples one
// past each edge; the edges themselves are comp_fz_*_edge). A null pixel
// outside the range is no error here (cfitsio's fits_read_imgnull flags
// it, comp_bad/i16_gzip1_w4_nullout.fits): its destination value is
// clipped to the range, which no read returns, and the first such stored
// value is returned so that an unmasked read can refuse.
func narrowTileInts(vals []byte, w int, dst []byte, bp bitpix.BITPIX, null []bool) (*int64, error) {
	n := len(vals) / w
	v := make([]int64, n)
	for i := range v {
		switch w {
		case 1:
			v[i] = int64(vals[i])
		case 2:
			v[i] = int64(bigendian.Int16(vals[i*2:]))
		case 4:
			v[i] = int64(bigendian.Int32(vals[i*4:]))
		default:
			v[i] = bigendian.Int64(vals[i*8:])
		}
	}
	return storeTileInts(v, dst, bp, null, false)
}

// hcompPlioTile decodes one HCOMPRESS_1 or PLIO_1 tile of an integer image
// into dst (big-endian ZBITPIX integers), flagging null pixels in null, as
// cfitsio's imcomp_decompress_tile does. Both decode to C ints
// (tiledatatype TINT at every ZBITPIX, 64 included), and the null value is
// compared with those ints uncast (comp_fz_i16small_hcomp_blank65536,
// _plio_blank65536). Then, per algorithm:
//
//   - HCOMPRESS_1: the reconstruction is exact in 64 bits
//     (fits_hdecompress64); a value outside the int32 range is refused,
//     because cfitsio's (int) cast wraps it (a saturated 32-bit pixel reads
//     as about -2^31 in cfitsio and funpack; comp_bad/i32sat_hcomp_s2,
//     i64_hcomp_2p40). Lossy HCOMPRESS reconstructs pixels slightly
//     outside the ZBITPIX range: for a 16-bit image they are clipped to
//     it, as funpack and cfitsio's TSHORT read do (it resets the overflow
//     status for HCOMPRESS_1 after fffi4i2) (comp_fz_i16sat_hcomp_s2,
//     u16sat_hcomp_s2; cfitsio's TDOUBLE read returns them unclipped, its
//     TUSHORT read of the uint16 image fails with status 412); for an 8-bit image they are refused, as
//     funpack refuses them (status 412; comp_bad/u8sat_hcomp_s2).
//   - PLIO_1: a value the ZBITPIX type cannot hold is refused (cfitsio's
//     read as the image type and funpack fail with 412;
//     comp_bad/u16_plio_bzero0, i32_plio_zbitpix8). With plioU16 (BZERO
//     32768 on ZBITPIX 16) the null value is compared with the PLIO value
//     and the stored int16 is that value - 32768
//     (comp_bad/i32_plio_zbitpix16_bzero32768 when it does not fit).
//
// A null pixel the ZBITPIX type cannot hold is no error; its stored value
// is returned (the first one) so that an unmasked read can refuse.
func hcompPlioTile(algo compress.Algorithm, dec compress.Decoder, src []byte, n int, dst []byte, bp bitpix.BITPIX, plioU16, haveBlank bool, blank int64, null []bool) (*int64, error) {
	vals := make([]int64, n)
	if algo == compress.HCOMPRESS1 {
		buf := make([]byte, n*8)
		if err := dec.Decode(src, buf, n, 8); err != nil {
			return nil, err
		}
		for i := range vals {
			v := bigendian.Int64(buf[i*8:])
			if v < math.MinInt32 || v > math.MaxInt32 {
				return nil, fmt.Errorf("HCOMPRESS_1 tile value %d is outside the int32 range of HCOMPRESS_1 tile integers (cfitsio's fits_hdecompress64 casts it to a C int, wrapping it); not supported", v)
			}
			vals[i] = v
		}
	} else {
		buf := make([]byte, n*4)
		if err := dec.Decode(src, buf, n, 4); err != nil {
			return nil, err
		}
		for i := range vals {
			vals[i] = int64(bigendian.Int32(buf[i*4:]))
		}
	}
	if haveBlank {
		for i, v := range vals {
			if v == blank {
				null[i] = true
			}
		}
	}
	if plioU16 {
		for i := range vals {
			vals[i] -= 32768
		}
	}
	return storeTileInts(vals, dst, bp, null, algo == compress.HCOMPRESS1 && bp == bitpix.Int16)
}

// storeTileInts writes integer tile values as the big-endian ZBITPIX
// integers of an integer image. A non-null value outside the ZBITPIX range
// is an error (cfitsio: status 412), or, with clip, is clipped to the range
// (16-bit HCOMPRESS_1, see hcompPlioTile). A null value outside the range
// is clipped (no read returns it) and the first one is returned.
func storeTileInts(vals []int64, dst []byte, bp bitpix.BITPIX, null []bool, clip bool) (*int64, error) {
	var lo, hi int64
	switch bp {
	case bitpix.Int8:
		lo, hi = 0, math.MaxUint8
	case bitpix.Int16:
		lo, hi = math.MinInt16, math.MaxInt16
	case bitpix.Int32:
		lo, hi = math.MinInt32, math.MaxInt32
	case bitpix.Int64:
		lo, hi = math.MinInt64, math.MaxInt64
	default:
		return nil, fmt.Errorf("ZBITPIX %d tile conversion is not supported", int(bp))
	}
	var overflow *int64
	es := bp.Size()
	for i, v := range vals {
		if v < lo || v > hi {
			if !null[i] && !clip {
				return nil, fmt.Errorf("tile value %d does not fit ZBITPIX %d (cfitsio: status 412, numerical overflow)", v, int(bp))
			}
			if null[i] && overflow == nil {
				ov := v
				overflow = &ov
			}
			v = min(max(v, lo), hi)
		}
		switch es {
		case 1:
			dst[i] = byte(v)
		case 2:
			bigendian.PutInt16(dst[i*2:], int16(v))
		case 4:
			bigendian.PutInt32(dst[i*4:], int32(v))
		default:
			bigendian.PutInt64(dst[i*8:], v)
		}
	}
	return overflow, nil
}

// markNaN flags NaN values of a float tile as null.
func markNaN(buf []byte, bp bitpix.BITPIX, null []bool) bool {
	any := false
	for i := range null {
		var nan bool
		if bp == bitpix.Float32 {
			f := bigendian.Float32(buf[i*4:])
			nan = f != f
		} else {
			nan = math.IsNaN(bigendian.Float64(buf[i*8:]))
		}
		if nan {
			null[i], any = true, true
		}
	}
	return any
}

// markBlank flags integer pixels whose stored value equals blank.
func markBlank(buf []byte, bp bitpix.BITPIX, blank int64, null []bool) bool {
	any := false
	for i := range null {
		var v int64
		switch bp {
		case bitpix.Int8:
			v = int64(buf[i])
		case bitpix.Int16:
			v = int64(bigendian.Int16(buf[i*2:]))
		case bitpix.Int32:
			v = int64(bigendian.Int32(buf[i*4:]))
		case bitpix.Int64:
			v = bigendian.Int64(buf[i*8:])
		}
		if v == blank {
			null[i], any = true, true
		}
	}
	return any
}

// applyFloatScale applies BSCALE/BZERO in place to a float tile (used for
// unquantized tiles of an image whose quantized tiles had the scaling
// folded in, so every tile leaves decompressAllTiles in physical units).
func applyFloatScale(buf []byte, bp bitpix.BITPIX, bscale, bzero float64) {
	n := len(buf) / bp.Size()
	for i := 0; i < n; i++ {
		if bp == bitpix.Float32 {
			f := float64(bigendian.Float32(buf[i*4:]))
			bigendian.PutFloat32(buf[i*4:], float32(float64(f*bscale)+bzero))
		} else {
			f := bigendian.Float64(buf[i*8:])
			bigendian.PutFloat64(buf[i*8:], float64(f*bscale)+bzero)
		}
	}
}

// copyTileNulls places a tile's null flags into the image-wide array.
func copyTileNulls(img, tile []bool, shape []int64, region tileRegion) {
	n := len(shape)
	idx := make([]int64, n)
	rowLen := region.size[0]
	var src int64
	for {
		var off int64
		stride := int64(1)
		for d := 0; d < n; d++ {
			off += (region.origin[d] + idx[d]) * stride
			stride *= shape[d]
		}
		copy(img[off:off+rowLen], tile[src:src+rowLen])
		src += rowLen
		d := 1
		for d < n {
			idx[d]++
			if idx[d] < region.size[d] {
				break
			}
			idx[d] = 0
			d++
		}
		if d >= n {
			break
		}
	}
}

// tileRegion describes one tile's coverage in the logical image grid.
type tileRegion struct {
	origin []int64 // lower corner in image coordinates (per-axis, 0-based)
	size   []int64 // extent along each axis (may be less than meta.Tile for edge tiles)
}

// tileBox computes the origin and extent of the tile at multi-index idx.
func tileBox(shape, tile, idx []int64) tileRegion {
	n := len(shape)
	origin := make([]int64, n)
	size := make([]int64, n)
	for d := 0; d < n; d++ {
		origin[d] = idx[d] * tile[d]
		end := origin[d] + tile[d]
		if end > shape[d] {
			end = shape[d]
		}
		size[d] = end - origin[d]
	}
	return tileRegion{origin: origin, size: size}
}

// copyTileIntoImage places the decoded bytes of one tile into the full
// image buffer at the correct location. Handles arbitrary NAXIS via a
// per-axis offset walk; the copy is done one "row" (innermost axis) at a
// time.
func copyTileIntoImage(imgBuf, tileBuf []byte, shape []int64, region tileRegion, elemSize int) {
	n := len(shape)
	// Stride per axis in the full image, in BYTES.
	imgStride := make([]int64, n)
	var s int64 = int64(elemSize)
	for d := 0; d < n; d++ {
		imgStride[d] = s
		s *= shape[d]
	}
	// Per-tile stride (in bytes) along each axis.
	tileStride := make([]int64, n)
	s = int64(elemSize)
	for d := 0; d < n; d++ {
		tileStride[d] = s
		s *= region.size[d]
	}
	// Walk every "row" (innermost axis) of the tile and copy it to the
	// matching offset in the image buffer.
	// The innermost axis (axis 0 in FITS / first axis stored fastest on
	// disk) is the row. A "row" is region.size[0] pixels long.
	rowBytes := region.size[0] * int64(elemSize)
	idx := make([]int64, n)
	for {
		// Compute source offset inside tileBuf (row start at axes 1..n-1).
		var tileOff int64
		for d := 1; d < n; d++ {
			tileOff += idx[d] * tileStride[d]
		}
		// Compute destination offset inside imgBuf (row start).
		var imgOff int64 = region.origin[0] * imgStride[0]
		for d := 1; d < n; d++ {
			imgOff += (region.origin[d] + idx[d]) * imgStride[d]
		}
		copy(imgBuf[imgOff:imgOff+rowBytes], tileBuf[tileOff:tileOff+rowBytes])
		// Advance index along axes 1..n-1.
		d := 1
		for d < n {
			idx[d]++
			if idx[d] < region.size[d] {
				break
			}
			idx[d] = 0
			d++
		}
		if d == n {
			break
		}
		if n == 1 {
			break // 1D: one row total
		}
	}
}

// findCompColumn returns the 1-based index of a column by name, or 0 if
// absent. Case-insensitive comparison — astropy writes the column names
// as upper-case.
func findCompColumn(cols []Column, name string) int {
	for _, c := range cols {
		if equalFold(c.Name, name) {
			return c.Index
		}
	}
	return 0
}

// equalFold is a tiny case-insensitive ASCII compare to avoid importing
// strings for a single call.
func equalFold(a, b string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := 0; i < len(a); i++ {
		ca, cb := a[i], b[i]
		if 'A' <= ca && ca <= 'Z' {
			ca += 'a' - 'A'
		}
		if 'A' <= cb && cb <= 'Z' {
			cb += 'a' - 'A'
		}
		if ca != cb {
			return false
		}
	}
	return true
}

// readVarBytes reads a VLA column and returns one []byte per row holding
// the row's elements as big-endian bytes (the on-disk heap layout). 1PB is
// the byte stream of RICE/GZIP/HCOMPRESS, 1PI the PLIO_1 token stream.
// UNCOMPRESSED_DATA holds the tile's values in the column type of ZBITPIX:
// 1PB/1PI/1PJ/1PE/1PD are each covered by a corpus file. 1PK (ZBITPIX 64)
// has no sample: no writer puts 64-bit tiles in UNCOMPRESSED_DATA.
func readVarBytes(t *BinaryTableHDU, col int, nrows int64) ([][]byte, error) {
	cols, err := t.Columns()
	if err != nil {
		return nil, err
	}
	switch elemCode := byte(cols[col-1].bin.VarType); elemCode {
	case 'B':
		return varRowsBytes[uint8](t, col, nrows, 1, func(b []byte, v uint8) { b[0] = v })
	case 'I':
		return varRowsBytes[int16](t, col, nrows, 2, func(b []byte, v int16) { bigendian.PutInt16(b, v) })
	case 'J':
		return varRowsBytes[int32](t, col, nrows, 4, func(b []byte, v int32) { bigendian.PutInt32(b, v) })
	case 'K':
		return varRowsBytes[int64](t, col, nrows, 8, func(b []byte, v int64) { bigendian.PutInt64(b, v) })
	case 'E':
		return varRowsBytes[float32](t, col, nrows, 4, func(b []byte, v float32) { bigendian.PutFloat32(b, v) })
	case 'D':
		return varRowsBytes[float64](t, col, nrows, 8, func(b []byte, v float64) { bigendian.PutFloat64(b, v) })
	default:
		return nil, fmt.Errorf("fits: compressed column %d has unsupported VLA element type %c", col, elemCode)
	}
}

func varRowsBytes[T Numeric](t *BinaryTableHDU, col int, nrows int64, size int, put func([]byte, T)) ([][]byte, error) {
	vc, err := ReadVarColumn[T](t, col)
	if err != nil {
		return nil, err
	}
	out := make([][]byte, nrows)
	for i := int64(0); i < nrows; i++ {
		row := vc.At(int(i))
		buf := make([]byte, len(row)*size)
		for k, v := range row {
			put(buf[k*size:], v)
		}
		out[i] = buf
	}
	return out, nil
}

// dequantizeTile converts one tile of quantized big-endian int32 values
// (ints) into floats of the image BITPIX written to buf, using the tile's
// scale/zero and dither method (cfitsio unquantize_i4r4 / _i4r8). Values
// equal to the tile's null value become NaN and are flagged in null;
// tileRow is the 1-based tile number and zdither0 the ZDITHER0 keyword.
// Reports whether any pixel was null.
func dequantizeTile(ints, buf []byte, bp bitpix.BITPIX, nelem int, zscale, zzero float64, blank int64, blankSet bool, dither compress.DitherMethod, tileRow, zdither0 int, null []bool) bool {
	q := make([]int32, nelem)
	for i := range q {
		q[i] = bigendian.Int32(ints[i*4:])
	}
	anyNull := false
	isNull := func(i int) bool {
		if blankSet && int64(q[i]) == blank {
			null[i], anyNull = true, true
			return true
		}
		return false
	}
	switch bp {
	case bitpix.Float32:
		out := make([]float32, nelem)
		compress.DequantizeFloat32(q, out, zscale, zzero, dither, tileRow, zdither0)
		for i := 0; i < nelem; i++ {
			if isNull(i) {
				bigendian.PutFloat32(buf[i*4:], float32(math.NaN()))
			} else {
				bigendian.PutFloat32(buf[i*4:], out[i])
			}
		}
	case bitpix.Float64:
		out := make([]float64, nelem)
		compress.Dequantize(q, out, zscale, zzero, dither, tileRow, zdither0)
		for i := 0; i < nelem; i++ {
			if isNull(i) {
				bigendian.PutFloat64(buf[i*8:], math.NaN())
			} else {
				bigendian.PutFloat64(buf[i*8:], out[i])
			}
		}
	}
	return anyNull
}
