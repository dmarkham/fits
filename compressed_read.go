package fits

import (
	"context"
	"fmt"
	"math"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// ReadPixelsCompressed reads and decompresses every tile of h, returning
// a flat slice of pixels in row-major order matching the logical image
// shape (Shape()). BSCALE/BZERO from the Z-prefixed keywords are applied
// as with a regular ImageHDU.
//
// This is the compressed counterpart to ReadPixels[T] for regular
// *ImageHDU. Users who open a file via fits.Open do not need to call
// this directly — they can continue to use the uniform ReadPixels[T]
// form and the library dispatches to the compressed path automatically.
func ReadPixelsCompressed[T Numeric](h *CompressedImageHDU) ([]T, error) {
	return ReadPixelsCompressedContext[T](context.Background(), h)
}

// ReadPixelsCompressedContext is ReadPixelsCompressed with an explicit
// context.
func ReadPixelsCompressedContext[T Numeric](ctx context.Context, h *CompressedImageHDU) ([]T, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	meta, err := h.metadata()
	if err != nil {
		return nil, err
	}
	bp := meta.Bitpix
	if !canConvertInto[T](bp) {
		var zero T
		return nil, &TypeMismatchError{
			Requested: fmt.Sprintf("%T", zero),
			BITPIX:    int(bp),
			Lossy:     true,
		}
	}
	// Decompress every tile into a flat byte buffer with the image layout.
	raw, err := decompressAllTiles(h, meta)
	if err != nil {
		return nil, err
	}
	// npix = product of the logical shape.
	var npix int64 = 1
	for _, n := range meta.Shape {
		npix *= n
	}
	out := make([]T, npix)
	bscale := 1.0
	bzero := 0.0
	if meta.BScale != 0 {
		bscale = meta.BScale
	}
	bzero = meta.BZero
	applyBScale := !(bscale == 1.0 && bzero == 0.0)
	if err := decodePixels(raw, out, bp, applyBScale, bscale, bzero); err != nil {
		return nil, err
	}
	return out, nil
}

// decompressAllTiles walks the tile grid, decompresses each tile via the
// selected algorithm, optionally applies per-tile quantization for float
// images, and stitches the tiles back into a full image byte buffer.
//
// The returned buffer has length npix * elemSize bytes, in big-endian
// FITS on-disk order so it can be fed directly into the same decodePixels
// path used by regular ImageHDU reads.
func decompressAllTiles(h *CompressedImageHDU, meta *compressedMetadata) ([]byte, error) {
	elemSize := meta.Bitpix.Size()
	var npix int64 = 1
	for _, n := range meta.Shape {
		npix *= n
	}
	out := make([]byte, npix*int64(elemSize))

	// Compute the tile grid. tileCounts[i] = ceil(shape[i] / tile[i]).
	tileCounts := make([]int64, meta.Naxis)
	totalTiles := int64(1)
	for i := 0; i < meta.Naxis; i++ {
		tc := (meta.Shape[i] + meta.Tile[i] - 1) / meta.Tile[i]
		tileCounts[i] = tc
		totalTiles *= tc
	}

	// Resolve the column indices for the per-row compressed payloads.
	cols, err := h.tbl.Columns()
	if err != nil {
		return nil, fmt.Errorf("fits: CompressedImageHDU: parse columns: %w", err)
	}
	compIdx := findCompColumn(cols, "COMPRESSED_DATA")
	uncompIdx := findCompColumn(cols, "UNCOMPRESSED_DATA")
	gzipIdx := findCompColumn(cols, "GZIP_COMPRESSED_DATA")
	zscaleIdx := findCompColumn(cols, "ZSCALE")
	zzeroIdx := findCompColumn(cols, "ZZERO")
	if compIdx == 0 {
		return nil, fmt.Errorf("fits: CompressedImageHDU: missing COMPRESSED_DATA column")
	}

	nrows := h.tbl.NumRows()
	if nrows != totalTiles {
		return nil, fmt.Errorf("fits: CompressedImageHDU: tile count mismatch: %d rows vs %d tiles", nrows, totalTiles)
	}

	// Read per-tile payloads via the variable-length column reader. For
	// files where astropy stores uncompressed or gzip fallback data on a
	// per-tile basis, we also need those columns.
	compPayloads, err := readVarBytes(h.tbl, compIdx, nrows)
	if err != nil {
		return nil, fmt.Errorf("fits: CompressedImageHDU: read COMPRESSED_DATA: %w", err)
	}
	var uncompPayloads, gzipPayloads [][]byte
	if uncompIdx != 0 {
		uncompPayloads, _ = readVarBytes(h.tbl, uncompIdx, nrows)
	}
	if gzipIdx != 0 {
		gzipPayloads, _ = readVarBytes(h.tbl, gzipIdx, nrows)
	}
	// Per-tile ZSCALE/ZZERO columns for float data.
	var zscaleVals, zzeroVals []float64
	if zscaleIdx != 0 && zzeroIdx != 0 {
		zscaleVals, _ = ReadColumn[float64](h.tbl, zscaleIdx)
		zzeroVals, _ = ReadColumn[float64](h.tbl, zzeroIdx)
	}

	// Decoder for the primary algorithm.
	decoder, err := compress.Select(meta.CmpType, meta.CmpParams)
	if err != nil {
		return nil, fmt.Errorf("fits: CompressedImageHDU: %w", err)
	}
	gzipDecoder, _ := compress.Select(compress.GZIP1, nil)

	// Iterate tiles. For each, decompress into a per-tile buffer sized to
	// that tile's pixel count, then copy into the correct region of the
	// full-image buffer.
	tileIdx := make([]int64, meta.Naxis)
	for t := int64(0); t < totalTiles; t++ {
		// Determine this tile's box in the logical image.
		box := tileBox(meta.Shape, meta.Tile, tileIdx)
		tileNpix := int64(1)
		for _, n := range box.size {
			tileNpix *= n
		}
		tileBytes := make([]byte, tileNpix*int64(elemSize))

		// Select which payload and decoder to use.
		var payload []byte
		var dec compress.Decoder
		if uncompPayloads != nil && len(uncompPayloads[t]) > 0 {
			payload = uncompPayloads[t]
			dec, _ = compress.Select(compress.NoCompress, nil)
		} else if len(compPayloads[t]) > 0 {
			payload = compPayloads[t]
			dec = decoder
		} else if gzipPayloads != nil && len(gzipPayloads[t]) > 0 {
			payload = gzipPayloads[t]
			dec = gzipDecoder
		} else {
			return nil, fmt.Errorf("fits: CompressedImageHDU: tile %d has no payload", t)
		}

		if err := dec.Decode(payload, tileBytes, int(tileNpix), elemSize); err != nil {
			return nil, fmt.Errorf("fits: CompressedImageHDU: tile %d: %w", t, err)
		}

		// For float BITPIX: apply per-tile ZSCALE/ZZERO dequantization.
		// The integer bytes in tileBytes become float bytes at the same
		// elemSize: ZBITPIX=-32 with BYTEPIX=4 means the underlying
		// decompressor produces int32 values which we reinterpret as
		// quantized integers, multiply by ZSCALE, add ZZERO, and write
		// back as float32 big-endian.
		if bp := meta.Bitpix; bp.IsFloat() {
			var zscale, zzero float64
			if zscaleVals != nil && int(t) < len(zscaleVals) {
				zscale = zscaleVals[t]
				zzero = zzeroVals[t]
			} else if meta.HasZScale {
				zscale = meta.ZScaleFixed
				zzero = meta.ZZeroFixed
			} else {
				zscale = 1
				zzero = 0
			}
			// Apply dither (and BLANK → NaN) per ZQUANTIZ.
			dither := compress.ParseDitherMethod(meta.Quantize)
			dequantizeTile(tileBytes, bp, int(tileNpix), zscale, zzero, meta.Blank, meta.BlankSet, dither, int(t+1), int(meta.Dither0))
		}

		// Copy tileBytes into the correct region of the output.
		copyTileIntoImage(out, tileBytes, meta.Shape, box, elemSize)

		// Advance multi-dim tile counter.
		for d := 0; d < meta.Naxis; d++ {
			tileIdx[d]++
			if tileIdx[d] < tileCounts[d] {
				break
			}
			tileIdx[d] = 0
		}
	}
	return out, nil
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

// readVarBytes reads a VLA column and returns one []byte per row containing
// the raw on-disk tile payload. Supports both 1PB (byte VLA) and 1PI
// (int16 VLA) element types — PLIO_1 uses 1PI to store its 16-bit token
// stream, while RICE/GZIP/NOCOMPRESS use 1PB.
func readVarBytes(t *BinaryTableHDU, col int, nrows int64) ([][]byte, error) {
	cols, _ := t.Columns()
	elemCode := byte(cols[col-1].bin.VarType)
	out := make([][]byte, nrows)
	switch elemCode {
	case 'B':
		vc, err := ReadVarColumn[uint8](t, col)
		if err != nil {
			return nil, err
		}
		for i := int64(0); i < nrows; i++ {
			row := vc.At(int(i))
			out[i] = append([]byte(nil), row...)
		}
	case 'I':
		vc, err := ReadVarColumn[int16](t, col)
		if err != nil {
			return nil, err
		}
		for i := int64(0); i < nrows; i++ {
			row := vc.At(int(i))
			buf := make([]byte, len(row)*2)
			for k, v := range row {
				buf[k*2] = byte(uint16(v) >> 8)
				buf[k*2+1] = byte(uint16(v))
			}
			out[i] = buf
		}
	default:
		return nil, fmt.Errorf("fits: compressed column %d has unsupported VLA element type %c", col, elemCode)
	}
	return out, nil
}

// dequantizeTile reinterprets a tile of quantized integer bytes as
// floating-point values using the per-tile ZSCALE/ZZERO parameters, and
// writes the results back into buf as big-endian floats. For BITPIX=-32
// the intermediate integer is read from 4 bytes per pixel; for -64 from
// 8 bytes. NULL pixels (value == ZBLANK) are written as NaN.
//
// When the image was compressed with SUBTRACTIVE_DITHER_1 or
// SUBTRACTIVE_DITHER_2 quantization, the per-pixel random dither offset
// is applied during reconstruction. tileRow is the 1-based tile index
// within the binary table; zdither0 is the ZDITHER0 seed from the
// header (typically 1).
func dequantizeTile(buf []byte, bp bitpix.BITPIX, nelem int, zscale, zzero float64, blank int64, blankSet bool, dither compress.DitherMethod, tileRow, zdither0 int) {
	// Decode the big-endian integer tile bytes into an int32 scratch
	// slice, dequantize via the dither-aware path, and write the float
	// result back into buf. Dither state must thread pixel-by-pixel
	// which is why we go through int32 intermediate storage rather than
	// inlining as in the original (non-dither) version.
	switch bp {
	case bitpix.Float32:
		ints := make([]int32, nelem)
		for i := 0; i < nelem; i++ {
			ints[i] = bigendian.Int32(buf[i*4:])
		}
		out := make([]float32, nelem)
		compress.DequantizeFloat32(ints, out, zscale, zzero, dither, tileRow, zdither0)
		for i := 0; i < nelem; i++ {
			if blankSet && int64(ints[i]) == blank {
				bigendian.PutFloat32(buf[i*4:], float32(math.NaN()))
			} else {
				bigendian.PutFloat32(buf[i*4:], out[i])
			}
		}
	case bitpix.Float64:
		ints := make([]int32, nelem)
		for i := 0; i < nelem; i++ {
			// cfitsio stores BITPIX=-64 quantized tiles as int32 still
			// (it's the same on-disk integer width after quantization —
			// the 64-bit BITPIX is recovered via ZSCALE/ZZERO to float64).
			ints[i] = int32(bigendian.Int64(buf[i*8:]))
		}
		out := make([]float64, nelem)
		compress.Dequantize(ints, out, zscale, zzero, dither, tileRow, zdither0)
		for i := 0; i < nelem; i++ {
			if blankSet && int64(ints[i]) == blank {
				bigendian.PutFloat64(buf[i*8:], math.NaN())
			} else {
				bigendian.PutFloat64(buf[i*8:], out[i])
			}
		}
	}
}
