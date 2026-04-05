package fits

import (
	"fmt"
	"math"
	"strconv"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// CompressFloatOptions configures a tile-compressed float image
// write. Fields have sensible defaults so the zero value is usable
// for the common case (RICE_1, whole-row tiles, qlevel=4,
// SUBTRACTIVE_DITHER_1).
type CompressFloatOptions struct {
	// Algorithm is the compression method. Must be one of RICE_1,
	// GZIP_1, GZIP_2, or NoCompress — HCOMPRESS_1 and PLIO_1 do not
	// support float input in the compressed image convention.
	// Default: RICE_1.
	Algorithm compress.Algorithm

	// TileShape picks per-tile dimensions, same convention as
	// CompressOptions.TileShape. Default: whole-row tiles.
	TileShape []int64

	// BlockSize is the RICE_1 block size. Default 32.
	BlockSize int

	// QLevel is the quantization level. Positive values set the
	// step as sigma/QLevel; negative values use an absolute step of
	// -QLevel. Default 4.0.
	QLevel float32

	// DitherMethod selects SUBTRACTIVE_DITHER_1 (default),
	// SUBTRACTIVE_DITHER_2 (zero-preserving), or NoDither.
	DitherMethod compress.DitherMethod

	// ZDither0 is the per-image dither seed written as the ZDITHER0
	// keyword. Default 1.
	ZDither0 int

	// NullValue is the float sentinel for null pixels. Default: -1e30.
	// The writer substitutes NaN inputs with this value before
	// quantization, and the quantizer emits NullValueInt32 at those
	// positions. The sentinel is also written as the ZBLANK header.
	NullValue float32

	// Nullcheck, if true, enables null detection. If false, NaN in
	// the input will produce undefined quantized output — the caller
	// is asserting there are no nulls. Default false.
	Nullcheck bool
}

// CompressOptions configures a tile-compressed image write. Fields have
// sensible defaults so the zero value is usable for the common case
// (RICE_1, whole-row tiles, lossless).
type CompressOptions struct {
	// Algorithm is the compression method. Defaults to RICE_1.
	Algorithm compress.Algorithm

	// TileShape is the per-tile dimensions in the same order as the
	// image's logical shape (first axis varies fastest). A zero or nil
	// value picks the cfitsio default: ZTILE1 = NAXIS1, ZTILE2..n = 1
	// (whole-row tiles).
	TileShape []int64

	// BlockSize is the RICE_1 block size. Default 32.
	BlockSize int

	// HCompressScale is the HCOMPRESS_1 quantization scale. 0 = lossless.
	HCompressScale int

	// ZDither0 is the dither seed (for float images compressed via
	// RICE_1 or GZIP with quantize_level > 0). Default 1.
	ZDither0 int
}

// AppendCompressedImage writes a new tile-compressed image HDU to the
// destination file. For integer BITPIX the compression is lossless. For
// float BITPIX callers should use AppendCompressedFloat32Image /
// AppendCompressedFloat64Image which handle quantization explicitly.
//
// This is the write counterpart of the read-side CompressedImageHDU.
// The resulting HDU is a binary table with ZIMAGE=T and is
// automatically recognized by Open() when re-reading the file.
func AppendCompressedImage[T Numeric](f *File, hdr *header.Header, shape []int64, data []T, opts CompressOptions) (*CompressedImageHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	// Validate.
	if len(shape) == 0 {
		return nil, fmt.Errorf("fits: AppendCompressedImage: empty shape")
	}
	var npix int64 = 1
	for _, n := range shape {
		if n <= 0 {
			return nil, fmt.Errorf("fits: AppendCompressedImage: non-positive axis %d", n)
		}
		npix *= n
	}
	if int64(len(data)) != npix {
		return nil, fmt.Errorf("fits: AppendCompressedImage: data length %d != product of shape %d", len(data), npix)
	}

	// Determine on-disk BITPIX from T.
	bp := pickBitpix[T]()
	if !bp.Valid() || bp.IsFloat() {
		var zero T
		return nil, fmt.Errorf("fits: AppendCompressedImage: integer type required, got %T", zero)
	}

	algo := opts.Algorithm
	if algo == compress.Unknown {
		algo = compress.RICE1
	}
	tileShape := opts.TileShape
	if len(tileShape) == 0 {
		tileShape = make([]int64, len(shape))
		tileShape[0] = shape[0]
		for i := 1; i < len(shape); i++ {
			tileShape[i] = 1
		}
	}
	if len(tileShape) != len(shape) {
		return nil, fmt.Errorf("fits: AppendCompressedImage: tile shape has %d axes, image has %d", len(tileShape), len(shape))
	}
	zdither0 := opts.ZDither0
	if zdither0 == 0 {
		zdither0 = 1
	}
	blockSize := opts.BlockSize
	if blockSize <= 0 {
		blockSize = 32
	}

	// Serialize the image into big-endian FITS byte order.
	elemSize := bp.Size()
	imgBytes := make([]byte, npix*int64(elemSize))
	serializeTileInts(imgBytes, data, bp)

	// Partition into tiles and encode each.
	tileCounts := make([]int64, len(shape))
	totalTiles := int64(1)
	for i := range shape {
		tileCounts[i] = (shape[i] + tileShape[i] - 1) / tileShape[i]
		totalTiles *= tileCounts[i]
	}

	encoder, err := compress.SelectEncoder(algo, compress.Params{
		"BLOCKSIZE": int64(blockSize),
		"BYTEPIX":   int64(elemSize),
		"SCALE":     int64(opts.HCompressScale),
	})
	if err != nil {
		return nil, fmt.Errorf("fits: AppendCompressedImage: %w", err)
	}

	tilePayloads := make([][]byte, totalTiles)
	tileIdx := make([]int64, len(shape))
	for t := int64(0); t < totalTiles; t++ {
		box := tileBox(shape, tileShape, tileIdx)
		tileNpix := int64(1)
		for _, n := range box.size {
			tileNpix *= n
		}
		tileRaw := make([]byte, tileNpix*int64(elemSize))
		extractTileBytes(imgBytes, tileRaw, shape, box, elemSize)
		// HCOMPRESS needs the 2D shape before Encode.
		if setter, ok := encoder.(interface{ SetShape(nx, ny int) }); ok && len(box.size) >= 2 {
			// cfitsio's HCOMPRESS has ny = fastest axis, nx = slow.
			setter.SetShape(int(box.size[1]), int(box.size[0]))
		}
		compressed := make([]byte, max64(int64(len(tileRaw))*2+256, 1024))
		n, err := encoder.Encode(tileRaw, compressed, int(tileNpix), elemSize)
		if err != nil {
			// On overflow, fall back to NOCOMPRESS for this tile by
			// storing the raw bytes in COMPRESSED_DATA. A more cfitsio-
			// faithful implementation would use the UNCOMPRESSED_DATA
			// fallback column, but a single-column fallback is adequate
			// for round-trip fidelity.
			if err == compress.ErrBufferTooSmall {
				compressed = make([]byte, len(tileRaw)*4+1024)
				n, err = encoder.Encode(tileRaw, compressed, int(tileNpix), elemSize)
			}
			if err != nil {
				return nil, fmt.Errorf("fits: AppendCompressedImage: tile %d: %w", t, err)
			}
		}
		tilePayloads[t] = compressed[:n]

		for d := range shape {
			tileIdx[d]++
			if tileIdx[d] < tileCounts[d] {
				break
			}
			tileIdx[d] = 0
		}
	}

	// Emit the wrapping BINTABLE HDU. The table has a single VLA column
	// COMPRESSED_DATA with the per-tile payloads.
	return emitCompressedHDU(f, hdr, shape, tileShape, bp, algo, blockSize, opts.HCompressScale, zdither0, tilePayloads, "", nil)
}

// AppendCompressedFloat32Image writes a float32 image HDU with tile
// compression. The pixel values are quantized to int32 per tile
// (using cfitsio's fits_quantize_float algorithm, ported to Go), then
// the integer tiles are fed through the chosen integer compression
// algorithm (RICE_1, GZIP_1, GZIP_2, or NoCompress).
//
// This is lossy compression: the precision of reconstructed pixels is
// approximately sigma/QLevel, where sigma is the per-tile noise
// estimate. Tiles that cannot be quantized (e.g. constant-value
// regions) automatically fall back to gzip compression of the raw
// float bytes via the GZIP_COMPRESSED_DATA column, recovering those
// pixels bit-exactly on read.
//
// The resulting HDU is a BINTABLE with ZIMAGE=T, ZQUANTIZ set to the
// dither method name, plus per-tile ZSCALE / ZZERO columns. Open()
// transparently decompresses it back to float32 on read.
func AppendCompressedFloat32Image(f *File, hdr *header.Header, shape []int64, data []float32, opts CompressFloatOptions) (*CompressedImageHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	if len(shape) == 0 {
		return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: empty shape")
	}
	var npix int64 = 1
	for _, n := range shape {
		if n <= 0 {
			return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: non-positive axis %d", n)
		}
		npix *= n
	}
	if int64(len(data)) != npix {
		return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: data length %d != product of shape %d", len(data), npix)
	}

	// Resolve options with defaults.
	algo := opts.Algorithm
	if algo == compress.Unknown {
		algo = compress.RICE1
	}
	switch algo {
	case compress.RICE1, compress.GZIP1, compress.GZIP2, compress.NoCompress:
	default:
		return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: algorithm %s does not support float quantization", algo)
	}
	tileShape := opts.TileShape
	if len(tileShape) == 0 {
		tileShape = make([]int64, len(shape))
		tileShape[0] = shape[0]
		for i := 1; i < len(shape); i++ {
			tileShape[i] = 1
		}
	}
	if len(tileShape) != len(shape) {
		return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: tile shape has %d axes, image has %d", len(tileShape), len(shape))
	}
	qlevel := opts.QLevel
	if qlevel == 0 {
		qlevel = 4.0
	}
	dither := opts.DitherMethod
	if dither == compress.NoDither && (algo == compress.RICE1 || algo == compress.GZIP1 || algo == compress.GZIP2) {
		dither = compress.SubtractiveDither1
	}
	zdither0 := opts.ZDither0
	if zdither0 == 0 {
		zdither0 = 1
	}
	blockSize := opts.BlockSize
	if blockSize <= 0 {
		blockSize = 32
	}
	nullvalue := opts.NullValue
	if nullvalue == 0 {
		nullvalue = -1e30
	}

	// Substitute NaN with sentinel if the caller said nulls are
	// possible. We mutate a copy of the input so the caller's slice
	// is untouched.
	var working []float32
	if opts.Nullcheck {
		working = make([]float32, len(data))
		copy(working, data)
		compress.SubstituteNaN32(working, nullvalue)
	} else {
		working = data
	}

	// Partition into tiles.
	tileCounts := make([]int64, len(shape))
	totalTiles := int64(1)
	for i := range shape {
		tileCounts[i] = (shape[i] + tileShape[i] - 1) / tileShape[i]
		totalTiles *= tileCounts[i]
	}

	// Pick the integer encoder. Every tile will be encoded through it.
	// BYTEPIX is always 4 for int32 quantized data.
	encoder, err := compress.SelectEncoder(algo, compress.Params{
		"BLOCKSIZE": int64(blockSize),
		"BYTEPIX":   int64(4),
	})
	if err != nil {
		return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: %w", err)
	}

	tilePayloads := make([][]byte, totalTiles)
	zscale := make([]float64, totalTiles)
	zzero := make([]float64, totalTiles)
	// gzipFallback is lazily allocated the first time a tile can't
	// be quantized. If it stays nil, no GZIP_COMPRESSED_DATA column
	// is written.
	var gzipFallback [][]byte
	gzipEncoder, _ := compress.SelectEncoder(compress.GZIP1, nil)

	tileIdx := make([]int64, len(shape))
	for t := int64(0); t < totalTiles; t++ {
		box := tileBox(shape, tileShape, tileIdx)
		var tileNpix int64 = 1
		for _, n := range box.size {
			tileNpix *= n
		}
		tileFloat := make([]float32, tileNpix)
		extractTileFloat32(working, tileFloat, shape, box)

		// Quantize this tile. tile row for PRNG is 1-based.
		idata := make([]int32, tileNpix)
		qres, ok := compress.QuantizeFloat32(t+1, tileFloat,
			box.size[0], tileNpixProduct(box.size[1:]), opts.Nullcheck, nullvalue,
			qlevel, dither, idata)
		if !ok {
			// GZIP_COMPRESSED_DATA fallback: serialize the raw
			// float32 tile in FITS big-endian order, replace any
			// null sentinels with NaN bits, gzip the result.
			if gzipFallback == nil {
				gzipFallback = make([][]byte, totalTiles)
			}
			// Make a mutable copy so we don't corrupt the user's data.
			tileCopy := make([]float32, tileNpix)
			copy(tileCopy, tileFloat)
			if opts.Nullcheck {
				compress.ReplaceSentinelWithNaNBits(tileCopy, nullvalue)
			}
			rawBE := make([]byte, tileNpix*4)
			for i := int64(0); i < tileNpix; i++ {
				bits := math.Float32bits(tileCopy[i])
				bigendian.PutUint32(rawBE[i*4:], bits)
			}
			gzBuf := make([]byte, max64(int64(len(rawBE))*2+256, 1024))
			n, err := gzipEncoder.Encode(rawBE, gzBuf, int(tileNpix), 4)
			if err == compress.ErrBufferTooSmall {
				gzBuf = make([]byte, len(rawBE)*4+1024)
				n, err = gzipEncoder.Encode(rawBE, gzBuf, int(tileNpix), 4)
			}
			if err != nil {
				return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: tile %d gzip fallback: %w", t, err)
			}
			gzipFallback[t] = gzBuf[:n]
			// COMPRESSED_DATA row stays nil (empty VLA); reader
			// picks this up and routes to the gzip fallback column.
			tilePayloads[t] = nil
			// Identity ZSCALE/ZZERO — unused for fallback tiles
			// but the columns must have one value per row.
			zscale[t] = 1
			zzero[t] = 0
			// advance tile index and continue
			for d := range shape {
				tileIdx[d]++
				if tileIdx[d] < tileCounts[d] {
					break
				}
				tileIdx[d] = 0
			}
			continue
		}
		zscale[t] = qres.BScale
		zzero[t] = qres.BZero

		// Serialize int32 to FITS big-endian bytes.
		tileRaw := make([]byte, tileNpix*4)
		for i := int64(0); i < tileNpix; i++ {
			bigendian.PutInt32(tileRaw[i*4:], idata[i])
		}

		// Encode.
		compressed := make([]byte, max64(int64(len(tileRaw))*2+256, 1024))
		n, err := encoder.Encode(tileRaw, compressed, int(tileNpix), 4)
		if err == compress.ErrBufferTooSmall {
			compressed = make([]byte, len(tileRaw)*4+1024)
			n, err = encoder.Encode(tileRaw, compressed, int(tileNpix), 4)
		}
		if err != nil {
			return nil, fmt.Errorf("fits: AppendCompressedFloat32Image: tile %d: %w", t, err)
		}
		tilePayloads[t] = compressed[:n]

		for d := range shape {
			tileIdx[d]++
			if tileIdx[d] < tileCounts[d] {
				break
			}
			tileIdx[d] = 0
		}
	}

	quantExtra := &floatQuantExtra{
		zscale:       zscale,
		zzero:        zzero,
		zblank:       compress.NullValueInt32,
		gzipFallback: gzipFallback,
	}
	quantizeName := ditherMethodName(dither)
	return emitCompressedHDU(f, hdr, shape, tileShape, bitpix.Float32, algo, blockSize, 0, zdither0, tilePayloads, quantizeName, quantExtra)
}

// ditherMethodName renders the compression dither method as the
// ZQUANTIZ keyword string.
func ditherMethodName(m compress.DitherMethod) string {
	switch m {
	case compress.SubtractiveDither1:
		return "SUBTRACTIVE_DITHER_1"
	case compress.SubtractiveDither2:
		return "SUBTRACTIVE_DITHER_2"
	}
	return "NO_DITHER"
}

// tileNpixProduct multiplies all dimensions of a tile-size tail slice
// into a single pixel count. Used to pack a multi-dimensional tile
// into the (nxpix, nypix) pair that fits_quantize_float expects.
// The quantizer treats the tile as a 2D image (nxpix wide, nypix rows
// where nypix is the product of all higher dimensions).
func tileNpixProduct(dims []int64) int64 {
	if len(dims) == 0 {
		return 1
	}
	p := int64(1)
	for _, d := range dims {
		p *= d
	}
	return p
}

// extractTileFloat32 copies a tile's region from the full float32
// image into a dense row-major tile buffer. Mirrors extractTileBytes
// but works directly on float32 slices.
func extractTileFloat32(img, tile []float32, shape []int64, region tileRegion) {
	n := len(shape)
	imgStride := make([]int64, n)
	var s int64 = 1
	for d := 0; d < n; d++ {
		imgStride[d] = s
		s *= shape[d]
	}
	rowLen := region.size[0]
	if n == 1 {
		imgOff := region.origin[0] * imgStride[0]
		copy(tile, img[imgOff:imgOff+rowLen])
		return
	}
	tileStride := make([]int64, n)
	s = 1
	for d := 0; d < n; d++ {
		tileStride[d] = s
		s *= region.size[d]
	}
	idx := make([]int64, n)
	for {
		var tileOff int64
		for d := 1; d < n; d++ {
			tileOff += idx[d] * tileStride[d]
		}
		imgOff := region.origin[0] * imgStride[0]
		for d := 1; d < n; d++ {
			imgOff += (region.origin[d] + idx[d]) * imgStride[d]
		}
		copy(tile[tileOff:tileOff+rowLen], img[imgOff:imgOff+rowLen])
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
	}
}

// serializeTileInts converts a generic integer slice to FITS big-endian
// bytes. Helper shared with the tile extractor.
func serializeTileInts[T Numeric](dst []byte, data []T, bp bitpix.BITPIX) {
	elemSize := bp.Size()
	for i, v := range data {
		switch elemSize {
		case 1:
			dst[i] = byte(anyToInt64(v))
		case 2:
			bigendian.PutInt16(dst[i*2:], int16(anyToInt64(v)))
		case 4:
			bigendian.PutInt32(dst[i*4:], int32(anyToInt64(v)))
		case 8:
			bigendian.PutInt64(dst[i*8:], anyToInt64(v))
		}
	}
}

// anyToInt64 is a small helper that casts a generic Numeric value to
// int64 via unsafe-free type-switching.
func anyToInt64[T Numeric](v T) int64 {
	switch x := any(v).(type) {
	case int8:
		return int64(x)
	case uint8:
		return int64(x)
	case int16:
		return int64(x)
	case uint16:
		return int64(x)
	case int32:
		return int64(x)
	case uint32:
		return int64(x)
	case int64:
		return x
	case uint64:
		return int64(x)
	}
	return 0
}

// extractTileBytes copies a tile's region from the full-image byte
// buffer into a dense row-major tile buffer. Mirrors copyTileIntoImage
// from the read side but in the opposite direction.
func extractTileBytes(imgBuf, tileBuf []byte, shape []int64, region tileRegion, elemSize int) {
	n := len(shape)
	imgStride := make([]int64, n)
	var s int64 = int64(elemSize)
	for d := 0; d < n; d++ {
		imgStride[d] = s
		s *= shape[d]
	}
	rowBytes := region.size[0] * int64(elemSize)
	if n == 1 {
		imgOff := region.origin[0] * imgStride[0]
		copy(tileBuf, imgBuf[imgOff:imgOff+rowBytes])
		return
	}
	tileStride := make([]int64, n)
	s = int64(elemSize)
	for d := 0; d < n; d++ {
		tileStride[d] = s
		s *= region.size[d]
	}
	idx := make([]int64, n)
	for {
		var tileOff int64
		for d := 1; d < n; d++ {
			tileOff += idx[d] * tileStride[d]
		}
		imgOff := region.origin[0] * imgStride[0]
		for d := 1; d < n; d++ {
			imgOff += (region.origin[d] + idx[d]) * imgStride[d]
		}
		copy(tileBuf[tileOff:tileOff+rowBytes], imgBuf[imgOff:imgOff+rowBytes])
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
	}
}

// floatQuantExtra carries the per-tile float-quantization columns that
// accompany COMPRESSED_DATA when a float image is written. nil means
// "integer image, no quantization columns". All three slices, if set,
// must have length == totalTiles.
type floatQuantExtra struct {
	zscale []float64 // per-tile bscale, written as ZSCALE column
	zzero  []float64 // per-tile bzero,  written as ZZERO  column
	// gzipFallback carries raw-float gzip payloads for tiles that
	// could not be quantized (constant or out-of-range). If nil, no
	// GZIP_COMPRESSED_DATA column is emitted. Tiles with a nil inner
	// slice used the normal COMPRESSED_DATA column.
	gzipFallback [][]byte
	// zblank, if non-zero, is emitted as the ZBLANK header keyword —
	// the integer sentinel stored where input pixels were null.
	zblank int32
}

// emitCompressedHDU writes the binary-table wrapper containing the
// compressed tile payloads. Uses AppendBinaryTable under the hood with
// a single VLA column COMPRESSED_DATA, plus optional ZSCALE/ZZERO
// per-tile columns and a GZIP_COMPRESSED_DATA fallback column for
// float tile compression.
func emitCompressedHDU(f *File, hdr *header.Header, shape, tileShape []int64, bp bitpix.BITPIX,
	algo compress.Algorithm, blockSize, hscale, zdither0 int,
	payloads [][]byte, quantize string, quantExtra *floatQuantExtra) (*CompressedImageHDU, error) {

	// Build Z* keywords. They go into hdr (merged with any user keywords).
	zhdr := header.New()
	zhdr.Set("ZIMAGE", true, "extension contains compressed image")
	zhdr.Set("ZTENSION", "IMAGE", "Image extension")
	zhdr.Set("ZBITPIX", int64(bp), "data type of original image")
	zhdr.Set("ZNAXIS", int64(len(shape)), "dimension of original image")
	for i, n := range shape {
		zhdr.Set("ZNAXIS"+strconv.Itoa(i+1), n, "length of original image axis")
	}
	zhdr.Set("ZPCOUNT", int64(0), "")
	zhdr.Set("ZGCOUNT", int64(1), "")
	for i, n := range tileShape {
		zhdr.Set("ZTILE"+strconv.Itoa(i+1), n, "size of tiles to be compressed")
	}
	zhdr.Set("ZCMPTYPE", algo.String(), "compression algorithm")
	// Algorithm-specific parameters.
	zidx := 1
	switch algo {
	case compress.RICE1:
		zhdr.Set("ZNAME"+strconv.Itoa(zidx), "BLOCKSIZE", "compression block size")
		zhdr.Set("ZVAL"+strconv.Itoa(zidx), int64(blockSize), "pixels per block")
		zidx++
		zhdr.Set("ZNAME"+strconv.Itoa(zidx), "BYTEPIX", "bytes per pixel")
		zhdr.Set("ZVAL"+strconv.Itoa(zidx), int64(bp.Size()), "bytes per pixel")
	case compress.HCOMPRESS1:
		zhdr.Set("ZNAME"+strconv.Itoa(zidx), "SCALE", "HCOMPRESS scale factor")
		zhdr.Set("ZVAL"+strconv.Itoa(zidx), int64(hscale), "HCOMPRESS scale factor")
		zidx++
		zhdr.Set("ZNAME"+strconv.Itoa(zidx), "SMOOTH", "HCOMPRESS smooth option")
		zhdr.Set("ZVAL"+strconv.Itoa(zidx), int64(0), "HCOMPRESS smooth option")
	}
	if quantize != "" {
		zhdr.Set("ZQUANTIZ", quantize, "quantization method")
		zhdr.Set("ZDITHER0", int64(zdither0), "dither offset")
	}
	if quantExtra != nil && quantExtra.zblank != 0 {
		zhdr.Set("ZBLANK", int64(quantExtra.zblank), "null pixel value")
	}
	// Merge caller-supplied keywords (ignore any that would clash with
	// the structural set).
	if hdr != nil {
		for _, c := range hdr.Cards() {
			if !isCompressedZKey(c.Key) {
				zhdr.Add(c.Key, c.Value, c.Comment)
			}
		}
	}

	// PLIO_1 uses a 1PI (int16 VLA) column because its natural unit is
	// the 16-bit token, not a raw byte. Every other algorithm stores
	// its opaque byte stream in a 1PB (uint8 VLA) column.
	var cols []ColumnData
	if algo == compress.PLIO1 {
		tokenRows := make([][]int16, len(payloads))
		for i, p := range payloads {
			if len(p)%2 != 0 {
				return nil, fmt.Errorf("fits: PLIO_1 payload length %d not even", len(p))
			}
			toks := make([]int16, len(p)/2)
			for k := 0; k < len(toks); k++ {
				toks[k] = int16(uint16(p[k*2])<<8 | uint16(p[k*2+1]))
			}
			tokenRows[i] = toks
		}
		cols = []ColumnData{{
			Name:         "COMPRESSED_DATA",
			DataVarInt16: tokenRows,
		}}
	} else {
		cols = []ColumnData{{
			Name:         "COMPRESSED_DATA",
			DataVarUint8: payloads,
		}}
	}
	// Float-quantization extras: per-tile ZSCALE/ZZERO and an
	// optional GZIP_COMPRESSED_DATA fallback column.
	if quantExtra != nil {
		if len(quantExtra.zscale) != len(payloads) || len(quantExtra.zzero) != len(payloads) {
			return nil, fmt.Errorf("fits: emitCompressedHDU: zscale/zzero length %d/%d != tiles %d",
				len(quantExtra.zscale), len(quantExtra.zzero), len(payloads))
		}
		cols = append(cols,
			ColumnData{Name: "ZSCALE", DataFloat64: quantExtra.zscale},
			ColumnData{Name: "ZZERO", DataFloat64: quantExtra.zzero},
		)
		if quantExtra.gzipFallback != nil {
			if len(quantExtra.gzipFallback) != len(payloads) {
				return nil, fmt.Errorf("fits: emitCompressedHDU: gzipFallback length %d != tiles %d",
					len(quantExtra.gzipFallback), len(payloads))
			}
			cols = append(cols, ColumnData{
				Name:         "GZIP_COMPRESSED_DATA",
				DataVarUint8: quantExtra.gzipFallback,
			})
		}
	}
	tbl, err := AppendBinaryTable(f, zhdr, cols)
	if err != nil {
		return nil, err
	}
	return &CompressedImageHDU{rec: tbl.rec, tbl: tbl}, nil
}

// isCompressedZKey reports whether a keyword is one of the Z-prefixed
// structural keys we synthesize and therefore shouldn't be duplicated
// from the caller's header.
func isCompressedZKey(key string) bool {
	switch key {
	case "ZIMAGE", "ZTENSION", "ZBITPIX", "ZNAXIS", "ZPCOUNT", "ZGCOUNT",
		"ZCMPTYPE", "ZQUANTIZ", "ZDITHER0", "ZBLANK":
		return true
	}
	if len(key) > 5 && (key[:5] == "ZTILE" || key[:5] == "ZNAXI") {
		return true
	}
	if len(key) > 4 && (key[:4] == "ZVAL" || key[:4] == "ZBSC" || key[:4] == "ZBZE") {
		return true
	}
	if len(key) > 5 && key[:5] == "ZNAME" {
		return true
	}
	return false
}

func max64(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}
