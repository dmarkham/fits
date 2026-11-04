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

// The writers in this file produce the same compressed images as
// cfitsio's own writer (fits_img_compress, the path fpack uses): for
// every case in compress/testdata/cref/golden/img (and img64, BITPIX -64)
// the RICE_1,
// HCOMPRESS_1 and PLIO_1 tiles are byte-identical to cfitsio's, GZIP
// tiles inflate to identical bytes, and ZSCALE/ZZERO and the compression
// keywords are equal (TestCompressedWriterMatchesCfitsio). Every
// algorithm x pixel type x dither x tiling the writers accept is also
// checked end to end with fitsverify, the cfitsio library, funpack and
// astropy (TestCompressedWriterInterop). Inputs cfitsio refuses are
// refused here too, with an error. Where this writer deliberately
// differs from cfitsio it says so on the option concerned: HCOMPRESS_1
// tiles with a thin edge are refused instead of enlarged, PLIO_1 uint16
// and 64-bit integers are refused, a caller BLANK on float input and a
// caller SIMPLE that is not logical are refused, no default EXTNAME is
// added, caller compression keywords and duplicate keywords are dropped,
// and a real pixel equal to FLOATNULLVALUE stays a value.

// CompressFloatOptions configures a tile-compressed float32 or float64
// image write. The zero value is usable: RICE_1, whole-row tiles, quantize
// level 4, SUBTRACTIVE_DITHER_1 with ZDITHER0 = 1 -- the same defaults as
// cfitsio/fpack (except that the dither seed is fixed, not clock-based).
type CompressFloatOptions struct {
	// Algorithm is the compression method: RICE_1, GZIP_1 or GZIP_2.
	// Default: RICE_1. HCOMPRESS_1 and PLIO_1 are refused for float
	// input by this writer, and NOCOMPRESS is never written (cfitsio and
	// funpack cannot read it).
	Algorithm compress.Algorithm

	// TileShape picks per-tile dimensions, same convention as
	// CompressOptions.TileShape. Default: whole-row tiles.
	TileShape []int64

	// BlockSize is the RICE_1 block size. Default 32 (what cfitsio
	// writes).
	BlockSize int

	// QLevel is the quantization level. Positive values set the step
	// as sigma/QLevel (sigma = cfitsio's per-tile noise estimate);
	// negative values use an absolute step of -QLevel. Default 4, or 16
	// when NoDither is set (cfitsio's defaults). NaN and -Inf are refused
	// with an error (no quantization step: cfitsio accepts them and every
	// pixel then reads back NaN); +Inf gives a zero step, so every tile is
	// stored losslessly (GZIP_COMPRESSED_DATA), as in cfitsio.
	QLevel float32

	// DitherMethod selects SUBTRACTIVE_DITHER_1 or SUBTRACTIVE_DITHER_2
	// (exact zeros preserved). compress.NoDither is the zero value of
	// this field and therefore means "default" (SUBTRACTIVE_DITHER_1);
	// use NoDither below to write without dithering.
	DitherMethod compress.DitherMethod

	// NoDither writes ZQUANTIZ = 'NO_DITHER': plain quantization with no
	// random offsets (cfitsio's quantize method NO_DITHER). It cannot be
	// combined with a DitherMethod.
	NoDither bool

	// ZDither0 is the dither seed written as the ZDITHER0 keyword, 1 to
	// 10000. Default 1. Tile i (1-based) is dithered from row
	// i + ZDither0 - 1 of cfitsio's random sequence, as in cfitsio.
	// Ignored (and no ZDITHER0 written) with NoDither.
	ZDither0 int

	// Pixels are taken the way cfitsio's writer reads them (fnan in
	// fffr4r4): NaN, +Inf and -Inf are always stored as nulls -- the
	// ZBLANK integer in quantized tiles, NaN in unquantizable
	// (GZIP_COMPRESSED_DATA) tiles -- and every reader gives NaN back;
	// +-0 and denormals (all-zero exponent) are stored as +0, also in
	// GZIP_COMPRESSED_DATA tiles. When Nullcheck is true, pixels equal
	// to NullValue are nulls too (and come back as NaN), compared with the
	// raw value as fits_write_imgnull compares them (a denormal NullValue
	// nulls the pixels equal to it). NullValue 0 means the default, -1e30,
	// so 0 itself cannot be the null value. NullValue is for the float32
	// writer; the float64 writer refuses a non-zero NullValue (it takes
	// NullValue64) and the float32 writer a non-zero NullValue64.
	//
	// Whatever null value the quantizer is handed internally, a real pixel
	// is never turned into a null: with Nullcheck off a real -1e30 (the
	// default NullValue) in a tile with NaNs stays a value (golden
	// rice_f32_nan_realnull_t40x5).
	//
	// One deliberate difference from fits_img_compress: in a tile that has
	// nulls, it also treats a real pixel equal to its internal sentinel
	// FLOATNULLVALUE (-9.11912e-36) as a null (golden
	// rice_f32_floatnullvalue_t40x5); this writer keeps it as a value. In
	// a GZIP_COMPRESSED_DATA (lossless) tile cfitsio's fits_read_img and
	// astropy read it back exactly; in a quantized tile it is quantized
	// like any value, and the tiles are the ones cfitsio writes when the
	// image is written with fits_write_imgnull and null value -1e30
	// (golden imgw/, TestCompressedWriterFloatNullValueQuantizedMatchesCfitsio).
	// funpack writes its output with FLOATNULLVALUE as the null value, so
	// that one value comes back NaN from funpack
	// (TestFloatNullValuePixelStaysValue).
	NullValue float32
	Nullcheck bool

	// NullValue64 is the null value of AppendCompressedFloat64Image (with
	// Nullcheck; 0 means -1e30). Proven against fits_write_imgnull (golden
	// imgw/rice_f64_nullval_m999_t40x5).
	NullValue64 float64
}

// CompressOptions configures a tile-compressed integer image write.
// The zero value is usable: RICE_1, cfitsio's default tiling, lossless.
type CompressOptions struct {
	// Algorithm is the compression method: RICE_1 (default), GZIP_1,
	// GZIP_2, HCOMPRESS_1 or PLIO_1. NOCOMPRESS is refused: cfitsio and
	// funpack cannot read it. PLIO_1 takes pixels in 0..2^24-1 and is
	// refused for uint16 (cfitsio and astropy disagree on that storage)
	// and, in practice, uint32 (stored as v-2^31, as in cfitsio).
	Algorithm compress.Algorithm

	// TileShape is the per-tile dimensions in the same order as the
	// image's logical shape (first axis varies fastest). A dimension
	// longer than its image axis is written as given (ZTILEn > NAXISn;
	// the tiles along that axis are then the whole axis), as cfitsio
	// does. Nil or empty
	// picks cfitsio's default: whole-row tiles (ZTILE1 = NAXIS1, the
	// rest 1), and for HCOMPRESS_1 NAXIS1 x (NAXIS2 if it is at most 30,
	// else 16 rows or the first of 24, 20, 30, 28, 26, 22, 18, 14 that
	// leaves no last tile of fewer than 4 rows, else 17).
	//
	// HCOMPRESS_1 needs an image of two or more axes and 2-D tiles over
	// the first two axes (ZTILE3.. = 1), each at least 4 pixels, with no
	// edge tile thinner than 4 pixels; other shapes are refused (cfitsio
	// would silently enlarge the tile instead). The default rule is
	// checked against cfitsio for every NAXIS2 from 4 to 2000 and the
	// 26/22/18/17-row fallbacks (TestHcompressDefaultTileMatchesCfitsio).
	TileShape []int64

	// BlockSize is the RICE_1 block size. Default 32 (what cfitsio
	// writes).
	BlockSize int

	// HCompressScale is the absolute HCOMPRESS_1 digitization scale
	// (cfitsio: fits_set_hcomp_scale with a negative value). 0 is
	// lossless. Written as ZVAL1 = -scale like cfitsio. Negative values
	// are refused (cfitsio's noise-relative scale is not supported).
	// Lossy decoding can overshoot the pixel type's range near its
	// limits (e.g. 256 for uint8); funpack then fails with "datatype
	// conversion overflow", exactly as for cfitsio-written files.
	HCompressScale int

	// ZDither0 is unused for integer images (kept for compatibility).
	ZDither0 int
}

// intStorage describes how an integer Go type is stored in a compressed
// image, exactly as cfitsio stores it (imcomp_init_table and
// imcomp_convert_tile_t* in imcompress.c).
type intStorage struct {
	bp       bitpix.BITPIX // ZBITPIX
	bzero    float64       // BZERO written to the header (0 = none)
	unsigned bool          // the type is stored with an offset
}

// AppendCompressedImage writes a new tile-compressed image HDU holding
// integer pixels (uint8, int8, int16, uint16, int32, uint32). The
// compression is lossless except for HCOMPRESS_1 with HCompressScale > 1.
// Float32 and float64 images use AppendCompressedFloat32Image and
// AppendCompressedFloat64Image (passing them here is an error). int64/uint64
// are refused: cfitsio cannot write compressed BITPIX 64 images (fpack:
// "Bad image datatype").
//
// Unsigned and signed-byte types are stored the way cfitsio stores them:
// uint16 as ZBITPIX 16 with BZERO = 32768, uint32 as ZBITPIX 32 with
// BZERO = 2147483648, int8 as ZBITPIX 8 with BZERO = -128 (BSCALE = 1).
// A caller header may not carry any other BSCALE/BZERO: cfitsio refuses
// scaled writes to compressed images, and this writer stores the given
// values unscaled.
//
// hdr may be nil or a header copied from an uncompressed image. Its
// structural keywords are replaced; SIMPLE, EXTEND, BLOCKED, CHECKSUM and
// DATASUM are stored as ZSIMPLE, ZEXTEND, ZBLOCKED, ZHECKSUM and
// ZDATASUM and every other card is copied, as cfitsio's
// imcomp_copy_img2comp does (golden img/rice_i16_hdr). Unlike cfitsio no
// default EXTNAME = 'COMPRESSED_IMAGE' is added, and (also unlike
// cfitsio, which copies them) compression keywords already in hdr
// (ZBITPIX, ZCMPTYPE, ZTILEn, ...) are dropped and a keyword is never
// written twice: the first card wins, so a header with both EXTEND and
// ZEXTEND yields one ZEXTEND (fitsverify-clean, TestCompressedWriterInterop
// i16_rice_caller_zkeys).
//
// Caller cards are copied with header.AddCard. A card read from a file is
// written as its bytes. An edited Cards() copy (Key/Type/Value/Comment
// changed) is written from its fields under the upper-cased keyword, and
// is left out like any caller card for a keyword the writer computes
// itself. An unedited lowercase-keyword (non-HIERARCH) card makes the
// writer return an error wrapping header.ErrInvalidKeyword: cfitsio and
// astropy read such a card differently.
//
// The resulting HDU is a binary table with ZIMAGE = T that Open()
// recognizes as a *CompressedImageHDU.
func AppendCompressedImage[T Numeric](f *File, hdr *header.Header, shape []int64, data []T, opts CompressOptions) (*CompressedImageHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	npix, err := checkCompressedShape("AppendCompressedImage", shape, len(data))
	if err != nil {
		return nil, err
	}
	var st intStorage
	var zero T
	switch any(zero).(type) {
	case uint8:
		st = intStorage{bp: bitpix.Int8}
	case int8:
		st = intStorage{bp: bitpix.Int8, bzero: -128, unsigned: true}
	case int16:
		st = intStorage{bp: bitpix.Int16}
	case uint16:
		st = intStorage{bp: bitpix.Int16, bzero: 32768, unsigned: true}
	case int32:
		st = intStorage{bp: bitpix.Int32}
	case uint32:
		st = intStorage{bp: bitpix.Int32, bzero: 2147483648, unsigned: true}
	case int64, uint64:
		return nil, fmt.Errorf("fits: AppendCompressedImage: %T pixels: cfitsio cannot write or read tile-compressed 64-bit integer images", zero)
	case float64:
		return nil, fmt.Errorf("fits: AppendCompressedImage: float64 pixels: use AppendCompressedFloat64Image")
	default:
		return nil, fmt.Errorf("fits: AppendCompressedImage: integer type required, got %T (use AppendCompressedFloat32Image for float32)", zero)
	}
	if err := checkCallerScaling("AppendCompressedImage", hdr, st.bzero, false); err != nil {
		return nil, err
	}

	algo := opts.Algorithm
	if algo == compress.Unknown {
		algo = compress.RICE1
	}
	switch algo {
	case compress.RICE1, compress.GZIP1, compress.GZIP2, compress.HCOMPRESS1, compress.PLIO1:
	case compress.NoCompress:
		return nil, errNoCompressWrite("AppendCompressedImage")
	default:
		return nil, fmt.Errorf("fits: AppendCompressedImage: unsupported algorithm %s", algo)
	}
	if opts.HCompressScale < 0 {
		return nil, fmt.Errorf("fits: AppendCompressedImage: HCompressScale %d < 0 (only absolute scales are supported)", opts.HCompressScale)
	}
	tileShape, err := resolveTileShape("AppendCompressedImage", algo, shape, opts.TileShape)
	if err != nil {
		return nil, err
	}
	blockSize := opts.BlockSize
	if blockSize <= 0 {
		blockSize = 32
	}

	// Serialize into the stored (offset) integers, big-endian.
	elemSize := st.bp.Size()
	encElem := elemSize
	if _, isU16 := any(zero).(uint16); isU16 && algo == compress.PLIO1 {
		// cfitsio stores PLIO_1 uint16 tiles as the raw value and its
		// reader subtracts 32768 from BZERO for PLIO only; astropy 8.0.2
		// applies BZERO as written and reads every pixel 32768 too high,
		// even in cfitsio's own files. No layout reads the same in both.
		return nil, fmt.Errorf("fits: AppendCompressedImage: PLIO_1 with uint16 pixels is not written: cfitsio and astropy decode that storage differently (use int16, int32 or uint8 masks)")
	}
	imgBytes := make([]byte, npix*int64(encElem))
	for i, v := range data {
		x := anyToInt64(v)
		switch {
		case elemSize == 1:
			if st.unsigned { // int8: flip the sign bit (v + 128)
				imgBytes[i] = byte(x) ^ 0x80
			} else {
				imgBytes[i] = byte(x)
			}
		case elemSize == 2:
			if st.unsigned { // uint16: v - 32768
				bigendian.PutInt16(imgBytes[i*2:], int16(uint16(x)^0x8000))
			} else {
				bigendian.PutInt16(imgBytes[i*2:], int16(x))
			}
		default:
			if st.unsigned { // uint32: v - 2^31
				bigendian.PutInt32(imgBytes[i*4:], int32(uint32(x)^0x80000000))
			} else {
				bigendian.PutInt32(imgBytes[i*4:], int32(x))
			}
		}
	}

	tileCounts, totalTiles := countTiles(shape, tileShape)
	encoder, err := compress.SelectEncoder(algo, compress.Params{
		"BLOCKSIZE": int64(blockSize),
		"BYTEPIX":   int64(encElem),
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
		tileRaw := make([]byte, tileNpix*int64(encElem))
		extractTileBytes(imgBytes, tileRaw, shape, box, encElem)
		if setter, ok := encoder.(interface{ SetShape(nx, ny int) }); ok {
			// cfitsio's HCOMPRESS has ny = fastest axis, nx = slow.
			setter.SetShape(int(box.size[1]), int(box.size[0]))
		}
		payload, err := encodeTile(encoder, tileRaw, int(tileNpix), encElem)
		if err != nil {
			return nil, fmt.Errorf("fits: AppendCompressedImage: tile %d: %w", t+1, err)
		}
		tilePayloads[t] = payload
		nextTile(tileIdx, tileCounts)
	}

	return emitCompressedHDU(f, hdr, shape, tileShape, compressedHDUParams{
		bp: st.bp, algo: algo, blockSize: blockSize, bytePix: elemSize,
		hscale: opts.HCompressScale, bzero: st.bzero,
	}, tilePayloads, nil)
}

// encodeTile runs one tile through the encoder, growing the output
// buffer until it fits.
func encodeTile(enc compress.Encoder, raw []byte, nelem, elemSize int) ([]byte, error) {
	size := max64(int64(len(raw))*2+256, 1024)
	for {
		buf := make([]byte, size)
		n, err := enc.Encode(raw, buf, nelem, elemSize)
		if err == nil {
			return buf[:n], nil
		}
		if err != compress.ErrBufferTooSmall || size > int64(len(raw))*64+1<<20 {
			return nil, err
		}
		size *= 4
	}
}

// maxCompressDim is cfitsio's MAX_COMPRESS_DIM: fits_img_compress refuses
// an image of more axes (status 212, golden img/rice_i16_7d), and funpack
// and fitsverify reject such a file.
const maxCompressDim = 6

func checkCompressedShape(fn string, shape []int64, ndata int) (int64, error) {
	if len(shape) == 0 {
		return 0, fmt.Errorf("fits: %s: empty shape", fn)
	}
	if len(shape) > maxCompressDim {
		return 0, fmt.Errorf("fits: %s: %d axes; tile compression takes at most %d (cfitsio MAX_COMPRESS_DIM)", fn, len(shape), maxCompressDim)
	}
	var npix int64 = 1
	for _, n := range shape {
		if n <= 0 {
			return 0, fmt.Errorf("fits: %s: non-positive axis %d", fn, n)
		}
		npix *= n
	}
	if int64(ndata) != npix {
		return 0, fmt.Errorf("fits: %s: %w: data length %d != product of shape %d", fn, ErrShapeMismatch, ndata, npix)
	}
	return npix, nil
}

// checkCallerScaling refuses caller BSCALE/BZERO other than the storage
// convention of the pixel type (BZERO = bzero, BSCALE = 1), a BLANK for a
// floating-point image, and a SIMPLE that is not a logical value.
//
// BLANK: the FITS standard (4.0, sect. 4.4.2.5) allows it only for integer
// BITPIX. cfitsio's writer copies a caller BLANK into the ZBITPIX -32
// header, and with no ZBLANK (a tile without NaNs) cfitsio, funpack and
// astropy then all take BLANK as the null value of the quantized integers:
// real pixels whose integer equals it decode as NaN (golden
// rice_f32_blank_hdr, TestFloatBLANKNullsRealPixels). Refused, not dropped:
// the caller's header claims a null value the data cannot carry.
//
// SIMPLE: written as ZSIMPLE, which must be T or F; a value the writer
// would have to replace is refused rather than silently rewritten.
func checkCallerScaling(fn string, hdr *header.Header, bzero float64, isFloat bool) error {
	if hdr == nil {
		return nil
	}
	if isFloat && hdr.Has("BLANK") {
		return fmt.Errorf("fits: %s: BLANK is not allowed in a floating-point image header (FITS 4.0 sect. 4.4.2.5): readers would take it as the null value of the quantized integers and null real pixels; remove it (NaN marks nulls in float images)", fn)
	}
	if hdr.Has("SIMPLE") {
		if _, err := hdr.Bool("SIMPLE"); err != nil {
			return fmt.Errorf("fits: %s: SIMPLE in the header is not a logical value (%v); it is written as ZSIMPLE, which must be T or F", fn, err)
		}
	}
	if hdr.Has("BSCALE") {
		v, err := hdr.Float("BSCALE")
		if err != nil || v != 1 {
			return fmt.Errorf("fits: %s: BSCALE in the header is not supported for compressed writes (cfitsio refuses scaled compressed writes)", fn)
		}
	}
	if hdr.Has("BZERO") {
		v, err := hdr.Float("BZERO")
		if err != nil || v != bzero {
			return fmt.Errorf("fits: %s: BZERO in the header must be %g for this pixel type (got %v); scaled compressed writes are not supported", fn, bzero, v)
		}
	}
	return nil
}

func errNoCompressWrite(fn string) error {
	return fmt.Errorf("fits: %s: NOCOMPRESS is not written: cfitsio and funpack cannot read ZCMPTYPE = 'NOCOMPRESS' (\"unknown compression algorithm\"); use GZIP_1 or RICE_1", fn)
}

// resolveTileShape applies cfitsio's default tiling and validates the
// tile shape for the algorithm.
func resolveTileShape(fn string, algo compress.Algorithm, shape, req []int64) ([]int64, error) {
	n := len(shape)
	if algo == compress.HCOMPRESS1 && n < 2 {
		return nil, fmt.Errorf("fits: %s: HCOMPRESS_1 cannot be used with 1-dimensional images", fn)
	}
	tile := make([]int64, n)
	if len(req) == 0 {
		tile[0] = shape[0]
		for i := 1; i < n; i++ {
			tile[i] = 1
		}
		if algo == compress.HCOMPRESS1 {
			tile[1] = hcompressDefaultRows(shape[1])
		}
	} else {
		if len(req) != n {
			return nil, fmt.Errorf("fits: %s: tile shape has %d axes, image has %d", fn, len(req), n)
		}
		for i, v := range req {
			if v <= 0 {
				return nil, fmt.Errorf("fits: %s: non-positive tile axis %d", fn, v)
			}
			// A tile longer than the axis is kept as requested (ZTILEn >
			// NAXISn), as cfitsio writes it (golden rice_i16_11x8_t100x3).
			tile[i] = v
		}
	}
	if algo == compress.HCOMPRESS1 {
		for i := 2; i < n; i++ {
			if tile[i] != 1 {
				return nil, fmt.Errorf("fits: %s: HCOMPRESS_1 tiles must be 2-D over the first two axes (ZTILE%d = %d)", fn, i+1, tile[i])
			}
		}
		for i := 0; i < 2; i++ {
			if tile[i] < 4 {
				return nil, fmt.Errorf("fits: %s: HCOMPRESS_1 minimum tile dimension is 4 pixels (ZTILE%d = %d)", fn, i+1, tile[i])
			}
			if r := shape[i] % tile[i]; r > 0 && r < 4 {
				return nil, fmt.Errorf("fits: %s: HCOMPRESS_1 last tile along axis %d would have %d < 4 pixels; choose another tile size", fn, i+1, r)
			}
		}
	}
	return tile, nil
}

// hcompressDefaultRows is cfitsio's default HCOMPRESS tile height
// (imcomp_init_table): the whole image up to 30 rows, else 16 rows or
// the first alternative that leaves the last tile at least 4 rows.
func hcompressDefaultRows(naxis2 int64) int64 {
	if naxis2 <= 30 {
		return naxis2
	}
	for _, h := range []int64{16, 24, 20, 30, 28, 26, 22, 18, 14} {
		if naxis2%h == 0 || naxis2%h > 3 {
			return h
		}
	}
	return 17
}

func countTiles(shape, tileShape []int64) ([]int64, int64) {
	tileCounts := make([]int64, len(shape))
	total := int64(1)
	for i := range shape {
		tileCounts[i] = (shape[i] + tileShape[i] - 1) / tileShape[i]
		total *= tileCounts[i]
	}
	return tileCounts, total
}

func nextTile(idx, counts []int64) {
	for d := range idx {
		idx[d]++
		if idx[d] < counts[d] {
			return
		}
		idx[d] = 0
	}
}

// AppendCompressedFloat32Image writes a float32 image HDU with tile
// compression. Each tile is quantized to int32 with a pure-Go port of
// cfitsio's fits_quantize_float (bit-exact, compress/testdata/cref), then
// compressed with RICE_1, GZIP_1 or GZIP_2, exactly as cfitsio does.
//
// This is lossy: the reconstructed pixels differ from the input by up to
// about half a quantization step (sigma/QLevel). Tiles that cannot be
// quantized (constant tiles, too few pixels, a range too wide for int32)
// are stored losslessly as gzipped raw floats in a GZIP_COMPRESSED_DATA
// column, as cfitsio does. So are tiles on which cfitsio's quantizer
// arithmetic overflows (a zero point beyond int64 steps, a pixel beyond
// int32 or on the null sentinel; see compress.QuantizeFloat32), which
// cfitsio writes as wrong pixels (TestQuantizerOverflowStoredLosslessly). NaN and +-Inf pixels are stored as nulls and
// read back as NaN; +-0 and denormals are stored as +0 (see
// CompressFloatOptions).
//
// The resulting HDU is a BINTABLE with ZIMAGE = T, ZQUANTIZ naming the
// dither method, and per-tile ZSCALE / ZZERO columns.
func AppendCompressedFloat32Image(f *File, hdr *header.Header, shape []int64, data []float32, opts CompressFloatOptions) (*CompressedImageHDU, error) {
	return appendCompressedFloat(f, "AppendCompressedFloat32Image", hdr, shape, data, opts)
}

// AppendCompressedFloat64Image is AppendCompressedFloat32Image for float64
// pixels: a ZBITPIX -64 image, each tile quantized to int32 with a port of
// cfitsio's fits_quantize_double, exactly as fits_img_compress (fpack)
// writes a BITPIX -64 image (golden compress/testdata/cref/golden/img64,
// TestCompressedWriterFloat64MatchesCfitsio). Unquantizable tiles are
// stored losslessly as gzipped raw float64 in GZIP_COMPRESSED_DATA. Pixels
// are classified as for float32 (NaN and +-Inf are nulls; +-0 and
// denormals are stored as +0); with Nullcheck the null value is
// NullValue64 (0 means -1e30), not NullValue. A real pixel equal to
// cfitsio's DOUBLENULLVALUE (-9.1191291391491e-36) in a tile with nulls
// stays a value (fits_img_compress nulls it, as it does FLOATNULLVALUE).
func AppendCompressedFloat64Image(f *File, hdr *header.Header, shape []int64, data []float64, opts CompressFloatOptions) (*CompressedImageHDU, error) {
	return appendCompressedFloat(f, "AppendCompressedFloat64Image", hdr, shape, data, opts)
}

// floatPixel is the pixel type of the quantizing writer.
type floatPixel interface{ float32 | float64 }

// appendCompressedFloat is the float32 / float64 writer: the C code
// (imcomp_convert_tile_tfloat / _tdouble, fits_quantize_float / _double,
// imcomp_float2nan / imcomp_double2nan) is the same with the pixel width
// changed.
func appendCompressedFloat[T floatPixel](f *File, fn string, hdr *header.Header, shape []int64, data []T, opts CompressFloatOptions) (*CompressedImageHDU, error) {
	var zero T
	_, is64 := any(zero).(float64)
	width := int64(4)
	if is64 {
		width = 8
	}
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	if _, err := checkCompressedShape(fn, shape, len(data)); err != nil {
		return nil, err
	}
	if err := checkCallerScaling(fn, hdr, 0, true); err != nil {
		return nil, err
	}
	algo := opts.Algorithm
	if algo == compress.Unknown {
		algo = compress.RICE1
	}
	switch algo {
	case compress.RICE1, compress.GZIP1, compress.GZIP2:
	case compress.NoCompress:
		return nil, errNoCompressWrite(fn)
	default:
		return nil, fmt.Errorf("fits: %s: algorithm %s is not supported for float input (use RICE_1, GZIP_1 or GZIP_2)", fn, algo)
	}
	tileShape, err := resolveTileShape(fn, algo, shape, opts.TileShape)
	if err != nil {
		return nil, err
	}
	dither := opts.DitherMethod
	switch {
	case opts.NoDither && dither != compress.NoDither:
		return nil, fmt.Errorf("fits: %s: NoDither and DitherMethod %d are mutually exclusive", fn, dither)
	case opts.NoDither:
	case dither == compress.NoDither:
		dither = compress.SubtractiveDither1
	case dither != compress.SubtractiveDither1 && dither != compress.SubtractiveDither2:
		return nil, fmt.Errorf("fits: %s: unknown DitherMethod %d", fn, dither)
	}
	qlevel := opts.QLevel
	if qlevel != qlevel || math.IsInf(float64(qlevel), -1) {
		// NaN: delta NaN; -Inf: delta +Inf and zero point 0*Inf = NaN.
		// cfitsio does not check (its range tests are false for NaN) and
		// writes tiles that read back all NaN; refuse instead.
		return nil, fmt.Errorf("fits: %s: QLevel %v has no quantization step (NaN or -Inf)", fn, qlevel)
	}
	if qlevel == 0 {
		if opts.NoDither {
			qlevel = 16
		} else {
			qlevel = 4
		}
	}
	zdither0 := opts.ZDither0
	if zdither0 == 0 {
		zdither0 = 1
	}
	if zdither0 < 1 || zdither0 > compress.NRandom {
		return nil, fmt.Errorf("fits: %s: ZDither0 %d out of range 1..%d", fn, zdither0, compress.NRandom)
	}
	blockSize := opts.BlockSize
	if blockSize <= 0 {
		blockSize = 32
	}
	var nullvalue T
	if is64 {
		if opts.NullValue != 0 {
			return nil, fmt.Errorf("fits: %s: NullValue %v is the float32 null value; the float64 writer takes NullValue64", fn, opts.NullValue)
		}
		nullvalue = T(opts.NullValue64)
	} else {
		if opts.NullValue64 != 0 {
			return nil, fmt.Errorf("fits: %s: NullValue64 %v is the float64 null value; the float32 writer takes NullValue", fn, opts.NullValue64)
		}
		nullvalue = T(opts.NullValue)
	}
	if nullvalue == 0 {
		nullvalue = T(-1e30)
	}

	tileCounts, totalTiles := countTiles(shape, tileShape)
	encoder, err := compress.SelectEncoder(algo, compress.Params{
		"BLOCKSIZE": int64(blockSize),
		"BYTEPIX":   int64(4),
	})
	if err != nil {
		return nil, fmt.Errorf("fits: %s: %w", fn, err)
	}

	tilePayloads := make([][]byte, totalTiles)
	zscale := make([]float64, totalTiles)
	zzero := make([]float64, totalTiles)
	// gzipFallback is allocated the first time a tile can't be
	// quantized; if it stays nil, no GZIP_COMPRESSED_DATA column is
	// written (as cfitsio adds the column on first use).
	var gzipFallback [][]byte
	gzipEncoder, _ := compress.SelectEncoder(compress.GZIP1, nil)
	anyNulls := false

	tileIdx := make([]int64, len(shape))
	for t := int64(0); t < totalTiles; t++ {
		box := tileBox(shape, tileShape, tileIdx)
		var tileNpix int64 = 1
		for _, n := range box.size {
			tileNpix *= n
		}
		tileFloat := make([]T, tileNpix)
		extractTile(data, tileFloat, shape, box)

		// Classify pixels the way cfitsio reads a float tile before
		// compressing it (fits_read_subset with a null value, i.e.
		// fffr4r4's fnan test): an all-ones exponent (NaN, +Inf, -Inf)
		// is a null; an all-zero exponent (+-0 and denormals) is read as
		// +0. Proven by goldens rice_f32_inf_* and
		// gzip2_f32_inf_fallback_t40x5. On request, pixels equal to
		// NullValue are nulls, compared with the raw value as
		// fits_write_imgnull does (so a denormal NullValue nulls its
		// pixels; golden imgw/rice_f32_nullval_denorm_t40x5).
		isNull := make([]bool, tileNpix)
		nulls := 0
		for i, v := range tileFloat {
			if opts.Nullcheck && v == nullvalue {
				isNull[i] = true
				nulls++
				continue
			}
			switch floatClass(v) {
			case classNull:
				isNull[i] = true
				nulls++
			case classZero:
				tileFloat[i] = 0
			}
		}
		// cfitsio quantizes with nullcheck only for tiles that contain
		// nulls, after replacing them with a sentinel. Its value cannot
		// reach the output (the quantizer skips null pixels), so pick
		// one that no real pixel of this tile equals.
		quantIn := tileFloat
		sentinel := nullvalue
		if nulls > 0 {
			anyNulls = true
			sentinel = absentSentinel(tileFloat, isNull, nullvalue)
			quantIn = make([]T, tileNpix)
			for i, v := range tileFloat {
				if isNull[i] {
					quantIn[i] = sentinel
				} else {
					quantIn[i] = v
				}
			}
		}

		// cfitsio imcomp_convert_tile_tfloat: irow = row + dither_seed - 1
		// (row is the 1-based tile number); 0 means "do not dither".
		var irow int64
		if !opts.NoDither {
			irow = t + 1 + int64(zdither0) - 1
		}
		idata := make([]int32, tileNpix)
		var qres compress.QuantizeResult
		var ok bool
		switch q := any(quantIn).(type) {
		case []float32:
			qres, ok = compress.QuantizeFloat32(irow, q,
				box.size[0], tileNpixProduct(box.size[1:]), nulls > 0, any(sentinel).(float32),
				qlevel, dither, idata)
		case []float64:
			qres, ok = compress.QuantizeFloat64(irow, q,
				box.size[0], tileNpixProduct(box.size[1:]), nulls > 0, any(sentinel).(float64),
				qlevel, dither, idata)
		}
		if !ok {
			// GZIP_COMPRESSED_DATA fallback: the raw float32/float64 tile
			// in FITS big-endian order with nulls as the all-ones NaN
			// pattern (cfitsio imcomp_float2nan / imcomp_double2nan),
			// gzipped. ZSCALE/ZZERO stay 0 for this row, as cfitsio
			// leaves them.
			if gzipFallback == nil {
				gzipFallback = make([][]byte, totalTiles)
			}
			rawBE := make([]byte, tileNpix*width)
			for i := int64(0); i < tileNpix; i++ {
				putFloatBE(rawBE[i*width:], tileFloat[i], isNull[i])
			}
			gz, err := encodeTile(gzipEncoder, rawBE, int(tileNpix), int(width))
			if err != nil {
				return nil, fmt.Errorf("fits: %s: tile %d gzip fallback: %w", fn, t+1, err)
			}
			gzipFallback[t] = gz
			tilePayloads[t] = nil // empty COMPRESSED_DATA cell
			nextTile(tileIdx, tileCounts)
			continue
		}
		zscale[t] = qres.BScale
		zzero[t] = qres.BZero

		tileRaw := make([]byte, tileNpix*4)
		for i := int64(0); i < tileNpix; i++ {
			bigendian.PutInt32(tileRaw[i*4:], idata[i])
		}
		payload, err := encodeTile(encoder, tileRaw, int(tileNpix), 4)
		if err != nil {
			return nil, fmt.Errorf("fits: %s: tile %d: %w", fn, t+1, err)
		}
		tilePayloads[t] = payload
		nextTile(tileIdx, tileCounts)
	}

	quantExtra := &floatQuantExtra{
		zscale:       zscale,
		zzero:        zzero,
		gzipFallback: gzipFallback,
	}
	if anyNulls {
		quantExtra.zblank = compress.NullValueInt32
	}
	params := compressedHDUParams{
		bp: floatBitpix(is64), algo: algo, blockSize: blockSize, bytePix: 4,
		quantize: ditherMethodName(dither, opts.NoDither),
	}
	if !opts.NoDither {
		params.zdither0 = zdither0
	}
	return emitCompressedHDU(f, hdr, shape, tileShape, params, tilePayloads, quantExtra)
}

// absentSentinel returns a value equal to no non-null pixel of tile,
// starting from want and stepping down.
func absentSentinel[T floatPixel](tile []T, isNull []bool, want T) T {
	s := want
	for {
		clash := false
		for i, v := range tile {
			if !isNull[i] && v == s {
				clash = true
				break
			}
		}
		if !clash {
			return s
		}
		switch p := any(&s).(type) {
		case *float32:
			*p = math.Nextafter32(*p, float32(math.Inf(-1)))
		case *float64:
			*p = math.Nextafter(*p, math.Inf(-1))
		}
	}
}

// Pixel classes of cfitsio's fnan / dnan tests (fitsio2.h), applied when
// fits_img_compress reads a tile with a null value: an all-ones exponent
// (NaN, +-Inf) is a null, an all-zero exponent (+-0, denormals) is read
// as 0.
const (
	classValue = iota
	classNull
	classZero
)

func floatClass[T floatPixel](v T) int {
	var exp, mask uint64
	switch x := any(v).(type) {
	case float32:
		mask = 0x7F800000
		exp = uint64(math.Float32bits(x)) & mask
	case float64:
		mask = 0x7FF0000000000000
		exp = math.Float64bits(x) & mask
	}
	switch exp {
	case mask:
		return classNull
	case 0:
		return classZero
	}
	return classValue
}

// putFloatBE writes v big-endian (4 or 8 bytes), or the all-ones NaN
// pattern for a null.
func putFloatBE[T floatPixel](b []byte, v T, null bool) {
	switch x := any(v).(type) {
	case float32:
		bits := math.Float32bits(x)
		if null {
			bits = 0xFFFFFFFF
		}
		bigendian.PutUint32(b, bits)
	case float64:
		bits := math.Float64bits(x)
		if null {
			bits = 0xFFFFFFFFFFFFFFFF
		}
		bigendian.PutUint64(b, bits)
	}
}

func floatBitpix(is64 bool) bitpix.BITPIX {
	if is64 {
		return bitpix.Float64
	}
	return bitpix.Float32
}

// ditherMethodName renders the compression dither method as the
// ZQUANTIZ keyword string.
func ditherMethodName(m compress.DitherMethod, noDither bool) string {
	if noDither {
		return "NO_DITHER"
	}
	switch m {
	case compress.SubtractiveDither2:
		return "SUBTRACTIVE_DITHER_2"
	}
	return "SUBTRACTIVE_DITHER_1"
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

// extractTile copies a tile's region from the full float32
// image into a dense row-major tile buffer. Mirrors extractTileBytes
// but works directly on float32 slices.
func extractTile[T any](img, tile []T, shape []int64, region tileRegion) {
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
// "integer image, no quantization columns". All slices, if set, must
// have length == totalTiles.
type floatQuantExtra struct {
	zscale []float64 // per-tile bscale, written as ZSCALE column
	zzero  []float64 // per-tile bzero,  written as ZZERO  column
	// gzipFallback carries raw-float gzip payloads for tiles that
	// could not be quantized. If nil, no GZIP_COMPRESSED_DATA column is
	// emitted. Tiles with a nil inner slice used COMPRESSED_DATA.
	gzipFallback [][]byte
	// zblank, if non-zero, is emitted as the ZBLANK header keyword —
	// the integer stored where input pixels were null.
	zblank int32
}

// compressedHDUParams are the header-level parameters of a compressed
// image HDU.
type compressedHDUParams struct {
	bp        bitpix.BITPIX // ZBITPIX
	algo      compress.Algorithm
	blockSize int     // RICE_1 ZVAL1
	bytePix   int     // RICE_1 ZVAL2
	hscale    int     // HCOMPRESS_1 absolute scale
	bzero     float64 // BZERO for unsigned storage (0 = none)
	quantize  string  // ZQUANTIZ ("" = integer image)
	zdither0  int     // ZDITHER0 (0 = not written)
}

// emitCompressedHDU writes the binary-table wrapper containing the
// compressed tile payloads: a COMPRESSED_DATA VLA column (1PI for
// PLIO_1, 1PB otherwise), plus ZSCALE/ZZERO and an optional
// GZIP_COMPRESSED_DATA column for float tile compression.
func emitCompressedHDU(f *File, hdr *header.Header, shape, tileShape []int64, p compressedHDUParams,
	payloads [][]byte, quantExtra *floatQuantExtra) (*CompressedImageHDU, error) {

	// Build Z* keywords. They go into hdr (merged with any user keywords).
	zb := newHeaderBuilder()
	zb.set("ZIMAGE", true, "extension contains compressed image")
	simple, primary := callerSimple(hdr)
	if primary {
		// A primary-HDU header: cfitsio records SIMPLE as ZSIMPLE and
		// writes no ZTENSION (imcomp_copy_img2comp).
		zb.set("ZSIMPLE", simple, "file does conform to FITS standard")
	} else {
		zb.set("ZTENSION", "IMAGE", "Image extension")
	}
	zb.set("ZBITPIX", int64(p.bp), "data type of original image")
	zb.set("ZNAXIS", int64(len(shape)), "dimension of original image")
	for i, n := range shape {
		zb.set("ZNAXIS"+strconv.Itoa(i+1), n, "length of original image axis")
	}
	if !primary { // an extension has PCOUNT/GCOUNT; a primary HDU does not
		zb.set("ZPCOUNT", int64(0), "")
		zb.set("ZGCOUNT", int64(1), "")
	}
	for i, n := range tileShape {
		zb.set("ZTILE"+strconv.Itoa(i+1), n, "size of tiles to be compressed")
	}
	zb.set("ZCMPTYPE", p.algo.String(), "compression algorithm")
	// Algorithm-specific parameters, as cfitsio writes them.
	switch p.algo {
	case compress.RICE1:
		zb.set("ZNAME1", "BLOCKSIZE", "compression block size")
		zb.set("ZVAL1", int64(p.blockSize), "pixels per block")
		zb.set("ZNAME2", "BYTEPIX", "bytes per pixel (1, 2, 4, or 8)")
		zb.set("ZVAL2", int64(p.bytePix), "bytes per pixel (1, 2, 4, or 8)")
	case compress.HCOMPRESS1:
		// cfitsio records the requested scale; a negative value means
		// "absolute", which is what HCompressScale is.
		zb.set("ZNAME1", "SCALE", "HCOMPRESS scale factor")
		// -0.0 would be written as "-0."; cfitsio writes 0.
		zval1 := 0.0
		if p.hscale != 0 {
			zval1 = -float64(p.hscale)
		}
		zb.set("ZVAL1", zval1, "HCOMPRESS scale factor")
		zb.set("ZNAME2", "SMOOTH", "HCOMPRESS smooth option")
		zb.set("ZVAL2", int64(0), "HCOMPRESS smooth option")
	}
	if p.quantize != "" {
		zb.set("ZQUANTIZ", p.quantize, "quantization method")
		if p.zdither0 != 0 {
			zb.set("ZDITHER0", int64(p.zdither0), "dithering offset when quantizing floats")
		}
	}
	if quantExtra != nil && quantExtra.zblank != 0 {
		zb.set("ZBLANK", int64(quantExtra.zblank), "null value in the compressed integer array")
	}
	if p.bzero != 0 {
		// The cards cfitsio writes (integer literals, its comments); the
		// header encoder refuses a comment that does not fit the card.
		comment := map[float64]string{
			-128:       "offset data range to that of signed byte",
			32768:      "offset data range to that of unsigned short",
			2147483648: "offset data range to that of unsigned long",
		}[p.bzero]
		zb.set("BZERO", int64(p.bzero), comment)
		zb.set("BSCALE", int64(1), "default scaling factor")
	}
	// Merge caller-supplied keywords (ignore any that would clash with
	// the structural set; BSCALE/BZERO were validated by the caller).
	// Keywords that describe the uncompressed HDU and are illegal or
	// wrong in a BINTABLE are renamed the way cfitsio's
	// imcomp_copy_img2comp renames them: EXTEND -> ZEXTEND, BLOCKED ->
	// ZBLOCKED, CHECKSUM -> ZHECKSUM, DATASUM -> ZDATASUM (SIMPLE ->
	// ZSIMPLE above); XTENSION is replaced by ZTENSION. The structural
	// image keywords (BITPIX, NAXISn, PCOUNT, GCOUNT) are dropped by
	// AppendBinaryTable. Proven by golden img/rice_i16_hdr.
	if hdr != nil {
		for _, c := range hdr.Cards() {
			key := c.Key
			switch key {
			case "SIMPLE", "XTENSION", "BSCALE", "BZERO":
				continue
			case "EXTEND":
				key = "ZEXTEND"
			case "BLOCKED":
				key = "ZBLOCKED"
			case "CHECKSUM":
				key = "ZHECKSUM"
			case "DATASUM":
				key = "ZDATASUM"
			}
			if isCompressedZKey(key) {
				continue
			}
			// No duplicate keywords (fitsverify warns on them): the first
			// card for a key wins, e.g. a caller header with both EXTEND
			// and ZEXTEND keeps one ZEXTEND. Commentary cards repeat.
			if key != "COMMENT" && key != "HISTORY" && key != "" && zb.h.Has(key) {
				continue
			}
			// Through the builder (core-header rule): a card read from a
			// file keeps its bytes; a renamed one gets only its key field
			// rewritten, as cfitsio overwrites the name in the card; a card
			// the encoder cannot represent fails the call.
			if key != c.Key {
				c.Key = key
				if c.Raw != ([header.CardWidth]byte{}) {
					copy(c.Raw[:8], fmt.Sprintf("%-8s", key))
				}
			}
			zb.add(c)
		}
	}

	// PLIO_1 uses a 1PI (int16 VLA) column because its natural unit is
	// the 16-bit token, not a raw byte. Every other algorithm stores
	// its opaque byte stream in a 1PB (uint8 VLA) column.
	var cols []ColumnData
	if p.algo == compress.PLIO1 {
		tokenRows := make([][]int16, len(payloads))
		for i, b := range payloads {
			if len(b)%2 != 0 {
				return nil, fmt.Errorf("fits: PLIO_1 payload length %d not even", len(b))
			}
			toks := make([]int16, len(b)/2)
			for k := 0; k < len(toks); k++ {
				toks[k] = int16(uint16(b[k*2])<<8 | uint16(b[k*2+1]))
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
	if zb.err != nil {
		return nil, fmt.Errorf("fits: compressed image header: %w", zb.err)
	}
	tbl, err := AppendBinaryTable(f, zb.h, cols)
	if err != nil {
		return nil, err
	}
	return &CompressedImageHDU{rec: tbl.rec, tbl: tbl}, nil
}

// callerSimple returns the caller header's SIMPLE value, if present
// (checkCallerScaling has refused a SIMPLE that is not logical).
func callerSimple(hdr *header.Header) (bool, bool) {
	if hdr == nil || !hdr.Has("SIMPLE") {
		return false, false
	}
	v, _ := hdr.Bool("SIMPLE")
	return v, true
}

// isCompressedZKey reports whether a keyword is one of the Z-prefixed
// structural keys we synthesize and therefore shouldn't be duplicated
// from the caller's header.
func isCompressedZKey(key string) bool {
	switch key {
	case "ZIMAGE", "ZTENSION", "ZSIMPLE", "ZBITPIX", "ZNAXIS", "ZPCOUNT", "ZGCOUNT",
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
