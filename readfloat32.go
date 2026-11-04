package fits

import (
	"context"
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// Siril's float conversion constants (Siril src/core/siril.h:76-77):
// INV_USHRT_MAX_SINGLE and INV_UCHAR_MAX_SINGLE. They are decimal
// literals, not 1/65535 and 1/255: the 8-bit one is 1 ULP below
// float32(1/255), so Siril maps 255 to 0.99999994, not 1.
var (
	sirilInvU16 float32 = .000015259022
	sirilInvU8  float32 = .0039215686
)

// ReadFloat32 reads every image pixel as float32 for processing pipelines
// (stacking, calibration, debayering, ...). It accepts *ImageHDU and
// *CompressedImageHDU. Pixels are returned in FITS order (NAXIS1 fastest).
//
// Integer BITPIX is mapped onto [0, 1]; the mapping increases with the
// physical value (raw·BSCALE + BZERO), so the physical minimum of the
// BITPIX range goes to 0 and the maximum to (about) 1:
//
//   - BITPIX=8 with BSCALE=1, BZERO=0, and BITPIX=16 with BSCALE=1 and
//     BZERO=0 or 32768 (the camera formats): out = float32(u) · c in
//     float32 arithmetic, where u is the unsigned sample (raw for BITPIX=8,
//     raw+32768 for BITPIX=16) and c is Siril's constant .0039215686f or
//     .000015259022f. This is bit-identical to Siril's loader
//     (ushort_to_float_range / uchar_to_float_range, and readfits with
//     force_float), proven over every 8- and 16-bit value against Siril's
//     own code compiled from source (testdata/imaging/sirilref). Note that
//     Siril's 8-bit constant maps 255 to 0.99999994.
//   - Every other integer case (BITPIX=32 and 64, BZERO=2^31 / 2^63
//     unsigned, signed bytes BZERO=-128, any other BSCALE/BZERO): out is the
//     float32 nearest to u / (rawmax − rawmin), where u = raw − rawmin for
//     BSCALE > 0 and rawmax − raw for BSCALE < 0 (the map is independent of
//     the size of BSCALE and of BZERO). This is NOT what Siril does: Siril
//     returns 32-bit integers unnormalized unless DATAMAX > 10, rejects
//     BITPIX=64, and multiplies other scaled data by 1/255 or 1/65535.
//   - Pixels whose stored value equals BLANK become NaN (Siril does not
//     mask BLANK; it normalizes the sentinel like any other value). For a
//     tile-compressed image the null value is the ZBLANK keyword, else the
//     BLANK keyword, of the compressed HDU, and a per-tile ZBLANK column
//     takes precedence for its tile, as cfitsio reads it (proven against
//     cfitsio's null flags, testdata/imaging/cfnull).
//   - BSCALE = 0 is an error. BSCALE / BZERO / BLANK are read by the
//     scaling-card rule of every image read (scalingCards; ImageHDU.BSCALE):
//     the last card of each keyword is applied (the first on a
//     tile-compressed image), as cfitsio and astropy apply it, and it must
//     be a valid value card ("= " in columns 9-10, upper-case keyword, a
//     FITS number cfitsio converts alike; BLANK an integer in the BITPIX
//     range), with no HIERARCH BSCALE / BZERO / BLANK, else an error
//     wrapping ErrInvalidScaling: e.g. a string, including a quoted number
//     such as '10.0'; "2.0.0"; a C hex float such as 0x1p3; 1E400, for
//     which cfitsio's fits_read_img silently uses the default (strtod's
//     value for 0x1p3). D exponents (1.0D1) are valid and read. The result
//     is the same on a fresh Open and after Header().
//
// For a tile-compressed image BSCALE and BZERO are the keywords of the
// compressed HDU header, as cfitsio reads them (the non-standard
// ZBSCALE/ZBZERO keywords are ignored, as cfitsio ignores them).
//
// Float BITPIX (-32, -64) is returned as the physical value
// float32(raw·BSCALE + BZERO) (evaluated in float64), with no
// normalization and no clamping: values outside [0, 1], ±Inf and NaN pass
// through. A finite physical value too large for float32 is an error
// rather than a silent +Inf. Siril additionally divides float data by
// 65535 when its maximum exceeds 10; ReadFloat32 does not.
//
// Differences from Siril for files Siril reads differently are pinned
// one by one in TestReadFloat32_SirilConvention, e.g. a BITPIX=16 file with
// no BZERO card: ReadFloat32 treats it as signed (+32768), Siril as
// unsigned without offset.
func ReadFloat32(h HDU) ([]float32, error) {
	var (
		raw           []byte
		bp            bitpix.BITPIX
		bscale, bzero float64
		blank         int64
		haveBlank     bool
		null          []bool // per-pixel null flags (tile-compressed images)
		err           error
	)
	switch img := h.(type) {
	case *ImageHDU:
		bp = img.bitpixType()
		cards, cerr := img.rec.imageCards()
		if cerr != nil {
			return nil, fmt.Errorf("fits.ReadFloat32: %w", cerr)
		}
		bscale, bzero, blank, haveBlank = cards.bscale, cards.bzero, cards.blank, cards.haveBlank
		n := img.pixelCount()
		if n == 0 {
			return []float32{}, nil
		}
		raw, err = readRawBytes(img, n)
	case *CompressedImageHDU:
		meta, merr := img.pixelMetadata()
		if merr != nil {
			return nil, merr
		}
		bp = meta.Bitpix
		// Scaling: BSCALE and BZERO of the compressed HDU header by the
		// scaling-card rule (metadata: scalingCards with the first card
		// winning, as cfitsio's imcomp_get_compressed_image_par and astropy
		// read it; an explicit BSCALE = 0 is kept so it is rejected below).
		// Null pixels are the tile decoder's null flags, cfitsio's rule:
		// the tile's ZBLANK column entry, else the ZBLANK keyword, else the
		// BLANK keyword.
		bscale, bzero = meta.BScale, meta.BZero
		if bp == bitpix.Int64 {
			// Not supported by ReadFloat32 (no test proves this path;
			// ReadPixelsCompressed reads astropy's ZBITPIX 64 GZIP and
			// NOCOMPRESS images, TestCompressedInt64Read).
			return nil, fmt.Errorf("fits.ReadFloat32: ZBITPIX 64 tile-compressed images are not supported")
		}
		// decompressAllTiles returns integer tiles unscaled (nulls are
		// matched against blank below) and float tiles as ReadPixels
		// takes them: quantized tiles already in physical units (the
		// image BSCALE/BZERO folded into each tile's ZSCALE/ZZERO, as
		// cfitsio imcomp_decompress_tile), unquantized tiles scaled by the
		// whole-image ZSCALE/ZZERO keywords in place of BSCALE/BZERO.
		raw, null, err = decompressAllTiles(context.Background(), img, meta)
		switch {
		case bp.IsFloat() && meta.ZScaleCols:
			bscale, bzero = 1, 0
		case bp.IsFloat() && meta.HasFloatScale:
			bscale, bzero = meta.FloatScale, meta.FloatZero
		}
	default:
		return nil, fmt.Errorf("fits.ReadFloat32: unsupported HDU type %T (need *ImageHDU or *CompressedImageHDU)", h)
	}
	if err != nil {
		return nil, err
	}
	// bscale and bzero are finite here: scalingCards refuses a card that
	// is not a finite number, and the header parser refuses NaN, Inf and
	// out-of-range spellings of ZSCALE / ZZERO (hdr.Float then fails and
	// gives 0, which is refused here as BSCALE = 0).
	if bscale == 0 {
		return nil, fmt.Errorf("fits.ReadFloat32: invalid scaling BSCALE=%g BZERO=%g", bscale, bzero)
	}
	if bp.IsFloat() {
		return float32Physical(raw, bp, bscale, bzero)
	}
	return normalizeInts(raw, bp, bscale, bzero, blank, haveBlank, null)
}

// strictFITSInteger parses a FITS integer value (FITS 4.0 §4.2.3:
// optional sign and decimal digits) that fits in an int64. ParseInt with
// base 10 accepts exactly that syntax (no underscores, no prefix).
func strictFITSInteger(tok string) (int64, bool) {
	n, err := strconv.ParseInt(tok, 10, 64)
	return n, err == nil
}

// strictFITSNumber parses a FITS integer or real value (FITS 4.0 §4.2.3,
// §4.2.4: optional sign, digits with an optional decimal point, optional
// E or D exponent) that is finite as a float64. A lower-case e exponent
// is accepted too: cfitsio 4.7.0 and astropy 8.0.2 both read 1.0e1 as 10.
// A lower-case d is not: cfitsio drops such a card (reads 0) while
// astropy reads 10.
func strictFITSNumber(tok string) (float64, bool) {
	t := tok
	if t != "" && (t[0] == '+' || t[0] == '-') {
		t = t[1:]
	}
	// Mantissa: digits with at most one decimal point and at least one
	// digit. Exponent: an optionally signed digit string. Both are checked
	// here because ParseFloat is looser (hex floats, NaN / Inf spellings,
	// underscores between digits: "1E1_0" is 1e10 to ParseFloat). An empty
	// exponent ("1E", "1E+") is left to ParseFloat, which refuses it, and
	// so is a value too large for a float64 (1E400: ErrRange).
	mant := t
	if i := strings.IndexAny(t, "EeD"); i >= 0 {
		mant = t[:i]
		if strings.Trim(strings.TrimLeft(t[i+1:], "+-"), "0123456789") != "" {
			return 0, false
		}
	}
	intPart, frac, _ := strings.Cut(mant, ".")
	if strings.Trim(intPart, "0123456789") != "" || strings.Trim(frac, "0123456789") != "" ||
		(intPart == "" && frac == "") {
		return 0, false
	}
	v, err := strconv.ParseFloat(strings.Replace(tok, "D", "E", 1), 64)
	if err != nil {
		return 0, false
	}
	return v, true
}

// float32Physical decodes float BITPIX to float32 physical values.
func float32Physical(raw []byte, bp bitpix.BITPIX, bscale, bzero float64) ([]float32, error) {
	scaled := bscale != 1 || bzero != 0
	size := bp.Size()
	out := make([]float32, len(raw)/size)
	for i := range out {
		var v float64
		if bp == bitpix.Float32 {
			f := bigendian.Float32(raw[i*4:])
			if !scaled {
				out[i] = f
				continue
			}
			v = float64(f)
		} else {
			v = bigendian.Float64(raw[i*8:])
		}
		if scaled {
			v = float64(v*bscale) + bzero
		}
		f := float32(v)
		if math.IsInf(float64(f), 0) && !math.IsInf(v, 0) && !math.IsNaN(v) {
			return nil, fmt.Errorf("fits.ReadFloat32: pixel %d: physical value %g overflows float32", i, v)
		}
		out[i] = f
	}
	return out, nil
}

// normalizeInts maps integer BITPIX samples onto [0, 1] (see ReadFloat32).
func normalizeInts(raw []byte, bp bitpix.BITPIX, bscale, bzero float64, blank int64, haveBlank bool, null []bool) ([]float32, error) {
	size := bp.Size()
	out := make([]float32, len(raw)/size)
	nan := float32(math.NaN())
	sirilRule := bscale == 1 && ((bp == bitpix.Int8 && bzero == 0) ||
		(bp == bitpix.Int16 && (bzero == 0 || bzero == 32768)))
	invert := bscale < 0

	// rawAt returns the stored integer of pixel i.
	rawAt := func(i int) int64 {
		switch bp {
		case bitpix.Int8:
			return int64(raw[i])
		case bitpix.Int16:
			return int64(bigendian.Int16(raw[i*2:]))
		case bitpix.Int32:
			return int64(bigendian.Int32(raw[i*4:]))
		default:
			return bigendian.Int64(raw[i*8:])
		}
	}
	var rawMin int64
	var span float64 // rawmax - rawmin, exact for BITPIX <= 32
	switch bp {
	case bitpix.Int8:
		rawMin, span = 0, 255
	case bitpix.Int16:
		rawMin, span = math.MinInt16, 65535
	case bitpix.Int32:
		rawMin, span = math.MinInt32, 4294967295
	case bitpix.Int64:
		rawMin = math.MinInt64
	default:
		return nil, fmt.Errorf("fits.ReadFloat32: unsupported BITPIX %d", int(bp))
	}

	var spanRat, uRat *big.Rat
	if bp == bitpix.Int64 {
		spanRat = new(big.Rat).SetInt(new(big.Int).SetUint64(math.MaxUint64))
		uRat = new(big.Rat)
	}
	for i := range out {
		r := rawAt(i)
		if (haveBlank && r == blank) || (null != nil && null[i]) {
			out[i] = nan
			continue
		}
		switch {
		case sirilRule && bp == bitpix.Int8:
			out[i] = float32(float32(r) * sirilInvU8)
		case sirilRule:
			out[i] = float32(float32(r+32768) * sirilInvU16)
		case bp == bitpix.Int64:
			// u = raw - rawmin as an exact uint64 (two's complement).
			u := uint64(r) ^ (1 << 63)
			if invert {
				u = math.MaxUint64 - u
			}
			uRat.SetFrac(new(big.Int).SetUint64(u), big.NewInt(1))
			q := new(big.Rat).Quo(uRat, spanRat)
			out[i], _ = q.Float32()
		default:
			u := float64(r - rawMin) // exact: < 2^32
			if invert {
				u = span - u
			}
			out[i] = divRoundFloat32(u, span)
		}
	}
	return out, nil
}

// divRoundFloat32 returns the float32 nearest to a/b for exactly
// representable doubles a >= 0, b > 0. The double quotient is turned into
// a round-to-odd result with its exact remainder (FMA), which then rounds
// to float32 correctly — no double-rounding error. For BITPIX 32 the naive
// float32(a/b) is wrong for exactly two values, a = 4294966655 and
// 4294967167 with b = 2^32-1 (an exhaustive search, find_i32_round.py);
// testdata/imaging/i32_round*.fits hold them.
func divRoundFloat32(a, b float64) float32 {
	q := a / b
	if r := math.FMA(-q, b, a); r != 0 && math.Float64bits(q)&1 == 0 {
		if r > 0 {
			q = math.Nextafter(q, math.Inf(1))
		} else {
			q = math.Nextafter(q, math.Inf(-1))
		}
	}
	return float32(q)
}
