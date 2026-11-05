package fits

import (
	"fmt"
	"math"
	"math/big"
	"math/bits"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
)

// This file holds the pixel value model shared by every image read and
// write path:
//
//	physical = BZERO + BSCALE * stored
//
// Floating-point targets follow cfitsio's arithmetic exactly (one double
// expression, then a single conversion to the target), so float32/float64
// reads equal cfitsio's TFLOAT/TDOUBLE reads bit for bit, except that a
// value beyond the float32 range reads as +-Inf into float32 where cfitsio
// returns +-FLT_MAX with NUM_OVERFLOW. Integer targets
// are computed exactly in 128-bit integer arithmetic whenever BSCALE and
// BZERO are integers, so BITPIX 64 data and the uint64 BZERO = 2^63
// convention never pass through float64.

// i128 is a two's-complement 128-bit integer (hi:lo).
type i128 struct {
	hi int64
	lo uint64
}

func i128From64(v int64) i128   { return i128{hi: v >> 63, lo: uint64(v)} }
func i128FromU64(v uint64) i128 { return i128{lo: v} }

func (a i128) add(b i128) i128 {
	lo, c := bits.Add64(a.lo, b.lo, 0)
	return i128{hi: a.hi + b.hi + int64(c), lo: lo}
}

func (a i128) neg() i128 {
	lo, b := bits.Sub64(0, a.lo, 0)
	return i128{hi: -a.hi - int64(b), lo: lo}
}

func (a i128) sub(b i128) i128 { return a.add(b.neg()) }

// cmp returns -1, 0, +1.
func (a i128) cmp(b i128) int {
	switch {
	case a.hi < b.hi:
		return -1
	case a.hi > b.hi:
		return 1
	case a.lo < b.lo:
		return -1
	case a.lo > b.lo:
		return 1
	}
	return 0
}

// mul64 is the exact product of two int64 values.
func mul64(a, b int64) i128 {
	neg := (a < 0) != (b < 0)
	ua, ub := uint64(a), uint64(b)
	if a < 0 {
		ua = -ua
	}
	if b < 0 {
		ub = -ub
	}
	hi, lo := bits.Mul64(ua, ub)
	r := i128{hi: int64(hi), lo: lo}
	if neg {
		r = r.neg()
	}
	return r
}

// i128FromFloat converts an integral float64 with |f| < 2^126 exactly.
func i128FromFloat(f float64) (i128, bool) {
	if f != math.Trunc(f) || math.Abs(f) >= 0x1p126 {
		return i128{}, false
	}
	bf := new(big.Float).SetFloat64(f)
	bi, _ := bf.Int(nil)
	neg := bi.Sign() < 0
	bi.Abs(bi)
	lo := new(big.Int).And(bi, new(big.Int).SetUint64(math.MaxUint64)).Uint64()
	hi := new(big.Int).Rsh(bi, 64).Uint64()
	r := i128{hi: int64(hi), lo: lo}
	if neg {
		r = r.neg()
	}
	return r, true
}

// storedRange is the raw on-disk integer range of an integer BITPIX
// (BITPIX 8 is unsigned 0..255).
func storedRange(bp bitpix.BITPIX) (lo, hi int64) {
	switch bp {
	case bitpix.Int8:
		return 0, math.MaxUint8
	case bitpix.Int16:
		return math.MinInt16, math.MaxInt16
	case bitpix.Int32:
		return math.MinInt32, math.MaxInt32
	case bitpix.Int64:
		return math.MinInt64, math.MaxInt64
	}
	return 0, 0
}

// intRange returns T's value range when T is an integer type.
func intRange[T Numeric]() (lo, hi i128, isInt, unsigned bool) {
	var z T
	switch any(z).(type) {
	case int8:
		return i128From64(math.MinInt8), i128From64(math.MaxInt8), true, false
	case uint8:
		return i128{}, i128From64(math.MaxUint8), true, true
	case int16:
		return i128From64(math.MinInt16), i128From64(math.MaxInt16), true, false
	case uint16:
		return i128{}, i128From64(math.MaxUint16), true, true
	case int32:
		return i128From64(math.MinInt32), i128From64(math.MaxInt32), true, false
	case uint32:
		return i128{}, i128From64(math.MaxUint32), true, true
	case int64:
		return i128From64(math.MinInt64), i128From64(math.MaxInt64), true, false
	case uint64:
		return i128{}, i128FromU64(math.MaxUint64), true, true
	}
	// Named types (~int16 etc.) are classified by behaviour.
	z--
	if z > 0 { // unsigned integer kinds wrap
		return i128{}, i128FromU64(uint64(z)), true, true
	}
	return i128{}, i128{}, false, false
}

// scaling is an image's BSCALE/BZERO as cfitsio uses them (doubles), plus
// their exact integer form when both are integral.
type scaling struct {
	bscale, bzero float64
	identity      bool // BSCALE == 1 && BZERO == 0
	exact         bool // BSCALE and BZERO integral (and |BSCALE| < 2^63)
	mul           int64
	off           i128
}

func newScaling(bscale, bzero float64) scaling {
	s := scaling{bscale: bscale, bzero: bzero, identity: bscale == 1 && bzero == 0}
	if bscale == math.Trunc(bscale) && math.Abs(bscale) < 0x1p63 {
		if off, ok := i128FromFloat(bzero); ok {
			s.exact, s.mul, s.off = true, int64(bscale), off
		}
	}
	return s
}

// physInt is the exact physical value of a stored integer (s.exact only).
func (s scaling) physInt(raw int64) i128 {
	if s.mul == 1 {
		return i128From64(raw).add(s.off)
	}
	return mul64(raw, s.mul).add(s.off)
}

// isU64Convention reports the BITPIX 64 unsigned convention (BZERO 2^63).
func (s scaling) isU64Convention(bp bitpix.BITPIX) bool {
	return bp == bitpix.Int64 && s.bscale == 1 && s.bzero == 0x1p63
}

// physRange is the exact physical range over every stored value
// (s.exact only).
func (s scaling) physRange(bp bitpix.BITPIX) (lo, hi i128) {
	rlo, rhi := storedRange(bp)
	a, b := s.physInt(rlo), s.physInt(rhi)
	if a.cmp(b) > 0 {
		a, b = b, a
	}
	return a, b
}

// canHold reports whether integer target T holds every physical value of
// an integer-BITPIX image exactly. Float targets always "hold" (with the
// rounding of a float conversion, as in cfitsio).
func canHold[T Numeric](bp bitpix.BITPIX, s scaling) bool {
	tlo, thi, isInt, _ := intRange[T]()
	if !isInt {
		return true
	}
	if bp.IsFloat() || !s.exact {
		return false
	}
	lo, hi := s.physRange(bp)
	return lo.cmp(tlo) >= 0 && hi.cmp(thi) <= 0
}

// --- keyword access ---

// imageCards holds the BSCALE, BZERO and BLANK values of an image by the
// scaling-card rule (scalingCards).
type imageCards struct {
	bscale, bzero float64
	blank         int64
	haveBlank     bool
	count         [3]int // cards per keyword (cardKeys order)
}

// cardKeys are the keywords scalingCards reads, in its slot order.
var cardKeys = [3]string{header.KeyBscale, header.KeyBzero, header.KeyBlankV}

// scalingCards is the one BSCALE / BZERO / BLANK rule of this package. Every
// image read and write uses it: ReadPixels, ReadPixelsMasked, CanConvert,
// ImageHDU.BSCALE / BZERO, ReadFloat32, ReadFrame, AsImage, image.Decode /
// DecodeConfig, OverwritePixels, WriteImage (the caller's cards) and
// tile-compressed images (the compressed HDU's cards). cards are the
// header's 80-byte card images: the stored bytes of a header that has not
// been parsed, or the encoded cards of a parsed (possibly edited) one,
// i.e. the bytes a Flush would write (encodedCards), so the answer is the
// same before and after Header() and before and after a Flush. bp is the
// image's BITPIX (ZBITPIX for a tile-compressed image).
//
//   - Which card: a card belongs to a keyword when its name (columns 1 up
//     to the first blank or '=', cfitsio's ffgknm) is that keyword in any
//     case. Of several, the LAST is used: cfitsio 4.7.0 (ffgphd) and
//     astropy 8.0.2 both use it on a plain image (keywords/dup_bscale,
//     dup_bzero, dup_blank, dup_bscale_bad_first; imaging/card_*_dup*).
//     With firstWins (a tile-compressed image) the FIRST is used, because
//     there cfitsio (imcomp_get_compressed_image_par: ffgky) and astropy
//     both use the first (imaging/comp_card_i16_bscale_dup_2d.fits.fz).
//   - That card must be a valid value card: the keyword upper case in
//     columns 1-8, "= " in columns 9-10, and a value that is a FITS
//     integer or real (strictFITSNumber; for BLANK an integer,
//     strictFITSInteger) which cfitsio converts to the same value (ffc2dd /
//     ffc2jj: no strtod range error, so no subnormal or underflowing
//     BSCALE). A BLANK must also lie in BITPIX's stored range (cfitsio
//     compares a BITPIX 32 BLANK as a 32-bit value, astropy as the full
//     integer) and must not be 1234554321 (cfitsio's NULL_UNDEFINED: no
//     BLANK; astropy applies it). Anything else is an error wrapping
//     ErrInvalidScaling: the references read every such card differently
//     from each other, refuse the file, or silently drop it.
//   - A HIERARCH card whose name is BSCALE, BZERO or BLANK (HIERARCH
//     BSCALE = 2.0) is an error: cfitsio applies it as scaling, astropy
//     ignores it (keywords/bscale_hierarch).
//   - BLANK is read on integer BITPIX only (it is undefined for float data,
//     FITS 4.0 §4.4.2.4); on float BITPIX a BLANK card, valid or not, is
//     not read, as both references ignore it.
//
// Proven per sample against cfitsio's and astropy's readings
// (testdata/images/golden/keywords_rule.json, TestScalingKeywordsRule;
// testdata/imaging/readfloat32_expected.json refs,
// TestReadFloat32_ScalingCardVariants). The errors are not prefixed; the
// caller names itself.
func scalingCards(cards func(yield func([]byte) bool), bp bitpix.BITPIX, firstWins bool) (imageCards, error) {
	withBlank := !bp.IsFloat()
	var pick [3][]byte
	var hierarch []byte
	out := imageCards{bscale: 1}
	cards(func(raw []byte) bool {
		if string(raw[:8]) == "END     " {
			return false
		}
		t := strings.TrimRight(string(raw), " ")
		isH := strings.HasPrefix(t, "HIERARCH ")
		name := cardKeyName(t)
		if isH {
			name = cfKeyName(t)
		}
		for i, k := range cardKeys {
			if !strings.EqualFold(name, k) || (i == keyBlank && !withBlank) {
				continue
			}
			if isH {
				if hierarch == nil {
					hierarch = raw
				}
				continue
			}
			out.count[i]++
			if !firstWins || pick[i] == nil {
				pick[i] = raw
			}
		}
		return true
	})
	if hierarch != nil {
		return imageCards{}, fmt.Errorf("%w: %q: a HIERARCH card named like a scaling keyword (cfitsio applies it, astropy ignores it on a plain image): %w",
			ErrInvalidScaling, strings.TrimRight(string(hierarch), " "), header.ErrInvalidKeyword)
	}
	for i, b := range pick {
		if b == nil {
			continue
		}
		key := cardKeys[i]
		if pad := fmt.Sprintf("%-8s", key); string(b[:8]) != pad || b[8] != '=' || b[9] != ' ' {
			err := fmt.Errorf("%w: %s card %q is not a value card (keyword in columns 1-8, \"= \" in columns 9-10)",
				ErrInvalidScaling, key, strings.TrimRight(string(b), " "))
			if string(b[:8]) != pad && strings.EqualFold(string(b[:8]), pad) {
				// A lower-case keyword is not a FITS keyword (FITS 4.0
				// sect. 4.1.2.1), as the writers report it.
				err = fmt.Errorf("%w: %w", err, header.ErrInvalidKeyword)
			}
			return imageCards{}, err
		}
		// A quoted string ('10.0') is a character string (FITS 4.0
		// §4.2.1), never a number: the strict parsers refuse the quote.
		v := string(b[10:])
		if j := strings.IndexByte(v, '/'); j >= 0 {
			v = v[:j]
		}
		text := strings.TrimSpace(v)
		if i == keyBlank {
			n, ok := strictFITSInteger(text)
			if c, cok := cfC2JJ(text); !ok || !cok || c != n {
				return imageCards{}, fmt.Errorf("%w: %s = %s is not an integer", ErrInvalidScaling, key, text)
			}
			if lo, hi := storedRange(bp); n < lo || n > hi {
				return imageCards{}, fmt.Errorf("%w: %s = %d is outside the BITPIX %d range [%d, %d] (cfitsio and astropy compare it differently)",
					ErrInvalidScaling, key, n, int(bp), lo, hi)
			}
			if n == nullUndefined {
				return imageCards{}, fmt.Errorf("%w: %s = %d is cfitsio's NULL_UNDEFINED, which cfitsio reads as no BLANK and astropy as the null value",
					ErrInvalidScaling, key, n)
			}
			out.blank, out.haveBlank = n, true
			continue
		}
		f, ok := strictFITSNumber(text)
		if !ok {
			return imageCards{}, fmt.Errorf("%w: %s = %s is not a finite FITS number", ErrInvalidScaling, key, text)
		}
		if c, cok := cfC2DD(text); !cok || c != f {
			return imageCards{}, fmt.Errorf("%w: %s = %s: cfitsio's strtod rejects it (range error) and would use the default, astropy reads the value",
				ErrInvalidScaling, key, text)
		}
		if i == keyBscale {
			out.bscale = f
		} else {
			out.bzero = f
		}
	}
	return out, nil
}

// cardKeyName is cfitsio's ffgknm as far as it matters here: columns 1
// up to the first blank or '=' ("BZERO=10.0" names BZERO; a HIERARCH card
// names HIERARCH, see cfKeyName). ffgknm also stops at 74 characters,
// which no name compared with BSCALE / BZERO / BLANK reaches.
func cardKeyName(t string) string {
	if i := strings.IndexAny(t, " ="); i >= 0 {
		return t[:i]
	}
	return t
}

// headerCards returns the card images scalingCards reads for r: the
// parsed header's encoded cards when it has been loaded (and possibly
// edited in memory), else the stored cards, so a malformed unrelated card
// does not stop a read that never parses the header.
func (r *hduRecord) headerCards() func(yield func([]byte) bool) {
	if r.parsed != nil {
		return encodedCards(r.parsed, nil)
	}
	return rawCards(r.rawHeader)
}

// imageCards returns the image's BSCALE / BZERO / BLANK by scalingCards.
func (r *hduRecord) imageCards() (imageCards, error) {
	return scalingCards(r.headerCards(), r.bitpix, false)
}

// imageKeywords returns the scaling and null value of the image
// (scalingCards).
func (r *hduRecord) imageKeywords() (s scaling, blank int64, hasBlank bool, err error) {
	c, err := r.imageCards()
	if err != nil {
		return scaling{}, 0, false, err
	}
	return newScaling(c.bscale, c.bzero), c.blank, c.haveBlank, nil
}
