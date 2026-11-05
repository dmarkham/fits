package fits

import (
	"errors"
	"fmt"
	"math"
	"math/big"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/tform"
)

// ErrNullValue is returned by the non-masked table readers when a cell is
// undefined (an integer equal to TNULLn, an undefined logical, an ASCII
// field equal to TNULLn) and the requested Go type has no way to say so.
// Float targets receive NaN instead. Use the *Masked readers to get the
// null flags explicitly.
var ErrNullValue = errors.New("fits: undefined (null) table value")

// tval is one decoded numeric table element after TSCAL/TZERO: an exact
// signed or unsigned integer, or a float64.
type tval struct {
	kind tvalKind
	i    int64
	u    uint64
	f    float64
	null bool
	// noInt: no integer target can hold the value, although f (its
	// nearest float64) may be integral: an integer outside the int64 and
	// uint64 ranges, or an ASCII decimal with a fraction.
	noInt bool
	// On the float64 scaling path f is the float64 arithmetic result
	// (what float targets get) and the exact value, which integer
	// targets need, is stored*TSCAL + TZERO with the exact literals of
	// sc: src says where the stored value is (srcInt: si, srcFloat: sf,
	// srcRat: sr).
	src tvalSrc
	si  int64
	sf  float64
	sr  *big.Rat
	sc  *colScaling
	// On the exact paths a float value is the exact physical value
	// correctly rounded, and float32 targets need the float32 nearest
	// that exact value (not the float32 of f, which would round twice):
	// xr holds it (an integer outside int64/uint64, an ASCII decimal with
	// a fraction under TSCAL 1 and an integer TZERO), or exactF says it is
	// exact() (an E/D value under TSCAL 1 and an integer TZERO).
	xr     *big.Rat
	exactF bool
}

type tvalKind uint8

const (
	tvInt tvalKind = iota
	tvUint
	tvFloat
)

type tvalSrc uint8

const (
	srcNone tvalSrc = iota
	srcInt
	srcFloat
	srcRat
)

// colScaling describes how stored values map to physical values for one
// column: phys = stored*TSCAL + TZERO, with the TSCALn/TZEROn literals
// read exactly.
type colScaling struct {
	scale, zero float64
	identity    bool     // TSCAL 1, TZERO 0
	intZero     bool     // TSCAL 1 and the TZERO literal is an integer that fits int64
	zeroI       int64    // that integer
	zero2p63    bool     // TSCAL 1 and TZERO exactly 2^63 (uint64 convention on K)
	bigZero     *big.Int // TSCAL 1 and an integer TZERO outside int64 (incl. 2^63 off K)
	scaleR      *big.Rat // the exact literals (float64 path, integer targets)
	zeroR       *big.Rat
	ambiguous   bool // TZERO is a float literal near, not at, 2^63: integer targets refused
	zeroIsF64   bool // the TZERO literal is exactly the float64 zero
	hasNull     bool
	nullInt     int64
}

func columnScaling(c Column) colScaling {
	s := colScaling{scale: c.Scale, zero: c.Zero, hasNull: c.hasNull, nullInt: c.nullInt, ambiguous: c.zeroNear2p63}
	if s.scale == 0 {
		s.scale = 1
	}
	s.scaleR, s.zeroR = c.scaleLit, c.zeroLit
	if s.scaleR == nil {
		s.scaleR = new(big.Rat).SetFloat64(s.scale)
	}
	if s.zeroR == nil {
		s.zeroR = new(big.Rat).SetFloat64(s.zero)
	}
	if s.scaleR.Cmp(ratOne) != 0 || !s.zeroR.IsInt() || s.ambiguous {
		return s // float64 path
	}
	z := s.zeroR.Num()
	s.zeroIsF64 = new(big.Rat).SetFloat64(s.zero).Cmp(s.zeroR) == 0
	switch {
	case z.Sign() == 0:
		s.identity = true
	case z.IsInt64():
		s.intZero, s.zeroI = true, z.Int64()
	default:
		s.zero2p63 = z.Cmp(big2p63) == 0
		s.bigZero = new(big.Int).Set(z)
	}
	return s
}

var (
	ratOne  = big.NewRat(1, 1)
	big2p63 = new(big.Int).Lsh(big.NewInt(1), 63)
)

// decodeInt applies null detection and scaling to a stored integer.
// is64 says whether the stored type is K (the only type where the uint64
// TZERO = 2^63 convention applies).
func (s *colScaling) decodeInt(raw int64, is64 bool) tval {
	if s.hasNull && raw == s.nullInt {
		return tval{null: true}
	}
	switch {
	case s.identity:
		return tval{kind: tvInt, i: raw}
	case s.zero2p63 && is64:
		return tval{kind: tvUint, u: uint64(raw) ^ (1 << 63)}
	case s.intZero:
		z := s.zeroI
		if z >= 0 && raw > math.MaxInt64-z {
			// The true value is in (2^63-1, 2^64): exact as a uint64.
			return tval{kind: tvUint, u: uint64(raw) + uint64(z)}
		}
		if z < 0 && raw < math.MinInt64-z {
			// Below -2^63: only a float can hold it, rounded once from
			// the exact sum (float64(raw)+zero would round twice).
			return bigTval(new(big.Int).Add(big.NewInt(raw), big.NewInt(z)))
		}
		return tval{kind: tvInt, i: raw + z}
	case s.bigZero != nil:
		return bigTval(new(big.Int).Add(big.NewInt(raw), s.bigZero))
	}
	return tval{kind: tvFloat, f: float64(float64(raw)*s.scale) + s.zero, src: srcInt, si: raw, sc: s}
}

// bigTval returns an exact integer as a tval: int64 or uint64 when it
// fits, otherwise its correctly rounded float64 marked noInt.
func bigTval(b *big.Int) tval {
	switch {
	case b.IsInt64():
		return tval{kind: tvInt, i: b.Int64()}
	case b.IsUint64():
		return tval{kind: tvUint, u: b.Uint64()}
	}
	r := new(big.Rat).SetInt(b)
	f, _ := r.Float64()
	return tval{kind: tvFloat, f: f, noInt: true, xr: r}
}

// ratTval returns an exact non-integer (or out-of-range) value on an exact
// path: its correctly rounded float64, never an integer target.
func ratTval(r *big.Rat) tval {
	f, _ := r.Float64()
	return tval{kind: tvFloat, f: f, noInt: true, xr: r}
}

// decodeFloat applies scaling to a stored float; NaN is the null value.
// ±Inf is a value (not null). Under TSCAL 1 and an integer TZERO the
// value is stored + TZERO correctly rounded (one IEEE addition when the
// TZERO literal is a double; the exact sum otherwise, since adding
// TZERO's double would round twice); otherwise it is the float64
// arithmetic float64(stored)*TSCAL + TZERO.
func (s *colScaling) decodeFloat(v float64) tval {
	if math.IsNaN(v) {
		return tval{kind: tvFloat, f: v, null: true}
	}
	if s.identity {
		return tval{kind: tvFloat, f: v}
	}
	t := tval{kind: tvFloat, src: srcFloat, sf: v, sc: s}
	if s.intZero || s.bigZero != nil {
		t.exactF = true
		if s.zeroIsF64 || math.IsInf(v, 0) {
			t.f = v + s.zero
		} else {
			t.f, _ = t.exact().Float64()
		}
		return t
	}
	t.f = float64(v*s.scale) + s.zero
	return t
}

// asFloat32 returns v for a float32 target: the float32 nearest the exact
// value on the exact paths, the float32 of the float64 arithmetic result
// on the float64 path.
func (v tval) asFloat32() float32 {
	switch v.kind {
	case tvInt:
		return float32(v.i)
	case tvUint:
		return float32(v.u)
	}
	r := v.xr
	if r == nil && v.exactF {
		r = v.exact()
	}
	if r != nil {
		f, _ := r.Float32()
		return f
	}
	return float32(v.f)
}

// exact returns the exact physical value of a float64-path element, or
// nil when it has none (a non-finite stored value, or an ambiguous TZERO).
func (v tval) exact() *big.Rat {
	if v.sc == nil || v.sc.ambiguous {
		return nil
	}
	var r *big.Rat
	switch v.src {
	case srcInt:
		r = new(big.Rat).SetInt64(v.si)
	case srcFloat:
		if math.IsInf(v.sf, 0) || math.IsNaN(v.sf) {
			return nil
		}
		r = new(big.Rat).SetFloat64(v.sf)
	case srcRat:
		r = new(big.Rat).Set(v.sr)
	default:
		return nil
	}
	r.Mul(r, v.sc.scaleR)
	return r.Add(r, v.sc.zeroR)
}

// decodeElem decodes one stored element of a numeric binary type.
func (s *colScaling) decodeElem(t tform.BinaryType, src []byte) tval {
	switch t {
	case tform.BinUint8:
		return s.decodeInt(int64(src[0]), false)
	case tform.BinInt16:
		return s.decodeInt(int64(bigendian.Int16(src)), false)
	case tform.BinInt32:
		return s.decodeInt(int64(bigendian.Int32(src)), false)
	case tform.BinInt64:
		return s.decodeInt(bigendian.Int64(src), true)
	case tform.BinFloat32:
		return s.decodeFloat(float64(bigendian.Float32(src)))
	case tform.BinFloat64:
		return s.decodeFloat(bigendian.Float64(src))
	}
	panic("fits: decodeElem on non-numeric type " + string(rune(t)))
}

func isNumericBinType(t tform.BinaryType) bool {
	switch t {
	case tform.BinUint8, tform.BinInt16, tform.BinInt32, tform.BinInt64,
		tform.BinFloat32, tform.BinFloat64:
		return true
	}
	return false
}

func (v tval) String() string {
	switch v.kind {
	case tvInt:
		return fmt.Sprint(v.i)
	case tvUint:
		return fmt.Sprint(v.u)
	}
	return fmt.Sprint(v.f)
}

func (v tval) asFloat() float64 {
	switch v.kind {
	case tvInt:
		return float64(v.i)
	case tvUint:
		return float64(v.u)
	}
	return v.f
}

// asInt returns v as an int64 if it is an integer in [lo, hi]. On the
// float64 scaling path the exact value decides, not the float64.
func (v tval) asInt(lo, hi int64) (int64, bool) {
	if v.noInt || (v.sc != nil && v.sc.ambiguous) {
		return 0, false
	}
	if v.src != srcNone {
		r := v.exact()
		if r == nil || !r.IsInt() || !r.Num().IsInt64() {
			return 0, false
		}
		n := r.Num().Int64()
		return n, n >= lo && n <= hi
	}
	switch v.kind {
	case tvInt:
		return v.i, v.i >= lo && v.i <= hi
	case tvUint:
		return int64(v.u), v.u <= uint64(hi) && hi >= 0
	}
	f := v.f
	if f != math.Trunc(f) || math.IsInf(f, 0) || math.IsNaN(f) {
		return 0, false
	}
	// float64(hi) may round up (2^63 for MaxInt64), so compare with >= on
	// the exclusive bound float64(hi)+1 (2^63 for MaxInt64, exactly).
	if f < float64(lo) || f >= float64(hi)+1 {
		return 0, false
	}
	return int64(f), true
}

// asUint returns v as a uint64 if it is an integer in [0, hi]. On the
// float64 scaling path the exact value decides, not the float64.
func (v tval) asUint(hi uint64) (uint64, bool) {
	if v.noInt || (v.sc != nil && v.sc.ambiguous) {
		return 0, false
	}
	if v.src != srcNone {
		r := v.exact()
		if r == nil || !r.IsInt() || !r.Num().IsUint64() {
			return 0, false
		}
		n := r.Num().Uint64()
		return n, n <= hi
	}
	switch v.kind {
	case tvInt:
		return uint64(v.i), v.i >= 0 && uint64(v.i) <= hi
	case tvUint:
		return v.u, v.u <= hi
	}
	f := v.f
	if f != math.Trunc(f) || math.IsInf(f, 0) || math.IsNaN(f) || f < 0 {
		return 0, false
	}
	if hi == math.MaxUint64 {
		if f >= 18446744073709551616.0 {
			return 0, false
		}
	} else if f > float64(hi) {
		return 0, false
	}
	return uint64(f), true
}
