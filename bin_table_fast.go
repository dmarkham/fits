package fits

import (
	"math"
	"reflect"
	"unsafe"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/tform"
)

// Fast paths for the numeric table readers (ReadColumn, ReadVectorColumn,
// ReadVarColumn and their Masked forms; ReadVarColumn[uint8] is how the
// tile-compressed reader fetches each tile's bytes). The general decoder (colScaling.decodeElem,
// numTo) builds a tval per element and, for integer targets under a TSCAL
// other than 1, an exact big.Rat; that costs several times the read
// itself. cellFast handles the common cases with the same arithmetic and
// returns false for every element it does not handle, which the caller
// then decodes the general way, so results (values, null flags and
// errors) are unchanged:
//
//   - integer stored types (B, I, J, K): a null (TNULLn) element is left
//     to the general path. With TSCAL and TZERO literals that are both
//     integers in int64 (TSCAL 1 and TZERO 0 included), the physical value
//     stored*TSCAL+TZERO is computed exactly in 128 bits; an integer
//     target gets it when it is in the target's range, a float target
//     gets its single conversion to float64 / float32 when TSCAL is 1 (as
//     decodeInt: one rounding of the exact sum). A float target under any
//     other scaling without a big TZERO gets decodeInt's float64
//     expression float64(stored)*TSCAL + TZERO, converted once.
//   - float stored types (E, D) under TSCAL 1 and TZERO 0: a float target
//     gets the stored value (NaN, the null value, goes to the general
//     path).
type fastDec struct {
	ok        bool
	intStored bool
	// integer targets: exact affine map with int64 m, z
	affine bool
	m, z   int64
	// float targets from integer stored values
	sumExact bool // TSCAL 1, integer TZERO in int64: float(raw+z)
	formula  bool // float64(float64(raw)*scale) + zero
	scale    float64
	zero     float64
	hasNull  bool
	nullInt  int64
	kind     reflect.Kind
	lo, hi   int64 // integer target range (hi = MaxInt64 for uint64)
}

func newFastDec(sc *colScaling, bt tform.BinaryType, kind reflect.Kind) fastDec {
	d := fastDec{kind: kind, hasNull: sc.hasNull, nullInt: sc.nullInt, scale: sc.scale, zero: sc.zero}
	isFloat := kind == reflect.Float32 || kind == reflect.Float64
	switch bt {
	case tform.BinUint8, tform.BinInt16, tform.BinInt32, tform.BinInt64:
		d.intStored = true
	case tform.BinFloat32, tform.BinFloat64:
		// Float stored values: only the unscaled float-target case.
		d.ok = isFloat && sc.identity
		return d
	default:
		return d
	}
	if isFloat {
		switch {
		case sc.identity:
			d.sumExact = true
		case sc.intZero:
			d.sumExact, d.z = true, sc.zeroI
		case sc.bigZero == nil:
			d.formula = true
		default:
			return d
		}
		d.m = 1
		d.ok = true
		return d
	}
	if sc.ambiguous || !sc.scaleR.IsInt() || !sc.zeroR.IsInt() ||
		!sc.scaleR.Num().IsInt64() || !sc.zeroR.Num().IsInt64() {
		return d
	}
	d.affine, d.m, d.z = true, sc.scaleR.Num().Int64(), sc.zeroR.Num().Int64()
	switch kind {
	case reflect.Int8:
		d.lo, d.hi = math.MinInt8, math.MaxInt8
	case reflect.Int16:
		d.lo, d.hi = math.MinInt16, math.MaxInt16
	case reflect.Int32:
		d.lo, d.hi = math.MinInt32, math.MaxInt32
	case reflect.Int64:
		d.lo, d.hi = math.MinInt64, math.MaxInt64
	case reflect.Uint8:
		d.lo, d.hi = 0, math.MaxUint8
	case reflect.Uint16:
		d.lo, d.hi = 0, math.MaxUint16
	case reflect.Uint32:
		d.lo, d.hi = 0, math.MaxUint32
	case reflect.Uint64:
		d.lo, d.hi = 0, math.MaxInt64 // larger values: general path
	default:
		return d
	}
	d.ok = true
	return d
}

func storedBinInt(bt tform.BinaryType, b []byte) int64 {
	switch bt {
	case tform.BinUint8:
		return int64(b[0])
	case tform.BinInt16:
		return int64(bigendian.Int16(b))
	case tform.BinInt32:
		return int64(bigendian.Int32(b))
	}
	return bigendian.Int64(b)
}

// cellFast decodes one element at src into *dst. It returns false, and
// leaves *dst alone, when the element needs the general path.
func cellFast[T any](d *fastDec, bt tform.BinaryType, src []byte, dst *T) bool {
	p := unsafe.Pointer(dst)
	if !d.intStored {
		var f float64
		if bt == tform.BinFloat32 {
			f = float64(bigendian.Float32(src))
		} else {
			f = bigendian.Float64(src)
		}
		if f != f { // NaN: null
			return false
		}
		if d.kind == reflect.Float64 {
			*(*float64)(p) = f
		} else {
			*(*float32)(p) = float32(f)
		}
		return true
	}
	raw := storedBinInt(bt, src)
	if d.hasNull && raw == d.nullInt {
		return false
	}
	if d.formula {
		f := float64(float64(raw)*d.scale) + d.zero
		if d.kind == reflect.Float64 {
			*(*float64)(p) = f
		} else {
			*(*float32)(p) = float32(f)
		}
		return true
	}
	// Exact raw*m + z in int64.
	var v int64
	if d.m == 1 {
		s := raw + d.z
		if (d.z >= 0) != (s >= raw) { // overflow
			return false
		}
		v = s
	} else {
		x := mul64(raw, d.m).add(i128From64(d.z))
		if x.hi != int64(x.lo)>>63 { // outside int64
			return false
		}
		v = int64(x.lo)
	}
	switch d.kind {
	case reflect.Float64:
		*(*float64)(p) = float64(v)
		return true
	case reflect.Float32:
		*(*float32)(p) = float32(v)
		return true
	}
	if v < d.lo || v > d.hi {
		return false // the general path reports the mismatch
	}
	switch d.kind {
	case reflect.Int8:
		*(*int8)(p) = int8(v)
	case reflect.Int16:
		*(*int16)(p) = int16(v)
	case reflect.Int32:
		*(*int32)(p) = int32(v)
	case reflect.Int64:
		*(*int64)(p) = v
	case reflect.Uint8:
		*(*uint8)(p) = uint8(v)
	case reflect.Uint16:
		*(*uint16)(p) = uint16(v)
	case reflect.Uint32:
		*(*uint32)(p) = uint32(v)
	case reflect.Uint64:
		*(*uint64)(p) = uint64(v)
	}
	return true
}
