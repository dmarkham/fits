package fits

import (
	"fmt"
	"math"

	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/bitpix"
)

// pixConv converts stored pixels of one BITPIX/scaling to target type T.
// Callers must have checked canHold[T] first; for integer T that guarantees
// every physical value is an exact integer inside T's range.
type pixConv[T Numeric] struct {
	bp       bitpix.BITPIX
	s        scaling
	size     int
	tInt     bool // T is an integer type
	unsigned bool // T is an unsigned integer type
	u64conv  bool // BITPIX 64 with BSCALE 1, BZERO 2^63
}

func newPixConv[T Numeric](bp bitpix.BITPIX, s scaling) *pixConv[T] {
	_, _, isInt, unsigned := intRange[T]()
	return &pixConv[T]{bp: bp, s: s, size: bp.Size(), tInt: isInt, unsigned: unsigned, u64conv: s.isU64Convention(bp)}
}

// storedInt reads the stored integer at b (integer BITPIX only).
func storedInt(b []byte, bp bitpix.BITPIX) int64 {
	switch bp {
	case bitpix.Int8:
		return int64(b[0]) // BITPIX 8 is unsigned on disk
	case bitpix.Int16:
		return int64(bigendian.Int16(b))
	case bitpix.Int32:
		return int64(bigendian.Int32(b))
	}
	return bigendian.Int64(b)
}

// storedFloat reads the stored IEEE value at b (float BITPIX only).
func storedFloat(b []byte, bp bitpix.BITPIX) float64 {
	if bp == bitpix.Float32 {
		return float64(bigendian.Float32(b))
	}
	return bigendian.Float64(b)
}

// one converts the stored element at b.
//
// Float targets reproduce cfitsio's conversions (getcol*.c): no arithmetic
// when BSCALE = 1 and BZERO = 0 (a direct integer-to-float conversion,
// rounded once), the exact unsigned value for the uint64 convention, and
// otherwise stored*BSCALE + BZERO in float64 followed by one conversion.
// Integer targets use exact integer arithmetic.
func (c *pixConv[T]) one(b []byte) T {
	if c.bp.IsFloat() {
		f := storedFloat(b, c.bp)
		if !c.s.identity {
			f = float64(f*c.s.bscale) + c.s.bzero
		}
		return T(f)
	}
	r := storedInt(b, c.bp)
	if c.tInt {
		v := c.s.physInt(r)
		if c.unsigned {
			return T(v.lo)
		}
		return T(int64(v.lo))
	}
	switch {
	case c.s.identity:
		return T(r)
	case c.u64conv:
		return T(uint64(r) ^ (1 << 63))
	}
	return T(float64(float64(r)*c.s.bscale) + c.s.bzero)
}

func (c *pixConv[T]) all(raw []byte, out []T) {
	if c.s.identity {
		// No scaling: one() is a single conversion of the stored value
		// (for an integer T, canHold has checked that every stored value
		// fits), done here in one loop per BITPIX.
		switch c.bp {
		case bitpix.Int8:
			for i := range out {
				out[i] = T(raw[i])
			}
		case bitpix.Int16:
			for i := range out {
				out[i] = T(bigendian.Int16(raw[i*2:]))
			}
		case bitpix.Int32:
			for i := range out {
				out[i] = T(bigendian.Int32(raw[i*4:]))
			}
		case bitpix.Int64:
			for i := range out {
				out[i] = T(bigendian.Int64(raw[i*8:]))
			}
		case bitpix.Float32:
			for i := range out {
				// Through float64 exactly as one() (a NaN's quiet bit
				// included).
				out[i] = T(float64(bigendian.Float32(raw[i*4:])))
			}
		case bitpix.Float64:
			for i := range out {
				out[i] = T(bigendian.Float64(raw[i*8:]))
			}
		}
		return
	}
	for i := range out {
		out[i] = c.one(raw[i*c.size:])
	}
}

// valid reports whether the stored element is defined: integer BITPIX
// pixels equal to BLANK and float pixels that are NaN are undefined.
func (c *pixConv[T]) valid(b []byte, blank int64, hasBlank bool) bool {
	if c.bp.IsFloat() {
		return !math.IsNaN(storedFloat(b, c.bp))
	}
	return !hasBlank || storedInt(b, c.bp) != blank
}

// pixEnc converts values of type T to the stored representation of an
// image with BITPIX bp and scaling s (the inverse of pixConv), following
// cfitsio's writers (putcol*.c): stored = (physical - BZERO) / BSCALE,
// rounded half away from zero for integer BITPIX. When T is an integer and
// BSCALE = 1 with an integral BZERO (which includes the
// int8/uint16/uint32/uint64 conventions) the subtraction is exact integer
// arithmetic. A value that cannot be stored (see put) is reported, never
// stored as a clamped or infinite value: the writers refuse the whole
// write with ErrNumOverflow.
type pixEnc[T Numeric] struct {
	bp       bitpix.BITPIX
	s        scaling
	size     int
	native   bool // T's bit layout is the stored layout and no scaling
	exactInt bool
	unsigned bool
	lo, hi   int64
}

func newPixEnc[T Numeric](bp bitpix.BITPIX, s scaling) *pixEnc[T] {
	_, _, isInt, unsigned := intRange[T]()
	e := &pixEnc[T]{bp: bp, s: s, size: bp.Size(), unsigned: unsigned}
	e.lo, e.hi = storedRange(bp)
	if s.identity {
		var z T
		switch any(z).(type) {
		case uint8:
			e.native = bp == bitpix.Int8
		case int16:
			e.native = bp == bitpix.Int16
		case int32:
			e.native = bp == bitpix.Int32
		case int64:
			e.native = bp == bitpix.Int64
		case float32:
			e.native = bp == bitpix.Float32
		case float64:
			e.native = bp == bitpix.Float64
		}
	}
	e.exactInt = isInt && bp.IsInt() && s.exact && s.mul == 1
	return e
}

// put stores v at dst and reports whether it could be stored. It cannot
// (and dst is left as is) when, for integer BITPIX, the stored value lies
// outside the stored range (by more than cfitsio's 0.49 rounding margin
// in the double arithmetic, which is where cfitsio reports NUM_OVERFLOW;
// or exactly, on the integer path), and when, for float BITPIX, a finite
// v scales to a value that is not finite in the stored type (cfitsio
// stores the infinity silently).
func (e *pixEnc[T]) put(dst []byte, v T) bool {
	if e.bp.IsFloat() {
		if e.native {
			if e.bp == bitpix.Float32 {
				bigendian.PutFloat32(dst, float32(v)) // keeps NaN payloads
			} else {
				bigendian.PutFloat64(dst, float64(v))
			}
			return true
		}
		in := float64(v)
		f := (in - e.s.bzero) / e.s.bscale
		finite := !math.IsInf(in, 0) && !math.IsNaN(in)
		if e.bp == bitpix.Float32 {
			f32 := float32(f)
			if finite && math.IsInf(float64(f32), 0) {
				return false
			}
			bigendian.PutFloat32(dst, f32)
		} else {
			if finite && math.IsInf(f, 0) {
				return false
			}
			bigendian.PutFloat64(dst, f)
		}
		return true
	}
	var raw int64
	switch {
	case e.native:
		raw = int64(v)
	case e.exactInt:
		var x i128
		if e.unsigned {
			x = i128FromU64(uint64(v))
		} else {
			x = i128From64(int64(v))
		}
		x = x.sub(e.s.off)
		if x.cmp(i128From64(e.lo)) < 0 || x.cmp(i128From64(e.hi)) > 0 {
			return false
		}
		raw = int64(x.lo)
	default:
		var ok bool
		if raw, ok = roundStored((float64(v)-e.s.bzero)/e.s.bscale, e.lo, e.hi); !ok {
			return false
		}
	}
	switch e.bp {
	case bitpix.Int8:
		dst[0] = byte(raw)
	case bitpix.Int16:
		bigendian.PutInt16(dst, int16(raw))
	case bitpix.Int32:
		bigendian.PutInt32(dst, int32(raw))
	case bitpix.Int64:
		bigendian.PutInt64(dst, raw)
	}
	return true
}

// roundStored is cfitsio's double-to-integer store: d rounded half away
// from zero; ok is false where cfitsio reports NUM_OVERFLOW (d beyond the
// range by more than 0.49, where cfitsio would store the range limit).
// d is never NaN: it is (v - BZERO) / BSCALE for an integer v, and the
// writers refuse a BSCALE of 0 and non-finite BSCALE/BZERO
// (ErrInvalidScaling), so d is finite or an overflow to +-Inf.
func roundStored(d float64, lo, hi int64) (int64, bool) {
	flo, fhi := float64(lo), float64(hi)
	if d < flo-0.49 || d > fhi+0.49 {
		return 0, false
	}
	if d >= 0 {
		d = math.Trunc(d + 0.5)
	} else {
		d = math.Trunc(d - 0.5)
	}
	// For BITPIX 64, float64(hi) is 2^63, which int64 cannot hold: d may
	// be exactly 2^63 here (cfitsio's own range test, d > DLONGLONG_MAX,
	// lets it through to an out-of-range cast that stores MinInt64 without
	// reporting an overflow; we refuse it). For BITPIX 8/16/32 float64(hi)
	// is exact and d <= hi after rounding. No lower test is needed: d >=
	// flo-0.49 rounds to at least flo, and float64(lo) is exact for every
	// BITPIX (the next float64 below -2^63 is 2048 away).
	if d > fhi || d == fhi && hi == math.MaxInt64 {
		return 0, false
	}
	return int64(d), true
}

// encode stores data. bad is the index of the first value that cannot be
// stored (-1 if none) and n the number of such values; buf is then not a
// valid data unit and must not be written.
func (e *pixEnc[T]) encode(data []T) (buf []byte, bad, n int) {
	buf = make([]byte, len(data)*e.size)
	bad = -1
	if e.native {
		// T's layout is the stored one and there is no scaling: every
		// value can be stored, with the bytes put would store.
		e.encodeNative(buf, data)
		return buf, bad, 0
	}
	for i, v := range data {
		if !e.put(buf[i*e.size:], v) {
			if n == 0 {
				bad = i
			}
			n++
		}
	}
	return buf, bad, n
}

// encodeNative is put's e.native case over all of data, one loop per
// BITPIX (put's per-element dispatch was most of the cost of a write).
func (e *pixEnc[T]) encodeNative(buf []byte, data []T) {
	switch e.bp {
	case bitpix.Int8:
		for i, v := range data {
			buf[i] = byte(v)
		}
	case bitpix.Int16:
		for i, v := range data {
			bigendian.PutInt16(buf[i*2:], int16(v))
		}
	case bitpix.Int32:
		for i, v := range data {
			bigendian.PutInt32(buf[i*4:], int32(v))
		}
	case bitpix.Int64:
		for i, v := range data {
			bigendian.PutInt64(buf[i*8:], int64(v))
		}
	case bitpix.Float32:
		for i, v := range data {
			bigendian.PutFloat32(buf[i*4:], float32(v)) // keeps NaN payloads
		}
	case bitpix.Float64:
		for i, v := range data {
			bigendian.PutFloat64(buf[i*8:], float64(v))
		}
	}
}

// encodeOrOverflow stores data, or returns an error wrapping
// ErrNumOverflow that names the first value that cannot be stored.
func encodeOrOverflow[T Numeric](bp bitpix.BITPIX, s scaling, data []T) ([]byte, error) {
	buf, bad, n := newPixEnc[T](bp, s).encode(data)
	if n > 0 {
		return nil, fmt.Errorf("%w: %d of %d values cannot be stored with BITPIX %d, BSCALE %v, BZERO %v; the first is element %d = %v",
			ErrNumOverflow, n, len(data), int(bp), s.bscale, s.bzero, bad, data[bad])
	}
	return buf, nil
}
