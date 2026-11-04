// Package crefgen regenerates the synthetic inputs of the cfitsio
// compression harness (compress/testdata/cref/cref.c) so they need not be
// committed. Every function is a line-for-line port of the C generator of
// the same name; the harness records the sha256 of each input it fed to
// libcfitsio, and the tests must reproduce that hash (In*SHA256) before
// they use an input, so a port that drifts by one bit fails loudly.
//
// Nothing here is truth: the inputs are just inputs. The truth is what
// cfitsio did with them (the committed goldens).
package crefgen

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"math"
)

// LCG is cref.c's lcg_t.
type LCG struct{ state uint64 }

const (
	lcgMul = 6364136223846793005
	lcgInc = 1442695040888963407
)

// NewLCG is lcg_init.
func NewLCG(seed uint64) *LCG { return &LCG{state: seed*lcgMul + lcgInc} }

// U32 is lcg_u32.
func (r *LCG) U32() uint32 {
	r.state = r.state*lcgMul + lcgInc
	return uint32(r.state >> 33)
}

// Uniform is lcg_uniform: [0, 1).
func (r *LCG) Uniform() float64 { return float64(r.U32()) / 4294967296.0 }

// twoPi is cref.c's 2.0 * 3.14159265358979323846 (exact: a power-of-two
// multiple of the double nearest pi).
const twoPi float64 = 2.0 * 3.14159265358979323846

// Gauss is lcg_gauss (Box-Muller, one sample per call). The explicit
// float64 conversions keep each product rounded on its own (no fused
// multiply-add), as gcc -O2 compiles the C.
func (r *LCG) Gauss() float64 {
	u1 := r.Uniform()
	u2 := r.Uniform()
	if u1 < 1e-12 {
		u1 = 1e-12
	}
	return float64(math.Sqrt(float64(-2.0*math.Log(u1)))) * math.Cos(float64(twoPi*u2))
}

// Range is lcg_range: uniform in [lo, hi]; the high word is drawn first.
func (r *LCG) Range(lo, hi int64) int64 {
	span := uint64(hi-lo) + 1
	h := uint64(r.U32())
	l := uint64(r.U32())
	return lo + int64(((h<<32)|l)%span)
}

func llround(x float64) int64 { return int64(math.Round(x)) }

func clamp(v, lo, hi int64) int64 {
	if v < lo {
		return lo
	}
	if v > hi {
		return hi
	}
	return v
}

// f32 rounds a double to float32 (the C cast).
func f32(x float64) float32 { return float32(x) }

var (
	nan32         = math.Float32frombits(0x7fc00000) // C NAN
	nullF float32 = -1.0e30
	// FloatNullValue is cfitsio's FLOATNULLVALUE (fitsio2.h).
	FloatNullValue float32 = -9.11912e-36
	// DenormNull is cref.c's DENORM_NULL, the denormal null value of
	// golden imgw/rice_f32_nullval_denorm_t40x5.
	DenormNull float32 = 1.04488e-40
)

// ---- gen cases (cref gen) ----

func fillGauss(seed uint64, n int64, f func(i int64, v float64) float32) []float32 {
	r := NewLCG(seed)
	out := make([]float32, n)
	for i := range out {
		v := float64(100.0 + float64(5.0*r.Gauss()))
		out[i] = f(int64(i), v)
	}
	return out
}

func fillRamp(n int64, lo, hi float64) []float32 {
	out := make([]float32, n)
	for i := int64(0); i < n; i++ {
		out[i] = f32(lo + float64(float64(i)*(hi-lo))/float64(n-1))
	}
	out[0] = f32(lo)
	out[n-1] = f32(hi)
	return out
}

func fillPatternExtreme(n int64, m float32) []float32 {
	out := make([]float32, n)
	for i := range out {
		if i%4 >= 2 {
			out[i] = 1
		}
	}
	out[33] = m
	return out
}

// Quantizer-limit constants of cref.c (Q_MARGIN, M_MARGIN, Q_HALF, M_HALF).
const (
	QMargin float32 = 0x1.ffffaap-1
	MMargin float32 = 0x1.0c615cp+32
	QHalf   float32 = 0x1.fffffap-1
	MHalf   float32 = 0x1.0c6132p+31
	// Q_MARGIN_FAR, M_MARGIN_FAR: a margin ratio whose centered extremes
	// stay clear of the sentinels (golden margin_far_64x1).
	QMarginFar float32 = 0x1.0000b8p+0
	MMarginFar float32 = 0x1.0c606ep+32
)

func fillRuns3(n int64) []float32 {
	r := NewLCG(47)
	out := make([]float32, n)
	for i := int64(0); i < n; i++ {
		v := float64(100.0 + float64(5.0*r.Gauss()))
		k := i % 23
		switch {
		case k >= 4 && k <= 6:
			out[i] = 97.25
		case k >= 14 && k <= 18:
			out[i] = 103.5
		default:
			out[i] = f32(v)
		}
	}
	return out
}

// Fill is the gen-case generator named fill (cref.c case_t.fill_name)
// for n pixels.
func Fill(fill string, n int64) ([]float32, error) {
	switch fill {
	case "fill_gauss_100_5":
		return fillGauss(42, n, func(_ int64, v float64) float32 { return f32(v) }), nil
	case "fill_gauss_with_nulls":
		return fillGauss(43, n, func(i int64, v float64) float32 {
			if i%37 == 0 {
				return nullF
			}
			return f32(v)
		}), nil
	case "fill_gauss_with_real_nans":
		return fillGauss(46, n, func(i int64, v float64) float32 {
			if i%41 == 0 {
				return nan32
			}
			return f32(v)
		}), nil
	case "fill_gauss_with_zeros":
		return fillGauss(44, n, func(i int64, v float64) float32 {
			if i%31 == 0 {
				return 0
			}
			return f32(v)
		}), nil
	case "fill_gradient_plus_noise":
		r := NewLCG(45)
		out := make([]float32, n)
		for i := int64(0); i < n; i++ {
			x, y := float64(i%64), float64(i/64)
			g := float64(float64(50.0+float64(0.5*x)) + float64(0.25*y))
			out[i] = f32(g + float64(2.0*r.Gauss()))
		}
		return out, nil
	case "fill_gauss_with_runs3":
		return fillRuns3(n), nil
	case "fill_runs3_with_nulls":
		out := fillRuns3(n)
		for i := range out {
			if i%29 == 11 {
				out[i] = nullF
			}
		}
		return out, nil
	case "fill_runs3_stale":
		r := NewLCG(48)
		out := make([]float32, n)
		for i := int64(0); i < n; i++ {
			row, col := i/200, i%200
			s := 5.0
			if row != 0 {
				s = 10.0
			}
			v := float64(100.0 + float64(s*r.Gauss()))
			if row == 1 && col%8 >= 2 && col%8 <= 4 {
				v = 100.0 + float64(col/8)
			}
			out[i] = f32(v)
		}
		return out, nil
	case "fill_constant":
		out := make([]float32, n)
		for i := range out {
			out[i] = 42.5
		}
		return out, nil
	case "fill_tiny":
		vals := [6]float32{1.0, 2.0, 1.5, 2.5, 1.25, 2.75}
		out := make([]float32, n)
		for i := range out {
			out[i] = vals[i%6]
		}
		return out, nil
	case "fill_tiny_4x5":
		out := make([]float32, n)
		for i := range out {
			out[i] = 10.0 + float32(i%3)
		}
		out[13] = 50
		out[17] = -7
		return out, nil
	case "fill_range_fsub":
		return fillRamp(n, -(math.Ldexp(1, -24) - 13*math.Ldexp(1, -32)), 1-math.Ldexp(1, -24)), nil
	case "fill_wide_center":
		return fillRamp(n, 0, 3.0e9), nil
	case "fill_center_fsum":
		return fillRamp(n, -120, 2147483520.0), nil
	case "fill_margin_pos":
		return fillPatternExtreme(n, MMargin), nil
	case "fill_halfmargin_pos":
		return fillPatternExtreme(n, MHalf), nil
	case "fill_margin_far":
		return fillPatternExtreme(n, MMarginFar), nil
	case "fill_n5_one":
		return fillLit(n, litN5One), nil
	case "fill_n5_two":
		return fillLit(n, litN5Two), nil
	case "fill_n5_zero":
		return fillLit(n, litN5Zero), nil
	case "fill_n1_short":
		return fillLit(n, litN1Short), nil
	case "fill_one_good":
		return fillLit(n, litOneGood), nil
	case "fill_clip_tie":
		out := make([]float32, n)
		for i := range out {
			if i%27 < 13 {
				out[i] = 26
			}
		}
		return out, nil
	case "fill_abs_nulls":
		return fillLit(n, litAbsNulls), nil
	case "fill_all_null":
		return fillLit(n, []float32{nullF}), nil
	case "fill_all_nan":
		return fillLit(n, []float32{nan32}), nil
	}
	return nil, fmt.Errorf("crefgen: unknown generator %q", fill)
}

// ---- encode cases (cref encode: rice / hcomp / plio) ----

func bpMin(bp int) int64 {
	switch bp {
	case 1:
		return 0
	case 2:
		return -32768
	}
	return math.MinInt32
}

func bpMax(bp int) int64 {
	switch bp {
	case 1:
		return 255
	case 2:
		return 32767
	}
	return math.MaxInt32
}

// IntKind is fill_int_kind: kind 0 extremes, 1 jumps, 2 smooth, 3 mixed,
// other zeros, for pixel width bp bytes.
func IntKind(n int64, bp, kind int, seed uint64) []int64 {
	o := make([]int64, n)
	r := NewLCG(seed)
	switch kind {
	case 0:
		p1 := [4]int64{0, 127, 128, 255}
		pn := [4]int64{0, bpMax(bp), bpMin(bp), -1}
		for i := range o {
			if bp == 1 {
				o[i] = p1[i%4]
			} else {
				o[i] = pn[i%4]
			}
		}
	case 1:
		for i := range o {
			o[i] = r.Range(bpMin(bp), bpMax(bp))
		}
	case 2:
		mid := int64(100000)
		if bp == 1 {
			mid = 100
		} else if bp == 2 {
			mid = 1000
		}
		for i := range o {
			o[i] = clamp(mid+int64(i%97)+llround(float64(3.0*r.Gauss())), bpMin(bp), bpMax(bp))
		}
	case 3:
		for i := range o {
			switch (i / 32) % 3 {
			case 0:
				if bp == 1 {
					o[i] = 7
				} else {
					o[i] = -5
				}
			case 1:
				base := int64(-20)
				if bp == 1 {
					base = 50
				}
				o[i] = clamp(base+llround(float64(2.0*r.Gauss())), bpMin(bp), bpMax(bp))
			default:
				o[i] = r.Range(bpMin(bp), bpMax(bp))
			}
		}
	}
	return o
}

// NoiseImage is gen_noise_image (HCOMPRESS kind 5).
func NoiseImage(n int64, bp int, seed uint64) []int64 {
	r := NewLCG(seed)
	mid, sig := 3.0e6, 4.0e4
	switch bp {
	case 1:
		mid, sig = 120, 10
	case 2:
		mid, sig = 1000, 50
	}
	o := make([]int64, n)
	for i := range o {
		o[i] = clamp(llround(float64(mid+float64(sig*r.Gauss()))+float64(i%40)), bpMin(bp), bpMax(bp))
	}
	return o
}

// Plio is gen_plio.
func Plio(n int64, kind int, seed uint64) []int64 {
	r := NewLCG(seed)
	o := make([]int64, n)
	big := [6]int64{16777215, 0, 4096, 1, 16777215, 8000000}
	dv := [14]int64{1, 4096, 8192, 4096, 8193, 4096, 8191, 4096, 1, 4097, 1, 4098, 4097, 1}
	for i := int64(0); i < n; i++ {
		switch kind {
		case 0:
			o[i] = [5]int64{0, 1, 2, 0, 7}[(i/13)%5]
		case 1:
			switch {
			case i < 5000:
				o[i] = 0
			case i < 9500:
				o[i] = 3
			case i < 9600:
				o[i] = 4100
			default:
				o[i] = 1
			}
		case 2:
			o[i] = big[(i/9)%6]
		case 3:
			if i%7 == 3 {
				o[i] = 1 + i%11
			}
		case 4:
		case 6:
			seg, pos := i/3, i%3
			ln := 1 + seg%3
			switch {
			case seg%5 == 4:
				o[i] = 0
			case pos < ln:
				o[i] = dv[seg%14]
			default:
				o[i] = dv[(seg+1)%14]
			}
		default:
			o[i] = r.Range(0, 16777215)
		}
	}
	return o
}

// ---- img cases (cref encode: whole images through fits_img_compress) ----

// TypeRange is type_range.
func TypeRange(t string) (lo, hi int64) {
	switch t {
	case "u8":
		return 0, 255
	case "s8":
		return -128, 127
	case "i16":
		return -32768, 32767
	case "u16":
		return 0, 65535
	case "i32":
		return math.MinInt32, math.MaxInt32
	}
	return 0, 4294967295
}

// TypeWidth is type_width.
func TypeWidth(t string) int {
	switch t {
	case "u8", "s8":
		return 1
	case "i16", "u16":
		return 2
	}
	return 4
}

// ImgInt is gen_img_int.
func ImgInt(n, nx int64, t string, kind int, seed uint64) []int64 {
	r := NewLCG(seed)
	lo, hi := TypeRange(t)
	mid := lo + (hi-lo)/2
	o := make([]int64, n)
	for i := int64(0); i < n; i++ {
		row, col := i/nx, i%nx
		if kind == 1 {
			top := hi
			if top > 16777215 {
				top = 16777215
			}
			vals := [6]int64{0, 1, 2, top, top/2 + 1, 0}
			if (row/3)%2 == 0 && col > 5 && col < 30 {
				o[i] = vals[(col/4)%6]
			}
			continue
		}
		if kind == 2 {
			sig := 8.0
			if hi-lo > 70000 {
				sig = 4.0e4
			} else if hi-lo > 300 {
				sig = 50
			}
			base := float64(mid)
			if hi-lo > 70000 {
				base = float64(mid) / 1000.0
			}
			if lo == 0 && hi-lo > 70000 {
				base = 3.0e6
			}
			o[i] = clamp(llround(float64(base+float64(sig*r.Gauss()))+float64(col)), lo, hi)
			continue
		}
		switch (row / 5) % 5 {
		case 0:
			o[i] = clamp(lo+(hi-lo)/3+col*3+row, lo, hi)
		case 1:
			o[i] = r.Range(lo, hi)
		case 2:
			o[i] = [4]int64{lo, hi, mid, lo + 1}[i%4]
		case 3:
			o[i] = clamp(mid+llround(float64(4.0*r.Gauss())), lo, hi)
		default:
			o[i] = mid
		}
	}
	return o
}

// ImgFloat is gen_img_float.
func ImgFloat(n, nx int64, kind int, seed uint64) []float32 {
	r := NewLCG(seed)
	switch kind {
	case 9:
		o := make([]float32, n)
		copy(o, fillPatternExtreme(nx, 0))
		o[5] = nan32
		copy(o[nx:], fillPatternExtreme(n-nx, MMargin))
		return o
	case 10:
		return fillRamp(n, 0, 3.0e9)
	case 11:
		return fillPatternExtreme(n, MHalf)
	case 19:
		o := make([]float32, n)
		for i := range o {
			o[i] = 100
		}
		return o
	case 20:
		return fillRamp(n, 1.0e10, 1.0e10+4294966784.0)
	case 21:
		o := make([]float32, n)
		copy(o, fillRamp(nx, -0x1.b0a9fep+23, 0x1.51b436p+32))
		for i := nx; i < n; i++ {
			o[i] = float32(i - nx)
		}
		o[nx+5] = nan32
		return o
	}
	o := make([]float32, n)
	for i := int64(0); i < n; i++ {
		row, col := i/nx, i%nx
		v := float64(float64(float64(100.0+float64(0.3*float64(col)))+float64(0.2*float64(row))) + float64(5.0*r.Gauss()))
		switch kind {
		case 1:
			if i%23 == 0 {
				o[i] = 0
			} else {
				o[i] = f32(v)
			}
		case 2:
			if i%31 == 0 {
				o[i] = nan32
			} else {
				o[i] = f32(v)
			}
		case 3:
			if row < 10 {
				o[i] = 42.5
			} else {
				o[i] = f32(v)
			}
		case 4, 6:
			if row < 5 {
				add := 0.0
				if i%2 == 1 {
					add = 1.0e20
				}
				o[i] = f32(add + v)
			} else {
				o[i] = f32(v)
			}
		default:
			o[i] = f32(v)
		}
		if kind == 12 {
			if i%37 == 5 {
				o[i] = nan32
			}
			if i%41 == 7 {
				o[i] = -1.0e30
			}
		}
		if kind == 13 {
			if i%200 == 3 {
				o[i] = nan32
			}
			if i%200 == 10 {
				o[i] = FloatNullValue
			}
		}
		if kind == 16 {
			if i%23 == 0 {
				o[i] = 0
			}
			if i%31 == 0 {
				o[i] = nan32
			}
		}
		if kind == 17 || kind == 18 {
			if i%200 == 3 {
				o[i] = nan32
			}
			if i%50 == 7 {
				if kind == 17 {
					o[i] = -999
				} else {
					o[i] = DenormNull
				}
			}
			if kind == 18 && i%50 == 9 {
				o[i] = -3.0e-42
			}
		}
		if kind == 7 || kind == 8 {
			if i == 3 {
				o[i] = nan32
			}
			if i == 20 {
				o[i] = -1.0e30
			}
			if kind == 8 && i == 10 {
				o[i] = FloatNullValue
			}
		}
		if kind == 5 || kind == 6 {
			switch {
			case i%31 == 0:
				if (i/31)%2 == 1 {
					o[i] = float32(math.Inf(1))
				} else {
					o[i] = float32(math.Inf(-1))
				}
			case i%37 == 0:
				o[i] = nan32
			case i%17 == 0:
				if (i/17)%2 == 1 {
					o[i] = 1.0e-40
				} else {
					o[i] = -3.0e-42
				}
			case i%13 == 0:
				o[i] = float32(math.Copysign(0, -1))
			case i%19 == 0:
				if (i/19)%2 == 1 {
					o[i] = math.Float32frombits(0x00800000) // FLT_MIN
				} else {
					o[i] = -math.Float32frombits(0x00800000)
				}
			}
		}
	}
	return o
}

// ---- raw byte layouts the harness hashed ----

// LEFloat32 is a float32 array as little-endian bytes.
func LEFloat32(v []float32) []byte {
	b := make([]byte, 4*len(v))
	for i, x := range v {
		binary.LittleEndian.PutUint32(b[4*i:], math.Float32bits(x))
	}
	return b
}

// LEInts encodes integer pixels little-endian, width bytes each (two's
// complement truncation, as cref.c put_le).
func LEInts(v []int64, width int) []byte {
	b := make([]byte, width*len(v))
	for i, x := range v {
		for k := 0; k < width; k++ {
			b[width*i+k] = byte(uint64(x) >> (8 * k))
		}
	}
	return b
}

// SHA256 is the hex sha256 of b.
func SHA256(b []byte) string {
	s := sha256.Sum256(b)
	return hex.EncodeToString(s[:])
}

// Check returns an error unless sha256(b) is want (the hash cref.c
// recorded for the input it gave cfitsio).
func Check(what string, b []byte, want string) error {
	if len(want) != 64 {
		return fmt.Errorf("crefgen: %s: golden has no input sha256 (%q)", what, want)
	}
	if got := SHA256(b); got != want {
		return fmt.Errorf("crefgen: %s: regenerated input sha256 %s, cfitsio's input %s: the Go port of the cref.c generator drifted", what, got[:16], want[:16])
	}
	return nil
}

// Literal gen inputs (cref.c lit_*).
var (
	litN5One = []float32{3, 1, 4, 1, 5, 9, 2, 6, 5,
		2, 7, 1, 8, 2, 8, 1, 8, 2,
		1, 7, 2, 3, 5, 5, 5, 8, 4}
	litN5Two = []float32{3, 1, 4, 1, 5, 9, 2, 6, 5, 3,
		2, 7, 1, 8, 2, 8, 1, 8, 2, 8,
		3, 8, 1, 2, 5, 5, 5, 9, 4, 6}
	litN5Zero = []float32{3, 1, 4, 1, 5, 9, 2, 6, 5,
		2, 7, 1, 8, 2, 8, 1, 8, 2,
		1, 7, 2, -1, 3, 5, 6, 8, 4}
	litN1Short = []float32{1, 4, 2, 8, 3,
		5, 1, 9, 2, 6,
		-1, 7, -1, -1, 2}
	litOneGood  = []float32{-1, -1, -1, -1, 7, -1, -1, -1, -1}
	litAbsNulls = []float32{12, 15, 11, 18, 14, 16, 13, 17, 19, 12,
		nullF, 3, 14, 16, 12, 15, 18, 11, 13, 17,
		16, 12, 15, 11, nullF, 14, 18, 13, 12, 16,
		14, 17, 12, 15, 13, 11, 16, 18, 12, 14}
)

// fillLit is cref.c's fill_lit: lit repeated over n pixels.
func fillLit(n int64, lit []float32) []float32 {
	out := make([]float32, n)
	for i := range out {
		out[i] = lit[i%len(lit)]
	}
	return out
}

// DoubleNullValue is cfitsio's DOUBLENULLVALUE (fitsio.h), the null
// sentinel fits_img_compress uses for float64 tiles.
const DoubleNullValue = -9.1191291391491e-36

// ImgDouble is gen_img_double (golden/img64).
func ImgDouble(n, nx int64, kind int, seed uint64) []float64 {
	r := NewLCG(seed)
	nan := math.Float64frombits(0x7ff8000000000000) // C NAN
	o := make([]float64, n)
	for i := int64(0); i < n; i++ {
		row, col := i/nx, i%nx
		v := float64(float64(float64(100.0+float64(0.3*float64(col)))+float64(0.2*float64(row))) + float64(5.0*r.noise4()))
		switch kind {
		case 1:
			if i%23 == 0 {
				o[i] = 0
			} else {
				o[i] = v
			}
		case 2:
			if i%31 == 0 {
				o[i] = nan
			} else {
				o[i] = v
			}
		case 3:
			if row < 10 {
				o[i] = 42.5
			} else {
				o[i] = v
			}
		case 4, 6:
			if row < 5 {
				add := 0.0
				if i%2 == 1 {
					add = 1.0e20
				}
				o[i] = add + v
			} else {
				o[i] = v
			}
		case 19:
			o[i] = float64(1.0e10 + float64(float64(i)*1.0e-6))
		case 15:
			o[i] = float64(1.0e39 + float64(1.0e33*float64(v-100.0)))
		case 14:
			d := float64(float64(float64(v-100.0)-float64(0.3*float64(col))) - float64(0.2*float64(row)))
			o[i] = float64(float64(1.0e8+float64(1.0e-4*float64(col))) + float64(0.01*d))
		default:
			o[i] = v
		}
		if kind == 7 || kind == 8 {
			if i == 3 {
				o[i] = nan
			}
			if i == 20 {
				o[i] = -1.0e30
			}
			if kind == 8 && i == 10 {
				o[i] = DoubleNullValue
			}
		}
		if kind == 16 {
			if i%23 == 0 {
				o[i] = 0
			}
			if i%31 == 0 {
				o[i] = nan
			}
		}
		if kind == 17 {
			if i%200 == 3 {
				o[i] = nan
			}
			if i%50 == 7 {
				o[i] = -999
			}
		}
		if kind == 5 || kind == 6 {
			switch {
			case i%31 == 0:
				if (i/31)%2 == 1 {
					o[i] = math.Inf(1)
				} else {
					o[i] = math.Inf(-1)
				}
			case i%37 == 0:
				o[i] = nan
			case i%17 == 0:
				if (i/17)%2 == 1 {
					o[i] = 1.0e-310
				} else {
					o[i] = -4.0e-320
				}
			case i%13 == 0:
				o[i] = math.Copysign(0, -1)
			case i%19 == 0:
				if (i/19)%2 == 1 {
					o[i] = 0x1p-1022
				} else {
					o[i] = -0x1p-1022
				}
			}
		}
	}
	return o
}

// LEFloat64 encodes float64 pixels little-endian (the bytes cref.c hashes).
func LEFloat64(v []float64) []byte {
	b := make([]byte, 8*len(v))
	for i, x := range v {
		binary.LittleEndian.PutUint64(b[8*i:], math.Float64bits(x))
	}
	return b
}

// noise4 is cref.c's lcg_noise4: a sum of 4 uniforms, centered, scaled to
// unit variance (exact IEEE + - *, no libm).
func (r *LCG) noise4() float64 {
	s := float64(float64(float64(r.Uniform()+r.Uniform())+r.Uniform()) + r.Uniform())
	return float64((s - 2.0) * 1.7320508075688772)
}
