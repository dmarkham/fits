package stretch

import (
	"errors"
	"fmt"
	"math"
)

// GHTParams holds the parameters of Siril's Generalized Hyperbolic
// Stretch (`ght -D= -B= -LP= -SP= -HP=`), normal (not inverse) mode.
type GHTParams struct {
	D  float64 // Stretch intensity, user scale, 0 <= D <= 10; Siril uses expm1(D) internally
	B  float64 // Local intensity, -5 <= B <= 15: 0 = exponential, -1 = logarithmic, >0 or <0 = power family. |B| < float32(0.001) = 0.00100000005 counts as 0, compared in float64 as siril-cli does (core/command.c:2328), so B = ±0.001 itself counts as 0 (Siril's GUI compares in float32 and keeps ±0.001)
	LP float64 // Low protection point, 0 <= LP <= SP
	SP float64 // Symmetry point, 0 <= SP <= 1
	HP float64 // High protection point, SP <= HP <= 1; 0 means Siril's default 1 unless HPZero is set
	// HPZero makes HP = 0 mean 0, as siril-cli's -HP=0 does (valid only
	// with LP = SP = 0; GHT is then the clamp alone). HP must be 0 when
	// HPZero is set. The zero value of GHTParams leaves HP at Siril's
	// default 1, as a ght command without -HP= does.
	HPZero bool
	BP     float64 // Unused by GHT (black point of Siril's linstretch; see Linear)
}

// hp is the high protection point GHT uses.
func (p GHTParams) hp() float64 {
	if p.HP == 0 && !p.HPZero {
		return 1
	}
	return p.HP
}

// ghtCoeffs mirrors the float fields of Siril's ght_compute_params used by
// STRETCH_PAYNE_NORMAL (filters/ght.h).
type ghtCoeffs struct {
	b1             float32
	a2, b2, c2, d2 float32
	e2             float32
	a3, b3, c3, d3 float32
	e3             float32
	a4, b4         float32
}

// GHT applies the Generalized Hyperbolic Stretch to one channel exactly as
// Siril's ght command does for a 32-bit image with the default
// independent-channel colour model (apply_linked_ght_to_fbuf_indep):
// each pixel is clamped to [0, 1], 0 stays 0, and the rest go through
// Siril's piecewise GHT (filters/ght.c GHTsetup/GHT, STRETCH_PAYNE_NORMAL)
// with the parameters converted to float32 as process_ght_args does
// (D -> (float)expm1(D)). Output is not clamped (Siril does not clamp it).
//
// Parameters outside the range Siril accepts are an error wrapping
// ErrGHTParams, and the pixels are left unchanged: siril-cli's
// process_ght_args (core/command.c:2363-2387) refuses D outside (0, 10],
// B outside [0, 15], SP outside [0, 1], LP outside [0, SP] and HP outside
// [SP, 1]. Two GUI extensions are accepted because Siril's GUI runs the
// same code on them: B in [-5, 0) (gui/uifiles/histogram_dialog.ui:77)
// and D = 0, the identity after the clamp. NaN or infinite parameters are
// refused.
func GHT(pixels []float32, params GHTParams) error {
	if err := params.validate(); err != nil {
		return err
	}
	g := prepareGHT(params)
	for i, v := range pixels {
		x := fminf(1, fmaxf(0, v))
		if x == 0 {
			pixels[i] = 0
			continue
		}
		pixels[i] = g.eval(x, &libm32)
	}
	return nil
}

// ErrGHTParams is wrapped by the error GHT returns for parameters outside
// the accepted range.
var ErrGHTParams = errors.New("stretch: GHT parameters out of range")

func (p GHTParams) validate() error {
	HP := p.hp()
	in := func(v, lo, hi float64) bool { return v >= lo && v <= hi } // false for NaN
	switch {
	case !in(p.D, 0, 10):
		return fmt.Errorf("%w: D = %v, want 0 <= D <= 10", ErrGHTParams, p.D)
	case !in(p.B, -5, 15):
		return fmt.Errorf("%w: B = %v, want -5 <= B <= 15", ErrGHTParams, p.B)
	case !in(p.SP, 0, 1):
		return fmt.Errorf("%w: SP = %v, want 0 <= SP <= 1", ErrGHTParams, p.SP)
	case !in(p.LP, 0, p.SP):
		return fmt.Errorf("%w: LP = %v, want 0 <= LP <= SP = %v", ErrGHTParams, p.LP, p.SP)
	case p.HPZero && p.HP != 0:
		return fmt.Errorf("%w: HPZero set with HP = %v", ErrGHTParams, p.HP)
	case !in(HP, p.SP, 1):
		return fmt.Errorf("%w: HP = %v, want SP = %v <= HP <= 1", ErrGHTParams, HP, p.SP)
	}
	return nil
}

// ghtPrepared is GHTParams converted the way process_ght_args does
// (command.c:2235: |B| < 1e-3f -> 0 compared in double, D -> (float)expm1(D),
// everything to float32) plus the GHTsetup coefficients.
type ghtPrepared struct {
	b, d, lp, sp, hp float32
	c                ghtCoeffs
}

func prepareGHT(params GHTParams) ghtPrepared {
	B := params.B
	if math.Abs(B) < float64(float32(1e-3)) {
		B = 0
	}
	HP := params.hp()
	g := ghtPrepared{b: float32(B), d: float32(math.Expm1(params.D)),
		lp: float32(params.LP), sp: float32(params.SP), hp: float32(HP)}
	g.c = ghtSetup(g.b, g.d, g.lp, g.sp, g.hp)
	return g
}

func (g *ghtPrepared) eval(x float32, m *float32Libm) float32 {
	return ghtEval(x, g.b, g.d, g.lp, g.sp, g.hp, &g.c, m)
}

// float32 libm used by ghtEval. Siril calls glibc's logf/powf/expf, which
// are not correctly rounded (glibc 2.44 on x86-64 disagrees with the
// correctly rounded value on about 0.06% of powf/expf and 0.9% of logf
// arguments, testdata/sref/libm_probe.c); Go has no float32 versions, so
// these round the float64 results. A GHT pixel can therefore differ from
// Siril's where glibc's result is one ULP away; the tests accept that only
// at the pixels each golden lists as libm-sensitive (ght_test.go ghtLibm).
type float32Libm struct {
	logf func(float32) float32
	powf func(float32, float32) float32
	expf func(float32) float32
}

var libm32 = float32Libm{logf: logf, powf: powf, expf: expf}

// ghtSetup is Siril's GHTsetup for STRETCH_PAYNE_NORMAL (filters/ght.c:34).
// All inputs and stored coefficients are float32; log1p/pow/exp are the C
// double functions, so a float32 expression that meets one of them is
// promoted to float64 for the rest of that C expression, exactly as C's
// usual arithmetic conversions do, and rounded to float32 on assignment.
// Products that feed an addition are wrapped in an explicit conversion,
// which the Go spec guarantees is never fused into an FMA (Siril's x86-64
// build has no FMA either).
func ghtSetup(B, D, LP, SP, HP float32) ghtCoeffs {
	var c ghtCoeffs
	if D == 0 {
		return c
	}
	f64 := func(x float32) float64 { return float64(x) }
	switch {
	case B == -1:
		qlp := float32(-1.0 * math.Log1p(f64(float32(D*(SP-LP)))))
		q0 := qlp - D*LP/(1+float32(D*(SP-LP)))
		qwp := float32(math.Log1p(f64(float32(D * (HP - SP)))))
		q1 := qwp + float32(D*(1-HP))/(1+float32(D*(HP-SP)))
		q := 1 / (q1 - q0)
		c.b1 = (1 + float32(D*(SP-LP))) / (D * q)
		c.a2 = float32((-q0) * q)
		c.b2 = -q
		c.c2 = 1 + float32(D*SP)
		c.d2 = -D
		c.a3 = float32((-q0) * q)
		c.b3 = q
		c.c3 = 1 - float32(D*SP)
		c.d3 = D
		c.a4 = float32((qwp - q0 - float32(D*HP)/(1+float32(D*(HP-SP)))) * q)
		c.b4 = float32(q*D) / (1 + float32(D*(HP-SP)))
	case B < 0:
		B = -B
		base1 := f64(1 + float32(float32(D*B)*(SP-LP)))
		base2 := f64(1 + float32(float32(D*B)*(HP-SP)))
		qlp := float32((1.0 - math.Pow(base1, f64((B-1)/B))) / f64(B-1))
		q0 := float32(f64(qlp) - float64(f64(float32(D*LP))*math.Pow(base1, f64(-1/B))))
		qwp := float32((math.Pow(base2, f64((B-1)/B)) - 1.0) / f64(B-1))
		q1 := float32(f64(qwp) + float64(f64(float32(D*(1-HP)))*math.Pow(base2, f64(-1/B))))
		q := 1 / (q1 - q0)
		c.b1 = float32(f64(D) * math.Pow(base1, f64(-1/B)) * f64(q))
		c.a2 = float32((1/(B-1) - q0) * q)
		c.b2 = -q / (B - 1)
		c.c2 = 1 + float32(float32(D*B)*SP)
		c.d2 = float32((-D) * B)
		c.e2 = (B - 1) / B
		c.a3 = float32((-1/(B-1) - q0) * q)
		c.b3 = q / (B - 1)
		c.c3 = 1 - float32(float32(D*B)*SP)
		c.d3 = float32(D * B)
		c.e3 = (B - 1) / B
		c.a4 = float32((f64(qwp-q0) - float64(f64(float32(D*HP))*math.Pow(base2, f64(-1/B)))) * f64(q))
		c.b4 = float32(f64(D) * math.Pow(base2, f64(-1/B)) * f64(q))
	case B == 0:
		e1 := math.Exp(f64(float32((-D) * (SP - LP))))
		e2 := math.Exp(f64(float32((-D) * (HP - SP))))
		qlp := float32(e1)
		q0 := float32(f64(qlp) - float64(f64(float32(D*LP))*e1))
		qwp := float32(2.0 - e2)
		q1 := float32(f64(qwp) + float64(f64(float32(D*(1-HP)))*e2))
		q := 1 / (q1 - q0)
		c.b1 = float32(f64(D) * e1 * f64(q))
		c.a2 = float32(-q0 * q)
		c.b2 = q
		c.c2 = float32((-D) * SP)
		c.d2 = D
		c.a3 = float32((2 - q0) * q)
		c.b3 = -q
		c.c3 = float32(D * SP)
		c.d3 = -D
		c.a4 = float32((f64(qwp-q0) - float64(f64(float32(D*HP))*e2)) * f64(q))
		c.b4 = float32(f64(D) * e2 * f64(q))
	default: // B > 0
		base1 := f64(1 + float32(float32(D*B)*(SP-LP)))
		base2 := f64(1 + float32(float32(D*B)*(HP-SP)))
		e := f64(-(1 + B) / B)
		qlp := float32(math.Pow(base1, f64(-1/B)))
		q0 := float32(f64(qlp) - float64(f64(float32(D*LP))*math.Pow(base1, e)))
		qwp := float32(2.0 - math.Pow(base2, f64(-1/B)))
		q1 := float32(f64(qwp) + float64(f64(float32(D*(1-HP)))*math.Pow(base2, e)))
		q := 1 / (q1 - q0)
		c.b1 = float32(f64(D) * math.Pow(base1, e) * f64(q))
		c.a2 = float32(-q0 * q)
		c.b2 = q
		c.c2 = 1 + float32(float32(D*B)*SP)
		c.d2 = float32((-D) * B)
		c.e2 = -1 / B
		c.a3 = float32((2 - q0) * q)
		c.b3 = -q
		c.c3 = 1 - float32(float32(D*B)*SP)
		c.d3 = float32(D * B)
		c.e3 = -1 / B
		eB := f64(-(B + 1) / B)
		c.a4 = float32((f64(qwp-q0) - float64(f64(float32(D*HP))*math.Pow(base2, eB))) * f64(q))
		c.b4 = float32((f64(D) * math.Pow(base2, eB)) * f64(q))
	}
	return c
}

// ghtEval is Siril's float GHT() for STRETCH_PAYNE_NORMAL with BP = 0
// (filters/ght.c), using float32 logf/powf/expf.
func ghtEval(in, B, D, LP, SP, HP float32, c *ghtCoeffs, m *float32Libm) float32 {
	in = fmaxf(0, in) // (in - 0) / (1 - 0)
	if D == 0 {
		return in
	}
	var res1, res2 float32
	switch {
	case B == -1:
		res1 = c.a2 + float32(c.b2*m.logf(c.c2+float32(c.d2*in)))
		res2 = c.a3 + float32(c.b3*m.logf(c.c3+float32(c.d3*in)))
	case B < 0 || B > 0:
		res1 = c.a2 + float32(c.b2*m.powf(c.c2+float32(c.d2*in), c.e2))
		res2 = c.a3 + float32(c.b3*m.powf(c.c3+float32(c.d3*in), c.e3))
	default:
		res1 = c.a2 + float32(c.b2*m.expf(c.c2+float32(c.d2*in)))
		res2 = c.a3 + float32(c.b3*m.expf(c.c3+float32(c.d3*in)))
	}
	switch {
	case in < LP:
		return float32(c.b1 * in)
	case in < SP:
		return res1
	case in < HP:
		return res2
	default:
		return c.a4 + float32(c.b4*in)
	}
}
