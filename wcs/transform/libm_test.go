package transform

import (
	"math"
	"testing"
)

// libmGolden is the mpmath truth (200 bits, rounded once) for sin, cos,
// asin and acos at the arguments the projections use, written by
// testdata/wref/gen_libm.py.
const libmGolden = "testdata/wref/libm/trig.json"

type libmTruth struct {
	Generator        string
	Mpmath           string
	Prec             int
	SinCos, AsinAcos [][3]float64
}

func loadLibm(t *testing.T) libmTruth {
	t.Helper()
	var meta struct {
		Generator string `json:"generator"`
		Mpmath    string `json:"mpmath"`
		Prec      int    `json:"precision_bits"`
		SinCos    int    `json:"sincos"`
		AsinAcos  int    `json:"asinacos"`
	}
	vals := loadF64Pack(t, libmGolden, &meta)
	g := libmTruth{Generator: meta.Generator, Mpmath: meta.Mpmath, Prec: meta.Prec,
		SinCos: f64packRows3(t, &vals, meta.SinCos), AsinAcos: f64packRows3(t, &vals, meta.AsinAcos)}
	if len(vals) != 0 {
		t.Fatalf("%s: %d doubles left over after the rows the metadata lists", libmGolden, len(vals))
	}
	return g
}

// ulpErr returns |got - want| in units of the last place of want (an
// exact zero want counts any non-zero got as +Inf).
func ulpErr(got, want float64) float64 {
	if got == want {
		return 0
	}
	if want == 0 || got != got {
		return math.Inf(1)
	}
	aw := math.Abs(want)
	return math.Abs(got-want) / (math.Nextafter(aw, math.Inf(1)) - aw)
}

// TestLibmAccuracy bounds the libm calls behind wcstrig.go against mpmath
// (200 bits, rounded once) at every committed argument: Go's math.Sin and
// math.Cos within 2 ULP (2^-100 absolute next to their zeros) and
// libAsin/libAcos within 1 ULP (the measured maxima are logged). These
// are the premises of the value budgets against wcslib (which calls
// glibc, itself within 1 ULP): a few ULP in
// a sine or cosine stays a few ULP in a well-conditioned projection. The
// test also logs how far math.Asin/math.Acos are off on the same
// arguments; they are not used because next to +-1 they lose up to 2.5e7
// ULP (see libAsin).
func TestLibmAccuracy(t *testing.T) {
	g := loadLibm(t)
	if g.Prec < 200 || len(g.SinCos) < 2000 || len(g.AsinAcos) < 900 {
		t.Fatalf("golden too small: %d bits, %d sin/cos, %d asin/acos", g.Prec, len(g.SinCos), len(g.AsinAcos))
	}
	worst := map[string]float64{}
	check := func(name string, budget, arg, got, want float64) {
		e := ulpErr(got, want)
		if e > worst[name] {
			worst[name] = e
		}
		if e > budget {
			t.Errorf("%s(%v) = %v, exact %v: %.3g ULP > %g", name, arg, got, want, e, budget)
		}
	}
	// sin and cos: ULP of the result, except that near a zero of the
	// function (an argument next to a multiple of pi/2) the error is
	// measured against 2^-100 absolute: Go reduces the argument with a
	// three-part pi/4, which leaves an absolute error near 1e-31 there
	// (7 ULP of cos(fl(pi/2)) = 6.1e-17); wcstrig.go returns exact zeros
	// at multiples of 90 degrees, and every other caller adds such a value
	// to, or compares it with, quantities of order 1.
	trig := func(name string, arg, got, want float64) {
		e := ulpErr(got, want)
		if a := math.Abs(got-want) / 0x1p-100; a < e {
			e = a
		}
		if e > worst[name] {
			worst[name] = e
		}
		if e > 2 {
			t.Errorf("%s(%v) = %v, exact %v: %.3g ULP (or units of 2^-100) > 2", name, arg, got, want, e)
		}
	}
	for _, r := range g.SinCos {
		s, c := math.Sincos(r[0])
		trig("math.Sincos sin", r[0], s, r[1])
		trig("math.Sincos cos", r[0], c, r[2])
		trig("math.Sin", r[0], math.Sin(r[0]), r[1])
		trig("math.Cos", r[0], math.Cos(r[0]), r[2])
	}
	for _, r := range g.AsinAcos {
		check("libAsin", 1, r[0], libAsin(r[0]), r[1])
		check("libAcos", 1, r[0], libAcos(r[0]), r[2])
		check("math.Asin (unused)", math.Inf(1), r[0], math.Asin(r[0]), r[1])
		check("math.Acos (unused)", math.Inf(1), r[0], math.Acos(r[0]), r[2])
	}
	for _, k := range []string{"math.Sincos sin", "math.Sincos cos", "math.Sin", "math.Cos", "libAsin", "libAcos", "math.Asin (unused)", "math.Acos (unused)"} {
		t.Logf("%-20s worst %.3g ULP", k, worst[k])
	}
	// Outside [-1, 1] asin/acos are NaN, as in C.
	for _, x := range []float64{1 + 0x1p-52, -1 - 0x1p-52, 2, math.Inf(1)} {
		if a, c := libAsin(x), libAcos(x); !math.IsNaN(a) || !math.IsNaN(c) {
			t.Errorf("libAsin/libAcos(%v) = %v, %v, want NaN", x, a, c)
		}
	}
	// Exact at the ends and at 0 (the projections rely on these).
	for _, c := range [][2]float64{{libAsin(1), math.Pi / 2}, {libAsin(-1), -math.Pi / 2}, {libAsin(0), 0},
		{libAcos(1), 0}, {libAcos(-1), math.Pi}, {libAcos(0), math.Pi / 2}} {
		if c[0] != c[1] {
			t.Errorf("libAsin/libAcos end value %v, want %v", c[0], c[1])
		}
	}
}

// libmAtan2Golden is the mpmath truth for atan2 at the arguments the
// projections pass it (behind atand and atan2d, and qscs2x's direct
// atan2), written by testdata/wref/gen_libm.py.
const libmAtan2Golden = "testdata/wref/libm/atan2.json"

// TestAtan2Accuracy bounds Go's math.Atan2 and math.Atan (atand) within 1
// ULP of mpmath (200 bits, rounded once) at every committed argument, the
// premise of the value budgets against wcslib, whose glibc atan2 is
// within 1 ULP too. The measured maxima are logged.
func TestAtan2Accuracy(t *testing.T) {
	var meta struct {
		Prec  int `json:"precision_bits"`
		Atan2 int `json:"atan2"`
	}
	vals := loadF64Pack(t, libmAtan2Golden, &meta)
	rows := f64packRows3(t, &vals, meta.Atan2)
	if len(vals) != 0 {
		t.Fatalf("%s: %d doubles left over", libmAtan2Golden, len(vals))
	}
	if meta.Prec < 200 || len(rows) < 2500 {
		t.Fatalf("golden too small: %d bits, %d rows", meta.Prec, len(rows))
	}
	worst, worstAtan := 0.0, 0.0
	for _, r := range rows {
		got := math.Atan2(r[0], r[1])
		e := ulpErr(got, r[2])
		worst = math.Max(worst, e)
		if e > 1 {
			t.Errorf("math.Atan2(%v, %v) = %v, exact %v: %.3g ULP", r[0], r[1], got, r[2], e)
		}
		if r[1] == 1 {
			got := math.Atan(r[0])
			e := ulpErr(got, r[2])
			worstAtan = math.Max(worstAtan, e)
			if e > 1 {
				t.Errorf("math.Atan(%v) = %v, exact %v: %.3g ULP", r[0], got, r[2], e)
			}
		}
	}
	t.Logf("%d atan2 arguments: math.Atan2 worst %.3g ULP, math.Atan (x = 1 rows) worst %.3g ULP", len(rows), worst, worstAtan)
}
