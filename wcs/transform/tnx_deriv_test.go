package transform

import (
	"math"
	"path/filepath"
	"testing"
)

// TestTNXBasisDerivatives: the Jacobian TNX.Inverse uses comes from
// tnxBasisWithDeriv. Truth: the closed forms of T_n, P_n and x^n and their
// derivatives for n <= 5 (Abramowitz & Stegun table 22.3), on a grid over
// [-1.2, 1.2]. Budget 1e-12 absolute: the recurrences and the closed forms
// are a handful of float64 operations on values of magnitude <= ~100.
func TestTNXBasisDerivatives(t *testing.T) {
	type poly struct{ f, d func(x float64) float64 }
	pw := func(n int) poly {
		return poly{
			func(x float64) float64 { return math.Pow(x, float64(n)) },
			func(x float64) float64 {
				if n == 0 {
					return 0
				}
				return float64(n) * math.Pow(x, float64(n-1))
			}}
	}
	one := poly{func(float64) float64 { return 1 }, func(float64) float64 { return 0 }}
	lin := poly{func(x float64) float64 { return x }, func(float64) float64 { return 1 }}
	closed := map[int][]poly{
		1: { // Chebyshev T_n
			one, lin,
			{func(x float64) float64 { return 2*x*x - 1 }, func(x float64) float64 { return 4 * x }},
			{func(x float64) float64 { return 4*x*x*x - 3*x }, func(x float64) float64 { return 12*x*x - 3 }},
			{func(x float64) float64 { return 8*x*x*x*x - 8*x*x + 1 }, func(x float64) float64 { return 32*x*x*x - 16*x }},
			{func(x float64) float64 { return 16*math.Pow(x, 5) - 20*x*x*x + 5*x },
				func(x float64) float64 { return 80*x*x*x*x - 60*x*x + 5 }},
		},
		2: { // Legendre P_n
			one, lin,
			{func(x float64) float64 { return (3*x*x - 1) / 2 }, func(x float64) float64 { return 3 * x }},
			{func(x float64) float64 { return (5*x*x*x - 3*x) / 2 }, func(x float64) float64 { return (15*x*x - 3) / 2 }},
			{func(x float64) float64 { return (35*x*x*x*x - 30*x*x + 3) / 8 },
				func(x float64) float64 { return (140*x*x*x - 60*x) / 8 }},
			{func(x float64) float64 { return (63*math.Pow(x, 5) - 70*x*x*x + 15*x) / 8 },
				func(x float64) float64 { return (315*x*x*x*x - 210*x*x + 15) / 8 }},
		},
		3: {pw(0), pw(1), pw(2), pw(3), pw(4), pw(5)}, // simple polynomial
	}
	for typ, ps := range closed {
		for k := 0; k <= 24; k++ {
			x := -1.2 + 0.1*float64(k)
			b, bp := tnxBasisWithDeriv(typ, x, len(ps))
			for n, p := range ps {
				if math.Abs(b[n]-p.f(x)) > 1e-12 {
					t.Errorf("type %d n=%d x=%.2f: basis %v, closed form %v", typ, n, x, b[n], p.f(x))
				}
				if math.Abs(bp[n]-p.d(x)) > 1e-12 {
					t.Errorf("type %d n=%d x=%.2f: derivative %v, closed form %v", typ, n, x, bp[n], p.d(x))
				}
			}
		}
	}
}

// TestTNXSurfaceDerivativesRealHeaders: evalWithDeriv's partials, including
// the domain-normalisation chain factor and the cross-term layout, agree
// with a central difference of eval on the committed headers' surfaces
// (real Chebyshev half, polynomial half, polynomial full; synthetic
// Legendre full). Step h = 1e-4 of the domain width on cubic surfaces:
// truncation ~h^2 relative, rounding ~1e-16/h; budget 1e-6 relative to
// the larger partial.
func TestTNXSurfaceDerivativesRealHeaders(t *testing.T) {
	for _, name := range []string{"noao_ct4m_tnx1", "registry_tnxdoc_tnx3", "kpno_kp4m_mosaic1_ccd1", "syn_leg_full"} {
		h, _ := readCardFile(t, filepath.Join(tnxHeaderDir, name+".hdr"))
		tnx, err := ParseTNX(h)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		for si, s := range []tnxSurface{tnx.Lon, tnx.Lat} {
			hx, hy := 1e-4*(s.XMax-s.XMin), 1e-4*(s.YMax-s.YMin)
			for i := 0; i <= 4; i++ {
				for j := 0; j <= 4; j++ {
					x := s.XMin + (s.XMax-s.XMin)*float64(i)/4
					y := s.YMin + (s.YMax-s.YMin)*float64(j)/4
					_, dx, dy := s.evalWithDeriv(x, y)
					fdx := (s.eval(x+hx, y) - s.eval(x-hx, y)) / (2 * hx)
					fdy := (s.eval(x, y+hy) - s.eval(x, y-hy)) / (2 * hy)
					scale := math.Max(math.Max(math.Abs(fdx), math.Abs(fdy)), 1e-12)
					if math.Abs(dx-fdx) > 1e-6*scale || math.Abs(dy-fdy) > 1e-6*scale {
						t.Errorf("%s surface %d at (%g, %g): d/dx %v vs central difference %v, d/dy %v vs %v",
							name, si, x, y, dx, fdx, dy, fdy)
					}
				}
			}
		}
	}
}
