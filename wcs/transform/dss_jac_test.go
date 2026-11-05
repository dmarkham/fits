package transform

import (
	"math"
	"testing"
)

// TestDSSJacobian: the DSS inverse is Newton iteration on the plate
// polynomial, and its Jacobian (dssPoly.jac) only steers the iteration: a
// wrong partial still converges to the same root, so no golden can catch
// it (a mutant that drops the dv*a[4] term of d/dm1 survived every header
// case). The partials are pinned here against central differences of the
// forward polynomial dssPoly.eval, whose values agree with wcslib's dssset
// / tpd5 on every DSS golden (TestHdrGoldens), on a grid out to 6000 px.
// Step h = 1e-2 px; budget 1e-6 + 1e-9 |grad f|: the round-off of the
// difference is ~eps |f| / h (<= 1e-6 for |f| <= 1e8) and the truncation
// error (h^2 / 6 times the third derivative) stays below 1e-9 |grad f| on
// these plates (measured worst 4.5e-5 at |grad f| = 7.2e5 on the strongly
// non-linear synthetic degree-5 plate, 2.6e-10 on the real one). Dropping
// the dv*a[4] term is off by 2.35 at (-6000, -6000) on the real plate.
func TestDSSJacobian(t *testing.T) {
	for _, name := range []string{"real_dss_plate_solution", "dss_all_terms_deg5"} {
		g := loadHdrGolden(t, name)
		tr, err := hdrTransform(g)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		if tr.dss == nil {
			t.Fatalf("%s: no DSS polynomial", name)
		}
		p := tr.dss
		const h = 1e-2
		worst := 0.0
		n := 0
		for _, m1 := range []float64{-6000, -2500, -300, 0, 17.25, 450, 3100, 6000} {
			for _, m2 := range []float64{-6000, -1200, -40, 0, 9.5, 800, 2900, 6000} {
				for j := 0; j < 2; j++ {
					val, d1, d2 := p.jac(j, m1, m2)
					if val != p.eval(j, m1, m2) {
						t.Errorf("%s axis %d at (%g, %g): jac value %v != eval %v", name, j, m1, m2, val, p.eval(j, m1, m2))
					}
					n1 := (p.eval(j, m1+h, m2) - p.eval(j, m1-h, m2)) / (2 * h)
					n2 := (p.eval(j, m1, m2+h) - p.eval(j, m1, m2-h)) / (2 * h)
					e := math.Max(math.Abs(d1-n1), math.Abs(d2-n2))
					tol := 1e-6 + 1e-9*math.Hypot(d1, d2)
					worst = math.Max(worst, e/tol)
					n++
					if !(e <= tol) {
						t.Errorf("%s axis %d at (%g, %g): partials (%.12g, %.12g), central differences (%.12g, %.12g)", name, j, m1, m2, d1, d2, n1, n2)
					}
				}
			}
		}
		t.Logf("%s: %d Jacobians, worst |analytic - central difference| = %.3g of the budget", name, n, worst)
	}
}
