package transform

import (
	"encoding/json"
	"math"
	"os"
	"testing"
)

// wcstrigGuardGolden is wcslib 8.10's own asind/acosd at arguments just
// inside, on and outside +-(1 + WCSTRIG_TOL), written by `wref trig`
// (testdata/wref/wref.c write_trig).
const wcstrigGuardGolden = "testdata/wref/libm/wcstrig_guards.json"

// TestWcstrigGuards pins the out-of-range tolerance of asind and acosd
// (wcstrigTol, wcslib's WCSTRIG_TOL = 1e-10), branch by branch, against
// wcslib: an argument up to 1e-10 past +-1 is treated as +-1 (an exact
// +-90, 0 or 180), one further out is a NaN. No projection reaches these
// arguments (the verifier toggled the tolerance at run time over all 61
// wref cases, inverse [-400, 400]^2 step 0.5 and the forward grids, and
// saw no output change; every caller's argument is within a few ULP of
// [-1, 1] or far outside it), so the tests through the projections could
// not tell 1e-10 from 1e-3: these rows can. Guard outputs must be
// bit-identical; inside (-1, 1) the result is asin/acos * R2D, where
// libAsin/libAcos may differ from glibc's asin/acos by 1 ULP, so 2 ULP of
// the result are allowed there.
func TestWcstrigGuards(t *testing.T) {
	data, err := os.ReadFile(wcstrigGuardGolden)
	if err != nil {
		t.Fatalf("golden missing: %v (run make -C testdata/wref wref-goldens)", err)
	}
	var g struct {
		WcslibVersion string        `json:"wcslib_version"`
		Tol           float64       `json:"wcstrig_tol"`
		Rows          [][3]*float64 `json:"rows"`
	}
	if err := json.Unmarshal(data, &g); err != nil {
		t.Fatal(err)
	}
	if g.WcslibVersion != wrefWcslibVersion || g.Tol != wcstrigTol || len(g.Rows) < 50 {
		t.Fatalf("golden: wcslib %q, WCSTRIG_TOL %g (ours %g), %d rows", g.WcslibVersion, g.Tol, wcstrigTol, len(g.Rows))
	}
	var nguard, nnan int
	for _, r := range g.Rows {
		v := *r[0]
		for i, f := range []struct {
			name string
			fn   func(float64) float64
		}{{"asind", asind}, {"acosd", acosd}} {
			got, want := f.fn(v), r[1+i]
			switch {
			case want == nil:
				nnan++
				if !math.IsNaN(got) {
					t.Errorf("%s(%.17g) = %.17g, wcslib NaN (argument beyond 1 + WCSTRIG_TOL)", f.name, v, got)
				}
			case math.Abs(v) >= 1:
				nguard++
				if math.Float64bits(got) != math.Float64bits(*want) {
					t.Errorf("%s(%.17g) = %.17g, wcslib %.17g (guard branch)", f.name, v, got, *want)
				}
			default:
				w := *want
				// 1 ULP of asin/acos (libAsin/libAcos vs glibc) becomes up
				// to 2 ULP of the degrees after the rounded * R2D:
				// acosd(0.5) is 59.999999999999993 here, 60.000000000000007
				// in wcslib (fl(pi/3) is 1 ULP apart in the two libms).
				ulp := math.Nextafter(math.Abs(w), math.Inf(1)) - math.Abs(w)
				if math.Abs(got-w) > 2*ulp {
					t.Errorf("%s(%.17g) = %.17g, wcslib %.17g", f.name, v, got, w)
				}
			}
		}
	}
	if nguard < 20 || nnan < 20 {
		t.Errorf("golden exercises %d guard results and %d NaN results; want both branches at both ends", nguard, nnan)
	}
}
