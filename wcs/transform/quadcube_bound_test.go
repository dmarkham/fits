package transform

import (
	"math"
	"testing"
)

// TestQuadcubeFaceBound pins the face-coordinate check that TSC and QSC
// share (quadcubeFaceBound), which no (phi, theta) reaches through the
// projections (see there), so TestWrefProjections cannot see its
// threshold. The rule is wcslib 8.10 prj.c tscs2x / qscs2x, verbatim:
//
//	if (fabs(xf) > 1.0) {
//	  if (fabs(xf) > 1.0+tol) { istat = 1; ... PRJERR_BAD_WORLD_SET }
//	  xf = copysign(1.0, xf);
//	}
//	(the same for yf), with tol = 1.0e-12.
func TestQuadcubeFaceBound(t *testing.T) {
	if quadcubeFaceTol != 1.0e-12 {
		t.Fatalf("quadcubeFaceTol = %g, wcslib tscs2x/qscs2x tol = 1.0e-12", quadcubeFaceTol)
	}
	for _, c := range []struct {
		xf, yf, wx, wy float64
		ok             bool
	}{
		{0.25, -0.5, 0.25, -0.5, true},
		{1, -1, 1, -1, true},
		{math.Nextafter(1, 2), 0, 1, 0, true},
		{1 + 9e-13, 0, 1, 0, true},
		{-1 - 9e-13, 0.5, -1, 0.5, true},
		{0.5, 1 + 9e-13, 0.5, 1, true},
		{1 + 1.1e-12, 0, 1, 0, false},
		{-1 - 1.1e-12, 0, -1, 0, false},
		{0, 1 + 1.1e-12, 0, 1, false},
		{0, -1 - 1e-6, 0, -1, false},
		{1 + 1e-3, 1 + 1e-13, 1, 1, false},
	} {
		x, y, ok := quadcubeFaceBound(c.xf, c.yf)
		if ok != c.ok || (ok && (x != c.wx || y != c.wy)) {
			t.Errorf("quadcubeFaceBound(%.17g, %.17g) = (%.17g, %.17g, %v), wcslib rule gives (%g, %g, %v)",
				c.xf, c.yf, x, y, ok, c.wx, c.wy, c.ok)
		}
	}
}
