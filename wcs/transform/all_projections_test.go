package transform

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// TestAllProjectionsRoundTrip is a SELF-CONSISTENCY check, not a truth
// test: for every implemented projection it builds a Transform from a
// hand-made header and requires pixel -> sky -> pixel to come back to the
// starting pixel, near the reference pixel (offsets of 3-5 px at
// 0.01 deg/px). Agreement with wcslib is proven elsewhere:
// TestWrefProjections (every code, default and non-default PV, both
// directions, grids past every boundary) and TestHdrGoldens (full header
// pipeline, one case per projection).
//
// Tolerances are the round-trip floor of each method: 1e-9 px for
// closed-form directions, 1e-7..1e-6 px where one direction iterates (MOL,
// PCO, ZPN, SZP), and 2e-2 px for CSC, whose forward and inverse
// polynomials are independent fits that do not invert each other exactly
// (wcslib's own CSC drifts the same way).
func TestAllProjectionsRoundTrip(t *testing.T) {
	type pv = map[wcs.PVKey]float64
	pvKey := func(axis, index int, val float64) pv {
		return pv{{Axis: axis, Index: index}: val}
	}
	conicPV := pvKey(2, 1, 45.0)
	zpnPV := pvKey(2, 1, 1.0)
	bonPV := pvKey(2, 1, 30.0)

	// offsets: small in-domain offsets from CRPIX-1. 3-5 pixels at
	// 0.01°/pixel = 0.03-0.05°, well inside every projection's valid range.
	smallOffsets := [][2]float64{{0, 0}, {3, 4}, {-2, 5}, {4, -3}, {-5, -4}}

	cases := []struct {
		code      string
		crval1    float64
		crval2    float64
		pv        pv
		tolerance float64
	}{
		// --- Zenithal (§5.1) ---
		{code: "AZP", crval1: 0, crval2: 0, pv: pvKey(2, 1, 2), tolerance: 1e-7},
		{code: "SZP", crval1: 0, crval2: 0, pv: pvKey(2, 1, 2), tolerance: 1e-6},
		{code: "TAN", crval1: 45, crval2: 30, tolerance: 1e-9},
		{code: "STG", crval1: 90, crval2: 45, tolerance: 1e-9},
		{code: "SIN", crval1: 0, crval2: 0, tolerance: 1e-9},
		{code: "ARC", crval1: 200, crval2: 10, tolerance: 1e-9},
		{code: "ZPN", crval1: 0, crval2: 0, pv: zpnPV, tolerance: 1e-7},
		{code: "ZEA", crval1: 60, crval2: -30, tolerance: 1e-9},
		{code: "AIR", crval1: 30, crval2: 20, tolerance: 1e-7},
		{code: "AIR", crval1: 30, crval2: 20, pv: pvKey(2, 1, 45), tolerance: 1e-7},
		// --- Cylindrical (§5.2) ---
		{code: "CYP", crval1: 0, crval2: 0, pv: pvKey(2, 1, 1), tolerance: 1e-8},
		{code: "CEA", crval1: 0, crval2: 0, pv: pvKey(2, 1, 1), tolerance: 1e-9},
		{code: "CAR", crval1: 0, crval2: 0, tolerance: 1e-9},
		{code: "MER", crval1: 0, crval2: 0, tolerance: 1e-9},
		// --- Pseudo-cylindrical (§5.3) ---
		{code: "SFL", crval1: 0, crval2: 0, tolerance: 1e-9},
		{code: "PAR", crval1: 0, crval2: 0, tolerance: 1e-9},
		{code: "MOL", crval1: 0, crval2: 0, tolerance: 1e-8},
		{code: "AIT", crval1: 0, crval2: 0, tolerance: 1e-9},
		// --- Conic (§5.4) ---
		{code: "COP", crval1: 0, crval2: 45, pv: conicPV, tolerance: 1e-8},
		{code: "COE", crval1: 0, crval2: 45, pv: conicPV, tolerance: 1e-8},
		{code: "COD", crval1: 0, crval2: 45, pv: conicPV, tolerance: 1e-8},
		{code: "COO", crval1: 0, crval2: 45, pv: conicPV, tolerance: 1e-8},
		// --- Polyconic (§5.5) ---
		{code: "BON", crval1: 0, crval2: 30, pv: bonPV, tolerance: 1e-8},
		{code: "PCO", crval1: 0, crval2: 0, tolerance: 1e-6},
		// --- Quad-cube (§5.6) ---
		{code: "TSC", crval1: 0, crval2: 0, tolerance: 1e-9},
		// CSC's forward/inverse polynomials are independent fits, not
		// exact inverses of each other.
		{code: "CSC", crval1: 0, crval2: 0, tolerance: 2e-2},
		{code: "QSC", crval1: 0, crval2: 0, tolerance: 1e-9},
		// --- HEALPix ---
		{code: "HPX", crval1: 0, crval2: 0, tolerance: 1e-9},
		{code: "XPH", crval1: 0, crval2: 90, tolerance: 1e-9},
	}

	const scale = 0.01
	const crpix1, crpix2 = 200.0, 200.0

	for _, c := range cases {
		t.Run(c.code, func(t *testing.T) {
			w := &wcs.Header{
				NAxis:    2,
				CType:    []string{"RA---" + c.code, "DEC--" + c.code},
				CRVal:    []float64{c.crval1, c.crval2},
				CRPix:    []float64{crpix1, crpix2},
				CDelt:    []float64{-scale, scale},
				PC:       [][]float64{{1, 0}, {0, 1}},
				PV:       c.pv,
				LonAxis:  1,
				LatAxis:  2,
				AxisType: []string{"RA", "DEC"},
				ProjCode: []string{c.code, c.code},
			}
			tr, err := New(w)
			if err != nil {
				t.Fatalf("New: %v", err)
			}
			// Use CRPIX - 1 + offset to place test pixels relative to the
			// 0-based reference.
			base0 := crpix1 - 1
			base1 := crpix2 - 1
			for _, off := range smallOffsets {
				px := base0 + off[0]
				py := base1 + off[1]
				a, d, err := tr.PixelToSky(px, py)
				if err != nil {
					t.Errorf("pixel (%v,%v): PixelToSky: %v", px, py, err)
					continue
				}
				p1, p2, err := tr.SkyToPixel(a, d)
				if err != nil {
					t.Errorf("pixel (%v,%v) → sky (%v,%v): SkyToPixel: %v", px, py, a, d, err)
					continue
				}
				if math.Abs(p1-px) > c.tolerance || math.Abs(p2-py) > c.tolerance {
					t.Errorf("pixel (%v,%v) → sky (%v,%v) → pixel (%v,%v); err (%v,%v), tol %v",
						px, py, a, d, p1, p2, p1-px, p2-py, c.tolerance)
				}
			}
		})
	}
}
