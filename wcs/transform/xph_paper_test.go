package transform

import (
	"math"
	"testing"
)

// xphPaperGolden is XPH evaluated from the equations of Calabretta & Lowe
// 2013 (PASA 30, e050) in mpmath, by testdata/wref/gen_xph_paper.py: a
// truth for XPH that does not come from wcslib, so a fault shared by
// wcslib's xphs2x/xphx2s and our port of them would show here.
const xphPaperGolden = "testdata/wref/paper/xph_paper.json"

type xphPaperTable struct {
	Forward, Inverse [][4]float64
}

func loadXPHPaper(t *testing.T) xphPaperTable {
	t.Helper()
	var meta struct {
		Source  string `json:"source"`
		Forward int    `json:"forward"`
		Inverse int    `json:"inverse"`
	}
	vals := loadF64Pack(t, xphPaperGolden, &meta)
	g := xphPaperTable{Forward: f64packRows4(t, &vals, meta.Forward), Inverse: f64packRows4(t, &vals, meta.Inverse)}
	if len(vals) != 0 {
		t.Fatalf("%s: %d doubles left over after the rows the metadata lists", xphPaperGolden, len(vals))
	}
	return g
}

// TestXPHPaperEquations holds the XPH port to the paper's projection
// equations (2)-(9) and deprojection equations (10)-(17), both directions,
// including the near-pole points where the paper's own small-angle forms
// (18)/(19) take over, and the zone boundary theta = asin(2/3) +- 1e-7 deg.
// Budget: wrefTolDeg (1e-11 deg) with the comparison of the wcslib goldens
// plus the conditioning allowance of an exact truth (wrefDiffConditioned).
// The paper's mpmath values are exact to the double
// they are rounded to; the port's (18)/(19) branch differs from the exact
// formula by under 2e-12 deg (truncation of the small-angle series at
// sigma < 1e-4). Where the port (like wcslib) evaluates (9),
// sigma = sqrt(3 (1 - |sin theta|)), in double next to a pole, the
// rounding of sin theta alone moves the result by up to
// xphSigmaRoundingDeg, which is allowed on top; see there.
func TestXPHPaperEquations(t *testing.T) {
	g := loadXPHPaper(t)
	if len(g.Forward) < 1200 || len(g.Inverse) < 1200 {
		t.Fatalf("golden too small: %d forward, %d inverse", len(g.Forward), len(g.Inverse))
	}
	proj, err := Select("XPH", nil, 2)
	if err != nil {
		t.Fatal(err)
	}
	d := proj.(degAdapter)
	for _, inverse := range []bool{false, true} {
		rows := g.Forward
		if inverse {
			rows = g.Inverse
		}
		var nbad int
		var worst float64
		for _, r := range rows {
			var o0, o1 float64
			var ok bool
			if inverse {
				o0, o1, ok = d.inverseDeg(r[0], r[1])
			} else {
				o0, o1, ok = d.forwardDeg(r[0], r[1])
			}
			if !ok {
				nbad++
				if nbad <= 10 {
					t.Errorf("inverse=%v (%v, %v): rejected, paper gives (%.15g, %.15g)", inverse, r[0], r[1], r[2], r[3])
				}
				continue
			}
			diff := wrefDiffConditioned(inverse, r[2], r[3], o0, o1)
			if !inverse {
				diff = math.Max(0, diff-xphSigmaRoundingDeg(r[1]))
			}
			worst = math.Max(worst, diff)
			if !(diff <= wrefTolDeg) {
				nbad++
				if nbad <= 10 {
					t.Errorf("inverse=%v (%.15g, %.15g): got (%.15g, %.15g), paper (%.15g, %.15g), diff %.3g deg",
						inverse, r[0], r[1], o0, o1, r[2], r[3], diff)
				}
			}
		}
		t.Logf("inverse=%v: %d points, worst %.3g deg (budget %g)", inverse, len(rows), worst, wrefTolDeg)
	}
}

// xphPaperSmallSigma is the paper's switch to the small-angle forms
// (18)/(19) (Calabretta & Lowe 2013: sigma < 1e-4), and xphPaperSmallTheta
// the native latitude where sigma = 1e-4 by eq. (19),
// sigma = (90 - theta) sqrt(3/2) pi/180, computed here independently of
// healpix.go.
const xphPaperSmallSigma = 1e-4

var xphPaperSmallTheta = 90 - xphPaperSmallSigma/(math.Sqrt(1.5)*math.Pi/180)

// xphSigmaRoundingDeg bounds what the double-precision evaluation of
// equation (9) can contribute at native latitude theta (deg): the polar
// zones' sigma = sqrt(3 (1 - |sin theta|)), with |sin theta| rounded to
// within 2^-53, is off by up to 3 * 2^-53 / (2 sigma), and x, y move by
// up to 45 deg * that * (kappa (1 + 1)) < 2 * 45 * 3 * 2^-53 / (2 sigma).
// It is zero where the paper's small-angle form (19) applies instead,
// sigma < 1e-4 (Calabretta & Lowe 2013, eqs. 18/19; xphPaperSmallSigma,
// converted to a latitude here from the paper's own constant, NOT from the
// port's xphW3/xphW5, so a port that moves its switch is not granted the
// allowance where the paper uses (19)); wcslib, and so the port, uses (19)
// only in the north. It is also zero in the equatorial zone. Next to the
// SOUTH pole it is large (7e-11 deg at theta = -89.99, 5e-6 deg at
// -89.9999999): there wcslib keeps (9), and the measured difference from
// the exact equations is up to 4.08e-7 deg (theta = -89.9999997; a wcslib
// precision limit the port shares, pinned by TestXPHPolarPrecision).
// The north-pole points with sigma < 1e-4 get no allowance, so a port that
// dropped (19) (3.96e-8 deg off at 89.9999) or moved its switch fails.
func xphSigmaRoundingDeg(theta float64) float64 {
	if math.Abs(math.Sin(theta*math.Pi/180)) <= 2.0/3.0 || theta >= xphPaperSmallTheta || math.Abs(theta) == 90 {
		return 0
	}
	// sigma from the identity 1 - |sin theta| = 2 sin^2((90 - |theta|)/2),
	// which does not cancel.
	sigma := math.Sqrt(6) * math.Sin((90-math.Abs(theta))*math.Pi/360)
	return 2 * 45 * 3 * 0x1p-53 / (2 * sigma)
}

// Measured worst differences of the port from the exact XPH equations in
// each polar zone (|sin theta| > 2/3), over every row of the paper golden.
// wcslib 8.10, which the port follows operation for operation
// (TestWrefProjections), uses the small-angle form (19) in the forward
// direction only next to the
// NORTH pole (theta >= 89.9953); next to the south pole xphs2x keeps (9),
// sigma = sqrt(3 (1 - |sin theta|)), whose double-precision evaluation
// limits the result: 1 - |sin theta| loses its digits as theta -> -90, and
// at theta = -(90 - 3e-7) deg sin theta rounds to -1, so the point lands on
// the butterfly's outer corner, 4.08e-7 deg from its true image. The
// inverse uses (18) next to both poles and is exact to ~1e-14 deg there.
// Away from the switch the north forward still evaluates (9) (2.6e-11 deg
// at theta = 89.99, inside xphSigmaRoundingDeg; 2.4e-11 at 89.994 when
// the port's sin was correctly rounded: Go's math.Sin is one ULP off the
// correctly rounded sin(89.99 deg), which 1 - sin(theta) magnifies, so
// the documented figure is 3e-11).
//
// These are wcslib's precision, not a port fault; the port keeps them so
// that it stays identical to wcslib. TestXPHPolarPrecision fails if the
// measured worst of any zone grows past its documented figure, so the
// figures quoted in docs/VALIDATION.md stay true.
var xphPolarWorstDeg = map[string]float64{
	"forward/north": 3e-11,
	"inverse/north": 1e-13,
	"forward/south": 4.1e-7,
	"inverse/south": 1e-13,
}

func TestXPHPolarPrecision(t *testing.T) {
	g := loadXPHPaper(t)
	proj, err := Select("XPH", nil, 2)
	if err != nil {
		t.Fatal(err)
	}
	d := proj.(degAdapter)
	worst := map[string]float64{}
	where := map[string]float64{}
	for _, inverse := range []bool{false, true} {
		rows, dir := g.Forward, "forward"
		if inverse {
			rows, dir = g.Inverse, "inverse"
		}
		for _, r := range rows {
			theta := r[1]
			var o0, o1 float64
			var ok bool
			if inverse {
				theta = r[3]
				o0, o1, ok = d.inverseDeg(r[0], r[1])
			} else {
				o0, o1, ok = d.forwardDeg(r[0], r[1])
			}
			if !ok || math.Abs(math.Sin(theta*math.Pi/180)) <= 2.0/3.0 {
				continue
			}
			zone := dir + "/north"
			if theta < 0 {
				zone = dir + "/south"
			}
			if diff := wrefDiffConditioned(inverse, r[2], r[3], o0, o1); diff > worst[zone] {
				worst[zone], where[zone] = diff, theta
			}
		}
	}
	for _, zone := range []string{"forward/north", "inverse/north", "forward/south", "inverse/south"} {
		t.Logf("%s polar zone: worst %.3g deg at theta = %.10g (documented %.3g)", zone, worst[zone], where[zone], xphPolarWorstDeg[zone])
		if worst[zone] > xphPolarWorstDeg[zone] {
			t.Errorf("%s polar zone: worst difference from the exact equations %.3g deg (theta = %.10g) exceeds the documented %.3g",
				zone, worst[zone], where[zone], xphPolarWorstDeg[zone])
		}
	}
}
