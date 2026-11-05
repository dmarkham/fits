package transform

import (
	"fmt"
	"math"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// hpxPaperGolden is HPX evaluated from the equations of Calabretta &
// Roukema 2007 (MNRAS 381, 865, section 6) in mpmath at 60 digits, by
// testdata/wref/gen_hpx_paper.py: a truth for HPX that does not come from
// wcslib, so a fault shared by wcslib's hpxs2x/hpxx2s and our port would
// show here. Cases (H, K) = (4, 3) and (5, 5).
const hpxPaperGolden = "testdata/wref/paper/hpx_paper.json"

// Measured worst differences from the exact HPX equations next to the
// poles (TestHPXPaperEquations logs them; docs/VALIDATION.md quotes them).
// The port is identical to wcslib there (TestWrefProjections), so these
// are wcslib 8.10's own precision, set by double-precision rounding in
// the polar-zone formulas, not by a port fault. The test fails if the
// measured worst grows past them, so the documented figures stay true.
const (
	// forward (x, y), at theta = +-(90 - 3e-7) deg: 1 - |sin theta| =
	// 1.4e-17 is below half an ULP of 1, so sin theta rounds to 1, sigma
	// to 0, and the point lands on the facet apex, (180 / H) sigma_exact
	// away from its true image.
	hpxPolarFwdWorstDeg = 3.0e-7
	// inverse theta, at sigma ~ 1e-8 (|y| 1e-8 * 180 / H short of the
	// apex): t = 1 - sigma^2 / K rounds to (or next to) 1 and theta comes
	// back as +-90 or a few ULP-steps of asin below it.
	hpxPolarInvWorstDeg = 4.7e-7
)

// hpxPolarRoundingDeg bounds what double precision can contribute to the
// port's (and wcslib's) polar-zone result at native latitude theta (deg),
// for H facets in longitude and K in latitude:
//
//   - forward: sigma = sqrt(K (1 - |sin theta|)) with |sin theta| rounded
//     to within 2^-53 is off by K 2^-53 / (2 sigma); x and y move by
//     (180 / H) times that (|phi - phi_c| <= 180 / H). Twice that is
//     allowed.
//   - inverse: theta = asin(t), t = 1 - sigma^2 / K rounded to within
//     2^-53, and d theta = dt / cos(theta) with cos(theta) =
//     sigma sqrt(2 / K) (rad). Twice that is allowed. When sigma^2 / K
//     drops below 2^-53, t rounds to 1 and theta comes back as exactly
//     90: the difference is then 90 - |theta| itself, which this bound
//     still covers (sigma sqrt(2/K) <= 2^-53 / (sigma sqrt(2/K)) there).
//
// sigma is computed from the exact theta with the identity
// 1 - |sin theta| = 2 sin^2((90 - |theta|) / 2), which does not cancel.
// Zero in the equatorial zone.
func hpxPolarRoundingDeg(inverse bool, theta, H, K float64) float64 {
	if math.Abs(math.Sin(theta*math.Pi/180)) <= (K-1)/K || math.Abs(theta) == 90 {
		return 0
	}
	sigma := math.Sqrt(2*K) * math.Sin((90-math.Abs(theta))*math.Pi/360)
	if inverse {
		return 2 * 0x1p-53 / (sigma * math.Sqrt(2/K)) * 180 / math.Pi
	}
	return 2 * (180 / H) * K * 0x1p-53 / (2 * sigma)
}

type hpxPaperCase struct {
	H, K             float64
	Forward, Inverse [][4]float64
}

func loadHPXPaper(t *testing.T) []hpxPaperCase {
	t.Helper()
	var meta struct {
		Source string `json:"source"`
		Cases  []struct {
			H       float64 `json:"h"`
			K       float64 `json:"k"`
			Forward int     `json:"forward"`
			Inverse int     `json:"inverse"`
		} `json:"cases"`
	}
	vals := loadF64Pack(t, hpxPaperGolden, &meta)
	var cases []hpxPaperCase
	for _, c := range meta.Cases {
		cases = append(cases, hpxPaperCase{H: c.H, K: c.K,
			Forward: f64packRows4(t, &vals, c.Forward),
			Inverse: f64packRows4(t, &vals, c.Inverse)})
	}
	if len(vals) != 0 {
		t.Fatalf("%s: %d doubles left over after the rows the metadata lists", hpxPaperGolden, len(vals))
	}
	if len(cases) != 2 || len(cases[0].Forward) < 1300 || len(cases[0].Inverse) < 1150 || len(cases[1].Forward) < 400 {
		t.Fatalf("golden too small: %d cases", len(cases))
	}
	return cases
}

// hpxPaperTolDeg is TestHPXPaperEquations' base budget: two ULP of a
// coordinate in [128, 256) deg (measured worst 1.24e-14 deg beyond the
// polar allowance, both cases, both directions). Until round 2 of the proof
// audit it was wrefTolDeg (1e-11 deg), which let a port that scales the
// forward sigma by 1 + 1e-13 (up to 4.5e-12 deg off) or the inverse
// sigma^2 / K by 1 + 1e-14 pass.
const hpxPaperTolDeg = 5e-14

// TestHPXPaperEquations holds the HPX port to the paper's equations, both
// directions. The budget is hpxPaperTolDeg (5e-14 deg) with the comparison
// of the wcslib goldens plus the latitude conditioning allowance of an
// exact truth (wrefDiffConditioned). In the polar zones
// hpxPolarRoundingDeg is allowed on top: it is the double-precision limit
// of the formulas wcslib uses (and the port reproduces bit for bit); the
// measured worst near-pole differences are hpxPolarFwdWorstDeg and
// hpxPolarInvWorstDeg, both required to hold.
func TestHPXPaperEquations(t *testing.T) {
	for _, c := range loadHPXPaper(t) {
		t.Run(fmt.Sprintf("H%g_K%g", c.H, c.K), func(t *testing.T) {
			proj, err := Select("HPX", map[wcs.PVKey]float64{{Axis: 2, Index: 1}: c.H, {Axis: 2, Index: 2}: c.K}, 2)
			if err != nil {
				t.Fatal(err)
			}
			d := proj.(degAdapter)
			for _, inverse := range []bool{false, true} {
				rows := c.Forward
				if inverse {
					rows = c.Inverse
				}
				var nbad int
				var worst, worstPolar, worstRatio float64
				for _, r := range rows {
					var o0, o1 float64
					var ok bool
					theta := r[1]
					if inverse {
						o0, o1, ok = d.inverseDeg(r[0], r[1])
						theta = r[3]
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
					raw := wrefDiffConditioned(inverse, r[2], r[3], o0, o1)
					allow := hpxPolarRoundingDeg(inverse, theta, c.H, c.K)
					if allow > 0 {
						worstPolar = math.Max(worstPolar, raw)
						worstRatio = math.Max(worstRatio, raw/(allow+hpxPaperTolDeg))
					}
					diff := math.Max(0, raw-allow)
					worst = math.Max(worst, diff)
					if !(diff <= hpxPaperTolDeg) {
						nbad++
						if nbad <= 10 {
							t.Errorf("inverse=%v (%.17g, %.17g): got (%.17g, %.17g), paper (%.17g, %.17g), diff %.3g deg (rounding allowance %.3g)",
								inverse, r[0], r[1], o0, o1, r[2], r[3], raw, allow)
						}
					}
				}
				limit := hpxPolarFwdWorstDeg
				if inverse {
					limit = hpxPolarInvWorstDeg
				}
				if worstPolar > limit {
					t.Errorf("inverse=%v: worst polar-zone difference %.3g deg exceeds the documented %.3g", inverse, worstPolar, limit)
				}
				t.Logf("inverse=%v: %d points, worst beyond allowance %.3g deg (budget %g); polar zones worst %.3g deg, %.2f of the rounding bound",
					inverse, len(rows), worst, hpxPaperTolDeg, worstPolar, worstRatio)
			}
		})
	}
}
