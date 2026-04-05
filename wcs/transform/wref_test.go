package transform

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// wcslib reference validation. These tests load the golden JSON files
// produced by testdata/wref/wref.c (which links libwcs) and diff the
// Go port output against the real wcslib implementation.
//
// Goldens are checked into the repo so CI doesn't need libwcs. To
// regenerate them locally:
//
//	cd wcs/transform/testdata/wref && make clean && make goldens
//
// If the golden directory is missing the tests skip rather than fail,
// so a fresh clone without the fixtures is still green.

const wrefGoldenDir = "testdata/wref/golden"

// wrefCase mirrors the schema wref.c writes. All angles are in degrees,
// per wcslib's native convention.
type wrefCase struct {
	Name    string    `json:"name"`
	Code    string    `json:"code"`
	PV      []float64 `json:"pv"`
	NPV     int       `json:"npv"`
	Forward []struct {
		Phi   float64 `json:"phi"`
		Theta float64 `json:"theta"`
		X     float64 `json:"x"`
		Y     float64 `json:"y"`
		Stat  int     `json:"stat"`
	} `json:"forward"`
	Inverse []struct {
		X     float64 `json:"x"`
		Y     float64 `json:"y"`
		Phi   float64 `json:"phi"`
		Theta float64 `json:"theta"`
		Stat  int     `json:"stat"`
	} `json:"inverse"`
	ForwardStatus int `json:"forward_status"`
	InverseStatus int `json:"inverse_status"`
}

func loadWrefCase(t *testing.T, name string) (*wrefCase, bool) {
	t.Helper()
	path := filepath.Join(wrefGoldenDir, name+".json")
	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, true
		}
		t.Fatalf("read %s: %v", path, err)
	}
	var c wrefCase
	if err := json.Unmarshal(data, &c); err != nil {
		t.Fatalf("parse %s: %v", path, err)
	}
	return &c, false
}

// wrefCases lists every golden file. Keep in sync with wref.c.
var wrefCases = []string{
	"szp_standard", "szp_infinity", "szp_tilted",
	"azp_standard", "azp_tilted",
	"zpn_linear", "zpn_quadratic", "zpn_cubic", "zpn_quintic",
	"tsc_allfaces", "csc_allfaces", "qsc_allfaces",
	"mol_standard", "pco_standard",
}

// buildProjectionFromWref constructs a Go Projection equivalent to a
// wcslib prjprm initialized via prjset(code, pv[]). wcslib's PV indexing
// starts at 0 and is attached to the latitude axis in our wcs package
// convention (latAxis=2, PVKey{Axis:2, Index:i}).
func buildProjectionFromWref(t *testing.T, code string, pv []float64) Projection {
	t.Helper()
	pvMap := make(map[wcs.PVKey]float64)
	for i, v := range pv {
		pvMap[wcs.PVKey{Axis: 2, Index: i}] = v
	}
	proj, err := Select(code, pvMap, 2)
	if err != nil {
		t.Fatalf("Select(%s): %v", code, err)
	}
	return proj
}

// TestWrefProjections diffs the Go port against wcslib on every
// golden case: both forward and inverse directions, skipping points
// where wcslib flagged the pixel as bad (stat != 0).
//
// Tolerances: 1e-9 degrees (~1e-11 radians), well above the
// conversion ULPs but tight enough to catch any algorithmic drift.
func TestWrefProjections(t *testing.T) {
	anyLoaded := false
	for _, name := range wrefCases {
		t.Run(name, func(t *testing.T) {
			c, skip := loadWrefCase(t, name)
			if skip {
				t.Skipf("fixture not found — run `make -C wcs/transform/testdata/wref goldens`")
				return
			}
			anyLoaded = true

			proj := buildProjectionFromWref(t, c.Code, c.PV)

			// Per-projection tolerance in output degrees. The
			// default 1e-9° (~3.5 µas) is below the noise floor of
			// our bit-exact ports. CSC relaxes to 1e-4° because
			// its polynomial evaluation is intentionally single-
			// precision (wcslib uses `float` throughout cscx2s /
			// cscs2x) and different compilers produce different
			// ULP-level rounding. The projection's own spec-
			// documented residual is ~0.02° so 1e-4° is still
			// 200× below the format's inherent precision.
			tolDeg := 1e-9
			if c.Code == "CSC" {
				tolDeg = 1e-4
			}

			// Forward: (phi, theta) -> (x, y)
			for i, row := range c.Forward {
				if row.Stat != 0 {
					continue // wcslib flagged this point — skip
				}
				phi := row.Phi * math.Pi / 180
				theta := row.Theta * math.Pi / 180
				xRad, yRad, ok := proj.Forward(phi, theta)
				if !ok {
					t.Errorf("forward[%d] (phi=%g, theta=%g): Go returned ok=false, wcslib ok",
						i, row.Phi, row.Theta)
					continue
				}
				xDeg := xRad * 180 / math.Pi
				yDeg := yRad * 180 / math.Pi
				if math.Abs(xDeg-row.X) > tolDeg || math.Abs(yDeg-row.Y) > tolDeg {
					t.Errorf("forward[%d] (phi=%g, theta=%g): got (%.15g, %.15g), want (%.15g, %.15g), diff (%g, %g)",
						i, row.Phi, row.Theta, xDeg, yDeg, row.X, row.Y,
						xDeg-row.X, yDeg-row.Y)
				}
			}

			// Inverse: (x, y) -> (phi, theta)
			for i, row := range c.Inverse {
				if row.Stat != 0 {
					continue
				}
				xRad := row.X * math.Pi / 180
				yRad := row.Y * math.Pi / 180
				phiRad, thetaRad, ok := proj.Inverse(xRad, yRad)
				if !ok {
					t.Errorf("inverse[%d] (x=%g, y=%g): Go returned ok=false, wcslib ok",
						i, row.X, row.Y)
					continue
				}
				phiDeg := phiRad * 180 / math.Pi
				thetaDeg := thetaRad * 180 / math.Pi
				// Phi wraps mod 360 — normalize both to [-180, 180] before comparing.
				dPhi := normalizeDegPi(phiDeg - row.Phi)
				if math.Abs(dPhi) > tolDeg || math.Abs(thetaDeg-row.Theta) > tolDeg {
					t.Errorf("inverse[%d] (x=%g, y=%g): got (phi=%.15g, theta=%.15g), want (%.15g, %.15g), diff (%g, %g)",
						i, row.X, row.Y, phiDeg, thetaDeg, row.Phi, row.Theta,
						dPhi, thetaDeg-row.Theta)
				}
			}
		})
	}
	if !anyLoaded {
		t.Skip("no wref fixtures — regenerate via `make -C wcs/transform/testdata/wref goldens`")
	}
}

// normalizeDegPi normalizes a degree-valued angle difference into
// (-180, 180]. Used when comparing longitudes, which are defined mod
// 360 so that 359° and -1° should diff as zero.
func normalizeDegPi(d float64) float64 {
	for d > 180 {
		d -= 360
	}
	for d <= -180 {
		d += 360
	}
	return d
}
