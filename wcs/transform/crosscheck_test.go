package transform

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// TestCrossCheckFixtures compares this library's pixel→sky mapping
// against pre-computed reference values from wcslib / astropy.wcs for a
// set of real observatory headers.
//
// The reference values live in testdata/wcs_crosscheck.json. That file
// is meant to be generated once by a Python helper (see
// testdata/generate_wcs_crosscheck.py) using astropy.wcs on a set of
// real FITS headers from HST, JWST, Chandra, Gaia, Planck, WMAP, COBE,
// and 2MASS. The Go test reads the file and asserts our library agrees
// with each tabulated (pixel, sky) pair to within tolerance.
//
// When the JSON file is absent, the test is skipped — this allows
// development to proceed without the Python toolchain installed. The
// v1.1 release gate requires the file to be present and the test to
// pass.
func TestCrossCheckFixtures(t *testing.T) {
	path := filepath.Join("..", "..", "testdata", "wcs_crosscheck.json")
	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			t.Skip("testdata/wcs_crosscheck.json not present; skipping wcslib cross-check")
		}
		t.Fatal(err)
	}
	var fixtures []CrossCheckFixture
	if err := json.Unmarshal(data, &fixtures); err != nil {
		t.Fatal(err)
	}
	if len(fixtures) == 0 {
		t.Skip("no cross-check fixtures declared")
	}
	for _, f := range fixtures {
		f := f
		t.Run(f.Name, func(t *testing.T) {
			pc := f.PC
			if pc == nil && f.CD == nil {
				// Default identity PC when neither matrix is in the JSON.
				pc = [][]float64{{1, 0}, {0, 1}}
			}
			w := &wcs.Header{
				NAxis:    2,
				CType:    []string{f.CType1, f.CType2},
				CRVal:    []float64{f.CRVal1, f.CRVal2},
				CRPix:    []float64{f.CRPix1, f.CRPix2},
				CDelt:    []float64{f.CDelt1, f.CDelt2},
				PC:       pc,
				CD:       f.CD,
				PV:       unmarshalPVMap(f.PV),
				LonAxis:  1,
				LatAxis:  2,
				AxisType: []string{splitAxisType(f.CType1), splitAxisType(f.CType2)},
				ProjCode: []string{splitProjCode(f.CType1), splitProjCode(f.CType2)},
			}
			tr, err := New(w)
			if err != nil {
				t.Fatalf("New: %v", err)
			}
			tol := f.Tolerance
			if tol == 0 {
				tol = 1e-9
			}
			for _, tv := range f.TestVectors {
				a, d, err := tr.PixelToSky(tv.Pixel[0], tv.Pixel[1])
				if err != nil {
					t.Errorf("PixelToSky(%v): %v", tv.Pixel, err)
					continue
				}
				if math.Abs(a-tv.Sky[0]) > tol || math.Abs(d-tv.Sky[1]) > tol {
					t.Errorf("pixel %v: got (%v,%v), want (%v,%v), tol %v",
						tv.Pixel, a, d, tv.Sky[0], tv.Sky[1], tol)
				}
			}
		})
	}
}

// CrossCheckFixture is the JSON-serializable schema for one test fixture.
// The Python reference generator emits a list of these.
type CrossCheckFixture struct {
	Name      string      `json:"name"`
	Source    string      `json:"source"` // mission / instrument, for documentation
	CType1    string      `json:"ctype1"`
	CType2    string      `json:"ctype2"`
	CRVal1    float64     `json:"crval1"`
	CRVal2    float64     `json:"crval2"`
	CRPix1    float64     `json:"crpix1"`
	CRPix2    float64     `json:"crpix2"`
	CDelt1    float64     `json:"cdelt1"`
	CDelt2    float64     `json:"cdelt2"`
	PC        [][]float64 `json:"pc,omitempty"`
	CD        [][]float64 `json:"cd,omitempty"`
	PV        []PVJSON    `json:"pv,omitempty"`
	Tolerance float64     `json:"tolerance"`
	TestVectors []CrossCheckVector `json:"test_vectors"`
}

// CrossCheckVector is one (pixel, sky) pair produced by the reference.
type CrossCheckVector struct {
	Pixel [2]float64 `json:"pixel"`
	Sky   [2]float64 `json:"sky"`
}

// PVJSON is the JSON serialization of one PV key-value pair.
type PVJSON struct {
	Axis  int     `json:"axis"`
	Index int     `json:"index"`
	Value float64 `json:"value"`
}

// unmarshalPVMap converts a slice of PVJSON into the map format used by
// *wcs.Header.
func unmarshalPVMap(pvs []PVJSON) map[wcs.PVKey]float64 {
	if len(pvs) == 0 {
		return nil
	}
	m := make(map[wcs.PVKey]float64, len(pvs))
	for _, p := range pvs {
		m[wcs.PVKey{Axis: p.Axis, Index: p.Index}] = p.Value
	}
	return m
}

// splitAxisType / splitProjCode duplicate the logic in wcs.splitCType
// because it is unexported in the wcs package. For a CTYPE like
// "RA---TAN" returns ("RA", "TAN"); for a linear CTYPE returns (raw, "").
func splitAxisType(ctype string) string {
	if len(ctype) >= 8 && ctype[len(ctype)-4] == '-' {
		return trimRightDashes(ctype[:len(ctype)-3])
	}
	return ctype
}

func splitProjCode(ctype string) string {
	if len(ctype) >= 8 && ctype[len(ctype)-4] == '-' {
		return ctype[len(ctype)-3:]
	}
	return ""
}

func trimRightDashes(s string) string {
	for len(s) > 0 && s[len(s)-1] == '-' {
		s = s[:len(s)-1]
	}
	return s
}
