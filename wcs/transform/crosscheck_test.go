package transform

import (
	"encoding/json"
	"math"
	"os"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCrossCheckFixtures compares this library with astropy.wcs (astropy
// 8.0.2, bundling wcslib 8.6) through the FULL public path: the header
// card images stored in testdata/wcs_crosscheck.json go through
// header.ParseCards -> wcs.Parse -> New (the parser is under test too),
// and both directions are checked: pixel -> sky against
// WCS.all_pix2world and sky -> pixel against WCS.all_world2pix.
//
// The JSON is produced by testdata/generate_wcs_crosscheck.py (see its
// docstring; TestCrossCheckFixturesReproduce re-runs it). A missing or
// empty file is a failure, never a skip.
//
// Tolerances (per fixture, from the JSON): 1e-8 arcsec pixel -> sky and
// 1e-9 px sky -> pixel, the same budgets as the wcslib header goldens
// (hdr_golden_test.go) for the same reasons; the smallest bug the WCS
// audit found was ~0.01 arcsec.
func TestCrossCheckFixtures(t *testing.T) {
	path := testutil.TestdataPath(t, "wcs_crosscheck.json")
	testutil.RequireFile(t, path)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var doc crossCheckDoc
	if err := json.Unmarshal(data, &doc); err != nil {
		t.Fatal(err)
	}
	if len(doc.Fixtures) == 0 {
		t.Fatalf("%s has no fixtures", path)
	}
	t.Logf("truth: astropy %s, wcslib %s, numpy %s", doc.Astropy, doc.Wcslib, doc.Numpy)
	for _, f := range doc.Fixtures {
		t.Run(f.Name, func(t *testing.T) {
			if len(f.P2S) == 0 || len(f.S2P) == 0 {
				t.Fatalf("fixture %s has no vectors", f.Name)
			}
			w, err := hdrParse(f.Cards, "")
			if err != nil {
				t.Fatalf("wcs.Parse: %v", err)
			}
			tr, err := New(w)
			if err != nil {
				t.Fatalf("New: %v", err)
			}
			var maxAS, maxPx float64
			for _, v := range f.P2S {
				a, d, err := tr.PixelToSky(v.Pixel[0], v.Pixel[1])
				if err != nil {
					t.Errorf("PixelToSky(%v): %v", v.Pixel, err)
					continue
				}
				sep := angSepArcsec(a, d, v.Sky[0], v.Sky[1])
				maxAS = math.Max(maxAS, sep)
				if !(sep <= f.TolArcsec) {
					t.Errorf("pixel %v: got (%.15g, %.15g), astropy (%.15g, %.15g): %.3g arcsec > %g",
						v.Pixel, a, d, v.Sky[0], v.Sky[1], sep, f.TolArcsec)
				}
			}
			for _, v := range f.S2P {
				x, y, err := tr.SkyToPixel(v.Sky[0], v.Sky[1])
				if err != nil {
					t.Errorf("SkyToPixel(%v): %v", v.Sky, err)
					continue
				}
				d := math.Max(math.Abs(x-v.Pixel[0]), math.Abs(y-v.Pixel[1]))
				maxPx = math.Max(maxPx, d)
				if !(d <= f.TolPx) {
					t.Errorf("sky %v: got pixel (%.15g, %.15g), astropy (%.15g, %.15g): %.3g px > %g",
						v.Sky, x, y, v.Pixel[0], v.Pixel[1], d, f.TolPx)
				}
			}
			t.Logf("%s: %d p2s max %.3g arcsec, %d s2p max %.3g px", f.Source, len(f.P2S), maxAS, len(f.S2P), maxPx)
		})
	}
}

// crossCheckDoc is the JSON written by testdata/generate_wcs_crosscheck.py.
type crossCheckDoc struct {
	Generator string              `json:"generator"`
	Astropy   string              `json:"astropy"`
	Wcslib    string              `json:"wcslib"`
	Numpy     string              `json:"numpy"`
	Fixtures  []CrossCheckFixture `json:"fixtures"`
}

// CrossCheckFixture is one header with astropy's answers in both
// directions.
type CrossCheckFixture struct {
	Name      string             `json:"name"`
	Source    string             `json:"source"`
	Cards     []string           `json:"cards"` // 80-column card images, END excluded
	TolArcsec float64            `json:"tol_arcsec"`
	TolPx     float64            `json:"tol_px"`
	P2S       []CrossCheckVector `json:"p2s"`
	S2P       []CrossCheckVector `json:"s2p"`
}

// CrossCheckVector is one (pixel, sky) pair (0-based pixels, degrees).
type CrossCheckVector struct {
	Pixel [2]float64 `json:"pixel"`
	Sky   [2]float64 `json:"sky"`
}
