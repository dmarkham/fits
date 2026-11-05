package transform

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits"
	"github.com/dmarkham/fits/wcs"
)

// Projection proofs through the full header pipeline.
//
// The header cases are ordinary TestHdrGoldens cases (hdrcases/*.case ->
// golden/hdr/*.json, wcslib 8.10 via hdr2sky); this file adds what is
// specific to the projection proof:
//
//   - prj_*.case   (gen_prj_cases.py): every projection code at CRVAL2 = 0
//     and at CRVAL2 != theta_0, with pixel grids running past the
//     projection boundary, so off-projection rejection is compared with
//     wcslib point by point; TestHdrProjectionCoverage enforces that.
//   - sample_*.case (gen_sample_cases.py): real headers - Calabretta's
//     Paper II example maps (1904-66_*, every projection incl. CSC),
//     HPXcvt HPX/XPH image headers, and MSX GLON-CAR / ROSAT GLON-AIT
//     survey headers with LONPOLE = 0.
//   - TestHPXcvtPixelCentres: an independent, data-level check of a real
//     HPXcvt image against healpy.

// prjHdrKnownBad lists projection-proof header cases that disagree with
// wcslib because of an open bug OUTSIDE the projection code; they are
// merged into hdrKnownBad (hdr_golden_test.go), so they stop failing the
// suite but start failing it the moment they agree (the fix must then
// delete the entry here). Empty since the WCS pipeline work
// (proof/wcs-pipeline) landed: explicit LONPOLE on non-zenithal
// projections, the PV1_0..2 native fiducial point and latitude-first
// headers (sample_msx/rosat_*, pvfid_*, prjswap_*) all agree with wcslib.
var prjHdrKnownBad = map[string]map[string]string{}

func init() {
	for name, dirs := range prjHdrKnownBad {
		if hdrKnownBad[name] == nil {
			hdrKnownBad[name] = map[string]string{}
		}
		for d, ref := range dirs {
			hdrKnownBad[name][d] = ref
		}
	}
}

// TestHdrProjectionCoverage requires, for every projection code we claim,
// a full-pipeline wcslib case at CRVAL2 = 0 and one at CRVAL2 != theta_0,
// and that the pixel grids of the bounded projections cross the boundary
// (wcslib rejected at least one pixel). It also requires the real-sample
// cases the projection proof depends on.
func TestHdrProjectionCoverage(t *testing.T) {
	type cov struct{ c0, cOther, rejects bool }
	got := map[string]*cov{}
	for _, name := range hdrCaseNames(t) {
		if !strings.HasPrefix(name, "prj_") {
			continue
		}
		g := loadHdrGolden(t, name)
		if len(g.WCS.Ctype) < 2 || len(g.WCS.Crval) < 2 {
			t.Fatalf("%s: golden has no celestial ctype/crval", name)
		}
		code := g.WCS.Ctype[0][len(g.WCS.Ctype[0])-3:]
		c := got[code]
		if c == nil {
			c = &cov{}
			got[code] = c
		}
		tr, err := hdrTransform(g)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		theta0 := radToDeg(tr.proj.Theta0())
		switch crval2 := g.WCS.Crval[1]; {
		case crval2 == 0:
			c.c0 = true
		case math.Abs(crval2-theta0) > 1e-9:
			c.cOther = true
		}
		for _, p := range g.P2S {
			if p.Status != 0 {
				c.rejects = true
			}
		}
	}
	// These map a hemisphere or more onto an unbounded plane, so a finite
	// pixel grid need not contain an off-projection pixel (their
	// rejections are exercised on the sky side, and in TestWrefProjections).
	unbounded := map[string]bool{"STG": true, "AIR": true, "TAN": true, "AZP": true, "SZP": true}
	for code := range allProjectionCodes {
		c := got[code]
		switch {
		case c == nil:
			t.Errorf("%s: no prj_* full-pipeline wcslib case", code)
		case !c.c0:
			t.Errorf("%s: no full-pipeline case at CRVAL2 = 0", code)
		case !c.cOther:
			t.Errorf("%s: no full-pipeline case at CRVAL2 != theta_0", code)
		case !c.rejects && !unbounded[code]:
			t.Errorf("%s: no full-pipeline case crosses the projection boundary (wcslib rejected nothing)", code)
		}
	}
	for _, name := range []string{
		"sample_1904-66_csc", "sample_hpxcvt_hpx_gal", "sample_hpxcvt_xph_n_gal",
		"sample_msx_glon_car_lonpole0", "sample_rosat_glon_ait_lonpole0",
	} {
		if _, err := os.Stat(filepath.Join(hdrCaseDir, name+".case")); err != nil {
			t.Errorf("required real-sample case %s missing: %v", name, err)
		}
	}
}

// TestHPXcvtPixelCentres checks a real HPX image written by wcslib's
// HPXcvt (testdata/wref/samples/hpxcvt_hpx_gal_n8.fits, NSIDE = 8, made by
// gen_sample_cases.py) end to end through this library: fits.Open ->
// wcs.Parse -> New -> PixelToSky at every filled pixel centre. Each pixel
// holds the RING index of the HEALPix pixel HPXcvt placed there, so its
// centre must be healpy.pix2ang(8, index) (samples/*.healpy.json, healpy
// 1.20.1) - an oracle independent of wcslib.
//
// Budget 1e-5 deg: HPXcvt writes CDELT and PC to 8 significant digits
// (CDELT = 7.9549513 for 45/(8 sqrt 2) deg), which alone moves pixel
// centres up to 2.6e-6 deg (the same figure wcslib gives through astropy).
// A pixel in the wrong facet would be off by degrees (pixel size ~8 deg).
//
// HPXcvt fills the duplicated edge pixels of the HPX layout; 56 of them sit
// exactly on the (staircase) projection boundary, numerically just outside
// it, and wcslib's own x2s rejects them (TestHdrGoldens/
// sample_hpxcvt_hpx_gal compares every pixel centre with wcslib, accept/
// reject included). Here a pixel we reject must be one wcslib rejects.
func TestHPXcvtPixelCentres(t *testing.T) {
	const dir = "testdata/wref/samples"
	f, err := fits.Open(filepath.Join(dir, "hpxcvt_hpx_gal_n8.fits"))
	if err != nil {
		t.Fatalf("open committed HPXcvt sample: %v", err)
	}
	defer f.Close()
	img, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	w, err := wcs.Parse(img.Header())
	if err != nil {
		t.Fatal(err)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	pix, err := fits.ReadPixels[float32](img)
	if err != nil {
		t.Fatal(err)
	}
	raw, err := os.ReadFile(filepath.Join(dir, "hpxcvt_hpx_gal_n8.healpy.json"))
	if err != nil {
		t.Fatalf("healpy golden missing: %v", err)
	}
	var g struct {
		Source string       `json:"source"`
		Nside  int          `json:"nside"`
		Rows   [][5]float64 `json:"rows"`
	}
	if err := json.Unmarshal(raw, &g); err != nil {
		t.Fatal(err)
	}
	wcslibRejects := map[[2]int]bool{}
	for _, p := range loadHdrGolden(t, "sample_hpxcvt_hpx_gal").P2S {
		if p.Status != 0 && p.X == math.Trunc(p.X) && p.Y == math.Trunc(p.Y) {
			wcslibRejects[[2]int{int(p.X), int(p.Y)}] = true
		}
	}
	const nx, ny = 40, 40
	if len(pix) != nx*ny || len(g.Rows) == 0 {
		t.Fatalf("sample has %d pixels, healpy golden %d rows", len(pix), len(g.Rows))
	}
	seen := map[int]bool{}
	var maxSep float64
	var nrej int
	for _, r := range g.Rows {
		x, y, idx := int(r[0]), int(r[1]), int(r[2])
		if v := pix[(y-1)*nx+(x-1)]; int(v) != idx {
			t.Fatalf("pixel (%d,%d): image value %v, golden index %d (golden not made from this file)", x, y, v, idx)
		}
		seen[idx] = true
		lon, lat, err := tr.PixelToSky(float64(x-1), float64(y-1))
		if err != nil {
			if !wcslibRejects[[2]int{x, y}] {
				t.Errorf("pixel (%d,%d) HEALPix %d: rejected (%v) but wcslib accepts it", x, y, idx, err)
			}
			nrej++
			continue
		}
		sep := angSepArcsec(lon, lat, r[3], r[4]) / 3600
		maxSep = math.Max(maxSep, sep)
		if sep > 1e-5 {
			t.Errorf("pixel (%d,%d) HEALPix %d: got (%.9f, %.9f), healpy (%.9f, %.9f): %.3g deg", x, y, idx, lon, lat, r[3], r[4], sep)
		}
	}
	if len(seen) != 12*g.Nside*g.Nside {
		t.Errorf("image holds %d distinct HEALPix pixels, want %d", len(seen), 12*g.Nside*g.Nside)
	}
	t.Logf("%d filled pixels (%d distinct HEALPix indices), %d on-boundary pixels rejected (as by wcslib), max |sky - healpy| = %.3g deg (%s)",
		len(g.Rows), len(seen), nrej, maxSep, g.Source)
}
