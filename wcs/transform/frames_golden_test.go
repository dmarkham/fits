package transform

import (
	"encoding/json"
	"fmt"
	"os"
	"testing"
)

// Sky-frame truth tests against astropy.coordinates and ERFA.
//
// testdata/frames/frames_golden.json is written by
// testdata/frames/gen_frames.py (astropy 8.0.2, pyerfa 2.0.1.5 / ERFA 2.0.1,
// pinned in testdata/frames/PROVENANCE.md). Every conversion in frames.go
// is compared point by point against it; a missing golden or a golden case
// with no Go function is a FAILURE.

const framesGolden = "testdata/frames/frames_golden.json"

// framesTolArcsec is the per-point angular budget. Every conversion is a
// product of at most a few 3x3 rotations (plus the FK4 E-term vector) in
// float64 on both sides, so the expected disagreement is a few ULP of a
// radian. Measured worst case over all 31 cases x 144 points: 4.1e-10
// arcsec (logged by the test). The budget, 1e-8 arcsec, leaves ~25x
// headroom and is six orders of magnitude below the smallest real defect
// in scope (the 0.02-0.03 arcsec ICRS/FK5 frame bias).
const framesTolArcsec = 1e-8

// framesGoldenFile: the 144 input directions are the same for every case,
// so they are stored once ("points") and each case holds its outputs in
// the same order.
type framesGoldenFile struct {
	Astropy string       `json:"astropy"`
	PyErfa  string       `json:"pyerfa"`
	Points  [][2]float64 `json:"points"`
	Cases   []struct {
		Name   string             `json:"name"`
		Params map[string]float64 `json:"params"`
		Out    [][2]float64       `json:"out"`
	} `json:"cases"`
}

// framesPoints is the number of input directions gen_frames.py writes.
const framesPoints = 144

// angSepArcsec (Vincenty form) is shared with hdr_golden_test.go.

func framesConversion(name string, p map[string]float64) (func(lon, lat float64) (float64, float64), bool) {
	switch name {
	case "icrs_to_galactic":
		return EquatorialToGalactic, true
	case "galactic_to_icrs":
		return GalacticToEquatorial, true
	case "fk5_to_galactic":
		return FK5ToGalactic, true
	case "galactic_to_fk5":
		return GalacticToFK5, true
	case "galactic_to_supergalactic":
		return GalacticToSupergalactic, true
	case "supergalactic_to_galactic":
		return SupergalacticToGalactic, true
	case "icrs_to_supergalactic":
		return EquatorialToSupergalactic, true
	case "supergalactic_to_icrs":
		return SupergalacticToEquatorial, true
	case "icrs_to_ecliptic":
		jd := p["jd"]
		return func(a, d float64) (float64, float64) { return EquatorialToEcliptic(a, d, jd) }, true
	case "ecliptic_to_icrs":
		jd := p["jd"]
		return func(a, d float64) (float64, float64) { return EclipticToEquatorial(a, d, jd) }, true
	case "precess_fk5":
		e1, e2 := p["equinox1"], p["equinox2"]
		return func(a, d float64) (float64, float64) { return PrecessFK5(a, d, e1, e2) }, true
	case "fk4_to_fk5":
		ob := p["obstime_byear"]
		return func(a, d float64) (float64, float64) { return FK4ToFK5(a, d, ob) }, true
	case "fk5_to_fk4":
		ob := p["obstime_byear"]
		return func(a, d float64) (float64, float64) { return FK5ToFK4(a, d, ob) }, true
	}
	return nil, false
}

func loadFramesGolden(t *testing.T) framesGoldenFile {
	t.Helper()
	raw, err := os.ReadFile(framesGolden)
	if err != nil {
		t.Fatalf("committed frames golden missing (regenerate with gen_frames.py, see testdata/frames/PROVENANCE.md): %v", err)
	}
	var g framesGoldenFile
	if err := json.Unmarshal(raw, &g); err != nil {
		t.Fatalf("decode %s: %v", framesGolden, err)
	}
	if g.Astropy != "8.0.2" || g.PyErfa != "2.0.1.5" {
		t.Fatalf("%s produced by astropy %s / pyerfa %s; PROVENANCE.md pins 8.0.2 / 2.0.1.5", framesGolden, g.Astropy, g.PyErfa)
	}
	if len(g.Cases) == 0 {
		t.Fatalf("%s has no cases", framesGolden)
	}
	if len(g.Points) != framesPoints {
		t.Fatalf("%s has %d input points, want %d", framesGolden, len(g.Points), framesPoints)
	}
	return g
}

// TestFramesAstropyGolden compares every frame conversion against astropy
// / ERFA. Output longitudes must also lie in [0, 360) and latitudes in
// [-90, 90] as documented.
func TestFramesAstropyGolden(t *testing.T) {
	g := loadFramesGolden(t)
	worst := 0.0
	for _, c := range g.Cases {
		fn, ok := framesConversion(c.Name, c.Params)
		if !ok {
			t.Errorf("golden case %q has no Go conversion", c.Name)
			continue
		}
		if len(c.Out) != len(g.Points) {
			t.Fatalf("case %s: %d outputs for %d input points", c.Name, len(c.Out), len(g.Points))
		}
		caseWorst, bad := 0.0, 0
		for i, in := range g.Points {
			lon, lat := fn(in[0], in[1])
			if !(lon >= 0 && lon < 360) || !(lat >= -90 && lat <= 90) {
				t.Errorf("%s %v: output (%v, %v) outside [0,360) x [-90,90]", c.Name, in, lon, lat)
			}
			sep := angSepArcsec(lon, lat, c.Out[i][0], c.Out[i][1])
			if sep > caseWorst {
				caseWorst = sep
			}
			if !(sep <= framesTolArcsec) {
				bad++
				if bad <= 3 {
					t.Errorf("%s %v %v: got (%.12f, %.12f), astropy/ERFA (%.12f, %.12f): %.3g arcsec",
						c.Name, c.Params, in, lon, lat, c.Out[i][0], c.Out[i][1], sep)
				}
			}
		}
		if bad > 3 {
			t.Errorf("%s %v: %d/%d points over %.0e arcsec, worst %.3g arcsec", c.Name, c.Params, bad, len(g.Points), framesTolArcsec, caseWorst)
		}
		if caseWorst > worst {
			worst = caseWorst
		}
	}
	t.Logf("worst disagreement over %d cases: %.3g arcsec (budget %.0e)", len(g.Cases), worst, framesTolArcsec)
}

// framesRequiredCases lists every (conversion, parameters) case the
// committed golden must hold, with its point count, so a golden with
// cases deleted fails here even when gen_frames.py cannot be run.
var framesRequiredCases = func() map[string]int {
	m := map[string]int{}
	for _, n := range []string{"icrs_to_galactic", "galactic_to_icrs", "fk5_to_galactic", "galactic_to_fk5",
		"galactic_to_supergalactic", "supergalactic_to_galactic", "icrs_to_supergalactic", "supergalactic_to_icrs"} {
		m[n+" map[]"] = framesPoints
	}
	for _, jd := range []string{"2.451545e+06", "2.43328242345905e+06", "2.4606765e+06", "2.4880695e+06", "2.41502031352e+06"} {
		m["icrs_to_ecliptic map[jd:"+jd+"]"] = framesPoints
		m["ecliptic_to_icrs map[jd:"+jd+"]"] = framesPoints
	}
	for _, e := range []string{"1950 equinox2:2000", "2000 equinox2:1950", "2000 equinox2:2050", "1900 equinox2:2100", "2025.5 equinox2:1987.25"} {
		m["precess_fk5 map[equinox1:"+e+"]"] = framesPoints
	}
	for _, ob := range []string{"1950", "1900", "1984", "2000"} {
		m["fk4_to_fk5 map[obstime_byear:"+ob+"]"] = framesPoints
		m["fk5_to_fk4 map[obstime_byear:"+ob+"]"] = framesPoints
	}
	return m
}()

// TestFramesGoldenComplete: the golden holds exactly the required cases,
// each once and with its full point set.
func TestFramesGoldenComplete(t *testing.T) {
	g := loadFramesGolden(t)
	seen := map[string]int{}
	for _, c := range g.Cases {
		key := c.Name + " " + fmt.Sprint(c.Params)
		if _, ok := framesRequiredCases[key]; !ok {
			t.Errorf("golden case %q is not in framesRequiredCases", key)
		}
		if seen[key]++; seen[key] > 1 {
			t.Errorf("golden case %q appears twice", key)
		}
		if want := framesRequiredCases[key]; want != 0 && len(c.Out) != want {
			t.Errorf("golden case %q has %d points, want %d", key, len(c.Out), want)
		}
	}
	for key := range framesRequiredCases {
		if seen[key] == 0 {
			t.Errorf("golden is missing required case %q", key)
		}
	}
}

// TestVec3ToLonLatWrap: a direction a hair below longitude 0 has
// atan2 = -tiny, and -tiny + 360 rounds to exactly 360 in float64. The
// documented range is [0, 360), so it must come back as 0.
func TestVec3ToLonLatWrap(t *testing.T) {
	for _, y := range []float64{-1e-20, -1e-300, -5e-324} {
		lon, lat := vec3ToLonLat([3]float64{1, y, 0})
		if lon != 0 || lat != 0 {
			t.Errorf("vec3ToLonLat(1, %g, 0) = (%v, %v), want (0, 0)", y, lon, lat)
		}
	}
	// Through a public conversion: the galactic direction (360 - tiny, 0)
	// round-trips to a longitude in [0, 360).
	if lon, _ := GalacticToSupergalactic(SupergalacticToGalactic(-1e-17, 0)); !(lon >= 0 && lon < 360) {
		t.Errorf("round trip of SGL -1e-17 gave lon %v", lon)
	}
}
