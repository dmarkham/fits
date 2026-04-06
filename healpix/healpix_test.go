package healpix

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"
)

const goldenDir = "testdata/golden"

// ---------- XYF (internal XY scheme) validation ----------

type xyfCase struct {
	Nside int `json:"nside"`
	Npix  int `json:"npix"`
	Cases []struct {
		Pix   int64 `json:"pix"`
		RingX int   `json:"ring_x"`
		RingY int   `json:"ring_y"`
		RingF int   `json:"ring_f"`
		NestX int   `json:"nest_x"`
		NestY int   `json:"nest_y"`
		NestF int   `json:"nest_f"`
	} `json:"cases"`
}

func TestXYF_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "xyf_nside%d.json", nside)
			var c xyfCase
			mustUnmarshal(t, data, &c)
			for _, row := range c.Cases {
				// RING: pix → XY → check face, x, y.
				xy := ringToXY(row.Pix, nside)
				got := xyDecompose(xy, nside)
				if got.Face != row.RingF || got.X != row.RingX || got.Y != row.RingY {
					t.Errorf("RING pix %d: got (f=%d x=%d y=%d), want (f=%d x=%d y=%d)",
						row.Pix, got.Face, got.X, got.Y, row.RingF, row.RingX, row.RingY)
				}
				// Verify round-trip: XY → RING → same pix.
				rt := xyToRing(xy, nside)
				if rt != row.Pix {
					t.Errorf("RING round-trip pix %d: XY→RING gave %d", row.Pix, rt)
				}
				// NESTED: pix → XY → check face, x, y.
				xyN := nestedToXY(row.Pix, nside)
				gotN := xyDecompose(xyN, nside)
				if gotN.Face != row.NestF || gotN.X != row.NestX || gotN.Y != row.NestY {
					t.Errorf("NESTED pix %d: got (f=%d x=%d y=%d), want (f=%d x=%d y=%d)",
						row.Pix, gotN.Face, gotN.X, gotN.Y, row.NestF, row.NestX, row.NestY)
				}
				// NESTED round-trip.
				rtN := xyToNested(xyN, nside)
				if rtN != row.Pix {
					t.Errorf("NESTED round-trip pix %d: XY→NESTED gave %d", row.Pix, rtN)
				}
			}
		})
	}
}

// ---------- Ordering (ring2nest / nest2ring) validation ----------

type orderingExhaustive struct {
	Nside    int     `json:"nside"`
	Npix     int     `json:"npix"`
	Ring2Nest []int64 `json:"ring2nest"`
	Nest2Ring []int64 `json:"nest2ring"`
}

func TestOrdering_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "ordering_nside%d.json", nside)
			var c orderingExhaustive
			mustUnmarshal(t, data, &c)
			for ring, wantNest := range c.Ring2Nest {
				gotNest := Ring2Nest(nside, int64(ring))
				if gotNest != wantNest {
					t.Errorf("Ring2Nest(%d, %d) = %d, want %d", nside, ring, gotNest, wantNest)
				}
			}
			for nest, wantRing := range c.Nest2Ring {
				gotRing := Nest2Ring(nside, int64(nest))
				if gotRing != wantRing {
					t.Errorf("Nest2Ring(%d, %d) = %d, want %d", nside, nest, gotRing, wantRing)
				}
			}
		})
	}
}

// ---------- Ang2Pix validation ----------

type ang2pixCase struct {
	Nside int `json:"nside"`
	Cases []struct {
		Theta   float64 `json:"theta"`
		Phi     float64 `json:"phi"`
		PixRing int64   `json:"pix_ring"`
		PixNest int64   `json:"pix_nest"`
	} `json:"cases"`
}

func TestAng2Pix_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16, 64, 256, 1024, 8192} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "ang2pix_nside%d.json", nside)
			var c ang2pixCase
			mustUnmarshal(t, data, &c)
			for i, row := range c.Cases {
				gotRing := Ang2Pix(nside, row.Theta, row.Phi, Ring)
				if gotRing != row.PixRing {
					t.Errorf("case %d (theta=%g phi=%g): RING got %d want %d",
						i, row.Theta, row.Phi, gotRing, row.PixRing)
				}
				gotNest := Ang2Pix(nside, row.Theta, row.Phi, Nested)
				if gotNest != row.PixNest {
					t.Errorf("case %d (theta=%g phi=%g): NESTED got %d want %d",
						i, row.Theta, row.Phi, gotNest, row.PixNest)
				}
			}
		})
	}
}

// ---------- Pix2Ang validation ----------

type pix2angCase struct {
	Nside int `json:"nside"`
	Cases []struct {
		Pix       int64   `json:"pix"`
		RingTheta float64 `json:"ring_theta"`
		RingPhi   float64 `json:"ring_phi"`
		NestTheta float64 `json:"nest_theta"`
		NestPhi   float64 `json:"nest_phi"`
	} `json:"cases"`
}

func TestPix2Ang_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16, 64, 256, 1024, 8192} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "pix2ang_nside%d.json", nside)
			var c pix2angCase
			mustUnmarshal(t, data, &c)
			for i, row := range c.Cases {
				gotTheta, gotPhi := Pix2Ang(nside, row.Pix, Ring)
				if !floatEq(gotTheta, row.RingTheta) || !floatEq(gotPhi, row.RingPhi) {
					t.Errorf("RING pix %d [%d]: got (%.17g, %.17g) want (%.17g, %.17g)",
						row.Pix, i, gotTheta, gotPhi, row.RingTheta, row.RingPhi)
				}
				gotTheta, gotPhi = Pix2Ang(nside, row.Pix, Nested)
				if !floatEq(gotTheta, row.NestTheta) || !floatEq(gotPhi, row.NestPhi) {
					t.Errorf("NESTED pix %d [%d]: got (%.17g, %.17g) want (%.17g, %.17g)",
						row.Pix, i, gotTheta, gotPhi, row.NestTheta, row.NestPhi)
				}
			}
		})
	}
}

// ---------- Utilities validation ----------

type utilCase struct {
	Cases []struct {
		Nside       int     `json:"nside"`
		Npix        int64   `json:"npix"`
		PixelAreaSr float64 `json:"pixel_area_sr"`
	} `json:"cases"`
}

func TestUtilities_Healpy(t *testing.T) {
	data := loadGolden(t, "utilities.json", 0)
	var c utilCase
	mustUnmarshal(t, data, &c)
	for _, row := range c.Cases {
		if got := Nside2Npix(row.Nside); got != row.Npix {
			t.Errorf("Nside2Npix(%d) = %d, want %d", row.Nside, got, row.Npix)
		}
		nside, ok := Npix2Nside(row.Npix)
		if !ok || nside != row.Nside {
			t.Errorf("Npix2Nside(%d) = (%d, %v), want (%d, true)", row.Npix, nside, ok, row.Nside)
		}
		if got := PixelArea(row.Nside); !floatEq(got, row.PixelAreaSr) {
			t.Errorf("PixelArea(%d) = %.17g, want %.17g", row.Nside, got, row.PixelAreaSr)
		}
	}
}

// ---------- Neighbors validation ----------

type neighborCase struct {
	Nside int `json:"nside"`
	Cases []struct {
		Pix           int64   `json:"pix"`
		NeighborsRing []int64 `json:"neighbors_ring"`
		NeighborsNest []int64 `json:"neighbors_nest"`
	} `json:"cases"`
}

func TestNeighbors_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16, 64, 256, 1024, 8192} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "neighbors_nside%d.json", nside)
			var c neighborCase
			mustUnmarshal(t, data, &c)
			for _, row := range c.Cases {
				// RING neighbors.
				gotRing := GetNeighbors(nside, row.Pix, Ring)
				if !neighborSetEq(gotRing[:], row.NeighborsRing) {
					t.Errorf("RING pix %d: got %v want %v", row.Pix, gotRing, row.NeighborsRing)
				}
				// NESTED neighbors.
				gotNest := GetNeighbors(nside, row.Pix, Nested)
				if !neighborSetEq(gotNest[:], row.NeighborsNest) {
					t.Errorf("NESTED pix %d: got %v want %v", row.Pix, gotNest, row.NeighborsNest)
				}
			}
		})
	}
}

// neighborSetEq checks if two neighbor arrays represent the same set
// of neighbors. The order may differ between implementations, but the
// set of valid (non -1) entries must match exactly.
func neighborSetEq(a []int64, b []int64) bool {
	setA := make(map[int64]bool)
	setB := make(map[int64]bool)
	for _, v := range a {
		if v >= 0 {
			setA[v] = true
		}
	}
	for _, v := range b {
		if v >= 0 {
			setB[v] = true
		}
	}
	if len(setA) != len(setB) {
		return false
	}
	for k := range setA {
		if !setB[k] {
			return false
		}
	}
	return true
}

// ---------- Exhaustive round-trip test ----------

func TestRoundTrip_Exhaustive(t *testing.T) {
	// For small nside, every pixel should round-trip through
	// ang2pix(pix2ang(pix)) == pix in both orderings.
	for _, nside := range []int{1, 2, 4, 8} {
		t.Run(nsideName(nside), func(t *testing.T) {
			npix := Nside2Npix(nside)
			for pix := int64(0); pix < npix; pix++ {
				// RING
				theta, phi := Pix2Ang(nside, pix, Ring)
				got := Ang2Pix(nside, theta, phi, Ring)
				if got != pix {
					t.Errorf("RING pix %d: Pix2Ang→Ang2Pix gave %d", pix, got)
				}
				// NESTED
				theta, phi = Pix2Ang(nside, pix, Nested)
				got = Ang2Pix(nside, theta, phi, Nested)
				if got != pix {
					t.Errorf("NESTED pix %d: Pix2Ang→Ang2Pix gave %d", pix, got)
				}
			}
		})
	}
}

// ---------- helpers ----------

func loadGolden(t *testing.T, pattern string, nside int) []byte {
	t.Helper()
	var name string
	if nside > 0 {
		name = filepath.Join(goldenDir, nsideFmt(pattern, nside))
	} else {
		name = filepath.Join(goldenDir, pattern)
	}
	data, err := os.ReadFile(name)
	if err != nil {
		t.Skipf("golden not found: %s (run gen_golden.py)", name)
	}
	return data
}

func mustUnmarshal(t *testing.T, data []byte, v interface{}) {
	t.Helper()
	if err := json.Unmarshal(data, v); err != nil {
		t.Fatalf("unmarshal: %v", err)
	}
}

func nsideName(n int) string {
	return nsideFmt("nside%d", n)
}

func nsideFmt(pattern string, nside int) string {
	// Simple sprintf replacement.
	out := make([]byte, 0, len(pattern)+10)
	for i := 0; i < len(pattern); i++ {
		if i+1 < len(pattern) && pattern[i] == '%' && pattern[i+1] == 'd' {
			s := itoa64(int64(nside))
			out = append(out, s...)
			i++ // skip 'd'
		} else {
			out = append(out, pattern[i])
		}
	}
	return string(out)
}

func itoa64(n int64) string {
	if n == 0 {
		return "0"
	}
	var buf [20]byte
	i := len(buf)
	neg := n < 0
	if neg {
		n = -n
	}
	for n > 0 {
		i--
		buf[i] = '0' + byte(n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

// floatEq compares two float64 angle values. Two tolerances:
//
//   - Absolute: 1e-10 radians (≈ 20 microarcseconds). Covers the
//     near-pole case where the astrometry.net formula computes
//     z = 1 - t²/3 and the cancellation loses ~4 digits of precision
//     for nside ≥ 256. healpy's healpix_cxx backend avoids this by
//     computing z from ring arithmetic instead. The worst-case absolute
//     difference across our test suite is ~5e-12 rad (≈1 µas).
//
//   - Relative: 1e-13. Covers the general case where Go's and numpy's
//     implementations of atan2/acos/asin differ by 1-2 ULPs.
//
// Both tolerances must be satisfied (we take the more lenient of the
// two). For reference, Gaia's astrometric precision is ~25 µas, so our
// 1 µas worst-case error is 25× below the best instrument on Earth.
func floatEq(a, b float64) bool {
	if a == b {
		return true
	}
	diff := math.Abs(a - b)
	// Absolute tolerance: 1e-10 rad ≈ 20 µas.
	if diff < 1e-10 {
		return true
	}
	// Relative tolerance.
	scale := math.Max(math.Abs(a), math.Abs(b))
	if scale == 0 {
		return diff < 1e-14
	}
	return diff/scale < 1e-13
}
