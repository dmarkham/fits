package healpix

import (
	"encoding/json"
	"fmt"
	"math"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

const goldenDir = "testdata/golden"

// ---------- XYF (internal XY scheme) validation ----------

type xyfCase struct {
	Nside int  `json:"nside"`
	Npix  int  `json:"npix"`
	RingX i64s `json:"ring_x"`
	RingY i64s `json:"ring_y"`
	RingF i64s `json:"ring_f"`
	NestX i64s `json:"nest_x"`
	NestY i64s `json:"nest_y"`
	NestF i64s `json:"nest_f"`
}

// TestXYF_Healpy: pix2xyf of every pixel (both orderings) equals
// healpy.pix2xyf, and xyf2pix inverts healpy's (x, y, face) exactly.
func TestXYF_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "xyf_nside%d.json", nside)
			var c xyfCase
			mustUnmarshal(t, data, &c)
			if c.Npix != 12*nside*nside || c.RingX.Len() != c.Npix {
				t.Fatalf("golden has %d values for npix %d", c.RingX.Len(), c.Npix)
			}
			for _, sc := range []struct {
				name    string
				scheme  Scheme
				x, y, f i64s
			}{{"RING", Ring, c.RingX, c.RingY, c.RingF}, {"NESTED", Nested, c.NestX, c.NestY, c.NestF}} {
				gx, gy, gf := make([]int64, c.Npix), make([]int64, c.Npix), make([]int64, c.Npix)
				for p := range gx {
					x, y, f := pix2xyf(nside, int64(p), sc.scheme)
					gx[p], gy[p], gf[p] = int64(x), int64(y), int64(f)
				}
				wx := sc.x.Resolve(t, sc.name+" pix2xyf x", gx)
				wy := sc.y.Resolve(t, sc.name+" pix2xyf y", gy)
				wf := sc.f.Resolve(t, sc.name+" pix2xyf face", gf)
				for p := range wx {
					if rt := xyf2pix(nside, int(wx[p]), int(wy[p]), int(wf[p]), sc.scheme); rt != int64(p) {
						t.Errorf("%s xyf2pix(%d,%d,%d) = %d, want %d", sc.name, wx[p], wy[p], wf[p], rt, p)
					}
				}
			}
		})
	}
}

// ---------- Ordering (ring2nest / nest2ring) validation ----------

type orderingExhaustive struct {
	Nside     int  `json:"nside"`
	Npix      int  `json:"npix"`
	Ring2Nest i64s `json:"ring2nest"`
	Nest2Ring i64s `json:"nest2ring"`
}

func TestOrdering_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16} {
		t.Run(nsideName(nside), func(t *testing.T) {
			data := loadGolden(t, "ordering_nside%d.json", nside)
			var c orderingExhaustive
			mustUnmarshal(t, data, &c)
			if c.Npix != 12*nside*nside {
				t.Fatalf("golden npix %d", c.Npix)
			}
			r2n, n2r := make([]int64, c.Npix), make([]int64, c.Npix)
			for p := range r2n {
				r2n[p] = Ring2Nest(nside, int64(p))
				n2r[p] = Nest2Ring(nside, int64(p))
			}
			c.Ring2Nest.Check(t, "Ring2Nest", r2n)
			c.Nest2Ring.Check(t, "Nest2Ring", n2r)
		})
	}
}

// TestOrderingSampled_Healpy loads ordering_sampled_nside{64,256,1024,8192}
// (ring -> nest pairs from healpy.ring2nest) and checks both directions.
func TestOrderingSampled_Healpy(t *testing.T) {
	for _, nside := range []int{64, 256, 1024, 8192} {
		t.Run(nsideName(nside), func(t *testing.T) {
			var c struct {
				Nside int `json:"nside"`
				Cases []struct {
					Ring int64 `json:"ring"`
					Nest int64 `json:"nest"`
				} `json:"cases"`
			}
			mustUnmarshal(t, loadGolden(t, "ordering_sampled_nside%d.json", nside), &c)
			if c.Nside != nside || len(c.Cases) < 100 {
				t.Fatalf("golden nside=%d with %d cases; want nside=%d and >= 100 cases", c.Nside, len(c.Cases), nside)
			}
			for _, row := range c.Cases {
				if got := Ring2Nest(nside, row.Ring); got != row.Nest {
					t.Errorf("Ring2Nest(%d, %d) = %d, want %d", nside, row.Ring, got, row.Nest)
				}
				if got := Nest2Ring(nside, row.Nest); got != row.Ring {
					t.Errorf("Nest2Ring(%d, %d) = %d, want %d", nside, row.Nest, got, row.Ring)
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
				if !angEq(gotTheta, row.RingTheta) || !angEq(gotPhi, row.RingPhi) {
					t.Errorf("RING pix %d [%d]: got (%.17g, %.17g) want (%.17g, %.17g) (%d, %d ULP)",
						row.Pix, i, gotTheta, gotPhi, row.RingTheta, row.RingPhi,
						ulpDiff(gotTheta, row.RingTheta), ulpDiff(gotPhi, row.RingPhi))
				}
				gotTheta, gotPhi = Pix2Ang(nside, row.Pix, Nested)
				if !angEq(gotTheta, row.NestTheta) || !angEq(gotPhi, row.NestPhi) {
					t.Errorf("NESTED pix %d [%d]: got (%.17g, %.17g) want (%.17g, %.17g) (%d, %d ULP)",
						row.Pix, i, gotTheta, gotPhi, row.NestTheta, row.NestPhi,
						ulpDiff(gotTheta, row.NestTheta), ulpDiff(gotPhi, row.NestPhi))
				}
			}
		})
	}
}

// ---------- Utilities validation ----------

type utilCase struct {
	Cases []struct {
		Nside         int     `json:"nside"`
		Npix          int64   `json:"npix"`
		PixelAreaSr   float64 `json:"pixel_area_sr"`
		PixelResolRad float64 `json:"pixel_resol_rad"`
		MaxPixRad     float64 `json:"max_pixrad"`
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
		checkUtilities(t, row.Nside, row.PixelAreaSr, row.PixelResolRad, row.MaxPixRad)
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
				// Ordered comparison: slot k must hold healpy's k-th
				// neighbour (SW, W, NW, N, NE, E, SE, S), -1 in place.
				gotRing := GetNeighbors(nside, row.Pix, Ring)
				if !neighborsEq(gotRing, row.NeighborsRing) {
					t.Errorf("RING pix %d: got %v want %v", row.Pix, gotRing, row.NeighborsRing)
				}
				gotNest := GetNeighbors(nside, row.Pix, Nested)
				if !neighborsEq(gotNest, row.NeighborsNest) {
					t.Errorf("NESTED pix %d: got %v want %v", row.Pix, gotNest, row.NeighborsNest)
				}
			}
		})
	}
}

// neighborsEq compares an ordered neighbour array with healpy's
// get_all_neighbours output slot by slot (including the -1 slots).
func neighborsEq(got [8]int64, want []int64) bool {
	if len(want) != 8 {
		return false
	}
	for k := range got {
		if got[k] != want[k] {
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
		name = filepath.Join(goldenDir, fmt.Sprintf(pattern, nside))
	} else {
		name = filepath.Join(goldenDir, pattern)
	}
	data, err := testutil.ReadGolden(name)
	if err != nil {
		t.Fatalf("committed golden missing (this is a failure, not a skip): %s: %v (regenerate with gen_golden.py, see testdata/PROVENANCE.md)", name, err)
	}
	return data
}

func mustUnmarshal(t *testing.T, data []byte, v any) {
	t.Helper()
	if err := json.Unmarshal(data, v); err != nil {
		t.Fatalf("unmarshal: %v", err)
	}
}

func nsideName(n int) string {
	return fmt.Sprintf("nside%d", n)
}

// angleULPBudget is the largest difference, in units in the last place,
// allowed between a Go angle (theta or phi, radians) and healpy's.
//
// Both sides evaluate the same healpix_cxx formulas in IEEE double, so the
// only source of difference is the final libm call: healpy's acos/atan2
// come from glibc (correctly rounded in practice) while Go's math.Acos and
// math.Atan2 are accurate to about 1 ULP. Measured maximum over every
// committed golden is reported by TestAngleULPBudget_Measured; the budget
// is set to that maximum, so any extra error fails.
const angleULPBudget = 1

// vecAbsBudget is the absolute tolerance for unit-vector components. A
// component that is mathematically ~0 (an equatorial z, or x at phi=pi/2)
// has no meaningful ULP scale, so components are compared absolutely in
// units of the double epsilon at 1.0 (2^-52).
const vecAbsBudget = 2 * 0x1p-52

// ulpDiff returns the number of representable doubles between a and b
// (0 when equal; math.MaxUint64 when exactly one is NaN or signs differ
// and neither is zero).
func ulpDiff(a, b float64) uint64 {
	if a == b {
		return 0
	}
	if math.IsNaN(a) || math.IsNaN(b) {
		return math.MaxUint64
	}
	ia, ib := orderedBits(a), orderedBits(b)
	if ia > ib {
		return uint64(ia - ib)
	}
	return uint64(ib - ia)
}

// orderedBits maps a double onto a monotonic int64 line so that
// subtracting two of them counts the doubles between.
func orderedBits(f float64) int64 {
	b := int64(math.Float64bits(f))
	if b < 0 {
		return math.MinInt64 - b
	}
	return b
}

func angEq(got, want float64) bool { return ulpDiff(got, want) <= angleULPBudget }

func vecEq(got, want float64) bool { return math.Abs(got-want) <= vecAbsBudget }
