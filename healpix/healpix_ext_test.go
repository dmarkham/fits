package healpix

import (
	"fmt"
	"math"
	"math/rand"
	"strings"
	"testing"
	"time"
)

// Tests in this file check claims that previously had no golden:
// MaxPixRad/PixelResol, IsNsidePow2/Npix2Nside, Vec2Pix/Pix2Vec, ordered
// neighbours at face edges, non-power-of-2 RING nside, nside 2^20 and
// 2^29, and the O(1) RING decomposition. Every expected value comes from
// healpy 1.20.1 via testdata/gen_golden.py (see testdata/PROVENANCE.md).

// maxPixRadAbsBudget: healpy computes max_pixrad as atan2(|A×B|, A·B)
// from rounded Cartesian components of two nearly coincident unit vectors,
// so its result carries an absolute (not relative) error of order 2^-53 —
// e.g. Go's and glibc's cos(π/128) differ in the last bit, which moves
// healpy's formula by 27 ULP of the 0.0167 result at Nside 64. Go's
// MaxPixRad is checked to 2 ULP against the exact (mpmath) value instead;
// against healpy only the absolute 2^-52 budget is meaningful.
const maxPixRadAbsBudget = 0x1p-52

// checkUtilities asserts one utilities golden row.
func checkUtilities(t *testing.T, nside int, area, resol, maxRad float64) {
	t.Helper()
	// healpy: nside2pixarea = 4*pi/npix and nside2resol = sqrt(that); the
	// Go expressions are the same IEEE operations, so they must be equal.
	if got := PixelArea(nside); got != area {
		t.Errorf("PixelArea(%d) = %.17g, want %.17g", nside, got, area)
	}
	if got := PixelResol(nside); got != resol {
		t.Errorf("PixelResol(%d) = %.17g, want %.17g", nside, got, resol)
	}
	if got := MaxPixRad(nside); math.Abs(got-maxRad) > maxPixRadAbsBudget {
		t.Errorf("MaxPixRad(%d) = %.17g, healpy %.17g (|diff| %.3g > %.3g, rel err %.2e)",
			nside, got, maxRad, math.Abs(got-maxRad), maxPixRadAbsBudget, math.Abs(got-maxRad)/maxRad)
	}
}

func TestUtilitiesExt_Healpy(t *testing.T) {
	var g struct {
		Cases []struct {
			Nside          int     `json:"nside"`
			Npix           int64   `json:"npix"`
			PixelAreaSr    float64 `json:"pixel_area_sr"`
			PixelResolRad  float64 `json:"pixel_resol_rad"`
			MaxPixRad      float64 `json:"max_pixrad"`
			MaxPixRadExact float64 `json:"max_pixrad_exact"`
		} `json:"cases"`
		IsNsideOK []struct {
			Nside  int  `json:"nside"`
			NestOK bool `json:"nest_ok"`
			RingOK bool `json:"ring_ok"`
		} `json:"isnsideok"`
		Npix2Nside []struct {
			Npix  int64 `json:"npix"`
			Valid bool  `json:"valid"`
			Nside int   `json:"nside"`
		} `json:"npix2nside"`
	}
	mustUnmarshal(t, loadGolden(t, "utilities_ext.json", 0), &g)
	if len(g.Cases) < 30 || len(g.IsNsideOK) < 100 || len(g.Npix2Nside) < 10 {
		t.Fatalf("utilities_ext.json too small: %d/%d/%d", len(g.Cases), len(g.IsNsideOK), len(g.Npix2Nside))
	}
	for _, c := range g.Cases {
		if got := Nside2Npix(c.Nside); got != c.Npix {
			t.Errorf("Nside2Npix(%d) = %d, want %d", c.Nside, got, c.Npix)
		}
		checkUtilities(t, c.Nside, c.PixelAreaSr, c.PixelResolRad, c.MaxPixRad)
		if c.MaxPixRadExact == 0 {
			t.Fatalf("nside %d: golden lacks max_pixrad_exact", c.Nside)
		}
		// 2 ULP: the measured maximum over all 38 Nside values (sin, sqrt
		// and asin each round once; nothing cancels).
		if got := MaxPixRad(c.Nside); ulpDiff(got, c.MaxPixRadExact) > 2 {
			t.Errorf("MaxPixRad(%d) = %.17g, exact %.17g (%d ULP > 2)", c.Nside, got, c.MaxPixRadExact, ulpDiff(got, c.MaxPixRadExact))
		}
	}
	for _, c := range g.IsNsideOK {
		// healpy isnsideok(nest=True) is "power of 2 and <= 2^29";
		// IsNsidePow2 is the power-of-2 part, NsideOK the whole check.
		pow2 := c.Nside > 0 && c.Nside&(c.Nside-1) == 0
		if got := IsNsidePow2(c.Nside); got != pow2 {
			t.Errorf("IsNsidePow2(%d) = %v, want %v", c.Nside, got, pow2)
		}
		if got := NsideOK(c.Nside, Nested); got != c.NestOK {
			t.Errorf("NsideOK(%d, Nested) = %v, healpy isnsideok(nest=True) = %v", c.Nside, got, c.NestOK)
		}
		if got := NsideOK(c.Nside, Ring); got != c.RingOK {
			t.Errorf("NsideOK(%d, Ring) = %v, healpy isnsideok(nest=False) = %v", c.Nside, got, c.RingOK)
		}
	}
	// Two healpy answers are wrong, and Npix2Nside deliberately differs:
	//   - isnpixok(0) is True (it accepts nside 0); 0 pixels is not a map.
	//   - isnpixok(12*2^58 + 12) is True with nside 2^29 because healpy
	//     converts npix to double, which drops the +12; 2^58+1 is not a
	//     square, so the count is invalid.
	//   - isnpixok(12*(2^29+1)^2) is True with nside 2^29+1: healpy has no
	//     upper limit, but nside > 2^29 is outside what healpix_cxx (and
	//     this package, MaxNside) can index.
	healpyWrong := map[int64]bool{0: true, 12<<58 + 12: true, 12 * (1<<29 + 1) * (1<<29 + 1): true}
	for npix := range healpyWrong {
		found := false
		for _, c := range g.Npix2Nside {
			found = found || c.Npix == npix
		}
		if !found {
			t.Errorf("npix %d listed in healpyWrong but missing from utilities_ext.json", npix)
		}
	}
	// The 2^29 limit itself: the largest valid count is accepted.
	if ns, ok := Npix2Nside(12 << 58); !ok || ns != MaxNside {
		t.Errorf("Npix2Nside(12*2^58) = (%d, %v), want (2^29, true)", ns, ok)
	}
	for _, c := range g.Npix2Nside {
		ns, ok := Npix2Nside(c.Npix)
		if healpyWrong[c.Npix] {
			if !c.Valid {
				t.Errorf("npix %d: healpy now rejects it; drop it from healpyWrong", c.Npix)
			}
			if ok {
				t.Errorf("Npix2Nside(%d) = (%d, true), want (0, false)", c.Npix, ns)
			}
			continue
		}
		if ok != c.Valid || (ok && ns != c.Nside) {
			t.Errorf("Npix2Nside(%d) = (%d, %v), healpy = (%d, %v)", c.Npix, ns, ok, c.Nside, c.Valid)
		}
	}
}

func TestAng2PixRandom_Healpy(t *testing.T) {
	var g struct {
		Theta   f64s `json:"theta"`
		Phi     f64s `json:"phi"`
		Results []struct {
			Nside   int  `json:"nside"`
			PixRing i64s `json:"pix_ring"`
			PixNest i64s `json:"pix_nest"`
		} `json:"results"`
	}
	mustUnmarshal(t, loadGolden(t, "ang2pix_random.json", 0), &g)
	if len(g.Theta) < 2000 || len(g.Phi) != len(g.Theta) || len(g.Results) < 21 {
		t.Fatalf("ang2pix_random.json too small: %d directions, %d nsides", len(g.Theta), len(g.Results))
	}
	for _, r := range g.Results {
		t.Run(nsideName(r.Nside), func(t *testing.T) {
			for _, sc := range []struct {
				name   string
				scheme Scheme
				want   i64s
			}{{"RING", Ring, r.PixRing}, {"NESTED", Nested, r.PixNest}} {
				if sc.scheme == Nested && !sc.want.Present() {
					continue
				}
				got := make([]int64, len(g.Theta))
				for i := range got {
					got[i] = Ang2Pix(r.Nside, g.Theta[i], g.Phi[i], sc.scheme)
				}
				sc.want.Check(t, sc.name+" Ang2Pix", got)
			}
		})
	}
}

// TestAng2PixPhiWrap: phi outside [0, 2π) (negative, 2π, large) is
// reduced exactly as healpy does.
func TestAng2PixPhiWrap(t *testing.T) {
	var g struct {
		Theta   []float64 `json:"theta"`
		Phi     []float64 `json:"phi"`
		Results []struct {
			Nside   int     `json:"nside"`
			PixRing []int64 `json:"pix_ring"`
			PixNest []int64 `json:"pix_nest"`
		} `json:"results"`
	}
	mustUnmarshal(t, loadGolden(t, "ang2pix_phiwrap.json", 0), &g)
	if len(g.Theta) < 60 || len(g.Results) < 6 {
		t.Fatalf("ang2pix_phiwrap.json too small")
	}
	for _, r := range g.Results {
		for i := range g.Theta {
			if got := Ang2Pix(r.Nside, g.Theta[i], g.Phi[i], Ring); got != r.PixRing[i] {
				t.Errorf("nside %d RING theta=%g phi=%g: got %d want %d", r.Nside, g.Theta[i], g.Phi[i], got, r.PixRing[i])
			}
			if r.PixNest != nil {
				if got := Ang2Pix(r.Nside, g.Theta[i], g.Phi[i], Nested); got != r.PixNest[i] {
					t.Errorf("nside %d NESTED theta=%g phi=%g: got %d want %d", r.Nside, g.Theta[i], g.Phi[i], got, r.PixNest[i])
				}
			}
		}
	}
}

func TestVec2Pix_Healpy(t *testing.T) {
	var g struct {
		VX      f64s `json:"vx"`
		VY      f64s `json:"vy"`
		VZ      f64s `json:"vz"`
		Results []struct {
			Nside   int  `json:"nside"`
			PixRing i64s `json:"pix_ring"`
			PixNest i64s `json:"pix_nest"`
		} `json:"results"`
	}
	mustUnmarshal(t, loadGolden(t, "vec2pix.json", 0), &g)
	if len(g.VX) < 1000 || len(g.VY) != len(g.VX) || len(g.VZ) != len(g.VX) || len(g.Results) < 21 {
		t.Fatalf("vec2pix.json too small: %d vectors, %d nsides", len(g.VX), len(g.Results))
	}
	for _, r := range g.Results {
		t.Run(nsideName(r.Nside), func(t *testing.T) {
			for _, sc := range []struct {
				name   string
				scheme Scheme
				want   i64s
			}{{"RING", Ring, r.PixRing}, {"NESTED", Nested, r.PixNest}} {
				if sc.scheme == Nested && !sc.want.Present() {
					continue
				}
				got := make([]int64, len(g.VX))
				for i := range got {
					got[i] = Vec2Pix(r.Nside, g.VX[i], g.VY[i], g.VZ[i], sc.scheme)
				}
				sc.want.Check(t, sc.name+" Vec2Pix", got)
			}
		})
	}
}

func TestPix2Vec_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16, 64, 256, 1024, 8192} {
		t.Run(nsideName(nside), func(t *testing.T) {
			var g struct {
				Pix  []int64 `json:"pix"`
				Ring [3]f64s `json:"ring"`
				Nest [3]f64s `json:"nest"`
			}
			mustUnmarshal(t, loadGolden(t, "pix2vec_nside%d.json", nside), &g)
			for i, p := range g.Pix {
				x, y, z := Pix2Vec(nside, p, Ring)
				if !vecEq(x, g.Ring[0][i]) || !vecEq(y, g.Ring[1][i]) || !vecEq(z, g.Ring[2][i]) {
					t.Errorf("RING pix %d: got (%.17g,%.17g,%.17g) want (%.17g,%.17g,%.17g)",
						p, x, y, z, g.Ring[0][i], g.Ring[1][i], g.Ring[2][i])
				}
				x, y, z = Pix2Vec(nside, p, Nested)
				if !vecEq(x, g.Nest[0][i]) || !vecEq(y, g.Nest[1][i]) || !vecEq(z, g.Nest[2][i]) {
					t.Errorf("NESTED pix %d: got (%.17g,%.17g,%.17g) want (%.17g,%.17g,%.17g)",
						p, x, y, z, g.Nest[0][i], g.Nest[1][i], g.Nest[2][i])
				}
			}
		})
	}
}

type pixfuncScheme struct {
	Theta      f64s `json:"theta"`
	Phi        f64s `json:"phi"`
	VX         f64s `json:"vx"`
	VY         f64s `json:"vy"`
	VZ         f64s `json:"vz"`
	X          i64s `json:"x"`
	Y          i64s `json:"y"`
	F          i64s `json:"f"`
	Neighbours i64s `json:"neighbours"` // 8 per pixel, SW,W,NW,N,NE,E,SE,S
}

// TestPixFuncs_Healpy covers non-power-of-2 RING nside (exhaustive for the
// small ones) and nside 2^20 and 2^29 in both orderings: pix2ang, pix2vec,
// pix2xyf/xyf2pix, ring2nest/nest2ring and ordered neighbours.
func TestPixFuncs_Healpy(t *testing.T) {
	for _, nside := range []int{3, 5, 6, 7, 12, 100, 1000, 59049, 3 << 27, 1<<29 - 1, 1 << 20, 1 << 29} {
		t.Run(nsideName(nside), func(t *testing.T) {
			var g struct {
				Nside     int            `json:"nside"`
				Pix       []int64        `json:"pix"`
				Ring      pixfuncScheme  `json:"ring"`
				Nest      *pixfuncScheme `json:"nest"`
				Ring2Nest i64s           `json:"ring2nest"`
				Nest2Ring i64s           `json:"nest2ring"`
			}
			mustUnmarshal(t, loadGolden(t, "pixfuncs_nside%d.json", nside), &g)
			if g.Nside != nside || len(g.Pix) == 0 {
				t.Fatalf("bad golden: nside %d, %d pixels", g.Nside, len(g.Pix))
			}
			if IsNsidePow2(nside) != (g.Nest != nil) {
				t.Fatalf("golden nest section presence %v disagrees with IsNsidePow2", g.Nest != nil)
			}
			n := len(g.Pix)
			check := func(name string, sc Scheme, w *pixfuncScheme) {
				if len(w.Theta) != n || len(w.Phi) != n || len(w.VX) != n || len(w.VY) != n || len(w.VZ) != n {
					t.Fatalf("%s: golden float arrays do not have %d values", name, n)
				}
				bad := 0
				fail := func(format string, args ...any) {
					if bad++; bad <= 8 {
						t.Errorf(name+" "+format, args...)
					}
				}
				gx, gy, gf := make([]int64, n), make([]int64, n), make([]int64, n)
				nbs := make([]int64, 0, 8*n)
				for i, p := range g.Pix {
					th, ph := Pix2Ang(nside, p, sc)
					if !angEq(th, w.Theta[i]) || !angEq(ph, w.Phi[i]) {
						fail("pix %d: Pix2Ang (%.17g, %.17g) want (%.17g, %.17g) (%d, %d ULP)",
							p, th, ph, w.Theta[i], w.Phi[i], ulpDiff(th, w.Theta[i]), ulpDiff(ph, w.Phi[i]))
					}
					x, y, z := Pix2Vec(nside, p, sc)
					if !vecEq(x, w.VX[i]) || !vecEq(y, w.VY[i]) || !vecEq(z, w.VZ[i]) {
						fail("pix %d: Pix2Vec (%.17g,%.17g,%.17g) want (%.17g,%.17g,%.17g)",
							p, x, y, z, w.VX[i], w.VY[i], w.VZ[i])
					}
					ix, iy, f := pix2xyf(nside, p, sc)
					gx[i], gy[i], gf[i] = int64(ix), int64(iy), int64(f)
					nb := GetNeighbors(nside, p, sc)
					nbs = append(nbs, nb[:]...)
				}
				if bad > 8 {
					t.Errorf("%s: %d mismatches in total", name, bad)
				}
				wx := w.X.Resolve(t, name+" pix2xyf x", gx)
				wy := w.Y.Resolve(t, name+" pix2xyf y", gy)
				wf := w.F.Resolve(t, name+" pix2xyf face", gf)
				for i, p := range g.Pix {
					if back := xyf2pix(nside, int(wx[i]), int(wy[i]), int(wf[i]), sc); back != p {
						t.Errorf("%s xyf2pix(%d,%d,%d) = %d want %d", name, wx[i], wy[i], wf[i], back, p)
					}
				}
				w.Neighbours.Check(t, name+" GetNeighbors (8 per pixel)", nbs)
			}
			check("RING", Ring, &g.Ring)
			if g.Nest != nil {
				check("NESTED", Nested, g.Nest)
				r2n, n2r := make([]int64, n), make([]int64, n)
				for i, p := range g.Pix {
					r2n[i], n2r[i] = Ring2Nest(nside, p), Nest2Ring(nside, p)
				}
				g.Ring2Nest.Check(t, "Ring2Nest", r2n)
				g.Nest2Ring.Check(t, "Nest2Ring", n2r)
			}
		})
	}
}

// TestNeighboursOrdered_Healpy checks every face-edge and face-corner
// pixel (where neighbours cross faces) against healpy, slot by slot.
func TestNeighboursOrdered_Healpy(t *testing.T) {
	for _, nside := range []int{1, 2, 4, 8, 16, 64, 8192} {
		t.Run(nsideName(nside), func(t *testing.T) {
			var g struct {
				Pix  []int64 `json:"pix"`
				Ring i64s    `json:"ring"` // 8 per pixel
				Nest i64s    `json:"nest"`
			}
			mustUnmarshal(t, loadGolden(t, "neighbours_ordered_nside%d.json", nside), &g)
			if len(g.Pix) == 0 {
				t.Fatalf("no pixels")
			}
			for _, sc := range []struct {
				name   string
				scheme Scheme
				want   i64s
			}{{"RING", Ring, g.Ring}, {"NESTED", Nested, g.Nest}} {
				got := make([]int64, 0, 8*len(g.Pix))
				for _, p := range g.Pix {
					nb := GetNeighbors(nside, p, sc.scheme)
					got = append(got, nb[:]...)
				}
				sc.want.Check(t, sc.name+" GetNeighbors (8 per pixel)", got)
			}
		})
	}
}

// TestAngleULPBudget_Measured makes the angle tolerance exactly as tight
// as the truth allows: the largest ULP difference over every committed
// pix2ang golden must equal angleULPBudget (a larger value fails the
// per-pixel tests; a smaller one means the budget is looser than needed).
func TestAngleULPBudget_Measured(t *testing.T) {
	var maxULP uint64
	var where string
	note := func(got, want float64, ctx string) {
		if d := ulpDiff(got, want); d > maxULP {
			maxULP, where = d, ctx
		}
	}
	for _, nside := range []int{1, 2, 4, 8, 16, 64, 256, 1024, 8192} {
		var c pix2angCase
		mustUnmarshal(t, loadGolden(t, "pix2ang_nside%d.json", nside), &c)
		for _, r := range c.Cases {
			th, ph := Pix2Ang(nside, r.Pix, Ring)
			note(th, r.RingTheta, fmt.Sprintf("nside %d RING pix %d theta", nside, r.Pix))
			note(ph, r.RingPhi, fmt.Sprintf("nside %d RING pix %d phi", nside, r.Pix))
			th, ph = Pix2Ang(nside, r.Pix, Nested)
			note(th, r.NestTheta, fmt.Sprintf("nside %d NESTED pix %d theta", nside, r.Pix))
			note(ph, r.NestPhi, fmt.Sprintf("nside %d NESTED pix %d phi", nside, r.Pix))
		}
	}
	t.Logf("max angle difference vs healpy: %d ULP (%s)", maxULP, where)
	if maxULP != angleULPBudget {
		t.Errorf("measured max %d ULP != angleULPBudget %d: set the budget to the measured value", maxULP, angleULPBudget)
	}
}

// TestRingDecomposition_ConstantTime guards the O(1) RING pixel
// decomposition: 200000 random RING pixels at nside 2^29 must convert in
// well under 2 s. An O(nside) ring walk needs ~10^8 iterations per
// equatorial or south-cap pixel, i.e. days for this loop, so the bound
// has a margin of many orders of magnitude and is not timing-flaky.
func TestRingDecomposition_ConstantTime(t *testing.T) {
	const nside = 1 << 29
	npix := Nside2Npix(nside)
	r := rand.New(rand.NewSource(1))
	deadline := time.Now().Add(2 * time.Second)
	const n = 200000
	for i := 0; i < n; i++ {
		p := r.Int63n(npix)
		x, y, f := pix2xyf(nside, p, Ring)
		if back := xyf2pix(nside, x, y, f, Ring); back != p {
			t.Fatalf("RING pix %d -> xyf (%d,%d,%d) -> %d", p, x, y, f, back)
		}
		if i%1000 == 0 && time.Now().After(deadline) {
			t.Fatalf("only %d of %d RING decompositions at nside 2^29 finished in 2 s: decomposition is not O(1)", i, n)
		}
	}
}

// TestInvalidInputPanics: invalid input is rejected loudly, by the
// package's own validation (the panic message is checked, so a Go runtime
// panic such as a negative shift does not count), rather than returning a
// wrong pixel.
func TestInvalidInputPanics(t *testing.T) {
	cases := []struct {
		name string
		f    func()
		want string
	}{
		{"Ang2Pix NESTED nside=3", func() { Ang2Pix(3, 1, 1, Nested) }, "NESTED ordering requires a power-of-2 nside, got 3"},
		{"Pix2Ang NESTED nside=3", func() { Pix2Ang(3, 0, Nested) }, "NESTED ordering requires a power-of-2 nside"},
		{"Pix2Vec NESTED nside=6", func() { Pix2Vec(6, 0, Nested) }, "NESTED ordering requires a power-of-2 nside"},
		{"Vec2Pix NESTED nside=12", func() { Vec2Pix(12, 1, 0, 0, Nested) }, "NESTED ordering requires a power-of-2 nside"},
		{"Ring2Nest nside=3", func() { Ring2Nest(3, 0) }, "NESTED ordering requires a power-of-2 nside"},
		{"Nest2Ring nside=3", func() { Nest2Ring(3, 0) }, "NESTED ordering requires a power-of-2 nside"},
		{"GetNeighbors NESTED nside=3", func() { GetNeighbors(3, 0, Nested) }, "NESTED ordering requires a power-of-2 nside"},
		{"nside 0", func() { Ang2Pix(0, 1, 1, Ring) }, "nside 0 out of range"},
		{"nside -4", func() { Pix2Ang(-4, 0, Ring) }, "nside -4 out of range"},
		{"nside 2^29+1", func() { Pix2Ang(1<<29+1, 0, Ring) }, "out of range [1, 2^29]"},
		{"nside 2^30", func() { Ang2Pix(1<<30, 1, 1, Ring) }, "out of range [1, 2^29]"},
		{"unknown scheme", func() { Ang2Pix(4, 1, 1, Scheme(7)) }, "unknown scheme 7"},
		{"pix -1", func() { Pix2Ang(4, -1, Ring) }, "pixel -1 out of range"},
		{"pix npix", func() { Pix2Ang(4, 192, Nested) }, "pixel 192 out of range"},
		{"neighbours pix npix", func() { GetNeighbors(4, 192, Ring) }, "pixel 192 out of range"},
		{"theta < 0", func() { Ang2Pix(4, -0.1, 1, Ring) }, "outside [0, pi]"},
		{"theta > pi", func() { Ang2Pix(4, math.Pi+1e-9, 1, Ring) }, "outside [0, pi]"},
		{"theta NaN", func() { Ang2Pix(4, math.NaN(), 1, Ring) }, "outside [0, pi]"},
		{"phi NaN", func() { Ang2Pix(4, 1, math.NaN(), Nested) }, "is not finite"},
		{"phi +Inf", func() { Ang2Pix(4, 1, math.Inf(1), Ring) }, "is not finite"},
		{"phi -Inf", func() { Ang2Pix(4, 1, math.Inf(-1), Nested) }, "is not finite"},
		{"zero vector", func() { Vec2Pix(4, 0, 0, 0, Ring) }, "zero vector"},
		{"vector +Inf x", func() { Vec2Pix(4, math.Inf(1), 0, 0, Ring) }, "is not finite"},
		{"vector -Inf z", func() { Vec2Pix(4, 0, 0, math.Inf(-1), Nested) }, "is not finite"},
		{"vector NaN y", func() { Vec2Pix(4, 1, math.NaN(), 0, Ring) }, "is not finite"},
		{"vector NaN only", func() { Vec2Pix(4, 0, 0, math.NaN(), Ring) }, "is not finite"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			defer func() {
				r := recover()
				if r == nil {
					t.Fatalf("%s: no panic", c.name)
				}
				msg := fmt.Sprint(r)
				if !strings.HasPrefix(msg, "healpix.") || !strings.Contains(msg, c.want) {
					t.Errorf("%s: panic %q, want the package's message containing %q", c.name, msg, c.want)
				}
			}()
			c.f()
		})
	}
}

// TestVec2PixExtremeLength: vectors whose squared length overflows or
// underflows. Truth: healpy.vec2pix of the same direction at unit scale
// (the golden vectors are unit-scale vectors times 2^k, exactly). The
// golden also records healpy's answer for the scaled vector, which is
// wrong when the length underflows; the test checks that such cases exist
// so the scaling is really exercised.
func TestVec2PixExtremeLength(t *testing.T) {
	var g struct {
		Cases []struct {
			Nside      int        `json:"nside"`
			V          [3]float64 `json:"v"`
			Ring       int64      `json:"ring"`
			Nest       int64      `json:"nest"`
			HealpyRing int64      `json:"healpy_ring"`
			HealpyNest int64      `json:"healpy_nest"`
		} `json:"cases"`
	}
	mustUnmarshal(t, loadGolden(t, "vec2pix_extreme.json", 0), &g)
	healpyWrong := 0
	for _, c := range g.Cases {
		if c.HealpyRing != c.Ring || c.HealpyNest != c.Nest {
			healpyWrong++
		}
		if got := Vec2Pix(c.Nside, c.V[0], c.V[1], c.V[2], Ring); got != c.Ring {
			t.Errorf("Vec2Pix(%d, %v, RING) = %d, want %d", c.Nside, c.V, got, c.Ring)
		}
		if got := Vec2Pix(c.Nside, c.V[0], c.V[1], c.V[2], Nested); got != c.Nest {
			t.Errorf("Vec2Pix(%d, %v, NESTED) = %d, want %d", c.Nside, c.V, got, c.Nest)
		}
	}
	if len(g.Cases) < 80 || healpyWrong == 0 {
		t.Fatalf("vec2pix_extreme.json: %d cases, %d where healpy fails; expected both to be non-trivial", len(g.Cases), healpyWrong)
	}
	// The scales just inside and outside 2^±511 (where the squares start
	// to overflow or lose bits) must be present for both signs.
	seen := map[int]bool{}
	for _, c := range g.Cases {
		_, e := math.Frexp(math.Max(math.Abs(c.V[0]), math.Max(math.Abs(c.V[1]), math.Abs(c.V[2]))))
		seen[e] = true
	}
	for _, band := range [][2]int{{502, 515}, {516, 560}, {561, 600}, {-600, -561}, {-560, -516}, {-515, -500}} {
		ok := false
		for e := band[0]; e <= band[1]; e++ {
			ok = ok || seen[e]
		}
		if !ok {
			t.Errorf("vec2pix_extreme.json: no vector with largest component in 2^[%d, %d]", band[0]-1, band[1])
		}
	}
}

// TestVec2PixScaledPolar: Vec2Pix depends only on the direction. Every
// polar-boundary vector of ang2pix_polar.json is scaled by 2^k, which is
// exact, for k from -1060 to 1020, and must land in exactly the pixel of
// the unscaled vector (the rescaling makes the two computations
// identical). The unscaled pixels are checked against the golden
// (healpix_cxx with correctly rounded atan2, equal to healpy at unit
// scale) as in TestAng2PixPolar_Healpy, boundary-ambiguous differences
// allowed and counted. healpy itself gives a different pixel for hundreds of them
// already at k = -490..-500 (squares of the small x, y components lose
// bits in the subnormal range while the largest component is still above
// 2^-500), and Vec2Pix did too before it rescaled every vector.
func TestVec2PixScaledPolar(t *testing.T) {
	g := loadPolar(t)
	n := len(g.VX)
	type want struct {
		nside  int
		scheme Scheme
		pix    []int64
	}
	var wants []want
	var tally pixTally
	for _, r := range g.Results {
		for _, sc := range []struct {
			scheme Scheme
			vec    i64s
		}{{Ring, r.VecRing}, {Nested, r.VecNest}} {
			if sc.scheme == Nested && !sc.vec.Present() {
				continue
			}
			ref := explicitValues(t, fmt.Sprintf("nside %d unscaled Vec2Pix", r.Nside), sc.vec)
			unscaled := make([]int64, n)
			for i := range unscaled {
				x, y, z := g.VX[i], g.VY[i], g.VZ[i]
				unscaled[i] = Vec2Pix(r.Nside, x, y, z, sc.scheme)
				tally.check(t, unscaled[i], ref[i],
					func() map[int64]bool { return vec2pixAlts(r.Nside, x, y, z, sc.scheme) },
					func() string {
						return fmt.Sprintf("Vec2Pix(%d, %.17g, %.17g, %.17g, %v) [%d]", r.Nside, x, y, z, sc.scheme, i)
					})
			}
			// The expected pixels of the scaled vectors: those of the
			// unscaled ones.
			wants = append(wants, want{r.Nside, sc.scheme, unscaled})
		}
	}
	tally.done(t, "unscaled vectors", polarAmbiguousMax)
	bad, checked := 0, 0
	// Every exponent from -512 to -470, where healpy's own unscaled formula
	// starts to move pixels (2^-482 at Nside 2^29: vec2pix_scaled_healpy),
	// so a Vec2Pix that stops rescaling anywhere below 2^-482 fails here.
	ks := []int{-1060, -1000, -700, -600, -540, -400, -1, 1, 400, 480, 490, 495, 499, 500, 505, 510, 512, 540, 600, 1000, 1020}
	for k := -512; k <= -470; k++ {
		ks = append(ks, k)
	}
	for _, k := range ks {
		for i := 0; i < n; i++ {
			v := [3]float64{math.Ldexp(g.VX[i], k), math.Ldexp(g.VY[i], k), math.Ldexp(g.VZ[i], k)}
			if math.Ldexp(v[0], -k) != g.VX[i] || math.Ldexp(v[1], -k) != g.VY[i] || math.Ldexp(v[2], -k) != g.VZ[i] {
				continue // scaling not exact (subnormal): a different direction
			}
			for _, w := range wants {
				checked++
				if got := Vec2Pix(w.nside, v[0], v[1], v[2], w.scheme); got != w.pix[i] {
					if bad++; bad <= 6 {
						t.Errorf("Vec2Pix(%d, 2^%d * (%.17g,%.17g,%.17g), %v) = %d, want %d", w.nside, k, g.VX[i], g.VY[i], g.VZ[i], w.scheme, got, w.pix[i])
					}
				}
			}
		}
	}
	if bad > 6 {
		t.Errorf("... %d mismatches in total", bad)
	}
	if checked < 200000 {
		t.Fatalf("only %d scaled vectors checked", checked)
	}
}

// TestAng2PixPolar_Healpy pins the near-pole sin(theta) branch of Ang2Pix
// (theta < 0.01 or > pi - 0.01) and Vec2Pix (|z| > 0.99), and the
// rounding of cos/sin/atan2, at the places where they decide the pixel:
// directions on, and up to 2 ULP either side of, polar-cap ring and pixel
// boundaries solved in 50-digit arithmetic at Nside 8192, 2^20, 2^29-1
// and 2^29 (gen_golden.py polar_directions), where 1 - cos(theta) has
// cancelled most of its bits.
//
// Expected values are the real healpix_cxx loc2pix (testdata/hpref) fed
// correctly rounded cos/sin/atan2. healpy's own answers are in the golden
// too: the generator asserts the harness reproduces healpy when it calls
// glibc itself, and lists every answer where healpy differs because glibc
// is not correctly rounded (none for the committed inputs, so here the
// expected values are healpy's answers).
//
// Every input lies within 2 ULP of a boundary by construction, so this is
// where Go's math.Cos/Sin/Atan2 (not correctly rounded) can pick the other
// side: a different pixel is accepted only at a boundary-ambiguous point
// (ambiguity_test.go), counted and logged, and on at most
// polarAmbiguousMax (6%) of the pixels per Nside.
func TestAng2PixPolar_Healpy(t *testing.T) {
	g := loadPolar(t)
	n := len(g.Theta)
	var all pixTally
	for _, r := range g.Results {
		t.Run(nsideName(r.Nside), func(t *testing.T) {
			if len(r.DiffsRing)+len(r.DiffsNest) != 0 {
				t.Fatalf("golden lists healpy/glibc rounding differences %v %v; this test assumes none", r.DiffsRing, r.DiffsNest)
			}
			var tally pixTally
			for _, sc := range []struct {
				name     string
				scheme   Scheme
				ang, vec i64s
			}{{"RING", Ring, r.AngRing, r.VecRing}, {"NESTED", Nested, r.AngNest, r.VecNest}} {
				if sc.scheme == Nested && !sc.ang.Present() {
					continue
				}
				wa := explicitValues(t, sc.name+" Ang2Pix", sc.ang)
				wv := explicitValues(t, sc.name+" Vec2Pix", sc.vec)
				if len(wa) != n || len(wv) != n {
					t.Fatalf("%s: %d/%d expected pixels for %d directions", sc.name, len(wa), len(wv), n)
				}
				gg := newGeom("test", r.Nside, sc.scheme)
				for i := 0; i < n; i++ {
					th, ph := g.Theta[i], g.Phi[i]
					x, y, z := g.VX[i], g.VY[i], g.VZ[i]
					// Exact on the correctly rounded libm values the
					// reference was computed from.
					if got := gg.ang2pixCS(th, ph, g.CosCR[i], g.SinCR[i], sc.scheme); got != wa[i] {
						t.Errorf("%s ang2pixCS(%d, %.17g, %.17g, correctly rounded cos/sin) [%d] = %d, want %d", sc.name, r.Nside, th, ph, i, got, wa[i])
					}
					m := math.Max(math.Abs(x), math.Max(math.Abs(y), math.Abs(z)))
					if got := gg.vec2pixPhi(x, y, z, m, g.Atan2CR[i], sc.scheme); got != wv[i] {
						t.Errorf("%s vec2pixPhi(%d, %.17g, %.17g, %.17g, correctly rounded atan2) [%d] = %d, want %d", sc.name, r.Nside, x, y, z, i, got, wv[i])
					}
					tally.check(t, Ang2Pix(r.Nside, th, ph, sc.scheme), wa[i],
						func() map[int64]bool { return ang2pixAlts(r.Nside, th, ph, sc.scheme) },
						func() string { return fmt.Sprintf("%s Ang2Pix(%d, %.17g, %.17g) [%d]", sc.name, r.Nside, th, ph, i) })
					tally.check(t, Vec2Pix(r.Nside, x, y, z, sc.scheme), wv[i],
						func() map[int64]bool { return vec2pixAlts(r.Nside, x, y, z, sc.scheme) },
						func() string {
							return fmt.Sprintf("%s Vec2Pix(%d, %.17g, %.17g, %.17g) [%d]", sc.name, r.Nside, x, y, z, i)
						})
				}
			}
			tally.done(t, nsideName(r.Nside), polarAmbiguousMax)
			all.checked += tally.checked
			all.ambiguous += tally.ambiguous
		})
	}
	all.done(t, "all Nside", polarAmbiguousMax)
}

// polarGolden is ang2pix_polar.json (gen_golden.py write_ang2pix_polar).
// CosCR, SinCR and Atan2CR are the correctly rounded cos(theta),
// sin(theta) and atan2(vy, vx) (mpmath) the expected pixels were computed
// from.
type polarGolden struct {
	Theta   f64s `json:"theta"`
	Phi     f64s `json:"phi"`
	VX      f64s `json:"vx"`
	VY      f64s `json:"vy"`
	VZ      f64s `json:"vz"`
	CosCR   f64s `json:"cos_cr"`
	SinCR   f64s `json:"sin_cr"`
	Atan2CR f64s `json:"atan2_cr"`
	Results []struct {
		Nside     int     `json:"nside"`
		AngRing   i64s    `json:"ang_ring"`
		AngNest   i64s    `json:"ang_nest"`
		VecRing   i64s    `json:"vec_ring"`
		VecNest   i64s    `json:"vec_nest"`
		DiffsRing [][]any `json:"healpy_differs_ring"`
		DiffsNest [][]any `json:"healpy_differs_nest"`
	} `json:"results"`
}

func loadPolar(t *testing.T) polarGolden {
	t.Helper()
	var g polarGolden
	mustUnmarshal(t, loadGolden(t, "ang2pix_polar.json", 0), &g)
	n := len(g.Theta)
	if n < 3000 || len(g.Phi) != n || len(g.VX) != n || len(g.VY) != n || len(g.VZ) != n || len(g.Results) != 4 ||
		len(g.CosCR) != n || len(g.SinCR) != n || len(g.Atan2CR) != n {
		t.Fatalf("ang2pix_polar.json: %d directions (%d phi, %d/%d/%d vectors, %d/%d/%d cos/sin/atan2), %d nsides",
			n, len(g.Phi), len(g.VX), len(g.VY), len(g.VZ), len(g.CosCR), len(g.SinCR), len(g.Atan2CR), len(g.Results))
	}
	return g
}

// TestBranchSwitch_Healpy pins the points where Ang2Pix switches from the
// cos(theta) to the sin(theta) formula (theta < 0.01, theta > 3.14159-0.01)
// and Vec2Pix does (|z| > 0.99). Every golden input is one where the two
// formulas give different pixels (gen_golden.py write_branch_switch, chosen
// with the forced-branch ops of testdata/hpref), and the inputs fill the
// bands a moved threshold would flip: [0.0099, 0.01), (3.13059, 3.13159],
// (3.13159, pi-0.01], |z| == 0.99 and (0.99, 0.9901]. Expected values are
// healpix_cxx fed correctly rounded cos/sin/atan2, and equal healpy's own
// answers (the generator drops any input where they differ).
//
// The branch choice and the arithmetic after the libm calls are checked
// exactly: ang2pixCS / vec2pixPhi fed the same correctly rounded values
// (golden cos_cr, sin_cr, atan2_cr) must give the expected pixel. Ang2Pix
// and Vec2Pix themselves (Go's math.Cos/Sin/Atan2) may differ from it only
// at boundary-ambiguous points (ambiguity_test.go): every input here sits
// on a pixel boundary by construction, and the cos(theta) formula just
// outside theta = 0.01 is where one ULP of cos moves a ring (measured:
// 12.2%, at most 20% allowed).
func TestBranchSwitch_Healpy(t *testing.T) {
	type rec struct {
		Nside   int       `json:"nside"`
		Band    string    `json:"band"`
		Theta   float64   `json:"theta"`
		Phi     float64   `json:"phi"`
		V       []float64 `json:"v"`
		CosCR   *float64  `json:"cos_cr"`
		SinCR   *float64  `json:"sin_cr"`
		Atan2CR *float64  `json:"atan2_cr"`
		Ring    int64     `json:"ring"`
		Nest    *int64    `json:"nest"`
	}
	var g struct {
		Ang []rec `json:"ang"`
		Vec []rec `json:"vec"`
	}
	mustUnmarshal(t, loadGolden(t, "branch_switch.json", 0), &g)
	region := func(r rec) string {
		if r.V != nil {
			switch z := math.Abs(r.V[2]); {
			case z == 0.99:
				return "|z|==0.99"
			case z > 0.99 && z <= 0.9901:
				return "|z| in (0.99,0.9901]"
			}
			return ""
		}
		switch th := r.Theta; {
		case th >= 0.0099 && th < 0.01:
			return "theta in [0.0099,0.01)"
		case th > 3.14159-0.011 && th <= 3.14159-0.01:
			return "theta in (3.13059,3.13159]"
		case th > 3.14159-0.01 && th <= math.Pi-0.01:
			return "theta in (3.13159,pi-0.01]"
		}
		return ""
	}
	count := map[string]map[int]int{}
	for _, r := range append(append([]rec{}, g.Ang...), g.Vec...) {
		if k := region(r); k != "" {
			if count[k] == nil {
				count[k] = map[int]int{}
			}
			count[k][r.Nside]++
		}
	}
	for _, k := range []string{"theta in [0.0099,0.01)", "theta in (3.13059,3.13159]", "theta in (3.13159,pi-0.01]", "|z|==0.99", "|z| in (0.99,0.9901]"} {
		for _, n := range []int{8192, 1 << 20, 1<<29 - 1, 1 << 29} {
			if count[k][n] < 4 {
				t.Fatalf("branch_switch.json: %d inputs with %s at Nside %d, want >= 4", count[k][n], k, n)
			}
		}
	}
	bad := 0
	fail := func(format string, args ...any) {
		if bad++; bad <= 8 {
			t.Errorf(format, args...)
		}
	}
	var tally pixTally
	for _, r := range g.Ang {
		if r.CosCR == nil || r.SinCR == nil {
			t.Fatalf("branch_switch.json: Ang2Pix record without cos_cr/sin_cr")
		}
		for _, sc := range []struct {
			scheme Scheme
			want   *int64
		}{{Ring, &r.Ring}, {Nested, r.Nest}} {
			if sc.want == nil {
				continue
			}
			gg := newGeom("test", r.Nside, sc.scheme)
			if got := gg.ang2pixCS(r.Theta, r.Phi, *r.CosCR, *r.SinCR, sc.scheme); got != *sc.want {
				fail("ang2pixCS(%d, %.17g, %.17g, correctly rounded cos/sin, %v) [%s] = %d, want %d", r.Nside, r.Theta, r.Phi, sc.scheme, r.Band, got, *sc.want)
			}
			tally.check(t, Ang2Pix(r.Nside, r.Theta, r.Phi, sc.scheme), *sc.want,
				func() map[int64]bool { return ang2pixAlts(r.Nside, r.Theta, r.Phi, sc.scheme) },
				func() string {
					return fmt.Sprintf("Ang2Pix(%d, %.17g, %.17g, %v) [%s]", r.Nside, r.Theta, r.Phi, sc.scheme, r.Band)
				})
		}
	}
	for _, r := range g.Vec {
		if r.Atan2CR == nil {
			t.Fatalf("branch_switch.json: Vec2Pix record without atan2_cr")
		}
		x, y, z := r.V[0], r.V[1], r.V[2]
		for _, sc := range []struct {
			scheme Scheme
			want   *int64
		}{{Ring, &r.Ring}, {Nested, r.Nest}} {
			if sc.want == nil {
				continue
			}
			gg := newGeom("test", r.Nside, sc.scheme)
			m := math.Max(math.Abs(x), math.Max(math.Abs(y), math.Abs(z)))
			if got := gg.vec2pixPhi(x, y, z, m, *r.Atan2CR, sc.scheme); got != *sc.want {
				fail("vec2pixPhi(%d, (%.17g, %.17g, %.17g), correctly rounded atan2, %v) [%s] = %d, want %d", r.Nside, x, y, z, sc.scheme, r.Band, got, *sc.want)
			}
			tally.check(t, Vec2Pix(r.Nside, x, y, z, sc.scheme), *sc.want,
				func() map[int64]bool { return vec2pixAlts(r.Nside, x, y, z, sc.scheme) },
				func() string {
					return fmt.Sprintf("Vec2Pix(%d, (%.17g, %.17g, %.17g), %v) [%s]", r.Nside, x, y, z, sc.scheme, r.Band)
				})
		}
	}
	if bad > 8 {
		t.Errorf("... %d mismatches in total", bad)
	}
	tally.done(t, "branch_switch", 0.20)
}

// TestAng2PixTwoThird_Healpy pins loc2pix's equatorial/polar split at
// z = cos(theta) exactly 2/3: theta = 0.8410686705679303 is the only
// double near acos(2/3) whose (correctly rounded) cos is the double
// nearest 2/3, and there 'za <= twothird' (equatorial formula) and
// 'za < twothird' (polar formula) give different pixels for most phi.
// Also the doubles 1 and 2 ULP either side, and the same around
// acos(-2/3), at phi = k pi/32 and every power-of-2 Nside 2^0..2^29,
// RING and NESTED. Truth: healpy 1.20.1 ang2pix (gen_golden.py
// write_ang2pix_twothird keeps only thetas whose glibc cos is correctly
// rounded).
//
// The split itself is checked exactly: ang2pixCS fed the correctly
// rounded cos (golden cos_cr; 2/3 exactly for 64 inputs) must give
// healpy's pixel. Ang2Pix with Go's math.Cos may differ from it only at
// boundary-ambiguous points (ambiguity_test.go), which these all are by
// construction where cos(theta) is within an ULP or two of 2/3 (measured:
// 2.4%, at most 5% allowed).
func TestAng2PixTwoThird_Healpy(t *testing.T) {
	var g struct {
		Theta   f64s `json:"theta"`
		Phi     f64s `json:"phi"`
		CosCR   f64s `json:"cos_cr"`
		Results []struct {
			Nside int  `json:"nside"`
			Ring  i64s `json:"ring"`
			Nest  i64s `json:"nest"`
		} `json:"results"`
	}
	mustUnmarshal(t, loadGolden(t, "ang2pix_twothird.json", 0), &g)
	n := len(g.Theta)
	if n < 512 || len(g.Phi) != n || len(g.CosCR) != n || len(g.Results) != 30 {
		t.Fatalf("ang2pix_twothird.json: %d theta, %d phi, %d cos, %d nsides", n, len(g.Phi), len(g.CosCR), len(g.Results))
	}
	onBoundary := 0
	for _, c := range g.CosCR {
		if c == twothird {
			onBoundary++
		}
	}
	if onBoundary < 64 {
		t.Fatalf("ang2pix_twothird.json: only %d inputs with correctly rounded cos(theta) == 2/3", onBoundary)
	}
	var tally pixTally
	for _, r := range g.Results {
		for _, sc := range []struct {
			name   string
			scheme Scheme
			want   i64s
		}{{"RING", Ring, r.Ring}, {"NESTED", Nested, r.Nest}} {
			want := explicitValues(t, fmt.Sprintf("nside %d %s", r.Nside, sc.name), sc.want)
			if len(want) != n {
				t.Fatalf("nside %d %s: %d expected pixels for %d inputs", r.Nside, sc.name, len(want), n)
			}
			gg := newGeom("test", r.Nside, sc.scheme)
			bad := 0
			for i := 0; i < n; i++ {
				th, ph := g.Theta[i], g.Phi[i]
				if got := gg.ang2pixCS(th, ph, g.CosCR[i], 0, sc.scheme); got != want[i] {
					if bad++; bad <= 4 {
						t.Errorf("nside %d %s: ang2pixCS(%.17g, %.17g, correctly rounded cos %.17g) = %d, healpy %d", r.Nside, sc.name, th, ph, g.CosCR[i], got, want[i])
					}
				}
				tally.check(t, Ang2Pix(r.Nside, th, ph, sc.scheme), want[i],
					func() map[int64]bool { return ang2pixAlts(r.Nside, th, ph, sc.scheme) },
					func() string { return fmt.Sprintf("Ang2Pix(%d, %.17g, %.17g, %s) [%d]", r.Nside, th, ph, sc.name, i) })
			}
		}
	}
	tally.done(t, "twothird", 0.05)
}

// TestVec2PixScaledPolar_HealpyMoves pins the Vec2Pix doc's statement
// about healpy: healpy.vec2pix on the polar-boundary vectors of
// ang2pix_polar.json scaled by 2^k (exact; every k from -490 to -480,
// and -495, -500) returns the same pixels as unscaled at k = -480 and
// -481, first moves one at k = -482 (Nside 2^29 - 1 RING), and moves tens
// of them at k <= -490, at Nside 2^29 - 1 and 2^29. Truth: healpy 1.20.1, its
// scaled answers recorded where they differ from its unscaled ones
// (gen_golden.py write_vec2pix_scaled_healpy). At every such vector
// Vec2Pix keeps the unscaled pixel, so it differs from healpy there.
func TestVec2PixScaledPolar_HealpyMoves(t *testing.T) {
	var g struct {
		Cases []struct {
			Nside  int       `json:"nside"`
			Nest   bool      `json:"nest"`
			K      int       `json:"k"`
			Exact  int       `json:"exact"`
			Differ [][]int64 `json:"healpy_differs"`
		} `json:"cases"`
	}
	mustUnmarshal(t, loadGolden(t, "vec2pix_scaled_healpy.json", 0), &g)
	p := loadPolar(t)
	if len(g.Cases) != 39 {
		t.Fatalf("vec2pix_scaled_healpy.json: %d cases, want 39", len(g.Cases))
	}
	firstMove := false
	coincide, moved := 0, 0
	for _, c := range g.Cases {
		scheme := Ring
		if c.Nest {
			scheme = Nested
		}
		exact := 0
		for i := range p.VX {
			if math.Ldexp(math.Ldexp(p.VX[i], c.K), -c.K) == p.VX[i] && math.Ldexp(math.Ldexp(p.VY[i], c.K), -c.K) == p.VY[i] &&
				math.Ldexp(math.Ldexp(p.VZ[i], c.K), -c.K) == p.VZ[i] {
				exact++
			}
		}
		if exact != c.Exact {
			t.Errorf("nside %d %v 2^%d: %d vectors scale exactly, golden %d", c.Nside, scheme, c.K, exact, c.Exact)
		}
		switch {
		case c.K >= -481 && len(c.Differ) != 0:
			t.Errorf("nside %d %v 2^%d: healpy moves %d vectors, the doc says none", c.Nside, scheme, c.K, len(c.Differ))
		case c.K == -482 && c.Nside == 1<<29-1 && !c.Nest:
			firstMove = len(c.Differ) > 0
		case c.K <= -490 && len(c.Differ) < 10:
			t.Errorf("nside %d %v 2^%d: healpy moves only %d vectors, the doc says tens", c.Nside, scheme, c.K, len(c.Differ))
		}
		for _, d := range c.Differ {
			i := int(d[0])
			want := Vec2Pix(c.Nside, p.VX[i], p.VY[i], p.VZ[i], scheme)
			got := Vec2Pix(c.Nside, math.Ldexp(p.VX[i], c.K), math.Ldexp(p.VY[i], c.K), math.Ldexp(p.VZ[i], c.K), scheme)
			if got != want {
				t.Errorf("nside %d %v vector %d * 2^%d: Vec2Pix %d, unscaled %d", c.Nside, scheme, i, c.K, got, want)
			}
			// healpy's scaled answer is a wrong pixel for the direction,
			// except that a boundary-ambiguous vector (which these
			// boundary vectors can be) has two right ones: then Vec2Pix
			// may coincide with it, counted. The pixel healpix_cxx gives
			// with correctly rounded inputs (vec2pixPhi on the golden's
			// atan2) must differ from it always.
			m := math.Max(math.Abs(p.VX[i]), math.Max(math.Abs(p.VY[i]), math.Abs(p.VZ[i])))
			if cr := newGeom("test", c.Nside, scheme).vec2pixPhi(p.VX[i], p.VY[i], p.VZ[i], m, p.Atan2CR[i], scheme); cr == d[1] {
				t.Errorf("nside %d %v vector %d * 2^%d: healpy's scaled pixel %d equals the correctly rounded unscaled one", c.Nside, scheme, i, c.K, d[1])
			}
			if got == d[1] {
				if !vec2pixAlts(c.Nside, p.VX[i], p.VY[i], p.VZ[i], scheme)[d[1]] {
					t.Errorf("nside %d %v vector %d * 2^%d: Vec2Pix %d equals healpy's scaled pixel and the vector is not boundary-ambiguous", c.Nside, scheme, i, c.K, got)
				}
				coincide++
			}
			moved++
		}
	}
	if !firstMove {
		t.Errorf("healpy moves no vector at 2^-482 (Nside 2^29-1 RING); the doc says it first moves there")
	}
	t.Logf("%d healpy-moved vectors: Vec2Pix equals healpy's scaled pixel at %d boundary-ambiguous ones", moved, coincide)
	if 10*coincide > moved {
		t.Errorf("Vec2Pix equals healpy's scaled pixel at %d of %d moved vectors, more than 10%%", coincide, moved)
	}
}

// TestVec2PixTinyVector_Healpy pins the Vec2Pix doc's example: healpy
// 1.20.1 vec2pix of (1e-200, 0, 1e-200), whose squares underflow to 0,
// is not the pixel of (1, 0, 1), the same direction: RING gives 0 at
// every Nside > 1 and NESTED another wrong pixel (255 for 171 at Nside
// 16); at Nside 1 both are pixel 0 anyway. Vec2Pix gives healpy's
// unit-vector pixel for both. Truth: vec2pix_scaled_healpy.json "tiny".
func TestVec2PixTinyVector_Healpy(t *testing.T) {
	var g struct {
		Tiny []struct {
			Nside int   `json:"nside"`
			Nest  bool  `json:"nest"`
			Unit  int64 `json:"unit"`
			Tiny  int64 `json:"tiny"`
		} `json:"tiny"`
	}
	mustUnmarshal(t, loadGolden(t, "vec2pix_scaled_healpy.json", 0), &g)
	if len(g.Tiny) != 10 {
		t.Fatalf("vec2pix_scaled_healpy.json: %d tiny cases, want 10", len(g.Tiny))
	}
	for _, c := range g.Tiny {
		scheme := Ring
		if c.Nest {
			scheme = Nested
		}
		switch {
		case c.Nside == 1 && c.Tiny != c.Unit:
			t.Errorf("nside 1 %v: healpy tiny %d, unit %d; the doc says they agree", scheme, c.Tiny, c.Unit)
		case c.Nside > 1 && !c.Nest && c.Tiny != 0:
			t.Errorf("nside %d RING: healpy tiny %d; the doc says 0", c.Nside, c.Tiny)
		case c.Nside > 1 && c.Tiny == c.Unit:
			t.Errorf("nside %d %v: healpy tiny == unit %d; the doc says healpy is wrong there", c.Nside, scheme, c.Unit)
		}
		if got := Vec2Pix(c.Nside, 1, 0, 1, scheme); got != c.Unit {
			t.Errorf("Vec2Pix(%d, 1, 0, 1, %v) = %d, healpy %d", c.Nside, scheme, got, c.Unit)
		}
		if got := Vec2Pix(c.Nside, 1e-200, 0, 1e-200, scheme); got != c.Unit {
			t.Errorf("Vec2Pix(%d, 1e-200, 0, 1e-200, %v) = %d, want %d (healpy's (1, 0, 1))", c.Nside, scheme, got, c.Unit)
		}
	}
}
