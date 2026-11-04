package transform

import (
	"math"
	"testing"
)

// TestQuadcubeFaceCentreBand pins QSC's forward small-angle switch
// (quadcube.go qscZecoSmall: qscs2x uses zeco = (p^2 + t^2)/2 instead of
// 1 - zeta when 1 - zeta < 1e-8, within ~0.0081 deg of a face centre)
// against wcslib 8.10, and holds TSC and CSC to the same standard next to
// the face centres.
//
// The two zeco formulas agree there to a relative ~1e-8 (1 - zeta loses
// the last bits of zeta), so a port that moves the switch is wrong only
// inside a thin ring around each face centre and by at most ~6e-11 deg:
// round 2 of the proof audit moved it to 0.5e-8 and wrote it as <=, and
// both survived every test, because no golden point had zeco between
// 0.5e-8 and 1e-8. wref.c centre_fwd now puts forward points at
// zeco = 1e-8 (1 +- {5e-4, 1e-3, 1e-2, 1e-1}), 0.5e-8 and 2e-8 (zeco as
// qscs2x computes it) on all six faces, in eight directions each; their
// wcslib images join the inverse set as round trips.
//
// On these rows the port runs wcslib's own operations (same face, same
// branch, the same degree arithmetic and constants), so both plane
// coordinates are held to 4 ULP of wcslib's value, and the round-trip
// inverse to 4 ULP of theta and of phi, each plus the row's own libm
// sensitivity (libm_sens_test.go: the change under a 4-ULP move of every
// libm result, which is what separates Go's math functions from glibc's).
// That allowance matters only where QSC magnifies the last bit of atan2:
// next to the centres of faces 2, 3 and 4 zeco is built from
// atan2(sin phi, cos phi) -+ pi/2 or -+ pi, a cancellation worth 5e4 ULP
// of the plane coordinate, and next to the pole face (zeco from 1 - zeta
// near theta = 90) up to 4e7 ULP; the log counts the rows that needed
// it. A moved switch changes x or y by ~1e-11 deg, more than that
// allowance: the switch moved to 0.5e-8, 2e-8 or written as <= still
// fails (checked by hand when the allowance was introduced). The test
// log gives the measured worst case per code.
func TestQuadcubeFaceCentreBand(t *testing.T) {
	ulp := func(v float64) float64 {
		a := math.Abs(v)
		return math.Nextafter(a, math.Inf(1)) - a
	}
	ulps := func(o, w float64) float64 {
		if o == w {
			return 0
		}
		return math.Abs(o-w) / math.Max(ulp(w), math.SmallestNonzeroFloat64)
	}
	for _, name := range []string{"qsc_default", "tsc_default", "csc_default"} {
		g := loadWref(t, name)
		proj, err := wrefSelect(g)
		if err != nil {
			t.Fatal(err)
		}
		// zeco below 2.5e-8 (selection only; the truth is wcslib's output).
		nearCentre := func(phi, theta float64) bool {
			sp, cp := math.Sincos(phi * math.Pi / 180)
			st, ct := math.Sincos(theta * math.Pi / 180)
			l, m, n := ct*cp, ct*sp, st
			z := math.Max(math.Max(math.Max(l, -l), math.Max(m, -m)), math.Max(n, -n))
			return 1-z < 2.5e-8
		}
		var cond condTally
		ulpDist := func(a0, a1, b0, b1 float64) float64 { return math.Max(ulps(b0, a0), ulps(b1, a1)) }
		sens := func(inverse bool, in0, in1 float64) func() float64 {
			return func() float64 {
				return libmSensitivity(func() (float64, float64, bool) { return wrefEvalDeg(proj, inverse, in0, in1) }, ulpDist)
			}
		}
		band := map[[2]uint64]bool{}
		var nf, ni int
		var below, above int
		var worstF, worstI float64
		for _, r := range g.Forward {
			if !nearCentre(*r[0], *r[1]) {
				continue
			}
			if *r[4] != 0 || r[2] == nil || r[3] == nil {
				t.Errorf("%s forward (%.17g, %.17g): wcslib rejected a face-centre point", name, *r[0], *r[1])
				continue
			}
			if z := qscZecoWcslib(*r[0], *r[1]); z < qscZecoSmall {
				below++
			} else {
				above++
			}
			o0, o1, ok := wrefEvalDeg(proj, false, *r[0], *r[1])
			if !ok {
				t.Errorf("%s forward (%.17g, %.17g): rejected, wcslib gives (%.17g, %.17g)", name, *r[0], *r[1], *r[2], *r[3])
				continue
			}
			nf++
			u := math.Max(ulps(o0, *r[2]), ulps(o1, *r[3]))
			worstF = math.Max(worstF, u)
			cond.rows++
			if !cond.valueOK(u, 4, sens(false, *r[0], *r[1])) {
				t.Errorf("%s forward (%.17g, %.17g): got (%.17g, %.17g), wcslib (%.17g, %.17g): %.3g ULP, budget 4 + libm sensitivity %.3g",
					name, *r[0], *r[1], o0, o1, *r[2], *r[3], u, sens(false, *r[0], *r[1])())
			}
			band[[2]uint64{math.Float64bits(*r[2]), math.Float64bits(*r[3])}] = true
		}
		for _, r := range g.Inverse {
			if !band[[2]uint64{math.Float64bits(*r[0]), math.Float64bits(*r[1])}] {
				continue
			}
			if *r[4] != 0 || r[2] == nil || r[3] == nil {
				t.Errorf("%s inverse (%.17g, %.17g): wcslib rejected a round trip", name, *r[0], *r[1])
				continue
			}
			o0, o1, ok := wrefEvalDeg(proj, true, *r[0], *r[1])
			if !ok {
				t.Errorf("%s inverse (%.17g, %.17g): rejected, wcslib gives (%.17g, %.17g)", name, *r[0], *r[1], *r[2], *r[3])
				continue
			}
			ni++
			u := math.Max(ulps(o0, *r[2]), ulps(o1, *r[3]))
			worstI = math.Max(worstI, u)
			cond.rows++
			if !cond.valueOK(u, 4, sens(true, *r[0], *r[1])) {
				t.Errorf("%s inverse (%.17g, %.17g): got (%.17g, %.17g), wcslib (%.17g, %.17g): %.3g ULP, budget 4 + libm sensitivity %.3g",
					name, *r[0], *r[1], o0, o1, *r[2], *r[3], u, sens(true, *r[0], *r[1])())
			}
		}
		// 6 faces x 8 directions x 10 zeco values (wref.c centre_fwd), 5 of
		// the 10 below the switch; plus extra_fwd's face-centre points.
		if nf < 480 || ni < 480 || below < 240 || above < 240 {
			t.Errorf("%s: %d forward / %d inverse face-centre rows (%d with zeco < 1e-8, %d above), want >= 480 / 480 (240 / 240)",
				name, nf, ni, below, above)
		}
		t.Logf("%s: %d forward rows (zeco < 1e-8: %d), worst %.3g ULP; %d round trips, worst %.3g ULP", name, nf, below, worstF, ni, worstI)
		cond.log(t, name)
		if 10*cond.conditioned > cond.rows {
			t.Errorf("%s: %d of %d rows needed their libm sensitivity, more than 10%%", name, cond.conditioned, cond.rows)
		}
	}
}

// qscZecoWcslib is 1 - zeta as qscs2x computes it (wref.c qsc_zeco),
// used to count the probes on each side of the switch.
func qscZecoWcslib(phi, theta float64) float64 {
	sinphi, cosphi := sincosd(phi)
	sinthe, costhe := sincosd(theta)
	_, zeta := quadcubeFace(float64(costhe*cosphi), float64(costhe*sinphi), sinthe)
	return 1.0 - zeta
}

// TestQSCZecoSwitchTie proves that qscs2x's zeco < 1e-8 can never be
// tied, so < and <= decide alike (the <= mutant survives for this reason,
// not for want of a probe). Below the switch zeta = 1 - zeco is in
// (1 - 1e-8, 1]: a double in [1/2, 1], so 1 - zeta is computed exactly
// (Sterbenz) and is an integer multiple of 2^-53 (zeta's ULP in [1/2, 1)).
// The double 1e-8 is not such a multiple.
func TestQSCZecoSwitchTie(t *testing.T) {
	k := math.Ldexp(qscZecoSmall, 53)
	if k == math.Trunc(k) {
		t.Fatalf("1e-8 * 2^53 = %v is an integer: 1 - zeta can equal the threshold and needs a wcslib probe", k)
	}
	// The doubles zeta next to 1 - 1e-8: 1 - zeta takes the values
	// either side of 1e-8 and never 1e-8 itself.
	var below, above bool
	zu, zd := 1-qscZecoSmall, 1-qscZecoSmall
	for i := 0; i < 1000; i++ {
		for _, z := range []float64{zu, zd} {
			c := 1.0 - z
			if c == qscZecoSmall {
				t.Errorf("zeta = %.17g gives zeco = 1e-8 exactly", z)
			}
			below = below || c < qscZecoSmall
			above = above || c > qscZecoSmall
		}
		zu, zd = math.Nextafter(zu, 2), math.Nextafter(zd, 0)
	}
	if !below || !above {
		t.Errorf("1 - zeta does not cross 1e-8 (below %v, above %v)", below, above)
	}
}

// TestQSCInverseZetaRange checks the premises of qscInverseZeta's proof
// that qscx2s's zeta < -1 rejection is unreachable (so the qsc_zeta_tol
// mutants, which widen its 1e-12 tolerance, are equivalent): every face
// coordinate quadcubeFaceXY returns is within [-1, 1], and for every
// such (xf, yf) zeco is in [0, 0.5] (zeta in [0.5, 1]), |omega| <= 1 and
// tau in [1, 2] (to rounding). Checked on a dense grid of the unfolded
// cube (0.05 deg plus every face and cross edge and the doubles either
// side of it) and on a dense face-coordinate grid with its edges,
// corners and centre lines.
func TestQSCInverseZetaRange(t *testing.T) {
	edges := []float64{}
	for k := -9; k <= 9; k++ {
		e := 45.0 * float64(k)
		edges = append(edges, e, math.Nextafter(e, -1e9), math.Nextafter(e, 1e9))
	}
	var nface int
	plane := func(x, y float64) {
		_, xf, yf, ok := quadcubeFaceXY(x, y)
		if !ok {
			return
		}
		nface++
		if math.Abs(xf) > 1 || math.Abs(yf) > 1 {
			t.Errorf("quadcubeFaceXY(%.17g, %.17g) = (%.17g, %.17g): face coordinate outside [-1, 1]", x, y, xf, yf)
		}
	}
	for x := -320.0; x <= 320; x += 0.05 {
		for _, y := range edges {
			plane(x, y)
			plane(y, x)
		}
	}
	for x := -320.0; x <= 320; x += 0.5 {
		for y := -140.0; y <= 140; y += 0.5 {
			plane(x, y)
		}
	}
	if nface < 100000 {
		t.Fatalf("only %d plane points on the cube", nface)
	}
	fc := []float64{-1, math.Nextafter(-1, 0), -0.5, math.Copysign(0, -1), 0, 1e-300, 0.5, math.Nextafter(1, 0), 1}
	for i := -1000; i <= 1000; i++ {
		fc = append(fc, float64(i)/1000)
	}
	var n int
	for _, xf := range fc {
		for _, yf := range fc {
			_, omega, tau, zeta, zeco := qscInverseZeta(xf, yf)
			n++
			if !(zeco >= 0 && zeco <= 0.5 && zeta >= 0.5 && zeta <= 1 && math.Abs(omega) <= 1+1e-15 && tau >= 1 && tau <= 2+1e-15) {
				t.Errorf("qscInverseZeta(%.17g, %.17g): omega %.17g tau %.17g zeta %.17g zeco %.17g outside the proven ranges", xf, yf, omega, tau, zeta, zeco)
			}
		}
	}
	t.Logf("%d plane points, %d face-coordinate pairs: zeta always in [0.5, 1]", nface, n)
}
