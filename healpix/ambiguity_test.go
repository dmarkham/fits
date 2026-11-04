package healpix

import (
	"math"
	"testing"
)

// Boundary ambiguity.
//
// Ang2Pix and Vec2Pix call Go's math.Cos, math.Sin and math.Atan2, which
// are not correctly rounded; healpy calls glibc's, and the polar and
// switch-over expected values are healpix_cxx fed correctly rounded ones.
// The three can differ in the last bit, and for a direction within about
// an ULP of a pixel boundary that last bit decides the pixel. Such a
// direction has no meaningful "correct" side, so a test may accept a
// pixel that differs from the reference there, and only there:
//
// A point is boundary-ambiguous when the pixel Go computes changes under a
// nudge of at most nudgeULP ULP of
//
//   - the inputs: theta and phi (Ang2Pix), or vx, vy, vz (Vec2Pix), each
//     independently, or
//   - the libm results the pixel is computed from: cos(theta), sin(theta)
//     (Ang2Pix, through ang2pixCS) and atan2(vy, vx) (Vec2Pix, through
//     vec2pixPhi, combined with the input nudges).
//
// The second kind is needed next to the poles: there an ULP of theta is
// far below an ULP of cos(theta) (at theta = 0.01, 1 ULP of z is 6500 ULP
// of theta), so no input nudge of a few ULP reproduces a 1-ULP difference
// in cos(theta), which is exactly what separates math.Cos from glibc.
// The nudged evaluations run the production code after the libm calls
// (ang2pixCS, vec2pixPhi), so a changed branch or threshold there changes
// the alternatives too.
//
// A mismatch is accepted only when the reference pixel is one of the
// pixels the nudges produce; anything else fails. Each test counts the
// accepted mismatches, logs the count and bounds it. The same tests also
// run the production code after the libm calls (ang2pixCS, vec2pixPhi) on
// the correctly rounded cos/sin/atan2 values the reference was computed
// from (stored in the golden) and require the reference pixel exactly
// there, so the allowance never covers the arithmetic or the branch
// choices themselves, only the last bit of Go's libm.
const nudgeULP = 4

// polarAmbiguousMax bounds the share of ang2pix_polar.json pixels that may
// differ from the reference at boundary-ambiguous points. Every one of
// those directions lies within 2 ULP of a boundary by construction;
// measured: 2.7% overall, 4.2% at worst (Nside 2^29).
const polarAmbiguousMax = 0.06

// nudged returns x moved by k ULP (k < 0 toward -Inf).
func nudged(x float64, k int) float64 {
	for ; k > 0; k-- {
		x = math.Nextafter(x, math.Inf(1))
	}
	for ; k < 0; k++ {
		x = math.Nextafter(x, math.Inf(-1))
	}
	return x
}

// ang2pixAlts returns the pixels Ang2Pix can give for (theta, phi) under
// the nudges described above (its own unnudged answer included).
func ang2pixAlts(nside int, theta, phi float64, scheme Scheme) map[int64]bool {
	alts := map[int64]bool{}
	for dt := -nudgeULP; dt <= nudgeULP; dt++ {
		tt := nudged(theta, dt)
		if !(tt >= 0 && tt <= pi) {
			continue
		}
		for dp := -nudgeULP; dp <= nudgeULP; dp++ {
			alts[Ang2Pix(nside, tt, nudged(phi, dp), scheme)] = true
		}
	}
	g := newGeom("ang2pixAlts", nside, scheme)
	c, sn := math.Cos(theta), math.Sin(theta)
	for dc := -nudgeULP; dc <= nudgeULP; dc++ {
		cc := nudged(c, dc)
		if cc > 1 || cc < -1 {
			continue
		}
		for ds := -nudgeULP; ds <= nudgeULP; ds++ {
			ss := nudged(sn, ds)
			if ss < 0 {
				continue
			}
			for dp := -nudgeULP; dp <= nudgeULP; dp++ {
				alts[g.ang2pixCS(theta, nudged(phi, dp), cc, ss, scheme)] = true
			}
		}
	}
	return alts
}

// vec2pixAlts is ang2pixAlts for Vec2Pix: every input nudge, each with
// atan2 nudged too.
func vec2pixAlts(nside int, vx, vy, vz float64, scheme Scheme) map[int64]bool {
	alts := map[int64]bool{}
	g := newGeom("vec2pixAlts", nside, scheme)
	for dx := -nudgeULP; dx <= nudgeULP; dx++ {
		for dy := -nudgeULP; dy <= nudgeULP; dy++ {
			for dz := -nudgeULP; dz <= nudgeULP; dz++ {
				x, y, z := nudged(vx, dx), nudged(vy, dy), nudged(vz, dz)
				m := math.Max(math.Abs(x), math.Max(math.Abs(y), math.Abs(z)))
				if m == 0 {
					continue
				}
				alts[Vec2Pix(nside, x, y, z, scheme)] = true
				if x == 0 && y == 0 {
					continue
				}
				phi := math.Atan2(y, x)
				for dp := -nudgeULP; dp <= nudgeULP; dp++ {
					alts[g.vec2pixPhi(x, y, z, m, nudged(phi, dp), scheme)] = true
				}
			}
		}
	}
	return alts
}

// pixTally counts the reference mismatches a test accepted as
// boundary-ambiguous.
type pixTally struct {
	checked, ambiguous, bad int
}

// check compares one pixel with the reference: equal, or a mismatch that
// is accepted when want is among alts() (computed only on a mismatch),
// else an error (the first few are reported).
func (p *pixTally) check(t testing.TB, got, want int64, alts func() map[int64]bool, describe func() string) {
	t.Helper()
	p.checked++
	if got == want {
		return
	}
	if alts()[want] {
		p.ambiguous++
		return
	}
	if p.bad++; p.bad <= 8 {
		t.Errorf("%s = %d, reference %d, not boundary-ambiguous (no nudge of up to %d ULP gives the reference)", describe(), got, want, nudgeULP)
	}
}

// done reports the totals and fails when more than maxFrac of the checked
// points needed the ambiguity allowance.
func (p *pixTally) done(t testing.TB, where string, maxFrac float64) {
	t.Helper()
	if p.bad > 8 {
		t.Errorf("%s: ... %d non-ambiguous mismatches in total", where, p.bad)
	}
	t.Logf("%s: %d pixels checked, %d differ from the reference at boundary-ambiguous points (%.3g%%)",
		where, p.checked, p.ambiguous, 100*float64(p.ambiguous)/math.Max(1, float64(p.checked)))
	if float64(p.ambiguous) > maxFrac*float64(p.checked) {
		t.Errorf("%s: %d of %d pixels differ at boundary-ambiguous points, more than the %.3g%% this test allows",
			where, p.ambiguous, p.checked, 100*maxFrac)
	}
}

// explicitValues returns the golden's explicit list (a digest cannot tell
// which elements differ, so ambiguity-tolerant tests need explicit lists).
func explicitValues(t testing.TB, where string, a i64s) []int64 {
	t.Helper()
	if a.explicit == nil {
		t.Fatalf("%s: golden array is missing or a digest; the ambiguity check needs explicit values (gen_golden.py)", where)
	}
	return a.explicit
}

// TestBoundaryAmbiguityDetects: the nudge sets really contain both sides
// of a boundary and nothing far away. At theta = 0.8410686705679303
// (cos(theta) within an ULP of 2/3, the equatorial/polar split) the
// nudged cosines fall on both sides of the split, which gives different
// pixels for most phi; a direction at a pixel centre has a single
// alternative.
func TestBoundaryAmbiguityDetects(t *testing.T) {
	split := 0
	for k := 0; k < 64; k++ {
		if len(ang2pixAlts(16, 0.8410686705679303, float64(k)*math.Pi/32, Ring)) >= 2 {
			split++
		}
	}
	if split < 16 {
		t.Errorf("z = 2/3: nudges reach both sides of the split for only %d of 64 phi", split)
	}
	th, ph := Pix2Ang(16, 1000, Ring)
	if a := ang2pixAlts(16, th, ph, Ring); len(a) != 1 || !a[1000] {
		t.Errorf("pixel centre: nudges give %v, want only 1000", a)
	}
	x, y, z := Pix2Vec(16, 1000, Nested)
	if a := vec2pixAlts(16, x, y, z, Nested); len(a) != 1 || !a[1000] {
		t.Errorf("pixel centre (vector): nudges give %v, want only 1000", a)
	}
}
