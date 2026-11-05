package transform

import (
	"bytes"
	"compress/zlib"
	"encoding/binary"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// Bare-projection wcslib truth tests.
//
// testdata/wref/wref.c links wcslib 8.10 and, for every case it defines,
// runs prjset + prjs2x + prjx2s over fixed grids that run well past every
// projection's boundary, writing testdata/wref/golden/<case>.json (metadata)
// and <case>.bin (every input and wcslib output as exact IEEE-754 bits,
// zlib-compressed; layout in wref.c). This
// test builds the same projection through Select (the PVi_m the case
// names, PV2_m on the latitude axis; every other parameter left to our
// defaults) and requires, for every point in both directions:
//
//   - wcslib rejected the point (stat != 0) <=> we return ok=false;
//   - otherwise the coordinates agree within the budget below.
//
// Missing, stale or orphan goldens are FAILURES, never skips. To
// regenerate after changing wref.c:
//
//	make -C wcs/transform/testdata/wref wref-goldens
//
// and review `git diff`. TestWrefGoldensReproduce rebuilds wref and
// requires it to reproduce every committed golden byte for byte.

const (
	wrefGoldenDir     = "testdata/wref/golden"
	wrefPrjsetStatus  = "prjset_status" // golden/prjset_status.json, not a projection case
	wrefWcslibVersion = "8.10"

	// Coordinate budget, degrees, for every case and direction: 1e-11 deg
	// = 36 nano-arcsec. The Go code is a line-for-line port of wcslib's
	// prj.c (same formulas, same iterations, same stopping rules, the same
	// single-precision arithmetic for CSC, wcslib's wcstrig structure over
	// Go's math functions, wcstrig.go), so the only differences left are
	// last-bit differences between Go's libm and glibc's; a row the
	// projection makes ill-conditioned may exceed the budget only by its
	// own libm sensitivity (libm_sens_test.go; the test log prints each
	// case's maximum and the rows that used the allowance).
	// The smallest real bug the audit found was ~3e-6 deg (0.01 arcsec).
	wrefTolDeg = 1e-11
)

// wrefGolden is one bare-projection case as the tests use it: every row
// [in0, in1, out0, out1, stat], outputs nil where wcslib rejected the point
// or returned a non-finite (or UNDEFINED) coordinate with stat 0.
type wrefGolden struct {
	Name          string
	Code          string
	WcslibVersion string
	PV            [][2]float64
	Prjset        struct {
		Phi0   *float64
		Theta0 *float64
	}
	Forward [][5]*float64
	Inverse [][5]*float64
}

// wrefMeta is golden/<case>.json, the metadata wref.c writes next to the
// binary payload golden/<case>.bin (format: see the header of wref.c).
type wrefMeta struct {
	Name          string       `json:"name"`
	Code          string       `json:"code"`
	WcslibVersion string       `json:"wcslib_version"`
	ZlibVersion   string       `json:"zlib_version"`
	PV            [][2]float64 `json:"pv"`
	Prjset        struct {
		Phi0   *float64 `json:"phi0"`
		Theta0 *float64 `json:"theta0"`
	} `json:"prjset"`
	Forward      int    `json:"forward"`
	Inverse      int    `json:"inverse"`
	Roundtrip    int    `json:"roundtrip"`
	Payload      string `json:"payload"`
	PayloadBytes int    `json:"payload_bytes"`
}

// wcslibUndefined is wcsmath.h's UNDEFINED (987654321.0e99), which wcslib
// can leave in an output it does not set.
const wcslibUndefined = 987654321.0e99

func wrefNames(t *testing.T) []string {
	t.Helper()
	files, err := filepath.Glob(filepath.Join(wrefGoldenDir, "*.json"))
	if err != nil {
		t.Fatal(err)
	}
	if len(files) == 0 {
		t.Fatalf("no bare-projection goldens in %s: run make -C wcs/transform/testdata/wref wref-goldens", wrefGoldenDir)
	}
	var names []string
	have := map[string]bool{}
	for _, f := range files {
		if n := strings.TrimSuffix(filepath.Base(f), ".json"); n != wrefPrjsetStatus {
			names = append(names, n)
			have[n] = true
		}
	}
	bins, _ := filepath.Glob(filepath.Join(wrefGoldenDir, "*.bin"))
	for _, b := range bins {
		if n := strings.TrimSuffix(filepath.Base(b), ".bin"); !have[n] {
			t.Errorf("orphan payload %s: no %s.json", b, n)
		}
	}
	sort.Strings(names)
	// The glob alone would let a deleted case go unnoticed (the tests would
	// just run over fewer cases), so the set must be exactly the one
	// written down in prj_manifest_test.go.
	for _, n := range setDiff(wrefCaseManifest, names) {
		t.Errorf("bare-projection golden %s is MISSING (listed in wrefCaseManifest)", n)
	}
	for _, n := range setDiff(names, wrefCaseManifest) {
		t.Errorf("bare-projection golden %s is not in wrefCaseManifest (prj_manifest_test.go)", n)
	}
	return names
}

// wrefPayload reads golden/<case>.json and decompresses its payload.
func wrefPayload(t *testing.T, dir, name string) (*wrefMeta, []byte) {
	t.Helper()
	data, err := os.ReadFile(filepath.Join(dir, name+".json"))
	if err != nil {
		t.Fatalf("golden %s missing: %v", name, err)
	}
	var m wrefMeta
	if err := json.Unmarshal(data, &m); err != nil {
		t.Fatalf("parse golden %s: %v", name, err)
	}
	if m.Name != name || m.Payload != name+".bin" {
		t.Fatalf("golden %s: metadata names %q / payload %q", name, m.Name, m.Payload)
	}
	z, err := os.ReadFile(filepath.Join(dir, m.Payload))
	if err != nil {
		t.Fatalf("golden %s: payload missing: %v", name, err)
	}
	zr, err := zlib.NewReader(bytes.NewReader(z))
	if err != nil {
		t.Fatalf("golden %s: payload: %v", name, err)
	}
	raw, err := io.ReadAll(zr) // zlib verifies its Adler-32 checksum at EOF
	if err != nil {
		t.Fatalf("golden %s: payload: %v", name, err)
	}
	if len(raw) != m.PayloadBytes {
		t.Fatalf("golden %s: payload is %d bytes, metadata says %d", name, len(raw), m.PayloadBytes)
	}
	return &m, raw
}

func loadWref(t *testing.T, name string) *wrefGolden {
	t.Helper()
	m, raw := wrefPayload(t, wrefGoldenDir, name)
	if m.WcslibVersion != wrefWcslibVersion {
		t.Fatalf("%s generated by wcslib %q, PROVENANCE pins %q", name, m.WcslibVersion, wrefWcslibVersion)
	}
	g := &wrefGolden{Name: m.Name, Code: m.Code, WcslibVersion: m.WcslibVersion, PV: m.PV}
	g.Prjset.Phi0, g.Prjset.Theta0 = m.Prjset.Phi0, m.Prjset.Theta0
	var err error
	g.Forward, g.Inverse, err = wrefDecode(raw)
	if err != nil {
		t.Fatalf("golden %s: %v", name, err)
	}
	if len(g.Forward) != m.Forward || len(g.Inverse) != m.Inverse+m.Roundtrip {
		t.Fatalf("golden %s: %d forward / %d inverse rows, metadata says %d / %d + %d",
			name, len(g.Forward), len(g.Inverse), m.Forward, m.Inverse, m.Roundtrip)
	}
	if len(g.Forward) == 0 || len(g.Inverse) == 0 {
		t.Fatalf("%s has no points", name)
	}
	return g
}

// wrefDecode parses a payload (the layout in wref.c's header) into rows.
// The round-trip inverse inputs are rebuilt from the forward rows exactly
// as wref.c appends them: the (x, y) of every forward row with stat 0 and
// finite outputs, in order.
func wrefDecode(raw []byte) (fwd, inv [][5]*float64, err error) {
	p := 0
	need := func(n int) bool {
		if p+n > len(raw) {
			err = fmt.Errorf("payload truncated at byte %d (need %d more)", p, n)
			return false
		}
		return true
	}
	u32 := func() int {
		if !need(4) {
			return 0
		}
		v := binary.BigEndian.Uint32(raw[p:])
		p += 4
		return int(v)
	}
	f64 := func() float64 {
		if !need(8) {
			return 0
		}
		v := math.Float64frombits(binary.BigEndian.Uint64(raw[p:]))
		p += 8
		return v
	}
	stats := func(n int) []byte {
		if !need(n) {
			return nil
		}
		s := raw[p : p+n]
		p += n
		return s
	}
	val := func(v float64) *float64 { return &v }
	out := func(v float64) *float64 {
		if math.IsNaN(v) || math.IsInf(v, 0) || v == wcslibUndefined {
			return nil
		}
		return &v
	}
	nf, nown, nrt := u32(), u32(), u32()
	if err != nil {
		return nil, nil, err
	}
	fwd = make([][5]*float64, nf)
	for k := range fwd {
		fwd[k][0], fwd[k][1] = val(f64()), val(f64())
	}
	fs := stats(nf)
	if err != nil {
		return nil, nil, err
	}
	var rt [][2]float64
	for k := range fwd {
		fwd[k][4] = val(float64(fs[k]))
		if fs[k] != 0 {
			continue
		}
		x, y := f64(), f64()
		fwd[k][2], fwd[k][3] = out(x), out(y)
		if !math.IsNaN(x) && !math.IsInf(x, 0) && !math.IsNaN(y) && !math.IsInf(y, 0) {
			rt = append(rt, [2]float64{x, y})
		}
	}
	if len(rt) != nrt {
		return nil, nil, fmt.Errorf("%d round-trip inputs from the forward rows, payload header says %d", len(rt), nrt)
	}
	inv = make([][5]*float64, nown+nrt)
	for k := 0; k < nown; k++ {
		inv[k][0], inv[k][1] = val(f64()), val(f64())
	}
	for k, r := range rt {
		inv[nown+k][0], inv[nown+k][1] = val(r[0]), val(r[1])
	}
	is := stats(nown + nrt)
	if err != nil {
		return nil, nil, err
	}
	for k := range inv {
		inv[k][4] = val(float64(is[k]))
		if is[k] == 0 {
			inv[k][2], inv[k][3] = out(f64()), out(f64())
		}
	}
	if err != nil {
		return nil, nil, err
	}
	if p != len(raw) {
		return nil, nil, fmt.Errorf("%d trailing payload bytes", len(raw)-p)
	}
	return fwd, inv, nil
}

// wrefSelect builds the Go projection exactly as a header would: only the
// PV2_m the case names are present.
func wrefSelect(g *wrefGolden) (Projection, error) {
	pv := map[wcs.PVKey]float64{}
	for _, p := range g.PV {
		pv[wcs.PVKey{Axis: 2, Index: int(p[0])}] = p[1]
	}
	return Select(g.Code, pv, 2)
}

// wrefTol returns the coordinate budget for one case and direction.
func wrefTol(*wrefGolden, bool) float64 { return wrefTolDeg }

// wrefFwdScale is the factor by which a forward row's difference is
// divided before it is compared with the budget: 1 while both coordinates
// are within 1000 deg (so the budget is exactly wrefTolDeg there, which
// covers every point of every projection's map: none is wider than 720
// deg), and |coordinate| / 1000 beyond, i.e. a relative budget of 1e-14
// (45 ULP). Only points next to a forward divergence get there: the
// boundary probes in wref.c come within 1e-7 deg of the AZP/TAN/SZP
// horizon, COP/COO's divergent pole and CYP's poles, where wcslib returns
// |x|, |y| of 1e5 to 1e11 deg and a few ULP of pow/tan, magnified by the
// divergence, is worth more than 1e-11 deg. 1e-14 relative is still five
// orders of magnitude below any real projection bug the audit found
// (relative 1e-9 and up).
func wrefFwdScale(w0, w1 float64) float64 {
	return math.Max(1, math.Max(math.Abs(w0), math.Abs(w1))/1000)
}

// wrefDiff compares one golden row with our result, in degrees. For
// inverse rows the longitude difference is reduced mod 360 and measured as
// arc length along the parallel (dphi * cos(theta)), so it is an angle on
// the sky like the latitude difference; at the poles longitude is
// undefined and only latitude counts. The cos(theta) weight matters only
// near a pole: the QSC/TSC probes 0.005 deg from the south pole (wref.c
// extra_inv) differ in phi by 2e-11 deg between libm's atan2 and Go's,
// which is 2e-15 deg on the sky. At |theta| <= 80 deg it relaxes the
// longitude budget by at most a factor 5.8. The latitude difference is
// taken as it is: with the cancellation-free asin/acos of wcstrig.go the
// QSC/TSC/CSC latitudes next to the poles stay within the budget without
// a conditioning allowance; rows elsewhere that need their libm
// sensitivity (e.g. MOL's inverse) are counted in the log. (An earlier version
// allowed 4 ULP of sin(theta) / cos(theta) there, quoting a TSC
// difference of 1.04e-11 deg that the current port no longer has.)
//
// Forward rows are compared in (x, y) with the absolute budget, scaled
// up only where the coordinate itself is huge (wrefFwdScale): next to a
// forward divergence (AZP/AIR/COO/TAN horizons, which the boundary
// probes in wref.c approach to 1e-7 deg) wcslib returns x, y of 1e5 to
// 1e11 deg, and there libm's last bit is worth more than 1e-11 deg.
func wrefDiff(inverse bool, want0, want1, got0, got1 float64) float64 {
	if !inverse {
		return math.Max(math.Abs(got0-want0), math.Abs(got1-want1)) / wrefFwdScale(want0, want1)
	}
	d := math.Abs(got1 - want1)
	if math.Abs(want1) < 90-1e-9 {
		dphi := math.Mod(got0-want0, 360)
		if dphi > 180 {
			dphi -= 360
		} else if dphi < -180 {
			dphi += 360
		}
		d = math.Max(d, math.Abs(dphi)*math.Cos(want1*math.Pi/180))
	}
	return d
}

// wrefPoleAllowDeg caps wrefDiffConditioned's near-pole latitude
// allowance; it is 3000 times smaller than the smallest real projection
// bug the audit found (3e-6 deg).
const wrefPoleAllowDeg = 1e-9

// wrefDiffConditioned is wrefDiff for a truth that is EXACT rather than
// wcslib's own double-precision answer (the paper goldens): next to a pole
// the port's latitude, theta = asind(n) with n = sin(theta) ~ 1, is
// ill-conditioned, so a 1-ULP change in n moves theta by eps / cos(theta)
// rad, which no double-precision evaluation can avoid. The latitude
// difference is reduced by 4 ULP of n (4 * 2^-52 / cos(theta) rad, in
// degrees, capped at wrefPoleAllowDeg) before it is compared: 5e-14 deg
// at theta = 0, 3e-13 deg at |theta| = 80, 2.9e-10 deg at |theta| = 89.99
// (the XPH paper rows there differ from the exact inverse by 1.74e-11).
func wrefDiffConditioned(inverse bool, want0, want1, got0, got1 float64) float64 {
	if !inverse {
		return wrefDiff(inverse, want0, want1, got0, got1)
	}
	dt := math.Abs(got1 - want1)
	if c := math.Cos(want1 * math.Pi / 180); c > 0 {
		dt = math.Max(0, dt-math.Min(4*0x1p-52/c*180/math.Pi, wrefPoleAllowDeg))
	}
	if math.Abs(want1) < 90-1e-9 {
		dphi := math.Mod(got0-want0, 360)
		if dphi > 180 {
			dphi -= 360
		} else if dphi < -180 {
			dphi += 360
		}
		return math.Max(dt, math.Abs(dphi)*math.Cos(want1*math.Pi/180))
	}
	return dt
}

func TestWrefProjections(t *testing.T) {
	names := wrefNames(t)

	// Every projection code we claim must have at least one wcslib case.
	covered := map[string]bool{}
	for _, n := range names {
		covered[loadWref(t, n).Code] = true
	}
	for code := range allProjectionCodes {
		if !covered[code] {
			t.Errorf("projection %s has no wcslib bare-projection golden", code)
		}
	}

	wrefCondTotal = condTally{}
	defer func() {
		// The libm allowances (libm_sens_test.go) must stay the exception:
		// they are used only next to a pole, a horizon or a divergence
		// (AIR, AZP, HPX, MOL, QSC, STG, XPH); see the log for the counts.
		c := wrefCondTotal
		t.Logf("all wref cases: %d rows, %d accepted with their libm sensitivity, %d boundary-ambiguous decisions", c.rows, c.conditioned, c.ambiguous)
		if c.rows < 50000 {
			return // a -run subset: the shares are only meaningful over all cases
		}
		if float64(c.conditioned) > 0.01*float64(c.rows) || float64(c.ambiguous) > 1e-4*float64(c.rows) {
			t.Errorf("libm allowances used on %d (conditioned) and %d (ambiguous) of %d rows: more than 1%% / 0.01%%", c.conditioned, c.ambiguous, c.rows)
		}
	}()
	for _, name := range names {
		t.Run(name, func(t *testing.T) {
			g := loadWref(t, name)
			proj, err := wrefSelect(g)
			if err != nil {
				t.Fatalf("Select(%s, %v): %v (wcslib prjset accepted it)", g.Code, g.PV, err)
			}
			if _, isDeg := proj.(degAdapter); !isDeg {
				t.Errorf("Select(%s) is not a degree-unit port of wcslib (degAdapter)", g.Code)
			}
			if g.Prjset.Theta0 != nil {
				if d := math.Abs(radToDeg(proj.Theta0()) - *g.Prjset.Theta0); d > 1e-12 {
					t.Errorf("Theta0 = %.15g deg, wcslib prjset theta0 = %.15g", radToDeg(proj.Theta0()), *g.Prjset.Theta0)
				}
			}
			wrefCheckLatAxis(t, g, proj)
			for _, inverse := range []bool{false, true} {
				wrefCheckRows(t, g, inverse, func(in0, in1 float64) (float64, float64, bool) {
					return wrefEvalDeg(proj, inverse, in0, in1)
				})
				if d, isDeg := proj.(degAdapter); isDeg {
					wrefCheckAdapter(t, g, d, inverse)
				}
			}
		})
	}
}

// wrefCheckLatAxis requires the projection parameters to be read from the
// latitude axis, whichever it is: the same case built with its PV on axis
// 1 and latAxis = 1 must return exactly what the axis-2 build returns at
// every golden row (which the caller compares with wcslib), and PV cards
// on the longitude axis (here values that would change every
// parameterised projection) must be ignored.
func wrefCheckLatAxis(t *testing.T, g *wrefGolden, proj Projection) {
	t.Helper()
	pv1 := map[wcs.PVKey]float64{}
	mixed := map[wcs.PVKey]float64{}
	for _, p := range g.PV {
		pv1[wcs.PVKey{Axis: 1, Index: int(p[0])}] = p[1]
		mixed[wcs.PVKey{Axis: 2, Index: int(p[0])}] = p[1]
	}
	for m := 1; m <= 5; m++ {
		mixed[wcs.PVKey{Axis: 1, Index: m}] = 0.37 * float64(m+1)
	}
	for _, alt := range []struct {
		what string
		pv   map[wcs.PVKey]float64
		lat  int
	}{{"PV on axis 1, latAxis = 1", pv1, 1}, {"extra PV1_m on the longitude axis", mixed, 2}} {
		p2, err := Select(g.Code, alt.pv, alt.lat)
		if err != nil {
			t.Errorf("%s: Select: %v", alt.what, err)
			continue
		}
		var nbad int
		for _, inverse := range []bool{false, true} {
			rows := g.Forward
			if inverse {
				rows = g.Inverse
			}
			for _, r := range rows {
				a0, a1, aok := wrefEvalDeg(proj, inverse, *r[0], *r[1])
				b0, b1, bok := wrefEvalDeg(p2, inverse, *r[0], *r[1])
				same := aok == bok && (!aok || (math.Float64bits(a0) == math.Float64bits(b0) && math.Float64bits(a1) == math.Float64bits(b1)))
				if !same {
					nbad++
					if nbad <= 3 {
						t.Errorf("%s: (%g, %g) inverse=%v: got (%v, %v, %v), latitude-axis-2 build (%v, %v, %v)",
							alt.what, *r[0], *r[1], inverse, bok, b0, b1, aok, a0, a1)
					}
				}
			}
		}
	}
}

// wrefEvalDeg evaluates one golden input through the degree-unit port
// (every code has one since round 2 of the proof audit, when TSC and QSC
// stopped being radian-native ports), else through the radian interface.
func wrefEvalDeg(proj Projection, inverse bool, in0, in1 float64) (float64, float64, bool) {
	if d, ok := proj.(degAdapter); ok {
		if inverse {
			return d.inverseDeg(in0, in1)
		}
		return d.forwardDeg(in0, in1)
	}
	var o0, o1 float64
	var ok bool
	if inverse {
		o0, o1, ok = proj.Inverse(degToRad(in0), degToRad(in1))
	} else {
		o0, o1, ok = proj.Forward(degToRad(in0), degToRad(in1))
	}
	return radToDeg(o0), radToDeg(o1), ok
}

// wrefCheckAdapter checks the public radian interface (Forward/Inverse)
// over the same rows: it must be exactly the degree port evaluated at the
// degree value the radian input actually represents, plus the deg<->rad
// conversion of the result. The degree port itself is what is compared
// with wcslib; this proves the adapter adds nothing but unit conversion.
// (Comparing the radian API with wcslib directly is not meaningful at
// inputs that do not survive deg->rad->deg exactly and sit on an exact
// boundary, e.g. ZEA at R = 2 r0 or an HPX facet seam: there a 1-ULP input
// change legitimately moves the answer by up to ~1e-6 deg or across the
// seam.)
func wrefCheckAdapter(t *testing.T, g *wrefGolden, d degAdapter, inverse bool) {
	t.Helper()
	rows := g.Forward
	if inverse {
		rows = g.Inverse
	}
	const tol = 1e-11 // deg; the result's rad->deg rounding at |v| <= 1000 (relative 1e-14 beyond, wrefFwdScale)
	var nbad int
	for _, r := range rows {
		in0, in1 := degToRad(*r[0]), degToRad(*r[1])
		var o0, o1, e0, e1 float64
		var ok, eok bool
		if inverse {
			o0, o1, ok = d.Inverse(in0, in1)
			e0, e1, eok = d.inverseDeg(radToDeg(in0), radToDeg(in1))
		} else {
			o0, o1, ok = d.Forward(in0, in1)
			e0, e1, eok = d.forwardDeg(radToDeg(in0), radToDeg(in1))
		}
		o0, o1 = radToDeg(o0), radToDeg(o1)
		scale := 1.0
		if !inverse {
			scale = wrefFwdScale(e0, e1)
		}
		if ok != eok || (ok && math.Max(math.Abs(o0-e0), math.Abs(o1-e1))/scale > tol) {
			nbad++
			if nbad <= 5 {
				t.Errorf("radian API (%g, %g): got (%v, %.15g, %.15g), degree port (%v, %.15g, %.15g)",
					*r[0], *r[1], ok, o0, o1, eok, e0, e1)
			}
		}
	}
}

func wrefCheckRows(t *testing.T, g *wrefGolden, inverse bool, eval func(in0, in1 float64) (float64, float64, bool)) {
	t.Helper()
	dir, rows := "forward", g.Forward
	if inverse {
		dir, rows = "inverse", g.Inverse
	}
	tol := wrefTol(g, inverse)
	var nbad, nrej, nok int
	var maxErr float64
	var cond condTally
	report := func(format string, args ...any) {
		nbad++
		if nbad <= 12 {
			t.Errorf("%s: "+format, append([]any{dir}, args...)...)
		}
	}
	for _, r := range rows {
		in0, in1, stat := *r[0], *r[1], *r[4]
		o0, o1, ok := eval(in0, in1)
		cond.rows++
		if stat != 0 {
			nrej++
			if ok {
				if decisionAmbiguous(eval, in0, in1) {
					cond.ambiguous++
					continue
				}
				report("(%g, %g): wcslib rejects (stat %g), we returned (%.15g, %.15g)", in0, in1, stat, o0, o1)
			}
			continue
		}
		if r[2] == nil || r[3] == nil {
			// wcslib reports success but its output is NaN or infinite
			// (an out-of-range asin/atan its bounds check does not catch:
			// CYP with |eta mu/sqrt(eta^2+1)| > 1, PCO far off the map,
			// COO at the apex of a southern cone).
			// There is no coordinate to agree with; we must refuse.
			nrej++
			if ok {
				report("(%g, %g): wcslib output is non-finite, we returned (%.15g, %.15g)", in0, in1, o0, o1)
			}
			continue
		}
		w0, w1 := *r[2], *r[3]
		if !ok {
			if decisionAmbiguous(eval, in0, in1) {
				cond.ambiguous++
				continue
			}
			report("(%g, %g): wcslib gives (%.15g, %.15g), we rejected", in0, in1, w0, w1)
			continue
		}
		nok++
		d := wrefDiff(inverse, w0, w1, o0, o1)
		maxErr = math.Max(maxErr, d)
		sens := func() float64 {
			return libmSensitivity(func() (float64, float64, bool) { return eval(in0, in1) },
				func(a0, a1, b0, b1 float64) float64 { return wrefDiff(inverse, a0, a1, b0, b1) })
		}
		if !cond.valueOK(d, tol, sens) {
			report("(%.15g, %.15g): got (%.15g, %.15g), wcslib (%.15g, %.15g), diff %.3g deg > %g + libm sensitivity %.3g",
				in0, in1, o0, o1, w0, w1, d, tol, sens())
		}
	}
	if nbad > 12 {
		t.Errorf("%s: ... %d mismatches in total", dir, nbad)
	}
	if nbad == 0 {
		t.Logf("%s: %d points agree (%d rejected by both), max diff %.3g deg (budget %g)",
			dir, nok+nrej, nrej, maxErr, tol)
	}
	cond.log(t, dir)
	wrefCondTotal.add(cond)
}

// wrefCondTotal sums the allowance counts over every wref case and
// direction; TestWrefConditioningRare bounds them.
var wrefCondTotal condTally

func (c *condTally) add(o condTally) {
	c.rows += o.rows
	c.conditioned += o.conditioned
	c.ambiguous += o.ambiguous
	c.worstRatio = math.Max(c.worstRatio, o.worstRatio)
}

// TestSelectRejectsWhatPrjsetRejects requires Select to refuse exactly the
// parameter sets wcslib's prjset refuses (golden/prjset_status.json, from
// wref.c): a header wcslib cannot set up must be an error here too, never
// a projection built from guessed parameters (e.g. BON without PV2_1 used
// to become SFL, CEA with lambda = 0 or 1.5 and AZP with mu = -1 used to
// be accepted).
func TestSelectRejectsWhatPrjsetRejects(t *testing.T) {
	data, err := os.ReadFile(filepath.Join(wrefGoldenDir, wrefPrjsetStatus+".json"))
	if err != nil {
		t.Fatalf("golden missing: %v", err)
	}
	var g struct {
		WcslibVersion string `json:"wcslib_version"`
		Cases         []struct {
			Name   string       `json:"name"`
			Code   string       `json:"code"`
			PV     [][2]float64 `json:"pv"`
			Status int          `json:"prjset_status"`
		} `json:"cases"`
	}
	if err := json.Unmarshal(data, &g); err != nil {
		t.Fatal(err)
	}
	if g.WcslibVersion != wrefWcslibVersion || len(g.Cases) == 0 {
		t.Fatalf("bad golden: version %q, %d cases", g.WcslibVersion, len(g.Cases))
	}
	for _, c := range g.Cases {
		pv := map[wcs.PVKey]float64{}
		for _, p := range c.PV {
			pv[wcs.PVKey{Axis: 2, Index: int(p[0])}] = p[1]
		}
		_, err := Select(c.Code, pv, 2)
		if (err != nil) != (c.Status != 0) {
			t.Errorf("%s: Select(%s, %v) error = %v, wcslib prjset status = %d", c.Name, c.Code, c.PV, err, c.Status)
		}
	}
}
