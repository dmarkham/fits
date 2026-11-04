package transform

import (
	"encoding/json"
	"math"
	"path/filepath"
	"regexp"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TPV headers with SIP keywords: the deliberate difference from wcslib.
//
// A header whose celestial CTYPEs are RA---TPV / DEC--TPV (no "-SIP")
// and that also carries SIP keywords is read with TPV only: the SIP set
// is ignored, as astropy 8.0.2 ignores it. wcslib 8.10 applies both (SIP
// as the prior distortion, TPV as the sequent one). The real case is the
// ATLAS difference stamp of PR #4 (hdrcases_differ/real_atlas_k2_tpv_sip),
// whose SIP set is a copy of its TPV correction: applying both moves every
// position on the stamp by 13-14 px (24-26 arcsec), and astropy's reading
// is the one that matches the image. With an explicit "-SIP" suffix
// (RA---TPV-SIP) both apply, as in wcslib; on TAN, SIN and every other
// projection SIP keywords apply with or without the suffix, as in wcslib.

const (
	hdrDifferDir   = "testdata/wref/hdrcases_differ"
	astropyHdrTool = "testdata/wref/astropy_hdr.py"
)

var (
	wcslibDiffersRe = regexp.MustCompile(`(?m)^# WCSLIB-DIFFERS-ARCSEC: (\S+) (\S+)$`)
	astropyTolASRe  = regexp.MustCompile(`(?m)^# ASTROPY-TOL-ARCSEC: (\S+)$`)
	astropyTolPxRe  = regexp.MustCompile(`(?m)^# ASTROPY-TOL-PX: (\S+)$`)
)

// astropyHdr is the JSON astropy_hdr.py prints for a case file.
type astropyHdr struct {
	Astropy  string       `json:"astropy"`
	Wcslib   string       `json:"wcslib"`
	Numpy    string       `json:"numpy"`
	Warnings []string     `json:"warnings"`
	Log      []string     `json:"log"`
	SIP      bool         `json:"sip"`
	P2S      [][4]float64 `json:"p2s"` // x, y (1-based), lng, lat
	S2P      [][4]float64 `json:"s2p"` // lng, lat, x, y (1-based)
}

// astropyHdrReading is astropy 8.0.2's reading of a case file (recorded
// in tier 1, run live in tier 2).
func astropyHdrReading(t *testing.T, caseFile string) *astropyHdr {
	t.Helper()
	out := testutil.AstropyScript(t, astropyHdrTool, caseFile)
	var a astropyHdr
	if err := json.Unmarshal(out, &a); err != nil {
		t.Fatalf("astropy_hdr.py %s: %v", caseFile, err)
	}
	if a.Astropy != "8.0.2" {
		t.Fatalf("astropy %s, the truth is pinned to 8.0.2", a.Astropy)
	}
	return &a
}

// astropyDroppedSIP fails unless astropy says it removed the SIP set.
func astropyDroppedSIP(t *testing.T, a *astropyHdr) {
	t.Helper()
	const msg = "Removed redundant SIP distortion parameters because CTYPE explicitly specifies TPV distortions"
	found := false
	for _, w := range a.Warnings {
		found = found || w == msg
	}
	if !found || a.SIP {
		t.Fatalf("astropy did not drop the SIP set (sip %v, warnings %q)", a.SIP, a.Warnings)
	}
}

// directiveFloat is the float of a "# NAME: value" directive line.
func directiveFloat(t *testing.T, re *regexp.Regexp, text string, n int) []float64 {
	t.Helper()
	m := re.FindStringSubmatch(text)
	if m == nil || len(m) != n+1 {
		t.Fatalf("case has no %q line", re.String())
	}
	var v []float64
	for _, s := range m[1:] {
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			t.Fatal(err)
		}
		v = append(v, f)
	}
	return v
}

// samePoints fails unless two readings were made on the same pixel list.
func samePoints(t *testing.T, a, b *hdrGolden) {
	t.Helper()
	if len(a.P2S) != len(b.P2S) {
		t.Fatalf("%d and %d pixel points", len(a.P2S), len(b.P2S))
	}
	for i := range a.P2S {
		if a.P2S[i].X != b.P2S[i].X || a.P2S[i].Y != b.P2S[i].Y {
			t.Fatalf("pixel point %d: (%g, %g) and (%g, %g)", i, a.P2S[i].X, a.P2S[i].Y, b.P2S[i].X, b.P2S[i].Y)
		}
	}
}

// wcslibOnCase is the built hdr2sky's reading of a case file that wcslib
// must accept, every pixel point computed.
func wcslibOnCase(t *testing.T, caseFile, text string) *hdrGolden {
	t.Helper()
	r := hdr2sky(t, caseFile)
	if r.Exit != 0 {
		t.Fatalf("hdr2sky refused %s (exit %d): %s", caseFile, r.Exit, r.Stderr)
	}
	g, err := decodeHdrGolden(r.Stdout, text)
	if err != nil {
		t.Fatalf("hdr2sky output: %v", err)
	}
	if g.WcslibVersion != hdrWcslibVersion {
		t.Fatalf("hdr2sky linked wcslib %s, want %s", g.WcslibVersion, hdrWcslibVersion)
	}
	for _, p := range g.P2S {
		if p.Status != 0 {
			t.Fatalf("wcslib rejected pixel (%g, %g) (status %d); every point must be computed", p.X, p.Y, p.Status)
		}
	}
	return g
}

// caseTransform builds our transform from a case file's cards.
func caseTransform(t *testing.T, f string) (*Transform, string) {
	t.Helper()
	cards, alt, text := readCaseCards(t, f)
	w, err := hdrParse(cards, alt)
	if err != nil {
		t.Fatalf("wcs.ParseAlt: %v", err)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatalf("transform.New: %v", err)
	}
	return tr, text
}

// TestHdrDeliberateDifferences: headers wcslib 8.10 ACCEPTS and computes
// differently from us, on purpose (hdrcases_differ/). Today every case is
// a RA---TPV header with SIP keywords and no "-SIP" suffix: we ignore the
// SIP set, as astropy 8.0.2 does. For each case, with "# SAME-AS: <name>"
// naming the hdrcases case whose wcslib golden is the answer (the same
// header without the SIP cards, same points):
//
//   - our transform of the case reproduces that golden in both
//     directions, to its tolerances;
//   - astropy 8.0.2's reading of the case (astropy_hdr.py: all_pix2world /
//     all_world2pix, a recorded reading) says it dropped the SIP set and
//     agrees with the golden to ASTROPY-TOL-ARCSEC, and with our
//     sky->pixel to ASTROPY-TOL-PX;
//   - wcslib's reading of the case as written (hdr2sky, recorded) computes
//     every point and lies WCSLIB-DIFFERS-ARCSEC lo..hi from the golden at
//     every point: the difference is real and its size is pinned.
func TestHdrDeliberateDifferences(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(hdrDifferDir, "*.case"))
	if len(files) == 0 {
		t.Fatalf("no cases in %s", hdrDifferDir)
	}
	sort.Strings(files)
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".case")
		t.Run(name, func(t *testing.T) {
			tr, text := caseTransform(t, f)
			same := sameAsRe.FindStringSubmatch(text)
			if same == nil {
				t.Fatalf("%s has no '# SAME-AS: <case>' line", f)
			}
			g := loadHdrGolden(t, same[1])
			if !tr.hdr.SIPIgnored {
				t.Errorf("wcs.Parse did not report the SIP set as ignored (SIPIgnored)")
			}
			tolAS, tolPx := hdrDefaultTolArcsec, hdrDefaultTolPx
			if g.TolArcsec != nil {
				tolAS = *g.TolArcsec
			}
			if g.TolPx != nil {
				tolPx = *g.TolPx
			}
			p := hdrCheckP2S(g, tr, tolAS)
			s := hdrCheckS2P(g, tr, tolPx)
			for _, m := range append(p.mismatches, s.mismatches...) {
				t.Error(m)
			}
			t.Logf("ours = wcslib on %s: p2s max %.3g arcsec, s2p max %.3g px over %d+%d points", same[1], p.maxErr, s.maxErr, p.n, s.n)

			// astropy 8.0.2 on the case as written.
			a := astropyHdrReading(t, f)
			astropyDroppedSIP(t, a)
			aTolAS := directiveFloat(t, astropyTolASRe, text, 1)[0]
			aTolPx := directiveFloat(t, astropyTolPxRe, text, 1)[0]
			if len(a.P2S) != len(g.P2S) || len(a.S2P) != len(a.P2S) {
				t.Fatalf("astropy gave %d/%d points, the golden has %d", len(a.P2S), len(a.S2P), len(g.P2S))
			}
			var maxA, maxOurs, maxPx float64
			for i, q := range a.P2S {
				gp := g.P2S[i]
				if q[0] != gp.X || q[1] != gp.Y || gp.Status != 0 {
					t.Fatalf("astropy point %d (%g, %g), golden (%g, %g) status %d", i, q[0], q[1], gp.X, gp.Y, gp.Status)
				}
				sep := angSepArcsec(q[2], q[3], *gp.Lng, *gp.Lat)
				maxA = math.Max(maxA, sep)
				if !(sep <= aTolAS) {
					t.Errorf("pixel (%g, %g): astropy (%.15g, %.15g), wcslib on %s (%.15g, %.15g): %.3g arcsec > %g", q[0], q[1], q[2], q[3], same[1], *gp.Lng, *gp.Lat, sep, aTolAS)
				}
				lng, lat, err := hdrPixelToSky(g, tr, q[0], q[1])
				if err != nil {
					t.Errorf("pixel (%g, %g): %v", q[0], q[1], err)
					continue
				}
				sep = angSepArcsec(q[2], q[3], lng, lat)
				maxOurs = math.Max(maxOurs, sep)
				if !(sep <= aTolAS) {
					t.Errorf("pixel (%g, %g): ours (%.15g, %.15g), astropy (%.15g, %.15g): %.3g arcsec > %g", q[0], q[1], lng, lat, q[2], q[3], sep, aTolAS)
				}
			}
			for _, q := range a.S2P {
				x, y, err := hdrSkyToPixel(g, tr, q[0], q[1])
				if err != nil {
					t.Errorf("sky (%.15g, %.15g): %v", q[0], q[1], err)
					continue
				}
				d := math.Max(math.Abs(x-q[2]), math.Abs(y-q[3]))
				maxPx = math.Max(maxPx, d)
				if !(d <= aTolPx) {
					t.Errorf("sky (%.15g, %.15g): ours (%.15g, %.15g), astropy (%.15g, %.15g): %.3g px > %g", q[0], q[1], x, y, q[2], q[3], d, aTolPx)
				}
			}
			t.Logf("astropy %s (wcslib %s) on the case: dropped the SIP set; vs wcslib on %s max %.3g arcsec, vs ours p2s max %.3g arcsec, s2p max %.3g px",
				a.Astropy, a.Wcslib, same[1], maxA, maxOurs, maxPx)

			// wcslib 8.10 on the case as written: both distortions.
			lohi := directiveFloat(t, wcslibDiffersRe, text, 2)
			wg := wcslibOnCase(t, f, text)
			samePoints(t, wg, g)
			lo, hi := math.Inf(1), 0.0
			for i, wp := range wg.P2S {
				sep := angSepArcsec(*wp.Lng, *wp.Lat, *g.P2S[i].Lng, *g.P2S[i].Lat)
				lo, hi = math.Min(lo, sep), math.Max(hi, sep)
			}
			if !(lo >= lohi[0] && hi <= lohi[1]) {
				t.Errorf("wcslib on the case as written is %.6g..%.6g arcsec from %s, the case pins %g..%g", lo, hi, same[1], lohi[0], lohi[1])
			}
			t.Logf("wcslib %s on the case as written (SIP and TPV both applied): %.6g..%.6g arcsec from %s (pinned %g..%g)", wg.WcslibVersion, lo, hi, same[1], lohi[0], lohi[1])
		})
	}
}

// atlasCorners are the four corner pixels of the 100x100 ATLAS stamp.
var atlasCorners = [][2]float64{{1, 1}, {100, 1}, {1, 100}, {100, 100}}

// p2sAt is a golden's pixel->sky output at pixel (x, y).
func p2sAt(t *testing.T, g *hdrGolden, x, y float64) (float64, float64) {
	t.Helper()
	for _, p := range g.P2S {
		if p.X == x && p.Y == y && p.Status == 0 {
			return *p.Lng, *p.Lat
		}
	}
	t.Fatalf("golden has no accepted point (%g, %g)", x, y)
	return 0, 0
}

// atlasLinearTAN is the ATLAS header's linear part alone (RA---TAN, the
// same CRPIX/CRVAL/CD, no PV and no SIP card): its sky->pixel turns a
// distorted position into the pixel offset the distortion put there.
func atlasLinearTAN(t *testing.T) (*Transform, *hdrGolden) {
	t.Helper()
	gS := loadHdrGolden(t, "real_atlas_k2_tansip")
	var cards []string
	for _, c := range gS.Cards {
		if regexp.MustCompile(`^(A|B|AP|BP)_(ORDER|[0-9]_[0-9])[ =]`).MatchString(c) {
			continue
		}
		cards = append(cards, strings.Replace(c, "-TAN-SIP'", "-TAN'    ", 1))
	}
	if len(cards) != 12 {
		t.Fatalf("linear ATLAS header has %d cards, want 12 (2 COMMENT, CTYPE, CRPIX, CRVAL, CD): %q", len(cards), cards)
	}
	w, err := hdrParse(cards, "")
	if err != nil {
		t.Fatal(err)
	}
	if w.SIP != nil || len(w.PV) != 0 || w.CelestialProjCode() != "TAN" || w.DistCode[0] != "" {
		t.Fatalf("linear ATLAS header is not plain TAN: %v", w)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	return tr, gS
}

// TestATLASSIPCopiesTPV: the SIP set of the real ATLAS header is a copy of
// its TPV correction. wcslib 8.10's reading of the header with the SIP
// cards removed (real_atlas_k2_tpv, TPV only) and of the TAN-SIP header
// made from the same CRPIX/CRVAL/CD and the SIP set alone
// (real_atlas_k2_tansip, no PV) agree at every point of the stamp grid
// (12 x 12 incl. the corners, and the centre) to 2.4e-7 arcsec
// (measured 2.31e-7, i.e. 1.2e-7 px at 1.86 arcsec/px). At the corners
// the pixel offset each correction makes (its position read back through
// the header's linear TAN part) is (-12.28, 6.06) px at (1, 1), Daryl
// Manning's figure, -12.03..-13.14 px and 5.63..6.33 px over the four
// corners (pinned to 5e-6 px), and the two corrections agree there to
// 1.5e-7 px (measured 1.2e-7).
func TestATLASSIPCopiesTPV(t *testing.T) {
	gT := loadHdrGolden(t, "real_atlas_k2_tpv")
	gS := loadHdrGolden(t, "real_atlas_k2_tansip")
	samePoints(t, gT, gS)
	const budgetAS = 2.4e-7
	var maxSep float64
	for i, p := range gT.P2S {
		q := gS.P2S[i]
		if p.Status != 0 || q.Status != 0 {
			t.Fatalf("wcslib rejected pixel (%g, %g)", p.X, p.Y)
		}
		maxSep = math.Max(maxSep, angSepArcsec(*p.Lng, *p.Lat, *q.Lng, *q.Lat))
	}
	if !(maxSep <= budgetAS) {
		t.Errorf("wcslib TPV-only and SIP-only readings of the ATLAS header differ by up to %.3g arcsec > %g: the SIP set is not a copy of TPV", maxSep, budgetAS)
	}
	t.Logf("wcslib: TPV only vs SIP only, %d points: max %.3g arcsec", len(gT.P2S), maxSep)

	lin, _ := atlasLinearTAN(t)
	offset := func(g *hdrGolden, x, y float64) (float64, float64) {
		a, d := p2sAt(t, g, x, y)
		px, py, err := lin.SkyToPixel(a, d)
		if err != nil {
			t.Fatal(err)
		}
		return px + 1 - x, py + 1 - y
	}
	// Pinned corner offsets (px), measured on the wcslib goldens through
	// the linear TAN inverse (itself checked against wcslib to 1e-9 px).
	want := map[[2]float64][2]float64{
		{1, 1}:     {-12.278215, 6.060993},
		{100, 1}:   {-13.143110, 6.325089},
		{1, 100}:   {-12.025883, 5.626913},
		{100, 100}: {-12.883547, 5.876910},
	}
	for _, c := range atlasCorners {
		tx, ty := offset(gT, c[0], c[1])
		sx, sy := offset(gS, c[0], c[1])
		w := want[c]
		t.Logf("corner (%g, %g): TPV offset (%.6f, %.6f) px, SIP offset (%.6f, %.6f) px, |TPV - SIP| %.2g px",
			c[0], c[1], tx, ty, sx, sy, math.Hypot(tx-sx, ty-sy))
		if d := math.Hypot(tx-sx, ty-sy); !(d <= 1.5e-7) {
			t.Errorf("corner (%g, %g): TPV and SIP offsets %.3g px apart", c[0], c[1], d)
		}
		if math.Abs(tx-w[0]) > 5e-6 || math.Abs(ty-w[1]) > 5e-6 {
			t.Errorf("corner (%g, %g): TPV offset (%.6f, %.6f) px, pinned (%.6f, %.6f)", c[0], c[1], tx, ty, w[0], w[1])
		}
	}
	// Daryl Manning's figure for the corner (1, 1), to its 2 decimals.
	if tx, ty := offset(gT, 1, 1); math.Abs(tx+12.28) > 0.005 || math.Abs(ty-6.06) > 0.005 {
		t.Errorf("corner (1, 1): TPV offset (%.4f, %.4f) px, reported (-12.28, 6.06)", tx, ty)
	}
}

// TestATLASTPVPlusRedundantSIP: the real ATLAS header as written
// (RA---TPV, TPV PLUS its SIP copy, no -SIP). wcslib 8.10 applies both
// corrections, astropy 8.0.2 drops the SIP set, and the two readings are
// 24.45-26.84 arcsec apart over the stamp grid (pinned 24.4..26.9),
// i.e. wcslib's position lands 13.28-14.59 px (pinned 13.2..14.6) from
// the pixel through the TPV solution: the correction counted twice (at
// the corner (1, 1) wcslib's position is the TPV solution's pixel
// (-11.28, 7.06), the pixel moved by a second (-12.28, 6.06) px).
// astropy's reading equals wcslib's reading of the TPV-only header
// (real_atlas_k2_tpv) to 1e-8 arcsec (measured 3.6e-10), and ours of the
// header as written equals astropy's to 1e-8 arcsec (measured 3.7e-10),
// so ours is not doubled.
func TestATLASTPVPlusRedundantSIP(t *testing.T) {
	f := filepath.Join(hdrDifferDir, "real_atlas_k2_tpv_sip.case")
	tr, text := caseTransform(t, f)
	gT := loadHdrGolden(t, "real_atlas_k2_tpv")
	trT, err := hdrTransform(gT)
	if err != nil {
		t.Fatal(err)
	}
	a := astropyHdrReading(t, f)
	astropyDroppedSIP(t, a)
	wg := wcslibOnCase(t, f, text)
	samePoints(t, wg, gT)
	if len(a.P2S) != len(wg.P2S) {
		t.Fatalf("astropy %d points, wcslib %d", len(a.P2S), len(wg.P2S))
	}
	const tolAS = 1e-8
	type span struct{ lo, hi float64 }
	asSep, pxSep := span{math.Inf(1), 0}, span{math.Inf(1), 0}
	var maxAT, maxOurs float64
	for i, q := range a.P2S {
		x, y := q[0], q[1]
		w := wg.P2S[i]
		if w.X != x || w.Y != y {
			t.Fatalf("point %d: astropy (%g, %g), wcslib (%g, %g)", i, x, y, w.X, w.Y)
		}
		// 2b: wcslib (both) vs astropy (TPV only).
		s := angSepArcsec(*w.Lng, *w.Lat, q[2], q[3])
		asSep = span{math.Min(asSep.lo, s), math.Max(asSep.hi, s)}
		px, py, err := hdrSkyToPixel(gT, trT, *w.Lng, *w.Lat)
		if err != nil {
			t.Fatal(err)
		}
		d := math.Hypot(px-x, py-y)
		pxSep = span{math.Min(pxSep.lo, d), math.Max(pxSep.hi, d)}
		// 2c: astropy on the full header = wcslib on the TPV-only header.
		tl, tb := p2sAt(t, gT, x, y)
		maxAT = math.Max(maxAT, angSepArcsec(q[2], q[3], tl, tb))
		// Ours on the full header = astropy.
		ol, ob, err := hdrPixelToSky(gT, tr, x, y)
		if err != nil {
			t.Fatal(err)
		}
		so := angSepArcsec(ol, ob, q[2], q[3])
		maxOurs = math.Max(maxOurs, so)
		if !(so <= tolAS) {
			t.Errorf("pixel (%g, %g): ours (%.15g, %.15g), astropy (%.15g, %.15g): %.3g arcsec > %g", x, y, ol, ob, q[2], q[3], so, tolAS)
		}
		if sw := angSepArcsec(ol, ob, *w.Lng, *w.Lat); !(sw >= 24.4) {
			t.Errorf("pixel (%g, %g): ours is %.3g arcsec from wcslib's doubled reading; want the TPV-only answer (>= 24.4 arcsec away)", x, y, sw)
		}
	}
	t.Logf("2b: wcslib (TPV + SIP) vs astropy (TPV), %d points: %.4f..%.4f arcsec, %.4f..%.4f px", len(a.P2S), asSep.lo, asSep.hi, pxSep.lo, pxSep.hi)
	t.Logf("2c: astropy on the header as written vs wcslib on the TPV-only header: max %.3g arcsec", maxAT)
	t.Logf("ours on the header as written vs astropy: max %.3g arcsec", maxOurs)
	if !(asSep.lo >= 24.4 && asSep.hi <= 26.9) {
		t.Errorf("wcslib's doubled reading is %.4f..%.4f arcsec from astropy's, pinned 24.4..26.9", asSep.lo, asSep.hi)
	}
	if !(pxSep.lo >= 13.2 && pxSep.hi <= 14.6) {
		t.Errorf("wcslib's doubled reading is %.4f..%.4f px from the pixel through TPV, pinned 13.2..14.6", pxSep.lo, pxSep.hi)
	}
	if !(maxAT <= tolAS) {
		t.Errorf("astropy on the header as written is %.3g arcsec from wcslib on the TPV-only header > %g", maxAT, tolAS)
	}
	// At a corner: the doubled position sits a second (-12.28, 6.06) px
	// correction away (wcslib applies SIP to the pixel offset, then TPV).
	px, py, _ := hdrSkyToPixel(gT, trT, *wg.P2S[0].Lng, *wg.P2S[0].Lat)
	t.Logf("corner (1, 1): wcslib's doubled position is at pixel (%.4f, %.4f) of the TPV solution", px, py)
	if math.Abs(px+11.278215) > 1e-5 || math.Abs(py-7.060993) > 1e-5 {
		t.Errorf("corner (1, 1): wcslib's doubled position at TPV pixel (%.6f, %.6f), pinned (-11.278215, 7.060993)", px, py)
	}
}

// TestTPVSIPSuffixAppliesBoth: RA---TPV-SIP declares both distortions, and
// both apply, SIP then TPV, as in wcslib 8.10, which accepts the CTYPE
// (its golden tpv_sip_suffix_atlas, checked point for point by
// TestHdrGoldens). The case is the real ATLAS TPV solution plus a SIP set
// that is not a copy of it: wcslib's answer (and ours) must move off the
// TPV-only answer (real_atlas_k2_tpv) by exactly the SIP correction: read
// back through the TPV solution, wcslib's position is the pixel plus
// (A(u, v), B(u, v)) to 1e-8 px (measured 1.1e-10 px; SIP is applied to
// the pixel offset before the matrix and TPV), 2.33-3.0 px on the stamp.
func TestTPVSIPSuffixAppliesBoth(t *testing.T) {
	gB := loadHdrGolden(t, "tpv_sip_suffix_atlas")
	gT := loadHdrGolden(t, "real_atlas_k2_tpv")
	samePoints(t, gB, gT)
	trB, err := hdrTransform(gB)
	if err != nil {
		t.Fatal(err)
	}
	if trB.sip == nil || trB.tpv == nil {
		t.Fatalf("RA---TPV-SIP: sip %v, tpv %v; both must apply", trB.sip != nil, trB.tpv != nil)
	}
	trT, err := hdrTransform(gT)
	if err != nil {
		t.Fatal(err)
	}
	p := hdrCheckP2S(gB, trB, hdrDefaultTolArcsec)
	for _, m := range p.mismatches {
		t.Error(m)
	}
	var minShift, worst = math.Inf(1), 0.0
	for _, q := range gB.P2S {
		if q.Status != 0 {
			t.Fatalf("wcslib rejected (%g, %g)", q.X, q.Y)
		}
		px, py, err := hdrSkyToPixel(gT, trT, *q.Lng, *q.Lat)
		if err != nil {
			t.Fatal(err)
		}
		u, v := q.X-gB.crpix(t, 0), q.Y-gB.crpix(t, 1)
		du, dv := trB.sip.Forward(u, v)
		du, dv = du-u, dv-v
		shift := math.Hypot(px-q.X, py-q.Y)
		minShift = math.Min(minShift, shift)
		d := math.Hypot(px-q.X-du, py-q.Y-dv)
		worst = math.Max(worst, d)
		if !(d <= 1e-8) {
			t.Errorf("pixel (%g, %g): wcslib moved (%.4f, %.4f) px off TPV, the SIP correction is (%.4f, %.4f)", q.X, q.Y, px-q.X, py-q.Y, du, dv)
		}
	}
	if !(minShift >= 2) {
		t.Errorf("wcslib's TPV-SIP answer is only %.3g px from TPV alone somewhere: the SIP set is not applied", minShift)
	}
	t.Logf("wcslib RA---TPV-SIP vs TPV alone: shift >= %.3f px, = the SIP vector to %.2g px; ours = wcslib to %.3g arcsec", minShift, worst, p.maxErr)
}

// crpix is CRPIX of pixel axis i+1 from a golden's cards.
func (g *hdrGolden) crpix(t *testing.T, i int) float64 {
	t.Helper()
	w, err := hdrParse(g.Cards, g.Alt)
	if err != nil {
		t.Fatal(err)
	}
	return w.CRPix[i]
}
