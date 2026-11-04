package transform

import (
	"bytes"
	"crypto/sha256"
	"encoding/base64"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"errors"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
	"github.com/dmarkham/fits/wcs"
)

// Full-pipeline wcslib truth tests.
//
// Every case file in testdata/wref/hdrcases/<name>.case holds raw FITS
// header card text plus pixel and sky points. testdata/wref/hdr2sky.c runs
// that header through wcslib 8.10 (wcspih -> wcsset -> wcsp2s / wcss2p) and
// writes testdata/wref/golden/hdr/<name>.json (committed gzipped as
// <name>.json.gz, read with testutil.ReadGolden). This test parses the SAME
// cards with header.ParseCards -> wcs.ParseAlt -> New and requires, for
// every point, in both directions:
//
//   - wcslib rejected the point  => we return an error (off-projection
//     points must never come back as coordinates);
//   - wcslib accepted the point  => we succeed and agree within tolerance.
//
// To add a case: drop a .case file in hdrcases/, run
//
//	make -C wcs/transform/testdata/wref hdr-goldens
//
// and commit both files. Nothing here needs editing unless the case
// exposes a known, not-yet-fixed bug (see hdrKnownBad).
//
// Missing or stale goldens are a test FAILURE, never a skip.

const (
	hdrCaseDir   = "testdata/wref/hdrcases"
	hdrGoldenDir = "testdata/wref/golden/hdr"

	// The wcslib version the committed goldens were produced with (pinned
	// in testdata/wref/PROVENANCE.md).
	hdrWcslibVersion = "8.10"

	// Default pixel->sky budget: 1e-8 arcsec (~4.8e-14 rad, a few hundred
	// ULP of a radian-scale angle). Both sides evaluate the same
	// closed-form trigonometry in float64; the smoke cases measure at most
	// 2.2e-10 arcsec (TAN 768 px off-axis, TAN-SIP 1.7e-10), so this leaves
	// ~45x headroom for libm differences while the smallest real bug in the
	// audit (~0.01 arcsec) is six orders of magnitude above it. A case whose
	// math is iterative on one side may set TOL_ARCSEC in its case file,
	// and must say why in a comment there.
	hdrDefaultTolArcsec = 1e-8
	// Default sky->pixel budget: 1e-9 px. Smoke cases measure at most
	// 3.7e-11 px (pixel coords ~1e3, ULP ~1e-13). Override: TOL_PX.
	hdrDefaultTolPx = 1e-9
)

// hdrKnownBad marks case directions that are KNOWN to disagree with
// wcslib because of an open, audited bug. Keyed by case name, then by
// direction ("p2s" or "s2p"); the value is the audit reference.
//
// A known-bad direction that unexpectedly AGREES with wcslib fails the
// test, so a fix cannot land without deleting its entry here, and an
// entry for a case that no longer exists also fails.
var hdrKnownBad = map[string]map[string]string{}

// hdrS2PBeyondWcslib lists cases where wcslib's sky->pixel iteration fails
// on points that its own pixel->sky accepted (a documented difference, see
// New): there a sky point that came from an accepted pixel ("from": "p2s")
// and that wcslib rejects must instead give back that source pixel, to
// the case's TOL_PX. The truth is wcslib's forward direction; at least
// minPoints such points must be present so the case keeps testing it.
var hdrS2PBeyondWcslib = map[string]struct {
	minPoints int
	why       string
}{
	"dss_far_outside_plate": {3, "dis.c disx2p fails to converge 6000 px outside a degree-5 plate solution; Newton converges"},
}

// hdrGolden is one decoded golden/hdr/<name>.json together with the
// inputs regenerated from its case file. The committed file (format
// "hdr2sky-2", written by hdr2sky.c) holds only what the case file does not
// already say: wcslib's derived parameters, digests of the case text, of
// the cards and of the input coordinates, and per point wcslib's status,
// stat and output coordinates as exact big-endian binary64 (base64). The
// input pixel / sky coordinates and the cards are rebuilt from the case
// file exactly as hdr2sky builds them, and the digests prove the rebuild
// is bit-identical to what wcslib was given.
type hdrGolden struct {
	WcslibVersion   string
	Case            string // case file text
	Note            string
	Alt             string
	TolArcsec       *float64
	TolPx           *float64
	Cards           []string // 80-column cards as handed to wcspih (without harness_injected)
	HarnessInjected []string
	WCS             hdrGoldenWCS
	P2S             []hdrP2S
	S2P             []hdrS2P
}

type hdrGoldenWCS struct {
	Naxis   int        `json:"naxis"`
	Lng     int        `json:"lng"`
	Lat     int        `json:"lat"`
	Ctype   []string   `json:"ctype"`
	Crval   []float64  `json:"crval"`
	Lonpole float64    `json:"lonpole"`
	Latpole float64    `json:"latpole"`
	CelRef  []*float64 `json:"cel_ref"`
}

// hdrP2S is one pixel -> sky point: input (X, Y) 1-based FITS pixels,
// wcslib's status/stat and, when Status == 0, its (Lng, Lat) in degrees.
type hdrP2S struct {
	X, Y         float64
	Status, Stat int
	Lng, Lat     *float64
}

// hdrS2P is one sky -> pixel point; From is "case" for an explicit SKY /
// SKYGRID point and "p2s" for the world coordinate of an accepted pixel
// point.
type hdrS2P struct {
	From         string
	Lng, Lat     float64
	Status, Stat int
	X, Y         *float64
}

// hdrGoldenFile is the committed JSON (format "hdr2sky-2").
type hdrGoldenFile struct {
	Format          string       `json:"format"`
	WcslibVersion   string       `json:"wcslib_version"`
	CaseSHA256      string       `json:"case_sha256"`
	CardsSHA256     string       `json:"cards_sha256"`
	NCards          int          `json:"ncards"`
	Note            string       `json:"note"`
	Alt             string       `json:"alt"`
	TolArcsec       *float64     `json:"tol_arcsec"`
	TolPx           *float64     `json:"tol_px"`
	HarnessInjected []string     `json:"harness_injected"`
	WCS             hdrGoldenWCS `json:"wcs"`
	NPix            int          `json:"npix"`
	NSkyCase        int          `json:"nsky_case"`
	InputsSHA256    string       `json:"inputs_sha256"`
	Sample          struct {
		P2S []struct {
			I      int      `json:"i"`
			X      float64  `json:"x"`
			Y      float64  `json:"y"`
			Status int      `json:"status"`
			Lng    *float64 `json:"lng"`
			Lat    *float64 `json:"lat"`
		} `json:"p2s"`
		S2P []struct {
			I      int      `json:"i"`
			Lng    float64  `json:"lng"`
			Lat    float64  `json:"lat"`
			Status int      `json:"status"`
			X      *float64 `json:"x"`
			Y      *float64 `json:"y"`
		} `json:"s2p"`
	} `json:"sample"`
	P2S string `json:"p2s"`
	S2P string `json:"s2p"`
}

// hdrCaseInput is what hdr2sky.c reads from a case file.
type hdrCaseInput struct {
	pixX, pixY []float64 // GRID / PIX points, in directive order
	skyA, skyD []float64 // SKYGRID / SKY points, in directive order
	alt        string
	cards      []string // 80 columns each
}

// hdrLin is hdr2sky.c's lin(): n evenly spaced values v0..v1, evaluated
// with the same operations in the same order (no product feeds an
// addition directly, so no FMA can change it).
func hdrLin(v0, v1 float64, n, i int) float64 {
	if n == 1 {
		return v0
	}
	return v0 + (v1-v0)*float64(i)/float64(n-1)
}

// parseHdrCase parses a case file exactly as hdr2sky.c does: directives
// before HEADER ('#' comments and blank lines skipped), then one card per
// line (empty lines skipped, padded to 80 columns, stopping at END).
func parseHdrCase(text string) (*hdrCaseInput, error) {
	c := &hdrCaseInput{}
	inHdr := false
	nums := func(f []string, n int) ([]float64, error) {
		if len(f) != n {
			return nil, fmt.Errorf("want %d values, got %q", n, f)
		}
		v := make([]float64, n)
		for k, x := range f {
			var err error
			if v[k], err = strconv.ParseFloat(x, 64); err != nil {
				return nil, err
			}
		}
		return v, nil
	}
	for _, line := range strings.Split(text, "\n") {
		if k := strings.IndexByte(line, '\r'); k >= 0 {
			line = line[:k]
		}
		if inHdr {
			if line == "" {
				continue
			}
			if len(line) > header.CardWidth {
				return nil, fmt.Errorf("header card longer than 80 chars: %q", line)
			}
			card := fmt.Sprintf("%-80s", line)
			if strings.HasPrefix(card, "END     ") {
				break
			}
			c.cards = append(c.cards, card)
			continue
		}
		p := strings.TrimLeft(line, " \t")
		if p == "" || p[0] == '#' {
			continue
		}
		f := strings.Fields(p)
		switch {
		case strings.HasPrefix(p, "HEADER"):
			inHdr = true
		case f[0] == "GRID" || f[0] == "SKYGRID":
			v, err := nums(f[1:], 6)
			if err != nil {
				return nil, fmt.Errorf("%s: %v", f[0], err)
			}
			na, nb := int(v[2]), int(v[5])
			for j := 0; j < nb; j++ {
				for i := 0; i < na; i++ {
					a, b := hdrLin(v[0], v[1], na, i), hdrLin(v[3], v[4], nb, j)
					if f[0] == "GRID" {
						c.pixX, c.pixY = append(c.pixX, a), append(c.pixY, b)
					} else {
						c.skyA, c.skyD = append(c.skyA, a), append(c.skyD, b)
					}
				}
			}
		case f[0] == "PIX" || f[0] == "SKY":
			v, err := nums(f[1:], 2)
			if err != nil {
				return nil, fmt.Errorf("%s: %v", f[0], err)
			}
			if f[0] == "PIX" {
				c.pixX, c.pixY = append(c.pixX, v[0]), append(c.pixY, v[1])
			} else {
				c.skyA, c.skyD = append(c.skyA, v[0]), append(c.skyD, v[1])
			}
		case f[0] == "ALT" && len(f) > 1:
			c.alt = f[1][:1]
		case f[0] == "TOL_ARCSEC", f[0] == "TOL_PX", f[0] == "NOTE":
		default:
			return nil, fmt.Errorf("unknown directive %q", p)
		}
	}
	if len(c.cards) == 0 {
		return nil, errors.New("case file has no HEADER cards")
	}
	return c, nil
}

func sha256Hex(b []byte) string {
	h := sha256.Sum256(b)
	return hex.EncodeToString(h[:])
}

// beFloats is the canonical byte encoding of the digests: big-endian
// IEEE-754 binary64.
func beFloats(v ...float64) []byte {
	b := make([]byte, 8*len(v))
	for i, x := range v {
		binary.BigEndian.PutUint64(b[8*i:], math.Float64bits(x))
	}
	return b
}

// hdrRecord is one packed point: status, stat and, when status == 0, two
// binary64 coordinates.
type hdrRecord struct {
	status, stat int
	a, b         float64
}

func decodeHdrRecords(s string, n int) ([]hdrRecord, error) {
	b, err := base64.StdEncoding.DecodeString(s)
	if err != nil {
		return nil, err
	}
	out := make([]hdrRecord, 0, n)
	for len(b) > 0 {
		if len(b) < 2 {
			return nil, errors.New("truncated record")
		}
		r := hdrRecord{status: int(b[0]), stat: int(b[1])}
		b = b[2:]
		if r.status == 0 {
			if len(b) < 16 {
				return nil, errors.New("truncated coordinates")
			}
			r.a = math.Float64frombits(binary.BigEndian.Uint64(b))
			r.b = math.Float64frombits(binary.BigEndian.Uint64(b[8:]))
			b = b[16:]
		}
		out = append(out, r)
	}
	if len(out) != n {
		return nil, fmt.Errorf("%d records, want %d", len(out), n)
	}
	return out, nil
}

// sameFloat compares a sample value with a decoded one bit for bit (both
// are exact: %.17g round-trips binary64). hdr2sky prints a non-finite
// value (wcslib can return status 0 with NaN, e.g. sky->pixel through a
// singular CD) as null.
func sameFloat(a *float64, b float64) bool {
	if a == nil {
		return math.IsNaN(b) || math.IsInf(b, 0)
	}
	return math.Float64bits(*a) == math.Float64bits(b)
}

// decodeHdrGolden checks a golden against its case text and rebuilds the
// full point lists. Every digest, count and sample must match.
func decodeHdrGolden(data []byte, caseText string) (*hdrGolden, error) {
	var f hdrGoldenFile
	if err := json.Unmarshal(data, &f); err != nil {
		return nil, err
	}
	if f.Format != "hdr2sky-2" {
		return nil, fmt.Errorf("golden format %q, want hdr2sky-2", f.Format)
	}
	if got := sha256Hex([]byte(caseText)); got != f.CaseSHA256 {
		return nil, fmt.Errorf("stale: it was generated from a different case file (case sha256 %s, golden %s)", got, f.CaseSHA256)
	}
	in, err := parseHdrCase(caseText)
	if err != nil {
		return nil, fmt.Errorf("case file: %v", err)
	}
	if len(in.cards) != f.NCards || sha256Hex([]byte(strings.Join(in.cards, ""))) != f.CardsSHA256 {
		return nil, fmt.Errorf("cards rebuilt from the case (%d) differ from the %d hdr2sky gave wcspih", len(in.cards), f.NCards)
	}
	if alt := f.Alt; alt != in.alt {
		return nil, fmt.Errorf("ALT %q in the golden, %q in the case", alt, in.alt)
	}
	if len(in.pixX) != f.NPix || len(in.skyA) != f.NSkyCase {
		return nil, fmt.Errorf("case gives %d pixel / %d sky points, golden %d / %d", len(in.pixX), len(in.skyA), f.NPix, f.NSkyCase)
	}
	var inBytes []byte
	for i := range in.pixX {
		inBytes = append(inBytes, beFloats(in.pixX[i], in.pixY[i])...)
	}
	for i := range in.skyA {
		inBytes = append(inBytes, beFloats(in.skyA[i], in.skyD[i])...)
	}
	if got := sha256Hex(inBytes); got != f.InputsSHA256 {
		return nil, fmt.Errorf("input coordinates rebuilt from the case differ from hdr2sky's (sha256 %s, golden %s)", got, f.InputsSHA256)
	}
	p2s, err := decodeHdrRecords(f.P2S, f.NPix)
	if err != nil {
		return nil, fmt.Errorf("p2s records: %v", err)
	}
	g := &hdrGolden{
		WcslibVersion: f.WcslibVersion, Case: caseText, Note: f.Note, Alt: f.Alt,
		TolArcsec: f.TolArcsec, TolPx: f.TolPx, Cards: in.cards,
		HarnessInjected: f.HarnessInjected, WCS: f.WCS,
	}
	skyA, skyD := append([]float64{}, in.skyA...), append([]float64{}, in.skyD...)
	for i, r := range p2s {
		p := hdrP2S{X: in.pixX[i], Y: in.pixY[i], Status: r.status, Stat: r.stat}
		if r.status == 0 {
			lng, lat := r.a, r.b
			p.Lng, p.Lat = &lng, &lat
			skyA, skyD = append(skyA, lng), append(skyD, lat)
		}
		g.P2S = append(g.P2S, p)
	}
	s2p, err := decodeHdrRecords(f.S2P, len(skyA))
	if err != nil {
		return nil, fmt.Errorf("s2p records: %v", err)
	}
	for i, r := range s2p {
		s := hdrS2P{From: "case", Lng: skyA[i], Lat: skyD[i], Status: r.status, Stat: r.stat}
		if i >= f.NSkyCase {
			s.From = "p2s"
		}
		if r.status == 0 {
			x, y := r.a, r.b
			s.X, s.Y = &x, &y
		}
		g.S2P = append(g.S2P, s)
	}
	for _, sp := range f.Sample.P2S {
		if sp.I < 0 || sp.I >= len(g.P2S) {
			return nil, fmt.Errorf("sample p2s index %d out of range", sp.I)
		}
		p := g.P2S[sp.I]
		if p.X != sp.X || p.Y != sp.Y || p.Status != sp.Status || (p.Status == 0 && !(sameFloat(sp.Lng, *p.Lng) && sameFloat(sp.Lat, *p.Lat))) {
			return nil, fmt.Errorf("sample p2s[%d] disagrees with the packed records", sp.I)
		}
	}
	for _, ss := range f.Sample.S2P {
		if ss.I < 0 || ss.I >= len(g.S2P) {
			return nil, fmt.Errorf("sample s2p index %d out of range", ss.I)
		}
		s := g.S2P[ss.I]
		if s.Lng != ss.Lng || s.Lat != ss.Lat || s.Status != ss.Status || (s.Status == 0 && !(sameFloat(ss.X, *s.X) && sameFloat(ss.Y, *s.Y))) {
			return nil, fmt.Errorf("sample s2p[%d] disagrees with the packed records", ss.I)
		}
	}
	return g, nil
}

// hdrCaseNames returns the sorted case names in hdrCaseDir. An empty or
// missing directory is fatal.
func hdrCaseNames(t *testing.T) []string {
	t.Helper()
	files, err := filepath.Glob(filepath.Join(hdrCaseDir, "*.case"))
	if err != nil {
		t.Fatal(err)
	}
	if len(files) == 0 {
		t.Fatalf("no case files in %s: the wcslib header-pipeline truth set is missing", hdrCaseDir)
	}
	names := make([]string, len(files))
	for i, f := range files {
		names[i] = strings.TrimSuffix(filepath.Base(f), ".case")
	}
	sort.Strings(names)
	return names
}

// loadHdrGolden reads golden/hdr/<name>.json and checks it was generated
// from the current case file by the pinned wcslib version.
func loadHdrGolden(t *testing.T, name string) *hdrGolden {
	t.Helper()
	caseText, err := os.ReadFile(filepath.Join(hdrCaseDir, name+".case"))
	if err != nil {
		t.Fatalf("read case: %v", err)
	}
	path := filepath.Join(hdrGoldenDir, name+".json")
	data, err := testutil.ReadGolden(path)
	if err != nil {
		t.Fatalf("golden for case %q missing (%v): run make -C wcs/transform/testdata/wref hdr-goldens and commit it", name, err)
	}
	g, err := decodeHdrGolden(data, string(caseText))
	if err != nil {
		t.Fatalf("%s: %v; regenerate with make hdr-goldens", path, err)
	}
	if g.WcslibVersion != hdrWcslibVersion {
		t.Fatalf("%s generated by wcslib %q, PROVENANCE pins %q", path, g.WcslibVersion, hdrWcslibVersion)
	}
	if len(g.Cards) == 0 || len(g.P2S)+len(g.S2P) == 0 {
		t.Fatalf("%s has no cards or no points", path)
	}
	return g
}

// hdrTransform builds a Transform from the golden's card text through the
// public parsing path: header.ParseCards -> wcs.ParseAlt -> New.
func hdrTransform(g *hdrGolden) (*Transform, error) {
	var buf bytes.Buffer
	for _, c := range g.Cards {
		if len(c) != header.CardWidth {
			return nil, fmt.Errorf("golden card is %d bytes, want %d: %q", len(c), header.CardWidth, c)
		}
		buf.WriteString(c)
	}
	buf.WriteString(fmt.Sprintf("%-80s", "END"))
	for buf.Len()%2880 != 0 {
		buf.WriteByte(' ')
	}
	cards, _, err := header.ParseCards(buf.Bytes())
	if err != nil {
		return nil, fmt.Errorf("header.ParseCards: %w", err)
	}
	w, err := wcs.ParseAlt(header.FromCards(cards), g.Alt)
	if err != nil {
		return nil, fmt.Errorf("wcs.ParseAlt: %w", err)
	}
	tr, err := New(w)
	if err != nil {
		return nil, fmt.Errorf("transform.New: %w", err)
	}
	return tr, nil
}

// angSepArcsec is the great-circle separation of two (lng, lat) positions
// in degrees, returned in arcsec (Vincenty form; stable at all distances
// and through RA wrap and the poles).
func angSepArcsec(a1, d1, a2, d2 float64) float64 {
	a1, d1, a2, d2 = degToRad(a1), degToRad(d1), degToRad(a2), degToRad(d2)
	da := a2 - a1
	sd1, cd1 := math.Sincos(d1)
	sd2, cd2 := math.Sincos(d2)
	sda, cda := math.Sincos(da)
	num := math.Hypot(cd2*sda, cd1*sd2-sd1*cd2*cda)
	den := sd1*sd2 + cd1*cd2*cda
	return radToDeg(math.Atan2(num, den)) * 3600
}

// hdrCheck is the result of one direction of one case.
type hdrCheck struct {
	n, rejected int
	beyond      int // hdrCheckS2PRoundTrip: wcslib-rejected points checked against their source pixel
	maxErr      float64
	mismatches  []string
	cond        condTally // rows within budget only with their libm sensitivity
}

func hdrCheckP2S(g *hdrGolden, tr *Transform, tol float64) hdrCheck {
	var r hdrCheck
	for _, p := range g.P2S {
		r.n++
		// Golden pixels are 1-based FITS coordinates; PixelToSky is 0-based
		// and takes the celestial pixel axes in axis order (hdrPixelToSky).
		a, d, err := hdrPixelToSky(g, tr, p.X, p.Y)
		if p.Status != 0 {
			r.rejected++
			if err == nil {
				r.mismatches = append(r.mismatches, fmt.Sprintf(
					"pixel (%g, %g): wcslib rejects (status %d stat %d) but we returned (%.12g, %.12g)",
					p.X, p.Y, p.Status, p.Stat, a, d))
			}
			continue
		}
		if err != nil {
			r.mismatches = append(r.mismatches, fmt.Sprintf(
				"pixel (%g, %g): wcslib gives (%.12g, %.12g) but we rejected: %v", p.X, p.Y, *p.Lng, *p.Lat, err))
			continue
		}
		sep := angSepArcsec(*p.Lng, *p.Lat, a, d)
		r.maxErr = math.Max(r.maxErr, sep)
		if !(sep <= tol) {
			r.mismatches = append(r.mismatches, fmt.Sprintf(
				"pixel (%g, %g): got (%.15g, %.15g), wcslib (%.15g, %.15g): %.3g arcsec > %g",
				p.X, p.Y, a, d, *p.Lng, *p.Lat, sep, tol))
		}
	}
	return r
}

func hdrCheckS2P(g *hdrGolden, tr *Transform, tol float64) hdrCheck {
	return hdrCheckS2PRoundTrip(g, tr, tol, false)
}

// hdrCheckS2PRoundTrip is hdrCheckS2P; with beyond set, a point wcslib
// rejects that came from an accepted pixel point must give back that
// pixel (hdrS2PBeyondWcslib); r.beyond counts them.
func hdrCheckS2PRoundTrip(g *hdrGolden, tr *Transform, tol float64, beyond bool) hdrCheck {
	var r hdrCheck
	var src []hdrP2S // accepted pixel points, in the order hdr2sky fed them to s2p
	for _, p := range g.P2S {
		if p.Status == 0 {
			src = append(src, p)
		}
	}
	k := 0
	for _, s := range g.S2P {
		r.n++
		x, y, err := hdrSkyToPixel(g, tr, s.Lng, s.Lat) // 1-based, golden order
		var from *hdrP2S
		if s.From == "p2s" {
			if k < len(src) {
				from = &src[k]
			}
			k++
		}
		if s.Status != 0 && beyond && from != nil {
			r.beyond++
			if err != nil {
				r.mismatches = append(r.mismatches, fmt.Sprintf(
					"sky (%.12g, %.12g) from pixel (%g, %g): wcslib's s2p fails here, and we rejected too: %v", s.Lng, s.Lat, from.X, from.Y, err))
				continue
			}
			dpx := math.Max(math.Abs(x-from.X), math.Abs(y-from.Y))
			r.maxErr = math.Max(r.maxErr, dpx)
			if !(dpx <= tol) {
				r.mismatches = append(r.mismatches, fmt.Sprintf(
					"sky (%.12g, %.12g) from pixel (%g, %g): got pixel (%.15g, %.15g), %.3g px > %g", s.Lng, s.Lat, from.X, from.Y, x, y, dpx, tol))
			}
			continue
		}
		if s.Status != 0 {
			r.rejected++
			if err == nil {
				r.mismatches = append(r.mismatches, fmt.Sprintf(
					"sky (%.12g, %.12g) [%s]: wcslib rejects (status %d stat %d) but we returned pixel (%.12g, %.12g)",
					s.Lng, s.Lat, s.From, s.Status, s.Stat, x, y))
			}
			continue
		}
		if err != nil {
			r.mismatches = append(r.mismatches, fmt.Sprintf(
				"sky (%.12g, %.12g) [%s]: wcslib gives pixel (%.12g, %.12g) but we rejected: %v",
				s.Lng, s.Lat, s.From, *s.X, *s.Y, err))
			continue
		}
		dpx := math.Max(math.Abs(x-*s.X), math.Abs(y-*s.Y))
		r.maxErr = math.Max(r.maxErr, dpx)
		// Beyond the pixel budget only by the point's own libm
		// sensitivity (libm_sens_test.go): next to a divergence of the
		// projection (COP's pole, where pixels run to 1e7 and more) the
		// last bits of Go's and glibc's sin/cos/atan2 differ by more than
		// 1e-7 px. Counted in r.conditioned.
		lng, lat := s.Lng, s.Lat
		sens := func() float64 {
			return libmSensitivity(func() (float64, float64, bool) {
				x, y, err := hdrSkyToPixel(g, tr, lng, lat)
				return x, y, err == nil
			}, func(a0, a1, b0, b1 float64) float64 { return math.Max(math.Abs(a0-b0), math.Abs(a1-b1)) })
		}
		if !r.cond.valueOK(dpx, tol, sens) {
			r.mismatches = append(r.mismatches, fmt.Sprintf(
				"sky (%.12g, %.12g) [%s]: got pixel (%.15g, %.15g), wcslib (%.15g, %.15g): %.3g px > %g + libm sensitivity %.3g",
				s.Lng, s.Lat, s.From, x, y, *s.X, *s.Y, dpx, tol, sens()))
		}
	}
	return r
}

// hdrReport fails or logs one direction, honouring hdrKnownBad.
func hdrReport(t *testing.T, caseName, dir, unit string, r hdrCheck, buildErr error) {
	t.Helper()
	if buildErr != nil {
		r.mismatches = []string{fmt.Sprintf("wcslib accepts this header but we cannot build a transform: %v", buildErr)}
	}
	ref, known := hdrKnownBad[caseName][dir]
	if known {
		if len(r.mismatches) == 0 {
			t.Errorf("KNOWN-BAD case %s/%s now AGREES with wcslib (%d points, max err %.3g %s): "+
				"the bug is fixed, delete its hdrKnownBad entry (%s)", caseName, dir, r.n, r.maxErr, unit, ref)
			return
		}
		t.Logf("KNOWN-BAD (still failing, as recorded): %s", ref)
		t.Logf("  %d of %d points disagree; first: %s", len(r.mismatches), r.n, r.mismatches[0])
		return
	}
	const show = 15
	for i, m := range r.mismatches {
		if i == show {
			t.Errorf("... and %d more", len(r.mismatches)-show)
			break
		}
		t.Error(m)
	}
	if len(r.mismatches) == 0 {
		t.Logf("%d points (%d rejected by wcslib and by us), max err %.3g %s", r.n, r.rejected, r.maxErr, unit)
		if r.cond.conditioned > 0 {
			t.Logf("%d points within budget only with their libm sensitivity (worst used %.2g of it)", r.cond.conditioned, r.cond.worstRatio)
		}
	}
}

func TestHdrGoldens(t *testing.T) {
	names := hdrCaseNames(t)

	// Every golden has a case file, every known-bad entry names a case.
	have := map[string]bool{}
	for _, n := range names {
		have[n] = true
	}
	goldens := testutil.GlobGoldens(filepath.Join(hdrGoldenDir, "*.json"))
	for _, gp := range goldens {
		if n := strings.TrimSuffix(filepath.Base(gp), ".json"); !have[n] {
			t.Errorf("orphan golden %s has no %s/%s.case", gp, hdrCaseDir, n)
		}
	}
	for n := range hdrS2PBeyondWcslib {
		if !have[n] {
			t.Errorf("hdrS2PBeyondWcslib entry %q names no case file", n)
		}
	}
	for n, dirs := range hdrKnownBad {
		if !have[n] {
			t.Errorf("hdrKnownBad entry %q names no case file", n)
		}
		for d := range dirs {
			if d != "p2s" && d != "s2p" {
				t.Errorf("hdrKnownBad[%q] has bad direction %q", n, d)
			}
		}
	}

	for _, name := range names {
		t.Run(name, func(t *testing.T) {
			g := loadHdrGolden(t, name)
			tolAS, tolPx := hdrDefaultTolArcsec, hdrDefaultTolPx
			if g.TolArcsec != nil {
				tolAS = *g.TolArcsec
			}
			if g.TolPx != nil {
				tolPx = *g.TolPx
			}
			if len(g.HarnessInjected) > 0 {
				t.Logf("wcslib was given extra cards %q (hdr2sky SIP-without-AP/BP workaround); we parse the original header", g.HarnessInjected)
			}
			tr, buildErr := hdrTransform(g)
			t.Run("p2s", func(t *testing.T) {
				var r hdrCheck
				if buildErr == nil {
					r = hdrCheckP2S(g, tr, tolAS)
				}
				hdrReport(t, name, "p2s", "arcsec", r, buildErr)
			})
			t.Run("s2p", func(t *testing.T) {
				var r hdrCheck
				bw, beyond := hdrS2PBeyondWcslib[name]
				if buildErr == nil {
					r = hdrCheckS2PRoundTrip(g, tr, tolPx, beyond)
				}
				if beyond && r.beyond < bw.minPoints {
					t.Errorf("%s: %d sky points where wcslib's s2p fails, want >= %d (%s)", name, r.beyond, bw.minPoints, bw.why)
				} else if beyond {
					t.Logf("%d points where wcslib's s2p fails (%s) give back their source pixel", r.beyond, bw.why)
				}
				hdrReport(t, name, "s2p", "px", r, buildErr)
			})
		})
	}
}
