package transform

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"regexp"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
	"github.com/dmarkham/fits/wcs"
)

// Pipeline-level checks on the wcslib header goldens beyond the
// coordinate comparison in TestHdrGoldens: what wcslib DERIVED from the
// header (WCS dimension, celestial axis indices, LONPOLE/LATPOLE after
// the native-pole solution), which headers wcslib REFUSES, and the
// integrity of the real archive headers.

// TestHdrDerivedParameters: for every golden, wcs.ParseAlt must find the
// same number of WCS axes (WCSAXESa / NAXIS / highest keyword index) and
// the same celestial axes as wcslib's wcspih, and New must derive the same
// LONPOLE and LATPOLE as wcslib's celset (wcsprm.cel.ref[2], ref[3]).
// Tolerance 1e-9 deg: both sides run the same closed-form native-pole
// solution in float64.
func TestHdrDerivedParameters(t *testing.T) {
	for _, name := range hdrCaseNames(t) {
		t.Run(name, func(t *testing.T) {
			g := loadHdrGolden(t, name)
			w, err := hdrParse(g.Cards, g.Alt)
			if err != nil {
				t.Fatalf("parse: %v", err)
			}
			if w.NAxis != g.WCS.Naxis {
				t.Errorf("NAxis = %d, wcslib naxis = %d", w.NAxis, g.WCS.Naxis)
			}
			if w.LonAxis-1 != g.WCS.Lng || w.LatAxis-1 != g.WCS.Lat {
				t.Errorf("celestial axes (lon %d, lat %d), wcslib (lng %d, lat %d) (0-based)", w.LonAxis-1, w.LatAxis-1, g.WCS.Lng, g.WCS.Lat)
			}
			tr, err := New(w)
			if err != nil {
				// Every committed hdrcase is one wcslib accepts; refusing
				// it is a failure here too, not a reason to skip.
				t.Fatalf("New: %v (wcslib accepts this header)", err)
			}
			if g.WCS.CelRef[2] == nil || g.WCS.CelRef[3] == nil {
				t.Fatal("golden has no cel_ref")
			}
			if d := math.Abs(tr.lonpole - *g.WCS.CelRef[2]); !(d <= 1e-9) {
				t.Errorf("LONPOLE = %.15g, wcslib %.15g", tr.lonpole, *g.WCS.CelRef[2])
			}
			if d := math.Abs(tr.latpole - *g.WCS.CelRef[3]); !(d <= 1e-9) {
				t.Errorf("LATPOLE = %.15g, wcslib %.15g", tr.latpole, *g.WCS.CelRef[3])
			}
		})
	}
}

// hdrParse runs card text through header.ParseCards -> wcs.ParseAlt.
func hdrParse(cards []string, alt string) (*wcs.Header, error) {
	var buf bytes.Buffer
	for _, c := range cards {
		if len(c) > header.CardWidth {
			return nil, fmt.Errorf("card longer than 80: %q", c)
		}
		buf.WriteString(fmt.Sprintf("%-80s", c))
	}
	buf.WriteString(fmt.Sprintf("%-80s", "END"))
	for buf.Len()%2880 != 0 {
		buf.WriteByte(' ')
	}
	parsed, _, err := header.ParseCards(buf.Bytes())
	if err != nil {
		return nil, err
	}
	return wcs.ParseAlt(header.FromCards(parsed), alt)
}

// readCaseCards returns the HEADER cards and the ALT letter of a case
// file (same rules as hdr2sky.c: blank lines skipped, END stops).
func readCaseCards(t *testing.T, path string) (cards []string, alt string, text string) {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	text = string(b)
	in := false
	for _, line := range strings.Split(text, "\n") {
		line = strings.TrimRight(line, "\r")
		if in {
			if line == "" {
				continue
			}
			if strings.HasPrefix(line, "END") && strings.TrimSpace(line) == "END" {
				break
			}
			cards = append(cards, line)
			continue
		}
		s := strings.TrimSpace(line)
		if s == "HEADER" {
			in = true
		} else if strings.HasPrefix(s, "ALT ") {
			alt = strings.TrimSpace(s[4:])
		}
	}
	if len(cards) == 0 {
		t.Fatalf("%s has no HEADER cards", path)
	}
	return cards, alt, text
}

const hdrRejectDir = "testdata/wref/hdrcases_reject"

var expectWcslibRe = regexp.MustCompile(`(?m)^# EXPECT-WCSLIB-ERROR: (.+)$`)

// TestHdrRejectCases: headers that wcslib 8.10 REFUSES (wcsset fails)
// must also be refused by wcs.ParseAlt or transform.New, never turned into
// coordinates. Each case in hdrcases_reject/ names the wcslib error text
// it expects; the wcslib side is the built hdr2sky's reading of the case
// (recorded in tier 1, run live in tier 2: testutil.Run).
func TestHdrRejectCases(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(hdrRejectDir, "*.case"))
	if len(files) == 0 {
		t.Fatalf("no reject cases in %s", hdrRejectDir)
	}
	sort.Strings(files)
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".case")
		t.Run(name, func(t *testing.T) {
			cards, alt, text := readCaseCards(t, f)
			m := expectWcslibRe.FindStringSubmatch(text)
			if m == nil {
				t.Fatalf("%s has no '# EXPECT-WCSLIB-ERROR: ...' line", f)
			}
			w, err := hdrParse(cards, alt)
			if err == nil {
				_, err = New(w)
			}
			if err == nil {
				t.Errorf("wcslib refuses this header (%q) but we built a transform", m[1])
			} else {
				t.Logf("we refuse: %v", err)
			}

			r := hdr2sky(t, f)
			out := append(r.Stdout, r.Stderr...)
			if r.Exit == 0 {
				t.Fatalf("hdr2sky ACCEPTED %s: the case no longer documents a wcslib rejection", name)
			}
			if !strings.Contains(string(out), m[1]) {
				t.Errorf("hdr2sky failed, but without %q:\n%s", m[1], out)
			}
		})
	}
}

const hdrStricterDir = "testdata/wref/hdrcases_stricter"

var expectGoRe = regexp.MustCompile(`(?m)^# EXPECT-GO-ERROR: (.+)$`)

// TestHdrStricterThanWcslib: headers wcslib 8.10 ACCEPTS (and computes
// coordinates for) that we refuse on purpose because the header is
// ambiguous (rule 7). The wcslib side proves the refusal really is a
// policy choice and not wcslib behaviour: the built hdr2sky must succeed
// on the case exactly as written (no harness-injected cards) and compute
// pixel->sky points. The Go side must refuse with the quoted error text.
// Interop as in TestHdrRejectCases.
func TestHdrStricterThanWcslib(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(hdrStricterDir, "*.case"))
	if len(files) == 0 {
		t.Fatalf("no cases in %s", hdrStricterDir)
	}
	sort.Strings(files)
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".case")
		t.Run(name, func(t *testing.T) {
			cards, alt, text := readCaseCards(t, f)
			m := expectGoRe.FindStringSubmatch(text)
			if m == nil {
				t.Fatalf("%s has no '# EXPECT-GO-ERROR: ...' line", f)
			}
			w, err := hdrParse(cards, alt)
			if err == nil {
				_, err = New(w)
			}
			if err == nil {
				t.Errorf("we built a transform; this header is refused by policy (%q)", m[1])
			} else if !strings.Contains(err.Error(), m[1]) {
				t.Errorf("refused, but error %q does not mention %q", err, m[1])
			}

			r := hdr2sky(t, f)
			if r.Exit != 0 {
				t.Fatalf("hdr2sky REFUSED %s (exit status %d): the case no longer documents a header wcslib accepts\n%s", name, r.Exit, r.Stderr)
			}
			g, err := decodeHdrGolden(r.Stdout, text)
			if err != nil {
				t.Fatalf("hdr2sky output: %v", err)
			}
			if len(g.HarnessInjected) != 0 {
				t.Fatalf("hdr2sky injected %q; the case must carry every card wcslib needs", g.HarnessInjected)
			}
			ok := 0
			for _, p := range g.P2S {
				if p.Status == 0 {
					ok++
				}
			}
			if ok == 0 || ok != len(g.P2S) {
				t.Fatalf("wcslib computed %d of %d pixel->sky points; want all", ok, len(g.P2S))
			}
		})
	}
}

const hdrLenientDir = "testdata/wref/hdrcases_lenient"

var sameAsRe = regexp.MustCompile(`(?m)^# SAME-AS: (\S+)$`)

// TestHdrLenientThanWcslib: headers wcslib 8.10 REFUSES that we accept on
// purpose, because wcslib's refusal is an artefact (card order, or an
// inverse polynomial that is only ever a starting guess) and an
// equivalent header has a wcslib answer. Each case in hdrcases_lenient/
// names the wcslib error text and, with "# SAME-AS: <name>", the
// hdrcases golden whose answer it must reproduce: both directions,
// including which points are rejected, to that golden's tolerances. The
// wcslib side (the built hdr2sky must FAIL on the case as written, with
// the quoted text) is interop as in TestHdrRejectCases.
func TestHdrLenientThanWcslib(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(hdrLenientDir, "*.case"))
	if len(files) == 0 {
		t.Fatalf("no cases in %s", hdrLenientDir)
	}
	sort.Strings(files)
	for _, f := range files {
		name := strings.TrimSuffix(filepath.Base(f), ".case")
		t.Run(name, func(t *testing.T) {
			cards, alt, text := readCaseCards(t, f)
			m := expectWcslibRe.FindStringSubmatch(text)
			same := sameAsRe.FindStringSubmatch(text)
			if m == nil || same == nil {
				t.Fatalf("%s needs '# EXPECT-WCSLIB-ERROR: ...' and '# SAME-AS: <case>' lines", f)
			}
			g := loadHdrGolden(t, same[1])
			w, err := hdrParse(cards, alt)
			var tr *Transform
			if err == nil {
				tr, err = New(w)
			}
			if err != nil {
				t.Fatalf("we refuse this header (%v); it must give the answer of %s", err, same[1])
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
			for _, mm := range append(p.mismatches, s.mismatches...) {
				t.Error(mm)
			}
			if p.n == 0 || s.n == 0 {
				t.Fatalf("golden %s has no points", same[1])
			}
			t.Logf("same as %s: p2s max %.3g arcsec, s2p max %.3g px over %d+%d points", same[1], p.maxErr, s.maxErr, p.n, s.n)

			r := hdr2sky(t, f)
			out := append(r.Stdout, r.Stderr...)
			if r.Exit == 0 {
				t.Fatalf("hdr2sky ACCEPTED %s: the case no longer documents a wcslib refusal", name)
			}
			if !strings.Contains(string(out), m[1]) {
				t.Errorf("hdr2sky failed, but without %q:\n%s", m[1], out)
			}
		})
	}
}

// toolHdr2sky is the wcslib header harness testdata/wref/hdr2sky (`make
// -C wcs/transform/testdata/wref hdr2sky` builds it for tier 2).
var toolHdr2sky = testutil.Tool{Name: "hdr2sky", Harness: "wcs/transform/testdata/wref/hdr2sky", Sources: []string{"wcs/transform/testdata/wref/hdr2sky.c"}}

// hdr2sky is wcslib's reading of a header case file.
func hdr2sky(t *testing.T, caseFile string) testutil.Result {
	t.Helper()
	return testutil.Run(t, toolHdr2sky, testutil.Call{Args: []string{caseFile}, Inputs: []string{caseFile}})
}

// TestTPVMissingPV1DefaultsToOne: the TPV registry says a missing PV1_1 /
// PV2_1 is 1. Golden tpv_explicit_pv1 is wcslib's answer for the header
// WITH explicit PV1_1 = PV2_1 = 1; we delete those two cards and must get
// the same coordinates (wcslib itself would treat the edited header as
// PVi_1 = 0, which is why its golden is made from the explicit form).
func TestTPVMissingPV1DefaultsToOne(t *testing.T) {
	g := loadHdrGolden(t, "tpv_explicit_pv1")
	var cards []string
	removed := 0
	for _, c := range g.Cards {
		if strings.HasPrefix(c, "PV1_1   =") || strings.HasPrefix(c, "PV2_1   =") {
			removed++
			continue
		}
		cards = append(cards, c)
	}
	if removed != 2 {
		t.Fatalf("expected to remove PV1_1 and PV2_1, removed %d cards", removed)
	}
	w, err := hdrParse(cards, "")
	if err != nil {
		t.Fatal(err)
	}
	tr, err := New(w)
	if err != nil {
		t.Fatal(err)
	}
	p := hdrCheckP2S(g, tr, hdrDefaultTolArcsec)
	s := hdrCheckS2P(g, tr, *g.TolPx)
	for _, m := range append(p.mismatches, s.mismatches...) {
		t.Error(m)
	}
	t.Logf("PV1_1/PV2_1 deleted: p2s max %.3g arcsec, s2p max %.3g px over %d+%d points", p.maxErr, s.maxErr, p.n, s.n)
}

var cardsShaRe = regexp.MustCompile(`(?m)^# header cards sha256: ([0-9a-f]{64})$`)

// TestRealHdrCasesIntact: every real_*.case (in hdrcases and
// hdrcases_differ) was written by
// realhdr/make_cases.py, which records the sha256 of the verbatim header
// cards it copied from the archive file (whose own sha256 is pinned in the
// script). The cards must still hash to that value. This catches an
// accidental edit of the cards; it cannot catch a deliberate edit of both
// the cards and the hash line (the pin against the archive bytes is the
// source sha256 in make_cases.py, which needs the downloaded sources).
func TestRealHdrCasesIntact(t *testing.T) {
	files, _ := filepath.Glob(filepath.Join(hdrCaseDir, "real_*.case"))
	differ, _ := filepath.Glob(filepath.Join(hdrDifferDir, "real_*.case"))
	files = append(files, differ...)
	if len(files) < 13 {
		t.Fatalf("expected >= 13 real header cases, found %d", len(files))
	}
	for _, f := range files {
		b, err := os.ReadFile(f)
		if err != nil {
			t.Fatal(err)
		}
		text := string(b)
		m := cardsShaRe.FindStringSubmatch(text)
		if m == nil {
			t.Errorf("%s: no '# header cards sha256' line", f)
			continue
		}
		i := strings.Index(text, "\nHEADER\n")
		if i < 0 {
			t.Errorf("%s: no HEADER line", f)
			continue
		}
		sum := sha256.Sum256([]byte(text[i+len("\nHEADER\n"):]))
		if got := hex.EncodeToString(sum[:]); got != m[1] {
			t.Errorf("%s: header cards sha256 %s, recorded %s (edited by hand? regenerate with realhdr/make_cases.py)", f, got, m[1])
		}
	}
}
