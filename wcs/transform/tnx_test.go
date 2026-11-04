package transform

import (
	"fmt"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
)

// tnxCards builds a parsed header from 80-column card text (raw bytes, as
// read from a file), so string values go through the header layer's
// trailing-blank trimming exactly as on disk.
func tnxCards(t *testing.T, lines ...string) *header.Header {
	t.Helper()
	var buf []byte
	for _, l := range append(lines, "END") {
		if len(l) > 80 {
			t.Fatalf("card too long: %q", l)
		}
		buf = append(buf, (l + strings.Repeat(" ", 80-len(l)))...)
	}
	cards, _, err := header.ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	return header.FromCards(cards)
}

var tnxBase = []string{
	"NAXIS   =                    2",
	"NAXIS1  =                  100",
	"NAXIS2  =                  100",
	"CTYPE1  = 'RA---TNX'",
	"CTYPE2  = 'DEC--TNX'",
	"CRVAL1  =                150.0",
	"CRVAL2  =                  2.0",
	"CRPIX1  =                 50.0",
	"CRPIX2  =                 50.0",
	"CD1_1   =              -7.0E-5",
	"CD1_2   =                  0.0",
	"CD2_1   =                  0.0",
	"CD2_2   =               7.0E-5",
}

func tnxHeader(t *testing.T, wat1, wat2 string, extra ...string) *header.Header {
	lines := append(append([]string(nil), tnxBase...), extra...)
	if wat1 != "" {
		lines = append(lines, "WAT1_001= '"+wat1+"'")
	}
	if wat2 != "" {
		lines = append(lines, "WAT2_001= '"+wat2+"'")
	}
	return tnxCards(t, lines...)
}

const (
	tnxGoodLng = `wtype=tnx axtype=ra lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`
	tnxGoodLat = `wtype=tnx axtype=dec latcor = "3. 2. 2. 1. -1 1 -1 1 0 -1E-4 0 3E-4"`
)

// TestTNXRejects: every malformed or unsupported TNX header is an error
// from NewTNX (never an uncorrected or guessed transform).
func TestTNXRejects(t *testing.T) {
	if _, err := NewTNX(tnxHeader(t, tnxGoodLng, tnxGoodLat)); err != nil {
		t.Fatalf("control header rejected: %v", err)
	}
	cases := []struct {
		name, wat1, wat2 string
		extra            []string
		want             string
	}{
		{"zpx wtype", strings.Replace(tnxGoodLng, "wtype=tnx", "wtype=zpx", 1), tnxGoodLat, nil, "wtype=zpx"},
		{"too few coeffs", `wtype=tnx lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4"`, tnxGoodLat, nil, "needs 4 coefficients, got 3"},
		{"too many coeffs", `wtype=tnx lngcor = "3. 2. 2. 0. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "needs 3 coefficients, got 4"},
		{"type 4", `wtype=tnx lngcor = "4. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "unknown surface type 4"},
		{"cross 3", `wtype=tnx lngcor = "3. 2. 2. 3. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "unknown cross-term code 3"},
		{"fractional order", `wtype=tnx lngcor = "3. 2.5 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "must be an integer"},
		{"order 0", `wtype=tnx lngcor = "3. 0. 2. 1. -1 1 -1 1"`, tnxGoodLat, nil, "orders must be >= 1"},
		{"empty domain", `wtype=tnx lngcor = "1. 2. 2. 1. 1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "empty domain"},
		{"bad number", `wtype=tnx lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 x 0"`, tnxGoodLat, nil, "bad number"},
		{"unclosed quote", `wtype=tnx lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0`, tnxGoodLat, nil, "no closing quote"},
		{"no equals", `wtype=tnx lngcor "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "has no '='"},
		{"ro", `ro=1. lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "ro=1."},
		// 57.3 is within 1e-5 of 180/pi: only a tight tolerance refuses it.
		{"ro near default", `ro=57.3 lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "ro=57.3"},
		{"ro in WAT2", tnxGoodLng, `ro=57.3 latcor = "3. 2. 2. 1. -1 1 -1 1 0 -1E-4 0 3E-4"`, nil, "WAT2 ro=57.3"},
		// WCSTools computes a zpx WAT2 as if it were tnx; it is refused here.
		{"zpx wtype in WAT2", tnxGoodLng, strings.Replace(tnxGoodLat, "wtype=tnx", "wtype=zpx", 1), nil, "WAT2 declares wtype=zpx"},
		{"crota without cd", tnxGoodLng, tnxGoodLat, nil, "CROTA2"},
		{"crota1 without cd", tnxGoodLng, tnxGoodLat, nil, "CROTA1"},
		// PV1_0 = 1 with PV1_1/PV1_2 moves the fiducial point (wcslib
		// prjoff); WCSTools ignores PVi_m on TNX. Refused like every
		// longitude-axis PVi_m (TestTNXRefusedHeaders has the WCSTools
		// answers for PV1_1/PV1_2 and PV1_3).
		{"pv1_0 fiducial offset", tnxGoodLng, tnxGoodLat,
			[]string{"PV1_0   =                  1.0", "PV1_1   =                  0.0", "PV1_2   =                 89.0"}, "PV1_0, PV1_1, PV1_2"},
		{"pv1_4 latpole", tnxGoodLng, tnxGoodLat, []string{"PV1_4   =                 10.0"}, "PV1_4"},
		{"duplicate attribute", `wtype=tnx wtype=tnx lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, `attribute "wtype" given twice`},
		{"duplicate lngcor", `wtype=tnx lngcor = "3. 1. 1. 0. -1 1 -1 1 0" `, tnxGoodLat,
			[]string{`WAT1_002= 'lngcor = "3. 1. 1. 0. -1 1 -1 1 1E-4"'`}, `attribute "lngcor" given twice`},
		{"nan coefficient", `wtype=tnx lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 NaN 0"`, tnxGoodLat, nil, `bad number "NaN"`},
		{"inf coefficient", `wtype=tnx lngcor = "3. 2. 2. 1. -1 1 -1 1 0 1E-4 2E-4 +Inf"`, tnxGoodLat, nil, `bad number "+Inf"`},
		{"inf domain", `wtype=tnx lngcor = "3. 2. 2. 1. -Inf 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, `bad number "-Inf"`},
		{"huge order", `wtype=tnx lngcor = "3. 2000. 2. 1. -1 1 -1 1 0 1E-4 2E-4 0"`, tnxGoodLat, nil, "must be an integer"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			var h *header.Header
			if strings.HasPrefix(c.name, "crota") {
				var lines []string
				for _, l := range tnxBase {
					if !strings.HasPrefix(l, "CD") {
						lines = append(lines, l)
					}
				}
				lines = append(lines, "CDELT1  =              -7.0E-5", "CDELT2  =               7.0E-5",
					c.want+"  =                 12.0", "WAT1_001= '"+c.wat1+"'", "WAT2_001= '"+c.wat2+"'")
				h = tnxCards(t, lines...)
			} else {
				h = tnxHeader(t, c.wat1, c.wat2, c.extra...)
			}
			_, err := NewTNX(h)
			if err == nil {
				t.Fatalf("NewTNX accepted it")
			}
			if !strings.Contains(err.Error(), c.want) {
				t.Fatalf("error %q does not mention %q", err, c.want)
			}
		})
	}
	// Not a TNX header at all.
	var tan []string
	for _, l := range tnxBase {
		tan = append(tan, strings.Replace(l, "TNX", "TAN", 1))
	}
	h := tnxCards(t, append(tan, "WAT1_001= '"+tnxGoodLng+"'")...)
	if _, err := NewTNX(h); err == nil {
		t.Error("NewTNX accepted a TAN header")
	}
	// A mixed pair (RA with GLAT) is not one celestial frame.
	var mixed []string
	for _, l := range tnxBase {
		mixed = append(mixed, strings.Replace(l, "'DEC--TNX'", "'GLAT-TNX'", 1))
	}
	if _, err := NewTNX(tnxCards(t, append(mixed, "WAT1_001= '"+tnxGoodLng+"'")...)); err == nil || !strings.Contains(err.Error(), "not a matching celestial pair") {
		t.Errorf("NewTNX on RA---TNX/GLAT-TNX: err = %v, want a 'not a matching celestial pair' error", err)
	}
	// DEC on axis 1, RA on axis 2: WCSTools handles it, NewTNX does not
	// (inherited from New) and must say so rather than mis-map the axes.
	swapped := tnxCards(t,
		"NAXIS   =                    2", "NAXIS1  =                 2048", "NAXIS2  =                 2048",
		"CTYPE1  = 'DEC--TNX'", "CTYPE2  = 'RA---TNX'",
		"CRVAL1  =           2.20583333", "CRVAL2  =         150.11916667",
		"CRPIX1  =               1024.5", "CRPIX2  =              1030.25",
		"CD1_1   =              3.9E-06", "CD1_2   =             7.25E-05",
		"CD2_1   =             -7.2E-05", "CD2_2   =              4.1E-06",
		"WAT1_001= '"+tnxGoodLat+"'", "WAT2_001= '"+tnxGoodLng+"'")
	if _, err := NewTNX(swapped); err == nil || !strings.Contains(err.Error(), "celestial axes must be 1,2") {
		t.Errorf("NewTNX on DEC/RA axis order: err = %v, want a 'celestial axes must be 1,2' error", err)
	}
}

// TestConcatWATRefusesWhatWCSToolsDrops: WCSTools hgetm stops at the
// first missing WATn_k and cuts the joined string at 1999 characters, and
// reads no key number above 499, so the cards after a gap or past the
// limit would be dropped without a word. NewTNX must refuse each of
// these headers (rule 7), naming the cards it would otherwise lose.
func TestConcatWATRefusesWhatWCSToolsDrops(t *testing.T) {
	lng := tnxGoodLng + " "
	lat := `wtype=tnx axtype=dec latcor = "3. 2. 2. 1. -1 1 -1 1 0 -1E-4 0 3E-4"`
	// Long WAT1 for the length limit: a 10x10 full-cross Legendre lngcor
	// whose joined string is over 1999 characters.
	var long strings.Builder
	long.WriteString(`wtype=tnx axtype=ra lngcor = "2. 10. 10. 1. -0.09 0.09 -0.09 0.09`)
	for i := range 100 {
		fmt.Fprintf(&long, " %.15e", 1e-6*float64(i%7-3))
	}
	long.WriteString(`"`)
	if long.Len() <= 1999 {
		t.Fatalf("long WAT1 is only %d characters", long.Len())
	}
	chunks := func(axis int, s string) []string {
		var out []string
		for i, n := 0, 1; i < len(s); i, n = i+68, n+1 {
			out = append(out, fmt.Sprintf("WAT%d_%03d= '%s'", axis, n, s[i:min(i+68, len(s))]))
		}
		return out
	}
	join := func(parts ...[]string) []string {
		lines := append([]string(nil), tnxBase...)
		for _, p := range parts {
			lines = append(lines, p...)
		}
		return lines
	}
	w1 := chunks(1, lng+strings.Repeat(" ", 70)+` latcor = "3. 2. 2. 1. -1 1 -1 1 0 -1E-4 0 3E-4"`)
	if len(w1) < 3 {
		t.Fatalf("need >= 3 WAT1 chunks, got %d", len(w1))
	}
	cases := []struct {
		name  string
		lines []string
		want  string
	}{
		{"gap", join(w1[:1], w1[2:], chunks(2, lat)), "WAT1_003"},
		{"gap before last", join(w1[:len(w1)-2], w1[len(w1)-1:], chunks(2, lat)), fmt.Sprintf("WAT1_%03d", len(w1))},
		{"mixed spelling", join(w1[:2], []string{"WAT1_3  = '" + w1[2][11:]}, chunks(2, lat)), "WAT1_3"},
		{"two first-card spellings", join(w1, []string{"WAT1_1  = '" + w1[0][11:]}, chunks(2, lat)), "first card in more than one spelling (WAT1_1, WAT1_001)"},
		{"no first card", join([]string{"WAT1_002= '" + lng + "'"}, chunks(2, lat)), "WAT1_002"},
		{"gap in WAT2", join(chunks(1, lng), chunks(2, lat)[:1], []string{"WAT2_003= 'axtype=dec'"}), "WAT2_003"},
		{"key number 500", join(chunks(1, lng), chunks(2, lat), func() []string {
			var l []string
			for i := 2; i <= 500; i++ {
				l = append(l, fmt.Sprintf("WAT2_%03d= ''", i))
			}
			return l
		}()), "WAT2_500"},
		{"over 1999 characters", join(chunks(1, long.String()), chunks(2, lat)), "reads at most 1999"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			_, err := NewTNX(tnxCards(t, c.lines...))
			if err == nil {
				t.Fatalf("NewTNX accepted it")
			}
			if !strings.Contains(err.Error(), c.want) {
				t.Fatalf("error %q does not mention %q", err, c.want)
			}
			t.Logf("refused: %v", err)
		})
	}
	// Control: the same cards without the gap are accepted.
	if _, err := NewTNX(tnxCards(t, join(w1, chunks(2, lat))...)); err != nil {
		t.Fatalf("control header without a gap: %v", err)
	}
	// The boundary: the committed syn_leg_wat_1999 header (WAT1 exactly
	// 1999 characters, proven against WCSTools by its golden) is accepted;
	// one more blank makes it 2000 characters, which WCSTools cuts to
	// "longpole=175." and NewTNX must refuse.
	h1999, text := readCardFile(t, filepath.Join(tnxHeaderDir, "syn_leg_wat_1999.hdr"))
	if s, err := concatWAT(h1999, "WAT1"); err != nil || len(s) != 1999 || !strings.HasSuffix(s, " longpole=175.5") {
		t.Fatalf("syn_leg_wat_1999 WAT1: %d characters, err %v", len(s), err)
	}
	lines := strings.Split(strings.TrimRight(string(text), "\n"), "\n")
	for i, l := range lines {
		if strings.HasPrefix(l, "WAT1_030= '") {
			lines[i] = "WAT1_030= ' " + l[len("WAT1_030= '"):]
		}
	}
	if _, err := NewTNX(tnxCards(t, lines[:len(lines)-1]...)); err == nil || !strings.Contains(err.Error(), "is 2000 characters") {
		t.Fatalf("2000-character WAT1: err = %v, want a length refusal", err)
	}
	// A WAT card continued with CONTINUE: WCSTools reads only the first
	// card, so the continuation would be lost.
	h := tnxCards(t, join([]string{
		"WAT1_001= 'wtype=tnx axtype=ra lngcor = \"3. 2. 2. 1. -1 1 -1 1 0 1E-4&'",
		"CONTINUE  ' 2E-4 0\"'",
	}, chunks(2, lat))...)
	if _, err := NewTNX(h); err == nil || !strings.Contains(err.Error(), "CONTINUE") {
		t.Fatalf("WAT card with CONTINUE: err = %v, want a CONTINUE error", err)
	}
}

// TestConcatWATKeepsChunkBlanks: the header layer trims trailing blanks
// from string values, but a WAT chunk boundary after a blank separates two
// numbers. Parsed from raw card bytes, the blank must survive.
func TestConcatWATKeepsChunkBlanks(t *testing.T) {
	h := tnxCards(t,
		"WAT1_001= 'wtype=tnx axtype=ra lngcor = \"3. 2. 2. 1. -1 1 -1 1 0. '",
		"WAT1_002= '5 0.25 0.125\"'",
	)
	if v, _ := h.String("WAT1_001"); strings.HasSuffix(v, " ") {
		t.Log("header layer no longer trims string values; concatWAT is still correct")
	}
	tnx, err := ParseTNX(h)
	if err != nil {
		t.Fatal(err)
	}
	want := []float64{0, 5, 0.25, 0.125}
	for i, c := range tnx.Lon.Coeffs {
		if c != want[i] {
			t.Fatalf("coeffs %v, want %v (chunks joined without the trailing blank)", tnx.Lon.Coeffs, want)
		}
	}
	// IRAF also emits the _1 and _01 key spellings that WCSTools accepts.
	h = tnxCards(t,
		"WAT1_1  = 'wtype=tnx lngcor = \"3. 1. 1. 0. -1 1 -1 1 '",
		"WAT1_2  = '7E-5\"'",
	)
	if tnx, err = ParseTNX(h); err != nil || len(tnx.Lon.Coeffs) != 1 || tnx.Lon.Coeffs[0] != 7e-5 {
		t.Fatalf("WAT1_1 spelling: %+v, %v", tnx, err)
	}
	h = tnxCards(t,
		"WAT1_01 = 'wtype=tnx lngcor = \"3. 1. 1. 0. -1 1 -1 1 '",
		"WAT1_02 = '6E-5\"'",
	)
	if tnx, err = ParseTNX(h); err != nil || len(tnx.Lon.Coeffs) != 1 || tnx.Lon.Coeffs[0] != 6e-5 {
		t.Fatalf("WAT1_01 spelling: %+v, %v", tnx, err)
	}
}

// TestWATCardStringDoubledQuote: inside a FITS string value a doubled
// single quote stands for one quote (FITS 4.0 section 4.2.1.1); reading the
// raw card must keep it as one quote, and keep the trailing blanks.
func TestWATCardStringDoubledQuote(t *testing.T) {
	h := tnxCards(t, "WAT1_001= 'a''b c '")
	c, ok := h.Get("WAT1_001")
	if !ok {
		t.Fatal("WAT1_001 missing")
	}
	got, err := watCardString(c)
	if err != nil {
		t.Fatal(err)
	}
	if got != "a'b c " {
		t.Fatalf("watCardString = %q, want %q", got, "a'b c ")
	}
}

// TestTNXInverseNonConvergence: Inverse must report !ok, not a value, when
// Newton cannot converge (here a correction whose Jacobian vanishes).
func TestTNXInverseNonConvergence(t *testing.T) {
	// lngcor = -xi makes xi' = 0 for every xi: singular Jacobian.
	tnx := &TNX{Lon: tnxSurface{Type: 3, XOrder: 2, YOrder: 1, Cross: 1,
		XMin: -1, XMax: 1, YMin: -1, YMax: 1, Coeffs: []float64{0, -1}, present: true}}
	if _, _, ok := tnx.Inverse(0.1, 0.1); ok {
		t.Fatal("Inverse reported ok for a singular correction")
	}
	// Strongly nonlinear: xi' = xi + 10 xi^2 has no real solution for
	// xi' < -1/40; Newton wanders and must give up.
	tnx = &TNX{Lon: tnxSurface{Type: 3, XOrder: 3, YOrder: 1, Cross: 1,
		XMin: -1, XMax: 1, YMin: -1, YMax: 1, Coeffs: []float64{0, 0, 10}, present: true}}
	if _, _, ok := tnx.Inverse(-1, 0); ok {
		t.Fatal("Inverse reported ok where Forward has no preimage")
	}
}

// TestTNXRealHeaderDetails pins facts about the committed real headers
// that the golden comparison relies on, so a change to the sample files
// cannot quietly remove the coverage.
func TestTNXRealHeaderDetails(t *testing.T) {
	for name, want := range map[string]struct {
		typ, cross, ncoef int
		split             string
	}{
		// CTIO 4m Mosaic II (registry tnx.fits): Chebyshev, half cross,
		// chunk boundary inside a mantissa ("0.305" | "7102603054315").
		"noao_ct4m_tnx1": {1, 2, 10, "0.305'"},
		// Registry tnx doc figure 1: polynomial, half cross, chunk ends
		// in a blank that the header layer trims.
		"registry_tnxdoc_tnx3": {3, 2, 10, "-0.1387962673564234 '"},
		// KPNO 4m Mosaic-1 InstCal ccd1 (NOIRLab archive): polynomial,
		// FULL cross terms (16 coefficients), chunk WAT1_003 ends in a
		// blank, and WAT1_006/007 come after WAT2_005 (interleaved).
		"kpno_kp4m_mosaic1_ccd1": {3, 1, 16, "-0.05218531830803252 '"},
	} {
		h, text := readCardFile(t, filepath.Join(tnxHeaderDir, name+".hdr"))
		if !strings.Contains(string(text), want.split) {
			t.Errorf("%s: sample no longer has the chunk split %q", name, want.split)
		}
		if name == "kpno_kp4m_mosaic1_ccd1" {
			w25 := strings.Index(string(text), "\nWAT2_005=")
			w16 := strings.Index(string(text), "\nWAT1_006=")
			if w25 < 0 || w16 < 0 || w16 < w25 {
				t.Errorf("%s: sample no longer has WAT1_006 after WAT2_005 (interleaved WAT cards)", name)
			}
		}
		tnx, err := ParseTNX(h)
		if err != nil {
			t.Fatalf("%s: %v", name, err)
		}
		for _, s := range []tnxSurface{tnx.Lon, tnx.Lat} {
			if s.Type != want.typ || s.Cross != want.cross || s.XOrder != 4 || s.YOrder != 4 || len(s.Coeffs) != want.ncoef {
				t.Errorf("%s: surface %+v", name, s)
			}
		}
	}
}
