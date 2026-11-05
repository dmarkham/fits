package transform

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// TNX truth tests against WCSTools 3.9.7.
//
// Each testdata/tnx/headers/<name>.hdr (FITS cards, one per line) has a
// testdata/tnx/golden/<name>.json written by testdata/tnx/gen_golden.py
// from WCSTools xy2sky / sky2xy, whose tnxpos.c is the C translation of
// IRAF's TNX driver. The test parses the SAME cards with
// header.ParseCards -> NewTNX and requires, for every point:
//
//   - p2s: PixelToSky(pixel) agrees with xy2sky;
//   - s2p: SkyToPixel(xy2sky sky) returns the original pixel (inverse of
//     the WCSTools forward map), and agrees with sky2xy within WCSTools'
//     own inverse convergence.
//
// A missing golden, an orphan golden, or a golden whose header hash does
// not match its .hdr is a FAILURE.

const (
	tnxHeaderDir = "testdata/tnx/headers"
	tnxGoldenDir = "testdata/tnx/golden"

	// p2s budget. xy2sky prints 12 decimals of a degree (rounding <= 5e-13
	// deg per coordinate, <= 2.5e-9 arcsec on the sky); both sides evaluate
	// the same closed-form polynomial + gnomonic formulas in float64.
	// Measured worst case over all headers: 2.5e-9 arcsec, i.e. the print
	// rounding. Budget 1e-8 arcsec; the smallest real TNX bug in scope
	// (type-3 normalisation, LONPOLE default) is > 1e2 arcsec.
	tnxTolArcsec = 1e-8
	// s2p budget against the original pixel: the golden sky is rounded to
	// 5e-13 deg, < 1e-8 px at these plate scales (measured worst 7.5e-9
	// px). Budget 1e-7 px.
	tnxTolPx = 1e-7
	// s2p budget against sky2xy: WCSTools stops its Newton iteration when
	// the step and residual are below 2.8e-8 deg (tnxpos.c), so sky2xy
	// itself is only good to ~1e-5 px here (measured max 1.4e-5 px vs the
	// pixel it came from).
	tnxTolWCSToolsPx = 1e-4
)

type tnxGolden struct {
	WCSTools     string            `json:"wcstools"`
	Header       string            `json:"header"`
	HeaderSHA256 string            `json:"header_sha256"`
	P2S          [][4]float64      `json:"p2s"`
	S2P          [][4]tnxJSONFloat `json:"s2p"`
	// NoAnswer is the xy2sky line of a refused header WCSTools does not
	// read as celestial (tnxRefusedNoAnswer); P2S and S2P are then empty.
	NoAnswer string `json:"wcstools_no_answer"`
}

// tnxJSONFloat is a golden number that may be null: sky2xy prints "nan"
// where its own inverse fails (only on refused headers), which
// gen_golden.py writes as null. It decodes to NaN, which fails every
// tolerance check of an accepted case.
type tnxJSONFloat float64

func (f *tnxJSONFloat) UnmarshalJSON(b []byte) error {
	if string(b) == "null" {
		*f = tnxJSONFloat(math.NaN())
		return nil
	}
	var v float64
	if err := json.Unmarshal(b, &v); err != nil {
		return err
	}
	*f = tnxJSONFloat(v)
	return nil
}

// readCardFile reads a one-card-per-line header file and returns the
// parsed header, exactly as a FITS reader would see it (80-column cards).
func readCardFile(t *testing.T, path string) (*header.Header, []byte) {
	t.Helper()
	text, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	var buf []byte
	for _, line := range strings.Split(strings.TrimRight(string(text), "\n"), "\n") {
		if len(line) > header.CardWidth {
			t.Fatalf("%s: card longer than 80 columns: %q", path, line)
		}
		buf = append(buf, line...)
		buf = append(buf, strings.Repeat(" ", header.CardWidth-len(line))...)
	}
	cards, _, err := header.ParseCards(buf)
	if err != nil {
		t.Fatalf("%s: ParseCards: %v", path, err)
	}
	return header.FromCards(cards), text
}

// tnxRequiredCases lists every TNX case the committed goldens must hold,
// with its point count (gen_golden.py: a 7x7 grid + 7 extra pixels = 56),
// so a deleted case or a thinned golden fails the default go test even
// when gen_golden.py cannot be run. The syn_* names are exactly the
// headers make_synthetic.py writes (TestTNXSyntheticHeadersReproduce).
var tnxRequiredCases = map[string]int{
	// Real headers.
	"kpno_kp4m_mosaic1_ccd1": 56, // KPNO 4 m Mosaic-1: polynomial, full cross, interleaved WAT
	"noao_ct4m_tnx1":         56, // CTIO 4 m Mosaic II: Chebyshev, half cross
	"registry_tnxdoc_tnx3":   56, // registry TNX doc: polynomial, half cross
	// Synthetic headers (make_synthetic.py).
	"syn_cheb_none_pc":              56, // Chebyshev, no cross terms, PC + CDELT
	"syn_leg_full":                  56, // Legendre, full cross terms
	"syn_leg_lonpole_card":          56, // LONPOLE card beats the TNX default 180
	"syn_leg_wat_1999":              56, // Legendre 9x9, joined WAT1 exactly 1999 characters (the WCSTools limit), longpole is its last token
	"syn_poly_cor_fallback":         56, // lngcor only in WAT2, latcor only in WAT1 (cross-axis fallback)
	"syn_poly_full_watlongpole":     56, // longpole only as a WAT attribute
	"syn_poly_ecliptic":             56, // ELON-TNX / ELAT-TNX
	"syn_poly_galactic":             56, // GLON-TNX / GLAT-TNX
	"syn_poly_helioecliptic":        56, // HLON-TNX / HLAT-TNX
	"syn_poly_supergalactic":        56, // SLON-TNX / SLAT-TNX
	"syn_poly_half_pole":            56, // CRVAL2 = +90, TNX default LONPOLE 180
	"syn_poly_lonpole_card_and_wat": 56, // LONPOLE card beats WAT longpole
	"syn_poly_wat_precedence":       56, // WAT1 longpole/lngcor and WAT2 latcor win
	"syn_tan_no_wat":                56, // RA---TNX with no WAT cards: plain TAN
	"syn_poly_pv2_ignored":          56, // PV2_1/PV2_2: ignored by WCSTools and by wcslib TAN
	"syn_poly_cunit_deg":            56, // CUNIT1/2 = 'deg' accepted
	"syn_poly_lonpole_360":          56, // LONPOLE = 360, the largest value WCSTools reads
	"syn_poly_lonpole_neg400":       56, // LONPOLE = -400 used as given
	"syn_poly_wat_nondigit_keys":    56, // WAT1_A / WAT2_1X are not WAT cards and are ignored
	"syn_poly_projp1":               56, // PROJP1 (wcs.Parse: PV2_1, latitude axis): ignored by both
	"syn_tan_no_wat_projp1":         56, // PROJP1 with no WAT cards: plain TAN in both
	"syn_poly_pc001001":             56, // CDELT + PC00i00j (WCSTools reads PC001001 when PC1_1 is absent)
	"syn_poly_pc_no_pc2_2":          56, // CDELT + PC1_1, PC1_2, PC2_1 (PC2_2 defaults to 1 in both)
	"syn_poly_cd_offdiag_only":      56, // CD1_2 / CD2_1 only (CD1_1, CD2_2 default 0 in both)
	"syn_poly_cd_and_cdelt":         56, // CD beside CDELT: both use CD and ignore CDELT
	"syn_poly_no_scale":             56, // no CD/CDELT/PC: 1 deg/pixel in both
	"syn_poly_cd_and_cd001001":      56, // CD001001.. beside equal CD1_1..: CD in both
	"syn_poly_cd_and_cdelt1":        56, // CD beside CDELT1 alone: CD in both (the CDELT pair rule needs no CD)
	"syn_poly_lngcor_only_pv2":      56, // lngcor only + PV2_m: tnxinit succeeds, PV2_m ignored by both
}

// tnxRefusedCases lists the headers of testdata/tnx/refused: headers that
// WCSTools (the TNX truth) reads one way and the wcslib conventions of
// wcs.Parse / New another, or where WCSTools ignores a card. NewTNX must
// refuse each with an error naming the cause; the value is a substring
// of that error. Each refused/<name>.json is the WCSTools answer, so a
// regression that accepts the header is reported with its disagreement.
var tnxRefusedCases = map[string]string{
	"pv1_fiducial":      "PV1_1",
	"pv1_3_lonpole":     "PV1_3",
	"cd_and_pc":         "both CD and PC",
	"lonpole_361":       "LONPOLE = 361",
	"lonpole_400":       "LONPOLE = 400",
	"cunit_arcmin":      "CUNIT1 = \"arcmin\"",
	"wat_duplicate_key": "WAT1_002 appears 2 times",
	"cdelt1_only":       "CDELT1 without CDELT2",
	"cdelt2_only":       "CDELT2 without CDELT1",
	"pc_no_cdelt":       "PC matrix without CDELT1",
	"pc_no_pc1_1":       "PC matrix without PC1_1",
	"cd001001_only":     "CD001001 without CD1_1",
	"nowat_pv2":         "no TNX correction surface and PV2_1",
	"watnocor_pv2":      "no TNX correction surface and PV2_1",
	"tnx_sip_suffix":    "TNX header carries SIP keywords",
	"tnx_sip_nosuffix":  "TNX header carries SIP keywords",
	"ctype_lowercase":   `CTYPE1 = "RA---tnx" is not upper case`,
}

// tnxRefusedNoAnswer are the refused cases WCSTools does not read as
// celestial at all: xy2sky prints asterisks for every point, so the
// golden records that line (wcstools_no_answer) and no points.
var tnxRefusedNoAnswer = map[string]bool{
	"ctype_lowercase": true,
}

func tnxCaseNames(t *testing.T) []string {
	t.Helper()
	hdrs, _ := filepath.Glob(filepath.Join(tnxHeaderDir, "*.hdr"))
	golds, _ := filepath.Glob(filepath.Join(tnxGoldenDir, "*.json"))
	have := map[string]bool{}
	for _, g := range golds {
		have[strings.TrimSuffix(filepath.Base(g), ".json")] = true
	}
	onDisk := map[string]bool{}
	for _, h := range hdrs {
		n := strings.TrimSuffix(filepath.Base(h), ".hdr")
		onDisk[n] = true
		if !have[n] {
			t.Errorf("header %s has no golden (run gen_golden.py, see testdata/tnx/PROVENANCE.md)", n)
		}
		if _, ok := tnxRequiredCases[n]; !ok {
			t.Errorf("header %s is not in tnxRequiredCases; add it there", n)
		}
		delete(have, n)
	}
	for n := range have {
		t.Errorf("orphan golden %s.json with no header", n)
	}
	var names []string
	for n := range tnxRequiredCases {
		if !onDisk[n] {
			t.Errorf("required TNX case %s: headers/%s.hdr is missing", n, n)
			continue
		}
		names = append(names, n)
	}
	sort.Strings(names)
	return names
}

// TestTNXRefusedHeaders: every refused/<name>.hdr (make_synthetic.py)
// must be refused by NewTNX with the error of tnxRefusedCases, and every
// listed case must have its header and WCSTools golden. If NewTNX accepts
// one, the test reports how far its answer is from WCSTools.
func TestTNXRefusedHeaders(t *testing.T) {
	const dir = "testdata/tnx/refused"
	hdrs, _ := filepath.Glob(filepath.Join(dir, "*.hdr"))
	onDisk := map[string]bool{}
	for _, p := range hdrs {
		n := strings.TrimSuffix(filepath.Base(p), ".hdr")
		onDisk[n] = true
		if _, ok := tnxRefusedCases[n]; !ok {
			t.Errorf("refused/%s.hdr is not in tnxRefusedCases", n)
		}
	}
	golds, _ := filepath.Glob(filepath.Join(dir, "*.json"))
	for _, g := range golds {
		if n := strings.TrimSuffix(filepath.Base(g), ".json"); !onDisk[n] {
			t.Errorf("orphan refused golden %s.json", n)
		}
	}
	names := make([]string, 0, len(tnxRefusedCases))
	for n := range tnxRefusedCases {
		names = append(names, n)
	}
	sort.Strings(names)
	for _, name := range names {
		t.Run(name, func(t *testing.T) {
			if !onDisk[name] {
				t.Fatalf("refused/%s.hdr is missing", name)
			}
			h, text := readCardFile(t, filepath.Join(dir, name+".hdr"))
			raw, err := os.ReadFile(filepath.Join(dir, name+".json"))
			if err != nil {
				t.Fatalf("golden: %v", err)
			}
			var g tnxGolden
			if err := json.Unmarshal(raw, &g); err != nil {
				t.Fatalf("decode golden: %v", err)
			}
			sum := sha256.Sum256(text)
			if g.HeaderSHA256 != hex.EncodeToString(sum[:]) {
				t.Fatalf("golden was made from a different %s.hdr; regenerate", name)
			}
			if tnxRefusedNoAnswer[name] {
				if !strings.HasPrefix(g.NoAnswer, "*") || len(g.P2S) != 0 || len(g.S2P) != 0 {
					t.Fatalf("golden should record WCSTools' asterisks and no points, has %q and %d/%d points", g.NoAnswer, len(g.P2S), len(g.S2P))
				}
			} else if len(g.P2S) != 56 || g.NoAnswer != "" {
				t.Fatalf("golden has %d p2s points (want 56) and no-answer %q", len(g.P2S), g.NoAnswer)
			}
			tr, err := NewTNX(h)
			if err == nil && g.NoAnswer != "" {
				a, d, _ := tr.PixelToSky(0, 0)
				t.Fatalf("NewTNX accepted it (pixel (1, 1) -> %.10f, %.10f), but WCSTools gives no answer (xy2sky: %q)", a, d, g.NoAnswer)
			}
			if err == nil {
				worst, at := 0.0, g.P2S[0]
				for _, p := range g.P2S {
					a, d, err := tr.PixelToSky(p[0]-1, p[1]-1)
					if err != nil {
						continue
					}
					if s := angSepArcsec(a, d, p[2], p[3]); s > worst {
						worst, at = s, p
					}
				}
				t.Fatalf("NewTNX accepted it; worst disagreement with WCSTools %.6g arcsec at pixel (%v, %v)", worst, at[0], at[1])
			}
			if want := tnxRefusedCases[name]; !strings.Contains(err.Error(), want) {
				t.Fatalf("error %q does not mention %q", err, want)
			}
		})
	}
}

func TestTNXWCSToolsGolden(t *testing.T) {
	for _, name := range tnxCaseNames(t) {
		t.Run(name, func(t *testing.T) {
			h, text := readCardFile(t, filepath.Join(tnxHeaderDir, name+".hdr"))
			raw, err := os.ReadFile(filepath.Join(tnxGoldenDir, name+".json"))
			if err != nil {
				t.Fatalf("golden: %v", err)
			}
			var g tnxGolden
			if err := json.Unmarshal(raw, &g); err != nil {
				t.Fatalf("decode golden: %v", err)
			}
			sum := sha256.Sum256(text)
			if g.HeaderSHA256 != hex.EncodeToString(sum[:]) {
				t.Fatalf("golden was made from a different %s.hdr (sha256 %s, file %x); regenerate", name, g.HeaderSHA256, sum)
			}
			if !strings.Contains(g.WCSTools, "WCSTools 3.9.7,") {
				t.Fatalf("golden from %q, PROVENANCE.md pins WCSTools 3.9.7", g.WCSTools)
			}
			if want := tnxRequiredCases[name]; len(g.P2S) != want || len(g.S2P) != want {
				t.Fatalf("golden has %d p2s / %d s2p points, tnxRequiredCases wants %d of each", len(g.P2S), len(g.S2P), want)
			}
			tr, err := NewTNX(h)
			if err != nil {
				t.Fatalf("NewTNX: %v", err)
			}
			var worstAs, worstPx, worstW float64
			bad := 0
			for i, p := range g.P2S {
				a, d, err := tr.PixelToSky(p[0]-1, p[1]-1)
				if err != nil {
					t.Errorf("pixel (%v, %v): %v", p[0], p[1], err)
					continue
				}
				sep := angSepArcsec(a, d, p[2], p[3])
				worstAs = math.Max(worstAs, sep)
				if !(sep <= tnxTolArcsec) {
					if bad++; bad <= 3 {
						t.Errorf("p2s pixel (%v, %v): got (%.12f, %.12f), WCSTools (%.12f, %.12f): %.3g arcsec",
							p[0], p[1], a, d, p[2], p[3], sep)
					}
				}
				gs := g.S2P[i]
				s := [4]float64{float64(gs[0]), float64(gs[1]), float64(gs[2]), float64(gs[3])}
				x, y, err := tr.SkyToPixel(s[0], s[1])
				if err != nil {
					t.Errorf("sky (%v, %v): %v", s[0], s[1], err)
					continue
				}
				x, y = x+1, y+1
				e := math.Max(math.Abs(x-p[0]), math.Abs(y-p[1]))
				w := math.Max(math.Abs(x-s[2]), math.Abs(y-s[3]))
				worstPx, worstW = math.Max(worstPx, e), math.Max(worstW, w)
				if !(e <= tnxTolPx) || !(w <= tnxTolWCSToolsPx) {
					if bad++; bad <= 3 {
						t.Errorf("s2p sky (%.12f, %.12f): got pixel (%.10f, %.10f), source pixel (%v, %v), sky2xy (%.10f, %.10f)",
							s[0], s[1], x, y, p[0], p[1], s[2], s[3])
					}
				}
			}
			if bad > 3 {
				t.Errorf("%d bad points of %d", bad, 2*len(g.P2S))
			}
			t.Logf("%d points: p2s worst %.3g arcsec; s2p worst %.3g px vs source pixel, %.3g px vs sky2xy",
				len(g.P2S), worstAs, worstPx, worstW)
		})
	}
}

// TestNewRejectsTNX: New sees only *wcs.Header, which has no WAT strings,
// so it must refuse a TNX header rather than return an uncorrected TAN.
func TestNewRejectsTNX(t *testing.T) {
	h, _ := readCardFile(t, filepath.Join(tnxHeaderDir, "noao_ct4m_tnx1.hdr"))
	w, err := wcs.Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if tr, err := New(w); err == nil {
		t.Fatalf("New accepted a TNX header (Code %q); it must reject it and leave TNX to NewTNX", tr.Code())
	}
}
