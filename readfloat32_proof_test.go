package fits

import (
	"bytes"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Proof tests for ReadFloat32 / ReadFrame / ReadFrameRGB against
// independent truth (testdata/imaging/PROVENANCE.md):
//   - readfloat32_expected.json: astropy-read stored values mapped by the
//     documented contract with exact rational arithmetic (gen_expected.py);
//   - sirilref/golden: Siril's own FITS->float code (sirilref.c, verbatim
//     Siril functions linked against cfitsio).

type rf32Expected struct {
	Bitpix int      `json:"bitpix"`
	Rule   string   `json:"rule"`
	Shape  []int64  `json:"shape"`
	Npix   int      `json:"npix"`
	Blank  *int64   `json:"blank"`
	Bits   pixArray `json:"bits"` // float32 bits, null = NaN; a digest past 1024 values
	// InvalidCard names the BSCALE or BZERO card whose value is not a
	// finite FITS number (rule "error" from gen_expected.py bad_scaling).
	InvalidCard string `json:"invalid_card"`
}

type sirilGolden struct {
	File        string   `json:"file"`
	Bitpix      int      `json:"bitpix"`
	Normalized  bool     `json:"normalized"`
	Status      int      `json:"status"`
	Unsupported bool     `json:"unsupported"`
	Bits        pixArray `json:"bits"` // float32 bits; a digest past 1024 values
}

func imagingPath(name string) string { return filepath.Join("testdata", "imaging", name) }

func mustReadFile(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("committed test file missing (this is a failure, not a skip): %v", err)
	}
	return b
}

func loadRF32Expected(t *testing.T) map[string]rf32Expected {
	t.Helper()
	var g struct {
		Files map[string]rf32Expected `json:"files"`
	}
	if err := json.Unmarshal(mustReadFile(t, imagingPath("readfloat32_expected.json")), &g); err != nil {
		t.Fatal(err)
	}
	if len(g.Files) < 40 {
		t.Fatalf("readfloat32_expected.json lists %d files, want >= 40", len(g.Files))
	}
	return g.Files
}

func loadSiril(t *testing.T, name string) sirilGolden {
	t.Helper()
	var g sirilGolden
	if err := json.Unmarshal(mustReadFile(t, imagingPath(filepath.Join("sirilref", "golden", name+".json"))), &g); err != nil {
		t.Fatal(err)
	}
	return g
}

// readFloat32File opens path and runs ReadFloat32 on the frame image
// (primary, or HDU 1 behind a NAXIS=0 primary).
func readFloat32File(t *testing.T, path string) ([]float32, error) {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	var img HDU
	prim, err := f.Primary()
	if err != nil {
		t.Fatalf("%s: %v", path, err)
	}
	img = prim
	if prim.NAXIS() == 0 {
		if img, err = f.HDU(1); err != nil {
			t.Fatalf("%s: HDU 1: %v", path, err)
		}
	}
	return ReadFloat32(img)
}

func f32bitsEq(got float32, want *int64) bool {
	if want == nil {
		return math.IsNaN(float64(got))
	}
	return int64(math.Float32bits(got)) == *want
}

// TestReadFloat32_ExpectedMatrix: every BITPIX x BSCALE/BZERO x BLANK x
// NaN sample (plain and fpack-compressed) against the exact expected bits.
func TestReadFloat32_ExpectedMatrix(t *testing.T) {
	exp := loadRF32Expected(t)
	names := make([]string, 0, len(exp))
	for n := range exp {
		names = append(names, n)
	}
	sort.Strings(names)
	rules := map[string]int{}
	for _, name := range names {
		e := exp[name]
		rules[e.Rule]++
		t.Run(name, func(t *testing.T) {
			got, err := readFloat32File(t, imagingPath(name))
			if e.Rule == "error" {
				if err == nil {
					t.Fatalf("ReadFloat32 returned %v, want an error (BSCALE=0, a physical value that overflows float32, or an integer image with a ZBLANK column)", got)
				}
				return
			}
			if err != nil {
				t.Fatalf("ReadFloat32: %v", err)
			}
			if len(got) != e.Npix {
				t.Fatalf("got %d pixels, want %d", len(got), e.Npix)
			}
			if !e.Bits.Present() {
				return // all-value files: checked against Siril below
			}
			e.Bits.Check(t, "ReadFloat32 (rule "+e.Rule+")", f32Canon(got))
		})
	}
	for _, r := range []string{"siril8", "siril16", "exact", "float", "error"} {
		if rules[r] == 0 {
			t.Errorf("no sample exercises rule %q", r)
		}
	}
}

func formatF32(f float32) string { return strconv.FormatFloat(float64(f), 'g', -1, 32) }

// TestReadFloat32_SirilConvention compares ReadFloat32 with Siril's own
// FITS->float loader (readfits with force_float) on every sample, and pins
// each documented difference.
func TestReadFloat32_SirilConvention(t *testing.T) {
	exp := loadRF32Expected(t)
	// relation of ReadFloat32 to Siril, per sample:
	//   equal      bit-identical at every pixel
	//   blank      bit-identical except at BLANK pixels, where ReadFloat32
	//              gives NaN and Siril normalizes the sentinel value
	//   differs    documented difference (reason given)
	//   unsupported Siril refuses the file (BITPIX=64)
	//   error      ReadFloat32 refuses the file (BSCALE=0)
	type rel struct{ kind, why string }
	relation := map[string]rel{
		"u8_all.fits":                    {"equal", ""},
		"u16_all.fits":                   {"equal", ""},
		"i16_all.fits":                   {"equal", ""},
		"asi533_u16_crop.fits":           {"equal", ""},
		"comp_asi533_u16_crop.fits.fz":   {"equal", ""},
		"u16_orient.fits":                {"equal", ""},
		"u16_rgb.fits":                   {"equal", ""},
		"u8_blank.fits":                  {"blank", ""},
		"comp_u8_blank.fits.fz":          {"blank", ""},
		"u16_blank.fits":                 {"blank", ""},
		"comp_u16_blank.fits.fz":         {"blank", ""},
		"f32_unit.fits":                  {"equal", ""},
		"f32_bscale.fits":                {"equal", ""},
		"i16_nobzero.fits":               {"differs", "no BZERO card: Siril treats BITPIX=16 as unsigned without the 32768 offset (manage_bitpix); ReadFloat32 treats it as signed"},
		"i16_bscale0.fits":               {"error", ""},
		"i16_bscale_neg.fits":            {"differs", "BSCALE=-1: Siril multiplies cfitsio's physical value by 1/65535 (negative results); ReadFloat32 maps the physical range onto [0,1]"},
		"i16_bscale_frac.fits":           {"differs", "BSCALE=0.5: Siril scales the physical value by 1/65535; ReadFloat32 maps the physical range"},
		"i8_signed.fits":                 {"differs", "BZERO=-128: Siril multiplies the physical value (-128..127) by 1/255; ReadFloat32 maps the physical range"},
		"u8_bscale.fits":                 {"differs", "BSCALE=2 BZERO=10: Siril scales the physical value by 1/255"},
		"i32.fits":                       {"differs", "BITPIX=32: Siril returns the raw values (normalizes only with DATAMAX > 10); ReadFloat32 maps the full range"},
		"u32.fits":                       {"differs", "BITPIX=32: as i32.fits"},
		"i32_blank_bscale.fits":          {"differs", "BITPIX=32: as i32.fits"},
		"comp_i32_blank_bscale.fits.fz":  {"differs", "BITPIX=32: as i32.fits"},
		"i32_bscale_neg.fits":            {"differs", "BITPIX=32: as i32.fits"},
		"i32_round.fits":                 {"differs", "BITPIX=32: as i32.fits"},
		"i32_round_neg.fits":             {"differs", "BITPIX=32: as i32.fits"},
		"f32_wide.fits":                  {"differs", "float data with maximum > 10: Siril divides by 65535; ReadFloat32 passes floats through"},
		"f64.fits":                       {"equal", ""},
		"f64_wide.fits":                  {"siril-error", "cfitsio reports NUM_OVERFLOW (412) converting the +Inf double to float, so Siril cannot load the file; ReadFloat32 keeps +Inf"},
		"f64_overflow.fits":              {"error", ""},
		"f32_overflow.fits":              {"error", ""},
		"u8_bscale_2d.fits":              {"differs", "BSCALE=2 BZERO=10: as u8_bscale.fits"},
		"u16_blank_2d.fits":              {"blank", ""},
		"u16_noblankhit_2d.fits":         {"equal", ""},
		"u16_bscale2_2d.fits":            {"differs", "BSCALE=2: Siril scales the physical value by 1/65535; ReadFloat32 maps the physical range"},
		"f32_2d.fits":                    {"differs", "float data with maximum > 10 (1e30): Siril divides by 65535"},
		"f32_bscale_2d.fits":             {"equal", ""},
		"f64_2d.fits":                    {"equal", ""},
		"f64_exact_2d.fits":              {"differs", "float data with maximum > 10 (65535): Siril divides by 65535"},
		"i64.fits":                       {"unsupported", ""},
		"u64.fits":                       {"unsupported", ""},
		"u64_intcard.fits":               {"unsupported", ""},
		"i64_blank.fits":                 {"unsupported", ""},
		"i64_round.fits":                 {"unsupported", ""},
		"i64_round_neg.fits":             {"unsupported", ""},
		"comp_u8_zbscale.fits.fz":        {"blank", ""},
		"comp_i8_signed.fits.fz":         {"differs", "BZERO=-128: as i8_signed.fits"},
		"comp_i16_bscale0.fits.fz":       {"error", ""},
		"comp_i16_bscale_neg.fits.fz":    {"differs", "BSCALE=-1: as i16_bscale_neg.fits"},
		"comp_i32_zblank_bscale.fits.fz": {"differs", "BITPIX=32: as i32.fits"},
		"comp_u16_zblank.fits.fz":        {"blank", ""},
		"comp_u16_zblank_blank.fits.fz":  {"blank", ""},
		"comp_u16_zblank_col.fits.fz":    {"blank", ""}, // per-tile ZBLANK column: cfitsio null flags
		"comp_u16_rgb.fits.fz":           {"equal", ""},
		"f32_bscale0_2d.fits":            {"error", ""},
		"f32_bscale_inexact_2d.fits":     {"equal", ""},
		"f32_bzero_2d.fits":              {"differs", "float data with maximum > 10 (physical 110): Siril divides by 65535"},
		"f64_bzero_2d.fits":              {"differs", "float data with maximum > 10 (physical 997): Siril divides by 65535"},
		"f32_bzero_inexact_2d.fits":      {"differs", "float data with maximum > 10 (physical 12): Siril divides by 65535"},
		"i16_bscale_str_2d.fits":         {"error", "BSCALE = 'abc': cfitsio (so Siril) silently uses BSCALE 1; ReadFloat32 refuses the invalid card"},
		"i16_bscale_bad_2d.fits":         {"error", "BSCALE = 2.0.0: as i16_bscale_str_2d.fits"},
		"f32_bzero_huge_2d.fits":         {"error", "BZERO = 1E400 (not a finite double): cfitsio uses BZERO 0; ReadFloat32 refuses it"},
		"frame_naxis1_primary.fits":      {"equal", ""},
		"f32_bzero_qnum_2d.fits":         {"error", "BZERO = '10.0' (a string): cfitsio (so Siril) silently uses BZERO 0; ReadFloat32 refuses it"},
		"i16_bscale_qnum_2d.fits":        {"error", "BSCALE = '0.5' (a string): cfitsio uses BSCALE 1; ReadFloat32 refuses it"},
		"f32_bzero_hex_2d.fits":          {"error", "BZERO = 0x1p3 (not a FITS number): cfitsio's strtod reads 8; ReadFloat32 refuses it"},
		"i16_dexp_2d.fits":               {"differs", "BSCALE=0.5 BZERO=-3 (D exponents): as i16_bscale_frac.fits"},
		"f32_bzero_dexp_2d.fits":         {"differs", "float data with maximum > 10 (physical 110): as f32_bzero_2d.fits"},
		// card_* (BSCALE / BZERO / BLANK card variants): Siril reads them
		// through cfitsio, so where ReadFloat32 refuses the card Siril
		// gets cfitsio's silent reading; the accepted float ones have a
		// physical maximum > 10, which Siril divides by 65535.
		"card_f32_bzero_dup_2d.fits":            {"differs", "float maximum > 10 (physical 120): as f32_bzero_2d.fits"},
		"card_f32_bzero_dupbadfirst_2d.fits":    {"differs", "float maximum > 10 (physical 120): as f32_bzero_2d.fits"},
		"card_f32_bzero_int_2d.fits":            {"differs", "float maximum > 10 (physical 110): as f32_bzero_2d.fits"},
		"card_f32_bzero_plus_2d.fits":           {"differs", "float maximum > 10 (physical 110): as f32_bzero_2d.fits"},
		"card_f32_bzero_dotlead_2d.fits":        {"differs", "float maximum > 10 (physical 105): as f32_bzero_2d.fits"},
		"card_f32_bzero_dottrail_2d.fits":       {"differs", "float maximum > 10 (physical 110): as f32_bzero_2d.fits"},
		"card_f32_bzero_lowere_2d.fits":         {"differs", "float maximum > 10 (physical 110): as f32_bzero_2d.fits"},
		"card_i16_bscale_dup_2d.fits":           {"differs", "BSCALE -1 (the last card): as i16_bscale_neg.fits"},
		"comp_card_i16_bscale_dup_2d.fits.fz":   {"differs", "tile-compressed, BSCALE 2.0 then -1.0: cfitsio, astropy and ReadFloat32 apply the first (2.0); Siril multiplies by 1/65535"},
		"comp_card_i16_bscale_str.fits.fz":      {"error", "tile-compressed BSCALE = 'abc': cfitsio (so Siril) uses BSCALE 1, astropy fails"},
		"comp_card_i16_bscale_hierarch.fits.fz": {"error", "tile-compressed HIERARCH BSCALE = -1.0: cfitsio (so Siril) and astropy apply it on this header, astropy ignores it on a plain image; the scaling-card rule refuses it everywhere"},
		"card_i16_blank_dup_2d.fits":            {"differs", "BSCALE=2: Siril multiplies the physical value by 1/65535; ReadFloat32 maps the range onto [0,1]"},
		"card_f32_bzero_dupbadlast_2d.fits":     {"error", "last BZERO = 'abc': cfitsio (so Siril) uses BZERO 0; astropy fails"},
		"card_f32_bzero_noeqsp_2d.fits":         {"error", "BZERO   =10.0: cfitsio (so Siril) reads 10, astropy 0"},
		"card_f32_bzero_col10eq_2d.fits":        {"error", "'=' in column 10: cfitsio reads 10, astropy 0"},
		"card_f32_bzero_col11eq_2d.fits":        {"error", "'=' in column 11: cfitsio reads 10, astropy 0"},
		"card_f32_bzero_spacein_2d.fits":        {"error", "BZERO = 1 0.0: cfitsio reads strtod's prefix 1, astropy fails"},
		"card_f32_bzero_trail_2d.fits":          {"error", "BZERO = 10.0 x: cfitsio reads 10, astropy fails"},
		"card_f32_bzero_lowerd_2d.fits":         {"error", "BZERO = 1.0d1: cfitsio drops it (0), astropy reads 10"},
		"card_f32_bzero_lowkey_2d.fits":         {"error", "bzero (lower case): cfitsio drops it (0), astropy reads 10"},
		"card_f32_bzero_empty_2d.fits":          {"error", "BZERO with no value: cfitsio uses 0, astropy fails"},
		"card_i16_blank_float_2d.fits":          {"error", "BLANK = 100.0: cfitsio and astropy ignore the null marker"},
		"card_i16_blank_str_2d.fits":            {"error", "BLANK = 'abc': cfitsio and astropy ignore the null marker"},
		"card_i16_blank_noeq_2d.fits":           {"error", "BLANK   =100: cfitsio applies it, astropy does not"},
		"card_f32_bzero_eqcol6_2d.fits":         {"error", "BZERO=10.0: cfitsio reads 10, astropy 0"},
		"card_f32_bzero_afterend_2d.fits":       {"differs", "float maximum > 10 (physical 110; the BZERO card after END is not read): as f32_bzero_2d.fits"},
		"card_f32_blank_str_2d.fits":            {"differs", "float maximum > 10 (100; BLANK on float data ignored): as f32_bzero_2d.fits"},
	}
	if len(relation) != len(exp) {
		t.Fatalf("relation table covers %d samples, expected json has %d", len(relation), len(exp))
	}
	for name, r := range relation {
		t.Run(name, func(t *testing.T) {
			e, ok := exp[name]
			if !ok {
				t.Fatalf("no expected entry")
			}
			s := loadSiril(t, name)
			if r.kind == "unsupported" {
				if !s.Unsupported {
					t.Fatalf("Siril read the file; relation table says unsupported")
				}
				return
			}
			got, err := readFloat32File(t, imagingPath(name))
			if r.kind == "error" {
				if err == nil {
					t.Fatalf("want error")
				}
				return
			}
			if r.kind == "siril-error" {
				if s.Status == 0 {
					t.Fatalf("Siril harness read the file; relation table says it fails (%s)", r.why)
				}
				if err != nil {
					t.Fatalf("ReadFloat32: %v", err)
				}
				return
			}
			if s.Unsupported || s.Status != 0 {
				t.Fatalf("Siril harness: unsupported=%v status=%d", s.Unsupported, s.Status)
			}
			if err != nil {
				t.Fatal(err)
			}
			if s.Bits.Digested() {
				// Large samples (all 65536 16-bit values, the camera
				// crop): the golden is an exact digest, so only an
				// identical result can be checked; each is "equal".
				if r.kind != "equal" {
					t.Fatalf("Siril golden is a digest; only relation \"equal\" can be checked against it, table says %q", r.kind)
				}
				s.Bits.Check(t, "ReadFloat32 vs Siril", f32Bits(got, false))
				return
			}
			sb := s.Bits.Explicit(t)
			if len(sb) != len(got) {
				t.Fatalf("Siril %d pixels, ReadFloat32 %d", len(sb), len(got))
			}
			var eb []*int64
			if e.Bits.Present() {
				eb = e.Bits.Explicit(t)
			}
			diff, blankDiff := 0, 0
			for i := range got {
				isBlank := eb != nil && eb[i] == nil // a null pixel (BLANK, ZBLANK or a ZBLANK column)
				sbits := uint32(*sb[i])
				same := math.Float32bits(got[i]) == sbits ||
					(math.IsNaN(float64(got[i])) && math.IsNaN(float64(math.Float32frombits(sbits))))
				if same {
					continue
				}
				if isBlank && math.IsNaN(float64(got[i])) {
					blankDiff++
					continue
				}
				if diff++; diff <= 4 && r.kind != "differs" {
					t.Errorf("pixel %d: ReadFloat32 0x%08x, Siril 0x%08x", i, math.Float32bits(got[i]), sbits)
				}
			}
			switch r.kind {
			case "equal":
				if diff != 0 || blankDiff != 0 {
					t.Errorf("%d pixels differ from Siril (%d BLANK)", diff, blankDiff)
				}
			case "blank":
				if diff != 0 || blankDiff == 0 {
					t.Errorf("want differences only at BLANK pixels: %d other, %d BLANK", diff, blankDiff)
				}
			case "differs":
				if diff == 0 {
					t.Errorf("documented difference (%s) not observed: ReadFloat32 equals Siril; update the doc", r.why)
				}
			}
		})
	}
}

// TestReadFrame_RealCamera reads the real ZWO ASI533MC Pro crop (BITPIX=16
// BZERO=32768) through ReadFrame and checks size, row-major layout and
// every pixel against Siril's loader.
func TestReadFrame_RealCamera(t *testing.T) {
	for _, name := range []string{"asi533_u16_crop.fits", "comp_asi533_u16_crop.fits.fz"} {
		t.Run(name, func(t *testing.T) {
			pix, w, h, err := ReadFrame(imagingPath(name))
			if err != nil {
				t.Fatal(err)
			}
			if w != 128 || h != 128 || len(pix) != w*h {
				t.Fatalf("got %dx%d (%d pixels), want 128x128", w, h, len(pix))
			}
			// pix[y*w+x] is FITS pixel (x, y): the same order as Siril's.
			loadSiril(t, name).Bits.Check(t, "ReadFrame vs Siril", f32Bits(pix, false))
		})
	}
}

// TestReadFrame_RowOrientation: u16_orient.fits stores physical value
// 1000*(4y+x) at 0-based FITS pixel (x, y) (x along NAXIS1). ReadFrame
// must return it at pix[y*width+x]: no flip, whatever ROWORDER says.
func TestReadFrame_RowOrientation(t *testing.T) {
	exp := loadRF32Expected(t)["u16_orient.fits"]
	pix, w, h, err := ReadFrame(imagingPath("u16_orient.fits"))
	if err != nil {
		t.Fatal(err)
	}
	if w != 4 || h != 3 {
		t.Fatalf("dims %dx%d, want 4x3", w, h)
	}
	raw, err := readPhysical(t, imagingPath("u16_orient.fits"))
	if err != nil {
		t.Fatal(err)
	}
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			if raw[y*w+x] != float64(1000*(4*y+x)) {
				t.Fatalf("fixture: pixel (%d,%d) holds %v, want %d", x, y, raw[y*w+x], 1000*(4*y+x))
			}
			if !f32bitsEq(pix[y*w+x], exp.Bits.Explicit(t)[y*w+x]) {
				t.Errorf("pixel (%d,%d): got %v", x, y, pix[y*w+x])
			}
		}
	}
}

// TestReadFrame_PrimaryNAXIS1UsesHDU1: the documented HDU rule is "the
// primary if it has NAXIS >= 2, otherwise HDU 1". frame_naxis1_primary.fits
// has a 1-D primary (5 pixels) in front of an IMAGE extension holding the
// u16_orient.fits pixels, so ReadFrame must return those (4x3, expected
// bits from readfloat32_expected.json for u16_orient.fits), not an error
// about the primary's NAXIS=1. astropy confirms the extension's pixels
// equal u16_orient's.
func TestReadFrame_PrimaryNAXIS1UsesHDU1(t *testing.T) {
	exp := loadRF32Expected(t)["u16_orient.fits"]
	pix, w, h, err := ReadFrame(imagingPath("frame_naxis1_primary.fits"))
	if err != nil {
		t.Fatalf("ReadFrame: %v (want HDU 1's 4x3 image)", err)
	}
	want := exp.Bits.Explicit(t)
	if w != 4 || h != 3 || len(pix) != 12 || len(want) != 12 {
		t.Fatalf("got %dx%d (%d pixels), want 4x3 (expected %d)", w, h, len(pix), len(want))
	}
	exp.Bits.Check(t, "ReadFrame HDU 1", f32Canon(pix))
	// The primary itself is read by ReadFloat32 (TestReadFloat32_ExpectedMatrix).
	if e := loadRF32Expected(t)["frame_naxis1_primary.fits"]; e.Npix != 5 {
		t.Fatalf("expected entry for the 1-D primary has %d pixels, want 5", e.Npix)
	}
}

// TestReadFloat32_InvalidScalingCards: a BSCALE or BZERO card that is
// present but not a finite number (a string, '2.0.0', 1E400) is an error
// from ReadFloat32, ReadFrame, AsImage and image.DecodeConfig, on the
// plain path as on the tile-compressed one: returning the data unscaled
// (what cfitsio 4.7.0 fits_read_img does, status 0) would be wrong data
// without an error. Truth: the documented contract, with the samples
// classified from their raw cards by gen_expected.py (rule "error",
// kind "reject").
//
// Each sample is read twice: on a fresh Open (the stored cards are read
// without the header parser) and after Header() has parsed the header;
// both must refuse it, naming the card. The quoted numeric strings
// ('10.0', '0.5') are character-string values, not numbers (FITS 4.0
// 4.2.1); 0x1p3 is a C99 hex float that strtod (so cfitsio) reads as 8
// but is not a FITS integer or real (4.2.3, 4.2.4). cfitsio 4.7.0
// fits_read_img returns status 0 with the default (or strtod) scaling for
// all of them (PROVENANCE.md).
func TestReadFloat32_InvalidScalingCards(t *testing.T) {
	exp := loadRF32Expected(t)
	for _, name := range []string{"i16_bscale_str_2d.fits", "i16_bscale_bad_2d.fits", "f32_bzero_huge_2d.fits",
		"f32_bzero_qnum_2d.fits", "i16_bscale_qnum_2d.fits", "f32_bzero_hex_2d.fits"} {
		t.Run(name, func(t *testing.T) {
			if exp[name].Rule != "error" {
				t.Fatalf("expected rule %q, want error", exp[name].Rule)
			}
			card := exp[name].InvalidCard
			if card == "" {
				t.Fatalf("expected entry names no invalid card")
			}
			for _, parsed := range []bool{false, true} {
				f, err := Open(imagingPath(name))
				if err != nil {
					t.Fatal(err)
				}
				p, err := f.Primary()
				if err != nil {
					t.Fatal(err)
				}
				if parsed {
					p.Header()
				}
				got, err := ReadFloat32(p)
				if err == nil || !strings.Contains(err.Error(), card) || strings.Count(err.Error(), "fits.ReadFloat32") != 1 {
					t.Errorf("parsed=%v: ReadFloat32 = %v, %v; want one fits.ReadFloat32 error naming %s", parsed, got, err, card)
				}
				img, err := p.AsImage()
				if err == nil || !strings.Contains(err.Error(), card) || strings.Contains(err.Error(), "ReadFloat32") {
					t.Errorf("parsed=%v: AsImage = %T, %v; want an AsImage error naming %s", parsed, img, err, card)
				}
				f.Close()
			}
			if got, _, _, err := ReadFrame(imagingPath(name)); err == nil {
				t.Errorf("ReadFrame = %v, nil; want an error", got)
			}
			f, err := Open(imagingPath(name))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			p, err := f.Primary()
			if err != nil {
				t.Fatal(err)
			}
			if img, err := p.AsImage(); err == nil {
				t.Errorf("AsImage = %T, nil; want an error", img)
			}
			if cfg, err := DecodeConfig(bytes.NewReader(mustReadFile(t, imagingPath(name)))); err == nil {
				t.Errorf("DecodeConfig = %+v, nil; want an error", cfg)
			}
		})
	}
}

// TestReadFloat32_DExponentScaling: BSCALE / BZERO written with a D
// exponent (5.0D-1, -3.0D0, 1.0D1) are valid FITS reals (FITS 4.0 4.2.4)
// and are applied, on a fresh Open (raw-card parser) and after Header()
// alike. Truth: gen_expected.py (astropy reads 0.5 / -3 / 10, and so does
// cfitsio 4.7.0 fits_read_img, PROVENANCE.md).
func TestReadFloat32_DExponentScaling(t *testing.T) {
	exp := loadRF32Expected(t)
	for _, name := range []string{"i16_dexp_2d.fits", "f32_bzero_dexp_2d.fits"} {
		t.Run(name, func(t *testing.T) {
			e := exp[name]
			if e.Rule == "error" || !e.Bits.Present() {
				t.Fatalf("expected entry %+v has no pixels", e)
			}
			for _, parsed := range []bool{false, true} {
				f, err := Open(imagingPath(name))
				if err != nil {
					t.Fatal(err)
				}
				p, err := f.Primary()
				if err != nil {
					t.Fatal(err)
				}
				if parsed {
					p.Header()
				}
				got, err := ReadFloat32(p)
				f.Close()
				if err != nil {
					t.Fatalf("parsed=%v: ReadFloat32: %v", parsed, err)
				}
				e.Bits.Check(t, fmt.Sprintf("parsed=%v ReadFloat32 (rule %s)", parsed, e.Rule), f32Canon(got))
			}
		})
	}
}

// TestReadFloat32_CompressedInvalidScaling: the compressed path's
// the scaling-card rule rejects a BSCALE / BZERO card that is not a number. The
// sample is comp_i16_bscale_neg.fits.fz with its BSCALE card value replaced
// by a string (header block rewritten in place, same length; the
// compressed data are untouched).
func TestReadFloat32_CompressedInvalidScaling(t *testing.T) {
	src := mustReadFile(t, imagingPath("comp_i16_bscale_neg.fits.fz"))
	for _, c := range []struct{ key, val string }{{"BSCALE", "'abc'"}, {"BZERO", "1.0.0"}} {
		t.Run(c.key, func(t *testing.T) {
			b := append([]byte(nil), src...)
			found := false
			for off := 2880; off+80 <= len(b); off += 80 {
				card := string(b[off : off+80])
				if strings.HasPrefix(card, "END ") {
					break
				}
				if strings.HasPrefix(card, fmt.Sprintf("%-8s= ", c.key)) {
					copy(b[off:off+80], fmt.Sprintf("%-8s= %20s", c.key, c.val)+strings.Repeat(" ", 50))
					found = true
					break
				}
			}
			if !found {
				t.Fatalf("no %s card in the compressed header", c.key)
			}
			path := filepath.Join(t.TempDir(), "bad.fits.fz")
			if err := os.WriteFile(path, b, 0o644); err != nil {
				t.Fatal(err)
			}
			got, err := readFloat32File(t, path)
			if err == nil || !strings.Contains(err.Error(), c.key) {
				t.Fatalf("ReadFloat32 = %v, %v; want an error naming %s", got, err, c.key)
			}
		})
	}
}

func readPhysical(t *testing.T, path string) ([]float64, error) {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	p, err := f.Primary()
	if err != nil {
		return nil, err
	}
	return ReadPixels[float64](p)
}

// TestReadFrameRGB_Planes: the planar u16 RGB sample splits into R, G, B
// planes of the expected values (plane k = pixels k*w*h .. (k+1)*w*h-1).
// comp_u16_rgb.fits.fz (fpack RICE_1 of the same cube) goes through the
// tile-compressed HDU-1 path; its expected values are cfitsio's own pixels
// and null flags (cfnull golden, gen_expected.py), and Siril reads it
// bit-identically (TestReadFloat32_SirilConvention, relation "equal").
func TestReadFrameRGB_Planes(t *testing.T) {
	for _, name := range []string{"u16_rgb.fits", "comp_u16_rgb.fits.fz"} {
		t.Run(name, func(t *testing.T) {
			exp := loadRF32Expected(t)[name]
			want := exp.Bits.Explicit(t)
			if len(want) != 60 {
				t.Fatalf("expected entry has %d values, want 60", len(want))
			}
			r, g, b, w, h, err := ReadFrameRGB(imagingPath(name))
			if err != nil {
				t.Fatal(err)
			}
			if w != 5 || h != 4 {
				t.Fatalf("dims %dx%d, want 5x4", w, h)
			}
			n := w * h
			for k, plane := range [][]float32{r, g, b} {
				for i := 0; i < n; i++ {
					if !f32bitsEq(plane[i], want[k*n+i]) {
						t.Errorf("plane %d pixel %d: got %v", k, i, plane[i])
					}
				}
			}
		})
	}
}

// TestReadFloat32_CompressedFloat_Cfitsio: for tile-compressed float
// images ReadFloat32 passes floats through, so it must equal float32 of
// cfitsio's own decompressed value (testdata/cfref goldens), NaN where
// cfitsio flags a null. Covers the astropy fixtures (quantized
// RICE/GZIP/HCOMPRESS; NO_DITHER, SUBTRACTIVE_DITHER_1 and _2) and every
// float sample of the fpack corpus (comp_fz_*: BITPIX -32 and -64, lossless
// and quantized, BSCALE/BZERO, ZSCALE/ZZERO keywords, ZBLANK column,
// NaN nulls, multi-tile cubes).
//
// One documented exception (testdata/COMP_FZ_PROVENANCE.md, "Fallback
// tiles"): cfitsio's fits_read_img returns a GZIP_COMPRESSED_DATA /
// UNCOMPRESSED_DATA fallback tile of a BSCALE/BZERO image unscaled, while
// funpack, astropy and this package scale it. For those image rows
// (compressedFallbackRows) the truth is cfitsio's read of the funpack twin
// (comp_fzun_*), exactly as TestCompressedCorpusRead does for ReadPixels.
func TestReadFloat32_CompressedFloat_Cfitsio(t *testing.T) {
	compressedFallbackRows := map[string][]int{
		"comp_fz_f32flat_rice_scaled.fits":   {5},
		"comp_fz_f32flat_uncomp_scaled.fits": {5},
		"comp_fz_f64flat_rice_scaled.fits":   {5},
		"comp_fz_f64flat_uncomp_scaled.fits": {5},
		// a NaN in a tile of an image without ZBLANK: fits_read_img does
		// not null-check it (see fzCase.rawRows)
		"comp_fz_f32_rice_zqnone.fits":        {10},
		"comp_fz_f32_rice_zqnone_scaled.fits": {10},
		"comp_fz_f64_rice_zqnone.fits":        {5, 15},
	}
	firstImage := func(f *testutil.CfrefFile) *testutil.CfrefImage {
		for i := range f.HDUs {
			if f.HDUs[i].Image != nil && f.HDUs[i].Image.Npix > 0 {
				return f.HDUs[i].Image
			}
		}
		return nil
	}
	names := []string{
		"comp_rice_f32_dither1.fits", "comp_rice_f32_dither2.fits", "comp_rice_f32_nodither.fits",
		"comp_gzip1_f32.fits", "comp_gzip2_f32.fits", "comp_hcompress_f32.fits",
	}
	for _, pat := range []string{"comp_fz_*f32*.fits", "comp_fz_*f64*.fits"} {
		m, _ := filepath.Glob(testutil.TestdataPath(t, pat))
		if len(m) == 0 {
			t.Fatalf("no committed samples match %s", pat)
		}
		for _, p := range m {
			names = append(names, filepath.Base(p))
		}
	}
	for _, name := range names {
		t.Run(name, func(t *testing.T) {
			want := testutil.LoadCfrefGolden(t, name)
			got, err := readFloat32File(t, testutil.TestdataPath(t, name))
			if err != nil {
				t.Fatal(err)
			}
			img := firstImage(want)
			if img == nil || img.Phys.Len() != len(got) {
				t.Fatalf("cfitsio golden has no matching image (got %d pixels)", len(got))
			}
			// got == float32(cfitsio TDOUBLE) is what is checked below; for a
			// digested golden, resolve that rounding against our reading.
			var phys []testutil.NullFloat
			if img.Bitpix == -64 && img.Phys.Digested() {
				// A digested -64 golden has no float32 digest: resolve the
				// TDOUBLE values against the library's float64 reading.
				phys = img.Phys.Resolve(t, name+" cfitsio TDOUBLE", func() []testutil.NullFloat { return ourPhysView(t, name, 1) })
			} else {
				phys = img.PhysAsFloat32(t, name+" cfitsio TDOUBLE as float32", func() []testutil.NullFloat {
					out := make([]testutil.NullFloat, len(got))
					for i, v := range got {
						if !math.IsNaN(float64(v)) {
							out[i] = testutil.NullFloat{V: float64(v), Valid: true}
						}
					}
					return out
				})
			}
			if rows, ok := compressedFallbackRows[name]; ok {
				twin := firstImage(testutil.LoadCfrefGolden(t, twinFile(t, strings.TrimSuffix(strings.TrimPrefix(name, "comp_fz_"), ".fits"))))
				if twin == nil || twin.Phys.Len() != len(phys) {
					t.Fatalf("funpack twin golden has no matching image")
				}
				w := int(img.Naxes[0])
				for _, r := range rows {
					differ := false
					for i := r * w; i < (r+1)*w; i++ {
						differ = differ || phys[i] != twin.Phys.Values(t)[i]
						phys[i] = twin.Phys.Values(t)[i]
					}
					if !differ {
						t.Fatalf("row %d: cfitsio now agrees with funpack; drop it from compressedFallbackRows", r)
					}
				}
			}
			bad := 0
			for i, p := range phys {
				ok := (!p.Valid && math.IsNaN(float64(got[i]))) || (p.Valid && got[i] == float32(p.V))
				if !ok {
					if bad++; bad <= 5 {
						t.Errorf("pixel %d: ReadFloat32 %v, cfitsio %+v", i, got[i], p)
					}
				}
			}
			if bad > 5 {
				t.Errorf("%d of %d pixels differ", bad, len(got))
			}
		})
	}
}
