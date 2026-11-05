package wcs

import (
	"strings"
	"testing"
)

// TestParseSIPOrderValidation (PR #3 review): A_ORDER/B_ORDER/AP_ORDER/
// BP_ORDER must be integers in 0..MaxSIPOrder. A negative order used to
// panic (make with a negative length) and a huge one allocated an
// order^2 grid; both must now be ordinary errors from Parse.
func TestParseSIPOrderValidation(t *testing.T) {
	base := func(extra ...any) []any {
		return append([]any{
			"NAXIS", int64(2),
			"CTYPE1", "RA---TAN-SIP", "CTYPE2", "DEC--TAN-SIP",
			"CDELT1", -1e-4, "CDELT2", 1e-4,
		}, extra...)
	}
	for _, tc := range []struct {
		name  string
		pairs []any
		want  string
	}{
		{"A_ORDER negative", base("A_ORDER", int64(-1), "B_ORDER", int64(2)), "A_ORDER = -1"},
		{"B_ORDER negative", base("A_ORDER", int64(2), "B_ORDER", int64(-3)), "B_ORDER = -3"},
		{"A_ORDER too large", base("A_ORDER", int64(10), "B_ORDER", int64(2)), "A_ORDER = 10"},
		{"A_ORDER huge", base("A_ORDER", int64(1<<40), "B_ORDER", int64(2)), "outside 0..9"},
		{"AP_ORDER negative", base("A_ORDER", int64(2), "B_ORDER", int64(2), "AP_ORDER", int64(-1), "BP_ORDER", int64(1)), "AP_ORDER = -1"},
		{"A_ORDER not integer", base("A_ORDER", 2.5, "B_ORDER", int64(2)), "A_ORDER"},
		{"coefficient beyond order", base("A_ORDER", int64(2), "B_ORDER", int64(2), "A_3_0", 1e-9), "A_3_0"},
		// "-SIP" with SIP keywords but a missing A_ORDER/B_ORDER used to be
		// read as undistorted, dropping the declared polynomial. Refusing
		// is our policy; wcslib 8.10 accepts these when AP/BP are present
		// (transform hdrcases_stricter/sip_a_without_b,
		// sip_coeffs_without_orders, TestHdrStricterThanWcslib).
		{"A without B_ORDER", base("A_ORDER", int64(2), "A_2_0", 1e-4), "B_ORDER is missing"},
		{"B without A_ORDER", base("B_ORDER", int64(2), "B_0_2", 1e-4), "A_ORDER is missing"},
		{"coefficient without any order", base("A_2_0", 1e-4), "A_ORDER is missing"},
		{"AP only", base("AP_ORDER", int64(2), "AP_2_0", 1e-4), "A_ORDER is missing"},
		// A coefficient that is not a number is an error, not a silent 0.
		{"A coefficient string", base("A_ORDER", int64(2), "B_ORDER", int64(2), "A_2_0", "1e-5"), "A_2_0 has a non-numeric value"},
		{"B coefficient logical", base("A_ORDER", int64(2), "B_ORDER", int64(2), "B_0_2", true), "B_0_2 has a non-numeric value"},
		{"AP coefficient string", base("A_ORDER", int64(2), "B_ORDER", int64(2), "AP_ORDER", int64(1), "BP_ORDER", int64(1), "AP_1_0", "x"), "AP_1_0 has a non-numeric value"},
		{"BP coefficient logical", base("A_ORDER", int64(2), "B_ORDER", int64(2), "AP_ORDER", int64(1), "BP_ORDER", int64(1), "BP_0_1", false), "BP_0_1 has a non-numeric value"},
		// wcslib's lexer takes single-digit SIP indices only (wcspih.l
		// <SIP2>{Z1}_{Z1}) and silently discards A_01_1, B_0_02, ...;
		// applying them as A_1_1 differs from wcslib and dropping them
		// loses a coefficient, so they are refused (transform
		// hdrcases_stricter/sip_index_leading_zero).
		{"A index zero-padded", base("A_ORDER", int64(2), "B_ORDER", int64(2), "A_01_1", 1e-5), "A_01_1 is not a SIP coefficient keyword"},
		{"B index zero-padded", base("A_ORDER", int64(2), "B_ORDER", int64(2), "B_0_02", 1e-5), "B_0_02 is not a SIP coefficient keyword"},
		{"A index two digits", base("A_ORDER", int64(9), "B_ORDER", int64(2), "A_1_10", 0.0), "A_1_10 is not a SIP coefficient keyword"},
		{"AP index zero-padded", base("A_ORDER", int64(2), "B_ORDER", int64(2), "AP_ORDER", int64(1), "BP_ORDER", int64(1), "AP_00_1", 1e-5), "AP_00_1 is not a SIP coefficient keyword"},
		{"BP index zero-padded", base("A_ORDER", int64(2), "B_ORDER", int64(2), "AP_ORDER", int64(1), "BP_ORDER", int64(1), "BP_1_00", 1e-5), "BP_1_00 is not a SIP coefficient keyword"},
		{"zero-padded index only", base("A_01_1", 1e-5), "A_01_1 is not a SIP coefficient keyword"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			defer func() {
				if r := recover(); r != nil {
					t.Fatalf("Parse panicked: %v", r)
				}
			}()
			_, err := Parse(newHeader(tc.pairs...))
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("Parse error = %v, want one mentioning %q", err, tc.want)
			}
		})
	}
	// Order 0 and 9 are valid; a zero coefficient beyond the order is
	// harmless and accepted.
	for _, ok := range [][]any{
		base("A_ORDER", int64(0), "B_ORDER", int64(0)),
		base("A_ORDER", int64(9), "B_ORDER", int64(9), "A_9_0", 1e-30),
		base("A_ORDER", int64(2), "B_ORDER", int64(2), "A_3_0", 0.0),
	} {
		w, err := Parse(newHeader(ok...))
		if err != nil || w.SIP == nil {
			t.Errorf("valid SIP rejected: %v", err)
		}
	}
	// "-SIP" with no SIP keyword at all is undistorted (wcslib agrees:
	// hdrcases/sip_suffix_no_coeffs).
	if w, err := Parse(newHeader(base()...)); err != nil || w.SIP != nil {
		t.Errorf("-SIP without SIP keywords: got %+v, err = %v; want SIP nil, no error", w, err)
	}
}

// TestParseWrongTypeIsAnError: a WCS keyword present with a value of the
// wrong type is reported, not silently replaced by its default (which
// would put the image somewhere else on the sky without a word).
func TestParseWrongTypeIsAnError(t *testing.T) {
	for _, tc := range [][]any{
		{"CRVAL1", "abc"},
		{"CDELT2", true},
		{"CTYPE1", 3.0},
		{"PC1_2", "x"},
		{"CD2_1", "x"},
		{"PV2_1", "x"},
		{"LONPOLE", "180"},
		{"WCSAXES", 2.0},
	} {
		pairs := []any{"NAXIS", int64(2), "CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN"}
		pairs = append(pairs, tc...)
		_, err := Parse(newHeader(pairs...))
		if err == nil || !strings.Contains(err.Error(), tc[0].(string)) {
			t.Errorf("%s = %v: Parse error = %v, want an error naming the keyword", tc[0], tc[1], err)
		}
	}
}

// TestParseNAxisRules pins detectNAxis (Paper I §2.2): WCSAXESa wins;
// otherwise max(NAXIS, highest axis index on this description's WCS
// keywords). The wcslib-backed versions are hdr goldens
// axes_wcsaxes3_naxis2, axes_index_above_naxis and alt_a_celestial
// (TestHdrDerivedParameters compares NAxis with wcslib's naxis).
func TestParseNAxisRules(t *testing.T) {
	for _, tc := range []struct {
		name  string
		alt   string
		pairs []any
		want  int
	}{
		{"NAXIS only", "", []any{"NAXIS", int64(2)}, 2},
		{"WCSAXES above NAXIS", "", []any{"NAXIS", int64(2), "WCSAXES", int64(4)}, 4},
		{"WCSAXES below NAXIS", "", []any{"NAXIS", int64(3), "WCSAXES", int64(2), "CTYPE3", "FREQ"}, 2},
		{"index above NAXIS", "", []any{"NAXIS", int64(2), "CRVAL3", 1.0}, 3},
		{"PC index above NAXIS", "", []any{"NAXIS", int64(2), "PC4_1", 0.0}, 4},
		{"no NAXIS", "", []any{"CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN"}, 2},
		{"alt uses WCSAXESA, not WCSAXES", "A", []any{"NAXIS", int64(2), "WCSAXES", int64(5), "WCSAXESA", int64(3), "CTYPE1A", "X"}, 3},
		{"alt without WCSAXESA ignores WCSAXES", "A", []any{"NAXIS", int64(2), "WCSAXES", int64(5), "CTYPE1A", "X"}, 2},
		{"primary ignores alt-only index", "", []any{"NAXIS", int64(2), "CTYPE4A", "X"}, 2},
	} {
		t.Run(tc.name, func(t *testing.T) {
			w, err := ParseAlt(newHeader(tc.pairs...), tc.alt)
			if err != nil {
				t.Fatal(err)
			}
			if w.NAxis != tc.want {
				t.Fatalf("NAxis = %d, want %d", w.NAxis, tc.want)
			}
		})
	}
}

// TestParseCelestialPairRules: two longitude axes (or two latitude axes)
// do not form a pair; the header must not silently pick one (wcslib:
// "Inconsistent projection types", hdr reject case two_lon_axes).
func TestParseCelestialPairRules(t *testing.T) {
	w, err := Parse(newHeader("NAXIS", int64(3), "CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN", "CTYPE3", "GLON-TAN"))
	if err != nil {
		t.Fatal(err)
	}
	if w.IsCelestial() || w.CelestialError == "" {
		t.Fatalf("RA + DEC + GLON: IsCelestial=%v CelestialError=%q", w.IsCelestial(), w.CelestialError)
	}
	for _, ct := range [][2]string{{"GLON-CAR", "GLAT-CAR"}, {"ELON-TAN", "ELAT-TAN"}, {"SLON-TAN", "SLAT-TAN"}, {"HLON-TAN", "HLAT-TAN"}, {"DEC--TAN", "RA---TAN"}} {
		w, err := Parse(newHeader("NAXIS", int64(2), "CTYPE1", ct[0], "CTYPE2", ct[1]))
		if err != nil || !w.IsCelestial() {
			t.Errorf("%v: not celestial (%v)", ct, err)
		}
	}
}

// TestParseUnsupportedDistortionListed: Paper IV / HST lookup-table
// distortion keywords and SIP on an alternate description are recorded so
// transform.New can refuse them.
func TestParseUnsupportedDistortionListed(t *testing.T) {
	w, err := Parse(newHeader("NAXIS", int64(2), "CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN",
		"CPDIS1", "LOOKUP", "D2IMDIS2", "LOOKUP", "CPDIS1A", "LOOKUP"))
	if err != nil {
		t.Fatal(err)
	}
	if got := strings.Join(w.UnsupportedDistortion, ","); got != "CPDIS1,D2IMDIS2" {
		t.Errorf("primary UnsupportedDistortion = %q, want CPDIS1,D2IMDIS2", got)
	}
	w, err = ParseAlt(newHeader("NAXIS", int64(2), "CTYPE1O", "RA---TAN-SIP", "CTYPE2O", "DEC--TAN-SIP",
		"A_ORDER", int64(2), "B_ORDER", int64(2)), "O")
	if err != nil {
		t.Fatal(err)
	}
	if w.SIP != nil || len(w.UnsupportedDistortion) != 1 || !strings.Contains(w.UnsupportedDistortion[0], "SIP") {
		t.Errorf("alt O -SIP: SIP=%v UnsupportedDistortion=%q", w.SIP, w.UnsupportedDistortion)
	}
}

// TestParseDefunctKeywords: PC00i00j / CD00i00j / PROJPn are read as
// PCi_j / CDi_j / PVi_n of the latitude axis (their wcslib goldens are
// hdrcases/legacy_*.case and real_ps1_pc00i00j.case); an alternate
// description does not see them; PROJPn without a latitude axis and a
// defunct card that disagrees with its modern twin are errors.
func TestParseDefunctKeywords(t *testing.T) {
	w, err := Parse(newHeader("NAXIS", int64(2), "CTYPE1", "RA---AZP", "CTYPE2", "DEC--AZP",
		"PC001002", 0.25, "CD002001", -0.5, "PROJP1", 2.0, "PC1_2", 0.25))
	if err != nil {
		t.Fatal(err)
	}
	if !w.HasPC || w.PC[0][1] != 0.25 || w.CD == nil || w.CD[1][0] != -0.5 || w.PV[PVKey{Axis: 2, Index: 1}] != 2 {
		t.Errorf("PC %v CD %v PV %v", w.PC, w.CD, w.PV)
	}
	a, err := ParseAlt(newHeader("NAXIS", int64(2), "CTYPE1A", "RA---TAN", "CTYPE2A", "DEC--TAN", "PC001002", 0.25), "A")
	if err != nil || a.HasPC {
		t.Errorf("alternate description picked up PC001002: %v %v", a, err)
	}
	for _, h := range []struct {
		name  string
		pairs []any
	}{
		{"projp-no-lat", []any{"NAXIS", int64(1), "CTYPE1", "FREQ", "PROJP1", 2.0}},
		{"pc-conflict", []any{"NAXIS", int64(2), "CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN", "PC1_1", 1.0, "PC001001", 0.9}},
		{"cd-conflict", []any{"NAXIS", int64(2), "CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN", "CD2_2", 1.0, "CD002002", 0.9}},
		{"projp-conflict", []any{"NAXIS", int64(2), "CTYPE1", "RA---AZP", "CTYPE2", "DEC--AZP", "PV2_1", 1.0, "PROJP1", 2.0}},
		{"projp-string", []any{"NAXIS", int64(2), "CTYPE1", "RA---AZP", "CTYPE2", "DEC--AZP", "PROJP1", "2"}},
	} {
		if _, err := Parse(newHeader(h.pairs...)); err == nil {
			t.Errorf("%s: Parse accepted it", h.name)
		}
	}
}
