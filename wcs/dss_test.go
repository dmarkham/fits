package wcs

import (
	"strings"
	"testing"
)

// dssCards is a complete DSS plate solution (the plate keywords of
// transform's hdrcases/dss_all_terms_deg5, first three AMD terms).
func dssCards(extra ...any) []any {
	return append([]any{
		"NAXIS", int64(2), "NAXIS1", int64(400), "NAXIS2", int64(300),
		"CNPIX1", int64(9120), "CNPIX2", int64(4410),
		"XPIXELSZ", 25.28445, "YPIXELSZ", 25.28445,
		"PPO3", 177500.0, "PPO6", 177500.0,
		"PLTRAH", int64(5), "PLTRAM", int64(34), "PLTRAS", 31.94,
		"PLTDECSN", "-", "PLTDECD", int64(22), "PLTDECM", int64(0), "PLTDECS", 52.2,
		"AMDX1", 67.165566279, "AMDX2", 0.018513587312, "AMDX3", -43.815032712,
		"AMDY1", 67.167031235, "AMDY2", -0.020129061902, "AMDY3", 21.304011803,
	}, extra...)
}

// TestParseDSSRefusals: DSS headers that wcslib 8.10 cannot read the way
// their author meant are errors. Partial plate keywords make wcslib read
// an uninitialised stack value (wcspih.l dsstmp[20] is never initialised),
// so no wcslib golden can exist; a non-zero AMDX14 makes wcslib 8.10
// crash (dssset's error message passes its arguments in the wrong order:
// hdr2sky segfaults on it). The other refusals have wcslib cases in
// transform's hdrcases_reject / hdrcases_stricter (dss_*).
func TestParseDSSRefusals(t *testing.T) {
	drop := func(key string) []any {
		c := dssCards()
		for i := 0; i < len(c); i += 2 {
			if c[i] == key {
				return append(append([]any{}, c[:i]...), c[i+2:]...)
			}
		}
		t.Fatalf("no %s in dssCards", key)
		return nil
	}
	for _, tc := range []struct {
		name  string
		pairs []any
		want  string
	}{
		{"PPO3 missing", drop("PPO3"), "PPO3 is missing"},
		{"CNPIX2 missing", drop("CNPIX2"), "CNPIX2 is missing"},
		{"PLTDECSN missing", drop("PLTDECSN"), "PLTDECSN is missing"},
		{"PLTDECS missing", drop("PLTDECS"), "PLTDECS is missing"},
		{"plate keyword alone", []any{"NAXIS", int64(2), "PLTRAH", int64(5)}, "is missing"},
		{"PLTDECSN not a sign", dssCards("PLTDECSN", "S"), "want '+' or '-'"},
		{"PLTRAS not numeric", dssCards("PLTRAS", "31.94"), "PLTRAS has a non-numeric value"},
		{"AMDX14 non-zero", dssCards("AMDX14", 1e-12), "AMDX14 = 1e-12"},
		{"AMDY20 non-zero", dssCards("AMDY20", -2e-9), "AMDY20 = -2e-09"},
		{"AMDX0", dssCards("AMDX0", 1.0), "AMDX0 is not a DSS coefficient keyword"},
		{"AMDY100", dssCards("AMDY100", 0.0), "AMDY100 is not a DSS coefficient keyword"},
		{"AMDX5 not numeric", dssCards("AMDX5", true), "AMDX5 has a non-numeric value"},
		{"zero plate scale", dssCards("AMDX1", 0.0, "AMDX2", 0.0), "plate scale"},
		{"DSS with SIP", dssCards("A_ORDER", int64(2), "B_ORDER", int64(2), "A_2_0", 1e-6), "both DSS plate-solution and SIP keywords"},
		{"DSS on 3 axes", dssCards("NAXIS", int64(3), "NAXIS3", int64(4)), "DSS plate solution on a 3-axis WCS"},
		{"DSS with NAXIS 0", []any{"NAXIS", int64(0), "AMDX1", 67.0}, "no WCS axes"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			_, err := Parse(newHeader(tc.pairs...))
			if err == nil || !strings.Contains(err.Error(), tc.want) {
				t.Errorf("Parse error = %v, want one mentioning %q", err, tc.want)
			}
		})
	}
}

// TestParseDSSTranslation: the zero AMD terms above 13 are skipped (as in
// wcslib), the translation replaces the header's own CTYPE / CD / CROTA,
// an alternate description is left alone, and PLATEID alone is not DSS.
func TestParseDSSTranslation(t *testing.T) {
	h := newHeader(dssCards("AMDX14", 0.0, "AMDY20", 0.0,
		"CTYPE1", "GLON-CAR", "CTYPE2", "GLAT-CAR", "CD1_1", 1.0, "CD2_2", 1.0, "CROTA2", 10.0,
		"CTYPE1A", "RA---TAN", "CTYPE2A", "DEC--TAN", "CRVAL1A", 1.0, "CDELT1A", 2.0)...)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.DSS == nil || w.CType[0] != "RA---TAN" || w.CType[1] != "DEC--TAN" || w.CD != nil || !w.HasPC || w.LonAxis != 1 || w.LatAxis != 2 {
		t.Fatalf("DSS translation not applied: %+v", w)
	}
	if w.CRVal[1] >= 0 {
		t.Errorf("PLTDECSN '-' gave CRVAL2 = %v", w.CRVal[1])
	}
	a, err := ParseAlt(h, "A")
	if err != nil {
		t.Fatal(err)
	}
	if a.DSS != nil || a.CRVal[0] != 1 || a.CDelt[0] != 2 {
		t.Errorf("alternate A changed by DSS: %+v", a)
	}
	w, err = Parse(newHeader("NAXIS", int64(2), "PLATEID", "0025", "CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN"))
	if err != nil || w.DSS != nil {
		t.Errorf("PLATEID alone: DSS = %v, err = %v; wcslib does not treat it as DSS", w.DSS, err)
	}
}
