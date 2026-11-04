package wcs

import (
	"math"
	"testing"

	"github.com/dmarkham/fits/header"
)

// newHeader builds a *header.Header from a flat list of (key, value) pairs
// for test convenience.
func newHeader(pairs ...any) *header.Header {
	h := header.New()
	for i := 0; i < len(pairs); i += 2 {
		k := pairs[i].(string)
		v := pairs[i+1]
		if err := h.Set(k, v, ""); err != nil {
			panic(err)
		}
	}
	return h
}

// TestParseMinimal: NAXIS=2 with defaults only. Everything should fall to
// sane defaults per Paper I.
func TestParseMinimal(t *testing.T) {
	h := newHeader("NAXIS", int64(2), "NAXIS1", int64(100), "NAXIS2", int64(100))
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.NAxis != 2 {
		t.Fatalf("NAxis=%d", w.NAxis)
	}
	if w.CDelt[0] != 1 || w.CDelt[1] != 1 {
		t.Fatalf("CDELT default not 1: %v", w.CDelt)
	}
	if w.CRPix[0] != 0 || w.CRPix[1] != 0 {
		t.Fatalf("CRPIX default not 0: %v", w.CRPix)
	}
	// PC should default to identity.
	if len(w.PC) != 2 || w.PC[0][0] != 1 || w.PC[1][1] != 1 || w.PC[0][1] != 0 {
		t.Fatalf("PC default not identity: %v", w.PC)
	}
	if w.CD != nil {
		t.Fatalf("CD should be nil when not present")
	}
	if w.IsCelestial() {
		t.Fatalf("minimal header should not be celestial")
	}
}

// TestParseTANHeader parses a typical HST-style TAN header with all
// standard keywords.
func TestParseTANHeader(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"NAXIS1", int64(512),
		"NAXIS2", int64(512),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CUNIT1", "deg",
		"CUNIT2", "deg",
		"CRVAL1", 180.0,
		"CRVAL2", 30.0,
		"CRPIX1", 256.5,
		"CRPIX2", 256.5,
		"CDELT1", -0.0002777778,
		"CDELT2", 0.0002777778,
		"RADESYS", "ICRS",
		"EQUINOX", 2000.0,
		"DATE-OBS", "2025-04-05T12:34:56",
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() {
		t.Fatalf("not celestial")
	}
	if w.LonAxis != 1 || w.LatAxis != 2 {
		t.Fatalf("axes: lon=%d lat=%d", w.LonAxis, w.LatAxis)
	}
	if w.AxisType[0] != "RA" || w.AxisType[1] != "DEC" {
		t.Fatalf("axis types: %v", w.AxisType)
	}
	if w.ProjCode[0] != "TAN" || w.ProjCode[1] != "TAN" {
		t.Fatalf("proj codes: %v", w.ProjCode)
	}
	if w.CelestialProjCode() != "TAN" {
		t.Fatalf("celestial code")
	}
	if w.CRVal[0] != 180.0 || w.CRVal[1] != 30.0 {
		t.Fatalf("CRVAL: %v", w.CRVal)
	}
	if math.Abs(w.CDelt[0]+0.0002777778) > 1e-15 {
		t.Fatalf("CDELT1: %v", w.CDelt[0])
	}
	if w.RadeSys != "ICRS" {
		t.Fatalf("RADESYS: %q", w.RadeSys)
	}
	if w.Equinox != 2000.0 {
		t.Fatalf("EQUINOX: %v", w.Equinox)
	}
	if w.DateObs != "2025-04-05T12:34:56" {
		t.Fatalf("DATE-OBS: %q", w.DateObs)
	}
}

// TestParseCDMatrix: CD matrix form instead of PC+CDELT.
func TestParseCDMatrix(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CRVAL1", 10.0,
		"CRVAL2", 20.0,
		"CRPIX1", 100.0,
		"CRPIX2", 100.0,
		"CD1_1", 1e-5,
		"CD1_2", 2e-6,
		"CD2_1", 3e-6,
		"CD2_2", -1e-5,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.CD == nil {
		t.Fatalf("CD matrix not populated")
	}
	if w.CD[0][0] != 1e-5 || w.CD[0][1] != 2e-6 || w.CD[1][0] != 3e-6 || w.CD[1][1] != -1e-5 {
		t.Fatalf("CD matrix: %v", w.CD)
	}
}

// TestParsePCAndCDBothAccepted: PC and CD both present is accepted, as
// wcslib does (hdr golden matrix_pc_and_cd shows wcslib uses PC + CDELT);
// both matrices are kept and HasPC tells transform.New to prefer PC.
func TestParsePCAndCDBothAccepted(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"PC1_1", 1.0,
		"CD1_1", 1e-5,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatalf("PC + CD must be accepted (wcslib accepts it): %v", err)
	}
	if !w.HasPC || w.CD == nil || w.CD[0][0] != 1e-5 {
		t.Fatalf("HasPC=%v CD=%v", w.HasPC, w.CD)
	}
}

// TestParsePCMatrix: explicit PC matrix (non-identity).
func TestParsePCMatrix(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"CDELT1", 0.01,
		"CDELT2", 0.01,
		"PC1_1", 0.707,
		"PC1_2", -0.707,
		"PC2_1", 0.707,
		"PC2_2", 0.707,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(w.PC[0][0]-0.707) > 1e-15 {
		t.Fatalf("PC[0][0]: %v", w.PC[0][0])
	}
	if math.Abs(w.PC[0][1]+0.707) > 1e-15 {
		t.Fatalf("PC[0][1]: %v", w.PC[0][1])
	}
}

// TestParsePV: PVi_m projection parameters.
func TestParsePV(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---SZP",
		"CTYPE2", "DEC--SZP",
		"PV2_1", 0.5,
		"PV2_2", 1.5,
		"PV2_3", 2.5,
	)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if w.PV[PVKey{Axis: 2, Index: 1}] != 0.5 {
		t.Fatalf("PV2_1: %v", w.PV[PVKey{Axis: 2, Index: 1}])
	}
	if w.PV[PVKey{Axis: 2, Index: 3}] != 2.5 {
		t.Fatalf("PV2_3: %v", w.PV[PVKey{Axis: 2, Index: 3}])
	}
}

// TestSplitCTypeVariants covers the quirky CTYPE formats seen in the wild.
func TestSplitCTypeVariants(t *testing.T) {
	cases := []struct {
		in       string
		axisType string
		projCode string
		distCode string
	}{
		{"RA---TAN", "RA", "TAN", ""},
		{"DEC--TAN", "DEC", "TAN", ""},
		{"GLON-CAR", "GLON", "CAR", ""},
		{"GLAT-CAR", "GLAT", "CAR", ""},
		{"ELON-AIT", "ELON", "AIT", ""},
		{"RA---SIN", "RA", "SIN", ""},
		{"DEC--SIN", "DEC", "SIN", ""},
		{"RA---TAN-SIP", "RA", "TAN", "SIP"},
		{"DEC--TAN-SIP", "DEC", "TAN", "SIP"},
		{"ra---tan-sip", "RA", "TAN", "SIP"},
		{"RA---TPV", "RA", "TPV", ""},
		{"FREQ", "FREQ", "", ""},
		{"WAVE", "WAVE", "", ""},
		{"", "", "", ""},
	}
	for _, c := range cases {
		at, pc, dc := splitCType(c.in)
		if at != c.axisType || pc != c.projCode || dc != c.distCode {
			t.Errorf("splitCType(%q) = (%q,%q,%q), want (%q,%q,%q)",
				c.in, at, pc, dc, c.axisType, c.projCode, c.distCode)
		}
	}
}

// sipCards is a TAN WCS with a quadratic SIP distortion and its inverse.
var sipCards = []any{
	"NAXIS", int64(2),
	"CRVAL1", 101.5, "CRVAL2", 28.9,
	"CRPIX1", 2640.5, "CRPIX2", 2640.5,
	"CD1_1", 4.5e-05, "CD1_2", 0.0, "CD2_1", 0.0, "CD2_2", -4.5e-05,
	"A_ORDER", int64(2), "A_2_0", -5.4e-08, "A_1_1", -1.0e-07,
	"B_ORDER", int64(2), "B_0_2", -7.1e-08,
	"AP_ORDER", int64(2), "AP_2_0", 5.4e-08,
	"BP_ORDER", int64(2), "BP_0_2", 7.1e-08,
}

// TestParseSIPSuffix: "RA---TAN-SIP" is a celestial TAN axis with a SIP
// distortion, and the SIP coefficients are read into Header.SIP.
func TestParseSIPSuffix(t *testing.T) {
	h := newHeader(append([]any{"CTYPE1", "RA---TAN-SIP", "CTYPE2", "DEC--TAN-SIP"}, sipCards...)...)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() {
		t.Fatalf("RA---TAN-SIP not celestial: AxisType=%v", w.AxisType)
	}
	if got := w.CelestialProjCode(); got != "TAN" {
		t.Errorf("CelestialProjCode() = %q, want TAN", got)
	}
	if w.DistCode[0] != "SIP" || w.DistCode[1] != "SIP" {
		t.Errorf("DistCode = %v, want [SIP SIP]", w.DistCode)
	}
	if w.SIP == nil {
		t.Fatal("SIP coefficients not parsed")
	}
	if w.SIP.A.Order != 2 || w.SIP.A.Coeffs[2][0] != -5.4e-08 || w.SIP.A.Coeffs[1][1] != -1.0e-07 {
		t.Errorf("A = %+v", w.SIP.A)
	}
	if w.SIP.B.Coeffs[0][2] != -7.1e-08 {
		t.Errorf("B_0_2 = %v", w.SIP.B.Coeffs[0][2])
	}
	if w.SIP.AP == nil || w.SIP.AP.Coeffs[2][0] != 5.4e-08 || w.SIP.BP == nil || w.SIP.BP.Coeffs[0][2] != 7.1e-08 {
		t.Errorf("AP/BP = %+v / %+v", w.SIP.AP, w.SIP.BP)
	}
}

// TestParseSIPKeysWithoutSuffix: SIP keywords are read for the primary
// description of a TAN or SIN (any non-TPV) header whatever CTYPE says,
// as wcslib 8.10 wcspih does (it sets the SIP prior distortion from the
// A_p_q / B_p_q keywords alone, wcspih.l <SIP2>/<SIP3>). The wcslib truth
// for the transform is in the header-pipeline goldens sip_nosuffix_tan and
// sip_nosuffix_sin (TestHdrGoldens); this test pins the parse. A RA---TPV
// header is the exception (TestParseTPVIgnoresSIPWithoutSuffix). SIP has
// no alternate form: an alternate description gets no SIP (golden
// sip_nosuffix_alt_a).
func TestParseSIPKeysWithoutSuffix(t *testing.T) {
	for _, ct := range [][2]string{{"RA---TAN", "DEC--TAN"}, {"RA---SIN", "DEC--SIN"}, {"RA---ZEA", "DEC--ZEA"}} {
		h := newHeader(append([]any{"CTYPE1", ct[0], "CTYPE2", ct[1]}, sipCards...)...)
		w, err := Parse(h)
		if err != nil {
			t.Fatal(err)
		}
		if w.SIP == nil || w.SIP.A == nil || w.SIP.B == nil {
			t.Errorf("%s: SIP keywords without a -SIP suffix not read (wcslib applies them)", ct[0])
		}
		if w.DistCode[0] != "" {
			t.Errorf("%s: DistCode = %q, want empty", ct[0], w.DistCode[0])
		}
		if w.SIPIgnored {
			t.Errorf("%s: SIPIgnored set; the SIP keywords apply", ct[0])
		}
	}
	// Alternate description: no SIP, and no error.
	h := newHeader(append([]any{"CTYPE1", "RA---TAN", "CTYPE2", "DEC--TAN", "CTYPE1A", "RA---TAN", "CTYPE2A", "DEC--TAN"}, sipCards...)...)
	w, err := ParseAlt(h, "A")
	if err != nil {
		t.Fatal(err)
	}
	if w.SIP != nil || len(w.UnsupportedDistortion) != 0 {
		t.Errorf("alternate A: SIP = %+v, unsupported %v; want undistorted (wcslib)", w.SIP, w.UnsupportedDistortion)
	}
}

// TestParseTPVIgnoresSIPWithoutSuffix: on a header whose two celestial
// axes are RA---TPV / DEC--TPV without a "-SIP" suffix the SIP keywords
// are not read (Header.SIP nil, SIPIgnored set), as astropy 8.0.2 drops
// them ("Removed redundant SIP distortion parameters because CTYPE
// explicitly specifies TPV distortions"); wcslib 8.10 applies them on top
// of TPV, which on the real ATLAS header (hdrcases_differ/
// real_atlas_k2_tpv_sip, TestATLASTPVPlusRedundantSIP) moves every
// position by 13-14 px. The SIP set is ignored whatever it holds, so an
// order-5 A/B/AP/BP set with linear AP/BP terms (the ATLAS form) and even
// a set this package would refuse on a TAN header (no B_ORDER, a
// coefficient above its order, degree 10) is no error. "-SIP" written on
// both axes (RA---TPV-SIP) declares the SIP set and it is read; "-SIP" on
// one axis only is read too (transform.New refuses the inconsistent
// distortion codes).
func TestParseTPVIgnoresSIPWithoutSuffix(t *testing.T) {
	atlasSIP := []any{
		"A_ORDER", int64(5), "A_2_0", -8.592785425686e-09, "A_5_0", -2.254046336565e-19,
		"B_ORDER", int64(5), "B_0_2", 1.011952116223e-08, "B_0_5", -1.783421432784e-19,
		"AP_ORDER", int64(5), "AP_1_0", 5.011756729401e-06, "AP_0_5", 5.858779290584e-22,
		"BP_ORDER", int64(5), "BP_0_1", 4.615336221203e-06, "BP_5_0", 1.780884109469e-21,
	}
	malformed := []any{"A_ORDER", int64(2), "A_3_0", 1e-8, "A_5_5", 0.0} // no B_ORDER, A_3_0 above A_ORDER, degree 10
	tpv := []any{"NAXIS", int64(2), "CTYPE1", "RA---TPV", "CTYPE2", "DEC--TPV",
		"CRVAL1", 272.84895, "CRVAL2", 8.28673, "CRPIX1", -3799.0, "CRPIX2", 1852.0,
		"CD1_1", -0.0005172389232442, "CD2_2", 0.0005171961669467, "PV1_1", 1.0, "PV2_1", 1.0, "PV1_7", -6.465683213014e-04}
	for name, set := range map[string][]any{"atlas order 5": atlasSIP, "sipCards": sipCards[2:], "malformed": malformed} {
		w, err := Parse(newHeader(append(append([]any{}, tpv...), set...)...))
		if err != nil {
			t.Errorf("RA---TPV + %s SIP set: %v; the SIP keywords are ignored, their presence is no error", name, err)
			continue
		}
		if w.SIP != nil || !w.SIPIgnored {
			t.Errorf("RA---TPV + %s SIP set: SIP = %+v, SIPIgnored = %v; want nil, true (astropy drops it)", name, w.SIP, w.SIPIgnored)
		}
		if w.CelestialProjCode() != "TPV" || w.PV[PVKey{1, 7}] != -6.465683213014e-04 {
			t.Errorf("RA---TPV + %s SIP set: projection %q, PV1_7 %v: the TPV solution must stay", name, w.CelestialProjCode(), w.PV[PVKey{1, 7}])
		}
	}
	// No SIP keyword: nothing ignored.
	if w, err := Parse(newHeader(tpv...)); err != nil || w.SIPIgnored || w.SIP != nil {
		t.Errorf("RA---TPV without SIP keywords: SIP %+v, SIPIgnored %v, err %v", w.SIP, w.SIPIgnored, err)
	}
	// Explicit -SIP: the SIP set is declared and read.
	for _, ct := range [][2]string{{"RA---TPV-SIP", "DEC--TPV-SIP"}, {"RA---TPV-SIP", "DEC--TPV"}, {"RA---TPV", "DEC--TPV-SIP"}} {
		h := newHeader(append(append(append([]any{}, tpv...), "CTYPE1", ct[0], "CTYPE2", ct[1]), atlasSIP...)...)
		w, err := Parse(h)
		if err != nil {
			t.Fatalf("%s / %s: %v", ct[0], ct[1], err)
		}
		if w.SIP == nil || w.SIPIgnored || w.SIP.A.Order != 5 || w.SIP.AP.Coeffs[1][0] != 5.011756729401e-06 {
			t.Errorf("%s / %s: SIP = %+v, SIPIgnored %v; want the order-5 set read", ct[0], ct[1], w.SIP, w.SIPIgnored)
		}
	}
	// The explicit suffix does not relax the SIP checks.
	h := newHeader(append(append(append([]any{}, tpv...), "CTYPE1", "RA---TPV-SIP", "CTYPE2", "DEC--TPV-SIP"), malformed...)...)
	if _, err := Parse(h); err == nil {
		t.Error("RA---TPV-SIP with a malformed SIP set parsed; it declares SIP, so the SIP checks apply")
	}
}

// TestParseSIPSuffixWithoutCoefficients: a "-SIP" CTYPE with no A_ORDER /
// B_ORDER parses as celestial with no SIP polynomial (as astropy.wcs).
func TestParseSIPSuffixWithoutCoefficients(t *testing.T) {
	h := newHeader("NAXIS", int64(2), "CTYPE1", "RA---TAN-SIP", "CTYPE2", "DEC--TAN-SIP",
		"CRVAL1", 10.0, "CRVAL2", 20.0, "CDELT1", -1e-3, "CDELT2", 1e-3)
	w, err := Parse(h)
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() || w.CelestialProjCode() != "TAN" {
		t.Fatalf("not a celestial TAN header: %v", w)
	}
	if w.DistCode[0] != "SIP" {
		t.Errorf("DistCode = %v", w.DistCode)
	}
	if w.SIP != nil {
		t.Errorf("SIP = %+v, want nil", w.SIP)
	}
}

// TestParseRADECSYSFallback: the pre-1998 RADECSYS spelling should still be
// recognized.
func TestParseRADECSYSFallback(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"RADECSYS", "FK5",
	)
	w, _ := Parse(h)
	if w.RadeSys != "FK5" {
		t.Fatalf("RADECSYS fallback: %q", w.RadeSys)
	}
}

// TestParseEpochFallback: pre-1988 EPOCH is a synonym for EQUINOX.
func TestParseEpochFallback(t *testing.T) {
	h := newHeader(
		"NAXIS", int64(2),
		"CTYPE1", "RA---TAN",
		"CTYPE2", "DEC--TAN",
		"EPOCH", 1950.0,
	)
	w, _ := Parse(h)
	if w.Equinox != 1950.0 {
		t.Fatalf("EPOCH fallback: %v", w.Equinox)
	}
}
