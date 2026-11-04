package transform

import (
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Header cases required BY NAME. The ledger (docs/VALIDATION.md) cites the
// keyword-rule cases through globs (cunit_*.case, ...), which any one file
// satisfies; several rules are proven by a single case (tpv_radial_cd is the
// only case that pins the TPV r-term Jacobian through the pipeline,
// celpole_car_eq_latpole_tie the only one for the LATPOLE tie, ...). So
// every hand-written case is listed here with what it pins, and
// TestHdrRequiredCases fails when a listed case (or its golden) is
// missing, or when a case on disk is not listed. The generated cases
// (prj_*, prjswap_*, pvfid_*, sample_*, hiorder_*) are required by their
// generators' reproduce tests, which diff the whole generated set.

// hdrRequiredCases: hdrcases/<name>.case with golden/hdr/<name>.json,
// checked against wcslib by TestHdrGoldens.
var hdrRequiredCases = map[string]string{
	"dss_all_terms_deg5":           "DSS plate solution with every AMD term 1..13 (degree-5 TPD form)",
	"dss_far_outside_plate":        "DSS degree 5 at 6000 px: wcslib's s2p fails, ours must return the source pixel",
	"dss_modern_cutout":            "DSS plate keywords replace an approximate TAN CD WCS (wcspih_final)",
	"real_dss_plate_solution":      "real archive header (realhdr/make_cases.py)",
	"celpole_z0_latpole95":         "celset z = 0 (theta0 = 0, |phip - phi0| = 90): LATPOLE 95 clamped to 90",
	"celpole_z0_latpolem95":        "celset z = 0: LATPOLE -95 clamped to -90",
	"lonpole_999_default":          "LONPOLE = 999 means the default phip (wcslib celset)",
	"alias_gls_c30":                "alias GLS: RA---GLS, CRVAL2 = 30",
	"alias_ncp_c60":                "alias NCP: RA---NCP, CRVAL2 = 60",
	"alias_ncp_cm25_pv":            "alias NCP: DEC = -25, PV2_1/PV2_2 overwritten",
	"alt_a_cd_lonpole":             "alternate WCS A, CDi_jA and LONPOLEA on CAR",
	"alt_a_celestial":              "alternate WCS A with LONPOLEa, PVi_ma, WCSAXESa",
	"alt_a_no_leak":                "alternate WCS A, primary LONPOLE/LATPOLE/PC/CD/PV must not leak",
	"alt_a_no_wcsaxesa":            "WCSAXES of the primary does not apply to alternate A",
	"axes_cube_freq_ra_dec":        "celestial axes at 2,3 behind a FREQ axis",
	"axes_index_above_naxis":       "WCS axis count from the highest keyword index",
	"axes_swapped":                 "swapped celestial axes (DEC on axis 1)",
	"axes_swapped_pc_cdelt":        "swapped celestial axes, PC + unequal CDELT",
	"axes_wcsaxes3_naxis2":         "WCSAXES > NAXIS",
	"car_allsky":                   "CAR all-sky, default LONPOLE, grid past the map edges",
	"celpole_car_eq_latpole_tie":   "LATPOLE equidistant from both native-pole solutions",
	"celpole_car_native_pole_sky":  "sky points exactly at both native poles of CAR (wcslib sphs2x longitude-origin convention)",
	"celpole_car_eq_south":         "native pole at the celestial south pole (celset z = 0, latp < 0)",
	"celpole_latp1_tol_out":        "celset: nearer native-pole latitude 90 + 1e-9 passed over (tolerance 1e-10)",
	"celpole_latp2_tol_out":        "celset: nearer second native-pole latitude -90 - 1e-9 passed over",
	"celpole_slz_within_tol":       "celset: slz = 1 + 5e-11 clamped to 1 (tolerance 1e-10); sky grid hits the native pole (wcslib sphs2x origin rule)",
	"celpole_z_above_tol":          "celset: z = 1.75e-10 just above the 1e-10 tolerance, general branch",
	"celpole_z_within_tol":         "celset: z = 1.6e-11 within the tolerance, pole shortcut branch",
	"crota2_tan":                   "CROTA2 on TAN",
	"crota2_unequal_cdelt":         "CROTA2 with unequal CDELT, CROTA1 ignored",
	"crota2_with_cd":               "CD wins over CROTA2",
	"crota2_with_pc":               "PC wins over CROTA2",
	"cunit_arcmin":                 "CUNIT arcmin",
	"cunit_arcsec":                 "CUNIT arcsec with PC + CDELT",
	"cunit_arcsec_cd":              "CUNIT arcsec with CD",
	"cunit_mas":                    "CUNIT mas",
	"cunit_mixed_arcmin_deg":       "per-axis CUNIT, arcmin longitude with blank (deg) latitude",
	"cunit_mixed_deg_arcsec":       "per-axis CUNIT, deg longitude with arcsec latitude",
	"cunit_rad":                    "CUNIT rad",
	"frame_elon_tan":               "ELON/ELAT",
	"frame_glon_car_lp0":           "GLON-CAR, LONPOLE = 0",
	"frame_hlon_zea":               "helioecliptic HLON/HLAT with ZEA",
	"frame_slon_tan":               "supergalactic SLON/SLAT with TAN",
	"latpole_car_two_solutions":    "LATPOLE choosing between two native-pole solutions",
	"latpole_coe_two_solutions":    "conic with LONPOLE and LATPOLE both explicit",
	"legacy_cd00i00j":              "defunct CD00i00j read as CDi_j",
	"legacy_cd00i00j_cdelt":        "defunct CD00i00j with CDELT (CDELT ignored)",
	"legacy_pc00i00j_rot":          "defunct PC00i00j rotation read as PCi_j",
	"legacy_projp_azp":             "defunct PROJP1/PROJP2 as PV2_1/PV2_2 on AZP",
	"legacy_projp_zpn":             "defunct PROJP1/PROJP3 on ZPN",
	"lonpole_ait_gal_offset":       "AIT with CRVAL2 = 20, default LONPOLE/LATPOLE",
	"lonpole_car_dec30_lp0":        "LONPOLE = 0 on a cylindrical projection, CRVAL2 = 30",
	"lonpole_car_decm30_default":   "default LONPOLE on CAR with CRVAL2 < theta_0",
	"lonpole_cea_lp90":             "LONPOLE = 90 on CEA, pole latitude from LATPOLE only",
	"lonpole_glon_ait_lp0":         "GLON-AIT all-sky with LONPOLE = 0",
	"lonpole_mol_dec20_lp0":        "MOL with LONPOLE = 0",
	"lonpole_tan_lp150":            "zenithal TAN with non-default LONPOLE",
	"lonpole_tan_southpole":        "zenithal TAN at CRVAL2 = -89.5, default LONPOLE",
	"lonpole_tan_dec40_lp0":        "explicit non-default LONPOLE = 0 on TAN at CRVAL2 = 40 (default 180) honoured",
	"pv1off_car_pv1_1_only":        "PV1_0 = 1 with only PV1_1: no fiducial offset (wcslib prjoff needs both)",
	"pv1off_car_pv1_2_only":        "PV1_0 = 1 with only PV1_2: no fiducial offset (wcslib prjoff needs both)",
	"matrix_cd_diag_missing":       "CD matrix with a missing diagonal element defaults to 0",
	"matrix_cd_partial":            "partial CD matrix",
	"matrix_pc_and_cd":             "PC and CD both present, PC wins",
	"pvlon_ait_phi0_offset":        "PVi_0..2 offset on a non-cylindrical projection (phi_0 matters)",
	"pvlon_car_offset":             "PVi_0..2 on the longitude axis, fiducial offset applied",
	"pvlon_car_offset_phi0_only":   "PV1_0 = 1 with PV1_1 only: no offset, phi_0 kept for the pole",
	"pvlon_car_offset_theta0_only": "PV1_0 = 1 with PV1_2 only: no offset, theta_0 kept for the pole",
	"pvlon_car_pv3_pv4":            "PVi_3 / PVi_4 on the longitude axis override LONPOLE / LATPOLE",
	"pvlon_car_theta0_tol":         "theta_0 just above 90 within the celset tolerance is clamped",
	"pvlon_tan_phi0":               "PV1_1 (phi_0) without offset on TAN",
	"pvlon_tan_theta0":             "PVi_2 (theta_0) on the longitude axis without offset",
	"real_anet_sip3":               "real archive header (realhdr/make_cases.py)",
	"real_anet_sip3_wide":          "real archive header (realhdr/make_cases.py)",
	"real_anet_sip4":               "real archive header (realhdr/make_cases.py)",
	"real_hst_acs_wfc_sip4":        "real archive header (realhdr/make_cases.py)",
	"real_hst_wfc3_uvis_sip4":      "real archive header (realhdr/make_cases.py)",
	"real_jwst_nircam_sip4":        "real archive header (realhdr/make_cases.py)",
	"real_jwst_nircam_sip5":        "real archive header (realhdr/make_cases.py)",
	"real_msx_glon_car_lp0":        "real archive header (realhdr/make_cases.py)",
	"real_ps1_pc00i00j":            "real archive header (realhdr/make_cases.py)",
	"real_ptf_tpv4":                "real archive header (realhdr/make_cases.py)",
	"real_atlas_k2_tpv":            "real ATLAS header (PR #4) without its SIP cards: the TPV answer of real_atlas_k2_tpv_sip",
	"real_atlas_k2_tansip":         "real ATLAS header's SIP set alone as TAN-SIP: equals TPV (TestATLASSIPCopiesTPV)",
	"real_rosat_glon_ait_lp0":      "real archive header (realhdr/make_cases.py)",
	"real_spitzer_irac_cryo_sip2":  "real archive header (realhdr/make_cases.py)",
	"real_spitzer_irac_warm_sip5":  "real archive header (realhdr/make_cases.py)",
	"real_ztf_tpv4":                "real archive header (realhdr/make_cases.py)",
	"sip4_no_apbp":                 "SIP order 4 without AP/BP (Newton inverse)",
	"sip4_rough_apbp":              "SIP AP/BP as starting guess only",
	"sip_nosuffix_alt_a":           "SIP keywords do not reach an alternate description without -SIP",
	"sip_nosuffix_primary_of_alt":  "primary of a header with an alternate: SIP keywords applied",
	"sip_nosuffix_sin":             "SIP keywords without the -SIP suffix are applied (SIN)",
	"sip_nosuffix_tan":             "SIP keywords without the -SIP suffix are applied (TAN)",
	"sip_cube_radec_freq":          "SIP on a 3-axis RA,DEC,FREQ cube (celestial pixel axes 1,2)",
	"sip_on_sin":                   "-SIP on SIN",
	"sip_swapped_nosuffix":         "SIP with DEC on pixel axis 1, no -SIP suffix (A stays on pixel axis 1)",
	"sip_swapped_suffix":           "SIP with DEC on pixel axis 1, -SIP suffix (A stays on pixel axis 1)",
	"sip_pc_cdelt":                 "SIP with PC + CDELT",
	"sip_suffix_no_coeffs":         "-SIP CTYPE without coefficients is undistorted",
	"tan_cd_rot":                   "TAN + rotated CD, no distortion, default LONPOLE (180)",
	"tan_sip3_apbp":                "TAN-SIP order 3 with AP/BP (astropy sip.fits)",
	"tpv_explicit_pv1":             "TPV with explicit PV1_1 = PV2_1 = 1 (registry default)",
	"tpv_pc_cdelt":                 "TPV with PC + CDELT, polynomial applied before CDELT (wcslib)",
	"tpv_radial_cd":                "TPV radial terms r, r^3, r^5, r^7 on both axes (CD form)",
	"tpv_sip_nosuffix_tpvonly":     "hdrcases_differ/tpv_sip_nosuffix without its SIP cards: the answer it must give",
	"tpv_sip_suffix":               "RA---TPV-SIP applies SIP then TPV",
	"tpv_sip_suffix_atlas":         "RA---TPV-SIP on the real ATLAS TPV solution with a non-redundant SIP: both apply",
}

// hdrRequiredReject: hdrcases_reject/<name>.case (TestHdrRejectCases), with
// the wcslib error text each one expects.
var hdrRequiredReject = map[string]string{
	"dss_with_sip":              "No distortion type, yet DPja keyvalues are present",
	"dss_zero_scale":            "PCi_ja matrix is singular",
	"celpole_latp_tol_out":      "No valid solution for latp",
	"celpole_slz_tol_out":       "is required for these values of phip, phi0, and theta0",
	"celpole_z0_lat30":          "lat0 == 0 is required for |phip - phi0| = 90 and theta0 == 0",
	"sip_zero_coeff_deg10":      "Invalid powers (5, 5) for SIP on axis 1",
	"sip_ap_beyond_order":       "Degree must not exceed 9, got -1",
	"sip_inverse_coeff_deg10":   "Invalid powers (9, 1) for SIP on axis 2",
	"cunit_unknown":             "Mismatched units type 'plane angle'",
	"cunit_upper":               "Invalid symbol in EXPON context",
	"lonpole_car_dec30_lp180":   "No valid solution for latp",
	"lonpole_glon_car_bm30_lp0": "Ill-conditioned coordinate transformation parameter",
	"lonpole_glon_ait_bm20_lp0": "Ill-conditioned coordinate transformation parameter",
	"ncp_equator":               "NCP blows up on the equator",
	"proj_mismatch":             "Inconsistent projection types",
	"pvlon_ait_theta0_95":       "theta0 > 90",
	"pvlon_axis2_pv5":           "Unrecognized coordinate transformation parameter",
	"pvlon_car_pv1_5":           "Unrecognized coordinate transformation parameter",
	"pvlon_car_theta0_90p5":     "theta0 > 90",
	"pvlon_car_theta0_tol_out":  "theta0 > 90",
	"tpv_pv1_40":                "Invalid parameter number (40) for TPV",
	"tpv_pv2_40":                "Invalid parameter number (40) for TPV",
	"two_lat_axes":              "Inconsistent projection types",
	"two_lon_axes":              "Inconsistent projection types",
}

// hdrRequiredStricter: hdrcases_stricter/<name>.case
// (TestHdrStricterThanWcslib), with the Go error text each one expects.
var hdrRequiredStricter = map[string]string{
	"dss_amd_leading_zero":      "AMDX07 is not a DSS coefficient keyword",
	"dss_naxis3":                "DSS plate solution on a 3-axis WCS",
	"dss_pltdecsn_unpadded":     "PLTDECSN = '-' without a trailing blank",
	"couple_lat_row":            "couples pixel axis 3 into celestial world axis 2",
	"couple_lng_row":            "couples pixel axis 3 into celestial world axis 1",
	"legacy_pc_conflict":        "disagree",
	"legacy_projp_conflict":     "disagree",
	"ps_numeric":                "PS1_0 has a non-string value",
	"singular_cd":               "linear transform is singular",
	"sip_a_without_b":           "B_ORDER is missing",
	"sip_a_without_b_nosuffix":  "B_ORDER is missing",
	"sip_coeffs_without_orders": "A_ORDER is missing",
	"sip_index_leading_zero":    "A_01_1 is not a SIP coefficient keyword",
	"sip_nosuffix_cube_axes23":  "SIP is defined on pixel axes 1,2 only",
	"tpv_cunit2_arcsec":         "TPV coefficients are defined in degrees",
	"tnx_sip_nosuffix":          "TNX header carries SIP keywords",
	"sip_beyond_order_nonzero":  "A_3_0 = 1e-08 has degree 3 > A_ORDER = 2",
	"sip_order10":               "A_ORDER = 10 is outside 0..9",
	"sip_neg_order":             "A_ORDER = -1 is outside 0..9",
	"sip_order_float":           "A_ORDER = 2 is not an integer",
	"sip_coef_string":           "A_0_2 has a non-numeric value",
	"sip_mixed_suffix":          "inconsistent distortion codes",
	"sip_alt_suffix":            "SIP (A_ORDER/B_ORDER have no alternate-description form",
	"sip_unknown_dist_code":     "unsupported distortion \"XYZ\" in CTYPE",
	"tpu_alias":                 "unknown projection code: \"TPU\"",
	"tnx_sip_suffix":            "TNX header carries SIP keywords",
}

// hdrRequiredLenient: hdrcases_lenient/<name>.case (TestHdrLenientThanWcslib),
// with the wcslib error text each one expects.
var hdrRequiredLenient = map[string]string{
	"dss_plateid_last": "No distortion type, yet DQia keyvalues are present",
	"sip_ap_one_axis":  "Degree must not exceed 9, got -1",
	"sip_no_apbp":      "Degree must not exceed 9, got -1",
}

// hdrRequiredDiffer: hdrcases_differ/<name>.case
// (TestHdrDeliberateDifferences), with what each one pins.
var hdrRequiredDiffer = map[string]string{
	"real_atlas_k2_tpv_sip": "real ATLAS RA---TPV + redundant SIP (PR #4): SIP ignored as astropy; wcslib applies both, 24-26 arcsec off",
	"real_ptf_tpv_plus_sip": "real PTF TPV and SIP solutions merged, RA---TPV: SIP ignored as astropy; wcslib applies both",
	"tpv_sip_nosuffix":      "RA---TPV + a non-redundant SIP set, no -SIP: still ignored (astropy); wcslib applies both",
}

// hdrGeneratedPrefixes are the case-name prefixes written by generator
// scripts (gen_prj_cases.py, gen_sample_cases.py, gen_highorder_cases.py).
var hdrGeneratedPrefixes = []string{"prj_", "prjswap_", "pvfid_", "sample_", "hiorder_"}

func TestHdrRequiredCases(t *testing.T) {
	for _, set := range []struct {
		dir    string
		cases  map[string]string
		golden bool
	}{
		{hdrCaseDir, hdrRequiredCases, true},
		{hdrRejectDir, hdrRequiredReject, false},
		{hdrStricterDir, hdrRequiredStricter, false},
		{hdrLenientDir, hdrRequiredLenient, false},
		{hdrDifferDir, hdrRequiredDiffer, false},
	} {
		files, err := filepath.Glob(filepath.Join(set.dir, "*.case"))
		if err != nil {
			t.Fatal(err)
		}
		onDisk := map[string]bool{}
		for _, f := range files {
			onDisk[strings.TrimSuffix(filepath.Base(f), ".case")] = true
		}
		var names []string
		for n := range set.cases {
			names = append(names, n)
		}
		sort.Strings(names)
		for _, n := range names {
			if !onDisk[n] {
				t.Errorf("required case %s/%s.case (%s) is missing", set.dir, n, set.cases[n])
				continue
			}
			if set.golden {
				if !testutil.GoldenExists(filepath.Join(hdrGoldenDir, n+".json")) {
					t.Errorf("required case %s has no golden %s/%s.json(.gz)", n, hdrGoldenDir, n)
				}
			}
		}
	next:
		for n := range onDisk {
			if _, ok := set.cases[n]; ok {
				continue
			}
			for _, p := range hdrGeneratedPrefixes {
				if set.golden && strings.HasPrefix(n, p) {
					continue next
				}
			}
			t.Errorf("%s/%s.case is not listed in hdr_required_test.go: add it with the rule it pins", set.dir, n)
		}
	}
}
