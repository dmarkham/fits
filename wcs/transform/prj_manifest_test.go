package transform

import (
	"bytes"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// The projection proof's committed truth, by name. Globbing the testdata
// directories alone would let a deleted case or golden go unnoticed (the
// tests just run over fewer files), so the expected set is written down
// here and TestProjectionProofManifest fails - without needing a C
// compiler, Python or wcslib - on any file that is missing or that is
// present but not listed. Adding a case means adding it here too.
//
// TestPrjCaseGeneratorReproduce and TestSampleGeneratorsReproduce
// (interop) close the loop on the other side: re-running the generators
// must write exactly these case files, byte for byte.

// wrefCaseManifest is every bare-projection golden (golden/<name>.json, from
// wref.c) the projection proof depends on.
var wrefCaseManifest = []string{
	"air_default", "air_thetab45", "ait_default", "arc_default",
	"azp_default", "azp_mu1p5_gammam25", "azp_mu2", "azp_mu2_gamma30",
	"azp_mum3_gamma20", "azp_mu1_gamma45", "bon_t0", "bon_t30", "bon_tm45", "car_default",
	"cea_default", "cea_lambda0p5", "cod_ta45", "cod_ta45_eta10",
	"cod_tam30_eta15", "coe_ta45", "coe_ta45_eta10", "coe_tam30_eta15",
	"coo_ta45", "coo_ta45_eta10", "coo_tam30_eta15", "cop_ta45",
	"cop_ta45_eta10", "cop_tam30_eta15", "csc_default", "cyp_default",
	"cyp_gall", "cyp_mu0_lambda1", "cyp_mu2_lambda0p5", "hpx_default",
	"hpx_h3_k2", "hpx_h4_k3", "hpx_h5_k5", "hpx_h6_k1", "mer_default",
	"mol_default", "par_default", "pco_default", "qsc_default", "sfl_default",
	"sin_default", "sin_ncp60", "sin_slant", "stg_default", "szp_default",
	"szp_mu0_phic30_thetac70", "szp_mu1p5_phic90_thetac45",
	"szp_mu2_phic180_thetac60", "tan_default", "tsc_default", "xph_default",
	"zea_default", "zpn_cubic", "zpn_cubic_inflect", "zpn_linear",
	"zpn_pv29", "zpn_quadratic", "zpn_quintic",
}

// prjHdrCaseManifest is every synthetic full-pipeline case gen_prj_cases.py
// writes (hdrcases/prj_*.case, pvfid_*.case and prjswap_*.case).
var prjHdrCaseManifest = []string{
	"prj_air_c0", "prj_air_c30_tb45", "prj_ait_c0", "prj_ait_cm40",
	"prj_arc_c0", "prj_arc_cm35", "prj_azp_c0", "prj_azp_c45_mu2_g30",
	"prj_bon_c0_t30", "prj_bon_c50_tm45", "prj_car_c0", "prj_car_c30",
	"prj_cea_c0", "prj_cea_c30_l0p5", "prj_cod_c0_ta45",
	"prj_cod_c60_ta45_eta10", "prj_coe_c0_ta45", "prj_coe_cm20_tam30_eta15",
	"prj_coo_c0_ta45", "prj_coo_cm20_tam30_eta15", "prj_cop_c0_ta45",
	"prj_cop_cm20_tam30_eta15", "prj_csc_c0", "prj_csc_c30", "prj_cyp_c0",
	"prj_cyp_c30_gall", "prj_hpx_c0", "prj_hpx_c0_h3k2", "prj_hpx_c0_h6k1",
	"prj_hpx_c30", "prj_mer_c0", "prj_mer_cm40", "prj_mol_c0", "prj_mol_c30",
	"prj_par_c0", "prj_par_cm40", "prj_pco_c0", "prj_pco_c30", "prj_qsc_c0",
	"prj_qsc_c30", "prj_sfl_c0", "prj_sfl_c30", "prj_sin_c0",
	"prj_sin_c50_slant", "prj_stg_c0", "prj_stg_cm35", "prj_szp_c0",
	"prj_szp_cm40_mu1p5", "prj_tan_c0", "prj_tan_c60", "prj_tsc_c0",
	"prj_tsc_c30", "prj_xph_c0", "prj_xph_cm50", "prj_zea_c0", "prj_zea_cm35",
	"prj_zpn_c0_cubic", "prj_zpn_c0_deg12", "prj_zpn_c70_quintic", "pvfid_car_pv1_off",
	"pvfid_car_pv1_t20", "pvfid_tan_pv1_off",
	"prjswap_azp_mu2_g30", "prjswap_cea_l0p5", "prjswap_cop_tam30_eta15",
	"prjswap_zpn_cubic",
}

// sampleHdrCaseManifest is every real-header case gen_sample_cases.py
// writes (hdrcases/sample_*.case).
var sampleHdrCaseManifest = []string{
	"sample_1904-66_air", "sample_1904-66_ait", "sample_1904-66_arc",
	"sample_1904-66_azp", "sample_1904-66_bon", "sample_1904-66_car",
	"sample_1904-66_cea", "sample_1904-66_cod", "sample_1904-66_coe",
	"sample_1904-66_coo", "sample_1904-66_cop", "sample_1904-66_csc",
	"sample_1904-66_cyp", "sample_1904-66_hpx", "sample_1904-66_mer",
	"sample_1904-66_mol", "sample_1904-66_ncp", "sample_1904-66_par",
	"sample_1904-66_pco", "sample_1904-66_qsc", "sample_1904-66_sfl",
	"sample_1904-66_sin", "sample_1904-66_stg", "sample_1904-66_szp",
	"sample_1904-66_tan", "sample_1904-66_tsc", "sample_1904-66_zea",
	"sample_1904-66_zpn", "sample_hpxcvt_hpx_equ", "sample_hpxcvt_hpx_gal",
	"sample_hpxcvt_hpx_gal_q1", "sample_hpxcvt_xph_n_gal",
	"sample_hpxcvt_xph_s_gal", "sample_msx_glon_car_lonpole0",
	"sample_rosat_glon_ait_lonpole0",
}

// projectionSampleFiles are the committed sample files under
// testdata/wref/samples the projection proof reads (TestHPXcvtPixelCentres).
var projectionSampleFiles = []string{"hpxcvt_hpx_gal_n8.fits", "hpxcvt_hpx_gal_n8.healpy.json"}

const (
	wrefDir       = "testdata/wref"
	wrefSampleDir = "testdata/wref/samples"
)

// globNames returns the base names (without ext) matching dir/pattern.
func globNames(t *testing.T, dir, pattern, ext string) []string {
	t.Helper()
	files, err := filepath.Glob(filepath.Join(dir, pattern))
	if err != nil {
		t.Fatal(err)
	}
	var names []string
	for _, f := range files {
		names = append(names, strings.TrimSuffix(filepath.Base(f), ext))
	}
	sort.Strings(names)
	return names
}

// setDiff returns the names in a but not in b.
func setDiff(a, b []string) []string {
	in := map[string]bool{}
	for _, n := range b {
		in[n] = true
	}
	var out []string
	for _, n := range a {
		if !in[n] {
			out = append(out, n)
		}
	}
	return out
}

// TestProjectionProofManifest requires every bare-projection golden, every
// projection header case with its wcslib golden, and every sample file the
// projection proof names, and refuses unlisted ones; and it requires a real
// header (Calabretta's 1904-66 example map, or an HPXcvt image for XPH)
// for every projection code we claim.
func TestProjectionProofManifest(t *testing.T) {
	var problems []string
	check := func(what string, want, have []string) {
		for _, n := range setDiff(want, have) {
			problems = append(problems, what+": "+n+" is MISSING")
		}
		for _, n := range setDiff(have, want) {
			problems = append(problems, what+": "+n+" is present but not in the manifest (add it to prj_manifest_test.go)")
		}
	}
	wref := setDiff(globNames(t, wrefGoldenDir, "*.json", ".json"), []string{wrefPrjsetStatus})
	check("bare-projection golden "+wrefGoldenDir, wrefCaseManifest, wref)
	if _, err := os.Stat(filepath.Join(wrefGoldenDir, wrefPrjsetStatus+".json")); err != nil {
		problems = append(problems, "prjset status golden is MISSING: "+err.Error())
	}

	prj := append(globNames(t, hdrCaseDir, "prj_*.case", ".case"), globNames(t, hdrCaseDir, "pvfid_*.case", ".case")...)
	prj = append(prj, globNames(t, hdrCaseDir, "prjswap_*.case", ".case")...)
	sort.Strings(prj)
	check("hdr case "+hdrCaseDir, prjHdrCaseManifest, prj)
	check("hdr case "+hdrCaseDir, sampleHdrCaseManifest, globNames(t, hdrCaseDir, "sample_*.case", ".case"))
	for _, n := range append(append([]string{}, prjHdrCaseManifest...), sampleHdrCaseManifest...) {
		if p := filepath.Join(hdrGoldenDir, n+".json"); !testutil.GoldenExists(p) {
			problems = append(problems, "wcslib golden for hdr case "+n+" is MISSING: "+p+"(.gz)")
		}
	}
	for _, f := range projectionSampleFiles {
		if _, err := os.Stat(filepath.Join(wrefSampleDir, f)); err != nil {
			problems = append(problems, "sample "+f+" is MISSING: "+err.Error())
		}
	}

	listed := map[string]bool{}
	for _, n := range sampleHdrCaseManifest {
		listed[n] = true
	}
	var codes []string
	for code := range allProjectionCodes {
		codes = append(codes, code)
	}
	sort.Strings(codes)
	for _, code := range codes {
		want := []string{"sample_1904-66_" + strings.ToLower(code)}
		if code == "XPH" {
			// Calabretta's Paper II set predates XPH; HPXcvt writes real XPH images.
			want = []string{"sample_hpxcvt_xph_n_gal", "sample_hpxcvt_xph_s_gal"}
		}
		for _, n := range want {
			if !listed[n] {
				problems = append(problems, "projection "+code+": real-header case "+n+" is not in the manifest")
			}
		}
	}
	if len(problems) > 0 {
		t.Fatalf("projection proof manifest: %d problem(s):\n  %s", len(problems), strings.Join(problems, "\n  "))
	}
}

// diffGenerated compares a generator's output directory with the committed
// one: same file names (matching pattern) and the same bytes.
func diffGenerated(t *testing.T, genDir, committedDir, pattern string) {
	t.Helper()
	got := globNames(t, genDir, pattern, "")
	have := globNames(t, committedDir, pattern, "")
	for _, n := range setDiff(got, have) {
		t.Errorf("%s: generator writes %s but it is not committed", committedDir, n)
	}
	for _, n := range setDiff(have, got) {
		t.Errorf("%s: %s is committed but the generator no longer writes it", committedDir, n)
	}
	for _, n := range got {
		g, err := os.ReadFile(filepath.Join(genDir, n))
		if err != nil {
			t.Fatal(err)
		}
		c, err := os.ReadFile(filepath.Join(committedDir, n))
		if err != nil {
			continue // reported above
		}
		if !bytes.Equal(g, c) {
			t.Errorf("%s/%s: regenerated file differs from the committed one (re-run the generator and review the diff)", committedDir, n)
		}
	}
}
