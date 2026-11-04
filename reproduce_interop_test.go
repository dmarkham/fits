//go:build interop

package fits

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// The *Reproduce* tests of docs/VALIDATION.md rule 2 for the golden sets
// of the root package: each re-runs the committed generator and requires
// the committed bytes back (tier 2, -tags interop: a missing reference
// tool fails the test).

// TestCfnullGoldensReproduce: testdata/imaging/cfnull (libcfitsio
// fits_read_pixnull of every tile-compressed imaging sample).
func TestCfnullGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	testutil.MakeVerify(t, "testdata/imaging/cfnull", []string{"cfnull"}, []string{"verify"})
}

// TestSirilrefGoldensReproduce: testdata/imaging/sirilref (Siril 1.4.0's
// FITS-to-float functions compiled verbatim), and the check that those
// functions are still byte-identical to the Siril source.
func TestSirilrefGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	testutil.MakeVerify(t, "testdata/imaging/sirilref", []string{"sirilref"}, []string{"verify"})
	src := testutil.ExternalSource(t, testutil.EnvSirilSrc, filepath.Join("src", "io", "image_format_fits.c"), "Siril 1.4.0 source tree")
	testutil.MakeVerify(t, "testdata/imaging/sirilref", nil, []string{"verbatim"}, "SIRIL_SRC="+src)
}

// TestImagesCorpusReproduce: testdata/images (gen.py: cfitsio-written and
// astropy-written samples, astropy's real files, framing and refused
// cases, every golden) regenerates byte for byte.
func TestImagesCorpusReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "astropy", "numpy")
	testutil.MakeVerify(t, "testdata/images", []string{"imgref", "cfref"}, []string{"verify"}, "PYTHON="+py)
}

// TestImagingExpectedReproduce: testdata/imaging/gen_expected.py
// (readfloat32_expected.json, stdimage_expected.json from astropy and the
// cfnull goldens) regenerates byte for byte in a copy of the directory.
func TestImagingExpectedReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "astropy", "numpy")
	dir := t.TempDir()
	testutil.CopyTree(t, filepath.Join("testdata", "imaging"), dir)
	testutil.RunIn(t, dir, py, "gen_expected.py")
	testutil.SameFiles(t, dir, filepath.Join("testdata", "imaging"),
		[]string{"readfloat32_expected.json", "stdimage_expected.json"})
}

// TestCompressedCorpusReproduce: the compressed read corpus regenerates
// byte for byte from its generators in a scratch copy: testdata/gen_fpack.sh
// (gen_fpack.c + libcfitsio, fpack/funpack 1.7.0, astropy's real files;
// every comp_fz_*, comp_fzsrc_* and comp_bad/* file; the twins are not
// committed, their goldens are checked by TestCompressedGoldensReproduce) and
// testdata/generate_compressed.py (astropy; every comp_<name> and
// comp_raw_<name> fixture). A hand-edited sample fails here.
func TestCompressedCorpusReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	for _, tool := range []string{"cc", "fpack", "funpack", "sh"} {
		testutil.RequireTool(t, tool)
	}
	py := testutil.PythonWith(t, "astropy", "numpy")
	dir := t.TempDir()
	if err := os.MkdirAll(filepath.Join(dir, "comp_bad"), 0o755); err != nil {
		t.Fatal(err)
	}
	for _, f := range []string{"gen_fpack.sh", "gen_fpack.c", "gen_twins.sh", "generate_compressed.py"} {
		b, err := os.ReadFile(filepath.Join("testdata", f))
		if err != nil {
			t.Fatalf("committed generator missing: %v", err)
		}
		if err := os.WriteFile(filepath.Join(dir, f), b, 0o644); err != nil {
			t.Fatal(err)
		}
	}
	cmd := exec.Command("sh", "gen_fpack.sh")
	cmd.Dir = dir
	cmd.Env = append(os.Environ(), "TMPDIR="+t.TempDir(), "FITS_ASTROPY_PYTHON="+py)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("gen_fpack.sh: %v\n%s", err, out)
	}
	testutil.RunIn(t, dir, py, "generate_compressed.py")

	var want, got []string
	for _, pat := range []string{"comp_*.fits", filepath.Join("comp_bad", "*.fits")} {
		w, _ := filepath.Glob(filepath.Join("testdata", pat))
		g, _ := filepath.Glob(filepath.Join(dir, pat))
		for _, p := range w {
			rel, _ := filepath.Rel("testdata", p)
			want = append(want, rel)
		}
		for _, p := range g {
			rel, _ := filepath.Rel(dir, p)
			got = append(got, rel)
		}
	}
	if len(want) < 200 {
		t.Fatalf("only %d committed compressed samples", len(want))
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Errorf("regenerated file set differs from the committed one:\n got  %v\n want %v", got, want)
	}
	testutil.SameFiles(t, dir, "testdata", want)
}

// TestTableSamplesReproduce: testdata/tables `make verify` regenerates the
// table SAMPLES (gen_tables.c with libcfitsio, gen_gaia.py and
// gen_theap.py with astropy) and their cfitsio goldens into a scratch
// directory and byte-compares both with the committed files, so a
// hand-edited sample (e.g. one bit of gaia_dr3_sample.fits with its golden
// regenerated to match) fails here. TestTableGoldensReproduce only
// re-dumps the committed samples.
func TestTableSamplesReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "astropy", "numpy")
	testutil.MakeVerify(t, "testdata/tables", []string{"gen_tables.out", "cf_tfloat.out"}, []string{"verify"}, "PYTHON="+py)
}

// TestImagingSamplesReproduce: the imaging samples regenerate byte for
// byte in a scratch copy with every generated file deleted first:
// gen_samples.py (numpy, raw FITS bytes), make_compressed.sh (fpack 1.7.0
// -C) and zblank_variants.py. The only input kept is
// asi533_u16_crop.fits, the real camera crop (crop_camera.py needs the
// original download; see PROVENANCE.md). Missing, extra or changed files
// all fail.
func TestImagingSamplesReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "numpy")
	testutil.RequireTool(t, "fpack")
	testutil.RequireTool(t, "sh")
	src := filepath.Join("testdata", "imaging")
	dir := t.TempDir()
	testutil.CopyTree(t, src, dir)
	ents, err := os.ReadDir(dir)
	if err != nil {
		t.Fatal(err)
	}
	removed := 0
	for _, e := range ents {
		n := e.Name()
		if e.IsDir() || n == "asi533_u16_crop.fits" || !(strings.HasSuffix(n, ".fits") || strings.HasSuffix(n, ".fz")) {
			continue
		}
		if err := os.Remove(filepath.Join(dir, n)); err != nil {
			t.Fatal(err)
		}
		removed++
	}
	if removed < 50 {
		t.Fatalf("only %d generated imaging samples found", removed)
	}
	testutil.RunIn(t, dir, py, "gen_samples.py")
	cmd := exec.Command("sh", "make_compressed.sh")
	cmd.Dir = dir
	cmd.Env = append(os.Environ(), "PYTHON="+py)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("make_compressed.sh: %v\n%s", err, out)
	}
	testutil.SameFiles(t, dir, src, nil)
}
