//go:build interop

package stretch

import (
	"bytes"
	"os"
	"os/exec"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestSirilGoldensReproduce re-runs gen_siril_goldens.py (the verbatim
// Siril 1.4.0 oracle in testdata/sref, built here with make) into a scratch
// directory and requires every file it writes -- the inputs, the 18 preset
// goldens and every golden/sref case -- to equal the committed one byte
// for byte, and nothing else of its set to be committed. The generator
// itself refuses to write the preset goldens unless the oracle reproduces
// the original siril-cli digests. Needs cc, c++, glib2, OpenCV 5 and numpy
// (tier 2: a failure without them).
func TestSirilGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "numpy")
	if _, err := exec.LookPath("make"); err != nil {
		t.Fatalf("make not found: %v", err)
	}
	if out, err := exec.Command("make", "-C", "testdata/sref").CombinedOutput(); err != nil {
		t.Fatalf("cannot build the Siril oracle (needs cc, c++, glib2, opencv5): %v\n%s", err, out)
	}
	dir := t.TempDir()
	testutil.RunIn(t, "testdata", py, "gen_siril_goldens.py", "--out", dir)
	sameGoldenSet(t, dir, sirilOwned)
}

// TestCLAHECv2GoldensReproduce: gen_clahe_cv2.py (opencv-python-headless
// 5.0.0.93) rewrites its goldens and inputs byte for byte.
func TestCLAHECv2GoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "cv2", "numpy")
	dir := t.TempDir()
	testutil.RunIn(t, "testdata", py, "gen_clahe_cv2.py", "--out", dir)
	sameGoldenSet(t, dir, func(rel string) bool {
		return !sirilOwned(rel) && rel != "siril_cli_sha256.json" && rel != "siril_cli_args.json"
	})
}

// sirilCLIPath is the siril-cli the CLI checks run: $FITS_SIRIL_CLI (no
// default; testdata/build_siril_cli.sh builds one).
func sirilCLIPath(t *testing.T) string {
	t.Helper()
	p := os.Getenv(testutil.EnvSirilCLI)
	if p == "" {
		t.Fatalf("%s not set (a siril-cli 1.4.0 executable; testdata/build_siril_cli.sh builds one)", testutil.EnvSirilCLI)
	}
	if _, err := exec.LookPath(p); err != nil {
		t.Fatalf("%s=%s is not an executable: %v", testutil.EnvSirilCLI, p, err)
	}
	return p
}

// TestSirilCLIReproduce runs a real siril-cli 1.4.0 (gen_golden.py): the
// 18 preset outputs must have the digests of golden/siril_cli_sha256.json,
// every oracle golden that names a siril command must come out of
// siril-cli with the golden's digest (146 of them, the 7 inf_* cases
// included: siril-cli's FITS loader rescales float data whose maximum
// exceeds 10, so gen_golden.py loads ±7 markers and turns them into ±Inf
// with siril-cli's pixel math, checked bit for bit), and the
// accept/refuse/pixel table golden/siril_cli_args.json must be what
// siril-cli produces now.
func TestSirilCLIReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "numpy", "astropy")
	cli := sirilCLIPath(t)
	cmd := exec.Command(py, "gen_golden.py", cli)
	cmd.Dir = "testdata"
	cmd.Env = append(os.Environ(), "TMPDIR="+t.TempDir())
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("gen_golden.py %s: %v\n%s", cli, err, out)
	}
	if !bytes.Contains(out, []byte("146 oracle goldens with a siril command, all compared bit for bit (7 with +-Inf input)")) {
		t.Fatalf("gen_golden.py compared a different set of oracle goldens:\n%s", out)
	}
	if !bytes.Contains(out, []byte("siril-cli agrees with every committed golden.")) {
		t.Fatalf("gen_golden.py did not report agreement:\n%s", out)
	}
}

// TestSirilCLIDigestsFromHistory recomputes golden/siril_cli_sha256.json
// from the siril-cli pixel files in git history (siril_cli_digests.py
// --check). Needs git with the history and numpy.
func TestSirilCLIDigestsFromHistory(t *testing.T) {
	py := testutil.PythonWith(t, "numpy")
	if _, err := exec.LookPath("git"); err != nil {
		t.Fatalf("git not found: %v", err)
	}
	var d sirilCLIDigests
	readGolden(t, "siril_cli_sha256.json", &d)
	if err := exec.Command("git", "cat-file", "-e", d.Input.Blob).Run(); err != nil {
		t.Fatalf("git history with blob %s (commit c17580b) not available: %v", d.Input.Blob, err)
	}
	testutil.RunIn(t, "testdata", py, "siril_cli_digests.py", "--check")
}
