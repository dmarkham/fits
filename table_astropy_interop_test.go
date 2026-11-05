//go:build interop

package fits

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestAstropyGoldensReproduce re-runs astropy_dump.py --golden on every
// sample and requires the committed golden byte for byte (tier 2: needs
// astropy, a failure without it).
func TestAstropyGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.AstropyPython(t)
	script := testutil.TestdataPath(t, "tables", "astropy_dump.py")
	paths := []string{testutil.TestdataPath(t, "tables", "astropy_only", "astropy_special.fits")}
	for _, s := range tableSamples(t) {
		paths = append(paths, s.path)
	}
	gold, _ := filepath.Glob(testutil.TestdataPath(t, "tables", "golden_astropy", "*.json"))
	if len(gold) != len(paths) {
		t.Errorf("golden_astropy has %d files for %d samples", len(gold), len(paths))
	}
	for _, p := range paths {
		live, err := exec.Command(py, script, "--golden", p).Output()
		if err != nil {
			t.Fatalf("%s: astropy_dump.py --golden: %v", p, err)
		}
		gp := testutil.TestdataPath(t, "tables", "golden_astropy", filepath.Base(p)+".json")
		testutil.RequireFile(t, gp)
		want, err := os.ReadFile(gp)
		if err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(live, want) {
			t.Errorf("%s: live astropy dump differs from %s (regenerate with make -C testdata/tables astropy-goldens and review)", filepath.Base(p), filepath.Base(gp))
		}
	}
}
