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

// TestCfmovnamReproducesGoldens rebuilds dupext.fits, names.fits and their
// goldens with the committed cfitsio harnesses and requires byte-identical
// output (tier 2: a failure when they are not built).
func TestCfmovnamReproducesGoldens(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	for _, c := range []struct{ golden, maker string }{
		{"dupext.json", "mkdupext"},
		{"names.json", "mknames"},
	} {
		t.Run(c.golden, func(t *testing.T) {
			g, names, src := loadCfmovnamGolden(t, c.golden)
			var bins []string
			for _, b := range []string{c.maker, "cfmovnam"} {
				bin := testutil.TestdataPath(t, "cfmovnam", b)
				if st, err := os.Stat(bin); err != nil || st.Mode()&0o111 == 0 {
					t.Fatalf("cfitsio harness %s not built; build it with: make -C testdata/cfmovnam (needs libcfitsio)", bin)
				}
				bins = append(bins, bin)
			}
			out := filepath.Join(t.TempDir(), g.File)
			if b, err := exec.Command(bins[0], out).CombinedOutput(); err != nil {
				t.Fatalf("%s: %v: %s", c.maker, err, b)
			}
			if !bytes.Equal(mustRead(t, out), mustRead(t, src)) {
				t.Errorf("%s output differs from the committed %s", c.maker, src)
			}
			js, err := exec.Command(bins[1], append([]string{src}, names...)...).Output()
			if err != nil {
				t.Fatalf("cfmovnam: %v", err)
			}
			gp := testutil.TestdataPath(t, "cfmovnam", "golden", c.golden)
			if !bytes.Equal(js, mustRead(t, gp)) {
				t.Errorf("cfmovnam output differs from the committed golden:\n%s\nwant\n%s", js, mustRead(t, gp))
			}
		})
	}
}
