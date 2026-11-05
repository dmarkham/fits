//go:build interop

package fits

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCompressedGoldensReproduce (docs/VALIDATION.md rule 2) regenerates
// the compressed read corpus in a scratch directory from its committed
// generators (testdata/gen_fpack.sh with gen_fpack.c, libcfitsio,
// fpack/funpack 1.7.0 and astropy's real files; testdata/
// generate_compressed.py with astropy 8.0.2), and then dumps every
// regenerated sample with the cfref harness (testdata/cfref/cfref.c,
// libcfitsio 4.7.0). It requires, byte for byte:
//
//   - the twin map testdata/comp_fzun_same.txt (which funpack twins were
//     byte-identical to an earlier twin and are not committed);
//   - every committed cfitsio golden of the corpus:
//     testdata/cfref/golden/comp_*.json (the samples and sources),
//     testdata/cfref/golden/twins/*.json (the funpack/astropy twins, which
//     gen_twins.sh writes and are not committed) and
//     testdata/comp_bad/golden/*.json (the refused layouts, which
//     testdata/cfref's `make verify` does not cover);
//   - no committed golden without a regenerated sample;
//   - every pool file (testdata/cfref/golden/pool, the mid-size arrays the
//     cfref/golden goldens reference) the regenerated corpus writes.
//
// TestCompressedCorpusReproduce (reproduce_interop_test.go) checks the samples
// themselves; together they reproduce everything the read corpus tests
// use. Tier 2 (-tags interop): a missing tool fails the test.
func TestCompressedGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	for _, tool := range []string{"cc", "fpack", "funpack", "sh", "make"} {
		testutil.RequireTool(t, tool)
	}
	py := testutil.PythonWith(t, "astropy", "numpy")
	if out, err := exec.Command("make", "-C", filepath.Join("testdata", "cfref"), "cfref").CombinedOutput(); err != nil {
		t.Fatalf("cannot build testdata/cfref/cfref (needs cc and libcfitsio): %v\n%s", err, out)
	}
	cfref, err := filepath.Abs(filepath.Join("testdata", "cfref", "cfref"))
	if err != nil {
		t.Fatal(err)
	}

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
	testutil.SameFiles(t, dir, "testdata", []string{"comp_fzun_same.txt"})

	// Every golden is made with "dump -pool" (mid-size arrays and header
	// cards stored once in cfref/golden/pool), comp_bad's too.
	pool := filepath.Join(t.TempDir(), "pool")
	if err := os.MkdirAll(pool, 0o755); err != nil {
		t.Fatal(err)
	}
	dump := func(sample string, pooled bool) []byte {
		t.Helper()
		args := []string{"dump", filepath.Base(sample)}
		if pooled {
			args = []string{"dump", "-pool", pool, filepath.Base(sample)}
		}
		c := exec.Command(cfref, args...)
		c.Dir = filepath.Dir(sample)
		var out, stderr bytes.Buffer
		c.Stdout, c.Stderr = &out, &stderr
		if err := c.Run(); err != nil {
			t.Fatalf("cfref dump %s: %v\n%s", sample, err, stderr.String())
		}
		return out.Bytes()
	}
	checked := 0
	for _, set := range []struct{ samples, golden, suffix string }{
		{"comp_*.fits", filepath.Join("cfref", "golden"), ".json"},
		// the funpack/astropy twins gen_twins.sh wrote (not committed)
		{filepath.Join("cfref", "twins.out", "comp_fzun_*.fits"), filepath.Join("cfref", "golden", "twins"), ".json"},
		{filepath.Join("comp_bad", "*.fits"), filepath.Join("comp_bad", "golden"), ".json"},
	} {
		regen, _ := filepath.Glob(filepath.Join(dir, set.samples))
		have := map[string]bool{}
		for _, s := range regen {
			name := filepath.Base(s) + set.suffix
			have[name] = true
			want, err := os.ReadFile(filepath.Join("testdata", set.golden, name))
			if err != nil {
				t.Errorf("committed golden missing for regenerated %s: %v", filepath.Base(s), err)
				continue
			}
			if got := dump(s, true); !bytes.Equal(got, want) {
				t.Errorf("%s: cfref dump of the regenerated sample differs from the committed golden", filepath.Join(set.golden, name))
			}
			checked++
		}
		committed, _ := filepath.Glob(filepath.Join("testdata", set.golden, "comp_*"+set.suffix))
		if set.golden != filepath.Join("cfref", "golden") {
			committed, _ = filepath.Glob(filepath.Join("testdata", set.golden, "*"+set.suffix))
		}
		for _, g := range committed {
			if !have[filepath.Base(g)] {
				t.Errorf("%s: committed golden without a regenerated sample", g)
			}
		}
	}
	// Every pool file the regenerated corpus wrote is committed byte for
	// byte (the goldens above reference only these).
	regenPool, _ := filepath.Glob(filepath.Join(pool, "*.json"))
	for _, p := range regenPool {
		got, _ := os.ReadFile(p)
		want, err := testutil.ReadGolden(filepath.Join("testdata", "cfref", "golden", "pool", filepath.Base(p)))
		if err != nil || !bytes.Equal(got, want) {
			t.Errorf("pool file %s: regenerated bytes differ from the committed file (%v)", filepath.Base(p), err)
		}
	}
	if len(regenPool) < 50 {
		t.Errorf("only %d pool files regenerated", len(regenPool))
	}
	if checked < 200 {
		t.Fatalf("only %d goldens checked", checked)
	}
	t.Logf("%d compressed-corpus goldens reproduce byte for byte (%s)", checked, strings.Join([]string{"cfref/golden/comp_*", "comp_bad/golden"}, ", "))
}
