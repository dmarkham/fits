//go:build interop

package stats

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"slices"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// EnvSirilSrc points at an unpacked Siril 1.4.0 source tree; it has no
// default.
const EnvSirilSrc = testutil.EnvSirilSrc

// TestGoldensReproduce rebuilds both truth generators from source and
// checks that they reproduce the committed golden files byte-for-byte,
// and that the Siril excerpts are verbatim copies of the Siril source.
// It needs a C/C++ compiler, make, a Python with numpy+astropy
// (FITS_ASTROPY_PYTHON) and the Siril source (FITS_SIRIL_SRC); without
// them it fails (tier 2, -tags interop).
func TestGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	makeBin := testutil.RequireTool(t, "make")
	testutil.RequireTool(t, "cc")
	testutil.RequireTool(t, "c++")
	py := testutil.AstropyPython(t)

	srefDir, err := filepath.Abs("testdata/sref")
	if err != nil {
		t.Fatal(err)
	}
	tmp := t.TempDir()
	bin := filepath.Join(tmp, "sref")
	run := func(cmd *exec.Cmd) {
		t.Helper()
		out, err := cmd.CombinedOutput()
		if err != nil {
			t.Fatalf("%v: %v\n%s", cmd.Args, err, out)
		}
	}
	run(exec.Command(makeBin, "-C", srefDir, "BIN="+bin, bin))

	t.Run("sref", func(t *testing.T) {
		out := filepath.Join(tmp, "sref-golden")
		if err := os.MkdirAll(out, 0o755); err != nil {
			t.Fatal(err)
		}
		run(exec.Command(bin, "gen", out))
		sameFile(t, filepath.Join(out, "siril_stats_golden.json"), sirilGoldenPath)
		onlyFiles(t, out, filepath.Dir(sirilGoldenPath))
	})

	t.Run("gen_golden.py", func(t *testing.T) {
		out := filepath.Join(tmp, "py-golden")
		cmd := exec.Command(py, "testdata/gen_golden.py", out)
		cmd.Env = append(os.Environ(), "SREF_BIN="+bin)
		run(cmd)
		sameFile(t, filepath.Join(out, "stats_golden.json"), statsGoldenPath)
		sameFile(t, filepath.Join(out, "mean_rounding.json"), meanRoundingGoldenPath)
		sameFile(t, filepath.Join(out, "histogram_ranges.json"), histogramRangesGoldenPath)
		sameFile(t, filepath.Join(out, "histogram_edges.json"), histogramEdgesGoldenPath)
		onlyFiles(t, out, filepath.Dir(statsGoldenPath))
	})

	t.Run("siril_excerpt_verbatim", func(t *testing.T) {
		src := testutil.ExternalSource(t, EnvSirilSrc, filepath.Join("src", "rt", "rt_algo.cc"), "Siril 1.4.0 source tree")
		run(exec.Command(py, filepath.Join(srefDir, "extract_siril.py"), "--check", src))
	})
}

// onlyFiles: the committed golden directory holds exactly the files the
// generator wrote to out. A stray committed file would be an unchecked
// golden, and a generated file with no committed copy an unproven one.
func onlyFiles(t *testing.T, out, committed string) {
	t.Helper()
	names := func(dir string) []string {
		ents, err := os.ReadDir(dir)
		if err != nil {
			t.Fatal(err)
		}
		var n []string
		for _, e := range ents {
			n = append(n, e.Name())
		}
		return n
	}
	if got, want := names(committed), names(out); !slices.Equal(got, want) {
		t.Errorf("%s holds %v; the generator writes %v", committed, got, want)
	}
}

func sameFile(t *testing.T, got, want string) {
	t.Helper()
	// testutil.ReadGolden: a golden committed as <name>.gz is compared
	// with the generator's <name>.gz after decompression.
	a, err := testutil.ReadGolden(got)
	if err != nil {
		t.Fatal(err)
	}
	b, err := testutil.ReadGolden(want)
	if err != nil {
		t.Fatalf("committed golden %s missing: %v", want, err)
	}
	if !bytes.Equal(a, b) {
		t.Errorf("regenerated %s differs from committed %s (%d vs %d bytes)", got, want, len(a), len(b))
	}
}
