//go:build interop

package healpix

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestHealpixGoldensReproduce: healpix/testdata/gen_golden.py (healpy
// 1.20.1, mpmath 1.4.1, and the healpix_cxx 8.3.0 harness hpref)
// regenerates every file of golden/ byte for byte, and no other file.
func TestHealpixGoldensReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	hcx := testutil.ExternalSource(t, testutil.EnvHealpixCxx, "healpix_base.cc", "healpix_cxx 8.3.0 source")
	testutil.MakeVerify(t, "testdata/hpref", []string{"hpref"}, nil, "HEALPIX_CXX="+hcx)
	py := testutil.PythonWith(t, "healpy", "mpmath", "numpy")
	out := t.TempDir()
	testutil.RunIn(t, "testdata", py, "gen_golden.py", out)
	testutil.SameFiles(t, out, filepath.Join("testdata", "golden"), nil)
}

// TestHealpixMapSamplesReproduce: testdata/maps/gen_map.py (healpy 1.20.1
// write_map) regenerates the RING and NESTED map samples and the healpy
// truth JSON byte for byte in a scratch copy with those files deleted
// first.
func TestHealpixMapSamplesReproduce(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	py := testutil.PythonWith(t, "healpy", "numpy")
	src := filepath.Join("testdata", "maps")
	dir := t.TempDir()
	testutil.CopyTree(t, src, dir)
	gen := []string{"healpy_map_nside16.json", "healpy_map_nside16_nest.fits", "healpy_map_nside16_ring.fits"}
	for _, n := range gen {
		if err := os.Remove(filepath.Join(dir, n)); err != nil {
			t.Fatal(err)
		}
	}
	testutil.RunIn(t, dir, py, "gen_map.py")
	testutil.SameFiles(t, dir, src, nil)
}
