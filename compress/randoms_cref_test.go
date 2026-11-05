package compress

import (
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestFitsRandomsMatchCfitsio compares the whole dither table with the
// one libcfitsio builds (fits_init_randoms -> fits_rand_value), dumped by
// `cref gen` into golden/fits_rand_value.json (each float32 widened to
// double, as an exact digest: sha256 over all N_RANDOM values plus a
// sample). Every one of the N_RANDOM float32 values must be bit-identical: the dither offsets of every
// SUBTRACTIVE_DITHER_1/2 pixel are read from this table, so a single
// wrong entry corrupts the pixels that use it.
//
// It also pins a few entries as literals (copied from the cfitsio dump)
// so a regenerated golden that drifted would be noticed in review. The
// way the table is walked is proven by writer parity goldens (see
// write_compressed_cfitsio_test.go): nextrand wrapping at N_RANDOM inside
// a tile by rice_f32_wrap_100x100 (one 10000-pixel tile, ZDITHER0=9000);
// the per-tile start row (row + ZDITHER0 - 1) % N_RANDOM reaching 0 by
// rice_f32_sd1_z10000; and iseed wrapping 9999 -> 0 inside a tile, in the
// no-null and the null branch of the quantizer, by
// rice_f32_wrapin_z10000 and rice_f32_wrapin_nan_z10000 (one
// 10000-pixel tile starting at row 9999).
func TestFitsRandomsMatchCfitsio(t *testing.T) {
	path := filepath.Join(cgoldenDir, "fits_rand_value.json")
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("committed golden missing (failure, not a skip): %v", err)
	}
	var g struct {
		NRandom int                                 `json:"n_random"`
		Values  testutil.Values[testutil.NullFloat] `json:"values"`
	}
	if err := json.Unmarshal(b, &g); err != nil {
		t.Fatalf("parse %s: %v", path, err)
	}
	if g.NRandom != NRandom || g.Values.Len() != NRandom {
		t.Fatalf("%s: %d values (n_random %d), want %d", path, g.Values.Len(), g.NRandom, NRandom)
	}
	got := FitsRandoms()
	if len(got) != NRandom {
		t.Fatalf("len(FitsRandoms) = %d, want %d", len(got), NRandom)
	}
	cand := make([]testutil.NullFloat, NRandom)
	for i, v := range got {
		cand[i] = testutil.NullFloat{V: float64(v), Valid: true}
	}
	// exact: any single differing entry fails (sha256), samples say where
	resolved := g.Values.Resolve(t, "FitsRandoms vs libcfitsio fits_rand_value", func() []testutil.NullFloat { return cand })
	want := make([]float32, NRandom)
	for i, v := range resolved {
		want[i] = float32(v.V)
	}

	// Literal spot checks: first, the start of the second "row", the
	// last entry, and the entry the wrap returns to.
	for _, c := range pinnedRandoms {
		if math.Float32bits(want[c.i]) != c.bits {
			t.Errorf("golden fits_rand_value[%d] = %#x, pinned %#x (golden regenerated with a different cfitsio?)", c.i, math.Float32bits(want[c.i]), c.bits)
		}
	}

}

// pinnedRandoms are fits_rand_value entries from libcfitsio 4.7.0,
// as raw float32 bits.
var pinnedRandoms = []struct {
	i    int
	bits uint32
}{
	{0, 0x37034e00},
	{1, 0x3e06b1d8},
	{2, 0x3f416f5a},
	{499, 0x3f27123b},
	{500, 0x3f1dda2e},
	{5000, 0x3e25c663},
	{9998, 0x3f310015},
	{9999, 0x3ef8d164},
}
