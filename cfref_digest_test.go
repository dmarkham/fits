package fits

import (
	"math"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCfrefDigestCatchesOneChange: a golden array stored as a digest
// (testdata/cfref/cfref.c jdigest) still proves every element. The
// library's decode of real M13 (90000 pixels) and of the 10000-row iter_a
// table equals the digest; changing any single element of the result (one
// float32 ULP at an index the sample does not cover, a sampled element, a
// null flag) makes the comparison fail, and the failure names what
// differs.
func TestCfrefDigestCatchesOneChange(t *testing.T) {
	const name = "comp_fz_real_m13_hcomp.fits"
	g := testutil.LoadCfrefGolden(t, name)
	im := g.HDUs[1].Image
	if im == nil || !im.Phys.Digested() || !im.I64.Digested() || im.Phys.Len() != 90000 {
		t.Fatalf("%s: expected digested 90000-pixel phys and i64 arrays", name)
	}
	got := ourPhysView(t, name, 1)
	if d := im.Phys.Diff(got); d != "" {
		t.Fatalf("library decode differs from the cfitsio digest: %s", d)
	}
	sampled := im.Phys.Sampled()
	unsampled := -1
	for k := 45000; k < len(got); k++ {
		if _, ok := sampled[k]; !ok && got[k].Valid {
			unsampled = k
			break
		}
	}
	if unsampled < 0 {
		t.Fatal("no unsampled valid pixel")
	}
	mutate := func(k int, f func(*testutil.NullFloat)) []testutil.NullFloat {
		m := append([]testutil.NullFloat(nil), got...)
		f(&m[k])
		return m
	}
	ulp := func(p *testutil.NullFloat) { p.V = float64(math.Nextafter32(float32(p.V), float32(math.Inf(1)))) }
	cases := []struct {
		name string
		v    []testutil.NullFloat
		want []string
	}{
		{"one ULP at an unsampled pixel", mutate(unsampled, ulp), []string{"sha256", "unsampled index"}},
		{"sampled pixel", mutate(0, ulp), []string{"[0] = ", "sha256"}},
		{"valid pixel made null", mutate(unsampled, func(p *testutil.NullFloat) { *p = testutil.NullFloat{} }), []string{"null runs", "sha256"}},
		{"sign flipped", mutate(unsampled, func(p *testutil.NullFloat) { p.V = -p.V }), []string{"sha256"}},
		{"one pixel missing", got[:len(got)-1], []string{"89999 values, cfitsio 90000"}},
	}
	for _, c := range cases {
		d := im.Phys.Diff(c.v)
		for _, w := range c.want {
			if !strings.Contains(d, w) {
				t.Errorf("%s: diff %q does not contain %q", c.name, d, w)
			}
		}
	}

	// A digested table column (iter_a: 10000 rows), checked through
	// Resolved: one changed cell fails.
	ga := testutil.LoadCfrefGolden(t, "iter_a.fit")
	col := ga.HDUs[1].Table.Columns[0]
	if !col.Digested() {
		t.Fatal("iter_a column 1 is not digested")
	}
	f, err := Open(testutil.TestdataPath(t, "iter_a.fit"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	hdu, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	bt := hdu.(*BinaryTableHDU)
	vals, ok, err := ReadVectorColumnMasked[int64](bt, 1)
	if err != nil {
		t.Fatal(err)
	}
	i64 := testutil.NullInt64s(vals, ok)
	if d := col.Flat().I64.Diff(i64); d != "" {
		t.Fatalf("iter_a Counts differs from the cfitsio digest: %s", d)
	}
	i64[5000].V++
	if d := col.Flat().I64.Diff(i64); !strings.Contains(d, "sha256") {
		t.Errorf("changed row 5000 not detected: %q", d)
	}
}
