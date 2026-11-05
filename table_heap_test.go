package fits

import (
	"fmt"
	"os"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestVLAHeapBoundaries pins the heap limits of FITS 4.0 7.3.5 at their
// boundaries: the heap starts THEAP bytes into the data unit (default
// NAXIS1*NAXIS2), never inside the fixed rows or past the data unit, and
// ends at NAXIS1*NAXIS2 + PCOUNT; a descriptor whose payload runs even one
// byte past it would read the block padding or the next HDU as data
// (cfitsio 4.7.0 does exactly that: end_plus1 reads 1792, 2304). The
// samples are testdata/tables/heap (gen_heap.py). fitsverify's committed
// verdict (heap/fitsverify.txt) agrees wherever it checks the rule; it
// measures descriptors against PCOUNT rather than the heap after a THEAP
// gap and does not check a THEAP past the data unit, which are pinned as
// the two places where it does not.
func TestVLAHeapBoundaries(t *testing.T) {
	verdicts := map[string]string{}
	vpath := testutil.TestdataPath(t, "tables", "heap", "fitsverify.txt")
	testutil.RequireFile(t, vpath)
	b, err := os.ReadFile(vpath)
	if err != nil {
		t.Fatal(err)
	}
	for _, line := range strings.Split(strings.TrimSpace(string(b)), "\n") {
		ok := strings.HasPrefix(line, "verification OK: ")
		name := strings.Fields(strings.TrimPrefix(strings.TrimPrefix(line, "verification OK: "), "verification FAILED: "))[0]
		name = strings.TrimSuffix(name, ",")
		verdicts[name] = map[bool]string{true: "ok", false: "error"}[ok]
	}
	cases := []struct {
		name    string
		wantErr string  // "" = reads
		want    []int32 // the row's values when it reads
		fv      string  // fitsverify's verdict
	}{
		{"end_exact", "", []int32{7, 9}, "ok"},
		{"end_plus1", "lies outside the heap", nil, "error"},
		{"off_eq", "lies outside the heap", nil, "error"},
		{"theap_end", "", []int32{}, "ok"},
		{"theap_end_plus1", "THEAP = 17 outside", nil, "ok"},      // fitsverify does not check THEAP past the data unit
		{"theap_fixed_minus1", "THEAP = 7 outside", nil, "error"}, // fitsverify: "Data fill area invalid"
		{"theap_gap_exact", "", []int32{9}, "ok"},
		{"theap_gap_plus1", "lies outside the heap", nil, "ok"}, // fitsverify measures against PCOUNT, not PCOUNT minus the gap
		// A negative offset (int32 -4; fitsverify reads it as 4294967292,
		// past the heap): without the offset check the reader would panic
		// (no THEAP) or return the gap byte 7 (THEAP 12).
		{"neg_off", "invalid VLA descriptor (1, -4)", nil, "error"},
		{"neg_off_gap", "invalid VLA descriptor (1, -4)", nil, "error"},
		// THEAP = 12.0: fitsverify "value = 12.0 is not an integer". Read
		// as the default THEAP it would return the gap byte 7.
		{"theap_float", "THEAP = 12 is not an integer", nil, "error"},
	}
	if len(verdicts) != len(cases) {
		t.Fatalf("fitsverify.txt has %d verdicts, %d cases", len(verdicts), len(cases))
	}
	for _, c := range cases {
		file := "heap_" + c.name + ".fits"
		if verdicts[file] != c.fv {
			t.Errorf("%s: fitsverify says %q, the test expects %q", file, verdicts[file], c.fv)
		}
		path := testutil.TestdataPath(t, "tables", "heap", file)
		testutil.RequireFile(t, path)
		bt := openTable(t, path).(*BinaryTableHDU)
		v, err := ReadVarColumn[int32](bt, 1)
		if c.wantErr != "" {
			if err == nil || !strings.Contains(err.Error(), c.wantErr) {
				t.Errorf("%s: err = %v, want %q", file, err, c.wantErr)
			}
			continue
		}
		if err != nil {
			t.Errorf("%s: %v", file, err)
			continue
		}
		if got := v.At(0); fmt.Sprint(got) != fmt.Sprint(c.want) {
			t.Errorf("%s: row 0 = %v, want %v", file, got, c.want)
		}
	}
}
