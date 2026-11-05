package fits

import (
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCfrefGoldenHDUAddresses is the smoke test for the shared cfitsio
// truth harness (testdata/cfref): the scanner's HDU byte layout for
// testprog.std must equal what cfitsio's fits_get_hduaddrll reports.
// Per-area tests assert the rest of the golden (cards, pixels, columns).
func TestCfrefGoldenHDUAddresses(t *testing.T) {
	const name = "testprog.std"
	want := testutil.LoadCfrefGolden(t, name)
	f, err := Open(testutil.TestdataPath(t, name))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	if f.NumHDU() != want.NHDU {
		t.Fatalf("NumHDU = %d, cfitsio = %d", f.NumHDU(), want.NHDU)
	}
	for i, h := range want.HDUs {
		r := f.hdus[i]
		if r.headerStart != h.HeaderStart || r.dataStart != h.DataStart || r.paddedEnd != h.DataEnd {
			t.Errorf("HDU %d: header/data/end = %d/%d/%d, cfitsio = %d/%d/%d",
				i, r.headerStart, r.dataStart, r.paddedEnd, h.HeaderStart, h.DataStart, h.DataEnd)
		}
	}
}
