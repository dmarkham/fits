//go:build interop

package fits

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCfcopyReproducesGoldens: the committed cfitsio goldens are exactly
// what testdata/cfcopy/cfcopy writes now (reproducibility, and a check that
// the scenario list here matches the harness).
func TestCfcopyReproducesGoldens(t *testing.T) {
	testutil.RequireGoldenPlatform(t)
	bin := testutil.TestdataPath(t, "cfcopy", "cfcopy")
	if st, err := os.Stat(bin); err != nil || st.Mode()&0o111 == 0 {
		t.Fatalf("cfitsio copy harness %s not built; build it with: make -C testdata/cfcopy (needs libcfitsio)", bin)
	}
	// The reserved-space input is itself cfitsio output (see the Makefile).
	res := filepath.Join(t.TempDir(), "reserve.fits")
	if b, err := exec.Command(bin, "-r", "40", res, testutil.TestdataPath(t, "testprog.std"), "0", "2", "6").CombinedOutput(); err != nil {
		t.Fatalf("cfcopy -r 40: %v: %s", err, b)
	}
	resWant := testutil.TestdataPath(t, "cfcopy", "reserve.fits")
	testutil.RequireFile(t, resWant)
	if !bytes.Equal(mustRead(t, res), mustRead(t, resWant)) {
		t.Error("cfcopy -r 40 output differs from the committed reserve.fits")
	}
	mk := testutil.TestdataPath(t, "cfcopy", "mkcomments")
	if st, err := os.Stat(mk); err != nil || st.Mode()&0o111 == 0 {
		t.Fatalf("cfitsio harness %s not built; build it with: make -C testdata/cfcopy (needs libcfitsio)", mk)
	}
	com := filepath.Join(t.TempDir(), "comments.fits")
	if b, err := exec.Command(mk, com).CombinedOutput(); err != nil {
		t.Fatalf("mkcomments: %v: %s", err, b)
	}
	comWant := testutil.TestdataPath(t, "cfcopy", "comments.fits")
	testutil.RequireFile(t, comWant)
	if !bytes.Equal(mustRead(t, com), mustRead(t, comWant)) {
		t.Error("mkcomments output differs from the committed comments.fits")
	}
	gs, _ := filepath.Glob(testutil.TestdataPath(t, "cfcopy", "golden", "*.fits"))
	if len(gs) != len(editCopyScenarios) {
		t.Errorf("%d goldens committed, %d scenarios here", len(gs), len(editCopyScenarios))
	}
	for _, sc := range editCopyScenarios {
		t.Run(sc.name, func(t *testing.T) {
			out := filepath.Join(t.TempDir(), sc.name+".fits")
			args := []string{out, testutil.TestdataPath(t, sc.sample)}
			for _, i := range sc.order {
				args = append(args, strconv.Itoa(i))
			}
			if b, err := exec.Command(bin, args...).CombinedOutput(); err != nil {
				t.Fatalf("cfcopy: %v: %s", err, b)
			}
			if !bytes.Equal(mustRead(t, out), mustRead(t, cfcopyGolden(t, sc.name))) {
				t.Fatal("cfcopy output differs from the committed golden")
			}
		})
	}
}
