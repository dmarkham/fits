// Package testutil holds shared helpers for tests that check the fits
// library against independent truth sources: the cfitsio harness in
// testdata/cfref, the fitsverify/fpack/funpack binaries and astropy.
//
// It is an ordinary (non-_test) package so tests in any package of this
// module can import it. Nothing outside tests should import it.
//
// The tests run in two tiers (docs/VALIDATION.md, "Test tiers"):
//
//   - Tier 1, plain `go test ./...`, needs nothing but Go and the files
//     committed in the repository. Every check compares the library with
//     committed goldens, samples and recorded tool readings. Nothing here
//     that tier 1 compiles starts another program: a reference tool's
//     reading of a file the library wrote is replayed from
//     testdata/recorded (Run, record.go), found by the sha256 of the
//     exact bytes, so tier 1 also requires the library to write the very
//     bytes the tools verified.
//   - Tier 2, `go test -tags interop ./...`, compiles interop.go and
//     record_live.go: the same Run calls start the reference programs live
//     and compare their reading with the recorded one (FITS_RECORD=1
//     records it), and the *Reproduce* tests regenerate the goldens. A
//     tool that is missing there is a test failure naming the tool.
//
// A committed golden or sample that is missing is a failure (t.Fatal) in
// both tiers, never a skip.
package testutil

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"sync"
	"testing"
)

// EnvAstropyPython names a Python interpreter that has astropy installed.
const EnvAstropyPython = "FITS_ASTROPY_PYTHON"

const modulePath = "github.com/dmarkham/fits"

// External source trees used only by the golden generators (tier 2). No
// test has a default location for them: unset is a failure naming the
// variable (ExternalSource).
const (
	// EnvSirilSrc is an unpacked Siril 1.4.0 source tree (the directory
	// holding src/), e.g. the release tarball siril-1.4.0.tar.bz2.
	EnvSirilSrc = "FITS_SIRIL_SRC"
	// EnvHealpixCxx is a healpix_cxx 8.3.0 source directory (the one
	// holding healpix_base.cc), e.g. Siril's subprojects/healpix_cxx.
	EnvHealpixCxx = "FITS_HEALPIX_CXX"
	// EnvSirilCLI is a siril-cli 1.4.0 executable.
	EnvSirilCLI = "FITS_SIRIL_CLI"
)

var (
	rootOnce sync.Once
	rootDir  string
	rootErr  error
)

// RepoRoot returns the absolute path of the module root (the directory
// holding go.mod for github.com/dmarkham/fits), found by walking up from
// the test's working directory.
func RepoRoot(t testing.TB) string {
	t.Helper()
	rootOnce.Do(func() {
		dir, err := os.Getwd()
		if err != nil {
			rootErr = err
			return
		}
		for {
			b, err := os.ReadFile(filepath.Join(dir, "go.mod"))
			if err == nil && isModule(b) {
				rootDir = dir
				return
			}
			parent := filepath.Dir(dir)
			if parent == dir {
				rootErr = fmt.Errorf("no go.mod for module %s above the working directory", modulePath)
				return
			}
			dir = parent
		}
	})
	if rootErr != nil {
		t.Fatalf("testutil.RepoRoot: %v", rootErr)
	}
	return rootDir
}

func isModule(gomod []byte) bool {
	for _, line := range bytes.Split(gomod, []byte("\n")) {
		f := bytes.Fields(line)
		if len(f) == 2 && string(f[0]) == "module" && string(f[1]) == modulePath {
			return true
		}
	}
	return false
}

// TestdataPath joins elem onto <repo>/testdata.
func TestdataPath(t testing.TB, elem ...string) string {
	t.Helper()
	return filepath.Join(append([]string{RepoRoot(t), "testdata"}, elem...)...)
}

// RequireFile fails the test when a committed sample or golden is missing.
func RequireFile(t testing.TB, path string) {
	t.Helper()
	st, err := os.Stat(path)
	if err != nil {
		t.Fatalf("committed test file missing (this is a failure, not a skip): %v", err)
	}
	if st.IsDir() {
		t.Fatalf("committed test file %s is a directory", path)
	}
}
