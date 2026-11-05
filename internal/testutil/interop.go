//go:build interop

package testutil

import (
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// Tier 2 helpers (-tags interop): everything here starts a reference
// program. Opting in to tier 2 means every tool is required: a missing
// one is a test failure that names it, never a skip.

// RequireTool returns the path of an external program (fitsverify, fpack,
// funpack, ...) found on PATH, failing the test when it is missing.
func RequireTool(t testing.TB, name string) string {
	t.Helper()
	p, err := exec.LookPath(name)
	if err != nil {
		t.Fatalf("tier 2 (-tags interop) needs %s: not found on PATH (%v); %s", name, err, toolHint(name))
	}
	return p
}

// toolHint says where a missing tool comes from.
func toolHint(name string) string {
	switch name {
	case "fpack", "funpack":
		return "install cfitsio (Arch: extra/cfitsio ships fpack and funpack)"
	case "fitsverify":
		return "install HEASARC fitsverify (packaged separately from cfitsio)"
	case "xy2sky", "sky2xy":
		return "install WCSTools 3.9.7 (the TNX truth source, see wcs/transform/testdata/tnx/PROVENANCE.md)"
	case "HPXcvt":
		return "install wcslib 8.10 (Arch: extra/wcslib ships HPXcvt)"
	case "cc", "c++", "make":
		return "install a C/C++ toolchain (gcc, make) to build the reference harnesses"
	case "git":
		return "install git"
	case "python3":
		return "set FITS_ASTROPY_PYTHON to the repo venv interpreter"
	}
	return "install it"
}

// Harness returns the absolute path of a built reference harness (rel is
// repository-relative, e.g. "testdata/cfref/cfref"), failing the test
// when it has not been built.
func Harness(t testing.TB, rel string) string {
	t.Helper()
	bin := filepath.Join(RepoRoot(t), filepath.FromSlash(rel))
	st, err := os.Stat(bin)
	if err != nil || st.IsDir() || st.Mode()&0o111 == 0 {
		t.Fatalf("tier 2 (-tags interop) needs the reference harness %s: not built; build it with: make -C %s", bin, filepath.Dir(rel))
	}
	return bin
}

// CfrefBinary is the built cfitsio harness testdata/cfref/cfref.
func CfrefBinary(t testing.TB) string {
	t.Helper()
	return Harness(t, ToolCfref.Harness)
}

// AstropyPython returns the FITS_ASTROPY_PYTHON interpreter after checking
// that it imports astropy and numpy. Unset or broken fails the test.
func AstropyPython(t testing.TB) string {
	t.Helper()
	py := os.Getenv(EnvAstropyPython)
	if py == "" {
		t.Fatalf("tier 2 (-tags interop) needs %s: unset; point it at a Python with astropy (e.g. the repo venv: %s=$PWD/venv/bin/python)", EnvAstropyPython, EnvAstropyPython)
	}
	if _, err := exec.LookPath(py); err != nil {
		t.Fatalf("%s=%q is not an executable: %v", EnvAstropyPython, py, err)
	}
	if v := pythonVersions(t); v == "" {
		t.Fatalf("%s=%q cannot import astropy/numpy", EnvAstropyPython, py)
	}
	return py
}

var (
	pyMu   sync.Mutex
	pyVers = map[string]string{}
)

// pythonVersions is "astropy <v>, numpy <v>" of the FITS_ASTROPY_PYTHON
// interpreter, "" if it cannot import them (cached per interpreter).
func pythonVersions(t testing.TB) string {
	t.Helper()
	py := os.Getenv(EnvAstropyPython)
	if py == "" {
		return ""
	}
	pyMu.Lock()
	defer pyMu.Unlock()
	if v, ok := pyVers[py]; ok {
		return v
	}
	out, err := exec.Command(py, "-c", "import astropy, numpy; print('astropy', astropy.__version__ + ',', 'numpy', numpy.__version__)").CombinedOutput()
	v := strings.TrimSpace(string(out))
	if err != nil {
		t.Logf("%s=%s: %v: %s", EnvAstropyPython, py, err, v)
		v = ""
	}
	pyVers[py] = v
	return v
}

// PythonWith returns AstropyPython(t) after checking that it can import
// every module in mods (healpy, cv2, mpmath, ...).
func PythonWith(t testing.TB, mods ...string) string {
	t.Helper()
	py := AstropyPython(t)
	if len(mods) > 0 {
		if out, err := exec.Command(py, "-c", "import "+strings.Join(mods, ", ")).CombinedOutput(); err != nil {
			t.Fatalf("%s=%s cannot import %s: %v\n%s", EnvAstropyPython, py, strings.Join(mods, ", "), err, out)
		}
	}
	return py
}

// ExternalSource returns the directory named by env after checking that
// probe (a path relative to it) exists; unset or wrong fails the test.
func ExternalSource(t testing.TB, env, probe, what string) string {
	t.Helper()
	dir := os.Getenv(env)
	if dir == "" {
		t.Fatalf("tier 2 (-tags interop) needs %s: unset (%s)", env, what)
	}
	if _, err := os.Stat(filepath.Join(dir, probe)); err != nil {
		t.Fatalf("%s=%s does not hold %s (%s): %v", env, dir, probe, what, err)
	}
	return dir
}

// RequireGoldenPlatform skips a byte-for-byte golden regeneration test on
// platforms where the reference libraries are not expected to reproduce
// the committed bytes. The goldens are defined by the reference tools
// built for amd64, where no compiler fuses x*y+z into an FMA. On arm64
// (Apple Silicon, Homebrew builds) clang contracts multiply-adds inside
// cfitsio, wcslib and WCSTools: measured on macOS 27 / arm64, wcslib 8.10
// gives identical accept/reject decisions on every committed case, sky
// positions within 0.009 arcsec (CSC, which computes in float) and
// mostly within a few ULP, and cfitsio's noise sigma differs by a few ULP.
// The library itself is FMA-proof, and the tests that compare it with the
// committed goldens run on every platform (tier 1); only the regeneration
// of the goldens (tier 2) is pinned to amd64. This is a platform skip.
func RequireGoldenPlatform(t testing.TB) {
	t.Helper()
	if !ReferencePlatform() {
		t.Skipf("PLATFORM SKIP: golden regeneration is byte-exact on %s only, this is %s/%s (see testutil.RequireGoldenPlatform)", referencePlatform, runtime.GOOS, runtime.GOARCH)
	}
}

// referencePlatform is where the committed goldens and tool readings are
// defined: the reference tools as built for glibc on amd64. Elsewhere the
// same tools differ in ways that are theirs, not this library's: clang on
// arm64 contracts multiply-adds inside them, macOS libc's strtod refuses
// some values at the normal/subnormal boundary that glibc accepts (cfitsio
// then falls back to BSCALE = 1), and the malloc fill used to pin
// cfitsio's unwritten buffers is 0xAA on macOS, 0xA5 on glibc.
const referencePlatform = "linux/amd64"

// ReferencePlatform reports whether this is the reference platform.
func ReferencePlatform() bool { return runtime.GOOS+"/"+runtime.GOARCH == referencePlatform }

// MakeVerify builds a harness with `make -C dir <build...>` and then runs
// `make -C dir <targets...> <vars...>`, whose recipes regenerate the
// goldens into a scratch directory and byte-compare them with the
// committed ones (the harness Makefiles' verify targets).
func MakeVerify(t testing.TB, dir string, build, targets []string, vars ...string) {
	t.Helper()
	mk := RequireTool(t, "make")
	if len(build) > 0 {
		args := append([]string{"-C", dir}, build...)
		if out, err := exec.Command(mk, append(args, vars...)...).CombinedOutput(); err != nil {
			t.Fatalf("cannot build the harness in %s (needs cc and the reference library): %v\n%s", dir, err, out)
		}
	}
	if len(targets) == 0 {
		return
	}
	args := append([]string{"-C", dir}, targets...)
	out, err := exec.Command(mk, append(args, vars...)...).CombinedOutput()
	if err != nil {
		t.Fatalf("make -C %s %s: the committed goldens do not reproduce: %v\n%s", dir, strings.Join(targets, " "), err, out)
	}
	t.Logf("make -C %s %s:\n%s", dir, strings.Join(targets, " "), lastLines(out, 4))
}

// RunIn runs name args... with working directory dir; a non-zero exit is a
// test failure (the generator itself refused or crashed).
func RunIn(t testing.TB, dir, name string, args ...string) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("%s %s (in %s): %v\n%s", name, strings.Join(args, " "), dir, err, lastLines(out, 20))
	}
}
