//go:build interop

package testutil

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

// Tier 2 is opted into: a missing tool is a failure naming it, never a skip.
func TestMissingToolFails(t *testing.T) {
	const missing = "fits-testutil-no-such-tool"
	r := observe(t, func(tb testing.TB) { RequireTool(tb, missing) })
	if !r.failed || r.skipped || !strings.Contains(r.msg, missing) {
		t.Fatalf("want failure naming the tool, got skipped=%v failed=%v %q", r.skipped, r.failed, r.msg)
	}
	r = observe(t, func(tb testing.TB) { Harness(tb, "testdata/no-such-harness/h") })
	if !r.failed || r.skipped || !strings.Contains(r.msg, "make -C testdata/no-such-harness") {
		t.Fatalf("want failure naming the harness and its build, got skipped=%v failed=%v %q", r.skipped, r.failed, r.msg)
	}
	for _, env := range []string{EnvSirilSrc, EnvHealpixCxx, EnvSirilCLI} {
		t.Setenv(env, "")
		r = observe(t, func(tb testing.TB) { ExternalSource(tb, env, "x", "test") })
		if !r.failed || r.skipped || !strings.Contains(r.msg, env) {
			t.Fatalf("unset %s: want failure naming it, got skipped=%v failed=%v %q", env, r.skipped, r.failed, r.msg)
		}
	}
}

func TestAstropyPythonBrokenFails(t *testing.T) {
	// Unset, missing or unable to import astropy: a failure naming the
	// variable.
	t.Setenv(EnvAstropyPython, "")
	r := observe(t, func(tb testing.TB) { AstropyPython(tb) })
	if !r.failed || !strings.Contains(r.msg, EnvAstropyPython) {
		t.Fatalf("unset %s: want failure naming it, got skipped=%v %q", EnvAstropyPython, r.skipped, r.msg)
	}
	t.Setenv(EnvAstropyPython, filepath.Join(t.TempDir(), "no-python-here"))
	r = observe(t, func(tb testing.TB) { AstropyPython(tb) })
	if !r.failed {
		t.Fatalf("want failure for missing explicit interpreter, got skipped=%v %q", r.skipped, r.msg)
	}
	// An explicit interpreter whose probe fails (/bin/false exits 1) must
	// also fail rather than skip.
	if _, err := os.Stat("/bin/false"); err == nil {
		t.Setenv(EnvAstropyPython, "/bin/false")
		r = observe(t, func(tb testing.TB) { AstropyPython(tb) })
		if !r.failed {
			t.Fatalf("want failure for interpreter that cannot import astropy, got skipped=%v %q", r.skipped, r.msg)
		}
	}
}

func TestAstropyPythonLive(t *testing.T) {
	py := AstropyPython(t)
	if py == "" {
		t.Fatal("empty interpreter path")
	}
}

// The committed goldens must be exactly what the harness produces now
// (reproducibility; also catches a libcfitsio change on this machine).
var poolRefFullRe = regexp.MustCompile(`\{"n": [0-9]+, "sha256": "([0-9a-f]{64})", "pool": "([a-z0-9]+)"\}`)

// They are made with "dump -pool" (cfref.c jarr), so the pool files the
// dump writes must be the committed ones too, and the golden decoded
// through the pool must equal the live explicit dump.
func TestCfrefDumpReproducesGoldens(t *testing.T) {
	RequireGoldenPlatform(t)
	bin := CfrefBinary(t)
	pool := t.TempDir()
	for _, s := range samples(t) {
		name := filepath.Base(s)
		t.Run(name, func(t *testing.T) {
			got, err := exec.Command(bin, "dump", "-pool", pool, s).Output()
			if err != nil {
				t.Fatalf("cfref dump %s: %v", s, err)
			}
			want, err := os.ReadFile(CfrefGoldenPath(t, name))
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(got, want) {
				t.Fatalf("harness output differs from committed golden %s; regenerate with make -C testdata/cfref goldens and review the diff", CfrefGoldenPath(t, name))
			}
			// The typed decoder must accept the live output too, and the
			// golden with every pool reference replaced by its pool file
			// must be the live explicit dump byte for byte.
			if cf := CfrefDump(t, s); cf.File != name {
				t.Fatalf("decoded file field %q", cf.File)
			}
			live, err := exec.Command(bin, "dump", s).Output()
			if err != nil {
				t.Fatal(err)
			}
			expanded := poolRefFullRe.ReplaceAllFunc(want, func(m []byte) []byte {
				sm := poolRefFullRe.FindSubmatch(m)
				b, err := ReadGolden(filepath.Join(cfrefPoolDir(), string(sm[2])+"-"+string(sm[1])+".json"))
				if err != nil {
					t.Fatal(err)
				}
				return bytes.TrimSuffix(b, []byte("\n"))
			})
			expanded = expandPooledCards(t, expanded)
			if !bytes.Equal(expanded, live) {
				t.Fatalf("golden with its pool files inlined differs from the explicit dump")
			}
		})
	}
	// The pool also holds arrays only the twin goldens (golden/twins, of
	// the uncommitted decompressed twins) reference; those are reproduced
	// by TestCompressedGoldensReproduce, and TestCfrefPoolComplete allows
	// no pool file that no golden references.
	files, _ := filepath.Glob(filepath.Join(pool, "*.json"))
	committed := GlobGoldens(filepath.Join(cfrefPoolDir(), "*.json"))
	if len(files) > len(committed) || len(files) < len(committed)/2 {
		t.Errorf("dump wrote %d pool files, %d are committed", len(files), len(committed))
	}
	for _, f := range files {
		got, _ := os.ReadFile(f)
		want, err := ReadGolden(filepath.Join(cfrefPoolDir(), filepath.Base(f)))
		if err != nil || !bytes.Equal(got, want) {
			t.Errorf("pool file %s differs from the committed one (%v)", filepath.Base(f), err)
		}
	}
}

func TestRequireToolLive(t *testing.T) {
	for _, n := range []string{"fitsverify", "fpack", "funpack"} {
		if p := RequireTool(t, n); p == "" {
			t.Fatalf("%s: empty path", n)
		}
	}
}

// The pure-Go git index reader (TestRepositorySizeBudget in tier 1) lists
// exactly what `git ls-files` lists.
func TestTrackedFilesMatchesGit(t *testing.T) {
	git := RequireTool(t, "git")
	root := RepoRoot(t)
	out, err := exec.Command(git, "-C", root, "ls-files", "-z").Output()
	if err != nil {
		t.Fatalf("git ls-files: %v", err)
	}
	want := strings.Split(strings.TrimRight(string(out), "\x00"), "\x00")
	got, ok, err := TrackedFiles(root)
	if err != nil || !ok {
		t.Fatalf("TrackedFiles: ok=%v err=%v", ok, err)
	}
	if strings.Join(got, "\n") != strings.Join(want, "\n") {
		t.Fatalf("TrackedFiles lists %d files, git ls-files %d", len(got), len(want))
	}
}
