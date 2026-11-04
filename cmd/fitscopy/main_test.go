package main

import (
	"bytes"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Every committed sample copies byte-for-byte with exit status 0.
func TestCopyEverySample(t *testing.T) {
	var samples []string
	for _, pat := range []string{"*.fits", "*.fit", "*.std"} {
		m, _ := filepath.Glob(testutil.TestdataPath(t, pat))
		samples = append(samples, m...)
	}
	sort.Strings(samples)
	if len(samples) < 43 {
		t.Fatalf("found %d samples, want at least 43", len(samples))
	}
	for _, in := range samples {
		t.Run(filepath.Base(in), func(t *testing.T) {
			out := filepath.Join(t.TempDir(), "out.fits")
			var errb bytes.Buffer
			if code := run([]string{in, out}, &errb); code != 0 {
				t.Fatalf("exit %d: %s", code, errb.String())
			}
			a, _ := os.ReadFile(in)
			b, err := os.ReadFile(out)
			if err != nil || !bytes.Equal(a, b) {
				t.Fatalf("output differs from input (%v)", err)
			}
		})
	}
}

// Copying a file onto itself (under any spelling) must not truncate it.
// Before, os.Create truncated the input before CopyTo read it.
func TestCopyOntoItselfRefused(t *testing.T) {
	src := testutil.TestdataPath(t, "iter_a.fit")
	orig, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	dir := t.TempDir()
	p := filepath.Join(dir, "x.fits")
	if err := os.WriteFile(p, orig, 0o644); err != nil {
		t.Fatal(err)
	}
	link := filepath.Join(dir, "link.fits")
	if err := os.Symlink(p, link); err != nil {
		t.Fatal(err)
	}
	for _, out := range []string{p, filepath.Join(dir, ".", "x.fits"), link} {
		var errb bytes.Buffer
		if code := run([]string{p, out}, &errb); code != 2 {
			t.Errorf("fitscopy %s %s: exit %d (%s), want 2", p, out, code, errb.String())
		}
		if b, _ := os.ReadFile(p); !bytes.Equal(b, orig) {
			t.Fatalf("fitscopy onto itself (%s) changed the file", out)
		}
	}
}

func TestCopyUsageAndErrors(t *testing.T) {
	var errb bytes.Buffer
	if code := run([]string{"one"}, &errb); code != 2 || !strings.Contains(errb.String(), "usage") {
		t.Fatalf("one arg: exit %d, %q", code, errb.String())
	}
	errb.Reset()
	dir := t.TempDir()
	if code := run([]string{filepath.Join(dir, "missing.fits"), filepath.Join(dir, "o.fits")}, &errb); code != 1 {
		t.Fatalf("missing input: exit %d", code)
	}
	// Not FITS: exit 1, nothing claimed identical.
	bad := filepath.Join(dir, "bad.fits")
	os.WriteFile(bad, bytes.Repeat([]byte("x"), 2880), 0o644)
	errb.Reset()
	if code := run([]string{bad, filepath.Join(dir, "o2.fits")}, &errb); code != 1 {
		t.Fatalf("non-FITS input: exit %d (%s)", code, errb.String())
	}
}

// Inputs whose CopyTo output is not the input (bytes after the last HDU,
// a last data unit without its final padding; see (*File).CopyTo) are
// reported as different (exit 3), not as copied.
func TestCopyReportsInputsCopyToChanges(t *testing.T) {
	read := func(elem ...string) []byte {
		p := testutil.TestdataPath(t, elem...)
		testutil.RequireFile(t, p)
		b, err := os.ReadFile(p)
		if err != nil {
			t.Fatal(err)
		}
		return b
	}
	prog := read("testprog.std")
	f32 := read("images", "cfwrite", "f32_plain.fits")
	for name, in := range map[string][]byte{
		"trailing junk": append(append([]byte(nil), prog...), []byte("thirteen junk")...),
		"unpadded":      f32[:2924],
	} {
		dir := t.TempDir()
		src := filepath.Join(dir, "in.fits")
		if err := os.WriteFile(src, in, 0o644); err != nil {
			t.Fatal(err)
		}
		var errb bytes.Buffer
		if code := run([]string{src, filepath.Join(dir, "out.fits")}, &errb); code != 3 || !strings.Contains(errb.String(), "differs") {
			t.Errorf("%s: exit %d, stderr %q; want 3 and a difference warning", name, code, errb.String())
		}
	}
}

// When the written file differs from the input, fitscopy says so and exits
// 3, whether one byte changed or the output is short. The copy step is
// replaced to produce such outputs.
func TestCopyDifferenceReported(t *testing.T) {
	defer func() { copyFn = copyFile }()
	src := testutil.TestdataPath(t, "vari.fits")
	testutil.RequireFile(t, src)
	orig, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	for name, bad := range map[string]func([]byte) []byte{
		"one byte": func(b []byte) []byte { b[len(b)/2] ^= 1; return b },
		"short":    func(b []byte) []byte { return b[:len(b)-2880] },
	} {
		copyFn = func(in, out string) error {
			if err := copyFile(in, out); err != nil {
				return err
			}
			b, err := os.ReadFile(out)
			if err != nil {
				return err
			}
			return os.WriteFile(out, bad(b), 0o644)
		}
		var errb bytes.Buffer
		code := run([]string{src, filepath.Join(t.TempDir(), "o.fits")}, &errb)
		if code != 3 || !strings.Contains(errb.String(), "differs") {
			t.Errorf("%s: exit %d, stderr %q; want 3 and a difference warning", name, code, errb.String())
		}
	}
	if b, _ := os.ReadFile(src); !bytes.Equal(b, orig) {
		t.Fatal("input changed")
	}
}
