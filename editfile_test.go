package fits

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
)

// writeThreeHDUFile writes a primary + two image extensions (EXTNAME KEEP
// and JUNK) and returns its path and bytes.
func writeThreeHDUFile(t *testing.T) (string, []byte) {
	t.Helper()
	path := filepath.Join(t.TempDir(), "multi.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	h1 := header.New()
	h1.Set("EXTNAME", "KEEP", "")
	if _, err := AppendImage(f, h1, []int64{4}, []float32{1, 2, 3, 4}); err != nil {
		t.Fatal(err)
	}
	h2 := header.New()
	h2.Set("EXTNAME", "JUNK", "")
	if _, err := AppendImage(f, h2, []int64{2}, []float32{9, 9}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return path, b
}

// assertUntouched checks that path still holds exactly orig and that no
// other file was left in its directory.
func assertUntouched(t *testing.T, path string, orig []byte) {
	t.Helper()
	after, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(after, orig) {
		t.Fatalf("source changed (len %d -> %d, first difference at byte %d)", len(orig), len(after), firstDiff(after, orig))
	}
	ents, err := os.ReadDir(filepath.Dir(path))
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range ents {
		if e.Name() != filepath.Base(path) {
			t.Fatalf("stray file after failed EditFile: %s", e.Name())
		}
	}
}

func TestEditFileDropHDU(t *testing.T) {
	path, orig := writeThreeHDUFile(t)

	// Drop the JUNK HDU through EditFile.
	err := EditFile(path, func(in *File, out *Writer) error {
		for i := range in.NumHDU() {
			h, err := in.HDU(i)
			if err != nil {
				return err
			}
			if name, _ := h.Header().String("EXTNAME"); name == "JUNK" {
				if err := out.SkipHDU(h); err != nil {
					return err
				}
				continue
			}
			if err := out.CopyHDU(h); err != nil {
				return err
			}
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	// The result is exactly the first two HDUs of the original, byte for
	// byte (nothing changes position, so nothing is converted).
	in, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer in.Close()
	if in.NumHDU() != 2 {
		t.Fatalf("NumHDU after drop = %d", in.NumHDU())
	}
	got, _ := os.ReadFile(path)
	src, err := OpenReader(bytes.NewReader(orig))
	if err != nil {
		t.Fatal(err)
	}
	if want := orig[:src.hdus[1].paddedEnd]; !bytes.Equal(got, want) {
		t.Fatalf("output is not the first two HDUs of the input byte for byte (first diff at %d)", firstDiff(got, want))
	}
	h, _ := in.HDU(1)
	if name, _ := h.Header().String("EXTNAME"); name != "KEEP" {
		t.Fatalf("HDU 1 EXTNAME = %q", name)
	}
}

// A callback error leaves the source byte-identical and no temp file.
// (The earlier version of this test compared only the length, and passed
// against an EditFile that overwrote a data byte of the source.)
func TestEditFileFailureLeavesSourceIntact(t *testing.T) {
	path, orig := writeThreeHDUFile(t)
	err := EditFile(path, func(in *File, out *Writer) error {
		h, _ := in.HDU(0)
		if err := out.CopyHDU(h); err != nil {
			return err
		}
		if err := out.CopyHDU(&ImageHDU{rec: &hduRecord{}}); err == nil {
			return fmt.Errorf("CopyHDU accepted an HDU with no backing file")
		}
		return errTestForce
	})
	if !errors.Is(err, errTestForce) {
		t.Fatalf("EditFile returned %v, want the callback's error", err)
	}
	assertUntouched(t, path, orig)
}

// A callback panic is returned as an error; the source is byte-identical.
func TestEditFilePanicLeavesSourceIntact(t *testing.T) {
	path, orig := writeThreeHDUFile(t)
	err := EditFile(path, func(in *File, out *Writer) error {
		h, _ := in.HDU(0)
		out.CopyHDU(h)
		panic("boom")
	})
	if err == nil {
		t.Fatal("EditFile returned nil after a panic")
	}
	assertUntouched(t, path, orig)
}

// A callback that returns nil without writing any HDU must not replace the
// source with an empty file.
func TestEditFileNoHDURefused(t *testing.T) {
	path, orig := writeThreeHDUFile(t)
	if err := EditFile(path, func(in *File, out *Writer) error { return nil }); err == nil {
		t.Fatal("EditFile with no HDU written returned nil")
	}
	assertUntouched(t, path, orig)
}

// CopyHDUWithHeader refuses a header describing a different data size
// (the old doc said the library does not check; passing the bytes through
// would corrupt every later HDU).
func TestEditFileCopyHDUWithHeaderRejectsSizeChange(t *testing.T) {
	path, orig := writeThreeHDUFile(t)
	err := EditFile(path, func(in *File, out *Writer) error {
		p, _ := in.HDU(0)
		if err := out.CopyHDU(p); err != nil {
			return err
		}
		h, _ := in.HDU(1)
		hdr := h.Header().Clone()
		hdr.Set("NAXIS1", int64(5), "")
		return out.CopyHDUWithHeader(h, hdr)
	})
	if err == nil {
		t.Fatal("CopyHDUWithHeader accepted a header with a different data size")
	}
	assertUntouched(t, path, orig)
}

// CopyHDUWithHeader refuses a header that keeps the data size but changes
// what kind of HDU the bytes are (image vs table, binary vs ASCII table,
// compressed vs plain): the data would be reinterpreted.
func TestEditFileCopyHDUWithHeaderRejectsKindChange(t *testing.T) {
	findKind := func(t *testing.T, f *File, k hduKind) int {
		for i, r := range f.hdus {
			if r.kind == k {
				return i
			}
		}
		t.Fatalf("no HDU of kind %d", k)
		return -1
	}
	cases := []struct {
		name, sample string
		kind         hduKind
		key          string
		val          any
	}{
		{"image-over-bintable", "testprog.std", kindBinTable, "XTENSION", "IMAGE"},
		{"bintable-over-ascii", "testprog.std", kindASCIITable, "XTENSION", "BINTABLE"},
		{"ascii-over-bintable", "testprog.std", kindBinTable, "XTENSION", "TABLE"},
		{"bintable-over-image", "testprog.std", kindImageExt, "XTENSION", "BINTABLE"},
		{"uncompress-zimage", "comp_rice_i16.fits", kindCompressed, "ZIMAGE", false},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			path, orig := copySample(t, c.sample)
			err := EditFile(path, func(in *File, out *Writer) error {
				k := findKind(t, in, c.kind)
				for i := 0; i < k; i++ {
					h, _ := in.HDU(i)
					if err := out.CopyHDU(h); err != nil {
						return err
					}
				}
				h, _ := in.HDU(k)
				hdr := h.Header().Clone()
				if err := hdr.Set(c.key, c.val, ""); err != nil {
					return err
				}
				if c.kind == kindImageExt {
					// A table header needs TFIELDS (else the refusal is
					// the missing keyword, not the kind change).
					if err := hdr.Set("TFIELDS", int64(0), ""); err != nil {
						return err
					}
				}
				return out.CopyHDUWithHeader(h, hdr)
			})
			if err == nil || !strings.Contains(err.Error(), "changes the HDU kind") {
				t.Fatalf("CopyHDUWithHeader error = %v, want \"changes the HDU kind\"", err)
			}
			assertUntouched(t, path, orig)
		})
	}
}

// CopyHDUWithHeader with a compatible header: the new card is there and the
// data bytes are the original's.
func TestEditFileCopyHDUWithHeader(t *testing.T) {
	path, orig := writeThreeHDUFile(t)
	err := EditFile(path, func(in *File, out *Writer) error {
		for i := range in.NumHDU() {
			h, _ := in.HDU(i)
			if i != 1 {
				if err := out.CopyHDU(h); err != nil {
					return err
				}
				continue
			}
			hdr := h.Header().Clone()
			hdr.Set("OBSERVER", "copy test", "added")
			if err := out.CopyHDUWithHeader(h, hdr); err != nil {
				return err
			}
		}
		return nil
	})
	if err != nil {
		t.Fatal(err)
	}
	src, _ := OpenReader(bytes.NewReader(orig))
	f, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	got, _ := os.ReadFile(path)
	for i := range src.hdus {
		s, g := src.hdus[i], f.hdus[i]
		if !bytes.Equal(orig[s.dataStart:s.paddedEnd], got[g.dataStart:g.paddedEnd]) {
			t.Errorf("HDU %d data bytes changed", i)
		}
	}
	h, _ := f.HDU(1)
	if v, err := h.Header().String("OBSERVER"); err != nil || v != "copy test" {
		t.Fatalf("OBSERVER = %q, %v", v, err)
	}
}

// The output keeps the source's permission bits exactly. 0600 checks the
// file is not widened; 0666 and 0606 check it is not narrowed by the umask
// (O_CREATE applies the umask, so only an explicit chmod gives these).
func TestEditFileKeepsPermissions(t *testing.T) {
	for _, mode := range []os.FileMode{0o600, 0o640, 0o666, 0o606} {
		t.Run(mode.String(), func(t *testing.T) {
			path, _ := writeThreeHDUFile(t)
			if err := os.Chmod(path, mode); err != nil {
				t.Fatal(err)
			}
			if err := EditFile(path, copyAll); err != nil {
				t.Fatal(err)
			}
			st, err := os.Stat(path)
			if err != nil {
				t.Fatal(err)
			}
			if st.Mode().Perm() != mode {
				t.Fatalf("mode after EditFile = %v, want %v", st.Mode().Perm(), mode)
			}
		})
	}
}

// A symbolic link is followed: the target file is replaced and the link is
// still a link.
func TestEditFileFollowsSymlink(t *testing.T) {
	path, orig := writeThreeHDUFile(t)
	link := filepath.Join(t.TempDir(), "link.fits")
	if err := os.Symlink(path, link); err != nil {
		t.Fatal(err)
	}
	err := EditFile(link, func(in *File, out *Writer) error {
		h, _ := in.HDU(0)
		return out.CopyHDU(h)
	})
	if err != nil {
		t.Fatal(err)
	}
	if st, err := os.Lstat(link); err != nil || st.Mode()&os.ModeSymlink == 0 {
		t.Fatalf("link replaced by a regular file (%v)", err)
	}
	got, _ := os.ReadFile(path)
	src, _ := OpenReader(bytes.NewReader(orig))
	if !bytes.Equal(got, orig[:src.hdus[0].paddedEnd]) {
		t.Fatal("link target does not hold the edited file")
	}
}

func copyAll(in *File, out *Writer) error {
	for i := range in.NumHDU() {
		h, err := in.HDU(i)
		if err != nil {
			return err
		}
		if err := out.CopyHDU(h); err != nil {
			return err
		}
	}
	return nil
}

// TestEditFileCrashAtEveryStep injects a crash after each step of EditFile
// (after the callback, after fsync, after rename, after the directory
// fsync). After recovery (by OpenForEdit, and separately by the next
// EditFile, see crashRecoverers) the file must be byte-identical to
// the original or to the fault-free result, monotonically, with no files
// left behind.
func TestEditFileCrashAtEveryStep(t *testing.T) {
	dropJunk := func(in *File, out *Writer) error {
		for i := 0; i < 2; i++ {
			h, _ := in.HDU(i)
			if err := out.CopyHDU(h); err != nil {
				return err
			}
		}
		return nil
	}
	path, pre := writeThreeHDUFile(t)
	var steps []string
	flushFaultHook = func(s string) error { steps = append(steps, s); return nil }
	err := EditFile(path, dropJunk)
	flushFaultHook = nil
	if err != nil {
		t.Fatal(err)
	}
	post, _ := os.ReadFile(path)
	if len(steps) == 0 || bytes.Equal(pre, post) {
		t.Fatalf("steps %v; changed %v", steps, !bytes.Equal(pre, post))
	}
	for _, rv := range crashRecoverers {
		t.Run(rv.name, func(t *testing.T) {
			reachedPost := false
			for k, step := range steps {
				dir := t.TempDir()
				p := filepath.Join(dir, "multi.fits")
				if err := os.WriteFile(p, pre, 0o644); err != nil {
					t.Fatal(err)
				}
				n := 0
				flushFaultHook = func(s string) error {
					if n == k {
						return errInjectedCrash
					}
					n++
					return nil
				}
				err := EditFile(p, dropJunk)
				flushFaultHook = nil
				if !errors.Is(err, errInjectedCrash) {
					t.Fatalf("step %s: EditFile = %v", step, err)
				}
				if err := rv.run(p); err != nil {
					t.Fatalf("crash after %s: recovery by %s: %v", step, rv.name, err)
				}
				got, _ := os.ReadFile(p)
				switch {
				case bytes.Equal(got, post):
					reachedPost = true
				case bytes.Equal(got, pre):
					if reachedPost {
						t.Errorf("crash after %s: completed edit lost", step)
					}
				default:
					t.Errorf("crash after %s: file is neither the original nor the result", step)
				}
				if s := sidecars(t, dir, "multi.fits"); len(s) != 0 {
					t.Errorf("crash after %s: left behind %v", step, s)
				}
			}
			if !reachedPost {
				t.Error("no crash point kept the completed edit")
			}
		})
	}
}

var errTestForce = &testError{"forced test failure"}

type testError struct{ msg string }

func (e *testError) Error() string { return e.msg }

// CopyHDUWithHeader re-encodes the caller's header, so it must refuse what
// the encoder cannot represent exactly, like Flush (lossyHeaderCases): the
// EditFile fails and the source stays byte-identical, or the output reads
// back in cfitsio with exactly what was set.
func TestEditFileCopyHDUWithHeaderRefusesLossyEncoding(t *testing.T) {
	for _, tc := range lossyHeaderCases {
		t.Run(tc.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			err := EditFile(path, func(in *File, out *Writer) error {
				for i := range in.NumHDU() {
					h, err := in.HDU(i)
					if err != nil {
						return err
					}
					if i != tc.hdu {
						if err := out.CopyHDU(h); err != nil {
							return err
						}
						continue
					}
					hdr := h.Header().Clone()
					if err := tc.mutate(hdr); err != nil {
						if tc.outcome != refusedBySet {
							t.Errorf("Set/Add refused the edit (%v); outcome %d expected", err, tc.outcome)
						}
						return err
					}
					if tc.outcome == refusedBySet {
						t.Errorf("Set/Add accepted the edit; it must refuse it")
					}
					if err := out.CopyHDUWithHeader(h, hdr); err != nil {
						if tc.outcome != refusedByWriter {
							t.Errorf("CopyHDUWithHeader refused (%v); outcome %d expected", err, tc.outcome)
						}
						return err
					}
					if tc.outcome == refusedByWriter {
						t.Errorf("CopyHDUWithHeader accepted the edit; it must refuse it")
					}
				}
				return nil
			})
			if err != nil {
				t.Logf("EditFile refused: %v", err)
				assertUntouched(t, path, orig)
				return
			}
			if tc.outcome != written {
				t.Fatalf("EditFile succeeded; outcome %d expected", tc.outcome)
			}
			tc.check(t, path)
		})
	}
}
