//go:build linux || darwin

package fits

import (
	"bytes"
	"errors"
	"os"
	"os/user"
	"path/filepath"
	"strconv"
	"syscall"
	"testing"
)

// The resizing Flush and EditFile replace the file with a temp file and a
// rename. What a rename alone would lose is either carried over or
// refused, never dropped silently (Linux and macOS):
//
//   - another hard link: refused with ErrHardLinked, every link left
//     byte-identical;
//   - the group (and owner): carried over; a chown the process may not make
//     is refused with ErrCannotPreserve, the file untouched;
//   - extended attributes: carried over, with POSIX ACLs on Linux
//     (replace_preserve_linux_interop_test.go) and the resource fork and Finder
//     info on macOS (replace_preserve_darwin_interop_test.go).
//
// Truth: the kernel's extended attributes read back with listxattr(2) /
// getxattr(2) before and after (tier 1: replace_xattr_linux_test.go,
// replace_xattr_darwin_test.go); getfattr / getfacl (Linux), xattr / stat /
// ls -le (macOS) in tier 2 (replace_preserve_{linux,darwin}_interop_test.go);
// byte equality for refusals.

// replacers are the two operations that replace the file.
var replacers = []struct {
	name string
	run  func(t *testing.T, path string) error
}{
	{"FlushResize", func(t *testing.T, path string) error {
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		growHDU1(t, f)
		return f.Flush()
	}},
	{"EditFile", func(t *testing.T, path string) error { return EditFile(path, copyAll) }},
}

func TestReplaceRefusesHardLinks(t *testing.T) {
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			other := filepath.Join(filepath.Dir(path), "other-name.fits")
			if err := os.Link(path, other); err != nil {
				t.Fatal(err)
			}
			err := r.run(t, path)
			if !errors.Is(err, ErrHardLinked) {
				t.Fatalf("%s on a file with two links: err = %v, want ErrHardLinked", r.name, err)
			}
			for _, p := range []string{path, other} {
				if !bytes.Equal(mustRead(t, p), orig) {
					t.Fatalf("%s changed after the refusal", p)
				}
			}
			a, _ := os.Stat(path)
			b, _ := os.Stat(other)
			if !os.SameFile(a, b) {
				t.Fatal("the two names no longer share an inode")
			}
			if s := sidecars(t, filepath.Dir(path), "testprog.std"); len(s) != 1 || s[0] != "other-name.fits" {
				t.Fatalf("left behind: %v", s)
			}
		})
	}
}

// supplementaryGroup returns a group of this user other than its primary
// one, or skips when there is none (root in a container, a CI user): a
// property of the account, not a missing tool.
func supplementaryGroup(t *testing.T) int {
	t.Helper()
	u, err := user.Current()
	if err != nil {
		t.Fatal(err)
	}
	gids, err := u.GroupIds()
	if err != nil {
		t.Fatal(err)
	}
	for _, g := range gids {
		if g != u.Gid {
			if n, err := strconv.Atoi(g); err == nil {
				return n
			}
		}
	}
	t.Skipf("ENVIRONMENT SKIP: the test needs a user with a supplementary group (to chgrp the file to it); %s has none", u.Username)
	return -1
}

// A chown the process is not allowed to make (another user's file) is
// refused, with the file untouched and no temp file left. The EPERM is
// injected: a non-root test cannot own a file it may not chown.
func TestReplaceRefusesWhenOwnerCannotBeKept(t *testing.T) {
	group := supplementaryGroup(t)
	chownFile = func(*os.File, int, int) error { return syscall.EPERM }
	defer func() { chownFile = (*os.File).Chown }()
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			if err := os.Chown(path, -1, group); err != nil {
				t.Fatal(err)
			}
			err := r.run(t, path)
			if !errors.Is(err, ErrCannotPreserve) {
				t.Fatalf("err = %v, want ErrCannotPreserve", err)
			}
			if !bytes.Equal(mustRead(t, path), orig) {
				t.Fatal("file changed")
			}
			if s := sidecars(t, filepath.Dir(path), "testprog.std"); len(s) != 0 {
				t.Fatalf("left behind: %v", s)
			}
		})
	}
}
