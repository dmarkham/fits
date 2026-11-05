//go:build darwin

package fits

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"syscall"
	"testing"

	"golang.org/x/sys/unix"
)

// Tier 1 twins of the /usr/bin/xattr tests (replace_preserve_darwin_interop_test.go,
// tier 2): the same setups and checks with the kernel's own interface,
// setxattr(2) to set and listxattr(2) / getxattr(2) to read back (the
// calls /usr/bin/xattr makes), so plain `go test` proves the replacements
// keep the group and every extended attribute on macOS without the xattr
// and stat programs. The ACL gap (TestReplaceDropsACLDarwin) stays tier 2:
// an ACL is set through libc's acl_set_file, not a system call Go reaches.

func setxattrDarwin(t *testing.T, path, name string, value []byte) {
	t.Helper()
	if err := unix.Setxattr(path, name, value, 0); err != nil {
		if errors.Is(err, unix.ENOTSUP) {
			t.Skipf("FILESYSTEM SKIP: %s does not support extended attributes: %v", filepath.Dir(path), err)
		}
		t.Fatalf("setxattr %s %s: %v", path, name, err)
	}
}

// xattrsDarwin is every extended attribute of path, name -> value.
func xattrsDarwin(t *testing.T, path string) map[string][]byte {
	t.Helper()
	size, err := unix.Listxattr(path, nil)
	if err != nil {
		t.Fatalf("listxattr %s: %v", path, err)
	}
	buf := make([]byte, size)
	n, err := unix.Listxattr(path, buf)
	if err != nil {
		t.Fatalf("listxattr %s: %v", path, err)
	}
	out := map[string][]byte{}
	for _, name := range strings.Split(strings.TrimRight(string(buf[:n]), "\x00"), "\x00") {
		if name == "" {
			continue
		}
		vs, err := unix.Getxattr(path, name, nil)
		if err != nil {
			t.Fatalf("getxattr %s %s: %v", path, name, err)
		}
		v := make([]byte, vs)
		vn, err := unix.Getxattr(path, name, v)
		if err != nil {
			t.Fatalf("getxattr %s %s: %v", path, name, err)
		}
		out[name] = v[:vn]
	}
	return out
}

func xattrDiff(got, want map[string][]byte) string {
	var d []string
	for k, v := range want {
		if g, ok := got[k]; !ok {
			d = append(d, "missing "+k)
		} else if !bytes.Equal(g, v) {
			d = append(d, "changed "+k)
		}
	}
	for k := range got {
		if _, ok := want[k]; !ok {
			d = append(d, "extra "+k)
		}
	}
	sort.Strings(d)
	return strings.Join(d, ", ")
}

// TestReplaceKeepsGroupAndXattrsDarwinKernel is
// TestReplaceKeepsGroupAndXattrsDarwin with the kernel as the truth: a
// user attribute, the Finder info, the quarantine flag and a 70000-byte
// resource fork, and the group, are kept byte for byte by both replacers.
func TestReplaceKeepsGroupAndXattrsDarwinKernel(t *testing.T) {
	group := supplementaryGroup(t)
	fork := make([]byte, 70000)
	for i := range fork {
		fork[i] = byte(i*7 + i>>8)
	}
	finder := make([]byte, 32)
	copy(finder, "TEXTFITS") // type and creator codes
	set := map[string][]byte{
		"com.example.fits.test":  []byte("kept"),
		"com.apple.FinderInfo":   finder,
		"com.apple.quarantine":   []byte("0081;5f000000;fits-test;"),
		"com.apple.ResourceFork": fork,
	}
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			if err := os.Chown(path, -1, group); err != nil {
				t.Fatal(err)
			}
			for _, name := range []string{"com.example.fits.test", "com.apple.FinderInfo", "com.apple.quarantine", "com.apple.ResourceFork"} {
				setxattrDarwin(t, path, name, set[name])
			}
			want := xattrsDarwin(t, path)
			for name, v := range set {
				if name != "com.apple.quarantine" && !bytes.Equal(want[name], v) {
					// (macOS may normalise the quarantine value; the
					// comparison below is with what it reports.)
					t.Fatalf("setup: %s reads back as %d bytes, set %d", name, len(want[name]), len(v))
				}
				if _, ok := want[name]; !ok {
					t.Fatalf("setup: %s not set", name)
				}
			}
			if err := r.run(t, path); err != nil {
				t.Fatal(err)
			}
			if r.name == "FlushResize" && len(mustRead(t, path)) == len(orig) {
				t.Fatal("Flush did not resize")
			}
			var st syscall.Stat_t
			if err := syscall.Stat(path, &st); err != nil {
				t.Fatal(err)
			}
			if int(st.Gid) != group {
				t.Errorf("group after %s = %d, want %d", r.name, st.Gid, group)
			}
			if d := xattrDiff(xattrsDarwin(t, path), want); d != "" {
				t.Errorf("extended attributes after %s: %s", r.name, d)
			}
		})
	}
}

// TestCopyXattrsMakesDestinationEqualDarwinKernel is
// TestCopyXattrsMakesDestinationEqualDarwin with the kernel as the truth:
// copyXattrs removes what only the destination has, gives a shared name
// the source's value and copies an empty attribute.
func TestCopyXattrsMakesDestinationEqualDarwinKernel(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "src")
	dst := filepath.Join(dir, "dst")
	for _, p := range []string{src, dst} {
		if err := os.WriteFile(p, []byte("x"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	setxattrDarwin(t, src, "com.example.both", []byte("from-src"))
	setxattrDarwin(t, dst, "com.example.both", []byte("from-dst"))
	setxattrDarwin(t, dst, "com.example.dst-only", []byte("drop me"))
	setxattrDarwin(t, src, "com.example.empty", []byte{})
	if err := copyXattrs(src, dst); err != nil {
		t.Fatal(err)
	}
	want := xattrsDarwin(t, src)
	if _, ok := want["com.example.empty"]; !ok {
		t.Fatalf("setup: the empty attribute is missing: %v", want)
	}
	if d := xattrDiff(xattrsDarwin(t, dst), want); d != "" {
		t.Fatalf("destination attributes differ from the source's: %s", d)
	}
}
