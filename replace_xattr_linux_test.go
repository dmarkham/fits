//go:build linux

package fits

import (
	"bytes"
	"encoding/binary"
	"errors"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"syscall"
	"testing"
)

// Tier 1 twins of the getfattr / getfacl tests (replace_preserve_linux_interop_test.go,
// tier 2): the same setups and checks with the kernel's own interface, so
// plain `go test` proves the replacements keep extended attributes and
// POSIX ACLs without the attr / acl programs. The attributes are set with
// setxattr(2) (an ACL in the kernel's posix_acl_xattr form, the bytes
// setfacl writes) and read back with listxattr(2) / getxattr(2), the calls
// getfattr and getfacl make; the truth is what the kernel reports for the
// file before and after.

// POSIX ACL xattr encoding (linux/posix_acl_xattr.h): a little-endian
// version 2 header, then entries {tag u16, perm u16, id u32}, ordered by
// tag and id.
const (
	aclUserObj  = 0x01
	aclUser     = 0x02
	aclGroupObj = 0x04
	aclGroup    = 0x08
	aclMask     = 0x10
	aclOther    = 0x20
	aclNoID     = 0xFFFFFFFF
	aclNobody   = 65534
)

type aclEntry struct {
	tag, perm uint16
	id        uint32
}

func aclXattr(entries ...aclEntry) []byte {
	b := binary.LittleEndian.AppendUint32(nil, 2)
	for _, e := range entries {
		b = binary.LittleEndian.AppendUint16(b, e.tag)
		b = binary.LittleEndian.AppendUint16(b, e.perm)
		b = binary.LittleEndian.AppendUint32(b, e.id)
	}
	return b
}

// setxattr sets an attribute; a file system without (that kind of)
// extended attributes skips the test: a property of the file system the
// temp directory is on, not a missing tool.
func setxattr(t *testing.T, path, name string, value []byte) {
	t.Helper()
	if err := syscall.Setxattr(path, name, value, 0); err != nil {
		if errors.Is(err, syscall.ENOTSUP) || errors.Is(err, syscall.EOPNOTSUPP) {
			t.Skipf("FILESYSTEM SKIP: %s does not support %s on %s: %v", filepath.Dir(path), name, path, err)
		}
		t.Fatalf("setxattr %s %s: %v", path, name, err)
	}
}

// xattrs is every extended attribute of path, name -> value, as
// listxattr(2) and getxattr(2) report them.
func xattrs(t *testing.T, path string) map[string][]byte {
	t.Helper()
	size, err := syscall.Listxattr(path, nil)
	if err != nil {
		t.Fatalf("listxattr %s: %v", path, err)
	}
	buf := make([]byte, size)
	n, err := syscall.Listxattr(path, buf)
	if err != nil {
		t.Fatalf("listxattr %s: %v", path, err)
	}
	out := map[string][]byte{}
	for _, name := range strings.Split(strings.TrimRight(string(buf[:n]), "\x00"), "\x00") {
		if name == "" {
			continue
		}
		vs, err := syscall.Getxattr(path, name, nil)
		if err != nil {
			t.Fatalf("getxattr %s %s: %v", path, name, err)
		}
		v := make([]byte, vs)
		vn, err := syscall.Getxattr(path, name, v)
		if err != nil {
			t.Fatalf("getxattr %s %s: %v", path, name, err)
		}
		out[name] = v[:vn]
	}
	return out
}

func sameXattrs(a, b map[string][]byte) bool {
	if len(a) != len(b) {
		return false
	}
	for k, v := range a {
		if w, ok := b[k]; !ok || !bytes.Equal(v, w) {
			return false
		}
	}
	return true
}

func xattrNames(m map[string][]byte) []string {
	var n []string
	for k := range m {
		n = append(n, k)
	}
	sort.Strings(n)
	return n
}

// TestReplaceDropsInheritedACLKernel is TestReplaceDropsInheritedACL with
// the kernel as the truth: a file without an extended ACL in a directory
// whose default ACL (u:nobody:rwx) a new file inherits; after the
// replacement the file has exactly its old extended attributes (none of
// the directory's ACL) and its old mode.
func TestReplaceDropsInheritedACLKernel(t *testing.T) {
	def := aclXattr(
		aclEntry{aclUserObj, 7, aclNoID}, aclEntry{aclUser, 7, aclNobody},
		aclEntry{aclGroupObj, 5, aclNoID}, aclEntry{aclMask, 7, aclNoID}, aclEntry{aclOther, 5, aclNoID})
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			dir := filepath.Dir(path)
			setxattr(t, dir, "system.posix_acl_default", def)
			// The hazard is real: a file created there now inherits it.
			pf := filepath.Join(dir, "probe.fits")
			if err := os.WriteFile(pf, nil, 0o644); err != nil {
				t.Fatal(err)
			}
			if _, ok := xattrs(t, pf)["system.posix_acl_access"]; !ok {
				t.Fatalf("setup: a new file in a directory with a default ACL has no access ACL (%v)", xattrNames(xattrs(t, pf)))
			}
			if err := os.Remove(pf); err != nil {
				t.Fatal(err)
			}
			want := xattrs(t, path)
			if _, ok := want["system.posix_acl_access"]; ok {
				t.Fatalf("setup: the file already has an access ACL")
			}
			before, err := os.Stat(path)
			if err != nil {
				t.Fatal(err)
			}
			if err := r.run(t, path); err != nil {
				t.Fatal(err)
			}
			if r.name == "FlushResize" && len(mustRead(t, path)) == len(orig) {
				t.Fatal("Flush did not resize")
			}
			if got := xattrs(t, path); !sameXattrs(got, want) {
				t.Errorf("extended attributes after %s: %v, want %v (the original's)", r.name, xattrNames(got), xattrNames(want))
			}
			after, err := os.Stat(path)
			if err != nil {
				t.Fatal(err)
			}
			if after.Mode() != before.Mode() {
				t.Errorf("mode after %s = %v, want %v", r.name, after.Mode(), before.Mode())
			}
		})
	}
}

// TestReplaceKeepsGroupXattrsAndACLKernel is TestReplaceKeepsGroupXattrsAndACL
// with the kernel as the truth: a file in a supplementary group with a
// user attribute and an access ACL (u:nobody:r, g:<group>:rw) keeps the
// group and every attribute byte for byte.
func TestReplaceKeepsGroupXattrsAndACLKernel(t *testing.T) {
	group := supplementaryGroup(t)
	acl := aclXattr(
		aclEntry{aclUserObj, 6, aclNoID}, aclEntry{aclUser, 4, aclNobody},
		aclEntry{aclGroupObj, 4, aclNoID}, aclEntry{aclGroup, 6, uint32(group)},
		aclEntry{aclMask, 6, aclNoID}, aclEntry{aclOther, 4, aclNoID})
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			if err := os.Chown(path, -1, group); err != nil {
				t.Fatal(err)
			}
			setxattr(t, path, "user.fits.test", []byte("kept"))
			setxattr(t, path, "system.posix_acl_access", acl)
			want := xattrs(t, path)
			if string(want["user.fits.test"]) != "kept" || want["system.posix_acl_access"] == nil {
				t.Fatalf("setup did not take: %v", xattrNames(want))
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
			if got := xattrs(t, path); !sameXattrs(got, want) {
				t.Errorf("extended attributes after %s: %v, want %v, byte for byte", r.name, xattrNames(got), xattrNames(want))
			}
		})
	}
}
