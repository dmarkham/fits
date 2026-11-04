//go:build interop && darwin

package fits

import (
	"encoding/hex"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"
)

// macOS: the resizing Flush and EditFile give the new file the old one's
// group and every listed extended attribute, byte for byte: a user
// attribute, the Finder info, the quarantine flag and a resource fork
// larger than 64 KiB (com.apple.ResourceFork, which is stored outside the
// inline attribute area). Truth: /usr/bin/xattr -lx (every attribute name
// with its value in hex) and stat -f %g before and after.
func TestReplaceKeepsGroupAndXattrsDarwin(t *testing.T) {
	requireTools(t, "xattr", "stat")
	group := supplementaryGroup(t)
	fork := make([]byte, 70000)
	for i := range fork {
		fork[i] = byte(i*7 + i>>8)
	}
	finder := make([]byte, 32)
	copy(finder, "TEXTFITS") // type and creator codes
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			if err := os.Chown(path, -1, group); err != nil {
				t.Fatal(err)
			}
			toolOutput(t, "xattr", "-w", "com.example.fits.test", "kept", path)
			toolOutput(t, "xattr", "-wx", "com.apple.FinderInfo", hex.EncodeToString(finder), path)
			toolOutput(t, "xattr", "-w", "com.apple.quarantine", "0081;5f000000;fits-test;", path)
			toolOutput(t, "xattr", "-wx", "com.apple.ResourceFork", hex.EncodeToString(fork), path)
			wantAttrs := toolOutput(t, "xattr", "-lx", path)
			for _, name := range []string{"com.example.fits.test", "com.apple.FinderInfo", "com.apple.quarantine", "com.apple.ResourceFork"} {
				if !strings.Contains(wantAttrs, name) {
					t.Fatalf("setup: %s not set:\n%s", name, wantAttrs)
				}
			}
			if err := r.run(t, path); err != nil {
				t.Fatal(err)
			}
			if r.name == "FlushResize" && len(mustRead(t, path)) == len(orig) {
				t.Fatal("Flush did not resize")
			}
			if got := strings.TrimSpace(toolOutput(t, "stat", "-f", "%g", path)); got != strconv.Itoa(group) {
				t.Errorf("group after %s = %s, want %d", r.name, got, group)
			}
			if got := toolOutput(t, "xattr", "-lx", path); got != wantAttrs {
				t.Errorf("extended attributes after %s:\n%s\nwant\n%s", r.name, got, wantAttrs)
			}
		})
	}
}

// copyXattrs also removes what the destination has and the source lacks,
// and gives an attribute the source's value when both have it. Truth:
// /usr/bin/xattr -lx of the destination equals that of the source.
func TestCopyXattrsMakesDestinationEqualDarwin(t *testing.T) {
	requireTools(t, "xattr")
	dir := t.TempDir()
	src := filepath.Join(dir, "src")
	dst := filepath.Join(dir, "dst")
	for _, p := range []string{src, dst} {
		if err := os.WriteFile(p, []byte("x"), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	toolOutput(t, "xattr", "-w", "com.example.both", "from-src", src)
	toolOutput(t, "xattr", "-w", "com.example.both", "from-dst", dst)
	toolOutput(t, "xattr", "-w", "com.example.dst-only", "drop me", dst)
	toolOutput(t, "xattr", "-w", "com.example.empty", "", src)
	if err := copyXattrs(src, dst); err != nil {
		t.Fatal(err)
	}
	want := toolOutput(t, "xattr", "-lx", src)
	if got := toolOutput(t, "xattr", "-lx", dst); got != want {
		t.Fatalf("destination attributes:\n%s\nwant (the source's)\n%s", got, want)
	}
	if !strings.Contains(want, "com.example.empty") {
		t.Fatalf("setup: the empty attribute is missing:\n%s", want)
	}
}

// Known gap, pinned so that the docs stay true: an ACL is NOT carried over
// on macOS (xattr_darwin.go). The rewritten file has no explicit ACL
// entries. Truth: ls -le (the ACL listing) before and after.
func TestReplaceDropsACLDarwin(t *testing.T) {
	requireTools(t, "chmod", "ls")
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, _ := copySample(t, "testprog.std")
			// An allow entry: a deny entry such as "deny delete" would make
			// macOS refuse the rename that replaces the file (the rewrite
			// then fails with EACCES, which is not what this test pins).
			toolOutput(t, "chmod", "+a", "everyone allow read", path)
			if acl := toolOutput(t, "ls", "-le", path); !strings.Contains(acl, "allow read") {
				t.Fatalf("setup: no ACL:\n%s", acl)
			}
			if err := r.run(t, path); err != nil {
				t.Fatal(err)
			}
			if acl := toolOutput(t, "ls", "-le", path); strings.Contains(acl, "allow read") {
				t.Fatalf("the ACL was carried over after %s; update xattr_darwin.go, the EditFile/Flush docs and docs/VALIDATION.md:\n%s", r.name, acl)
			}
		})
	}
}
