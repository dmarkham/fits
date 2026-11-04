//go:build interop && linux

package fits

import (
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"syscall"
	"testing"
)

// A file without an extended ACL in a directory with a default ACL: the
// temp file inherits ACL entries from the directory that the original
// never had, and the replacement must not keep them (nor any other
// extended attribute the temp file got that the original lacks). Truth:
// getfacl / getfattr and the mode bits of the file before and after.
func TestReplaceDropsInheritedACL(t *testing.T) {
	requireTools(t, "setfacl", "getfacl", "getfattr")
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			toolOutput(t, "setfacl", "-d", "-m", "u:nobody:rwx", filepath.Dir(path))
			toolOutput(t, "setfacl", "-b", path)
			wantACL := toolOutput(t, "getfacl", "-n", path)
			wantAttrs := toolOutput(t, "getfattr", "-d", "-m", "-", path)
			before, err := os.Stat(path)
			if err != nil {
				t.Fatal(err)
			}
			if strings.Contains(wantACL, "user:65534") {
				t.Fatalf("setup: the file has the directory's ACL entries:\n%s", wantACL)
			}
			if err := r.run(t, path); err != nil {
				t.Fatal(err)
			}
			if r.name == "FlushResize" && len(mustRead(t, path)) == len(orig) {
				t.Fatal("Flush did not resize")
			}
			if got := toolOutput(t, "getfacl", "-n", path); got != wantACL {
				t.Errorf("ACL after %s:\n%s\nwant (the original's)\n%s", r.name, got, wantACL)
			}
			if got := toolOutput(t, "getfattr", "-d", "-m", "-", path); got != wantAttrs {
				t.Errorf("extended attributes after %s:\n%s\nwant\n%s", r.name, got, wantAttrs)
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

func TestReplaceKeepsGroupXattrsAndACL(t *testing.T) {
	requireTools(t, "setfattr", "getfattr", "setfacl", "getfacl")
	// A supplementary group of this user other than the primary one, so the
	// temp file (created with the primary group) needs a chown.
	group := supplementaryGroup(t)
	for _, r := range replacers {
		t.Run(r.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			if err := os.Chown(path, -1, group); err != nil {
				t.Fatal(err)
			}
			toolOutput(t, "setfattr", "-n", "user.fits.test", "-v", "kept", path)
			toolOutput(t, "setfacl", "-m", "u:nobody:r,g:"+strconv.Itoa(group)+":rw", path)
			wantAttrs := toolOutput(t, "getfattr", "-d", "-m", "-", path)
			wantACL := toolOutput(t, "getfacl", "-n", path)
			if !strings.Contains(wantAttrs, "user.fits.test") || !strings.Contains(wantACL, "user:65534:r--") {
				t.Fatalf("setup did not take:\n%s\n%s", wantAttrs, wantACL)
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
			if got := toolOutput(t, "getfattr", "-d", "-m", "-", path); got != wantAttrs {
				t.Errorf("extended attributes after %s:\n%s\nwant\n%s", r.name, got, wantAttrs)
			}
			if got := toolOutput(t, "getfacl", "-n", path); got != wantACL {
				t.Errorf("ACL after %s:\n%s\nwant\n%s", r.name, got, wantACL)
			}
		})
	}
}
