//go:build unix

package fits

import (
	"fmt"
	"os"
	"syscall"
)

// checkReplaceable refuses a temp-file-and-rename replacement of target
// that would break hard links: the rename gives the path a new inode, so
// the other names would keep the old content (ErrHardLinked).
func checkReplaceable(target string, st os.FileInfo) error {
	if s, ok := st.Sys().(*syscall.Stat_t); ok && s.Nlink > 1 {
		return fmt.Errorf("%w: %s has %d links", ErrHardLinked, target, s.Nlink)
	}
	return nil
}

// chownFile is (*os.File).Chown; tests replace it to simulate a chown the
// process is not allowed to make.
var chownFile = (*os.File).Chown

// carryOver gives tmp (at tmpPath) the owner, group and extended
// attributes (on Linux they hold POSIX ACLs; on macOS the resource fork
// and Finder info, but not the ACL: xattr_darwin.go) of the file it will
// replace. Stat_t differs per system (Nlink is uint16 on darwin, uint64
// on linux/amd64); only Nlink, Uid and Gid are used. If one
// cannot be carried over, it returns ErrCannotPreserve and the caller
// removes the temp file, so the original is never replaced by a file that
// silently lost them.
func carryOver(target string, st os.FileInfo, tmp *os.File, tmpPath string) error {
	if s, ok := st.Sys().(*syscall.Stat_t); ok {
		ts, err := tmp.Stat()
		if err != nil {
			return err
		}
		if t, ok := ts.Sys().(*syscall.Stat_t); ok && (t.Uid != s.Uid || t.Gid != s.Gid) {
			if err := chownFile(tmp, int(s.Uid), int(s.Gid)); err != nil {
				return fmt.Errorf("%w: owner %d:%d of %s: %v", ErrCannotPreserve, s.Uid, s.Gid, target, err)
			}
		}
	}
	return copyXattrs(target, tmpPath)
}
