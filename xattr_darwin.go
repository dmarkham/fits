//go:build darwin

package fits

import (
	"bytes"
	"errors"
	"fmt"

	"golang.org/x/sys/unix"
)

// copyXattrs makes dst's extended attributes those of src (macOS): it
// removes every attribute dst has that src lacks, then copies src's. On
// macOS the resource fork (com.apple.ResourceFork), the Finder info
// (com.apple.FinderInfo), the quarantine flag and any user-defined
// attribute are extended attributes and are carried over byte for byte.
// The syscall package has no xattr calls on darwin; golang.org/x/sys/unix
// supplies them. Every call uses XATTR_NOFOLLOW (the L* variants): src and
// dst are the resolved target and the temp file, never links, and a link
// swapped in meanwhile must not redirect the copy.
//
// Not carried over:
//
//   - ACLs. macOS keeps a file's ACL in the protected attribute
//     com.apple.system.Security, which listxattr does not list and only
//     root may read (xnu bsd/vfs/vfs_syscalls.c, xattr_entitlement_check);
//     carrying it needs acl_get_file / acl_set_file, which this package
//     does not call. The rewritten file has no explicit ACL: it has
//     whatever ACL entries it inherits from its directory, like any new
//     file there. Documented in the EditFile and Flush docs.
//   - attributes the kernel assigns and a process cannot set or remove
//     (kernelManagedXattr): the new file gets the kernel's own value.
func copyXattrs(src, dst string) error {
	names, err := listXattrs(src)
	if err != nil {
		return fmt.Errorf("%w: list extended attributes of %s: %v", ErrCannotPreserve, src, err)
	}
	have, err := listXattrs(dst)
	if err != nil {
		return fmt.Errorf("%w: list extended attributes of %s: %v", ErrCannotPreserve, dst, err)
	}
	keep := map[string]bool{}
	for _, name := range names {
		keep[name] = true
	}
	for _, name := range have {
		if keep[name] || kernelManagedXattr[name] {
			continue
		}
		if err := unix.Lremovexattr(dst, name); err != nil {
			return fmt.Errorf("%w: remove extended attribute %s that %s does not have: %v", ErrCannotPreserve, name, src, err)
		}
	}
	for _, name := range names {
		if kernelManagedXattr[name] {
			continue
		}
		val, err := getXattr(src, name)
		if err != nil {
			return fmt.Errorf("%w: read extended attribute %s of %s: %v", ErrCannotPreserve, name, src, err)
		}
		if err := unix.Lsetxattr(dst, name, val, 0); err != nil {
			return fmt.Errorf("%w: extended attribute %s of %s: %v", ErrCannotPreserve, name, src, err)
		}
	}
	return nil
}

// kernelManagedXattr names attributes that macOS sets itself on files a
// process creates and that the process may not change:
// com.apple.provenance (macOS 13 and later) records which app created the
// file.
var kernelManagedXattr = map[string]bool{"com.apple.provenance": true}

func listXattrs(path string) ([]string, error) {
	for {
		n, err := unix.Llistxattr(path, nil)
		if errors.Is(err, unix.ENOTSUP) {
			return nil, nil // the file system has none
		}
		if err != nil {
			return nil, err
		}
		if n == 0 {
			return nil, nil
		}
		buf := make([]byte, n)
		m, err := unix.Llistxattr(path, buf)
		if errors.Is(err, unix.ERANGE) {
			continue // grew in between
		}
		if err != nil {
			return nil, err
		}
		var names []string
		for _, b := range bytes.Split(buf[:m], []byte{0}) {
			if len(b) > 0 {
				names = append(names, string(b))
			}
		}
		return names, nil
	}
}

func getXattr(path, name string) ([]byte, error) {
	for {
		n, err := unix.Lgetxattr(path, name, nil)
		if err != nil {
			return nil, err
		}
		buf := make([]byte, n)
		if n == 0 {
			return buf, nil
		}
		m, err := unix.Lgetxattr(path, name, buf)
		if errors.Is(err, unix.ERANGE) {
			continue
		}
		if err != nil {
			return nil, err
		}
		return buf[:m], nil
	}
}
