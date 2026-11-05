//go:build linux

package fits

import (
	"bytes"
	"errors"
	"fmt"
	"syscall"
)

// copyXattrs makes dst's extended attributes (user.*, and the
// system.posix_acl_* attributes that hold ACLs) those of src: it removes
// every attribute dst has that src lacks (such as the ACL a new file
// inherits from a default ACL on its directory), then copies src's.
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
		if keep[name] {
			continue
		}
		if err := syscall.Removexattr(dst, name); err != nil {
			return fmt.Errorf("%w: remove extended attribute %s that %s does not have: %v", ErrCannotPreserve, name, src, err)
		}
	}
	for _, name := range names {
		val, err := getXattr(src, name)
		if err != nil {
			return fmt.Errorf("%w: read extended attribute %s of %s: %v", ErrCannotPreserve, name, src, err)
		}
		if err := syscall.Setxattr(dst, name, val, 0); err != nil {
			return fmt.Errorf("%w: extended attribute %s of %s: %v", ErrCannotPreserve, name, src, err)
		}
	}
	return nil
}

func listXattrs(path string) ([]string, error) {
	for {
		n, err := syscall.Listxattr(path, nil)
		if errors.Is(err, syscall.ENOTSUP) {
			return nil, nil // the file system has none
		}
		if err != nil {
			return nil, err
		}
		if n == 0 {
			return nil, nil
		}
		buf := make([]byte, n)
		m, err := syscall.Listxattr(path, buf)
		if errors.Is(err, syscall.ERANGE) {
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
		n, err := syscall.Getxattr(path, name, nil)
		if err != nil {
			return nil, err
		}
		buf := make([]byte, n)
		m, err := syscall.Getxattr(path, name, buf)
		if errors.Is(err, syscall.ERANGE) {
			continue
		}
		if err != nil {
			return nil, err
		}
		return buf[:m], nil
	}
}
