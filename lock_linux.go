//go:build linux

package fits

import (
	"errors"
	"os"
	"syscall"
)

// Linux side of the locks described in lock_ofd.go. The edit lock is a
// flock(2); flock and OFD record locks do not interact on Linux, so the
// read/write lock covers the whole file.

// Linux fcntl commands for open file description locks (<fcntl.h>; the
// syscall package does not name them). They are the same on every Linux
// architecture.
const (
	fOFDSetlk  = 37
	fOFDSetlkw = 38
)

// The read/write lock range: the whole file (length 0 = to the end of
// any file).
const (
	rwLockStart = 0
	rwLockLen   = 0
)

// editLockOpenFlag is the open flag for a descriptor that will carry the
// edit lock: flock needs no access mode, so the caller's is kept.
func editLockOpenFlag(flag int) int { return flag }

// lockForEdit takes a non-blocking exclusive flock(2) on f's open file
// description. It is released when f is closed. A lock held by another
// editor (another OpenForEdit or EditFile of the same file, in this or
// another process) makes it fail with ErrLocked.
func lockForEdit(f *os.File) error {
	for {
		err := syscall.Flock(int(f.Fd()), syscall.LOCK_EX|syscall.LOCK_NB)
		if errors.Is(err, syscall.EINTR) {
			continue
		}
		return editLockError(f, err)
	}
}

// noRecordLocks reports whether err says the file system does not support
// record locks.
func noRecordLocks(err error) bool {
	return errors.Is(err, syscall.ENOLCK) || errors.Is(err, syscall.EINVAL) || errors.Is(err, syscall.EOPNOTSUPP)
}
