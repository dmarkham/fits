//go:build darwin

package fits

import (
	"errors"
	"os"
	"syscall"
)

// macOS side of the locks described in lock_ofd.go.
//
// Open file description locks. XNU implements F_OFD_SETLK / F_OFD_SETLKW
// (bsd/kern/kern_descrip.c, fcntl_nocancel) with the same meaning as on
// Linux: the lock is owned by the open file description, so two
// descriptors from two os.Open calls conflict even in one process, and
// closing one descriptor does not release locks held through another. The
// kernel and libc's fcntl wrapper (libsyscall fcntl-base.c, which passes
// the struct flock pointer through for these commands) have them since
// macOS 11 (xnu-7195), the oldest macOS Go supports; the SDK header
// <sys/fcntl.h> marked them private until macOS 14 (xnu-10002), which is
// why the syscall package does not name them. Classic POSIX locks
// (F_SETLK) are not used: they belong to the process and closing ANY
// descriptor of the file in the process (an os.ReadFile of it, or the
// caller closing the *os.File it gave OpenReader) silently drops every
// one of them.
//
// No flock. On XNU flock(2) and fcntl record locks share one lock list per
// vnode and lf_getblock (bsd/kern/kern_lockf.c) treats any overlapping
// lock of another owner as conflicting, whatever its flavour; flock always
// covers the whole file. A flock edit lock would therefore conflict with
// every reader's shared record lock: an open reader would stop an editor
// from opening (ErrLocked) and a reader would wait for the editor's
// Close. So the edit lock is an exclusive OFD record lock too, on the one
// byte at editLockOffset, and the read/write lock covers
// [0, editLockOffset); the two ranges never overlap, so they never
// conflict, exactly like flock and OFD locks on Linux.
//
// A record lock for writing needs a descriptor open for writing
// (EBADF otherwise), so the descriptor EditFile holds the edit lock on is
// opened read-write here (editLockOpenFlag): on macOS EditFile needs write
// permission on the file, not only on its directory.
//
// A file system without record locks (fcntl fails with ENOLCK, EINVAL,
// ENOTSUP or EOPNOTSUPP; some network file systems): as on Linux a
// read-only Open goes ahead unlocked and every in-place write is refused.
// Unlike Linux, where the edit lock is a flock, the edit lock fails too,
// so OpenForEdit, EditFile and Create refuse such a file with the fcntl
// error instead of running without one editor at a time.

// XNU fcntl commands (<sys/fcntl.h>).
const (
	fOFDSetlk  = 90
	fOFDSetlkw = 91
)

// editLockOffset is the byte the edit lock covers: 2^62, far past the end
// of any FITS file (record locks may lie beyond the end of a file).
const editLockOffset = 1 << 62

// The read/write lock range: everything before the edit lock byte.
const (
	rwLockStart = 0
	rwLockLen   = editLockOffset
)

// editLockOpenFlag makes a descriptor that will carry the edit lock
// read-write (see above).
func editLockOpenFlag(flag int) int {
	return flag&^(os.O_RDONLY|os.O_WRONLY|os.O_RDWR) | os.O_RDWR
}

// lockForEdit takes a non-blocking exclusive OFD lock on the edit lock
// byte of f's open file description. It is released when the last
// descriptor of that description is closed. A lock held by another editor
// (another OpenForEdit or EditFile of the same file, in this or another
// process) makes it fail with ErrLocked. It does not go through ofdLock,
// which tests replace to simulate a file system without record locks for
// the read/write lock.
func lockForEdit(f *os.File) error {
	return editLockError(f, fcntlOFDRange(f, syscall.F_WRLCK, fOFDSetlk, editLockOffset, 1))
}

// noRecordLocks reports whether err says the file system does not support
// record locks. On macOS ENOTSUP (45) and EOPNOTSUPP (102) are distinct.
func noRecordLocks(err error) bool {
	return errors.Is(err, syscall.ENOLCK) || errors.Is(err, syscall.EINVAL) ||
		errors.Is(err, syscall.ENOTSUP) || errors.Is(err, syscall.EOPNOTSUPP)
}
