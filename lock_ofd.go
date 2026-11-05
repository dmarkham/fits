//go:build linux || darwin

package fits

import (
	"errors"
	"fmt"
	"io"
	"os"
	"syscall"
)

// Two independent locks guard a file this package edits (Linux and macOS,
// local file systems):
//
//   - the EDIT LOCK, exclusive, taken without waiting by OpenForEdit,
//     EditFile and Create for the whole edit, so there is one editor at a
//     time (ErrLocked). Linux: flock(2). macOS: an OFD record lock on one
//     byte far past any FITS data (lock_darwin.go says why it cannot be
//     flock there).
//   - the READ/WRITE LOCK, an open file description (OFD) record lock
//     (fcntl F_OFD_SETLK / F_OFD_SETLKW) over the file's data range
//     (rwLockStart, rwLockLen). A read-only Open (and OpenReader of an
//     *os.File) holds it shared until Close; every write this package
//     makes inside the file (in-place Flush, OverwritePixels, an append in
//     ModeEdit, journal rollback) holds it exclusive while it runs. A
//     write is refused with ErrOpenForReading while a reader holds the
//     file, and a reader that opens while a write runs waits for it to
//     finish.
//
// Every lock belongs to an open file description (one os.Open), not to the
// process, so two *os.File of the same file in ONE process conflict exactly
// as they would in two processes, and closing one descriptor of the file
// never releases a lock held through another. (Classic POSIX fcntl locks,
// F_SETLK, would break both: they belong to the process, and closing any
// descriptor of the file drops all of them.) A lock is released by its
// unlock call or when the last descriptor of its open file description is
// closed, including by the death of the process.
//
// The edit lock and the read/write lock never conflict with each other, so
// a reader does not keep an editor from opening, only from changing bytes
// in place. All of them are advisory: programs that do not use this
// package are not stopped.

// editSupported reports whether OpenForEdit and EditFile can run here.
func editSupported() error { return nil }

// ofdLock is the fcntl call on the read/write lock range; tests replace it
// to simulate a file system without record locks.
var ofdLock = fcntlOFDLock

func fcntlOFDLock(f *os.File, typ int16, cmd int) error {
	return fcntlOFDRange(f, typ, cmd, rwLockStart, rwLockLen)
}

// fcntlOFDRange sets or clears (typ F_UNLCK) an OFD lock on bytes
// [start, start+n) of f (n == 0: to the end of any file), retrying on
// EINTR.
func fcntlOFDRange(f *os.File, typ int16, cmd int, start, n int64) error {
	lk := syscall.Flock_t{Type: typ, Whence: io.SeekStart, Start: start, Len: n}
	for {
		err := syscall.FcntlFlock(f.Fd(), cmd, &lk)
		if errors.Is(err, syscall.EINTR) {
			continue
		}
		return err
	}
}

// lockWaitHook, when non-nil, is called with "READ" or "WRITE" and the
// file's name when a waiting lock request (lockForRead, lockForWrite with
// wait) finds the lock held by another open file description and is about
// to block. Tests use it to see that a call waited.
var lockWaitHook func(kind, name string)

// ofdLockWait takes the lock, waiting for it while another open file
// description holds a conflicting one. It first tries without waiting so
// a wait is visible to lockWaitHook; the result is the same as a single
// F_OFD_SETLKW.
func ofdLockWait(f *os.File, typ int16) error {
	err := ofdLock(f, typ, fOFDSetlk)
	if !isLockConflict(err) {
		return err
	}
	if lockWaitHook != nil {
		kind := "READ"
		if typ == syscall.F_WRLCK {
			kind = "WRITE"
		}
		lockWaitHook(kind, f.Name())
	}
	return ofdLock(f, typ, fOFDSetlkw)
}

// isLockConflict reports whether a non-waiting lock request failed because
// another open file description holds a conflicting lock.
func isLockConflict(err error) bool {
	return errors.Is(err, syscall.EAGAIN) || errors.Is(err, syscall.EACCES)
}

// lockForRead takes the shared read lock on f, waiting while a write of
// this package holds the file. It is released by unlockRead or when f is
// closed. On a file system without record locks the read goes ahead
// unlocked: lockForWrite fails there too, so no in-place write of this
// package can run on the file.
func lockForRead(f *os.File) error {
	err := ofdLockWait(f, syscall.F_RDLCK)
	switch {
	case err == nil:
		return nil
	case noRecordLocks(err):
		return nil
	default:
		return fmt.Errorf("fits: lock %s for reading: %w", f.Name(), err)
	}
}

// unlockRead releases the lock lockForRead took.
func unlockRead(f *os.File) { _ = ofdLock(f, syscall.F_UNLCK, fOFDSetlk) }

// lockForWrite takes the exclusive write lock on f without waiting: while
// a reader holds the file it fails with ErrOpenForReading. wait makes it
// wait instead (journal rollback, which no reader can be using: a
// read-only Open refuses while a journal is pending).
func lockForWrite(f *os.File, wait bool) error {
	var err error
	if wait {
		err = ofdLockWait(f, syscall.F_WRLCK)
	} else {
		err = ofdLock(f, syscall.F_WRLCK, fOFDSetlk)
	}
	switch {
	case err == nil:
		return nil
	case isLockConflict(err):
		return fmt.Errorf("%w (%s)", ErrOpenForReading, f.Name())
	default:
		// Includes a file system without record locks: an in-place write
		// is refused there, since readers cannot be excluded.
		return fmt.Errorf("fits: lock %s for writing (needed to keep readers out of an in-place write): %w", f.Name(), err)
	}
}

// unlockWrite releases the lock lockForWrite took.
func unlockWrite(f *os.File) { _ = ofdLock(f, syscall.F_UNLCK, fOFDSetlk) }

// lockForCreate is the lock Create holds until Close: the edit lock, so
// no editor opens the file while it is being written and Create does not
// rewrite a file an editor has.
func lockForCreate(f *os.File) error { return lockForEdit(f) }

// truncateForCreate empties f for Create, holding the write lock while it
// does so it never cuts a file a reader of this package has open
// (ErrOpenForReading). On a file system without record locks no reader
// can have locked it, and the file is truncated.
func truncateForCreate(f *os.File) error {
	err := lockForWrite(f, false)
	switch {
	case err == nil:
		defer unlockWrite(f)
	case noRecordLocks(err):
	default:
		return err
	}
	return f.Truncate(0)
}

// editLockError maps the result of the non-waiting edit lock request.
func editLockError(f *os.File, err error) error {
	switch {
	case err == nil:
		return nil
	case isLockConflict(err), errors.Is(err, syscall.EWOULDBLOCK):
		return fmt.Errorf("%w (%s)", ErrLocked, f.Name())
	default:
		return fmt.Errorf("fits: lock %s for editing: %w", f.Name(), err)
	}
}
