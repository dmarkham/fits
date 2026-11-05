package fits

import (
	"fmt"
	"os"
)

// lockOpenHook, when non-nil, runs in openLocked between the open and the
// edit lock. Tests use it to replace the file at that moment.
var lockOpenHook func()

// openLocked opens name with flag and takes the edit lock on it
// (lockForEdit). A descriptor opened just before another editor's resizing
// Flush or EditFile renamed a new file over name refers to the replaced
// inode, whose lock that editor no longer holds: locking it would succeed
// while the file at name is still being edited (or let this editor write
// to a deleted inode). So after the lock is taken, the descriptor must
// still be the file name refers to; if it is not, the open is retried on
// the new file (TestEditLockFollowsReplacedFile). On macOS the descriptor
// is opened read-write whatever flag asks (editLockOpenFlag,
// lock_darwin.go).
func openLocked(name string, flag int) (*os.File, error) {
	if err := editSupported(); err != nil {
		return nil, err
	}
	return openWithLock(name, editLockOpenFlag(flag), 0, lockForEdit)
}

// openWithLock is openLocked with the lock function and the permission
// bits for a file that flag creates.
func openWithLock(name string, flag int, perm os.FileMode, lock func(*os.File) error) (*os.File, error) {
	const attempts = 100
	for range attempts {
		fh, err := os.OpenFile(name, flag, perm)
		if err != nil {
			return nil, err
		}
		if lockOpenHook != nil {
			lockOpenHook()
		}
		if err := lock(fh); err != nil {
			fh.Close()
			return nil, err
		}
		a, err := fh.Stat()
		if err != nil {
			fh.Close()
			return nil, err
		}
		b, err := os.Stat(name)
		if err != nil {
			fh.Close()
			return nil, err
		}
		if os.SameFile(a, b) {
			return fh, nil
		}
		fh.Close() // replaced between the open and the lock: try the new file
	}
	return nil, fmt.Errorf("fits: %s was replaced %d times while it was being locked for editing", name, attempts)
}
