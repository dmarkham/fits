//go:build !linux && !darwin

package fits

import (
	"fmt"
	"os"
	"runtime"
)

// editSupported refuses editing outside Linux and macOS (Windows, the
// BSDs, ...). The edit lock, the reader/writer lock that keeps an
// in-place write from changing bytes a read-only Open is reading, and
// carrying extended attributes over to a rewritten file are implemented
// for Linux and macOS only (lock_ofd.go); without them an edit could
// change a file under a reader, run beside a second editor, or drop
// metadata, all without an error.
func editSupported() error {
	return fmt.Errorf("%w (GOOS %s)", ErrEditUnsupported, runtime.GOOS)
}

// The functions below are never reached (OpenForEdit and EditFile stop at
// editSupported); a read-only Open takes no lock here.

func editLockOpenFlag(flag int) int { return flag }

func lockForEdit(*os.File) error { return editSupported() }

func lockForRead(*os.File) error { return nil }

func unlockRead(*os.File) {}

func lockForWrite(*os.File, bool) error { return editSupported() }

func unlockWrite(*os.File) {}

// Create takes no lock here (no editor or reader of this package locks a
// file on this system) and just empties the file.

func lockForCreate(*os.File) error { return nil }

func truncateForCreate(f *os.File) error { return f.Truncate(0) }
