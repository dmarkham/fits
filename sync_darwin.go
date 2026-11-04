//go:build darwin

package fits

import (
	"errors"
	"os"
	"syscall"
)

// syncFile makes f's data and metadata durable (macOS). fsync(2) on macOS
// only hands the data to the drive, which may keep it in its volatile
// cache; fcntl F_FULLFSYNC also asks the drive to flush that cache to
// stable storage, which is what the journal and the temp-file rename need
// across a power loss. (*os.File).Sync issues F_FULLFSYNC on darwin
// (internal/poll/fd_fsync_darwin.go) and itself falls back to fsync when
// the file system answers ENOTSUP (SMB mounts, golang.org/issue/64215).
// XNU passes F_FULLFSYNC to the file system as an ioctl, and a file system
// that does not implement it may also answer ENOTTY or EINVAL
// (bsd/kern/kern_descrip.c treats all three as "not supported" for
// F_BARRIERFSYNC); for those this falls back to plain fsync too, the best
// that file system offers.
func syncFile(f *os.File) error {
	err := f.Sync()
	if errors.Is(err, syscall.ENOTTY) || errors.Is(err, syscall.EINVAL) {
		return plainFsync(f)
	}
	return err
}

func plainFsync(f *os.File) error {
	rc, err := f.SyscallConn()
	if err != nil {
		return err
	}
	var serr error
	if err := rc.Control(func(fd uintptr) {
		for {
			serr = syscall.Fsync(int(fd))
			if !errors.Is(serr, syscall.EINTR) {
				return
			}
		}
	}); err != nil {
		return err
	}
	return serr
}
