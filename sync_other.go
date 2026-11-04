//go:build !darwin

package fits

import "os"

// syncFile makes f's data and metadata durable: (*os.File).Sync, which is
// fsync(2) on Linux. (macOS needs more: sync_darwin.go.)
func syncFile(f *os.File) error { return f.Sync() }
