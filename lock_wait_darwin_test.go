//go:build darwin

package fits

import "testing"

// kernelShowsWaiter: macOS has no /proc/locks or other listing of blocked
// lock requests, so the evidence of a wait is the kernel's EAGAIN to the
// non-waiting request that lockWaitHook reported.
func kernelShowsWaiter(*testing.T, string, string, <-chan struct{}) {}
