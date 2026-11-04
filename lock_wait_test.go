//go:build linux || darwin

package fits

import (
	"os"
	"sync"
	"testing"
	"time"
)

// lockWaitWatch records the waits lockWaitHook reports while a test runs.
type lockWaitWatch struct {
	mu     sync.Mutex
	events []lockWaitEvent
}

type lockWaitEvent struct {
	kind string
	st   os.FileInfo // the file the waiting descriptor names
}

// watchLockWaits installs lockWaitHook for the rest of the test. Install
// it before starting the call that is expected to wait.
func watchLockWaits(t *testing.T) *lockWaitWatch {
	t.Helper()
	w := &lockWaitWatch{}
	lockWaitHook = func(kind, name string) {
		st, err := os.Stat(name)
		if err != nil {
			return
		}
		w.mu.Lock()
		w.events = append(w.events, lockWaitEvent{kind, st})
		w.mu.Unlock()
	}
	t.Cleanup(func() { lockWaitHook = nil })
	return w
}

func (w *lockWaitWatch) seen(path, kind string) bool {
	st, err := os.Stat(path)
	if err != nil {
		return false
	}
	w.mu.Lock()
	defer w.mu.Unlock()
	for _, e := range w.events {
		if e.kind == kind && os.SameFile(e.st, st) {
			return true
		}
	}
	return false
}

// await waits until a lock request of kind ("READ" or "WRITE") on the
// file at path found the lock held and went to wait for it, or done
// reports a result first (the call did not wait). It returns false in the
// second case. The evidence is the kernel's refusal of the non-waiting
// request (lockWaitHook) and, on Linux, the blocked request in /proc/locks
// (kernelShowsWaiter).
func (w *lockWaitWatch) await(t *testing.T, path, kind string, done <-chan struct{}) bool {
	t.Helper()
	deadline := time.Now().Add(20 * time.Second)
	for time.Now().Before(deadline) {
		select {
		case <-done:
			return false
		default:
		}
		if w.seen(path, kind) {
			kernelShowsWaiter(t, path, kind, done)
			return true
		}
		time.Sleep(time.Millisecond)
	}
	t.Fatalf("no %s lock request on %s waited within 20 s", kind, path)
	return false
}
