//go:build darwin

package fits

import (
	"os"
	"syscall"
	"testing"
)

// Design check behind lock_darwin.go: on XNU, flock(2) and fcntl record
// locks share one lock list, so an exclusive flock conflicts with a shared
// OFD record lock held through another descriptor (which is why the edit
// lock is not a flock on macOS). Logged, not asserted: the edit lock does
// not depend on the answer, and TestLocksAcrossProcesses proves the
// guarantee itself (a reader does not stop an editor from opening).
func TestDarwinFlockConflictsWithRecordLock(t *testing.T) {
	path, _ := copySample(t, "testprog.std")
	r, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	if err := fcntlOFDLock(r, syscall.F_RDLCK, fOFDSetlk); err != nil {
		t.Fatal(err)
	}
	e, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer e.Close()
	err = syscall.Flock(int(e.Fd()), syscall.LOCK_EX|syscall.LOCK_NB)
	t.Logf("flock(LOCK_EX|LOCK_NB) beside a shared OFD record lock: %v (the XNU source predicts %v)", err, syscall.EWOULDBLOCK)
}
