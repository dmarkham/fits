//go:build linux

package fits

import (
	"bufio"
	"fmt"
	"os"
	"strings"
	"syscall"
	"testing"
	"time"
)

// procLocksWaiter reports whether /proc/locks lists a blocked lock request
// ("->") of kind ("READ" or "WRITE") on the file at path: a descriptor of
// this process waiting in F_OFD_SETLKW.
func procLocksWaiter(t *testing.T, path, kind string) bool {
	t.Helper()
	st, err := os.Stat(path)
	if err != nil {
		t.Fatal(err)
	}
	ino := fmt.Sprintf(":%d", st.Sys().(*syscall.Stat_t).Ino)
	fh, err := os.Open("/proc/locks")
	if err != nil {
		t.Fatalf("/proc/locks: %v", err)
	}
	defer fh.Close()
	sc := bufio.NewScanner(fh)
	for sc.Scan() {
		f := strings.Fields(sc.Text())
		if len(f) < 7 || f[1] != "->" {
			continue
		}
		// "N: -> OFDLCK ADVISORY READ -1 MAJ:MIN:INO START END"
		if f[2] == "OFDLCK" && f[4] == kind && strings.HasSuffix(f[6], ino) {
			return true
		}
	}
	return false
}

// kernelShowsWaiter waits until /proc/locks lists the blocked request of
// kind on path (Linux): the call that reported its wait through
// lockWaitHook really is blocked in the kernel. done reporting first
// means the request was granted without blocking after all, which fails
// the test.
func kernelShowsWaiter(t *testing.T, path, kind string, done <-chan struct{}) {
	t.Helper()
	deadline := time.Now().Add(20 * time.Second)
	for time.Now().Before(deadline) {
		if procLocksWaiter(t, path, kind) {
			return
		}
		select {
		case <-done:
			t.Fatalf("the %s request on %s reported a wait but returned without one showing in /proc/locks", kind, path)
		default:
		}
		time.Sleep(time.Millisecond)
	}
	t.Fatalf("no blocked %s lock request on %s in /proc/locks within 20 s", kind, path)
}
