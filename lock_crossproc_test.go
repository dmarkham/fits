//go:build linux || darwin

package fits

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"os/exec"
	"strings"
	"testing"
	"time"
)

// The lock guarantees across PROCESSES (the other lock tests run both
// sides in one process, which the locks treat the same way). A child
// process, the test binary itself re-run with FITS_LOCK_HELPER set
// (TestLockHelperProcess), opens, edits and reads the same file on
// command. Truth: the documented error of each refused call
// (ErrLocked, ErrOpenForReading) and byte equality of the file after a
// refusal; for the wait, the child's report that its read lock request
// was refused by the kernel and that it then read the written header.

const lockHelperEnv = "FITS_LOCK_HELPER"

// lockResult names the outcome of a call for the protocol.
func lockResult(err error) string {
	switch {
	case err == nil:
		return "ok"
	case errors.Is(err, ErrLocked):
		return "ErrLocked"
	case errors.Is(err, ErrOpenForReading):
		return "ErrOpenForReading"
	case errors.Is(err, ErrPendingJournal):
		return "ErrPendingJournal"
	default:
		return "error: " + strings.ReplaceAll(err.Error(), "\n", " ")
	}
}

// TestLockHelperProcess is the child side; it does nothing in a normal
// test run. Commands arrive one per line on stdin; every reply is one line
// on stdout starting with "R: ".
//
//	edit           OpenForEdit and keep the *File
//	close-edit     Close it
//	read           Open and keep the *File
//	close-read     Close it
//	editfile       EditFile with a pass-through copy
//	create         Create (truncates the file if it is allowed)
//	flush V        OpenForEdit, set INPLACE = V in HDU 2 (same block
//	               count), Flush, Close
//	read-wait      Open, reporting "R: waiting" when its read lock request
//	               is refused and it waits; then "R: ok INPLACE=<value>"
//	exit
func TestLockHelperProcess(t *testing.T) {
	if os.Getenv(lockHelperEnv) == "" {
		return
	}
	path := os.Getenv(lockHelperEnv)
	reply := func(s string) { fmt.Printf("R: %s\n", s) }
	var editor, reader *File
	sc := bufio.NewScanner(os.Stdin)
	for sc.Scan() {
		cmd, arg, _ := strings.Cut(sc.Text(), " ")
		switch cmd {
		case "edit":
			f, err := OpenForEdit(path)
			if err == nil {
				editor = f
			}
			reply(lockResult(err))
		case "close-edit":
			reply(lockResult(editor.Close()))
			editor = nil
		case "read":
			f, err := Open(path)
			if err == nil {
				reader = f
			}
			reply(lockResult(err))
		case "close-read":
			reply(lockResult(reader.Close()))
			reader = nil
		case "editfile":
			reply(lockResult(EditFile(path, copyAll)))
		case "create":
			f, err := Create(path)
			if err == nil {
				f.Close()
			}
			reply(lockResult(err))
		case "flush":
			f, err := OpenForEdit(path)
			if err != nil {
				reply(lockResult(err))
				continue
			}
			h, err := f.HDU(2)
			if err == nil {
				err = h.Header().Set("INPLACE", arg, "in-place edit")
			}
			if err == nil {
				err = f.Flush()
			}
			f.Close()
			reply(lockResult(err))
		case "read-wait":
			lockWaitHook = func(kind, _ string) {
				if kind == "READ" {
					reply("waiting")
				}
			}
			f, err := Open(path)
			lockWaitHook = nil
			if err != nil {
				reply(lockResult(err))
				continue
			}
			h, err := f.HDU(2)
			v := ""
			if err == nil {
				v, _ = h.Header().String("INPLACE")
			}
			f.Close()
			reply("ok INPLACE=" + v)
		case "exit":
			return
		default:
			reply("unknown command " + cmd)
		}
	}
}

// lockPeer is the parent's handle on a helper process.
type lockPeer struct {
	t     *testing.T
	cmd   *exec.Cmd
	in    io.WriteCloser
	lines chan string
}

func startLockPeer(t *testing.T, path string) *lockPeer {
	t.Helper()
	cmd := exec.Command(os.Args[0], "-test.run=^TestLockHelperProcess$", "-test.count=1")
	cmd.Env = append(os.Environ(), lockHelperEnv+"="+path)
	cmd.Stderr = os.Stderr
	in, err := cmd.StdinPipe()
	if err != nil {
		t.Fatal(err)
	}
	out, err := cmd.StdoutPipe()
	if err != nil {
		t.Fatal(err)
	}
	if err := cmd.Start(); err != nil {
		t.Fatal(err)
	}
	p := &lockPeer{t: t, cmd: cmd, in: in, lines: make(chan string, 16)}
	go func() {
		sc := bufio.NewScanner(out)
		for sc.Scan() {
			if s, ok := strings.CutPrefix(sc.Text(), "R: "); ok {
				p.lines <- s
			}
		}
		close(p.lines)
	}()
	t.Cleanup(func() {
		fmt.Fprintln(in, "exit")
		in.Close()
		if err := cmd.Wait(); err != nil {
			t.Errorf("helper process: %v", err)
		}
	})
	return p
}

// send writes one command without waiting for its reply.
func (p *lockPeer) send(command string) {
	p.t.Helper()
	if _, err := fmt.Fprintln(p.in, command); err != nil {
		p.t.Fatalf("helper process: send %q: %v", command, err)
	}
}

// next returns the helper's next reply.
func (p *lockPeer) next() string {
	p.t.Helper()
	select {
	case s, ok := <-p.lines:
		if !ok {
			p.t.Fatal("helper process exited")
		}
		return s
	case <-time.After(20 * time.Second):
		p.t.Fatal("helper process: no reply within 20 s")
	}
	return ""
}

// do runs command in the helper and checks its reply.
func (p *lockPeer) do(command, want string) {
	p.t.Helper()
	p.send(command)
	if got := p.next(); got != want {
		p.t.Fatalf("helper %q: %s, want %s", command, got, want)
	}
}

func TestLocksAcrossProcesses(t *testing.T) {
	// One editor at a time, both ways round; the lock moves with the file
	// across a resizing Flush.
	t.Run("editor-vs-editor", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		peer := startLockPeer(t, path)
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		peer.do("edit", "ErrLocked")
		peer.do("editfile", "ErrLocked")
		peer.do("create", "ErrLocked")
		if !bytes.Equal(mustRead(t, path), orig) {
			t.Fatal("a refused editor in the other process changed the file")
		}
		growHDU1(t, f)
		if err := f.Flush(); err != nil {
			t.Fatal(err)
		}
		peer.do("edit", "ErrLocked")
		resized := mustRead(t, path)
		f.Close()

		peer.do("edit", "ok")
		if _, err := OpenForEdit(path); !errors.Is(err, ErrLocked) {
			t.Fatalf("OpenForEdit while the other process edits: %v, want ErrLocked", err)
		}
		if err := EditFile(path, copyAll); !errors.Is(err, ErrLocked) {
			t.Fatalf("EditFile while the other process edits: %v, want ErrLocked", err)
		}
		if _, err := Create(path); !errors.Is(err, ErrLocked) {
			t.Fatalf("Create while the other process edits: %v, want ErrLocked", err)
		}
		if !bytes.Equal(mustRead(t, path), resized) {
			t.Fatal("a refused editor changed the file")
		}
		peer.do("close-edit", "ok")
		g, err := OpenForEdit(path)
		if err != nil {
			t.Fatalf("OpenForEdit after the other process closed its editor: %v", err)
		}
		g.Close()
	})

	// A reader in one process refuses the in-place writes of an editor in
	// the other, without keeping the editor from opening; both ways round.
	t.Run("reader-vs-writer", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		peer := startLockPeer(t, path)
		peer.do("read", "ok")
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatalf("OpenForEdit while the other process reads: %v", err)
		}
		inPlaceHDU2(t, f)
		if err := f.Flush(); !errors.Is(err, ErrOpenForReading) {
			t.Fatalf("in-place Flush while the other process reads: %v, want ErrOpenForReading", err)
		}
		f.Close()
		if _, err := Create(path); !errors.Is(err, ErrOpenForReading) {
			t.Fatalf("Create while the other process reads: %v, want ErrOpenForReading", err)
		}
		if !bytes.Equal(mustRead(t, path), orig) {
			t.Fatal("a refused write changed the file")
		}
		peer.do("close-read", "ok")
		f, err = OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		inPlaceHDU2(t, f)
		if err := f.Flush(); err != nil {
			t.Fatalf("in-place Flush after the other process closed its reader: %v", err)
		}
		f.Close()

		r, err := Open(path)
		if err != nil {
			t.Fatal(err)
		}
		before := mustRead(t, path)
		peer.do("flush child-1", "ErrOpenForReading")
		if !bytes.Equal(mustRead(t, path), before) {
			t.Fatal("the other process's refused Flush changed the file")
		}
		r.Close()
		peer.do("flush child-2", "ok")
		checkInplace(t, path, "child-2")
	})

	// A reader that opens while an in-place write of the other process
	// runs waits for it, then reads the new header.
	t.Run("reader-waits-for-writer", func(t *testing.T) {
		path, _ := copySample(t, "testprog.std")
		peer := startLockPeer(t, path)
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		inPlaceHDU2(t, f)
		flushFaultHook = func(step string) error {
			if step == "header-written" {
				peer.send("read-wait")
				if got := peer.next(); got != "waiting" {
					t.Errorf("the other process's Open during the write: %s, want it to wait", got)
				}
			}
			return nil
		}
		err = f.Flush()
		flushFaultHook = nil
		if err != nil {
			t.Fatal(err)
		}
		if got := peer.next(); got != "ok INPLACE=same block count" {
			t.Fatalf("the other process's Open after the write: %s, want the written header", got)
		}
	})

	// Closing another descriptor of the file in this process (an
	// os.ReadFile, the caller's *os.File after OpenReader) releases none of
	// the locks this process holds through its other descriptors. Classic
	// POSIX record locks would all be dropped by such a close.
	t.Run("other-close-keeps-locks", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		peer := startLockPeer(t, path)
		closeOthers := func() {
			t.Helper()
			mustRead(t, path)
			fh, err := os.Open(path)
			if err != nil {
				t.Fatal(err)
			}
			r, err := OpenReader(fh)
			if err != nil {
				t.Fatal(err)
			}
			r.Close()
			fh.Close()
		}
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		closeOthers()
		peer.do("edit", "ErrLocked")
		f.Close()

		r, err := Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		closeOthers()
		peer.do("flush child", "ErrOpenForReading")
		if !bytes.Equal(mustRead(t, path), orig) {
			t.Fatal("the other process's refused Flush changed the file")
		}
	})
}

// checkInplace checks that HDU 2 of the file at path has INPLACE = want.
func checkInplace(t *testing.T, path, want string) {
	t.Helper()
	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	h, err := r.HDU(2)
	if err != nil {
		t.Fatal(err)
	}
	if v, _ := h.Header().String("INPLACE"); v != want {
		t.Fatalf("INPLACE = %q, want %q", v, want)
	}
}
