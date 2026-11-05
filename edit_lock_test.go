//go:build linux || darwin

package fits

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"syscall"
	"testing"
)

// Two editors of one file are refused: while an OpenForEdit *File is open
// (also after a resizing Flush replaced the file) or an EditFile runs, a
// second OpenForEdit or EditFile of the same path fails with ErrLocked and
// changes nothing, and in particular does not roll back the first
// editor's pending journal. The edit lock (flock(2) on Linux, an OFD
// record lock on macOS) belongs to an open file description, so a second
// open in this process conflicts like one in another process would
// (TestLocksAcrossProcesses runs the same refusals across processes).
func TestEditLockRefusesConcurrentEditors(t *testing.T) {
	second := []struct {
		name string
		run  func(path string) error
	}{
		{"OpenForEdit", func(p string) error {
			f, err := OpenForEdit(p)
			if err == nil {
				f.Close()
			}
			return err
		}},
		{"EditFile", func(p string) error { return EditFile(p, copyAll) }},
	}
	for _, s := range second {
		t.Run("open/"+s.name, func(t *testing.T) {
			path, orig := copySample(t, "testprog.std")
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			if err := s.run(path); !errors.Is(err, ErrLocked) {
				t.Fatalf("%s while an OpenForEdit is open: %v, want ErrLocked", s.name, err)
			}
			// After a resizing Flush the file is a new inode; the lock moved
			// with it.
			growHDU1(t, f)
			if err := f.Flush(); err != nil {
				t.Fatal(err)
			}
			if err := s.run(path); !errors.Is(err, ErrLocked) {
				t.Fatalf("%s after the first editor's resizing Flush: %v, want ErrLocked", s.name, err)
			}
			f.Close()
			if err := s.run(path); err != nil {
				t.Fatalf("%s after the first editor closed: %v", s.name, err)
			}
			_ = orig
		})
		t.Run("editfile/"+s.name, func(t *testing.T) {
			path, _ := copySample(t, "testprog.std")
			var inner error
			if err := EditFile(path, func(in *File, out *Writer) error {
				inner = s.run(path)
				return copyAll(in, out)
			}); err != nil {
				t.Fatal(err)
			}
			if !errors.Is(inner, ErrLocked) {
				t.Fatalf("%s during an EditFile: %v, want ErrLocked", s.name, inner)
			}
		})
	}

	// The first editor's journal is not rolled back by a refused second
	// editor: a crash-like stop after the header write leaves a journal;
	// while the first *File is open, a second OpenForEdit is refused and
	// the new header stays on disk.
	t.Run("journal-kept", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		inPlaceHDU2(t, f)
		flushFaultHook = func(s string) error {
			if s == "header-written" {
				return errInjectedCrash
			}
			return nil
		}
		err = f.Flush()
		flushFaultHook = nil
		if !errors.Is(err, errInjectedCrash) {
			t.Fatal(err)
		}
		mid := mustRead(t, path)
		if bytes.Equal(mid, orig) {
			t.Fatal("stop point is before the header write")
		}
		if _, err := OpenForEdit(path); !errors.Is(err, ErrLocked) {
			t.Fatalf("second OpenForEdit: %v, want ErrLocked", err)
		}
		if !bytes.Equal(mustRead(t, path), mid) {
			t.Fatal("the refused second editor rolled the first editor's journal back")
		}
	})
}

// A descriptor opened just before another editor's rename refers to the
// replaced inode, whose lock nobody holds any more. OpenForEdit and
// EditFile must not take the edit lock on it: after the flock they check
// that the path still names the locked inode and otherwise open the new
// file. lockOpenHook runs the other editor between the open and the flock.
// Truth: the edit lock of the editor still running (ErrLocked), and the
// bytes at the path after the second editor's Flush.
func TestEditLockFollowsReplacedFile(t *testing.T) {
	defer func() { lockOpenHook = nil }()
	second := []struct {
		name string
		run  func(path string) error
	}{
		{"OpenForEdit", func(p string) error {
			f, err := OpenForEdit(p)
			if err == nil {
				f.Close()
			}
			return err
		}},
		{"EditFile", func(p string) error { return EditFile(p, copyAll) }},
	}
	for _, s := range second {
		// The first editor resizes the file with Flush and is still open.
		t.Run("resized-and-open/"+s.name, func(t *testing.T) {
			path, _ := copySample(t, "testprog.std")
			var first *File
			lockOpenHook = func() {
				lockOpenHook = nil
				f, err := OpenForEdit(path)
				if err != nil {
					t.Fatal(err)
				}
				growHDU1(t, f)
				if err := f.Flush(); err != nil {
					t.Fatal(err)
				}
				first = f
			}
			err := s.run(path)
			if first == nil {
				t.Fatal("hook did not run")
			}
			defer first.Close()
			if !errors.Is(err, ErrLocked) {
				t.Fatalf("%s opened before the first editor's resizing Flush and locked after it: %v, want ErrLocked", s.name, err)
			}
		})
		// An EditFile has replaced the file and returned; a third editor
		// now holds the new file.
		t.Run("editfile-then-editor/"+s.name, func(t *testing.T) {
			path, _ := copySample(t, "testprog.std")
			var third *File
			lockOpenHook = func() {
				lockOpenHook = nil
				if err := EditFile(path, copyAll); err != nil {
					t.Fatal(err)
				}
				f, err := OpenForEdit(path)
				if err != nil {
					t.Fatal(err)
				}
				third = f
			}
			err := s.run(path)
			if third == nil {
				t.Fatal("hook did not run")
			}
			defer third.Close()
			if !errors.Is(err, ErrLocked) {
				t.Fatalf("%s locked the replaced inode while another editor holds the file: %v, want ErrLocked", s.name, err)
			}
		})
	}
	// With nobody holding the new file, the second editor opens it and its
	// edits reach the path (not a deleted inode).
	t.Run("retry-edits-the-file", func(t *testing.T) {
		path, _ := copySample(t, "testprog.std")
		lockOpenHook = func() {
			lockOpenHook = nil
			if err := EditFile(path, copyAll); err != nil {
				t.Fatal(err)
			}
		}
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		inPlaceHDU2(t, f)
		if err := f.Flush(); err != nil {
			t.Fatal(err)
		}
		f.Close()
		r, err := Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		h, err := r.HDU(2)
		if err != nil {
			t.Fatal(err)
		}
		if v, _ := h.Header().String("INPLACE"); v != "same block count" {
			t.Fatalf("the flushed edit is not in the file at the path (INPLACE = %q): it went to the replaced inode", v)
		}
	})
}

// An in-place write (in-place Flush, OverwritePixels, the appends) while
// the file is open for reading through this package is refused with
// ErrOpenForReading and writes nothing, so a reader never sees a header or
// data unit that is partly old and partly new; once the readers are
// closed the same call succeeds and gives the fault-free result. A
// resizing Flush replaces the file through a rename, which leaves an open
// reader on the old, complete inode, so it is allowed. Truth: byte
// equality with the original file and with the fault-free result.
func TestInPlaceWriteRefusedWhileReaderOpen(t *testing.T) {
	flushOp := inPlaceOp{name: "Flush", sample: "testprog.std",
		run: func(t *testing.T, f *File) error {
			inPlaceHDU2(t, f)
			return f.Flush()
		}}
	openers := []struct {
		name string
		open func(t *testing.T, path string) (*File, func())
	}{
		{"Open", func(t *testing.T, path string) (*File, func()) {
			r, err := Open(path)
			if err != nil {
				t.Fatal(err)
			}
			return r, func() { r.Close() }
		}},
		{"OpenReader", func(t *testing.T, path string) (*File, func()) {
			fh, err := os.Open(path)
			if err != nil {
				t.Fatal(err)
			}
			r, err := OpenReader(fh)
			if err != nil {
				t.Fatal(err)
			}
			// Close releases the lock and leaves the caller's file open.
			return r, func() {
				r.Close()
				if _, err := fh.Stat(); err != nil {
					t.Errorf("OpenReader's Close closed the caller's file: %v", err)
				}
				fh.Close()
			}
		}},
	}
	for _, op := range append([]inPlaceOp{flushOp}, inPlaceOps...) {
		for _, o := range openers {
			t.Run(op.name+"/"+o.name, func(t *testing.T) {
				_, ppath, _ := setupOp(t, op)
				if err := runOp(t, op, ppath); err != nil {
					t.Fatal(err)
				}
				post := mustRead(t, ppath)

				dir, path, pre := setupOp(t, op)
				r, closeReader := o.open(t, path)
				f, err := OpenForEdit(path)
				if err != nil {
					t.Fatalf("OpenForEdit while a reader is open: %v", err)
				}
				defer f.Close()
				nhdu := f.NumHDU()
				if err := op.run(t, f); !errors.Is(err, ErrOpenForReading) {
					t.Fatalf("%s while the file is open for reading: %v, want ErrOpenForReading", op.name, err)
				}
				if got := mustRead(t, path); !bytes.Equal(got, pre) {
					t.Fatalf("refused %s changed the file (first diff at %d)", op.name, firstDiff(got, pre))
				}
				if s := sidecars(t, dir, op.file()); len(s) != 0 {
					t.Fatalf("refused %s left %v", op.name, s)
				}
				if f.NumHDU() != nhdu {
					t.Fatalf("refused %s changed NumHDU %d -> %d", op.name, nhdu, f.NumHDU())
				}
				if _, err := r.HDU(r.NumHDU() - 1); err != nil {
					t.Fatalf("reader: %v", err)
				}
				closeReader()
				if err := op.run(t, f); err != nil {
					t.Fatalf("%s after the reader closed: %v", op.name, err)
				}
				if err := f.Close(); err != nil {
					t.Fatal(err)
				}
				if got := mustRead(t, path); !bytes.Equal(got, post) {
					t.Fatalf("%s after the reader closed differs from the fault-free result (first diff at %d)", op.name, firstDiff(got, post))
				}
			})
		}
	}

	t.Run("resize-leaves-reader-on-old-file", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		r, err := Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer r.Close()
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		growHDU1(t, f)
		if err := f.Flush(); err != nil {
			t.Fatalf("resizing Flush with a reader open: %v", err)
		}
		f.Close()
		if bytes.Equal(mustRead(t, path), orig) {
			t.Fatal("Flush did not change the file")
		}
		var buf bytes.Buffer
		if err := r.CopyTo(&buf); err != nil {
			t.Fatal(err)
		}
		if !bytes.Equal(buf.Bytes(), orig) {
			t.Fatalf("the reader opened before the resize no longer reads the original bytes (first diff at %d)", firstDiff(buf.Bytes(), orig))
		}
	})
}

// On a file system without record locks (fcntl fails with ENOLCK, EINVAL,
// ENOTSUP or EOPNOTSUPP, which differ on macOS) a read-only Open still
// works, and every in-place write is refused with nothing written, since
// readers cannot be kept out. The fcntl failure of the read/write lock is
// injected (the edit lock, a flock on Linux and a separate record lock on
// macOS, is left working).
func TestNoRecordLocksRefusesInPlaceWrites(t *testing.T) {
	defer func() { ofdLock = fcntlOFDLock }()
	errnos := []syscall.Errno{syscall.ENOLCK, syscall.EINVAL, syscall.ENOTSUP}
	if syscall.EOPNOTSUPP != syscall.ENOTSUP {
		errnos = append(errnos, syscall.EOPNOTSUPP)
	}
	for _, errno := range errnos {
		t.Run(errno.Error(), func(t *testing.T) {
			ofdLock = func(*os.File, int16, int) error { return errno }
			path, orig := copySample(t, "testprog.std")
			r, err := Open(path)
			if err != nil {
				t.Fatalf("Open without record locks: %v", err)
			}
			r.Close()
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			inPlaceHDU2(t, f)
			if err := f.Flush(); err == nil || !errors.Is(err, errno) {
				t.Fatalf("in-place Flush without record locks: %v, want an error wrapping %v", err, errno)
			}
			if !bytes.Equal(mustRead(t, path), orig) {
				t.Fatal("file changed")
			}
			if s := sidecars(t, filepath.Dir(path), "testprog.std"); len(s) != 0 {
				t.Fatalf("left %v", s)
			}
		})
	}
}
