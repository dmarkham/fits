//go:build linux || darwin

package fits

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"testing"
)

type openResult struct {
	f   *File
	err error
}

// Open takes the read lock first and only then looks for a journal, and it
// waits for a running in-place write instead of failing (the fits.go
// comment and lockForRead's doc). Here Open is called while an in-place
// Flush holds the write lock with its journal on disk:
//
//   - the write completes: Open must wait, then succeed and read the new
//     header. Checking the journal before taking the lock would fail with
//     ErrPendingJournal at once (the journal of a write still running);
//     taking the lock without waiting would fail at once too;
//   - the writer dies at that point (injected crash, which drops the write
//     lock as a process death does): Open must wait, then fail with
//     ErrPendingJournal, not read the half-written file.
//
// Truth: the order the doc states; the kernel refused the non-waiting
// request (lockWaitHook) and, on Linux, /proc/locks shows the blocked one.
func TestOpenWaitsForInPlaceWriteThenChecksJournal(t *testing.T) {
	for _, crash := range []bool{false, true} {
		t.Run(fmt.Sprintf("crash=%v", crash), func(t *testing.T) {
			watch := watchLockWaits(t)
			path, _ := copySample(t, "testprog.std")
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			inPlaceHDU2(t, f)
			res := make(chan openResult, 1)
			done := make(chan struct{})
			waited := false
			flushFaultHook = func(step string) error {
				if step != "header-written" {
					return nil
				}
				go func() {
					r, err := Open(path)
					res <- openResult{r, err}
					close(done)
				}()
				waited = watch.await(t, path, "READ", done)
				if crash {
					return errors.New("crash")
				}
				return nil
			}
			ferr := f.Flush()
			flushFaultHook = nil
			got := <-res
			if got.f != nil {
				defer got.f.Close()
			}
			if !waited {
				t.Fatalf("Open during the in-place write returned at once (err %v); it must wait for the write", got.err)
			}
			if crash {
				if ferr == nil {
					t.Fatal("the injected crash did not stop Flush")
				}
				if !errors.Is(got.err, ErrPendingJournal) {
					t.Fatalf("Open after the writer died mid-write: %v, want ErrPendingJournal", got.err)
				}
				return
			}
			if ferr != nil {
				t.Fatal(ferr)
			}
			if got.err != nil {
				t.Fatalf("Open that waited for the write: %v", got.err)
			}
			h, err := got.f.HDU(2)
			if err != nil {
				t.Fatal(err)
			}
			if v, _ := h.Header().String("INPLACE"); v != "same block count" {
				t.Fatalf("Open that waited reads INPLACE = %q, want the written header", v)
			}
		})
	}
}

// Close of an OpenReader *File releases its read lock while the caller's
// *os.File stays open, and so does an OpenReader that fails after taking
// the lock (the scan of a non-FITS file): an in-place write afterwards
// succeeds. Without the release the caller's descriptor would keep every
// in-place write refused with ErrOpenForReading. Truth: the documented
// contract; byte equality of the written file with the fault-free result.
func TestOpenReaderReleasesLockWithCallerFileOpen(t *testing.T) {
	_, post := func() (string, []byte) {
		p, _ := copySample(t, "testprog.std")
		f, err := OpenForEdit(p)
		if err != nil {
			t.Fatal(err)
		}
		inPlaceHDU2(t, f)
		if err := f.Flush(); err != nil {
			t.Fatal(err)
		}
		f.Close()
		return p, mustRead(t, p)
	}()
	flushWithCallerFileOpen := func(t *testing.T, path string, fh *os.File) {
		t.Helper()
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		inPlaceHDU2(t, f)
		if err := f.Flush(); err != nil {
			t.Fatalf("in-place Flush with the caller's file still open after OpenReader: %v", err)
		}
		if !bytes.Equal(mustRead(t, path), post) {
			t.Fatal("Flush result differs from the fault-free result")
		}
		if _, err := fh.Stat(); err != nil {
			t.Fatalf("the caller's file was closed: %v", err)
		}
	}
	t.Run("Close", func(t *testing.T) {
		path, _ := copySample(t, "testprog.std")
		fh, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer fh.Close()
		r, err := OpenReader(fh)
		if err != nil {
			t.Fatal(err)
		}
		if err := r.Close(); err != nil {
			t.Fatal(err)
		}
		flushWithCallerFileOpen(t, path, fh)
	})
	t.Run("failed-scan", func(t *testing.T) {
		path, orig := copySample(t, "testprog.std")
		if err := os.WriteFile(path, bytes.Repeat([]byte("x"), 2880), 0o644); err != nil {
			t.Fatal(err)
		}
		fh, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		defer fh.Close()
		if r, err := OpenReader(fh); err == nil {
			r.Close()
			t.Fatal("OpenReader accepted a non-FITS file")
		}
		// Same inode, now a FITS file.
		if err := os.WriteFile(path, orig, 0o644); err != nil {
			t.Fatal(err)
		}
		flushWithCallerFileOpen(t, path, fh)
	})
}

// The temp file EditFile renames over the path is locked before the
// rename, so the path is never unlocked while the edit runs: an
// OpenForEdit or EditFile between the rename and the directory fsync
// fails with ErrLocked. Truth: the EditFile doc (flock(2) contract).
func TestEditFileKeepsPathLockedAfterRename(t *testing.T) {
	path, _ := copySample(t, "testprog.std")
	var during []error
	inHook := false
	flushFaultHook = func(step string) error {
		if step != "renamed" || inHook {
			return nil // the calls below must not recurse when they succeed
		}
		inHook = true
		defer func() { inHook = false }()
		f, err := OpenForEdit(path)
		if err == nil {
			f.Close()
		}
		during = append(during, err)
		during = append(during, EditFile(path, copyAll))
		return nil
	}
	err := EditFile(path, copyAll)
	flushFaultHook = nil
	if err != nil {
		t.Fatal(err)
	}
	if len(during) != 2 {
		t.Fatalf("the renamed step ran %d times", len(during)/2)
	}
	for i, name := range []string{"OpenForEdit", "EditFile"} {
		if !errors.Is(during[i], ErrLocked) {
			t.Errorf("%s after EditFile's rename, before it returned: %v, want ErrLocked", name, during[i])
		}
	}
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("OpenForEdit after EditFile returned: %v", err)
	}
	f.Close()
}

// Journal rollback takes the write lock and waits for it: a reader that
// took its read lock just as the writer died (and is about to see the
// journal and refuse) only delays the recovery, it does not make it fail.
// Truth: the rollbackJournal / lockForWrite doc; byte equality of the
// recovered file with the original.
func TestRollbackWaitsForReaderLock(t *testing.T) {
	watch := watchLockWaits(t)
	path, orig := copySample(t, "testprog.std")
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	inPlaceHDU2(t, f)
	flushFaultHook = func(step string) error {
		if step == "header-written" {
			return errors.New("crash")
		}
		return nil
	}
	ferr := f.Flush()
	flushFaultHook = nil
	f.Close()
	if ferr == nil {
		t.Fatal("the injected crash did not stop Flush")
	}
	if _, err := os.Stat(journalPath(path)); err != nil {
		t.Fatalf("no journal after the crash: %v", err)
	}
	reader, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	if err := lockForRead(reader); err != nil {
		t.Fatal(err)
	}
	res := make(chan error, 1)
	done := make(chan struct{})
	go func() {
		e, err := OpenForEdit(path)
		if err == nil {
			e.Close()
		}
		res <- err
		close(done)
	}()
	waited := watch.await(t, path, "WRITE", done)
	reader.Close()
	if err := <-res; err != nil || !waited {
		t.Fatalf("OpenForEdit recovery while a reader held the lock: err %v, waited %v; it must wait and recover", err, waited)
	}
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("recovery did not restore the original bytes")
	}
}

// OpenReader of an *os.File counts a journal next to its name only while
// that name still refers to the file it has open: after the path was
// replaced by another file with a pending journal, the old descriptor is
// read without a false ErrPendingJournal, while Open of the path refuses.
// Truth: the OpenReader doc; the pixels of the file the descriptor holds.
func TestOpenReaderIgnoresJournalOfReplacedPath(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "a.fits")
	if err := createU8(t, path, 1); err != nil {
		t.Fatal(err)
	}
	fh, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	other := filepath.Join(dir, "b.fits")
	if err := createU8(t, other, 2); err != nil {
		t.Fatal(err)
	}
	if err := os.Rename(other, path); err != nil {
		t.Fatal(err)
	}
	st, _ := os.Stat(path)
	if err := os.WriteFile(journalPath(path), encodeJournal(undoJournal{FileSize: st.Size()}), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := Open(path); !errors.Is(err, ErrPendingJournal) {
		t.Fatalf("Open of the path with its journal: %v, want ErrPendingJournal", err)
	}
	r, err := OpenReader(fh)
	if err != nil {
		t.Fatalf("OpenReader of the replaced file: %v, want no error (the journal belongs to the new file)", err)
	}
	defer r.Close()
	p, err := r.Primary()
	if err != nil {
		t.Fatal(err)
	}
	px, err := ReadPixels[uint8](p)
	if err != nil || !bytes.Equal(px, bytes.Repeat([]byte{1}, 8)) {
		t.Fatalf("pixels %v (%v), want the replaced file's 1s", px, err)
	}
}
