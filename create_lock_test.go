package fits

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"testing"
)

// createU8 writes a 4x2 uint8 primary with every pixel v through Create.
func createU8(t *testing.T, path string, v uint8) error {
	t.Helper()
	f, err := Create(path)
	if err != nil {
		return err
	}
	if _, err := WriteImage(f, nil, []int64{4, 2}, bytes.Repeat([]byte{v}, 8)); err != nil {
		f.Close()
		t.Fatal(err)
	}
	return f.Close()
}

func readU8(t *testing.T, path string) []uint8 {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer f.Close()
	p, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	px, err := ReadPixels[uint8](p)
	if err != nil {
		t.Fatal(err)
	}
	return px
}

// Create must not rewrite a file that another edit of this package owns:
//
//   - with a pending undo journal (an interrupted OverwritePixels), Create
//     used to truncate and rewrite the file, and the next OpenForEdit
//     "recovered" by writing the stale journal's bytes into the new file
//     (the journal records only the file size, so a new file of the same
//     size passed), which silently undid the new data. Now Create fails
//     with ErrPendingJournal, the file and journal untouched, and recovery
//     gives back the bytes from before the interrupted write;
//   - while an OpenForEdit editor holds the file, Create rewrote it under
//     the editor, whose index and later in-place writes assume the old
//     layout. Now it fails with ErrLocked, file untouched;
//   - while Open has the file, truncating it would tear the reader's view.
//     Now it fails with ErrOpenForReading, file untouched;
//   - while Create's *File is open no editor can open the file
//     (ErrLocked).
//
// Truth: byte equality with the file before the refused Create; the
// pixels the next reader sees.
func TestCreateRespectsJournalAndLocks(t *testing.T) {
	t.Run("pending-journal", func(t *testing.T) {
		path := filepath.Join(t.TempDir(), "c.fits")
		if err := createU8(t, path, 1); err != nil {
			t.Fatal(err)
		}
		pre := mustRead(t, path)
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		p, err := f.Primary()
		if err != nil {
			t.Fatal(err)
		}
		flushFaultHook = func(step string) error {
			if step == "data-synced" {
				return errors.New("crash")
			}
			return nil
		}
		err = OverwritePixels(p, bytes.Repeat([]byte{2}, 8))
		flushFaultHook = nil
		f.Close()
		if err == nil {
			t.Fatal("the injected crash did not stop OverwritePixels")
		}
		crashed := mustRead(t, path)
		journal := mustRead(t, journalPath(path))

		err = createU8(t, path, 3)
		if !errors.Is(err, ErrPendingJournal) {
			px := fmt.Sprint(readU8RecoveringT(t, path))
			t.Fatalf("Create over a file with a pending journal: %v, want ErrPendingJournal (after the next OpenForEdit the pixels are %s; Create wrote 3s)", err, px)
		}
		if !bytes.Equal(mustRead(t, path), crashed) || !bytes.Equal(mustRead(t, journalPath(path)), journal) {
			t.Fatal("refused Create changed the file or its journal")
		}
		e, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		e.Close()
		if !bytes.Equal(mustRead(t, path), pre) {
			t.Fatal("recovery after the refused Create did not restore the file from before the interrupted write")
		}
		if err := createU8(t, path, 3); err != nil {
			t.Fatalf("Create after recovery: %v", err)
		}
		if got := readU8(t, path); !bytes.Equal(got, bytes.Repeat([]byte{3}, 8)) {
			t.Fatalf("pixels %v after Create, want 3s", got)
		}
	})
	t.Run("editor-open", func(t *testing.T) {
		path := filepath.Join(t.TempDir(), "c.fits")
		if err := createU8(t, path, 1); err != nil {
			t.Fatal(err)
		}
		pre := mustRead(t, path)
		e, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		err = createU8(t, path, 3)
		if !errors.Is(err, ErrLocked) {
			e.Close()
			t.Fatalf("Create while an editor has the file: %v, want ErrLocked", err)
		}
		if !bytes.Equal(mustRead(t, path), pre) {
			t.Fatal("refused Create changed the file")
		}
		e.Close()
		if err := createU8(t, path, 3); err != nil {
			t.Fatalf("Create after the editor closed: %v", err)
		}
	})
	t.Run("reader-open", func(t *testing.T) {
		path := filepath.Join(t.TempDir(), "c.fits")
		if err := createU8(t, path, 1); err != nil {
			t.Fatal(err)
		}
		pre := mustRead(t, path)
		r, err := Open(path)
		if err != nil {
			t.Fatal(err)
		}
		err = createU8(t, path, 3)
		if !errors.Is(err, ErrOpenForReading) {
			r.Close()
			t.Fatalf("Create while the file is open for reading: %v, want ErrOpenForReading", err)
		}
		if !bytes.Equal(mustRead(t, path), pre) {
			t.Fatal("refused Create changed the file")
		}
		r.Close()
		if err := createU8(t, path, 3); err != nil {
			t.Fatalf("Create after the reader closed: %v", err)
		}
	})
	t.Run("create-open-blocks-editors", func(t *testing.T) {
		path := filepath.Join(t.TempDir(), "c.fits")
		f, err := Create(path)
		if err != nil {
			t.Fatal(err)
		}
		if _, err := WriteImage(f, nil, []int64{2}, []uint8{1, 2}); err != nil {
			t.Fatal(err)
		}
		if _, err := OpenForEdit(path); !errors.Is(err, ErrLocked) {
			t.Fatalf("OpenForEdit while Create has the file: %v, want ErrLocked", err)
		}
		if err := f.Close(); err != nil {
			t.Fatal(err)
		}
		e, err := OpenForEdit(path)
		if err != nil {
			t.Fatalf("OpenForEdit after Create closed: %v", err)
		}
		e.Close()
	})
	t.Run("replaced-before-lock", func(t *testing.T) {
		// The path is renamed over (as a resizing Flush does) between
		// Create's open and its lock: Create must write the file now at the
		// path, not the replaced inode.
		dir := t.TempDir()
		path := filepath.Join(dir, "c.fits")
		if err := createU8(t, path, 1); err != nil {
			t.Fatal(err)
		}
		other := filepath.Join(dir, "other.fits")
		if err := createU8(t, other, 9); err != nil {
			t.Fatal(err)
		}
		lockOpenHook = func() {
			lockOpenHook = nil
			if err := os.Rename(other, path); err != nil {
				t.Error(err)
			}
		}
		err := createU8(t, path, 3)
		lockOpenHook = nil
		if err != nil {
			t.Fatal(err)
		}
		if got := readU8(t, path); !bytes.Equal(got, bytes.Repeat([]byte{3}, 8)) {
			t.Fatalf("pixels at the path after Create %v, want 3s (Create wrote the replaced inode)", got)
		}
	})
}

// readU8RecoveringT runs the recovery the next editor would run and reads
// the pixels, for a failure message.
func readU8RecoveringT(t *testing.T, path string) []uint8 {
	t.Helper()
	e, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("OpenForEdit: %v", err)
	}
	e.Close()
	return readU8(t, path)
}
