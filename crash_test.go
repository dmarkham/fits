package fits

import (
	"os"
	"path/filepath"
	"testing"
)

// TestJournalRecoverStaleIsNoOp verifies that if a stale journal file exists
// but the main file is at its pre-shift length (i.e. the shift never
// started), OpenForEdit deletes the journal and opens normally.
func TestJournalRecoverStaleIsNoOp(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "j1.fits")
	f, _ := Create(path)
	WriteImage(f, nil, []int64{4}, []float32{1, 2, 3, 4})
	f.Close()
	origSize, _ := os.Stat(path)

	// Write a fake journal pointing at the current (pre-shift) state.
	err := writeJournal(journalPath(path), journalEntry{
		HDUIndex:     0,
		OldHeaderEnd: 2880,
		NewHeaderEnd: 5760,
		FileSize:     origSize.Size(),
	})
	if err != nil {
		t.Fatal(err)
	}

	// OpenForEdit must succeed and delete the journal.
	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("OpenForEdit: %v", err)
	}
	fe.Close()

	if _, err := os.Stat(journalPath(path)); !os.IsNotExist(err) {
		t.Fatalf("stale journal not removed")
	}

	// File size unchanged.
	fi, _ := os.Stat(path)
	if fi.Size() != origSize.Size() {
		t.Fatalf("file size changed: %d → %d", origSize.Size(), fi.Size())
	}
}

// TestJournalRecoverRollback verifies that if a shift began (file size
// differs from the journal's recorded pre-shift size), the recovery
// conservatively truncates back to the recorded size and deletes the
// journal. This preserves the user's data at the cost of losing the
// in-flight edit.
func TestJournalRecoverRollback(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "j2.fits")
	f, _ := Create(path)
	WriteImage(f, nil, []int64{4}, []float32{1, 2, 3, 4})
	f.Close()
	preSize, _ := os.Stat(path)

	// Simulate a crash mid-shift: write a journal claiming pre-size == preSize,
	// then extend the file by 2880 bytes (like a half-done shift).
	fh, _ := os.OpenFile(path, os.O_RDWR, 0)
	fh.Seek(0, 2)
	fh.Write(make([]byte, 2880))
	fh.Close()

	err := writeJournal(journalPath(path), journalEntry{
		HDUIndex:     0,
		OldHeaderEnd: 2880,
		NewHeaderEnd: 5760,
		FileSize:     preSize.Size(),
	})
	if err != nil {
		t.Fatal(err)
	}

	// Recovery should truncate back to preSize.
	fe, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("OpenForEdit: %v", err)
	}
	fe.Close()

	fi, _ := os.Stat(path)
	if fi.Size() != preSize.Size() {
		t.Fatalf("rollback failed: %d != %d", fi.Size(), preSize.Size())
	}
	if _, err := os.Stat(journalPath(path)); !os.IsNotExist(err) {
		t.Fatalf("journal not cleaned up")
	}
}

// TestEditFileFailureRemovesTempFile verifies that when the EditFile
// callback returns an error, the temp file is cleaned up and the source is
// untouched.
func TestEditFileFailureRemovesTempFile(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "ef.fits")
	f, _ := Create(path)
	WriteImage(f, nil, []int64{4}, []float32{1, 2, 3, 4})
	f.Close()

	origContents, _ := os.ReadFile(path)

	EditFile(path, func(in *File, out *Writer) error {
		return errTestForce // defined in editfile_test.go
	})

	// Source unchanged.
	after, _ := os.ReadFile(path)
	if len(after) != len(origContents) {
		t.Fatalf("source changed: %d → %d", len(origContents), len(after))
	}
	// No stray temp files.
	entries, _ := os.ReadDir(dir)
	for _, e := range entries {
		if e.Name() != "ef.fits" {
			t.Fatalf("stray file after failed EditFile: %s", e.Name())
		}
	}
}
