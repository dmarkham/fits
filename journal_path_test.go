//go:build linux || darwin

package fits

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// crashInPlaceFlush runs an in-place Flush through path and stops it (an
// injected crash) after the header write, leaving the undo journal.
func crashInPlaceFlush(t *testing.T, path string) {
	t.Helper()
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	inPlaceHDU2(t, f)
	flushFaultHook = func(step string) error {
		if step == "header-written" {
			return errInjectedCrash
		}
		return nil
	}
	err = f.Flush()
	flushFaultHook = nil
	if !errors.Is(err, errInjectedCrash) {
		t.Fatalf("Flush: %v, want the injected crash", err)
	}
}

// checkJournalSeenAndRecovered: every spelling in readers must see the
// pending journal (Open fails with ErrPendingJournal), and an OpenForEdit
// through recoverVia must roll it back to orig.
func checkJournalSeenAndRecovered(t *testing.T, readers []string, recoverVia, real string, orig []byte) {
	t.Helper()
	if bytes.Equal(mustRead(t, real), orig) {
		t.Fatal("the crash point is before the header write")
	}
	for _, p := range readers {
		if r, err := Open(p); !errors.Is(err, ErrPendingJournal) {
			if r != nil {
				r.Close()
			}
			t.Errorf("Open(%s) with the journal pending: %v, want ErrPendingJournal", p, err)
		}
	}
	f, err := OpenForEdit(recoverVia)
	if err != nil {
		t.Fatalf("OpenForEdit(%s): %v", recoverVia, err)
	}
	f.Close()
	if !bytes.Equal(mustRead(t, real), orig) {
		t.Fatalf("OpenForEdit(%s) did not roll the journal back", recoverVia)
	}
	if s := sidecars(t, filepath.Dir(real), filepath.Base(real)); len(s) != 0 {
		t.Fatalf("left %v", s)
	}
}

// The journal is named after the file with symbolic links resolved, so a
// crash through one spelling of the path is seen and recovered through
// any other: through a symlinked directory and through the real one, both
// ways round, through a link to the file in another directory, and
// through the fully resolved spelling of the temp dir
// (on macOS t.TempDir is under /var, a link to /private/var). Truth: the
// documented ErrPendingJournal and byte equality with the original.
func TestJournalFoundThroughSymlinkedDir(t *testing.T) {
	setup := func(t *testing.T) (realPath, linkPath, resolved string, orig []byte) {
		realPath, orig = copySample(t, "testprog.std")
		link := filepath.Join(t.TempDir(), "link")
		if err := os.Symlink(filepath.Dir(realPath), link); err != nil {
			t.Fatal(err)
		}
		linkPath = filepath.Join(link, filepath.Base(realPath))
		resolved, err := filepath.EvalSymlinks(realPath)
		if err != nil {
			t.Fatal(err)
		}
		return realPath, linkPath, resolved, orig
	}
	t.Run("crash-via-link", func(t *testing.T) {
		realPath, linkPath, resolved, orig := setup(t)
		crashInPlaceFlush(t, linkPath)
		checkJournalSeenAndRecovered(t, []string{realPath, linkPath, resolved}, realPath, realPath, orig)
	})
	t.Run("crash-via-real", func(t *testing.T) {
		realPath, linkPath, resolved, orig := setup(t)
		crashInPlaceFlush(t, realPath)
		checkJournalSeenAndRecovered(t, []string{realPath, linkPath, resolved}, linkPath, realPath, orig)
	})
	// A link to the file itself from another directory: a journal named
	// after the link would sit in that other directory, where Open of the
	// real path does not look.
	t.Run("crash-via-file-link", func(t *testing.T) {
		realPath, linkPath, resolved, orig := setup(t)
		alias := filepath.Join(t.TempDir(), "alias.fits")
		if err := os.Symlink(realPath, alias); err != nil {
			t.Fatal(err)
		}
		crashInPlaceFlush(t, alias)
		checkJournalSeenAndRecovered(t, []string{realPath, linkPath, resolved, alias}, realPath, realPath, orig)
		if s := sidecars(t, filepath.Dir(alias), "alias.fits"); len(s) != 0 {
			t.Fatalf("left %v next to the link", s)
		}
	})
	t.Run("crash-via-resolved", func(t *testing.T) {
		realPath, linkPath, resolved, orig := setup(t)
		crashInPlaceFlush(t, resolved)
		checkJournalSeenAndRecovered(t, []string{realPath, linkPath, resolved}, linkPath, realPath, orig)
	})
}

// On a case-insensitive file system (the macOS default) two spellings
// that differ in case name one file, and their journal names name one
// journal: a crash through one spelling is seen and recovered through the
// other. Skipped on a case-sensitive file system (Linux, case-sensitive
// APFS), where the two spellings are two files. Truth: the documented
// ErrPendingJournal and byte equality with the original.
func TestJournalFoundWithOtherCase(t *testing.T) {
	dir := t.TempDir()
	orig := mustRead(t, func() string { p, _ := copySample(t, "testprog.std"); return p }())
	lower := filepath.Join(dir, "case.fits")
	if err := os.WriteFile(lower, orig, 0o644); err != nil {
		t.Fatal(err)
	}
	upper := filepath.Join(dir, "CASE.FITS")
	a, err := os.Stat(lower)
	if err != nil {
		t.Fatal(err)
	}
	b, err := os.Stat(upper)
	if err != nil || !os.SameFile(a, b) {
		t.Skipf("the file system of %s is case-sensitive (CASE.FITS is not case.fits); the property is about case-insensitive file systems", dir)
	}
	mixed := filepath.Join(dir, "Case.Fits")
	crashInPlaceFlush(t, upper)
	if _, err := os.Stat(journalPath(lower)); err != nil {
		t.Fatalf("the journal written as %s is not found as %s: %v", journalPath(upper), journalPath(lower), err)
	}
	checkJournalSeenAndRecovered(t, []string{lower, upper, mixed}, mixed, lower, orig)
	if names, _ := os.ReadDir(dir); len(names) != 1 || !strings.EqualFold(names[0].Name(), "case.fits") {
		t.Fatalf("directory holds %v, want the one file", names)
	}
}
