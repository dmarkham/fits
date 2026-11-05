package fits

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Journal recovery edge cases. The crash-at-every-step proofs for Flush
// are in flush_shift_test.go; these cover journals that recovery must
// reject or discard, and a crash during recovery itself.

// copySample copies a committed testdata sample into a temp dir and
// returns the copy's path and the original bytes.
func copySample(t *testing.T, sample string) (string, []byte) {
	t.Helper()
	src := testutil.TestdataPath(t, sample)
	testutil.RequireFile(t, src)
	orig, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	path := filepath.Join(t.TempDir(), filepath.Base(sample))
	if err := os.WriteFile(path, orig, 0o644); err != nil {
		t.Fatal(err)
	}
	return path, orig
}

func mustRead(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	return b
}

// A version 1 journal (written by the earlier tail-shift code, which saved
// no pre-image) cannot be recovered: OpenForEdit must refuse and leave the
// file and the journal alone instead of guessing (the old recovery
// truncated the file, which left shifted bytes in place).
func TestRecoverLegacyJournalRefuses(t *testing.T) {
	path, orig := copySample(t, "iter_a.fit")
	v1 := make([]byte, 48)
	copy(v1, journalMagic[:])
	v1[8] = 1
	binary.LittleEndian.PutUint64(v1[24:], 2880)
	binary.LittleEndian.PutUint64(v1[32:], 5760)
	binary.LittleEndian.PutUint64(v1[40:], uint64(len(orig)))
	if err := os.WriteFile(journalPath(path), v1, 0o644); err != nil {
		t.Fatal(err)
	}
	_, err := OpenForEdit(path)
	if !errors.Is(err, ErrLegacyJournal) {
		t.Fatalf("OpenForEdit with a v1 journal: err = %v, want ErrLegacyJournal", err)
	}
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("file modified")
	}
	if !bytes.Equal(mustRead(t, journalPath(path)), v1) {
		t.Fatal("legacy journal modified or removed")
	}
}

// A journal that fails its CRC was torn while being written, which only
// happens before the FITS file is touched: it is removed and the file is
// left byte-identical.
func TestRecoverTornJournalRemoved(t *testing.T) {
	path, orig := copySample(t, "iter_a.fit")
	j := encodeJournal(undoJournal{FileSize: int64(len(orig)), Entries: []undoEntry{{Offset: 0, Data: bytes.Repeat([]byte("X"), 2880)}}})
	flipped := append([]byte{}, j...)
	flipped[100] ^= 1
	for _, torn := range [][]byte{j[:len(j)-1], j[:40], j[:9], flipped} {
		if err := os.WriteFile(journalPath(path), torn, 0o644); err != nil {
			t.Fatal(err)
		}
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatalf("OpenForEdit with torn journal (%d bytes): %v", len(torn), err)
		}
		f.Close()
		if !bytes.Equal(mustRead(t, path), orig) {
			t.Fatal("torn journal was applied to the file")
		}
		if _, err := os.Stat(journalPath(path)); !os.IsNotExist(err) {
			t.Fatal("torn journal not removed")
		}
	}
}

// A valid journal written for a file of a different size does not belong
// to this file's current state: recovery refuses rather than overwrite.
func TestRecoverJournalSizeMismatchRefuses(t *testing.T) {
	path, orig := copySample(t, "iter_a.fit")
	j := encodeJournal(undoJournal{FileSize: int64(len(orig)) + 2880, Entries: []undoEntry{{Offset: 0, Data: bytes.Repeat([]byte("X"), 2880)}}})
	if err := os.WriteFile(journalPath(path), j, 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := OpenForEdit(path); err == nil {
		t.Fatal("OpenForEdit accepted a journal written for a different file size")
	}
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("file modified")
	}
}

// A journal whose flags byte has a bit this release does not know was
// written by a later release, whose rollback this one cannot do: OpenForEdit
// and EditFile refuse it (CRC valid, so it is not a torn journal) and leave
// the file and the journal as they are. Every known flag combination is
// accepted. Truth: the journal format (encodeJournal's layout comment) and
// byte equality.
func TestRecoverJournalUnknownFlagsRefuses(t *testing.T) {
	for _, flags := range []byte{0x02, 0x80, 0x03} {
		for _, rv := range crashRecoverers {
			t.Run(fmt.Sprintf("%#x/%s", flags, rv.name), func(t *testing.T) {
				path, orig := copySample(t, "iter_a.fit")
				j := encodeJournal(undoJournal{FileSize: int64(len(orig)), Entries: []undoEntry{{Offset: 0, Data: bytes.Repeat([]byte("X"), 2880)}}})
				j[9] = flags
				binary.LittleEndian.PutUint32(j[len(j)-4:], crc32.ChecksumIEEE(j[:len(j)-4]))
				if err := os.WriteFile(journalPath(path), j, 0o644); err != nil {
					t.Fatal(err)
				}
				err := rv.run(path)
				if err == nil || !strings.Contains(err.Error(), "unknown journal flags") {
					t.Fatalf("%s with journal flags %#x: %v, want an unknown-flags error", rv.name, flags, err)
				}
				if !bytes.Equal(mustRead(t, path), orig) {
					t.Fatal("file modified")
				}
				if !bytes.Equal(mustRead(t, journalPath(path)), j) {
					t.Fatal("journal removed or changed")
				}
			})
		}
	}
	for _, truncate := range []bool{false, true} {
		path, orig := copySample(t, "iter_a.fit")
		saved := append([]byte(nil), orig[:2880]...)
		bad := bytes.Repeat([]byte("X"), 2880)
		if err := os.WriteFile(path, append(bad, orig[2880:]...), 0o644); err != nil {
			t.Fatal(err)
		}
		j := encodeJournal(undoJournal{FileSize: int64(len(orig)), Truncate: truncate, Entries: []undoEntry{{Offset: 0, Data: saved}}})
		if err := os.WriteFile(journalPath(path), j, 0o644); err != nil {
			t.Fatal(err)
		}
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatalf("known flags (truncate=%v): %v", truncate, err)
		}
		f.Close()
		if !bytes.Equal(mustRead(t, path), orig) {
			t.Fatalf("known flags (truncate=%v): not rolled back", truncate)
		}
	}
}

// A crash during recovery (after the saved bytes are written back, before
// the journal is removed) is recovered by the next OpenForEdit.
func TestRecoveryCrashIsRecoverable(t *testing.T) {
	path, orig := copySample(t, "testprog.std")
	flushFaultHook = func(step string) error {
		if step == "header-written" {
			return errInjectedCrash
		}
		return nil
	}
	err := editAndFlush(t, path, inPlaceHDU2)
	flushFaultHook = nil
	if !errors.Is(err, errInjectedCrash) {
		t.Fatalf("Flush: %v, want injected crash", err)
	}
	if bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("crash point is before the header write; test does not exercise rollback")
	}
	flushFaultHook = func(step string) error {
		if step == "recovery-restored" {
			return errInjectedCrash
		}
		return nil
	}
	_, err = OpenForEdit(path)
	flushFaultHook = nil
	if !errors.Is(err, errInjectedCrash) {
		t.Fatalf("first recovery: %v, want injected crash", err)
	}
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatalf("second recovery: %v", err)
	}
	f.Close()
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("file after two recoveries differs from the original")
	}
	if s := sidecars(t, filepath.Dir(path), "testprog.std"); len(s) != 0 {
		t.Fatalf("left behind: %v", s)
	}
}

// A read-only Open of a file whose Flush crashed half way (one of two
// headers written) refuses with ErrPendingJournal, also through a symbolic
// link, and with a legacy journal. After OpenForEdit rolls the journal
// back, Open works and reads the original header values.
func TestOpenRefusesPendingJournal(t *testing.T) {
	path, orig := copySample(t, "testprog.std")
	link := filepath.Join(t.TempDir(), "link.fits")
	if err := os.Symlink(path, link); err != nil {
		t.Fatal(err)
	}
	flushFaultHook = func(step string) error {
		if step == "header-written" {
			return errInjectedCrash
		}
		return nil
	}
	err := editAndFlush(t, path, func(t *testing.T, f *File) {
		inPlaceHDU2(t, f)
		if err := hduHeader(t, f, 3).Set("EXTNAME", "renamed-ascii", ""); err != nil {
			t.Fatal(err)
		}
	})
	flushFaultHook = nil
	if !errors.Is(err, errInjectedCrash) {
		t.Fatalf("Flush: %v, want injected crash", err)
	}
	if bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("crash point is before the first header write")
	}
	for _, p := range []string{path, link} {
		if f, err := Open(p); !errors.Is(err, ErrPendingJournal) {
			if f != nil {
				f.Close()
			}
			t.Fatalf("Open(%s) with a pending journal: err = %v, want ErrPendingJournal", p, err)
		}
	}
	f, err := OpenForEdit(link)
	if err != nil {
		t.Fatal(err)
	}
	f.Close()
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("recovery did not restore the original bytes")
	}
	r, err := Open(path)
	if err != nil {
		t.Fatalf("Open after recovery: %v", err)
	}
	defer r.Close()
	h3 := hduHeader(t, r, 3)
	if v, err := h3.String("EXTNAME"); err != nil || v != "Test-ASCII" {
		t.Fatalf("HDU 3 EXTNAME after recovery = %q, %v; want the original \"Test-ASCII\"", v, err)
	}

	// A legacy journal is pending too.
	path2, _ := copySample(t, "iter_a.fit")
	v1 := make([]byte, 48)
	copy(v1, journalMagic[:])
	v1[8] = 1
	if err := os.WriteFile(journalPath(path2), v1, 0o644); err != nil {
		t.Fatal(err)
	}
	if f, err := Open(path2); !errors.Is(err, ErrPendingJournal) {
		if f != nil {
			f.Close()
		}
		t.Fatalf("Open with a legacy journal: err = %v, want ErrPendingJournal", err)
	}
}

// TestFsyncOrder checks where every fsync happens relative to the commit
// steps. The crash-at-every-step tests model a process death, where data
// already written survives without fsync; this test covers the other half
// of the power-loss argument: the journal (and its directory entry) is
// durable before the FITS file is touched, the new headers are durable
// before the journal is removed, the temp file is durable before the
// rename, and the directory is synced after each rename or removal. It
// checks the calls, not what a kernel or disk does with them.
func TestFsyncOrder(t *testing.T) {
	var trace []string
	record := func() {
		trace = nil
		flushFaultHook = func(s string) error { trace = append(trace, s); return nil }
		syncHook = func(role string) { trace = append(trace, "fsync "+role) }
	}
	stop := func() { flushFaultHook, syncHook = nil, nil }
	defer stop()
	check := func(t *testing.T, want []string) {
		t.Helper()
		if strings.Join(trace, ", ") != strings.Join(want, ", ") {
			t.Errorf("fsync/step order:\n got  %v\n want %v", trace, want)
		}
	}
	journalTrace := func(headers int) []string {
		w := []string{"fsync journal", "fsync dir", "journal-written"}
		for range headers {
			w = append(w, "header-written")
		}
		return append(w, "fsync target", "data-synced", "fsync dir", "journal-removed")
	}
	renameTrace := []string{"temp-written", "fsync temp", "temp-synced", "renamed", "fsync dir", "dir-synced"}

	for _, sc := range shiftScenarios {
		t.Run("Flush/"+sc.name, func(t *testing.T) {
			path, _ := copySample(t, sc.sample)
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			sc.mutate(t, f)
			record()
			err = f.Flush()
			stop()
			f.Close()
			if err != nil {
				t.Fatal(err)
			}
			switch sc.name {
			case "inplace":
				check(t, journalTrace(1))
			case "inplace-two-hdus":
				check(t, journalTrace(2))
			default: // every other scenario resizes a header
				check(t, renameTrace)
			}
		})
	}
	// OverwritePixels and the edit-mode appends use the same journal
	// protocol around their data write.
	for _, op := range inPlaceOps {
		t.Run(op.name, func(t *testing.T) {
			_, path, _ := setupOp(t, op)
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			record()
			err = op.run(t, f)
			stop()
			f.Close()
			if err != nil {
				t.Fatal(err)
			}
			mid := "hdu-written"
			if strings.HasPrefix(op.name, "OverwritePixels") {
				mid = "pixels-written"
			}
			check(t, []string{"fsync journal", "fsync dir", "journal-written", mid, "fsync target", "data-synced", "fsync dir", "journal-removed"})
		})
	}
	t.Run("EditFile", func(t *testing.T) {
		path, _ := writeThreeHDUFile(t)
		record()
		err := EditFile(path, copyAll)
		stop()
		if err != nil {
			t.Fatal(err)
		}
		check(t, renameTrace)
	})
	t.Run("recovery", func(t *testing.T) {
		path, _ := copySample(t, "testprog.std")
		flushFaultHook = func(s string) error {
			if s == "header-written" {
				return errInjectedCrash
			}
			return nil
		}
		if err := editAndFlush(t, path, inPlaceHDU2); !errors.Is(err, errInjectedCrash) {
			t.Fatalf("Flush: %v, want injected crash", err)
		}
		record()
		f, err := OpenForEdit(path)
		stop()
		if err != nil {
			t.Fatal(err)
		}
		f.Close()
		check(t, []string{"fsync target", "recovery-restored", "fsync dir"})
	})
}

// OpenReader given the *os.File of a file with a pending journal refuses
// like Open (the journal is found next to the file the descriptor names,
// when that path is the same file). Readers that are not an *os.File carry
// no path and are not checked: a bytes.Reader of the same bytes opens.
func TestOpenReaderRefusesPendingJournal(t *testing.T) {
	path, orig := copySample(t, "testprog.std")
	flushFaultHook = func(step string) error {
		if step == "header-written" {
			return errInjectedCrash
		}
		return nil
	}
	err := editAndFlush(t, path, inPlaceHDU2)
	flushFaultHook = nil
	if !errors.Is(err, errInjectedCrash) {
		t.Fatalf("Flush: %v, want injected crash", err)
	}
	fh, err := os.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fh.Close()
	if f, err := OpenReader(fh); !errors.Is(err, ErrPendingJournal) {
		if f != nil {
			f.Close()
		}
		t.Fatalf("OpenReader(*os.File) with a pending journal: err = %v, want ErrPendingJournal", err)
	}
	if f, err := OpenReader(bytes.NewReader(mustRead(t, path))); err != nil {
		t.Fatalf("OpenReader(bytes.Reader): %v (no path, not checked)", err)
	} else {
		f.Close()
	}
	// A descriptor whose name is now another file is not judged by that
	// file's journal.
	other, _ := copySample(t, "testprog.std")
	if err := os.Rename(path, other); err != nil {
		t.Fatal(err)
	}
	_ = orig
	if f, err := OpenReader(fh); err != nil {
		t.Fatalf("OpenReader of a descriptor whose path was replaced: %v", err)
	} else {
		f.Close()
	}
}
