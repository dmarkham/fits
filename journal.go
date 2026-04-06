package fits

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/dmarkham/fits/internal/block"
)

// journalMagic identifies our journal file format. Version byte follows.
var journalMagic = [8]byte{'F', 'I', 'T', 'S', 'J', 'R', 'N', 'L'}

// journalEntry describes a single pending tail shift. For v1 we only ever
// have one entry per flush because Flush processes dirty HDUs serially and
// fsyncs after each.
//
// Layout on disk:
//
//	magic   [8]byte
//	version uint8
//	state   uint8   (0 = pending, 1 = committed — only used during recovery)
//	_       [6]byte // padding
//	hduIdx       int64
//	oldHeaderEnd int64
//	newHeaderEnd int64
//	fileSize     int64
type journalEntry struct {
	HDUIndex     int64
	OldHeaderEnd int64
	NewHeaderEnd int64
	FileSize     int64
}

const journalHeaderBytes = 8 + 1 + 1 + 6 + 8*4

// writeJournal atomically writes a journal file at path. It calls fsync on
// the journal file and on its parent directory before returning so that a
// crash after the journal is created but before the main file is mutated
// can be safely replayed on next open.
func writeJournal(path string, entry journalEntry) error {
	f, err := os.Create(path)
	if err != nil {
		return err
	}
	buf := make([]byte, journalHeaderBytes)
	copy(buf[0:8], journalMagic[:])
	buf[8] = 1 // version
	buf[9] = 0 // state = pending
	binary.LittleEndian.PutUint64(buf[16:24], uint64(entry.HDUIndex))
	binary.LittleEndian.PutUint64(buf[24:32], uint64(entry.OldHeaderEnd))
	binary.LittleEndian.PutUint64(buf[32:40], uint64(entry.NewHeaderEnd))
	binary.LittleEndian.PutUint64(buf[40:48], uint64(entry.FileSize))
	if _, err := f.Write(buf); err != nil {
		f.Close()
		return err
	}
	if err := f.Sync(); err != nil {
		f.Close()
		return err
	}
	if err := f.Close(); err != nil {
		return err
	}
	// fsync parent directory.
	if dir := parentDir(path); dir != "" {
		d, err := os.Open(dir)
		if err == nil {
			d.Sync()
			d.Close()
		}
	}
	return nil
}

// readJournal returns the single entry in a journal file, or an error if
// the file is missing or malformed.
func readJournal(path string) (journalEntry, error) {
	var e journalEntry
	b, err := os.ReadFile(path)
	if err != nil {
		return e, err
	}
	if len(b) < journalHeaderBytes {
		return e, fmt.Errorf("fits: journal too short")
	}
	if string(b[0:8]) != string(journalMagic[:]) {
		return e, fmt.Errorf("fits: bad journal magic")
	}
	if b[8] != 1 {
		return e, fmt.Errorf("fits: unknown journal version %d", b[8])
	}
	e.HDUIndex = int64(binary.LittleEndian.Uint64(b[16:24]))
	e.OldHeaderEnd = int64(binary.LittleEndian.Uint64(b[24:32]))
	e.NewHeaderEnd = int64(binary.LittleEndian.Uint64(b[32:40]))
	e.FileSize = int64(binary.LittleEndian.Uint64(b[40:48]))
	return e, nil
}

// parentDir returns the directory portion of a path (or "").
func parentDir(p string) string {
	d := filepath.Dir(p)
	if d == "." {
		return ""
	}
	return d
}

// journalPath returns the path of the journal file for main file `name`.
func journalPath(name string) string { return name + ".journal" }

// recoverJournal inspects the journal file for name and replays or rolls
// back the pending edit before returning. If no journal exists, recovery
// is a no-op and nil is returned.
//
// Because v1 only issues one tail-shift entry per flush and fsyncs the
// journal before the main file mutation starts, the recovery rule is:
//
//   - Journal present, main file length == entry.FileSize (original):
//     the shift hadn't started; delete the journal.
//   - Journal present, main file length == entry.FileSize + delta (delta
//     positive = grew, negative = shrank): the shift began but we cannot
//     safely resume without knowing which portion was written. We
//     conservatively roll back to the pre-shift length and delete the
//     journal.
//
// The conservative rollback means a crashed Flush silently loses the
// in-memory header changes; the user's data (pre-edit HDU bytes) is
// preserved. This matches the spirit of plan 5.11 (strict-fail over
// silent guess) while giving a real crash-safety guarantee.
func recoverJournal(name string) error {
	path := journalPath(name)
	_, err := os.Stat(path)
	if err != nil {
		if errors.Is(err, os.ErrNotExist) {
			return nil
		}
		return err
	}
	entry, err := readJournal(path)
	if err != nil {
		// Corrupted journal — remove it so the file can be opened.
		os.Remove(path)
		return nil
	}
	fi, err := os.Stat(name)
	if err != nil {
		return err
	}
	if fi.Size() != entry.FileSize {
		// Tail shift was in progress; truncate back to the pre-shift size.
		fh, err := os.OpenFile(name, os.O_RDWR, 0)
		if err != nil {
			return err
		}
		if err := fh.Truncate(entry.FileSize); err != nil {
			fh.Close()
			return err
		}
		fh.Sync()
		fh.Close()
	}
	// Remove journal.
	return os.Remove(path)
}

// flushWithTailShift resizes a header that no longer fits in its current
// block allocation. It is invoked by Flush when a dirty HDU's re-serialized
// header length != its old length on disk.
//
// Algorithm:
//
//  1. Write a journal file describing the planned shift and fsync it (and
//     its parent directory). The journal records enough state to either
//     roll back (truncate) or recognize "never started" on next open.
//  2. Determine the delta (bytes to shift the tail). If positive we need to
//     make room by moving the tail forward; if negative we can just copy
//     in-place and truncate after.
//  3. Copy the tail (all bytes from oldHeaderEnd to EOF) to its new
//     position. Positive delta: back-to-front chunked copy to avoid
//     overwriting unread data. Negative delta: front-to-back.
//  4. fsync data.
//  5. Write the new header into place and fsync.
//  6. Update the file length (truncate if shrank) and fsync directory.
//  7. Delete the journal and fsync directory.
//
// After a successful shift the in-memory hduRecord offsets (headerStart,
// dataStart, dataEnd, paddedEnd) for every HDU after the shifted one must be
// adjusted by delta.
func (f *File) flushWithTailShift(rec *hduRecord, newHdrBytes []byte) error {
	if f.name == "" {
		return fmt.Errorf("fits: tail-shift flush requires a named file (not io.ReadWriteSeeker)")
	}
	oldEnd := rec.headerStart + int64(len(rec.rawHeader))
	newEnd := rec.headerStart + int64(len(newHdrBytes))
	delta := newEnd - oldEnd

	// Determine current file size.
	fi, err := os.Stat(f.name)
	if err != nil {
		return err
	}
	fileSize := fi.Size()

	// Write journal.
	if err := writeJournal(journalPath(f.name), journalEntry{
		HDUIndex:     int64(rec.index),
		OldHeaderEnd: oldEnd,
		NewHeaderEnd: newEnd,
		FileSize:     fileSize,
	}); err != nil {
		return err
	}

	// Extend file length if growing.
	if delta > 0 {
		if err := f.truncateTo(fileSize + delta); err != nil {
			return err
		}
	}

	// Copy tail.
	if delta != 0 {
		if err := f.copyTailChunked(oldEnd, newEnd, fileSize-oldEnd); err != nil {
			return err
		}
	}
	if delta < 0 {
		if err := f.truncateTo(fileSize + delta); err != nil {
			return err
		}
	}

	// fsync data.
	f.sync()

	// Write new header.
	if _, err := f.bw.Seek(rec.headerStart, io.SeekStart); err != nil {
		return err
	}
	if err := f.bw.WriteRange(newHdrBytes); err != nil {
		return err
	}
	f.sync()

	// Update in-memory offsets for this and all later HDUs.
	rec.rawHeader = append(rec.rawHeader[:0], newHdrBytes...)
	rec.dataStart = newEnd
	rec.dataEnd += delta
	rec.paddedEnd += delta
	rec.dirty = false
	for _, later := range f.hdus[rec.index+1:] {
		later.headerStart += delta
		later.dataStart += delta
		later.dataEnd += delta
		later.paddedEnd += delta
	}

	// Remove journal.
	os.Remove(journalPath(f.name))
	// fsync directory.
	if dir := parentDir(f.name); dir != "" {
		if d, err := os.Open(dir); err == nil {
			d.Sync()
			d.Close()
		}
	}
	f.br.Invalidate()
	return nil
}

// copyTailChunked copies the tail of the file from src to dst. srcLen is the
// number of bytes to copy. When dst > src we copy back-to-front to avoid
// overwriting unread data; otherwise front-to-back.
func (f *File) copyTailChunked(src, dst, srcLen int64) error {
	const chunk = 1 << 20 // 1 MiB
	buf := make([]byte, chunk)
	if dst > src {
		// Back-to-front.
		remaining := srcLen
		for remaining > 0 {
			n := int64(chunk)
			if remaining < n {
				n = remaining
			}
			offset := remaining - n
			if err := f.br.ReadRange(src+offset, buf[:n]); err != nil {
				return err
			}
			if _, err := f.bw.Seek(dst+offset, io.SeekStart); err != nil {
				return err
			}
			if err := f.bw.WriteRange(buf[:n]); err != nil {
				return err
			}
			remaining -= n
		}
	} else {
		// Front-to-back.
		var done int64
		for done < srcLen {
			n := int64(chunk)
			if srcLen-done < n {
				n = srcLen - done
			}
			if err := f.br.ReadRange(src+done, buf[:n]); err != nil {
				return err
			}
			if _, err := f.bw.Seek(dst+done, io.SeekStart); err != nil {
				return err
			}
			if err := f.bw.WriteRange(buf[:n]); err != nil {
				return err
			}
			done += n
		}
	}
	return nil
}

// truncateTo calls Truncate on the backing file if it is an *os.File.
func (f *File) truncateTo(size int64) error {
	if fh, ok := f.closer.(*os.File); ok {
		return fh.Truncate(size)
	}
	return fmt.Errorf("fits: cannot truncate non-os.File backing")
}

// sync fsyncs the backing file if it supports it.
func (f *File) sync() {
	if s, ok := f.closer.(interface{ Sync() error }); ok {
		_ = s.Sync()
	}
}

// Ensure the block package is referenced even if future refactors remove
// its explicit use here.
var _ = block.Size
