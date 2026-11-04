package fits

import (
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"hash/crc32"
	"io"
	"math"
	"os"
	"path/filepath"

	"github.com/dmarkham/fits/internal/block"
)

// Crash safety for in-place edits (OpenForEdit + Flush, OverwritePixels,
// appends in ModeEdit) and EditFile.
//
// Two mechanisms are used:
//
//   - Header rewrite with no size change, OverwritePixels, and an append at
//     the end of a ModeEdit file: an UNDO JOURNAL (journaledWrite). Before
//     the first byte of the FITS file is touched, the current on-disk bytes
//     of every region about to be overwritten (the headers, the data unit,
//     or the bytes from the append point to the end of the file) are
//     written to <name>.journal with the file size and a CRC-32 and fsynced
//     (file and directory). The new bytes are then written and fsynced, and
//     the journal is removed: that removal is the commit point. For a write
//     that can grow the file the journal is flagged to cut the file back to
//     its old size on rollback. If the process dies before it, OpenForEdit finds
//     the journal and writes the saved bytes back, so the file is exactly
//     what it was before the call. A journal that fails its CRC was torn
//     while being written, which can only happen before the FITS file is
//     touched, so it is simply removed. While a journal exists, a
//     read-only Open of the file fails with ErrPendingJournal instead of
//     reading headers that recovery would roll back. The write holds an
//     exclusive OFD lock from before the journal is written until after it
//     is removed, and a read-only Open holds a shared one until Close
//     (lock_ofd.go), so no reader of this package is open while bytes
//     change in place: the write is refused (ErrOpenForReading) instead.
//
//   - Header resize (a header gains or loses a 2880-byte block, so every
//     later byte moves): TEMP FILE + RENAME. The complete new file is
//     written to <name>.fits-tmp next to the target, fsynced and renamed
//     over the target (the commit point), then the directory is fsynced. A
//     crash before the rename leaves the original untouched (the stale
//     temp file is removed by the next OpenForEdit or EditFile); a crash
//     after it leaves the complete new file.
//
// Either way, once the next OpenForEdit or EditFile has run recovery, a
// crash leaves the file byte-identical to its state before the Flush or to
// the fully edited state, never anything in between. Until then a
// read-only Open refuses with ErrPendingJournal if a journal is pending.
//
// What the tests prove: for a process crash, by injecting one after every
// step (flushFaultHook) and comparing bytes. For a power loss the argument
// also needs every fsync to be issued before the step that relies on it;
// TestFsyncOrder checks that order (syncHook), but no test cuts power, so
// whether the kernel and disk honour fsync is not tested. On macOS every
// one of these syncs is F_FULLFSYNC (syncFile, sync_darwin.go), since a
// plain fsync there leaves the data in the drive's volatile cache.
//
// The resize path replaces the directory entry: a file with other hard
// links is refused (ErrHardLinked), and the new file gets the permission
// bits, owner and group (Unix) and extended attributes including POSIX
// ACLs (Linux) or the resource fork and Finder info (macOS, where an ACL
// is not carried over: xattr_darwin.go) of the old one, or the Flush is
// refused (ErrCannotPreserve) with the file untouched. A symbolic link is
// followed: the file it points to is replaced and the link is kept.
//
// Path spelling: the journal and temp file are named after the target
// with symbolic links resolved (resolveTarget), so every path that leads
// to the file through links (macOS /tmp and /var are links to /private/...)
// finds the same journal (TestJournalFoundThroughSymlinkedDir). On a
// case-insensitive file system (the macOS default) names that differ only
// in case are one file, and so are their journal names
// (TestJournalFoundWithOtherCase). A second hard link to the file is a
// different name with a different journal name: an in-place write through
// one name is not seen as pending through the other.

// journalMagic identifies our journal file format. Version byte follows.
var journalMagic = [8]byte{'F', 'I', 'T', 'S', 'J', 'R', 'N', 'L'}

// journalVersion 2 is the undo journal. Version 1 (the tail-shift journal
// written by earlier releases) saved no pre-image and cannot be recovered;
// see ErrLegacyJournal.
const journalVersion = 2

// ErrLegacyJournal is returned by OpenForEdit when <name>.journal was
// written by an earlier release of this package, whose header-resize
// journal did not save the bytes it overwrote. An interrupted Flush of that
// release can leave the FITS file half shifted, and it cannot be repaired
// automatically. Restore the file from a backup (or inspect it), then
// delete the journal file.
var ErrLegacyJournal = errors.New("fits: journal from an earlier release found; the file may be inconsistent and cannot be recovered automatically")

// ErrPendingJournal is returned by Open when <name>.journal exists: a Flush
// was interrupted and the file may hold some of its new headers. A
// read-only Open cannot roll the file back, so it refuses rather than show
// a half-applied edit. OpenForEdit or EditFile on the same path rolls the
// journal back (or reports ErrLegacyJournal); Open works again afterwards.
var ErrPendingJournal = errors.New("fits: an interrupted Flush left an undo journal; open the file with OpenForEdit or EditFile once to roll it back")

// ErrHardLinked is returned by a resizing Flush and by EditFile when the
// file has more than one hard link: replacing it through a temp file and a
// rename would leave the other links with the old content.
var ErrHardLinked = errors.New("fits: the file has other hard links, which replacing it would break")

// ErrCannotPreserve is returned by a resizing Flush and by EditFile when
// the new file cannot be given the old one's owner, group or extended
// attributes (POSIX ACLs included); the original is left untouched.
var ErrCannotPreserve = errors.New("fits: cannot carry the file's owner or extended attributes over to the rewritten file")

// ErrLocked is returned by OpenForEdit and EditFile when another editor
// (an OpenForEdit *File not yet closed, or an EditFile in progress) holds
// the file, in this process or another. The lock is advisory (flock(2) on
// Linux, an open file description record lock on macOS; see lock_ofd.go):
// programs that do not use this package are not stopped by it.
var ErrLocked = errors.New("fits: the file is being edited by another OpenForEdit or EditFile")

// ErrOpenForReading is returned by an in-place write (Flush without a
// header resize, OverwritePixels, an append in ModeEdit) while the file is
// open for reading through this package (Open, or OpenReader of an
// *os.File, not yet closed, in this process or another): changing bytes in
// place could hand that reader a header or data unit that is partly old
// and partly new. Nothing is written; close the readers and retry.
var ErrOpenForReading = errors.New("fits: the file is open for reading (Open/OpenReader); an in-place write is refused until the readers are closed")

// ErrEditUnsupported is returned by OpenForEdit and EditFile on systems
// other than Linux and macOS (Windows, the BSDs, ...), where this package
// cannot lock the file against other editors and readers or carry its
// extended attributes over to a rewritten file.
var ErrEditUnsupported = errors.New("fits: editing (OpenForEdit, EditFile) is supported on Linux and macOS only")

// checkNoPendingJournal returns ErrPendingJournal if name (after following
// symbolic links) has a journal next to it. A path that cannot be resolved
// is left for the caller's open to report.
func checkNoPendingJournal(name string) error {
	target, err := resolveTarget(name)
	if err != nil {
		return nil
	}
	jp := journalPath(target)
	if _, err := os.Lstat(jp); err == nil {
		return fmt.Errorf("%w (%s)", ErrPendingJournal, jp)
	}
	return nil
}

// checkDescriptorJournal is checkNoPendingJournal for an open file: the
// journal next to fh.Name() counts only if that path still names the file
// fh has open.
func checkDescriptorJournal(fh *os.File) error {
	a, err := fh.Stat()
	if err != nil {
		return nil
	}
	b, err := os.Stat(fh.Name())
	if err != nil || !os.SameFile(a, b) {
		return nil
	}
	return checkNoPendingJournal(fh.Name())
}

// tempSuffix names the temp file used by a resizing Flush and by EditFile.
const tempSuffix = ".fits-tmp"

// journalPath returns the path of the journal file for main file `name`.
func journalPath(name string) string { return name + ".journal" }

// undoEntry is one saved region of the FITS file.
type undoEntry struct {
	Offset int64
	Data   []byte
}

// undoJournal is the decoded content of a version 2 journal.
type undoJournal struct {
	FileSize int64 // size of the FITS file when the journal was written
	// Truncate is set when the protected write may grow the file (an
	// append, or OverwritePixels on a data unit that lacks its final
	// padding): rollback then also cuts the file back to FileSize.
	Truncate bool
	Entries  []undoEntry
}

// journalFlagTruncate is bit 0 of the flags byte.
const journalFlagTruncate = 1

// Journal layout (all integers little-endian):
//
//	magic    [8]byte "FITSJRNL"
//	version  uint8   2
//	flags    uint8   bit 0: truncate to fileSize on rollback (0 in
//	                 journals written before appends were journaled)
//	_        [6]byte
//	fileSize int64
//	count    int64
//	count x { offset int64; length int64; data [length]byte }
//	crc      uint32  CRC-32 (IEEE) of every preceding byte
func encodeJournal(j undoJournal) []byte {
	var b bytes.Buffer
	b.Write(journalMagic[:])
	b.WriteByte(journalVersion)
	var flags byte
	if j.Truncate {
		flags |= journalFlagTruncate
	}
	b.WriteByte(flags)
	b.Write(make([]byte, 6))
	var n [8]byte
	put := func(v int64) {
		binary.LittleEndian.PutUint64(n[:], uint64(v))
		b.Write(n[:])
	}
	put(j.FileSize)
	put(int64(len(j.Entries)))
	for _, e := range j.Entries {
		put(e.Offset)
		put(int64(len(e.Data)))
		b.Write(e.Data)
	}
	var c [4]byte
	binary.LittleEndian.PutUint32(c[:], crc32.ChecksumIEEE(b.Bytes()))
	b.Write(c[:])
	return b.Bytes()
}

// errTornJournal marks a journal that is incomplete or fails its CRC.
var errTornJournal = errors.New("fits: torn journal")

func decodeJournal(b []byte) (undoJournal, error) {
	var j undoJournal
	if len(b) < 9 || !bytes.Equal(b[:8], journalMagic[:]) {
		return j, errTornJournal
	}
	if b[8] != journalVersion {
		if b[8] == 1 {
			return j, ErrLegacyJournal
		}
		return j, fmt.Errorf("fits: unknown journal version %d", b[8])
	}
	if len(b) < 32+4 {
		return j, errTornJournal
	}
	body, sum := b[:len(b)-4], binary.LittleEndian.Uint32(b[len(b)-4:])
	if crc32.ChecksumIEEE(body) != sum {
		return j, errTornJournal
	}
	if flags := body[9]; flags&^journalFlagTruncate != 0 {
		return j, fmt.Errorf("fits: unknown journal flags %#x", flags)
	} else {
		j.Truncate = flags&journalFlagTruncate != 0
	}
	get := func(off int) int64 { return int64(binary.LittleEndian.Uint64(body[off : off+8])) }
	j.FileSize = get(16)
	count := get(24)
	off := 32
	for i := int64(0); i < count; i++ {
		if off+16 > len(body) {
			return j, errTornJournal
		}
		e := undoEntry{Offset: get(off)}
		n := get(off + 8)
		off += 16
		if n < 0 || int64(off)+n > int64(len(body)) {
			return j, errTornJournal
		}
		e.Data = body[off : off+int(n)]
		off += int(n)
		j.Entries = append(j.Entries, e)
	}
	if off != len(body) {
		return j, errTornJournal
	}
	return j, nil
}

// writeFileDurable creates path with content and fsyncs it and its
// directory.
func writeFileDurable(path string, content []byte) error {
	f, err := os.OpenFile(path, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0o644)
	if err != nil {
		return err
	}
	if _, err := f.Write(content); err != nil {
		f.Close()
		return err
	}
	if err := fsync(f, "journal"); err != nil {
		f.Close()
		return err
	}
	if err := f.Close(); err != nil {
		return err
	}
	return syncDir(path)
}

// removeDurable removes path (absent is fine) and fsyncs its directory.
func removeDurable(path string) error {
	if err := os.Remove(path); err != nil && !errors.Is(err, os.ErrNotExist) {
		return err
	}
	return syncDir(path)
}

// syncDir fsyncs the directory holding path so a create, rename or remove
// in it is durable.
func syncDir(path string) error {
	d, err := os.Open(filepath.Dir(path))
	if err != nil {
		return err
	}
	defer d.Close()
	return fsync(d, "dir")
}

// syncHook, when non-nil, is called after every successful fsync with the
// role of the synced file ("journal", "target", "temp" or "dir"). Tests
// record it next to the flushFaultHook steps to check that each fsync
// happens before the step that depends on it.
var syncHook func(role string)

// fsync syncs f (syncFile: fsync on Linux, F_FULLFSYNC on macOS) and
// reports it to syncHook.
func fsync(f *os.File, role string) error {
	if err := syncFile(f); err != nil {
		return err
	}
	if syncHook != nil {
		syncHook(role)
	}
	return nil
}

// resolveTarget follows symbolic links so temp files and renames act on the
// real file and keep the link.
func resolveTarget(name string) (string, error) {
	return filepath.EvalSymlinks(name)
}

// recoverEdit puts name back into a consistent state after a crashed Flush
// or EditFile: it removes a stale temp file and rolls back an undo journal.
// With nothing to recover it does nothing.
func recoverEdit(name string) error {
	target, err := resolveTarget(name)
	if err != nil {
		return err
	}
	if err := removeStaleTemp(target); err != nil {
		return err
	}
	return rollbackJournal(target)
}

// removeStaleTemp removes target's temp file. A temp file is never
// authoritative: the rename over the target is the commit point, so a temp
// file that still exists was never committed.
func removeStaleTemp(target string) error {
	if _, err := os.Lstat(target + tempSuffix); err != nil {
		return nil
	}
	return removeDurable(target + tempSuffix)
}

// rollbackJournal restores the regions saved in target's journal, if any,
// then removes the journal. Restoring is idempotent, so a crash during
// recovery is recovered by running it again.
func rollbackJournal(target string) error {
	jp := journalPath(target)
	b, err := os.ReadFile(jp)
	if errors.Is(err, os.ErrNotExist) {
		return nil
	}
	if err != nil {
		return err
	}
	j, err := decodeJournal(b)
	if errors.Is(err, errTornJournal) {
		// The journal is fsynced before the FITS file is touched, so a torn
		// journal means the FITS file was never modified.
		return removeDurable(jp)
	}
	if err != nil {
		return fmt.Errorf("%w (journal %s)", err, jp)
	}
	fh, err := os.OpenFile(target, os.O_RDWR, 0)
	if err != nil {
		return err
	}
	if err := lockForWrite(fh, true); err != nil {
		fh.Close()
		return err
	}
	st, err := fh.Stat()
	if err != nil {
		fh.Close()
		return err
	}
	// The protected write never shrinks the file, and grows it only when
	// the journal says so.
	if st.Size() < j.FileSize || (st.Size() > j.FileSize && !j.Truncate) {
		fh.Close()
		return fmt.Errorf("fits: journal %s was written for a %d-byte file but %s is %d bytes; not restoring (check the file before removing the journal)", jp, j.FileSize, target, st.Size())
	}
	for _, e := range j.Entries {
		if _, err := fh.WriteAt(e.Data, e.Offset); err != nil {
			fh.Close()
			return err
		}
	}
	if j.Truncate && st.Size() > j.FileSize {
		if err := fh.Truncate(j.FileSize); err != nil {
			fh.Close()
			return err
		}
	}
	if err := fsync(fh, "target"); err != nil {
		fh.Close()
		return err
	}
	if err := fh.Close(); err != nil {
		return err
	}
	if err := faultPoint("recovery-restored"); err != nil {
		return err
	}
	return removeDurable(jp)
}

// pendingHeader is one header Flush will write.
type pendingHeader struct {
	rec      *hduRecord
	newBytes []byte
}

// overwriteHeadersJournaled writes headers whose length did not change,
// protected by the undo journal.
func (f *File) overwriteHeadersJournaled(plan []pendingHeader) error {
	var regions []undoRegion
	for _, p := range plan {
		regions = append(regions, undoRegion{p.rec.headerStart, int64(len(p.newBytes))})
	}
	err := f.journaledWrite(regions, false, func(fh *os.File) error {
		for _, p := range plan {
			if err := ioFault("write"); err != nil {
				return err
			}
			if _, err := fh.WriteAt(p.newBytes, p.rec.headerStart); err != nil {
				return err
			}
			if err := faultPoint("header-written"); err != nil {
				return err
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	for _, p := range plan {
		p.rec.adoptHeader(p.newBytes)
	}
	return nil
}

// undoRegion is a byte range [Off, Off+Len) of the FITS file that a
// journaled write may change.
type undoRegion struct{ Off, Len int64 }

// journaledWrite runs write (which changes only the given regions of the
// file, and grows the file only if grows is set) under the undo journal:
//
//  1. the current bytes of every region (clipped to the file size) are
//     saved to <name>.journal with the file size, and the journal and its
//     directory are fsynced ("journal-written");
//  2. write runs;
//  3. the file is fsynced ("data-synced");
//  4. the journal is removed and the directory fsynced ("journal-removed"):
//     the commit point.
//
// A process crash before step 4 leaves the journal, and the next
// OpenForEdit or EditFile puts back the saved bytes (and cuts the file
// back to its old size when grows is set), so the file is byte-identical
// to its state before the call. A real error from write, the fsync or the
// removal rolls the file back the same way before it is returned. An
// injected crash (faultPoint) returns at once with no cleanup.
func (f *File) journaledWrite(regions []undoRegion, grows bool, write func(fh *os.File) error) error {
	if f.name == "" {
		return fmt.Errorf("fits: a journaled write requires a file opened with OpenForEdit")
	}
	target, err := resolveTarget(f.name)
	if err != nil {
		return err
	}
	fh, ok := f.closer.(*os.File)
	if !ok {
		return fmt.Errorf("fits: a journaled write requires a file opened with OpenForEdit")
	}
	// No reader may be using the bytes this write changes; the lock is
	// held until the journal is gone (a process death drops it).
	if err := lockForWrite(fh, false); err != nil {
		return err
	}
	defer unlockWrite(fh)
	st, err := fh.Stat()
	if err != nil {
		return err
	}
	j := undoJournal{FileSize: st.Size(), Truncate: grows}
	for _, r := range regions {
		n := min(r.Len, st.Size()-r.Off)
		if n <= 0 {
			continue // past the end: undone by the truncation
		}
		old := make([]byte, n)
		if _, err := fh.ReadAt(old, r.Off); err != nil {
			return fmt.Errorf("fits: read bytes %d..%d for the journal: %w", r.Off, r.Off+n, err)
		}
		j.Entries = append(j.Entries, undoEntry{Offset: r.Off, Data: old})
	}
	jp := journalPath(target)
	if err := writeFileDurable(jp, encodeJournal(j)); err != nil {
		os.Remove(jp)
		return err
	}
	if err := faultPoint("journal-written"); err != nil {
		return err
	}
	// From here on a real error rolls the file back from the journal.
	fail := func(err error) error {
		if f.br != nil {
			f.br.Invalidate()
		}
		// rollbackJournal takes the write lock through its own descriptor,
		// which conflicts with this one's.
		unlockWrite(fh)
		if rerr := rollbackJournal(target); rerr != nil {
			return fmt.Errorf("%w; rollback from journal failed (%v): the journal %s is kept and the next OpenForEdit retries it", err, rerr, jp)
		}
		return err
	}
	if err := write(fh); err != nil {
		if isFaultStop(err) {
			return err
		}
		return fail(err)
	}
	if err := ioFault("fsync"); err != nil {
		return fail(err)
	}
	if err := fsync(fh, "target"); err != nil {
		return fail(err)
	}
	if err := faultPoint("data-synced"); err != nil {
		return err
	}
	if err := ioFault("remove-journal"); err != nil {
		return fail(err)
	}
	if err := removeDurable(jp); err != nil {
		// The new bytes are on disk but the journal may still exist, so the
		// next OpenForEdit would roll them back. Roll back now so the on-disk
		// state matches the error we report.
		return fail(err)
	}
	if err := faultPoint("journal-removed"); err != nil {
		return err
	}
	if f.br != nil {
		f.br.Invalidate()
	}
	return nil
}

// journaledAppend runs emit, which writes one complete HDU at the end of
// the file, under the undo journal when the file is in ModeEdit (see
// journaledWrite: a crash or error leaves the file exactly as it was,
// including any trailing bytes the new HDU overwrites). In ModeCreate it
// just runs emit. start is where the HDU will be written.
func (f *File) journaledAppend(start int64, emit func() (*hduRecord, error)) (*hduRecord, error) {
	if f.mode != ModeEdit {
		return emit()
	}
	var rec *hduRecord
	err := f.journaledWrite([]undoRegion{{start, math.MaxInt64 - start}}, true, func(*os.File) error {
		if err := ioFault("write"); err != nil {
			return err
		}
		r, err := emit()
		if err != nil {
			return err
		}
		rec = r
		return faultPoint("hdu-written")
	})
	if err != nil {
		if rec != nil && len(f.hdus) > 0 && f.hdus[len(f.hdus)-1] == rec {
			f.hdus = f.hdus[:len(f.hdus)-1] // rolled back: the HDU is not in the file
		}
		return nil, err
	}
	return rec, nil
}

// rewriteWithHeaders writes the whole file again, with the planned headers,
// to a temp file and renames it over the target.
func (f *File) rewriteWithHeaders(plan []pendingHeader) error {
	target, err := resolveTarget(f.name)
	if err != nil {
		return err
	}
	st, err := os.Stat(target)
	if err != nil {
		return err
	}
	if err := checkReplaceable(target, st); err != nil {
		return err
	}
	newHdr := map[*hduRecord]pendingHeader{}
	for _, p := range plan {
		newHdr[p.rec] = p
	}
	tmpPath := target + tempSuffix
	if err := os.Remove(tmpPath); err != nil && !errors.Is(err, os.ErrNotExist) {
		return err
	}
	tmp, err := os.OpenFile(tmpPath, os.O_RDWR|os.O_CREATE|os.O_EXCL, st.Mode().Perm())
	if err != nil {
		return err
	}
	cleanup := func(err error) error {
		tmp.Close()
		os.Remove(tmpPath)
		return err
	}
	if err := tmp.Chmod(st.Mode().Perm()); err != nil {
		return cleanup(err)
	}
	if err := carryOver(target, st, tmp, tmpPath); err != nil {
		return cleanup(err)
	}
	// The temp file becomes the file after the rename: lock it now so the
	// edit lock moves with it, with no moment when the path is unlocked.
	if err := lockForEdit(tmp); err != nil {
		return cleanup(err)
	}
	// Write every HDU and compute the new layout.
	type layout struct{ headerStart, dataStart, dataEnd, paddedEnd int64 }
	lay := make([]layout, len(f.hdus))
	var pos int64
	buf := make([]byte, 1<<20)
	// copyRange copies the stored bytes [from, to) of the file.
	copyRange := func(from, to int64, what string) error {
		for done := from; done < to; {
			n := min(int64(len(buf)), to-done)
			if err := f.br.ReadRange(done, buf[:n]); err != nil {
				return fmt.Errorf("fits: read %s: %w", what, err)
			}
			if _, err := tmp.Write(buf[:n]); err != nil {
				return err
			}
			done += n
		}
		return nil
	}
	// Only the bytes the file has are copied: a last data unit that lacks
	// its final padding stays without it, and the bytes after the last
	// HDU (TrailingBytes) follow it unchanged, as after an in-place Flush.
	total := st.Size() // now: appends in ModeEdit grow the file after Open
	for i, rec := range f.hdus {
		hdr := rec.rawHeader
		if p, ok := newHdr[rec]; ok {
			hdr = p.newBytes
		}
		if _, err := tmp.Write(hdr); err != nil {
			return cleanup(err)
		}
		dataLen := rec.paddedEnd - rec.dataStart
		ds := pos + int64(len(hdr))
		lay[i] = layout{pos, ds, ds + (rec.dataEnd - rec.dataStart), ds + dataLen}
		if err := copyRange(rec.dataStart, min(rec.paddedEnd, total), fmt.Sprintf("HDU %d data", rec.index)); err != nil {
			return cleanup(err)
		}
		pos = lay[i].paddedEnd
	}
	if n := len(f.hdus); n > 0 {
		if err := copyRange(f.hdus[n-1].paddedEnd, total, "the bytes after the last HDU"); err != nil {
			return cleanup(err)
		}
	}
	if err := faultPoint("temp-written"); err != nil {
		tmp.Close() // a process death closes its descriptors (and drops the lock)
		return err
	}
	if err := fsync(tmp, "temp"); err != nil {
		return cleanup(err)
	}
	if err := faultPoint("temp-synced"); err != nil {
		tmp.Close() // a process death closes its descriptors (and drops the lock)
		return err
	}
	if err := os.Rename(tmpPath, target); err != nil {
		return cleanup(err)
	}
	// Committed. Errors after this point are reported, but the file on disk
	// is the complete edited file.
	if err := faultPoint("renamed"); err != nil {
		tmp.Close() // a process death closes its descriptors (and drops the lock)
		return err
	}
	if err := syncDir(target); err != nil {
		tmp.Close()
		return fmt.Errorf("fits: header resize committed but directory fsync failed: %w", err)
	}
	if err := faultPoint("dir-synced"); err != nil {
		tmp.Close() // a process death closes its descriptors (and drops the lock)
		return err
	}
	// The old descriptor refers to the replaced inode: switch to the temp
	// file's descriptor (the new file, holding the edit lock) and adopt the
	// new layout.
	fh := tmp
	if _, err := fh.Seek(0, io.SeekStart); err != nil {
		fh.Close()
		return fmt.Errorf("fits: header resize committed but the new file cannot be used: %w", err)
	}
	br, err := block.NewReader(fh)
	if err != nil {
		fh.Close()
		return err
	}
	bw, err := block.NewWriter(fh)
	if err != nil {
		fh.Close()
		return err
	}
	if f.closer != nil {
		f.closer.Close()
	}
	f.closer, f.rs, f.rws, f.br, f.bw = fh, fh, fh, br, bw
	for i, rec := range f.hdus {
		if p, ok := newHdr[rec]; ok {
			rec.adoptHeader(p.newBytes)
		}
		rec.headerStart, rec.dataStart, rec.dataEnd, rec.paddedEnd = lay[i].headerStart, lay[i].dataStart, lay[i].dataEnd, lay[i].paddedEnd
	}
	return nil
}

// adoptHeader records that newBytes is now the header on disk. The
// structural keywords cannot have changed (planHeader checks them); the
// EXTNAME HDUByName compares is read from the header itself
// (currentExtname), so nothing else in the record follows the header.
func (r *hduRecord) adoptHeader(newBytes []byte) {
	r.rawHeader = append([]byte(nil), newBytes...)
	r.dirty = false
}

// flushFaultHook, when non-nil, is called at every named step of Flush, of
// journal recovery and of EditFile. Tests set it to simulate a crash at that
// step: a non-nil return makes the operation stop at once with no cleanup,
// leaving the on-disk state exactly as a process death at that point would.
var flushFaultHook func(step string) error

func faultPoint(step string) error {
	if flushFaultHook == nil {
		return nil
	}
	if err := flushFaultHook(step); err != nil {
		return faultStop{err}
	}
	return nil
}

// faultStop marks an injected crash: journaledWrite returns it without
// rolling back, as a process death would.
type faultStop struct{ err error }

func (e faultStop) Error() string { return e.err.Error() }
func (e faultStop) Unwrap() error { return e.err }

func isFaultStop(err error) bool {
	var fs faultStop
	return errors.As(err, &fs)
}

// ioFaultHook, when non-nil, is called before each I/O step of a journaled
// write ("write", "fsync", "remove-journal"). Tests return an error to
// simulate a real I/O failure (a full disk, an EIO): unlike an injected
// crash it must be rolled back from the journal before it is returned.
var ioFaultHook func(op string) error

func ioFault(op string) error {
	if ioFaultHook == nil {
		return nil
	}
	return ioFaultHook(op)
}
