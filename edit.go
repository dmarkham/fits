package fits

import (
	"context"
	"errors"
	"fmt"
	"math"
	"os"
	"reflect"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/block"
)

// OpenForEdit opens a FITS file for in-place mutation.
//
// The returned *File supports header keyword update/add/delete through the
// Header() surface (persisted by Flush or CloseAndFlush; Close alone
// discards them), same-shape pixel overwrite (OverwritePixels), and HDU
// append at EOF (AppendImage/AppendBinaryTable). Header edits must not
// change the structural keywords that size or type the data (SIMPLE,
// XTENSION, BITPIX, NAXIS, NAXISn, PCOUNT, GCOUNT, TFIELDS, ZIMAGE):
// Flush refuses them. Structural mutations (resize, insert-in-middle,
// delete) must be done through EditFile.
//
// Editing is supported on Linux and macOS; elsewhere (Windows, the BSDs)
// OpenForEdit and EditFile fail with ErrEditUnsupported, because the
// locks below and the carry-over of extended attributes are implemented
// for Linux and macOS only (lock_ofd.go, lock_linux.go, lock_darwin.go).
//
// One editor at a time: OpenForEdit takes an exclusive advisory lock
// (flock(2) on Linux, an open file description record lock on macOS) on
// the file that is held until Close and moves with the file across a
// resizing Flush; a second OpenForEdit or EditFile of the same
// file, in this process or another, fails with ErrLocked and touches
// nothing (TestEditLockRefusesConcurrentEditors), also when it opened the
// file just before the first editor's rename and locked it just after
// (TestEditLockFollowsReplacedFile).
//
// Readers: a read-only Open (or OpenReader of an *os.File) of the file,
// in this process or another, keeps the editor's in-place writes (a Flush
// that resizes no header, OverwritePixels, the appends) from running: they
// fail with ErrOpenForReading and write nothing until the readers are
// closed (TestInPlaceWriteRefusedWhileReaderOpen). A resizing Flush is
// allowed, since it replaces the file and leaves the readers on the old,
// complete one. Both locks are advisory: programs that do not use this
// package are not stopped by them.
//
// Before opening, OpenForEdit recovers from a crashed earlier Flush or
// EditFile on the same file: an undo journal (<name>.journal) is rolled
// back, restoring the exact bytes the file had before that Flush, and a
// stale temp file (<name>.fits-tmp) is removed. A journal left by an
// earlier release of this package cannot be recovered and makes
// OpenForEdit fail with ErrLegacyJournal. See journal.go for the protocol.
//
// Every change made through the returned *File is crash safe: Flush,
// OverwritePixels and the appends (AppendImage / WriteImage,
// AppendBinaryTable, AppendASCIITable, AppendCompressedImage,
// AppendCompressedFloat32Image) each run under the undo journal, so a
// process crash during one of them is rolled back to the bytes the file
// had before that call by the next OpenForEdit or EditFile
// (TestInPlaceWriteCrashAtEveryStep), and a real I/O error is rolled back
// before it is returned (TestInPlaceWriteIOErrorRollsBack). Each call is
// its own commit: an append that has returned survives a later crash.
func OpenForEdit(name string) (*File, error) {
	return OpenForEditContext(context.Background(), name)
}

// OpenForEditContext is OpenForEdit with an explicit context, checked
// before anything is touched (no journal recovery runs under a cancelled
// context) and once per HDU of the scan (TestContextFunctions).
func OpenForEditContext(ctx context.Context, name string) (*File, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	// Lock first: recovery must not roll back the journal of an editor
	// that is still running.
	fh, err := openLocked(name, os.O_RDWR)
	if err != nil {
		return nil, err
	}
	if err := recoverEdit(name); err != nil {
		fh.Close()
		return nil, fmt.Errorf("fits: recovering %s from an interrupted edit: %w", name, err)
	}
	fi, err := fh.Stat()
	if err != nil {
		fh.Close()
		return nil, err
	}
	if fi.Size() == 0 {
		fh.Close()
		return nil, ErrEmptyFile
	}
	f := &File{mode: ModeEdit, name: name, closer: fh, rs: fh, rws: fh}
	br, err := block.NewReader(fh)
	if err != nil {
		fh.Close()
		return nil, err
	}
	f.br = br
	bw, err := block.NewWriter(fh)
	if err != nil {
		fh.Close()
		return nil, err
	}
	f.bw = bw
	if err := f.scanHDUs(ctx); err != nil {
		fh.Close()
		return nil, err
	}
	return f, nil
}

// Flush persists pending header mutations made through Header(). Only
// meaningful in ModeEdit; a no-op in ModeCreate (writes go straight to
// disk) and ErrReadOnly for read-only files.
//
// A header whose cards are unchanged is not rewritten, so a Flush with no
// mutation leaves the file byte-identical. For each changed header Flush
// re-encodes the cards and checks the encoding before anything is written:
// it re-parses the new bytes and fails (writing nothing) if any card's key,
// value (including its type: a float must not come back as an integer) or
// comment would not read back as set, if a card is lost (a keyword named
// END), if a HIERARCH card's path is itself a standard keyword name
// ("HIERARCH BSCALE", "HIERARCH END": cfitsio and astropy read it as that
// keyword, TestFlushRefusesHierarchAliasOfStandardKeyword), or if a
// structural keyword changed. Then:
//
//   - If every changed header keeps its block count, the headers are
//     overwritten in place under an undo journal.
//   - If any header needs more or fewer 2880-byte blocks, the whole file is
//     rewritten to a temp file next to it and renamed over it. This costs a
//     full copy of the file. Every byte after the changed headers is
//     copied as stored, as an in-place overwrite leaves it: a last data
//     unit that lacks its final padding stays so, and the bytes after the
//     last HDU (TrailingBytes) are kept (TestFlushKeepsBytesAfterLastHDU).
//
// All changed headers of one Flush are committed together. If the process
// dies at any point, the next OpenForEdit or EditFile of the file leaves it
// byte-identical either to its state before Flush or to the complete
// result; until that recovery has run, Open refuses with ErrPendingJournal
// when a journal is pending. Power loss relies on the fsync order (see
// journal.go).
func (f *File) Flush() error {
	if f.mode == ModeRead {
		return ErrReadOnly
	}
	if f.mode == ModeCreate {
		// Nothing to flush for ModeCreate — writes go straight to disk.
		return nil
	}
	if f.name == "" {
		return fmt.Errorf("fits: Flush requires a file opened with OpenForEdit")
	}

	var plan []pendingHeader
	resize := false
	for _, rec := range f.hdus {
		if rec.parsed == nil {
			continue
		}
		p, changed, err := planHeader(rec)
		if err != nil {
			return err
		}
		if !changed {
			continue
		}
		if len(p.newBytes) != len(rec.rawHeader) {
			resize = true
		}
		plan = append(plan, p)
	}
	if len(plan) == 0 {
		return nil
	}
	if resize {
		return f.rewriteWithHeaders(plan)
	}
	return f.overwriteHeadersJournaled(plan)
}

// planHeader decides whether rec's parsed header differs from what is on
// disk and, if so, encodes and validates its replacement.
func planHeader(rec *hduRecord) (pendingHeader, bool, error) {
	onDisk, _, err := header.ParseCards(rec.rawHeader)
	if err != nil {
		return pendingHeader{}, false, fmt.Errorf("fits: re-parse HDU %d header: %w", rec.index, err)
	}
	cur := rec.parsed.Cards()
	if sameCards(onDisk, cur) {
		return pendingHeader{}, false, nil
	}
	newBytes, err := encodeVerified(rec.parsed, rec.index)
	if err != nil {
		return pendingHeader{}, false, err
	}
	// Compare with the structure the on-disk header declares (records made
	// by the writers do not fill in every structural field).
	oldRec, err := extractStructural(rec.rawHeader, rec.index)
	if err != nil {
		return pendingHeader{}, false, fmt.Errorf("fits: HDU %d: on-disk header: %w", rec.index, err)
	}
	newRec, err := structuralOfEdited(newBytes, rec.rawHeader, rec.index)
	if err != nil {
		return pendingHeader{}, false, fmt.Errorf("fits: HDU %d: edited header: %w; structural changes need EditFile", rec.index, err)
	}
	if err := sameStructure(oldRec, newRec); err != nil {
		return pendingHeader{}, false, fmt.Errorf("fits: HDU %d: %w; structural changes need EditFile", rec.index, err)
	}
	return pendingHeader{rec: rec, newBytes: newBytes}, true, nil
}

// encodeVerified encodes h and parses the result back. Unless every card
// reads back with the key, value and comment that was set, it returns an
// error instead of the bytes: header.Encode cannot represent every card
// (commentary text longer than a card is split into more cards), and
// writing such a header would silently change or lose data. CONTINUE long
// strings and integer-valued floats read back exactly and are written;
// NaN and a keyword named END never get here (header.Set refuses them).
// The read-back is this package's parser; because cfitsio names HIERARCH
// cards differently (hierarchAlias), a HIERARCH card whose path is a
// standard keyword name is refused too, with an error wrapping
// header.ErrInvalidKeyword. idx names the HDU in the error.
func encodeVerified(h *header.Header, idx int) ([]byte, error) {
	cur := h.Cards()
	b, err := header.Encode(h)
	if err != nil {
		return nil, fmt.Errorf("fits: re-encode HDU %d header: %w", idx, err)
	}
	back, _, err := header.ParseCards(b)
	if err != nil {
		return nil, fmt.Errorf("fits: HDU %d: re-encoded header does not parse: %w", idx, err)
	}
	if len(back) != len(cur) {
		return nil, fmt.Errorf("fits: HDU %d: re-encoded header has %d cards, want %d; not writing", idx, len(back), len(cur))
	}
	for i := range cur {
		if err := refuseHierarchAlias(back[i].Key); err != nil {
			return nil, fmt.Errorf("fits: HDU %d card %d (%s): %w", idx, i+1, cur[i].Key, err)
		}
		if !sameCard(back[i], cur[i]) {
			return nil, fmt.Errorf("fits: HDU %d card %d (%s): the encoder cannot represent it (reads back as key %q value %v (%T) comment %q, set as key %q value %v (%T) comment %q); not writing",
				idx, i+1, cur[i].Key, back[i].Key, back[i].Value, back[i].Value, back[i].Comment, cur[i].Key, cur[i].Value, cur[i].Value, cur[i].Comment)
		}
	}
	return b, nil
}

// hierarchAlias reports whether key is a HIERARCH key whose path is itself
// a standard keyword name (at most 8 characters of A-Z, 0-9, '-', '_'),
// and returns that name. cfitsio's ffgknm names a HIERARCH card by the
// text between "HIERARCH " and '=' (getkey.c), so for cfitsio such a card
// is that keyword ("HIERARCH BSCALE" scales the data, "HIERARCH END" ends
// the header), and astropy names it the same; this package's parser keys
// it "HIERARCH BSCALE" and gives it no meaning.
// refuseHierarchAlias is the one rule every header write applies (Flush,
// and every writer through headerBuilder.add): a HIERARCH alias of a
// standard keyword (hierarchAlias) is refused with an error wrapping
// header.ErrInvalidKeyword.
func refuseHierarchAlias(key string) error {
	short, ok := hierarchAlias(key)
	if !ok {
		return nil
	}
	return fmt.Errorf("cfitsio names the card %q %s, the standard keyword, while this package keys it %q (and astropy ignores it as scaling), so the readers would disagree on it (cfitsio applies HIERARCH BSCALE as scaling and ends the header at HIERARCH END); set %s instead; not writing: %w",
		key, short, key, short, header.ErrInvalidKeyword)
}

func hierarchAlias(key string) (string, bool) {
	p, ok := strings.CutPrefix(key, "HIERARCH ")
	if !ok || p == "" || len(p) > 8 {
		return "", false
	}
	for i := 0; i < len(p); i++ {
		b := p[i]
		if !(b >= 'A' && b <= 'Z' || b >= 'a' && b <= 'z' || b >= '0' && b <= '9' || b == '-' || b == '_') {
			return "", false
		}
	}
	return strings.ToUpper(p), true
}

// sameCards reports whether two card lists carry the same keys, values
// and comments (the raw bytes are not compared).
func sameCards(a, b []header.Card) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !sameCard(a[i], b[i]) {
			return false
		}
	}
	return true
}

func sameCard(a, b header.Card) bool {
	if a.Key != b.Key || strings.TrimRight(a.Comment, " ") != strings.TrimRight(b.Comment, " ") {
		return false
	}
	return sameValue(a.Value, b.Value)
}

// sameValue compares card values. Trailing blanks in a string are not
// significant (FITS Standard 4.0, 4.2.1.1). Anything else must read back
// with the same type and value: a float that the encoder writes without a
// decimal point reads back as an integer, which other readers (astropy,
// cfitsio's keyword type) see as a different type, so it is not the same
// value. NaN never compares equal (FITS cannot represent it). Other floats
// compare by their bits, so a -0.0 that the encoder writes as 0.0 is not
// the same value (TestFlushEditMatrixMatchesCfitsio float-negative-zero).
func sameValue(a, b any) bool {
	if sa, ok := a.(string); ok {
		sb, ok := b.(string)
		return ok && strings.TrimRight(sa, " ") == strings.TrimRight(sb, " ")
	}
	if fa, ok := a.(float64); ok {
		fb, ok := b.(float64)
		return ok && !math.IsNaN(fa) && math.Float64bits(fa) == math.Float64bits(fb)
	}
	return reflect.DeepEqual(a, b)
}

// sameStructure reports an error if newRec's structural keywords differ
// from rec's.
func sameStructure(rec, newRec *hduRecord) error {
	switch {
	case rec.simple != newRec.simple || rec.xtension != newRec.xtension:
		return fmt.Errorf("SIMPLE/XTENSION changed")
	case rec.bitpix != newRec.bitpix:
		return fmt.Errorf("BITPIX changed from %d to %d", rec.bitpix, newRec.bitpix)
	case rec.naxis != newRec.naxis || !reflect.DeepEqual(rec.shape, newRec.shape):
		return fmt.Errorf("NAXIS/NAXISn changed from %v to %v", rec.shape, newRec.shape)
	case rec.pcount != newRec.pcount || rec.gcount != newRec.gcount:
		return fmt.Errorf("PCOUNT/GCOUNT changed")
	case rec.tfields != newRec.tfields:
		return fmt.Errorf("TFIELDS changed")
	case rec.zimage != newRec.zimage:
		return fmt.Errorf("ZIMAGE changed")
	}
	return nil
}

// structuralOfEdited is extractStructural for a header that was edited
// through header.Header, with onDisk the header it replaces. A
// *DuplicateKeywordError there has one more cause than in a file: when
// the header on disk already repeats that keyword with one value
// (accepted at Open; QFITS-written kdtree headers do), Set or Delete,
// which act on the first card of a repeated keyword, changed that card
// only. The error says so in that case, since the user never saw two
// cards; a repeat the edit itself introduced (Add) is reported as is.
func structuralOfEdited(buf, onDisk []byte, hduIdx int) (*hduRecord, error) {
	rec, err := extractStructural(buf, hduIdx)
	var dup *DuplicateKeywordError
	if errors.As(err, &dup) && cardCount(onDisk, dup.Keyword) > 1 {
		return nil, fmt.Errorf("%w (the header on disk repeats it with one value; Set and Delete change its first card only)", err)
	}
	return rec, err
}

// cardCount is the number of cards named key in raw header bytes.
func cardCount(raw []byte, key string) int {
	n := 0
	for off := 0; off+header.CardWidth <= len(raw); off += header.CardWidth {
		if strings.TrimRight(string(raw[off:off+8]), " ") == key {
			n++
		}
	}
	return n
}

// CloseAndFlush flushes pending header mutations (ModeEdit; see Flush) and
// then closes the file. Close alone never flushes. If Flush fails, its
// error is returned and the file is left open (the caller may retry or
// Close).
func (f *File) CloseAndFlush() error {
	if f.mode == ModeEdit || f.mode == ModeCreate {
		if err := f.Flush(); err != nil {
			return err
		}
	}
	return f.Close()
}

// bytesEqual is a local helper: bytes.Equal without the import here.
func bytesEqual(a, b []byte) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// currentExtname is the EXTNAME value the HDU's header holds now: the
// parsed header's, once Header() has loaded it (it may carry edits not yet
// flushed), otherwise the header on disk. With more than one EXTNAME card
// the first one counts, as for cfitsio's fits_movnam_hdu, which reads it
// with fits_read_key_str (TestHDUByNameDuplicateEXTNAMEMatchesCfitsio).
// An HDU whose first EXTNAME card has no string value, and an HDU without
// one, have no name ("").
func (r *hduRecord) currentExtname() string {
	var cards []header.Card
	if r.parsed != nil {
		cards = r.parsed.Cards()
	} else {
		c, _, err := header.ParseCards(r.rawHeader)
		if err != nil {
			return ""
		}
		cards = c
	}
	for _, c := range cards {
		if c.Key == header.KeyExtname {
			v, _ := c.Value.(string)
			return v
		}
	}
	return ""
}
