package fits

import (
	"bytes"
	"context"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
	"github.com/dmarkham/fits/internal/block"
)

// EditFile provides a streaming-rebuild surface for structural changes:
// dropping, reordering, duplicating or inserting HDUs, and replacing a
// header while keeping its data bytes.
//
// It first recovers path from any crashed earlier Flush or EditFile (as
// OpenForEdit does), opens it read-only, and passes the caller an (in,
// out) pair. The output goes to a temp file (<path>.fits-tmp) next to the
// target, created with the source's permission bits. When the callback
// returns nil, the temp file is fsynced and atomically renamed over the
// source and the directory is fsynced. When the callback returns an error
// or panics, or writes no HDU at all, the temp file is removed and the
// source is left byte-identical.
//
// EditFile holds the same edit lock as OpenForEdit for the whole call: it
// fails with ErrLocked while another editor has the file, and another
// editor is refused while it runs. Like OpenForEdit it is supported on
// Linux and macOS (ErrEditUnsupported elsewhere). On macOS the edit lock
// is a record lock for writing, which needs the file open for writing, so
// EditFile needs write permission on the file (on Linux, permission to
// write the directory is enough).
//
// Crash safety: the rename is the commit point. A process crash before it
// leaves the source untouched (the stale temp file is removed by the next
// OpenForEdit or EditFile); a crash after it leaves the complete output.
// Power loss additionally relies on the fsyncs above, whose order is
// checked by the tests but which no test exercises with a real power cut.
// A symbolic link at path is followed and kept. The rename replaces the
// directory entry, so a file with other hard links is refused
// (ErrHardLinked), and the new file is given the source's permission bits
// and, on Unix, its owner and group, and its extended attributes: on Linux
// POSIX ACLs included, on macOS the resource fork, Finder info and other
// attributes but NOT an ACL (macOS keeps it outside the listed
// attributes; the new file has only the ACL entries it inherits from its
// directory, see xattr_darwin.go); if one of them cannot be set the edit
// is refused with ErrCannotPreserve and the source is left as it was
// (TestReplaceKeepsGroupXattrsAndACL, TestReplaceKeepsGroupAndXattrsDarwin,
// TestReplaceRefusesHardLinks, TestReplaceRefusesWhenOwnerCannotBeKept).
func EditFile(path string, fn func(in *File, out *Writer) error) error {
	return EditFileContext(context.Background(), path, fn)
}

// EditFileContext is EditFile with an explicit context. It is checked
// before anything is touched (no recovery of an interrupted edit, no temp
// file, fn not called) and once per HDU while the source is scanned; fn
// itself should check ctx for long edits (TestContextFunctions).
func EditFileContext(ctx context.Context, path string, fn func(in *File, out *Writer) error) error {
	if err := ctx.Err(); err != nil {
		return err
	}
	// Hold the edit lock for the whole edit, taken before recovery so a
	// running editor's journal is never rolled back under it.
	lk, err := openLocked(path, os.O_RDONLY)
	if err != nil {
		return err
	}
	defer lk.Close()
	if err := recoverEdit(path); err != nil {
		return fmt.Errorf("fits: recovering %s from an interrupted edit: %w", path, err)
	}
	target, err := resolveTarget(path)
	if err != nil {
		return err
	}
	in, err := OpenContext(ctx, target)
	if err != nil {
		return err
	}
	defer in.Close()
	st, err := os.Stat(target)
	if err != nil {
		return err
	}
	if err := checkReplaceable(target, st); err != nil {
		return err
	}

	tmpPath := target + tempSuffix
	tmpFile, err := os.OpenFile(tmpPath, os.O_RDWR|os.O_CREATE|os.O_EXCL, st.Mode().Perm())
	if err != nil {
		return err
	}
	// cleanupTmp is called on every error path.
	cleanupTmp := func() {
		tmpFile.Close()
		os.Remove(tmpPath)
	}
	if err := tmpFile.Chmod(st.Mode().Perm()); err != nil {
		cleanupTmp()
		return err
	}
	if err := carryOver(target, st, tmpFile, tmpPath); err != nil {
		cleanupTmp()
		return err
	}
	// The temp file becomes the file at the rename: lock it now so the
	// path is never unlocked while this edit runs.
	if err := lockForEdit(tmpFile); err != nil {
		cleanupTmp()
		return err
	}
	bw, err := block.NewWriter(tmpFile)
	if err != nil {
		cleanupTmp()
		return err
	}
	w := &Writer{out: bw}

	// Invoke the callback. Any error or panic rolls back.
	runErr := func() (ret error) {
		defer func() {
			if r := recover(); r != nil {
				ret = fmt.Errorf("fits: EditFile callback panic: %v", r)
			}
		}()
		return fn(in, w)
	}()
	if runErr == nil && w.hduIndex == 0 {
		runErr = fmt.Errorf("fits: EditFile: the callback wrote no HDU; refusing to replace %s with an empty file", path)
	}
	if runErr != nil {
		cleanupTmp()
		return runErr
	}
	if err := faultPoint("temp-written"); err != nil {
		tmpFile.Close() // a process death closes its descriptors (and drops the lock)
		return err
	}
	if err := fsync(tmpFile, "temp"); err != nil {
		cleanupTmp()
		return err
	}
	if err := faultPoint("temp-synced"); err != nil {
		tmpFile.Close() // a process death closes its descriptors (and drops the lock)
		return err
	}
	if err := os.Rename(tmpPath, target); err != nil {
		cleanupTmp()
		return err
	}
	// Committed. The temp file's descriptor (the new file, holding the edit
	// lock) is closed last.
	defer tmpFile.Close()
	if err := faultPoint("renamed"); err != nil {
		return err
	}
	// fsync the containing directory so the rename is durable.
	if err := syncDir(target); err != nil {
		return fmt.Errorf("fits: EditFile committed but directory fsync failed: %w", err)
	}
	if err := faultPoint("dir-synced"); err != nil {
		return err
	}
	return tmpFile.Close()
}

// Writer is the streaming output side of EditFile. Callers emit HDUs to it,
// in output order, with CopyHDU (pass-through), CopyHDUWithHeader (new
// header, same data bytes) and WriteImage(out, ...) (a new image HDU).
// SkipHDU records that an input HDU is left out (it writes nothing).
// Binary tables cannot be created through a Writer (AppendBinaryTable
// takes a *File).
//
// The first HDU of a FITS file must be a primary array (SIMPLE) and every
// later one an extension (XTENSION), so an HDU whose position changes is
// converted exactly as cfitsio's fits_copy_hdu converts it:
//
//   - A primary array written at position 1 or later becomes an IMAGE
//     extension: SIMPLE is replaced by XTENSION = 'IMAGE', PCOUNT = 0 and
//     GCOUNT = 1 are inserted after the NAXISn cards, and EXTEND and the
//     two standard "FITS (Flexible Image Transport System) format ..."
//     COMMENT cards are dropped. A random-groups primary (GROUPS = T)
//     cannot be an extension and is refused.
//   - An IMAGE extension written at position 0 becomes the primary array:
//     XTENSION is replaced by SIMPLE = T, EXTEND = T and the two standard
//     COMMENT cards are inserted after the NAXISn cards, and PCOUNT and
//     GCOUNT are dropped.
//   - A table (including a tile-compressed image) written at position 0 is
//     preceded by an empty primary array (BITPIX = 8, NAXIS = 0).
//
// A converted header is rebuilt from its keywords: reserved blank cards
// between the last keyword and END in the input are not carried over (nor
// does fits_copy_hdu copy them), and no new reserve is added, so the header
// can take fewer blocks than cfitsio's, which sometimes ends with blank
// fill. Headers that keep their form are copied byte-for-byte, reserved
// space included.
type Writer struct {
	out      *block.Writer
	hduIndex int // number of HDUs written so far
}

// recordOf returns the backing record of any HDU type, including
// tile-compressed images.
func recordOf(h HDU) *hduRecord {
	if c, ok := h.(*CompressedImageHDU); ok {
		return c.rec
	}
	return hduRecordOf(h)
}

// CopyHDU writes h to the output: its on-disk header bytes (converted
// between primary and extension form if its position changes, see Writer)
// and its data bytes verbatim. In-memory header changes not yet flushed are
// not copied; use CopyHDUWithHeader for those.
func (w *Writer) CopyHDU(h HDU) error {
	rec := recordOf(h)
	if rec == nil || rec.file == nil {
		return fmt.Errorf("fits: CopyHDU: HDU has no backing record")
	}
	return w.copyFromRecord(rec, rec.rawHeader)
}

// CopyHDUWithHeader writes h to the output with its header replaced by hdr
// and its data bytes passed through verbatim. hdr must describe the same
// data: if the data size implied by hdr's structural keywords (BITPIX,
// NAXISn, PCOUNT, GCOUNT) differs from h's, or hdr changes the HDU kind
// (image or which table type), it is refused. hdr is converted between
// primary and extension form like CopyHDU. A header the encoder cannot
// represent exactly is refused, as in Flush: every card must read back
// with the key, value (and value type) and comment that was set. Today
// commentary (COMMENT/HISTORY) text longer than one card is refused by
// this check, and NaN and a keyword named END are refused earlier by
// header.Set/Add; EditFile then fails and leaves the source
// byte-identical. CONTINUE long strings and integer-valued floats are
// written (the float with a decimal point) and read back by cfitsio
// (TestEditFileCopyHDUWithHeaderRefusesLossyEncoding pins each outcome).
func (w *Writer) CopyHDUWithHeader(h HDU, hdr *header.Header) error {
	rec := recordOf(h)
	if rec == nil || rec.file == nil {
		return fmt.Errorf("fits: CopyHDUWithHeader: HDU has no backing record")
	}
	b, err := encodeVerified(hdr, rec.index)
	if err != nil {
		return fmt.Errorf("fits: CopyHDUWithHeader: %w", err)
	}
	newRec, err := extractStructural(b, rec.index)
	if err != nil {
		return fmt.Errorf("fits: CopyHDUWithHeader: %w", err)
	}
	if got, want := dataSize(newRec), dataSize(rec); got != want {
		return fmt.Errorf("fits: CopyHDUWithHeader: new header describes %d data bytes, HDU %d has %d", got, rec.index, want)
	}
	isImage := func(r *hduRecord) bool { return r.xtension == "" || r.xtension == "IMAGE" }
	if isImage(newRec) != isImage(rec) || (!isImage(rec) && newRec.xtension != rec.xtension) || newRec.zimage != rec.zimage {
		return fmt.Errorf("fits: CopyHDUWithHeader: new header changes the HDU kind")
	}
	return w.copyFromRecord(rec, b)
}

// SkipHDU records that h is left out of the output. It writes nothing; any
// HDU, including the primary, may be skipped (the next HDU written becomes
// the primary, converted as described on Writer).
func (w *Writer) SkipHDU(h HDU) error {
	if recordOf(h) == nil {
		return fmt.Errorf("fits: SkipHDU: HDU has no backing record")
	}
	return nil
}

// copyFromRecord writes hdrBytes (converted for the output position) and
// the data of rec to the output.
func (w *Writer) copyFromRecord(rec *hduRecord, hdrBytes []byte) error {
	hdrBytes, needDummy, err := convertForPosition(hdrBytes, w.hduIndex)
	if err != nil {
		return err
	}
	if needDummy {
		if err := w.out.WriteRange(emptyPrimaryHeader()); err != nil {
			return err
		}
		w.hduIndex++
	}
	if err := w.out.WriteRange(hdrBytes); err != nil {
		return err
	}

	// Data: pass through byte-for-byte (including trailing pad).
	dataLen := rec.paddedEnd - rec.dataStart
	if dataLen > 0 {
		const chunk = 1 << 20 // 1 MiB
		buf := make([]byte, chunk)
		var done int64
		for done < dataLen {
			n := min(int64(chunk), dataLen-done)
			if err := rec.file.br.ReadRange(rec.dataStart+done, buf[:n]); err != nil {
				return err
			}
			if _, err := w.out.Write(buf[:n]); err != nil {
				return err
			}
			done += n
		}
	}
	w.hduIndex++
	return nil
}

// Standard COMMENT cards written by cfitsio (putkey.c ffphprll, edithdu.c
// ffcphd). Converted headers use the same text so they match cfitsio's
// fits_copy_hdu output card for card.
const (
	cardFITSComment1 = "COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronomy"
	cardFITSComment2 = "COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H"
)

// fixedCard formats a fixed-format value card: keyword in columns 1-8,
// "= " in 9-10, the value right-justified to column 30 (a string
// left-justified from column 11 and padded to at least 8 characters),
// then " / comment".
func fixedCard(key, value string, isString bool, comment string) string {
	var v string
	if isString {
		s := value
		for len(s) < 8 {
			s += " "
		}
		v = fmt.Sprintf("%-20s", "'"+s+"'")
	} else {
		v = fmt.Sprintf("%20s", value)
	}
	return fmt.Sprintf("%-8s= %s / %s", key, v, comment)
}

// emptyPrimaryHeader is the header cfitsio writes for an empty primary
// array (fits_create_img with BITPIX = 8, NAXIS = 0).
func emptyPrimaryHeader() []byte {
	return packHeader([]string{
		fixedCard("SIMPLE", "T", false, "file does conform to FITS standard"),
		fixedCard("BITPIX", "8", false, "number of bits per data pixel"),
		fixedCard("NAXIS", "0", false, "number of data axes"),
		fixedCard("EXTEND", "T", false, "FITS dataset may contain extensions"),
		cardFITSComment1,
		cardFITSComment2,
	})
}

// packHeader pads each card to 80 columns, appends END and pads to a whole
// number of 2880-byte blocks with blanks.
func packHeader(cards []string) []byte {
	var b bytes.Buffer
	for _, c := range append(cards, "END") {
		b.WriteString(c)
		b.Write(bytes.Repeat([]byte(" "), header.CardWidth-len(c)))
	}
	if r := b.Len() % block.Size; r != 0 {
		b.Write(bytes.Repeat([]byte(" "), block.Size-r))
	}
	return b.Bytes()
}

// splitCards returns the 80-column cards before END.
func splitCards(hdr []byte) ([]string, error) {
	var cards []string
	for off := 0; off+header.CardWidth <= len(hdr); off += header.CardWidth {
		c := string(hdr[off : off+header.CardWidth])
		if strings.TrimRight(c, " ") == "END" {
			return cards, nil
		}
		cards = append(cards, c)
	}
	return nil, fmt.Errorf("fits: header has no END card")
}

func cardKey(c string) string { return strings.TrimRight(c[:8], " ") }

// cardValue returns the value field of a fixed-format value card (from
// column 11 up to a comment slash outside quotes), trimmed.
func cardValue(c string) string {
	if len(c) < 10 || c[8:10] != "= " {
		return ""
	}
	v := c[10:]
	inStr := false
	for i := 0; i < len(v); i++ {
		switch v[i] {
		case '\'':
			inStr = !inStr
		case '/':
			if !inStr {
				return strings.TrimSpace(v[:i])
			}
		}
	}
	return strings.TrimSpace(v)
}

// convertForPosition converts hdr between primary and extension form for
// output position pos (see Writer). needDummy reports that an empty primary
// array must be written first.
func convertForPosition(hdr []byte, pos int) (out []byte, needDummy bool, err error) {
	if len(hdr) < header.CardWidth {
		return nil, false, fmt.Errorf("fits: header too short")
	}
	first := cardKey(string(hdr[:header.CardWidth]))
	isPrimary := first == "SIMPLE"
	if !isPrimary && first != "XTENSION" {
		return nil, false, fmt.Errorf("fits: header starts with %q, not SIMPLE or XTENSION", first)
	}
	if (pos == 0) == isPrimary {
		return hdr, false, nil // already in the right form
	}
	cards, err := splitCards(hdr)
	if err != nil {
		return nil, false, err
	}
	if pos == 0 && strings.Trim(cardValue(cards[0]), "' ") != "IMAGE" {
		// A table cannot be the primary HDU: put an empty primary in front.
		return hdr, true, nil
	}
	// Trailing blank cards (reserved space) are not keywords.
	for len(cards) > 0 && strings.TrimSpace(cards[len(cards)-1]) == "" {
		cards = cards[:len(cards)-1]
	}
	// The mandatory cards: SIMPLE/XTENSION, BITPIX, NAXIS, NAXIS1..n.
	if len(cards) < 3 || cardKey(cards[1]) != "BITPIX" || cardKey(cards[2]) != "NAXIS" {
		return nil, false, fmt.Errorf("fits: header does not start with the mandatory BITPIX, NAXIS cards")
	}
	naxis, err := strconv.Atoi(cardValue(cards[2]))
	if err != nil || naxis < 0 || naxis > 999 || len(cards) < 3+naxis {
		return nil, false, fmt.Errorf("fits: bad NAXIS card %q", strings.TrimRight(cards[2], " "))
	}
	for i := 1; i <= naxis; i++ {
		if cardKey(cards[2+i]) != "NAXIS"+strconv.Itoa(i) {
			return nil, false, fmt.Errorf("fits: card %d is %q, want NAXIS%d", 3+i, cardKey(cards[2+i]), i)
		}
	}
	mandatory, rest := cards[1:3+naxis], cards[3+naxis:]
	var outCards []string
	if isPrimary {
		// Primary array -> IMAGE extension.
		for _, c := range rest {
			if cardKey(c) == "GROUPS" && cardValue(c) == "T" {
				return nil, false, fmt.Errorf("fits: a random-groups primary array cannot be written as an extension")
			}
		}
		outCards = append(outCards, fixedCard("XTENSION", "IMAGE", true, "IMAGE extension"))
		outCards = append(outCards, mandatory...)
		outCards = append(outCards,
			fixedCard("PCOUNT", "0", false, "number of random group parameters"),
			fixedCard("GCOUNT", "1", false, "number of random groups"))
		for _, c := range rest {
			// The same prefixes cfitsio compares (58 and 47 characters).
			if cardKey(c) == "EXTEND" || strings.HasPrefix(c, cardFITSComment1[:58]) || strings.HasPrefix(c, cardFITSComment2[:47]) {
				continue
			}
			outCards = append(outCards, c)
		}
	} else {
		// IMAGE extension -> primary array.
		outCards = append(outCards, fixedCard("SIMPLE", "T", false, "file does conform to FITS standard"))
		outCards = append(outCards, mandatory...)
		outCards = append(outCards,
			fixedCard("EXTEND", "T", false, "FITS dataset may contain extensions"),
			cardFITSComment1, cardFITSComment2)
		for _, c := range rest {
			if k := cardKey(c); k == "PCOUNT" || k == "GCOUNT" {
				continue
			}
			outCards = append(outCards, c)
		}
	}
	for i, c := range outCards {
		outCards[i] = strings.TrimRight(c, " ")
	}
	return packHeader(outCards), false, nil
}

// writeImage on *Writer: append a new image HDU. Implements WriteTarget.
func (w *Writer) writeImage(bp bitpix.BITPIX, hdr *header.Header, shape []int64, data anyData) (*ImageHDU, error) {
	rec, err := emitImageHDU(w.out, bp, hdr, shape, data, w.hduIndex == 0, w.hduIndex)
	if err != nil {
		return nil, err
	}
	w.hduIndex++
	return &ImageHDU{rec: rec}, nil
}
