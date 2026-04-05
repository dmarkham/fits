package fits

import (
	"context"
	"fmt"
	"io"
	"os"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/block"
)

// OpenForEdit opens a FITS file for in-place mutation.
//
// The returned *File supports header keyword update/add/delete (via the
// Header() surface), same-shape pixel overwrite (OverwritePixels), and HDU
// append at EOF (AppendImage/AppendBinaryTable). Structural mutations
// (resize, insert-in-middle, delete) must be done through EditFile.
//
// On Open we check for a stale journal file from a crashed previous edit;
// if present it is either replayed or rolled back before returning. See
// plan decision 5.8 and implementation step 27.
func OpenForEdit(name string) (*File, error) {
	return OpenForEditContext(context.Background(), name)
}

// OpenForEditContext is OpenForEdit with an explicit context.
func OpenForEditContext(ctx context.Context, name string) (*File, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	// Check for a stale journal from a crashed previous edit.
	if err := recoverJournal(name); err != nil {
		return nil, fmt.Errorf("fits: journal recovery: %w", err)
	}

	fh, err := os.OpenFile(name, os.O_RDWR, 0)
	if err != nil {
		return nil, err
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

// Flush persists any pending header mutations. Only meaningful in ModeEdit
// and ModeCreate. Returns ErrReadOnly for read-only files.
//
// Implementation strategy (plan decision 5.8):
//
//  1. For each dirty HDU, re-serialize the full header into bytes.
//  2. If new length == old length (block count unchanged), overwrite in
//     place and fsync. O(KB) per HDU.
//  3. If new length differs, walk the dirty HDUs in reverse order and shift
//     the file tail using a journal protocol (step 27). For v1 we
//     conservatively refuse any flush that would require a tail shift and
//     point the caller at EditFile. Journaled shift is implemented below
//     but gated behind explicit opt-in while the crash-safety tests mature.
func (f *File) Flush() error {
	if f.mode == ModeRead {
		return ErrReadOnly
	}
	if f.mode == ModeCreate {
		// Nothing to flush for ModeCreate — writes go straight to disk.
		return nil
	}

	// Walk every HDU with a loaded parsed header, re-encode, and compare.
	// If the encoding changed, rewrite. This subsumes an explicit dirty flag
	// and catches mutations performed directly on *Header returned from
	// Header() — no extra bookkeeping needed.
	for _, rec := range f.hdus {
		if rec.parsed == nil {
			continue
		}
		newBytes, err := header.Encode(rec.parsed)
		if err != nil {
			return fmt.Errorf("fits: re-encode HDU %d header: %w", rec.index, err)
		}
		if bytesEqual(newBytes, rec.rawHeader) {
			continue
		}
		oldLen := int64(len(rec.rawHeader))
		if int64(len(newBytes)) == oldLen {
			if err := f.flushInPlace(rec, newBytes); err != nil {
				return err
			}
			continue
		}
		// Header grew or shrank across block boundary — journaled tail shift.
		if err := f.flushWithTailShift(rec, newBytes); err != nil {
			return err
		}
	}
	// Sync the file to disk if we own it.
	if syncer, ok := f.closer.(interface{ Sync() error }); ok {
		if err := syncer.Sync(); err != nil {
			return err
		}
	}
	return nil
}

// flushInPlace rewrites a header whose new length matches the old. Fast path.
func (f *File) flushInPlace(rec *hduRecord, newBytes []byte) error {
	if _, err := f.bw.Seek(rec.headerStart, io.SeekStart); err != nil {
		return err
	}
	if err := f.bw.WriteRange(newBytes); err != nil {
		return err
	}
	// Update cached raw header.
	rec.rawHeader = append(rec.rawHeader[:0], newBytes...)
	rec.dirty = false
	// Reader cache must be invalidated — we just wrote through.
	f.br.Invalidate()
	return nil
}

// Close finalizes the file. In write modes it flushes pending header
// mutations and syncs. In all modes it releases any owned file handle.
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
