package fits

import (
	"context"
	"fmt"
	"os"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
	"github.com/dmarkham/fits/internal/block"
)

// EditFile provides a streaming-rebuild surface for structural changes
// (image resize, row/column insert/delete, HDU insert-in-middle, HDU delete,
// VLA heap reorganization).
//
// It opens src read-only, creates a temp file in the same directory, passes
// the caller an (in, writer) pair, fsyncs on success, and atomically
// renames the temp file over src. On any error (callback returns err or
// panic) the temp file is removed and src is untouched.
//
// Crash-safe by construction (temp + rename on same filesystem).
//
// See plan decision 5.8.
func EditFile(path string, fn func(in *File, out *Writer) error) error {
	return EditFileContext(context.Background(), path, fn)
}

// EditFileContext is EditFile with an explicit context.
func EditFileContext(ctx context.Context, path string, fn func(in *File, out *Writer) error) error {
	if err := ctx.Err(); err != nil {
		return err
	}
	in, err := OpenContext(ctx, path)
	if err != nil {
		return err
	}
	defer in.Close()

	tmpPath := fmt.Sprintf("%s.tmp-%d", path, os.Getpid())
	tmpFile, err := os.Create(tmpPath)
	if err != nil {
		return err
	}
	// cleanupTmp is called on every error path.
	cleanupTmp := func() {
		tmpFile.Close()
		os.Remove(tmpPath)
	}
	bw, err := block.NewWriter(tmpFile)
	if err != nil {
		cleanupTmp()
		return err
	}
	w := &Writer{out: bw, osFile: tmpFile, tmpPath: tmpPath, targetPath: path}

	// Invoke the callback. Any error or panic rolls back.
	runErr := func() (ret error) {
		defer func() {
			if r := recover(); r != nil {
				ret = fmt.Errorf("fits: EditFile callback panic: %v", r)
			}
		}()
		return fn(in, w)
	}()
	if runErr != nil {
		cleanupTmp()
		return runErr
	}
	if err := tmpFile.Sync(); err != nil {
		cleanupTmp()
		return err
	}
	if err := tmpFile.Close(); err != nil {
		os.Remove(tmpPath)
		return err
	}
	if err := os.Rename(tmpPath, path); err != nil {
		os.Remove(tmpPath)
		return err
	}
	// fsync containing directory so the rename is durable.
	if dir := parentDir(path); dir != "" {
		if d, derr := os.Open(dir); derr == nil {
			d.Sync()
			d.Close()
		}
	}
	return nil
}

// Writer is the streaming output side of EditFile. Users emit HDUs to it
// using CopyHDU (pass-through), CopyHDUWithHeader (new header, same data
// bytes), SkipHDU (omit from output), or AppendImage/AppendBinaryTable to
// add new HDUs.
//
// Writer tracks whether the caller has emitted a primary HDU and rejects
// SkipHDU for the primary (every valid FITS file must have a primary).
type Writer struct {
	out        *block.Writer
	osFile     *os.File
	tmpPath    string
	targetPath string

	primaryEmitted bool
	hduIndex       int
}

// CopyHDU writes h to the output verbatim — both header bytes and data
// bytes are passed through with no re-parse.
func (w *Writer) CopyHDU(h HDU) error {
	rec := hduRecordOf(h)
	if rec == nil {
		return fmt.Errorf("fits: CopyHDU: HDU has no backing record")
	}
	return w.copyFromRecord(rec, nil)
}

// CopyHDUWithHeader writes h to the output, replacing its header with hdr
// but passing the data bytes through verbatim. The caller's hdr must not
// include any mandatory structural keyword that would conflict with the
// data shape on disk — the library does not validate this; the output
// file's integrity is the caller's responsibility.
func (w *Writer) CopyHDUWithHeader(h HDU, hdr *header.Header) error {
	rec := hduRecordOf(h)
	if rec == nil {
		return fmt.Errorf("fits: CopyHDUWithHeader: HDU has no backing record")
	}
	return w.copyFromRecord(rec, hdr)
}

// SkipHDU omits h from the output. Rejects skipping the primary (HDU 0).
func (w *Writer) SkipHDU(h HDU) error {
	rec := hduRecordOf(h)
	if rec == nil {
		return fmt.Errorf("fits: SkipHDU: HDU has no backing record")
	}
	if rec.index == 0 {
		return fmt.Errorf("fits: SkipHDU: cannot skip the primary HDU")
	}
	return nil
}

// copyFromRecord writes the header and data of rec to the output, optionally
// replacing the header with replaceHdr.
func (w *Writer) copyFromRecord(rec *hduRecord, replaceHdr *header.Header) error {
	// Header.
	var hdrBytes []byte
	if replaceHdr != nil {
		b, err := header.Encode(replaceHdr)
		if err != nil {
			return err
		}
		hdrBytes = b
	} else {
		hdrBytes = rec.rawHeader
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
			n := int64(chunk)
			if dataLen-done < n {
				n = dataLen - done
			}
			if err := rec.file.br.ReadRange(rec.dataStart+done, buf[:n]); err != nil {
				return err
			}
			if _, err := w.out.Write(buf[:n]); err != nil {
				return err
			}
			done += n
		}
	}

	if w.hduIndex == 0 {
		w.primaryEmitted = true
	}
	w.hduIndex++
	return nil
}

// writeImage on *Writer: append a new image HDU. Implements WriteTarget.
func (w *Writer) writeImage(bp bitpix.BITPIX, hdr *header.Header, shape []int64, data anyData) (*ImageHDU, error) {
	isPrimary := w.hduIndex == 0
	rec, err := emitImageHDU(w.out, bp, hdr, shape, data, isPrimary, w.hduIndex)
	if err != nil {
		return nil, err
	}
	w.hduIndex++
	if isPrimary {
		w.primaryEmitted = true
	}
	return &ImageHDU{rec: rec}, nil
}
