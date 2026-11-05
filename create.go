package fits

import (
	"context"
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/dmarkham/fits/internal/block"
)

// Create opens a new file for writing, truncating any existing contents.
// The returned *File is in ModeCreate and accepts HDUs through WriteImage
// and the Append* functions. When the first HDU written is an extension
// (AppendBinaryTable, AppendASCIITable, AppendCompressedImage,
// AppendCompressedFloat32Image), an empty primary (BITPIX 16, NAXIS 0) is
// written in front of it, as cfitsio does. Each call writes its HDU completely (header,
// data and the zero padding to a 2880-byte block) before returning; Close
// releases the file handle. Nothing is crash-safe here: a crash while
// writing leaves a partial file.
//
// Create does not rewrite a file another edit of this package owns
// (TestCreateRespectsJournalAndLocks): a pending undo journal
// (<name>.journal, left by an interrupted Flush, OverwritePixels or
// append) makes it fail with ErrPendingJournal, so the next recovery
// cannot write stale bytes into the new file; an open OpenForEdit or
// EditFile makes it fail with ErrLocked; and on Linux and macOS a
// read-only Open (or OpenReader of an *os.File) of the file makes it fail
// with ErrOpenForReading. Each leaves the file untouched. On Linux and
// macOS the *File holds the edit lock until Close, so OpenForEdit and
// EditFile of the file fail with ErrLocked while it is being written; a
// read-only Open is not stopped and sees the HDUs written so far. Like
// the edit locks these are advisory: programs that do not use this
// package are not stopped.
func Create(name string) (*File, error) { return CreateContext(context.Background(), name) }

// CreateContext is Create with an explicit context, checked before the
// file is created (a cancelled context creates nothing).
func CreateContext(ctx context.Context, name string) (*File, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	// Lock (and check for a pending journal) before truncating: the file
	// may be one another edit of this package owns.
	fh, err := openWithLock(name, os.O_RDWR|os.O_CREATE, 0o666, lockForCreate)
	if err != nil {
		return nil, err
	}
	if err := checkDescriptorJournal(fh); err != nil {
		fh.Close()
		return nil, err
	}
	if err := truncateForCreate(fh); err != nil {
		fh.Close()
		return nil, err
	}
	f := &File{mode: ModeCreate, name: name, closer: fh, rs: fh, rws: fh}
	bw, err := block.NewWriter(fh)
	if err != nil {
		fh.Close()
		return nil, err
	}
	f.bw = bw
	// br is initialized lazily after Close-and-reopen, or never if the caller
	// only writes.
	return f, nil
}

// ErrStreamNotEmpty is returned by CreateWriter for a stream that already
// holds bytes.
var ErrStreamNotEmpty = errors.New("fits: CreateWriter needs an empty stream")

// CreateWriter opens a FITS writer over an io.ReadWriteSeeker. The stream
// must be empty: HDUs written after existing content would make the result
// invalid FITS, and the stream is not truncated, so a stream whose end is
// not at offset 0 is refused with ErrStreamNotEmpty and nothing is
// written. The caller retains ownership of rws; *File.Close does not close
// it.
func CreateWriter(rws io.ReadWriteSeeker) (*File, error) {
	end, err := rws.Seek(0, io.SeekEnd)
	if err != nil {
		return nil, err
	}
	if end != 0 {
		return nil, fmt.Errorf("%w (it holds %d bytes)", ErrStreamNotEmpty, end)
	}
	f := &File{mode: ModeCreate, rs: rws, rws: rws}
	bw, err := block.NewWriter(rws)
	if err != nil {
		return nil, err
	}
	f.bw = bw
	return f, nil
}

// writeAssertMode returns an error if the file is not in a write-capable mode.
func (f *File) writeAssertMode() error {
	if f.mode != ModeCreate && f.mode != ModeEdit {
		return ErrReadOnly
	}
	if f.bw == nil {
		return fmt.Errorf("fits: no writer available (mode=%v)", f.mode)
	}
	return nil
}
