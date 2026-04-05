package fits

import (
	"context"
	"fmt"
	"io"
	"os"

	"github.com/dmarkham/fits/internal/block"
)

// Create opens a new file for writing, truncating any existing contents.
// The returned *File is in ModeCreate and accepts HDU appends through the
// Append* family of functions. Callers must call Close() to finalize the
// trailing block pad.
func Create(name string) (*File, error) { return CreateContext(context.Background(), name) }

// CreateContext is Create with an explicit context.
func CreateContext(ctx context.Context, name string) (*File, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	fh, err := os.Create(name)
	if err != nil {
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

// CreateWriter opens a FITS writer over an existing io.ReadWriteSeeker. The
// caller retains ownership of rws and must not close it via *File.Close.
func CreateWriter(rws io.ReadWriteSeeker) (*File, error) {
	// Rewind and truncate if possible. We require the caller to pass an empty
	// stream or to have already truncated — FITS Create always starts at
	// offset zero.
	if _, err := rws.Seek(0, io.SeekStart); err != nil {
		return nil, err
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
