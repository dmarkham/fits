// Package block provides 2880-byte block I/O over io.ReadSeeker and
// io.WriteSeeker.
//
// All FITS files are a sequence of 2880-byte blocks (FITS Standard v4.0
// §3.3.1, §7.3.3). Headers and data sections are always padded up to a
// block boundary. This package is the single canonical place that knows
// about that constant.
//
// This package is INTERNAL. It is not part of the public API and its shape
// may change between minor releases of the parent module.
package block

import (
	"errors"
	"fmt"
	"io"
)

// Size is the FITS logical block size in bytes.
const Size = 2880

// ErrShortBlock indicates an attempt to read or write a partial final block.
var ErrShortBlock = errors.New("fits/internal/block: short block")

// Reader reads fixed-size 2880-byte blocks from an io.ReadSeeker.
//
// Reader keeps a tiny single-block cache so that a hot sequence of small
// reads inside one block avoids redundant seeks and syscalls. The OS page
// cache is the primary cache; this is a lightweight tier on top.
//
// Reader is not safe for concurrent use by multiple goroutines.
type Reader struct {
	r          io.ReadSeeker
	size       int64 // total stream length, -1 if unknown
	cacheIdx   int64 // block index currently held in cache, -1 if empty
	cache      [Size]byte
}

// NewReader returns a new block Reader over r. If r also implements io.Seeker
// (which io.ReadSeeker guarantees) the stream length is determined up front.
func NewReader(r io.ReadSeeker) (*Reader, error) {
	// Determine total length by seeking to end then back to 0.
	end, err := r.Seek(0, io.SeekEnd)
	if err != nil {
		return nil, fmt.Errorf("block: seek end: %w", err)
	}
	if _, err := r.Seek(0, io.SeekStart); err != nil {
		return nil, fmt.Errorf("block: seek start: %w", err)
	}
	return &Reader{r: r, size: end, cacheIdx: -1}, nil
}

// Size returns the total byte length of the underlying stream as observed at
// NewReader time.
func (r *Reader) Size() int64 { return r.size }

// BlockCount returns the number of complete 2880-byte blocks in the stream.
// If the stream length is not an exact multiple of 2880 the remainder is
// reported separately by the caller via Size() — the parser treats a partial
// trailing block as garbage (plan decision 5.11).
func (r *Reader) BlockCount() int64 { return r.size / Size }

// ReadBlock reads block index n (zero-based). It returns ErrShortBlock if the
// underlying stream does not have a full block at that position.
func (r *Reader) ReadBlock(n int64) (*[Size]byte, error) {
	if n < 0 {
		return nil, fmt.Errorf("block: negative block index %d", n)
	}
	if r.cacheIdx == n {
		return &r.cache, nil
	}
	off := n * Size
	if _, err := r.r.Seek(off, io.SeekStart); err != nil {
		return nil, fmt.Errorf("block: seek to block %d: %w", n, err)
	}
	if _, err := io.ReadFull(r.r, r.cache[:]); err != nil {
		r.cacheIdx = -1
		if errors.Is(err, io.ErrUnexpectedEOF) || errors.Is(err, io.EOF) {
			return nil, ErrShortBlock
		}
		return nil, fmt.Errorf("block: read block %d: %w", n, err)
	}
	r.cacheIdx = n
	return &r.cache, nil
}

// ReadRange reads len(dst) bytes starting at absolute byte offset off.
// The read may span any number of blocks. It bypasses the single-block cache
// and is intended for bulk reads of header-blocks or data sections.
func (r *Reader) ReadRange(off int64, dst []byte) error {
	if off < 0 {
		return fmt.Errorf("block: negative offset %d", off)
	}
	if _, err := r.r.Seek(off, io.SeekStart); err != nil {
		return fmt.Errorf("block: seek to %d: %w", off, err)
	}
	if _, err := io.ReadFull(r.r, dst); err != nil {
		if errors.Is(err, io.ErrUnexpectedEOF) || errors.Is(err, io.EOF) {
			return ErrShortBlock
		}
		return fmt.Errorf("block: read at %d: %w", off, err)
	}
	// Any read via ReadRange potentially overlaps the cached block and would
	// otherwise return stale data. Invalidate to be safe.
	r.cacheIdx = -1
	return nil
}

// Invalidate drops any cached state. Call after an external write to the
// underlying stream.
func (r *Reader) Invalidate() { r.cacheIdx = -1 }

// Writer writes fixed-size 2880-byte blocks to an io.WriteSeeker.
//
// Writer does not buffer writes. It is a thin wrapper that enforces block
// alignment and offers WriteBlock / WriteRange / PadToBlock helpers.
//
// Writer is not safe for concurrent use by multiple goroutines.
type Writer struct {
	w   io.WriteSeeker
	pos int64 // current seek position tracked locally
}

// NewWriter returns a new block Writer over w.
func NewWriter(w io.WriteSeeker) (*Writer, error) {
	pos, err := w.Seek(0, io.SeekCurrent)
	if err != nil {
		return nil, fmt.Errorf("block: seek current: %w", err)
	}
	return &Writer{w: w, pos: pos}, nil
}

// Pos returns the locally-tracked write offset.
func (w *Writer) Pos() int64 { return w.pos }

// Seek sets the write position in the underlying stream.
func (w *Writer) Seek(offset int64, whence int) (int64, error) {
	n, err := w.w.Seek(offset, whence)
	if err != nil {
		return n, err
	}
	w.pos = n
	return n, nil
}

// WriteBlock writes one full 2880-byte block at the current position.
func (w *Writer) WriteBlock(b *[Size]byte) error {
	n, err := w.w.Write(b[:])
	w.pos += int64(n)
	if err != nil {
		return fmt.Errorf("block: write block: %w", err)
	}
	if n != Size {
		return fmt.Errorf("block: short write %d/%d", n, Size)
	}
	return nil
}

// WriteRange writes the bytes in src at the current position. len(src) need
// not be block-aligned; the caller is responsible for subsequent padding.
func (w *Writer) WriteRange(src []byte) error {
	_, err := w.Write(src)
	return err
}

// Write implements io.Writer. It keeps Pos() in sync with the underlying
// stream so downstream consumers can compose block.Writer with any code
// that expects an io.Writer (e.g. generic pixel emitters).
func (w *Writer) Write(p []byte) (int, error) {
	n, err := w.w.Write(p)
	w.pos += int64(n)
	if err != nil {
		return n, fmt.Errorf("block: write: %w", err)
	}
	if n != len(p) {
		return n, fmt.Errorf("block: short write %d/%d", n, len(p))
	}
	return n, nil
}

// PadToBlock writes zero bytes (or, if space != 0, the byte `space`) to bring
// the current write position up to the next block boundary. Headers are
// padded with ASCII space (0x20) per FITS §4.1.1; data sections are padded
// with zero bytes per §5.
func (w *Writer) PadToBlock(pad byte) error {
	rem := w.pos % Size
	if rem == 0 {
		return nil
	}
	need := Size - int(rem)
	buf := make([]byte, need)
	if pad != 0 {
		for i := range buf {
			buf[i] = pad
		}
	}
	return w.WriteRange(buf)
}

// RoundUpBlocks returns n rounded up to the next multiple of Size.
func RoundUpBlocks(n int64) int64 {
	if n <= 0 {
		return 0
	}
	return ((n + Size - 1) / Size) * Size
}
