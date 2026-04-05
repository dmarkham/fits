package fits

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"io"
	"os"
	"sync"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
	"github.com/dmarkham/fits/internal/block"
)

// Mode is the access mode of an open *File.
type Mode int

const (
	// ModeRead opens the file read-only. Set/Add/Delete on returned headers
	// may succeed in memory but (*File).Flush returns ErrReadOnly.
	ModeRead Mode = iota
	// ModeEdit opens the file read-write for in-place mutations and
	// end-of-file HDU appends. Structural changes require EditFile.
	ModeEdit
	// ModeCreate opens a new file for writing. Existing contents are
	// truncated.
	ModeCreate
)

// File is a handle to an open FITS file. Not safe for concurrent use by
// multiple goroutines; see the package doc for the concurrency model.
type File struct {
	mode    Mode
	name    string
	closer  io.Closer         // non-nil if File owns the underlying stream
	rs      io.ReadSeeker     // read path
	rws     io.ReadWriteSeeker // set only in ModeEdit / ModeCreate
	br      *block.Reader
	bw      *block.Writer

	hdus []*hduRecord // ordered by file position
}

// hduRecord carries the minimal eager metadata extracted during the HDU scan
// plus a lazily-parsed full header (plan decision 5.2).
type hduRecord struct {
	index int
	kind  hduKind

	headerStart int64 // byte offset of first header block
	dataStart   int64 // byte offset of first data byte (dataStart == headerEnd)
	dataEnd     int64 // byte offset just past last data byte (before block pad)
	paddedEnd   int64 // byte offset of the start of the next HDU

	// Eager structural metadata.
	simple   bool   // primary HDU only
	xtension string // "IMAGE", "BINTABLE", "TABLE", or "" for primary
	bitpix   bitpix.BITPIX
	naxis    int
	shape    []int64 // NAXIS1..NAXISn
	pcount   int64
	gcount   int64
	extname  string
	extver   int64
	tfields  int64 // only for table HDUs

	// Raw header bytes (excluding the trailing END/padding) preserved for
	// lazy full parsing.
	rawHeader []byte

	// Lazy-parsed full header, cached behind once.
	once   sync.Once
	parsed *header.Header
	parseErr error

	// dirty indicates whether the in-memory header has been mutated through
	// the public API and must be re-serialized on Flush (ModeEdit only).
	dirty bool

	// Parent File reference for lazy-load callbacks.
	file *File
}

// hduKind is the coarse type of an HDU, derived from SIMPLE/XTENSION.
type hduKind int

const (
	kindPrimaryImage hduKind = iota
	kindImageExt
	kindBinTable
	kindASCIITable
	kindRandomGroups // not supported in v1
	kindCompressed   // detected tile-compressed, not supported in v1
	kindUnknown
)

// Open opens a FITS file by path for read-only access.
func Open(name string) (*File, error) { return OpenContext(context.Background(), name) }

// OpenContext is Open with explicit context.Context. The context is consulted
// during the initial HDU scan; if it is cancelled, the open fails early.
func OpenContext(ctx context.Context, name string) (*File, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	f, err := os.Open(name)
	if err != nil {
		return nil, err
	}
	fi, err := f.Stat()
	if err != nil {
		f.Close()
		return nil, err
	}
	if fi.Size() == 0 {
		f.Close()
		return nil, ErrEmptyFile
	}
	file := &File{mode: ModeRead, name: name, closer: f, rs: f}
	if err := file.initReader(ctx); err != nil {
		f.Close()
		return nil, err
	}
	return file, nil
}

// OpenReader opens a FITS file over an existing io.ReadSeeker. The caller
// retains ownership of r — Close() on the returned *File will NOT close r.
func OpenReader(r io.ReadSeeker) (*File, error) {
	file := &File{mode: ModeRead, rs: r}
	if err := file.initReader(context.Background()); err != nil {
		return nil, err
	}
	return file, nil
}

// OpenReadAll slurps the contents of an io.Reader into memory and returns a
// regular *File. This is the helper for stream sources that do not already
// provide random access (see plan decision 5.3).
func OpenReadAll(r io.Reader) (*File, error) {
	buf, err := io.ReadAll(r)
	if err != nil {
		return nil, fmt.Errorf("fits: slurp reader: %w", err)
	}
	if len(buf) == 0 {
		return nil, ErrEmptyFile
	}
	return OpenReader(bytes.NewReader(buf))
}

// initReader runs the HDU scan pass over the file's current content.
func (f *File) initReader(ctx context.Context) error {
	br, err := block.NewReader(f.rs)
	if err != nil {
		return err
	}
	f.br = br
	return f.scanHDUs(ctx)
}

// Close releases any resources owned by the *File. It is safe to call on a
// *File returned from OpenReader; in that case only internal state is
// released.
func (f *File) Close() error {
	if f.closer != nil {
		return f.closer.Close()
	}
	return nil
}

// Mode returns the access mode of the *File.
func (f *File) Mode() Mode { return f.mode }

// Name returns the file name supplied at Open, or "" for OpenReader/OpenReadAll.
func (f *File) Name() string { return f.name }

// NumHDU returns the total number of HDUs in the file.
func (f *File) NumHDU() int { return len(f.hdus) }

// HDU returns the i-th HDU (zero-indexed; HDU(0) is the primary). An out-of-
// range index returns an error.
func (f *File) HDU(i int) (HDU, error) {
	if i < 0 || i >= len(f.hdus) {
		return nil, fmt.Errorf("fits: HDU index %d out of range [0,%d)", i, len(f.hdus))
	}
	return makeHDU(f.hdus[i])
}

// HDUByName returns the first extension HDU whose EXTNAME matches name.
// EXTNAME comparison is case-sensitive per convention. Returns an error if
// no matching HDU is found.
func (f *File) HDUByName(name string) (HDU, error) {
	for _, r := range f.hdus {
		if r.extname == name {
			return makeHDU(r)
		}
	}
	return nil, fmt.Errorf("fits: no HDU with EXTNAME %q", name)
}

// Primary returns the primary HDU as an *ImageHDU. The primary HDU is always
// an image (though it may be a zero-axis placeholder, i.e. NAXIS=0).
func (f *File) Primary() (*ImageHDU, error) {
	if len(f.hdus) == 0 {
		return nil, errors.New("fits: file has no HDUs")
	}
	h, err := f.HDU(0)
	if err != nil {
		return nil, err
	}
	img, ok := h.(*ImageHDU)
	if !ok {
		return nil, fmt.Errorf("fits: primary HDU is not an image")
	}
	return img, nil
}
