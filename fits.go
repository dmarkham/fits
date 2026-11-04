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

// File is a handle to an open FITS file. It is not safe for concurrent use
// by multiple goroutines: reads share one seek position and block cache.
type File struct {
	mode   Mode
	name   string
	closer io.Closer          // non-nil if File owns the underlying stream
	rs     io.ReadSeeker      // read path
	rws    io.ReadWriteSeeker // set only in ModeEdit / ModeCreate
	br     *block.Reader
	bw     *block.Writer

	hdus []*hduRecord // ordered by file position

	// trailing is the number of bytes after the last HDU that do not form
	// an HDU (see scanHDUs).
	trailing int64

	// readLocked is the caller's *os.File that OpenReader took the read
	// lock on (lockForRead); Close releases it without closing the file.
	readLocked *os.File
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
	tfields  int64 // only for table HDUs

	// zimage is true if this HDU is a binary table with ZIMAGE=T,
	// indicating it holds a tile-compressed image per Pence et al. 2010.
	// Detected during the minimal-eager scan; used by classifyKind.
	zimage bool
	ztable bool // ZTABLE = T: tile-compressed table

	// Raw header bytes (excluding the trailing END/padding) preserved for
	// lazy full parsing.
	rawHeader []byte

	// Lazy-parsed full header, cached behind once.
	once     sync.Once
	parsed   *header.Header
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
	kindRandomGroups // legacy radio-interferometry format — not supported
	kindCompressed   // tile-compressed image, routed to CompressedImageHDU
	kindTiledTable   // tile-compressed table (ZTABLE = T) — not supported
	kindUnknown
)

// Open opens a FITS file by path for read-only access.
//
// Open does not repair a file left by an interrupted Flush. If an undo
// journal (<name>.journal) is pending, Open fails with ErrPendingJournal;
// open the file once with OpenForEdit or EditFile to roll it back.
//
// On Linux and macOS Open holds a shared lock on the file until Close, so an
// OpenForEdit editor of it (in this process or another) cannot change its
// bytes in place meanwhile (ErrOpenForReading, see OpenForEdit); if such a
// write is running when Open is called, Open waits for it to finish, then
// checks for a journal (TestInPlaceWriteRefusedWhileReaderOpen). The wait
// does not observe a context. Programs that do not use this package are
// not stopped. OpenReader does the same for an *os.File; OpenReadAll and
// other readers see only the bytes they are given and cannot be protected.
func Open(name string) (*File, error) { return OpenContext(context.Background(), name) }

// OpenContext is Open with explicit context.Context. The context is checked
// before the file is opened and once per HDU during the initial HDU scan;
// if it is cancelled, the open fails with its error (TestContextFunctions).
func OpenContext(ctx context.Context, name string) (*File, error) {
	if err := ctx.Err(); err != nil {
		return nil, err
	}
	f, err := os.Open(name)
	if err != nil {
		return nil, err
	}
	// Take the read lock first, then look for a journal: a write of this
	// package holds the write lock from before it writes its journal until
	// after it removes it, so once the lock is held no in-place write is
	// running or can start until Close.
	if err := lockForRead(f); err != nil {
		f.Close()
		return nil, err
	}
	if err := checkDescriptorJournal(f); err != nil {
		f.Close()
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
//
// When r is an *os.File, OpenReader takes the same read lock as Open
// (released by Close, which leaves r open), and when its name still
// refers to the same file a pending undo journal next to it makes
// OpenReader fail with ErrPendingJournal, as Open does
// (TestOpenReaderRefusesPendingJournal). Any other reader has no file to
// lock and no path to check: it is not protected from an interrupted or
// concurrent edit.
func OpenReader(r io.ReadSeeker) (*File, error) {
	file := &File{mode: ModeRead, rs: r}
	if fh, ok := r.(*os.File); ok {
		if err := lockForRead(fh); err != nil {
			return nil, err
		}
		if err := checkDescriptorJournal(fh); err != nil {
			unlockRead(fh)
			return nil, err
		}
		file.readLocked = fh
	}
	if err := file.initReader(context.Background()); err != nil {
		file.Close()
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
	if f.readLocked != nil {
		unlockRead(f.readLocked)
		f.readLocked = nil
	}
	if f.closer != nil {
		return f.closer.Close()
	}
	return nil
}

// Mode returns the access mode of the *File.
func (f *File) Mode() Mode { return f.mode }

// Name returns the file name supplied at Open, or "" for OpenReader/OpenReadAll.
func (f *File) Name() string { return f.name }

// TrailingBytes returns the number of bytes after the end of the last HDU
// that are not part of any HDU: a trailing zero block, junk text, or a
// partial block. Like cfitsio, Open ignores them rather than failing; this
// reports that they were there. It is 0 for a well-formed file (and for a
// file whose last data unit merely lacks its final padding).
func (f *File) TrailingBytes() int64 { return f.trailing }

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

// HDUByName returns the first HDU, in file order and including the
// primary, whose EXTNAME value equals name exactly: the comparison is
// case-sensitive (trailing blanks are not significant, as in every FITS
// string value). Returns an error if no HDU matches. The value compared is
// the one the HDU's Header() holds now, so in ModeEdit an EXTNAME set or
// deleted but not yet flushed counts, as it does for cfitsio
// (TestHDUByNameSeesUnflushedRename).
//
// A header with more than one EXTNAME card is named by the first one, as
// cfitsio's fits_movnam_hdu reads it (TestHDUByNameDuplicateEXTNAMEMatchesCfitsio,
// truth testdata/cfmovnam). An HDU whose first EXTNAME card has no string
// value has no name here.
//
// cfitsio's fits_movnam_hdu also scans from the primary HDU, but it
// compares case-insensitively, falls back to the HDUNAME keyword when an
// HDU has no EXTNAME, and compares a non-string EXTNAME by its value text.
// HDUByName does none of these, so where they decide cfitsio's answer
// HDUByName silently returns a different HDU (a later one whose EXTNAME
// matches exactly) or a "no HDU" error: on a file whose HDU 1 is 'sci' and
// HDU 2 'SCI', "SCI" is HDU 2 here and HDU 1 for cfitsio; with HDU 3
// carrying only HDUNAME = 'ERR' and HDU 4 EXTNAME = 'ERR', "ERR" is HDU 4
// here and HDU 3 for cfitsio; "5" names no HDU here and the HDU with
// EXTNAME = 5 for cfitsio (TestHDUByNameDifferencesFromCfitsio, truth
// testdata/cfmovnam/golden/names.json; TestHDUByNameFirstMatchCaseSensitive).
func (f *File) HDUByName(name string) (HDU, error) {
	for _, r := range f.hdus {
		if r.currentExtname() == name {
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
