package fits

import (
	"fmt"

	"github.com/dmarkham/fits/header"
)

// HDUType enumerates the broad HDU categories exposed by the library.
type HDUType int

const (
	// TypeImage identifies a primary HDU or an IMAGE extension.
	TypeImage HDUType = iota
	// TypeASCIITable identifies an ASCII TABLE extension (§7.2).
	TypeASCIITable
	// TypeBinaryTable identifies a BINTABLE extension (§7.3).
	TypeBinaryTable
)

func (t HDUType) String() string {
	switch t {
	case TypeImage:
		return "image"
	case TypeASCIITable:
		return "ascii_table"
	case TypeBinaryTable:
		return "binary_table"
	}
	return fmt.Sprintf("HDUType(%d)", int(t))
}

// HDU is the common interface implemented by every HDU kind.
//
// Concrete types are *ImageHDU, *ASCIITableHDU, *BinaryTableHDU, and
// *CompressedImageHDU. Narrow to the concrete type via type assertion.
type HDU interface {
	Type() HDUType
	// Header returns the full parsed header. Mutations to the returned
	// *header.Header persist in memory but are only written to disk when
	// (*File).Flush is called on a ModeEdit *File. See plan decision 5.8.
	Header() *header.Header
	// Index returns the 0-based HDU index in the parent file.
	Index() int
	// Compressed reports whether this HDU is a tile-compressed image
	// backed by a binary table with ZIMAGE=T. Always false for plain
	// ImageHDU / BinaryTableHDU / ASCIITableHDU; true only for
	// CompressedImageHDU.
	Compressed() bool
	// CompressionType returns the ZCMPTYPE value ("RICE_1", "GZIP_1",
	// etc.) for compressed HDUs, or "" for uncompressed HDUs.
	CompressionType() string
}

// makeHDU constructs the concrete HDU wrapper for a given record. It does
// not trigger the lazy full-header parse; that happens on first Header()
// call.
func makeHDU(r *hduRecord) (HDU, error) {
	switch r.kind {
	case kindPrimaryImage, kindImageExt:
		return &ImageHDU{rec: r}, nil
	case kindBinTable:
		return &BinaryTableHDU{rec: r}, nil
	case kindASCIITable:
		return &ASCIITableHDU{rec: r}, nil
	case kindCompressed:
		return &CompressedImageHDU{rec: r, tbl: &BinaryTableHDU{rec: r}}, nil
	case kindRandomGroups:
		return nil, ErrRandomGroups
	}
	return nil, fmt.Errorf("fits: unknown HDU kind %v at index %d", r.kind, r.index)
}

// loadHeader parses the retained raw header bytes into a full
// *header.Header on first call, caches the result, and returns the cache on
// subsequent calls. See plan decision 5.2.
func (r *hduRecord) loadHeader() (*header.Header, error) {
	r.once.Do(func() {
		cards, _, err := header.ParseCards(r.rawHeader)
		if err != nil {
			r.parseErr = err
			return
		}
		r.parsed = header.FromCards(cards)
	})
	return r.parsed, r.parseErr
}
