package fits

import (
	"math/big"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/tform"
)

// ColumnType enumerates the supported table column types.
type ColumnType int

const (
	ColByte ColumnType = iota
	ColInt16
	ColInt32
	ColInt64
	ColFloat32
	ColFloat64
	ColString
	ColLogical
	ColBit
	ColComplex64
	ColComplex128
	ColVarArray
)

// Column describes a single table column.
//
// For binary tables Repeat is the TFORM repeat count (number of elements,
// or of bits for X, or of characters for A) and Dim the parsed TDIMn. For
// ASCII tables Repeat is 1 and Width the TFORM field width.
type Column struct {
	Index   int    // 1-based per FITS convention
	Name    string // TTYPEn
	Unit    string // TUNITn
	TForm   string // raw TFORMn
	Repeat  int64
	Type    ColumnType
	Scale   float64 // TSCALn, default 1
	Zero    float64 // TZEROn, default 0
	Null    any     // TNULLn if set: int64 for binary tables, string for ASCII tables
	Display string  // TDISPn
	Dim     []int64 // TDIMn if set
	// Width is the ASCII-table field width (TFORM Aw/Iw/Fw.d/...), or the
	// substring width w of a binary "rAw" TFORM; 0 otherwise. (A TDIMn on
	// an rAw column sets the string width instead: ReadStringArrayColumn.)
	Width int

	// Derived on column-metadata pass.
	byteOffset int // byte offset within a table row
	byteSize   int // total byte size of this column's row cells
	bin        tform.BinaryForm
	ascii      tform.ASCIIForm

	// Integer null (binary tables, B/I/J/K and their P/Q forms only).
	hasNull bool
	nullInt int64
	// ASCII null string (TNULLn) and whether it is defined.
	hasNullStr bool
	nullStr    string
	// scaleLit/zeroLit: the TSCALn/TZEROn literals, exactly (nil when
	// the card is absent).
	scaleLit, zeroLit *big.Rat
	// zeroNear2p63: TZEROn is a float literal other than 2^63 whose
	// nearest double is 2^63 (see readScaleZero).
	zeroNear2p63 bool
}

// BinaryTableHDU represents a BINTABLE extension (§7.3).
type BinaryTableHDU struct {
	rec *hduRecord

	// Column metadata, lazily populated on first Columns()/ReadColumn call.
	columnsOnce bool
	columns     []Column
	columnsErr  error
	rowBytes    int
}

// Type returns TypeBinaryTable.
func (h *BinaryTableHDU) Type() HDUType { return TypeBinaryTable }

// Index returns the 0-based HDU index.
func (h *BinaryTableHDU) Index() int { return h.rec.index }

// Compressed returns false — plain binary tables are not compressed.
// (CompressedImageHDU is a separate type backed by a binary table with
// ZIMAGE=T, not a BinaryTableHDU.)
func (*BinaryTableHDU) Compressed() bool { return false }

// CompressionType returns the empty string for plain binary tables.
func (*BinaryTableHDU) CompressionType() string { return "" }

// Header returns the parsed header (lazy).
func (h *BinaryTableHDU) Header() *header.Header {
	hdr, err := h.rec.loadHeader()
	if err != nil {
		panic(err)
	}
	return hdr
}

// NumRows returns the number of rows from NAXIS2.
func (h *BinaryTableHDU) NumRows() int64 {
	if len(h.rec.shape) >= 2 {
		return h.rec.shape[1]
	}
	return 0
}
