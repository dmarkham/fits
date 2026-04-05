package fits

import "github.com/dmarkham/fits/header"

// ASCIITableHDU represents an ASCII TABLE extension (§7.2).
type ASCIITableHDU struct {
	rec *hduRecord

	columnsOnce bool
	columns     []Column
	columnsErr  error
}

// Type returns TypeASCIITable.
func (h *ASCIITableHDU) Type() HDUType { return TypeASCIITable }

// Index returns the 0-based HDU index.
func (h *ASCIITableHDU) Index() int { return h.rec.index }

// Header returns the parsed header (lazy).
func (h *ASCIITableHDU) Header() *header.Header {
	hdr, err := h.rec.loadHeader()
	if err != nil {
		panic(err)
	}
	return hdr
}

// NumRows returns the number of rows from NAXIS2.
func (h *ASCIITableHDU) NumRows() int64 {
	if len(h.rec.shape) >= 2 {
		return h.rec.shape[1]
	}
	return 0
}
