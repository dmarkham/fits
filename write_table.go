package fits

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bigendian"
	"github.com/dmarkham/fits/internal/block"
	"github.com/dmarkham/fits/internal/tform"
)

// ColumnData describes a single column to write into a new binary table.
// Exactly one of Data* fields may be set; the type selects the TFORM code.
type ColumnData struct {
	Name    string
	Unit    string
	Display string
	TForm   string // optional override; if empty, inferred from Data*
	Dim     []int64

	// One of the following (mutually exclusive).
	DataUint8   []uint8
	DataInt16   []int16
	DataInt32   []int32
	DataInt64   []int64
	DataFloat32 []float32
	DataFloat64 []float64
	DataString  []string // for Char columns; all entries must be the same length
	DataBool    []bool
}

// AppendBinaryTable appends a new BINTABLE HDU to the file. All columns must
// have the same row count. The TFORM code is inferred from the DataXxx slice
// that is populated (or from an explicit .TForm if set).
func AppendBinaryTable(f *File, hdr *header.Header, cols []ColumnData) (*BinaryTableHDU, error) {
	if err := f.writeAssertMode(); err != nil {
		return nil, err
	}
	if len(cols) == 0 {
		return nil, fmt.Errorf("fits: AppendBinaryTable: no columns")
	}
	nrows, rowBytes, colInfos, err := describeBinaryColumns(cols)
	if err != nil {
		return nil, err
	}
	// Seek to EOF.
	if f.mode == ModeEdit && len(f.hdus) > 0 {
		last := f.hdus[len(f.hdus)-1]
		if _, err := f.bw.Seek(last.paddedEnd, io.SeekStart); err != nil {
			return nil, err
		}
	} else if f.mode == ModeCreate {
		if _, err := f.bw.Seek(0, io.SeekEnd); err != nil {
			return nil, err
		}
	}
	start := f.bw.Pos()
	// Build mandatory header.
	h := header.New()
	h.Set(header.KeyXtension, "BINTABLE", "binary table extension")
	h.Set(header.KeyBitpix, int64(8), "")
	h.Set(header.KeyNaxis, int64(2), "")
	h.Set("NAXIS1", int64(rowBytes), "width of table row in bytes")
	h.Set("NAXIS2", int64(nrows), "number of rows")
	h.Set(header.KeyPcount, int64(0), "no heap in v1 writer")
	h.Set(header.KeyGcount, int64(1), "")
	h.Set(header.KeyTfields, int64(len(cols)), "number of columns")
	for i, ci := range colInfos {
		k := strconv.Itoa(i + 1)
		h.Set("TFORM"+k, ci.tform, "")
		if cols[i].Name != "" {
			h.Set("TTYPE"+k, cols[i].Name, "")
		}
		if cols[i].Unit != "" {
			h.Set("TUNIT"+k, cols[i].Unit, "")
		}
		if cols[i].Display != "" {
			h.Set("TDISP"+k, cols[i].Display, "")
		}
		if len(cols[i].Dim) > 0 {
			dims := make([]string, len(cols[i].Dim))
			for j, d := range cols[i].Dim {
				dims[j] = strconv.FormatInt(d, 10)
			}
			h.Set("TDIM"+k, "("+strings.Join(dims, ",")+")", "")
		}
	}
	// Append user-supplied header cards (skipping any structural collisions).
	if hdr != nil {
		for _, c := range hdr.Cards() {
			if isMandatoryTableKey(c.Key) {
				continue
			}
			h.Add(c.Key, c.Value, c.Comment)
		}
	}
	hdrBytes, err := header.Encode(h)
	if err != nil {
		return nil, err
	}
	if err := f.bw.WriteRange(hdrBytes); err != nil {
		return nil, err
	}
	dataStart := f.bw.Pos()

	// Emit row bytes.
	if err := emitBinaryRows(f.bw, cols, colInfos, nrows, rowBytes); err != nil {
		return nil, err
	}
	dataEnd := f.bw.Pos()
	if err := f.bw.PadToBlock(0); err != nil {
		return nil, err
	}
	paddedEnd := f.bw.Pos()

	rec := &hduRecord{
		index:       len(f.hdus),
		kind:        kindBinTable,
		headerStart: start,
		dataStart:   dataStart,
		dataEnd:     dataEnd,
		paddedEnd:   paddedEnd,
		xtension:    "BINTABLE",
		bitpix:      8,
		naxis:       2,
		shape:       []int64{rowBytes, nrows},
		gcount:      1,
		tfields:     int64(len(cols)),
		rawHeader:   hdrBytes,
		parsed:      h,
		file:        f,
	}
	rec.once.Do(func() {})
	f.hdus = append(f.hdus, rec)
	return &BinaryTableHDU{rec: rec}, nil
}

type binColInfo struct {
	tform string
	bin   tform.BinaryForm
	size  int64 // bytes per row for this column
	repeat int64
}

// describeBinaryColumns inspects each ColumnData, infers the TFORM, and
// validates that every column has the same number of rows.
func describeBinaryColumns(cols []ColumnData) (nrows, rowBytes int64, out []binColInfo, err error) {
	out = make([]binColInfo, len(cols))
	nrows = -1
	for i, c := range cols {
		n, tf, bf, sz, rpt, e := describeOneColumn(c)
		if e != nil {
			return 0, 0, nil, fmt.Errorf("column %d (%q): %w", i+1, c.Name, e)
		}
		if nrows < 0 {
			nrows = n
		} else if n != nrows {
			return 0, 0, nil, fmt.Errorf("column %d (%q): row count %d != %d", i+1, c.Name, n, nrows)
		}
		out[i] = binColInfo{tform: tf, bin: bf, size: sz, repeat: rpt}
		rowBytes += sz
	}
	if nrows < 0 {
		nrows = 0
	}
	return nrows, rowBytes, out, nil
}

// describeOneColumn looks at a single ColumnData and returns row count,
// TFORM string, parsed binary form, byte size per row, and repeat.
func describeOneColumn(c ColumnData) (int64, string, tform.BinaryForm, int64, int64, error) {
	var n int64
	var tf string
	switch {
	case c.DataUint8 != nil:
		n = int64(len(c.DataUint8))
		tf = "1B"
	case c.DataInt16 != nil:
		n = int64(len(c.DataInt16))
		tf = "1I"
	case c.DataInt32 != nil:
		n = int64(len(c.DataInt32))
		tf = "1J"
	case c.DataInt64 != nil:
		n = int64(len(c.DataInt64))
		tf = "1K"
	case c.DataFloat32 != nil:
		n = int64(len(c.DataFloat32))
		tf = "1E"
	case c.DataFloat64 != nil:
		n = int64(len(c.DataFloat64))
		tf = "1D"
	case c.DataString != nil:
		n = int64(len(c.DataString))
		maxLen := 0
		for _, s := range c.DataString {
			if len(s) > maxLen {
				maxLen = len(s)
			}
		}
		if maxLen == 0 {
			maxLen = 1
		}
		tf = strconv.Itoa(maxLen) + "A"
	case c.DataBool != nil:
		n = int64(len(c.DataBool))
		tf = "1L"
	default:
		return 0, "", tform.BinaryForm{}, 0, 0, fmt.Errorf("no data set")
	}
	if c.TForm != "" {
		tf = c.TForm
	}
	bf, err := tform.ParseBinary(tf)
	if err != nil {
		return 0, "", tform.BinaryForm{}, 0, 0, err
	}
	var colBytes int64
	switch bf.Type {
	case tform.BinChar:
		colBytes = bf.Repeat
	case tform.BinBit:
		colBytes = (bf.Repeat + 7) / 8
	default:
		colBytes = bf.Repeat * int64(bf.Type.ElementSize())
	}
	return n, tf, bf, colBytes, bf.Repeat, nil
}

// emitBinaryRows writes nrows rows of data to bw, with each row packed in
// the declared column order.
func emitBinaryRows(bw *block.Writer, cols []ColumnData, info []binColInfo, nrows, rowBytes int64) error {
	row := make([]byte, rowBytes)
	for r := int64(0); r < nrows; r++ {
		var off int
		for i, ci := range info {
			fillCell(row[off:off+int(ci.size)], cols[i], ci, r)
			off += int(ci.size)
		}
		if _, err := bw.Write(row); err != nil {
			return err
		}
	}
	return nil
}

// fillCell serializes the r-th row of one column into cell (exactly ci.size
// bytes).
func fillCell(cell []byte, col ColumnData, ci binColInfo, r int64) {
	switch ci.bin.Type {
	case tform.BinUint8:
		cell[0] = col.DataUint8[r]
	case tform.BinInt16:
		bigendian.PutInt16(cell, col.DataInt16[r])
	case tform.BinInt32:
		bigendian.PutInt32(cell, col.DataInt32[r])
	case tform.BinInt64:
		bigendian.PutInt64(cell, col.DataInt64[r])
	case tform.BinFloat32:
		bigendian.PutFloat32(cell, col.DataFloat32[r])
	case tform.BinFloat64:
		bigendian.PutFloat64(cell, col.DataFloat64[r])
	case tform.BinLogical:
		if col.DataBool[r] {
			cell[0] = 'T'
		} else {
			cell[0] = 'F'
		}
	case tform.BinChar:
		s := col.DataString[r]
		n := copy(cell, s)
		for i := n; i < len(cell); i++ {
			cell[i] = ' '
		}
	}
}

// isMandatoryTableKey reports whether key is a structural table keyword we
// synthesize from the column data.
func isMandatoryTableKey(key string) bool {
	switch key {
	case header.KeyXtension, header.KeyBitpix, header.KeyNaxis, header.KeyPcount,
		header.KeyGcount, header.KeyTfields, header.KeyEnd:
		return true
	}
	// NAXISn, TFORMn, TTYPEn, TUNITn, TDISPn, TDIMn, TSCALn, TZEROn, TNULLn, TBCOLn.
	prefixes := []string{"NAXIS", "TFORM", "TTYPE", "TUNIT", "TDISP", "TDIM", "TSCAL", "TZERO", "TNULL", "TBCOL"}
	for _, p := range prefixes {
		if strings.HasPrefix(key, p) {
			rest := key[len(p):]
			if rest == "" {
				continue
			}
			onlyDigits := true
			for i := 0; i < len(rest); i++ {
				if rest[i] < '0' || rest[i] > '9' {
					onlyDigits = false
					break
				}
			}
			if onlyDigits {
				return true
			}
		}
	}
	return false
}
