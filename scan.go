package fits

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/bitpix"
	"github.com/dmarkham/fits/internal/block"
)

// maxHeaderBlocks caps the number of 2880-byte blocks a single header may
// span before we give up and return ErrUnterminatedHeader. Real FITS files
// almost never exceed a handful of header blocks; 1024 (≈3 MB) is a generous
// safety limit against pathological or malicious inputs.
const maxHeaderBlocks = 1024

// scanHDUs walks the entire file, extracts the minimal eager structural
// metadata for each HDU, and records offsets into f.hdus. It does NOT parse
// full headers — that happens lazily via hduRecord.loadHeader (plan
// decision 5.2).
//
// Framing tolerance follows cfitsio (proven against fits_get_num_hdus /
// fits_get_hduaddrll on the samples in testdata/images/framing):
//
//   - After the last HDU, anything that does not begin a new extension —
//     fewer than 2880 bytes, or a block not starting with "XTENSION=" (a
//     zero block, junk text) — ends the HDU list. Those bytes are not part
//     of any HDU; their count is reported by (*File).TrailingBytes.
//   - The last data unit may lack its final block padding: the data bytes
//     are all present, and the missing pad reads as zeros.
//
// Deliberately stricter than cfitsio: a block that starts with "XTENSION="
// but has no END card before EOF is an error (ErrUnterminatedHeader), as is
// a data unit whose declared bytes run past EOF (ErrTruncatedData); cfitsio
// silently drops the former and fails later on the latter.
func (f *File) scanHDUs(ctx context.Context) error {
	total := f.br.Size()
	if total == 0 {
		return ErrEmptyFile
	}

	off := int64(0)
	idx := 0
	for off < total {
		if err := ctx.Err(); err != nil {
			return err
		}
		if idx > 0 && !f.startsExtension(off, total) {
			f.trailing = total - off
			break
		}
		// The first card must be SIMPLE (FITS 4.0 §4.4.1.1); check it
		// before looking for END, so a non-FITS file is ErrNotFITS rather
		// than an unterminated header.
		if idx == 0 && !f.startsPrimary(total) {
			return ErrNotFITS
		}
		// Scan header blocks until we find an END card or hit the safety cap.
		headerBytes, headerEnd, err := f.readHeaderAt(off, idx)
		if err != nil {
			return err
		}
		// Very first HDU: the first card must be SIMPLE = T.
		if idx == 0 {
			if len(headerBytes) < 10 || !bytes.HasPrefix(headerBytes, []byte("SIMPLE  =")) {
				return ErrNotFITS
			}
		}
		// Extract structural keywords via a restricted parse.
		rec, err := extractStructural(headerBytes, idx)
		if err != nil {
			return err
		}
		rec.index = idx
		rec.file = f
		rec.headerStart = off
		rec.dataStart = headerEnd
		rec.rawHeader = headerBytes

		// Compute data size per Eq. 2:
		//   |BITPIX|/8 * GCOUNT * (PCOUNT + NAXIS1*...*NAXISm)
		dataBytes := dataSize(rec)
		rec.dataEnd = rec.dataStart + dataBytes
		rec.paddedEnd = rec.dataStart + block.RoundUpBlocks(dataBytes)

		if rec.paddedEnd > total {
			if rec.dataEnd > total {
				// Declared data runs past EOF — strict failure per plan 5.11.
				return &ErrTruncatedData{HDU: idx, Expected: dataBytes, Got: total - rec.dataStart}
			}
			// Only the final zero padding is missing.
			f.br.ZeroFillTo(rec.paddedEnd)
		}

		f.hdus = append(f.hdus, rec)
		off = rec.paddedEnd
		idx++
	}
	return nil
}

// startsPrimary reports whether the file begins with a "SIMPLE  =" card.
func (f *File) startsPrimary(total int64) bool {
	n := min(total, block.Size)
	b := make([]byte, n)
	if err := f.br.ReadRange(0, b); err != nil {
		return false
	}
	return bytes.HasPrefix(b, []byte("SIMPLE  ="))
}

// startsExtension reports whether a complete block beginning with an
// XTENSION card starts at off.
func (f *File) startsExtension(off, total int64) bool {
	if total-off < block.Size || off%block.Size != 0 {
		return false
	}
	b, err := f.br.ReadBlock(off / block.Size)
	if err != nil {
		return false
	}
	return bytes.HasPrefix(b[:], []byte("XTENSION="))
}

// readHeaderAt reads header blocks starting at off until an END card is
// located. Returns the concatenated header bytes (always a multiple of 2880)
// and the byte offset just past the final header block.
func (f *File) readHeaderAt(off int64, hduIdx int) ([]byte, int64, error) {
	start := off
	var buf bytes.Buffer
	blocksRead := 0
	for {
		if blocksRead >= maxHeaderBlocks {
			return nil, 0, &ErrUnterminatedHeader{HDU: hduIdx, BytesScanned: int64(blocksRead) * block.Size}
		}
		blkIdx := off / block.Size
		b, err := f.br.ReadBlock(blkIdx)
		if errors.Is(err, block.ErrShortBlock) {
			// EOF before the END card.
			return nil, 0, &ErrUnterminatedHeader{HDU: hduIdx, BytesScanned: off - start}
		}
		if err != nil {
			return nil, 0, fmt.Errorf("fits: read header block at offset %d: %w", off, err)
		}
		buf.Write(b[:])
		blocksRead++
		off += block.Size
		if containsEND(b[:]) {
			return buf.Bytes(), off, nil
		}
	}
}

// containsEND reports whether the 2880-byte block holds an END card.
// The END card is exactly "END" followed by 77 spaces, aligned on an 80-byte
// boundary within the block.
func containsEND(b []byte) bool {
	for i := 0; i < len(b); i += header.CardWidth {
		if i+header.CardWidth > len(b) {
			return false
		}
		card := b[i : i+header.CardWidth]
		if card[0] == 'E' && card[1] == 'N' && card[2] == 'D' && (card[3] == ' ' || card[3] == 0) {
			// Verify columns 4..80 are spaces (§4.4.1).
			rest := card[3:]
			for _, c := range rest {
				if c != ' ' && c != 0 {
					// Not a valid END — continue scanning (defensive).
					goto next
				}
			}
			return true
		}
	next:
	}
	return false
}

// dataSize returns the number of data-area bytes declared by the HDU's
// eager structural metadata. Zero-NAXIS HDUs return zero (no data area).
// Random groups (FITS 4.0 §6, eq. 4) leave NAXIS1 = 0 out of the product.
func dataSize(r *hduRecord) int64 {
	if r.naxis == 0 {
		return 0
	}
	axes := r.shape
	if r.kind == kindRandomGroups {
		axes = axes[1:]
	}
	var prod int64 = 1
	for _, n := range axes {
		prod *= n
	}
	// (NAXIS1*...*NAXISm + PCOUNT) * GCOUNT * |BITPIX|/8
	bpb := int64(r.bitpix.Size())
	groups := r.gcount
	if groups == 0 {
		groups = 1
	}
	return (prod + r.pcount) * groups * bpb
}

// structuralKeys are the only cards the open-time scan decodes. Every other
// card is left to the lazy full parse, so one malformed non-structural card
// does not make the file unopenable (cfitsio opens such files too).
var structuralKeys = map[string]bool{
	header.KeySimple: true, header.KeyXtension: true, header.KeyBitpix: true,
	header.KeyNaxis: true, header.KeyPcount: true, header.KeyGcount: true,
	header.KeyTfields: true, "ZIMAGE": true, "ZTABLE": true, "GROUPS": true,
}

func isStructuralKey(k string) bool {
	if structuralKeys[k] {
		return true
	}
	if len(k) > 5 && k[:5] == header.KeyNaxis {
		_, err := strconv.Atoi(k[5:])
		return err == nil
	}
	return false
}

// extractStructural parses the headerBytes in restricted mode, pulling out
// only the keywords needed to compute data size and answer hot questions.
// It reports *ErrMissingRequiredKeyword for any missing mandatory keyword.
//
// Only structural cards are decoded (see structuralKeys); a malformed
// structural card is an error. NAXISn cards with n > NAXIS are ignored, as
// cfitsio ignores them; each NAXISn for n <= NAXIS must be present once.
// A second SIMPLE, XTENSION, BITPIX, NAXIS, NAXISn (n <= NAXIS), PCOUNT,
// GCOUNT or TFIELDS card is a *DuplicateKeywordError.
// A primary HDU with GROUPS = T and NAXIS1 = 0 is random groups (kind
// kindRandomGroups): sized per FITS 4.0 eq. 4 so later HDUs are found, and
// refused with ErrRandomGroups by (*File).HDU.
func extractStructural(buf []byte, hduIdx int) (*hduRecord, error) {
	if len(buf)%header.CardWidth != 0 {
		return nil, fmt.Errorf("fits: HDU %d header not card-aligned", hduIdx)
	}
	rec := &hduRecord{gcount: 1}
	nCards := len(buf) / header.CardWidth
	sawEND := false
	hadSimple := false
	hadXtension := false
	hadBitpix := false
	hadNaxis := false
	hadPcount, hadGcount, hadTfields := false, false, false
	groups := false
	axes := map[int]int64{}
	// Structural keywords seen, with their first value. A second card with
	// a different value is refused: cfitsio reads BITPIX/NAXIS/NAXISn by
	// position and an image's PCOUNT/GCOUNT last wins, so such a header
	// has no single reading. A repeat of the same value reads the same
	// either way and is accepted (astrometry.net 5200-series index files
	// written by QFITS repeat BITPIX..TFIELDS in every kdtree header).
	seen := map[string]any{}
	dupAxis := 0 // first NAXISn seen twice with different values (checked against NAXIS below)
	for i := range nCards {
		off := i * header.CardWidth
		raw := buf[off : off+header.CardWidth]
		if malformedEND(raw) {
			return nil, fmt.Errorf("fits: HDU %d card %d: %w", hduIdx, i, ErrMalformedEND)
		}
		key := strings.TrimRight(string(raw[:8]), " ")
		if key == header.KeyEnd {
			sawEND = true
			break
		}
		if !isStructuralKey(key) {
			continue
		}
		c, err := header.DecodeCard(raw)
		if err != nil {
			return nil, fmt.Errorf("fits: HDU %d card %d: %w", hduIdx, i, err)
		}
		switch c.Key {
		case header.KeySimple, header.KeyXtension, header.KeyBitpix, header.KeyNaxis,
			header.KeyPcount, header.KeyGcount, header.KeyTfields:
			v := structuralValue(c.Key, c.Value)
			if first, dup := seen[c.Key]; dup {
				if first != v {
					return nil, &DuplicateKeywordError{HDU: hduIdx, Keyword: c.Key}
				}
				continue
			}
			seen[c.Key] = v
		}
		switch c.Key {
		case header.KeySimple:
			if v, ok := c.Value.(bool); ok {
				rec.simple = v
				hadSimple = true
			}
		case header.KeyXtension:
			if v, ok := c.Value.(string); ok {
				rec.xtension = strings.ToUpper(strings.TrimSpace(v))
				hadXtension = true
			}
		case header.KeyBitpix:
			if v, ok := c.Value.(int64); ok {
				rec.bitpix = bitpix.BITPIX(v)
				hadBitpix = true
			}
		case header.KeyNaxis:
			if v, ok := c.Value.(int64); ok {
				if v < 0 || v > 999 {
					return nil, fmt.Errorf("fits: HDU %d: invalid NAXIS %d", hduIdx, v)
				}
				rec.naxis = int(v)
				hadNaxis = true
			}
		case header.KeyPcount:
			if v, ok := c.Value.(int64); ok {
				rec.pcount = v
				hadPcount = true
			}
		case header.KeyGcount:
			if v, ok := c.Value.(int64); ok {
				rec.gcount = v
				hadGcount = true
			}
		case header.KeyTfields:
			if v, ok := c.Value.(int64); ok {
				rec.tfields = v
				hadTfields = true
			}
		case "ZIMAGE":
			if v, ok := c.Value.(bool); ok && v {
				rec.zimage = true
			}
		case "ZTABLE":
			if v, ok := c.Value.(bool); ok && v {
				rec.ztable = true
			}
		case "GROUPS":
			if v, ok := c.Value.(bool); ok && v {
				groups = true
			}
		default:
			// NAXISn
			n, _ := strconv.Atoi(c.Key[len(header.KeyNaxis):])
			if v, ok := c.Value.(int64); ok && n >= 1 {
				if v < 0 {
					return nil, fmt.Errorf("fits: HDU %d: negative %s = %d", hduIdx, c.Key, v)
				}
				if first, dup := axes[n]; dup {
					if first != v && (dupAxis == 0 || n < dupAxis) {
						dupAxis = n
					}
				} else {
					axes[n] = v
				}
			}
		}
	}
	if !sawEND {
		return nil, &ErrUnterminatedHeader{HDU: hduIdx, BytesScanned: int64(len(buf))}
	}
	if hduIdx == 0 {
		if !hadSimple {
			return nil, &ErrMissingRequiredKeyword{HDU: 0, Keyword: header.KeySimple}
		}
	} else {
		if !hadXtension {
			return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: header.KeyXtension}
		}
	}
	if !hadBitpix {
		return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: header.KeyBitpix}
	}
	if !hadNaxis {
		return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: header.KeyNaxis}
	}
	// A table extension must carry PCOUNT, GCOUNT and TFIELDS (FITS 4.0
	// §7.2.1, §7.3.1); cfitsio will not move to it otherwise (NO_PCOUNT
	// 228, NO_GCOUNT 229, NO_TFIELDS 230). An IMAGE extension may omit
	// PCOUNT/GCOUNT: cfitsio defaults them to 0/1 and so do we.
	switch rec.xtension {
	case "BINTABLE", "A3DTABLE", "3DTABLE", "TABLE":
		switch {
		case !hadPcount:
			return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: header.KeyPcount}
		case !hadGcount:
			return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: header.KeyGcount}
		case !hadTfields:
			return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: header.KeyTfields}
		}
	}
	if !rec.bitpix.Valid() {
		return nil, fmt.Errorf("fits: HDU %d: invalid BITPIX %d", hduIdx, int(rec.bitpix))
	}
	// For every NAXISi 1..naxis there must be exactly one card; NAXISn
	// beyond NAXIS are not axes and are ignored.
	if dupAxis > 0 && dupAxis <= rec.naxis {
		return nil, &DuplicateKeywordError{HDU: hduIdx, Keyword: fmt.Sprintf("NAXIS%d", dupAxis)}
	}
	rec.shape = make([]int64, rec.naxis)
	for i := 1; i <= rec.naxis; i++ {
		v, ok := axes[i]
		if !ok {
			return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: fmt.Sprintf("NAXIS%d", i)}
		}
		rec.shape[i-1] = v
	}
	// Determine kind. Tile-compressed binary tables get reclassified
	// as kindCompressed by classifyKind so makeHDU returns a
	// CompressedImageHDU wrapper.
	rec.kind = classifyKind(rec)
	if hduIdx == 0 && groups && rec.naxis >= 1 && rec.shape[0] == 0 {
		rec.kind = kindRandomGroups
	}
	return rec, nil
}

// classifyKind maps the SIMPLE/XTENSION fields onto an hduKind value.
// Binary tables with ZIMAGE=T are reclassified as kindCompressed so the
// makeHDU factory returns a CompressedImageHDU wrapper.
func classifyKind(r *hduRecord) hduKind {
	if r.xtension == "" {
		return kindPrimaryImage
	}
	switch r.xtension {
	case "IMAGE", "IUEIMAGE":
		// IUEIMAGE: the legacy IUE name cfitsio reads as IMAGE_HDU.
		return kindImageExt
	case "BINTABLE", "A3DTABLE", "3DTABLE":
		// A3DTABLE / 3DTABLE: the legacy names cfitsio reads as BINARY_TBL.
		if r.zimage {
			return kindCompressed
		}
		if r.ztable {
			return kindTiledTable
		}
		return kindBinTable
	case "TABLE":
		return kindASCIITable
	}
	return kindUnknown
}

// structuralValue is the value extractStructural compares when a
// structural keyword repeats: the decoded value with its type (so 8 and
// 8.0 differ), XTENSION trimmed and upper-cased as it is read.
func structuralValue(key string, v any) string {
	if s, ok := v.(string); ok && key == header.KeyXtension {
		v = strings.ToUpper(strings.TrimSpace(s))
	}
	return fmt.Sprintf("%T %v", v, v)
}
