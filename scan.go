package fits

import (
	"bytes"
	"context"
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
			// Data section runs past EOF — strict failure per plan 5.11.
			return &ErrTruncatedData{HDU: idx, Expected: dataBytes, Got: total - rec.dataStart}
		}

		f.hdus = append(f.hdus, rec)
		off = rec.paddedEnd
		idx++
	}
	return nil
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
		if err != nil {
			return nil, 0, fmt.Errorf("fits: read header block at offset %d: %w", off, err)
		}
		buf.Write(b[:])
		blocksRead++
		off += block.Size
		if containsEND(b[:]) {
			_ = start
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
func dataSize(r *hduRecord) int64 {
	if r.naxis == 0 {
		return 0
	}
	var prod int64 = 1
	for _, n := range r.shape {
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

// extractStructural parses the headerBytes in restricted mode, pulling out
// only the keywords needed to compute data size and answer hot questions.
// It reports *ErrMissingRequiredKeyword for any missing mandatory keyword.
//
// The restricted parse still uses the card decoder but only inspects the
// handful of keywords it cares about, skipping CONTINUE joining and full
// keyword-index building. Full parsing happens lazily on first Header() call.
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
	axesSeen := map[int]bool{}
	for i := range nCards {
		off := i * header.CardWidth
		c, err := header.DecodeCard(buf[off : off+header.CardWidth])
		if err != nil {
			return nil, fmt.Errorf("fits: HDU %d card %d: %w", hduIdx, i, err)
		}
		if c.IsEnd() {
			sawEND = true
			break
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
				rec.naxis = int(v)
				hadNaxis = true
				rec.shape = make([]int64, rec.naxis)
			}
		case header.KeyPcount:
			if v, ok := c.Value.(int64); ok {
				rec.pcount = v
			}
		case header.KeyGcount:
			if v, ok := c.Value.(int64); ok {
				rec.gcount = v
			}
		case header.KeyTfields:
			if v, ok := c.Value.(int64); ok {
				rec.tfields = v
			}
		case "ZIMAGE":
			if v, ok := c.Value.(bool); ok && v {
				rec.zimage = true
			}
		case header.KeyExtname:
			if v, ok := c.Value.(string); ok {
				rec.extname = v
			}
		case header.KeyExtver:
			if v, ok := c.Value.(int64); ok {
				rec.extver = v
			}
		default:
			// NAXISn: only keys starting with "NAXIS" and having a numeric suffix.
			if strings.HasPrefix(c.Key, header.KeyNaxis) && len(c.Key) > len(header.KeyNaxis) {
				rest := c.Key[len(header.KeyNaxis):]
				n, convErr := strconv.Atoi(rest)
				if convErr == nil && n >= 1 {
					if v, ok := c.Value.(int64); ok {
						if n > len(rec.shape) {
							// NAXISn appeared before NAXIS or beyond NAXIS count;
							// grow shape if needed. FITS allows any ordering but
							// NAXIS must cover the maximum n.
							if rec.naxis < n {
								// Defer strict validation — we will re-check below.
								newShape := make([]int64, n)
								copy(newShape, rec.shape)
								rec.shape = newShape
							}
						}
						rec.shape[n-1] = v
						axesSeen[n] = true
					}
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
	if !rec.bitpix.Valid() {
		return nil, fmt.Errorf("fits: HDU %d: invalid BITPIX %d", hduIdx, int(rec.bitpix))
	}
	// For every NAXISi 1..naxis there must be a card.
	for i := 1; i <= rec.naxis; i++ {
		if !axesSeen[i] {
			return nil, &ErrMissingRequiredKeyword{HDU: hduIdx, Keyword: fmt.Sprintf("NAXIS%d", i)}
		}
	}
	// Determine kind.
	rec.kind = classifyKind(rec)
	if rec.kind == kindCompressed {
		// Still record the HDU — caller will see ErrCompressed on access.
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
	case "IMAGE":
		return kindImageExt
	case "BINTABLE", "A3DTABLE":
		if r.zimage {
			return kindCompressed
		}
		return kindBinTable
	case "TABLE":
		return kindASCIITable
	}
	return kindUnknown
}
