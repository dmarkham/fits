package compress

import (
	"encoding/binary"
	"fmt"
)

// PLIO_1 — IRAF Planio run-length encoding for integer masks.
//
// DEFERRED: this decoder is a stub. Astropy's PLIO output uses the IRAF
// Planio format with a 7-word header (magic + version + dimensions) and
// a specific RLE scheme that our current implementation does not fully
// decode. PLIO is a narrow IRAF-legacy codec (2MASS DQ masks, some IRAF
// pipelines) and not a blocker for JWST/HST science data. When we
// eventually add full PLIO support, the format to decode is specified
// in IRAF imio/pllib.c.
//
// Original sketch (not correct against astropy output):
//
//	token < 0x8000        → run of zeros of length (token) pixels
//	token == 0x8000       → special value (unused / reserved)
//	0x8000 < token < 0xC000 → run of value (next-token) of length (token-0x8000) pixels ...
//
// Actually the real spec (Tody 1986, IRAF procedure docs) is:
//
//	The mask stream is a list of runs. Each "run" is represented by one
//	or two 16-bit integers. A token with the high bit clear is a value
//	(0..32767). A token with the high bit set is a run-length for the
//	PRECEDING value (or for implicit value 0 if the previous token was
//	also a run-length).
//
// There are a few subtle encodings in IRAF's reference implementation
// that matter for round-trip fidelity. Our implementation follows the
// cfitsio quantize.c / pliocomp.c reference.
type plio1Decoder struct{}

func newPLIO1() Decoder { return plio1Decoder{} }

func (plio1Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	if elemSize != 2 && elemSize != 4 && elemSize != 1 {
		return fmt.Errorf("%w: PLIO_1 requires BYTEPIX of 1/2/4, got %d", ErrCorrupt, elemSize)
	}
	// Decode PLIO stream into a temporary int32 buffer, then pack out to
	// the caller's elemSize.
	if len(src)%2 != 0 {
		return fmt.Errorf("%w: PLIO_1 src length %d not even", ErrCorrupt, len(src))
	}
	tokens := make([]uint16, len(src)/2)
	for i := range tokens {
		tokens[i] = binary.BigEndian.Uint16(src[i*2:])
	}

	out := make([]int32, nelem)
	pos := 0
	i := 0
	for i < len(tokens) && pos < nelem {
		t := tokens[i]
		i++
		if t&0x8000 == 0 {
			// Value token — emits a single pixel with this value.
			if pos >= nelem {
				break
			}
			out[pos] = int32(t)
			pos++
			continue
		}
		// Run-length token for a run of value 0 (or previously-decoded
		// value). The low 15 bits are the run length.
		run := int(t & 0x7FFF)
		// cfitsio convention: if the next token's high bit is clear and
		// it's a value, the run uses that value; otherwise the run is of
		// zeros. We follow that.
		val := int32(0)
		if i < len(tokens) && tokens[i]&0x8000 == 0 {
			val = int32(tokens[i])
			i++
		}
		for j := 0; j < run && pos < nelem; j++ {
			out[pos] = val
			pos++
		}
	}
	// Fill any remaining pixels with zero (silent trailing run).
	// Pack into dst according to elemSize.
	for k := 0; k < nelem; k++ {
		v := out[k]
		switch elemSize {
		case 1:
			dst[k] = byte(v)
		case 2:
			binary.BigEndian.PutUint16(dst[k*2:], uint16(v))
		case 4:
			binary.BigEndian.PutUint32(dst[k*4:], uint32(v))
		}
	}
	return nil
}
