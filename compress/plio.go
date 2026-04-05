package compress

import (
	"encoding/binary"
	"fmt"
)

// PLIO_1 — IRAF Planio line-list encoding for integer masks.
//
// PLIO is a run-length scheme optimized for sparse integer mask data
// (IRAF DQ arrays, NIRCam quality flags, cosmic-ray masks). The
// compressed stream is a sequence of 16-bit tokens interpreted as
// opcodes over a running "pixel value" (PV) and a running "x position".
//
// Reference: cfitsio/pliocomp.c (pl_l2pi), itself an f2c translation of
// the IRAF imio/pllib.c SPP source (Tody 1986).
//
// # Stream format (int16 tokens, 1-based in pliocomp.c)
//
//	Word 1    : unused
//	Word 2    : header length in words (7 for the long form)
//	Word 3    : magic number -100 (long form) or length (short form > 0)
//	Word 4    : low 15 bits of list length (long form)
//	Word 5    : high bits of list length (long form)
//	Word 6..7 : reserved
//	Word 8+   : data tokens
//
// Total list length lllen = (word[5] << 15) + word[4] words
// (header + data, inclusive). Data starts at word llfirt = word[2] + 1
// for long form, or word 4 for short form.
//
// # Token format
//
// Each 16-bit data token splits into:
//
//	opcode = (word >> 12) & 0xf
//	data   = word & 0xfff
//
// Opcodes (direct dispatch from pliocomp.c; case numbers are
// `opcode + 1` because pliocomp.c uses a 1-based switch after
// `++sw0001`):
//
//	0  (case 1) → run of `data` zeros                     [L160 w/ opcode!=4,!=5]
//	1  (case 2) → set PV: low = data, next word = high    [L220]
//	2  (case 3) → PV += data (no pixel emitted)           [L230]
//	3  (case 4) → PV -= data (no pixel emitted)           [L240]
//	4  (case 5) → run of `data` pixels with value PV      [L160 w/ opcode==4]
//	5  (case 6) → run of `data` zeros with PV at the end  [L160 w/ opcode==5]
//	6  (case 7) → PV += data, emit single pixel           [L250]
//	7  (case 8) → PV -= data, emit single pixel           [L260]
//
// PV starts at 1, x starts at 1 (1-based). Each opcode advances x by
// `data` (run opcodes) or 1 (single-pixel opcodes 6/7) or 0 (PV-update
// opcodes 1/2/3). Output is 1D; the FITS-side layer reshapes to the
// tile's logical shape (row-major in the natural iteration order).

// plio1Decoder decodes PLIO_1 tiles. The Decode method reads the input
// bytes as int16 tokens in FITS big-endian and delegates to decodePLIO1.
type plio1Decoder struct{}

func newPLIO1() Decoder { return plio1Decoder{} }

func (plio1Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	if len(src)%2 != 0 {
		return fmt.Errorf("%w: PLIO_1 src length %d not even", ErrCorrupt, len(src))
	}
	tokens := make([]int16, len(src)/2)
	for i := range tokens {
		tokens[i] = int16(binary.BigEndian.Uint16(src[i*2:]))
	}
	return decodePLIO1(tokens, dst, nelem, elemSize)
}

// decodePLIO1 is the token-level decoder. Exported at package scope so
// unit tests can feed tokens directly without round-tripping through
// the byte serialization.
func decodePLIO1(tokens []int16, dst []byte, nelem, elemSize int) error {
	if len(tokens) < 4 {
		return fmt.Errorf("%w: PLIO_1 stream shorter than minimal header (%d words)", ErrCorrupt, len(tokens))
	}
	// Header: 1-based indexing in pliocomp.c → 0-based here.
	//
	//	tokens[0] = ll_src[1] (unused)
	//	tokens[1] = ll_src[2] = header length in words
	//	tokens[2] = ll_src[3] = -100 for long form (negative), or length for short form (positive)
	//	tokens[3] = ll_src[4] = low 15 bits of list length (long form)
	//	tokens[4] = ll_src[5] = high bits of list length (long form)
	var lllen int
	var dataStart int
	if tokens[2] > 0 {
		// Short form: tokens[2] is the list length directly, data starts at 0-based index 3.
		lllen = int(tokens[2])
		dataStart = 3
	} else {
		// Long form: magic -100 in tokens[2], length in tokens[3..4].
		if len(tokens) < 7 {
			return fmt.Errorf("%w: PLIO_1 long-form header too short (%d words)", ErrCorrupt, len(tokens))
		}
		lllen = (int(tokens[4]) << 15) | (int(tokens[3]) & 0x7fff)
		dataStart = int(tokens[1]) // 1-based llfirt = ll_src[2]+1 → 0-based = tokens[1]
	}
	if lllen > len(tokens) {
		return fmt.Errorf("%w: PLIO_1 declared length %d > stream length %d", ErrCorrupt, lllen, len(tokens))
	}

	// Decode state — cfitsio pl_l2pi initializes pv=1, x1=1, op=1 (all 1-based).
	pv := int32(1)
	x := 1  // 1-based x position
	xe := nelem
	out := make([]int32, nelem)

	i := dataStart
	for i < lllen && x <= xe {
		word := uint16(tokens[i])
		i++
		opcode := int(word>>12) & 0xf
		data := int(word & 0xfff)

		switch opcode {
		case 0:
			// Run of `data` zeros (out is already zero-initialized).
			x2 := x + data - 1
			if x2 > xe {
				x2 = xe
			}
			x = x2 + 1
		case 1:
			// Set PV: low 12 bits = data, high bits in the next word.
			if i >= lllen {
				return fmt.Errorf("%w: PLIO_1 opcode 1 at end of stream", ErrCorrupt)
			}
			hi := int32(tokens[i])
			i++
			pv = (hi << 12) | int32(data)
		case 2:
			pv += int32(data)
		case 3:
			pv -= int32(data)
		case 4:
			// Run of `data` pixels with value PV.
			x2 := x + data - 1
			i2 := x2
			if i2 > xe {
				i2 = xe
			}
			for k := x; k <= i2; k++ {
				out[k-1] = pv
			}
			x = x2 + 1
		case 5:
			// Run of `data` zeros, with PV at the last pixel of the run
			// (if the run reaches its intended end within xe).
			x2 := x + data - 1
			if x2 >= 1 && x2 <= xe {
				out[x2-1] = pv
			}
			x = x2 + 1
		case 6:
			// PV += data, emit single pixel.
			pv += int32(data)
			if x >= 1 && x <= xe {
				out[x-1] = pv
			}
			x++
		case 7:
			// PV -= data, emit single pixel.
			pv -= int32(data)
			if x >= 1 && x <= xe {
				out[x-1] = pv
			}
			x++
		default:
			return fmt.Errorf("%w: PLIO_1 unknown opcode %d at word %d", ErrCorrupt, opcode, i-1)
		}
	}

	// Pack out[] into dst[] at the caller's element size.
	switch elemSize {
	case 1:
		for k := 0; k < nelem; k++ {
			dst[k] = byte(out[k])
		}
	case 2:
		for k := 0; k < nelem; k++ {
			binary.BigEndian.PutUint16(dst[k*2:], uint16(out[k]))
		}
	case 4:
		for k := 0; k < nelem; k++ {
			binary.BigEndian.PutUint32(dst[k*4:], uint32(out[k]))
		}
	default:
		return fmt.Errorf("%w: PLIO_1 unsupported elemSize %d", ErrCorrupt, elemSize)
	}
	return nil
}
