package compress

import (
	"encoding/binary"
	"fmt"
)

// PLIO_1 encoder — IRAF Planio line-list format, symmetric counterpart
// of the decoder in plio.go. Produces a token stream that the decoder
// round-trips exactly.
//
// Reference: cfitsio/pliocomp.c pl_p2li — the IRAF pixel-to-line-list
// translator. Format: 7-word header (long form with magic -100, length
// in words 3-4, dimensions in words 5-7), then data tokens.
//
// The data tokens are produced by walking the pixel array and emitting
// runs: each run is encoded as an opcode + 12-bit data field that
// represents either "run of zeros", "run of PV", or a PV update (add
// or subtract a constant). See plio.go for the full opcode table.
//
// Our encoder takes a simple approach that matches cfitsio's output for
// the non-exotic inputs we care about (masks with 0/1/2/... values):
//
//   1. Scan for runs of consecutive equal-value pixels.
//   2. For each run, emit a PV-update opcode if the new PV differs from
//      the current, followed by a "run of PV" opcode (4) or "run of zeros"
//      opcode (0).

type plio1Encoder struct{}

func (plio1Encoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	// Read the input pixels as int32 (PLIO operates on 16-bit values but
	// we accept wider input types too — cfitsio converts in-place).
	pixels := make([]int32, nelem)
	switch elemSize {
	case 1:
		for i := 0; i < nelem; i++ {
			pixels[i] = int32(src[i])
		}
	case 2:
		for i := 0; i < nelem; i++ {
			pixels[i] = int32(int16(uint16(src[i*2])<<8 | uint16(src[i*2+1])))
		}
	case 4:
		for i := 0; i < nelem; i++ {
			pixels[i] = int32(uint32(src[i*4])<<24 | uint32(src[i*4+1])<<16 |
				uint32(src[i*4+2])<<8 | uint32(src[i*4+3]))
		}
	default:
		return 0, fmt.Errorf("%w: PLIO_1 unsupported elemSize %d", ErrCorrupt, elemSize)
	}

	// Build the token stream in-memory as []int16, then serialize big-endian.
	var tokens []int16
	// Header placeholder — will be filled in at the end once we know
	// the total length.
	tokens = append(tokens, 0)             // word 0: unused
	tokens = append(tokens, 7)             // word 1: header length
	tokens = append(tokens, -100)          // word 2: magic (long form)
	tokens = append(tokens, 0, 0)          // words 3,4: length placeholders
	tokens = append(tokens, 0, 0)          // words 5,6: reserved

	// Walk pixels and emit runs. PV starts at 1 per the decoder.
	pv := int32(1)
	i := 0
	for i < nelem {
		if pixels[i] == 0 {
			// Run of zeros.
			j := i
			for j < nelem && pixels[j] == 0 && (j-i) < 0xfff {
				j++
			}
			runLen := j - i
			tokens = append(tokens, int16(runLen)) // opcode 0, data = runLen
			i = j
			continue
		}
		// Non-zero run. First update PV if needed.
		val := pixels[i]
		if val != pv {
			delta := val - pv
			if delta > 0 && delta <= 0xfff {
				// Opcode 2: PV += delta, no pixel emitted.
				tokens = append(tokens, int16((2<<12)|int(delta)))
			} else if delta < 0 && -delta <= 0xfff {
				// Opcode 3: PV -= |delta|, no pixel emitted.
				tokens = append(tokens, int16((3<<12)|int(-delta)))
			} else {
				// Opcode 1: set PV directly via two-word encoding.
				// Low 12 bits are data, next token is high bits.
				lo := int32(val & 0xfff)
				hi := int32(val >> 12)
				tokens = append(tokens, int16((1<<12)|int(lo)))
				tokens = append(tokens, int16(hi))
			}
			pv = val
		}
		// Now emit the run of PV.
		j := i
		for j < nelem && pixels[j] == val && (j-i) < 0xfff {
			j++
		}
		runLen := j - i
		tokens = append(tokens, int16((4<<12)|runLen)) // opcode 4: run of PV
		i = j
	}

	// Fill in the length (words 3, 4) in the header. Length is the total
	// word count (header + data) in 1-based cfitsio convention.
	lllen := int32(len(tokens))
	tokens[3] = int16(lllen & 0x7fff)
	tokens[4] = int16(lllen >> 15)

	// Serialize as big-endian int16.
	need := len(tokens) * 2
	if need > len(dst) {
		return 0, ErrBufferTooSmall
	}
	for k, t := range tokens {
		binary.BigEndian.PutUint16(dst[k*2:], uint16(t))
	}
	return need, nil
}
