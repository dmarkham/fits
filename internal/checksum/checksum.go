// Package checksum implements the FITS checksum algorithm defined in
// Appendix J of the FITS Standard v4.0.
//
// The FITS checksum is a 32-bit 1's-complement sum computed over the
// big-endian uint16 view of an HDU's 2880-byte-aligned bytes. Overflow bits
// are folded back into the low 16 bits so that every bit position is sampled
// evenly. The final 32-bit value is ASCII-encoded via Rob Seaman's encoding
// so that the CHECKSUM keyword value is printable and round-trippable.
//
// Port of cfitsio checksum.c (ffcsum / ffesum / ffdsum). Algorithm verbatim,
// idioms Go.
package checksum

import (
	"encoding/binary"
)

// Update folds the 2880-byte block b into the running 32-bit 1's-complement
// checksum sum and returns the new value. It is the caller's responsibility
// to ensure b is exactly 2880 bytes.
func Update(sum uint32, b []byte) uint32 {
	// Read each block as 1440 big-endian uint16 values and accumulate into two
	// 32-bit halves (hi: even-indexed shorts, lo: odd-indexed shorts). Fold
	// overflow back in afterwards.
	hi := uint32(sum >> 16)
	lo := uint32(sum & 0xFFFF)
	for i := 0; i < len(b); i += 4 {
		hi += uint32(binary.BigEndian.Uint16(b[i : i+2]))
		lo += uint32(binary.BigEndian.Uint16(b[i+2 : i+4]))
	}
	hicarry := hi >> 16
	locarry := lo >> 16
	for hicarry|locarry != 0 {
		hi = (hi & 0xFFFF) + locarry
		lo = (lo & 0xFFFF) + hicarry
		hicarry = hi >> 16
		locarry = lo >> 16
	}
	return (hi << 16) + lo
}

// Compute returns the 32-bit 1's-complement checksum of the supplied bytes.
// The length must be a multiple of 4 (i.e. the caller has already padded the
// input to the block boundary as FITS requires).
func Compute(data []byte) uint32 {
	var sum uint32
	return Update(sum, data)
}

// excludeBytes are the ASCII punctuation characters between the digits and
// letters that cfitsio/FITS skip when encoding.
var excludeBytes = [13]byte{
	0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f, 0x40,
	0x5b, 0x5c, 0x5d, 0x5e, 0x5f, 0x60,
}

// Encode converts a 32-bit 1's-complement checksum into its 16-character
// ASCII representation per Appendix J.
//
// If complement is true the function encodes the bit-wise complement of sum,
// which is the form written to the CHECKSUM keyword.
func Encode(sum uint32, complement bool) [16]byte {
	value := sum
	if complement {
		value = 0xFFFFFFFF - sum
	}
	masks := [4]uint32{0xFF000000, 0x00FF0000, 0x0000FF00, 0x000000FF}
	var asc [16]byte
	const offset = 0x30 // '0'

	for ii := range 4 {
		byteVal := int((value & masks[ii]) >> (24 - (8 * ii)))
		quotient := byteVal/4 + offset
		remainder := byteVal % 4
		var ch [4]int
		for jj := range 4 {
			ch[jj] = quotient
		}
		ch[0] += remainder

		for {
			check := 0
			for _, ex := range excludeBytes {
				for jj := 0; jj < 4; jj += 2 {
					if byte(ch[jj]) == ex || byte(ch[jj+1]) == ex {
						ch[jj]++
						ch[jj+1]--
						check++
					}
				}
			}
			if check == 0 {
				break
			}
		}

		for jj := range 4 {
			asc[4*jj+ii] = byte(ch[jj])
		}
	}

	// Shift the bytes 1 to the right (permutation).
	var out [16]byte
	for ii := range 16 {
		out[ii] = asc[(ii+15)%16]
	}
	return out
}

// Decode is the inverse of Encode. It returns the 32-bit checksum sum that
// produced the given 16-character ASCII encoding.
func Decode(ascii [16]byte, complement bool) uint32 {
	var cbuf [16]int
	for ii := range 16 {
		cbuf[ii] = int(ascii[(ii+1)%16]) - 0x30
	}
	var hi, lo uint32
	for ii := 0; ii < 16; ii += 4 {
		hi += uint32(cbuf[ii]<<8 + cbuf[ii+1])
		lo += uint32(cbuf[ii+2]<<8 + cbuf[ii+3])
	}
	hicarry := hi >> 16
	locarry := lo >> 16
	for hicarry|locarry != 0 {
		hi = (hi & 0xFFFF) + locarry
		lo = (lo & 0xFFFF) + hicarry
		hicarry = hi >> 16
		locarry = lo >> 16
	}
	sum := (hi << 16) + lo
	if complement {
		sum = 0xFFFFFFFF - sum
	}
	return sum
}
