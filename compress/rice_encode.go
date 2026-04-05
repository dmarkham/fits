package compress

import (
	"fmt"
)

// RICE_1 encoder — Rice coding with per-block k parameter. Symmetric
// counterpart of the RICE_1 decoder in rice.go; follows cfitsio/
// ricecomp.c's fits_rcomp / fits_rcomp_short / fits_rcomp_byte.
//
// The encoder:
//
//  1. Writes the first pixel as initBits raw bits.
//  2. For each block of BLOCKSIZE pixels: computes diffs relative to a
//     running "last" value, zigzag-encodes to non-negative, selects an
//     optimal fs parameter based on the block's average magnitude, then
//     writes the block with one of three encodings:
//       - fs = 0 + all-zero block: a single 0 in the fs field, no data
//       - fs >= fsmax: literal escape (stored value = fsmax+1), every
//         pixel written as bbits bits
//       - normal: stored value = fs+1, each pixel as top zeros + 1 + fs
//         low bits
//
// This port is straightforward except for the optimization cfitsio uses
// where the bit buffer lives in two local variables (`lbitbuffer`,
// `lbits_to_go`) inside the hot inner loop. We keep the same structure
// to simplify review against the C reference.

type rice1Encoder struct {
	blockSize int
}

func newRICE1Encoder(params Params) Encoder {
	bs := int(params.Get("BLOCKSIZE", 32))
	if bs <= 0 {
		bs = 32
	}
	return rice1Encoder{blockSize: bs}
}

func (r rice1Encoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	if nelem == 0 {
		return 0, nil
	}
	// fsbits / fsmax per BYTEPIX (from ricecomp.c):
	var fsbits, fsmax int
	switch elemSize {
	case 1:
		fsbits = 3
		fsmax = 6
	case 2:
		fsbits = 4
		fsmax = 14
	case 4:
		fsbits = 5
		fsmax = 25
	case 8:
		fsbits = 6
		fsmax = 57
	default:
		return 0, fmt.Errorf("%w: RICE_1 unsupported BYTEPIX %d", ErrCorrupt, elemSize)
	}
	bbits := 1 << uint(fsbits)

	// Read pixel k from src as a signed integer.
	readPixel := func(k int) int64 {
		off := k * elemSize
		switch elemSize {
		case 1:
			return int64(int8(src[off]))
		case 2:
			v := int16(uint16(src[off])<<8 | uint16(src[off+1]))
			return int64(v)
		case 4:
			v := int32(uint32(src[off])<<24 | uint32(src[off+1])<<16 |
				uint32(src[off+2])<<8 | uint32(src[off+3]))
			return int64(v)
		case 8:
			return int64(uint64(src[off])<<56 | uint64(src[off+1])<<48 |
				uint64(src[off+2])<<40 | uint64(src[off+3])<<32 |
				uint64(src[off+4])<<24 | uint64(src[off+5])<<16 |
				uint64(src[off+6])<<8 | uint64(src[off+7]))
		}
		return 0
	}

	bw := &bitWriter{dst: dst}
	initBits := elemSize * 8
	// Write first pixel as initBits bits (BITPIX=8 is unsigned on disk,
	// others signed — the bit pattern is identical either way).
	if err := bw.writeBits(uint64(readPixel(0)), initBits); err != nil {
		return 0, err
	}
	last := readPixel(0)

	diff := make([]uint64, r.blockSize)
	pos := 1 // we've output pixel 0's value; subsequent diffs start at pixel 1
	// cfitsio's loop iterates nelem total diffs including pixel 0's
	// "always zero" self-diff. We mirror that: the first block decodes
	// pos=0 as a zero diff, and we start the block loop from 0 as well.
	// The initial value has already been written — so the first block
	// will encode diff[0] = a[0] - last_init = 0 and then the rest.
	pos = 0

	for pos < nelem {
		blockLen := r.blockSize
		if pos+blockLen > nelem {
			blockLen = nelem - pos
		}
		var pixelsum float64
		for j := 0; j < blockLen; j++ {
			next := readPixel(pos + j)
			pd := next - last
			// Zigzag: (pd<<1) for pd>=0, (~(pd<<1)) for pd<0.
			var d uint64
			if pd >= 0 {
				d = uint64(pd << 1)
			} else {
				d = uint64(^(pd << 1))
			}
			diff[j] = d
			pixelsum += float64(d)
			last = next
		}
		// Compute fs from the block average.
		dpsum := (pixelsum - float64(blockLen/2) - 1) / float64(blockLen)
		if dpsum < 0 {
			dpsum = 0
		}
		psum := uint(dpsum) >> 1
		fs := 0
		for psum > 0 {
			fs++
			psum >>= 1
		}

		switch {
		case fs >= fsmax:
			// Literal block: fs field = fsmax+1, then bbits bits per diff.
			if err := bw.writeBits(uint64(fsmax+1), fsbits); err != nil {
				return 0, err
			}
			for j := 0; j < blockLen; j++ {
				if err := bw.writeBits(diff[j], bbits); err != nil {
					return 0, err
				}
			}
		case fs == 0 && pixelsum == 0:
			// All-zero block: fs field = 0, no data follows.
			if err := bw.writeBits(0, fsbits); err != nil {
				return 0, err
			}
		default:
			// Normal block: fs field = fs+1, then top zeros + 1 + fs low bits per diff.
			if err := bw.writeBits(uint64(fs+1), fsbits); err != nil {
				return 0, err
			}
			fsmask := uint64(1<<uint(fs)) - 1
			for j := 0; j < blockLen; j++ {
				v := diff[j]
				top := v >> uint(fs)
				// Write `top` zero bits then a 1.
				for top > 0 {
					if err := bw.writeBits(0, 1); err != nil {
						return 0, err
					}
					top--
				}
				if err := bw.writeBits(1, 1); err != nil {
					return 0, err
				}
				if fs > 0 {
					if err := bw.writeBits(v&fsmask, fs); err != nil {
						return 0, err
					}
				}
			}
		}
		pos += blockLen
	}
	return bw.finish(), nil
}

// bitWriter is the MSB-first-within-byte bit writer used by the RICE
// encoder. Bits accumulate in a small buffer and flush to dst a byte
// at a time.
type bitWriter struct {
	dst       []byte
	byteIdx   int
	bitBuf    uint64 // pending bits, left-justified at bit position bitLen-1
	bitLen    int    // number of pending bits
	overflow  bool
}

func (w *bitWriter) writeBits(v uint64, n int) error {
	w.bitBuf = (w.bitBuf << uint(n)) | (v & ((1 << uint(n)) - 1))
	w.bitLen += n
	for w.bitLen >= 8 {
		w.bitLen -= 8
		b := byte((w.bitBuf >> uint(w.bitLen)) & 0xff)
		if w.byteIdx >= len(w.dst) {
			w.overflow = true
			return ErrBufferTooSmall
		}
		w.dst[w.byteIdx] = b
		w.byteIdx++
	}
	return nil
}

// finish flushes any remaining bits (padded with zeros on the right) and
// returns the total number of bytes written.
func (w *bitWriter) finish() int {
	if w.bitLen > 0 {
		b := byte((w.bitBuf << uint(8-w.bitLen)) & 0xff)
		if w.byteIdx < len(w.dst) {
			w.dst[w.byteIdx] = b
			w.byteIdx++
		}
		w.bitLen = 0
	}
	return w.byteIdx
}
