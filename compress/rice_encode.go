package compress

import (
	"fmt"
)

// RICE_1 encoder — Rice coding with per-block k parameter. Symmetric
// counterpart of the RICE_1 decoder in rice.go; a port of cfitsio's
// ricecomp.c fits_rcomp (BYTEPIX 4), fits_rcomp_short (BYTEPIX 2) and
// fits_rcomp_byte (BYTEPIX 1). Its output is byte-identical to those
// functions in libcfitsio 4.7.0 (TestRiceEncodeMatchesCfitsio, golden
// streams from compress/testdata/cref).
//
// The encoder:
//
//  1. Writes the first pixel as BYTEPIX*8 raw bits.
//  2. For each block of BLOCKSIZE pixels: computes the difference to the
//     previous pixel IN THE PIXEL WIDTH (cfitsio's pdiff is a signed char
//     / short / int, so a jump of more than half the range wraps exactly
//     like the decoder's unsigned addition undoes it), zigzag-maps it to
//     a non-negative value, selects fs from the block's mean magnitude,
//     then writes the block with one of three encodings:
//       - fs = 0 + all-zero block: a single 0 in the fs field, no data
//       - fs >= fsmax: literal escape (stored value = fsmax+1), every
//         pixel written as bbits bits
//       - normal: stored value = fs+1, each pixel as top zeros + 1 + fs
//         low bits
//
// BYTEPIX 8 is refused: cfitsio has no 64-bit Rice coder, so there is no
// reference for such a stream (and cfitsio cannot write compressed
// BITPIX 64 images at all).

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
	default:
		return 0, fmt.Errorf("%w: RICE_1 encoder supports BYTEPIX 1, 2, 4 (as cfitsio does), not %d", ErrUnsupportedCompression, elemSize)
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
		default:
			v := int32(uint32(src[off])<<24 | uint32(src[off+1])<<16 |
				uint32(src[off+2])<<8 | uint32(src[off+3]))
			return int64(v)
		}
	}
	// wrap narrows a difference to the signed pixel width, as the C
	// assignment to cfitsio's `signed char` / `short` / `int` pdiff does.
	wrap := func(d int64) int64 {
		switch elemSize {
		case 1:
			return int64(int8(d))
		case 2:
			return int64(int16(d))
		default:
			return int64(int32(d))
		}
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
	// cfitsio's loop iterates nelem total diffs including pixel 0's
	// "always zero" self-diff. We mirror that: the initial value has
	// already been written, so the first block encodes
	// diff[0] = a[0] - a[0] = 0 and then the rest.
	pos := 0

	for pos < nelem {
		blockLen := r.blockSize
		if pos+blockLen > nelem {
			blockLen = nelem - pos
		}
		var pixelsum float64
		for j := 0; j < blockLen; j++ {
			next := readPixel(pos + j)
			pd := wrap(next - last)
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
	return bw.finish()
}

// bitWriter is the MSB-first-within-byte bit writer used by the RICE
// encoder. Bits accumulate in a small buffer and flush to dst a byte
// at a time.
type bitWriter struct {
	dst      []byte
	byteIdx  int
	bitBuf   uint64 // pending bits, left-justified at bit position bitLen-1
	bitLen   int    // number of pending bits
	overflow bool
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
// returns the total number of bytes written, or ErrBufferTooSmall when
// the final partial byte does not fit (never a truncated stream).
func (w *bitWriter) finish() (int, error) {
	if w.bitLen > 0 {
		b := byte((w.bitBuf << uint(8-w.bitLen)) & 0xff)
		if w.byteIdx >= len(w.dst) {
			w.overflow = true
			return 0, ErrBufferTooSmall
		}
		w.dst[w.byteIdx] = b
		w.byteIdx++
		w.bitLen = 0
	}
	return w.byteIdx, nil
}
