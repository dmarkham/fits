package compress

import "fmt"

// debugRice is a compile-time flag for decoder debugging. Flip to true
// and re-run the RICE unit tests to trace per-block fs values.
const debugRice = false

// RICE_1 — Rice coding with per-block k parameter.
//
// Reference: Pence 2000 (cfitsio ricecomp.c) and White & Greenfield 1998.
// Adapted to Go idiom. The algorithm:
//
//	1. Read the "last" pixel value from the stream as the initial value
//	   (elemSize * 8 bits).
//	2. For each block of BLOCKSIZE pixels (default 32):
//	   a. Read the fs parameter (k) from the stream: a fixed number of
//	      bits equal to fsbits. fsbits is 5 for BYTEPIX<=4, 6 for BYTEPIX=8.
//	   b. If fs < 0 (all bits set to 1), the block is encoded literally:
//	      read BLOCKSIZE * elemSize bits verbatim as raw pixel values.
//	   c. Otherwise, decode each pixel via the fs-split code:
//	      - Read a unary-coded prefix (count of 0 bits up to the first 1).
//	      - Read fs low-order bits as a fixed-width field.
//	      - Combine: value = (prefix << fs) | lowBits.
//	      - Map back from zigzag to signed delta.
//	      - Add to running last = last + delta.
//	3. Write each pixel as big-endian bytes into dst.
//
// This is stateful per-stream but resets between tiles (each tile row
// of the binary table is an independent RICE stream).
type rice1Decoder struct {
	blockSize int // pixels per block (BLOCKSIZE, default 32)
}

func newRICE1(params Params) Decoder {
	bs := int(params.Get("BLOCKSIZE", 32))
	if bs <= 0 {
		bs = 32
	}
	return rice1Decoder{blockSize: bs}
}

func (r rice1Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	if nelem == 0 {
		return nil
	}
	if len(dst) < nelem*elemSize {
		return fmt.Errorf("%w: RICE_1 dst length %d < needed %d", ErrCorrupt, len(dst), nelem*elemSize)
	}
	// fsbits and fsmax depend on BYTEPIX per Pence 2000.
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
		return fmt.Errorf("%w: RICE_1 unsupported BYTEPIX %d", ErrCorrupt, elemSize)
	}

	br := &bitReader{src: src}
	// Read the initial "last" value — elemSize * 8 bits. The initial
	// value is used as the baseline for the first delta but is NOT
	// written to the output directly. Instead, the first block's first
	// diff will be (a[0] - last_init) which, for cfitsio's encoder, is
	// exactly zero (the encoder initializes last = a[0] before
	// encoding). The decoder computes array[0] = last + first_diff,
	// which gives back a[0].
	initBits := elemSize * 8
	lastVal, err := br.readBits(initBits)
	if err != nil {
		return fmt.Errorf("%w: RICE_1 init: %v", ErrCorrupt, err)
	}
	last := int64(lastVal)
	// Sign-extend from the top bit for signed integer BITPIX.
	if initBits < 64 {
		shift := uint(64 - initBits)
		last = (last << shift) >> shift
	}

	writeOne := func(idx int, v int64) {
		off := idx * elemSize
		switch elemSize {
		case 1:
			dst[off] = byte(v)
		case 2:
			dst[off] = byte(v >> 8)
			dst[off+1] = byte(v)
		case 4:
			dst[off] = byte(v >> 24)
			dst[off+1] = byte(v >> 16)
			dst[off+2] = byte(v >> 8)
			dst[off+3] = byte(v)
		case 8:
			dst[off] = byte(v >> 56)
			dst[off+1] = byte(v >> 48)
			dst[off+2] = byte(v >> 40)
			dst[off+3] = byte(v >> 32)
			dst[off+4] = byte(v >> 24)
			dst[off+5] = byte(v >> 16)
			dst[off+6] = byte(v >> 8)
			dst[off+7] = byte(v)
		}
	}

	// cfitsio encodes ALL nelem pixels as deltas relative to the running
	// "last" value. The first diff decodes to (array[0] - last_init) and
	// is typically zero because the encoder initializes last = array[0]
	// before the encoding loop starts. We decode every position in the
	// block loop — nothing is written before pos=0.
	pos := 0
	for pos < nelem {
		// Read the block's fs parameter. cfitsio's encoder writes:
		//   0            → all-zero block: no diff data, all pixels equal last
		//   fsmax+1      → literal block: each diff stored raw
		//   1 .. fsmax   → normal block with actual fs = stored - 1
		fsU, err := br.readBits(fsbits)
		if err != nil {
			return fmt.Errorf("%w: RICE_1 fs read at pixel %d: %v", ErrCorrupt, pos, err)
		}
		fsRaw := int(fsU)
		if debugRice {
			fmt.Printf("block @ pos=%d fsRaw=%d\n", pos, fsRaw)
		}
		blockLen := r.blockSize
		if pos+blockLen > nelem {
			blockLen = nelem - pos
		}

		// Dispatch on the fs field encoding.
		switch {
		case fsRaw == 0:
			// All-zero block: every diff is zero, so every pixel equals
			// the running `last` value. No diff data to read.
			for i := 0; i < blockLen; i++ {
				writeOne(pos+i, last)
			}
		case fsRaw == fsmax+1:
			// Literal block: each diff stored as bbits = initBits bits.
			for i := 0; i < blockLen; i++ {
				v, err := br.readBits(initBits)
				if err != nil {
					return fmt.Errorf("%w: RICE_1 literal read: %v", ErrCorrupt, err)
				}
				// The written value is the raw zigzag-encoded diff.
				delta := unzigzag(v)
				last += delta
				writeOne(pos+i, last)
			}
		default:
			// Normal block: actual fs = fsRaw - 1.
			fs := fsRaw - 1
			// fs-coded block. The unary prefix length is unbounded in
			// general — for tile row boundaries on multi-dim tiles the
			// diff can be much larger than fsmax permits. The cap here
			// is just a sanity bound against corrupted input (8*buffer
			// bytes is the absolute maximum number of bits possible).
			unaryCap := len(src) * 8
			for i := 0; i < blockLen; i++ {
				top, err := br.readUnary(unaryCap)
				if err != nil {
					return fmt.Errorf("%w: RICE_1 unary read at pixel %d: %v", ErrCorrupt, pos+i, err)
				}
				var low uint64
				if fs > 0 {
					low, err = br.readBits(fs)
					if err != nil {
						return fmt.Errorf("%w: RICE_1 low read: %v", ErrCorrupt, err)
					}
				}
				zigzag := (uint64(top) << uint(fs)) | low
				delta := unzigzag(zigzag)
				last += delta
				writeOne(pos+i, last)
			}
		}
		pos += blockLen
	}
	return nil
}

// unzigzag maps a non-negative zigzag-encoded integer back to its signed
// delta. cfitsio uses the mapping:
//
//	0 → 0
//	1 → -1
//	2 → 1
//	3 → -2
//	4 → 2
//	...
//
// Equivalently: if zigzag is odd, -(zigzag+1)/2; if even, zigzag/2.
func unzigzag(z uint64) int64 {
	if z&1 != 0 {
		return -int64((z + 1) / 2)
	}
	return int64(z / 2)
}

// bitReader reads variable-width bit fields from a byte stream in MSB-first
// order. cfitsio's RICE encoder writes bits MSB-first within each byte,
// which is the convention we follow here.
type bitReader struct {
	src  []byte
	byteIdx int
	bitIdx  uint // 0..7, bits consumed from the current byte (MSB-first)
}

// readBits consumes n bits (0 <= n <= 64) and returns them right-aligned
// in the returned uint64.
func (br *bitReader) readBits(n int) (uint64, error) {
	if n == 0 {
		return 0, nil
	}
	if n > 64 {
		return 0, fmt.Errorf("bitReader: readBits n=%d > 64", n)
	}
	var v uint64
	for n > 0 {
		if br.byteIdx >= len(br.src) {
			return 0, fmt.Errorf("bitReader: unexpected EOF at byte %d bit %d", br.byteIdx, br.bitIdx)
		}
		// How many bits remain in the current byte?
		remain := 8 - int(br.bitIdx)
		take := remain
		if take > n {
			take = n
		}
		// Extract the top `take` bits from the current byte position.
		b := br.src[br.byteIdx]
		shift := uint(remain - take)
		chunk := uint64((b >> shift) & ((1 << uint(take)) - 1))
		v = (v << uint(take)) | chunk
		br.bitIdx += uint(take)
		if br.bitIdx == 8 {
			br.byteIdx++
			br.bitIdx = 0
		}
		n -= take
	}
	return v, nil
}

// readUnary reads a unary-coded value: counts consecutive 0 bits until
// the first 1 bit. Returns the count. Maximum allowed count is cap — if
// no 1 bit is seen within cap bits, returns an error.
func (br *bitReader) readUnary(cap int) (int, error) {
	count := 0
	for count <= cap {
		b, err := br.readBits(1)
		if err != nil {
			return 0, err
		}
		if b == 1 {
			return count, nil
		}
		count++
	}
	return 0, fmt.Errorf("bitReader: unary run exceeds cap=%d", cap)
}
