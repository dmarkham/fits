package compress

import (
	"encoding/binary"
	"fmt"
)

// HCOMPRESS_1 — H-transform compression with quadtree coding.
//
// Reference: cfitsio/fits_hdecompress.c (White 1992, "High-Performance
// Compression of Astronomical Images", Journal of ASP 104:630). This Go
// port follows the cfitsio reference function-for-function.
//
// # Stream format (as emitted by cfitsio's fits_hcompress):
//
//	bytes 0-1     : magic 0xDD 0x99
//	bytes 2-5     : nx (4-byte big-endian int, number of rows)
//	bytes 6-9     : ny (fastest axis — ny varies fastest in the output)
//	bytes 10-13   : scale (digitization scale factor)
//	bytes 14-21   : sumall (8-byte big-endian int, sum of all pixels)
//	byte 22       : nbitplanes[0] — top-level quadrant bit-plane count
//	byte 23       : nbitplanes[1] — X/Y quadrants (shared)
//	byte 24       : nbitplanes[2] — XY quadrant
//	bytes 25..    : 4 quadtree-coded quadrants
//	              : nybble-packed EOF marker
//	              : sign bits for non-zero pixels
//
// # Algorithm
//
//  1. qtree_decode for each of the 4 quadrants reconstructs the H-transform
//     coefficient array bit-plane by bit-plane, starting from the most
//     significant bit and working down.
//  2. Sign bits are then read for every non-zero coefficient.
//  3. a[0] is replaced by sumall (the DC level, which the encoder wrote
//     separately for precision).
//  4. undigitize multiplies every coefficient by scale (if scale > 1).
//  5. hinv applies the inverse H-transform via log2n butterfly passes,
//     producing the reconstructed image in row-major order.
//
// Note on axis order: cfitsio's HCOMPRESS convention has ny as the fastest
// varying dimension, opposite the FITS NAXIS1/NAXIS2 intuition. When
// astropy writes an HCOMPRESS tile it follows this same convention.

var hcompressMagic = [2]byte{0xDD, 0x99}

// hcompress1Decoder decodes HCOMPRESS_1 tiles.
type hcompress1Decoder struct {
	smooth int // 0 = no smoothing, else smooth during inversion
}

func newHCOMPRESS1(params Params) Decoder {
	return hcompress1Decoder{
		smooth: int(params.Get("SMOOTH", 0)),
	}
}

func (h hcompress1Decoder) Decode(src, dst []byte, nelem, elemSize int) error {
	// Allocate the working int32 array with enough headroom.
	work := make([]int32, nelem)
	nx, ny, scale, err := decodeHCOMPRESS(src, work, h.smooth)
	if err != nil {
		return err
	}
	if nx*ny != nelem {
		return fmt.Errorf("%w: HCOMPRESS_1 stream declares %dx%d=%d pixels, tile expects %d",
			ErrCorrupt, nx, ny, nx*ny, nelem)
	}
	_ = scale // only used internally by hinv

	// Pack work into dst at the caller's element size (big-endian FITS).
	switch elemSize {
	case 1:
		for i := 0; i < nelem; i++ {
			dst[i] = byte(work[i])
		}
	case 2:
		for i := 0; i < nelem; i++ {
			binary.BigEndian.PutUint16(dst[i*2:], uint16(int16(work[i])))
		}
	case 4:
		for i := 0; i < nelem; i++ {
			binary.BigEndian.PutUint32(dst[i*4:], uint32(work[i]))
		}
	case 8:
		for i := 0; i < nelem; i++ {
			binary.BigEndian.PutUint64(dst[i*8:], uint64(int64(work[i])))
		}
	default:
		return fmt.Errorf("%w: HCOMPRESS_1 unsupported elemSize %d", ErrCorrupt, elemSize)
	}
	return nil
}

// decodeHCOMPRESS is the top-level entry point: parses the stream header,
// quadtree-decodes the 4 bit-plane quadrants, reads sign bits, restores
// the sum-of-all DC level, undigitizes, and runs the inverse H-transform.
// Returns (nx, ny, scale).
func decodeHCOMPRESS(src []byte, a []int32, smooth int) (int, int, int, error) {
	r := &hcompressReader{src: src}
	// Magic.
	var magic [2]byte
	if err := r.readBytes(magic[:]); err != nil {
		return 0, 0, 0, fmt.Errorf("%w: HCOMPRESS_1: %v", ErrCorrupt, err)
	}
	if magic != hcompressMagic {
		return 0, 0, 0, fmt.Errorf("%w: HCOMPRESS_1 bad magic %#x %#x", ErrCorrupt, magic[0], magic[1])
	}
	// nx, ny, scale (4 bytes each, big-endian).
	nx, err := r.readInt32BE()
	if err != nil {
		return 0, 0, 0, err
	}
	ny, err := r.readInt32BE()
	if err != nil {
		return 0, 0, 0, err
	}
	scale, err := r.readInt32BE()
	if err != nil {
		return 0, 0, 0, err
	}
	// sumall (8 bytes).
	sumall, err := r.readInt64BE()
	if err != nil {
		return 0, 0, 0, err
	}
	// nbitplanes[3].
	var nbitplanes [3]byte
	if err := r.readBytes(nbitplanes[:]); err != nil {
		return 0, 0, 0, err
	}

	if int(nx)*int(ny) > len(a) {
		return 0, 0, 0, fmt.Errorf("%w: HCOMPRESS_1 pixel count %d exceeds buffer %d",
			ErrCorrupt, int(nx)*int(ny), len(a))
	}
	// a is already zero per Go make.
	if err := dodecodeHCOMPRESS(r, a, int(nx), int(ny), nbitplanes); err != nil {
		return 0, 0, 0, err
	}
	// Replace a[0] with sumall (the DC level).
	a[0] = int32(sumall)

	// Undigitize: multiply by scale.
	if scale > 1 {
		s := int32(scale)
		for i := range a[:int(nx)*int(ny)] {
			a[i] *= s
		}
	}
	// Inverse H-transform.
	if err := hinvHCOMPRESS(a, int(nx), int(ny), smooth, int(scale)); err != nil {
		return 0, 0, 0, err
	}
	return int(nx), int(ny), int(scale), nil
}

// hcompressReader is a combined byte + bit reader for the HCOMPRESS
// stream. Bit reading is restarted twice (once for the quadtree body,
// once for sign bits) so we keep byte cursor and bit buffer separate
// and expose a start() method to flush the bit buffer.
type hcompressReader struct {
	src     []byte
	byteIdx int
	bitBuf  int32 // pending bits (used by input_nbits/input_nybble/input_bit)
	bitsTG  int   // bits remaining in bitBuf
}

func (r *hcompressReader) readBytes(dst []byte) error {
	if r.byteIdx+len(dst) > len(r.src) {
		return fmt.Errorf("%w: HCOMPRESS_1 unexpected EOF", ErrCorrupt)
	}
	copy(dst, r.src[r.byteIdx:r.byteIdx+len(dst)])
	r.byteIdx += len(dst)
	return nil
}

func (r *hcompressReader) readInt32BE() (int32, error) {
	var b [4]byte
	if err := r.readBytes(b[:]); err != nil {
		return 0, err
	}
	return int32(binary.BigEndian.Uint32(b[:])), nil
}

func (r *hcompressReader) readInt64BE() (int64, error) {
	var b [8]byte
	if err := r.readBytes(b[:]); err != nil {
		return 0, err
	}
	return int64(binary.BigEndian.Uint64(b[:])), nil
}

// startBits resets the bit-reader state. Called once before the quadtree
// decoding and once before the sign-bit pass.
func (r *hcompressReader) startBits() {
	r.bitsTG = 0
	r.bitBuf = 0
}

// inputBit returns the next single bit from the stream.
func (r *hcompressReader) inputBit() int {
	if r.bitsTG == 0 {
		if r.byteIdx >= len(r.src) {
			return 0
		}
		r.bitBuf = int32(r.src[r.byteIdx])
		r.byteIdx++
		r.bitsTG = 8
	}
	r.bitsTG--
	return int((r.bitBuf >> uint(r.bitsTG)) & 1)
}

// inputNbits returns the next n bits (n <= 8) right-aligned.
func (r *hcompressReader) inputNbits(n int) int {
	if r.bitsTG < n {
		if r.byteIdx >= len(r.src) {
			return 0
		}
		r.bitBuf = (r.bitBuf << 8) | int32(r.src[r.byteIdx])
		r.byteIdx++
		r.bitsTG += 8
	}
	r.bitsTG -= n
	mask := int32((1 << uint(n)) - 1)
	return int((r.bitBuf >> uint(r.bitsTG)) & mask)
}

// inputNybble returns the next 4 bits.
func (r *hcompressReader) inputNybble() int {
	return r.inputNbits(4)
}

// inputNnybble reads n 4-bit nybbles into dst[0..n-1]. Mirrors the
// cfitsio bulk reader (with the same backspace trick for cases where
// bits_to_go==8 at entry).
func (r *hcompressReader) inputNnybble(n int, dst []byte) {
	if n == 1 {
		dst[0] = byte(r.inputNybble())
		return
	}
	if r.bitsTG == 8 {
		r.byteIdx--
		r.bitsTG = 0
	}
	shift1 := r.bitsTG + 4
	shift2 := r.bitsTG
	kk := 0
	if r.bitsTG == 0 {
		for ii := 0; ii < n/2; ii++ {
			r.bitBuf = (r.bitBuf << 8) | int32(r.src[r.byteIdx])
			r.byteIdx++
			dst[kk] = byte((r.bitBuf >> 4) & 15)
			dst[kk+1] = byte(r.bitBuf & 15)
			kk += 2
		}
	} else {
		for ii := 0; ii < n/2; ii++ {
			r.bitBuf = (r.bitBuf << 8) | int32(r.src[r.byteIdx])
			r.byteIdx++
			dst[kk] = byte((r.bitBuf >> uint(shift1)) & 15)
			dst[kk+1] = byte((r.bitBuf >> uint(shift2)) & 15)
			kk += 2
		}
	}
	if (n/2)*2 != n {
		dst[n-1] = byte(r.inputNybble())
	}
}

// inputHuffman is the fixed Huffman decoder for qtree nybbles. Codes:
//
//	3e → 0 (6 bits)
//	00 → 1 (3 bits)
//	01 → 2 (3 bits)
//	08 → 3 (4 bits)
//	02 → 4 (3 bits)
//	09 → 5 (4 bits)
//	1a → 6 (5 bits)
//	1b → 7 (5 bits)
//	03 → 8 (3 bits)
//	1c → 9 (5 bits)
//	0a → 10 (4 bits)
//	1d → 11 (5 bits)
//	0b → 12 (4 bits)
//	1e → 13 (5 bits)
//	3f → 14 (6 bits)
//	0c → 15 (4 bits)
func (r *hcompressReader) inputHuffman() int {
	c := r.inputNbits(3)
	if c < 4 {
		return 1 << uint(c) // 1, 2, 4, 8
	}
	c = r.inputBit() | (c << 1)
	if c < 13 {
		switch c {
		case 8:
			return 3
		case 9:
			return 5
		case 10:
			return 10
		case 11:
			return 12
		case 12:
			return 15
		}
	}
	c = r.inputBit() | (c << 1)
	if c < 31 {
		switch c {
		case 26:
			return 6
		case 27:
			return 7
		case 28:
			return 9
		case 29:
			return 11
		case 30:
			return 13
		}
	}
	c = r.inputBit() | (c << 1)
	if c == 62 {
		return 0
	}
	return 14
}

// dodecodeHCOMPRESS decodes the 4 H-transform quadrants into a[].
func dodecodeHCOMPRESS(r *hcompressReader, a []int32, nx, ny int, nbitplanes [3]byte) error {
	nel := nx * ny
	nx2 := (nx + 1) / 2
	ny2 := (ny + 1) / 2

	r.startBits()

	// Quadrant LL (top-left, nx2 x ny2 region, nbitplanes[0]).
	if err := qtreeDecode(r, a, 0, ny, nx2, ny2, int(nbitplanes[0])); err != nil {
		return err
	}
	// Quadrant LH (top-right, nx2 x ny/2, nbitplanes[1]).
	if err := qtreeDecode(r, a, ny2, ny, nx2, ny/2, int(nbitplanes[1])); err != nil {
		return err
	}
	// Quadrant HL (bottom-left, nx/2 x ny2, nbitplanes[1]).
	if err := qtreeDecode(r, a, ny*nx2, ny, nx/2, ny2, int(nbitplanes[1])); err != nil {
		return err
	}
	// Quadrant HH (bottom-right, nx/2 x ny/2, nbitplanes[2]).
	if err := qtreeDecode(r, a, ny*nx2+ny2, ny, nx/2, ny/2, int(nbitplanes[2])); err != nil {
		return err
	}
	// EOF nybble — must be 0.
	if r.inputNybble() != 0 {
		return fmt.Errorf("%w: HCOMPRESS_1 bad bitplane EOF nybble", ErrCorrupt)
	}
	// Sign bits. Re-start bit reader.
	r.startBits()
	for i := 0; i < nel; i++ {
		if a[i] != 0 {
			if r.inputBit() != 0 {
				a[i] = -a[i]
			}
		}
	}
	return nil
}

// qtreeDecode decodes one quadrant's bit planes via the quadtree coding.
// base is the index into a[] where this quadrant starts; n is the full
// row stride (= ny); nqx, nqy are the quadrant dimensions; nbitplanes is
// the number of bit planes to decode.
func qtreeDecode(r *hcompressReader, a []int32, base, n, nqx, nqy, nbitplanes int) error {
	if nbitplanes <= 0 {
		return nil
	}
	nqmax := nqx
	if nqy > nqmax {
		nqmax = nqy
	}
	log2n := 0
	for 1<<uint(log2n) < nqmax {
		log2n++
	}
	nqx2 := (nqx + 1) / 2
	nqy2 := (nqy + 1) / 2
	scratch := make([]byte, nqx2*nqy2)

	for bit := nbitplanes - 1; bit >= 0; bit-- {
		b := r.inputNybble()
		if b == 0 {
			// Direct bit-image representation: one nybble per 2x2 block.
			nNybbles := ((nqx + 1) / 2) * ((nqy + 1) / 2)
			r.inputNnybble(nNybbles, scratch)
			qtreeBitins(scratch, nqx, nqy, a[base:], n, bit)
		} else if b != 0xf {
			return fmt.Errorf("%w: HCOMPRESS_1 bad qtree format code %#x", ErrCorrupt, b)
		} else {
			// Quadtree expansion.
			scratch[0] = byte(r.inputHuffman())
			nx := 1
			ny := 1
			nfx := nqx
			nfy := nqy
			c := 1 << uint(log2n)
			for k := 1; k < log2n; k++ {
				c >>= 1
				nx <<= 1
				ny <<= 1
				if nfx <= c {
					nx--
				} else {
					nfx -= c
				}
				if nfy <= c {
					ny--
				} else {
					nfy -= c
				}
				qtreeExpand(r, scratch, nx, ny, scratch)
			}
			qtreeBitins(scratch, nqx, nqy, a[base:], n, bit)
		}
	}
	return nil
}

// qtreeExpand does one quadtree expansion step on a[(nx+1)/2, (ny+1)/2],
// writing the results into b[nx, ny] (which may be the same array).
func qtreeExpand(r *hcompressReader, a []byte, nx, ny int, b []byte) {
	// First copy a to b, expanding each 4-bit value into a 2x2 block.
	qtreeCopy(a, nx, ny, b, ny)
	// Then read new 4-bit values for each non-zero element.
	for i := nx*ny - 1; i >= 0; i-- {
		if b[i] != 0 {
			b[i] = byte(r.inputHuffman())
		}
	}
}

// qtreeCopy expands the (nx2, ny2) block of 4-bit values into a
// (nx, ny) bit image with each value expanded to its 2x2 constituent.
// a and b may alias — the copy is done end-first.
func qtreeCopy(a []byte, nx, ny int, b []byte, n int) {
	nx2 := (nx + 1) / 2
	ny2 := (ny + 1) / 2

	// Copy 4-bit values to b end-first so a and b may alias.
	k := ny2*(nx2-1) + ny2 - 1
	for i := nx2 - 1; i >= 0; i-- {
		s00 := 2 * (n*i + ny2 - 1)
		for j := ny2 - 1; j >= 0; j-- {
			b[s00] = a[k]
			k--
			s00 -= 2
		}
	}
	// Expand each 2x2 block. Bit layout:
	//	b[s00  ] = (a[k] >> 3) & 1  (top-left)
	//	b[s00+1] = (a[k] >> 2) & 1  (top-right)
	//	b[s10  ] = (a[k] >> 1) & 1  (bottom-left)
	//	b[s10+1] = (a[k]     ) & 1  (bottom-right)
	i := 0
	for ; i < nx-1; i += 2 {
		s00 := n * i
		s10 := s00 + n
		j := 0
		for ; j < ny-1; j += 2 {
			v := b[s00]
			b[s10+1] = v & 1
			b[s10] = (v >> 1) & 1
			b[s00+1] = (v >> 2) & 1
			b[s00] = (v >> 3) & 1
			s00 += 2
			s10 += 2
		}
		if j < ny {
			// Odd row length — only top-left and bottom-left.
			v := b[s00]
			b[s10] = (v >> 1) & 1
			b[s00] = (v >> 3) & 1
		}
	}
	if i < nx {
		// Odd column length — only top row.
		s00 := n * i
		j := 0
		for ; j < ny-1; j += 2 {
			v := b[s00]
			b[s00+1] = (v >> 2) & 1
			b[s00] = (v >> 3) & 1
			s00 += 2
		}
		if j < ny {
			v := b[s00]
			b[s00] = (v >> 3) & 1
		}
	}
}

// qtreeBitins inserts 4-bit values from a[(nx+1)/2, (ny+1)/2] into
// bitplane `bit` of b[nx, ny]. Each 4-bit value expands to 2x2 pixels.
// a and b MUST NOT alias.
func qtreeBitins(a []byte, nx, ny int, b []int32, n, bit int) {
	planeVal := int32(1) << uint(bit)
	k := 0
	i := 0
	for ; i < nx-1; i += 2 {
		s00 := n * i
		j := 0
		for ; j < ny-1; j += 2 {
			v := a[k]
			// 4 bits → 2x2 block at (s00, s00+1, s00+n, s00+n+1).
			if v&1 != 0 {
				b[s00+n+1] |= planeVal
			}
			if v&2 != 0 {
				b[s00+n] |= planeVal
			}
			if v&4 != 0 {
				b[s00+1] |= planeVal
			}
			if v&8 != 0 {
				b[s00] |= planeVal
			}
			s00 += 2
			k++
		}
		if j < ny {
			v := a[k]
			if v&2 != 0 {
				b[s00+n] |= planeVal
			}
			if v&8 != 0 {
				b[s00] |= planeVal
			}
			k++
		}
	}
	if i < nx {
		s00 := n * i
		j := 0
		for ; j < ny-1; j += 2 {
			v := a[k]
			if v&4 != 0 {
				b[s00+1] |= planeVal
			}
			if v&8 != 0 {
				b[s00] |= planeVal
			}
			s00 += 2
			k++
		}
		if j < ny {
			v := a[k]
			if v&8 != 0 {
				b[s00] |= planeVal
			}
			k++
		}
	}
}

// hinvHCOMPRESS applies the inverse H-transform to a[nx*ny]. This is
// a 2D butterfly operation performed log2n times, where each pass
// un-shuffles and un-combines the wavelet coefficients to reconstruct
// successively larger pixel blocks.
func hinvHCOMPRESS(a []int32, nx, ny, smooth, scale int) error {
	nmax := nx
	if ny > nmax {
		nmax = ny
	}
	log2n := 0
	for 1<<uint(log2n) < nmax {
		log2n++
	}

	// Temporary unshuffle buffer.
	tmp := make([]int32, (nmax+1)/2)

	shift := 1
	bit0 := int32(1) << uint(log2n-1)
	bit1 := bit0 << 1
	bit2 := bit0 << 2
	mask0 := -bit0
	mask1 := mask0 << 1
	mask2 := mask0 << 2
	prnd0 := bit0 >> 1
	prnd1 := bit1 >> 1
	prnd2 := bit2 >> 1
	nrnd0 := prnd0 - 1
	nrnd1 := prnd1 - 1
	nrnd2 := prnd2 - 1

	// Round a[0] to multiple of bit2.
	if a[0] >= 0 {
		a[0] = (a[0] + prnd2) & mask2
	} else {
		a[0] = (a[0] + nrnd2) & mask2
	}

	nxtop := 1
	nytop := 1
	nxf := nx
	nyf := ny
	c := 1 << uint(log2n)

	for k := log2n - 1; k >= 0; k-- {
		c >>= 1
		nxtop <<= 1
		nytop <<= 1
		if nxf <= c {
			nxtop--
		} else {
			nxf -= c
		}
		if nyf <= c {
			nytop--
		} else {
			nyf -= c
		}
		// Double shift and fix nrnd0 on last pass (prnd0 = 0 there).
		if k == 0 {
			nrnd0 = 0
			shift = 2
		}
		// Unshuffle in each dimension to interleave coefficients.
		for i := 0; i < nxtop; i++ {
			unshuffleInts(a[ny*i:], nytop, 1, tmp)
		}
		for j := 0; j < nytop; j++ {
			unshuffleInts(a[j:], nxtop, ny, tmp)
		}
		// (smoothing omitted — we use the standard no-smooth path,
		// matching cfitsio behavior when SMOOTH=0 in the header.)
		oddx := nxtop % 2
		oddy := nytop % 2

		var i int
		for i = 0; i < nxtop-oddx; i += 2 {
			s00 := ny * i
			s10 := s00 + ny
			for j := 0; j < nytop-oddy; j += 2 {
				h0 := a[s00]
				hx := a[s10]
				hy := a[s00+1]
				hc := a[s10+1]
				// Round hx, hy to multiple of bit1; hc to multiple of bit0.
				if hx >= 0 {
					hx = (hx + prnd1) & mask1
				} else {
					hx = (hx + nrnd1) & mask1
				}
				if hy >= 0 {
					hy = (hy + prnd1) & mask1
				} else {
					hy = (hy + nrnd1) & mask1
				}
				if hc >= 0 {
					hc = (hc + prnd0) & mask0
				} else {
					hc = (hc + nrnd0) & mask0
				}
				// Propagate bit0 of hc to hx, hy.
				lowbit0 := hc & bit0
				if hx >= 0 {
					hx -= lowbit0
				} else {
					hx += lowbit0
				}
				if hy >= 0 {
					hy -= lowbit0
				} else {
					hy += lowbit0
				}
				// Propagate bits 0 and 1 of hc, hx, hy to h0.
				lowbit1 := (hc ^ hx ^ hy) & bit1
				if h0 >= 0 {
					h0 = h0 + lowbit0 - lowbit1
				} else {
					if lowbit0 == 0 {
						h0 = h0 + lowbit1
					} else {
						h0 = h0 + (lowbit0 - lowbit1)
					}
				}
				// Divide sums by 2 (4 on last pass).
				a[s10+1] = (h0 + hx + hy + hc) >> uint(shift)
				a[s10] = (h0 + hx - hy - hc) >> uint(shift)
				a[s00+1] = (h0 - hx + hy - hc) >> uint(shift)
				a[s00] = (h0 - hx - hy + hc) >> uint(shift)
				s00 += 2
				s10 += 2
			}
			if oddy != 0 {
				// Last element in row — s00+1, s10+1 are off the edge.
				h0 := a[s00]
				hx := a[s10]
				if hx >= 0 {
					hx = (hx + prnd1) & mask1
				} else {
					hx = (hx + nrnd1) & mask1
				}
				lowbit1 := hx & bit1
				if h0 >= 0 {
					h0 -= lowbit1
				} else {
					h0 += lowbit1
				}
				a[s10] = (h0 + hx) >> uint(shift)
				a[s00] = (h0 - hx) >> uint(shift)
			}
		}
		if oddx != 0 {
			// Last row — s10 off edge.
			s00 := ny * i
			for j := 0; j < nytop-oddy; j += 2 {
				h0 := a[s00]
				hy := a[s00+1]
				if hy >= 0 {
					hy = (hy + prnd1) & mask1
				} else {
					hy = (hy + nrnd1) & mask1
				}
				lowbit1 := hy & bit1
				if h0 >= 0 {
					h0 -= lowbit1
				} else {
					h0 += lowbit1
				}
				a[s00+1] = (h0 + hy) >> uint(shift)
				a[s00] = (h0 - hy) >> uint(shift)
				s00 += 2
			}
			if oddy != 0 {
				// Corner element.
				a[s00] = a[s00] >> uint(shift)
			}
		}
		// Divide all masks and rounding values by 2.
		bit2 = bit1
		bit1 = bit0
		bit0 >>= 1
		mask1 = mask0
		mask0 >>= 1
		prnd1 = prnd0
		prnd0 >>= 1
		nrnd1 = nrnd0
		nrnd0 = prnd0 - 1
	}
	return nil
}

// unshuffleInts is the cfitsio unshuffle helper: moves the first half of
// an n-length stride-n2 array into the even positions, the second half
// into the odd positions. Used by the inverse H-transform to interleave
// wavelet coefficients back into pixel order.
func unshuffleInts(a []int32, n, n2 int, tmp []int32) {
	nhalf := (n + 1) >> 1
	// Copy second half of array to tmp.
	pt := 0
	p1 := n2 * nhalf
	for i := nhalf; i < n; i++ {
		tmp[pt] = a[p1]
		p1 += n2
		pt++
	}
	// Distribute first half of array to even elements.
	p2 := n2 * (nhalf - 1)
	p1 = (n2 * (nhalf - 1)) << 1
	for i := nhalf - 1; i >= 0; i-- {
		a[p1] = a[p2]
		p2 -= n2
		p1 -= (n2 + n2)
	}
	// Distribute second half (in tmp) to odd elements.
	pt = 0
	p1 = n2
	for i := 1; i < n; i += 2 {
		a[p1] = tmp[pt]
		p1 += (n2 + n2)
		pt++
	}
}
