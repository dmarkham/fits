package compress

import (
	"encoding/binary"
	"fmt"
)

// HCOMPRESS_1 encoder. Symmetric counterpart of the decoder in
// hcompress.go, ported from cfitsio/fits_hcompress.c.
//
// Pipeline:
//
//  1. htrans   — forward H-transform (log2n passes of 2x2 butterfly +
//                row/column shuffle). Produces signed wavelet coefficients.
//  2. digitize — divide by scale (0 = lossless, leaves coefficients
//                unchanged).
//  3. encode   — write header (magic + nx + ny + scale + sumall), pull
//                sign bits into a separate buffer, compute per-quadrant
//                nbitplanes[3], then call doencode.
//  4. doencode — for each quadrant, call qtreeEncodeHC which walks the
//                bit planes from MSB to LSB. For each plane, try
//                quadtree-Huffman coding; if that would exceed the
//                direct-bit-map size, write the raw bit map instead.
//                Append EOF nybble + sign bits.
//
// HCOMPRESS is 2D — the caller must set nx and ny via SetShape before
// calling Encode. nelem must equal nx*ny.
type hcompress1Encoder struct {
	scale  int
	smooth int
	nx, ny int // 2D shape; set by SetShape before Encode
}

// SetShape sets the 2D tile dimensions. HCOMPRESS is the only encoder
// that needs shape information (beyond nelem) because its transform is
// 2D. Callers use a type-assertion to get at this method.
func (h *hcompress1Encoder) SetShape(nx, ny int) {
	h.nx = nx
	h.ny = ny
}

func newHCOMPRESS1Encoder(params Params) Encoder {
	return &hcompress1Encoder{
		scale:  int(params.Get("SCALE", 0)),
		smooth: int(params.Get("SMOOTH", 0)),
	}
}

func (h *hcompress1Encoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	work := make([]int32, nelem)
	for i := 0; i < nelem; i++ {
		switch elemSize {
		case 1:
			work[i] = int32(int8(src[i]))
		case 2:
			work[i] = int32(int16(uint16(src[i*2])<<8 | uint16(src[i*2+1])))
		case 4:
			work[i] = int32(binary.BigEndian.Uint32(src[i*4:]))
		case 8:
			work[i] = int32(binary.BigEndian.Uint64(src[i*8:]))
		default:
			return 0, fmt.Errorf("%w: HCOMPRESS_1 unsupported elemSize %d", ErrCorrupt, elemSize)
		}
	}
	nx := h.nx
	ny := h.ny
	if nx == 0 || ny == 0 {
		// Fall back to treating it as a 1D strip (nx=1, ny=nelem) if the
		// caller did not set a 2D shape. The H-transform still works on
		// degenerate dimensions.
		nx = 1
		ny = nelem
	}
	if nx*ny != nelem {
		return 0, fmt.Errorf("%w: HCOMPRESS_1 nx*ny=%d != nelem=%d", ErrCorrupt, nx*ny, nelem)
	}
	if err := htransHC(work, nx, ny); err != nil {
		return 0, err
	}
	digitizeHC(work, nx, ny, h.scale)
	w := &hcompressWriter{dst: dst}
	if err := encodeHC(w, work, nx, ny, h.scale); err != nil {
		return 0, err
	}
	if w.overflow {
		return 0, ErrBufferTooSmall
	}
	return w.byteIdx, nil
}

// ------------------ htrans (forward H-transform) ------------------

// htransHC is the forward H-transform: log2n passes over the image,
// each pass reducing 2x2 blocks to 4 wavelet coefficients (h0 = sum,
// hx, hy = partial differences, hc = cross-term) followed by a row/
// column shuffle to group coefficients by order.
//
// In-place on a[nx*ny]. Port of cfitsio's htrans().
func htransHC(a []int32, nx, ny int) error {
	nmax := nx
	if ny > nmax {
		nmax = ny
	}
	log2n := 0
	for 1<<uint(log2n) < nmax {
		log2n++
	}
	tmp := make([]int32, (nmax+1)/2)

	shift := 0
	mask := int32(-2)
	mask2 := mask << 1
	prnd := int32(1)
	prnd2 := prnd << 1
	nrnd2 := prnd2 - 1

	nxtop := nx
	nytop := ny
	for k := 0; k < log2n; k++ {
		oddx := nxtop % 2
		oddy := nytop % 2
		var i int
		for i = 0; i < nxtop-oddx; i += 2 {
			s00 := i * ny
			s10 := s00 + ny
			for j := 0; j < nytop-oddy; j += 2 {
				h0 := (a[s10+1] + a[s10] + a[s00+1] + a[s00]) >> uint(shift)
				hx := (a[s10+1] + a[s10] - a[s00+1] - a[s00]) >> uint(shift)
				hy := (a[s10+1] - a[s10] + a[s00+1] - a[s00]) >> uint(shift)
				hc := (a[s10+1] - a[s10] - a[s00+1] + a[s00]) >> uint(shift)

				a[s10+1] = hc
				if hx >= 0 {
					a[s10] = (hx + prnd) & mask
				} else {
					a[s10] = hx & mask
				}
				if hy >= 0 {
					a[s00+1] = (hy + prnd) & mask
				} else {
					a[s00+1] = hy & mask
				}
				if h0 >= 0 {
					a[s00] = (h0 + prnd2) & mask2
				} else {
					a[s00] = (h0 + nrnd2) & mask2
				}
				s00 += 2
				s10 += 2
			}
			if oddy != 0 {
				h0 := (a[s10] + a[s00]) << uint(1-shift)
				hx := (a[s10] - a[s00]) << uint(1-shift)
				if hx >= 0 {
					a[s10] = (hx + prnd) & mask
				} else {
					a[s10] = hx & mask
				}
				if h0 >= 0 {
					a[s00] = (h0 + prnd2) & mask2
				} else {
					a[s00] = (h0 + nrnd2) & mask2
				}
			}
		}
		if oddx != 0 {
			s00 := i * ny
			for j := 0; j < nytop-oddy; j += 2 {
				h0 := (a[s00+1] + a[s00]) << uint(1-shift)
				hy := (a[s00+1] - a[s00]) << uint(1-shift)
				if hy >= 0 {
					a[s00+1] = (hy + prnd) & mask
				} else {
					a[s00+1] = hy & mask
				}
				if h0 >= 0 {
					a[s00] = (h0 + prnd2) & mask2
				} else {
					a[s00] = (h0 + nrnd2) & mask2
				}
				s00 += 2
			}
			if oddy != 0 {
				h0 := a[s00] << uint(2-shift)
				if h0 >= 0 {
					a[s00] = (h0 + prnd2) & mask2
				} else {
					a[s00] = (h0 + nrnd2) & mask2
				}
			}
		}
		// Row-wise and column-wise shuffle.
		for ii := 0; ii < nxtop; ii++ {
			shuffleInts(a[ny*ii:], nytop, 1, tmp)
		}
		for jj := 0; jj < nytop; jj++ {
			shuffleInts(a[jj:], nxtop, ny, tmp)
		}
		nxtop = (nxtop + 1) >> 1
		nytop = (nytop + 1) >> 1
		shift = 1
		mask = mask2
		prnd = prnd2
		mask2 = mask2 << 1
		prnd2 = prnd2 << 1
		nrnd2 = prnd2 - 1
	}
	return nil
}

// shuffleInts is the inverse of unshuffleInts (used by the decoder).
// It groups even-indexed elements in the first half and odd-indexed in
// the second half — moving interleaved coefficients back into sorted
// order.
func shuffleInts(a []int32, n, n2 int, tmp []int32) {
	nhalf := (n + 1) >> 1
	// Copy odd (second-half destined) elements to tmp.
	pt := 0
	p1 := n2
	for i := 1; i < n; i += 2 {
		tmp[pt] = a[p1]
		p1 += n2 + n2
		pt++
	}
	// Collapse even elements to the front.
	p2 := 0
	p1 = 0
	for i := 0; i < nhalf; i++ {
		a[p1] = a[p2]
		p1 += n2
		p2 += n2 + n2
	}
	// Append tmp (the odd elements) to the second half.
	pt = 0
	p1 = n2 * nhalf
	for i := nhalf; i < n; i++ {
		a[p1] = tmp[pt]
		p1 += n2
		pt++
	}
}

// ------------------ digitize ------------------

// digitizeHC divides every element by scale, rounding toward zero with
// a half-step bias so that both positive and negative values round
// symmetrically. scale <= 1 is a no-op (lossless HCOMPRESS).
func digitizeHC(a []int32, nx, ny, scale int) {
	if scale <= 1 {
		return
	}
	d := int32(scale / 2)
	s := int32(scale)
	for i := 0; i < nx*ny; i++ {
		if a[i] > 0 {
			a[i] = (a[i] + d) / s
		} else {
			a[i] = (a[i] - d) / s
		}
	}
}

// ------------------ encode ------------------

// encodeHC writes the HCOMPRESS header + quadtree body + sign bits to
// the output.
func encodeHC(w *hcompressWriter, a []int32, nx, ny, scale int) error {
	// Magic, nx, ny, scale.
	w.writeBytes(hcompressMagic[:])
	w.writeInt32BE(int32(nx))
	w.writeInt32BE(int32(ny))
	w.writeInt32BE(int32(scale))
	// sumall (8 bytes) — the DC level, stored separately.
	w.writeInt64BE(int64(a[0]))
	a[0] = 0

	// Build sign-bit buffer: one bit per non-zero pixel.
	nel := nx * ny
	signbits := make([]byte, (nel+7)/8)
	nsign := 0
	bitsToGo := 8
	for i := 0; i < nel; i++ {
		switch {
		case a[i] > 0:
			signbits[nsign] <<= 1
			bitsToGo--
		case a[i] < 0:
			signbits[nsign] <<= 1
			signbits[nsign] |= 1
			bitsToGo--
			a[i] = -a[i]
		}
		if bitsToGo == 0 {
			bitsToGo = 8
			nsign++
		}
	}
	if bitsToGo != 8 {
		signbits[nsign] <<= uint(bitsToGo)
		nsign++
	}

	// Per-quadrant max → nbitplanes[3].
	var vmax [3]int32
	nx2 := (nx + 1) / 2
	ny2 := (ny + 1) / 2
	jj := 0
	kk := 0
	for i := 0; i < nel; i++ {
		q := 0
		if jj >= ny2 {
			q++
		}
		if kk >= nx2 {
			q++
		}
		if vmax[q] < a[i] {
			vmax[q] = a[i]
		}
		jj++
		if jj >= ny {
			jj = 0
			kk++
		}
	}
	var nbitplanes [3]byte
	for q := 0; q < 3; q++ {
		for nbitplanes[q] = 0; vmax[q] > 0; vmax[q] >>= 1 {
			nbitplanes[q]++
		}
	}
	w.writeBytes(nbitplanes[:])

	// doencode: quadtree bit-plane coding for 4 quadrants.
	w.startBits()
	if err := qtreeEncodeHC(w, a, 0, ny, nx2, ny2, int(nbitplanes[0])); err != nil {
		return err
	}
	if err := qtreeEncodeHC(w, a, ny2, ny, nx2, ny/2, int(nbitplanes[1])); err != nil {
		return err
	}
	if err := qtreeEncodeHC(w, a, ny*nx2, ny, nx/2, ny2, int(nbitplanes[1])); err != nil {
		return err
	}
	if err := qtreeEncodeHC(w, a, ny*nx2+ny2, ny, nx/2, ny/2, int(nbitplanes[2])); err != nil {
		return err
	}
	// EOF nybble.
	w.outputNybble(0)
	w.flushBits()
	// Sign bits (as raw bytes, not through the bit writer).
	if nsign > 0 {
		w.writeBytes(signbits[:nsign])
	}
	return nil
}

// ------------------ qtree_encode ------------------

// Fixed Huffman code table (same values the decoder reads via
// inputHuffman), used to encode 4-bit quadtree nybbles.
var hcHuffCode = [16]uint16{
	0x3e, 0x00, 0x01, 0x08, 0x02, 0x09, 0x1a, 0x1b,
	0x03, 0x1c, 0x0a, 0x1d, 0x0b, 0x1e, 0x3f, 0x0c,
}
var hcHuffNBits = [16]int{
	6, 3, 3, 4, 3, 4, 5, 5,
	3, 5, 4, 5, 4, 5, 6, 4,
}

// qtreeEncodeHC encodes the coefficients of one quadrant, bit plane by
// bit plane from MSB down.
func qtreeEncodeHC(w *hcompressWriter, a []int32, base, n, nqx, nqy, nbitplanes int) error {
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
	bmax := (nqx2*nqy2 + 1) / 2
	scratch := make([]byte, 2*bmax)
	buffer := make([]byte, bmax)

	for bit := nbitplanes - 1; bit >= 0; bit-- {
		b := 0
		bitbuffer := 0
		bitsToGo3 := 0

		// First pass: copy bit plane of A to scratch.
		qtreeOnebitHC(a[base:], n, nqx, nqy, scratch, bit)
		nx := (nqx + 1) >> 1
		ny := (nqy + 1) >> 1

		expanded := bufcopyHC(scratch, nx*ny, buffer, &b, bmax, &bitbuffer, &bitsToGo3)
		if expanded {
			writeBdirectHC(w, a[base:], n, nqx, nqy, scratch, bit)
			continue
		}
		// log2n-1 reductions.
		done := false
		for k := 1; k < log2n; k++ {
			qtreeReduceHC(scratch, ny, nx, ny, scratch)
			nx = (nx + 1) >> 1
			ny = (ny + 1) >> 1
			expanded = bufcopyHC(scratch, nx*ny, buffer, &b, bmax, &bitbuffer, &bitsToGo3)
			if expanded {
				writeBdirectHC(w, a[base:], n, nqx, nqy, scratch, bit)
				done = true
				break
			}
		}
		if done {
			continue
		}
		// Write warning code 0xF, then the remaining bit buffer contents,
		// then the buffer in reverse order.
		w.outputNybble(0xF)
		if b == 0 {
			if bitsToGo3 > 0 {
				w.outputNbits(uint(bitbuffer&((1<<uint(bitsToGo3))-1)), bitsToGo3)
			} else {
				// No 1s in the plane at all — write a zero Huffman code.
				w.outputNbits(uint(hcHuffCode[0]), hcHuffNBits[0])
			}
		} else {
			if bitsToGo3 > 0 {
				w.outputNbits(uint(bitbuffer&((1<<uint(bitsToGo3))-1)), bitsToGo3)
			}
			for i := b - 1; i >= 0; i-- {
				w.outputNbits(uint(buffer[i]), 8)
			}
		}
	}
	return nil
}

// qtreeOnebitHC packs bit `bit` of A's coefficient quadrant into
// scratch, expanding each 2x2 block of A to a 4-bit value in one byte
// of scratch.
func qtreeOnebitHC(a []int32, n, nx, ny int, b []byte, bit int) {
	b0 := int32(1) << uint(bit)
	b1 := b0 << 1
	b2 := b0 << 2
	b3 := b0 << 3
	k := 0
	var i int
	for i = 0; i < nx-1; i += 2 {
		s00 := n * i
		s10 := s00 + n
		var j int
		for j = 0; j < ny-1; j += 2 {
			b[k] = byte(((a[s10+1] & b0) |
				((a[s10] << 1) & b1) |
				((a[s00+1] << 2) & b2) |
				((a[s00] << 3) & b3)) >> uint(bit))
			k++
			s00 += 2
			s10 += 2
		}
		if j < ny {
			b[k] = byte((((a[s10] << 1) & b1) |
				((a[s00] << 3) & b3)) >> uint(bit))
			k++
		}
	}
	if i < nx {
		s00 := n * i
		var j int
		for j = 0; j < ny-1; j += 2 {
			b[k] = byte((((a[s00+1] << 2) & b2) |
				((a[s00] << 3) & b3)) >> uint(bit))
			k++
			s00 += 2
		}
		if j < ny {
			b[k] = byte(((a[s00] << 3) & b3) >> uint(bit))
			k++
		}
	}
}

// qtreeReduceHC reduces one quadtree level: 2x2 blocks of scratch → one
// 4-bit value (non-zero flag per 2x2 sub-block).
func qtreeReduceHC(a []byte, n, nx, ny int, b []byte) {
	k := 0
	var i int
	for i = 0; i < nx-1; i += 2 {
		s00 := n * i
		s10 := s00 + n
		var j int
		for j = 0; j < ny-1; j += 2 {
			var v byte
			if a[s10+1] != 0 {
				v |= 1
			}
			if a[s10] != 0 {
				v |= 2
			}
			if a[s00+1] != 0 {
				v |= 4
			}
			if a[s00] != 0 {
				v |= 8
			}
			b[k] = v
			k++
			s00 += 2
			s10 += 2
		}
		if j < ny {
			var v byte
			if a[s10] != 0 {
				v |= 2
			}
			if a[s00] != 0 {
				v |= 8
			}
			b[k] = v
			k++
		}
	}
	if i < nx {
		s00 := n * i
		var j int
		for j = 0; j < ny-1; j += 2 {
			var v byte
			if a[s00+1] != 0 {
				v |= 4
			}
			if a[s00] != 0 {
				v |= 8
			}
			b[k] = v
			k++
			s00 += 2
		}
		if j < ny {
			var v byte
			if a[s00] != 0 {
				v |= 8
			}
			b[k] = v
			k++
		}
	}
}

// bufcopyHC copies non-zero 4-bit codes from a[] to buffer[], encoded
// via the Huffman table. Returns true if the buffer would overflow (the
// caller then falls back to the direct bit-map output path).
func bufcopyHC(a []byte, n int, buffer []byte, b *int, bmax int, bitbuffer, bitsToGo3 *int) bool {
	for i := 0; i < n; i++ {
		if a[i] != 0 {
			*bitbuffer |= int(hcHuffCode[a[i]]) << uint(*bitsToGo3)
			*bitsToGo3 += hcHuffNBits[a[i]]
			if *bitsToGo3 >= 8 {
				buffer[*b] = byte(*bitbuffer & 0xff)
				*b++
				if *b >= bmax {
					return true
				}
				*bitbuffer >>= 8
				*bitsToGo3 -= 8
			}
		}
	}
	return false
}

// writeBdirectHC writes the direct-bit-map representation of one bit
// plane: 0 nybble followed by one 4-bit nybble per 2x2 block.
func writeBdirectHC(w *hcompressWriter, a []int32, n, nqx, nqy int, scratch []byte, bit int) {
	w.outputNybble(0)
	qtreeOnebitHC(a, n, nqx, nqy, scratch, bit)
	nNybbles := ((nqx + 1) / 2) * ((nqy + 1) / 2)
	w.outputNnybble(nNybbles, scratch)
}

// ------------------ bit writer ------------------

// hcompressWriter is the HCOMPRESS bit writer. MSB-first within each
// byte, matching cfitsio's output_nbits / output_nybble / output_nnybble
// and the hcompressReader decoder.
type hcompressWriter struct {
	dst      []byte
	byteIdx  int
	bitBuf   int // pending bits, left-justified
	bitsToGo int // bits free in the top byte (starts at 8)
	overflow bool
}

func (w *hcompressWriter) writeBytes(b []byte) {
	if w.byteIdx+len(b) > len(w.dst) {
		w.overflow = true
		return
	}
	copy(w.dst[w.byteIdx:], b)
	w.byteIdx += len(b)
}

func (w *hcompressWriter) writeInt32BE(v int32) {
	var b [4]byte
	binary.BigEndian.PutUint32(b[:], uint32(v))
	w.writeBytes(b[:])
}

func (w *hcompressWriter) writeInt64BE(v int64) {
	var b [8]byte
	binary.BigEndian.PutUint64(b[:], uint64(v))
	w.writeBytes(b[:])
}

func (w *hcompressWriter) startBits() {
	w.bitBuf = 0
	w.bitsToGo = 8
}

func (w *hcompressWriter) outputNbits(bits uint, n int) {
	mask := (1 << uint(n)) - 1
	w.bitBuf = (w.bitBuf << uint(n)) | (int(bits) & mask)
	w.bitsToGo -= n
	if w.bitsToGo <= 0 {
		if w.byteIdx >= len(w.dst) {
			w.overflow = true
		} else {
			w.dst[w.byteIdx] = byte((w.bitBuf >> uint(-w.bitsToGo)) & 0xff)
			w.byteIdx++
		}
		w.bitsToGo += 8
	}
}

func (w *hcompressWriter) outputNybble(bits int) {
	w.outputNbits(uint(bits&0xf), 4)
}

func (w *hcompressWriter) outputNnybble(n int, array []byte) {
	for i := 0; i < n; i++ {
		w.outputNybble(int(array[i]))
	}
}

// flushBits writes any bits remaining in the bit buffer, padding with
// zeros on the right. Called at the end of the quadtree body before
// the raw sign-bit bytes.
func (w *hcompressWriter) flushBits() {
	if w.bitsToGo < 8 {
		if w.byteIdx < len(w.dst) {
			w.dst[w.byteIdx] = byte((w.bitBuf << uint(w.bitsToGo)) & 0xff)
			w.byteIdx++
		} else {
			w.overflow = true
		}
		w.bitsToGo = 8
		w.bitBuf = 0
	}
}
