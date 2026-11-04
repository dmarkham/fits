package compress

import (
	"encoding/binary"
	"fmt"
)

// PLIO_1 encoder — IRAF PLIO line-list format, symmetric counterpart of
// the decoder in plio.go. A line-for-line port of cfitsio's pliocomp.c
// pl_p2li (itself f2c output of the IRAF SPP routine), so the token
// stream is byte-identical to what cfitsio writes
// (TestPLIOEncodeMatchesCfitsio, golden line lists from
// compress/testdata/cref).
//
// Like cfitsio (imcompress.c imcomp_compress_tile), the encoder refuses
// any pixel outside [0, 2^24-1] with ErrPLIORange: pl_p2li silently
// clamps negatives to 0 and cannot represent larger values, so writing
// them would store different data.
//
// Input pixels: BYTEPIX 1 is the unsigned on-disk byte, 2 and 4 are
// signed big-endian integers (the caller converts unsigned types the way
// cfitsio does before calling).

type plio1Encoder struct{}

// plioMaxValue is the largest pixel value PLIO_1 can encode.
const plioMaxValue = 16777215

func (plio1Encoder) Encode(src, dst []byte, nelem, elemSize int) (int, error) {
	if len(src) < nelem*elemSize {
		return 0, fmt.Errorf("%w: PLIO_1 src length %d < %d", ErrCorrupt, len(src), nelem*elemSize)
	}
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
			pixels[i] = int32(binary.BigEndian.Uint32(src[i*4:]))
		}
	default:
		return 0, fmt.Errorf("%w: PLIO_1 unsupported elemSize %d", ErrUnsupportedCompression, elemSize)
	}
	for i, v := range pixels {
		if v < 0 || v > plioMaxValue {
			return 0, fmt.Errorf("%w: pixel %d = %d", ErrPLIORange, i, v)
		}
	}
	tokens := plP2LI(pixels)
	need := len(tokens) * 2
	if need > len(dst) {
		return 0, ErrBufferTooSmall
	}
	for k, t := range tokens {
		binary.BigEndian.PutUint16(dst[k*2:], uint16(t))
	}
	return need, nil
}

// plP2LI is cfitsio's pl_p2li(pxsrc, xs=1, lldst, npix): it converts a
// pixel array to a PLIO line list and returns the list (length = the C
// return value). The goto structure of the f2c original is kept as
// labelled control flow so the port can be checked against it line by
// line; ll and px are 1-based like the C after its "--lldst; --pxsrc".
func plP2LI(pixels []int32) []int16 {
	npix := len(pixels)
	if npix <= 0 {
		return nil
	}
	px := func(i int) int { return int(pixels[i-1]) } // pxsrc[i], 1-based
	// Worst case is a few words per pixel; grow on demand.
	ll := make([]int16, 8, 8+npix)
	set := func(i int, v int) {
		for len(ll) <= i {
			ll = append(ll, 0)
		}
		ll[i] = int16(v)
	}
	get := func(i int) int { return int(ll[i]) }
	maxi := func(a, b int) int {
		if a > b {
			return a
		}
		return b
	}
	mini := func(a, b int) int {
		if a < b {
			return a
		}
		return b
	}
	abs := func(a int) int {
		if a < 0 {
			return -a
		}
		return a
	}

	const xs = 1
	set(3, -100)
	set(2, 7)
	set(1, 0)
	set(6, 0)
	set(7, 0)
	xe := xs + npix - 1
	op := 8
	pv := maxi(0, px(xs))
	x1 := xs
	iz := xs
	hi := 1
	nv := 0
	for ip := xs; ip <= xe; ip++ {
		if ip < xe {
			nv = maxi(0, px(ip+1))
			if nv == pv {
				continue // L120
			}
			if pv == 0 {
				pv = nv
				x1 = ip + 1
				continue // L120
			}
			// goto L131
		} else if pv == 0 { // L130
			x1 = xe + 1
		}
		// L131
		np := ip - x1 + 1
		nz := x1 - iz
		emitted := false // true when the code reaches L91 directly
		if pv > 0 {
			dv := pv - hi
			if dv != 0 {
				hi = pv
				if abs(dv) > 4095 {
					set(op, (pv&4095)+4096)
					op++
					set(op, pv/4096)
					op++
				} else {
					if dv < 0 {
						set(op, -dv+12288)
					} else {
						set(op, dv+8192)
					}
					op++
					if np == 1 && nz == 0 {
						v := get(op - 1)
						set(op-1, v|16384)
						emitted = true // goto L91
					}
				}
			}
		}
		if !emitted {
			if nz > 0 {
				for nz > 0 {
					set(op, mini(4095, nz))
					op++
					nz -= 4095
				}
				if np == 1 && pv > 0 {
					set(op-1, get(op-1)+20481)
					emitted = true // goto L91
				}
			}
		}
		if !emitted {
			for np > 0 {
				set(op, mini(4095, np)+16384)
				op++
				np -= 4095
			}
		}
		// L91
		x1 = ip + 1
		iz = x1
		pv = nv
	}
	set(4, (op-1)%32768)
	set(5, (op-1)/32768)
	// The C list is lldst[1..op-1]; return it 0-based.
	return append([]int16(nil), ll[1:op]...)
}

// ErrPLIORange is returned when a PLIO_1 tile holds a pixel outside
// [0, 2^24-1], the range cfitsio's PLIO codec accepts.
var ErrPLIORange = fmt.Errorf("%w: data out of range for PLIO compression (0 - 2**24)", ErrUnsupportedCompression)
