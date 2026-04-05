package compress

import (
	"encoding/binary"
	"testing"
)

// TestRiceKnown100 decodes the first row (100 pixels) of the real
// astropy-generated comp_rice_i16.fits fixture. 4 blocks of 32/32/32/4
// pixels exercise the multi-block path and the partial final block.
func TestRiceKnown100(t *testing.T) {
	hexStr := "000019249249249249249249249244924924924924924924924924492492492492492492492492449240"
	compressed := make([]byte, len(hexStr)/2)
	for i := 0; i < len(compressed); i++ {
		var b byte
		for j := 0; j < 2; j++ {
			c := hexStr[i*2+j]
			var v byte
			if c >= '0' && c <= '9' {
				v = c - '0'
			} else {
				v = c - 'a' + 10
			}
			b = (b << 4) | v
		}
		compressed[i] = b
	}
	params := Params{"BLOCKSIZE": 32, "BYTEPIX": 2}
	dec := newRICE1(params)
	const nelem = 100
	const elemSize = 2
	dst := make([]byte, nelem*elemSize)
	if err := dec.Decode(compressed, dst, nelem, elemSize); err != nil {
		t.Fatalf("decode: %v", err)
	}
	for i := 0; i < nelem; i++ {
		got := int16((int16(dst[i*2])<<8) | int16(dst[i*2+1]))
		if got != int16(i) {
			t.Errorf("pixel %d: got %d want %d", i, got, i)
			if i > 5 {
				break
			}
		}
	}
}

// TestRiceKnownTiny decodes a tiny RICE_1 tile whose compressed bytes
// and expected output were captured from astropy.
//
// Input data: int16 array [10, 11, 12, 13, 14, 15, 16, 17] (8 pixels).
// Compressed bytes (6 bytes, captured from astropy): 00 0a 19 24 92 40
// Expected decoded: same 8 pixels.
func TestRiceKnownTiny(t *testing.T) {
	compressed := []byte{0x00, 0x0a, 0x19, 0x24, 0x92, 0x40}
	params := Params{"BLOCKSIZE": 32, "BYTEPIX": 2}
	dec := newRICE1(params)

	const nelem = 8
	const elemSize = 2
	dst := make([]byte, nelem*elemSize)
	if err := dec.Decode(compressed, dst, nelem, elemSize); err != nil {
		t.Fatalf("decode: %v", err)
	}
	got := make([]int16, nelem)
	for i := 0; i < nelem; i++ {
		got[i] = int16(binary.BigEndian.Uint16(dst[i*2:]))
	}
	want := []int16{10, 11, 12, 13, 14, 15, 16, 17}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("pixel %d: got %d want %d (full got=%v)", i, got[i], want[i], got)
		}
	}
}
