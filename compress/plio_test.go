package compress

import (
	"encoding/binary"
	"testing"
)

// TestPLIO1Tile0 decodes the actual tokens for tile 0 of
// testdata/comp_plio_mask.fits (captured via the fits package's column
// reader). Expected output is 50 zeros.
func TestPLIO1Tile0(t *testing.T) {
	tokens := []int16{0, 7, -100, 8, 0, 0, 0, 50}
	dst := make([]byte, 50*2)
	if err := decodePLIO1(tokens, dst, 50, 2); err != nil {
		t.Fatal(err)
	}
	for i := 0; i < 50; i++ {
		v := int16(binary.BigEndian.Uint16(dst[i*2:]))
		if v != 0 {
			t.Errorf("pixel %d: got %d, want 0 (full dst=%v)", i, v, dst)
			break
		}
	}
}

// TestPLIO1Tile30 decodes the tokens for tile 30: 5 zeros, 40 twos,
// 5 zeros. Mask was set to 2 for rows 30..39, cols 5..44.
func TestPLIO1Tile30(t *testing.T) {
	tokens := []int16{0, 7, -100, 11, 0, 0, 0, 8193, 5, 16424, 5}
	dst := make([]byte, 50*2)
	if err := decodePLIO1(tokens, dst, 50, 2); err != nil {
		t.Fatal(err)
	}
	got := make([]int16, 50)
	for i := 0; i < 50; i++ {
		got[i] = int16(binary.BigEndian.Uint16(dst[i*2:]))
	}
	want := make([]int16, 50)
	for i := 5; i < 45; i++ {
		want[i] = 2
	}
	for i := 0; i < 50; i++ {
		if got[i] != want[i] {
			t.Errorf("pixel %d: got %d want %d (full got=%v)", i, got[i], want[i], got)
			break
		}
	}
}

// TestPLIO1Tile10 decodes tile 10: 15 zeros, 20 ones, 15 zeros.
// Mask was set to 1 for rows 10..19, cols 15..34.
func TestPLIO1Tile10(t *testing.T) {
	tokens := []int16{0, 7, -100, 10, 0, 0, 0, 15, 16404, 15}
	dst := make([]byte, 50*2)
	if err := decodePLIO1(tokens, dst, 50, 2); err != nil {
		t.Fatal(err)
	}
	got := make([]int16, 50)
	for i := 0; i < 50; i++ {
		got[i] = int16(binary.BigEndian.Uint16(dst[i*2:]))
	}
	want := make([]int16, 50)
	for i := 15; i < 35; i++ {
		want[i] = 1
	}
	for i := 0; i < 50; i++ {
		if got[i] != want[i] {
			t.Errorf("pixel %d: got %d want %d (full got=%v)", i, got[i], want[i], got)
			break
		}
	}
}
