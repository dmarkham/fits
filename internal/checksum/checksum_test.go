package checksum

import (
	"bytes"
	"testing"
)

// TestAppendixJExample exercises the round-trip Encode→Decode identity over
// a spread of values. The worked example in Appendix J §J.3 uses an all-zero
// block and expects the "0000000000000000" encoding for sum=0 with
// complement=false; our encoder reproduces that.
func TestAppendixJZero(t *testing.T) {
	enc := Encode(0, false)
	if string(enc[:]) != "0000000000000000" {
		t.Fatalf("zero encoding %q", enc)
	}
	if Decode(enc, false) != 0 {
		t.Fatalf("zero decode")
	}
}

func TestEncodeDecodeRoundTrip(t *testing.T) {
	for _, sum := range []uint32{
		0,
		1,
		0xFFFFFFFF,
		0x12345678,
		0xDEADBEEF,
		0xCAFEBABE,
		0x55555555,
		0xAAAAAAAA,
	} {
		for _, compl := range []bool{false, true} {
			enc := Encode(sum, compl)
			if len(enc) != 16 {
				t.Fatalf("len %d", len(enc))
			}
			// Verify no punctuation characters appear in the output.
			for _, b := range enc {
				for _, ex := range excludeBytes {
					if b == ex {
						t.Fatalf("sum=%#x compl=%v produced excluded char %#x", sum, compl, b)
					}
				}
				if b < 0x30 || b > 0x7e {
					t.Fatalf("sum=%#x produced non-printable %#x", sum, b)
				}
			}
			if got := Decode(enc, compl); got != sum {
				t.Errorf("roundtrip sum=%#x compl=%v: got %#x", sum, compl, got)
			}
		}
	}
}

func TestUpdateAssociativity(t *testing.T) {
	// Two blocks summed in one pass must equal two Updates.
	block := make([]byte, 2880*2)
	for i := range block {
		block[i] = byte(i * 7)
	}
	full := Compute(block)
	piece := Update(0, block[:2880])
	piece = Update(piece, block[2880:])
	if piece != full {
		t.Fatalf("block-split inconsistent: %#x vs %#x", piece, full)
	}
}

// Sanity check: an all-zero FITS block has sum zero.
func TestZeroBlockSum(t *testing.T) {
	z := make([]byte, 2880)
	if Compute(z) != 0 {
		t.Fatalf("zero-block sum non-zero")
	}
	_ = bytes.Repeat // keep import if trimmed
}
