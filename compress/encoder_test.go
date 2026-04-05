package compress

import (
	"bytes"
	"encoding/binary"
	"testing"
)

// TestEncoderDecoderRoundTripNoCompress: trivial identity check.
func TestEncoderDecoderRoundTripNoCompress(t *testing.T) {
	src := makeInt16Bytes([]int16{1, 2, 3, 4, 5, 6, 7, 8})
	enc, _ := SelectEncoder(NoCompress, nil)
	dec, _ := Select(NoCompress, nil)
	roundTripCheck(t, src, 8, 2, enc, dec)
}

// TestEncoderDecoderRoundTripGZIP1: zlib encode → zlib decode identity
// for a sequence that has a lot of compression-friendly repetition.
func TestEncoderDecoderRoundTripGZIP1(t *testing.T) {
	src := makeInt16Bytes([]int16{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
	})
	enc, _ := SelectEncoder(GZIP1, nil)
	dec, _ := Select(GZIP1, nil)
	roundTripCheck(t, src, 30, 2, enc, dec)
}

// TestEncoderDecoderRoundTripGZIP2: byte-shuffle + zlib, for BYTEPIX=2.
func TestEncoderDecoderRoundTripGZIP2(t *testing.T) {
	src := makeInt16Bytes([]int16{0x0102, 0x0304, 0x0506, 0x0708, 0x0102, 0x0304, 0x0506, 0x0708})
	enc, _ := SelectEncoder(GZIP2, Params{"BYTEPIX": 2})
	dec, _ := Select(GZIP2, Params{"BYTEPIX": 2})
	roundTripCheck(t, src, 8, 2, enc, dec)
}

// TestEncoderDecoderRoundTripRICE1Int16: the critical RICE test —
// encode a sequential ramp (small diffs → high compression) and verify
// the decoder recovers it byte-exact.
func TestEncoderDecoderRoundTripRICE1Int16(t *testing.T) {
	pixels := make([]int16, 100)
	for i := range pixels {
		pixels[i] = int16(i)
	}
	src := makeInt16Bytes(pixels)
	enc, _ := SelectEncoder(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 2})
	dec, _ := Select(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 2})
	roundTripCheck(t, src, 100, 2, enc, dec)
}

// TestEncoderDecoderRoundTripRICE1Int32: same for int32 with a larger
// dynamic range (tests BYTEPIX=4 path).
func TestEncoderDecoderRoundTripRICE1Int32(t *testing.T) {
	pixels := make([]int32, 100)
	for i := range pixels {
		pixels[i] = int32(i) * 100
	}
	src := makeInt32Bytes(pixels)
	enc, _ := SelectEncoder(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 4})
	dec, _ := Select(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 4})
	roundTripCheck(t, src, 100, 4, enc, dec)
}

// TestEncoderDecoderRoundTripRICE1LargeDiffs: non-sequential pixels
// that force the decoder into larger unary runs.
func TestEncoderDecoderRoundTripRICE1LargeDiffs(t *testing.T) {
	// Generate a pattern with occasional large jumps.
	pixels := make([]int16, 64)
	for i := range pixels {
		if i%8 == 0 {
			pixels[i] = int16(1000 * (i / 8))
		} else {
			pixels[i] = pixels[i-1] + 1
		}
	}
	src := makeInt16Bytes(pixels)
	enc, _ := SelectEncoder(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 2})
	dec, _ := Select(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 2})
	roundTripCheck(t, src, 64, 2, enc, dec)
}

// TestEncoderDecoderRoundTripHCOMPRESS1: 2D tile through the forward
// H-transform + quadtree + Huffman, decoded back via the inverse.
func TestEncoderDecoderRoundTripHCOMPRESS1(t *testing.T) {
	// 8x8 image = 64 pixels, simple pattern.
	nx, ny := 8, 8
	pixels := make([]int32, nx*ny)
	for i := 0; i < nx; i++ {
		for j := 0; j < ny; j++ {
			pixels[i*ny+j] = int32(i*10 + j)
		}
	}
	src := makeInt32Bytes(pixels)
	enc := &hcompress1Encoder{scale: 0}
	enc.SetShape(nx, ny)
	dec, _ := Select(HCOMPRESS1, nil)
	dst := make([]byte, len(src)*2+100)
	n, err := enc.Encode(src, dst, nx*ny, 4)
	if err != nil {
		t.Fatalf("encode: %v", err)
	}
	round := make([]byte, len(src))
	if err := dec.Decode(dst[:n], round, nx*ny, 4); err != nil {
		t.Fatalf("decode: %v", err)
	}
	for i := 0; i < nx*ny; i++ {
		got := int32(binary.BigEndian.Uint32(round[i*4:]))
		if got != pixels[i] {
			t.Errorf("pixel %d: got %d want %d", i, got, pixels[i])
			if i > 5 {
				break
			}
		}
	}
}

// TestEncoderDecoderRoundTripPLIO1: mask data through the encoder and
// decoder. Uses the exact pattern from our astropy fixture test.
func TestEncoderDecoderRoundTripPLIO1(t *testing.T) {
	// Row 30 of the mask: 5 zeros, 40 twos, 5 zeros.
	pixels := make([]int16, 50)
	for i := 5; i < 45; i++ {
		pixels[i] = 2
	}
	src := makeInt16Bytes(pixels)
	enc, _ := SelectEncoder(PLIO1, nil)
	dec, _ := Select(PLIO1, nil)
	roundTripCheck(t, src, 50, 2, enc, dec)
}

// roundTripCheck runs a src buffer through Encode then Decode and
// byte-compares the result against the original.
func roundTripCheck(t *testing.T, src []byte, nelem, elemSize int, enc Encoder, dec Decoder) {
	t.Helper()
	dst := make([]byte, len(src)*4+100)
	n, err := enc.Encode(src, dst, nelem, elemSize)
	if err != nil {
		t.Fatalf("encode: %v", err)
	}
	round := make([]byte, len(src))
	if err := dec.Decode(dst[:n], round, nelem, elemSize); err != nil {
		t.Fatalf("decode: %v", err)
	}
	if !bytes.Equal(src, round) {
		t.Fatalf("round-trip mismatch: encoded %d→%d bytes, recovered != original", len(src), n)
	}
}

func makeInt16Bytes(p []int16) []byte {
	out := make([]byte, len(p)*2)
	for i, v := range p {
		binary.BigEndian.PutUint16(out[i*2:], uint16(v))
	}
	return out
}

func makeInt32Bytes(p []int32) []byte {
	out := make([]byte, len(p)*4)
	for i, v := range p {
		binary.BigEndian.PutUint32(out[i*4:], uint32(v))
	}
	return out
}
