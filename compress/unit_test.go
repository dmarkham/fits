package compress

import (
	"bytes"
	"compress/zlib"
	"errors"
	"testing"
)

// TestNOCOMPRESSIdentity verifies that NOCOMPRESS is a byte-copy identity.
func TestNOCOMPRESSIdentity(t *testing.T) {
	src := []byte{1, 2, 3, 4, 5, 6, 7, 8}
	dst := make([]byte, len(src))
	dec := newNoCompress()
	if err := dec.Decode(src, dst, 4, 2); err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(src, dst) {
		t.Fatalf("identity failed: %v vs %v", src, dst)
	}
}

// TestNOCOMPRESSShortSrc rejects a short source buffer.
func TestNOCOMPRESSShortSrc(t *testing.T) {
	dec := newNoCompress()
	err := dec.Decode([]byte{1, 2}, make([]byte, 8), 4, 2)
	if !errors.Is(err, ErrCorrupt) {
		t.Fatalf("expected ErrCorrupt, got %v", err)
	}
}

// TestGZIP1Roundtrip compresses a known sequence via stdlib zlib and
// verifies the decoder recovers it.
func TestGZIP1Roundtrip(t *testing.T) {
	original := []byte{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}
	var buf bytes.Buffer
	w := zlib.NewWriter(&buf)
	w.Write(original)
	w.Close()

	dec := newGZIP1()
	dst := make([]byte, len(original))
	if err := dec.Decode(buf.Bytes(), dst, 8, 2); err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(original, dst) {
		t.Fatalf("roundtrip failed: %v vs %v", original, dst)
	}
}

// TestGZIP2Roundtrip compresses a BYTEPIX=2 sequence after byte-shuffling,
// then verifies the decoder unshuffles correctly.
func TestGZIP2Roundtrip(t *testing.T) {
	// Pixels: [0x0102, 0x0304, 0x0506, 0x0708] as big-endian bytes.
	pixels := []byte{
		0x01, 0x02,
		0x03, 0x04,
		0x05, 0x06,
		0x07, 0x08,
	}
	// Shuffled: all first bytes then all second bytes.
	shuffled := []byte{
		0x01, 0x03, 0x05, 0x07,
		0x02, 0x04, 0x06, 0x08,
	}
	var buf bytes.Buffer
	w := zlib.NewWriter(&buf)
	w.Write(shuffled)
	w.Close()

	dec := newGZIP2(Params{"BYTEPIX": 2})
	dst := make([]byte, len(pixels))
	if err := dec.Decode(buf.Bytes(), dst, 4, 2); err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(pixels, dst) {
		t.Fatalf("GZIP_2 unshuffle failed: %v vs %v", pixels, dst)
	}
}

// TestParseAlgorithm covers every FITS standard ZCMPTYPE value plus
// whitespace handling and case-insensitivity.
func TestParseAlgorithm(t *testing.T) {
	cases := []struct {
		in   string
		want Algorithm
	}{
		{"RICE_1", RICE1},
		{"  RICE_1  ", RICE1},
		{"rice_1", RICE1},
		{"GZIP_1", GZIP1},
		{"GZIP_2", GZIP2},
		{"HCOMPRESS_1", HCOMPRESS1},
		{"PLIO_1", PLIO1},
		{"NOCOMPRESS", NoCompress},
		{"unknown", Unknown},
		{"", Unknown},
	}
	for _, c := range cases {
		if got := ParseAlgorithm(c.in); got != c.want {
			t.Errorf("ParseAlgorithm(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

// TestSelectUnknown returns ErrUnknownCompression for a bogus algorithm.
func TestSelectUnknown(t *testing.T) {
	_, err := Select(Unknown, nil)
	if !errors.Is(err, ErrUnknownCompression) {
		t.Fatalf("expected ErrUnknownCompression, got %v", err)
	}
}

// TestSelectHCompressStub verifies HCOMPRESS_1 returns
// ErrUnsupportedCompression (we deferred its implementation).
func TestSelectHCompressStub(t *testing.T) {
	_, err := Select(HCOMPRESS1, nil)
	if !errors.Is(err, ErrUnsupportedCompression) {
		t.Fatalf("expected ErrUnsupportedCompression, got %v", err)
	}
}
