package compress

import (
	"bytes"
	"compress/gzip"
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

// TestParseAlgorithm covers every FITS standard ZCMPTYPE value, the
// RICE_ONE alias, trailing blanks, and cfitsio's case-sensitive match.
func TestParseAlgorithm(t *testing.T) {
	cases := []struct {
		in   string
		want Algorithm
	}{
		{"RICE_1", RICE1},
		{"RICE_ONE", RICE1},
		{"RICE_1  ", RICE1},
		{"  RICE_1", Unknown}, // leading blanks are significant in FITS strings
		{"rice_1", Unknown},   // cfitsio: FSTRCMP is case-sensitive
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

// TestSelectHCompress verifies HCOMPRESS_1 now constructs successfully.
func TestSelectHCompress(t *testing.T) {
	d, err := Select(HCOMPRESS1, nil)
	if err != nil {
		t.Fatalf("Select(HCOMPRESS1) returned error: %v", err)
	}
	if d == nil {
		t.Fatal("Select(HCOMPRESS1) returned nil decoder")
	}
}

// TestRiceDecoderRefusesBytepix8: cfitsio has no 64-bit Rice coder
// (ricecomp.c has fits_rdecomp, fits_rdecomp_short and fits_rdecomp_byte
// only), so there is no reference stream for BYTEPIX 8. The decoder
// refuses it with ErrUnsupportedCompression, like the encoder
// (TestEncodersRefuseBytepix8), instead of decoding an unproven layout.
func TestRiceDecoderRefusesBytepix8(t *testing.T) {
	dec, err := Select(RICE1, Params{"BLOCKSIZE": 32, "BYTEPIX": 8})
	if err != nil {
		t.Fatal(err)
	}
	src := make([]byte, 64)
	err = dec.Decode(src, make([]byte, 8*4), 4, 8)
	if !errors.Is(err, ErrUnsupportedCompression) {
		t.Fatalf("Decode elemSize 8: err %v, want ErrUnsupportedCompression", err)
	}
}

// TestGZIPDecodersRefuseWrongLength: a GZIP_1/GZIP_2 tile that inflates
// (stdlib gzip writer, as cfitsio writes gzip) to anything but
// nelem*elemSize bytes is refused with ErrCorrupt: never truncated (the
// leftover bytes are another tile width, see InflateTile) and never
// zero-filled.
func TestGZIPDecodersRefuseWrongLength(t *testing.T) {
	gz := func(n int) []byte {
		var buf bytes.Buffer
		w := gzip.NewWriter(&buf)
		w.Write(bytes.Repeat([]byte{7}, n))
		w.Close()
		return buf.Bytes()
	}
	for _, algo := range []Algorithm{GZIP1, GZIP2} {
		dec, _ := Select(algo, nil)
		for _, n := range []int{15, 17, 32} { // want 16 (8 pixels x 2 bytes)
			if err := dec.Decode(gz(n), make([]byte, 16), 8, 2); !errors.Is(err, ErrCorrupt) {
				t.Errorf("%v: %d inflated bytes for 16: err %v, want ErrCorrupt", algo, n, err)
			}
		}
		if err := dec.Decode(gz(16), make([]byte, 16), 8, 2); err != nil {
			t.Errorf("%v: exact length: %v", algo, err)
		}
	}
}

// TestInflateTileWidths: InflateTile takes the value width from the
// inflated length as cfitsio does (1, 2, 4, 8 x nelem, imcompress.c:
// 6372-6417), un-shuffles GZIP_2 by that width, and refuses any other
// length. (The decode of real cfitsio/astropy files with such tiles is
// proven in the fits package, TestCompressedCorpusRead real_m13_gzip and
// the *_w<N> cases.)
func TestInflateTileWidths(t *testing.T) {
	gz := func(b []byte) []byte {
		var buf bytes.Buffer
		w := gzip.NewWriter(&buf)
		w.Write(b)
		w.Close()
		return buf.Bytes()
	}
	// Two 4-byte values 0x01020304, 0x05060708; shuffled: 01 05 02 06 03 07 04 08.
	plain := []byte{1, 2, 3, 4, 5, 6, 7, 8}
	shuf := []byte{1, 5, 2, 6, 3, 7, 4, 8}
	for _, c := range []struct {
		algo  Algorithm
		in    []byte
		nelem int
		width int
	}{
		{GZIP1, plain, 8, 1}, {GZIP1, plain, 4, 2}, {GZIP1, plain, 2, 4}, {GZIP1, plain, 1, 8},
		{GZIP2, shuf, 2, 4},
	} {
		got, w, err := InflateTile(c.algo, gz(c.in), c.nelem)
		if err != nil || w != c.width {
			t.Fatalf("%v nelem %d: width %d err %v, want %d", c.algo, c.nelem, w, err, c.width)
		}
		if c.width == 4 && !bytes.Equal(got, plain) {
			t.Errorf("%v: %v, want %v", c.algo, got, plain)
		}
	}
	for _, n := range []int{3, 5, 7} { // 8 bytes is not 1/2/4/8 x n
		if _, _, err := InflateTile(GZIP1, gz(plain), n); !errors.Is(err, ErrCorrupt) {
			t.Errorf("nelem %d: err %v, want ErrCorrupt", n, err)
		}
	}
	if _, _, err := InflateTile(GZIP1, gz(append(plain, 0)), 2); !errors.Is(err, ErrCorrupt) {
		t.Errorf("leftover byte: err %v, want ErrCorrupt", err)
	}
}
