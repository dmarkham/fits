package block

import (
	"bytes"
	"errors"
	"io"
	"testing"
)

// seekBuf is an in-memory io.ReadWriteSeeker for tests.
type seekBuf struct {
	data []byte
	pos  int64
}

func (s *seekBuf) Read(p []byte) (int, error) {
	if s.pos >= int64(len(s.data)) {
		return 0, io.EOF
	}
	n := copy(p, s.data[s.pos:])
	s.pos += int64(n)
	return n, nil
}

func (s *seekBuf) Write(p []byte) (int, error) {
	need := int(s.pos) + len(p)
	if need > len(s.data) {
		s.data = append(s.data, make([]byte, need-len(s.data))...)
	}
	n := copy(s.data[s.pos:], p)
	s.pos += int64(n)
	return n, nil
}

func (s *seekBuf) Seek(offset int64, whence int) (int64, error) {
	switch whence {
	case io.SeekStart:
		s.pos = offset
	case io.SeekCurrent:
		s.pos += offset
	case io.SeekEnd:
		s.pos = int64(len(s.data)) + offset
	}
	return s.pos, nil
}

func TestReadBlock(t *testing.T) {
	data := make([]byte, Size*3)
	for i := range data {
		data[i] = byte(i / Size) // block 0 = 0x00, block 1 = 0x01, block 2 = 0x02
	}
	r, err := NewReader(&seekBuf{data: data})
	if err != nil {
		t.Fatal(err)
	}
	if r.Size() != Size*3 {
		t.Fatalf("size %d", r.Size())
	}
	if r.BlockCount() != 3 {
		t.Fatalf("blocks %d", r.BlockCount())
	}
	for i := range int64(3) {
		b, err := r.ReadBlock(i)
		if err != nil {
			t.Fatalf("block %d: %v", i, err)
		}
		if b[0] != byte(i) || b[Size-1] != byte(i) {
			t.Fatalf("block %d contents wrong", i)
		}
	}
}

func TestReadBlockCache(t *testing.T) {
	data := bytes.Repeat([]byte{0xAB}, Size*2)
	r, _ := NewReader(&seekBuf{data: data})
	b1, _ := r.ReadBlock(0)
	b2, _ := r.ReadBlock(0)
	if b1 != b2 {
		t.Fatalf("cache did not return same pointer")
	}
}

func TestReadBlockShort(t *testing.T) {
	data := make([]byte, Size-1) // not a full block
	r, _ := NewReader(&seekBuf{data: data})
	if _, err := r.ReadBlock(0); !errors.Is(err, ErrShortBlock) {
		t.Fatalf("expected ErrShortBlock, got %v", err)
	}
}

func TestWriteAndPad(t *testing.T) {
	buf := &seekBuf{}
	w, err := NewWriter(buf)
	if err != nil {
		t.Fatal(err)
	}
	if err := w.WriteRange([]byte("SIMPLE  =")); err != nil {
		t.Fatal(err)
	}
	if err := w.PadToBlock(' '); err != nil {
		t.Fatal(err)
	}
	if len(buf.data) != Size {
		t.Fatalf("padded len %d != %d", len(buf.data), Size)
	}
	if buf.data[Size-1] != ' ' {
		t.Fatalf("tail not space: %v", buf.data[Size-1])
	}
}

func TestRoundUpBlocks(t *testing.T) {
	cases := []struct{ n, want int64 }{
		{0, 0},
		{1, Size},
		{Size, Size},
		{Size + 1, Size * 2},
		{Size * 5, Size * 5},
		{Size*5 + 10, Size * 6},
	}
	for _, c := range cases {
		if got := RoundUpBlocks(c.n); got != c.want {
			t.Errorf("RoundUpBlocks(%d)=%d want %d", c.n, got, c.want)
		}
	}
}
