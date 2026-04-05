package bigendian

import (
	"math"
	"testing"
)

func TestRoundTripInt(t *testing.T) {
	var buf [8]byte
	PutInt16(buf[:], -1)
	if Int16(buf[:]) != -1 {
		t.Fatalf("int16 roundtrip")
	}
	PutInt32(buf[:], -123456789)
	if Int32(buf[:]) != -123456789 {
		t.Fatalf("int32 roundtrip")
	}
	PutInt64(buf[:], math.MinInt64)
	if Int64(buf[:]) != math.MinInt64 {
		t.Fatalf("int64 roundtrip")
	}
}

func TestRoundTripFloat(t *testing.T) {
	var buf [8]byte
	PutFloat32(buf[:], float32(math.Pi))
	if Float32(buf[:]) != float32(math.Pi) {
		t.Fatalf("float32 roundtrip")
	}
	PutFloat64(buf[:], math.E)
	if Float64(buf[:]) != math.E {
		t.Fatalf("float64 roundtrip")
	}
}

func TestBigEndianLayout(t *testing.T) {
	var buf [4]byte
	PutInt32(buf[:], 0x01020304)
	if buf != [4]byte{0x01, 0x02, 0x03, 0x04} {
		t.Fatalf("big-endian layout wrong: % x", buf)
	}
}

func TestBulkSwap16(t *testing.T) {
	b := []byte{0x01, 0x02, 0x03, 0x04}
	BulkSwap16(b)
	if b[0] != 0x02 || b[1] != 0x01 || b[2] != 0x04 || b[3] != 0x03 {
		t.Fatalf("swap16 wrong: % x", b)
	}
}

func TestBulkSwap32(t *testing.T) {
	b := []byte{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08}
	BulkSwap32(b)
	want := []byte{0x04, 0x03, 0x02, 0x01, 0x08, 0x07, 0x06, 0x05}
	for i := range b {
		if b[i] != want[i] {
			t.Fatalf("swap32 wrong: % x", b)
		}
	}
}

func TestBulkSwap64(t *testing.T) {
	b := []byte{0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08}
	BulkSwap64(b)
	want := []byte{0x08, 0x07, 0x06, 0x05, 0x04, 0x03, 0x02, 0x01}
	for i := range b {
		if b[i] != want[i] {
			t.Fatalf("swap64 wrong: % x", b)
		}
	}
}
