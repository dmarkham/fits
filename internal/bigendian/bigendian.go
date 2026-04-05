// Package bigendian provides typed big-endian read/write helpers and bulk
// in-place byte-swapping for FITS data sections.
//
// FITS data is always stored in big-endian byte order regardless of the
// host architecture (FITS Standard v4.0 §5.2).
package bigendian

import (
	"encoding/binary"
	"math"
)

// Int16 decodes one big-endian int16 from b[0:2].
func Int16(b []byte) int16 { return int16(binary.BigEndian.Uint16(b)) }

// Int32 decodes one big-endian int32 from b[0:4].
func Int32(b []byte) int32 { return int32(binary.BigEndian.Uint32(b)) }

// Int64 decodes one big-endian int64 from b[0:8].
func Int64(b []byte) int64 { return int64(binary.BigEndian.Uint64(b)) }

// Uint16 decodes one big-endian uint16 from b[0:2].
func Uint16(b []byte) uint16 { return binary.BigEndian.Uint16(b) }

// Uint32 decodes one big-endian uint32 from b[0:4].
func Uint32(b []byte) uint32 { return binary.BigEndian.Uint32(b) }

// Uint64 decodes one big-endian uint64 from b[0:8].
func Uint64(b []byte) uint64 { return binary.BigEndian.Uint64(b) }

// Float32 decodes one IEEE-754 big-endian float32 from b[0:4].
func Float32(b []byte) float32 {
	return math.Float32frombits(binary.BigEndian.Uint32(b))
}

// Float64 decodes one IEEE-754 big-endian float64 from b[0:8].
func Float64(b []byte) float64 {
	return math.Float64frombits(binary.BigEndian.Uint64(b))
}

// PutInt16 encodes v as big-endian into b[0:2].
func PutInt16(b []byte, v int16) { binary.BigEndian.PutUint16(b, uint16(v)) }

// PutInt32 encodes v as big-endian into b[0:4].
func PutInt32(b []byte, v int32) { binary.BigEndian.PutUint32(b, uint32(v)) }

// PutInt64 encodes v as big-endian into b[0:8].
func PutInt64(b []byte, v int64) { binary.BigEndian.PutUint64(b, uint64(v)) }

// PutUint16 encodes v as big-endian into b[0:2].
func PutUint16(b []byte, v uint16) { binary.BigEndian.PutUint16(b, v) }

// PutUint32 encodes v as big-endian into b[0:4].
func PutUint32(b []byte, v uint32) { binary.BigEndian.PutUint32(b, v) }

// PutUint64 encodes v as big-endian into b[0:8].
func PutUint64(b []byte, v uint64) { binary.BigEndian.PutUint64(b, v) }

// PutFloat32 encodes v as big-endian IEEE-754 into b[0:4].
func PutFloat32(b []byte, v float32) {
	binary.BigEndian.PutUint32(b, math.Float32bits(v))
}

// PutFloat64 encodes v as big-endian IEEE-754 into b[0:8].
func PutFloat64(b []byte, v float64) {
	binary.BigEndian.PutUint64(b, math.Float64bits(v))
}

// BulkSwap16 swaps byte order for every 2-byte pair in b in place.
// len(b) must be a multiple of 2.
func BulkSwap16(b []byte) {
	for i := 0; i+1 < len(b); i += 2 {
		b[i], b[i+1] = b[i+1], b[i]
	}
}

// BulkSwap32 swaps byte order for every 4-byte group in b in place.
// len(b) must be a multiple of 4.
func BulkSwap32(b []byte) {
	for i := 0; i+3 < len(b); i += 4 {
		b[i], b[i+3] = b[i+3], b[i]
		b[i+1], b[i+2] = b[i+2], b[i+1]
	}
}

// BulkSwap64 swaps byte order for every 8-byte group in b in place.
// len(b) must be a multiple of 8.
func BulkSwap64(b []byte) {
	for i := 0; i+7 < len(b); i += 8 {
		b[i], b[i+7] = b[i+7], b[i]
		b[i+1], b[i+6] = b[i+6], b[i+1]
		b[i+2], b[i+5] = b[i+5], b[i+2]
		b[i+3], b[i+4] = b[i+4], b[i+3]
	}
}
