package testutil

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"fmt"
	"math"
	"strconv"
)

// FloatsDigest is a recorded pixel reading too large to commit value by
// value (Call.Compact): the exact sha256 of the values plus a
// deterministic sample of them, so a mismatch stays debuggable. The
// values are []*float64, nil meaning null; the canonical encoding is one
// byte per value (0 null, 1 value) followed, for a value, by its float64
// bits big-endian. A NaN is a value like any other (its bits count).
type FloatsDigest struct {
	N      int               `json:"n"`
	SHA256 string            `json:"sha256"`
	Sample map[string]string `json:"sample"` // index -> "null" or hex float64 bits
}

// DigestFloats digests v (every value; the sample holds the first and
// last values and every max(1, n/8)-th one).
func DigestFloats(v []*float64) FloatsDigest {
	h := sha256.New()
	var b [9]byte
	for _, p := range v {
		if p == nil {
			h.Write([]byte{0})
			continue
		}
		b[0] = 1
		binary.BigEndian.PutUint64(b[1:], math.Float64bits(*p))
		h.Write(b[:])
	}
	d := FloatsDigest{N: len(v), SHA256: hex.EncodeToString(h.Sum(nil)), Sample: map[string]string{}}
	step := len(v) / 8
	if step < 1 {
		step = 1
	}
	for i := 0; i < len(v); i += step {
		d.Sample[strconv.Itoa(i)] = floatToken(v[i])
	}
	if len(v) > 0 {
		d.Sample[strconv.Itoa(len(v)-1)] = floatToken(v[len(v)-1])
	}
	return d
}

func floatToken(p *float64) string {
	if p == nil {
		return "null"
	}
	return fmt.Sprintf("%016x", math.Float64bits(*p))
}

// Matches reports whether cand is exactly the digested reading; when it
// is not, why names the first sampled difference (or only the digest).
func (d FloatsDigest) Matches(cand []*float64) (ok bool, why string) {
	if len(cand) != d.N {
		return false, fmt.Sprintf("%d values, the recorded reading has %d", len(cand), d.N)
	}
	if got := DigestFloats(cand); got.SHA256 == d.SHA256 {
		return true, ""
	}
	for i := 0; i < len(cand); i++ {
		if w, ok := d.Sample[strconv.Itoa(i)]; ok && w != floatToken(cand[i]) {
			return false, fmt.Sprintf("value %d is %s, the recorded reading has %s", i, floatToken(cand[i]), w)
		}
	}
	return false, "the sha256 differs (no sampled value does)"
}
