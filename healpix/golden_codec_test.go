package healpix

import (
	"bytes"
	"crypto/sha256"
	"encoding/base64"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"sort"
	"strings"
	"testing"
)

// Compact exact golden encodings written by gen_golden.py (f64, i64;
// see the comment there):
//
//   - f64s: an explicit JSON number list, or "f64:" + base64 of the
//     big-endian float64 bytes of every value. Decodes to the exact values.
//   - i64s: an explicit JSON integer list, or a digest {"n", "sha256",
//     "sample"} of a longer array: sha256 over the big-endian int64 bytes
//     of every value plus [index, value] samples. Check compares a test's
//     own full result against either form exactly (a digest fails on any
//     differing value through the sha256, and names the sampled indices
//     that differ).

// f64s is an exact float64 array (explicit, or base64 big-endian bits).
type f64s []float64

func (a *f64s) UnmarshalJSON(b []byte) error {
	b = bytes.TrimSpace(b)
	if len(b) == 0 || b[0] != '"' {
		var v []float64
		if err := json.Unmarshal(b, &v); err != nil {
			return err
		}
		*a = v
		return nil
	}
	var s string
	if err := json.Unmarshal(b, &s); err != nil {
		return err
	}
	enc, ok := strings.CutPrefix(s, "f64:")
	if !ok {
		return fmt.Errorf("golden float array %.20q: want a list or an f64: string", s)
	}
	raw, err := base64.StdEncoding.DecodeString(enc)
	if err != nil || len(raw)%8 != 0 {
		return fmt.Errorf("golden f64 array: bad base64 (%d bytes, %v)", len(raw), err)
	}
	v := make([]float64, len(raw)/8)
	for i := range v {
		v[i] = math.Float64frombits(binary.BigEndian.Uint64(raw[8*i:]))
	}
	*a = v
	return nil
}

// i64s is an exact int64 array: explicit, or a digest.
type i64s struct {
	explicit []int64
	n        int
	sha      string
	sample   map[int]int64
}

func (a *i64s) UnmarshalJSON(b []byte) error {
	b = bytes.TrimSpace(b)
	*a = i64s{}
	if len(b) > 0 && b[0] == '[' {
		var v []int64
		if err := json.Unmarshal(b, &v); err != nil {
			return err
		}
		if v == nil {
			v = []int64{}
		}
		a.explicit = v
		return nil
	}
	var d struct {
		N      int        `json:"n"`
		SHA256 string     `json:"sha256"`
		Sample [][2]int64 `json:"sample"`
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&d); err != nil {
		return fmt.Errorf("golden int digest: %w", err)
	}
	if len(d.SHA256) != 64 || d.N <= 0 || len(d.Sample) == 0 {
		return fmt.Errorf("golden int digest: n %d, sha256 %q, %d samples", d.N, d.SHA256, len(d.Sample))
	}
	a.n, a.sha, a.sample = d.N, d.SHA256, map[int]int64{}
	last := -1
	for _, s := range d.Sample {
		k := int(s[0])
		if k <= last || k >= d.N {
			return fmt.Errorf("golden int digest: sample index %d (n %d)", k, d.N)
		}
		a.sample[k], last = s[1], k
	}
	return nil
}

// Present reports whether the field was in the golden.
func (a i64s) Present() bool { return a.explicit != nil || a.n > 0 }

// Len is the number of values in either form.
func (a i64s) Len() int {
	if a.explicit != nil {
		return len(a.explicit)
	}
	return a.n
}

func digestInt64s(v []int64) string {
	h := sha256.New()
	var buf [8]byte
	for _, x := range v {
		binary.BigEndian.PutUint64(buf[:], uint64(x))
		h.Write(buf[:])
	}
	return hex.EncodeToString(h.Sum(nil))
}

// Diff describes every difference between got and the golden ("" when
// equal): all differing elements of an explicit array (the first 8
// listed), or for a digest the length, every sampled value and the
// sha256 over all values.
func (a i64s) Diff(got []int64) string {
	if len(got) != a.Len() {
		return fmt.Sprintf("%d values, golden %d", len(got), a.Len())
	}
	var d []string
	if a.explicit != nil {
		n := 0
		for k := range got {
			if got[k] != a.explicit[k] {
				if n < 8 {
					d = append(d, fmt.Sprintf("[%d] = %d, golden %d", k, got[k], a.explicit[k]))
				}
				n++
			}
		}
		if n > 8 {
			d = append(d, fmt.Sprintf("... %d differing values in all", n))
		}
		return strings.Join(d, "; ")
	}
	keys := make([]int, 0, len(a.sample))
	for k := range a.sample {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		if got[k] != a.sample[k] && len(d) < 8 {
			d = append(d, fmt.Sprintf("[%d] = %d, golden %d", k, got[k], a.sample[k]))
		}
	}
	if h := digestInt64s(got); h != a.sha {
		msg := fmt.Sprintf("sha256 over all %d values %s, golden %s", len(got), h[:16], a.sha[:16])
		if len(d) == 0 {
			msg += fmt.Sprintf(" (the %d sampled values agree: the difference is at an unsampled index)", len(keys))
		}
		d = append(d, msg)
	} else if len(d) > 0 {
		d = append(d, "sha256 agrees but samples differ: the golden is inconsistent")
	}
	return strings.Join(d, "; ")
}

// Check reports any difference between got and the golden with t.Errorf
// (prefixed by where) and returns true when they are equal.
func (a i64s) Check(t testing.TB, where string, got []int64) bool {
	t.Helper()
	if !a.Present() {
		t.Errorf("%s: golden array missing", where)
		return false
	}
	if d := a.Diff(got); d != "" {
		t.Errorf("%s: %s", where, d)
		return false
	}
	return true
}

// Resolve returns the golden's values element by element for checks that
// need more than equality: the explicit values, or for a digest got
// itself once it is proven identical to the golden (else the test stops).
func (a i64s) Resolve(t testing.TB, where string, got []int64) []int64 {
	t.Helper()
	if a.explicit != nil {
		return a.explicit
	}
	if !a.Check(t, where, got) {
		t.FailNow()
	}
	return got
}

// i64sFromInts converts int results for Check.
func i64sFromInts(v []int) []int64 {
	out := make([]int64, len(v))
	for i, x := range v {
		out[i] = int64(x)
	}
	return out
}

// TestGoldenCodec pins the decoders against hand-made values: the f64
// base64 big-endian form (bit for bit incl. -0, NaN payload, subnormal)
// and the i64 digest's sha256 (a known vector) and sample checks.
func TestGoldenCodec(t *testing.T) {
	var f f64s
	// 1.0, -0.0, 2^-1074 and a NaN with payload 1, big-endian
	if err := json.Unmarshal([]byte(`"f64:P/AAAAAAAACAAAAAAAAAAAAAAAAAAAABf/AAAAAAAAE="`), &f); err != nil {
		t.Fatal(err)
	}
	want := []uint64{0x3ff0000000000000, 0x8000000000000000, 1, 0x7ff0000000000001}
	if len(f) != len(want) {
		t.Fatalf("decoded %d floats", len(f))
	}
	for i, w := range want {
		if math.Float64bits(f[i]) != w {
			t.Errorf("f64[%d] bits %#x, want %#x", i, math.Float64bits(f[i]), w)
		}
	}
	// sha256 of the 8 zero bytes of int64(0): a published vector.
	if h := digestInt64s([]int64{0}); h != "af5570f5a1810b7af78caf4bc70a660f0df51e42baf91d4de5b2328de0e83dfc" {
		t.Errorf("sha256(int64 0) = %s", h)
	}
	var d i64s
	vals := []int64{5, -1, 7}
	js := fmt.Sprintf(`{"n":3,"sha256":"%s","sample":[[0,5],[2,7]]}`, digestInt64s(vals))
	if err := json.Unmarshal([]byte(js), &d); err != nil {
		t.Fatal(err)
	}
	if s := d.Diff(vals); s != "" {
		t.Errorf("Diff(equal) = %q", s)
	}
	if s := d.Diff([]int64{5, -2, 7}); !strings.Contains(s, "unsampled") {
		t.Errorf("Diff(unsampled change) = %q", s)
	}
	if s := d.Diff([]int64{6, -1, 7}); !strings.Contains(s, "[0] = 6, golden 5") {
		t.Errorf("Diff(sampled change) = %q", s)
	}
}
