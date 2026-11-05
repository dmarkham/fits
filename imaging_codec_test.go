package fits

import (
	"bytes"
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"sort"
	"strings"
	"testing"
)

// pixArray is a testdata/imaging truth array (readfloat32_expected.json
// "bits", stdimage_expected.json "values", sirilref/golden "bits"): an
// explicit JSON list (null = a NaN / null pixel), or, past 1024 values,
// the exact digest written by gen_expected.py and sirilref.c: {"n",
// "sha256", "sample"}, sha256 over the big-endian int64 of every value
// (null as -1; the values are float32 bit patterns or pixel values, all
// >= 0) plus [index, value] samples. Check compares a test's full result
// with either form exactly: a digest fails on any differing value (sha256)
// and names the sampled indices that differ.
type pixArray struct {
	explicit []*int64
	n        int
	sha      string
	sample   map[int]int64
}

// pixNull is the canonical value of a null (NaN) pixel.
const pixNull = -1

func (a *pixArray) UnmarshalJSON(b []byte) error {
	b = bytes.TrimSpace(b)
	*a = pixArray{}
	switch {
	case string(b) == "null":
		return nil
	case len(b) > 0 && b[0] == '[':
		var v []*int64
		if err := json.Unmarshal(b, &v); err != nil {
			return err
		}
		if v == nil {
			v = []*int64{}
		}
		for _, p := range v {
			if p != nil && *p < 0 {
				return fmt.Errorf("imaging truth value %d < 0", *p)
			}
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
		return fmt.Errorf("imaging truth digest: %w", err)
	}
	if len(d.SHA256) != 64 || d.N <= 0 || len(d.Sample) == 0 {
		return fmt.Errorf("imaging truth digest: n %d, sha256 %q, %d samples", d.N, d.SHA256, len(d.Sample))
	}
	a.n, a.sha, a.sample = d.N, d.SHA256, map[int]int64{}
	last := -1
	for _, s := range d.Sample {
		if k := int(s[0]); k <= last || k >= d.N {
			return fmt.Errorf("imaging truth digest: sample index %d (n %d)", k, d.N)
		}
		a.sample[int(s[0])], last = s[1], int(s[0])
	}
	return nil
}

// Present reports whether the golden holds values (JSON null: none).
func (a pixArray) Present() bool { return a.explicit != nil || a.n > 0 }

// Digested reports whether the golden stores a digest.
func (a pixArray) Digested() bool { return a.n > 0 }

// Len is the number of values in either form.
func (a pixArray) Len() int {
	if a.explicit != nil {
		return len(a.explicit)
	}
	return a.n
}

// Explicit returns the explicit values (nil entries are nulls); a digest
// fails the test (compare it with Check instead).
func (a pixArray) Explicit(t testing.TB) []*int64 {
	t.Helper()
	if a.n > 0 {
		t.Fatalf("imaging truth array of %d values is a digest; compare it with Check", a.n)
	}
	return a.explicit
}

func pixDigest(v []int64) string {
	h := sha256.New()
	var buf [8]byte
	for _, x := range v {
		binary.BigEndian.PutUint64(buf[:], uint64(x))
		h.Write(buf[:])
	}
	return hex.EncodeToString(h.Sum(nil))
}

// Diff describes every difference between got (pixNull for a null) and
// the golden, "" when equal.
func (a pixArray) Diff(got []int64) string {
	if len(got) != a.Len() {
		return fmt.Sprintf("%d values, golden %d", len(got), a.Len())
	}
	show := func(v int64) string {
		if v == pixNull {
			return "null"
		}
		return fmt.Sprintf("%#x", v)
	}
	var d []string
	if a.explicit != nil {
		n := 0
		for k, p := range a.explicit {
			w := int64(pixNull)
			if p != nil {
				w = *p
			}
			if got[k] != w {
				if n < 8 {
					d = append(d, fmt.Sprintf("[%d] = %s, golden %s", k, show(got[k]), show(w)))
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
			d = append(d, fmt.Sprintf("[%d] = %s, golden %s", k, show(got[k]), show(a.sample[k])))
		}
	}
	if h := pixDigest(got); h != a.sha {
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

// Check reports any difference with t.Errorf and returns true when equal.
func (a pixArray) Check(t testing.TB, where string, got []int64) bool {
	t.Helper()
	if !a.Present() {
		t.Errorf("%s: golden has no values", where)
		return false
	}
	if d := a.Diff(got); d != "" {
		t.Errorf("%s: %s", where, d)
		return false
	}
	return true
}

// f32Canon maps ReadFloat32 output to the truth arrays' values: the
// float32 bit pattern, pixNull for any NaN (readfloat32_expected.json
// records every NaN pixel as null).
func f32Canon(v []float32) []int64 {
	out := make([]int64, len(v))
	for i, f := range v {
		if math.IsNaN(float64(f)) {
			out[i] = pixNull
		} else {
			out[i] = int64(math.Float32bits(f))
		}
	}
	return out
}

// f32Bits maps float32 values to their bit patterns, every NaN as
// 0x7fc00000 when canonNaN (stdimage_expected.json digests; Siril's bits
// are compared raw).
func f32Bits(v []float32, canonNaN bool) []int64 {
	out := make([]int64, len(v))
	for i, f := range v {
		b := math.Float32bits(f)
		if canonNaN && math.IsNaN(float64(f)) {
			b = 0x7fc00000
		}
		out[i] = int64(b)
	}
	return out
}

// TestImagingGoldenCodec pins pixArray on hand-made values: an explicit
// list with a null, a digest (sha256 of the big-endian int64 0 is the
// published sha256 of 8 zero bytes) and its sample/unsampled diffs.
func TestImagingGoldenCodec(t *testing.T) {
	var a pixArray
	if err := json.Unmarshal([]byte(`[1065353216,null,0]`), &a); err != nil {
		t.Fatal(err)
	}
	if d := a.Diff(f32Canon([]float32{1, float32(math.NaN()), 0})); d != "" {
		t.Errorf("explicit Diff = %q", d)
	}
	if h := pixDigest([]int64{0}); h != "af5570f5a1810b7af78caf4bc70a660f0df51e42baf91d4de5b2328de0e83dfc" {
		t.Errorf("sha256(int64 0) = %s", h)
	}
	vals := []int64{7, pixNull, 9}
	var d pixArray
	if err := json.Unmarshal([]byte(fmt.Sprintf(`{"n":3,"sha256":"%s","sample":[[0,7],[2,9]]}`, pixDigest(vals))), &d); err != nil {
		t.Fatal(err)
	}
	if s := d.Diff(vals); s != "" {
		t.Errorf("digest Diff(equal) = %q", s)
	}
	if s := d.Diff([]int64{7, 8, 9}); !strings.Contains(s, "unsampled") {
		t.Errorf("digest Diff(unsampled) = %q", s)
	}
	if s := d.Diff([]int64{6, pixNull, 9}); !strings.Contains(s, "[0] = 0x6, golden 0x7") {
		t.Errorf("digest Diff(sampled) = %q", s)
	}
}
