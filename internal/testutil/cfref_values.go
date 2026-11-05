package testutil

import (
	"bytes"
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"path/filepath"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"testing"
)

// A harness value array (an image's pixels, a flattened table column) is
// either explicit, or, past cfref.c's DIGEST_MIN elements, a digest: the
// element count, a sha256 over a canonical big-endian encoding of every
// element, every null as a run, and explicit values at a deterministic
// sample of indices (first/last 8, image corners, first min/max, 64
// strided). The encoding is specified in cfref.c (jdigest) and mirrored
// by the canon methods below. Check verifies a result against either form
// exactly: for a digest the hash must match (any single differing element
// or null flag fails) and every sampled value and null run is compared
// too, so a failure names the indices that differ.

// Bit is one TBIT element (0 or 1).
type Bit int

// Cplx is one TDBLCOMPLEX element [re, im].
type Cplx [2]Float

// Elem is a harness array element type.
type Elem interface {
	NullFloat | NullInt64 | NullUint64 | NullString | NullBool | Bit | Cplx
}

// canonNaN is the one NaN bit pattern the canonical encoding uses.
const canonNaN = 0x7ff8000000000000

func floatBits(v float64) uint64 {
	if math.IsNaN(v) {
		return canonNaN
	}
	return math.Float64bits(v)
}

// appendCanon appends the canonical encoding of e (cfref.c jdigest).
func appendCanon[T Elem](b []byte, e T) []byte {
	flag := func(valid bool) byte {
		if valid {
			return 1
		}
		return 0
	}
	switch v := any(e).(type) {
	case NullFloat:
		var u uint64
		if v.Valid {
			u = floatBits(v.V)
		}
		return binary.BigEndian.AppendUint64(append(b, flag(v.Valid)), u)
	case NullInt64:
		var u uint64
		if v.Valid {
			u = uint64(v.V)
		}
		return binary.BigEndian.AppendUint64(append(b, flag(v.Valid)), u)
	case NullUint64:
		var u uint64
		if v.Valid {
			u = v.V
		}
		return binary.BigEndian.AppendUint64(append(b, flag(v.Valid)), u)
	case NullString:
		s := ""
		if v.Valid {
			s = v.V
		}
		b = binary.BigEndian.AppendUint32(append(b, flag(v.Valid)), uint32(len(s)))
		return append(b, s...)
	case NullBool:
		return append(b, flag(v.Valid), flag(v.Valid && v.V))
	case Bit:
		return append(b, flag(v != 0))
	case Cplx:
		b = binary.BigEndian.AppendUint64(b, floatBits(float64(v[0])))
		return binary.BigEndian.AppendUint64(b, floatBits(float64(v[1])))
	}
	panic("unreachable")
}

func isNull[T Elem](e T) bool {
	switch v := any(e).(type) {
	case NullFloat:
		return !v.Valid
	case NullInt64:
		return !v.Valid
	case NullUint64:
		return !v.Valid
	case NullString:
		return !v.Valid
	case NullBool:
		return !v.Valid
	}
	return false
}

// same reports whether a and b have the same canonical encoding (floats
// bit for bit, every NaN alike, -0 != +0; nulls alike whatever the value).
func same[T Elem](a, b T) bool {
	var x, y [32]byte
	return bytes.Equal(appendCanon(x[:0], a), appendCanon(y[:0], b))
}

// show formats one element for a diff.
func show[T Elem](e T) string {
	switch v := any(e).(type) {
	case NullFloat:
		if !v.Valid {
			return "null"
		}
		return fmtFloat(v.V)
	case NullInt64:
		if !v.Valid {
			return "null"
		}
		return strconv.FormatInt(v.V, 10)
	case NullUint64:
		if !v.Valid {
			return "null"
		}
		return strconv.FormatUint(v.V, 10)
	case NullString:
		if !v.Valid {
			return "null"
		}
		return strconv.Quote(v.V)
	case NullBool:
		if !v.Valid {
			return "null"
		}
		return strconv.FormatBool(v.V)
	case Bit:
		return strconv.Itoa(int(v))
	case Cplx:
		return "(" + fmtFloat(float64(v[0])) + ", " + fmtFloat(float64(v[1])) + ")"
	}
	return fmt.Sprint(e)
}

func fmtFloat(v float64) string {
	return fmt.Sprintf("%s [%#016x]", strconv.FormatFloat(v, 'g', -1, 64), math.Float64bits(v))
}

// Digest is the stored form of a large harness array.
type Digest struct {
	N      int64                `json:"n"`
	SHA256 string               `json:"sha256"`
	Nulls  [][2]int64           `json:"nulls"`  // [start, count] runs, ascending
	Sample [][2]json.RawMessage `json:"sample"` // [index, value], ascending
}

// Values is a harness value array: explicit, or a Digest.
type Values[T Elem] struct {
	present  bool // false: the harness wrote null (the read failed)
	explicit []T
	dig      *Digest
	sample   map[int64]T
}

// UnmarshalJSON accepts an explicit array, a digest object, a pool
// reference, or null.
func (a *Values[T]) UnmarshalJSON(b []byte) error {
	b = bytes.TrimSpace(b)
	*a = Values[T]{}
	switch {
	case string(b) == "null":
		return nil
	case len(b) > 0 && b[0] == '[':
		var v []T
		if err := json.Unmarshal(b, &v); err != nil {
			return err
		}
		if v == nil {
			v = []T{}
		}
		a.present, a.explicit = true, v
		return nil
	}
	if bytes.Contains(b, []byte(`"pool"`)) {
		dec := json.NewDecoder(bytes.NewReader(b))
		dec.DisallowUnknownFields()
		var r poolRef
		if err := dec.Decode(&r); err != nil {
			return fmt.Errorf("testutil: pooled harness array: %w", err)
		}
		v, err := loadPooled[T](cfrefPoolDir(), r)
		if err != nil {
			return err
		}
		a.present, a.explicit = true, v
		return nil
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	var d Digest
	if err := dec.Decode(&d); err != nil {
		return fmt.Errorf("testutil: harness digest: %w", err)
	}
	if len(d.SHA256) != 64 {
		return fmt.Errorf("testutil: harness digest sha256 %q", d.SHA256)
	}
	a.sample = make(map[int64]T, len(d.Sample))
	last := int64(-1)
	for _, s := range d.Sample {
		k, err := strconv.ParseInt(string(s[0]), 10, 64)
		if err != nil || k <= last || k >= d.N {
			return fmt.Errorf("testutil: harness digest sample index %s (n %d)", s[0], d.N)
		}
		var v T
		if err := json.Unmarshal(s[1], &v); err != nil {
			return fmt.Errorf("testutil: harness digest sample %d: %w", k, err)
		}
		a.sample[k], last = v, k
	}
	end := int64(0)
	for _, r := range d.Nulls {
		if r[0] < end || r[1] <= 0 || r[0]+r[1] > d.N {
			return fmt.Errorf("testutil: harness digest null run %v (n %d)", r, d.N)
		}
		end = r[0] + r[1]
	}
	a.present, a.dig = true, &d
	return nil
}

// A mid-size harness array (cfref.c "dump -pool": 65..1024 values) may be
// stored once in testdata/cfref/golden/pool and referenced from the golden
// as {"n": N, "sha256": "<hex>", "pool": "<kind>"}; the pool file
// <kind>-<sha256>.json holds the explicit JSON array. Loading it checks
// the element kind, the count and the sha256 of the canonical encoding,
// so the values a test sees are exactly the ones cfitsio wrote; the array
// then behaves as an explicit one.

// poolRef is the stored form of a pooled array.
type poolRef struct {
	N      int64  `json:"n"`
	SHA256 string `json:"sha256"`
	Pool   string `json:"pool"`
}

// cfrefPoolDir is <repo>/testdata/cfref/golden/pool, found from this
// source file (UnmarshalJSON has no testing.TB to locate the repo with).
var cfrefPoolDir = sync.OnceValue(func() string {
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		return ""
	}
	return filepath.Join(filepath.Dir(file), "..", "..", "testdata", "cfref", "golden", "pool")
})

// poolKind is the cfref.c kind name of element type T.
func poolKind[T Elem]() string {
	var z T
	switch any(z).(type) {
	case NullFloat:
		return "f64"
	case NullInt64:
		return "i64"
	case NullUint64:
		return "u64"
	case NullString:
		return "str"
	case NullBool:
		return "log"
	case Bit:
		return "bit"
	case Cplx:
		return "cplx"
	}
	return ""
}

var poolCache sync.Map // pool file path -> []byte

// loadPooled returns the explicit values of a pooled array.
func loadPooled[T Elem](dir string, r poolRef) ([]T, error) {
	if r.Pool != poolKind[T]() {
		return nil, fmt.Errorf("testutil: pooled %s array where %s values are expected", r.Pool, poolKind[T]())
	}
	if len(r.SHA256) != 64 || strings.Trim(r.SHA256, "0123456789abcdef") != "" {
		return nil, fmt.Errorf("testutil: pooled array sha256 %q", r.SHA256)
	}
	p := filepath.Join(dir, r.Pool+"-"+r.SHA256+".json")
	var b []byte
	if c, ok := poolCache.Load(p); ok {
		b = c.([]byte)
	} else {
		var err error
		if b, err = ReadGolden(p); err != nil {
			return nil, fmt.Errorf("testutil: committed pool file missing (a failure, not a skip): %w", err)
		}
		poolCache.Store(p, b)
	}
	var v []T
	if err := json.Unmarshal(b, &v); err != nil {
		return nil, fmt.Errorf("testutil: pool file %s: %w", p, err)
	}
	if int64(len(v)) != r.N {
		return nil, fmt.Errorf("testutil: pool file %s holds %d values, the golden says %d", p, len(v), r.N)
	}
	if h := digestOf(v); h != r.SHA256 {
		return nil, fmt.Errorf("testutil: pool file %s hashes to %s, not its name", p, h)
	}
	return v, nil
}

// ExplicitValues wraps explicit values (e.g. a flattened table column).
func ExplicitValues[T Elem](v []T) Values[T] {
	if v == nil {
		v = []T{}
	}
	return Values[T]{present: true, explicit: v}
}

// Present is false when the harness wrote null (the cfitsio read failed;
// the HDU's errors say why) or the field was absent.
func (a Values[T]) Present() bool { return a.present }

// Digested reports whether the golden stores a digest, not every value.
func (a Values[T]) Digested() bool { return a.dig != nil }

// Len is the number of elements in either form.
func (a Values[T]) Len() int {
	if a.dig != nil {
		return int(a.dig.N)
	}
	return len(a.explicit)
}

// Values returns the explicit values; a digested array fails the test
// (compare it with Check / CheckFunc instead, which verify every element).
func (a Values[T]) Values(t testing.TB) []T {
	t.Helper()
	if a.dig != nil {
		t.Fatalf("testutil: harness array of %d values is stored as a digest; compare it with Check/CheckFunc", a.dig.N)
	}
	return a.explicit
}

// IsNull reports cfitsio's null flag for element k in either form.
func (a Values[T]) IsNull(k int) bool {
	if a.dig == nil {
		return isNull(a.explicit[k])
	}
	i := sort.Search(len(a.dig.Nulls), func(i int) bool { return a.dig.Nulls[i][0]+a.dig.Nulls[i][1] > int64(k) })
	return i < len(a.dig.Nulls) && a.dig.Nulls[i][0] <= int64(k)
}

// NullCount is the number of null elements.
func (a Values[T]) NullCount() int {
	n := 0
	if a.dig == nil {
		for _, e := range a.explicit {
			if isNull(e) {
				n++
			}
		}
		return n
	}
	for _, r := range a.dig.Nulls {
		n += int(r[1])
	}
	return n
}

// Sampled returns the explicit indices available in either form: every
// index for an explicit array, the sample indices for a digest.
func (a Values[T]) Sampled() map[int]T {
	out := map[int]T{}
	if a.dig == nil {
		for k, e := range a.explicit {
			out[k] = e
		}
		return out
	}
	for k, e := range a.sample {
		out[int(k)] = e
	}
	return out
}

// SHA256 is the hex digest of the canonical encoding (computed for an
// explicit array).
func (a Values[T]) SHA256() string {
	if a.dig != nil {
		return a.dig.SHA256
	}
	return digestOf(a.explicit)
}

func digestOf[T Elem](v []T) string {
	h := sha256.New()
	var buf []byte
	for _, e := range v {
		buf = appendCanon(buf[:0], e)
		h.Write(buf)
	}
	return hex.EncodeToString(h.Sum(nil))
}

func nullRuns[T Elem](v []T) [][2]int64 {
	var out [][2]int64
	for k := 0; k < len(v); {
		if !isNull(v[k]) {
			k++
			continue
		}
		s := k
		for k < len(v) && isNull(v[k]) {
			k++
		}
		out = append(out, [2]int64{int64(s), int64(k - s)})
	}
	return out
}

const maxDiffs = 12

// Diff compares got with the golden exactly and returns a readable
// description of every difference found ("" when equal). For an explicit
// golden every element is compared; for a digest the element count, the
// null runs, every sampled index and the sha256 of all elements.
func (a Values[T]) Diff(got []T) string {
	var d []string
	if !a.present {
		return "golden has no values (the cfitsio read failed; see the HDU errors)"
	}
	if len(got) != a.Len() {
		return fmt.Sprintf("%d values, cfitsio %d", len(got), a.Len())
	}
	if a.dig == nil {
		n := 0
		for k := range got {
			if !same(got[k], a.explicit[k]) {
				if n < maxDiffs {
					d = append(d, fmt.Sprintf("[%d] = %s, cfitsio %s", k, show(got[k]), show(a.explicit[k])))
				}
				n++
			}
		}
		if n > maxDiffs {
			d = append(d, fmt.Sprintf("... %d differing elements in all", n))
		}
		return strings.Join(d, "; ")
	}
	if g, w := nullRuns(got), a.dig.Nulls; fmt.Sprint(g) != fmt.Sprint(w) {
		d = append(d, fmt.Sprintf("null runs %s, cfitsio %s", clip(fmt.Sprint(g)), clip(fmt.Sprint(w))))
	}
	keys := make([]int64, 0, len(a.sample))
	for k := range a.sample {
		keys = append(keys, k)
	}
	sort.Slice(keys, func(i, j int) bool { return keys[i] < keys[j] })
	n := 0
	for _, k := range keys {
		if w := a.sample[k]; !same(got[k], w) {
			if n < maxDiffs {
				d = append(d, fmt.Sprintf("[%d] = %s, cfitsio %s", k, show(got[k]), show(w)))
			}
			n++
		}
	}
	if n > maxDiffs {
		d = append(d, fmt.Sprintf("... %d of %d sampled elements differ", n, len(keys)))
	}
	if h := digestOf(got); h != a.dig.SHA256 {
		msg := fmt.Sprintf("sha256 over all %d values %s, cfitsio %s", len(got), h[:16], a.dig.SHA256[:16])
		if len(d) == 0 {
			msg += fmt.Sprintf(" (the %d sampled values and the null runs agree: the difference is at an unsampled index)", len(keys))
		}
		d = append(d, msg)
	} else if len(d) > 0 {
		d = append(d, "sha256 agrees but samples differ: the golden is inconsistent")
	}
	return strings.Join(d, "; ")
}

func clip(s string) string {
	if len(s) > 200 {
		return s[:200] + "..."
	}
	return s
}

// Check compares got with the golden exactly (see Diff) and reports any
// difference with t.Errorf, prefixed by where. It returns true when equal.
func (a Values[T]) Check(t testing.TB, where string, got []T) bool {
	t.Helper()
	if d := a.Diff(got); d != "" {
		t.Errorf("%s: %s", where, d)
		return false
	}
	return true
}

// CheckFunc is Check for got[k] = at(k), k in [0, Len()).
func (a Values[T]) CheckFunc(t testing.TB, where string, at func(k int) T) bool {
	t.Helper()
	got := make([]T, a.Len())
	for k := range got {
		got[k] = at(k)
	}
	return a.Check(t, where, got)
}

// Same reports whether two harness arrays hold the same values (both
// cfitsio views; either form). It returns a description when they differ.
func (a Values[T]) Same(b Values[T]) string {
	switch {
	case a.present != b.present:
		return fmt.Sprintf("present %v, cfitsio %v", a.present, b.present)
	case !a.present:
		return ""
	case a.dig == nil:
		return b.Diff(a.explicit)
	case b.dig == nil:
		return "harness wrote " + a.Diff(b.explicit) // a's digest vs b's values
	}
	if a.dig.N != b.dig.N || a.dig.SHA256 != b.dig.SHA256 {
		return fmt.Sprintf("%d values sha256 %s, cfitsio %d values %s", a.dig.N, a.dig.SHA256[:16], b.dig.N, b.dig.SHA256[:16])
	}
	return ""
}

// NullFloats builds checkable values from a result and its validity mask
// (valid nil: every element valid).
func NullFloats(v []float64, valid []bool) []NullFloat {
	out := make([]NullFloat, len(v))
	for k, x := range v {
		ok := valid == nil || valid[k]
		if ok {
			out[k] = NullFloat{V: x, Valid: true}
		}
	}
	return out
}

// NullInt64s is NullFloats for int64 results.
func NullInt64s(v []int64, valid []bool) []NullInt64 {
	out := make([]NullInt64, len(v))
	for k, x := range v {
		if valid == nil || valid[k] {
			out[k] = NullInt64{V: x, Valid: true}
		}
	}
	return out
}

// NullUint64s is NullFloats for uint64 results.
func NullUint64s(v []uint64, valid []bool) []NullUint64 {
	out := make([]NullUint64, len(v))
	for k, x := range v {
		if valid == nil || valid[k] {
			out[k] = NullUint64{V: x, Valid: true}
		}
	}
	return out
}

// Resolve returns the golden's values element by element, for tests that
// need more than an equality check (per-element rules, relations between
// goldens). An explicit golden returns its own values and cand is not
// called. For a digest, cand (the test's own reading, mapped into
// cfitsio's terms) is verified against the golden exactly (Diff: count,
// null runs, every sample, sha256 over all elements); on any difference
// the test fails with the diff, otherwise cand is returned: it is then,
// element for element, the golden. cand nil, or returning nil, fails a
// digested golden (no way to recover its values).
func (a Values[T]) Resolve(t testing.TB, where string, cand func() []T) []T {
	t.Helper()
	if a.dig == nil {
		return a.explicit
	}
	if cand == nil {
		t.Fatalf("%s: golden is a digest of %d values and the test has no candidate reading to verify against it", where, a.dig.N)
	}
	got := cand()
	if got == nil {
		t.Fatalf("%s: golden is a digest of %d values and the candidate reading failed", where, a.dig.N)
	}
	if d := a.Diff(got); d != "" {
		t.Fatalf("%s: %s", where, d)
	}
	return got
}

// CellCandidates supplies, per value kind, the test's own reading of a
// whole column (all rows, flattened) for CfrefColumn.ResolveRows. Only
// the kinds a digested column carries are called.
type CellCandidates struct {
	Phys func() []NullFloat
	I64  func() []NullInt64
	U64  func() []NullUint64
	Raw  func() []NullInt64
	Str  func() []NullString
	Log  func() []NullBool
	Bits func() []Bit
	Cplx func() []Cplx
}

// ResolveRows returns the column's per-row cells in either form: the
// explicit rows, or, for a digested column, each kind resolved against
// its candidate (Values.Resolve: exact, else the test fails) and split
// into rows of PerRow elements.
func (c *CfrefColumn) ResolveRows(t testing.TB, where string, cand CellCandidates) []CfrefCell {
	t.Helper()
	if c.cells == nil {
		return c.rows
	}
	f := c.cells
	n := f.Len()
	if f.PerRow <= 0 || n%f.PerRow != 0 {
		t.Fatalf("%s: digested column with %d values, %d per row", where, n, f.PerRow)
	}
	rows := make([]CfrefCell, n/f.PerRow)
	split := func(set func(r, lo, hi int)) {
		for r := range rows {
			set(r, r*f.PerRow, (r+1)*f.PerRow)
		}
	}
	if f.Str.Present() {
		v := f.Str.Resolve(t, where+" str", cand.Str)
		split(func(r, lo, hi int) { rows[r].Str = v[lo:hi:hi] })
	}
	if f.Log.Present() {
		v := f.Log.Resolve(t, where+" log", cand.Log)
		split(func(r, lo, hi int) { rows[r].Log = v[lo:hi:hi] })
	}
	if f.Bits.Present() {
		v := f.Bits.Resolve(t, where+" bits", cand.Bits)
		split(func(r, lo, hi int) {
			for _, b := range v[lo:hi] {
				rows[r].Bits = append(rows[r].Bits, int(b))
			}
		})
	}
	if f.Cplx.Present() {
		v := f.Cplx.Resolve(t, where+" cplx", cand.Cplx)
		split(func(r, lo, hi int) {
			for _, x := range v[lo:hi] {
				rows[r].Cplx = append(rows[r].Cplx, [2]Float(x))
			}
		})
	}
	if f.Phys.Present() {
		v := f.Phys.Resolve(t, where+" phys", cand.Phys)
		split(func(r, lo, hi int) { rows[r].Phys = v[lo:hi:hi] })
	}
	if f.I64.Present() {
		v := f.I64.Resolve(t, where+" i64", cand.I64)
		split(func(r, lo, hi int) { rows[r].I64 = v[lo:hi:hi] })
	}
	if f.U64.Present() {
		v := f.U64.Resolve(t, where+" u64", cand.U64)
		split(func(r, lo, hi int) { rows[r].U64 = v[lo:hi:hi] })
	}
	if f.Raw.Present() {
		v := f.Raw.Resolve(t, where+" raw", cand.Raw)
		split(func(r, lo, hi int) { rows[r].Raw = v[lo:hi:hi] })
	}
	return rows
}

// Resolved returns a copy of the column with explicit rows (ResolveRows),
// so per-row code (RowCells) works on a digested column too.
func (c CfrefColumn) Resolved(t testing.TB, where string, cand CellCandidates) CfrefColumn {
	t.Helper()
	if c.cells == nil {
		return c
	}
	c.rows, c.cells = c.ResolveRows(t, where, cand), nil
	return c
}
