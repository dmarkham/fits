package stats

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"slices"
	"strconv"
	"strings"
	"testing"
)

// Tests in this file compare the Go functions bit-for-bit with Siril's own
// code: findMinMaxPercentile (src/rt/rt_algo.cc), histogram_median_float,
// quickmedian_float (src/algos/sorting.c), siril_stats_float_sd,
// siril_stats_float_mad and reassign_to_non_null_data_float
// (src/algos/statistics_float.c), copied verbatim into the C harness at
// testdata/sref and run on the inputs committed in the golden file, or on
// inputs both sides generate: a "formula" (formulaData below) or a "gen"
// spec (genValues32 / genValues64). A generated input of more than 64
// values carries the sha256 of its bits (in_sha256), checked before any
// result, and its FilterNonZero output a count and sha256.
//
// The generic (non-[]float32) Percentile path is compared with the same
// verbatim findMinMaxPercentile compiled with float replaced by double.
//
// Regenerate: make -C testdata/sref goldens  (see testdata/sref/PROVENANCE.md)

const sirilGoldenPath = "testdata/sref/golden/siril_stats_golden.json"

type sirilGolden struct {
	Percentiles []string       `json:"percentiles"`
	F32         []sirilF32Case `json:"f32"`
	F64         []sirilF64Case `json:"f64"`
}

type formulaSpec struct {
	Kind string `json:"kind"`
	Seed uint32 `json:"seed"`
}

// genSpec is a "gen" input of sref.cc (the generators it used to commit
// the output of). Lo and Hi are IEEE bit patterns in hex (float32 for the
// []float32 cases, float64 for "uniform" in the generic ones); IntLo and
// IntHi are the "ints" bounds.
type genSpec struct {
	Kind string          `json:"kind"`
	Seed uint32          `json:"seed"`
	Lo   json.RawMessage `json:"lo"`
	Hi   json.RawMessage `json:"hi"`
}

type sirilF32Case struct {
	Name     string       `json:"name"`
	N        int          `json:"n"`
	Input    *string      `json:"input"`
	Formula  *formulaSpec `json:"formula"`
	Gen      *genSpec     `json:"gen"`
	InSHA256 string       `json:"in_sha256"`
	// FilterNonZero truth for generated inputs: count and sha256.
	NonnullN      *int   `json:"nonnull_n"`
	NonnullSHA256 string `json:"nonnull_sha256"`
	SD            struct {
		Mean  string `json:"mean"`
		Stdev string `json:"stdev"`
	} `json:"sd"`
	Nonnull     *string           `json:"nonnull"`
	Pct         map[string]string `json:"pct"`
	Median      *string           `json:"median"`
	Mad         *string           `json:"mad"`
	Quickmedian *string           `json:"quickmedian"`
	NClean      int               `json:"nclean"`
}

type sirilF64Case struct {
	Name     string            `json:"name"`
	Dtype    string            `json:"dtype"`
	N        int               `json:"n"`
	Input    *string           `json:"input"`
	Formula  *formulaSpec      `json:"formula"`
	Gen      *genSpec          `json:"gen"`
	InSHA256 string            `json:"in_sha256"`
	Pct      map[string]string `json:"pct"`
	Median   string            `json:"median"`
	Mad      string            `json:"mad"`
}

// mixedState mirrors mixed_state() in sref.cc.
func mixedState(seed uint32) uint32 {
	x := seed * 2654435761
	if x == 0 {
		x = 1
	}
	for range 8 {
		xorshift32(&x)
	}
	return x
}

func specHex32(t *testing.T, raw json.RawMessage) float32 {
	t.Helper()
	var s string
	if err := json.Unmarshal(raw, &s); err != nil {
		t.Fatalf("gen bound %s: %v", raw, err)
	}
	return hex32(t, s)
}

func specHex64(t *testing.T, raw json.RawMessage) float64 {
	t.Helper()
	var s string
	if err := json.Unmarshal(raw, &s); err != nil {
		t.Fatalf("gen bound %s: %v", raw, err)
	}
	return hex64(t, s)
}

// genValues32 mirrors the float32 generators of sref.cc (gen_uniform,
// gen_gradient, gen_ties3, gen_adu_quarter and the zeroing / NaN
// variants). Every float32 operation is written with an explicit
// conversion, so it is rounded once as in C (and never fused).
func genValues32(t *testing.T, g genSpec, n int) []float32 {
	t.Helper()
	out := make([]float32, n)
	uniform := func(lo, hi float32) {
		x := mixedState(g.Seed)
		for i := range out {
			u := float32(float32(xorshift32(&x)>>8) * 5.9604644775390625e-8)
			out[i] = float32(lo + float32(float32(hi-lo)*u))
		}
	}
	gradient := func(lo, hi float32) {
		for i := range out {
			out[i] = float32(lo + float32(float32(float32(hi-lo)*float32(i))/float32(n-1)))
		}
	}
	switch g.Kind {
	case "uniform":
		uniform(specHex32(t, g.Lo), specHex32(t, g.Hi))
	case "uniform_nan9":
		uniform(specHex32(t, g.Lo), specHex32(t, g.Hi))
		for i := 0; i < n; i += 9 {
			out[i] = float32(math.NaN())
		}
	case "gradient":
		gradient(specHex32(t, g.Lo), specHex32(t, g.Hi))
	case "gradient_zero60":
		gradient(specHex32(t, g.Lo), specHex32(t, g.Hi))
		for i := 0; i < 60 && i < n; i++ {
			out[i] = 0
		}
	case "ties3":
		three := [3]float32{0.1, 0.2, 0.3}
		for i := range out {
			out[i] = three[(i*7)%3]
		}
	case "adu_quarter":
		x := mixedState(g.Seed)
		for i := range out {
			out[i] = float32(30000 + float32(float32(xorshift32(&x)%41)*0.25))
		}
	default:
		t.Fatalf("unknown gen kind %q", g.Kind)
	}
	return out
}

// genValues64 mirrors gen_ints and gen_f64 of sref.cc.
func genValues64(t *testing.T, g genSpec, n int) []float64 {
	t.Helper()
	out := make([]float64, n)
	x := mixedState(g.Seed)
	switch g.Kind {
	case "ints":
		var lo, hi int64
		if json.Unmarshal(g.Lo, &lo) != nil || json.Unmarshal(g.Hi, &hi) != nil {
			t.Fatalf("ints bounds %s %s", g.Lo, g.Hi)
		}
		span := uint64(hi-lo) + 1
		for i := range out {
			out[i] = float64(lo + int64(uint64(xorshift32(&x))%span))
		}
	case "uniform":
		lo, hi := specHex64(t, g.Lo), specHex64(t, g.Hi)
		for i := range out {
			a := xorshift32(&x)
			b := xorshift32(&x)
			u := float64(float64((uint64(a)<<32|uint64(b))>>11) * 0x1p-53)
			out[i] = float64(lo + float64(float64(hi-lo)*u))
		}
	default:
		t.Fatalf("unknown gen kind %q", g.Kind)
	}
	return out
}

// sha32 / sha64: sha256 over the big-endian IEEE bit patterns, as
// sref.cc's sha32 / sha64.
func sha32(d []float32) string {
	h := sha256.New()
	var b [4]byte
	for _, v := range d {
		binary.BigEndian.PutUint32(b[:], math.Float32bits(v))
		h.Write(b[:])
	}
	return hex.EncodeToString(h.Sum(nil))
}

func sha64(d []float64) string {
	h := sha256.New()
	var b [8]byte
	for _, v := range d {
		binary.BigEndian.PutUint64(b[:], math.Float64bits(v))
		h.Write(b[:])
	}
	return hex.EncodeToString(h.Sum(nil))
}

// checkInputDigest: a generated input of more than 64 values must carry
// the harness's digest, and Go must have generated the same bits.
func checkInputDigest(t *testing.T, name string, generated bool, n int, want, got func() string) {
	t.Helper()
	if !generated || n <= 64 {
		return
	}
	w := want()
	if w == "" {
		t.Fatalf("%s: generated input of %d values has no in_sha256", name, n)
	}
	if g := got(); g != w {
		t.Fatalf("%s: regenerated input sha256 %s, the harness used %s (the Go generator does not mirror sref.cc)", name, g, w)
	}
}

func loadSirilGolden(t *testing.T) sirilGolden {
	t.Helper()
	data, err := os.ReadFile(sirilGoldenPath)
	if err != nil {
		t.Fatalf("Siril golden %s missing (committed file; regenerate with `make -C testdata/sref goldens`): %v", sirilGoldenPath, err)
	}
	var g sirilGolden
	if err := json.Unmarshal(data, &g); err != nil {
		t.Fatal(err)
	}
	if len(g.F32) == 0 || len(g.F64) == 0 || len(g.Percentiles) == 0 {
		t.Fatalf("Siril golden %s has no cases", sirilGoldenPath)
	}
	checkRequiredSirilCases(t, g)
	return g
}

// sirilF32CaseCount is the number of []float32 golden inputs; siril.go's
// MeanStdevSiril doc cites it ("424 inputs from n=0 to n=2*10^7").
const sirilF32CaseCount = 424

// requiredSirilF32 and requiredSirilF64 are the golden cases that specific
// claims depend on (the inputs above 2^24 elements that kill a float32
// count, the degenerate sizes, the element types). A regenerated golden
// that drops one must fail here, not silently test less: the reproduce
// test only compares the generator with the committed file.
var requiredSirilF32 = map[string]int{
	"gradient_200": 200, "gradient_1000": 1000, "random_500": 500, "random_10000": 10000,
	"with_zeros": 300, "all_same": 100, "single": 1, "empty": 0,
	"n2": 2, "n3": 3, "n8": 8, "n9": 9, "n10": 10,
	"astro_u16_5000": 5000, "negative_301": 301, "ties_1001": 1001,
	"with_nan_200": 200, "mixed_zero_nan": 11, "all_zero": 7, "adu_offset_1000": 1000,
	"f_astro_65537": 65537, "f_uniform_16777219": 16777219,
	"f_uniform_20000000": 20000000, "f_astro_20000000": 20000000,
}

var requiredSirilF64 = map[string]struct {
	dtype string
	n     int
}{
	"int16_random_500":           {"int16", 500},
	"uint8_random_500":           {"uint8", 500},
	"int32_negative_200":         {"int32", 200},
	"uint16_5000":                {"uint16", 5000},
	"int64_100":                  {"int64", 100},
	"float64_random_1000":        {"float64", 1000},
	"float64_offset_777":         {"float64", 777},
	"float64_n2":                 {"float64", 2},
	"float64_n3":                 {"float64", 3},
	"float64_single":             {"float64", 1},
	"float64_all_same":           {"float64", 50},
	"f_float64_uniform_70001":    {"float64", 70001},
	"f_float64_uniform_20000000": {"float64", 20000000},
	"f_uint16_adu_20000000":      {"uint16", 20000000},
	"int8_3000":                  {"int8", 3000},
	"uint32_3000":                {"uint32", 3000},
	"uint64_3000":                {"uint64", 3000},
	"int64_wide_3000":            {"int64", 3000},
	"namedf32_5000":              {"namedf32", 5000},
}

func checkRequiredSirilCases(t *testing.T, g sirilGolden) {
	t.Helper()
	if len(g.F32) != sirilF32CaseCount {
		t.Fatalf("Siril golden has %d float32 cases; siril.go documents %d", len(g.F32), sirilF32CaseCount)
	}
	small := 0
	f32 := map[string]sirilF32Case{}
	for _, c := range g.F32 {
		f32[c.Name] = c
		if strings.HasPrefix(c.Name, "small_") {
			small++
		}
	}
	if small != 400 {
		t.Errorf("Siril golden has %d small_* cases, want 400", small)
	}
	for name, n := range requiredSirilF32 {
		c, ok := f32[name]
		if !ok {
			t.Errorf("required Siril golden case %q missing", name)
		} else if c.N != n {
			t.Errorf("Siril golden case %q has n=%d, want %d", name, c.N, n)
		}
	}
	f64 := map[string]sirilF64Case{}
	for _, c := range g.F64 {
		f64[c.Name] = c
	}
	for name, want := range requiredSirilF64 {
		c, ok := f64[name]
		if !ok {
			t.Errorf("required Siril golden case %q missing", name)
		} else if c.Dtype != want.dtype || c.N != want.n {
			t.Errorf("Siril golden case %q is %s n=%d, want %s n=%d", name, c.Dtype, c.N, want.dtype, want.n)
		}
	}
	if t.Failed() {
		t.FailNow()
	}
}

// xorshift32, formulaData and formulaData64 mirror formula() and
// formula64() in testdata/sref/sref.cc.
// Only exactly rounded IEEE operations are involved, so Go and C produce
// the same float32 values.
func xorshift32(x *uint32) uint32 {
	*x ^= *x << 13
	*x ^= *x >> 17
	*x ^= *x << 5
	return *x
}

func formulaData(t *testing.T, f formulaSpec, n int) []float32 {
	t.Helper()
	out, err := formulaValues(f, n)
	if err != nil {
		t.Fatal(err)
	}
	return out
}

// formulaValues is formulaData without a *testing.T (the stats golden
// loader also uses it).
func formulaValues(f formulaSpec, n int) ([]float32, error) {
	out := make([]float32, n)
	x := f.Seed
	for i := range out {
		b := xorshift32(&x)
		switch f.Kind {
		case "uniform24":
			out[i] = float32(b>>8) * 5.9604644775390625e-8
		case "astro_u16":
			adu := 2900 + (b & 0xff) + ((b >> 8) & 0xff) + ((b >> 16) & 0xff) + (b >> 24)
			if b%509 == 0 {
				adu += (b >> 7) % 60000
			}
			if adu > 65535 {
				adu = 65535
			}
			out[i] = float32(adu) / 65535
		case "sky_hot", "adu_f32":
			// gen_golden.py formula_f32 only (not in sref.cc).
			k := int32((b&0xff)+((b>>8)&0xff)+((b>>16)&0xff)+(b>>24)) - 510
			switch {
			case f.Kind == "adu_f32":
				out[i] = float32(30000 + float32(float32(k)*0x1p-8))
			case b%100 == 7:
				out[i] = float32(float32(0.9) + float32(float32((b>>16)&0xfff)*0x1p-16))
			default:
				out[i] = float32(float32(0.05) + float32(float32(k)*0x1p-17))
			}
		default:
			return nil, fmt.Errorf("unknown formula kind %q", f.Kind)
		}
	}
	return out, nil
}

// formulaData64 mirrors formula64() in sref.cc: integer-valued (or
// integer times a power of two) doubles for the generic path.
func formulaData64(t *testing.T, f formulaSpec, n int) []float64 {
	t.Helper()
	out, err := formulaValues64(f, n)
	if err != nil {
		t.Fatal(err)
	}
	return out
}

// formulaValues64 is formulaData64 without a *testing.T. The kinds
// "uniform53_1e9" (1e9 + a uniform53 value, rounded once) and "wide_i32"
// are not in sref.cc: they are stats golden inputs, mirrored by
// formula_f64 in testdata/gen_golden.py (as are the float32 kinds
// "sky_hot" and "adu_f32" by formula_f32).
func formulaValues64(f formulaSpec, n int) ([]float64, error) {
	out := make([]float64, n)
	if f.Kind == "uniform24" {
		d, err := formulaValues(f, n)
		if err != nil {
			return nil, err
		}
		for i, v := range d {
			out[i] = float64(v)
		}
		return out, nil
	}
	x := f.Seed
	for i := range out {
		b := xorshift32(&x)
		switch f.Kind {
		case "astro_adu":
			adu := 2900 + (b & 0xff) + ((b >> 8) & 0xff) + ((b >> 16) & 0xff) + (b >> 24)
			if b%509 == 0 {
				adu += (b >> 7) % 60000
			}
			if adu > 65535 {
				adu = 65535
			}
			out[i] = float64(adu)
		case "u16":
			out[i] = float64(b >> 16)
		case "i8":
			v := int(b & 0xff)
			if v >= 128 {
				v -= 256
			}
			out[i] = float64(v)
		case "u32":
			out[i] = float64(b)
		case "u64":
			var v uint64
			if b&1 != 0 {
				v = 1<<63 + uint64(b)<<11
			} else {
				v = uint64(b) << 20
			}
			out[i] = float64(v)
		case "i64":
			out[i] = float64(int64(int32(b)) << 30)
		case "uniform53":
			b2 := xorshift32(&x)
			m := (uint64(b)<<32 | uint64(b2)) >> 11
			out[i] = float64(m) * 0x1p-53
		case "wide_i32":
			// gen_golden.py formula_f64 only: int32(b) * 2^467, exact.
			out[i] = float64(int32(b)) * 0x1p467
		case "uniform53_1e9":
			b2 := xorshift32(&x)
			m := (uint64(b)<<32 | uint64(b2)) >> 11
			out[i] = 1e9 + float64(m)*0x1p-53
		default:
			return nil, fmt.Errorf("unknown formula64 kind %q", f.Kind)
		}
	}
	return out, nil
}

func parseHex32List(t *testing.T, s string) []float32 {
	t.Helper()
	fields := strings.Fields(s)
	out := make([]float32, len(fields))
	for i, h := range fields {
		u, err := strconv.ParseUint(h, 16, 32)
		if err != nil {
			t.Fatalf("bad hex %q: %v", h, err)
		}
		out[i] = math.Float32frombits(uint32(u))
	}
	return out
}

func hex32(t *testing.T, s string) float32 {
	t.Helper()
	u, err := strconv.ParseUint(s, 16, 32)
	if err != nil {
		t.Fatalf("bad hex %q: %v", s, err)
	}
	return math.Float32frombits(uint32(u))
}

func hex64(t *testing.T, s string) float64 {
	t.Helper()
	u, err := strconv.ParseUint(s, 16, 64)
	if err != nil {
		t.Fatalf("bad hex %q: %v", s, err)
	}
	return math.Float64frombits(u)
}

func (c sirilF32Case) data(t *testing.T) []float32 {
	t.Helper()
	switch {
	case c.Input != nil:
		d := parseHex32List(t, *c.Input)
		if len(d) != c.N {
			t.Fatalf("%s: input has %d values, n=%d", c.Name, len(d), c.N)
		}
		return d
	case c.Formula != nil, c.Gen != nil:
		var d []float32
		if c.Formula != nil {
			d = formulaData(t, *c.Formula, c.N)
		} else {
			d = genValues32(t, *c.Gen, c.N)
		}
		checkInputDigest(t, c.Name, true, c.N, func() string { return c.InSHA256 }, func() string { return sha32(d) })
		return d
	}
	t.Fatalf("%s: no input, formula or gen", c.Name)
	return nil
}

// same32 reports whether a and b are the same float32 (bit-identical, or
// both NaN: NaN payloads are not part of any contract).
func same32(a, b float32) bool {
	if a != a && b != b {
		return true
	}
	return math.Float32bits(a) == math.Float32bits(b)
}

func same64(a, b float64) bool {
	if a != a && b != b {
		return true
	}
	return math.Float64bits(a) == math.Float64bits(b)
}

// checkF32Result asserts that a float64 result returned for []float32 input
// is exactly the float32 Siril produced (no extra precision, no rounding
// difference).
func checkF32Result(t *testing.T, what string, got float64, want float32) {
	t.Helper()
	if float64(float32(got)) != got && !math.IsNaN(got) {
		t.Errorf("%s: result %.17g is not a float32 value", what, got)
	}
	if !same32(float32(got), want) {
		t.Errorf("%s: got %.9g (%08x) want Siril %.9g (%08x)", what,
			float32(got), math.Float32bits(float32(got)), want, math.Float32bits(want))
	}
}

// TestSirilOracle_Float32 checks every []float32 entry point that claims
// Siril equivalence against the verbatim Siril code, bit for bit.
func TestSirilOracle_Float32(t *testing.T) {
	g := loadSirilGolden(t)
	histo := make([]uint32, MaxHistoSize)
	sdMismatchVsGeneric := 0
	for _, c := range g.F32 {
		t.Run(c.Name, func(t *testing.T) {
			raw := c.data(t)

			t.Run("MeanStdevSiril", func(t *testing.T) {
				m, s := MeanStdevSiril(raw)
				wm, ws := hex32(t, c.SD.Mean), hex32(t, c.SD.Stdev)
				if !same32(m, wm) || !same32(s, ws) {
					t.Errorf("got (%.9g %08x, %.9g %08x) want Siril (%.9g %08x, %.9g %08x)",
						m, math.Float32bits(m), s, math.Float32bits(s),
						wm, math.Float32bits(wm), ws, math.Float32bits(ws))
				}
				if len(raw) >= 2 {
					if _, gs := MeanStdev(raw); float32(gs) != s {
						sdMismatchVsGeneric++
					}
				}
			})

			if c.NonnullN != nil {
				t.Run("FilterNonZero", func(t *testing.T) {
					got := FilterNonZero(raw)
					if *c.NonnullN == 0 {
						if len(got) != len(raw) || (len(raw) > 0 && &got[0] != &raw[0]) {
							t.Errorf("all-zero/NaN input: want the input slice back, got %d values", len(got))
						}
						return
					}
					if len(got) != *c.NonnullN || sha32(got) != c.NonnullSHA256 {
						t.Errorf("%d values, sha256 %s; want Siril's %d values, sha256 %s", len(got), sha32(got), *c.NonnullN, c.NonnullSHA256)
					}
				})
			} else if c.Input == nil && c.N > 64 {
				t.Errorf("generated input of %d values has no FilterNonZero digest", c.N)
			}
			if c.Nonnull != nil {
				t.Run("FilterNonZero", func(t *testing.T) {
					want := parseHex32List(t, *c.Nonnull)
					got := FilterNonZero(raw)
					if len(want) == 0 {
						// Siril computes no statistics here; FilterNonZero's
						// documented fallback is the input slice itself.
						if len(got) != len(raw) || (len(raw) > 0 && &got[0] != &raw[0]) {
							t.Errorf("all-zero/NaN input: want the input slice back, got %v", got)
						}
						return
					}
					if len(got) != len(want) {
						t.Fatalf("len %d want %d", len(got), len(want))
					}
					for i := range want {
						if math.Float32bits(got[i]) != math.Float32bits(want[i]) {
							t.Errorf("[%d] = %08x want %08x", i, math.Float32bits(got[i]), math.Float32bits(want[i]))
						}
					}
				})
			}

			if c.NClean == 0 {
				return
			}
			t.Run("Percentile", func(t *testing.T) {
				for _, ps := range g.Percentiles {
					p, err := strconv.ParseFloat(ps, 64)
					if err != nil {
						t.Fatal(err)
					}
					want := hex32(t, c.Pct[ps])
					checkF32Result(t, "Percentile p="+ps, Percentile(raw, p), want)
					checkF32Result(t, "PercentileBuf p="+ps, PercentileBuf(raw, p, histo), want)
				}
			})
			t.Run("Median", func(t *testing.T) {
				want := hex32(t, *c.Median)
				checkF32Result(t, "Median", Median(raw), want)
				checkF32Result(t, "MedianBuf", MedianBuf(raw, histo), want)
			})
			t.Run("MAD", func(t *testing.T) {
				wantMed, wantMad := hex32(t, *c.Median), hex32(t, *c.Mad)
				checkF32Result(t, "MAD", MAD(raw), wantMad)
				med, mad := MADWithMedian(raw)
				checkF32Result(t, "MADWithMedian median", med, wantMed)
				checkF32Result(t, "MADWithMedian mad", mad, wantMad)
				absDev := make([]float64, 0, len(raw))
				med, mad = MADWithMedianBuf(raw, histo, absDev)
				checkF32Result(t, "MADWithMedianBuf median", med, wantMed)
				checkF32Result(t, "MADWithMedianBuf mad", mad, wantMad)
			})
			t.Run("QuickSelect_vs_quickmedian_float", func(t *testing.T) {
				clean := make([]float32, 0, len(raw))
				for _, v := range raw {
					if v == v {
						clean = append(clean, v)
					}
				}
				n := len(clean)
				want := hex64(t, *c.Quickmedian)
				work := slices.Clone(clean)
				hi := QuickSelect(work, n/2)
				var got float64
				switch {
				case n%2 == 1:
					got = float64(hi)
				case n < 9:
					// sortnet_median_float: float sum, then /2.0 in double.
					lo := QuickSelect(work, n/2-1)
					got = float64(lo+hi) / 2
				default:
					// quickmedian_float: ((double)a + b) / 2.0.
					lo := QuickSelect(work, n/2-1)
					got = (float64(lo) + float64(hi)) / 2
				}
				if !same64(got, want) {
					t.Errorf("median from QuickSelect %.17g want quickmedian_float %.17g", got, want)
				}
			})
		})
	}
	// MeanStdevSiril's doc says it differs from MeanStdev at the ULP level
	// on ordinary inputs; the fixture set must show that.
	if sdMismatchVsGeneric == 0 {
		t.Errorf("MeanStdevSiril equalled float32(MeanStdev) on every case; the documented ULP-level difference is not exercised")
	}
	t.Logf("MeanStdevSiril != float32(MeanStdev) on %d of %d cases", sdMismatchVsGeneric, len(g.F32))
}

// namedF32 is a named type whose underlying type is float32. The package
// documents that such types take the float64 path.
type namedF32 float32

func convertSlice[T Numeric](t *testing.T, d []float64) []T {
	t.Helper()
	out := make([]T, len(d))
	for i, v := range d {
		out[i] = T(v)
		if float64(out[i]) != v {
			t.Fatalf("golden value %v not representable in the target type", v)
		}
	}
	return out
}

func genericCheck[T Numeric](t *testing.T, data []T, c sirilF64Case, ps []string) {
	t.Helper()
	histo := make([]uint32, MaxHistoSize)
	for _, s := range ps {
		p, err := strconv.ParseFloat(s, 64)
		if err != nil {
			t.Fatal(err)
		}
		want := hex64(t, c.Pct[s])
		if got := Percentile(data, p); !same64(got, want) {
			t.Errorf("Percentile p=%s: got %.17g want %.17g", s, got, want)
		}
		if got := PercentileBuf(data, p, histo); !same64(got, want) {
			t.Errorf("PercentileBuf p=%s: got %.17g want %.17g", s, got, want)
		}
	}
	wantMed, wantMad := hex64(t, c.Median), hex64(t, c.Mad)
	if got := Median(data); !same64(got, wantMed) {
		t.Errorf("Median: got %.17g want %.17g", got, wantMed)
	}
	if got := MAD(data); !same64(got, wantMad) {
		t.Errorf("MAD: got %.17g want %.17g", got, wantMad)
	}
	med, mad := MADWithMedianBuf(data, histo, make([]float64, 0, len(data)))
	if !same64(med, wantMed) || !same64(mad, wantMad) {
		t.Errorf("MADWithMedianBuf: got (%.17g, %.17g) want (%.17g, %.17g)", med, mad, wantMed, wantMad)
	}
}

// TestSirilOracle_Generic checks the float64 path (every element type
// except float32) against findMinMaxPercentile compiled at double precision.
func TestSirilOracle_Generic(t *testing.T) {
	g := loadSirilGolden(t)
	for _, c := range g.F64 {
		t.Run(c.Name, func(t *testing.T) {
			var d []float64
			switch {
			case c.Input != nil:
				for _, h := range strings.Fields(*c.Input) {
					d = append(d, hex64(t, h))
				}
			case c.Formula != nil:
				d = formulaData64(t, *c.Formula, c.N)
			case c.Gen != nil:
				d = genValues64(t, *c.Gen, c.N)
			default:
				t.Fatal("no input, formula or gen")
			}
			checkInputDigest(t, c.Name, c.Input == nil, len(d), func() string { return c.InSHA256 }, func() string { return sha64(d) })
			if len(d) != c.N {
				t.Fatalf("input has %d values, n=%d", len(d), c.N)
			}
			switch c.Dtype {
			case "float64":
				genericCheck(t, d, c, g.Percentiles)
			case "int16":
				genericCheck(t, convertSlice[int16](t, d), c, g.Percentiles)
			case "uint8":
				genericCheck(t, convertSlice[uint8](t, d), c, g.Percentiles)
			case "int32":
				genericCheck(t, convertSlice[int32](t, d), c, g.Percentiles)
			case "uint16":
				genericCheck(t, convertSlice[uint16](t, d), c, g.Percentiles)
			case "int64":
				genericCheck(t, convertSlice[int64](t, d), c, g.Percentiles)
			case "int8":
				genericCheck(t, convertSlice[int8](t, d), c, g.Percentiles)
			case "uint32":
				genericCheck(t, convertSlice[uint32](t, d), c, g.Percentiles)
			case "uint64":
				genericCheck(t, convertSlice[uint64](t, d), c, g.Percentiles)
			case "namedf32":
				// Underlying type float32, but not []float32: the
				// float64 path, compared with the float64 oracle.
				genericCheck(t, convertSlice[namedF32](t, d), c, g.Percentiles)
			default:
				t.Fatalf("unknown dtype %q", c.Dtype)
			}
		})
	}
}

// TestMeanStdevSiril_Edges pins the degenerate-input behaviour the doc
// lists, and checks that the Siril harness agrees for the "empty", "single"
// and "with_nan_200" golden cases.
func TestMeanStdevSiril_Edges(t *testing.T) {
	g := loadSirilGolden(t)
	byName := map[string]sirilF32Case{}
	for _, c := range g.F32 {
		byName[c.Name] = c
	}
	for _, name := range []string{"empty", "single", "with_nan_200"} {
		if _, ok := byName[name]; !ok {
			t.Fatalf("golden case %q missing", name)
		}
	}
	m, s := MeanStdevSiril(nil)
	if !math.IsNaN(float64(m)) || math.Float32bits(s) != 0x80000000 {
		t.Errorf("n=0: got (%v, %v [%08x]) want (NaN, -0)", m, s, math.Float32bits(s))
	}
	if ws := hex32(t, byName["empty"].SD.Stdev); math.Float32bits(ws) != 0x80000000 {
		t.Errorf("Siril n=0 stdev is %08x, doc says -0", math.Float32bits(ws))
	}
	m, s = MeanStdevSiril([]float32{3.14})
	if m != 3.14 || !math.IsNaN(float64(s)) {
		t.Errorf("n=1: got (%v, %v) want (3.14, NaN)", m, s)
	}
	if ws := hex32(t, byName["single"].SD.Stdev); !math.IsNaN(float64(ws)) {
		t.Errorf("Siril n=1 stdev is %v, doc says NaN", ws)
	}
	if ws := hex32(t, byName["with_nan_200"].SD.Stdev); !math.IsNaN(float64(ws)) {
		t.Errorf("Siril stdev with NaN input is %v, doc says NaN", ws)
	}
}
