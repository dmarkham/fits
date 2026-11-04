package stats

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"math"
	"math/big"
	"os"
	"strconv"
	"testing"
)

const meanRoundingGoldenPath = "testdata/golden/mean_rounding.json"

type meanRoundingGolden struct {
	Seed    uint32            `json:"seed"`
	Sets    int               `json:"sets"`
	Inexact int               `json:"inexact"`
	SHA256  string            `json:"sha256_be_f64"`
	Samples map[string]string `json:"samples"`
}

// meanRoundingSets mirrors mean_rounding_sets() in testdata/gen_golden.py:
// sets of 2..64 integers ±(2^e + k), e in 40..62, |k| < 2^23.
func meanRoundingSets(seed uint32, sets int) [][]int64 {
	out := make([][]int64, 0, sets)
	x := seed
	for range sets {
		n := 2 + int(xorshift32(&x)%63)
		e := 40 + xorshift32(&x)%23
		neg := xorshift32(&x)&1 != 0
		vals := make([]int64, n)
		for i := range vals {
			v := int64(1)<<e + int64(xorshift32(&x)>>8) - 1<<23
			if neg {
				v = -v
			}
			vals[i] = v
		}
		out = append(out, vals)
	}
	return out
}

// TestMeanStdevTwoPass_MeanCorrectlyRounded proves the MeanStdevTwoPass
// doc's claim that the mean is correctly rounded: on 20000 sets of large
// integers (beyond 2^53 for most, so the float64 conversion rounds and
// the sum carries far more bits than a float64) the mean, for []int64 and
// for the same values as []float64, is bit-identical to the exact mean of
// the float64 values correctly rounded (Python exact integer arithmetic,
// testdata/gen_golden.py). The inputs are regenerated here; the golden
// holds a sha256 of every expected mean (big-endian float64 bits) plus
// explicit samples. On a mismatch math/big locates the first wrong set
// (diagnostic only; the truth is the golden).
func TestMeanStdevTwoPass_MeanCorrectlyRounded(t *testing.T) {
	b, err := os.ReadFile(meanRoundingGoldenPath)
	if err != nil {
		t.Fatalf("golden %s missing (committed; regenerate with testdata/gen_golden.py): %v", meanRoundingGoldenPath, err)
	}
	var g meanRoundingGolden
	if err := json.Unmarshal(b, &g); err != nil {
		t.Fatal(err)
	}
	if g.Sets != 20000 || len(g.Samples) < 30 || g.Inexact < g.Sets/2 {
		t.Fatalf("golden has %d sets, %d samples, %d inexact means; want 20000 sets, >= 30 samples, mostly inexact", g.Sets, len(g.Samples), g.Inexact)
	}
	sets := meanRoundingSets(g.Seed, g.Sets)
	for name, conv := range map[string]func([]int64) float64{
		"int64": func(v []int64) float64 { m, _ := MeanStdevTwoPass(v); return m },
		"float64": func(v []int64) float64 {
			f := make([]float64, len(v))
			for i, x := range v {
				f[i] = float64(x)
			}
			m, _ := MeanStdevTwoPass(f)
			return m
		},
	} {
		t.Run(name, func(t *testing.T) {
			h := sha256.New()
			got := make([]float64, len(sets))
			var buf [8]byte
			for i, s := range sets {
				got[i] = conv(s)
				binary.BigEndian.PutUint64(buf[:], math.Float64bits(got[i]))
				h.Write(buf[:])
			}
			bad := 0
			for is, hx := range g.Samples {
				i, err := strconv.Atoi(is)
				if err != nil || i < 0 || i >= len(sets) {
					t.Fatalf("bad sample index %q", is)
				}
				raw, err := hex.DecodeString(hx)
				if err != nil || len(raw) != 8 {
					t.Fatalf("bad sample %q", hx)
				}
				want := math.Float64frombits(binary.BigEndian.Uint64(raw))
				if math.Float64bits(got[i]) != math.Float64bits(want) {
					bad++
					t.Errorf("set %d (n=%d): mean %.17g want %.17g (%d ULP)", i, len(sets[i]), got[i], want, ulpDiff(got[i], want))
				}
			}
			if hex.EncodeToString(h.Sum(nil)) != g.SHA256 {
				wrong, first := 0, -1
				for i, s := range sets {
					if math.Float64bits(got[i]) != math.Float64bits(exactMeanBig(s)) {
						wrong++
						if first < 0 {
							first = i
						}
					}
				}
				t.Errorf("sha256 of the %d means differs from the golden: %d of %d sets not correctly rounded (math/big), first set %d", len(sets), wrong, len(sets), first)
			}
		})
	}
}

// exactMeanBig is a diagnostic only: the exact mean of float64(v),
// correctly rounded, via math/big.
func exactMeanBig(v []int64) float64 {
	var s big.Int
	for _, x := range v {
		bf := new(big.Float).SetFloat64(float64(x))
		bi, _ := bf.Int(nil)
		s.Add(&s, bi)
	}
	f, _ := new(big.Rat).SetFrac(&s, big.NewInt(int64(len(v)))).Float64()
	return f
}
