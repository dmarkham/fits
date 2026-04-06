package stats

import (
	"encoding/json"
	"math"
	"os"
	"testing"
)

// Tests in this file validate against Siril's actual findMinMaxPercentile
// implementation (RawTherapee rt_algo.cc), compiled in a standalone C++
// harness at testdata/sref/. This is the authoritative oracle — Siril
// is the processor we're matching.
//
// Regenerate goldens:  make -C testdata/sref clean && make -C testdata/sref goldens

const sirilGoldenPath = "testdata/sref/golden/siril_stats_golden.json"

type sirilCase struct {
	Len         int                `json:"len"`
	Median      *float64           `json:"median"`
	Mad         *float64           `json:"mad"`
	Percentiles map[string]float64 `json:"percentiles"`
	SigmaClip   map[string]struct {
		Mean   float64 `json:"mean"`
		Stdev  float64 `json:"stdev"`
		Median float64 `json:"median"`
		NGood  int     `json:"ngood"`
	} `json:"sigma_clip"`
}

func loadSirilGoldens(t *testing.T) map[string]sirilCase {
	t.Helper()
	data, err := os.ReadFile(sirilGoldenPath)
	if err != nil {
		t.Skipf("Siril golden not found — run `make -C testdata/sref goldens`")
	}
	var cases map[string]sirilCase
	if err := json.Unmarshal(data, &cases); err != nil {
		t.Fatal(err)
	}
	return cases
}

// generateTestData produces the same deterministic arrays as sref.cc.
// Must match exactly — same formula, same seed, same values.
func generateTestData(name string) []float32 {
	switch name {
	case "gradient_200":
		return genGradient(200, 0.01, 0.30)
	case "gradient_1000":
		return genGradient(1000, 0.0, 1.0)
	case "random_500", "random_10000":
		return nil // C++ mt19937, can't reproduce in Go — skipped
	case "with_zeros":
		d := genGradient(300, 0.01, 0.25)
		for i := 0; i < 60; i++ {
			d[i] = 0
		}
		return d
	case "all_same":
		d := make([]float32, 100)
		for i := range d {
			d[i] = 0.42
		}
		return d
	case "single":
		return []float32{3.14}
	}
	return nil
}

func genGradient(n int, lo, hi float64) []float32 {
	// Must match sref.cc's gen_gradient: float32 arithmetic throughout.
	loF := float32(lo)
	rangeF := float32(hi) - loF
	out := make([]float32, n)
	for i := range out {
		out[i] = loF + rangeF*float32(i)/float32(n-1)
	}
	return out
}

// Note: the "random_500" and "random_10000" cases use C++ mt19937
// which we can't reproduce in Go. Those cases are skipped — the
// gradient and deterministic cases provide the core algorithm
// validation against Siril's findMinMaxPercentile.

func TestMedian_Siril(t *testing.T) {
	cases := loadSirilGoldens(t)
	for name, c := range cases {
		if c.Median == nil {
			continue
		}
		data := generateTestData(name)
		if data == nil {
			continue // random data we can't reproduce
		}
		t.Run(name, func(t *testing.T) {
			got := Median(data)
			if !sirilApprox(got, *c.Median, 1e-6) {
				t.Errorf("median: got %.17g want %.17g (diff %g)",
					got, *c.Median, got-*c.Median)
			}
		})
	}
}

func TestMAD_Siril(t *testing.T) {
	cases := loadSirilGoldens(t)
	for name, c := range cases {
		if c.Mad == nil {
			continue
		}
		data := generateTestData(name)
		if data == nil {
			continue
		}
		t.Run(name, func(t *testing.T) {
			got := MAD(data)
			if !sirilApprox(got, *c.Mad, 1e-6) {
				t.Errorf("mad: got %.17g want %.17g (diff %g)",
					got, *c.Mad, got-*c.Mad)
			}
		})
	}
}

func TestPercentile_Siril(t *testing.T) {
	cases := loadSirilGoldens(t)
	for name, c := range cases {
		if len(c.Percentiles) == 0 {
			continue
		}
		data := generateTestData(name)
		if data == nil {
			continue
		}
		t.Run(name, func(t *testing.T) {
			for pStr, want := range c.Percentiles {
				var p float64
				switch pStr {
				case "0.01":
					p = 0.01
				case "0.1":
					p = 0.1
				case "0.25":
					p = 0.25
				case "0.5":
					p = 0.5
				case "0.75":
					p = 0.75
				case "0.9":
					p = 0.9
				case "0.99":
					p = 0.99
				default:
					continue
				}
				got := Percentile(data, p)
				if !sirilApprox(got, want, 1e-6) {
					t.Errorf("p=%s: got %.17g want %.17g (diff %g)",
						pStr, got, want, got-want)
				}
			}
		})
	}
}

func TestSigmaClip_Siril(t *testing.T) {
	cases := loadSirilGoldens(t)
	for name, c := range cases {
		if len(c.SigmaClip) == 0 {
			continue
		}
		data := generateTestData(name)
		if data == nil {
			continue
		}
		t.Run(name, func(t *testing.T) {
			for scKey, expected := range c.SigmaClip {
				t.Run(scKey, func(t *testing.T) {
					var k float64
					var cf CenterFunc
					switch scKey {
					case "k2.0_mean":
						k, cf = 2, CenterMean
					case "k3.0_mean":
						k, cf = 3, CenterMean
					case "k5.0_mean":
						k, cf = 5, CenterMean
					case "k2.0_median":
						k, cf = 2, CenterMedian
					case "k3.0_median":
						k, cf = 3, CenterMedian
					case "k5.0_median":
						k, cf = 5, CenterMedian
					default:
						t.Skip("unknown key")
					}
					got := SigmaClip(data, k, k, 10, cf)

					if !sirilApprox(got.Mean, expected.Mean, 1e-6) {
						t.Errorf("mean: got %.17g want %.17g", got.Mean, expected.Mean)
					}
					if !sirilApprox(got.Stdev, expected.Stdev, 1e-6) {
						t.Errorf("stdev: got %.17g want %.17g", got.Stdev, expected.Stdev)
					}
					if !sirilApprox(got.Median, expected.Median, 1e-6) {
						t.Errorf("median: got %.17g want %.17g", got.Median, expected.Median)
					}
					if got.NGood != expected.NGood {
						t.Errorf("ngood: got %d want %d", got.NGood, expected.NGood)
					}
				})
			}
		})
	}
}

// sirilApprox compares with a tight tolerance — we're matching the
// SAME algorithm (findMinMaxPercentile), so the only difference is
// Go float64 vs C++ float arithmetic ordering.
func sirilApprox(a, b, tol float64) bool {
	if a == b {
		return true
	}
	diff := math.Abs(a - b)
	if diff < 1e-12 {
		return true
	}
	scale := math.Max(math.Abs(a), math.Abs(b))
	if scale == 0 {
		return diff < tol
	}
	return diff/scale < tol
}
