package stats

import (
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"errors"
	"math"
	"os"
	"slices"
	"testing"
)

const histogramRangesGoldenPath = "testdata/golden/histogram_ranges.json"

type histogramRangeCase struct {
	Values     []string `json:"values"` // big-endian float64 bits, hex
	NBins      int      `json:"nbins"`
	Counts     []int64  `json:"counts"`
	Edges      []string `json:"edges"`
	NumpyError string   `json:"numpy_error"`
}

func beHex64(t *testing.T, s string) float64 {
	t.Helper()
	b, err := hex.DecodeString(s)
	if err != nil || len(b) != 8 {
		t.Fatalf("bad float64 hex %q", s)
	}
	return math.Float64frombits(binary.BigEndian.Uint64(b))
}

// TestBuildHistogram_NumpyRanges runs numpy.histogram's behaviour at the
// limits of float64 (testdata/gen_golden.py HIST_RANGES): single values
// at 2^52 and 2^53 (the ±0.5 padding kept or lost), max-min overflowing,
// subnormal ranges with one, two and three bins, and integers at 2^60
// with just enough and one too many bins. Where numpy bins, counts and
// edges must be bit-identical; where it raises (ValueError "Too many bins
// for data range", or IndexError on an overflowing range with one bin),
// BuildHistogram must return ErrHistogramRange and a zero Histogram.
func TestBuildHistogram_NumpyRanges(t *testing.T) {
	b, err := os.ReadFile(histogramRangesGoldenPath)
	if err != nil {
		t.Fatalf("golden %s missing (committed; regenerate with testdata/gen_golden.py): %v", histogramRangesGoldenPath, err)
	}
	var cases []histogramRangeCase
	if err := json.Unmarshal(b, &cases); err != nil {
		t.Fatal(err)
	}
	raised, binned := 0, 0
	for i, c := range cases {
		vals := make([]float64, len(c.Values))
		for j, h := range c.Values {
			vals[j] = beHex64(t, h)
		}
		h, err := BuildHistogram(vals, c.NBins)
		if c.NumpyError != "" {
			raised++
			if !errors.Is(err, ErrHistogramRange) || h.NBins != 0 || h.Counts != nil || h.Edges != nil {
				t.Errorf("case %d %v nbins %d: numpy raises %q; got err %v, histogram %+v", i, vals, c.NBins, c.NumpyError, err, h)
			}
			continue
		}
		binned++
		if err != nil {
			t.Errorf("case %d %v nbins %d: numpy bins it; got %v", i, vals, c.NBins, err)
			continue
		}
		if len(h.Counts) != len(c.Counts) || len(h.Edges) != len(c.Edges) {
			t.Errorf("case %d: shape %d/%d want %d/%d", i, len(h.Counts), len(h.Edges), len(c.Counts), len(c.Edges))
			continue
		}
		for j := range c.Counts {
			if h.Counts[j] != c.Counts[j] {
				t.Errorf("case %d %v nbins %d: counts %v want %v", i, vals, c.NBins, h.Counts, c.Counts)
				break
			}
		}
		for j, e := range c.Edges {
			if w := beHex64(t, e); math.Float64bits(h.Edges[j]) != math.Float64bits(w) {
				t.Errorf("case %d %v nbins %d: edge %d %.17g want %.17g", i, vals, c.NBins, j, h.Edges[j], w)
			}
		}
	}
	if raised < 8 || binned < 6 {
		t.Errorf("golden has %d raising and %d binning cases; want >= 8 / 6", raised, binned)
	}
}

const histogramEdgesGoldenPath = "testdata/golden/histogram_edges.json"

type histogramEdgeCase struct {
	Values     []string `json:"values"` // committed input (else built from Edges)
	NBins      int      `json:"nbins"`
	Edges      []string `json:"edges"`
	Counts     []int64  `json:"counts"`
	Decrements int      `json:"decrements"`
	Increments int      `json:"increments"`
}

// edgeInput is gen_golden.py edge_input: [edges[0], edges[n]], then for
// each edge e: the float64 below it (not for edges[0]), e itself (interior
// edges only) and the float64 above it (not for edges[n]).
func edgeInput(edges []float64) []float64 {
	n := len(edges) - 1
	v := []float64{edges[0], edges[n]}
	for i, e := range edges {
		if i > 0 {
			v = append(v, math.Nextafter(e, math.Inf(-1)))
		}
		if i > 0 && i < n {
			v = append(v, e)
		}
		if i < n {
			v = append(v, math.Nextafter(e, math.Inf(1)))
		}
	}
	return v
}

// TestBuildHistogram_NumpyEdgeCorrections: values on numpy's bin edges
// and one float64 either side of them (testdata/gen_golden.py
// histogram_edges_golden), where numpy's computed index is one off and
// numpy corrects it downwards (x < edges[i]) or upwards (x >= edges[i+1]
// outside the last bin). Counts and edges must be bit-identical to
// numpy.histogram; the golden records how often numpy corrected, and
// both corrections must be exercised.
func TestBuildHistogram_NumpyEdgeCorrections(t *testing.T) {
	b, err := os.ReadFile(histogramEdgesGoldenPath)
	if err != nil {
		t.Fatalf("golden %s missing (committed; regenerate with testdata/gen_golden.py): %v", histogramEdgesGoldenPath, err)
	}
	var cases []histogramEdgeCase
	if err := json.Unmarshal(b, &cases); err != nil {
		t.Fatal(err)
	}
	dec, inc := 0, 0
	for i, c := range cases {
		edges := make([]float64, len(c.Edges))
		for j, h := range c.Edges {
			edges[j] = beHex64(t, h)
		}
		nbins := len(edges) - 1
		var vals []float64
		if c.Values != nil {
			if c.NBins != nbins {
				t.Fatalf("case %d: nbins %d but %d edges", i, c.NBins, len(edges))
			}
			for _, h := range c.Values {
				vals = append(vals, beHex64(t, h))
			}
		} else {
			vals = edgeInput(edges)
		}
		dec += c.Decrements
		inc += c.Increments
		h, err := BuildHistogram(vals, nbins)
		if err != nil {
			t.Errorf("case %d: %v", i, err)
			continue
		}
		for j := range edges {
			if math.Float64bits(h.Edges[j]) != math.Float64bits(edges[j]) {
				t.Errorf("case %d: edge %d %.17g want numpy %.17g", i, j, h.Edges[j], edges[j])
			}
		}
		if !slices.Equal(h.Counts, c.Counts) {
			t.Errorf("case %d (%d bins over [%.17g, %.17g], numpy corrected %d down, %d up): counts %v want numpy %v",
				i, nbins, edges[0], edges[nbins], c.Decrements, c.Increments, h.Counts, c.Counts)
		}
	}
	if dec < 50 || inc < 50 {
		t.Fatalf("golden exercises numpy's decrement %d and increment %d times; want both >= 50", dec, inc)
	}
}
