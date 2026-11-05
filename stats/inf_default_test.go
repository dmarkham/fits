package stats

import (
	"math"
	"reflect"
	"testing"
)

// TestSigmaClipDropsInfOnGenericPath: on the float64 (generic) path ±Inf
// is removed before clipping, as on the float32 path, so adding infinities
// to the input leaves the result unchanged.
func TestSigmaClipDropsInfOnGenericPath(t *testing.T) {
	finite := make([]float64, 0, 64)
	for i := 0; i < 60; i++ {
		finite = append(finite, float64(i%17)*0.25+1)
	}
	finite = append(finite, 40) // one outlier to clip
	withInf := append(append([]float64{math.Inf(1)}, finite...), math.Inf(-1), math.Inf(1))
	for _, c := range []struct {
		name   string
		center CenterFunc
	}{{"median", CenterMedian}, {"mean", CenterMean}} {
		want := SigmaClip(finite, 3, 3, 10, c.center)
		got := SigmaClip(withInf, 3, 3, 10, c.center)
		if want.NGood == 0 {
			t.Fatalf("%s: reference clip kept nothing", c.name)
		}
		if !reflect.DeepEqual(got, want) {
			t.Errorf("%s: with ±Inf %+v, without %+v", c.name, got, want)
		}
	}
}

// TestBuildHistogramDefaultBins: nbins <= 0 means 256 bins.
func TestBuildHistogramDefaultBins(t *testing.T) {
	data := []float64{0, 1, 2, 3, 4, 5, 6, 7, 8, 9.5}
	want, err := BuildHistogram(data, 256)
	if err != nil {
		t.Fatal(err)
	}
	for _, n := range []int{0, -1} {
		got, err := BuildHistogram(data, n)
		if err != nil {
			t.Fatal(err)
		}
		if got.NBins != 256 || !reflect.DeepEqual(got, want) {
			t.Errorf("nbins=%d: NBins %d, equal to 256-bin histogram: %v", n, got.NBins, reflect.DeepEqual(got, want))
		}
	}
}
