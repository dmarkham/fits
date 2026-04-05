package transform

import (
	"testing"

	"github.com/dmarkham/fits/wcs"
)

// makeTestTransform builds a Transform for the given projection code for
// benchmark use. Uses TAN defaults for the linear pipeline.
func makeTestTransform(code string, pv map[wcs.PVKey]float64) (*Transform, error) {
	w := &wcs.Header{
		NAxis:    2,
		CType:    []string{"RA---" + code, "DEC--" + code},
		CRVal:    []float64{0, 0},
		CRPix:    []float64{200, 200},
		CDelt:    []float64{-0.01, 0.01},
		PC:       [][]float64{{1, 0}, {0, 1}},
		PV:       pv,
		LonAxis:  1,
		LatAxis:  2,
		AxisType: []string{"RA", "DEC"},
		ProjCode: []string{code, code},
	}
	return New(w)
}

// BenchmarkPixelToSkyTAN measures single-call pixel→sky throughput for
// the TAN projection (the workhorse for optical astronomy). Target: sub-
// microsecond per call, which is within 2× of wcslib.
func BenchmarkPixelToSkyTAN(b *testing.B) {
	tr, err := makeTestTransform("TAN", nil)
	if err != nil {
		b.Fatal(err)
	}
	b.ResetTimer()
	for b.Loop() {
		_, _, _ = tr.PixelToSky(250, 250)
	}
}

// BenchmarkSkyToPixelTAN is the inverse direction.
func BenchmarkSkyToPixelTAN(b *testing.B) {
	tr, err := makeTestTransform("TAN", nil)
	if err != nil {
		b.Fatal(err)
	}
	b.ResetTimer()
	for b.Loop() {
		_, _, _ = tr.SkyToPixel(0.5, 0.5)
	}
}

// BenchmarkPixelToSky1Mpixel performs 1M forward transforms to measure
// bulk throughput. Target: ~1e7 calls/sec on modern hardware, giving
// ~100ns per call.
func BenchmarkPixelToSky1Mpixel(b *testing.B) {
	tr, _ := makeTestTransform("TAN", nil)
	b.ResetTimer()
	b.SetBytes(int64(1_000_000 * 16)) // 16 bytes per pixel (2 × float64)
	for b.Loop() {
		for i := 0; i < 1000; i++ {
			for j := 0; j < 1000; j++ {
				_, _, _ = tr.PixelToSky(float64(i), float64(j))
			}
		}
	}
}

// BenchmarkHPXForward measures the HEALPix forward path. HPX has a
// piecewise polar-cap fold so benchmarks against it stress the common
// all-sky projection case for CMB-era data products.
func BenchmarkHPXForward(b *testing.B) {
	tr, _ := makeTestTransform("HPX", nil)
	b.ResetTimer()
	for b.Loop() {
		_, _, _ = tr.PixelToSky(200, 200)
	}
}

// BenchmarkMOLInverse measures Mollweide, which uses Newton iteration
// in its forward path — this is the most expensive common projection.
func BenchmarkMOLForward(b *testing.B) {
	tr, _ := makeTestTransform("MOL", nil)
	b.ResetTimer()
	for b.Loop() {
		_, _, _ = tr.SkyToPixel(30, 20)
	}
}
