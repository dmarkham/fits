package healpix

import (
	"math"
	"testing"
)

// benchDirs returns 1024 fixed directions spread over the sphere (a
// Fibonacci lattice), as (theta, phi) and as unit vectors, so the
// benchmarks exercise the equatorial belt and both polar caps.
func benchDirs() (theta, phi, vx, vy, vz []float64) {
	const n = 1024
	golden := math.Pi * (3 - math.Sqrt(5))
	for i := 0; i < n; i++ {
		z := 1 - (2*float64(i)+1)/n
		th := math.Acos(z)
		ph := math.Mod(float64(i)*golden, 2*math.Pi)
		s := math.Sin(th)
		theta = append(theta, th)
		phi = append(phi, ph)
		vx = append(vx, float64(s*math.Cos(ph)))
		vy = append(vy, float64(s*math.Sin(ph)))
		vz = append(vz, z)
	}
	return
}

var benchSink int64

// BenchmarkAng2PixRing measures Ang2Pix at Nside 1024, RING ordering.
func BenchmarkAng2PixRing(b *testing.B) {
	theta, phi, _, _, _ := benchDirs()
	i := 0
	for b.Loop() {
		benchSink += Ang2Pix(1024, theta[i], phi[i], Ring)
		i = (i + 1) & 1023
	}
}

// BenchmarkVec2PixNested measures Vec2Pix at Nside 1024, NESTED ordering.
func BenchmarkVec2PixNested(b *testing.B) {
	_, _, vx, vy, vz := benchDirs()
	i := 0
	for b.Loop() {
		benchSink += Vec2Pix(1024, vx[i], vy[i], vz[i], Nested)
		i = (i + 1) & 1023
	}
}
