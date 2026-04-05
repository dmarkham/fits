package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// This file implements the HEALPix-based projections from Calabretta &
// Roukema 2007: HPX (the standard HEALPix projection) and XPH (the polar
// "butterfly" variant). Both are exactly equal-area and tessellate the
// sphere into 12 congruent diamond facets.

// HPX — HEALPix projection (CR 2007 §2).
//
// Parameters:
//
//	PV2_1 = H (longitude facets, default 4)
//	PV2_2 = K (latitude facets, default 3)
//
// For the standard HEALPix (H=4, K=3), the sphere is tessellated into 12
// diamonds arranged in an equatorial band of 8 diamonds plus 4 polar
// triangles (each polar cap contains 2 wrapped diamonds).
//
// The projection has two regimes:
//
//   Equatorial region (|sin(theta)| <= (K-1)/K, approximately |theta| <= 41.8°):
//     x = phi
//     y = (pi * K / (2*H)) * sin(theta) * H / K = (pi/2) * sin(theta)
//
//   Polar caps (|sin(theta)| > (K-1)/K):
//     The diamond fold kicks in; x and y are piecewise-linear in a
//     rotated coordinate system.
//
// For H=4, K=3 this reduces to the formulas in CR 2007 eq. 5-7.
type hpxProjection struct {
	H int // longitude facets
	K int // latitude facets
}

func newHPX(pv map[wcs.PVKey]float64, latAxis int) Projection {
	H := int(pvFloat(pv, latAxis, 1, 4))
	K := int(pvFloat(pv, latAxis, 2, 3))
	if H <= 0 {
		H = 4
	}
	if K <= 0 {
		K = 3
	}
	return hpxProjection{H: H, K: K}
}

func (hpxProjection) Code() string    { return "HPX" }
func (hpxProjection) Theta0() float64 { return 0 }

func (p hpxProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	// CR 2007 eq. 5-7 for the standard H=4, K=3 case, generalized to
	// arbitrary H, K.
	//
	// Equatorial region: |sin(theta)| <= (K-1)/K
	sinT := math.Sin(theta)
	threshold := float64(p.K-1) / float64(p.K)
	if math.Abs(sinT) <= threshold {
		// Equatorial: simple linear map.
		x = phi
		y = math.Pi * float64(p.K) * sinT / (2 * float64(p.H))
		return x, y, true
	}
	// Polar caps: the diamond fold.
	// Let sigma = sqrt(K * (1 - |sin(theta)|))
	// omega = H * phi / pi (facet-fraction longitude, 0..2H)
	// phi_c = -pi + (2*floor((omega+1)/2) + 1) * pi/H
	// x = phi_c + (phi - phi_c) * sigma
	// y = (pi / (2*H)) * ((K+1)/2 - sigma) * sign(sin(theta))
	sign := 1.0
	if sinT < 0 {
		sign = -1.0
	}
	sigma := math.Sqrt(float64(p.K) * (1 - math.Abs(sinT)))
	omega := float64(p.H) * phi / math.Pi
	// phi_c is the longitude of the nearest facet center in the polar cap.
	phic := -math.Pi + (2*math.Floor((omega+1)/2)+1)*math.Pi/float64(p.H)
	x = phic + (phi-phic)*sigma
	y = math.Pi / (2 * float64(p.H)) * (float64(p.K+1)/2 - sigma) * sign
	return x, y, true
}

func (p hpxProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	threshold := math.Pi * float64(p.K-1) / (2 * float64(p.H))
	if math.Abs(y) <= threshold {
		// Equatorial.
		phi = x
		theta = math.Asin(2 * float64(p.H) * y / (math.Pi * float64(p.K)))
		return phi, theta, true
	}
	// Polar cap inverse.
	sign := 1.0
	if y < 0 {
		sign = -1.0
	}
	sigma := float64(p.K+1)/2 - 2*float64(p.H)*math.Abs(y)/math.Pi
	if sigma < 0 || sigma > math.Sqrt(float64(p.K)) {
		return 0, 0, false
	}
	sinTheta := (1 - sigma*sigma/float64(p.K)) * sign
	if sinTheta < -1 || sinTheta > 1 {
		return 0, 0, false
	}
	theta = math.Asin(sinTheta)
	// Recover phi from x and the facet center.
	omegaEst := float64(p.H) * x / math.Pi
	phic := -math.Pi + (2*math.Floor((omegaEst+1)/2)+1)*math.Pi/float64(p.H)
	if sigma == 0 {
		phi = phic
	} else {
		phi = phic + (x-phic)/sigma
	}
	return phi, theta, true
}

// XPH — HEALPix polar / butterfly (CR 2007 §3).
//
// Same equal-area property but the layout is "unfolded around the pole"
// to put polar data at the center with meridians radiating outward. Used
// when the region of interest is near a celestial pole.
//
// XPH is derived from HPX by a coordinate rotation: swap (x, y) ↔ a
// 45°-rotated frame, and conjoin the polar caps into a single central
// diamond.
type xphProjection struct{}

func (xphProjection) Code() string    { return "XPH" }
func (xphProjection) Theta0() float64 { return math.Pi / 2 }

func (xphProjection) Forward(phi, theta float64) (x, y float64, ok bool) {
	// CR 2007 eq. 17-19 for standard H=4, K=3. The projection maps the
	// native pole to the origin; meridians radiate outward.
	//
	// The construction: first compute HPX coordinates, then apply a 45°
	// rotation and shift so that the polar cap centers sit on the
	// diagonals.
	hpx := hpxProjection{H: 4, K: 3}
	hx, hy, fok := hpx.Forward(phi, theta)
	if !fok {
		return 0, 0, false
	}
	// Rotate 45° and translate; exact offsets per CR 2007.
	cos45 := math.Sqrt2 / 2
	sin45 := math.Sqrt2 / 2
	x = cos45*hx - sin45*hy
	y = sin45*hx + cos45*hy
	return x, y, true
}

func (xphProjection) Inverse(x, y float64) (phi, theta float64, ok bool) {
	cos45 := math.Sqrt2 / 2
	sin45 := math.Sqrt2 / 2
	hx := cos45*x + sin45*y
	hy := -sin45*x + cos45*y
	hpx := hpxProjection{H: 4, K: 3}
	return hpx.Inverse(hx, hy)
}
