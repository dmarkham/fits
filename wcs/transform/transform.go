package transform

import (
	"errors"
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Transform maps between image pixel coordinates and celestial spherical
// coordinates for a 2-D celestial WCS.
//
// Higher-dimensional WCS with more than one celestial axis pair (rare in
// practice) is not supported in v1.
//
// Pixel coordinates are 0-BASED to match Go slice indexing: pixel (0, 0)
// is the center of the first image cell on disk, which corresponds to
// FITS pixel (1, 1). The library applies the +1 offset internally.
//
// Angles in the public API (alpha, delta) are in DEGREES. Internally the
// pipeline works in radians.
type Transform struct {
	hdr  *wcs.Header
	proj Projection

	// Linear forward transform: matrix M such that
	//    (x, y) = M * (pixel - CRPIX)
	// where (x, y) is the intermediate world coordinate in radians (the
	// library pre-scales by CDELT and by deg→rad so the projection step
	// can work in radians directly).
	//
	// lin[0] = {m00, m01}, lin[1] = {m10, m11}.
	lin    [2][2]float64
	linInv [2][2]float64

	// CRPIX for the two celestial axes, in 0-based pixel indexing (so the
	// "FITS 1-based → Go 0-based" offset of -1 is already baked in).
	crpixX float64
	crpixY float64

	// Celestial rotation parameters in radians.
	alphaP float64 // reference-point longitude (CRVAL of lon axis, radians)
	deltaP float64 // reference-point latitude (CRVAL of lat axis, radians)
	phiP   float64 // native longitude of the celestial pole (radians)
}

// New builds a Transform from a parsed *wcs.Header. Returns an error if
// the header is not a 2-D celestial WCS, if the linear transform is
// singular, or if the projection code is not yet implemented (see
// ErrUnsupportedProjection in projection.go).
func New(h *wcs.Header) (*Transform, error) {
	if h == nil {
		return nil, errors.New("wcs/transform: nil header")
	}
	if !h.IsCelestial() {
		return nil, errors.New("wcs/transform: header is not a celestial WCS")
	}
	if h.NAxis < 2 {
		return nil, fmt.Errorf("wcs/transform: need NAxis >= 2, got %d", h.NAxis)
	}
	// Today we require the celestial axes to be axis 1 and axis 2. Headers
	// with celestial axes at other positions (rare) can be handled later.
	if h.LonAxis != 1 || h.LatAxis != 2 {
		return nil, fmt.Errorf("wcs/transform: celestial axes must be 1,2 (got lon=%d lat=%d)",
			h.LonAxis, h.LatAxis)
	}

	code := h.CelestialProjCode()
	proj, err := Select(code, h.PV, h.LatAxis)
	if err != nil {
		return nil, err
	}

	t := &Transform{
		hdr:  h,
		proj: proj,
	}

	// Build the pixel-offset-to-radians linear matrix from either CD or
	// PC+CDELT, assuming CUNIT is "deg" or empty (treated as "deg"). CUNIT
	// of "rad" would skip the deg→rad scaling; other units (e.g. "arcsec")
	// are rare for celestial axes and not supported in v1.
	if err := t.buildLinear(h); err != nil {
		return nil, err
	}

	// Invert the 2x2 linear matrix for the sky-to-pixel direction.
	if err := t.invertLinear(); err != nil {
		return nil, err
	}

	// CRPIX: convert FITS 1-based to 0-based by subtracting 1.
	t.crpixX = h.CRPix[h.LonAxis-1] - 1
	t.crpixY = h.CRPix[h.LatAxis-1] - 1

	// Celestial rotation setup. For projections where the fiducial point
	// is the native pole (theta_0 = pi/2, all zenithal), (alphaP, deltaP) =
	// (CRVAL_lon, CRVAL_lat) directly. For other projections with
	// theta_0 = 0 (cylindrical, pseudo-cylindrical, HPX, PCO, TSC/CSC/QSC),
	// the fiducial point is on the native equator and (alphaP, deltaP) is
	// the celestial coordinate of the NATIVE POLE, which lies at a
	// distance (pi/2 - theta_0) from the reference point along the
	// meridian containing the reference.
	//
	// For conic projections (theta_0 = thetaA), the same general rule
	// applies: rotate the reference from native (0, theta_0) to the
	// native pole (any phi, pi/2) through angle (pi/2 - theta_0).
	alpha0 := degToRad(h.CRVal[h.LonAxis-1])
	delta0 := degToRad(h.CRVal[h.LatAxis-1])
	theta0 := proj.Theta0()
	// LONPOLE default per Paper II §2.4 eq. 8. For zenithal projections
	// (theta_0 = pi/2) the rule is phi_p = 0 if delta_0 >= theta_0 else
	// 180°; for non-zenithal (theta_0 = 0 or other small value) the
	// default is 180° regardless of delta_0, matching astropy.wcs /
	// wcslib behavior. We unify this into a single rule: default to 180°
	// unless zenithal with delta_0 == theta_0 (pole observation, where
	// LONPOLE = 0° is the conventional choice).
	lonPole := 180.0
	if h.LonPoleSet {
		lonPole = h.LonPole
	} else if theta0 == math.Pi/2 && delta0 >= theta0 {
		lonPole = 0
	}
	t.phiP = degToRad(lonPole)

	// Compute the celestial coordinates of the native pole given the
	// reference point (alpha0, delta0) at native (phi_0, theta_0) and the
	// native longitude of the celestial pole (phi_p = lonPole).
	// For zenithal projections theta_0 = pi/2 → native pole == reference.
	// For other projections theta_0 != pi/2 → apply Paper II §2.4 eq. 7.
	if theta0 == math.Pi/2 {
		t.alphaP = alpha0
		t.deltaP = delta0
	} else {
		// Non-zenithal: the fiducial point sits on the native equator
		// (theta_0 = 0) or on a small parallel (conics). The native pole
		// is (pi/2 - theta_0) "above" the reference along the celestial
		// meridian containing it.
		//
		// For the library default LONPOLE = 180° (which matches
		// astropy.wcs / wcslib for non-zenithal headers), the rotation
		// parameters are:
		//
		//	alphaP = alpha_0
		//	deltaP = delta_0 + (pi/2 - theta_0)
		//	phiP   = pi (i.e. LONPOLE = 180°)
		//
		// The pole is verified to place the fiducial point correctly
		// via Paper II eq. 2 by construction.
		t.alphaP = alpha0
		t.deltaP = delta0 + (math.Pi/2 - theta0)
	}

	return t, nil
}

// buildLinear constructs the pixel-offset → radians linear matrix from
// either CD or PC+CDELT, assuming CUNIT is degrees.
func (t *Transform) buildLinear(h *wcs.Header) error {
	// We need the 2x2 sub-matrix for the celestial axis pair (which we
	// already required to be axes 1 and 2). Convert the degrees-per-pixel
	// values to radians-per-pixel by multiplying each row by pi/180.
	deg2rad := math.Pi / 180
	if h.CD != nil {
		t.lin[0][0] = h.CD[0][0] * deg2rad
		t.lin[0][1] = h.CD[0][1] * deg2rad
		t.lin[1][0] = h.CD[1][0] * deg2rad
		t.lin[1][1] = h.CD[1][1] * deg2rad
		return nil
	}
	// PC form: lin = diag(CDELT) * PC, then × deg2rad.
	cd1 := h.CDelt[0] * deg2rad
	cd2 := h.CDelt[1] * deg2rad
	t.lin[0][0] = cd1 * h.PC[0][0]
	t.lin[0][1] = cd1 * h.PC[0][1]
	t.lin[1][0] = cd2 * h.PC[1][0]
	t.lin[1][1] = cd2 * h.PC[1][1]
	return nil
}

// invertLinear computes the inverse of the 2x2 linear matrix.
func (t *Transform) invertLinear() error {
	a, b := t.lin[0][0], t.lin[0][1]
	c, d := t.lin[1][0], t.lin[1][1]
	det := a*d - b*c
	if det == 0 || math.IsNaN(det) || math.IsInf(det, 0) {
		return fmt.Errorf("wcs/transform: linear transform is singular (det=%g)", det)
	}
	t.linInv[0][0] = d / det
	t.linInv[0][1] = -b / det
	t.linInv[1][0] = -c / det
	t.linInv[1][1] = a / det
	return nil
}


// PixelToSky maps 0-based pixel coordinates (p1, p2) to celestial
// coordinates (alpha, delta) in degrees. alpha is normalized to [0, 360).
func (t *Transform) PixelToSky(p1, p2 float64) (alpha, delta float64, err error) {
	// Step 1: pixel offset.
	dx := p1 - t.crpixX
	dy := p2 - t.crpixY

	// Step 2: linear to intermediate (radians).
	x := t.lin[0][0]*dx + t.lin[0][1]*dy
	y := t.lin[1][0]*dx + t.lin[1][1]*dy

	// Step 3: projection → native (phi, theta) in radians.
	phi, theta, ok := t.proj.Inverse(x, y)
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: pixel (%g, %g) outside %s projection domain", p1, p2, t.proj.Code())
	}

	// Step 4: native → celestial rotation.
	a, d := nativeToCelestial(phi, theta, t.alphaP, t.deltaP, t.phiP)
	return normalizeDeg(radToDeg(a)), radToDeg(d), nil
}

// SkyToPixel is the inverse of PixelToSky. Returns 0-based pixel
// coordinates.
func (t *Transform) SkyToPixel(alpha, delta float64) (p1, p2 float64, err error) {
	// Step 4 (inverse): celestial → native.
	phi, theta := celestialToNative(degToRad(alpha), degToRad(delta), t.alphaP, t.deltaP, t.phiP)

	// Step 3 (inverse): native → intermediate.
	x, y, ok := t.proj.Forward(phi, theta)
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: sky (%g, %g) outside %s projection domain", alpha, delta, t.proj.Code())
	}

	// Step 2 (inverse): intermediate → pixel offset.
	dx := t.linInv[0][0]*x + t.linInv[0][1]*y
	dy := t.linInv[1][0]*x + t.linInv[1][1]*y

	// Step 1 (inverse): add CRPIX.
	return dx + t.crpixX, dy + t.crpixY, nil
}

// Code returns the projection code in use (e.g. "TAN").
func (t *Transform) Code() string { return t.proj.Code() }
