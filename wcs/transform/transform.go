package transform

import (
	"errors"
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Transform maps between image pixel coordinates and celestial spherical
// coordinates for the celestial axis pair of a WCS.
//
// The pipeline follows wcslib 8.10 (wcsp2s / wcss2p) step for step and is
// checked against it, through header parsing, on every case in
// testdata/wref/hdrcases (real HST, JWST, Spitzer, ZTF, PTF,
// astrometry.net, ROSAT and MSX headers plus synthetic cases for each
// keyword rule):
//
//	pixel offset q = p - CRPIX                (the two celestial pixel axes)
//	SIP (prior)    q' = q + (A(q), B(q))
//	linear         m = M q'                   (M = PC, or CD, or from CROTA)
//	TPV (sequent)  m' = TPV(m)                (before CDELT, as in wcslib)
//	DSS (sequent)  m' = DSS(m)                (plate polynomial, ditto)
//	scale          (x, y) = CDELT * m'        (degrees; CUNIT converted)
//	projection     (phi, theta) = P^-1(x + x0, y + y0)
//	rotation       (alpha, delta) from (phi, theta) and the native pole
//
// Supported: any position of the celestial axes in an N-dimensional WCS
// (including swapped DEC/RA) as long as the matrix does not couple a
// non-celestial pixel axis into a celestial world axis; PCi_j, CDi_j,
// CROTAi; CUNIT deg, arcmin, arcsec, mas, rad; LONPOLE / LATPOLE and
// PVi_0..4 on the longitude axis (fiducial point and offset); SIP (read
// from the SIP keywords with or without the "-SIP" CTYPE suffix, as
// wcslib's wcspih does, except on a RA---TPV header without "-SIP",
// whose SIP set is ignored; applied together with TPV on RA---TPV-SIP);
// TPV; the DSS plate solution (wcs.Header.DSS: the plate
// keywords replace the header's linear WCS and the plate polynomial is
// applied, as wcslib does).
// The header forms New refuses (listed on New) are errors, not silently
// different results; every header in testdata/wref/hdrcases_reject, which
// wcslib 8.10 refuses, is refused here too (TestHdrRejectCases), and the
// headers in testdata/wref/hdrcases_stricter, which wcslib accepts, are
// refused on purpose (TestHdrStricterThanWcslib). Three accepted forms
// give a different answer from wcslib on purpose: a RA---TPV / DEC--TPV
// header (no "-SIP") that also carries SIP keywords is read with TPV
// alone, as astropy 8.0.2 reads it, where wcslib applies SIP and TPV
// both (wcs.Header.SIPIgnored; testdata/wref/hdrcases_differ,
// TestHdrDeliberateDifferences: on the real ATLAS header of PR #4 the SIP
// set is a copy of the TPV correction and wcslib's reading lands 13.3-14.6
// px from the image's, astropy's and ours); a TPV header that
// omits PV1_1 or PV2_1 takes the TPV registry default 1, where wcslib
// uses 0 (see ParseTPV; TestTPVMissingPV1DefaultsToOne); and sky->pixel
// on a DSS header far outside the plate cutout, where wcslib's
// separable iteration (disx2p) can fail to converge and report the point
// as invalid while Newton iteration converges to the pixel (seen 6000 px
// out on a degree-5 plate solution: testdata/wref/hdrcases/
// dss_far_outside_plate.case, where the returned pixel is checked against
// the source pixel of wcslib's own pixel->sky). Every other
// non-convergence is an error, not a value. Keywords on
// non-celestial axes do not enter the celestial answer and are not
// checked.
//
// Pixel coordinates are 0-BASED to match Go slice indexing: pixel (0, 0)
// is the center of the first image cell on disk, which corresponds to
// FITS pixel (1, 1). p1 and p2 are the coordinates along the two pixel
// axes that carry the celestial world axes, in increasing axis order (so
// for CTYPE1 = 'DEC--TAN', CTYPE2 = 'RA---TAN', p1 is still the column).
//
// Angles in the public API (alpha, delta) are in DEGREES: alpha is the
// longitude-axis coordinate and delta the latitude-axis coordinate,
// whatever their axis order. Internally the pipeline works in radians.
type Transform struct {
	hdr  *wcs.Header
	proj Projection

	// tpv is the TPV polynomial for RA---TPV / DEC--TPV headers, applied to
	// the intermediate coordinates between the matrix and CDELT. Nil for
	// every other projection.
	tpv *TPV

	// sip is the SIP polynomial, applied to the pixel offset from CRPIX
	// before the matrix whenever wcs.Parse read SIP keywords for the
	// primary description (with or without "-SIP" in CTYPE, as in wcslib;
	// not on RA---TPV without "-SIP"). Nil when
	// it carries none.
	sip *SIP

	// dss is the DSS plate polynomial (wcs.Header.DSS), a sequent
	// distortion applied between the matrix and CDELT like TPV, as in
	// wcslib. Nil without DSS keywords.
	dss *dssPoly

	// lin maps the pixel offsets (q1, q2) on the celestial pixel axes (in
	// increasing axis order) to the (longitude, latitude) intermediate
	// coordinates BEFORE the CDELT scale; scale holds that scale in
	// radians per unit (CDELT, CUNIT and deg->rad combined).
	lin    [2][2]float64
	linInv [2][2]float64
	scale  [2]float64

	// CRPIX for the two celestial pixel axes, 0-based.
	crpixX float64
	crpixY float64

	// Projection offset (wcslib prjprm x0, y0) in radians: non-zero only
	// when PVi_0 on the longitude axis requests a fiducial offset.
	x0, y0 float64

	// Celestial rotation parameters in radians: the celestial coordinates
	// of the native pole and the native longitude of the celestial pole.
	alphaP float64
	deltaP float64
	// sinDp, cosDp = sincos(deltaP), computed once (wcslib keeps them in
	// celprm euler[3], euler[4] for the same reason).
	sinDp, cosDp float64
	phiP         float64

	// eul is the same rotation as wcslib's celprm euler[] (degrees), used
	// by SkyToPixel (sphS2X, a port of wcslib sphs2x).
	eul eulerAngles

	// lonpole / latpole as derived (degrees), for diagnostics and tests.
	lonpole, latpole float64
}

// New builds a Transform from a parsed *wcs.Header. It returns an error,
// rather than a transform that would give wrong answers, when:
//
//   - the header has no single celestial longitude/latitude axis pair, or
//     the two celestial axes declare different projection or distortion
//     codes;
//   - the projection code is not supported (see ErrUnsupportedProjection
//     in projection.go), or CTYPE declares a distortion other than SIP;
//   - the header carries distortion keywords this package does not
//     implement (wcs.Header.UnsupportedDistortion: CPDIS/CQDIS/D2IM,
//     or SIP on an alternate description);
//   - the linear matrix couples a non-celestial pixel axis into a
//     celestial world axis, or is singular;
//   - a celestial CUNIT is not deg/arcmin/arcsec/mas/rad (blank = deg),
//     or a TPV header uses a celestial CUNIT other than deg;
//   - SIP is declared on celestial pixel axes other than 1 and 2;
//   - a PVi_m index is out of range: m > 4 on the longitude axis of a
//     non-TPV projection, m > 39 on either celestial axis of TPV (both
//     refused by wcslib);
//   - LONPOLE / LATPOLE / PVi_m admit no valid native pole (wcslib
//     celset "ill-conditioned coordinate transformation parameters").
//
// SIP keywords are applied with or without the "-SIP" CTYPE suffix, on
// any projection New supports, as in wcslib, except that wcs.Parse drops
// them from a RA---TPV / DEC--TPV header without "-SIP" (astropy's rule;
// wcs.Header.SIPIgnored), so New applies TPV alone there; RA---TPV-SIP
// gets SIP then TPV, as in wcslib. A
// "-SIP" CTYPE without SIP keywords is undistorted, as in wcslib. SIP
// keywords on a TNX / ZPX header are refused by wcs.Parse (so by New and
// NewTNX too), where wcslib applies SIP and ignores the WAT correction
// when the SIP cards come last (hdrcases_stricter/tnx_sip_*).
func New(h *wcs.Header) (*Transform, error) {
	if h == nil {
		return nil, errors.New("wcs/transform: nil header")
	}
	if h.CelestialError != "" {
		return nil, errors.New("wcs/transform: " + h.CelestialError)
	}
	if !h.IsCelestial() {
		return nil, errors.New("wcs/transform: header is not a celestial WCS")
	}
	if h.NAxis < 2 {
		return nil, fmt.Errorf("wcs/transform: need NAxis >= 2, got %d", h.NAxis)
	}
	if len(h.UnsupportedDistortion) > 0 {
		return nil, fmt.Errorf("wcs/transform: header uses distortions this package does not implement: %v", h.UnsupportedDistortion)
	}
	lng, lat := h.LonAxis-1, h.LatAxis-1
	if strAt(h.ProjCode, lng) != strAt(h.ProjCode, lat) {
		return nil, fmt.Errorf("wcs/transform: inconsistent projection codes %q (CTYPE%d) and %q (CTYPE%d)",
			strAt(h.ProjCode, lng), h.LonAxis, strAt(h.ProjCode, lat), h.LatAxis)
	}
	if strAt(h.DistCode, lng) != strAt(h.DistCode, lat) {
		return nil, fmt.Errorf("wcs/transform: inconsistent distortion codes %q (CTYPE%d) and %q (CTYPE%d)",
			strAt(h.DistCode, lng), h.LonAxis, strAt(h.DistCode, lat), h.LatAxis)
	}

	// Celestial pixel axes in increasing order.
	pix1, pix2 := lng, lat
	if pix1 > pix2 {
		pix1, pix2 = pix2, pix1
	}

	// SIP applies whenever wcs.Parse read SIP keywords (h.SIP), with or
	// without a "-SIP" suffix, as in wcslib; on RA---TPV without "-SIP"
	// it did not read them (see wcs.Header.SIP).
	switch dist := strAt(h.DistCode, lng); dist {
	case "", "SIP":
	default:
		return nil, fmt.Errorf("wcs/transform: unsupported distortion %q in CTYPE", dist)
	}
	sip := sipFromWCS(h.SIP)
	if sip != nil && (pix1 != 0 || pix2 != 1) {
		return nil, fmt.Errorf("wcs/transform: SIP on celestial pixel axes %d,%d (SIP is defined on pixel axes 1,2 only)", pix1+1, pix2+1)
	}

	code := h.CelestialProjCode()
	// PVi_m index ranges wcslib enforces rather than ignoring the card:
	// TPV defines PVi_0..39 on each celestial axis (wcslib tpvset "Invalid
	// parameter number"), any other projection PVi_0..4 on the longitude
	// axis (Paper II sect. 2.5; wcslib wcsset "Unrecognized coordinate
	// transformation parameter"). Latitude-axis projection parameters are
	// the projection's business (Select).
	for k := range h.PV {
		switch {
		case code == "TPV" && (k.Axis == h.LonAxis || k.Axis == h.LatAxis) && k.Index > 39:
			return nil, fmt.Errorf("wcs/transform: PV%d_%d: TPV defines PVi_0..39 only", k.Axis, k.Index)
		case code != "TPV" && k.Axis == h.LonAxis && k.Index > 4:
			return nil, fmt.Errorf("wcs/transform: PV%d_%d: only PVi_0..4 are defined on the celestial longitude axis", k.Axis, k.Index)
		}
	}
	pv := h.PV
	// AIPS aliases, translated as wcslib 8.10 wcsset does ("simple alias
	// translations", wcs.c): NCP is SIN with PV2_1 = 0 and PV2_2 =
	// cot(CRVAL of the latitude axis), overriding any PVi_1/PVi_2 given on
	// that axis; GLS is SFL with the fiducial offset phi0 = 0, theta0 =
	// CRVAL of the latitude axis (applied below).
	gls := false
	switch code {
	case "NCP":
		u, _ := angleUnit(strAt(h.CUnit, lat), h.LatAxis) // a bad unit is reported below
		d0 := float64(h.CRVal[lat] * u)
		if d0 == 0 {
			return nil, errors.New("wcs/transform: NCP projection with the reference point on the equator (cot(0) is infinite; wcslib: \"NCP blows up on the equator\")")
		}
		pv = make(map[wcs.PVKey]float64, len(h.PV)+2)
		for k, v := range h.PV {
			pv[k] = v
		}
		pv[wcs.PVKey{Axis: h.LatAxis, Index: 1}] = 0
		pv[wcs.PVKey{Axis: h.LatAxis, Index: 2}] = cosd(d0) / sind(d0)
		code = "SIN"
	case "GLS":
		code, gls = "SFL", true
	}
	var tpv *TPV
	if code == "TPV" {
		// TPV is TAN plus a polynomial distortion. Its PVi_m are the
		// polynomial coefficients, not TAN projection parameters (wcslib
		// moves them to DQia), so they are not passed to the projection
		// and do not set the fiducial point below.
		tpv = ParseTPV(h)
		code, pv = "TAN", nil
	}
	proj, err := Select(code, pv, h.LatAxis)
	if err != nil {
		return nil, err
	}

	t := &Transform{hdr: h, proj: proj, tpv: tpv, sip: sip}
	if h.DSS != nil {
		if tpv != nil || sip != nil || pix1 != 0 || pix2 != 1 {
			return nil, errors.New("wcs/transform: DSS plate solution combined with another distortion or on celestial axes other than 1,2")
		}
		t.dss = newDSSPoly(h.DSS)
	}

	// Celestial units -> degrees (wcslib wcsset/wcs_units).
	ulng, err := angleUnit(strAt(h.CUnit, lng), h.LonAxis)
	if err != nil {
		return nil, err
	}
	ulat, err := angleUnit(strAt(h.CUnit, lat), h.LatAxis)
	if err != nil {
		return nil, err
	}
	if tpv != nil && (ulng != 1 || ulat != 1) {
		return nil, fmt.Errorf("wcs/transform: TPV with celestial CUNIT %q/%q: TPV coefficients are defined in degrees", strAt(h.CUnit, lng), strAt(h.CUnit, lat))
	}
	if err := t.buildLinear(h, lng, lat, pix1, pix2, ulng, ulat); err != nil {
		return nil, err
	}
	if err := t.invertLinear(); err != nil {
		return nil, err
	}

	// CRPIX: convert FITS 1-based to 0-based by subtracting 1.
	t.crpixX = h.CRPix[pix1] - 1
	t.crpixY = h.CRPix[pix2] - 1

	// Fiducial point and native pole (wcslib wcsset -> celset).
	lng0 := float64(h.CRVal[lng] * ulng)
	lat0 := float64(h.CRVal[lat] * ulat)
	phi0, theta0 := 0.0, radToDeg(proj.Theta0())
	lonpole, lonpoleSet := h.LonPole, h.LonPoleSet
	latpole := 90.0
	if h.LatPoleSet {
		latpole = h.LatPole
	}
	if tpv == nil {
		// PVi_0..4 on the longitude axis (Paper II §2.5), which take
		// precedence over LONPOLE / LATPOLE.
		pvl := func(m int) (float64, bool) {
			v, ok := h.PV[wcs.PVKey{Axis: h.LonAxis, Index: m}]
			return v, ok
		}
		offset := false
		if v, ok := pvl(0); ok && v != 0 {
			offset = true
		}
		v1, set1 := pvl(1)
		v2, set2 := pvl(2)
		if set1 {
			phi0 = v1
		}
		if set2 {
			theta0 = v2
		}
		if v, ok := pvl(3); ok {
			lonpole, lonpoleSet = v, true
		}
		if v, ok := pvl(4); ok {
			latpole = v
		}
		if gls {
			// wcslib wcsset for GLS: offset = 1, phi0 = 0, theta0 = CRVAL.
			offset, set1, set2 = true, true, true
			phi0, theta0 = 0, lat0
		}
		if math.Abs(theta0) > 90 {
			if math.Abs(theta0) > 90+celTol {
				return nil, fmt.Errorf("wcs/transform: PV%d_2 = theta0 = %g is outside [-90, 90]", h.LonAxis, theta0)
			}
			theta0 = math.Copysign(90, theta0)
		}
		if offset && set1 && set2 {
			// (x, y) = (0, 0) moves to the fiducial point (wcslib prjoff).
			x0, y0, ok := proj.Forward(degToRad(phi0), degToRad(theta0))
			if !ok {
				return nil, fmt.Errorf("wcs/transform: fiducial point (phi0, theta0) = (%g, %g) is outside the %s projection", phi0, theta0, proj.Code())
			}
			t.x0, t.y0 = x0, y0
		}
	}
	lngp, latp, phip, err := celPole(lng0, lat0, phi0, theta0, lonpole, lonpoleSet, latpole)
	if err != nil {
		return nil, err
	}
	t.alphaP, t.deltaP, t.phiP = degToRad(lngp), degToRad(latp), degToRad(phip)
	t.sinDp, t.cosDp = lmSincos(t.deltaP)
	t.eul = newEulerAngles(lngp, latp, phip)
	t.lonpole, t.latpole = phip, latp
	return t, nil
}

// strAt / floatAt read a per-axis slice of a (possibly hand-built)
// wcs.Header, treating a missing entry as the keyword's default.
func strAt(s []string, i int) string {
	if i < len(s) {
		return s[i]
	}
	return ""
}

func floatAt(s []float64, i int, def float64) float64 {
	if i < len(s) {
		return s[i]
	}
	return def
}

// angleUnit returns the factor converting a celestial axis's CUNIT to
// degrees. FITS units are case sensitive; a blank CUNIT means degrees.
func angleUnit(u string, axis int) (float64, error) {
	switch u {
	case "", "deg":
		return 1, nil
	case "arcmin":
		return 1.0 / 60, nil
	case "arcsec":
		return 1.0 / 3600, nil
	case "mas":
		return 1.0 / 3600e3, nil
	case "rad":
		return 180 / math.Pi, nil
	}
	return 0, fmt.Errorf("wcs/transform: CUNIT%d = %q is not a supported angle unit (deg, arcmin, arcsec, mas, rad)", axis, u)
}

// buildLinear selects the linear transformation exactly as wcslib does
// (wcs.h altlin): PC when any PCi_j is present, else CD (CDELT ignored),
// else CROTA of the latitude axis with CDELT, else identity PC with
// CDELT. It keeps the 2x2 block mapping the celestial pixel axes to the
// celestial world axes and refuses a matrix that mixes any other pixel
// axis into them.
func (t *Transform) buildLinear(h *wcs.Header, lng, lat, pix1, pix2 int, ulng, ulat float64) error {
	n := h.NAxis
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	cdelt := [2]float64{h.CDelt[lng], h.CDelt[lat]}
	switch {
	case h.HasPC || (h.CD == nil && !h.HasCRota):
		for i := range m {
			copy(m[i], h.PC[i])
		}
	case h.CD != nil:
		for i := range m {
			copy(m[i], h.CD[i])
		}
		cdelt = [2]float64{1, 1}
	default:
		// AIPS CROTA (wcslib wcsset): rotation of the latitude axis.
		for i := range m {
			m[i][i] = 1
		}
		if h.CDelt[lng] == 0 {
			return errors.New("wcs/transform: CDELT of the longitude axis is zero (singular CROTA transformation)")
		}
		rho := floatAt(h.CRota, lat, 0)
		lambda := h.CDelt[lat] / h.CDelt[lng]
		s, c := sincosd(rho)
		m[lng][lng], m[lat][lat] = c, c
		m[lng][lat] = float64(-lambda * s)
		m[lat][lng] = s / lambda
	}
	for _, r := range []int{lng, lat} {
		for j := 0; j < n; j++ {
			if j != pix1 && j != pix2 && m[r][j] != 0 {
				return fmt.Errorf("wcs/transform: linear matrix couples pixel axis %d into celestial world axis %d; only the celestial pixel axes may feed the celestial world axes", j+1, r+1)
			}
		}
	}
	t.lin[0][0], t.lin[0][1] = m[lng][pix1], m[lng][pix2]
	t.lin[1][0], t.lin[1][1] = m[lat][pix1], m[lat][pix2]
	deg2rad := math.Pi / 180
	t.scale[0] = float64(cdelt[0] * ulng * deg2rad)
	t.scale[1] = float64(cdelt[1] * ulat * deg2rad)
	if t.scale[0] == 0 || t.scale[1] == 0 || math.IsNaN(t.scale[0]) || math.IsNaN(t.scale[1]) {
		return fmt.Errorf("wcs/transform: zero or invalid CDELT on a celestial axis (%g, %g)", cdelt[0], cdelt[1])
	}
	return nil
}

// invertLinear computes the inverse of the 2x2 linear matrix.
func (t *Transform) invertLinear() error {
	a, b := t.lin[0][0], t.lin[0][1]
	c, d := t.lin[1][0], t.lin[1][1]
	det := float64(a*d) - float64(b*c)
	if det == 0 || math.IsNaN(det) || math.IsInf(det, 0) {
		return fmt.Errorf("wcs/transform: linear transform is singular (det=%g)", det)
	}
	t.linInv[0][0] = d / det
	t.linInv[0][1] = -b / det
	t.linInv[1][0] = -c / det
	t.linInv[1][1] = a / det
	return nil
}

// PixelToSky maps 0-based pixel coordinates (p1, p2) on the celestial
// pixel axes (increasing axis order) to celestial coordinates (alpha,
// delta) in degrees. alpha is normalized to [0, 360).
func (t *Transform) PixelToSky(p1, p2 float64) (alpha, delta float64, err error) {
	dx := p1 - t.crpixX
	dy := p2 - t.crpixY

	if t.sip != nil {
		dx, dy = t.sip.Forward(dx, dy)
	}

	mx := float64(t.lin[0][0]*dx) + float64(t.lin[0][1]*dy)
	my := float64(t.lin[1][0]*dx) + float64(t.lin[1][1]*dy)

	if t.tpv != nil {
		mx, my = t.tpv.Forward(mx, my)
	}
	if t.dss != nil {
		mx, my = t.dss.Forward(mx, my)
	}

	x := float64(mx*t.scale[0]) + t.x0
	y := float64(my*t.scale[1]) + t.y0

	phi, theta, ok := t.proj.Inverse(x, y)
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: pixel (%g, %g) outside %s projection domain", p1, p2, t.proj.Code())
	}

	a, d := nativeToCelestial(phi, theta, t.alphaP, t.sinDp, t.cosDp, t.phiP)
	return normalizeDeg(radToDeg(a)), radToDeg(d), nil
}

// SkyToPixel is the inverse of PixelToSky. Returns 0-based pixel
// coordinates in the same axis order as PixelToSky takes them.
func (t *Transform) SkyToPixel(alpha, delta float64) (p1, p2 float64, err error) {
	phiD, thetaD := sphS2X(&t.eul, alpha, delta)

	x, y, ok := t.proj.Forward(degToRad(phiD), degToRad(thetaD))
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: sky (%g, %g) outside %s projection domain", alpha, delta, t.proj.Code())
	}
	mx := (x - t.x0) / t.scale[0]
	my := (y - t.y0) / t.scale[1]

	if t.tpv != nil {
		if mx, my, ok = t.tpv.Inverse(mx, my); !ok {
			return 0, 0, fmt.Errorf("wcs/transform: TPV inverse did not converge at sky (%g, %g)", alpha, delta)
		}
	}

	if t.dss != nil {
		if mx, my, ok = t.dss.Inverse(mx, my); !ok {
			return 0, 0, fmt.Errorf("wcs/transform: DSS plate polynomial inverse did not converge at sky (%g, %g)", alpha, delta)
		}
	}

	dx := float64(t.linInv[0][0]*mx) + float64(t.linInv[0][1]*my)
	dy := float64(t.linInv[1][0]*mx) + float64(t.linInv[1][1]*my)

	if t.sip != nil {
		if dx, dy, ok = t.sip.Inverse(dx, dy); !ok {
			return 0, 0, fmt.Errorf("wcs/transform: SIP inverse did not converge at sky (%g, %g)", alpha, delta)
		}
	}

	return dx + t.crpixX, dy + t.crpixY, nil
}

// Code returns the projection code in use (e.g. "TAN", or "TPV" for TAN
// with a TPV distortion).
func (t *Transform) Code() string {
	if t.tpv != nil {
		return "TPV"
	}
	return t.proj.Code()
}
