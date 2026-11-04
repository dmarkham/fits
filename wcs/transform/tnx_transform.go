package transform

import (
	"errors"
	"fmt"
	"sort"
	"strings"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// TNXTransform is the pixel <-> sky transform for an IRAF TNX header
// ("RA---TNX" / "DEC--TNX"): the linear CD (or PC + CDELT) step, the TNX
// polynomial correction from the WAT1_nnn / WAT2_nnn strings, then the
// gnomonic (TAN) projection and celestial rotation.
//
// It exists beside Transform because the WAT strings are not part of
// *wcs.Header: New(*wcs.Header) rejects a TNX header (it cannot see the
// correction and refuses to return an uncorrected TAN), and NewTNX reads
// the raw *header.Header instead.
//
// Pixel coordinates are 0-based, angles in degrees, exactly as Transform.
// Proven against WCSTools 3.9.7 (xy2sky / sky2xy) on real NOAO headers
// and synthetic headers covering every surface type and cross-term mode;
// see tnx_golden_test.go and testdata/tnx/PROVENANCE.md.
type TNXTransform struct {
	base *Transform // the TAN part (linear step, projection, rotation)
	tnx  *TNX
}

// NewTNX builds a TNXTransform from a raw FITS header whose celestial
// CTYPEs are a matching TNX pair on axes 1 and 2: RA---TNX / DEC--TNX,
// GLON-TNX / GLAT-TNX, ELON-TNX / ELAT-TNX, SLON-TNX / SLAT-TNX or
// HLON-TNX / HLAT-TNX (each proven by a golden; the result is in that
// frame, no frame conversion is done). A mixed pair such as RA / GLAT is
// an error.
//
// Conventions, each proven by a WCSTools 3.9.7 golden (testdata/tnx):
//   - linear step: the CD matrix if any CDi_j card is present (CDELT is
//     then ignored), else PC + CDELT with CDELT in degrees and both
//     CDELT1 and CDELT2 given (PC read only with PC1_1 or PC001001
//     present, as WCSTools does), else 1 deg/pixel with no PC;
//   - LONPOLE: the LONPOLE card if present (at most 360; negative values
//     are used as given), else a "longpole" attribute in the WAT strings,
//     else 180 deg (NOT the Paper II zenithal default of 0 deg when
//     CRVAL2 = +90); LATPOLE as wcslib (WCSTools agrees);
//   - PVi_m on the latitude axis are ignored when the header has a
//     correction surface (TAN has no parameters, WCSTools ignores them);
//     PROJPn (which wcs.Parse reads as latitude-axis PVi_n) is ignored
//     with or without a surface;
//   - lngcor is read from WAT1 (falling back to WAT2), latcor from WAT2
//     (falling back to WAT1); a header with neither is plain TAN.
//
// Returns an error, never a silently uncorrected transform, when: the
// CTYPEs are not a matching TNX pair or are not upper case (WCSTools does
// not read 'RA---tnx'); a WAT card sequence has a gap, mixed key
// spellings, a duplicated key, a key number above 499 or a CONTINUE card,
// or a joined WAT string is longer than 1999 characters (WCSTools, the
// truth source, silently drops or cuts those; see concatWAT); a WAT
// string declares a wtype other than tnx (e.g. zpx); a surface is
// malformed (unknown type, non-integer or < 1 order, unknown cross-term
// code, empty domain, wrong coefficient count); the WAT strings carry a
// non-default "ro"; the header uses CROTAi with no CD matrix (CROTA is
// not applied by the linear step); or the header is one that WCSTools
// and wcslib read differently or where WCSTools ignores a card
// (testdata/tnx/refused, each with its WCSTools answer): any PVi_m on
// the longitude axis; latitude-axis PVi_m with no correction surface
// (WCSTools then applies them as TPV); both CD and PC matrices; a
// CD00i00j card without its CDi_j; with no CD matrix, only one of
// CDELT1/CDELT2, or a PC matrix without CDELT1 or without PC1_1; SIP
// keywords, with or without a -SIP CTYPE suffix (wcs.Parse refuses SIP on
// a TNX header: wcslib would apply SIP + TAN and ignore the WAT correction
// when the SIP cards come last); LONPOLE above 360; a
// celestial CUNIT other than 'deg'; or a WATn_k key given twice.
func NewTNX(h *header.Header) (*TNXTransform, error) {
	if h == nil {
		return nil, errors.New("wcs/transform: nil header")
	}
	w, err := wcs.Parse(h)
	if err != nil {
		return nil, err
	}
	if !w.IsCelestial() {
		return nil, errors.New("wcs/transform: header is not a celestial WCS")
	}
	if w.ProjCode[w.LonAxis-1] != "TNX" || w.ProjCode[w.LatAxis-1] != "TNX" {
		return nil, fmt.Errorf("wcs/transform: NewTNX needs RA---TNX/DEC--TNX, got %q/%q",
			w.CType[w.LonAxis-1], w.CType[w.LatAxis-1])
	}
	pairs := map[string]string{"RA": "DEC", "GLON": "GLAT", "ELON": "ELAT", "SLON": "SLAT", "HLON": "HLAT"}
	if lon, lat := w.AxisType[w.LonAxis-1], w.AxisType[w.LatAxis-1]; pairs[lon] != lat {
		return nil, fmt.Errorf("wcs/transform: TNX axes %q/%q are not a matching celestial pair", w.CType[w.LonAxis-1], w.CType[w.LatAxis-1])
	}
	// WCSTools handles DEC on axis 1, but the WAT lngcor/latcor reading
	// here has only been proven for the lon-first order; refuse the
	// other rather than mis-assign the corrections.
	if w.LonAxis != 1 || w.LatAxis != 2 {
		return nil, fmt.Errorf("wcs/transform: TNX celestial axes must be 1,2 (got lon=%d lat=%d)", w.LonAxis, w.LatAxis)
	}
	// WCSTools recognises the TNX CTYPEs only in upper case (a lowercase
	// 'RA---tnx' is a LINEAR axis to it, and xy2sky prints asterisks),
	// so a lowercase CTYPE has no truth to be proven against.
	for _, ax := range []int{w.LonAxis, w.LatAxis} {
		key := fmt.Sprintf("CTYPE%d", ax)
		if v, err := h.String(key); err == nil && v != strings.ToUpper(v) {
			return nil, fmt.Errorf("wcs/transform: %s = %q is not upper case; WCSTools (the TNX truth source) does not read it as TNX; not supported", key, v)
		}
	}
	// SIP keywords (A_ORDER, A_p_q, ...): with the -SIP suffix WCSTools
	// applies SIP to the pixel offsets before TNX, without it WCSTools
	// ignores them, and wcslib 8.10 applies SIP with TAN and drops the WAT
	// correction. Two distortion solutions with no agreed meaning.
	if wcs.HasSIPKeywords(h) {
		return nil, errors.New("wcs/transform: TNX header carries SIP keywords: two distortion solutions (WCSTools applies SIP before TNX only with a -SIP CTYPE suffix and ignores it otherwise; wcslib applies SIP and drops the WAT correction); not supported")
	}
	t, err := ParseTNX(h)
	if err != nil {
		return nil, err
	}
	if t == nil {
		t = &TNX{}
	}
	// No correction surface at all (no WAT cards, or WAT strings with
	// neither lngcor nor latcor): WCSTools tnxinit fails and wcsinit
	// demotes the projection to TAN, whose TPV block then applies any
	// PVi_m (PV1_m default 0); wcslib TAN ignores latitude-axis PVi_m.
	// PROJPn is not a PVi_m card to WCSTools and stays accepted.
	if !t.Lon.present && !t.Lat.present {
		if k := tnxLatPVCard(h, w.LatAxis); k != "" {
			return nil, fmt.Errorf("wcs/transform: TNX header has no TNX correction surface and %s; WCSTools then applies PVi_m as a TPV polynomial, wcslib TAN ignores them; not supported", k)
		}
	}
	if err := tnxCheckLinear(h, w); err != nil {
		return nil, err
	}
	// Where the wcslib conventions of wcs.Parse / New and WCSTools (the
	// TNX truth) read a header differently, or WCSTools ignores a card,
	// refuse it (testdata/tnx/refused, each with its WCSTools answer).
	//
	// PVi_m on the longitude axis (wcslib: fiducial offset, phi0, theta0,
	// LONPOLE, LATPOLE): WCSTools ignores PVi_m on a TNX header.
	var pvLon []string
	for k := range w.PV {
		if k.Axis == w.LonAxis {
			pvLon = append(pvLon, fmt.Sprintf("PV%d_%d", k.Axis, k.Index))
		}
	}
	if len(pvLon) > 0 {
		sort.Strings(pvLon)
		return nil, fmt.Errorf("wcs/transform: TNX header has %s; WCSTools ignores PVi_m on TNX while wcslib applies them to the celestial rotation; not supported",
			strings.Join(pvLon, ", "))
	}
	// Both CD and PC (Paper I forbids it): WCSTools uses CD, wcslib PC.
	if w.CD != nil && w.HasPC {
		return nil, errors.New("wcs/transform: TNX header has both CD and PC matrices; WCSTools uses CD, wcslib PC; not supported")
	}
	// LONPOLE > 360: WCSTools tnxinit treats it as unset (falls back to
	// the WAT longpole or 180) and so ignores the card.
	if w.LonPoleSet && w.LonPole > 360 {
		return nil, fmt.Errorf("wcs/transform: TNX header has LONPOLE = %g; WCSTools ignores a LONPOLE above 360; not supported", w.LonPole)
	}
	// CUNIT: the TNX surfaces are in degrees and WCSTools ignores CUNIT
	// (CDELT is read as degrees), so only blank or 'deg' is accepted.
	for _, ax := range []int{w.LonAxis, w.LatAxis} {
		if u := w.CUnit[ax-1]; u != "" && u != "deg" {
			return nil, fmt.Errorf("wcs/transform: TNX header has CUNIT%d = %q; WCSTools ignores CUNIT and TNX is defined in degrees; only 'deg' is supported", ax, u)
		}
	}
	// CROTAi is not applied by the linear step (it reads CD or PC +
	// CDELT only); refuse rather than silently drop a rotation.
	if w.CD == nil {
		for _, k := range []string{"CROTA1", "CROTA2"} {
			if v, err := h.Float(k); err == nil && v != 0 {
				return nil, fmt.Errorf("wcs/transform: TNX header has %s = %g and no CD matrix; CROTA is not supported", k, v)
			}
		}
	}
	tan := *w
	tan.CType = append([]string(nil), w.CType...)
	tan.ProjCode = append([]string(nil), w.ProjCode...)
	for _, ax := range []int{w.LonAxis - 1, w.LatAxis - 1} {
		tan.ProjCode[ax] = "TAN"
		tan.CType[ax] = tan.CType[ax][:len(tan.CType[ax])-3] + "TAN"
	}
	if !tan.LonPoleSet {
		tan.LonPole, tan.LonPoleSet = 180, true
		if t.LonPoleSet {
			tan.LonPole = t.LonPole
		}
	}
	base, err := New(&tan)
	if err != nil {
		return nil, err
	}
	return &TNXTransform{base: base, tnx: t}, nil
}

// tnxLatPVCard returns the first PV<lat>_m card of the primary
// description (key as written), or "".
func tnxLatPVCard(h *header.Header, lat int) string {
	prefix := fmt.Sprintf("PV%d_", lat)
	for _, c := range h.Cards() {
		rest, ok := strings.CutPrefix(c.Key, prefix)
		if ok && rest != "" && strings.Trim(rest, "0123456789") == "" {
			return c.Key
		}
	}
	return ""
}

// tnxCheckLinear refuses the linear-step forms that WCSTools wcsinit (the
// TNX truth) and the wcslib conventions of wcs.Parse read differently
// (refused/cdelt*_only, pc_no_*, cd001001_only):
//   - a defunct CD00i00j card without its CDi_j: wcslib reads it as
//     CDi_j, wcsinit reads CDi_j only and falls back to CDELT;
//   - with no CD matrix, CDELT1 without CDELT2 (wcsinit: CDELT2 =
//     -CDELT1; wcslib: 1) or CDELT2 without CDELT1 (wcsinit then takes 1
//     deg/pixel and ignores CDELT2 and PC);
//   - with no CD matrix, a PC matrix without CDELT1 (wcsinit ignores PC)
//     or without PC1_1 / PC001001 (wcsinit reads PC only when that card
//     is present, else CROTA; wcslib applies the other PCi_j).
func tnxCheckLinear(h *header.Header, w *wcs.Header) error {
	for i := 1; i <= 2; i++ {
		for j := 1; j <= 2; j++ {
			old, mod := fmt.Sprintf("CD%03d%03d", i, j), fmt.Sprintf("CD%d_%d", i, j)
			if h.Has(old) && !h.Has(mod) {
				return fmt.Errorf("wcs/transform: TNX header has %s without %s; WCSTools reads only CDi_j, wcslib reads CD00i00j as CDi_j; not supported", old, mod)
			}
		}
	}
	if w.CD != nil {
		return nil
	}
	c1, c2 := h.Has("CDELT1"), h.Has("CDELT2")
	switch {
	case c1 && !c2:
		return errors.New("wcs/transform: TNX header has CDELT1 without CDELT2 and no CD matrix; WCSTools sets CDELT2 = -CDELT1, wcslib 1; not supported")
	case c2 && !c1:
		return errors.New("wcs/transform: TNX header has CDELT2 without CDELT1 and no CD matrix; WCSTools then ignores CDELT2 (1 deg/pixel), wcslib uses it; not supported")
	}
	if w.HasPC {
		if !c1 {
			return errors.New("wcs/transform: TNX header has a PC matrix without CDELT1 and no CD matrix; WCSTools ignores PC without CDELT1, wcslib applies it; not supported")
		}
		if !h.Has("PC1_1") && !h.Has("PC001001") {
			return errors.New("wcs/transform: TNX header has a PC matrix without PC1_1 (or PC001001); WCSTools then ignores the PC matrix, wcslib applies it; not supported")
		}
	}
	return nil
}

// Code returns "TNX".
func (t *TNXTransform) Code() string { return "TNX" }

// PixelToSky maps 0-based pixel coordinates to (alpha, delta) in degrees,
// alpha in [0, 360).
func (t *TNXTransform) PixelToSky(p1, p2 float64) (alpha, delta float64, err error) {
	b := t.base
	dx := p1 - b.crpixX
	dy := p2 - b.crpixY
	// lin is the matrix before CDELT; scale is CDELT * CUNIT in radians.
	x := float64((float64(b.lin[0][0]*dx) + float64(b.lin[0][1]*dy)) * b.scale[0])
	y := float64((float64(b.lin[1][0]*dx) + float64(b.lin[1][1]*dy)) * b.scale[1])
	// The TNX surfaces are defined on (xi, eta) in degrees.
	xd, yd := t.tnx.Forward(radToDeg(x), radToDeg(y))
	phi, theta, ok := b.proj.Inverse(degToRad(xd), degToRad(yd))
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: pixel (%g, %g) outside TNX projection domain", p1, p2)
	}
	a, d := nativeToCelestial(phi, theta, b.alphaP, b.sinDp, b.cosDp, b.phiP)
	return normalizeDeg(radToDeg(a)), radToDeg(d), nil
}

// SkyToPixel is the inverse of PixelToSky (Newton iteration on the TNX
// correction). Returns an error when the sky position is outside the TAN
// hemisphere or the iteration does not converge.
func (t *TNXTransform) SkyToPixel(alpha, delta float64) (p1, p2 float64, err error) {
	b := t.base
	phi, theta := celestialToNative(degToRad(alpha), degToRad(delta), b.alphaP, b.sinDp, b.cosDp, b.phiP)
	x, y, ok := b.proj.Forward(phi, theta)
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: sky (%g, %g) outside TNX projection domain", alpha, delta)
	}
	xd, yd, ok := t.tnx.Inverse(radToDeg(x), radToDeg(y))
	if !ok {
		return 0, 0, fmt.Errorf("wcs/transform: TNX inverse did not converge at sky (%g, %g)", alpha, delta)
	}
	x, y = degToRad(xd)/b.scale[0], degToRad(yd)/b.scale[1]
	dx := float64(b.linInv[0][0]*x) + float64(b.linInv[0][1]*y)
	dy := float64(b.linInv[1][0]*x) + float64(b.linInv[1][1]*y)
	return dx + b.crpixX, dy + b.crpixY, nil
}
