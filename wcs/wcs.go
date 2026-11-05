// Package wcs parses the FITS World Coordinate System keyword set defined
// in Greisen & Calabretta 2002 ("Representations of world coordinates in
// FITS", Paper I) and Calabretta & Greisen 2002 (Paper II) into a typed
// in-memory struct.
//
// This package is PARSING ONLY. It does not compute any projection math,
// celestial rotations, or coordinate conversions. The math lives in
// sibling package fits/wcs/transform, which is checked end to end (this
// parser included) against wcslib 8.10 on real and synthetic headers; see
// wcs/transform/testdata/wref/PROVENANCE.md.
//
// Keyword vocabulary read (each with the alternate-description suffix a
// for ParseAlt, where Paper I defines one):
//
//   - WCS dimension: WCSAXESa
//   - Per-axis: CTYPEia, CUNITia, CRVALia, CRPIXia, CDELTia, CROTAia
//   - Linear transform: PCi_ja, CDi_ja
//   - Projection parameters: PVi_ma, PSi_ma
//   - Celestial: LONPOLEa, LATPOLEa, RADESYSa (RADECSYS), EQUINOXa (EPOCH)
//   - Observation time: MJD-OBS, DATE-OBS
//   - SIP distortion (primary description only): A_ORDER, A_p_q, B_ORDER,
//     B_p_q, AP_ORDER, AP_p_q, BP_ORDER, BP_p_q
//   - DSS plate solution (primary description only): CNPIX1/2, PPO3,
//     PPO6, XPIXELSZ, YPIXELSZ, PLTRAH/M/S, PLTDECSN/D/M/S, AMDXm, AMDYm
//
// RADESYS, EQUINOX, MJD-OBS and DATE-OBS are parsed for the caller; the
// transform package does not use them (it does not convert between
// reference frames on its own).
//
// Axes are 1-based per FITS convention; the exported slices in Header are
// indexed 0-based (so CTYPE[0] corresponds to CTYPE1, CRVAL[1] to CRVAL2,
// etc.) because Go slices are 0-based.
package wcs

import (
	"errors"
	"fmt"
	"regexp"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
)

// PVKey identifies a PVi_m or PSi_m projection parameter: axis index i
// (1-based, matching the FITS keyword) and parameter index m.
type PVKey struct {
	Axis  int // i, 1-based
	Index int // m, 0-based per the FITS convention
}

// Header is the parsed WCS keyword set for one coordinate representation.
//
// The slices (CType, CUnit, CRVal, CRPix, CDelt, CRota, ProjCode,
// DistCode, AxisType) are indexed 0-based: CType[0] corresponds to the
// FITS keyword CTYPE1.
//
// PC and CD are row-major NAxis × NAxis matrices in the same 0-based
// indexing: PC[i][j] corresponds to PCi+1_j+1. PC is always populated
// (defaulting to the identity matrix). CD is nil unless the header has at
// least one CDi_j card; unspecified CDi_j then default to zero.
//
// Which matrix applies follows wcslib (wcs.h, wcsprm::altlin): PCi_j when
// any PCi_j card is present (HasPC), else CDi_j when any CDi_j card is
// present (CDELT and CROTA are then ignored), else CROTA of the latitude
// axis combined with CDELT (AIPS convention), else the identity PC with
// CDELT. A header carrying both PC and CD (forbidden by Paper I but seen
// in the wild) is accepted and PC wins, as in wcslib.
type Header struct {
	// NAxis is the number of WCS axes: WCSAXESa when present, otherwise
	// the larger of NAXIS and the highest axis index on any WCS keyword of
	// this description (Paper I §2.2). Keywords with an axis index above
	// NAxis are ignored.
	NAxis int

	// Per-axis keywords; each slice has length NAxis.
	CType []string
	CUnit []string
	CRVal []float64
	CRPix []float64
	CDelt []float64
	CRota []float64 // CROTAia in degrees (0 when absent)

	// Linear transform matrices. PC is always NAxis × NAxis (identity if
	// the header has no PC cards). CD is nil unless the header has CD
	// cards. HasPC / HasCRota record whether those keywords were present.
	PC       [][]float64
	CD       [][]float64
	HasPC    bool
	HasCRota bool

	// Projection parameters, keyed by (axis, m). A parameter absent from
	// the header is absent from the map.
	PV map[PVKey]float64
	PS map[PVKey]string

	// Global celestial keywords.
	LonPole float64 // LONPOLEa, degrees (valid when LonPoleSet)
	LatPole float64 // LATPOLEa, degrees (valid when LatPoleSet)
	RadeSys string  // RADESYSa (e.g. "ICRS", "FK5", "FK4"), or RADECSYS
	Equinox float64 // EQUINOXa, or EPOCH

	// Time of observation.
	DateObs string  // DATE-OBS
	MJDObs  float64 // MJD-OBS

	// Derived metadata, populated by Parse.

	// ProjCode holds the 3-letter projection code extracted from the tail
	// of each CTYPEi (e.g. "TAN" from "RA---TAN"). For axes whose CTYPE is
	// not of the "NNNN-PPP" form this is the empty string.
	ProjCode []string

	// DistCode holds the 3-letter distortion code appended to a celestial
	// CTYPEi per Paper IV (e.g. "SIP" from "RA---TAN-SIP"), or the empty
	// string when CTYPEi declares none.
	DistCode []string

	// SIP holds the SIP distortion coefficients of the PRIMARY
	// description whenever the header carries SIP keywords, with or
	// without the "-SIP" CTYPE suffix, as in wcslib 8.10 (wcspih makes
	// SIP the prior distortion from the A_p_q / B_p_q keywords alone),
	// with one exception:
	//
	//   - CTYPE RA---TPV / DEC--TPV (both celestial axes TPV, neither with
	//     "-SIP"): the SIP keywords are IGNORED (SIP nil, SIPIgnored set),
	//     as astropy 8.0.2 ignores them ("Removed redundant SIP distortion
	//     parameters because CTYPE explicitly specifies TPV distortions").
	//     Survey pipelines (ATLAS) write a SIP copy of the TPV correction
	//     next to it for SIP-only readers. This is a deliberate difference
	//     from wcslib 8.10, which applies both: on the real ATLAS header of
	//     PR #4 (wcs/transform/testdata/wref/hdrcases_differ/
	//     real_atlas_k2_tpv_sip.case) its SIP set equals its TPV
	//     correction to 1.3e-7 px ((-12.28, 6.06) px at the stamp corner),
	//     and wcslib's reading lands 13.3-14.6 px (24.5-26.8 arcsec) from
	//     the TPV-only one at every point of the stamp, the correction
	//     counted twice, while astropy's equals TPV alone to 4e-10 arcsec.
	//     The ignored cards are not read at all, so they are no error
	//     whatever they hold.
	//   - CTYPE RA---TPV-SIP (the suffix written out): the header declares
	//     both and gets both, SIP before the matrix and TPV after it, as
	//     wcslib applies them.
	//
	// "-SIP" with no SIP keyword at all is undistorted, as in wcslib. SIP
	// keywords without A_ORDER or B_ORDER are a Parse error, and so are SIP
	// keywords on a TNX or ZPX header (with or without "-SIP"): its WAT
	// strings carry a second distortion solution, and wcslib 8.10 keeps
	// one or the other by card order (SIP after WAT: SIP + TAN, the WAT
	// correction ignored), so we refuse it. SIP is nil for an
	// alternate description (SIP has no alternate form; wcslib leaves
	// alternates undistorted), and an alternate that declares "-SIP" is
	// listed in UnsupportedDistortion.
	SIP *SIP

	// SIPIgnored reports that the header carries SIP keywords that were
	// not read because the celestial CTYPEs declare TPV without "-SIP"
	// (see SIP).
	SIPIgnored bool

	// DSS is the Digitized Sky Survey plate solution of the PRIMARY
	// description when the header carries DSS plate keywords (see DSS),
	// nil otherwise. Parse then replaces CTYPE1/2, CRPIX, CRVAL, CDELT and
	// PC in this Header with the values wcslib 8.10 derives from the plate
	// solution (CD nil, CROTA unused), as wcslib does, and wcs/transform
	// applies the plate polynomial. Alternate descriptions are not
	// affected.
	DSS *DSS

	// UnsupportedDistortion lists distortion keywords present in the
	// header that this library does not implement (Paper IV CPDISja /
	// CQDISia lookup or polynomial records, HST D2IMDISa detector-to-image
	// tables, or SIP keywords referenced from an alternate description,
	// which SIP does not define). transform.New refuses a header with a
	// non-empty list rather than silently ignoring the distortion. A
	// caller that knowingly accepts the residual error may clear it.
	UnsupportedDistortion []string

	// AxisType holds the axis type of each CTYPEi: the part before the
	// projection code with the '-' padding removed (e.g. "RA", "DEC",
	// "GLON", "GLAT", "ELON", "ELAT"), or, for a CTYPE without a
	// projection code, the whole upper-cased value (e.g. "FREQ").
	AxisType []string

	// LonAxis and LatAxis are the 1-based indices of the longitude and
	// latitude axes, respectively, or 0 if a celestial axis pair was not
	// identified. Recognised pairs: RA/DEC, GLON/GLAT, ELON/ELAT,
	// SLON/SLAT, HLON/HLAT. More than one longitude (or latitude) axis
	// leaves both at 0 and records the problem in CelestialError.
	LonAxis int
	LatAxis int

	// CelestialError is non-empty when the CTYPEs name celestial axes that
	// do not form exactly one longitude/latitude pair.
	CelestialError string

	// LonPoleSet / LatPoleSet record whether the user supplied a value
	// explicitly; transform needs this to apply the default-by-projection
	// rules from Paper II §2.4 correctly.
	LonPoleSet bool
	LatPoleSet bool
}

// ErrConflictingMatrix is no longer returned. A header with both PCi_j and
// CDi_j is accepted and PCi_j takes precedence, as in wcslib. The value is
// kept so that existing references still compile.
var ErrConflictingMatrix = errors.New("fits/wcs: header contains both PC and CD matrices")

// Parse reads a *header.Header and returns the primary (unsuffixed) WCS
// description. Equivalent to ParseAlt(h, "").
func Parse(h *header.Header) (*Header, error) { return ParseAlt(h, "") }

// ParseAlt reads a *header.Header (the parsed FITS header of an image
// HDU) and returns a typed WCS Header for the given alternate WCS
// description. Alt must be the empty string (for the primary WCS) or a
// single uppercase letter A-Z (for alternate descriptions per Paper I
// §3.7). Every keyword that Paper I defines per description is searched
// with the alt suffix appended (e.g. CTYPE1A, CD1_1A, WCSAXESA,
// LONPOLEA); NAXIS, MJD-OBS, DATE-OBS, RADECSYS and EPOCH are shared.
//
// For the primary description the defunct keywords of the early Paper
// I/II drafts are read as wcslib 8.10 (WCSHDR_all) and astropy read them:
// PC00i00j and CD00i00j as PCi_j and CDi_j, PROJPn (n = 0..9) as PVi_n of
// the latitude axis. A defunct card whose modern twin is also present
// with a different value is an error, as is PROJPn without a latitude
// axis.
//
// Parse never panics on malformed input. A WCS keyword that is present
// with a value of the wrong type (e.g. CRVAL1 = 'abc'), or an invalid SIP
// order, is an error naming the keyword: the value is not silently
// replaced by a default. SIP keywords on a TNX / ZPX header are an error
// (see Header.SIP).
func ParseAlt(h *header.Header, alt string) (*Header, error) {
	if len(alt) > 1 {
		return nil, fmt.Errorf("fits/wcs: alt WCS identifier %q must be 0 or 1 character", alt)
	}
	if len(alt) == 1 {
		c := alt[0]
		if c < 'A' || c > 'Z' {
			return nil, fmt.Errorf("fits/wcs: alt WCS identifier %q must be A-Z", alt)
		}
	}
	return parseWithAlt(h, alt)
}

// getFloat returns (value, present, error) for a numeric keyword.
func getFloat(h *header.Header, key string) (float64, bool, error) {
	c, ok := h.Get(key)
	if !ok || c.Type == header.TypeEmpty {
		return 0, false, nil
	}
	if f, ok := floatValue(c.Value); ok {
		return f, true, nil
	}
	return 0, false, fmt.Errorf("fits/wcs: %s has a non-numeric value %v", key, c.Value)
}

// getString returns (value, present, error) for a string keyword.
func getString(h *header.Header, key string) (string, bool, error) {
	c, ok := h.Get(key)
	if !ok || c.Type == header.TypeEmpty {
		return "", false, nil
	}
	if s, ok := c.Value.(string); ok {
		return strings.TrimSpace(s), true, nil
	}
	return "", false, fmt.Errorf("fits/wcs: %s has a non-string value %v", key, c.Value)
}

// parseWithAlt is the worker that does the actual parse.
func parseWithAlt(h *header.Header, alt string) (*Header, error) {
	w := &Header{
		PV: make(map[PVKey]float64),
		PS: make(map[PVKey]string),
	}

	// Determine the number of WCS axes.
	n, err := detectNAxis(h, alt)
	if err != nil {
		return nil, err
	}
	w.NAxis = n
	if w.NAxis <= 0 {
		if alt == "" && HasDSSKeywords(h) {
			return nil, errors.New("fits/wcs: DSS plate solution keywords on a header with no WCS axes (NAXIS = 0)")
		}
		return w, nil
	}

	// Pre-allocate per-axis slices.
	w.CType = make([]string, w.NAxis)
	w.CUnit = make([]string, w.NAxis)
	w.CRVal = make([]float64, w.NAxis)
	w.CRPix = make([]float64, w.NAxis)
	w.CDelt = make([]float64, w.NAxis)
	w.CRota = make([]float64, w.NAxis)
	w.ProjCode = make([]string, w.NAxis)
	w.DistCode = make([]string, w.NAxis)
	w.AxisType = make([]string, w.NAxis)

	// Defaults per Paper I:
	//   CRPIXi = 0, CRVALi = 0, CDELTi = 1, CTYPEi = "" (linear), CUNITi = ""
	// PC defaults to the identity matrix.
	for i := 1; i <= w.NAxis; i++ {
		w.CDelt[i-1] = 1.0
	}

	// Per-axis scalar keywords, with the alt suffix appended.
	for i := 1; i <= w.NAxis; i++ {
		var err error
		str := func(base string, dst *string) {
			if err != nil {
				return
			}
			if v, ok, e := getString(h, ax(base, i)+alt); e != nil {
				err = e
			} else if ok {
				*dst = v
			}
		}
		num := func(base string, dst *float64) bool {
			if err != nil {
				return false
			}
			v, ok, e := getFloat(h, ax(base, i)+alt)
			if e != nil {
				err = e
				return false
			}
			if ok {
				*dst = v
			}
			return ok
		}
		str("CTYPE", &w.CType[i-1])
		str("CUNIT", &w.CUnit[i-1])
		num("CRVAL", &w.CRVal[i-1])
		num("CRPIX", &w.CRPix[i-1])
		num("CDELT", &w.CDelt[i-1])
		if num("CROTA", &w.CRota[i-1]) {
			w.HasCRota = true
		}
		if err != nil {
			return nil, err
		}
	}

	// PC / CD / PV / PS scan. Each keyword has its alt suffix appended to
	// the "i_j" or "i_m" part — so PC1_1A matches for alt="A", PC1_1 for
	// alt="".
	var haveCD bool
	pc := identityMatrix(w.NAxis)
	var cd [][]float64
	var legacy []legacyCard
	modernPC, modernCD := map[[2]int]float64{}, map[[2]int]float64{}
	for _, c := range h.Cards() {
		if c.Key == "" || c.Type == header.TypeEmpty {
			continue
		}
		key := stripAltSuffix(c.Key, alt)
		if key == "" {
			continue // didn't match the requested alt
		}
		if alt == "" {
			// Defunct forms from the early drafts of WCS Papers I/II, which
			// wcslib 8.10 wcspih (WCSHDR_PC00i00j / CD00i00j / PROJPn, all
			// in WCSHDR_all) and astropy (relax=True) read for the primary
			// description: PC00i00j and CD00i00j are PCi_j and CDi_j,
			// PROJPn (n <= 9) is PVi_n of the latitude axis.
			if p, i, j, ok := parseLegacyMatrixKey(key); ok {
				f, ok := floatValue(c.Value)
				if !ok {
					return nil, fmt.Errorf("fits/wcs: %s has a non-numeric value %v", c.Key, c.Value)
				}
				legacy = append(legacy, legacyCard{key: c.Key, p: p, i: i, j: j, v: f})
				continue
			}
			if n, ok := parseProjPKey(key); ok {
				f, ok := floatValue(c.Value)
				if !ok {
					return nil, fmt.Errorf("fits/wcs: %s has a non-numeric value %v", c.Key, c.Value)
				}
				legacy = append(legacy, legacyCard{key: c.Key, p: "PROJP", j: n, v: f})
				continue
			}
		}
		for _, p := range []string{"PC", "CD", "PV", "PS"} {
			i, j, ok := parseMatrixIndex(key, p)
			if !ok {
				continue
			}
			if p == "PS" {
				s, ok := c.Value.(string)
				if !ok {
					return nil, fmt.Errorf("fits/wcs: %s has a non-string value %v", c.Key, c.Value)
				}
				w.PS[PVKey{Axis: i, Index: j}] = s
				continue
			}
			f, ok := floatValue(c.Value)
			if !ok {
				return nil, fmt.Errorf("fits/wcs: %s has a non-numeric value %v", c.Key, c.Value)
			}
			switch p {
			case "PV":
				w.PV[PVKey{Axis: i, Index: j}] = f
			case "PC":
				modernPC[[2]int{i, j}] = f
				if i <= w.NAxis && j >= 1 && j <= w.NAxis {
					pc[i-1][j-1] = f
					w.HasPC = true
				}
			case "CD":
				modernCD[[2]int{i, j}] = f
				if i <= w.NAxis && j >= 1 && j <= w.NAxis {
					if cd == nil {
						cd = zeroMatrix(w.NAxis)
					}
					cd[i-1][j-1] = f
					haveCD = true
				}
			}
		}
	}
	var projP []legacyCard
	for _, l := range legacy {
		switch l.p {
		case "PROJP":
			projP = append(projP, l) // needs the latitude axis, below
		case "PC", "CD":
			modern := modernPC
			if l.p == "CD" {
				modern = modernCD
			}
			if v, dup := modern[[2]int{l.i, l.j}]; dup && v != l.v {
				return nil, fmt.Errorf("fits/wcs: %s = %v and %s%d_%d = %v disagree", l.key, l.v, l.p, l.i, l.j, v)
			}
			if l.i > w.NAxis || l.j < 1 || l.j > w.NAxis {
				continue
			}
			if l.p == "PC" {
				pc[l.i-1][l.j-1] = l.v
				w.HasPC = true
			} else {
				if cd == nil {
					cd = zeroMatrix(w.NAxis)
				}
				cd[l.i-1][l.j-1] = l.v
				haveCD = true
			}
		}
	}
	w.PC = pc
	if haveCD {
		w.CD = cd
	}

	// Global celestial keywords (with alt suffix where applicable).
	if v, ok, err := getFloat(h, "LONPOLE"+alt); err != nil {
		return nil, err
	} else if ok {
		w.LonPole, w.LonPoleSet = v, true
	}
	if v, ok, err := getFloat(h, "LATPOLE"+alt); err != nil {
		return nil, err
	} else if ok {
		w.LatPole, w.LatPoleSet = v, true
	}
	if v, ok, err := getString(h, "RADESYS"+alt); err != nil {
		return nil, err
	} else if ok {
		w.RadeSys = v
	} else if v, ok, err := getString(h, "RADECSYS"); err != nil {
		return nil, err
	} else if ok {
		// RADECSYS is the pre-1998 spelling (no alt variant was defined).
		w.RadeSys = v
	}
	if v, ok, err := getFloat(h, "EQUINOX"+alt); err != nil {
		return nil, err
	} else if ok {
		w.Equinox = v
	} else if v, ok, err := getFloat(h, "EPOCH"); err != nil {
		return nil, err
	} else if ok {
		// EPOCH is the pre-1988 spelling for EQUINOX (no alt variant).
		w.Equinox = v
	}
	if v, ok, err := getString(h, "DATE-OBS"); err != nil {
		return nil, err
	} else if ok {
		w.DateObs = v
	}
	if v, ok, err := getFloat(h, "MJD-OBS"); err != nil {
		return nil, err
	} else if ok {
		w.MJDObs = v
	}

	// DSS plate solution (primary only, as in wcslib 8.10 wcspih_final):
	// it replaces the header's linear WCS before anything is derived from
	// it.
	if alt == "" {
		d, err := parseDSS(h)
		if err != nil {
			return nil, err
		}
		if d != nil {
			if HasSIPKeywords(h) {
				return nil, errors.New("fits/wcs: header has both DSS plate-solution and SIP keywords (wcslib 8.10 refuses it: one translated distortion per header)")
			}
			if w.NAxis != 2 {
				return nil, fmt.Errorf("fits/wcs: DSS plate solution on a %d-axis WCS; it is defined for a 2-axis plate image only", w.NAxis)
			}
			w.applyDSS(d)
		}
	}

	// Derive axis type, projection and distortion codes from each CTYPE.
	for i := range w.CType {
		w.AxisType[i], w.ProjCode[i], w.DistCode[i] = splitCType(w.CType[i])
	}
	// Identify the longitude/latitude axis pair (1-based indices).
	var lons, lats []int
	for i, a := range w.AxisType {
		switch a {
		case "RA", "GLON", "ELON", "SLON", "HLON":
			lons = append(lons, i+1)
		case "DEC", "GLAT", "ELAT", "SLAT", "HLAT":
			lats = append(lats, i+1)
		}
	}
	switch {
	case len(lons) == 1 && len(lats) == 1:
		w.LonAxis, w.LatAxis = lons[0], lats[0]
	case len(lons)+len(lats) > 0:
		w.CelestialError = fmt.Sprintf("celestial CTYPEs %q do not form one longitude/latitude pair", w.CType)
	}
	if len(projP) > 0 {
		if w.LatAxis == 0 {
			return nil, fmt.Errorf("fits/wcs: %s needs a celestial latitude axis (PROJPn is PVi_n of that axis)", projP[0].key)
		}
		for _, l := range projP {
			k := PVKey{Axis: w.LatAxis, Index: l.j}
			if v, dup := w.PV[k]; dup && v != l.v {
				return nil, fmt.Errorf("fits/wcs: %s = %v and PV%d_%d = %v disagree", l.key, l.v, k.Axis, k.Index, v)
			}
			w.PV[k] = l.v
		}
	}

	// Distortion keywords we do not implement: refuse rather than ignore.
	for _, c := range h.Cards() {
		if m := unsupportedDistRe.FindStringSubmatch(c.Key); m != nil && m[3] == alt {
			w.UnsupportedDistortion = append(w.UnsupportedDistortion, c.Key)
		}
	}

	// SIP, as wcslib 8.10 wcspih reads it: the A_p_q / B_p_q keywords
	// alone make SIP the prior distortion of the PRIMARY description,
	// whatever its CTYPEs say (no "-SIP" suffix is needed) - except on a
	// RA---TPV / DEC--TPV header without "-SIP", whose SIP set is ignored
	// as astropy ignores it (see Header.SIP). SIP has no
	// alternate-description form, so an alternate description is never
	// distorted by it; one that declares "-SIP" is refused.
	if HasSIPKeywords(h) && alt == "" && w.tpvWithoutSIPSuffix() {
		w.SIPIgnored = true
	} else if HasSIPKeywords(h) {
		// An IRAF TNX / ZPX header carries its own distortion in the WAT
		// strings. wcslib 8.10 keeps whichever of SIP and WAT its lexer met
		// last (SIP cards after the WAT cards: SIP with TAN, the WAT
		// correction ignored; the other order: wcspih fails), so the
		// header has two solutions and an order-dependent meaning. Refused
		// (hdrcases_stricter/tnx_sip_*).
		if alt == "" && w.IsCelestial() {
			for _, ax := range []int{w.LonAxis, w.LatAxis} {
				if c := w.ProjCode[ax-1]; c == "TNX" || c == "ZPX" {
					return nil, fmt.Errorf("fits/wcs: %s header carries SIP keywords (CTYPE%d = %q): two distortion solutions, refused", c, ax, w.CType[ax-1])
				}
			}
		}
		if alt == "" {
			sip, err := ParseSIP(h)
			if err != nil {
				return nil, err
			}
			w.SIP = sip
		} else if w.IsCelestial() && w.DistCode[w.LonAxis-1] == "SIP" {
			w.UnsupportedDistortion = append(w.UnsupportedDistortion,
				fmt.Sprintf("SIP (A_ORDER/B_ORDER have no alternate-description form; CTYPE%d%s declares -SIP)", w.LonAxis, alt))
		}
	}

	return w, nil
}

// tpvWithoutSIPSuffix reports whether both celestial axes declare TPV
// and neither declares the "-SIP" distortion code: the form whose SIP
// keywords astropy 8.0.2 drops (wcs.WCS._fix_pre2012_scamp_tpv: both
// CTYPEs end in "-TPV").
func (w *Header) tpvWithoutSIPSuffix() bool {
	if !w.IsCelestial() {
		return false
	}
	for _, ax := range []int{w.LonAxis, w.LatAxis} {
		if w.ProjCode[ax-1] != "TPV" || w.DistCode[ax-1] != "" {
			return false
		}
	}
	return true
}

// unsupportedDistRe matches Paper IV / HST distortion keywords:
// CPDISja, CQDISia, D2IMDISa (group 3 is the alt letter, "" for primary).
var unsupportedDistRe = regexp.MustCompile(`^(CPDIS|CQDIS|D2IMDIS)([0-9]+)([A-Z]?)$`)

// IsCelestial reports whether the header describes a celestial coordinate
// pair (longitude + latitude axis both identified).
func (w *Header) IsCelestial() bool {
	return w.LonAxis > 0 && w.LatAxis > 0
}

// CelestialProjCode returns the projection code of the celestial
// longitude axis. It panics if IsCelestial() is false. It does not check
// that the latitude axis declares the same code; transform.New rejects a
// header whose two celestial axes disagree, as wcslib does.
func (w *Header) CelestialProjCode() string {
	if !w.IsCelestial() {
		panic("wcs: CelestialProjCode called on non-celestial header")
	}
	return w.ProjCode[w.LonAxis-1]
}

// ax formats a per-axis keyword name "BASE" + i.
func ax(base string, i int) string { return base + strconv.Itoa(i) }

// stripAltSuffix returns the keyword name with the alt suffix removed if
// it matches, or "" if the key does not belong to the requested alt
// description. For alt="", keys ending in a letter A-Z are rejected (they
// belong to alternate descriptions); keys without a trailing letter are
// returned unchanged.
func stripAltSuffix(key, alt string) string {
	// Only keys matching known WCS prefixes can carry an alt suffix.
	prefixes := []string{"PC", "CD", "PV", "PS"}
	matched := false
	for _, p := range prefixes {
		if strings.HasPrefix(key, p) {
			matched = true
			break
		}
	}
	if !matched {
		return key
	}
	lastChar := key[len(key)-1]
	keyHasAlt := lastChar >= 'A' && lastChar <= 'Z'
	if alt == "" {
		if keyHasAlt {
			return "" // primary WCS skips alt-suffixed keys
		}
		return key
	}
	if !keyHasAlt || lastChar != alt[0] {
		return "" // wrong alt
	}
	return key[:len(key)-1]
}

// axisKeyRe / matrixKeyRe / paramKeyRe match the keywords whose axis
// index counts toward the WCS dimension (the last group is the alt
// suffix).
var (
	axisKeyRe   = regexp.MustCompile(`^(?:CTYPE|CUNIT|CRVAL|CRPIX|CDELT|CROTA)([0-9]+)([A-Z]?)$`)
	matrixKeyRe = regexp.MustCompile(`^(?:PC|CD)([0-9]+)_([0-9]+)([A-Z]?)$`)
	paramKeyRe  = regexp.MustCompile(`^(?:PV|PS)([0-9]+)_([0-9]+)([A-Z]?)$`)
)

// detectNAxis determines the WCS dimension of description alt: WCSAXESa
// when present; otherwise the larger of NAXIS and the highest axis index
// on any CTYPE/CUNIT/CRVAL/CRPIX/CDELT/CROTA, PCi_j/CDi_j (both indices)
// or PVi_m/PSi_m (axis index) keyword of that description (Paper I §2.2).
func detectNAxis(h *header.Header, alt string) (int, error) {
	if c, ok := h.Get("WCSAXES" + alt); ok && c.Type != header.TypeEmpty {
		v, ok := c.Value.(int64)
		if !ok || v < 0 || v > 99 {
			return 0, fmt.Errorf("fits/wcs: WCSAXES%s = %v is not an integer in 0..99", alt, c.Value)
		}
		return int(v), nil
	}
	max := 0
	if c, ok := h.Get("NAXIS"); ok {
		if v, ok := c.Value.(int64); ok && v > 0 && v <= 999 {
			max = int(v)
		}
	}
	bump := func(s string) {
		if n, err := strconv.Atoi(s); err == nil && n > max && n <= 99 {
			max = n
		}
	}
	for _, c := range h.Cards() {
		if m := axisKeyRe.FindStringSubmatch(c.Key); m != nil && m[2] == alt {
			bump(m[1])
		} else if m := matrixKeyRe.FindStringSubmatch(c.Key); m != nil && m[3] == alt {
			bump(m[1])
			bump(m[2])
		} else if m := paramKeyRe.FindStringSubmatch(c.Key); m != nil && m[3] == alt {
			bump(m[1])
		}
	}
	return max, nil
}

// legacyCard is a defunct-form WCS card (PC00i00j, CD00i00j, PROJPn).
type legacyCard struct {
	key  string
	p    string // "PC", "CD" or "PROJP"
	i, j int    // matrix indices, or j = n for PROJPn
	v    float64
}

// parseLegacyMatrixKey matches "PC" or "CD" followed by two 3-digit
// indices (PC001002 is PC1_2), the form wcslib's WCSHDR_PC00i00j /
// WCSHDR_CD00i00j accept.
func parseLegacyMatrixKey(key string) (string, int, int, bool) {
	if len(key) != 8 || (key[:2] != "PC" && key[:2] != "CD") {
		return "", 0, 0, false
	}
	for _, r := range key[2:] {
		if r < '0' || r > '9' {
			return "", 0, 0, false
		}
	}
	i, _ := strconv.Atoi(key[2:5])
	j, _ := strconv.Atoi(key[5:8])
	if i < 1 || j < 1 {
		return "", 0, 0, false
	}
	return key[:2], i, j, true
}

// parseProjPKey matches PROJPn with one digit n (wcslib WCSHDR_PROJPn).
func parseProjPKey(key string) (int, bool) {
	if len(key) != 6 || key[:5] != "PROJP" || key[5] < '0' || key[5] > '9' {
		return 0, false
	}
	return int(key[5] - '0'), true
}

// parseMatrixIndex parses a keyword of the form "PREFIX<i>_<j>" and returns
// (i, j, true) if the pattern matches. For PV/PS keywords the separator is
// still "_" per Paper I.
func parseMatrixIndex(key, prefix string) (int, int, bool) {
	if !strings.HasPrefix(key, prefix) {
		return 0, 0, false
	}
	rest := key[len(prefix):]
	us := strings.IndexByte(rest, '_')
	if us <= 0 || us == len(rest)-1 {
		return 0, 0, false
	}
	i, err1 := strconv.Atoi(rest[:us])
	j, err2 := strconv.Atoi(rest[us+1:])
	if err1 != nil || err2 != nil || i < 1 || j < 0 {
		return 0, 0, false
	}
	return i, j, true
}

// identityMatrix returns an n×n identity matrix.
func identityMatrix(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
		m[i][i] = 1
	}
	return m
}

// zeroMatrix returns an n×n zero matrix.
func zeroMatrix(n int) [][]float64 {
	m := make([][]float64, n)
	for i := range m {
		m[i] = make([]float64, n)
	}
	return m
}

// floatValue extracts a float64 from a Card.Value. Integer values are
// promoted to float64. Anything else returns ok=false.
func floatValue(v any) (float64, bool) {
	switch x := v.(type) {
	case float64:
		return x, true
	case int64:
		return float64(x), true
	}
	return 0, false
}

// splitCType parses a CTYPE value into an axis type, projection code and
// distortion code.
//
// CTYPEi for celestial axes follows the pattern "NNNNN-TTT" where NNNNN is
// a 1-to-4-character axis type code (e.g. "RA", "DEC", "GLON") padded on
// the right with '-' characters to column 5, followed by the 3-character
// projection code. Examples: "RA---TAN", "DEC--TAN", "GLON-CAR". Paper IV
// allows a further "-DDD" distortion code, e.g. "RA---TAN-SIP".
//
// For non-celestial or free-form CTYPE values the entire (upper-cased)
// string becomes the axis type and ProjCode and DistCode are empty.
func splitCType(ctype string) (axisType, projCode, distCode string) {
	ctype = strings.TrimSpace(ctype)
	if len(ctype) >= 12 && ctype[len(ctype)-4] == '-' && ctype[len(ctype)-8] == '-' {
		distCode = strings.ToUpper(ctype[len(ctype)-3:])
		ctype = ctype[:len(ctype)-4]
	}
	if len(ctype) < 5 {
		return strings.ToUpper(ctype), "", ""
	}
	if len(ctype) >= 8 && ctype[len(ctype)-4] == '-' {
		prefix := strings.TrimRight(ctype[:len(ctype)-3], "-")
		projCode = ctype[len(ctype)-3:]
		return strings.ToUpper(prefix), strings.ToUpper(projCode), distCode
	}
	return strings.ToUpper(ctype), "", ""
}

// String returns a short human-readable summary used in diagnostic output.
func (w *Header) String() string {
	return fmt.Sprintf("wcs.Header{NAxis=%d Lon=%d Lat=%d Proj=%v RadeSys=%s Equinox=%g}",
		w.NAxis, w.LonAxis, w.LatAxis, w.ProjCode, w.RadeSys, w.Equinox)
}
