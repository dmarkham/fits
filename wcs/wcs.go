// Package wcs parses the FITS World Coordinate System keyword set defined
// in Greisen & Calabretta 2002 ("Representations of world coordinates in
// FITS", Paper I) into a typed in-memory struct.
//
// This package is PARSING ONLY. It does not compute any projection math,
// celestial rotations, or coordinate conversions. The math lives in
// sibling package fits/wcs/transform.
//
// Supported keyword vocabulary:
//
//   - Per-axis: CTYPEi, CUNITi, CRVALi, CRPIXi, CDELTi
//   - Linear transform: PCi_j, CDi_j (mutually exclusive per Paper I §2.1)
//   - Projection parameters: PVi_m, PSi_m
//   - Celestial: LONPOLE, LATPOLE, RADESYS, EQUINOX, EPOCH
//   - Observation time: MJD-OBS, DATE-OBS
//
// Axes are 1-based per FITS convention; the exported slices in Header are
// indexed 0-based (so CTYPE[0] corresponds to CTYPE1, CRVAL[1] to CRVAL2,
// etc.) because Go slices are 0-based. This trade-off matches
// astropy.wcs.WCS and is documented on each field.
package wcs

import (
	"errors"
	"fmt"
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
// The slices (CType, CUnit, CRVal, CRPix, CDelt, ProjCode, AxisType) are
// indexed 0-based: CType[0] corresponds to the FITS keyword CTYPE1.
//
// PC and CD are row-major NAxis × NAxis matrices in the same 0-based
// indexing: PC[i][j] corresponds to PCi+1_j+1. PC is always populated
// (defaulting to the identity matrix). CD is nil unless the header used
// the CDi_j form. The two forms are mutually exclusive per Paper I §2.1;
// if both are present in the header, Parse returns an error.
type Header struct {
	// NAxis is the number of WCS axes. It is the larger of WCSAXES (if
	// present) and the maximum axis index found on any per-axis keyword,
	// and falls back to the image's NAXIS if neither is set.
	NAxis int

	// Per-axis keywords; each slice has length NAxis.
	CType []string
	CUnit []string
	CRVal []float64
	CRPix []float64
	CDelt []float64

	// Linear transform matrices. PC is always NAxis × NAxis (identity if
	// the header has neither PC nor CD). CD is nil unless the header used
	// the CDi_j form.
	PC [][]float64
	CD [][]float64

	// Projection parameters. Absent entries are returned as zero value
	// with ok=false via Has/Get helpers below.
	PV map[PVKey]float64
	PS map[PVKey]string

	// Global celestial keywords.
	LonPole float64 // LONPOLE, degrees
	LatPole float64 // LATPOLE, degrees
	RadeSys string  // RADESYS (e.g. "ICRS", "FK5", "FK4")
	Equinox float64 // EQUINOX

	// Time of observation.
	DateObs string  // DATE-OBS
	MJDObs  float64 // MJD-OBS

	// Derived metadata, populated by Parse.

	// ProjCode holds the 3-letter projection code extracted from the tail
	// of each CTYPEi (e.g. "TAN" from "RA---TAN"). For non-celestial axes
	// this is the empty string.
	ProjCode []string

	// AxisType holds the axis kind extracted from the leading portion of
	// CTYPEi (e.g. "RA", "DEC", "GLON", "GLAT", "ELON", "ELAT"). Empty for
	// non-celestial axes.
	AxisType []string

	// LonAxis and LatAxis are the 1-based indices of the longitude and
	// latitude axes, respectively, or 0 if a celestial axis pair was not
	// identified.
	LonAxis int
	LatAxis int

	// lonPoleSet / latPoleSet record whether the user supplied a value
	// explicitly; transform needs this to apply the default-by-projection
	// rules from Paper II §2.4 correctly.
	LonPoleSet bool
	LatPoleSet bool
}

// ErrConflictingMatrix is returned when the header contains both PCi_j and
// CDi_j cards, which is forbidden by Paper I §2.1.
var ErrConflictingMatrix = errors.New("fits/wcs: header contains both PC and CD matrices")

// Parse reads a *header.Header and returns the primary (unsuffixed) WCS
// description. Equivalent to ParseAlt(h, "").
func Parse(h *header.Header) (*Header, error) { return ParseAlt(h, "") }

// ParseAlt reads a *header.Header (the parsed FITS header of an image
// HDU) and returns a typed WCS Header for the given alternate WCS
// description. Alt must be the empty string (for the primary WCS) or a
// single uppercase letter A-Z (for alternate descriptions per Paper I
// §3.7). Every per-axis and matrix keyword is searched with the alt
// suffix appended (e.g. CTYPE1A, CD1_1A) and falls back to no suffix
// only for WCSAXES and NAXIS (which are not duplicated per description).
//
// Parse never panics on malformed input; keywords that fail to parse as
// their expected type are silently replaced with their defaults, matching
// the behavior of astropy.wcs.WCS. Structural contradictions (PC+CD both
// present) return an error.
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

// parseWithAlt is the worker that does the actual parse.
func parseWithAlt(h *header.Header, alt string) (*Header, error) {
	w := &Header{
		PV: make(map[PVKey]float64),
		PS: make(map[PVKey]string),
	}

	// Determine the number of WCS axes.
	w.NAxis = detectNAxis(h)
	if w.NAxis <= 0 {
		return w, nil
	}

	// Pre-allocate per-axis slices.
	w.CType = make([]string, w.NAxis)
	w.CUnit = make([]string, w.NAxis)
	w.CRVal = make([]float64, w.NAxis)
	w.CRPix = make([]float64, w.NAxis)
	w.CDelt = make([]float64, w.NAxis)
	w.ProjCode = make([]string, w.NAxis)
	w.AxisType = make([]string, w.NAxis)

	// Defaults per Paper I:
	//   CRPIXi = 0, CRVALi = 0, CDELTi = 1, CTYPEi = "" (linear), CUNITi = ""
	// PC defaults to the identity matrix.
	for i := 1; i <= w.NAxis; i++ {
		w.CDelt[i-1] = 1.0
	}

	// Per-axis scalar keywords, with the alt suffix appended.
	for i := 1; i <= w.NAxis; i++ {
		if v, err := h.String(ax("CTYPE", i) + alt); err == nil {
			w.CType[i-1] = strings.TrimSpace(v)
		}
		if v, err := h.String(ax("CUNIT", i) + alt); err == nil {
			w.CUnit[i-1] = strings.TrimSpace(v)
		}
		if v, err := h.Float(ax("CRVAL", i) + alt); err == nil {
			w.CRVal[i-1] = v
		}
		if v, err := h.Float(ax("CRPIX", i) + alt); err == nil {
			w.CRPix[i-1] = v
		}
		if v, err := h.Float(ax("CDELT", i) + alt); err == nil {
			w.CDelt[i-1] = v
		}
	}

	// PC / CD / PV / PS scan. Each keyword has its alt suffix appended to
	// the "i_j" or "i_m" part — so PC1_1A matches for alt="A", PC1_1 for
	// alt="".
	var havePC, haveCD bool
	pc := identityMatrix(w.NAxis)
	var cd [][]float64
	for _, c := range h.Cards() {
		if c.Key == "" {
			continue
		}
		key := stripAltSuffix(c.Key, alt)
		if key == "" {
			continue // didn't match the requested alt
		}
		// PCi_j
		if strings.HasPrefix(key, "PC") {
			if i, j, ok := parseMatrixIndex(key, "PC"); ok {
				if i >= 1 && i <= w.NAxis && j >= 1 && j <= w.NAxis {
					if f, ok := floatValue(c.Value); ok {
						pc[i-1][j-1] = f
						havePC = true
					}
				}
			}
		}
		// CDi_j
		if strings.HasPrefix(key, "CD") {
			if i, j, ok := parseMatrixIndex(key, "CD"); ok {
				if i >= 1 && i <= w.NAxis && j >= 1 && j <= w.NAxis {
					if cd == nil {
						cd = zeroMatrix(w.NAxis)
					}
					if f, ok := floatValue(c.Value); ok {
						cd[i-1][j-1] = f
						haveCD = true
					}
				}
			}
		}
		// PVi_m
		if strings.HasPrefix(key, "PV") {
			if i, m, ok := parseMatrixIndex(key, "PV"); ok {
				if f, ok := floatValue(c.Value); ok {
					w.PV[PVKey{Axis: i, Index: m}] = f
				}
			}
		}
		// PSi_m
		if strings.HasPrefix(key, "PS") {
			if i, m, ok := parseMatrixIndex(key, "PS"); ok {
				if s, ok := c.Value.(string); ok {
					w.PS[PVKey{Axis: i, Index: m}] = s
				}
			}
		}
	}
	if havePC && haveCD {
		return nil, ErrConflictingMatrix
	}
	w.PC = pc
	if haveCD {
		w.CD = cd
	}

	// Global celestial keywords (with alt suffix where applicable).
	if v, err := h.Float("LONPOLE" + alt); err == nil {
		w.LonPole = v
		w.LonPoleSet = true
	}
	if v, err := h.Float("LATPOLE" + alt); err == nil {
		w.LatPole = v
		w.LatPoleSet = true
	}
	if v, err := h.String("RADESYS" + alt); err == nil {
		w.RadeSys = strings.TrimSpace(v)
	} else if v, err := h.String("RADECSYS"); err == nil {
		// RADECSYS is the pre-1998 spelling (no alt variant was defined).
		w.RadeSys = strings.TrimSpace(v)
	}
	if v, err := h.Float("EQUINOX" + alt); err == nil {
		w.Equinox = v
	} else if v, err := h.Float("EPOCH"); err == nil {
		// EPOCH is the pre-1988 spelling for EQUINOX (no alt variant).
		w.Equinox = v
	}
	if v, err := h.String("DATE-OBS"); err == nil {
		w.DateObs = strings.TrimSpace(v)
	}
	if v, err := h.Float("MJD-OBS"); err == nil {
		w.MJDObs = v
	}

	// Derive projection code and axis type from each CTYPE.
	for i := range w.CType {
		w.AxisType[i], w.ProjCode[i] = splitCType(w.CType[i])
	}
	// Identify the longitude/latitude axis pair (1-based indices).
	for i, a := range w.AxisType {
		switch a {
		case "RA", "GLON", "ELON", "SLON", "HLON":
			w.LonAxis = i + 1
		case "DEC", "GLAT", "ELAT", "SLAT", "HLAT":
			w.LatAxis = i + 1
		}
	}

	return w, nil
}

// IsCelestial reports whether the header describes a celestial coordinate
// pair (longitude + latitude axis both identified).
func (w *Header) IsCelestial() bool {
	return w.LonAxis > 0 && w.LatAxis > 0
}

// CelestialProjCode returns the common projection code for the celestial
// axis pair. It panics if IsCelestial() is false. If the longitude and
// latitude axes declare different projection codes (which is malformed
// input but seen in the wild), the longitude axis wins.
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

// detectNAxis determines the WCS axis count from WCSAXES, then NAXIS, then
// the maximum observed axis index on any CTYPEi keyword.
func detectNAxis(h *header.Header) int {
	if v, err := h.Int("WCSAXES"); err == nil && v > 0 {
		return int(v)
	}
	if v, err := h.Int("NAXIS"); err == nil && v > 0 {
		return int(v)
	}
	// Scan for the highest CTYPEi seen.
	max := 0
	for _, c := range h.Cards() {
		if !strings.HasPrefix(c.Key, "CTYPE") {
			continue
		}
		rest := c.Key[len("CTYPE"):]
		if n, err := strconv.Atoi(rest); err == nil && n > max {
			max = n
		}
	}
	return max
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
	if err1 != nil || err2 != nil || i < 1 {
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

// splitCType parses a CTYPE value into an axis type and projection code.
//
// CTYPEi for celestial axes follows the pattern "NNNNN-TTT" where NNNNN is
// a 1-to-4-character axis type code (e.g. "RA", "DEC", "GLON") padded on
// the right with '-' characters to column 5, followed by the 3-character
// projection code. Examples: "RA---TAN", "DEC--TAN", "GLON-CAR".
//
// For non-celestial or free-form CTYPE values the entire string becomes
// the axis type and ProjCode is empty.
func splitCType(ctype string) (axisType, projCode string) {
	ctype = strings.TrimSpace(ctype)
	if len(ctype) < 5 {
		return ctype, ""
	}
	// The hyphen-padded axis type ends at the first non-hyphen character of
	// the trailing 3-character projection code. Equivalently: strip trailing
	// characters after the last run of hyphens.
	// Safe bet: if position [len-4] == '-', it's a celestial CTYPE.
	if len(ctype) >= 8 && ctype[len(ctype)-4] == '-' {
		prefix := strings.TrimRight(ctype[:len(ctype)-3], "-")
		projCode = ctype[len(ctype)-3:]
		return strings.ToUpper(prefix), strings.ToUpper(projCode)
	}
	return strings.ToUpper(ctype), ""
}

// String returns a short human-readable summary used in diagnostic output.
func (w *Header) String() string {
	return fmt.Sprintf("wcs.Header{NAxis=%d Lon=%d Lat=%d Proj=%v RadeSys=%s Equinox=%g}",
		w.NAxis, w.LonAxis, w.LatAxis, w.ProjCode, w.RadeSys, w.Equinox)
}
