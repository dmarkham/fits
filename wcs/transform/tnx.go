package transform

import (
	"fmt"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
)

// TNX — IRAF polynomial distortion.
//
// TNX stores its distortion coefficients inside free-form text strings
// under the WAT0_001, WAT1_001, WAT2_001 (and sometimes WATi_002 ...)
// header keywords. The text is a Forth-like grammar listing the
// polynomial surface type, its domain, and its coefficients. IRAF emits
// TNX headers from MSCRED, DAOFIND, and some 2MASS pipelines; legacy but
// still encountered in survey archives.
//
// Format of the concatenated WATi_nnn strings (example, axis 1):
//
//	"wtype=tnx axtype=ra projection=tan lngcor = \"3. 4. 4. 2. ... <coeffs>\""
//
// The "3." is the surface function type:
//
//	1 = Chebyshev
//	2 = Legendre
//	3 = simple polynomial (power series)
//
// Followed by xorder, yorder, cross-term flag, x domain, y domain, and
// the flat coefficient list. The polynomial is evaluated in the same
// intermediate world coordinate space as TPV, and the forward/inverse
// pipeline mirrors TPV exactly.

// TNX holds a parsed TNX distortion for both axes.
type TNX struct {
	Lon tnxSurface
	Lat tnxSurface
}

// tnxSurface is one axis's polynomial surface.
type tnxSurface struct {
	Type    int // 1=Chebyshev, 2=Legendre, 3=simple polynomial
	XOrder  int
	YOrder  int
	Cross   int // 0 = none, 1 = full, 2 = half (diagonal only)
	XMin    float64
	XMax    float64
	YMin    float64
	YMax    float64
	Coeffs  []float64 // row-major, length depends on Cross
	present bool
}

// ParseTNX reads TNX coefficients from a *header.Header. Returns nil if
// no TNX keywords are present (WAT1_* or WAT2_* absent).
func ParseTNX(h *header.Header) (*TNX, error) {
	wat1 := concatWAT(h, "WAT1")
	wat2 := concatWAT(h, "WAT2")
	if wat1 == "" && wat2 == "" {
		return nil, nil
	}
	t := &TNX{}
	if wat1 != "" {
		s, err := parseTNXWat(wat1, "lngcor")
		if err != nil {
			return nil, fmt.Errorf("wcs/transform: TNX axis 1: %w", err)
		}
		t.Lon = s
	}
	if wat2 != "" {
		s, err := parseTNXWat(wat2, "latcor")
		if err != nil {
			return nil, fmt.Errorf("wcs/transform: TNX axis 2: %w", err)
		}
		t.Lat = s
	}
	return t, nil
}

// concatWAT concatenates all WATn_mmm strings for a given axis prefix,
// following IRAF's convention that the value may be split across
// WATn_001, WATn_002, ... in 68-character chunks to fit inside a FITS
// header card.
func concatWAT(h *header.Header, prefix string) string {
	var sb strings.Builder
	for i := 1; i < 100; i++ {
		key := prefix + "_" + fmt.Sprintf("%03d", i)
		v, err := h.String(key)
		if err != nil {
			break
		}
		sb.WriteString(v)
	}
	return sb.String()
}

// parseTNXWat extracts the polynomial surface from a concatenated WAT
// string. Looks for the "lngcor" or "latcor" subsection containing the
// coefficients.
func parseTNXWat(wat, tag string) (tnxSurface, error) {
	var s tnxSurface
	idx := strings.Index(wat, tag)
	if idx < 0 {
		return s, nil // no correction on this axis
	}
	// Find the start of the quoted string after "tag = ".
	rest := wat[idx+len(tag):]
	// Skip "= "
	eq := strings.IndexByte(rest, '"')
	if eq < 0 {
		return s, fmt.Errorf("TNX %s: no opening quote", tag)
	}
	rest = rest[eq+1:]
	end := strings.IndexByte(rest, '"')
	if end < 0 {
		return s, fmt.Errorf("TNX %s: no closing quote", tag)
	}
	body := rest[:end]
	fields := strings.Fields(body)
	if len(fields) < 8 {
		return s, fmt.Errorf("TNX %s: need at least 8 tokens, got %d", tag, len(fields))
	}
	var err error
	if s.Type, err = parseIntFloat(fields[0]); err != nil {
		return s, err
	}
	if s.XOrder, err = parseIntFloat(fields[1]); err != nil {
		return s, err
	}
	if s.YOrder, err = parseIntFloat(fields[2]); err != nil {
		return s, err
	}
	if s.Cross, err = parseIntFloat(fields[3]); err != nil {
		return s, err
	}
	if s.XMin, err = strconv.ParseFloat(fields[4], 64); err != nil {
		return s, err
	}
	if s.XMax, err = strconv.ParseFloat(fields[5], 64); err != nil {
		return s, err
	}
	if s.YMin, err = strconv.ParseFloat(fields[6], 64); err != nil {
		return s, err
	}
	if s.YMax, err = strconv.ParseFloat(fields[7], 64); err != nil {
		return s, err
	}
	s.Coeffs = make([]float64, 0, len(fields)-8)
	for _, tok := range fields[8:] {
		v, err := strconv.ParseFloat(tok, 64)
		if err != nil {
			return s, fmt.Errorf("TNX %s: bad coefficient %q", tag, tok)
		}
		s.Coeffs = append(s.Coeffs, v)
	}
	s.present = true
	return s, nil
}

// parseIntFloat accepts "3" or "3." and returns the integer value. IRAF
// sometimes writes integer fields with a trailing dot.
func parseIntFloat(tok string) (int, error) {
	tok = strings.TrimSuffix(tok, ".")
	return strconv.Atoi(tok)
}

// tnxTermAllowed reports whether the (i, j) term is present in the
// coefficient list given the surface's Cross flag. The three options
// mirror the IRAF convention: 0 = no cross terms, 1 = full cross, 2 =
// half cross (total degree < max(XOrder, YOrder)).
func (s tnxSurface) tnxTermAllowed(i, j int) bool {
	switch s.Cross {
	case 0:
		return !(i > 0 && j > 0)
	case 2:
		maxOrder := s.XOrder
		if s.YOrder > maxOrder {
			maxOrder = s.YOrder
		}
		return i+j < maxOrder
	}
	return true // Cross == 1: full
}

// eval computes the polynomial value at (x, y).
func (s tnxSurface) eval(x, y float64) float64 {
	if !s.present {
		return 0
	}
	xn := 2*(x-s.XMin)/(s.XMax-s.XMin) - 1
	yn := 2*(y-s.YMin)/(s.YMax-s.YMin) - 1
	xBasis, _ := tnxBasisWithDeriv(s.Type, xn, s.XOrder)
	yBasis, _ := tnxBasisWithDeriv(s.Type, yn, s.YOrder)

	var sum float64
	k := 0
	for j := 0; j < s.YOrder; j++ {
		for i := 0; i < s.XOrder; i++ {
			if !s.tnxTermAllowed(i, j) {
				continue
			}
			if k >= len(s.Coeffs) {
				return sum
			}
			sum += s.Coeffs[k] * xBasis[i] * yBasis[j]
			k++
		}
	}
	return sum
}

// evalWithDeriv computes the polynomial value and its partial
// derivatives ∂/∂x and ∂/∂y at (x, y). Uses the chain rule through the
// domain normalization: if xn = 2(x-XMin)/(XMax-XMin) - 1, then
// dxn/dx = 2/(XMax-XMin), so ∂P/∂x = ∂P/∂xn · 2/(XMax-XMin). The basis
// derivatives ∂b_i/∂xn are computed alongside b_i by differentiating the
// recurrence relation — see tnxBasisWithDeriv.
func (s tnxSurface) evalWithDeriv(x, y float64) (val, dvdx, dvdy float64) {
	if !s.present {
		return 0, 0, 0
	}
	scaleX := 2 / (s.XMax - s.XMin)
	scaleY := 2 / (s.YMax - s.YMin)
	xn := (x-s.XMin)*scaleX - 1
	yn := (y-s.YMin)*scaleY - 1

	xBasis, xBasisPrime := tnxBasisWithDeriv(s.Type, xn, s.XOrder)
	yBasis, yBasisPrime := tnxBasisWithDeriv(s.Type, yn, s.YOrder)

	k := 0
	for j := 0; j < s.YOrder; j++ {
		for i := 0; i < s.XOrder; i++ {
			if !s.tnxTermAllowed(i, j) {
				continue
			}
			if k >= len(s.Coeffs) {
				return val, dvdx * scaleX, dvdy * scaleY
			}
			c := s.Coeffs[k]
			val += c * xBasis[i] * yBasis[j]
			dvdx += c * xBasisPrime[i] * yBasis[j]
			dvdy += c * xBasis[i] * yBasisPrime[j]
			k++
		}
	}
	// Chain-rule back from normalized to raw coordinates.
	return val, dvdx * scaleX, dvdy * scaleY
}

// tnxBasisWithDeriv returns the basis values [b_0, ..., b_{order-1}]
// and their derivatives [b'_0, ..., b'_{order-1}] for the normalized
// coordinate xn under the given surface function type. Computed in a
// single forward pass by differentiating the basis recurrence.
func tnxBasisWithDeriv(surfaceType int, xn float64, order int) (b, bp []float64) {
	b = make([]float64, order)
	bp = make([]float64, order)
	if order == 0 {
		return
	}
	b[0] = 1
	bp[0] = 0
	if order > 1 {
		b[1] = xn
		bp[1] = 1
	}
	switch surfaceType {
	case 1:
		// Chebyshev T_n(x) = 2x·T_{n-1} - T_{n-2}.
		// Differentiating: T_n'(x) = 2·T_{n-1} + 2x·T_{n-1}' - T_{n-2}'.
		for n := 2; n < order; n++ {
			b[n] = 2*xn*b[n-1] - b[n-2]
			bp[n] = 2*b[n-1] + 2*xn*bp[n-1] - bp[n-2]
		}
	case 2:
		// Legendre: n·P_n = (2n-1)·x·P_{n-1} - (n-1)·P_{n-2}.
		// Differentiating: n·P_n' = (2n-1)·P_{n-1} + (2n-1)·x·P_{n-1}'
		//                          - (n-1)·P_{n-2}'.
		for n := 2; n < order; n++ {
			f := float64(n)
			b[n] = ((2*f-1)*xn*b[n-1] - (f-1)*b[n-2]) / f
			bp[n] = ((2*f-1)*b[n-1] + (2*f-1)*xn*bp[n-1] - (f-1)*bp[n-2]) / f
		}
	case 3:
		// Simple power series: x^n. b_n'(x) = n·x^(n-1).
		for n := 2; n < order; n++ {
			b[n] = xn * b[n-1]
			bp[n] = float64(n) * b[n-1] // n * x^(n-1)
		}
	}
	return
}

// Forward applies the TNX correction at intermediate-world coordinates
// (xi, eta). Output is (xi + lng_correction, eta + lat_correction).
func (t *TNX) Forward(xi, eta float64) (xiOut, etaOut float64) {
	dxi := t.Lon.eval(xi, eta)
	deta := t.Lat.eval(xi, eta)
	return xi + dxi, eta + deta
}

// Inverse iterates Newton's method against Forward, using the analytic
// Jacobian computed via evalWithDeriv. The forward map is
// F(xi, eta) = (xi + Lon(xi, eta), eta + Lat(xi, eta)), so the
// Jacobian is
//
//	J = I + [∂Lon/∂xi  ∂Lon/∂eta]
//	        [∂Lat/∂xi  ∂Lat/∂eta]
func (t *TNX) Inverse(xiPrime, etaPrime float64) (xi, eta float64, ok bool) {
	xi, eta = xiPrime, etaPrime
	for range 50 {
		lonVal, dLonDxi, dLonDeta := t.Lon.evalWithDeriv(xi, eta)
		latVal, dLatDxi, dLatDeta := t.Lat.evalWithDeriv(xi, eta)
		fx := xi + lonVal
		fy := eta + latVal
		rx := xiPrime - fx
		ry := etaPrime - fy
		if absFloat(rx) < 1e-13 && absFloat(ry) < 1e-13 {
			return xi, eta, true
		}
		// J = I + [dLonDxi dLonDeta; dLatDxi dLatDeta]
		j11 := 1 + dLonDxi
		j12 := dLonDeta
		j21 := dLatDxi
		j22 := 1 + dLatDeta
		det := j11*j22 - j12*j21
		if det == 0 {
			return 0, 0, false
		}
		xi += (j22*rx - j12*ry) / det
		eta += (-j21*rx + j11*ry) / det
	}
	return xi, eta, true
}

// absFloat is a tiny helper to avoid importing math for a single call.
func absFloat(x float64) float64 {
	if x < 0 {
		return -x
	}
	return x
}
