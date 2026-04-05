package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// TPV — TAN-PV polynomial distortion (the SCAMP / legacy ESO convention).
//
// TPV encodes a polynomial correction directly via PVi_m projection
// parameters on the intermediate world coordinate axes. Unlike SIP (which
// corrects in pixel space before the CD matrix), TPV corrects in
// intermediate world coordinate space AFTER the linear step and BEFORE
// the TAN projection.
//
// CTYPE1 = 'RA---TPV', CTYPE2 = 'DEC--TPV'
//
// The polynomial uses up to 39 coefficients per axis (degree 7). Paper
// IV defines the coefficient layout:
//
//	PV1_0  = constant
//	PV1_1  = xi
//	PV1_2  = eta
//	PV1_3  = r = sqrt(xi² + eta²)
//	PV1_4  = xi²
//	PV1_5  = xi*eta
//	PV1_6  = eta²
//	PV1_7  = xi³
//	PV1_8  = xi²*eta
//	...and so on up through degree 7
//
// where (xi, eta) are the intermediate world coordinates (in degrees).
// The r terms (PV1_3, PV1_11, PV1_23, PV1_39) are radial and make the
// polynomial non-polynomial in (xi, eta) — they provide a way to encode
// radial distortion efficiently.
//
// Forward is direct polynomial evaluation. Inverse requires Newton
// iteration against the forward polynomial.

// tpvTerms lists the (i, j, r-power) triples for each of the 40 TPV
// coefficient slots. Radial coefficient r is sqrt(xi² + eta²).
//
// This table mirrors the SCAMP / wcslib definition: slots 3, 11, 23, 39
// are radial-only (i.e. r, r³, r⁵, r⁷), the rest are polynomial in xi and
// eta alone.
var tpvTerms = [40]struct {
	xiPow  int
	etaPow int
	rPow   int
}{
	{0, 0, 0},  // PV1_0  = 1
	{1, 0, 0},  // PV1_1  = xi
	{0, 1, 0},  // PV1_2  = eta
	{0, 0, 1},  // PV1_3  = r
	{2, 0, 0},  // PV1_4  = xi²
	{1, 1, 0},  // PV1_5  = xi*eta
	{0, 2, 0},  // PV1_6  = eta²
	{3, 0, 0},  // PV1_7  = xi³
	{2, 1, 0},  // PV1_8  = xi²*eta
	{1, 2, 0},  // PV1_9  = xi*eta²
	{0, 3, 0},  // PV1_10 = eta³
	{0, 0, 3},  // PV1_11 = r³
	{4, 0, 0},  // PV1_12 = xi⁴
	{3, 1, 0},  // PV1_13 = xi³*eta
	{2, 2, 0},  // PV1_14 = xi²*eta²
	{1, 3, 0},  // PV1_15 = xi*eta³
	{0, 4, 0},  // PV1_16 = eta⁴
	{5, 0, 0},  // PV1_17 = xi⁵
	{4, 1, 0},  // PV1_18 = xi⁴*eta
	{3, 2, 0},  // PV1_19 = xi³*eta²
	{2, 3, 0},  // PV1_20 = xi²*eta³
	{1, 4, 0},  // PV1_21 = xi*eta⁴
	{0, 5, 0},  // PV1_22 = eta⁵
	{0, 0, 5},  // PV1_23 = r⁵
	{6, 0, 0},  // PV1_24 = xi⁶
	{5, 1, 0},  // PV1_25 = xi⁵*eta
	{4, 2, 0},  // PV1_26 = xi⁴*eta²
	{3, 3, 0},  // PV1_27 = xi³*eta³
	{2, 4, 0},  // PV1_28 = xi²*eta⁴
	{1, 5, 0},  // PV1_29 = xi*eta⁵
	{0, 6, 0},  // PV1_30 = eta⁶
	{7, 0, 0},  // PV1_31 = xi⁷
	{6, 1, 0},  // PV1_32 = xi⁶*eta
	{5, 2, 0},  // PV1_33 = xi⁵*eta²
	{4, 3, 0},  // PV1_34 = xi⁴*eta³
	{3, 4, 0},  // PV1_35 = xi³*eta⁴
	{2, 5, 0},  // PV1_36 = xi²*eta⁵
	{1, 6, 0},  // PV1_37 = xi*eta⁶
	{0, 7, 0},  // PV1_38 = eta⁷
	{0, 0, 7},  // PV1_39 = r⁷
}

// TPV holds a parsed TPV distortion: one coefficient array per axis.
type TPV struct {
	Axis1 [40]float64
	Axis2 [40]float64
	Lon   int // 1-based longitude axis index (conventionally 1)
	Lat   int // 1-based latitude axis index (conventionally 2)
}

// ParseTPV extracts TPV coefficients from a parsed WCS header. Returns
// nil if no TPV coefficients are present (i.e. no PV1_* or PV2_* values).
func ParseTPV(w *wcs.Header) *TPV {
	if w.LonAxis == 0 || w.LatAxis == 0 {
		return nil
	}
	t := &TPV{Lon: w.LonAxis, Lat: w.LatAxis}
	hasAny := false
	for i := 0; i < 40; i++ {
		if v, ok := w.PV[wcs.PVKey{Axis: w.LonAxis, Index: i}]; ok {
			t.Axis1[i] = v
			hasAny = true
		}
		if v, ok := w.PV[wcs.PVKey{Axis: w.LatAxis, Index: i}]; ok {
			t.Axis2[i] = v
			hasAny = true
		}
	}
	if !hasAny {
		return nil
	}
	// Set unit-diagonal defaults: if PV_1 is not explicitly set, assume 1
	// (identity). This matches the SCAMP convention.
	if _, ok := w.PV[wcs.PVKey{Axis: w.LonAxis, Index: 1}]; !ok {
		t.Axis1[1] = 1.0
	}
	if _, ok := w.PV[wcs.PVKey{Axis: w.LatAxis, Index: 1}]; !ok {
		t.Axis2[1] = 1.0
	}
	return t
}

// evalOne evaluates one axis's polynomial at intermediate coordinates
// (xi, eta), in the same units as the PV coefficients (typically degrees).
func evalTPV(coefs *[40]float64, xi, eta float64) float64 {
	r := math.Hypot(xi, eta)
	var sum float64
	// Pre-compute power tables up to 7.
	xiPow := [8]float64{1, xi, xi * xi, xi * xi * xi, 0, 0, 0, 0}
	etaPow := [8]float64{1, eta, eta * eta, eta * eta * eta, 0, 0, 0, 0}
	rPow := [8]float64{1, r, 0, 0, 0, 0, 0, 0}
	for k := 4; k < 8; k++ {
		xiPow[k] = xiPow[k-1] * xi
		etaPow[k] = etaPow[k-1] * eta
		rPow[k] = rPow[k-1] * r
	}
	for k, c := range coefs {
		if c == 0 {
			continue
		}
		term := tpvTerms[k]
		v := c * xiPow[term.xiPow] * etaPow[term.etaPow]
		if term.rPow > 0 {
			v *= rPow[term.rPow]
		}
		sum += v
	}
	return sum
}

// Forward applies the TPV polynomial to intermediate world coordinates.
// Input (xi, eta) is the pre-distortion position; output is the
// distortion-corrected position, both in degrees.
func (t *TPV) Forward(xi, eta float64) (xiOut, etaOut float64) {
	return evalTPV(&t.Axis1, xi, eta), evalTPV(&t.Axis2, xi, eta)
}

// Inverse applies the inverse TPV correction via Newton iteration
// against the forward polynomial.
func (t *TPV) Inverse(xiPrime, etaPrime float64) (xi, eta float64, ok bool) {
	// Initial guess: pass-through (TPV is close to identity for small
	// distortions).
	xi, eta = xiPrime, etaPrime
	for range 50 {
		fx, fy := t.Forward(xi, eta)
		rx := xiPrime - fx
		ry := etaPrime - fy
		if math.Abs(rx) < 1e-13 && math.Abs(ry) < 1e-13 {
			return xi, eta, true
		}
		const h = 1e-7
		fx1, fy1 := t.Forward(xi+h, eta)
		fx2, fy2 := t.Forward(xi, eta+h)
		j11 := (fx1 - fx) / h
		j12 := (fx2 - fx) / h
		j21 := (fy1 - fy) / h
		j22 := (fy2 - fy) / h
		det := j11*j22 - j12*j21
		if det == 0 {
			return 0, 0, false
		}
		xi += (j22*rx - j12*ry) / det
		eta += (-j21*rx + j11*ry) / det
	}
	return xi, eta, true
}
