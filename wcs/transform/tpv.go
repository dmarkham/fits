package transform

import (
	"math"

	"github.com/dmarkham/fits/wcs"
)

// TPV — TAN-PV polynomial distortion (the SCAMP / legacy ESO convention),
// defined by the FITS registry note
// https://fits.gsfc.nasa.gov/registry/tpvwcs/tpv.html ("the registry").
//
// CTYPE1 = 'RA---TPV', CTYPE2 = 'DEC--TPV'
//
// The registry defines TPV on CD-matrix intermediate coordinates
// (xi, eta) = CD (p - CRPIX), in degrees: (xi', eta') = f(xi, eta), then
// TAN. For a PC + CDELT header (which the registry does not define) we
// follow wcslib, which applies TPV as a sequent distortion BEFORE CDELT,
// i.e. on PC (p - CRPIX) (wcslib dis.h, "Historical idiosyncrasies": the
// coefficients must then be scaled for the change of units). With
// CDELT = 1, which is what astropy/wcslib write, both readings agree.
//
// The polynomial uses up to 40 coefficients per axis (PVi_0 .. PVi_39,
// degree 7), in the registry's layout:
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
// The r terms (PV1_3, PV1_11, PV1_23, PV1_39) are radial (odd powers of r
// only).
//
// The axis-2 polynomial uses the same layout with xi and eta exchanged:
// PV2_1 multiplies eta, PV2_2 multiplies xi, PV2_4 is eta², and so on
// (registry, wcslib tpv.c, SCAMP).
//
// Defaults: the registry says missing PV keywords default to 0 except
// PV1_1 and PV2_1, which default to 1, and we follow it. wcslib 8.10
// differs: when any PVi_m is present it defaults a missing PVi_1 to 0
// (verified with astropy 8.0.2 / wcslib 8.10). Real SCAMP/DECam/ZTF/PTF
// headers always write PVi_1, so this only matters for hand-made headers;
// TestTPVMissingPV1DefaultsToOne pins our choice.
//
// Forward is direct polynomial evaluation. Inverse is Newton iteration
// against the forward polynomial (see (*TPV).Inverse).

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
	{0, 0, 0}, // PV1_0  = 1
	{1, 0, 0}, // PV1_1  = xi
	{0, 1, 0}, // PV1_2  = eta
	{0, 0, 1}, // PV1_3  = r
	{2, 0, 0}, // PV1_4  = xi²
	{1, 1, 0}, // PV1_5  = xi*eta
	{0, 2, 0}, // PV1_6  = eta²
	{3, 0, 0}, // PV1_7  = xi³
	{2, 1, 0}, // PV1_8  = xi²*eta
	{1, 2, 0}, // PV1_9  = xi*eta²
	{0, 3, 0}, // PV1_10 = eta³
	{0, 0, 3}, // PV1_11 = r³
	{4, 0, 0}, // PV1_12 = xi⁴
	{3, 1, 0}, // PV1_13 = xi³*eta
	{2, 2, 0}, // PV1_14 = xi²*eta²
	{1, 3, 0}, // PV1_15 = xi*eta³
	{0, 4, 0}, // PV1_16 = eta⁴
	{5, 0, 0}, // PV1_17 = xi⁵
	{4, 1, 0}, // PV1_18 = xi⁴*eta
	{3, 2, 0}, // PV1_19 = xi³*eta²
	{2, 3, 0}, // PV1_20 = xi²*eta³
	{1, 4, 0}, // PV1_21 = xi*eta⁴
	{0, 5, 0}, // PV1_22 = eta⁵
	{0, 0, 5}, // PV1_23 = r⁵
	{6, 0, 0}, // PV1_24 = xi⁶
	{5, 1, 0}, // PV1_25 = xi⁵*eta
	{4, 2, 0}, // PV1_26 = xi⁴*eta²
	{3, 3, 0}, // PV1_27 = xi³*eta³
	{2, 4, 0}, // PV1_28 = xi²*eta⁴
	{1, 5, 0}, // PV1_29 = xi*eta⁵
	{0, 6, 0}, // PV1_30 = eta⁶
	{7, 0, 0}, // PV1_31 = xi⁷
	{6, 1, 0}, // PV1_32 = xi⁶*eta
	{5, 2, 0}, // PV1_33 = xi⁵*eta²
	{4, 3, 0}, // PV1_34 = xi⁴*eta³
	{3, 4, 0}, // PV1_35 = xi³*eta⁴
	{2, 5, 0}, // PV1_36 = xi²*eta⁵
	{1, 6, 0}, // PV1_37 = xi*eta⁶
	{0, 7, 0}, // PV1_38 = eta⁷
	{0, 0, 7}, // PV1_39 = r⁷
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
	// Registry defaults: PV1_1 and PV2_1 are 1 when absent (wcslib 8.10
	// uses 0 instead; see the file comment).
	if _, ok := w.PV[wcs.PVKey{Axis: w.LonAxis, Index: 1}]; !ok {
		t.Axis1[1] = 1.0
	}
	if _, ok := w.PV[wcs.PVKey{Axis: w.LatAxis, Index: 1}]; !ok {
		t.Axis2[1] = 1.0
	}
	return t
}

// tpvPowers computes the power tables for xi, eta, and r = sqrt(xi²+eta²)
// at indices 0..7. Every index is filled correctly; unlike a previous
// version that had a bug leaving rPow[2..7] = 0, this is usable for the
// radial terms (slots 3, 11, 23, 39 → powers 1, 3, 5, 7) and for
// analytic derivatives that need r^(p-2).
func tpvPowers(xi, eta float64) (xiPow, etaPow, rPow [8]float64) {
	r := math.Hypot(xi, eta)
	xiPow[0], etaPow[0], rPow[0] = 1, 1, 1
	for k := 1; k < 8; k++ {
		xiPow[k] = float64(xiPow[k-1] * xi)
		etaPow[k] = float64(etaPow[k-1] * eta)
		rPow[k] = float64(rPow[k-1] * r)
	}
	return
}

// evalTPV evaluates one axis's polynomial at (xi, eta), in the same
// units as the PV coefficients (typically degrees).
func evalTPV(coefs *[40]float64, xi, eta float64) float64 {
	xiPow, etaPow, rPow := tpvPowers(xi, eta)
	var sum float64
	for k, c := range coefs {
		if c == 0 {
			continue
		}
		term := tpvTerms[k]
		v := float64(c * xiPow[term.xiPow] * etaPow[term.etaPow])
		if term.rPow > 0 {
			v *= rPow[term.rPow]
		}
		sum += v
	}
	return sum
}

// evalTPVWithJacobian evaluates one axis's polynomial at (xi, eta) and
// its analytic partial derivatives ∂/∂xi and ∂/∂eta in a single pass.
//
// For a non-radial term c·xi^a·eta^b:
//
//	∂/∂xi  = c·a·xi^(a-1)·eta^b     (0 if a=0)
//	∂/∂eta = c·b·xi^a·eta^(b-1)     (0 if b=0)
//
// For a radial term c·r^p (TPV has only p=1, 3, 5, 7 with xi/eta
// exponents zero):
//
//	∂/∂xi  = c·p·r^(p-2)·xi
//	∂/∂eta = c·p·r^(p-2)·eta
//
// For p=1 the derivative is c·xi/r, singular at the origin. We skip
// the radial contribution when r is below a tiny threshold — the
// linear PV_1 / PV_2 coefficient dominates the Jacobian near the
// origin so this doesn't affect convergence.
func evalTPVWithJacobian(coefs *[40]float64, xi, eta float64) (val, dvdxi, dvdeta float64) {
	xiPow, etaPow, rPow := tpvPowers(xi, eta)
	r := rPow[1]
	const rMin = 1e-20
	for k, c := range coefs {
		if c == 0 {
			continue
		}
		term := tpvTerms[k]
		a, b, p := term.xiPow, term.etaPow, term.rPow

		if p == 0 {
			// Non-radial monomial.
			val += float64(c * xiPow[a] * etaPow[b])
			if a >= 1 {
				dvdxi += float64(c * float64(a) * xiPow[a-1] * etaPow[b])
			}
			if b >= 1 {
				dvdeta += float64(c * xiPow[a] * float64(b) * etaPow[b-1])
			}
		} else {
			// Radial term: c·r^p (TPV has a=b=0 for all radial slots).
			val += float64(c * rPow[p])
			if r > rMin {
				if p == 1 {
					dvdxi += c * xi / r
					dvdeta += c * eta / r
				} else {
					// p ≥ 3 odd: c·p·r^(p-2)·(xi|eta)
					rp2 := rPow[p-2]
					dvdxi += float64(c * float64(p) * rp2 * xi)
					dvdeta += float64(c * float64(p) * rp2 * eta)
				}
			}
		}
	}
	return
}

// Forward applies the TPV polynomial to intermediate world coordinates.
// Input (xi, eta) is the pre-distortion position; output is the
// distortion-corrected position, both in degrees.
func (t *TPV) Forward(xi, eta float64) (xiOut, etaOut float64) {
	// Axis 2 evaluates the shared term table with xi and eta exchanged.
	return evalTPV(&t.Axis1, xi, eta), evalTPV(&t.Axis2, eta, xi)
}

// Inverse inverts the TPV polynomial by Newton iteration from (xi',
// eta') with the analytic Jacobian (one polynomial walk per axis per
// iteration; finite differences would lose precision at small steps). See
// newton2 for the relative convergence test. ok is false when the
// iteration does not converge (no inverse near the start point, e.g.
// beyond a fold of the polynomial); the returned values are then
// meaningless.
func (t *TPV) Inverse(xiPrime, etaPrime float64) (xi, eta float64, ok bool) {
	return newton2(func(xi, eta float64) (fx, fy, j11, j12, j21, j22 float64) {
		fx, dfxdxi, dfxdeta := evalTPVWithJacobian(&t.Axis1, xi, eta)
		// Axis 2 swaps xi and eta, so its partials come back swapped too.
		fy, dfydeta, dfydxi := evalTPVWithJacobian(&t.Axis2, eta, xi)
		return fx, fy, dfxdxi, dfxdeta, dfydxi, dfydeta
	}, xiPrime, etaPrime, xiPrime, etaPrime, tpvNewtonFloor)
}

// tpvNewtonFloor is the coordinate size below which TPV convergence is
// tested absolutely (degrees in the usual CD form: 1e-10 deg = 0.36 uas).
const tpvNewtonFloor = 1e-10
