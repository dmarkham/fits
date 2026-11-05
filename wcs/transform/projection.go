// Package transform computes forward and inverse FITS World Coordinate
// System mappings between image pixel coordinates and celestial spherical
// coordinates.
//
// The coordinate pipeline follows Greisen & Calabretta 2002 Paper I §2.1:
//
//  1. Pixel offset from reference point:
//     (p_i - CRPIX_i)
//
//  2. Linear transform to intermediate world coordinates (applied as
//     either a CD matrix or a PC matrix with CDELT scaling):
//     q_i = sum_j CD_{ij} * (p_j - CRPIX_j)
//     or
//     q_i = CDELT_i * sum_j PC_{ij} * (p_j - CRPIX_j)
//
//  3. Projection from intermediate world coordinates to native spherical
//     coordinates (phi, theta) — this is the projection-specific step.
//
//  4. Spherical rotation from native to celestial coordinates (alpha,
//     delta) using LONPOLE / LATPOLE and the reference (CRVAL_lon,
//     CRVAL_lat).
//
// # Supported projections
//
// All 28 projection codes that wcslib 8.10 implements are supported: the
// 26 of Calabretta & Greisen 2002 (Paper II), HPX (Calabretta & Roukema
// 2007) and XPH (Calabretta & Lowe 2013).
//
//	Zenithal (§5.1):           AZP SZP TAN STG SIN ARC ZPN ZEA AIR
//	Cylindrical (§5.2):        CYP CEA CAR MER
//	Pseudo-cylindrical (§5.3): SFL PAR MOL AIT
//	Conic (§5.4):              COP COE COD COO
//	Polyconic (§5.5):          BON PCO
//	Quad-cube (§5.6):          TSC CSC QSC
//	HEALPix:                   HPX XPH
//
// New also accepts the AIPS aliases NCP and GLS in CTYPE and translates
// them as wcslib's wcsset does: NCP is SIN with PV2_1 = 0 and PV2_2 =
// cot(CRVAL of the latitude axis) (refused on the equator), GLS is SFL
// with the fiducial offset phi0 = 0, theta0 = CRVAL of the latitude axis
// (hdrcases/alias_*.case). Select itself knows only the 28 codes.
//
// The projection code is a port of wcslib's prj.c (wcslib 8.10), in
// wcslib's own degree units and with its exact trigonometric special cases
// (wcstrig.go), so that both the values AND the decision of which points
// lie on the projection match wcslib. Every code, with default and
// non-default PVi_m parameters, is checked against wcslib point by point in
// both directions (wref_test.go, goldens from testdata/wref/wref.c): on
// grids that run past every boundary, and at probe points 1e-7 deg either
// side of every accept/reject crossing on 80 inverse and 48 forward scan
// lines, at and around every pole image, and next to both poles, so each
// boundary is pinned to 1e-7 deg, not to the grid spacing. Also through the full
// header pipeline (hdr_golden_test.go, goldens from testdata/wref/hdr2sky.c).
// The sin/cos/tan/atan/atan2 underneath are Go's math functions (asin and
// acos are cancellation-free identities on top of them, wcstrig.go), which
// can differ from glibc's, which wcslib calls, in the last bit. Values
// agree with wcslib within 1e-11 deg (forward rows beyond 1000 deg: 1e-14
// relative); where the projection itself magnifies a last-bit libm
// difference (next to a pole, a horizon, a divergence or a cube-face
// centre) a row may exceed that only by its own libm sensitivity, the
// change of the output when every libm result moves by 4 ULP, and such
// rows are counted (0.2% of the wref rows). Accept/reject decisions agree
// with wcslib except at boundary-ambiguous points, where the port's own
// decision flips under a nudge of a few ULP (2 of 171178 wref rows, both
// on AZP's horizon). See libm_sens_test.go.
// XPH and HPX are also checked
// against the equations of their defining papers, independently of wcslib
// (xph_paper_test.go, hpx_paper_test.go). Next to the poles both inherit
// wcslib's double-precision limit (measured against the exact equations:
// HPX up to 3.0e-7 deg forward and 4.7e-7 deg inverse within ~1e-6 deg of
// a pole, XPH 4.1e-7 deg forward next to the south pole), which the port
// keeps so that it stays identical to wcslib.
// A point that is not on the projection is rejected (ok=false / an error),
// never returned as a coordinate; where wcslib itself returns a NaN or
// infinite coordinate with a success status (CYP with an unreachable
// latitude, PCO far off the map, COO at the apex of a southern cone), we
// reject the point.
//
// The projection functions here work at wcslib's default fiducial point
// (phi_0, theta_0). A header that moves it with PVi_0 = 1 and PVi_1/PVi_2
// on the LONGITUDE axis is handled by transform.New, which applies the
// fiducial offsets x0/y0 as wcslib's prjoff does; that is proven against
// wcslib 8.10 in both directions by the pvfid_* header cases
// (wcs/transform/testdata/wref/hdrcases/pvfid_*.case, TestHdrGoldens),
// including PV1_0 = 1 with only one of PV1_1/PV1_2 given, where wcslib
// applies no offset (pv1off_car_pv1_1_only, pv1off_car_pv1_2_only).
package transform

import (
	"errors"
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// Guards that no input reaches.
//
// The port keeps every check of wcslib 8.10 prj.c, so that its accept/
// reject rule is wcslib's. A few of them cannot be reached by any input,
// in wcslib or here; changing their tolerance changes no result, so no
// test can pin it through the projections (a mutant that widens it is an
// equivalent mutant). Each is listed with the reason, and where a direct
// pin is possible the test that does it. Proven:
//
//   - asind/acosd's WCSTRIG_TOL (wcstrig.go): callers' arguments are within
//     a few ULP of [-1, 1] or far outside it (the verifier toggled the
//     tolerance over all wref cases, inverse [-400, 400]^2 step 0.5 and the
//     forward grids: no output changed). Pinned branch by branch against
//     wcslib's own asind/acosd: TestWcstrigGuards.
//   - TSC/QSC forward face overshoot (quadcubeFaceBound): |xf| = |m| / zeta
//     <= 1 exactly for TSC; QSC exceeds 1 only by rounding. A wcslib scan
//     of 4.62e8 points per code found no rejection. Pinned directly:
//     TestQuadcubeFaceBound.
//   - prjbchk (the final [-180, 180] x [-90, 90] check, tolerance 1e-13 or
//     1e-12) where its inputs cannot leave the range by more than a few
//     ULP, or not at all. By construction: TAN, SIN (both branches), ZEA
//     and SZP (phi = atan2d, within [-180, 180]; theta from atan2d(r0, r),
//     asind, acosd or 90 - 2 asind(s) with 0 <= s <= 1), and AIT (s is
//     clamped to >= 0.5 just before, so x0 = 2 z^2 - 1 >= 0 up to rounding
//     and phi = 2 atan2d(y0, x0) passes 180 by a few 1e-14 deg at most).
//     TSC, CSC, QSC and AIR never even reach the clamp: their phi is
//     atan2d (|atan2d| <= fl(fl(pi) r2d) = 180 exactly), their theta asind
//     (<= fl(fl(pi/2) r2d) = 90) or, for AIR, 90 - 2 acosd(c) with c in
//     (0, 1] (>= -90 exactly) or 90 - 2 r w6 next to the pole
//     (TestPrjbchkRangeByConstruction). AZP's and XPH's prjbchk, by
//     contrast, DO decide points (AZP with mu = 1 and a tilt lands a few
//     1e-13 deg below -90 at its horizon; XPH lands 1e-14 .. 1e-7 deg past
//     +-180 next to its north pole on the x = 0 seam); both are pinned by
//     wcslib probes on either side of the tolerance (wref cases
//     azp_mu1_gamma45 and xph_default, wref.c xph_seam_inv).
//   - SIN's and SZP's root tolerance on sin(theta) (a root within 1e-13
//     of +-1 is clamped to +-1). Both solve
//     f(s) = a s^2 + 2 b s + c = 0 with a > 0 and, in exact arithmetic,
//     f(1) = r^2 (r^2 = x^2 + y^2 in radians) and f(-1) = |P - 2 Q|^2 >= 0
//     (P = (x, y), Q = the slant point; for SIN Q = (xi, eta)). A root
//     1 + e (e > 0) needs f'(1) = 2 (1 + P.Q) < 0, so -P.Q > 1, and a real
//     root then needs (1 + P.Q)^2 >= a r^2 > r^2 |Q|^2 >= (P.Q)^2
//     (a = 1 + |Q|^2, Cauchy-Schwarz), i.e. (s - 1)^2 > s^2 with
//     s = -P.Q > 1: impossible. Likewise a root -1 - e needs
//     m = P.Q - 2|Q|^2 > 1 and (m - 1)^2 >= a |P - 2Q|^2 > m^2 / |Q|^2 *
//     (1 + |Q|^2) > m^2: impossible. So no root lies outside [-1, 1]
//     except by rounding. That is a few ULP: next to +1 (the pole, where
//     r^2 >= 1e-10 keeps the root ~5e-11 below 1) the root is (-b + d)/a
//     with -b ~ |Q|^2 >= 0 and d >= 0, next to -1 (P ~ 2 Q) it is
//     (-b - d)/a with -b ~ -|Q|^2 <= 0, so neither cancels, and a 1e-13
//     and a 1e-6 tolerance clamp the same points (the 1e-6 mutants
//     survive for this reason).
//   - The comparison (< against <=) of four small-angle switches: SIN's
//     forward t = (90 - |theta|) d2r < 1e-5, AIR's forward
//     xi = d2r (90 - theta) / 2 < 1e-4, XPH's inverse
//     sigma = (90 - |eta|)/45 < 1e-4 and QSC's forward
//     zeco = 1 - zeta < 1e-8. No double input makes the computed value
//     equal the threshold, so < and <= decide alike
//     (TestSmallAngleSwitchTies, TestQSCZecoSwitchTie). The inverse
//     SIN/SZP and AIR switches can be tied and are probed at the tie
//     (wref.c tie_inv).
//   - QSC's inverse zeta tolerance: zeta = 1 - zeco is always in [0.5, 1],
//     never near -1 (qscInverseZeta; TestQSCInverseZetaRange).
//   - CSC's float32-underflow guards (chi^2, psi^2, |chi psi| > 1e-16):
//     the term they drop is below half an ULP of the partial sum it is
//     added to, or of the term added next (cscXfYf;
//     TestCSCUnderflowGuardPremises), for any threshold up to 1e-8.
//
// Measured, not proven:
//
//   - CSC's forward overshoot tolerance (1e-7 in face units): xf, yf reach
//     +-1 only at the face edges, and wcslib's own cscs2x rejects none of
//     6.5e8 points of a 0.01 deg grid; the float32 polynomial evaluated
//     at every float32 chi in [0.99, 1] times 4001 psi never exceeds 1.
//     No proof that rounding can never push it past 1 + 1e-7.
//   - AIR's bracketing failure (airProjection.x2s, k == 30): when the 30
//     halvings fail, the regula falsi that would follow starts with both
//     bracket ends below r, and its iterates rt = R(cos xi) are computed
//     as fl(-w1 fl(1 / cos xi)) (the log term is below half an ULP there)
//     for |w1| >= 0.4 (theta_b >= ~20 deg), monotone in cos xi, so
//     rt <= R(2^-30) < r, r - rt >= ULP(r) > 1e-12, and its own k == 100
//     check rejects the point too. For smaller |w1| this is measured only:
//     an exhaustive scan of 2.2e6 doubles r above R(2^-30) and 2e5
//     log-spaced r up to 1e7 R(2^-30) found no difference for either
//     committed theta_b. The decisions at r = 1e6..1e15 deg are in the
//     air_* goldens (wref.c extra_inv).
//
// Everything else that decides a result - every boundary, small-angle
// switch and parameter check - is pinned against wcslib:
// TestWrefProjections (boundaries to 1e-7 deg, every point to 1e-11 deg),
// TestPCOSmallAngleBand (PCO's 1e-4 deg switches, to 0.1% either side and
// at the threshold itself, 4 ULP), TestPoleSmallAngleBand (the switches
// next to the native pole: SIN forward and inverse, SZP, AIR and XPH, to
// 0.1%, SIN's forward to 0.05%, 4 ULP), TestQuadcubeFaceCentreBand (QSC's
// zeco < 1e-8 switch on all six faces, to 0.05%, 4 ULP; measured
// bit-identical), TestSelectRejectsWhatPrjsetRejects and the paper tests.
// The mutation sweep that checked this (every threshold moved by 0.1%,
// 2%, 20%, 2x and 10x either way, < against <=, every prjbchk tolerance
// to 1e-6) is recorded in docs/VALIDATION.md.

// ErrUnsupportedProjection is returned by Select for a code listed in
// allProjectionCodes with implemented=false. No such code exists today
// (all 28 codes are implemented); the sentinel is kept so that callers
// can distinguish "recognized but unsupported" from "unknown".
var ErrUnsupportedProjection = errors.New("wcs/transform: unsupported projection")

// ErrUnknownProjection is returned by Select when the 3-letter code is not
// recognized at all (i.e. it is not in the Paper II projection list).
var ErrUnknownProjection = errors.New("wcs/transform: unknown projection code")

// Projection is the interface implemented by every spherical projection.
//
// All angles are in RADIANS. Implementations work with FITS "native"
// longitude phi and "native" latitude theta on the unit sphere.
//
// Forward: given (phi, theta), return intermediate world coordinates
// (x, y) — also in radians. ok=false indicates the input falls outside
// the projection's valid domain (e.g. the far hemisphere for SIN).
//
// Inverse: given intermediate (x, y), return the native (phi, theta).
// ok=false indicates the point is not on the projection's image.
//
// Theta0 returns the native latitude of the projection's reference
// (fiducial) point, in radians, as wcslib's prjset assigns it: pi/2 for
// the zenithal projections and XPH, theta_a for the conics, 0 for the
// rest (BON included). transform.New uses it to place the native pole.
type Projection interface {
	Code() string
	Theta0() float64
	Forward(phi, theta float64) (x, y float64, ok bool)
	Inverse(x, y float64) (phi, theta float64, ok bool)
}

// allProjectionCodes enumerates every projection code defined by Paper II,
// HPX (Calabretta & Roukema 2007) and XPH (Calabretta & Lowe 2013), with a
// bit saying whether it is implemented. TestWrefProjections requires a
// wcslib golden for every code listed here.
var allProjectionCodes = map[string]bool{
	// Zenithal (§5.1)
	"AZP": true, "SZP": true, "TAN": true, "STG": true, "SIN": true,
	"ARC": true, "ZPN": true, "ZEA": true, "AIR": true,
	// Cylindrical (§5.2)
	"CYP": true, "CEA": true, "CAR": true, "MER": true,
	// Pseudo-cylindrical and conventional (§5.3)
	"SFL": true, "PAR": true, "MOL": true, "AIT": true,
	// Conic (§5.4)
	"COP": true, "COE": true, "COD": true, "COO": true,
	// Polyconic and pseudoconic (§5.5)
	"BON": true, "PCO": true,
	// Quad-cube (§5.6)
	"TSC": true, "CSC": true, "QSC": true,
	// HEALPix (Calabretta & Roukema 2007; Calabretta & Lowe 2013)
	"HPX": true, "XPH": true,
}

// Select returns the Projection for a 3-letter FITS projection code,
// configured from the projection-parameter map pv (PVi_m values keyed by
// (axis, index)); the parameters are read from latAxis, the 1-based
// latitude axis, as the standard requires. Absent parameters take the
// same defaults as wcslib's prjset.
//
// Returns ErrUnknownProjection for codes not in the FITS standard list,
// and an error for a parameter set wcslib's prjset also rejects (e.g. COP
// without PV2_1, CEA with lambda outside (0, 1], ZPN with every
// coefficient zero).
func Select(code string, pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	var p degProjection
	var err error
	switch code {
	// Zenithal (§5.1)
	case "AZP":
		p, err = newAZP(pv, latAxis)
	case "SZP":
		p, err = newSZP(pv, latAxis)
	case "TAN":
		p = tanProjection{}
	case "STG":
		p = stgProjection{}
	case "SIN":
		p = newSIN(pv, latAxis)
	case "ARC":
		p = arcProjection{}
	case "ZPN":
		p, err = newZPN(pv, latAxis)
	case "ZEA":
		p = zeaProjection{}
	case "AIR":
		p, err = newAIR(pv, latAxis)
	// Cylindrical (§5.2)
	case "CYP":
		p, err = newCYP(pv, latAxis)
	case "CEA":
		p, err = newCEA(pv, latAxis)
	case "CAR":
		p = carProjection{}
	case "MER":
		p = merProjection{}
	// Pseudo-cylindrical (§5.3)
	case "SFL":
		p = sflProjection{}
	case "PAR":
		p = parProjection{}
	case "MOL":
		p = molProjection{}
	case "AIT":
		p = aitProjection{}
	// Conic (§5.4)
	case "COP":
		p, err = newCOP(pv, latAxis)
	case "COE":
		p, err = newCOE(pv, latAxis)
	case "COD":
		p, err = newCOD(pv, latAxis)
	case "COO":
		p, err = newCOO(pv, latAxis)
	// Polyconic (§5.5)
	case "BON":
		p, err = newBON(pv, latAxis)
	case "PCO":
		p = pcoProjection{}
	// Quad-cube (§5.6)
	case "TSC":
		p = tscProjection{}
	case "CSC":
		p = cscProjection{}
	case "QSC":
		p = qscProjection{}
	// HEALPix
	case "HPX":
		p, err = newHPX(pv, latAxis)
	case "XPH":
		p = xphProjection{}
	default:
		if allProjectionCodes[code] {
			return nil, fmt.Errorf("%w: %q (recognized but not supported)", ErrUnsupportedProjection, code)
		}
		return nil, fmt.Errorf("%w: %q", ErrUnknownProjection, code)
	}
	if err != nil {
		return nil, err
	}
	return degAdapter{p}, nil
}

// degProjection is a projection ported from wcslib prj.c and evaluated in
// wcslib's native units: every angle and every plane coordinate in
// DEGREES (wcslib's r0 = 180/pi). ok=false marks a point wcslib's
// prjs2x/prjx2s gives a non-zero stat for.
type degProjection interface {
	Code() string
	// theta0 is the native latitude of the fiducial point, degrees (the
	// theta0 wcslib's prjset assigns).
	theta0() float64
	s2x(phi, theta float64) (x, y float64, ok bool)
	x2s(x, y float64) (phi, theta float64, ok bool)
}

// errBadParam is the error for a parameter set wcslib's prjset rejects.
func errBadParam(code, format string, args ...any) error {
	return fmt.Errorf("wcs/transform: %s: invalid projection parameters: %s", code, fmt.Sprintf(format, args...))
}

// degAdapter presents a degProjection through the radian Projection
// interface. It also refuses a non-finite result: wcslib returns NaN or
// infinity with a zero stat in a few places (see the package comment) and
// such a value must never escape as a coordinate.
type degAdapter struct{ p degProjection }

func (a degAdapter) Code() string { return a.p.Code() }

func (a degAdapter) Theta0() float64 {
	t := a.p.theta0()
	if t == 90 {
		return math.Pi / 2 // exact, transform.New tests for it
	}
	return degToRad(t)
}

// forwardDeg / inverseDeg are the degree-unit entry points used by the
// wcslib comparison tests; Forward / Inverse wrap them.
func (a degAdapter) forwardDeg(phi, theta float64) (x, y float64, ok bool) {
	x, y, ok = a.p.s2x(phi, theta)
	if !ok || !isFinite(x) || !isFinite(y) {
		return 0, 0, false
	}
	return x, y, true
}

func (a degAdapter) inverseDeg(x, y float64) (phi, theta float64, ok bool) {
	phi, theta, ok = a.p.x2s(x, y)
	if !ok || !isFinite(phi) || !isFinite(theta) {
		return 0, 0, false
	}
	return phi, theta, true
}

func (a degAdapter) Forward(phi, theta float64) (x, y float64, ok bool) {
	x, y, ok = a.forwardDeg(radToDeg(phi), radToDeg(theta))
	if !ok {
		return 0, 0, false
	}
	return degToRad(x), degToRad(y), true
}

func (a degAdapter) Inverse(x, y float64) (phi, theta float64, ok bool) {
	phi, theta, ok = a.inverseDeg(radToDeg(x), radToDeg(y))
	if !ok {
		return 0, 0, false
	}
	return degToRad(phi), degToRad(theta), true
}

func isFinite(v float64) bool { return !math.IsNaN(v) && !math.IsInf(v, 0) }

// pvSet reports whether PVaxis_index is present.
func pvSet(pv map[wcs.PVKey]float64, axis, index int) bool {
	if pv == nil {
		return false
	}
	_, ok := pv[wcs.PVKey{Axis: axis, Index: index}]
	return ok
}

// pvFloat returns the value of PVaxis_index from pv, or def if absent.
func pvFloat(pv map[wcs.PVKey]float64, axis, index int, def float64) float64 {
	if pv == nil {
		return def
	}
	if v, ok := pv[wcs.PVKey{Axis: axis, Index: index}]; ok {
		return v
	}
	return def
}

// degToRad converts degrees to radians.
func degToRad(d float64) float64 { return d * math.Pi / 180 }

// radToDeg converts radians to degrees.
func radToDeg(r float64) float64 { return r * 180 / math.Pi }

// normalizeDeg maps a degree value into [0, 360).
func normalizeDeg(d float64) float64 {
	d = math.Mod(d, 360)
	if d < 0 {
		d += 360
	}
	return d
}
