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
// Every projection from Paper II is either fully implemented or
// recognized and returning ErrUnsupportedProjection — never a silent
// wrong answer.
//
// Zenithal family (§5.1):   AZP SZP TAN STG SIN ARC ZPN ZEA
// Cylindrical (§5.2):       CYP CEA CAR MER
// Pseudo-cylindrical (§5.3): SFL PAR MOL AIT
// Conic (§5.4):             COP COE COD COO
// Polyconic (§5.5):         BON PCO
// Quad-cube (§5.6):         TSC CSC QSC
// HEALPix (CR 2007):        HPX XPH
package transform

import (
	"errors"
	"fmt"
	"math"

	"github.com/dmarkham/fits/wcs"
)

// ErrUnsupportedProjection is returned by Select when asked for a
// projection code the transform package does not yet implement.
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
// (fiducial) point, in radians. Zenithal projections have theta_0 = pi/2;
// most others have theta_0 = 0. Used by the transform constructor to
// apply the correct LONPOLE / LATPOLE defaults per Paper II §2.4.
type Projection interface {
	Code() string
	Theta0() float64
	Forward(phi, theta float64) (x, y float64, ok bool)
	Inverse(x, y float64) (phi, theta float64, ok bool)
}

// allProjectionCodes enumerates every code defined by Paper II (and the
// HEALPix extension paper) along with a bit indicating whether it is
// actually implemented by this package. Codes with implemented=false are
// still recognized — Select returns ErrUnsupportedProjection for them
// rather than treating them as unknown.
var allProjectionCodes = map[string]bool{
	// Zenithal (§5.1)
	"AZP": true, "SZP": true, "TAN": true, "STG": true, "SIN": true,
	"ARC": true, "ZPN": true, "ZEA": true,
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
	// HEALPix (Calabretta & Roukema 2007)
	"HPX": true, "XPH": true,
}

// Select returns the Projection for a 3-letter FITS projection code,
// configured from the projection-parameter map pv (PVi_m values keyed by
// (axis, index) — conventionally axis 2 for the latitude axis). The
// latAxis argument is the 1-based latitude axis index; most callers pass
// 2 (the common convention where CRVAL2 is declination).
//
// Returns ErrUnknownProjection for codes not in the FITS standard list.
// All 27 standard projections plus HPX/XPH are implemented, so
// ErrUnsupportedProjection is currently unreachable in normal use and
// exists only as a forward-compatibility sentinel.
func Select(code string, pv map[wcs.PVKey]float64, latAxis int) (Projection, error) {
	switch code {
	// Zenithal (§5.1)
	case "AZP":
		return newAZP(pv, latAxis), nil
	case "SZP":
		return newSZP(pv, latAxis), nil
	case "TAN":
		return tanProjection{}, nil
	case "STG":
		return stgProjection{}, nil
	case "SIN":
		return newSIN(pv, latAxis), nil
	case "ARC":
		return arcProjection{}, nil
	case "ZPN":
		return newZPN(pv, latAxis)
	case "ZEA":
		return zeaProjection{}, nil
	// Cylindrical (§5.2)
	case "CYP":
		return newCYP(pv, latAxis), nil
	case "CEA":
		return newCEA(pv, latAxis), nil
	case "CAR":
		return carProjection{}, nil
	case "MER":
		return merProjection{}, nil
	// Pseudo-cylindrical (§5.3)
	case "SFL":
		return sflProjection{}, nil
	case "PAR":
		return parProjection{}, nil
	case "MOL":
		return molProjection{}, nil
	case "AIT":
		return aitProjection{}, nil
	// Conic (§5.4)
	case "COP":
		return newCOP(pv, latAxis)
	case "COE":
		return newCOE(pv, latAxis)
	case "COD":
		return newCOD(pv, latAxis)
	case "COO":
		return newCOO(pv, latAxis)
	// Polyconic (§5.5)
	case "BON":
		return newBON(pv, latAxis), nil
	case "PCO":
		return pcoProjection{}, nil
	// Quad-cube (§5.6)
	case "TSC":
		return tscProjection{}, nil
	case "CSC":
		return cscProjection{}, nil
	case "QSC":
		return qscProjection{}, nil
	// HEALPix (CR 2007)
	case "HPX":
		return newHPX(pv, latAxis), nil
	case "XPH":
		return xphProjection{}, nil
	}
	if allProjectionCodes[code] {
		return nil, fmt.Errorf("%w: %q (known but not yet implemented)", ErrUnsupportedProjection, code)
	}
	return nil, fmt.Errorf("%w: %q", ErrUnknownProjection, code)
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

// normalizeRadPi maps a radian value into [-pi, pi].
func normalizeRadPi(r float64) float64 {
	r = math.Mod(r, 2*math.Pi)
	if r > math.Pi {
		r -= 2 * math.Pi
	} else if r < -math.Pi {
		r += 2 * math.Pi
	}
	return r
}
