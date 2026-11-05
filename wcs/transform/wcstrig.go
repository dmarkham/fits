package transform

import "math"

// Degree-based trigonometry with wcslib's exact special cases.
//
// The projection code in this package is a port of wcslib's prj.c, which
// works in degrees and relies on wcstrig.c returning EXACT values at
// multiples of 90 degrees (cosd(90) == 0, not 6.1e-17) and on its
// tolerant inverse functions. Whether a point is on or off a projection
// (TAN at theta = 0, AZP/SZP at theta = -90, a conic at |theta - theta_a| =
// 90, ...) is decided by those exact values, so reproducing wcslib's
// accept/reject decisions requires reproducing them. Each *d function below
// is a line-for-line port of the wcslib 8.10 wcstrig.c function of the
// same name. Where wcstrig.c calls C's sin/cos/tan/atan/atan2 they call
// Go's math.Sin/Cos/Sincos/Tan/Atan/Atan2; for asin and acos they call
// libAsin and libAcos below (cancellation-free identities on top of
// math.Atan2 and math.Atan), because math.Asin and math.Acos lose most of
// their bits next to +-1.
//
// None of this is correctly rounded. Go's libm and glibc's (which wcslib
// calls) differ in the last bit for a small fraction of arguments, so a
// value can differ from wcslib's by a few ULP, and by more where the
// projection itself is ill-conditioned (next to a pole, a horizon or a
// cube-face centre). Values are checked against wcslib within stated
// budgets (docs/VALIDATION.md), and accept/reject decisions are identical
// except at boundary-ambiguous points: points within a few ULP of an
// accept/reject boundary, where our own answer flips when the input is
// nudged by a few ULP, so neither answer is more right than the other.

// wpi, d2r and r2d are wcsmath.h's PI, D2R = PI/180.0 and R2D = 180.0/PI.
// They are variables, not constants, on purpose: C evaluates PI/180.0 in
// double arithmetic on the already-rounded PI, whereas a Go constant
// expression is evaluated exactly and rounded once, which can differ in
// the last bit. Ported code uses wpi (never math.Pi) in any expression
// wcslib writes with PI.
var (
	wpi = float64(math.Pi)
	d2r = wpi / 180.0
	r2d = 180.0 / wpi
)

const (
	// wcstrigTol is WCSTRIG_TOL (wcstrig.h): how far outside [-1, 1] the
	// argument of acosd/asind may stray and still be treated as +-1.
	wcstrigTol = 1e-10
)

// libm calls. Every sin, cos, tan, atan, atan2, log, exp and pow the
// projections and the celestial rotation evaluate goes through one of the
// lm* functions below, which return Go's math result unchanged unless a
// test has set libmPerturbULP. Then each result is moved by that many ULP
// (alternating in sign from call to call when libmPerturbAlt is set), so
// a test can measure how far a given output moves when the libm results
// it is built from differ in their last bits, which is all that separates
// this port from wcslib calling glibc (see libmSensitivity in the tests).
// Exact zeros and infinities are left alone. In production the cost is one
// load and a predicted branch per call.
var (
	libmPerturbULP int64
	libmPerturbAlt bool
	libmPerturbN   int64
)

func lmP(v float64) float64 {
	if libmPerturbULP == 0 {
		return v
	}
	return lmPerturb(v)
}

func lmPerturb(v float64) float64 {
	if v == 0 || math.IsInf(v, 0) || v != v {
		return v
	}
	k := libmPerturbULP
	if libmPerturbAlt {
		libmPerturbN++
		if libmPerturbN&1 == 1 {
			k = -k
		}
	}
	// Adding to the bit pattern moves |v| by k ULP (k > 0: away from 0).
	b := int64(math.Float64bits(v)) + k
	if r := math.Float64frombits(uint64(b)); r == r && !math.IsInf(r, 0) && (r < 0) == (v < 0) {
		return r
	}
	return v
}

func lmSin(x float64) float64           { return lmP(math.Sin(x)) }
func lmCos(x float64) float64           { return lmP(math.Cos(x)) }
func lmTan(x float64) float64           { return lmP(math.Tan(x)) }
func lmAtan(x float64) float64          { return lmP(math.Atan(x)) }
func lmAtan2(y, x float64) float64      { return lmP(math.Atan2(y, x)) }
func lmLog(x float64) float64           { return lmP(math.Log(x)) }
func lmExp(x float64) float64           { return lmP(math.Exp(x)) }
func lmPow(x, y float64) float64        { return lmP(math.Pow(x, y)) }
func lmSincos(x float64) (s, c float64) { s, c = math.Sincos(x); return lmP(s), lmP(c) }

func cosd(angle float64) float64 {
	if math.Mod(angle, 90.0) == 0.0 {
		switch abs(int(math.Floor(angle/90.0+0.5))) % 4 {
		case 0:
			return 1.0
		case 1:
			return 0.0
		case 2:
			return -1.0
		case 3:
			return 0.0
		}
	}
	return lmCos(float64(angle * d2r))
}

func sind(angle float64) float64 {
	if math.Mod(angle, 90.0) == 0.0 {
		switch abs(int(math.Floor(angle/90.0-0.5))) % 4 {
		case 0:
			return 1.0
		case 1:
			return 0.0
		case 2:
			return -1.0
		case 3:
			return 0.0
		}
	}
	return lmSin(float64(angle * d2r))
}

func sincosd(angle float64) (s, c float64) {
	if math.Mod(angle, 90.0) == 0.0 {
		switch abs(int(math.Floor(angle/90.0+0.5))) % 4 {
		case 0:
			return 0.0, 1.0
		case 1:
			if angle > 0.0 {
				return 1.0, 0.0
			}
			return -1.0, 0.0
		case 2:
			return 0.0, -1.0
		case 3:
			if angle > 0.0 {
				return -1.0, 0.0
			}
			return 1.0, 0.0
		}
	}
	return lmSincos(float64(angle * d2r))
}

func tand(angle float64) float64 {
	resid := math.Mod(angle, 360.0)
	switch {
	case resid == 0.0 || math.Abs(resid) == 180.0:
		return 0.0
	case resid == 45.0 || resid == 225.0:
		return 1.0
	case resid == -135.0 || resid == -315.0:
		return -1.0
	}
	return lmTan(float64(angle * d2r))
}

// libAsin and libAcos stand in for C's asin and acos. Go's math.Asin and
// math.Acos form 1 - x*x, which cancels near |x| = 1: against glibc 2.44
// math.Asin is off by up to 1024 ULP and math.Acos by up to 2.5e7 ULP
// (acos(0.9999999995) is wrong in its 10th significant digit), and near a
// pole every projection inverse takes an asin or acos of an argument next
// to +-1 (MOL even takes asin of a value built from another asin). These
// form the exact differences 1 - x and 1 + x instead (Sterbenz) and
// call math.Atan2 / math.Atan, which keeps them within a few ULP of the
// exact value everywhere (TestLibmAccuracy) at about the cost of
// math.Asin / math.Acos. NaN outside [-1, 1], like C.
func libAsin(x float64) float64 {
	return lmAtan2(x, math.Sqrt(float64((1-x)*(1+x))))
}

func libAcos(x float64) float64 {
	return float64(2 * lmAtan(math.Sqrt((1-x)/(1+x))))
}

func acosd(v float64) float64 {
	if v >= 1.0 {
		if v-1.0 < wcstrigTol {
			return 0.0
		}
	} else if v == 0.0 {
		return 90.0
	} else if v <= -1.0 {
		if v+1.0 > -wcstrigTol {
			return 180.0
		}
	}
	return float64(libAcos(v) * r2d)
}

func asind(v float64) float64 {
	if v <= -1.0 {
		if v+1.0 > -wcstrigTol {
			return -90.0
		}
	} else if v == 0.0 {
		return 0.0
	} else if v >= 1.0 {
		if v-1.0 < wcstrigTol {
			return 90.0
		}
	}
	return float64(libAsin(v) * r2d)
}

func atand(v float64) float64 {
	switch v {
	case -1.0:
		return -45.0
	case 0.0:
		return 0.0
	case 1.0:
		return 45.0
	}
	return float64(lmAtan(v) * r2d)
}

func atan2d(y, x float64) float64 {
	if y == 0.0 {
		if x >= 0.0 {
			return 0.0
		} else if x < 0.0 {
			return 180.0
		}
	} else if x == 0.0 {
		if y > 0.0 {
			return 90.0
		} else if y < 0.0 {
			return -90.0
		}
	}
	return float64(lmAtan2(y, x) * r2d)
}

func abs(i int) int {
	if i < 0 {
		return -i
	}
	return i
}

// prjbchk is wcslib's prjbchk for one point: native longitude outside
// [-180, 180] or latitude outside [-90, 90] by more than tol (degrees) is
// rejected; within tol it is clamped onto the boundary.
func prjbchk(tol float64, phi, theta *float64) bool {
	if *phi < -180.0 {
		if *phi < -180.0-tol {
			return false
		}
		*phi = -180.0
	} else if 180.0 < *phi {
		if 180.0+tol < *phi {
			return false
		}
		*phi = 180.0
	}
	if *theta < -90.0 {
		if *theta < -90.0-tol {
			return false
		}
		*theta = -90.0
	} else if 90.0 < *theta {
		if 90.0+tol < *theta {
			return false
		}
		*theta = 90.0
	}
	return true
}
