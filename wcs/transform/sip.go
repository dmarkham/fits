package transform

import (
	"math"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// SIP — Simple Imaging Polynomial distortion (Shupe et al. 2005).
//
// SIP adds a polynomial correction to the pixel offset from CRPIX, before
// the linear (CD or PC+CDELT) step of the WCS pipeline. HST, Spitzer, JWST
// and astrometry.net headers use it; real headers of each are in
// testdata/wref/hdrcases (real_*) and checked against wcslib. HST FLT
// headers also carry lookup-table distortions (CPDISj/D2IMDISj) that this
// package does not implement; New rejects such headers.
//
// # Header representation
//
// The SIP keywords are:
//
//	CTYPE1 = 'RA---TAN-SIP'
//	CTYPE2 = 'DEC--TAN-SIP'
//	A_ORDER = n       / polynomial order for the x direction
//	B_ORDER = n       / polynomial order for the y direction
//	A_i_j   = coefficient of u^i * v^j in the x correction
//	B_i_j   = coefficient of u^i * v^j in the y correction
//	AP_ORDER = n      / optional inverse polynomial order (x direction)
//	BP_ORDER = n      / optional inverse polynomial order (y direction)
//	AP_i_j  = inverse x-correction coefficient
//	BP_i_j  = inverse y-correction coefficient
//
// where (u, v) = (pixel - CRPIX).
//
// # Forward direction
//
// Given image pixel (u, v):
//
//	f(u, v) = sum over i, j of A_i_j * u^i * v^j
//	g(u, v) = sum over i, j of B_i_j * u^i * v^j
//	(u', v') = (u + f(u, v), v + g(u, v))
//
// The corrected (u', v') is then fed into the linear CD matrix step.
//
// # Inverse direction
//
// The forward polynomial is inverted by Newton iteration, started from
// AP/BP when present (they are only an approximate inverse), exactly as
// wcslib disx2p does. See (*SIP).Inverse.
//
// # Parsing
//
// wcs.Parse reads the SIP keywords into wcs.Header.SIP when the celestial
// CTYPEs end in "-SIP", and New applies them automatically. ParseSIP reads
// them directly from a *header.Header for callers that apply SIP
// themselves.

// SIPPoly holds a SIP polynomial as a 2D coefficient grid indexed by
// (i, j), where Coeffs[i][j] is the coefficient of u^i * v^j. Order is
// the maximum total degree (i + j <= Order).
type SIPPoly struct {
	Order  int
	Coeffs [][]float64 // [i][j], 0-indexed; entries beyond Order are zero
}

// eval returns the polynomial value at (u, v).
func (p *SIPPoly) eval(u, v float64) float64 {
	val, _, _ := p.evalWithDeriv(u, v)
	return val
}

// evalWithDeriv returns the polynomial value and its partial derivatives
// ∂/∂u and ∂/∂v at (u, v). For a monomial c·u^i·v^j:
//
//	∂/∂u = c·i·u^(i-1)·v^j      (zero if i=0)
//	∂/∂v = c·u^i·j·v^(j-1)      (zero if j=0)
//
// Both derivatives are computed in the same walk as the value for
// no extra allocation cost.
func (p *SIPPoly) evalWithDeriv(u, v float64) (val, dvdu, dvdv float64) {
	if p == nil || p.Order < 0 {
		return 0, 0, 0
	}
	// Power tables; uPow[k] = u^k, vPow[k] = v^k.
	uPow := make([]float64, p.Order+1)
	vPow := make([]float64, p.Order+1)
	uPow[0] = 1
	vPow[0] = 1
	for k := 1; k <= p.Order; k++ {
		uPow[k] = float64(uPow[k-1] * u)
		vPow[k] = float64(vPow[k-1] * v)
	}
	for i := 0; i <= p.Order; i++ {
		for j := 0; j+i <= p.Order; j++ {
			c := p.Coeffs[i][j]
			if c == 0 {
				continue
			}
			val += float64(c * uPow[i] * vPow[j])
			if i >= 1 {
				dvdu += float64(c * float64(i) * uPow[i-1] * vPow[j])
			}
			if j >= 1 {
				dvdv += float64(c * uPow[i] * float64(j) * vPow[j-1])
			}
		}
	}
	return
}

// SIP holds a parsed SIP distortion: forward (A, B) and optional inverse
// (AP, BP) polynomial pairs.
type SIP struct {
	A  *SIPPoly
	B  *SIPPoly
	AP *SIPPoly // may be nil (solved iteratively if absent)
	BP *SIPPoly // may be nil
}

// ParseSIP extracts SIP coefficients from the given image header; it does
// not look at CTYPE. It returns nil, nil only when the header carries no
// SIP keyword at all, and an error (see wcs.ParseSIP) when any SIP keyword
// is present but A_ORDER or B_ORDER is missing, for an order outside
// 0..wcs.MaxSIPOrder (including negative orders), a non-numeric
// coefficient, or a non-zero coefficient beyond its declared order.
//
// wcs.Parse already puts these coefficients in wcs.Header.SIP whenever the
// primary description carries SIP keywords (with or without "-SIP" in
// CTYPE, as wcslib does; not on a RA---TPV header without "-SIP", where
// they are ignored as astropy ignores them, wcs.Header.SIPIgnored), and
// New applies them; call ParseSIP only when applying SIP yourself.
func ParseSIP(h *header.Header) (*SIP, error) {
	s, err := wcs.ParseSIP(h)
	if err != nil {
		return nil, err
	}
	return sipFromWCS(s), nil
}

// sipFromWCS converts the parsed keywords into a *SIP, or nil.
func sipFromWCS(s *wcs.SIP) *SIP {
	if s == nil {
		return nil
	}
	poly := func(p *wcs.SIPPoly) *SIPPoly {
		if p == nil {
			return nil
		}
		return &SIPPoly{Order: p.Order, Coeffs: p.Coeffs}
	}
	return &SIP{A: poly(s.A), B: poly(s.B), AP: poly(s.AP), BP: poly(s.BP)}
}

// Forward applies the SIP forward correction to a pixel offset (u, v)
// from CRPIX. Returns the corrected (u', v') that should be fed into the
// linear CD step.
func (s *SIP) Forward(u, v float64) (uOut, vOut float64) {
	du := s.A.eval(u, v)
	dv := s.B.eval(u, v)
	return u + du, v + dv
}

// Inverse inverts the SIP forward map F(u, v) = (u + A(u, v), v + B(u, v))
// for a corrected offset (u', v'), as wcslib disx2p does: AP/BP, when
// present, give only the STARTING point (they are a fitted approximation,
// not an exact inverse), and Newton iteration on A/B with the analytic
// Jacobian
//
//	J = I + [∂A/∂u  ∂A/∂v]
//	        [∂B/∂u  ∂B/∂v]
//
// refines it to round-off (see newton2 for the convergence test). Without
// AP/BP the start is (u', v'). ok is false when the iteration does not
// converge (e.g. a point outside the region where the polynomial is
// invertible); the returned (u, v) are then meaningless.
func (s *SIP) Inverse(uPrime, vPrime float64) (u, v float64, ok bool) {
	u0, v0 := uPrime, vPrime
	if s.AP != nil && s.BP != nil {
		u0 += s.AP.eval(uPrime, vPrime)
		v0 += s.BP.eval(uPrime, vPrime)
	}
	return newton2(func(u, v float64) (fu, fv, j11, j12, j21, j22 float64) {
		aVal, dAdu, dAdv := s.A.evalWithDeriv(u, v)
		bVal, dBdu, dBdv := s.B.evalWithDeriv(u, v)
		return u + aVal, v + bVal, 1 + dAdu, dAdv, dBdu, 1 + dBdv
	}, uPrime, vPrime, u0, v0, 1)
}

// newton2 solves F(x, y) = (tx, ty) by Newton iteration from (x0, y0). f
// returns F and its Jacobian [[j11 j12] [j21 j22]].
//
// Convergence is RELATIVE to the size of the solution: the iteration
// stops when the Newton step is at most newtonRelTol * max(|x|, |y|,
// floor), i.e. a few tens of ULP of the answer (floor is the coordinate
// scale below which an absolute test applies: 1 px for SIP). A fixed
// absolute tolerance does not work here: 1e-13 px is already below one
// ULP of a 2640 px offset, so such a test can never pass on a real
// detector and used to be "passed" by returning ok after the last
// iteration regardless. If the step stops shrinking while already at the
// round-off level (newtonStallTol) that is also accepted. Anything else
// -- a singular or non-finite Jacobian, a non-finite iterate, or
// newtonMaxIter iterations without convergence -- returns ok = false.
func newton2(f func(x, y float64) (fx, fy, j11, j12, j21, j22 float64), tx, ty, x0, y0, floor float64) (x, y float64, ok bool) {
	x, y = x0, y0
	prev := math.Inf(1)
	for iter := 0; iter < newtonMaxIter; iter++ {
		fx, fy, j11, j12, j21, j22 := f(x, y)
		rx, ry := tx-fx, ty-fy
		det := float64(j11*j22) - float64(j12*j21)
		if det == 0 || math.IsNaN(det) || math.IsInf(det, 0) {
			return x, y, false
		}
		dx := (float64(j22*rx) - float64(j12*ry)) / det
		dy := (float64(-j21*rx) + float64(j11*ry)) / det
		x += dx
		y += dy
		if math.IsNaN(x) || math.IsNaN(y) || math.IsInf(x, 0) || math.IsInf(y, 0) {
			return x, y, false
		}
		step := math.Max(math.Abs(dx), math.Abs(dy))
		scale := math.Max(math.Max(math.Abs(x), math.Abs(y)), floor)
		if step <= float64(newtonRelTol*scale) {
			return x, y, true
		}
		if iter >= 2 && step >= prev && step <= float64(newtonStallTol*scale) {
			return x, y, true
		}
		prev = step
	}
	return x, y, false
}

const (
	newtonMaxIter  = 100
	newtonRelTol   = 1e-14 // ~45 ULP of the solution
	newtonStallTol = 1e-11 // round-off plateau accepted only below this
)
