package transform

import (
	"math"

	"github.com/dmarkham/fits/header"
)

// SIP — Simple Imaging Polynomial distortion (Shupe et al. 2005).
//
// SIP adds a polynomial correction between the linear CD matrix step and
// the projection step of the WCS pipeline. It is the dominant distortion
// convention for HST, Spitzer, and JWST image data, and is required to
// correctly map pixels to sky for any of those missions.
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
// If AP/BP are present, the inverse polynomial is applied directly.
// Otherwise we solve the forward polynomial via Newton iteration.
//
// # Parsing
//
// SIP coefficients are read from an image HDU's header outside of the
// standard WCS keyword set. Because they are image-specific and not part
// of Paper I, the wcs package does not parse them automatically. Instead,
// ParseSIP reads them directly from a *header.Header.

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
		uPow[k] = uPow[k-1] * u
		vPow[k] = vPow[k-1] * v
	}
	for i := 0; i <= p.Order; i++ {
		for j := 0; j+i <= p.Order; j++ {
			c := p.Coeffs[i][j]
			if c == 0 {
				continue
			}
			val += c * uPow[i] * vPow[j]
			if i >= 1 {
				dvdu += c * float64(i) * uPow[i-1] * vPow[j]
			}
			if j >= 1 {
				dvdv += c * uPow[i] * float64(j) * vPow[j-1]
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

// ParseSIP extracts SIP coefficients from the given image header. Returns
// nil, nil if the header does not declare SIP (i.e. CTYPE suffix is not
// -SIP, or A_ORDER / B_ORDER is absent). An error indicates a malformed
// SIP section.
func ParseSIP(h *header.Header) (*SIP, error) {
	aOrder, errA := h.Int("A_ORDER")
	bOrder, errB := h.Int("B_ORDER")
	if errA != nil || errB != nil {
		return nil, nil
	}
	sip := &SIP{
		A: readSIPPoly(h, "A", int(aOrder)),
		B: readSIPPoly(h, "B", int(bOrder)),
	}
	if apOrder, err := h.Int("AP_ORDER"); err == nil {
		if bpOrder, err := h.Int("BP_ORDER"); err == nil {
			sip.AP = readSIPPoly(h, "AP", int(apOrder))
			sip.BP = readSIPPoly(h, "BP", int(bpOrder))
		}
	}
	return sip, nil
}

// readSIPPoly reads a polynomial from the header given its prefix ("A",
// "B", "AP", "BP") and its declared order.
func readSIPPoly(h *header.Header, prefix string, order int) *SIPPoly {
	p := &SIPPoly{
		Order:  order,
		Coeffs: make([][]float64, order+1),
	}
	for i := 0; i <= order; i++ {
		p.Coeffs[i] = make([]float64, order+1)
		for j := 0; j+i <= order; j++ {
			key := prefix + "_" + itoa(i) + "_" + itoa(j)
			if v, err := h.Float(key); err == nil {
				p.Coeffs[i][j] = v
			}
		}
	}
	return p
}

// itoa is a tiny helper that avoids importing strconv just for this.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var buf [8]byte
	i := len(buf)
	neg := n < 0
	if neg {
		n = -n
	}
	for n > 0 {
		i--
		buf[i] = '0' + byte(n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}

// Forward applies the SIP forward correction to a pixel offset (u, v)
// from CRPIX. Returns the corrected (u', v') that should be fed into the
// linear CD step.
func (s *SIP) Forward(u, v float64) (uOut, vOut float64) {
	du := s.A.eval(u, v)
	dv := s.B.eval(u, v)
	return u + du, v + dv
}

// Inverse applies the SIP inverse correction. If AP/BP are present, it
// evaluates them directly. Otherwise it iterates Newton's method against
// the forward polynomial using the analytic Jacobian.
//
// The forward map is F(u, v) = (u + A(u, v), v + B(u, v)), so the
// Jacobian is
//
//	J = I + [∂A/∂u  ∂A/∂v]
//	        [∂B/∂u  ∂B/∂v]
//
// Both partial derivatives are computed in a single walk over the
// A/B coefficient tables, making each iteration cost the same as a
// forward evaluation (vs. 3× cost for the old numerical version).
func (s *SIP) Inverse(uPrime, vPrime float64) (u, v float64, ok bool) {
	if s.AP != nil && s.BP != nil {
		du := s.AP.eval(uPrime, vPrime)
		dv := s.BP.eval(uPrime, vPrime)
		return uPrime + du, vPrime + dv, true
	}
	// Newton iteration with analytic Jacobian.
	u, v = uPrime, vPrime
	for range 50 {
		aVal, dAdu, dAdv := s.A.evalWithDeriv(u, v)
		bVal, dBdu, dBdv := s.B.evalWithDeriv(u, v)
		fu := u + aVal
		fv := v + bVal
		ru := uPrime - fu
		rv := vPrime - fv
		if math.Abs(ru) < 1e-13 && math.Abs(rv) < 1e-13 {
			return u, v, true
		}
		j11 := 1 + dAdu
		j12 := dAdv
		j21 := dBdu
		j22 := 1 + dBdv
		det := j11*j22 - j12*j21
		if det == 0 {
			return 0, 0, false
		}
		u += (j22*ru - j12*rv) / det
		v += (-j21*ru + j11*rv) / det
	}
	return u, v, true
}
