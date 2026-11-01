package wcs

import (
	"strconv"

	"github.com/dmarkham/fits/header"
)

// SIPPoly holds one Simple Imaging Polynomial (Shupe et al. 2005) as a
// coefficient grid: Coeffs[i][j] is the coefficient of u^i · v^j, with
// i + j <= Order.
type SIPPoly struct {
	Order  int
	Coeffs [][]float64
}

// SIP holds the SIP distortion keywords of an image header: the forward
// polynomials A and B and, when present, the inverse polynomials AP and
// BP (nil otherwise).
type SIP struct {
	A, B   *SIPPoly
	AP, BP *SIPPoly
}

// ParseSIP reads the SIP keywords (A_ORDER, A_i_j, B_ORDER, B_i_j and the
// optional AP_/BP_ inverse) from h. It returns nil if A_ORDER or B_ORDER
// is absent. It does not look at CTYPE; Parse calls it only when the
// celestial CTYPEs carry the "-SIP" suffix.
func ParseSIP(h *header.Header) *SIP {
	aOrder, errA := h.Int("A_ORDER")
	bOrder, errB := h.Int("B_ORDER")
	if errA != nil || errB != nil {
		return nil
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
	return sip
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
			key := prefix + "_" + strconv.Itoa(i) + "_" + strconv.Itoa(j)
			if v, err := h.Float(key); err == nil {
				p.Coeffs[i][j] = v
			}
		}
	}
	return p
}
