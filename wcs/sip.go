package wcs

import (
	"fmt"
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

// MaxSIPOrder is the largest SIP polynomial order ParseSIP accepts. Real
// solutions stay well below it (HST, Spitzer and astrometry.net use 2-5);
// the bound keeps a corrupt order from forcing a huge allocation.
const MaxSIPOrder = 20

// ParseSIP reads the SIP keywords (A_ORDER, A_i_j, B_ORDER, B_i_j and the
// optional AP_/BP_ inverse) from h. It returns nil, nil if A_ORDER or
// B_ORDER is absent, and an error if any declared order is outside
// 0..MaxSIPOrder. It does not look at CTYPE; Parse calls it only when the
// celestial CTYPEs carry the "-SIP" suffix.
func ParseSIP(h *header.Header) (*SIP, error) {
	aOrder, errA := h.Int("A_ORDER")
	bOrder, errB := h.Int("B_ORDER")
	if errA != nil || errB != nil {
		return nil, nil
	}
	sip := &SIP{}
	var err error
	if sip.A, err = readSIPPoly(h, "A", aOrder); err != nil {
		return nil, err
	}
	if sip.B, err = readSIPPoly(h, "B", bOrder); err != nil {
		return nil, err
	}
	if apOrder, err := h.Int("AP_ORDER"); err == nil {
		if bpOrder, err := h.Int("BP_ORDER"); err == nil {
			if sip.AP, err = readSIPPoly(h, "AP", apOrder); err != nil {
				return nil, err
			}
			if sip.BP, err = readSIPPoly(h, "BP", bpOrder); err != nil {
				return nil, err
			}
		}
	}
	return sip, nil
}

// readSIPPoly reads a polynomial from the header given its prefix ("A",
// "B", "AP", "BP") and its declared order.
func readSIPPoly(h *header.Header, prefix string, declared int64) (*SIPPoly, error) {
	if declared < 0 || declared > MaxSIPOrder {
		return nil, fmt.Errorf("fits/wcs: %s_ORDER = %d outside 0..%d", prefix, declared, MaxSIPOrder)
	}
	order := int(declared)
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
	return p, nil
}
