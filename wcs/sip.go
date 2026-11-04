package wcs

import (
	"fmt"
	"regexp"
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
// polynomials A and B and, when both AP_ORDER and BP_ORDER are present,
// the inverse polynomials AP and BP (nil otherwise). AP/BP are only ever
// an approximation of the inverse; wcs/transform uses them as the
// starting point of an iterative inversion of A/B, as wcslib does.
type SIP struct {
	A, B   *SIPPoly
	AP, BP *SIPPoly
}

// MaxSIPOrder is the highest SIP polynomial order accepted. The SIP
// convention sets no limit; 9 is the highest degree wcslib's TPD
// implementation of SIP supports ("Degree must not exceed 9"), and no
// published SIP solution comes close to it.
const MaxSIPOrder = 9

// sipKeyRe matches SIP coefficient keywords A_p_q, B_p_q, AP_p_q, BP_p_q.
// It also matches multi-digit and zero-padded indices (A_01_1, B_1_10),
// which ParseSIP refuses: the SIP convention and wcslib's lexer (wcspih.l
// <SIP2>{Z1}_{Z1}) use single digits, and wcslib silently discards such
// a card.
var sipKeyRe = regexp.MustCompile(`^(A|B|AP|BP)_([0-9]+)_([0-9]+)$`)

// ParseSIP reads the SIP keywords (A_ORDER, A_i_j, B_ORDER, B_i_j and the
// optional AP_/BP_ inverse) from h. It returns nil, nil when h carries no
// SIP keyword at all (no *_ORDER and no A_/B_/AP_/BP_p_q card). It does not
// look at CTYPE; Parse calls it for the primary description whenever the
// header carries a SIP keyword, with or without the "-SIP" suffix, as
// wcslib does, except on a RA---TPV / DEC--TPV header without "-SIP",
// whose SIP keywords Parse ignores as astropy does (Header.SIPIgnored).
//
// It returns an error (never panics) when:
//
//   - a coefficient card has p + q > MaxSIPOrder (A_5_5, BP_9_1, ...),
//     even with value zero: wcslib 8.10 refuses it too ("Invalid
//     powers"; hdrcases_reject/sip_*_deg10);
//   - a coefficient keyword has an index that is not a single digit
//     (A_01_1, B_0_02, A_1_10): wcslib 8.10 discards such a card without
//     a message, so the header means something different to wcslib than
//     its author intended;
//   - any SIP keyword is present but A_ORDER or B_ORDER is missing;
//   - an *_ORDER keyword is not an integer in 0..MaxSIPOrder (10, -1, or
//     2.0 written as a real);
//   - a coefficient is not numeric;
//   - a NON-ZERO coefficient has i + j greater than its polynomial's
//     order (A_3_0 with A_ORDER = 2).
//
// Apart from the first, these are STRICTER than wcslib 8.10, which never
// reads A_ORDER/B_ORDER (sipset takes each degree from the coefficients
// present) and so accepts and applies such headers, reading a missing B
// polynomial as zero, applying coefficients above or without a declared
// order and skipping an unreadable coefficient. Each refusal has a
// wcslib-accepts case in testdata/wref/hdrcases_stricter
// (TestHdrStricterThanWcslib): sip_index_leading_zero, sip_a_without_b,
// sip_coeffs_without_orders, sip_order10, sip_neg_order, sip_order_float,
// sip_coef_string, sip_beyond_order_nonzero.
//
// The inverse AP/BP is kept only when both AP_ORDER and BP_ORDER are
// present. It is never required for correctness: transform uses it only
// as the starting guess of an iteration on A/B (as wcslib does), so an
// incomplete or missing inverse changes no result. wcslib 8.10 instead
// refuses a header whose AP or BP has no coefficient card ("Degree must
// not exceed 9, got -1"), including every SIP header without AP/BP; that
// is a deliberate difference (testdata/wref/hdrcases_lenient).
func ParseSIP(h *header.Header) (*SIP, error) {
	if !HasSIPKeywords(h) {
		return nil, nil
	}
	for _, c := range h.Cards() {
		m := sipKeyRe.FindStringSubmatch(c.Key)
		if m == nil {
			continue
		}
		if len(m[2]) != 1 || len(m[3]) != 1 {
			return nil, fmt.Errorf("fits/wcs: %s is not a SIP coefficient keyword (SIP indices are single digits; wcslib 8.10 ignores this card)", c.Key)
		}
		// wcslib 8.10 sipset refuses any SIP coefficient card of degree
		// p + q > 9 ("Invalid powers"), forward or inverse, whatever its
		// value or *_ORDER.
		if d := int(m[2][0]-'0') + int(m[3][0]-'0'); d > MaxSIPOrder {
			return nil, fmt.Errorf("fits/wcs: %s has degree %d > %d (wcslib 8.10 refuses it: Invalid powers)", c.Key, d, MaxSIPOrder)
		}
	}
	for _, k := range []string{"A_ORDER", "B_ORDER"} {
		if !h.Has(k) {
			return nil, fmt.Errorf("fits/wcs: SIP keywords present but %s is missing", k)
		}
	}
	order := func(key string) (int, error) {
		c, _ := h.Get(key)
		v, ok := c.Value.(int64)
		if !ok {
			return 0, fmt.Errorf("fits/wcs: %s = %v is not an integer", key, c.Value)
		}
		if v < 0 || v > MaxSIPOrder {
			return 0, fmt.Errorf("fits/wcs: %s = %d is outside 0..%d", key, v, MaxSIPOrder)
		}
		return int(v), nil
	}
	orders := map[string]int{}
	for _, p := range []string{"A", "B", "AP", "BP"} {
		if !h.Has(p + "_ORDER") {
			continue
		}
		n, err := order(p + "_ORDER")
		if err != nil {
			return nil, err
		}
		orders[p] = n
	}
	polys := map[string]*SIPPoly{}
	for p, n := range orders {
		polys[p] = newSIPPoly(n)
	}
	for _, c := range h.Cards() {
		m := sipKeyRe.FindStringSubmatch(c.Key)
		if m == nil {
			continue
		}
		p := polys[m[1]]
		if p == nil {
			continue // AP_p_q without AP_ORDER: no inverse polynomial declared
		}
		v, ok := floatValue(c.Value)
		if !ok {
			return nil, fmt.Errorf("fits/wcs: %s has a non-numeric value %v", c.Key, c.Value)
		}
		i, _ := strconv.Atoi(m[2])
		j, _ := strconv.Atoi(m[3])
		if i+j > p.Order {
			if v != 0 {
				return nil, fmt.Errorf("fits/wcs: %s = %g has degree %d > %s_ORDER = %d", c.Key, v, i+j, m[1], p.Order)
			}
			continue
		}
		p.Coeffs[i][j] = v
	}
	sip := &SIP{A: polys["A"], B: polys["B"]}
	if polys["AP"] != nil && polys["BP"] != nil {
		sip.AP, sip.BP = polys["AP"], polys["BP"]
	}
	return sip, nil
}

// HasSIPKeywords reports whether h carries any SIP keyword: A_ORDER,
// B_ORDER, AP_ORDER, BP_ORDER or an A_p_q / B_p_q / AP_p_q / BP_p_q
// coefficient.
func HasSIPKeywords(h *header.Header) bool {
	for _, c := range h.Cards() {
		switch c.Key {
		case "A_ORDER", "B_ORDER", "AP_ORDER", "BP_ORDER":
			return true
		}
		if sipKeyRe.MatchString(c.Key) {
			return true
		}
	}
	return false
}

func newSIPPoly(order int) *SIPPoly {
	p := &SIPPoly{Order: order, Coeffs: make([][]float64, order+1)}
	for i := range p.Coeffs {
		p.Coeffs[i] = make([]float64, order+1)
	}
	return p
}
