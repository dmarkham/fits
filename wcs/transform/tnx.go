package transform

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
)

// TNX — IRAF polynomial distortion of the gnomonic projection.
//
// TNX (FITS registry, https://fits.gsfc.nasa.gov/registry/tnx.html) stores
// its correction as free-form attribute strings in WAT1_nnn (longitude
// axis) and WAT2_nnn (latitude axis). IRAF (MSCRED, ccmap, ...) writes it;
// the real headers proven here are NOAO Mosaic II on the CTIO 4 m and
// Mosaic-1 on the KPNO 4 m (Chebyshev and plain polynomial surfaces);
// Legendre surfaces are proven on synthetic headers only, as no public
// real one was found (testdata/tnx/PROVENANCE.md). The WATn_nnn
// cards are joined in key-number order, wherever they sit in the header
// (the KPNO sample interleaves them with other keywords).
// The cards of one axis are joined VERBATIM — the string value of each card
// between its quotes, trailing blanks included, because IRAF cuts the
// attribute string every 68 characters with no regard for token boundaries
// (WCSTools hgetm does the same). Example (registry tnx doc, figure 1):
//
//	WAT1_003= '838772 0.01749134520424022 -0.01082784423020123 -0.1387962673564234 '
//	WAT1_004= '-4.307309762939804E-4 0.009069288008295441 0.002875265278754504 -0.0'
//
// The joined string is a list of "key = value" attributes, values quoted
// with double quotes when they contain blanks:
//
//	wtype=tnx axtype=ra lngcor = "3. 4. 4. 2. <xmin> <xmax> <ymin> <ymax> <coeffs>"
//
// A surface is: function type (1 = Chebyshev, 2 = Legendre, 3 = plain
// polynomial), x order, y order (number of terms, i.e. highest power + 1),
// cross terms (0 = none, 1 = full, 2 = half: total power < max order), the
// (xi, eta) domain, then the coefficients C(m, n) with m (xi power)
// varying fastest. Chebyshev and Legendre are evaluated on the domain
// normalised to [-1, 1]; the plain polynomial is evaluated on raw xi, eta
// (degrees). The corrected coordinates are xi + lngcor(xi, eta) and
// eta + latcor(xi, eta). All of this follows IRAF mwcs/wfgsurfit.x as
// translated in WCSTools tnxpos.c, and is proven against WCSTools 3.9.7
// (tnx_golden_test.go).
//
// Parsing is strict: anything WCSTools would read differently or not at
// all (wrong coefficient count, unknown function type or cross-term code,
// non-integer orders, empty domain, a non-TNX wtype such as zpx, a
// non-default "ro") is an error, never a guess.

// TNX holds a parsed TNX distortion for both axes.
type TNX struct {
	Lon tnxSurface // lngcor: correction added to xi
	Lat tnxSurface // latcor: correction added to eta

	// LonPole is the "longpole" WAT attribute (degrees) when present
	// (LonPoleSet). NewTNX uses it when the header has no LONPOLE card.
	LonPole    float64
	LonPoleSet bool
}

// tnxSurface is one axis's polynomial surface.
type tnxSurface struct {
	Type    int // 1=Chebyshev, 2=Legendre, 3=simple polynomial
	XOrder  int
	YOrder  int
	Cross   int // 0 = none, 1 = full, 2 = half
	XMin    float64
	XMax    float64
	YMin    float64
	YMax    float64
	Coeffs  []float64 // C(m, n), m fastest, terms filtered by Cross
	present bool
}

// ParseTNX reads the TNX correction from the WAT1_nnn / WAT2_nnn cards of
// a raw header. Returns (nil, nil) when neither axis has WAT cards. lngcor
// is taken from WAT1 (else WAT2) and latcor from WAT2 (else WAT1), as
// WCSTools tnxinit.
func ParseTNX(h *header.Header) (*TNX, error) {
	wat1, err := concatWAT(h, "WAT1")
	if err != nil {
		return nil, err
	}
	wat2, err := concatWAT(h, "WAT2")
	if err != nil {
		return nil, err
	}
	if wat1 == "" && wat2 == "" {
		return nil, nil
	}
	a1, err := parseWATAttrs(wat1)
	if err != nil {
		return nil, fmt.Errorf("wcs/transform: TNX WAT1: %w", err)
	}
	a2, err := parseWATAttrs(wat2)
	if err != nil {
		return nil, fmt.Errorf("wcs/transform: TNX WAT2: %w", err)
	}
	for i, a := range []map[string]string{a1, a2} {
		if wt, ok := a["wtype"]; ok && !strings.EqualFold(wt, "tnx") {
			return nil, fmt.Errorf("wcs/transform: TNX WAT%d declares wtype=%s; only tnx is supported", i+1, wt)
		}
		if ro, ok := a["ro"]; ok {
			v, err := strconv.ParseFloat(ro, 64)
			if err != nil || math.Abs(v-180/math.Pi) > 1e-12*(180/math.Pi) {
				return nil, fmt.Errorf("wcs/transform: TNX WAT%d ro=%s; only the default 180/pi is supported", i+1, ro)
			}
		}
	}
	t := &TNX{}
	for _, a := range []map[string]string{a1, a2} {
		if lp, ok := a["longpole"]; ok {
			v, err := strconv.ParseFloat(lp, 64)
			if err != nil {
				return nil, fmt.Errorf("wcs/transform: TNX longpole=%q: %w", lp, err)
			}
			t.LonPole, t.LonPoleSet = v, true
			break // WAT1 takes precedence (tnxinit)
		}
	}
	pick := func(tag string, first, second map[string]string) (string, bool) {
		if v, ok := first[tag]; ok {
			return v, true
		}
		v, ok := second[tag]
		return v, ok
	}
	if v, ok := pick("lngcor", a1, a2); ok {
		if t.Lon, err = parseTNXSurface(v); err != nil {
			return nil, fmt.Errorf("wcs/transform: TNX lngcor: %w", err)
		}
	}
	if v, ok := pick("latcor", a2, a1); ok {
		if t.Lat, err = parseTNXSurface(v); err != nil {
			return nil, fmt.Errorf("wcs/transform: TNX latcor: %w", err)
		}
	}
	return t, nil
}

// tnxMaxWAT is the longest joined WAT string NewTNX accepts: 1999
// characters. WCSTools 3.9.7 (the truth source) joins the cards into a
// 2000-byte buffer (tnxpos.c SZ_ATSTRING, hget.c hgetm) and silently cuts
// anything longer, so a longer string has no independent proof and is an
// error here. A golden at exactly 1999 characters (syn_leg_wat_1999)
// pins the boundary.
const tnxMaxWAT = 1999

// tnxMaxWATCards: hgetm reads WATn_1 .. WATn_499 and ignores higher key
// numbers; a header that has one is refused rather than read differently.
const tnxMaxWATCards = 499

// concatWAT joins prefix_001, prefix_002, ... (or the _1 / _01 spellings;
// the first card fixes the spelling) in key-number order,
// using each card's string value exactly as written between the quotes.
// The header layer trims trailing blanks from string values (FITS
// §4.2.1.1 says they are insignificant), but in a WAT chunk they separate
// tokens, so the value is re-read from the card's raw 80 bytes when
// available. Cards built in memory (no raw bytes) contribute their Value
// unchanged.
//
// It is an error, never a silent drop, when: the numbering has a gap (a
// WATn_k card with WATn_(k-1) missing; WCSTools stops at the gap and
// ignores the rest), a key appears more than once (WCSTools ksearch
// takes the first card and ignores the others), the cards mix spellings (WAT1_001 and WAT1_2,
// or WAT1_1 beside WAT1_001), a key
// number exceeds 499, a WAT card is continued with CONTINUE (WCSTools
// reads only the first card), or the joined string is longer than
// tnxMaxWAT.
func concatWAT(h *header.Header, prefix string) (string, error) {
	// Every card that names a WAT key of this axis, by its exact key.
	// A key whose suffix is not all digits (WAT1_A) is not a WAT card:
	// hgetm looks up numbered keys only.
	present := map[string]bool{}
	count := map[string]int{}
	for _, c := range h.Cards() {
		rest, ok := strings.CutPrefix(c.Key, prefix+"_")
		if !ok || rest == "" || strings.Trim(rest, "0123456789") != "" {
			continue
		}
		present[c.Key] = true
		count[c.Key]++
	}
	if len(present) == 0 {
		return "", nil
	}
	for _, k := range sortedKeys(present) {
		if count[k] > 1 {
			return "", fmt.Errorf("wcs/transform: %s appears %d times; WCSTools reads only the first, so the others would be lost; not supported", k, count[k])
		}
	}
	// The spelling is fixed by the first card. hgetm tries _1, _01, _001
	// in that order, but a header with more than one spelling of the
	// first card is refused here, so no input depends on that order.
	format := ""
	var firsts []string
	for _, f := range []string{"%s_%d", "%s_%02d", "%s_%03d"} {
		if k := fmt.Sprintf(f, prefix, 1); present[k] {
			format = f
			firsts = append(firsts, k)
		}
	}
	switch {
	case len(firsts) > 1:
		return "", fmt.Errorf("wcs/transform: %s has its first card in more than one spelling (%s); not supported", prefix, strings.Join(firsts, ", "))
	case format == "":
		keys := sortedKeys(present)
		return "", fmt.Errorf("wcs/transform: %s has no first card (%s_1/_01/_001) but has %s", prefix, prefix, strings.Join(keys, ", "))
	}
	var sb strings.Builder
	used := 0
	for i := 1; ; i++ {
		key := fmt.Sprintf(format, prefix, i)
		if !present[key] {
			break
		}
		if i > tnxMaxWATCards {
			return "", fmt.Errorf("wcs/transform: %s: WAT key numbers above %d are not read by WCSTools; not supported", key, tnxMaxWATCards)
		}
		c, _ := h.Get(key)
		if len(c.RawCards()) > 1 {
			return "", fmt.Errorf("wcs/transform: %s is continued with CONTINUE cards; not supported in a WAT string", key)
		}
		v, err := watCardString(c)
		if err != nil {
			return "", fmt.Errorf("wcs/transform: %s: %w", key, err)
		}
		sb.WriteString(v)
		used++
	}
	if used != len(present) {
		var stray []string
		for _, k := range sortedKeys(present) {
			n := 0
			fmt.Sscanf(k[len(prefix)+1:], "%d", &n)
			if n < 1 || n > used || fmt.Sprintf(format, prefix, n) != k {
				stray = append(stray, k)
			}
		}
		return "", fmt.Errorf("wcs/transform: %s..%s are followed by %s after a gap or in a different spelling (WCSTools would drop them); not supported",
			fmt.Sprintf(format, prefix, 1), fmt.Sprintf(format, prefix, used), strings.Join(stray, ", "))
	}
	if sb.Len() > tnxMaxWAT {
		return "", fmt.Errorf("wcs/transform: joined %s string is %d characters; WCSTools (the TNX truth source) reads at most %d, longer strings are not supported",
			prefix, sb.Len(), tnxMaxWAT)
	}
	return sb.String(), nil
}

func sortedKeys(m map[string]bool) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// watCardString returns a WAT card's string value without trimming.
func watCardString(c header.Card) (string, error) {
	if c.Type != header.TypeString {
		return "", fmt.Errorf("not a string card")
	}
	if c.Raw == ([header.CardWidth]byte{}) {
		s, _ := c.Value.(string)
		return s, nil
	}
	raw := string(c.Raw[:])
	if len(raw) < 10 || raw[8:10] != "= " {
		return "", fmt.Errorf("no value indicator in raw card")
	}
	body := raw[10:]
	q := strings.IndexByte(body, '\'')
	if q < 0 {
		return "", fmt.Errorf("no opening quote in raw card")
	}
	var sb strings.Builder
	for j := q + 1; j < len(body); j++ {
		if body[j] == '\'' {
			if j+1 < len(body) && body[j+1] == '\'' {
				sb.WriteByte('\'')
				j++
				continue
			}
			return sb.String(), nil
		}
		sb.WriteByte(body[j])
	}
	return "", fmt.Errorf("unterminated string in raw card")
}

// parseWATAttrs splits a joined WAT string into its key = value
// attributes. Values are either a double-quoted string (may contain
// blanks) or a single blank-free token.
func parseWATAttrs(s string) (map[string]string, error) {
	attrs := map[string]string{}
	i := 0
	skip := func() {
		for i < len(s) && s[i] == ' ' {
			i++
		}
	}
	for {
		skip()
		if i >= len(s) {
			return attrs, nil
		}
		k0 := i
		for i < len(s) && s[i] != ' ' && s[i] != '=' {
			i++
		}
		key := s[k0:i]
		skip()
		if i >= len(s) || s[i] != '=' {
			return nil, fmt.Errorf("attribute %q has no '='", key)
		}
		i++
		skip()
		var val string
		if i < len(s) && s[i] == '"' {
			end := strings.IndexByte(s[i+1:], '"')
			if end < 0 {
				return nil, fmt.Errorf("attribute %q: no closing quote", key)
			}
			val = s[i+1 : i+1+end]
			i += end + 2
		} else {
			v0 := i
			for i < len(s) && s[i] != ' ' {
				i++
			}
			val = s[v0:i]
		}
		if _, dup := attrs[key]; dup {
			return nil, fmt.Errorf("attribute %q given twice", key)
		}
		attrs[key] = val
	}
}

// tnxNCoeff is the number of coefficients a surface must carry (WCSTools
// wf_gsrestore).
func tnxNCoeff(xo, yo, cross int) int {
	switch cross {
	case 0:
		return xo + yo - 1
	case 2:
		o := min(xo, yo)
		return xo*yo - o*(o-1)/2
	}
	return xo * yo
}

// parseTNXSurface parses the quoted body of an lngcor / latcor attribute.
func parseTNXSurface(body string) (tnxSurface, error) {
	var s tnxSurface
	fields := strings.Fields(body)
	if len(fields) < 8 {
		return s, fmt.Errorf("need at least 8 numbers, got %d", len(fields))
	}
	nums := make([]float64, len(fields))
	for i, f := range fields {
		v, err := strconv.ParseFloat(f, 64)
		if err != nil || math.IsNaN(v) || math.IsInf(v, 0) {
			return s, fmt.Errorf("bad number %q", f)
		}
		nums[i] = v
	}
	ints := [4]int{}
	for i := range 4 {
		v := nums[i]
		if v != math.Trunc(v) || math.Abs(v) > 1000 {
			return s, fmt.Errorf("field %d (%q) must be an integer", i+1, fields[i])
		}
		ints[i] = int(v)
	}
	s.Type, s.XOrder, s.YOrder, s.Cross = ints[0], ints[1], ints[2], ints[3]
	s.XMin, s.XMax, s.YMin, s.YMax = nums[4], nums[5], nums[6], nums[7]
	switch {
	case s.Type < 1 || s.Type > 3:
		return s, fmt.Errorf("unknown surface type %d (1 Chebyshev, 2 Legendre, 3 polynomial)", s.Type)
	case s.XOrder < 1 || s.YOrder < 1:
		return s, fmt.Errorf("orders must be >= 1, got %d, %d", s.XOrder, s.YOrder)
	case s.Cross < 0 || s.Cross > 2:
		return s, fmt.Errorf("unknown cross-term code %d (0 none, 1 full, 2 half)", s.Cross)
	case !(s.XMax > s.XMin) || !(s.YMax > s.YMin):
		return s, fmt.Errorf("empty domain [%g, %g] x [%g, %g]", s.XMin, s.XMax, s.YMin, s.YMax)
	}
	if want := tnxNCoeff(s.XOrder, s.YOrder, s.Cross); len(nums)-8 != want {
		return s, fmt.Errorf("type %d orders %dx%d cross %d needs %d coefficients, got %d",
			s.Type, s.XOrder, s.YOrder, s.Cross, want, len(nums)-8)
	}
	s.Coeffs = nums[8:]
	s.present = true
	return s, nil
}

// tnxTermAllowed reports whether the (i, j) = (xi power, eta power) term
// is present given the surface's Cross flag (0 none, 1 full, 2 half:
// i + j < max(XOrder, YOrder)). This is the term set of WCSTools
// wf_gseval's shrinking inner loop.
func (s tnxSurface) tnxTermAllowed(i, j int) bool {
	switch s.Cross {
	case 0:
		return !(i > 0 && j > 0)
	case 2:
		return i+j < max(s.XOrder, s.YOrder)
	}
	return true
}

// normalize maps (x, y) to the basis argument: the [-1, 1] normalised
// domain for Chebyshev / Legendre (computed as WCSTools: (x + k1) * k2,
// k1 = -(max+min)/2, k2 = 2/(max-min)), raw (x, y) for the polynomial.
// It also returns d(xn)/dx and d(yn)/dy.
func (s tnxSurface) normalize(x, y float64) (xn, yn, sx, sy float64) {
	if s.Type == 3 {
		return x, y, 1, 1
	}
	sx = 2 / (s.XMax - s.XMin)
	sy = 2 / (s.YMax - s.YMin)
	return float64((x - (s.XMax+s.XMin)/2) * sx), float64((y - (s.YMax+s.YMin)/2) * sy), sx, sy
}

// eval computes the surface value at (x, y).
func (s tnxSurface) eval(x, y float64) float64 {
	v, _, _ := s.evalWithDeriv(x, y)
	return v
}

// evalWithDeriv computes the surface value and its partial derivatives
// with respect to raw x and y.
func (s tnxSurface) evalWithDeriv(x, y float64) (val, dvdx, dvdy float64) {
	if !s.present {
		return 0, 0, 0
	}
	xn, yn, sx, sy := s.normalize(x, y)
	xb, xbp := tnxBasisWithDeriv(s.Type, xn, s.XOrder)
	yb, ybp := tnxBasisWithDeriv(s.Type, yn, s.YOrder)
	k := 0
	for j := 0; j < s.YOrder; j++ {
		for i := 0; i < s.XOrder; i++ {
			if !s.tnxTermAllowed(i, j) {
				continue
			}
			c := s.Coeffs[k]
			k++
			val += float64(c * xb[i] * yb[j])
			dvdx += float64(c * xbp[i] * yb[j])
			dvdy += float64(c * xb[i] * ybp[j])
		}
	}
	return val, float64(dvdx * sx), float64(dvdy * sy)
}

// tnxBasisWithDeriv returns the basis values [b_0, ..., b_{order-1}] and
// their derivatives for argument xn under the given surface type.
func tnxBasisWithDeriv(surfaceType int, xn float64, order int) (b, bp []float64) {
	b = make([]float64, order)
	bp = make([]float64, order)
	if order == 0 {
		return
	}
	b[0] = 1
	if order > 1 {
		b[1] = xn
		bp[1] = 1
	}
	switch surfaceType {
	case 1:
		// Chebyshev T_n = 2x T_{n-1} - T_{n-2}.
		for n := 2; n < order; n++ {
			b[n] = float64(2*xn*b[n-1]) - b[n-2]
			bp[n] = float64(2*b[n-1]) + float64(2*xn*bp[n-1]) - bp[n-2]
		}
	case 2:
		// Legendre n P_n = (2n-1) x P_{n-1} - (n-1) P_{n-2}.
		for n := 2; n < order; n++ {
			f := float64(n)
			b[n] = (float64((float64(2*f)-1)*xn*b[n-1]) - float64((f-1)*b[n-2])) / f
			bp[n] = (float64((float64(2*f)-1)*b[n-1]) + float64((float64(2*f)-1)*xn*bp[n-1]) - float64((f-1)*bp[n-2])) / f
		}
	case 3:
		for n := 2; n < order; n++ {
			b[n] = float64(xn * b[n-1])
			bp[n] = float64(float64(n) * b[n-1])
		}
	}
	return
}

// Forward applies the TNX correction at intermediate world coordinates
// (xi, eta) in degrees: (xi + lngcor(xi, eta), eta + latcor(xi, eta)).
func (t *TNX) Forward(xi, eta float64) (xiOut, etaOut float64) {
	return xi + t.Lon.eval(xi, eta), eta + t.Lat.eval(xi, eta)
}

// Inverse solves Forward(xi, eta) = (xiPrime, etaPrime) by Newton's method
// with the analytic Jacobian, starting from (xiPrime, etaPrime). ok is
// false (and the result must not be used) if the iteration does not
// converge to 1e-13 deg in 50 steps or the Jacobian is singular.
func (t *TNX) Inverse(xiPrime, etaPrime float64) (xi, eta float64, ok bool) {
	xi, eta = xiPrime, etaPrime
	for range 50 {
		lonVal, dLonDxi, dLonDeta := t.Lon.evalWithDeriv(xi, eta)
		latVal, dLatDxi, dLatDeta := t.Lat.evalWithDeriv(xi, eta)
		rx := xiPrime - (xi + lonVal)
		ry := etaPrime - (eta + latVal)
		if math.Abs(rx) < 1e-13 && math.Abs(ry) < 1e-13 {
			return xi, eta, true
		}
		j11, j12 := 1+dLonDxi, dLonDeta
		j21, j22 := dLatDxi, 1+dLatDeta
		det := float64(j11*j22) - float64(j12*j21)
		if det == 0 || math.IsNaN(det) || math.IsInf(det, 0) {
			return 0, 0, false
		}
		xi += (float64(j22*rx) - float64(j12*ry)) / det
		eta += (float64(-j21*rx) + float64(j11*ry)) / det
	}
	return 0, 0, false
}
