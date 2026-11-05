package fits

import (
	"fmt"
	"math"
	"math/big"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
)

// This file holds the cfitsio-compatible card readers used by the
// scaling-card rule (scalingCards, scaling.go) and the checksum code:
//
//   - the keyword name is ffgknm's (cfKeyName): up to the first blank or
//     '=', the name after "HIERARCH " for a HIERARCH card;
//   - the value is ffpsvc's (cfValueString);
//   - cfC2DD is ffc2dd (C strtod after replacing the first uppercase 'D'
//     with 'E'; trailing junk, ERANGE, Inf and NaN fail) and cfC2JJ is
//     ffc2jj (C strtoll base 10). scalingCards accepts a BSCALE / BZERO /
//     BLANK value only when these agree with the strict FITS number
//     grammar, so no accepted card is read differently by cfitsio.
//
// cStrtod / cStrtoll are proven against glibc (TestCStrtodMatchesGlibc)
// and cfC2DD / cfC2JJ against cfitsio's own functions
// (TestC2DDAndC2JJMatchCfitsio).

// nullUndefined is cfitsio's NULL_UNDEFINED (fitsio2.h): cfitsio reads a
// BLANK equal to it as "no BLANK" (astropy does not), so scalingCards
// refuses it.
const nullUndefined = 1234554321

// Slots of cardKeys.
const (
	keyBscale = iota
	keyBzero
	keyBlank
)

// rawCards yields the 80-byte cards of raw header bytes.
func rawCards(hdr []byte) func(yield func([]byte) bool) {
	return func(yield func([]byte) bool) {
		for off := 0; off+header.CardWidth <= len(hdr); off += header.CardWidth {
			if !yield(hdr[off : off+header.CardWidth]) {
				return
			}
		}
	}
}

// encodedCards yields the first physical card of each card of h as it
// would be written (a card read from disk is its stored bytes), skipping
// those for which skip reports true. A card that cannot be encoded never
// reaches the disk; if it is a scaling keyword it is yielded as a card
// whose value is an unterminated string ("BZERO   = '"), which
// cfitsio-style conversion treats as unconvertible. (The value must start
// in column 11: a card shorter than 9 columns has no value for ffpsvc, and
// an empty value converts to 0.)
func encodedCards(h *header.Header, skip func(key string) bool) func(yield func([]byte) bool) {
	return func(yield func([]byte) bool) {
		for _, c := range h.Cards() {
			if skip != nil && skip(c.Key) {
				continue
			}
			enc, err := header.EncodeCardWithContinuation(c)
			if err != nil || len(enc) == 0 {
				if isScalingKey(c.Key) {
					bad := []byte(fmt.Sprintf("%-8s= '", strings.ToUpper(c.Key)))
					if !yield(bad) {
						return
					}
				}
				continue
			}
			if !yield(enc[0][:]) {
				return
			}
		}
	}
}

func isScalingKey(key string) bool {
	for _, n := range cardKeys {
		if strings.EqualFold(key, n) {
			return true
		}
	}
	return false
}

// cfKeyName is cfitsio's ffgknm on a card with trailing blanks removed.
func cfKeyName(card string) string {
	if strings.HasPrefix(card, "HIERARCH ") {
		eq := strings.IndexByte(card, '=')
		if eq < 0 {
			return "HIERARCH"
		}
		// eq >= 9: the "HIERARCH " prefix holds no '='.
		return strings.TrimRight(strings.TrimLeft(card[9:eq], " "), " ")
	}
	const maxName = 74 // FLEN_KEYWORD - 1
	for i := 0; i < len(card) && i < maxName; i++ {
		if card[i] == ' ' || card[i] == '=' {
			return card[:i]
		}
	}
	return card[:min(len(card), maxName)]
}

// cfValueString is cfitsio's ffpsvc value string of a card with trailing
// blanks removed: a quoted string keeps its quotes; ok is false where
// ffpsvc reports an error (a complex value without its closing ')').
func cfValueString(card string) (value string, ok bool) {
	const maxValue = 70 // FLEN_VALUE - 1
	n := len(card)
	var valpos int
	switch {
	case strings.HasPrefix(card, "HIERARCH "):
		eq := strings.IndexByte(card, '=')
		if eq < 0 {
			return "", true
		}
		valpos = eq + 1
	case n < 9 || strings.HasPrefix(card, "COMMENT ") || strings.HasPrefix(card, "HISTORY ") ||
		strings.HasPrefix(card, "END     ") || strings.HasPrefix(card, "CONTINUE") ||
		strings.HasPrefix(card, "        "):
		return "", true
	case n >= 10 && card[8:10] == "= ":
		valpos = 10
	default:
		eq := strings.IndexByte(card, '=')
		if eq < 0 {
			return "", true
		}
		valpos = eq + 1
	}
	i := valpos
	for i < n && card[i] == ' ' {
		i++
	}
	if i == n {
		return "", true // no value: legal, the value is undefined
	}
	switch card[i] {
	case '/':
		return "", true
	case '\'':
		var b strings.Builder
		b.WriteByte('\'')
		jj := 1
		for i++; i < n && jj < maxValue; i, jj = i+1, jj+1 {
			if card[i] == '\'' {
				if i+1 < n && card[i+1] == '\'' {
					b.WriteByte('\'')
					i++
					jj++
				} else {
					b.WriteByte('\'')
					return b.String(), true
				}
			}
			if jj < maxValue {
				b.WriteByte(card[i])
			}
		}
		s := b.String()
		if len(s) > maxValue-1 {
			s = s[:maxValue-1]
		}
		return s + "'", true // no closing quote: cfitsio closes it
	case '(':
		j := strings.IndexByte(card[i:], ')')
		if j < 0 || j >= maxValue {
			return "", false
		}
		return card[i : i+j+1], true
	}
	j := strings.IndexAny(card[i:], " /")
	if j < 0 {
		j = n - i
	}
	if j > maxValue {
		j = maxValue
	}
	return card[i : i+j], true
}

// cfC2DD is cfitsio's ffc2dd: true when the value string converts to a
// finite double with nothing after it.
func cfC2DD(v string) (float64, bool) {
	if strings.IndexByte(v, 'D') >= 0 {
		if len(v) > 72 {
			return 0, false
		}
		v = strings.Replace(v, "D", "E", 1)
	}
	f, n, erange := cStrtod(v)
	if n < len(v) && v[n] != ' ' {
		return 0, false
	}
	if erange || math.IsInf(f, 0) || math.IsNaN(f) {
		return 0, false
	}
	return f, true
}

// cfC2JJ is cfitsio's ffc2jj (strtoll base 10, nothing after the number,
// no overflow).
func cfC2JJ(v string) (int64, bool) {
	x, n, erange := cStrtoll(v)
	if n < len(v) && v[n] != ' ' {
		return 0, false
	}
	return x, !erange
}

// cStrtoll is C strtoll with base 10: optional white space and sign, then
// decimal digits. It returns the value (LLONG_MAX / LLONG_MIN on
// overflow), the bytes consumed (0 when there is no digit) and ERANGE.
func cStrtoll(v string) (x int64, n int, erange bool) {
	i := 0
	for i < len(v) && isCSpace(v[i]) {
		i++
	}
	neg := false
	if i < len(v) && (v[i] == '+' || v[i] == '-') {
		neg = v[i] == '-'
		i++
	}
	d := i
	for i < len(v) && isDigit(v[i]) {
		i++
	}
	if i == d {
		return 0, 0, false
	}
	digits := v[d:i]
	if neg {
		digits = "-" + digits
	}
	x, err := strconv.ParseInt(digits, 10, 64)
	if err != nil { // only ErrRange: x is the clamped value
		return x, i, true
	}
	return x, i, false
}

// cfAtof is C atof (strtod, junk and range errors ignored), as cfitsio's
// ffvcks reads DATASUM.
func cfAtof(v string) float64 {
	f, _, _ := cStrtod(v)
	return f
}

func isCSpace(c byte) bool {
	return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r'
}

func isDigit(c byte) bool { return c >= '0' && c <= '9' }

func isHexDigit(c byte) bool {
	return isDigit(c) || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F')
}

func hasPrefixFold(s, p string) bool { return len(s) >= len(p) && strings.EqualFold(s[:len(p)], p) }

// cStrtod is C99 strtod in the "C" locale as glibc implements it: the
// longest prefix of s that is a decimal or hexadecimal floating constant,
// "inf"/"infinity" or "nan"/"nan(chars)" (case-insensitive), after
// optional white space and sign. It returns the value (correctly rounded),
// the number of bytes consumed (0 when no conversion was performed) and
// whether glibc sets ERANGE: on overflow, and on underflow, i.e. a value
// that is tiny (below 2^-1022 after rounding to 53 bits with an unbounded
// exponent) and not exactly representable.
func cStrtod(s string) (v float64, n int, erange bool) {
	i := 0
	for i < len(s) && isCSpace(s[i]) {
		i++
	}
	neg := false
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		neg = s[i] == '-'
		i++
	}
	sign := func(f float64) float64 {
		if neg {
			return -f
		}
		return f
	}
	r := s[i:]
	switch {
	case hasPrefixFold(r, "infinity"):
		return sign(math.Inf(1)), i + 8, false
	case hasPrefixFold(r, "inf"):
		return sign(math.Inf(1)), i + 3, false
	case hasPrefixFold(r, "nan"):
		m := 3
		if len(r) > 3 && r[3] == '(' {
			j := 4
			for j < len(r) && (isDigit(r[j]) || r[j] == '_' || (r[j]|0x20 >= 'a' && r[j]|0x20 <= 'z')) {
				j++
			}
			if j < len(r) && r[j] == ')' {
				m = j + 1
			}
		}
		return math.NaN(), i + m, false
	}
	if len(r) >= 2 && r[0] == '0' && (r[1] == 'x' || r[1] == 'X') {
		j, digits, nonzero := 2, 0, false
		for ; j < len(r) && isHexDigit(r[j]); j++ {
			digits++
			nonzero = nonzero || r[j] != '0'
		}
		if j < len(r) && r[j] == '.' {
			for j++; j < len(r) && isHexDigit(r[j]); j++ {
				digits++
				nonzero = nonzero || r[j] != '0'
			}
		}
		if digits == 0 { // "0x" with no hex digit: just the "0"
			return sign(0), i + 1, false
		}
		lit := r[:j]
		exp := "p0"
		if j < len(r) && (r[j] == 'p' || r[j] == 'P') {
			k := j + 1
			if k < len(r) && (r[k] == '+' || r[k] == '-') {
				k++
			}
			d := k
			for k < len(r) && isDigit(r[k]) {
				k++
			}
			if k > d {
				exp, j = r[j:k], k
			}
		}
		f, er := parseCFloat(lit+exp, nonzero, true)
		return sign(f), i + j, er
	}
	j, digits, nonzero := 0, 0, false
	for ; j < len(r) && isDigit(r[j]); j++ {
		digits++
		nonzero = nonzero || r[j] != '0'
	}
	if j < len(r) && r[j] == '.' {
		for j++; j < len(r) && isDigit(r[j]); j++ {
			digits++
			nonzero = nonzero || r[j] != '0'
		}
	}
	if digits == 0 {
		return 0, 0, false
	}
	if j < len(r) && (r[j] == 'e' || r[j] == 'E') {
		k := j + 1
		if k < len(r) && (r[k] == '+' || r[k] == '-') {
			k++
		}
		d := k
		for k < len(r) && isDigit(r[k]) {
			k++
		}
		if k > d {
			j = k
		}
	}
	f, er := parseCFloat(r[:j], nonzero, false)
	return sign(f), i + j, er
}

// parseCFloat converts an unsigned decimal or hex ("0x...p...") literal
// and reports glibc's ERANGE.
func parseCFloat(lit string, nonzero, hex bool) (float64, bool) {
	f, err := strconv.ParseFloat(lit, 64)
	if err != nil {
		// The literal is well formed, so the only error is ErrRange
		// (overflow); ParseFloat then returns +Inf.
		return math.Inf(1), true
	}
	if f == 0 {
		return 0, nonzero
	}
	if f > 0x1p-1022 {
		return f, false
	}
	// At or below the smallest normal. glibc reports underflow when the
	// result is tiny and inexact; tiny means the exact value rounded to 53
	// bits with an unbounded exponent is below 2^-1022 (so a value just
	// below 2^-1022 that rounds up to it at 53 bits is not tiny).
	var exact *big.Rat
	if hex {
		bf, _, perr := new(big.Float).SetPrec(2048).Parse(lit, 0)
		if perr != nil {
			return f, true
		}
		exact, _ = bf.Rat(nil)
	} else {
		var ok bool
		if exact, ok = new(big.Rat).SetString(lit); !ok {
			return f, true
		}
	}
	r53 := new(big.Float).SetPrec(53).SetMode(big.ToNearestEven).SetRat(exact)
	tiny := r53.Cmp(big.NewFloat(0x1p-1022)) < 0
	return f, tiny && exact.Cmp(new(big.Rat).SetFloat64(f)) != 0
}

// cfUnquote is cfitsio's ffc2s on a value string: a quoted string without
// its quotes (a doubled quote read as one), trailing blanks dropped;
// anything else as is.
func cfUnquote(v string) string {
	if !strings.HasPrefix(v, "'") {
		return v
	}
	var b strings.Builder
	for i := 1; i < len(v); i++ {
		if v[i] == '\'' {
			if i+1 < len(v) && v[i+1] == '\'' {
				i++
			} else {
				break
			}
		}
		b.WriteByte(v[i])
	}
	return strings.TrimRight(b.String(), " ")
}

// cfFindCard is cfitsio's ffgcrd for a plain keyword name: starting at
// card index from (0-based, among the cards before END), it searches to
// END and then wraps to the top, comparing each card's ffgknm name in
// upper case. It returns the index and ffpsvc value string of the card
// found, and whether one was found.
func cfFindCard(hdr []byte, key string, from int) (pos int, value string, found bool) {
	var names []string
	var cards []string
	rawCards(hdr)(func(raw []byte) bool {
		card := strings.TrimRight(string(raw), " ")
		name := strings.ToUpper(cfKeyName(card))
		if name == header.KeyEnd {
			return false
		}
		names = append(names, name)
		cards = append(cards, card)
		return true
	})
	n := len(names)
	for k := range n {
		i := (from + k) % n
		if names[i] == key {
			v, _ := cfValueString(cards[i])
			return i, v, true
		}
	}
	return -1, "", false
}
