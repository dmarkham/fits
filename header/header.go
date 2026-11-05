package header

import (
	"errors"
	"fmt"
	"maps"
	"math"
	"math/big"
	"strconv"
	"strings"
)

// Header is an ordered collection of FITS header cards with
// case-insensitive keyword lookup.
//
// Insertion order is preserved. Lookups upper-case the name (FITS keywords
// are upper case, §4.1.2.1); a HIERARCH card is looked up by its full key,
// "HIERARCH ESO DET NAME". Get, Set and Delete act on the first card with
// a keyword; Add always appends, so commentary keywords (COMMENT, HISTORY,
// blank) and any duplicate keyword can repeat. Commentary cards are found
// by History, Comments and Delete, not by Get.
//
// Mutations through Set, Add and Delete stay in memory; writing them to
// disk is the caller's job (for example (*fits.File).Flush in edit mode).
type Header struct {
	cards []Card
	// idx maps canonical (uppercase) key to the index of the first matching
	// card in cards. Commentary keys are not tracked here (they are looked up
	// by scanning cards directly).
	idx map[string]int
}

// ErrKeyNotFound is returned by typed getters when a keyword is absent.
var ErrKeyNotFound = errors.New("fits/header: key not found")

// ErrWrongType is returned when a typed getter is called on a card whose
// value type does not match (including a TypeUnparsed card).
var ErrWrongType = errors.New("fits/header: wrong value type")

// ErrIntOverflow is returned by Int for an integer outside the int64
// range (a TypeBigInt card; cfitsio's integer read reports the same card
// as a conversion overflow). BigInt and Float read such cards.
var ErrIntOverflow = errors.New("fits/header: integer overflows int64")

// New returns an empty Header.
func New() *Header { return &Header{idx: make(map[string]int)} }

// FromCards returns a Header initialized with the given cards. Card order is
// preserved, and for a keyword that appears more than once lookups use the
// first card (as cfitsio's keyword search from the top of the header
// does). Cards are not validated here; Encode validates every card it has
// to write from its fields (a card built in memory, or a decoded card
// whose fields were edited, see Cards). Such a card is stored with its
// keyword upper-cased (and an edited copy without its old Raw), as AddCard
// stores it and as it is written, so Cards and Get report the keyword the
// file will carry. A decoded card that still matches its bytes keeps its
// decoded name, lower-case letters included.
func FromCards(cards []Card) *Header {
	h := &Header{cards: make([]Card, len(cards)), idx: make(map[string]int)}
	copy(h.cards, cards)
	for i := range h.cards {
		c := &h.cards[i]
		k := canonKey(c.Key)
		if k != c.Key && (c.Raw == ([CardWidth]byte{}) || !rawMatches(c)) {
			// Written from its fields under the upper-cased keyword:
			// store it so, as AddCard does, so that Key says what the
			// file will carry.
			c.Key, c.Raw = k, [CardWidth]byte{}
		}
		if _, dup := h.idx[k]; !dup && !isCommentaryKey(k) {
			h.idx[k] = i
		}
	}
	return h
}

// canonKey returns the canonical (uppercase) form of keyword name. HIERARCH
// paths are upper-cased too, so HIERARCH lookups are case-insensitive.
func canonKey(name string) string { return strings.ToUpper(name) }

func isCommentaryKey(k string) bool {
	return k == KeyComment || k == KeyHistory || k == ""
}

// Len returns the number of cards.
func (h *Header) Len() int { return len(h.cards) }

// Cards returns the underlying card slice in insertion order. The returned
// slice must not be modified in place; use Set/Add/Delete to mutate.
//
// To build a changed header from it, copy the slice, edit the copies and
// pass them to FromCards or AddCard. A copied card keeps its Raw bytes,
// but Encode writes those bytes only while the card's Key, Type, Value and
// Comment are still what they decode to: an edited copy is written from
// its fields (validated as Set validates them, its keyword upper-cased as
// Set and cfitsio write it, so an edited card that Set could not build is
// an error), never as its old bytes. The bytes count as edited for any
// difference: a change of key letter case, of the sign of a float zero,
// of either part of a complex, of a complex between integer and float.
// FromCards and AddCard store such a copy under its upper-cased keyword,
// so a writer that computes a structural keyword itself (NAXIS1, BITPIX,
// TFIELDS, ...) leaves an edited "naxis1" copy out as it leaves out an
// upper-case NAXIS1.
func (h *Header) Cards() []Card { return h.cards }

// Has reports whether a keyword is present (case-insensitive).
func (h *Header) Has(name string) bool {
	_, ok := h.idx[canonKey(name)]
	return ok
}

// Get returns the first card matching name, or ok=false.
func (h *Header) Get(name string) (Card, bool) {
	i, ok := h.idx[canonKey(name)]
	if !ok {
		return Card{}, false
	}
	return h.cards[i], true
}

// String returns the string-typed value of name.
func (h *Header) String(name string) (string, error) {
	c, ok := h.Get(name)
	if !ok {
		return "", fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	if c.Type != TypeString {
		return "", fmt.Errorf("%w: %s is %v", ErrWrongType, name, c.Type)
	}
	return c.Value.(string), nil
}

// Int returns the integer-typed value of name. A float-typed card is not
// converted (ErrWrongType); an integer outside int64 returns ErrIntOverflow.
func (h *Header) Int(name string) (int64, error) {
	c, ok := h.Get(name)
	if !ok {
		return 0, fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	switch c.Type {
	case TypeInt:
		return c.Value.(int64), nil
	case TypeBigInt:
		return 0, fmt.Errorf("%w: %s = %s", ErrIntOverflow, name, c.Value.(*big.Int))
	}
	return 0, fmt.Errorf("%w: %s is %v", ErrWrongType, name, c.Type)
}

// BigInt returns an integer-typed value of name exactly, whether or not it
// fits in int64 (TypeInt or TypeBigInt). The result is a new *big.Int.
func (h *Header) BigInt(name string) (*big.Int, error) {
	c, ok := h.Get(name)
	if !ok {
		return nil, fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	switch c.Type {
	case TypeInt:
		return big.NewInt(c.Value.(int64)), nil
	case TypeBigInt:
		return new(big.Int).Set(c.Value.(*big.Int)), nil
	}
	return nil, fmt.Errorf("%w: %s is %v", ErrWrongType, name, c.Type)
}

// Float returns the floating-point value of name. Integer-typed cards
// (TypeInt and TypeBigInt) are converted to the nearest float64, which is
// what cfitsio's fits_read_key_dbl returns for them; this is exact for
// integers up to 2^53 in magnitude. An integer card written "-0" gives -0,
// as cfitsio's strtod-based read does, while its Value is int64(0) (a
// copy of it edited to another value gives that value).
func (h *Header) Float(name string) (float64, error) {
	c, ok := h.Get(name)
	if !ok {
		return 0, fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	switch c.Type {
	case TypeFloat:
		return c.Value.(float64), nil
	case TypeInt:
		v := c.Value.(int64)
		if v == 0 && c.negZero {
			return math.Copysign(0, -1), nil
		}
		return float64(v), nil
	case TypeBigInt:
		f, _ := new(big.Float).SetInt(c.Value.(*big.Int)).Float64()
		return f, nil
	}
	return 0, fmt.Errorf("%w: %s is %v", ErrWrongType, name, c.Type)
}

// Bool returns the logical value of name.
func (h *Header) Bool(name string) (bool, error) {
	c, ok := h.Get(name)
	if !ok {
		return false, fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	if c.Type != TypeLogical {
		return false, fmt.Errorf("%w: %s is %v", ErrWrongType, name, c.Type)
	}
	return c.Value.(bool), nil
}

// Set replaces the value and comment of the first card named name, or
// appends a new card when there is none. For commentary keywords Set is
// equivalent to Add. The replaced card is re-encoded from its new value
// (its original bytes are dropped); other cards keep theirs.
//
// The keyword, value and comment are validated: the keyword must be a
// legal FITS keyword (at most 8 characters from A–Z 0–9 '-' '_' after
// upper-casing, or "HIERARCH " plus a path), floats must be finite, and
// text must be printable ASCII. An error leaves the header unchanged.
// Supported value types are listed at buildCard: Go integers (uint64 above
// MaxInt64 is kept exactly as TypeBigInt), *big.Int, floats, bool, string,
// Complex, RawValue and nil (undefined value).
func (h *Header) Set(name string, value any, comment string) error {
	k := canonKey(name)
	if isCommentaryKey(k) {
		return h.Add(name, value, comment)
	}
	c, err := buildCard(k, value, comment)
	if err != nil {
		return err
	}
	if i, ok := h.idx[k]; ok {
		h.cards[i] = c
		return nil
	}
	h.cards = append(h.cards, c)
	h.idx[k] = len(h.cards) - 1
	return nil
}

// Add appends a new card. For commentary keys this is the normal path; for
// other keys it permits duplicates (use with care — only HISTORY/COMMENT are
// expected to appear more than once). Commentary text longer than 72
// characters is kept as one card in memory and written as several
// 72-character cards (see EncodeCardWithContinuation), so it reads back as
// several cards. Validation is as for Set.
func (h *Header) Add(name string, value any, comment string) error {
	k := canonKey(name)
	c, err := buildCard(k, value, comment)
	if err != nil {
		return err
	}
	h.cards = append(h.cards, c)
	if !isCommentaryKey(k) {
		if _, dup := h.idx[k]; !dup {
			h.idx[k] = len(h.cards) - 1
		}
	}
	return nil
}

// AddCard appends c. A card decoded from a file (Raw non-zero, as
// returned by ParseCards or Header.Cards of a header read from disk) is
// kept as it is and written as the bytes it was read from, including its
// CONTINUE cards; this is how a header is copied into a new HDU, and it
// works for cards Set could not build (a keyword cfitsio reads as
// "DATE.OBS" or "KEYLONG10", a malformed value, a non-ASCII comment byte).
// A copy of a decoded card whose Key, Type, Value or Comment was changed
// no longer matches its bytes: it is treated as a card built in memory
// (its Raw dropped) and written from its fields, so the writers write
// what the header's getters return. A card built in memory
// (Raw zero) is validated as Set validates it, its keyword upper-cased,
// and an error leaves the header unchanged.
//
// A decoded card that still matches its bytes but whose keyword is not
// upper case ("bzero", "naxis1", "object"; a HIERARCH path is exempt) is
// refused with ErrInvalidKeyword. Such a card is illegal FITS (§4.1.2.1)
// and readers disagree on it: cfitsio's HDU parser does not apply a
// "bzero" card to the pixels while astropy does, and astropy takes
// "naxis1" for NAXIS1. No copy of its bytes reads the same everywhere,
// and a writer filtering structural keywords would not recognise it. To
// copy it, set the copy's Key to upper case (it is then written from its
// fields, under that keyword) or leave it out. A lower-case HIERARCH path
// is copied as its bytes: fitsverify accepts it and both readers find it
// by any case.
func (h *Header) AddCard(c Card) error {
	if c.Raw != ([CardWidth]byte{}) && !rawMatches(&c) {
		c.Raw = [CardWidth]byte{}
	}
	if c.Raw == ([CardWidth]byte{}) {
		c.Key = canonKey(c.Key)
		if c.Key == KeyEnd {
			return fmt.Errorf("%w: END is written by Encode, not set", ErrInvalidKeyword)
		}
		if err := validateCard(&c); err != nil {
			return err
		}
	} else if k := canonKey(c.Key); k != c.Key && !c.IsHierarch() {
		return fmt.Errorf("%w: %q is not upper case (FITS 4.0 §4.1.2.1), and readers disagree on such a card "+
			"(cfitsio's HDU parser matches BSCALE, BZERO, BLANK and the structural keywords case-sensitively, "+
			"astropy upper-cases every keyword), so its bytes are not copied: set the copy's Key to %q to have it "+
			"written from its fields, or leave the card out", ErrInvalidKeyword, c.Key, k)
	}
	h.cards = append(h.cards, c)
	if k := canonKey(c.Key); !isCommentaryKey(k) {
		if _, dup := h.idx[k]; !dup {
			h.idx[k] = len(h.cards) - 1
		}
	}
	return nil
}

// Delete removes the first card matching name. Returns ErrKeyNotFound if
// the card is absent.
func (h *Header) Delete(name string) error {
	k := canonKey(name)
	i, ok := h.idx[k]
	if !ok {
		// Commentary cards: delete the first occurrence by linear scan.
		if isCommentaryKey(k) {
			for j := range h.cards {
				if canonKey(h.cards[j].Key) == k {
					h.cards = append(h.cards[:j], h.cards[j+1:]...)
					h.rebuildIdx()
					return nil
				}
			}
		}
		return fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	h.cards = append(h.cards[:i], h.cards[i+1:]...)
	h.rebuildIdx()
	return nil
}

// Clone returns a copy of h whose cards can be changed with Set, Add and
// Delete without affecting h. Card values are shared, not copied: do not
// modify a *big.Int taken from a card in place.
func (h *Header) Clone() *Header {
	out := &Header{
		cards: make([]Card, len(h.cards)),
		idx:   make(map[string]int, len(h.idx)),
	}
	copy(out.cards, h.cards)
	maps.Copy(out.idx, h.idx)
	return out
}

// History returns the text of every HISTORY card in insertion order.
func (h *Header) History() []string {
	var out []string
	for i := range h.cards {
		if h.cards[i].Key == KeyHistory {
			out = append(out, h.cards[i].Comment)
		}
	}
	return out
}

// Comments returns the text of every COMMENT card in insertion order.
func (h *Header) Comments() []string {
	var out []string
	for i := range h.cards {
		if h.cards[i].Key == KeyComment {
			out = append(out, h.cards[i].Comment)
		}
	}
	return out
}

// NAXIS returns NAXIS (number of axes) — convenience.
func (h *Header) NAXIS() (int, error) {
	v, err := h.Int(KeyNaxis)
	return int(v), err
}

// NAXISn returns the integer value of the keyword "NAXIS" + i (i is the
// 1-based axis number).
func (h *Header) NAXISn(i int) (int64, error) {
	return h.Int(naxisKey(i))
}

func naxisKey(i int) string {
	return KeyNaxis + strconv.Itoa(i)
}

func (h *Header) rebuildIdx() {
	h.idx = make(map[string]int, len(h.cards))
	for i := range h.cards {
		k := canonKey(h.cards[i].Key)
		if isCommentaryKey(k) {
			continue
		}
		if _, dup := h.idx[k]; !dup {
			h.idx[k] = i
		}
	}
}

// buildCard constructs a validated Card from a Go value. Supported value
// types: string, bool, int, int8..int64, uint, uint8..uint64 (values above
// MaxInt64 become TypeBigInt and are written as exact decimal integers, as
// cfitsio writes uint64 keywords), *big.Int, float32, float64, Complex,
// RawValue (written verbatim: it carries a value read from a malformed
// card; its Text must be one token that reads back as the same unparsed
// text, so no blanks, no '/', no leading quote, no valid FITS value), and
// nil (an undefined value). A commentary key takes its text from the
// comment or from a string value, not both.
//
// It refuses, with an error: keywords that are not legal FITS keywords
// (more than 8 characters without "HIERARCH ", characters outside A–Z 0–9
// '-' '_'; the name is upper-cased first; HIERARCH paths with other than
// single blanks between words), END and CONTINUE, NaN and ±Inf, a RawValue
// that would not read back, and text outside printable ASCII.
func buildCard(key string, value any, comment string) (Card, error) {
	c := Card{Key: key, Comment: comment}
	if key == KeyEnd {
		return c, fmt.Errorf("%w: END is written by Encode, not set", ErrInvalidKeyword)
	}
	if key == "CONTINUE" {
		return c, fmt.Errorf("%w: CONTINUE cards are written by the encoder for long strings; set the string on its own keyword", ErrInvalidKeyword)
	}
	if isCommentaryKey(key) {
		// Commentary cards carry one text in Comment; a string value is
		// taken as the text when no comment is given. Both at once is an
		// error (one of them would be dropped).
		c.Type = TypeEmpty
		switch s := value.(type) {
		case nil:
		case string:
			if comment != "" && s != "" {
				return c, fmt.Errorf("fits/header: commentary key %q has one text; got both a string value %q and a comment %q", key, s, comment)
			}
			if comment == "" {
				c.Comment = s
			}
		default:
			return c, fmt.Errorf("fits/header: commentary key %q takes text, not %T", key, value)
		}
		return c, validateCard(&c)
	}
	switch v := value.(type) {
	case nil:
		c.Type = TypeEmpty
	case string:
		c.Type, c.Value = TypeString, v
	case bool:
		c.Type, c.Value = TypeLogical, v
	case int:
		c.Type, c.Value = TypeInt, int64(v)
	case int8:
		c.Type, c.Value = TypeInt, int64(v)
	case int16:
		c.Type, c.Value = TypeInt, int64(v)
	case int32:
		c.Type, c.Value = TypeInt, int64(v)
	case int64:
		c.Type, c.Value = TypeInt, v
	case uint:
		setUint(&c, uint64(v))
	case uint8:
		c.Type, c.Value = TypeInt, int64(v)
	case uint16:
		c.Type, c.Value = TypeInt, int64(v)
	case uint32:
		c.Type, c.Value = TypeInt, int64(v)
	case uint64:
		setUint(&c, v)
	case *big.Int:
		if v == nil {
			return c, fmt.Errorf("fits/header: nil *big.Int for key %q", key)
		}
		if v.IsInt64() {
			c.Type, c.Value = TypeInt, v.Int64()
		} else {
			c.Type, c.Value = TypeBigInt, new(big.Int).Set(v)
		}
	case float32:
		c.Type, c.Value = TypeFloat, float64(v)
	case float64:
		c.Type, c.Value = TypeFloat, v
	case Complex:
		c.Type, c.Value = TypeComplexFloat, v
	case RawValue:
		c.Type, c.Value = TypeUnparsed, v
	default:
		return c, fmt.Errorf("fits/header: unsupported value type %T for key %q", value, key)
	}
	return c, validateCard(&c)
}

func setUint(c *Card, v uint64) {
	if v <= math.MaxInt64 {
		c.Type, c.Value = TypeInt, int64(v)
		return
	}
	c.Type, c.Value = TypeBigInt, new(big.Int).SetUint64(v)
}
