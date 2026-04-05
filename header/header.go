package header

import (
	"errors"
	"fmt"
	"maps"
	"strconv"
	"strings"
)

// Header is an ordered, case-insensitive collection of FITS header cards.
//
// Insertion order is preserved. Lookups are case-insensitive (§4.1.2.1
// states all keywords are uppercase). Duplicate keys are allowed for
// commentary keywords (COMMENT, HISTORY, blank) which have "last wins"
// semantics only in the sense that Set() updates the first occurrence;
// Add() always appends.
//
// A Header is safe to mutate directly; persistence back to disk is the
// caller's responsibility (via *fits.File in edit mode).
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
// value type does not match.
var ErrWrongType = errors.New("fits/header: wrong value type")

// New returns an empty Header.
func New() *Header { return &Header{idx: make(map[string]int)} }

// FromCards returns a Header initialized with the given cards. Card order is
// preserved.
func FromCards(cards []Card) *Header {
	h := &Header{cards: make([]Card, len(cards)), idx: make(map[string]int)}
	copy(h.cards, cards)
	for i := range h.cards {
		k := canonKey(h.cards[i].Key)
		if _, dup := h.idx[k]; !dup && !isCommentaryKey(k) {
			h.idx[k] = i
		}
	}
	return h
}

// canonKey returns the canonical (uppercase) form of keyword name. HIERARCH
// keys are case-sensitive in their path segments per convention, but FITS
// itself is always uppercase.
func canonKey(name string) string { return strings.ToUpper(name) }

func isCommentaryKey(k string) bool {
	return k == KeyComment || k == KeyHistory || k == ""
}

// Len returns the number of cards.
func (h *Header) Len() int { return len(h.cards) }

// Cards returns the underlying card slice in insertion order. The returned
// slice must not be modified in place; use Set/Add/Delete to mutate.
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

// Int returns the integer-typed value of name. Float-typed cards that
// represent exact integers are not coerced.
func (h *Header) Int(name string) (int64, error) {
	c, ok := h.Get(name)
	if !ok {
		return 0, fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	if c.Type != TypeInt {
		return 0, fmt.Errorf("%w: %s is %v", ErrWrongType, name, c.Type)
	}
	return c.Value.(int64), nil
}

// Float returns the floating-point value of name. Integer-typed cards are
// promoted to float64 here (no loss for values that fit in a float64 mantissa
// exactly).
func (h *Header) Float(name string) (float64, error) {
	c, ok := h.Get(name)
	if !ok {
		return 0, fmt.Errorf("%w: %s", ErrKeyNotFound, name)
	}
	switch c.Type {
	case TypeFloat:
		return c.Value.(float64), nil
	case TypeInt:
		return float64(c.Value.(int64)), nil
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

// Set updates the value (and optionally comment) of name. If the keyword
// does not exist, Set appends a new card — equivalent to Add in that case.
// For commentary keywords Set is equivalent to Add.
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
		// Preserve Raw on update only when the raw bytes encode the same
		// value as the new card; otherwise zero Raw so Encode re-emits from
		// the typed value.
		c.Raw = [CardWidth]byte{}
		h.cards[i] = c
		return nil
	}
	h.cards = append(h.cards, c)
	h.idx[k] = len(h.cards) - 1
	return nil
}

// Add appends a new card. For commentary keys this is the normal path; for
// other keys it permits duplicates (use with care — only HISTORY/COMMENT are
// expected to appear more than once).
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

// Clone returns a deep copy of h that is safe to mutate independently.
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

// NAXISn returns NAXISi for 1-based axis index. NAXIS1..NAXIS999 are supported.
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

// buildCard constructs a Card from an arbitrary Go value. Supported value
// types: string, bool, int, int8..int64, uint, uint8..uint64, float32,
// float64, Complex, and nil (for commentary cards whose value is text).
func buildCard(key string, value any, comment string) (Card, error) {
	c := Card{Key: key, Comment: comment}
	if isCommentaryKey(key) {
		// Commentary cards carry text in Comment; value is ignored except if
		// it's a string in which case we adopt it as the comment text.
		c.Type = TypeEmpty
		if s, ok := value.(string); ok && comment == "" {
			c.Comment = s
		}
		return c, nil
	}
	switch v := value.(type) {
	case nil:
		c.Type = TypeEmpty
	case string:
		c.Type = TypeString
		c.Value = v
	case bool:
		c.Type = TypeLogical
		c.Value = v
	case int:
		c.Type = TypeInt
		c.Value = int64(v)
	case int8:
		c.Type = TypeInt
		c.Value = int64(v)
	case int16:
		c.Type = TypeInt
		c.Value = int64(v)
	case int32:
		c.Type = TypeInt
		c.Value = int64(v)
	case int64:
		c.Type = TypeInt
		c.Value = v
	case uint:
		c.Type = TypeInt
		c.Value = int64(v)
	case uint8:
		c.Type = TypeInt
		c.Value = int64(v)
	case uint16:
		c.Type = TypeInt
		c.Value = int64(v)
	case uint32:
		c.Type = TypeInt
		c.Value = int64(v)
	case uint64:
		c.Type = TypeInt
		c.Value = int64(v)
	case float32:
		c.Type = TypeFloat
		c.Value = float64(v)
	case float64:
		c.Type = TypeFloat
		c.Value = v
	case Complex:
		c.Type = TypeComplexFloat
		c.Value = v
	default:
		return c, fmt.Errorf("fits/header: unsupported value type %T for key %q", value, key)
	}
	return c, nil
}
