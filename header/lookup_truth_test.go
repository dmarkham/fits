package header

import (
	"bytes"
	"encoding/json"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// The Header doc claims case-insensitive keyword lookup for Get, Has, the
// typed getters, Set and Delete, and that each acts on the FIRST card with
// a keyword, as cfitsio's fits_read_keyword / fits_delete_key do (they
// upper-case the name and search from the top of the header). The truth
// is libcfitsio 4.7.0 itself, through the hdrref harness ("hdrref
// readkey" / "hdrref delkey", header/testdata/hdrref/hdrref.c), run on the
// header this package encodes.

type readkeyResult struct {
	Key     string  `json:"key"`
	Status  int     `json:"status"`
	Count   int     `json:"count"`
	Value   *string `json:"value"`
	Comment *string `json:"comment"`
}

func hdrrefReadkey(t *testing.T, path string, names ...string) map[string]readkeyResult {
	t.Helper()
	var res []readkeyResult
	if err := json.Unmarshal(runHdrref(t, append([]string{"readkey", path}, names...)...), &res); err != nil {
		t.Fatalf("decode hdrref readkey: %v", err)
	}
	out := map[string]readkeyResult{}
	for _, r := range res {
		out[r.Key] = r
	}
	return out
}

// lookupHeader is built ONLY with lower- and mixed-case names, so every
// call below exercises the case folding.
func lookupHeader(t *testing.T) *Header {
	t.Helper()
	h := New()
	for _, c := range []struct {
		k string
		v any
	}{{"simple", true}, {"bitpix", int64(8)}, {"naxis", int64(0)}} {
		if err := h.Set(c.k, c.v, ""); err != nil {
			t.Fatal(err)
		}
	}
	steps := []func() error{
		func() error { return h.Set("object", "M31", "first") },
		// Set on an existing keyword in another case replaces that card;
		// it must not append a second OBJECT.
		func() error { return h.Set("Object", "M32", "replaced") },
		// A new, non-commentary key through Add, then a duplicate of it:
		// lookups must keep finding the FIRST card.
		func() error { return h.Add("dupkey", int64(1), "first dup") },
		func() error { return h.Add("DupKey", int64(2), "second dup") },
		func() error { return h.Add("foo", int64(7), "") },
		func() error { return h.Set("flag", true, "") },
		func() error { return h.Set("ratio", 0.25, "") },
	}
	for i, s := range steps {
		if err := s(); err != nil {
			t.Fatalf("step %d: %v", i, err)
		}
	}
	return h
}

func writeHeader(t *testing.T, h *Header, name string) string {
	t.Helper()
	b, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func TestLookupCaseAndFirstCardMatchCfitsio(t *testing.T) {
	h := lookupHeader(t)

	// In memory, independent of cfitsio: one OBJECT card, two DUPKEY.
	n := map[string]int{}
	for _, c := range h.Cards() {
		n[c.Key]++
	}
	if n["OBJECT"] != 1 || n["DUPKEY"] != 2 {
		t.Fatalf("card counts OBJECT %d DUPKEY %d, want 1 and 2", n["OBJECT"], n["DUPKEY"])
	}
	for _, k := range []string{"object", "OBJECT", "dupkey", "Foo", "flag", "RATIO"} {
		if !h.Has(k) {
			t.Errorf("Has(%q) = false", k)
		}
	}
	if h.Has("nosuch") {
		t.Error(`Has("nosuch") = true`)
	}
	if v, err := h.String("oBjEcT"); err != nil || v != "M32" {
		t.Errorf(`String("oBjEcT") = %q, %v; want M32`, v, err)
	}
	if v, err := h.Int("dupkey"); err != nil || v != 1 {
		t.Errorf(`Int("dupkey") = %d, %v; want the first card's 1`, v, err)
	}
	if v, err := h.Int("foo"); err != nil || v != 7 {
		t.Errorf(`Int("foo") = %d, %v; want 7`, v, err)
	}
	if v, err := h.Bool("FLAG"); err != nil || !v {
		t.Errorf(`Bool("FLAG") = %v, %v`, v, err)
	}
	if v, err := h.Float("ratio"); err != nil || v != 0.25 {
		t.Errorf(`Float("ratio") = %v, %v`, v, err)
	}
	if c, ok := h.Get("Foo"); !ok || c.Key != "FOO" {
		t.Errorf(`Get("Foo") = %+v, %v`, c, ok)
	}

	// cfitsio reads the encoded header by the same lower-case names.
	path := writeHeader(t, h, "lookup.fits")
	cf := hdrrefReadkey(t, path, "object", "dupkey", "foo", "flag", "ratio", "nosuch")
	for _, k := range []string{"object", "dupkey", "foo", "flag", "ratio"} {
		r := cf[k]
		if r.Status != 0 || r.Value == nil {
			t.Fatalf("cfitsio fits_read_keyword(%q): status %d", k, r.Status)
		}
		c, ok := h.Get(k)
		if !ok {
			t.Fatalf("Get(%q) not found; cfitsio finds %q", k, *r.Value)
		}
		enc, err := EncodeCard(c)
		if err != nil {
			t.Fatal(err)
		}
		// The value field of the card Get returns is the one cfitsio's
		// first-match read returns.
		var vc, cm string
		if vc, cm, err = valueAndComment(enc); err != nil {
			t.Fatal(err)
		}
		if vc != *r.Value || cm != *r.Comment {
			t.Errorf("%s: Get gives value %q comment %q; cfitsio %q / %q", k, vc, cm, *r.Value, *r.Comment)
		}
	}
	if cf["object"].Count != 1 || cf["dupkey"].Count != 2 {
		t.Errorf("cfitsio sees OBJECT x%d, DUPKEY x%d; want 1 and 2", cf["object"].Count, cf["dupkey"].Count)
	}
	if cf["nosuch"].Status != 202 {
		t.Errorf("cfitsio status for a missing key %d, want 202 (KEY_NO_EXIST)", cf["nosuch"].Status)
	}
}

// valueAndComment splits an encoded card into the value text and comment
// the way cfitsio's fits_parse_value reports them (value trimmed, comment
// after the '/' with its leading blank dropped). Only used on simple
// one-card values built above.
func valueAndComment(card [CardWidth]byte) (string, string, error) {
	s := string(card[10:])
	val, cm := s, ""
	if i := strings.Index(s, " /"); i >= 0 && !strings.HasPrefix(strings.TrimSpace(s), "'") {
		val, cm = s[:i], s[i+2:]
	} else if strings.HasPrefix(strings.TrimSpace(s), "'") {
		// value is a quoted string: find its closing quote.
		t := strings.TrimLeft(s, " ")
		j := 1
		for j < len(t) {
			if t[j] == '\'' {
				if j+1 < len(t) && t[j+1] == '\'' {
					j += 2
					continue
				}
				break
			}
			j++
		}
		val = t[:j+1]
		rest := t[j+1:]
		if k := strings.Index(rest, "/"); k >= 0 {
			cm = rest[k+1:]
		}
	}
	cm = strings.TrimPrefix(cm, " ")
	return strings.TrimSpace(val), strings.TrimRight(cm, " "), nil
}

// Delete with a lower-case name removes the FIRST card of a duplicated
// keyword, as cfitsio's fits_delete_key does: the encoded header after
// Delete("dupkey") is byte-identical to cfitsio deleting "dupkey" from the
// header before.
func TestDeleteFirstCardMatchesCfitsio(t *testing.T) {
	h := lookupHeader(t)
	before := writeHeader(t, h, "before.fits")
	if err := h.Delete("dupkey"); err != nil {
		t.Fatal(err)
	}
	if v, err := h.Int("DUPKEY"); err != nil || v != 2 {
		t.Fatalf("after Delete, Int(DUPKEY) = %d, %v; want the second card's 2", v, err)
	}
	ours, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	cfOut := filepath.Join(t.TempDir(), "cfdel.fits")
	runHdrrefWriting(t, []string{cfOut}, "delkey", before, cfOut, "dupkey")
	cf, err := os.ReadFile(cfOut)
	if err != nil {
		t.Fatal(err)
	}
	if !bytes.Equal(ours, cf) {
		t.Fatalf("header after Delete(\"dupkey\") differs from cfitsio fits_delete_key: %s", firstCardDiff(ours, cf))
	}
	// Deleting the remaining card, then a third time: ErrKeyNotFound.
	if err := h.Delete("DupKey"); err != nil {
		t.Fatal(err)
	}
	if err := h.Delete("dupkey"); !errors.Is(err, ErrKeyNotFound) {
		t.Fatalf("Delete of an absent key: %v, want ErrKeyNotFound", err)
	}
}

// ErrKeyNotFound is what the doc says every typed getter (and Delete)
// returns for an absent keyword, and ErrWrongType for a present one of
// another type; callers branch on errors.Is, so both are pinned.
func TestTypedGettersErrKeyNotFound(t *testing.T) {
	h := lookupHeader(t)
	getters := map[string]func(string) error{
		"String": func(k string) error { _, err := h.String(k); return err },
		"Int":    func(k string) error { _, err := h.Int(k); return err },
		"BigInt": func(k string) error { _, err := h.BigInt(k); return err },
		"Float":  func(k string) error { _, err := h.Float(k); return err },
		"Bool":   func(k string) error { _, err := h.Bool(k); return err },
		"Delete": func(k string) error { return h.Delete(k) },
	}
	for name, get := range getters {
		if err := get("nosuch"); !errors.Is(err, ErrKeyNotFound) || errors.Is(err, ErrWrongType) {
			t.Errorf("%s(absent): %v, want ErrKeyNotFound", name, err)
		}
	}
	wrong := map[string]func() error{
		"String": func() error { _, err := h.String("foo"); return err },
		"Int":    func() error { _, err := h.Int("object"); return err },
		"BigInt": func() error { _, err := h.BigInt("flag"); return err },
		"Float":  func() error { _, err := h.Float("object"); return err },
		"Bool":   func() error { _, err := h.Bool("foo"); return err },
	}
	for name, get := range wrong {
		if err := get(); !errors.Is(err, ErrWrongType) || errors.Is(err, ErrKeyNotFound) {
			t.Errorf("%s(wrong type): %v, want ErrWrongType", name, err)
		}
	}
}
