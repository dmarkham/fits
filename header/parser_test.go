package header

import (
	"bytes"
	"testing"
)

// mkBlock builds a 2880-byte header block from an arbitrary count of
// 80-byte card strings; pads with all-space cards.
func mkBlock(cards ...string) []byte {
	const size = 2880
	out := make([]byte, 0, size)
	for _, s := range cards {
		c := make([]byte, CardWidth)
		for i := range c {
			c[i] = ' '
		}
		copy(c, s)
		out = append(out, c...)
	}
	// Pad with blank cards up to block boundary.
	for len(out) < size {
		blank := bytes.Repeat([]byte{' '}, CardWidth)
		out = append(out, blank...)
	}
	return out
}

func TestParseSimple(t *testing.T) {
	buf := mkBlock(
		"SIMPLE  =                    T / standard FITS",
		"BITPIX  =                    8",
		"NAXIS   =                    0",
		"END",
	)
	cards, endIdx, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	if endIdx != 3 {
		t.Fatalf("end idx %d", endIdx)
	}
	if len(cards) != 3 {
		t.Fatalf("card count %d", len(cards))
	}
	if cards[0].Key != "SIMPLE" || cards[0].Value.(bool) != true {
		t.Fatalf("SIMPLE: %+v", cards[0])
	}
	if cards[1].Value.(int64) != 8 {
		t.Fatalf("BITPIX: %+v", cards[1])
	}
}

func TestParseNoEnd(t *testing.T) {
	buf := mkBlock("SIMPLE  =                    T")
	_, _, err := ParseCards(buf)
	if err == nil {
		t.Fatal("expected ErrNoEnd")
	}
}

func TestParseContinue(t *testing.T) {
	// A long string that spans two cards using the CONTINUE convention.
	// Note: two-card payload shows joining of "hello..." + "world!" without
	// the '&'.
	buf := mkBlock(
		"LONGSTR = 'first part, with &'",
		"CONTINUE  'second part.'",
		"END",
	)
	cards, _, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	if len(cards) != 1 {
		t.Fatalf("expected 1 logical card, got %d", len(cards))
	}
	if cards[0].Value.(string) != "first part, with second part." {
		t.Fatalf("joined: %q", cards[0].Value)
	}
}

func TestEncodeRoundTripCard(t *testing.T) {
	// A card built in memory must encode to 80 bytes and parse back identical.
	h := New()
	h.Set("OBJECT", "NGC 5194", "galaxy")
	h.Set("NAXIS", int64(2), "two axes")
	h.Set("BSCALE", 1.5e-3, "scale")
	h.Set("SIMPLE", true, "")

	buf, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	if len(buf)%2880 != 0 {
		t.Fatalf("buf length %d not block-aligned", len(buf))
	}
	cards, _, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	// Insertion order preserved.
	if cards[0].Key != "OBJECT" || cards[0].Value.(string) != "NGC 5194" {
		t.Fatalf("OBJECT: %+v", cards[0])
	}
	if cards[1].Key != "NAXIS" || cards[1].Value.(int64) != 2 {
		t.Fatalf("NAXIS: %+v", cards[1])
	}
	if cards[2].Key != "BSCALE" {
		t.Fatalf("BSCALE: %+v", cards[2])
	}
	if cards[3].Key != "SIMPLE" || cards[3].Value.(bool) != true {
		t.Fatalf("SIMPLE: %+v", cards[3])
	}
}

func TestEncodeLongString(t *testing.T) {
	h := New()
	long := ""
	for range 200 {
		long += "x"
	}
	h.Set("LONG", long, "")
	buf, err := Encode(h)
	if err != nil {
		t.Fatal(err)
	}
	cards, _, err := ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	if len(cards) != 1 {
		t.Fatalf("expected 1 logical card, got %d", len(cards))
	}
	if cards[0].Value.(string) != long {
		t.Fatalf("roundtrip mismatch: got %q", cards[0].Value)
	}
}

func TestHeaderGettersSetters(t *testing.T) {
	h := New()
	h.Set("STR", "hello", "")
	h.Set("INT", int64(42), "")
	h.Set("FLT", 3.14, "")
	h.Set("BOOL", true, "")

	if s, _ := h.String("str"); s != "hello" {
		t.Fatalf("str: %q", s)
	}
	if i, _ := h.Int("INT"); i != 42 {
		t.Fatalf("int: %d", i)
	}
	if f, _ := h.Float("flt"); f != 3.14 {
		t.Fatalf("flt: %v", f)
	}
	if b, _ := h.Bool("BOOL"); !b {
		t.Fatalf("bool")
	}

	// Wrong type
	if _, err := h.Int("STR"); err == nil {
		t.Fatalf("wrong type should error")
	}
	// Missing
	if _, err := h.String("MISSING"); err == nil {
		t.Fatalf("missing should error")
	}
}

func TestHeaderDelete(t *testing.T) {
	h := New()
	h.Set("A", int64(1), "")
	h.Set("B", int64(2), "")
	h.Set("C", int64(3), "")
	if err := h.Delete("B"); err != nil {
		t.Fatal(err)
	}
	if h.Has("B") {
		t.Fatal("B still present")
	}
	// Re-lookup A and C still works.
	if v, _ := h.Int("A"); v != 1 {
		t.Fatal("A")
	}
	if v, _ := h.Int("C"); v != 3 {
		t.Fatal("C")
	}
}
