package fits

import (
	"encoding/binary"
	"errors"
	"math"
	"math/big"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// TZERO = 9.223372036854776E+18 is a float literal whose decimal value is
// 9223372036854776000, not 2^63, while its nearest double is 2^63. cfitsio
// and astropy compare the double and apply the uint64 convention; the
// literal says stored + 9223372036854776000. Neither reading can be
// proven to be the writer's, so integer targets are refused and float
// targets get the float64 arithmetic (identical under both readings).
func TestAmbiguousTZERONear2p63(t *testing.T) {
	row := make([]byte, 16)
	binary.BigEndian.PutUint64(row, uint64(5))
	binary.BigEndian.PutUint64(row[8:], 1<<63) // stored -2^63
	p := rawTableFile(t, binHeader(8, 2, 0, "TFIELDS =                    1", "TFORM1  = '1K'",
		"TZERO1  =  9.223372036854776E+18"), row)
	bt := openTable(t, p).(*BinaryTableHDU)
	if _, err := ReadColumn[uint64](bt, 1); !errors.Is(err, ErrTypeMismatch) || !strings.Contains(err.Error(), "near 2^63") {
		t.Errorf("uint64 read: err = %v, want ErrTypeMismatch naming the ambiguous TZERO", err)
	}
	if _, err := ReadColumn[int64](bt, 1); !errors.Is(err, ErrTypeMismatch) {
		t.Errorf("int64 read: err = %v, want ErrTypeMismatch", err)
	}
	f, err := ReadColumn[float64](bt, 1)
	if err != nil || f[0] != 5+9223372036854775808.0 || f[1] != 0 {
		t.Errorf("float64 read = %v, %v; want [2^63 0]", f, err)
	}
	// The exact literal 2^63, as an integer or as a float spelling of the
	// same decimal, is the uint64 convention.
	for _, lit := range []string{"9223372036854775808", "9223372036854775808.0", "9.223372036854775808E18"} {
		p := rawTableFile(t, binHeader(8, 2, 0, "TFIELDS =                    1", "TFORM1  = '1K'",
			"TZERO1  = "+strings.Repeat(" ", max(0, 20-len(lit)))+lit), row)
		u, err := ReadColumn[uint64](openTable(t, p).(*BinaryTableHDU), 1)
		if err != nil || u[0] != 1<<63+5 || u[1] != 0 {
			t.Errorf("TZERO1 = %s: uint64 read = %v, %v; want [2^63+5 0]", lit, u, err)
		}
	}
}

// The writer writes TZERO = 2^63 (DataUint64) as the exact integer literal
// cfitsio and astropy write, and refuses a hdr TSCALn/TZEROn literal it
// could not write back with the same decimal value (it holds TScal/TZero
// as float64).
func TestWriterScalingLiteralsExact(t *testing.T) {
	f, path := createTemp(t)
	if _, err := AppendBinaryTable(f, nil, []ColumnData{{Name: "U", DataUint64: []uint64{0, math.MaxUint64}}}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	bt := mustBinTable(t, r, 1)
	c, _ := bt.Header().Get("TZERO1")
	if c.Type != header.TypeBigInt || !strings.HasPrefix(strings.TrimSpace(string(c.Raw[10:])), "9223372036854775808") {
		t.Errorf("TZERO1 card %q (type %v), want the integer literal 9223372036854775808", c.Raw[:], c.Type)
	}
	if u, err := ReadColumn[uint64](bt, 1); err != nil || u[0] != 0 || u[1] != math.MaxUint64 {
		t.Errorf("read back %v, %v", u, err)
	}

	for _, tc := range []struct {
		card string
		ok   bool
	}{
		{"TZERO1  =  9223372036854775809", false}, // nearest double 2^63
		{"TZERO1  =  9223372036854775808", true},
		{"TSCAL1  =                  0.1", true},
		{"TSCAL1  = 0.10000000000000000001", false},
		{"TZERO1  =                 1000", true},
	} {
		cards, _, err := header.ParseCards([]byte(tcard(tc.card) + tcard("END") + strings.Repeat(" ", 2880-160)))
		if err != nil {
			t.Fatalf("%s: %v", tc.card, err)
		}
		hdr := header.FromCards(cards)
		g, _ := createTemp(t)
		_, err = AppendBinaryTable(g, hdr, []ColumnData{{Name: "K", DataInt64: []int64{1, 2}}})
		if (err == nil) != tc.ok {
			t.Errorf("hdr %s: err = %v, want ok=%v", tc.card, err, tc.ok)
		}
		g.Close()
	}
}

func createTemp(t *testing.T) (*File, string) {
	t.Helper()
	p := filepath.Join(t.TempDir(), "w.fits")
	f, err := Create(p)
	if err != nil {
		t.Fatal(err)
	}
	return f, p
}

// An integral TScal/TZero is written as its exact integer literal, so a
// reader that takes the literal exactly (FITS 4.0 7.3.2; this library)
// gets the float64 itself. That literal fits a card up to 70 characters
// (bytes 11-80 with no comment): 1e60 and 1e69 are written, and cfitsio
// reads back a value string whose exact decimal is the float64. Beyond
// that no exact decimal of the float64 fits a card (every float64 of that
// size is an integer of more than 70 digits), and a shorter literal such
// as 1E+100 is a different number to an exact reader, so the writer
// refuses it with an error saying so.
func TestWriterHugeScalingLiterals(t *testing.T) {
	for _, tc := range []struct {
		v     float64
		zero  bool
		fails bool
	}{
		{1e30, true, false}, {1e60, true, false}, {1e69, true, false}, {-1e68, false, false},
		{-1e69, true, true}, {1e100, true, true}, {-1e200, true, true}, {1e71, false, true},
	} {
		f, path := createTemp(t)
		v := tc.v
		col := ColumnData{Name: "D", DataFloat64: []float64{1.5, -2}}
		key := "TZERO1"
		if tc.zero {
			col.TZero = &v
		} else {
			col.TScal, key = &v, "TSCAL1"
		}
		_, err := AppendBinaryTable(f, nil, []ColumnData{col})
		f.Close()
		if tc.fails {
			if err == nil || !strings.Contains(err.Error(), "no exact decimal literal") {
				t.Errorf("%s = %g: err = %v, want a refusal naming the missing exact decimal literal", key, tc.v, err)
			}
			continue
		}
		if err != nil {
			t.Errorf("%s = %g: %v", key, tc.v, err)
			continue
		}
		cf := testutil.CfrefDump(t, path)
		c := cf.HDUs[1].Card(key)
		if c == nil {
			t.Fatalf("%s = %g: cfitsio sees no %s", key, tc.v, key)
		}
		lit, ok := new(big.Rat).SetString(strings.TrimSpace(string(c.Value)))
		if !ok || lit.Cmp(new(big.Rat).SetFloat64(tc.v)) != 0 {
			t.Errorf("%s = %g: cfitsio reads the value string %q, not the exact decimal of the float64", key, tc.v, c.Value)
		}
	}
}

// A non-integral TScal/TZero is written as the shortest decimal that
// reads back as the same float64 ("0.1", "0.3333333333333333"), not as
// the float64's exact decimal (0.1000000000000000055511151231257827...),
// which no card could hold for most values. A reader that takes the
// literal exactly, as this library does (FITS 4.0 4.2.4: the value is the
// decimal written), then scales by that decimal: stored 10 and 30 with
// TScal 0.1 are exactly 1 and 3 for integer targets, and stored 3 with
// TScal 1/3 is 0.9999999999999999, not an integer, while float targets
// get the float64 arithmetic (1, as cfitsio computes).
func TestWriterNonIntegralScalingLiteral(t *testing.T) {
	third := 1.0 / 3
	tenth := 0.1
	half := -2.5
	f, path := createTemp(t)
	if _, err := AppendBinaryTable(f, nil, []ColumnData{
		{Name: "K", DataInt64: []int64{10, 30}, TScal: &tenth},
		{Name: "T", DataInt64: []int64{3, 6}, TScal: &third, TZero: &half},
	}); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	bt := mustBinTable(t, r, 1)
	for key, want := range map[string]string{"TSCAL1": "0.1", "TSCAL2": "0.3333333333333333", "TZERO2": "-2.5"} {
		c, _ := bt.Header().Get(key)
		lit, ok := numericLiteral(c)
		if !ok || lit != want {
			t.Errorf("%s literal %q, want %q (card %q)", key, lit, want, c.Raw[:])
		}
		if v, err := strconv.ParseFloat(want, 64); err != nil || (key == "TSCAL1" && v != tenth) || (key == "TSCAL2" && v != third) {
			t.Errorf("%s: %q does not read back as the float64 written", key, want)
		}
	}
	if lit, _ := new(big.Rat).SetString("0.1"); lit.Cmp(new(big.Rat).SetFloat64(tenth)) == 0 {
		t.Errorf("the literal 0.1 is the float64 0.1 exactly; the test premise is wrong")
	}
	if v, err := ReadColumn[int64](bt, 1); err != nil || v[0] != 1 || v[1] != 3 {
		t.Errorf("TScal 0.1 as int64: %v %v, want [1 3] (10 * 0.1 and 30 * 0.1 exactly)", v, err)
	}
	if _, err := ReadColumn[int64](bt, 2); !errors.Is(err, ErrTypeMismatch) {
		t.Errorf("TScal 0.3333333333333333 as int64: err = %v, want ErrTypeMismatch (3 * 0.3333333333333333 - 2.5 is not an integer)", err)
	}
	// float64(...) keeps the expected product rounded on its own: without
	// it Go may fuse x*y-2.5 into one FMA (arm64), which the library, like
	// cfitsio's TDOUBLE read, does not do.
	if v, err := ReadColumn[float64](bt, 2); err != nil || v[0] != float64(float64(3)*third)-2.5 || v[1] != float64(float64(6)*third)-2.5 {
		t.Errorf("TScal 1/3 as float64: %v %v, want the float64 arithmetic", v, err)
	}
}
