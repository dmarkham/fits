package fits

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/header"
)

// decodedCards decodes 80-column card images the way a header read from a
// file is decoded. Each element of lines is one physical card; a long
// string is given as its first card followed by its CONTINUE cards.
func decodedCards(t *testing.T, lines ...string) []header.Card {
	t.Helper()
	var buf []byte
	for _, l := range lines {
		buf = append(buf, fmt.Sprintf("%-80s", l)...)
	}
	buf = append(buf, fmt.Sprintf("%-80s", "END")...)
	for len(buf)%2880 != 0 {
		buf = append(buf, ' ')
	}
	cards, _, err := header.ParseCards(buf)
	if err != nil {
		t.Fatal(err)
	}
	return cards
}

// userCardWriters runs every writer that takes caller header cards.
var userCardWriters = []struct {
	name  string
	write func(f *File, h *header.Header) error
}{
	{"WriteImage", func(f *File, h *header.Header) error {
		_, err := WriteImage[int16](f, h, []int64{2}, []int16{1, 2})
		return err
	}},
	{"AppendBinaryTable", func(f *File, h *header.Header) error {
		if _, err := WriteImage[int16](f, nil, []int64{1}, []int16{1}); err != nil {
			return fmt.Errorf("primary: %w", err)
		}
		_, err := AppendBinaryTable(f, h, []ColumnData{{Name: "A", DataInt32: []int32{1, 2}}})
		return err
	}},
	{"AppendCompressedImage", func(f *File, h *header.Header) error {
		if _, err := WriteImage[int16](f, nil, []int64{1}, []int16{1}); err != nil {
			return fmt.Errorf("primary: %w", err)
		}
		_, err := AppendCompressedImage[int16](f, h, []int64{4, 4}, make([]int16, 16), CompressOptions{})
		return err
	}},
}

// TestWritersKeepCardsReadFromAFile: header cards decoded from a file are
// written by every writer as the bytes they were read from, including
// cards this package cannot build itself (a keyword cfitsio reads as
// KEYLONG10 or DATE.OBS, a Latin-1 comment byte, a malformed value) and a
// CONTINUE'd string. The truth is the source card bytes. Before the fix
// AppendBinaryTable and AppendCompressedImage ignored the error of
// Header.Add and silently left these cards out, and WriteImage refused
// them.
func TestWritersKeepCardsReadFromAFile(t *testing.T) {
	lines := []string{
		"KEYLONG10= 5 / cfitsio names this KEYLONG10",
		"DATE.OBS= '2020-01-01'",
		"OBSERVER= 'x' / caf\xe9",
		"BADNUM  = 1.5E3D / malformed number",
		"LONGSTR = 'first part of a long string value that goes over the card&'",
		"CONTINUE  'and continues here' / joined comment",
		"GOOD    =                    7",
	}
	src := decodedCards(t, lines...)
	if len(src) != 6 {
		t.Fatalf("decoded %d cards from %d lines, want 6", len(src), len(lines))
	}
	var want [][]byte
	for i := range src {
		var b []byte
		for _, r := range src[i].RawCards() {
			b = append(b, r[:]...)
		}
		want = append(want, b)
	}
	for _, w := range userCardWriters {
		t.Run(w.name, func(t *testing.T) {
			p := filepath.Join(t.TempDir(), "out.fits")
			f, err := Create(p)
			if err != nil {
				t.Fatal(err)
			}
			werr := w.write(f, header.FromCards(src))
			if cerr := f.Close(); werr == nil && cerr != nil {
				werr = cerr
			}
			if werr != nil {
				t.Fatalf("%s refused cards read from a file: %v", w.name, werr)
			}
			out, err := os.ReadFile(p)
			if err != nil {
				t.Fatal(err)
			}
			for i, b := range want {
				if !bytes.Contains(out, b) {
					t.Errorf("card %d (%s) is not in the written file as its source bytes:\n%q", i, src[i].Key, b)
				}
			}
		})
	}
}

// TestWritersRefuseInvalidUserCards: cards built in memory that FITS
// cannot hold, and column metadata that is not printable ASCII, make every
// writer return an error. Before the fix AppendBinaryTable and
// AppendCompressedImage returned nil and wrote the HDU without the card.
func TestWritersRefuseInvalidUserCards(t *testing.T) {
	bad := []header.Card{
		{Key: "DATE.OBS", Type: header.TypeString, Value: "2020-01-01"},
		{Key: "KEYLONG10", Type: header.TypeInt, Value: int64(5)},
		{Key: "OBSERVER", Type: header.TypeString, Value: "x", Comment: "caf\xe9"},
		{Key: "RAW", Type: header.TypeUnparsed, Value: header.RawValue{Text: "5 / injected"}},
	}
	for _, c := range bad {
		for _, w := range userCardWriters {
			t.Run(c.Key+"/"+w.name, func(t *testing.T) {
				f, err := Create(filepath.Join(t.TempDir(), "out.fits"))
				if err != nil {
					t.Fatal(err)
				}
				defer f.Close()
				if err := w.write(f, header.FromCards([]header.Card{c})); err == nil {
					t.Fatalf("%s accepted in-memory card %+v", w.name, c)
				}
			})
		}
	}
	for _, col := range []ColumnData{
		{Name: "caf\xe9", DataInt32: []int32{1}},
		{Name: "A", Unit: "\xb5m", DataInt32: []int32{1}},
		{Name: "A", Display: "I\t4", DataInt32: []int32{1}},
	} {
		f, err := Create(filepath.Join(t.TempDir(), "out.fits"))
		if err != nil {
			t.Fatal(err)
		}
		if _, err := WriteImage[int16](f, nil, []int64{1}, []int16{1}); err != nil {
			t.Fatal(err)
		}
		if _, err := AppendBinaryTable(f, nil, []ColumnData{col}); err == nil {
			t.Errorf("AppendBinaryTable accepted column %q unit %q display %q", col.Name, col.Unit, col.Display)
		}
		f.Close()
	}
}
