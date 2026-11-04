package fits

import (
	"bytes"
	"errors"
	"fmt"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
)

// TestWritersRefuseHierarchAlias: a caller card HIERARCH <name> whose path
// is itself a standard keyword name (HIERARCH BSCALE, BZERO, BLANK, END,
// NAXIS1, EXTNAME, TSCAL1, ...) is refused by every writer, on a ModeCreate
// file and by an append to a ModeEdit file, with the rule Flush applies
// (hierarchAlias; TestFlushRefusesHierarchAliasOfStandardKeyword proves
// that cfitsio reads such a card as that keyword, e.g. applies HIERARCH
// BSCALE as scaling, while astropy ignores it): an error wrapping
// header.ErrInvalidKeyword, no HDU added, and an edited file left
// byte-identical. A HIERARCH card with a longer path is written.
func TestWritersRefuseHierarchAlias(t *testing.T) {
	writers := map[string]func(f *File, h *header.Header) error{
		"WriteImage": func(f *File, h *header.Header) error {
			_, err := WriteImage(f, h, []int64{3}, []int16{1, 2, 3})
			return err
		},
		"AppendImage": func(f *File, h *header.Header) error {
			_, err := AppendImage(f, h, []int64{3}, []float32{1, 2, 3})
			return err
		},
		"AppendBinaryTable": func(f *File, h *header.Header) error {
			_, err := AppendBinaryTable(f, h, []ColumnData{{Name: "A", DataInt32: []int32{1, 2}}})
			return err
		},
		"AppendASCIITable": func(f *File, h *header.Header) error {
			_, err := AppendASCIITable(f, h, []ASCIIColumnData{{Name: "A", TForm: "I4", DataInt64: []int64{1, 2}}})
			return err
		},
		"AppendCompressedImage": func(f *File, h *header.Header) error {
			_, err := AppendCompressedImage(f, h, []int64{4, 3}, []int16{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}, CompressOptions{})
			return err
		},
		"AppendCompressedFloat32Image": func(f *File, h *header.Header) error {
			_, err := AppendCompressedFloat32Image(f, h, []int64{4, 3}, make([]float32, 12), CompressFloatOptions{})
			return err
		},
		"AppendCompressedFloat64Image": func(f *File, h *header.Header) error {
			_, err := AppendCompressedFloat64Image(f, h, []int64{4, 3}, make([]float64, 12), CompressFloatOptions{})
			return err
		},
	}
	cards := []struct {
		key     string
		value   any
		refused bool
	}{
		{"HIERARCH BSCALE", 2.0, true},
		{"HIERARCH BZERO", 100.0, true},
		{"HIERARCH BLANK", int64(3), true},
		{"HIERARCH END", int64(1), true},
		{"HIERARCH NAXIS1", int64(8), true},
		{"HIERARCH EXTNAME", "SCI", true},
		{"HIERARCH TSCAL1", 2.0, true},
		{"HIERARCH ESO DET BSCALE", 2.0, false},
	}
	for wname, write := range writers {
		for _, c := range cards {
			for _, mode := range []string{"create", "edit"} {
				t.Run(fmt.Sprintf("%s/%s/%s", wname, strings.ReplaceAll(c.key, " ", "_"), mode), func(t *testing.T) {
					path := filepath.Join(t.TempDir(), "w.fits")
					f, err := Create(path)
					if err != nil {
						t.Fatal(err)
					}
					if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
						t.Fatal(err)
					}
					if mode == "edit" {
						if err := f.Close(); err != nil {
							t.Fatal(err)
						}
						if f, err = OpenForEdit(path); err != nil {
							t.Fatal(err)
						}
					}
					before := mustRead(t, path)
					nhdu := f.NumHDU()
					h := header.New()
					if err := h.Set(c.key, c.value, "caller card"); err != nil {
						t.Fatalf("Set: %v", err)
					}
					werr := write(f, h)
					if !c.refused {
						if werr != nil {
							t.Fatalf("%s: %v", wname, werr)
						}
						f.Close()
						return
					}
					if werr == nil || !errors.Is(werr, header.ErrInvalidKeyword) || !strings.Contains(werr.Error(), "HIERARCH") {
						f.Close()
						t.Fatalf("%s with %s = %v, want an error wrapping header.ErrInvalidKeyword naming the card", wname, c.key, werr)
					}
					if f.NumHDU() != nhdu {
						t.Errorf("NumHDU %d after the refused write, %d before", f.NumHDU(), nhdu)
					}
					if err := f.Close(); err != nil {
						t.Fatal(err)
					}
					if !bytes.Equal(mustRead(t, path), before) {
						t.Errorf("refused write changed the file")
					}
				})
			}
		}
	}
}
