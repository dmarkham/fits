package fits

import (
	"encoding/binary"
	"strings"
	"testing"
)

// String layouts whose meaning the references do not agree on are
// refused, never returned as one string cut at a NUL: a TDIM that covers
// fewer characters than the repeat count (cfitsio splits the cell into
// repeat/width strings, the standard leaves the rest undefined; cfitsio's
// own fits_write_tdim and fitsverify call it illegal), an rAw whose r is
// not a multiple of w, and a TDIM string array on a VLA column. A TDIM on
// an rAw column wins over the w (cfitsio and astropy agree; proven on
// cf_bintable_strarray.fits W8 / W42). The supported layouts are proven
// against cfitsio and astropy in TestTablesMatchCfitsio and
// TestTablesMatchAstropy (cf_bintable_strarray.fits,
// strarray_astropy.fits, cf_bintable_substr.fits).
func TestStringLayoutsRefused(t *testing.T) {
	cell := []byte("abcdefgh")
	for _, tc := range []struct {
		name  string
		cards []string
		want  string
	}{
		{"tdim short", []string{"TFORM1  = '8A'", "TDIM1   = '(4,1)'"}, "covers 4 of the 8 characters"},
		{"rAw not multiple", []string{"TFORM1  = '8A3'"}, "not a multiple"},
	} {
		p := rawTableFile(t, binHeader(8, 1, 0, append([]string{"TFIELDS =                    1"}, tc.cards...)...), cell)
		bt := openTable(t, p).(*BinaryTableHDU)
		if _, err := ReadStringArrayColumn(bt, 1); err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%s: ReadStringArrayColumn err = %v, want %q", tc.name, err, tc.want)
		}
		if _, err := ReadColumn[string](bt, 1); err == nil || !strings.Contains(err.Error(), tc.want) {
			t.Errorf("%s: ReadColumn[string] err = %v, want %q", tc.name, err, tc.want)
		}
	}
	// A zero-width TDIM string array (0A with TDIM (0,3)): w*k == repeat
	// holds (0 == 0), but there is no string width to split by; refused
	// rather than returned as three empty strings.
	p0 := rawTableFile(t, binHeader(8, 1, 0, "TFIELDS =                    2", "TFORM1  = '0A'", "TDIM1   = '(0,3)'", "TFORM2  = '8A'"), cell)
	bt0 := openTable(t, p0).(*BinaryTableHDU)
	if _, err := ReadStringArrayColumn(bt0, 1); err == nil || !strings.Contains(err.Error(), "TDIM [0 3] covers 0 of the 0 characters") {
		t.Errorf("0A TDIM (0,3): ReadStringArrayColumn err = %v", err)
	}
	if got, err := ReadColumn[string](bt0, 2); err != nil || got[0] != "abcdefgh" {
		t.Errorf("8A next to 0A: %q %v", got, err)
	}
	// TDIM (4,2) on 8A: two strings; TDIM (8) on 8A: one; the same on
	// 8A2, where the TDIM wins over the substring width.
	for _, tc := range []struct {
		tform, dim string
		want       []string
	}{{"8A", "(4,2)", []string{"abcd", "efgh"}}, {"8A", "(8)", []string{"abcdefgh"}},
		{"8A2", "(4,2)", []string{"abcd", "efgh"}}, {"8A2", "(8)", []string{"abcdefgh"}}} {
		p := rawTableFile(t, binHeader(8, 1, 0, "TFIELDS =                    1", "TFORM1  = '"+tc.tform+"'", "TDIM1   = '"+tc.dim+"'"), cell)
		got, err := ReadStringArrayColumn(openTable(t, p).(*BinaryTableHDU), 1)
		if err != nil || strings.Join(got[0], "|") != strings.Join(tc.want, "|") {
			t.Errorf("%s TDIM %s: %q %v, want %q", tc.tform, tc.dim, got, err, tc.want)
		}
	}
	// A VLA string column with a TDIM string array.
	row := make([]byte, 8)
	binary.BigEndian.PutUint32(row, 8)
	p := rawTableFile(t, binHeader(8, 1, 8, "TFIELDS =                    1", "TFORM1  = '1PA(8)'", "TDIM1   = '(4,2)'"), append(row, cell...))
	if _, err := ReadVarStringColumn(openTable(t, p).(*BinaryTableHDU), 1); err == nil || !strings.Contains(err.Error(), "string arrays in VLA columns") {
		t.Errorf("VLA TDIM string array: err = %v", err)
	}
	// The writer refuses a TDIM that would declare a string array.
	f, _ := createTemp(t)
	defer f.Close()
	if _, err := AppendBinaryTable(f, nil, []ColumnData{{Name: "S", Dim: []int64{5, 2}, DataString: []string{"ab   cdefg"}}}); err == nil || !strings.Contains(err.Error(), "string arrays is not supported") {
		t.Errorf("writer Dim [5 2] on DataString: err = %v", err)
	}
	if _, err := AppendBinaryTable(f, nil, []ColumnData{{Name: "S", Dim: []int64{10}, DataString: []string{"ab   cdefg"}}}); err != nil {
		t.Errorf("writer Dim [10] on DataString: %v", err)
	}
}
