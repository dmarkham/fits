package fits

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// fitsCards pads cards to 80 columns, appends END and pads the header to
// a whole 2880-byte block.
func fitsCards(cards ...string) []byte {
	var b strings.Builder
	for _, c := range append(cards, "END") {
		fmt.Fprintf(&b, "%-80s", c)
	}
	for b.Len()%2880 != 0 {
		b.WriteByte(' ')
	}
	return []byte(b.String())
}

// kdtreeHeader is a BINTABLE header shaped like the kdtree headers of the
// astrometry.net 5200-series index files (written by ESO-QFITS): the block
// BITPIX..TFIELDS appears twice, every repeat with the same value.
func kdtreeHeader(naxis1, naxis2 int, tform string, repeatValue func(key string, v int) int) []byte {
	block := func(rv func(string, int) int) []string {
		return []string{
			fmt.Sprintf("BITPIX  = %20d / 8-bits character format", rv("BITPIX", 8)),
			fmt.Sprintf("NAXIS   = %20d / Tables are 2-D char. array", rv("NAXIS", 2)),
			fmt.Sprintf("NAXIS1  = %20d / Bytes in row", rv("NAXIS1", naxis1)),
			fmt.Sprintf("NAXIS2  = %20d / No. of rows in table", rv("NAXIS2", naxis2)),
			fmt.Sprintf("PCOUNT  = %20d / Parameter count always 0", rv("PCOUNT", 0)),
			fmt.Sprintf("GCOUNT  = %20d / Group count always 1", rv("GCOUNT", 1)),
			fmt.Sprintf("TFIELDS = %20d / No. of col in table", rv("TFIELDS", 1)),
			fmt.Sprintf("TFORM1  = '%-8s' / Format of field", tform),
			"TTYPE1  = 'kdtree_header_codes' / Field label",
			"ORIGIN  = 'ESO-QFITS' / Written by QFITS",
		}
	}
	same := func(_ string, v int) int { return v }
	cards := []string{"XTENSION= 'BINTABLE'           / FITS Binary Table Extension"}
	cards = append(cards, block(same)...)
	cards = append(cards, block(repeatValue)...)
	cards = append(cards, "AN_FILE = 'CKDT    ' / This file is a code kdtree.")
	return fitsCards(cards...)
}

func writeKdtreeFile(t *testing.T, repeatValue func(string, int) int) string {
	t.Helper()
	var b []byte
	b = append(b, fitsCards(
		"SIMPLE  =                    T / Standard FITS file",
		"BITPIX  =                    8 / ASCII or bytes array",
		"NAXIS   =                    0 / Minimal header",
		"EXTEND  =                    T / There may be FITS ext",
	)...)
	// HDU 1: 3 rows of 4 bytes; HDU 2: empty table, as in the index files.
	b = append(b, kdtreeHeader(4, 3, "4A", repeatValue)...)
	data := make([]byte, 2880)
	copy(data, "abcdefghijkl")
	b = append(b, data...)
	b = append(b, kdtreeHeader(0, 0, "0A", repeatValue)...)
	path := filepath.Join(t.TempDir(), "kdtree.fits")
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return path
}

// TestDuplicateStructuralKeywordsSameValue: a structural keyword repeated
// with the same value has a single reading (cfitsio's first-card and
// last-card rules agree) and opens on every path, with the HDUs framed as
// declared. v1.2.2 refused these headers, which made every astrometry.net
// 5200-series index file unopenable.
func TestDuplicateStructuralKeywordsSameValue(t *testing.T) {
	path := writeKdtreeFile(t, func(_ string, v int) int { return v })
	for _, op := range openPaths {
		t.Run(op.name, func(t *testing.T) {
			f, err := op.open(t, path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			if n := f.NumHDU(); n != 3 {
				t.Fatalf("NumHDU = %d, want 3", n)
			}
			for i, want := range []struct{ naxis1, naxis2 int64 }{{4, 3}, {0, 0}} {
				h, err := f.HDU(i + 1)
				if err != nil {
					t.Fatal(err)
				}
				hdr := h.Header()
				n1, err1 := hdr.Int("NAXIS1")
				n2, err2 := hdr.Int("NAXIS2")
				if err1 != nil || err2 != nil || n1 != want.naxis1 || n2 != want.naxis2 {
					t.Errorf("HDU %d: NAXIS1 %d (%v), NAXIS2 %d (%v); want %d, %d", i+1, n1, err1, n2, err2, want.naxis1, want.naxis2)
				}
			}
		})
	}
}

// TestDuplicateStructuralKeywordsDifferentValueInTable: the same layout
// with one repeat changed is still refused, naming the keyword.
func TestDuplicateStructuralKeywordsDifferentValueInTable(t *testing.T) {
	for _, key := range []string{"BITPIX", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS"} {
		t.Run(key, func(t *testing.T) {
			path := writeKdtreeFile(t, func(k string, v int) int {
				if k == key {
					return v + 1
				}
				return v
			})
			for _, op := range openPaths {
				_, err := op.open(t, path)
				var e *DuplicateKeywordError
				if !errors.As(err, &e) || e.HDU != 1 || e.Keyword != key {
					t.Errorf("%s: err %v, want DuplicateKeywordError{HDU: 1, Keyword: %s}", op.name, err, key)
				}
			}
		})
	}
}
