package fits

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
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

// fitsData pads a data unit to a whole block.
func fitsData(b []byte) []byte {
	out := make([]byte, (len(b)+2879)/2880*2880)
	copy(out, b)
	return out
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
	b = append(b, fitsData([]byte("abcdefghijkl"))...)
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
	for _, key := range []string{"BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS"} {
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

// with returns base followed by extra, never sharing base's array.
func with(base []string, extra ...string) []string {
	return append(slices.Clone(base), extra...)
}

// TestDuplicateStructuralValueComparison pins how a repeat is compared:
// with its type (an integer 8 and a float 8.0 differ, as cfitsio's
// typed reads of them would), XTENSION after the trim and upper-casing
// its first card gets, and SIMPLE as a logical. The accepted cases open
// on every path and read as the first card says.
func TestDuplicateStructuralValueComparison(t *testing.T) {
	prim := []string{
		"SIMPLE  =                    T",
		"BITPIX  =                    8",
		"NAXIS   =                    0",
		"EXTEND  =                    T",
	}
	tbl := []string{
		"XTENSION= 'BINTABLE'",
		"BITPIX  =                    8",
		"NAXIS   =                    2",
		"NAXIS1  =                    4",
		"NAXIS2  =                    1",
		"PCOUNT  =                    0",
		"GCOUNT  =                    1",
		"TFIELDS =                    1",
		"TFORM1  = '4A      '",
	}
	img := []string{
		"SIMPLE  =                    T",
		"BITPIX  =                   16",
		"NAXIS   =                    1",
		"NAXIS1  =                    4",
	}
	cases := []struct {
		name    string
		primary []string
		ext     []string // nil: the primary image alone
		refused *DuplicateKeywordError
	}{
		{"BITPIX int then float", prim, with(tbl, "BITPIX  =                  8.0"), &DuplicateKeywordError{HDU: 1, Keyword: "BITPIX"}},
		{"BITPIX float then int", prim, with([]string{tbl[0], "BITPIX  =                  8.0"}, with(tbl[2:], tbl[1])...), &DuplicateKeywordError{HDU: 1, Keyword: "BITPIX"}},
		{"XTENSION lower-case padded", prim, with(tbl, "XTENSION= 'bintable  '"), nil},
		{"XTENSION different", prim, with(tbl, "XTENSION= 'IMAGE   '"), &DuplicateKeywordError{HDU: 1, Keyword: "XTENSION"}},
		{"NAXIS different", prim, with(tbl, "NAXIS   =                    3"), &DuplicateKeywordError{HDU: 1, Keyword: "NAXIS"}},
		{"SIMPLE repeated T", with(prim, "SIMPLE  =                    T"), tbl, nil},
		{"SIMPLE T then F", with(prim, "SIMPLE  =                    F"), tbl, &DuplicateKeywordError{HDU: 0, Keyword: "SIMPLE"}},
		{"image NAXIS and NAXIS1 repeated", with(img, "NAXIS   =                    1", "NAXIS1  =                    4"), nil, nil},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			var b []byte
			b = append(b, fitsCards(c.primary...)...)
			if c.ext != nil {
				b = append(b, fitsCards(c.ext...)...)
				b = append(b, fitsData([]byte("abcd"))...)
			} else {
				b = append(b, fitsData([]byte{0, 1, 0, 2, 0, 3, 0, 4})...)
			}
			path := filepath.Join(t.TempDir(), "dup.fits")
			if err := os.WriteFile(path, b, 0o644); err != nil {
				t.Fatal(err)
			}
			for _, op := range openPaths {
				f, err := op.open(t, path)
				if c.refused != nil {
					var e *DuplicateKeywordError
					if !errors.As(err, &e) || *e != *c.refused {
						t.Errorf("%s: err %v, want %+v", op.name, err, *c.refused)
					}
					continue
				}
				if err != nil {
					t.Fatalf("%s: %v", op.name, err)
				}
				if c.ext != nil {
					if n := f.NumHDU(); n != 2 {
						t.Fatalf("%s: NumHDU = %d, want 2", op.name, n)
					}
					h, err := f.HDU(1)
					if err != nil {
						t.Fatalf("%s: %v", op.name, err)
					}
					if bt, ok := h.(*BinaryTableHDU); !ok {
						t.Errorf("%s: HDU 1 is %T, want a 1-row binary table", op.name, h)
					} else if bt.NumRows() != 1 {
						t.Errorf("%s: HDU 1 has %d rows, want 1", op.name, bt.NumRows())
					}
				} else {
					h, err := f.HDU(0)
					if err != nil {
						t.Fatalf("%s: %v", op.name, err)
					}
					px, err := ReadPixels[int16](h.(*ImageHDU))
					if err != nil || !slices.Equal(px, []int16{1, 2, 3, 4}) {
						t.Errorf("%s: pixels %v, %v; want [1 2 3 4]", op.name, px, err)
					}
				}
				f.Close()
			}
		})
	}
}

// TestEditOfRepeatedStructuralKeywordNamesTheCause: Header.Set changes
// the first card of a repeated keyword, so a structural edit on a header
// that repeats it (accepted at Open) leaves two different values. Flush
// and CopyHDUWithHeader refuse, say why, and write nothing. A repeat the
// edit itself adds to a header that had one card gets no such hint.
func TestEditOfRepeatedStructuralKeywordNamesTheCause(t *testing.T) {
	path := writeKdtreeFile(t, func(_ string, v int) int { return v })
	before, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	wantParts := []string{`structural keyword "NAXIS2" appears more than once`, "first card", "EditFile"}
	check := func(path, what string, err error, parts []string, before []byte) {
		t.Helper()
		if !errors.Is(err, ErrDuplicateKeyword) {
			t.Fatalf("%s: err %v, want ErrDuplicateKeyword", what, err)
		}
		for _, p := range parts {
			if !strings.Contains(err.Error(), p) {
				t.Errorf("%s: err %q does not mention %q", what, err, p)
			}
		}
		after, rerr := os.ReadFile(path)
		if rerr != nil || !bytes.Equal(after, before) {
			t.Errorf("%s: file changed (%v)", what, rerr)
		}
	}
	noHint := func(what string, err error) {
		t.Helper()
		if strings.Contains(err.Error(), "first card") {
			t.Errorf("%s: err %q hints at a repeated keyword on disk; the header had one card", what, err)
		}
	}

	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	if err := h.Header().Set("NAXIS2", int64(2), ""); err != nil {
		t.Fatal(err)
	}
	check(path, "Flush", f.Flush(), wantParts, before)
	f.Close()

	err = EditFile(path, func(in *File, out *Writer) error {
		h, err := in.HDU(1)
		if err != nil {
			return err
		}
		hdr := h.Header().Clone()
		if err := hdr.Set("NAXIS2", int64(2), ""); err != nil {
			return err
		}
		return out.CopyHDUWithHeader(h, hdr)
	})
	check(path, "CopyHDUWithHeader", err, wantParts[:2], before)

	// The same edits with Add on a header that carries one BITPIX: the
	// repeat is the user's own, refused without the on-disk hint.
	plain := filepath.Join(t.TempDir(), "plain.fits")
	src, err := os.ReadFile(testutil.TestdataPath(t, "testprog.std"))
	if err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(plain, src, 0o644); err != nil {
		t.Fatal(err)
	}
	plainBefore, err := os.ReadFile(plain)
	if err != nil {
		t.Fatal(err)
	}
	f, err = OpenForEdit(plain)
	if err != nil {
		t.Fatal(err)
	}
	h, err = f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	if err := h.Header().Add("BITPIX", int64(16), ""); err != nil {
		t.Fatal(err)
	}
	err = f.Flush()
	check(plain, "Flush after Add", err, []string{`structural keyword "BITPIX" appears more than once`, "EditFile"}, plainBefore)
	noHint("Flush after Add", err)
	f.Close()
	err = EditFile(plain, func(in *File, out *Writer) error {
		h, err := in.HDU(1)
		if err != nil {
			return err
		}
		hdr := h.Header().Clone()
		if err := hdr.Add("BITPIX", int64(16), ""); err != nil {
			return err
		}
		return out.CopyHDUWithHeader(h, hdr)
	})
	check(plain, "CopyHDUWithHeader after Add", err, []string{`structural keyword "BITPIX" appears more than once`}, plainBefore)
	noHint("CopyHDUWithHeader after Add", err)
}

// qfitsRepeated are the structural keywords ESO-QFITS writes twice in
// the kdtree headers of the astrometry.net 5200-series index files.
var qfitsRepeated = []string{"BITPIX", "NAXIS", "NAXIS1", "NAXIS2", "PCOUNT", "GCOUNT", "TFIELDS"}

// anetIndexSamples are the real index files cut down by
// testdata/tables/extract_anet.py: every header verbatim, tables cut to
// their first 8 rows. kdtree lists the HDUs whose header repeats
// qfitsRepeated (the 'CKDT' and 'SKDT' kdtrees); index-4107 (2013,
// Tycho-2) has the same layout without repeats and is the control.
var anetIndexSamples = []struct {
	name   string
	nhdu   int
	kdtree []int
}{
	{"anet_5205-00.fits", 13, []int{2, 7}},
	{"anet_5206-38.fits", 13, []int{2, 7}},
	{"anet_4107.fits", 14, nil},
}

// TestAstrometryIndexFilesMatchCfitsio: the cut-down real index files
// open on every path and every HDU is framed (header start, data start,
// data end) exactly as cfitsio 4.7.0 fits_get_hduaddrll frames it, with
// the repeated keywords read as cfitsio reads them. The golden is also
// made to prove the fixture: cfitsio sees each of qfitsRepeated twice,
// with equal values, in the kdtree HDUs and once elsewhere. (The cells
// are compared in TestTablesMatchCfitsio and TestTablesMatchAstropy like
// every other table sample.)
func TestAstrometryIndexFilesMatchCfitsio(t *testing.T) {
	for _, s := range anetIndexSamples {
		path := testutil.TestdataPath(t, "tables", s.name)
		g := loadGoldenFile(t, testutil.TestdataPath(t, "tables", "golden", s.name+".json"))
		if g.NHDU != s.nhdu {
			t.Fatalf("%s: cfitsio golden has %d HDUs, the fixture is documented with %d", s.name, g.NHDU, s.nhdu)
		}
		// keys expected in HDU i: the primary has NAXIS = 0 and no
		// table keywords.
		keysOf := func(i int) []string {
			if i == 0 {
				return qfitsRepeated[:2]
			}
			return qfitsRepeated
		}
		for i, gh := range g.HDUs {
			want := 1
			if slices.Contains(s.kdtree, i) {
				want = 2
			}
			for _, k := range keysOf(i) {
				var vals []string
				for _, c := range gh.Cards {
					if string(c.Key) == k {
						vals = append(vals, string(c.Value))
					}
				}
				if len(vals) != want {
					t.Errorf("%s HDU %d: cfitsio sees %s %d times, the fixture is documented with %d", s.name, i, k, len(vals), want)
				}
				if len(vals) == 2 && vals[0] != vals[1] {
					t.Errorf("%s HDU %d: %s repeated with different values %v; the fixture is the same-value case", s.name, i, k, vals)
				}
			}
		}
		for _, op := range openPaths {
			t.Run(s.name+"/"+op.name, func(t *testing.T) {
				f, err := op.open(t, path)
				if err != nil {
					t.Fatalf("Open: %v (cfitsio opens it with %d HDUs)", err, g.NHDU)
				}
				defer f.Close()
				if f.NumHDU() != g.NHDU {
					t.Fatalf("NumHDU = %d, cfitsio %d", f.NumHDU(), g.NHDU)
				}
				for i, gh := range g.HDUs {
					r := f.hdus[i]
					if r.headerStart != gh.HeaderStart || r.dataStart != gh.DataStart || r.paddedEnd != gh.DataEnd {
						t.Errorf("HDU %d: header/data/end = %d/%d/%d, cfitsio %d/%d/%d",
							i, r.headerStart, r.dataStart, r.paddedEnd, gh.HeaderStart, gh.DataStart, gh.DataEnd)
					}
					h, err := f.HDU(i)
					if err != nil {
						t.Errorf("HDU %d: %v", i, err)
						continue
					}
					hdr := h.Header()
					for _, k := range keysOf(i) {
						var cf []int64
						for _, c := range gh.Cards {
							if string(c.Key) == k && c.Int != nil {
								cf = append(cf, *c.Int)
							}
						}
						n := 0
						for _, c := range hdr.Cards() {
							if c.Key == k {
								n++
							}
						}
						v, verr := hdr.Int(k)
						if n != len(cf) || verr != nil || len(cf) == 0 || v != cf[0] {
							t.Errorf("HDU %d: %s: %d cards, Int %d (%v); cfitsio %v", i, k, n, v, verr, cf)
						}
					}
					if bt, ok := h.(*BinaryTableHDU); ok && gh.Table != nil && bt.NumRows() != gh.Table.NRows {
						t.Errorf("HDU %d: NumRows %d, cfitsio %d", i, bt.NumRows(), gh.Table.NRows)
					}
				}
			})
		}
	}
}
