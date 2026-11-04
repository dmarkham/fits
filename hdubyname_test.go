package fits

import (
	"bytes"
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// fitsBytes packs HDUs given as header cards (END appended, blank-padded
// to 2880) and data bytes (zero-padded to 2880).
func fitsBytes(hdus ...struct {
	cards []string
	data  []byte
}) []byte {
	var b strings.Builder
	for _, h := range hdus {
		for _, c := range append(h.cards, "END") {
			fmt.Fprintf(&b, "%-80s", c)
		}
		for b.Len()%2880 != 0 {
			b.WriteByte(' ')
		}
		b.Write(h.data)
		for b.Len()%2880 != 0 {
			b.WriteByte(0)
		}
	}
	return []byte(b.String())
}

// TestHDUByNameFirstMatchCaseSensitive pins the HDUByName doc: the FIRST
// HDU in file order whose EXTNAME equals name, the primary included (as
// cfitsio's fits_movnam_hdu loop, which starts at HDU 1), compared
// case-sensitively (unlike fits_movnam_hdu's CASEINSEN compare), with
// trailing blanks not significant.
func TestHDUByNameFirstMatchCaseSensitive(t *testing.T) {
	img := func(first string, extname string, a, b byte) struct {
		cards []string
		data  []byte
	} {
		return struct {
			cards []string
			data  []byte
		}{[]string{first, "BITPIX  =                    8", "NAXIS   =                    1",
			"NAXIS1  =                    2", "PCOUNT  =                    0", "GCOUNT  =                    1",
			"EXTNAME = " + extname}, []byte{a, b}}
	}
	prim := img("SIMPLE  =                    T", "'SCI     '", 1, 2)
	prim.cards = append(prim.cards[:4], "EXTEND  =                    T", prim.cards[6])
	b := fitsBytes(prim,
		img("XTENSION= 'IMAGE   '", "'SCI     '", 3, 4),
		img("XTENSION= 'IMAGE   '", "'sci     '", 5, 6),
		img("XTENSION= 'IMAGE   '", "'ERR'", 7, 8),
	)
	p := filepath.Join(t.TempDir(), "names.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := Open(p)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	if f.NumHDU() != 4 {
		t.Fatalf("NumHDU %d", f.NumHDU())
	}
	for _, c := range []struct {
		name  string
		index int
		pix   []uint8
	}{
		{"SCI", 0, []uint8{1, 2}}, // the primary carries EXTNAME and comes first
		{"sci", 2, []uint8{5, 6}}, // case-sensitive: not the 'SCI' HDUs
		{"ERR", 3, []uint8{7, 8}},
	} {
		h, err := f.HDUByName(c.name)
		if err != nil {
			t.Fatalf("HDUByName(%q): %v", c.name, err)
		}
		im := h.(*ImageHDU)
		if h.Index() != c.index {
			t.Errorf("HDUByName(%q) = HDU %d, want %d", c.name, h.Index(), c.index)
		}
		px, err := ReadPixels[uint8](im)
		if err != nil || fmt.Sprint(px) != fmt.Sprint(c.pix) {
			t.Errorf("HDUByName(%q) pixels %v, %v; want %v", c.name, px, err, c.pix)
		}
	}
	for _, n := range []string{"Sci", "err", "SCI ", "NOSUCH"} {
		if n == "SCI " {
			// Trailing blanks are not significant in the card value, but
			// the name argument is compared as given.
			if _, err := f.HDUByName(n); err == nil {
				t.Errorf("HDUByName(%q) matched; the argument is compared verbatim", n)
			}
			continue
		}
		if _, err := f.HDUByName(n); err == nil {
			t.Errorf("HDUByName(%q) matched; want no HDU (case-sensitive)", n)
		}
	}
}

// cfmovnamGolden is testdata/cfmovnam/golden/<file>.json, written by
// cfmovnam (cfitsio fits_movnam_hdu).
type cfmovnamGolden struct {
	File string         `json:"file"`
	HDU  map[string]int `json:"hdu"`
}

func loadCfmovnamGolden(t *testing.T, golden string) (cfmovnamGolden, []string, string) {
	t.Helper()
	gp := testutil.TestdataPath(t, "cfmovnam", "golden", golden)
	testutil.RequireFile(t, gp)
	raw := mustRead(t, gp)
	var g cfmovnamGolden
	if err := json.Unmarshal(raw, &g); err != nil {
		t.Fatal(err)
	}
	// The names in golden order (the order cfmovnam printed them).
	dec := json.NewDecoder(bytes.NewReader(raw))
	var names []string
	depth := 0
	for {
		tok, err := dec.Token()
		if err != nil {
			break
		}
		switch v := tok.(type) {
		case json.Delim:
			if v == '{' {
				depth++
			} else if v == '}' {
				depth--
			}
		case string:
			if depth == 2 {
				names = append(names, v)
				dec.Token() // the value
			}
		}
	}
	if len(names) != len(g.HDU) || len(names) == 0 {
		t.Fatalf("golden %s: %d names in order, %d in the map", gp, len(names), len(g.HDU))
	}
	src := testutil.TestdataPath(t, "cfmovnam", g.File)
	testutil.RequireFile(t, src)
	return g, names, src
}

// With two EXTNAME cards in one header, HDUByName reads the FIRST one, as
// cfitsio's fits_movnam_hdu does (fits_read_key_str finds the first card).
// Truth: testdata/cfmovnam/golden/dupext.json, cfitsio's answer for every
// name on dupext.fits (written by cfitsio, see its PROVENANCE.md). Checked
// for a read-only file, and in ModeEdit before and after Header() has
// loaded each header (the two code paths of currentExtname).
func TestHDUByNameDuplicateEXTNAMEMatchesCfitsio(t *testing.T) {
	g, names, src := loadCfmovnamGolden(t, "dupext.json")
	check := func(t *testing.T, f *File) {
		t.Helper()
		for _, name := range names {
			want := g.HDU[name]
			h, err := f.HDUByName(name)
			got := -1
			if err == nil {
				got = recordOf(h).index
			}
			if got != want {
				t.Errorf("HDUByName(%q) = HDU %d (err %v), cfitsio fits_movnam_hdu: %d", name, got, err, want)
			}
		}
	}
	t.Run("Open", func(t *testing.T) {
		f, err := Open(src)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		check(t, f)
	})
	for _, loaded := range []bool{false, true} {
		t.Run(fmt.Sprintf("OpenForEdit/headers-loaded=%v", loaded), func(t *testing.T) {
			path, _ := copySample(t, filepath.Join("cfmovnam", g.File))
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			if loaded {
				for i := range f.NumHDU() {
					hduHeader(t, f, i)
				}
			}
			check(t, f)
		})
	}
}

// HDUByName differs from cfitsio's fits_movnam_hdu on purpose in three
// ways (the HDUByName doc): it compares case-sensitively, it does not fall
// back to HDUNAME, and an HDU whose first EXTNAME card has no string value
// has no name (cfitsio compares the value text, "5"). Where they differ
// it returns another HDU, or none. Truth: testdata/cfmovnam/golden/names.json,
// cfitsio's answer on the cfitsio-made names.fits (HDU 1 'sci', 2 'SCI',
// 3 HDUNAME 'ERR' only, 4 'ERR', 5 EXTNAME = 5 then 'NUM', 6 'NUM'). Each
// name's expected HDU is pinned here with whether cfitsio agrees, and the
// golden must say the same, so every difference the doc names is real and
// every agreement is cfitsio's. Checked for a read-only file and in
// ModeEdit with the headers loaded.
func TestHDUByNameDifferencesFromCfitsio(t *testing.T) {
	g, names, src := loadCfmovnamGolden(t, "names.json")
	want := map[string]struct {
		hdu int
		why string // "" where cfitsio gives the same HDU
	}{
		"SCI":  {2, "case-sensitive: HDU 1 is 'sci'"},
		"sci":  {1, ""},
		"Sci":  {-1, "case-sensitive: no HDU is 'Sci'"},
		"ERR":  {4, "no HDUNAME fallback: HDU 3 has only HDUNAME = 'ERR'"},
		"NUM":  {6, ""}, // HDU 5's first EXTNAME card (= 5) names it, not the second
		"5":    {-1, "a non-string first EXTNAME card is no name: HDU 5 has EXTNAME = 5"},
		"NONE": {-1, ""},
	}
	if len(names) != len(want) {
		t.Fatalf("golden has names %v; the test pins %d", names, len(want))
	}
	check := func(t *testing.T, f *File) {
		t.Helper()
		for _, name := range names {
			w, ok := want[name]
			if !ok {
				t.Fatalf("golden name %q is not pinned here", name)
			}
			cf := g.HDU[name]
			if (w.why == "") != (cf == w.hdu) {
				t.Errorf("%q: cfitsio fits_movnam_hdu gives HDU %d, the test pins %d with difference %q", name, cf, w.hdu, w.why)
			}
			h, err := f.HDUByName(name)
			got := -1
			if err == nil {
				got = recordOf(h).index
			}
			if got != w.hdu {
				t.Errorf("HDUByName(%q) = HDU %d (err %v), want %d (cfitsio: %d; %s)", name, got, err, w.hdu, cf, w.why)
			}
		}
	}
	t.Run("Open", func(t *testing.T) {
		f, err := Open(src)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		check(t, f)
	})
	t.Run("OpenForEdit/headers-loaded", func(t *testing.T) {
		path, _ := copySample(t, filepath.Join("cfmovnam", g.File))
		f, err := OpenForEdit(path)
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		for i := range f.NumHDU() {
			hduHeader(t, f, i)
		}
		check(t, f)
	})
}
