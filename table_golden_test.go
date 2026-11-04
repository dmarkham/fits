package fits

import (
	"fmt"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Table read tests against cfitsio. For every table HDU in every committed
// sample (the cfitsio test programs' files in testdata/ and everything in
// testdata/tables/), every column is read with the Go API and every cell
// is compared with what the system libcfitsio returned for it (the cfref
// goldens). Integers, logicals, bits and strings must match exactly;
// floats must match bit-for-bit (NaN = NaN) except where noted.

// tableSample is one sample file and its golden.
type tableSample struct {
	name   string // display name
	path   string // FITS file
	golden string // golden JSON
}

func tableSamples(t *testing.T) []tableSample {
	t.Helper()
	var out []tableSample
	for _, n := range []string{"testprog.std", "testf77.std", "iter_a.fit", "iter_b.fit", "iter_c.fit"} {
		out = append(out, tableSample{n, testutil.TestdataPath(t, n), testutil.CfrefGoldenPath(t, n)})
	}
	dir := testutil.TestdataPath(t, "tables")
	fits, _ := filepath.Glob(filepath.Join(dir, "*.fits"))
	sort.Strings(fits)
	if len(fits) != 39 {
		t.Fatalf("testdata/tables has %d samples, want the 39 committed ones (see PROVENANCE.md)", len(fits))
	}
	for _, p := range fits {
		b := filepath.Base(p)
		out = append(out, tableSample{"tables/" + b, p, filepath.Join(dir, "golden", b+".json")})
	}
	return out
}

func loadGoldenFile(t *testing.T, path string) *testutil.CfrefFile {
	t.Helper()
	testutil.RequireFile(t, path)
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	cf, err := testutil.DecodeCfref(b)
	if err != nil {
		t.Fatalf("decode %s: %v", path, err)
	}
	return cf
}

// Every sample in testdata/tables has a golden and vice versa.
func TestTableGoldensCoverSamples(t *testing.T) {
	dir := testutil.TestdataPath(t, "tables")
	fits, _ := filepath.Glob(filepath.Join(dir, "*.fits"))
	gold, _ := filepath.Glob(filepath.Join(dir, "golden", "*.json"))
	have := map[string]bool{}
	for _, f := range fits {
		have[filepath.Base(f)] = true
		testutil.RequireFile(t, filepath.Join(dir, "golden", filepath.Base(f)+".json"))
	}
	for _, g := range gold {
		if n := strings.TrimSuffix(filepath.Base(g), ".json"); !have[n] {
			t.Errorf("golden %s has no sample", filepath.Base(g))
		}
	}
	for _, g := range gold {
		cf := loadGoldenFile(t, g)
		if cf.Cfitsio != "4.7.0" {
			t.Errorf("%s: made with cfitsio %s, PROVENANCE.md pins 4.7.0", g, cf.Cfitsio)
		}
	}
}

// Known, deliberate differences from cfitsio, each pinned: the test checks
// that cfitsio still says what we claim it says, and that we say what the
// FITS standard (and astropy) say.
type asciiDivergence struct {
	sample string
	hdu    int
	col    string
	row    int
	want   float64 // our value; cfitsio reports null
}

var asciiTNULLPrefixDivergences = []asciiDivergence{
	// TNULL5 = '-' and the right-justified field "-1.234E-100" fills the
	// whole E11.3 window, so it starts with '-'. cfitsio's prefix match
	// makes it null; the standard pads TNULL to the field width, so it is
	// the number. astropy 8.0.2 also reads -1.234e-100.
	{"tables/cf_ascii.fits", 1, "EVAL", 5, -1.234e-100},
}

func TestTablesMatchCfitsio(t *testing.T) {
	for _, s := range tableSamples(t) {
		t.Run(s.name, func(t *testing.T) {
			want := loadGoldenFile(t, s.golden)
			f, err := Open(s.path)
			if err != nil {
				t.Fatalf("Open: %v", err)
			}
			defer f.Close()
			if f.NumHDU() != want.NHDU {
				t.Fatalf("NumHDU %d, cfitsio %d", f.NumHDU(), want.NHDU)
			}
			ntab := 0
			inexact := 0 // cells cfitsio itself does not read exactly
			for _, gh := range want.HDUs {
				if gh.Table == nil {
					continue
				}
				ntab++
				h, err := f.HDU(gh.Index)
				if err != nil {
					t.Fatal(err)
				}
				switch v := h.(type) {
				case *BinaryTableHDU:
					if gh.Type != "BINARY_TBL" {
						t.Fatalf("HDU %d: Go BinaryTableHDU, cfitsio %s", gh.Index, gh.Type)
					}
					checkBinaryTable(t, s.name, s.path, gh, v, &inexact)
				case *ASCIITableHDU:
					if gh.Type != "ASCII_TBL" {
						t.Fatalf("HDU %d: Go ASCIITableHDU, cfitsio %s", gh.Index, gh.Type)
					}
					checkASCIITable(t, s.name, gh, v, &inexact)
				default:
					t.Fatalf("HDU %d: Go %T, cfitsio %s", gh.Index, h, gh.Type)
				}
			}
			if ntab == 0 {
				t.Fatalf("no table HDUs in golden")
			}
			why, listed := cfitsioInexactSamples[s.name]
			switch {
			case listed && inexact == 0:
				t.Errorf("cfitsio now reads every cell exactly (%q no longer holds); drop the cfitsioInexactSamples entry", why)
			case !listed && inexact > 0:
				t.Errorf("cfitsio does not read %d cells exactly; the exact model or cfitsio changed", inexact)
			case listed:
				t.Logf("pinned: cfitsio not exact on %d cells (%s); checked against the exact values", inexact, why)
			}
		})
	}
}

func checkColumnMeta(t *testing.T, where string, c Column, gc testutil.CfrefColumn, nrows int64, gnrows int64) {
	t.Helper()
	if nrows != gnrows {
		t.Errorf("%s: NumRows %d, cfitsio %d", where, nrows, gnrows)
	}
	if gc.TType != nil && c.Name != string(*gc.TType) {
		t.Errorf("%s: Name %q, cfitsio TTYPE %q", where, c.Name, *gc.TType)
	}
	if gc.TUnit != nil && c.Unit != string(*gc.TUnit) {
		t.Errorf("%s: Unit %q, cfitsio %q", where, c.Unit, *gc.TUnit)
	}
	if c.Scale != float64(gc.TScal) || c.Zero != float64(gc.TZero) {
		t.Errorf("%s: TSCAL/TZERO %v/%v, cfitsio %v/%v", where, c.Scale, c.Zero, gc.TScal, gc.TZero)
	}
}

func checkBinaryTable(t *testing.T, sample, path string, gh testutil.CfrefHDU, h *BinaryTableHDU, inexact *int) {
	cols, err := h.Columns()
	if err != nil {
		t.Fatalf("HDU %d Columns: %v", gh.Index, err)
	}
	if len(cols) != len(gh.Table.Columns) {
		t.Fatalf("HDU %d: %d columns, cfitsio %d", gh.Index, len(cols), len(gh.Table.Columns))
	}
	for i, gc := range gh.Table.Columns {
		c := cols[i]
		where := sample + " HDU " + itoa(gh.Index) + " col " + itoa(gc.N) + " " + c.Name + " " + c.TForm
		// A large column is stored as digests: verify our own reading
		// against them exactly, then check it cell by cell as before.
		gc = gc.Resolved(t, where, binaryCandidates(h, c.Index))
		checkColumnMeta(t, where, c, gc, h.NumRows(), gh.Table.NRows)
		if gc.TNull != nil {
			if c.Null == nil {
				t.Errorf("%s: TNULL %s not parsed", where, *gc.TNull)
			} else if itoa64(c.Null.(int64)) != string(*gc.TNull) {
				t.Errorf("%s: TNULL %v, cfitsio %s", where, c.Null, *gc.TNull)
			}
		}
		if card := gh.Card("TDIM" + itoa(gc.N)); card != nil && !equalI64(c.Dim, gc.TDim) {
			t.Errorf("%s: Dim %v, cfitsio %v", where, c.Dim, gc.TDim)
		}
		if gc.Typecode < 0 {
			checkVarColumn(t, where, sample, path, gh, h, c, gc, inexact)
		} else {
			checkFixedColumn(t, where, path, gh, h, c, gc, inexact)
		}
	}
}

// cfitsio type codes (fitsio.h).
const (
	cfBit        = 1
	cfByte       = 11
	cfLogical    = 14
	cfString     = 16
	cfShort      = 21
	cfLong       = 41
	cfFloat      = 42
	cfLongLong   = 81
	cfDouble     = 82
	cfComplex    = 83
	cfDblComplex = 163
)

// binaryCandidates reads a fixed numeric column as cfitsio's TDOUBLE /
// TLONGLONG / TULONGLONG views (all rows, flattened) for
// CfrefColumn.Resolved; a read that fails returns nil (the test fails).
func binaryCandidates(h *BinaryTableHDU, col int) testutil.CellCandidates {
	return testutil.CellCandidates{
		Phys: func() []testutil.NullFloat {
			v, ok, err := ReadVectorColumnMasked[float64](h, col)
			if err != nil {
				return nil
			}
			return testutil.NullFloats(v, ok)
		},
		I64: func() []testutil.NullInt64 {
			v, ok, err := ReadVectorColumnMasked[int64](h, col)
			if err != nil {
				return nil
			}
			return testutil.NullInt64s(v, ok)
		},
		U64: func() []testutil.NullUint64 {
			v, ok, err := ReadVectorColumnMasked[uint64](h, col)
			if err != nil {
				return nil
			}
			return testutil.NullUint64s(v, ok)
		},
	}
}

func checkFixedColumn(t *testing.T, where, path string, gh testutil.CfrefHDU, h *BinaryTableHDU, c Column, gc testutil.CfrefColumn, inexact *int) {
	t.Helper()
	switch gc.Typecode {
	case cfString:
		// cfitsio splits a cell into repeat/width strings (width from a
		// TDIM string array or an rAw TFORM); ReadStringArrayColumn must
		// return exactly those, and ReadColumn[string] must refuse a
		// multi-string cell instead of returning it as one string.
		rows := gc.RowCells(t)
		arr, err := ReadStringArrayColumn(h, c.Index)
		if err != nil {
			t.Errorf("%s: ReadStringArrayColumn: %v", where, err)
			return
		}
		multi := false
		for r, cell := range rows {
			multi = multi || len(cell.Str) != 1
			if len(arr[r]) != len(cell.Str) {
				t.Errorf("%s row %d: %d strings, cfitsio %d", where, r, len(arr[r]), len(cell.Str))
				continue
			}
			for k, w := range cell.Str {
				// cfitsio: a string whose first byte is NUL is a null
				// string (V == ""), and an all-blank one comes back as " ".
				if want := strings.TrimRight(w.V, " "); arr[r][k] != want {
					t.Errorf("%s row %d string %d: %q, cfitsio %q", where, r, k, arr[r][k], want)
				}
			}
		}
		got, err := ReadColumn[string](h, c.Index)
		_, _, merr := ReadColumnMasked[string](h, c.Index)
		if multi {
			if err == nil || !strings.Contains(err.Error(), "ReadStringArrayColumn") || merr == nil {
				t.Errorf("%s: ReadColumn[string] on a multi-string cell: err = %v / %v, want an error naming ReadStringArrayColumn", where, err, merr)
			}
			// RowBytes is the raw row: split at cfitsio's width it gives
			// cfitsio's strings.
			for r, cell := range rows {
				row, err := h.RowBytes(int64(r))
				if err != nil {
					t.Errorf("%s row %d: RowBytes: %v", where, r, err)
					continue
				}
				raw := row[c.byteOffset : c.byteOffset+c.byteSize]
				for k, w := range cell.Str {
					if sub := cellString(raw[int64(k)*gc.Width : int64(k+1)*gc.Width]); sub != strings.TrimRight(w.V, " ") {
						t.Errorf("%s row %d string %d: RowBytes gives %q, cfitsio %q", where, r, k, sub, w.V)
					}
				}
			}
			if _, err := h.RowBytes(h.NumRows()); err == nil {
				t.Errorf("%s: RowBytes past the last row returned no error", where)
			}
			return
		}
		if err != nil || merr != nil {
			t.Errorf("%s: %v / %v", where, err, merr)
			return
		}
		for r, cell := range rows {
			if want := strings.TrimRight(cell.Str[0].V, " "); got[r] != want {
				t.Errorf("%s row %d: %q, cfitsio %q", where, r, got[r], want)
			}
		}
	case cfLogical:
		if strings.HasPrefix(where, badLogicalSample+" ") {
			// Deliberate difference: cfitsio reads the 0x01 bytes as true;
			// the library refuses them (TestInvalidLogicalByteIsAnError).
			return
		}
		vals, valid, err := ReadLogicalColumn(h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		checkLogical(t, where, vals, valid, flatLog(gc.RowCells(t)))
		if c.Repeat == 1 {
			v2, ok2, err := ReadColumnMasked[bool](h, c.Index)
			if err != nil {
				t.Errorf("%s ReadColumnMasked[bool]: %v", where, err)
			} else {
				checkLogical(t, where+" (ReadColumnMasked)", v2, ok2, flatLog(gc.RowCells(t)))
			}
		}
	case cfBit:
		bits, err := ReadBitColumn(h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		var want []int
		for _, cell := range gc.RowCells(t) {
			want = append(want, cell.Bits...)
		}
		checkBits(t, where, bits, want)
	case cfComplex, cfDblComplex:
		got, err := ReadComplexColumn[complex128](h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		var want [][2]testutil.Float
		for _, cell := range gc.RowCells(t) {
			want = append(want, cell.Cplx...)
		}
		checkComplex(t, where, got, want)
	default:
		var phys []testutil.NullFloat
		var i64 []testutil.NullInt64
		var u64 []testutil.NullUint64
		for _, cell := range gc.RowCells(t) {
			phys = append(phys, cell.Phys...)
			i64 = append(i64, cell.I64...)
			u64 = append(u64, cell.U64...)
		}
		got, valid, err := ReadVectorColumnMasked[float64](h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		exact, off := binaryExact(t, where, path, gh, gc, got, valid)
		*inexact += off
		if off == 0 {
			// cfitsio is exact here: its reads must match bit for bit.
			checkPhys(t, where, got, valid, phys)
			if len(i64) > 0 {
				gi, vi, err := ReadVectorColumnMasked[int64](h, c.Index)
				if err != nil {
					t.Errorf("%s int64: %v", where, err)
				} else {
					checkI64(t, where, gi, vi, i64)
				}
			}
			if len(u64) > 0 {
				gu, vu, err := ReadVectorColumnMasked[uint64](h, c.Index)
				if err != nil {
					t.Errorf("%s uint64: %v", where, err)
				} else {
					checkU64(t, where, gu, vu, u64)
				}
			}
		}
		checkAllTargets(t, where, exact, fixedReadersFor(h, c.Index, false))
		if c.Repeat == 1 {
			// The scalar readers must agree with the vector readers.
			checkAllTargets(t, where+" (ReadColumn)", exact, fixedReadersFor(h, c.Index, true))
		}
	}
}

func checkVarColumn(t *testing.T, where, sample, path string, gh testutil.CfrefHDU, h *BinaryTableHDU, c Column, gc testutil.CfrefColumn, inexact *int) {
	t.Helper()
	// Lengths from the descriptors.
	wantLen := make([]int64, len(gc.RowCells(t)))
	for r, cell := range gc.RowCells(t) {
		if cell.Desc == nil {
			t.Fatalf("%s row %d: no descriptor in golden", where, r)
		}
		wantLen[r] = cell.Desc[0]
	}
	checkLens := func(lens func(int) int) {
		for r := range wantLen {
			if int64(lens(r)) != wantLen[r] {
				t.Errorf("%s row %d: length %d, cfitsio descriptor %d", where, r, lens(r), wantLen[r])
			}
		}
	}
	switch -gc.Typecode {
	case cfString:
		got, err := ReadVarStringColumn(h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		for r, cell := range gc.RowCells(t) {
			want := ""
			if len(cell.Str) == 1 {
				want = cell.Str[0].V
			} else if wantLen[r] != 0 {
				t.Errorf("%s row %d: golden has no string for a %d-char cell", where, r, wantLen[r])
				continue
			}
			// cfitsio fails to read a zero-length string VLA (status 308,
			// recorded in the golden's errors); the value is "".
			if got[r] != want {
				t.Errorf("%s row %d: %q, cfitsio %q", where, r, got[r], want)
			}
		}
	case cfLogical:
		vals, valid, err := ReadVarLogicalColumn(h, c.Index)
		if sample == "tables/astropy_vla_logical_pre_fix.fits" {
			// Written by an old astropy with logical bytes that are not
			// 'T', 'F' or 0 (fitsverify: 3 errors). cfitsio reads them as
			// false; there is no right value, so we refuse.
			if err == nil || !strings.Contains(err.Error(), "invalid logical byte") {
				t.Errorf("%s: err = %v, want an invalid-logical-byte error", where, err)
			}
			return
		}
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		checkLens(func(r int) int { return len(vals.At(r)) })
		checkLogical(t, where, vals.values, valid.values, flatLog(gc.RowCells(t)))
	case cfBit:
		got, err := ReadVarBitColumn(h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		checkLens(func(r int) int { return len(got.At(r)) })
		var want []int
		for _, cell := range gc.RowCells(t) {
			want = append(want, cell.Bits...)
		}
		checkBits(t, where, got.values, want)
	case cfComplex, cfDblComplex:
		got, err := ReadVarComplexColumn[complex128](h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		checkLens(func(r int) int { return len(got.At(r)) })
		var want [][2]testutil.Float
		for _, cell := range gc.RowCells(t) {
			want = append(want, cell.Cplx...)
		}
		checkComplex(t, where, got.values, want)
	default:
		var phys []testutil.NullFloat
		var i64 []testutil.NullInt64
		var u64 []testutil.NullUint64
		for _, cell := range gc.RowCells(t) {
			phys = append(phys, cell.Phys...)
			i64 = append(i64, cell.I64...)
			u64 = append(u64, cell.U64...)
		}
		got, valid, err := ReadVarColumnMasked[float64](h, c.Index)
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		checkLens(func(r int) int { return len(got.At(r)) })
		exact, off := binaryExact(t, where, path, gh, gc, got.values, valid.values)
		*inexact += off
		if off == 0 {
			checkPhys(t, where, got.values, valid.values, phys)
			if len(i64) > 0 {
				gi, vi, err := ReadVarColumnMasked[int64](h, c.Index)
				if err != nil {
					t.Errorf("%s int64: %v", where, err)
				} else {
					checkI64(t, where, gi.values, vi.values, i64)
				}
			}
			if len(u64) > 0 {
				gu, vu, err := ReadVarColumnMasked[uint64](h, c.Index)
				if err != nil {
					t.Errorf("%s uint64: %v", where, err)
				} else {
					checkU64(t, where, gu.values, vu.values, u64)
				}
			}
		}
		checkAllTargets(t, where, exact, varReadersFor(h, c.Index))
	}
}

func checkASCIITable(t *testing.T, sample string, gh testutil.CfrefHDU, h *ASCIITableHDU, inexact *int) {
	cols, err := h.Columns()
	if err != nil {
		t.Fatalf("HDU %d Columns: %v", gh.Index, err)
	}
	if len(cols) != len(gh.Table.Columns) {
		t.Fatalf("HDU %d: %d columns, cfitsio %d", gh.Index, len(cols), len(gh.Table.Columns))
	}
	for i, gc := range gh.Table.Columns {
		c := cols[i]
		where := sample + " HDU " + itoa(gh.Index) + " col " + itoa(gc.N) + " " + c.Name + " " + c.TForm
		checkColumnMeta(t, where, c, gc, h.NumRows(), gh.Table.NRows)
		if gc.TBCol == nil || int64(c.byteOffset+1) != *gc.TBCol || int64(c.Width) != gc.Width {
			t.Errorf("%s: TBCOL/width %d/%d, cfitsio %v/%d", where, c.byteOffset+1, c.Width, gc.TBCol, gc.Width)
		}
		if gc.TNull != nil && (c.Null == nil || c.Null.(string) != string(*gc.TNull)) {
			t.Errorf("%s: TNULL %v, cfitsio %q", where, c.Null, *gc.TNull)
		}
		if gc.Typecode == cfString {
			got, valid, err := ReadASCIIColumnMasked[string](h, c.Index)
			if err != nil {
				t.Errorf("%s: %v", where, err)
				continue
			}
			for r, cell := range gc.RowCells(t) {
				w := cell.Str[0]
				// cfitsio returns " " for an all-blank field; we return "".
				if valid[r] != w.Valid || (w.Valid && got[r] != strings.TrimRight(w.V, " ")) {
					t.Errorf("%s row %d: %q valid=%v, cfitsio %q valid=%v", where, r, got[r], valid[r], w.V, w.Valid)
				}
			}
			continue
		}
		exact, off := asciiExact(t, where, gh, gc)
		*inexact += off
		for r, cell := range gc.RowCells(t) {
			cfNull := !cell.Phys[0].Valid
			if cfNull == exact[r].null {
				continue
			}
			diverges := false
			for _, d := range asciiTNULLPrefixDivergences {
				if d.sample == sample && d.hdu == gh.Index && d.col == c.Name && d.row == r {
					diverges = true
					if !cfNull || exact[r].null || exact[r].f64() != d.want {
						t.Errorf("%s row %d: divergence entry no longer matches (cfitsio null %v, standard %v)", where, r, cfNull, exact[r])
					}
				}
			}
			if !diverges {
				t.Errorf("%s row %d: cfitsio null=%v, FITS 7.2.2 TNULL rule null=%v", where, r, cfNull, exact[r].null)
			}
		}
		checkAllTargets(t, where, exact, asciiReadersFor(h, c.Index))
	}
}

func flatLog(rows []testutil.CfrefCell) []testutil.NullBool {
	var out []testutil.NullBool
	for _, cell := range rows {
		out = append(out, cell.Log...)
	}
	return out
}

func checkLogical(t *testing.T, where string, vals, valid []bool, want []testutil.NullBool) {
	t.Helper()
	if len(vals) != len(want) {
		t.Errorf("%s: %d logicals, cfitsio %d", where, len(vals), len(want))
		return
	}
	for i, w := range want {
		if valid[i] != w.Valid || (w.Valid && vals[i] != w.V) {
			t.Errorf("%s elem %d: %v valid=%v, cfitsio %v valid=%v", where, i, vals[i], valid[i], w.V, w.Valid)
		}
	}
}

func checkBits(t *testing.T, where string, got []bool, want []int) {
	t.Helper()
	if len(got) != len(want) {
		t.Errorf("%s: %d bits, cfitsio %d", where, len(got), len(want))
		return
	}
	for i := range want {
		if got[i] != (want[i] == 1) {
			t.Errorf("%s bit %d: %v, cfitsio %d", where, i, got[i], want[i])
		}
	}
}

func checkComplex(t *testing.T, where string, got []complex128, want [][2]testutil.Float) {
	t.Helper()
	if len(got) != len(want) {
		t.Errorf("%s: %d complex values, cfitsio %d", where, len(got), len(want))
		return
	}
	for i, w := range want {
		re := testutil.NullFloat{V: float64(w[0]), Valid: true}
		im := testutil.NullFloat{V: float64(w[1]), Valid: true}
		okRe := sameFloat(real(got[i]), re.V) || cfitsioFlushed(real(got[i]), true, re)
		okIm := sameFloat(imag(got[i]), im.V) || cfitsioFlushed(imag(got[i]), true, im)
		if okRe && okIm {
			continue
		}
		if !sameFloat(real(got[i]), float64(w[0])) || !sameFloat(imag(got[i]), float64(w[1])) {
			t.Errorf("%s elem %d: %v, cfitsio (%v,%v)", where, i, got[i], w[0], w[1])
		}
	}
}

// checkPhys compares physical values: bit-identical floats, NaN equals
// NaN, and the same null flags. cfitsio flags NaN floats as null.
func checkPhys(t *testing.T, where string, got []float64, valid []bool, want []testutil.NullFloat) {
	t.Helper()
	if len(got) != len(want) {
		t.Errorf("%s: %d values, cfitsio %d", where, len(got), len(want))
		return
	}
	bad := 0
	for i, w := range want {
		if cfitsioFlushed(got[i], valid[i], w) {
			continue
		}
		if valid[i] != w.Valid || (w.Valid && !sameFloat(got[i], w.V)) {
			if bad < 5 {
				t.Errorf("%s elem %d: %v valid=%v, cfitsio %v valid=%v", where, i, got[i], valid[i], w.V, w.Valid)
			}
			bad++
		}
	}
	if bad > 5 {
		t.Errorf("%s: %d mismatches in total", where, bad)
	}
}

func checkI64(t *testing.T, where string, got []int64, valid []bool, want []testutil.NullInt64) {
	t.Helper()
	if len(got) != len(want) {
		t.Errorf("%s: %d int64 values, cfitsio %d", where, len(got), len(want))
		return
	}
	for i, w := range want {
		if valid[i] != w.Valid || (w.Valid && got[i] != w.V) {
			t.Errorf("%s elem %d: int64 %d valid=%v, cfitsio %d valid=%v", where, i, got[i], valid[i], w.V, w.Valid)
		}
	}
}

func checkU64(t *testing.T, where string, got []uint64, valid []bool, want []testutil.NullUint64) {
	t.Helper()
	if len(got) != len(want) {
		t.Errorf("%s: %d uint64 values, cfitsio %d", where, len(got), len(want))
		return
	}
	for i, w := range want {
		if valid[i] != w.Valid || (w.Valid && got[i] != w.V) {
			t.Errorf("%s elem %d: uint64 %d valid=%v, cfitsio %d valid=%v", where, i, got[i], valid[i], w.V, w.Valid)
		}
	}
}

// cfitsioFlushed reports whether got differs from cfitsio's value only by
// cfitsio's null-checking float read (fits_read_colnull -> fnan()): it
// flags ±Inf as null and returns +0 for -0 and subnormals. We keep the
// IEEE value; TestTablesMatchAstropy proves it bit-for-bit against astropy.
func cfitsioFlushed(got float64, valid bool, w testutil.NullFloat) bool {
	if !w.Valid && valid && math.IsInf(got, 0) {
		return true
	}
	if w.Valid && valid && math.Float64bits(w.V) == 0 && got != 0 && math.Abs(got) < 0x1p-126 {
		return true // subnormal (float32 or float64 range) flushed to +0
	}
	if w.Valid && valid && math.Float64bits(w.V) == 0 && got == 0 && math.Signbit(got) {
		return true // -0 flushed to +0
	}
	return false
}

// sameFloat (images_truth_test.go): equal bits, or both NaN.

func equalI64(a, b []int64) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func itoa(i int) string { return strconv.Itoa(i) }

func itoa64(i int64) string { return strconv.FormatInt(i, 10) }

// badLogicalSample holds logical bytes that are neither 'T', 'F' nor 0.
const badLogicalSample = "tables/cf_bintable_badlogical.fits"

// TestInvalidLogicalByteIsAnError pins the documented deliberate
// difference "invalid logical bytes are an error" on the fixed-width
// paths. tables/cf_bintable_badlogical.fits (gen_tables.c) is a
// cfitsio-written table whose 1L column (row 0) and 3L column (row 1,
// element 1) carry the byte 0x01. cfitsio 4.7.0 reads both as true
// (golden); every reader of the library must refuse them with "invalid
// logical byte" rather than return a value or a null.
func TestInvalidLogicalByteIsAnError(t *testing.T) {
	g := loadGoldenFile(t, testutil.TestdataPath(t, "tables/golden/cf_bintable_badlogical.fits.json"))
	cols := g.HDUs[1].Table.Columns
	if len(cols) != 3 || cols[0].RowCells(t)[0].Log[0] != (testutil.NullBool{V: true, Valid: true}) ||
		cols[1].RowCells(t)[1].Log[1] != (testutil.NullBool{V: true, Valid: true}) {
		t.Fatalf("cfitsio golden does not read the patched bytes as true: %+v", cols)
	}
	f, err := Open(testutil.TestdataPath(t, badLogicalSample))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	bt := h.(*BinaryTableHDU)
	raw, err := bt.RowBytes(0)
	if err != nil || raw[0] != 0x01 {
		t.Fatalf("sample row 0 = %x %v, want the patched 0x01 first", raw, err)
	}
	wantErr := func(what string, err error) {
		t.Helper()
		if err == nil || !strings.Contains(err.Error(), "invalid logical byte") {
			t.Errorf("%s: err = %v, want an invalid-logical-byte error", what, err)
		}
	}
	_, err = ReadColumn[bool](bt, 1)
	wantErr("ReadColumn[bool] L1", err)
	_, _, err = ReadColumnMasked[bool](bt, 1)
	wantErr("ReadColumnMasked[bool] L1", err)
	_, _, err = ReadLogicalColumn(bt, 1)
	wantErr("ReadLogicalColumn L1", err)
	_, _, err = ReadLogicalColumn(bt, 2)
	wantErr("ReadLogicalColumn L3", err)
	// The untouched column still reads.
	if v, err := ReadColumn[int32](bt, 3); err != nil || fmt.Sprint(v) != "[1 2 3]" {
		t.Errorf("J1 = %v %v", v, err)
	}
}
