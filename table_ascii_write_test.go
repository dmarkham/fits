package fits

import (
	"encoding/json"
	"math"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// ASCII table writer: Go writes, fitsverify judges, cfitsio and astropy
// read back. Expected values are the Go inputs.

// Row 5 is negative zero (a Go literal -0.0 is +0).
var asciiWriteFloats = []float64{math.Pi, 1e-300, -math.MaxFloat64, 0.1, math.SmallestNonzeroFloat64, math.Copysign(0, -1), 123456789.125, 0}

// The exact field text the FITS formats call for (§7.2.5): Fw.d rounds to
// d decimals; Ew.d / Dw.d have one digit before the point, d after it and
// an E (Dw.d: D) exponent. The texts are C printf's "%.dE" / "%.df"
// (glibc, correctly rounded), right-justified, with E -> D for Dw.d.
var asciiWantFields = map[string][]string{
	"F": {"       1.250", "      -2.500", "       0.000", " 1000000.000", "      -0.001", "   99999.999", "       3.142", "      -2.718"},
	"E": {"  3.1415926535897931E+00", " 1.0000000000000000E-300", "-1.7976931348623157E+308", "  1.0000000000000001E-01",
		" 4.9406564584124654E-324", " -0.0000000000000000E+00", "  1.2345678912500000E+08", "  0.0000000000000000E+00"},
	"D": {"  3.14159265358979312D+00", " 1.00000000000000003D-300", "-1.79769313486231571D+308", "  1.00000000000000006D-01",
		" 4.94065645841246544D-324", " -0.00000000000000000D+00", "  1.23456789125000000D+08", "  0.00000000000000000D+00"},
	"ENUL":  {" 1.00000E+00", "NULL        ", " 2.50000E-10", "-3.00000E+00", "NULL        ", " 4.00000E+00", " 5.00000E+00", " 6.00000E+00"},
	"ID":    {"-9223372036854775808", " 9223372036854775807", "                   0", "                  -1", "    9007199254740992", "    9007199254740993", "*                   ", "                   7"},
	"SMALL": {"-99999", "999999", "     0", "     1", "    -1", "    12", "   345", "     6"},
	// Fw.0 / Ew.0 keep an explicit decimal point (printf "%#.0f" /
	// "%#.0E"; fitsverify reports a field without one as an error).
	"F0": {"   12.", "   -3.", "    0.", "99999.", "    0.", "    2.", "   -0.", "    7."},
	"E0": {"  1.E+01", " -3.E+00", "  0.E+00", "  1.E+05", "  4.E-01", "  2.E+00", " -5.E-01", "  7.E+00"},
	// A TNULL exactly as wide as the field.
	"N3": {"  1", "-12", "  0", "999", "-99", "  5", "***", "  3"},
}

var asciiFloat0 = []float64{12, -3, 0, 99999, 0.4, 2.5, -0.5, 7}

func asciiWriteCols() []ASCIIColumnData {
	return []ASCIIColumnData{
		{Name: "NAME", DataString: []string{"alpha", "", "  lead", "x", "trail  ", "NONE", "z", "12345678"}},
		{Name: "ID", Unit: "count", Null: "*", DataInt64: []int64{math.MinInt64, math.MaxInt64, 0, -1, 1 << 53, 1<<53 + 1, 42, 7},
			Valid: []bool{true, true, true, true, true, true, false, true}},
		{Name: "SMALL", TForm: "I6", DataInt64: []int64{-99999, 999999, 0, 1, -1, 12, 345, 6}},
		{Name: "F", TForm: "F12.3", Unit: "Jy", DataFloat64: []float64{1.25, -2.5, 0, 1e6, -0.001, 99999.999, 3.14159, -2.71828}},
		{Name: "E", DataFloat64: asciiWriteFloats},
		{Name: "D", TForm: "D25.17", DataFloat64: asciiWriteFloats},
		{Name: "ENUL", TForm: "E12.5", Null: "NULL", DataFloat64: []float64{1, math.NaN(), 2.5e-10, -3, math.Inf(1), 4, 5, 6},
			Valid: []bool{true, false, true, true, false, true, true, true}},
		{Name: "F0", TForm: "F6.0", DataFloat64: asciiFloat0},
		{Name: "E0", TForm: "E8.0", DataFloat64: asciiFloat0},
		{Name: "N3", TForm: "I3", Null: "***", DataInt64: []int64{1, -12, 0, 999, -99, 5, 0, 3},
			Valid: []bool{true, true, true, true, true, true, false, true}},
	}
}

func writeASCIIFile(t *testing.T) string {
	t.Helper()
	path := filepath.Join(t.TempDir(), "ascii.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	if _, err := AppendASCIITable(f, nil, asciiWriteCols()); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	return path
}

// ulpDiff is the distance in units of the last place between two finite
// floats; anything involving Inf or NaN that is not bit-identical is
// "infinitely" far (so an overflow to -Inf is never "1 ulp" from -MaxFloat64).
func ulpDiff(a, b float64) uint64 {
	if math.IsInf(a, 0) || math.IsInf(b, 0) || math.IsNaN(a) || math.IsNaN(b) {
		if math.Float64bits(a) == math.Float64bits(b) {
			return 0
		}
		return math.MaxUint64
	}
	ia, ib := int64(math.Float64bits(a)), int64(math.Float64bits(b))
	if ia < 0 {
		ia = math.MinInt64 - ia
	}
	if ib < 0 {
		ib = math.MinInt64 - ib
	}
	if ia > ib {
		return uint64(ia - ib)
	}
	return uint64(ib - ia)
}

func TestWriteASCIITableInteropCfitsio(t *testing.T) {
	path := writeASCIIFile(t)
	if r := testutil.Fitsverify(t, path); !r.OK() {
		t.Errorf("fitsverify: %d warnings, %d errors\n%s", r.Warnings, r.Errors, r.Report)
	}
	cf := testutil.CfrefDump(t, path)
	gh := cf.HDUs[1]
	if gh.Type != "ASCII_TBL" || len(gh.Errors) > 0 {
		t.Fatalf("cfitsio: type %s errors %v", gh.Type, gh.Errors)
	}
	cols := asciiWriteCols()
	// The documented layout: fields in column order, one blank between
	// neighbouring fields, NAXIS1 the end of the last field; cfitsio's
	// TBCOLn and NAXIS1 must say exactly that.
	end := int64(0)
	for i, gc := range gh.Table.Columns {
		want := end + 1
		if i > 0 {
			want++ // the blank after the previous field
		}
		if gc.TBCol == nil || *gc.TBCol != want {
			t.Errorf("%s: cfitsio TBCOL %v, want %d (one blank between fields)", cols[i].Name, gc.TBCol, want)
		} else {
			end = *gc.TBCol - 1 + gc.Width
		}
	}
	if gh.Table.RowLen != end {
		t.Errorf("NAXIS1 = %d, want %d (the end of the last field)", gh.Table.RowLen, end)
	}
	// The raw field text, as cfitsio returns the rows.
	if len(gh.Table.RawRows) != len(asciiWriteFloats) {
		t.Fatalf("cfitsio returned %d raw rows", len(gh.Table.RawRows))
	}
	for i, c := range cols {
		want, ok := asciiWantFields[c.Name]
		if !ok {
			continue
		}
		gc := gh.Table.Columns[i]
		for r, row := range gh.Table.RawRows {
			lo := int(*gc.TBCol) - 1
			if got := string(row)[lo : lo+int(gc.Width)]; got != want[r] {
				t.Errorf("%s row %d: field %q, want %q (TFORM %s)", c.Name, r, got, want[r], gc.TForm)
			}
		}
	}
	for i, c := range cols {
		gc := gh.Table.Columns[i]
		where := c.Name
		if gc.TType == nil || string(*gc.TType) != c.Name {
			t.Errorf("%s: TTYPE %v", where, gc.TType)
		}
		// TUNITn from Unit, and no TUNITn card for a column without one.
		switch {
		case c.Unit != "" && (gc.TUnit == nil || string(*gc.TUnit) != c.Unit):
			t.Errorf("%s: cfitsio TUNIT %v, wrote Unit %q", where, gc.TUnit, c.Unit)
		case c.Unit == "" && gc.TUnit != nil:
			t.Errorf("%s: TUNIT %q written although no Unit was given", where, *gc.TUnit)
		}
		for r, cell := range gc.RowCells(t) {
			null := c.Valid != nil && !c.Valid[r]
			switch {
			case c.DataString != nil:
				got := strings.TrimRight(cell.Str[0].V, " ")
				if got != strings.TrimRight(c.DataString[r], " ") {
					t.Errorf("%s row %d: cfitsio %q, wrote %q", where, r, got, c.DataString[r])
				}
			case c.DataInt64 != nil:
				w := c.DataInt64[r]
				p := cell.Phys[0]
				if p.Valid == null {
					t.Errorf("%s row %d: cfitsio valid=%v, wrote null=%v", where, r, p.Valid, null)
					continue
				}
				if null {
					continue
				}
				// cfitsio parses I fields through a double: compare its
				// TDOUBLE value with the nearest double of what we wrote,
				// and its TLONGLONG value exactly where the double is exact
				// (|v| <= 2^53). The extremes are proven exact by astropy.
				if p.V != float64(w) {
					t.Errorf("%s row %d: cfitsio %v, wrote %d", where, r, p.V, w)
				}
				if w >= -(1<<53) && w <= 1<<53 && (len(cell.I64) != 1 || cell.I64[0].V != w) {
					t.Errorf("%s row %d: cfitsio int64 %v, wrote %d", where, r, cell.I64, w)
				}
			default:
				w := c.DataFloat64[r]
				p := cell.Phys[0]
				if p.Valid == null {
					t.Errorf("%s row %d: cfitsio valid=%v, wrote null=%v", where, r, p.Valid, null)
					continue
				}
				if null {
					continue
				}
				form := c.TForm
				if form == "" {
					form = "E24.16"
				}
				want := w
				switch form {
				case "F12.3":
					want = roundDecimal(w, 3)
				case "E12.5":
					want = roundSig(w, 6)
				case "F6.0":
					want = roundDecimal(w, 0)
				case "E8.0":
					want = roundSig(w, 1)
				}
				// cfitsio accumulates the digits in a double and scales by
				// pow(10, exp), which is not correctly rounded for 17-digit
				// mantissas; its result may differ from the correctly
				// rounded value by at most 1 ulp. astropy (exact) is the
				// bit-exact check, in TestWriteASCIITableInteropAstropy.
				// Subnormals lose more: pow(10,-324) underflows, so the
				// smallest subnormal comes back as 0.
				d := ulpDiff(p.V, want)
				if want == math.SmallestNonzeroFloat64 && p.V == 0 {
					t.Logf("%s row %d: cfitsio reads the smallest subnormal as 0 (pow(10,-324) underflow)", where, r)
				} else if d > 1 {
					t.Errorf("%s row %d: cfitsio %v (%d ulp from %v)", where, r, p.V, d, want)
				} else if d == 1 {
					t.Logf("%s row %d: cfitsio %v is 1 ulp from %v (not correctly rounded)", where, r, p.V, want)
				}
			}
		}
	}
	// Our own reader: exact.
	f, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, _ := f.HDU(1)
	at := h.(*ASCIITableHDU)
	ids, ok, err := ReadASCIIColumnMasked[int64](at, 2)
	if err != nil {
		t.Fatal(err)
	}
	for r, w := range cols[1].DataInt64 {
		if ok[r] != cols[1].Valid[r] || (ok[r] && ids[r] != w) {
			t.Errorf("ID row %d: %d valid=%v, wrote %d", r, ids[r], ok[r], w)
		}
	}
	for _, ci := range []int{4, 5} {
		got, _, err := ReadASCIIColumnMasked[float64](at, ci+1)
		if err != nil {
			t.Fatal(err)
		}
		for r, w := range asciiWriteFloats {
			if math.Float64bits(got[r]) != math.Float64bits(w) {
				t.Errorf("col %d row %d: read %v, wrote %v", ci+1, r, got[r], w)
			}
		}
	}
}

func roundDecimal(v float64, d int) float64 {
	s := formatFixed(v, d)
	return parseF(s)
}

func roundSig(v float64, sig int) float64 {
	return parseF(formatE(v, sig-1))
}

func TestWriteASCIITableInteropAstropy(t *testing.T) {
	path := writeASCIIFile(t)
	out := testutil.AstropyScript(t, testutil.TestdataPath(t, "tables", "astropy_dump.py"), path)
	var all [][]astropyCol
	if err := json.Unmarshal(out, &all); err != nil {
		t.Fatal(err)
	}
	f, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	seen := map[string]bool{}
	for _, ac := range all[0] {
		seen[ac.Name] = true
		checkAstropyColumn(t, "ascii.fits", f, ac, nil)
	}
	for _, n := range []string{"ID", "SMALL", "F", "E", "D", "ENUL"} {
		if !seen[n] {
			t.Errorf("astropy did not report column %s", n)
		}
	}
	// The full-precision columns read back bit-exactly through astropy too
	// (checkAstropyColumn compares astropy with our reader; this pins our
	// reader to the written values).
	h, _ := f.HDU(1)
	got, _, _ := ReadASCIIColumnMasked[float64](h.(*ASCIITableHDU), 5)
	for r, w := range asciiWriteFloats {
		if math.Float64bits(got[r]) != math.Float64bits(w) {
			t.Errorf("E row %d: %v, wrote %v", r, got[r], w)
		}
	}
}

func TestWriteASCIITableRejects(t *testing.T) {
	cases := []struct {
		name string
		col  ASCIIColumnData
		want string
	}{
		{"NaN without null", ASCIIColumnData{DataFloat64: []float64{math.NaN()}}, "cannot be written"},
		{"int too wide", ASCIIColumnData{TForm: "I3", DataInt64: []int64{12345}}, "does not fit I3"},
		{"string too wide", ASCIIColumnData{TForm: "A3", DataString: []string{"abcd"}}, "does not fit A3"},
		{"float too wide", ASCIIColumnData{TForm: "F6.2", DataFloat64: []float64{12345.5}}, "does not fit F6.2"},
		{"wrong format kind", ASCIIColumnData{TForm: "I6", DataFloat64: []float64{1}}, "does not match"},
		{"null without TNULL", ASCIIColumnData{DataInt64: []int64{1}, Valid: []bool{false}}, "no Null"},
		{"value spells TNULL", ASCIIColumnData{TForm: "A4", Null: "NONE", DataString: []string{"NONE"}}, "equals TNULL"},
		{"non-ASCII", ASCIIColumnData{DataString: []string{"\x01"}}, "not printable ASCII"},
		{"+Inf", ASCIIColumnData{DataFloat64: []float64{math.Inf(1)}}, "cannot be written"},
		{"-Inf", ASCIIColumnData{DataFloat64: []float64{math.Inf(-1)}}, "cannot be written"},
		{"TNULL wider than field", ASCIIColumnData{TForm: "I3", Null: "NULLX", DataInt64: []int64{1}}, "wider than the field"},
		{"TNULL one wider than field", ASCIIColumnData{TForm: "I3", Null: "NULL", DataInt64: []int64{1}}, "wider than the field"},
		{"TNULL not printable", ASCIIColumnData{Null: "\x01", DataInt64: []int64{1}}, "not printable ASCII"},
		{"Valid length", ASCIIColumnData{Null: "*", DataInt64: []int64{1, 2}, Valid: []bool{true}}, "Valid has 1 flags, data 2"},
		// One case per remaining guard (each could be removed before
		// without a failing test).
		{"two Data fields", ASCIIColumnData{DataInt64: []int64{1}, DataFloat64: []float64{1}}, "2 Data fields set, want exactly one"},
		{"no Data field", ASCIIColumnData{Name: "A"}, "0 Data fields set, want exactly one"},
		{"TNULL control byte", ASCIIColumnData{Null: "N\x7f", DataInt64: []int64{1}}, `column 1 (""): TNULL: `},
		{"name not ASCII", ASCIIColumnData{Name: "caf\xc3\xa9", DataInt64: []int64{1}}, "TTYPE1"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "a.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			WriteImage(f, nil, []int64{}, []float32{})
			_, err = AppendASCIITable(f, nil, []ASCIIColumnData{c.col})
			wantErr(t, err, c.want)
		})
	}
	t.Run("row counts differ", func(t *testing.T) {
		f, err := Create(filepath.Join(t.TempDir(), "a.fits"))
		if err != nil {
			t.Fatal(err)
		}
		defer f.Close()
		WriteImage(f, nil, []int64{}, []float32{})
		_, err = AppendASCIITable(f, nil, []ASCIIColumnData{{DataInt64: []int64{1, 2}}, {DataFloat64: []float64{1}}})
		wantErr(t, err, "1 rows, earlier columns have 2")
	})
	for _, tc := range []struct {
		name string
		hdr  *header.Header
		cols []ASCIIColumnData
		want string
	}{
		{"no columns", nil, nil, "no columns"},
		{"1000 columns", nil, make([]ASCIIColumnData, 1000), "1000 columns, FITS allows 999"},
		{"hdr key beyond the columns", func() *header.Header { h := header.New(); h.Set("TTYPE2", "B", ""); return h }(),
			[]ASCIIColumnData{{DataInt64: []int64{1}}}, "hdr has TTYPE2 but the table has 1 columns"},
	} {
		t.Run(tc.name, func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "a.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			WriteImage(f, nil, []int64{}, []float32{})
			_, err = AppendASCIITable(f, tc.hdr, tc.cols)
			wantErr(t, err, tc.want)
		})
	}
	for _, key := range []string{"TTYPE1", "TUNIT1", "TNULL1"} {
		t.Run("hdr "+key+" conflicts", func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "a.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			WriteImage(f, nil, []int64{}, []float32{})
			hdr := header.New()
			hdr.Set(key, "X", "")
			_, err = AppendASCIITable(f, hdr, []ASCIIColumnData{{Name: "Y", Unit: "Y", Null: "Y", DataInt64: []int64{1}}})
			wantErr(t, err, "conflicts with column 1")
		})
	}
	// Column keywords the ASCII writer does not write must not be passed
	// through from hdr (they would describe values the writer did not
	// scale or lay out).
	for _, key := range []string{"TSCAL1", "TZERO1", "TDIM1", "TDISP1", "TFORM1"} {
		t.Run("hdr "+key, func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "a.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			WriteImage(f, nil, []int64{}, []float32{})
			hdr := header.New()
			var v any = 2.0
			if key == "TDIM1" || key == "TDISP1" || key == "TFORM1" {
				v = "(1)"
			}
			hdr.Set(key, v, "")
			_, err = AppendASCIITable(f, hdr, []ASCIIColumnData{{DataFloat64: []float64{1}}})
			wantErr(t, err, "hdr "+key+" is not supported")
		})
	}
}

func formatFixed(v float64, d int) string { return strconv.FormatFloat(v, 'f', d, 64) }
func formatE(v float64, d int) string     { return strconv.FormatFloat(v, 'E', d, 64) }
func parseF(s string) float64 {
	v, err := strconv.ParseFloat(s, 64)
	if err != nil {
		panic(err)
	}
	return v
}
