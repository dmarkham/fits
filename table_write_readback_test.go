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

// Write tests: the Go writer produces tables with every column type, the
// unsigned TZERO convention, TSCAL/TZERO/TNULL (from ColumnData and from
// a user header), nulls, vectors with TDIM, and P and Q VLAs of every
// element type. Each file must be
//   - accepted by fitsverify with 0 warnings and 0 errors,
//   - read back by cfitsio (the cfref harness) with exactly the values we
//     meant to write (the Go input slices are the expected values; cfitsio
//     is the independent reader), and
//   - read back by astropy with the same numbers (bit-exact).
// Our own reader is checked last, as a round trip.

func f64p(v float64) *float64 { return &v }
func i64p(v int64) *int64     { return &v }

// writeCase is one column and how cfitsio must report it.
type writeCase struct {
	col ColumnData
	// Expected per-element values, flattened over rows (VLA: rows
	// concatenated). Exactly one of these is set.
	phys  []float64 // physical value; NaN = null for floats
	nulls []bool    // optional: element is null (cfitsio null flag)
	i64   []int64   // also check cfitsio's TLONGLONG read
	u64   []uint64  // also check cfitsio's TULONGLONG read
	log   []int     // 1, 0, -1 = undefined
	bits  []bool
	str   []string
	cplx  []complex128
	lens  []int64 // VLA row lengths
}

func fixedWriteCases() []writeCase {
	nan := math.NaN()
	k1 := []int64{math.MinInt64, -1, 1<<62 + 1, 4295806720896741504, 6917528997577384320, math.MaxInt64}
	u64 := []uint64{0, 1, 1 << 63, 1<<63 - 1, 1<<64 - 2, 1<<64 - 1}
	bits13 := make([]bool, 6*13)
	for i := range bits13 {
		bits13[i] = (i*7+i/5)%3 == 0
	}
	return []writeCase{
		{col: ColumnData{Name: "L1", DataBool: []bool{true, false, false, true, false, true}, Valid: []bool{true, true, false, true, true, false}},
			log: []int{1, 0, -1, 1, 0, -1}},
		{col: ColumnData{Name: "L3", TForm: "3L", DataBool: make([]bool, 18)}, log: make([]int, 18)},
		{col: ColumnData{Name: "X1", DataBits: []bool{true, false, true, true, false, false}},
			bits: []bool{true, false, true, true, false, false}},
		{col: ColumnData{Name: "X13", TForm: "13X", DataBits: bits13}, bits: bits13},
		{col: ColumnData{Name: "B1", DataUint8: []uint8{0, 1, 127, 128, 254, 255}},
			phys: []float64{0, 1, 127, 128, 254, 255}, i64: []int64{0, 1, 127, 128, 254, 255}},
		{col: ColumnData{Name: "B4", Dim: []int64{2, 2}, DataUint8: seqU8(24)}, phys: seqF(24, func(i int) float64 { return float64(uint8(i * 11)) })},
		{col: ColumnData{Name: "I1", DataInt16: []int16{math.MinInt16, -1, 0, 1, 12345, math.MaxInt16}},
			i64: []int64{math.MinInt16, -1, 0, 1, 12345, math.MaxInt16}},
		{col: ColumnData{Name: "J1", DataInt32: []int32{math.MinInt32, -1, 0, 1, 123456789, math.MaxInt32}},
			i64: []int64{math.MinInt32, -1, 0, 1, 123456789, math.MaxInt32}},
		{col: ColumnData{Name: "K1", DataInt64: k1}, i64: k1},
		{col: ColumnData{Name: "A12", TForm: "12A", DataString: []string{"hello", "", "trailing  ", "exactly12chr", "a b c", " lead"}},
			str: []string{"hello", "", "trailing", "exactly12chr", "a b c", " lead"}},
		{col: ColumnData{Name: "E1", DataFloat32: []float32{1.5, -2.25, math.MaxFloat32, 1e-30, float32(nan), 0}},
			phys: []float64{1.5, -2.25, math.MaxFloat32, float64(float32(1e-30)), nan, 0}},
		{col: ColumnData{Name: "E3", Dim: []int64{3}, DataFloat32: []float32{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18}},
			phys: seqF(18, func(i int) float64 { return float64(i + 1) })},
		{col: ColumnData{Name: "D1", Display: "G25.16", Unit: "deg", DataFloat64: []float64{math.Pi, -1e300, 0.1, nan, 1.0 / 3, -7}},
			phys: []float64{math.Pi, -1e300, 0.1, nan, 1.0 / 3, -7}},
		// Valid=false on a float column writes NaN (the float null).
		{col: ColumnData{Name: "EN", DataFloat32: []float32{1, 2, 3, 4, 5, 6}, Valid: []bool{true, false, true, true, false, true}},
			phys: []float64{1, nan, 3, 4, nan, 6}},
		{col: ColumnData{Name: "DN", DataFloat64: []float64{1, 2, 3, 4, 5, 6}, Valid: []bool{false, true, true, true, true, false}},
			phys: []float64{nan, 2, 3, 4, 5, nan}},
		{col: ColumnData{Name: "D6", Dim: []int64{3, 2}, DataFloat64: seqF(36, func(i int) float64 { return float64(i)*1.25 - 7 })},
			phys: seqF(36, func(i int) float64 { return float64(i)*1.25 - 7 })},
		{col: ColumnData{Name: "C1", DataComplex64: []complex64{1 - 2i, 3 - 4i, 0.5 + 0.25i, -1e30 + 1e-30i, 0, 7}},
			cplx: []complex128{1 - 2i, 3 - 4i, 0.5 + 0.25i, complex(float64(float32(-1e30)), float64(float32(1e-30))), 0, 7}},
		{col: ColumnData{Name: "M2", TForm: "2M", DataComplex128: []complex128{math.Pi - math.E*1i, 1e300, 0.1 + 0.2i, 1i, 7, 8, 9, -10i, 1, 2, 3, 4}},
			cplx: []complex128{math.Pi - math.E*1i, 1e300, 0.1 + 0.2i, 1i, 7, 8, 9, -10i, 1, 2, 3, 4}},
		{col: ColumnData{Name: "S1", DataInt8: []int8{-128, -1, 0, 1, 100, 127}}, i64: []int64{-128, -1, 0, 1, 100, 127}},
		{col: ColumnData{Name: "U1", DataUint16: []uint16{0, 1, 32767, 32768, 40000, 65535}}, i64: []int64{0, 1, 32767, 32768, 40000, 65535}},
		{col: ColumnData{Name: "V1", DataUint32: []uint32{0, 1, 1<<31 - 1, 1 << 31, 3000000000, 1<<32 - 1}},
			i64: []int64{0, 1, 1<<31 - 1, 1 << 31, 3000000000, 1<<32 - 1}},
		{col: ColumnData{Name: "W1", DataUint64: u64}, u64: u64},
		// Stored values with TSCAL/TZERO: phys = stored*0.5 + 10.
		{col: ColumnData{Name: "JSC", TScal: f64p(0.5), TZero: f64p(10), DataInt32: []int32{0, 1, -30, 1999980, -20, 2}},
			phys: []float64{10, 10.5, -5, 1e6, 0, 11}},
		// TNULL with nulls; TSCAL 1.5 TZERO -3 given through the HEADER.
		{col: ColumnData{Name: "ISC", TNull: i64p(-32768), DataInt16: []int16{0, 2, 3, 5, -1, 0}, Valid: []bool{true, true, true, true, true, false}},
			phys: []float64{-3, 0, 1.5, 4.5, -4.5, 0}, nulls: []bool{false, false, false, false, false, true}},
		// TNULL on the top of the B range and TNULL = 0.
		{col: ColumnData{Name: "BNUL", TNull: i64p(255), DataUint8: []uint8{0, 1, 254, 0, 7, 3}, Valid: []bool{true, true, true, false, true, true}},
			i64: []int64{0, 1, 254, 0, 7, 3}, nulls: []bool{false, false, false, true, false, false}},
		{col: ColumnData{Name: "JNUL0", TNull: i64p(0), DataInt32: []int32{5, 0, -1, 0, 2147483647, -2147483648}, Valid: []bool{true, false, true, false, true, true}},
			i64: []int64{5, 0, -1, 0, 2147483647, -2147483648}, nulls: []bool{false, true, false, true, false, false}},
		{col: ColumnData{Name: "KNUL", TNull: i64p(math.MinInt64), DataInt64: []int64{-1, 0, 2, 3, 0, math.MaxInt64}, Valid: []bool{true, false, true, true, false, true}},
			i64: []int64{-1, 0, 2, 3, 0, math.MaxInt64}, nulls: []bool{false, true, false, false, true, false}},
	}
}

func varWriteCases(q bool) []writeCase {
	tf := func(code string) string {
		if q {
			return "1Q" + code
		}
		return ""
	}
	tfAlways := func(code string) string {
		if q {
			return "1Q" + code
		}
		return "1P" + code
	}
	return []writeCase{
		{col: ColumnData{Name: "VL", TForm: tf("L"), DataVarBool: [][]bool{{true, false, true}, {}, {true}, {false, false, true, true}},
			VarValid: [][]bool{{true, true, true}, {}, {true}, {true, false, true, true}}},
			log: []int{1, 0, 1, 1, 0, -1, 1, 1}, lens: []int64{3, 0, 1, 4}},
		{col: ColumnData{Name: "VX", TForm: tf("X"), DataVarBits: [][]bool{{true, false, false, true}, {false, true}, {}, {true, true, true, true, true, true, true, true, false, true, true}}},
			bits: []bool{true, false, false, true, false, true, true, true, true, true, true, true, true, true, false, true, true}, lens: []int64{4, 2, 0, 11}},
		{col: ColumnData{Name: "VB", TForm: tf("B"), TNull: i64p(77), DataVarUint8: [][]uint8{{1, 2, 3}, {}, {255}, {0, 9}},
			VarValid: [][]bool{{true, false, true}, {}, {true}, {true, true}}},
			phys: []float64{1, 0, 3, 255, 0, 9}, nulls: []bool{false, true, false, false, false, false}, lens: []int64{3, 0, 1, 2}},
		{col: ColumnData{Name: "VI", TForm: tf("I"), DataVarInt16: [][]int16{{math.MinInt16}, {1, -1}, {}, {math.MaxInt16, 0, 5}}},
			i64: []int64{math.MinInt16, 1, -1, math.MaxInt16, 0, 5}, lens: []int64{1, 2, 0, 3}},
		{col: ColumnData{Name: "VJ", TForm: tf("J"), DataVarInt32: [][]int32{{}, {math.MinInt32, math.MaxInt32}, {7}, {}}},
			i64: []int64{math.MinInt32, math.MaxInt32, 7}, lens: []int64{0, 2, 1, 0}},
		{col: ColumnData{Name: "VK", TForm: tf("K"), DataVarInt64: [][]int64{{1<<62 + 3}, {}, {math.MinInt64, math.MaxInt64}, {-5}}},
			i64: []int64{1<<62 + 3, math.MinInt64, math.MaxInt64, -5}, lens: []int64{1, 0, 2, 1}},
		{col: ColumnData{Name: "VA", TForm: tfAlways("A"), DataVarString: []string{"abc", "", "x", "hello world"}},
			str: []string{"abc", "", "x", "hello world"}, lens: []int64{3, 0, 1, 11}},
		{col: ColumnData{Name: "VE", TForm: tf("E"), DataVarFloat32: [][]float32{{1.5}, {}, {-0.25, 3}, {float32(math.NaN())}},
			VarValid: [][]bool{{true}, {}, {false, true}, {true}}},
			phys: []float64{1.5, math.NaN(), 3, math.NaN()}, lens: []int64{1, 0, 2, 1}},
		{col: ColumnData{Name: "VD", TForm: tf("D"), Unit: "km/s", DataVarFloat64: [][]float64{{math.Pi, math.E}, {1e-300}, {}, {0.1}}},
			phys: []float64{math.Pi, math.E, 1e-300, 0.1}, lens: []int64{2, 1, 0, 1}},
		{col: ColumnData{Name: "VC", TForm: tf("C"), DataVarComplex64: [][]complex64{{1 + 2i}, {}, {3 - 4i, 5}, {}}},
			cplx: []complex128{1 + 2i, 3 - 4i, 5}, lens: []int64{1, 0, 2, 0}},
		{col: ColumnData{Name: "VM", TForm: tf("M"), DataVarComplex128: [][]complex128{{}, {math.Pi + 1i}, {}, {-1e300 - 1e-300i, 2}}},
			cplx: []complex128{math.Pi + 1i, -1e300 - 1e-300i, 2}, lens: []int64{0, 1, 0, 2}},
	}
}

func seqU8(n int) []uint8 {
	out := make([]uint8, n)
	for i := range out {
		out[i] = uint8(i * 11)
	}
	return out
}

func seqF(n int, f func(int) float64) []float64 {
	out := make([]float64, n)
	for i := range out {
		out[i] = f(i)
	}
	return out
}

func writeTableFile(t *testing.T, path string, hdr *header.Header, cases []writeCase) {
	t.Helper()
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	cols := make([]ColumnData, len(cases))
	for i := range cases {
		cols[i] = cases[i].col
	}
	if _, err := AppendBinaryTable(f, hdr, cols); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
}

func TestWriteTableInteropCfitsio(t *testing.T) {
	dir := t.TempDir()
	fixed := fixedWriteCases()
	// TSCAL/TZERO for ISC come from the user header, like a header copied
	// from another file; they must not be dropped.
	hdr := header.New()
	iscIdx := -1
	for i, c := range fixed {
		if c.col.Name == "ISC" {
			iscIdx = i + 1
		}
	}
	hdr.Set("TSCAL"+strconv.Itoa(iscIdx), 1.5, "")
	hdr.Set("TZERO"+strconv.Itoa(iscIdx), int64(-3), "")
	hdr.Set("OBSERVER", "proof", "user card is kept")
	files := map[string][]writeCase{
		"fixed.fits": fixed,
		"vla_p.fits": varWriteCases(false),
		"vla_q.fits": varWriteCases(true),
	}
	for name, cases := range files {
		t.Run(name, func(t *testing.T) {
			path := filepath.Join(dir, name)
			var h *header.Header
			if name == "fixed.fits" {
				h = hdr
			}
			writeTableFile(t, path, h, cases)
			if r := testutil.Fitsverify(t, path); !r.OK() {
				t.Errorf("fitsverify: %d warnings, %d errors\n%s", r.Warnings, r.Errors, r.Report)
			}
			cf := testutil.CfrefDump(t, path)
			if cf.NHDU != 2 || cf.HDUs[1].Table == nil {
				t.Fatalf("cfitsio sees %d HDUs", cf.NHDU)
			}
			gh := cf.HDUs[1]
			if len(gh.Errors) > 0 {
				// The only cfitsio error tolerated is its failure to read a
				// zero-length string VLA (status 308), which the expected
				// values below account for.
				for _, e := range gh.Errors {
					if !strings.Contains(e, "TSTRING: status 308") {
						t.Errorf("cfitsio error: %s", e)
					}
				}
			}
			if name == "fixed.fits" {
				if c := gh.Card("OBSERVER"); c == nil || c.Str == nil || string(*c.Str) != "proof" {
					t.Errorf("user header card OBSERVER lost")
				}
			}
			for i, wc := range cases {
				gc := gh.Table.Columns[i]
				checkWritten(t, name+" "+wc.col.Name, wc, gc)
				checkWrittenLayout(t, name+" "+wc.col.Name, wc, gh, gc, name == "vla_q.fits")
			}
			if name != "fixed.fits" {
				// Every column is a VLA: a row is one descriptor per column,
				// 16 bytes for Q and 8 for P.
				per := int64(8)
				if name == "vla_q.fits" {
					per = 16
				}
				if gh.Table.RowLen != per*int64(len(cases)) {
					t.Errorf("NAXIS1 = %d, want %d (%d-byte descriptors)", gh.Table.RowLen, per*int64(len(cases)), per)
				}
			}
			// Our own reader agrees with what we wrote (round trip).
			f, err := Open(path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			tbl := mustBinTable(t, f, 1)
			inexact := 0
			checkBinaryTable(t, name, path, gh, tbl, &inexact)
			if inexact != 0 {
				t.Errorf("cfitsio does not read %d cells of our file exactly", inexact)
			}
		})
	}
}

// checkWrittenLayout checks the keywords cfitsio sees for a written
// column: TDIMn exactly when Dim was given (and equal to it), and for VLA
// columns the descriptor kind (TFORM 1Q... when asked for, else 1P...).
func checkWrittenLayout(t *testing.T, where string, wc writeCase, gh testutil.CfrefHDU, gc testutil.CfrefColumn, q bool) {
	t.Helper()
	card := gh.Card("TDIM" + itoa(gc.N))
	if wc.col.Dim != nil {
		if card == nil {
			t.Errorf("%s: no TDIM%d card written for Dim %v", where, gc.N, wc.col.Dim)
		} else if !equalI64(gc.TDim, wc.col.Dim) {
			t.Errorf("%s: cfitsio TDIM %v, wrote Dim %v", where, gc.TDim, wc.col.Dim)
		}
	} else if card != nil {
		t.Errorf("%s: TDIM%d written although no Dim was given", where, gc.N)
	}
	switch {
	case wc.col.Display != "" && (gc.TDisp == nil || string(*gc.TDisp) != wc.col.Display):
		t.Errorf("%s: cfitsio TDISP%d %v, wrote Display %q", where, gc.N, gc.TDisp, wc.col.Display)
	case wc.col.Display == "" && gc.TDisp != nil:
		t.Errorf("%s: TDISP%d %q written although no Display was given", where, gc.N, *gc.TDisp)
	}
	switch {
	case wc.col.Unit != "" && (gc.TUnit == nil || string(*gc.TUnit) != wc.col.Unit):
		t.Errorf("%s: cfitsio TUNIT%d %v, wrote Unit %q", where, gc.N, gc.TUnit, wc.col.Unit)
	case wc.col.Unit == "" && gc.TUnit != nil:
		t.Errorf("%s: TUNIT%d %q written although no Unit was given", where, gc.N, *gc.TUnit)
	}
	if wc.lens != nil {
		prefix := "1P"
		if q {
			prefix = "1Q"
		}
		if gc.Typecode >= 0 || !strings.HasPrefix(string(gc.TForm), prefix) {
			t.Errorf("%s: cfitsio TFORM %q typecode %d, want a %s VLA", where, gc.TForm, gc.Typecode, prefix)
		}
	}
}

func checkWritten(t *testing.T, where string, wc writeCase, gc testutil.CfrefColumn) {
	t.Helper()
	var phys []testutil.NullFloat
	var i64 []testutil.NullInt64
	var u64 []testutil.NullUint64
	var log []testutil.NullBool
	var bits []int
	var strs []string
	var cplx [][2]testutil.Float
	for r, cell := range gc.RowCells(t) {
		if wc.lens != nil {
			if cell.Desc == nil || cell.Desc[0] != wc.lens[r] {
				t.Errorf("%s row %d: cfitsio descriptor %v, want length %d", where, r, cell.Desc, wc.lens[r])
			}
		}
		phys = append(phys, cell.Phys...)
		i64 = append(i64, cell.I64...)
		u64 = append(u64, cell.U64...)
		log = append(log, cell.Log...)
		bits = append(bits, cell.Bits...)
		cplx = append(cplx, cell.Cplx...)
		for _, s := range cell.Str {
			strs = append(strs, s.V)
		}
		if cell.Str == nil && wc.str != nil {
			strs = append(strs, "") // zero-length string VLA: cfitsio status 308
		}
	}
	isNull := func(k int) bool { return wc.nulls != nil && wc.nulls[k] }
	switch {
	case wc.phys != nil:
		if len(phys) != len(wc.phys) {
			t.Fatalf("%s: cfitsio %d values, wrote %d", where, len(phys), len(wc.phys))
		}
		for k, w := range wc.phys {
			wantNull := isNull(k) || math.IsNaN(w)
			if phys[k].Valid == wantNull || (!wantNull && !sameFloat(phys[k].V, w)) {
				t.Errorf("%s elem %d: cfitsio %v (valid %v), wrote %v (null %v)", where, k, phys[k].V, phys[k].Valid, w, wantNull)
			}
		}
	case wc.log != nil:
		if len(log) != len(wc.log) {
			t.Fatalf("%s: cfitsio %d logicals, wrote %d", where, len(log), len(wc.log))
		}
		for k, w := range wc.log {
			if (w == -1) == log[k].Valid || (w != -1 && log[k].V != (w == 1)) {
				t.Errorf("%s elem %d: cfitsio %+v, wrote %d", where, k, log[k], w)
			}
		}
	case wc.bits != nil:
		if len(bits) != len(wc.bits) {
			t.Fatalf("%s: cfitsio %d bits, wrote %d", where, len(bits), len(wc.bits))
		}
		for k, w := range wc.bits {
			if (bits[k] == 1) != w {
				t.Errorf("%s bit %d: cfitsio %d, wrote %v", where, k, bits[k], w)
			}
		}
	case wc.str != nil:
		if len(strs) != len(wc.str) {
			t.Fatalf("%s: cfitsio %d strings, wrote %d", where, len(strs), len(wc.str))
		}
		for k, w := range wc.str {
			if strings.TrimRight(strs[k], " ") != w {
				t.Errorf("%s row %d: cfitsio %q, wrote %q", where, k, strs[k], w)
			}
		}
	case wc.cplx != nil:
		if len(cplx) != len(wc.cplx) {
			t.Fatalf("%s: cfitsio %d complex, wrote %d", where, len(cplx), len(wc.cplx))
		}
		for k, w := range wc.cplx {
			if !sameFloat(float64(cplx[k][0]), real(w)) || !sameFloat(float64(cplx[k][1]), imag(w)) {
				t.Errorf("%s elem %d: cfitsio %v, wrote %v", where, k, cplx[k], w)
			}
		}
	}
	if wc.i64 != nil {
		if len(i64) != len(wc.i64) {
			t.Fatalf("%s: cfitsio %d int64, wrote %d", where, len(i64), len(wc.i64))
		}
		for k, w := range wc.i64 {
			if i64[k].Valid == isNull(k) || (!isNull(k) && i64[k].V != w) {
				t.Errorf("%s elem %d: cfitsio int64 %+v, wrote %d (null %v)", where, k, i64[k], w, isNull(k))
			}
		}
	}
	if wc.u64 != nil {
		if len(u64) != len(wc.u64) {
			t.Fatalf("%s: cfitsio %d uint64 (eqtype %d), wrote %d", where, len(u64), gc.Eqtype, len(wc.u64))
		}
		for k, w := range wc.u64 {
			if !u64[k].Valid || u64[k].V != w {
				t.Errorf("%s elem %d: cfitsio uint64 %+v, wrote %d", where, k, u64[k], w)
			}
		}
	}
}

// astropy reads the same Go-written files back (second independent reader).
func TestWriteTableInteropAstropy(t *testing.T) {
	dir := t.TempDir()
	fixed := filepath.Join(dir, "fixed.fits")
	writeTableFile(t, fixed, nil, fixedWriteCasesNoHeaderScale())
	vla := filepath.Join(dir, "vla.fits")
	// astropy 8.0.2 rejects 1PX (see astropySkip), so the bit VLA is left out.
	var vcases []writeCase
	for _, c := range varWriteCases(false) {
		if c.col.Name != "VX" {
			vcases = append(vcases, c)
		}
	}
	writeTableFile(t, vla, nil, vcases)
	script := testutil.TestdataPath(t, "tables", "astropy_dump.py")
	out := testutil.AstropyScript(t, script, fixed, vla)
	var all [][]astropyCol
	if err := json.Unmarshal(out, &all); err != nil {
		t.Fatal(err)
	}
	for i, path := range []string{fixed, vla} {
		f, err := Open(path)
		if err != nil {
			t.Fatal(err)
		}
		if len(all[i]) == 0 {
			t.Fatalf("%s: astropy found no numeric columns", path)
		}
		for _, ac := range all[i] {
			checkAstropyColumn(t, filepath.Base(path), f, ac, nil)
		}
		f.Close()
	}
	// astropy also sees the uint64 column as uint64 with exact values.
	var sawW bool
	for _, ac := range all[0] {
		if ac.Name == "W1" {
			sawW = true
			if ac.Kind != "uint" {
				t.Errorf("astropy reads W1 as %s, want uint (TZERO 2^63 convention)", ac.Kind)
			}
		}
	}
	if !sawW {
		t.Errorf("astropy did not report column W1")
	}
}

// fixedWriteCasesNoHeaderScale is fixedWriteCases with ISC's scaling set
// on the column itself (the astropy test passes no header).
func fixedWriteCasesNoHeaderScale() []writeCase {
	cs := fixedWriteCases()
	for i := range cs {
		if cs[i].col.Name == "ISC" {
			cs[i].col.TScal = f64p(1.5)
			cs[i].col.TZero = f64p(-3)
		}
	}
	return cs
}

func TestWriteTableRejects(t *testing.T) {
	cases := []struct {
		name string
		col  ColumnData
		want string
	}{
		{"string longer than TFORM", ColumnData{TForm: "2A", DataString: []string{"hello"}}, "does not fit TFORM 2A"},
		{"TFORM type mismatch", ColumnData{TForm: "1J", DataFloat64: []float64{1.5}}, "does not match"},
		{"vector length not a multiple", ColumnData{TForm: "3E", DataFloat32: []float32{1, 2, 3, 4}}, "not a multiple"},
		{"null without TNULL", ColumnData{DataInt32: []int32{1, 2}, Valid: []bool{true, false}}, "no TNULL"},
		{"value equals TNULL", ColumnData{TNull: i64p(-1), DataInt32: []int32{-1}}, "equals TNULL"},
		{"TNULL on float", ColumnData{TNull: i64p(0), DataFloat32: []float32{1}}, "TNULL is only allowed"},
		{"TNULL out of stored range", ColumnData{TNull: i64p(300), DataUint8: []uint8{1}}, "does not fit"},
		{"TNULL 256 on B", ColumnData{TNull: i64p(256), DataUint8: []uint8{1}}, "does not fit"},
		{"TNULL -1 on B", ColumnData{TNull: i64p(-1), DataUint8: []uint8{1}}, "does not fit"},
		{"TNULL 32768 on I", ColumnData{TNull: i64p(32768), DataInt16: []int16{1}}, "does not fit"},
		{"TNULL -32769 on I", ColumnData{TNull: i64p(-32769), DataInt16: []int16{1}}, "does not fit"},
		{"TNULL 2^31 on J", ColumnData{TNull: i64p(1 << 31), DataInt32: []int32{1}}, "does not fit"},
		{"TNULL -2^31-1 on J", ColumnData{TNull: i64p(-1<<31 - 1), DataInt32: []int32{1}}, "does not fit"},
		{"TSCAL on logical", ColumnData{TScal: f64p(2), DataBool: []bool{true}}, "not allowed on L"},
		{"TZERO with uint16 data", ColumnData{TZero: f64p(1), DataUint16: []uint16{1}}, "cannot be combined"},
		{"two Data fields", ColumnData{DataInt32: []int32{1}, DataInt16: []int16{1}}, "want exactly one"},
		{"non-ASCII string", ColumnData{DataString: []string{"caf\xc3\xa9"}}, "not printable ASCII"},
		{"Dim bigger than repeat", ColumnData{TForm: "2E", Dim: []int64{3}, DataFloat32: []float32{1, 2}}, "more than the repeat"},
		{"null complex", ColumnData{DataComplex64: []complex64{1}, Valid: []bool{false}}, "no null value"},
		{"VLA TForm on fixed data", ColumnData{TForm: "1PE", DataFloat32: []float32{1}}, "does not match"},
		{"VLA row longer than TFORM max", ColumnData{TForm: "1PE(2)", DataVarFloat32: [][]float32{{1}, {1, 2, 3}}}, "max length 2 < longest row 3"},
		{"VarValid row length", ColumnData{DataVarInt32: [][]int32{{1, 2}, {3}}, VarValid: [][]bool{{true, true}, {true, false}}}, "VarValid row 1 has 2 flags, data 1"},
		{"VarValid row count", ColumnData{DataVarInt32: [][]int32{{1, 2}, {3}}, VarValid: [][]bool{{true, true}}}, "VarValid has 1 rows, data 2"},
		{"substring array TForm", ColumnData{TForm: "20A5", DataString: []string{"abc"}}, "substring arrays is not supported"},
		{"VLA null complex", ColumnData{DataVarComplex64: [][]complex64{{1}}, VarValid: [][]bool{{false}}}, "no null value"},
		// A, X (fixed and VLA) have no null value: a false flag must not be
		// dropped silently.
		{"null string", ColumnData{DataString: []string{"a", "b"}, Valid: []bool{true, false}}, "character columns have no null value"},
		{"null bit", ColumnData{DataBits: []bool{true, false}, Valid: []bool{true, false}}, "bit columns have no null value"},
		{"VLA null string", ColumnData{DataVarString: []string{"ab"}, VarValid: [][]bool{{false, true}}}, "character columns have no null value"},
		{"VLA non-ASCII string", ColumnData{DataVarString: []string{"caf\xc3\xa9"}}, "not printable ASCII"},
		{"Valid length", ColumnData{TNull: i64p(-1), DataInt32: []int32{1, 2}, Valid: []bool{true}}, "Valid has 1 flags, data 2"},
		{"TSCAL 0", ColumnData{TScal: f64p(0), DataInt32: []int32{1}}, "TSCAL = 0"},
		{"TSCAL NaN", ColumnData{TScal: f64p(math.NaN()), DataInt32: []int32{1}}, "TSCAL = NaN"},
		{"TSCAL Inf", ColumnData{TScal: f64p(math.Inf(1)), DataInt32: []int32{1}}, "TSCAL = +Inf"},
		{"TZERO NaN", ColumnData{TZero: f64p(math.NaN()), DataInt32: []int32{1}}, "TZERO = NaN"},
		{"TZERO Inf", ColumnData{TZero: f64p(math.Inf(-1)), DataInt32: []int32{1}}, "TZERO = -Inf"},
		{"negative Dim", ColumnData{TForm: "2E", Dim: []int64{-1, 2}, DataFloat32: []float32{1, 2}}, "negative TDIM axis"},
		// One case per remaining guard (each could be removed before
		// without a failing test; Valid on a VLA column then silently
		// dropped the null).
		{"VLA repeat 0", ColumnData{TForm: "0PE", DataVarFloat32: [][]float32{{1}}}, "needs a VLA repeat of 1"},
		{"Valid on a VLA column", ColumnData{DataVarInt32: [][]int32{{1}, {2}, {3}}, Valid: []bool{true, false, true}}, "use VarValid, not Valid"},
		{"TDIM on a VLA column", ColumnData{Dim: []int64{2}, DataVarFloat32: [][]float32{{1, 2}}}, "TDIM on a VLA column"},
		{"VarValid on a fixed column", ColumnData{DataInt32: []int32{1}, VarValid: [][]bool{{true}}}, "use Valid, not VarValid"},
		{"null complex128", ColumnData{DataComplex128: []complex128{1, 2}, Valid: []bool{true, false}}, "no null value"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "w.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			WriteImage(f, nil, []int64{}, []float32{})
			_, err = AppendBinaryTable(f, nil, []ColumnData{c.col})
			wantErr(t, err, c.want)
		})
	}
}

// Table-level writer errors: columns of different lengths, and hdr column
// keywords that contradict the ColumnData (the doc promises an error,
// never a silent choice of one of them).
func TestWriteTableRejectsTable(t *testing.T) {
	hdrWith := func(key string, v any) *header.Header {
		h := header.New()
		h.Set(key, v, "")
		return h
	}
	cases := []struct {
		name string
		hdr  *header.Header
		cols []ColumnData
		want string
	}{
		{"row counts differ", nil, []ColumnData{{DataInt32: []int32{1, 2}}, {DataInt16: []int16{1}}}, "1 rows, earlier columns have 2"},
		{"VLA row count differs", nil, []ColumnData{{DataInt32: []int32{1, 2}}, {DataVarInt16: [][]int16{{1}}}}, "1 rows, earlier columns have 2"},
		{"hdr TSCAL conflicts", hdrWith("TSCAL1", 2.0), []ColumnData{{TScal: f64p(3), DataInt32: []int32{1}}}, "conflicts with column 1"},
		{"hdr TZERO conflicts", hdrWith("TZERO1", 2.0), []ColumnData{{TZero: f64p(3), DataInt32: []int32{1}}}, "conflicts with column 1"},
		{"hdr TNULL conflicts", hdrWith("TNULL1", int64(7)), []ColumnData{{TNull: i64p(8), DataInt32: []int32{1}}}, "conflicts with column 1"},
		{"hdr TTYPE conflicts", hdrWith("TTYPE1", "A"), []ColumnData{{Name: "B", DataInt32: []int32{1}}}, "conflicts with column 1"},
		{"hdr TFORM conflicts", hdrWith("TFORM1", "1E"), []ColumnData{{TForm: "1J", DataInt32: []int32{1}}}, "conflicts with column 1"},
		{"hdr key beyond the columns", hdrWith("TSCAL2", 2.0), []ColumnData{{DataInt32: []int32{1}}}, "the table has 1 columns"},
		{"no columns", nil, nil, "no columns"},
		{"1000 columns", nil, make([]ColumnData, 1000), "1000 columns, FITS allows 999"},
		{"hdr TDIM unparsable", hdrWith("TDIM1", "(2,x)"), []ColumnData{{TForm: "2E", DataFloat32: []float32{1, 2}}}, "hdr TDIM1"},
		{"hdr TSCAL not a number", hdrWith("TSCAL1", "two"), []ColumnData{{DataInt32: []int32{1}}}, "hdr TSCAL1"},
		{"hdr TZERO not a number", hdrWith("TZERO1", true), []ColumnData{{DataInt32: []int32{1}}}, "hdr TZERO1"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			f, err := Create(filepath.Join(t.TempDir(), "w.fits"))
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			WriteImage(f, nil, []int64{}, []float32{})
			_, err = AppendBinaryTable(f, c.hdr, c.cols)
			wantErr(t, err, c.want)
		})
	}
}

// A 1P descriptor has 32-bit signed fields: a heap past 2^31-1 bytes (or a
// longer row) is an error, never a wrapped descriptor. Building a 2 GiB
// heap in a test is not practical, so the descriptor writer that every VLA
// row goes through is called directly at the boundary.
func TestPutDescriptorRange(t *testing.T) {
	cell := make([]byte, 16)
	if err := putDescriptor(cell, false, 'E', 3, math.MaxInt32-12, math.MaxInt32); err != nil {
		t.Fatalf("heap ending at 2^31-1: %v", err)
	}
	if got := [2]int32{int32(binaryBE32(cell)), int32(binaryBE32(cell[4:]))}; got != [2]int32{3, math.MaxInt32 - 12} {
		t.Errorf("P descriptor %v", got)
	}
	err := putDescriptor(cell, false, 'E', 3, math.MaxInt32-11, math.MaxInt32+1)
	wantErr(t, err, "exceeds the 1P descriptor range")
	err = putDescriptor(cell, false, 'X', math.MaxInt32+1, 0, 1<<28)
	wantErr(t, err, "exceeds the 1P descriptor range")
	if err := putDescriptor(cell, true, 'E', 3, 1<<40, 1<<40+12); err != nil {
		t.Fatalf("Q descriptor: %v", err)
	}
	if n, off := binaryBE64(cell), binaryBE64(cell[8:]); n != 3 || off != 1<<40 {
		t.Errorf("Q descriptor (%d, %d)", n, off)
	}
}

func binaryBE32(b []byte) uint32 {
	return uint32(b[0])<<24 | uint32(b[1])<<16 | uint32(b[2])<<8 | uint32(b[3])
}
func binaryBE64(b []byte) uint64 { return uint64(binaryBE32(b))<<32 | uint64(binaryBE32(b[4:])) }

// Complex cells keep ±Inf, -0, NaN and subnormals bit-for-bit. cfitsio is
// no truth here: its null-checking read flags Inf/NaN as null without
// returning a value (the cfref harness prints complex values without a null
// flag) and flushes subnormals to 0. astropy reads the bytes as IEEE
// values, so it is the reader: astropy must see exactly what we wrote, and
// our reader must agree with astropy, for C and M, fixed and VLA.
func TestComplexSpecialValuesAstropy(t *testing.T) {
	inf := math.Inf(1)
	sub32 := float64(math.SmallestNonzeroFloat32)
	sub64 := math.SmallestNonzeroFloat64
	negz := math.Copysign(0, -1)
	nan := math.NaN()
	c64 := []complex64{complex(float32(inf), float32(-inf)), complex(float32(negz), float32(sub32)),
		complex(float32(-sub32), float32(nan)), complex(float32(0x1p-140), float32(math.MaxFloat32))}
	c128 := []complex128{complex(inf, -inf), complex(negz, sub64), complex(-sub64, nan), complex(0x1p-1030, math.MaxFloat64)}
	cols := []ColumnData{
		{Name: "C1", DataComplex64: c64},
		{Name: "M1", DataComplex128: c128},
		{Name: "VC", DataVarComplex64: [][]complex64{c64[:2], {}, c64[2:], c64[3:]}},
		{Name: "VM", DataVarComplex128: [][]complex128{c128[:1], c128[1:3], {}, c128[3:]}},
	}
	want := map[string][][]complex128{
		"C1": {{complex128(c64[0])}, {complex128(c64[1])}, {complex128(c64[2])}, {complex128(c64[3])}},
		"M1": {{c128[0]}, {c128[1]}, {c128[2]}, {c128[3]}},
		"VC": {{complex128(c64[0]), complex128(c64[1])}, {}, {complex128(c64[2]), complex128(c64[3])}, {complex128(c64[3])}},
		"VM": {{c128[0]}, {c128[1], c128[2]}, {}, {c128[3]}},
	}
	path := filepath.Join(t.TempDir(), "cplx.fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []float32{}); err != nil {
		t.Fatal(err)
	}
	if _, err := AppendBinaryTable(f, nil, cols); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}
	out := testutil.AstropyScript(t, testutil.TestdataPath(t, "tables", "astropy_dump.py"), path)
	var all [][]astropyCol
	if err := json.Unmarshal(out, &all); err != nil {
		t.Fatal(err)
	}
	r, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer r.Close()
	seen := 0
	for _, ac := range all[0] {
		w, ok := want[ac.Name]
		if !ok {
			continue
		}
		seen++
		if ac.Kind != "complex" || len(ac.Rows) != len(w) {
			t.Fatalf("%s: astropy kind %s, %d rows", ac.Name, ac.Kind, len(ac.Rows))
		}
		for row, raw := range ac.Rows {
			got := astropyRow(t, raw, "complex")
			if len(got) != len(w[row]) {
				t.Errorf("%s row %d: astropy %d elements, wrote %d", ac.Name, row, len(got), len(w[row]))
				continue
			}
			for k, z := range w[row] {
				if s := fbits(real(z)) + "," + fbits(imag(z)); got[k] != s {
					t.Errorf("%s row %d elem %d: astropy %s, wrote %v (%s)", ac.Name, row, k, got[k], z, s)
				}
			}
		}
		checkAstropyColumn(t, "cplx.fits", r, ac, nil)
	}
	if seen != len(want) {
		t.Errorf("astropy reported %d of the %d complex columns", seen, len(want))
	}
}
