package fits

import (
	"bytes"
	"encoding/binary"
	"errors"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Loud-failure tests: everything the table code cannot represent must be
// an error, never zeros, a wrapped integer or a silently dropped value.

func tcard(s string) string { return s + strings.Repeat(" ", 80-len(s)) }

// rawTableFile builds primary + one extension from cards and data bytes.
func rawTableFile(t *testing.T, cards []string, data []byte) string {
	t.Helper()
	var b bytes.Buffer
	prim := []string{"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0", "EXTEND  =                    T", "END"}
	for _, c := range prim {
		b.WriteString(tcard(c))
	}
	for b.Len()%2880 != 0 {
		b.WriteByte(' ')
	}
	for _, c := range append(cards, "END") {
		b.WriteString(tcard(c))
	}
	for b.Len()%2880 != 0 {
		b.WriteByte(' ')
	}
	b.Write(data)
	for b.Len()%2880 != 0 {
		b.WriteByte(0)
	}
	p := filepath.Join(t.TempDir(), "t.fits")
	if err := os.WriteFile(p, b.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

func openTable(t *testing.T, path string) HDU {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	t.Cleanup(func() { f.Close() })
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	return h
}

func binHeader(naxis1, naxis2, pcount int, extra ...string) []string {
	c := []string{
		"XTENSION= 'BINTABLE'",
		"BITPIX  =                    8",
		"NAXIS   =                    2",
		"NAXIS1  = " + padNum(naxis1),
		"NAXIS2  = " + padNum(naxis2),
		"PCOUNT  = " + padNum(pcount),
		"GCOUNT  =                    1",
	}
	return append(c, extra...)
}

func padNum(n int) string {
	s := itoa(n)
	return strings.Repeat(" ", 20-len(s)) + s
}

func wantErr(t *testing.T, err error, sub string) {
	t.Helper()
	if err == nil {
		t.Fatalf("got no error, want one containing %q", sub)
	}
	if !strings.Contains(err.Error(), sub) {
		t.Fatalf("error %q does not contain %q", err, sub)
	}
}

func TestTableColumnsRejectBadKeywords(t *testing.T) {
	cases := []struct {
		name  string
		cards []string
		want  string
	}{
		{"TDIM bigger than repeat", []string{"TFIELDS =                    1", "TFORM1  = '2J'", "TDIM1   = '(3)'"}, "more than the repeat count"},
		{"malformed TDIM", []string{"TFIELDS =                    1", "TFORM1  = '2J'", "TDIM1   = '(a)'"}, "TDIM1"},
		{"malformed TFORM", []string{"TFIELDS =                    1", "TFORM1  = '2J5'"}, "TFORM1"},
		{"TSCAL zero", []string{"TFIELDS =                    1", "TFORM1  = '2J'", "TSCAL1  =                  0.0"}, "TSCAL1 = 0"},
		{"non-integer TNULL", []string{"TFIELDS =                    1", "TFORM1  = '2J'", "TNULL1  =                  1.5"}, "TNULL1"},
		{"widths != NAXIS1", []string{"TFIELDS =                    1", "TFORM1  = '1J'"}, "!= NAXIS1"},
		{"substring width > repeat", []string{"TFIELDS =                    1", "TFORM1  = '8A10'"}, "substring width exceeds repeat"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			p := rawTableFile(t, binHeader(8, 1, 0, c.cards...), make([]byte, 8))
			_, err := openTable(t, p).(*BinaryTableHDU).Columns()
			wantErr(t, err, c.want)
		})
	}
}

// Column layouts whose byte arithmetic overflows int64 must fail in
// Columns, not pass the NAXIS1 check after wrapping and then panic in a
// reader (8*(2^61+1) wraps to 8; the TDIM product 2*(2^62+1)*2 wraps to
// 4; four column widths adding up to 2^64+8 wrap to 8).
func TestColumnLayoutOverflow(t *testing.T) {
	for _, c := range []struct {
		name   string
		naxis1 int
		cards  []string
		want   string
	}{
		{"repeat*size", 8, []string{"TFIELDS =                    1", "TFORM1  = '2305843009213693953D'"}, `TFORM1 "2305843009213693953D": 2305843009213693953 elements of 8 bytes overflow`},
		{"TDIM product", 4, []string{"TFIELDS =                    1", "TFORM1  = '4A'", "TDIM1   = '(2,4611686018427387905,2)'"}, "TDIM1 (2,4611686018427387905,2): the element count overflows"},
		{"VLA TDIM product", 8, []string{"TFIELDS =                    1", "TFORM1  = '1PB(4)'", "TDIM1   = '(2,4611686018427387905,2)'"}, "TDIM1 (2,4611686018427387905,2): the element count overflows"},
		{"row width", 8, []string{"TFIELDS =                    4", "TFORM1  = '4611686018427387904B'", "TFORM2  = '4611686018427387904B'", "TFORM3  = '9223372036854775807B'", "TFORM4  = '9B'"}, "the column widths overflow"},
	} {
		t.Run(c.name, func(t *testing.T) {
			p := rawTableFile(t, binHeader(c.naxis1, 1, 0, c.cards...), make([]byte, c.naxis1))
			bt := openTable(t, p).(*BinaryTableHDU)
			_, err := bt.Columns()
			wantErr(t, err, c.want)
			if _, err := ReadVectorColumn[float64](bt, 1); err == nil {
				t.Errorf("ReadVectorColumn: no error")
			}
			if _, err := ReadStringArrayColumn(bt, 1); err == nil {
				t.Errorf("ReadStringArrayColumn: no error")
			}
		})
	}
}

func TestVarColumnRejectsBadHeap(t *testing.T) {
	// One 1PJ column, 1 row: descriptor (n=2, off=0), heap 8 bytes.
	row := []byte{0, 0, 0, 2, 0, 0, 0, 0}
	heap := []byte{0, 0, 0, 1, 0, 0, 0, 2}
	t.Run("ok", func(t *testing.T) {
		p := rawTableFile(t, binHeader(8, 1, 8, "TFIELDS =                    1", "TFORM1  = '1PJ(2)'"), append(append([]byte{}, row...), heap...))
		v, err := ReadVarColumn[int32](openTable(t, p).(*BinaryTableHDU), 1)
		if err != nil || len(v.At(0)) != 2 || v.At(0)[1] != 2 {
			t.Fatalf("v=%v err=%v", v, err)
		}
	})
	t.Run("descriptor beyond PCOUNT", func(t *testing.T) {
		p := rawTableFile(t, binHeader(8, 1, 4, "TFIELDS =                    1", "TFORM1  = '1PJ(2)'"), append(append([]byte{}, row...), heap...))
		_, err := ReadVarColumn[int32](openTable(t, p).(*BinaryTableHDU), 1)
		wantErr(t, err, "outside the heap")
	})
	t.Run("THEAP inside the rows", func(t *testing.T) {
		p := rawTableFile(t, binHeader(8, 1, 8, "TFIELDS =                    1", "TFORM1  = '1PJ(2)'", "THEAP   =                    4"), append(append([]byte{}, row...), heap...))
		_, err := ReadVarColumn[int32](openTable(t, p).(*BinaryTableHDU), 1)
		wantErr(t, err, "THEAP = 4 outside")
	})
	// 1Q descriptors are 64-bit: element counts and offsets whose byte
	// range overflows int64 must be errors, not panics (n*4 and off+bytes
	// wrap around).
	for _, d := range []struct {
		name   string
		n, off uint64
	}{
		{"Q count*size overflows", 1 << 62, 0},
		{"Q count*size wraps negative", 1 << 61, 0},
		{"Q offset+size overflows", 1, 1<<63 - 2},
		{"Q huge count", 1<<63 - 1, 0},
	} {
		t.Run(d.name, func(t *testing.T) {
			cell := make([]byte, 16)
			for i := 0; i < 8; i++ {
				cell[i] = byte(d.n >> (56 - 8*i))
				cell[8+i] = byte(d.off >> (56 - 8*i))
			}
			p := rawTableFile(t, binHeader(16, 1, 8, "TFIELDS =                    1", "TFORM1  = '1QJ(2)'"), append(cell, heap...))
			_, err := ReadVarColumn[int32](openTable(t, p).(*BinaryTableHDU), 1)
			wantErr(t, err, "outside the heap")
		})
	}
	t.Run("negative descriptor", func(t *testing.T) {
		bad := []byte{0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0}
		p := rawTableFile(t, binHeader(8, 1, 8, "TFIELDS =                    1", "TFORM1  = '1PJ(2)'"), append(bad, heap...))
		_, err := ReadVarColumn[int32](openTable(t, p).(*BinaryTableHDU), 1)
		wantErr(t, err, "invalid VLA descriptor")
	})
}

// A repeat count of 0 on a P/Q column ("0PJ", "0QD") occupies no bytes
// in the row (FITS 4.0 7.3.1: the field width is the repeat count times
// the element size) and every row is empty. No reference tool can be the
// truth here: cfitsio 4.7.0 sizes such columns 8/16 bytes when it writes
// NAXIS1 but 0 when it reads, and its harness crashes on them; astropy
// 8.0.2 cannot build the record array; fitsverify reads descriptors from
// the next column. The file below is laid out by the standard by hand.
func TestZeroRepeatVLA(t *testing.T) {
	rows := []byte{0, 0, 0, 1, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 2, 0xff, 0xff, 0xff, 0xfe, 0, 0, 0, 3, 0xff, 0xff, 0xff, 0xfd}
	p := rawTableFile(t, binHeader(8, 3, 0, "TFIELDS =                    4",
		"TFORM1  = '1J'", "TFORM2  = '0PJ'", "TFORM3  = '0QD'", "TFORM4  = '1J'"), rows)
	h := openTable(t, p).(*BinaryTableHDU)
	for _, col := range []int{2, 3} {
		v, err := ReadVarColumn[float64](h, col)
		if err != nil {
			t.Fatalf("col %d: %v", col, err)
		}
		if v.Len() != 3 {
			t.Fatalf("col %d: %d rows, want 3", col, v.Len())
		}
		for r := 0; r < 3; r++ {
			if len(v.At(r)) != 0 {
				t.Errorf("col %d row %d: %v, want empty", col, r, v.At(r))
			}
		}
	}
	a, err := ReadColumn[int32](h, 1)
	if err != nil || a[0] != 1 || a[1] != 2 || a[2] != 3 {
		t.Errorf("col 1: %v %v", a, err)
	}
	b, err := ReadColumn[int32](h, 4)
	if err != nil || b[0] != -1 || b[1] != -2 || b[2] != -3 {
		t.Errorf("col 4: %v %v (the 0P/0Q columns must take no bytes)", b, err)
	}
}

func TestASCIIRejectsBadFields(t *testing.T) {
	hdr := []string{
		"XTENSION= 'TABLE   '", "BITPIX  =                    8", "NAXIS   =                    2",
		"NAXIS1  =                    6", "NAXIS2  =                    1", "PCOUNT  =                    0",
		"GCOUNT  =                    1", "TFIELDS =                    1", "TBCOL1  =                    1",
	}
	t.Run("garbage in numeric field", func(t *testing.T) {
		p := rawTableFile(t, append(hdr, "TFORM1  = 'I6'"), []byte("  12x "))
		_, err := ReadASCIIColumn[int64](openTable(t, p).(*ASCIITableHDU), 1)
		wantErr(t, err, "cannot read a number")
	})
	t.Run("lowercase exponent", func(t *testing.T) {
		p := rawTableFile(t, append(hdr, "TFORM1  = 'E6.1'"), []byte(" 1.5e2"))
		_, err := ReadASCIIColumn[float64](openTable(t, p).(*ASCIITableHDU), 1)
		wantErr(t, err, "cannot read a number")
	})
	t.Run("field beyond NAXIS1", func(t *testing.T) {
		p := rawTableFile(t, append(hdr, "TFORM1  = 'I8'"), []byte("     12 "))
		_, err := openTable(t, p).(*ASCIITableHDU).Columns()
		wantErr(t, err, "does not fit NAXIS1")
	})
	t.Run("TBCOL 0", func(t *testing.T) {
		// TBCOL counts from 1 (FITS 4.0 7.2.1): 0 would put the field one
		// byte before the row (the previous row's last byte, or before the
		// data for row 0).
		h := append(append([]string{}, hdr...), "TFORM1  = 'I3'")
		h[8] = "TBCOL1  =                    0"
		p := rawTableFile(t, h, []byte("  12  "))
		_, err := openTable(t, p).(*ASCIITableHDU).Columns()
		wantErr(t, err, "TBCOL1 = 0 with width 3 does not fit NAXIS1 = 6")
	})
	t.Run("TBCOL near MaxInt64", func(t *testing.T) {
		// TBCOL-1+width must not wrap around int64 past the NAXIS1 check.
		h := append(append([]string{}, hdr...), "TFORM1  = 'I3'")
		h[8] = "TBCOL1  =  9223372036854775807"
		p := rawTableFile(t, h, []byte("  12  "))
		a := openTable(t, p).(*ASCIITableHDU)
		_, err := a.Columns()
		wantErr(t, err, "TBCOL1 = 9223372036854775807 with width 3 does not fit NAXIS1 = 6")
		if _, err := ReadASCIIColumn[int64](a, 1); err == nil {
			t.Errorf("ReadASCIIColumn: no error")
		}
	})
	t.Run("width near MaxInt64", func(t *testing.T) {
		h := append(append([]string{}, hdr...), "TFORM1  = 'I9223372036854775807'")
		h[8] = "TBCOL1  =                    2"
		p := rawTableFile(t, h, []byte("  12  "))
		_, err := openTable(t, p).(*ASCIITableHDU).Columns()
		wantErr(t, err, "TBCOL1 = 2 with width 9223372036854775807 does not fit NAXIS1 = 6")
	})
	t.Run("exponent out of range", func(t *testing.T) {
		// 7 significant exponent digits: refused, not parsed as 10^0 or
		// silently turned into Inf.
		h := append(append([]string{}, hdr...), "TFORM1  = 'E12.1'")
		h[3] = "NAXIS1  =                   12"
		p := rawTableFile(t, h, []byte("1.5E+1234567"))
		_, err := ReadASCIIColumn[float64](openTable(t, p).(*ASCIITableHDU), 1)
		wantErr(t, err, "exponent out of range")
	})
	t.Run("exponent overflow and underflow", func(t *testing.T) {
		// Within the 6-digit exponent limit, a value beyond the float64
		// range is +-Inf (IEEE overflow of the exact decimal) and a tiny
		// one 0; neither is an integer.
		h := append(append([]string{}, hdr...), "TFORM1  = 'E12.1'")
		h[3] = "NAXIS1  =                   12"
		for field, want := range map[string]float64{"   1.5E+400 ": math.Inf(1), "  -1.5E+400 ": math.Inf(-1), "   1.5E-400 ": 0} {
			p := rawTableFile(t, h, []byte(field))
			a := openTable(t, p).(*ASCIITableHDU)
			v, err := ReadASCIIColumn[float64](a, 1)
			if err != nil || v[0] != want {
				t.Errorf("%q: %v %v, want %v", field, v, err, want)
			}
			if _, err := ReadASCIIColumn[int64](a, 1); !errors.Is(err, ErrTypeMismatch) {
				t.Errorf("%q as int64: err = %v, want ErrTypeMismatch", field, err)
			}
		}
	})
	t.Run("null string read unmasked", func(t *testing.T) {
		p := rawTableFile(t, append(hdr, "TFORM1  = 'A6'", "TNULL1  = 'N/A'"), []byte("N/A   "))
		h := openTable(t, p).(*ASCIITableHDU)
		_, err := ReadASCIIColumn[string](h, 1)
		if !errors.Is(err, ErrNullValue) {
			t.Fatalf("err = %v, want ErrNullValue", err)
		}
		v, ok, err := ReadASCIIColumnMasked[string](h, 1)
		if err != nil || ok[0] || v[0] != "" {
			t.Fatalf("masked: %q %v %v", v, ok, err)
		}
	})
}

func TestTableReadExactnessErrors(t *testing.T) {
	f, err := Open(testutil.TestdataPath(t, "tables", "cf_bintable_fixed.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	tb := mustBinTable(t, f, 1)
	col := func(name string) int {
		i := tb.ColumnIndex(name)
		if i == 0 {
			t.Fatalf("no column %s", name)
		}
		return i
	}
	// uint16 values (TZERO 32768) do not fit int16: error, not a wrap.
	if _, err := ReadColumn[int16](tb, col("U1")); !errors.Is(err, ErrTypeMismatch) {
		t.Errorf("U1 as int16: err = %v, want ErrTypeMismatch", err)
	}
	if v, err := ReadColumn[uint16](tb, col("U1")); err != nil || v[5] != 65535 || v[3] != 32768 {
		t.Errorf("U1 as uint16: %v %v", v, err)
	}
	// Negative values never become unsigned (no wrap).
	for name, read := range map[string]func(int) error{
		"I1 as uint16": func(c int) error { _, err := ReadColumn[uint16](tb, c); return err },
		"I1 as uint64": func(c int) error { _, err := ReadColumn[uint64](tb, c); return err },
		"J1 as uint32": func(c int) error { _, err := ReadColumn[uint32](tb, c); return err },
		"K1 as uint64": func(c int) error { _, err := ReadColumn[uint64](tb, c); return err },
		"S1 as uint8":  func(c int) error { _, err := ReadColumn[uint8](tb, c); return err },
	} {
		if err := read(col(strings.Fields(name)[0])); !errors.Is(err, ErrTypeMismatch) {
			t.Errorf("%s: err = %v, want ErrTypeMismatch (column holds negative values)", name, err)
		}
	}
	// TSCAL 0.5: 10.5 is not an integer.
	if _, err := ReadColumn[int32](tb, col("JSC")); !errors.Is(err, ErrTypeMismatch) {
		t.Errorf("JSC as int32: err = %v, want ErrTypeMismatch", err)
	}
	// int64 extremes and 2^62+1 are exact.
	k, err := ReadColumn[int64](tb, col("K1"))
	if err != nil || k[0] != math.MinInt64 || k[2] != 1<<62+1 || k[3] != 4295806720896741504 || k[5] != math.MaxInt64 {
		t.Errorf("K1: %v %v", k, err)
	}
	// TNULL: integer targets need the masked reader; floats get NaN.
	if _, err := ReadColumn[int16](tb, col("INUL")); !errors.Is(err, ErrNullValue) {
		t.Errorf("INUL as int16: err = %v, want ErrNullValue", err)
	}
	fv, err := ReadColumn[float64](tb, col("INUL"))
	if err != nil || !math.IsNaN(fv[1]) || !math.IsNaN(fv[4]) || fv[5] != 32767 {
		t.Errorf("INUL as float64: %v %v", fv, err)
	}
	if _, err := ReadColumn[bool](tb, col("L1")); !errors.Is(err, ErrNullValue) {
		t.Errorf("L1 with undefined logicals as bool: err = %v, want ErrNullValue", err)
	}
	// Readers for other types refuse instead of returning zeros.
	for name, read := range map[string]func(int) error{
		"X13": func(c int) error { _, err := ReadVectorColumn[float64](tb, c); return err },
		"C2":  func(c int) error { _, err := ReadVectorColumn[float64](tb, c); return err },
		"M1":  func(c int) error { _, err := ReadColumn[float64](tb, c); return err },
		"L3":  func(c int) error { _, err := ReadVectorColumn[uint8](tb, c); return err },
		"A12": func(c int) error { _, err := ReadColumn[float64](tb, c); return err },
		"E3":  func(c int) error { _, err := ReadColumn[float64](tb, c); return err },
	} {
		if err := read(col(name)); err == nil {
			t.Errorf("%s: wrong reader returned no error", name)
		}
	}

	fv2, err := Open(testutil.TestdataPath(t, "tables", "cf_vla_p.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer fv2.Close()
	vt := mustBinTable(t, fv2, 1)
	for _, name := range []string{"VL", "VX", "VA", "VC", "VM"} {
		if _, err := ReadVarColumn[float64](vt, vt.ColumnIndex(name)); err == nil {
			t.Errorf("ReadVarColumn on %s returned no error", name)
		}
	}
	// VLA TNULL (77) on integer targets.
	if _, err := ReadVarColumn[int16](vt, vt.ColumnIndex("VI")); !errors.Is(err, ErrNullValue) {
		t.Errorf("VI as int16: err = %v, want ErrNullValue", err)
	}
	// 2^62+3 in a K VLA is exact as int64 but not representable in int32.
	kv, err := ReadVarColumn[int64](vt, vt.ColumnIndex("VK"))
	if err == nil {
		t.Errorf("VK has TNULL rows; want ErrNullValue, got %v", kv)
	}
	kvm, _, err := ReadVarColumnMasked[int64](vt, vt.ColumnIndex("VK"))
	if err != nil || kvm.At(3)[3] != 1<<62+3 {
		t.Errorf("VK row 3: %v %v", kvm.At(3), err)
	}
}

func TestComplexScalingRejected(t *testing.T) {
	// A literal whose nearest double is 1 still declares scaling.
	for _, key := range []string{"TSCAL1  =                  2.0", "TZERO1  =                  1.0",
		"TSCAL1  = 1.00000000000000000001"} {
		p := rawTableFile(t, binHeader(8, 1, 0, "TFIELDS =                    1", "TFORM1  = '1C'", key), make([]byte, 8))
		_, err := ReadComplexColumn[complex64](openTable(t, p).(*BinaryTableHDU), 1)
		wantErr(t, err, "TSCAL/TZERO on complex columns is not supported")
		if key == "TZERO1  =                  1.0" {
			// TSCAL = 1.0 / TZERO = 0.0 written out are no scaling.
			p := rawTableFile(t, binHeader(8, 1, 0, "TFIELDS =                    1", "TFORM1  = '1C'", "TSCAL1  =                  1.0", "TZERO1  =                  0.0"), make([]byte, 8))
			if _, err := ReadComplexColumn[complex64](openTable(t, p).(*BinaryTableHDU), 1); err != nil {
				t.Errorf("TSCAL 1.0 / TZERO 0.0 on 1C: %v", err)
			}
		}
		// The same for a complex VLA: one row of one element.
		row := []byte{0, 0, 0, 1, 0, 0, 0, 0}
		p = rawTableFile(t, binHeader(8, 1, 8, "TFIELDS =                    1", "TFORM1  = '1PC(1)'", key), append(row, make([]byte, 8)...))
		_, err = ReadVarComplexColumn[complex64](openTable(t, p).(*BinaryTableHDU), 1)
		wantErr(t, err, "TSCAL/TZERO on complex columns is not supported")
	}
}

// Why TSCAL/TZERO on complex columns is refused: FITS 4.0 7.3.2 scales
// the real part (TSCAL*re + TZERO) and only multiplies the imaginary part
// (TSCAL*im), but neither reference does that. The committed goldens of
// testdata/tables/refused/cf_bintable_cplx_scaled.fits (written by
// cfitsio with scaling off) record cfitsio 4.7.0 adding TZERO to the
// imaginary part too and astropy 8.0.2 dropping the imaginary part (it
// returns a float column). The standard's value is computed here from the
// stored IEEE values decoded from the sample bytes; the test pins that
// both references still disagree with it and that we refuse the columns.
func TestComplexScalingReferencesDisagree(t *testing.T) {
	path := testutil.TestdataPath(t, "tables", "refused", "cf_bintable_cplx_scaled.fits")
	cf := loadGoldenFile(t, path+".cfref.json")
	ag := loadAstropyGolden(t, path+".astropy.json")
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	gh := cf.HDUs[1]
	f, err := Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	bt := mustBinTable(t, f, 1)
	off := int64(0)
	for i, gc := range gh.Table.Columns {
		scale, zero := float64(gc.TScal), float64(gc.TZero)
		if scale == 1 || zero == 0 {
			t.Fatalf("%s: TSCAL %v / TZERO %v, the sample must scale both", *gc.TType, scale, zero)
		}
		es := typecodeBytes(gc.Typecode) / 2
		ac := ag.Columns[i]
		if ac.Kind != "float" || ac.Error != "" {
			t.Errorf("%s: astropy kind %q (error %q); it used to drop the imaginary part (kind float)", *gc.TType, ac.Kind, ac.Error)
		}
		for r, cell := range gc.RowCells(t) {
			p := raw[gh.DataStart+int64(r)*gh.Table.RowLen+off:]
			var re, im float64
			if es == 4 {
				re = float64(math.Float32frombits(binary.BigEndian.Uint32(p)))
				im = float64(math.Float32frombits(binary.BigEndian.Uint32(p[4:])))
			} else {
				re = math.Float64frombits(binary.BigEndian.Uint64(p))
				im = math.Float64frombits(binary.BigEndian.Uint64(p[8:]))
			}
			stdRe, stdIm := re*scale+zero, im*scale
			c := cell.Cplx[0]
			if float64(c[0]) != stdRe || float64(c[1]) != im*scale+zero || float64(c[1]) == stdIm {
				t.Errorf("%s row %d: cfitsio (%v, %v), standard (%v, %v): cfitsio no longer adds TZERO to the imaginary part", *gc.TType, r, c[0], c[1], stdRe, stdIm)
			}
			if a := astropyRow(t, ac.Rows[r], ac.Kind); len(a) != 1 || a[0] != fbits(stdRe) {
				t.Errorf("%s row %d: astropy %v, want the real part %v alone", *gc.TType, r, a, stdRe)
			}
		}
		off += 2 * es
	}
	if _, err := ReadComplexColumn[complex64](bt, 1); err == nil || !strings.Contains(err.Error(), "TSCAL/TZERO on complex columns is not supported") {
		t.Errorf("1C: err = %v", err)
	}
	if _, err := ReadComplexColumn[complex128](bt, 2); err == nil || !strings.Contains(err.Error(), "TSCAL/TZERO on complex columns is not supported") {
		t.Errorf("1M: err = %v", err)
	}
}
