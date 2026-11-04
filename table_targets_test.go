package fits

import (
	"encoding/binary"
	"errors"
	"fmt"
	"math"
	"math/big"
	"os"
	"reflect"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Every Go target type, on every numeric column of every sample.
//
// For each numeric cell the test builds the exact physical value from
// the cfitsio golden, the sample bytes and the FITS definitions alone,
// never from the Go reader. The physical value is stored*TSCAL + TZERO
// with both literals read as exact decimals from the header card text
// (FITS 4.0 7.3.2 / 7.2.5), evaluated in math/big:
//
//   - B/I/J/K with TSCAL 1 and an integer TZERO literal of any size (the
//     uint64 convention TZERO = 2^63 on K is one case of it): the exact
//     integer stored + TZERO, the stored integer being cfitsio's read
//     with scaling off ("raw"), or its TLONGLONG read when the column is
//     not scaled. Float targets get that integer correctly rounded.
//   - Any other scaled column (TSCAL != 1 or a fractional TZERO), integer
//     or float: the exact rational stored*TSCAL + TZERO decides what an
//     integer target gets (the integer, or ErrTypeMismatch when it is not
//     one or out of range). The stored value is cfitsio's raw read for
//     B/I/J/K and the IEEE value decoded from the sample bytes (at the
//     data_start and column widths cfitsio reports) for E/D. Float
//     targets get cfitsio's TDOUBLE read: the documented float64
//     arithmetic.
//   - Unscaled E/D: cfitsio's TDOUBLE read. Where cfitsio's
//     null-checking read flushes a value (+-Inf flagged null, -0 and
//     subnormals read as +0; see cfitsioFlushed) the IEEE value from our
//     float64 read is used: those cells are proven bit-exact against the
//     committed astropy goldens in TestTablesMatchAstropy, which also
//     requires a cell we call null to be NaN (or TNULL) in astropy.
//   - ASCII I/F/E/D fields: the field text (cfitsio's raw row bytes) read
//     as an exact decimal per FITS 4.0 7.2.5 (blanks ignored, implied
//     decimal point of Fw.d/Ew.d/Dw.d, E or D exponent, a field with no
//     digits is 0 as in cfitsio), null when it equals TNULLn padded to the
//     field width (7.2.2); then TSCAL/TZERO exactly as above. Float
//     targets on the scaled path get the field's float64 scaled in
//     float64 (documented). cfitsio's TDOUBLE read must agree with that
//     to within 2 ulp; it is not correctly rounded for long fields, which
//     is pinned per sample.
//
// The expected result for a target T then follows from the documented
// rules (bin_table_read.go): an integer T gets the value when every
// non-null cell of the column is an integer in T's range, otherwise the
// read is an ErrTypeMismatch error; float64 gets the value above,
// float32 the float32 conversion of that float64; nulls are valid=false
// in the masked readers and NaN in the unmasked float reads.
// cf_bintable_ranges.fits puts values on and one step past the range of
// every integer type, through all three arithmetic paths;
// cf_bintable_scaled_exact.fits and cf_ascii_scaled_exact.fits hold
// values whose exact physical value is not the float64 arithmetic result.

type exactCell struct {
	null bool
	// exact integer value, or nil when the value is v (a float64 that
	// may or may not be integral) or r
	i *big.Int
	v float64
	// r: the exact value on the float64 scaling path, where float
	// targets get v (the float64 arithmetic) and integer targets the
	// exact value; nil for a non-finite stored value.
	r *big.Rat
	// noInt: the exact value is not an integer even though v may be
	// (an ASCII decimal with a fraction that rounds to an integral double)
	noInt bool
	// x: the exact value on an exact path whose float targets get it
	// correctly rounded (v is then its float64; float32 targets get the
	// float32 nearest x, not the float32 of v). nil elsewhere.
	x *big.Rat
}

// integral returns the cell's value as an exact integer, if it is one.
func (c exactCell) integral() (*big.Int, bool) {
	if c.i != nil {
		return c.i, true
	}
	if c.r != nil {
		if !c.r.IsInt() {
			return nil, false
		}
		return new(big.Int).Set(c.r.Num()), true
	}
	if c.noInt || math.IsInf(c.v, 0) || math.IsNaN(c.v) || c.v != math.Trunc(c.v) {
		return nil, false
	}
	b, _ := new(big.Float).SetFloat64(c.v).Int(nil)
	return b, true
}

// f64 is the value a float64 target must receive.
func (c exactCell) f64() float64 {
	if c.i != nil {
		f, _ := new(big.Float).SetInt(c.i).Float64()
		return f
	}
	return c.v
}

// f32 is the value a float32 target must receive: the float32 nearest
// the exact value on the exact paths (an exact integer, or x), otherwise
// the float32 conversion of the float64 the column's arithmetic gives.
func (c exactCell) f32() float32 {
	switch {
	case c.i != nil:
		f, _ := new(big.Float).SetInt(c.i).Float32()
		return f
	case c.x != nil:
		f, _ := c.x.Float32()
		return f
	}
	return float32(c.v)
}

func (c exactCell) String() string {
	switch {
	case c.null:
		return "null"
	case c.i != nil:
		return c.i.String()
	case c.r != nil:
		return c.r.RatString() + " (float64 " + fmt.Sprint(c.v) + ")"
	}
	return fmt.Sprint(c.v)
}

// targetReaders reads one column as T: masked (values + validity, flat)
// and unmasked.
type targetReaders[T Numeric] struct {
	masked   func() ([]T, []bool, error)
	unmasked func() ([]T, error)
}

func kindRange(k reflect.Kind) (lo, hi *big.Int) {
	switch k {
	case reflect.Int8:
		return big.NewInt(math.MinInt8), big.NewInt(math.MaxInt8)
	case reflect.Int16:
		return big.NewInt(math.MinInt16), big.NewInt(math.MaxInt16)
	case reflect.Int32:
		return big.NewInt(math.MinInt32), big.NewInt(math.MaxInt32)
	case reflect.Int64:
		return big.NewInt(math.MinInt64), big.NewInt(math.MaxInt64)
	case reflect.Uint8:
		return big.NewInt(0), big.NewInt(math.MaxUint8)
	case reflect.Uint16:
		return big.NewInt(0), big.NewInt(math.MaxUint16)
	case reflect.Uint32:
		return big.NewInt(0), big.NewInt(math.MaxUint32)
	case reflect.Uint64:
		return big.NewInt(0), new(big.Int).SetUint64(math.MaxUint64)
	}
	panic("not an integer kind")
}

func toBig[T Numeric](x T) *big.Int {
	v := reflect.ValueOf(x)
	switch v.Kind() {
	case reflect.Int8, reflect.Int16, reflect.Int32, reflect.Int64:
		return big.NewInt(v.Int())
	}
	return new(big.Int).SetUint64(v.Uint())
}

func checkTarget[T Numeric](t *testing.T, where string, want []exactCell, rd targetReaders[T]) {
	t.Helper()
	kind := reflect.TypeFor[T]().Kind()
	where = where + " as " + kind.String()
	got, valid, err := rd.masked()
	if kind == reflect.Float32 || kind == reflect.Float64 {
		if err != nil {
			t.Errorf("%s: %v", where, err)
			return
		}
		if len(got) != len(want) {
			t.Errorf("%s: %d values, want %d", where, len(got), len(want))
			return
		}
		un, uerr := rd.unmasked()
		if uerr != nil {
			t.Errorf("%s (unmasked): %v", where, uerr)
			return
		}
		bad := 0
		for i, w := range want {
			var g, gu, e float64
			if kind == reflect.Float32 {
				g, gu = float64(any(got[i]).(float32)), float64(any(un[i]).(float32))
				e = float64(w.f32())
			} else {
				g, gu = any(got[i]).(float64), any(un[i]).(float64)
				e = w.f64()
			}
			switch {
			case w.null && (valid[i] || g != 0 || !math.IsNaN(gu)):
				t.Errorf("%s elem %d: null cell read as %v valid=%v, unmasked %v (want valid=false, 0, NaN)", where, i, g, valid[i], gu)
				bad++
			case !w.null && (!valid[i] || !sameFloat(g, e) || !sameFloat(gu, e)):
				t.Errorf("%s elem %d: %v valid=%v, unmasked %v, want %v (exact %s)", where, i, g, valid[i], gu, e, w)
				bad++
			}
			if bad >= 5 {
				return
			}
		}
		return
	}
	lo, hi := kindRange(kind)
	fits := true
	for _, w := range want {
		if w.null {
			continue
		}
		b, ok := w.integral()
		if !ok || b.Cmp(lo) < 0 || b.Cmp(hi) > 0 {
			fits = false
			break
		}
	}
	if !fits {
		if !errors.Is(err, ErrTypeMismatch) {
			t.Errorf("%s: err = %v, want ErrTypeMismatch (a value is not an integer in range); got %v", where, err, got)
		}
		return
	}
	if err != nil {
		t.Errorf("%s: %v", where, err)
		return
	}
	if len(got) != len(want) {
		t.Errorf("%s: %d values, want %d", where, len(got), len(want))
		return
	}
	hasNull := false
	for i, w := range want {
		if w.null {
			hasNull = true
			var zero T
			if valid[i] || got[i] != zero {
				t.Errorf("%s elem %d: null cell read as %v valid=%v", where, i, got[i], valid[i])
			}
			continue
		}
		b, _ := w.integral()
		if !valid[i] || toBig(got[i]).Cmp(b) != 0 {
			t.Errorf("%s elem %d: %v valid=%v, want %s", where, i, got[i], valid[i], b)
		}
	}
	if _, uerr := rd.unmasked(); hasNull && !errors.Is(uerr, ErrNullValue) {
		t.Errorf("%s (unmasked): err = %v, want ErrNullValue (column has nulls)", where, uerr)
	} else if !hasNull && uerr != nil {
		t.Errorf("%s (unmasked): %v", where, uerr)
	}
}

// checkAllTargets runs checkTarget for every Numeric Go type.
func checkAllTargets(t *testing.T, where string, want []exactCell, readers func(k reflect.Kind) any) {
	t.Helper()
	checkTarget(t, where, want, readers(reflect.Int8).(targetReaders[int8]))
	checkTarget(t, where, want, readers(reflect.Int16).(targetReaders[int16]))
	checkTarget(t, where, want, readers(reflect.Int32).(targetReaders[int32]))
	checkTarget(t, where, want, readers(reflect.Int64).(targetReaders[int64]))
	checkTarget(t, where, want, readers(reflect.Uint8).(targetReaders[uint8]))
	checkTarget(t, where, want, readers(reflect.Uint16).(targetReaders[uint16]))
	checkTarget(t, where, want, readers(reflect.Uint32).(targetReaders[uint32]))
	checkTarget(t, where, want, readers(reflect.Uint64).(targetReaders[uint64]))
	checkTarget(t, where, want, readers(reflect.Float32).(targetReaders[float32]))
	checkTarget(t, where, want, readers(reflect.Float64).(targetReaders[float64]))
}

func fixedReaders[T Numeric](h *BinaryTableHDU, col int, scalar bool) targetReaders[T] {
	if scalar {
		return targetReaders[T]{
			masked:   func() ([]T, []bool, error) { return ReadColumnMasked[T](h, col) },
			unmasked: func() ([]T, error) { return ReadColumn[T](h, col) },
		}
	}
	return targetReaders[T]{
		masked:   func() ([]T, []bool, error) { return ReadVectorColumnMasked[T](h, col) },
		unmasked: func() ([]T, error) { return ReadVectorColumn[T](h, col) },
	}
}

func fixedReadersFor(h *BinaryTableHDU, col int, scalar bool) func(reflect.Kind) any {
	return func(k reflect.Kind) any {
		switch k {
		case reflect.Int8:
			return fixedReaders[int8](h, col, scalar)
		case reflect.Int16:
			return fixedReaders[int16](h, col, scalar)
		case reflect.Int32:
			return fixedReaders[int32](h, col, scalar)
		case reflect.Int64:
			return fixedReaders[int64](h, col, scalar)
		case reflect.Uint8:
			return fixedReaders[uint8](h, col, scalar)
		case reflect.Uint16:
			return fixedReaders[uint16](h, col, scalar)
		case reflect.Uint32:
			return fixedReaders[uint32](h, col, scalar)
		case reflect.Uint64:
			return fixedReaders[uint64](h, col, scalar)
		case reflect.Float32:
			return fixedReaders[float32](h, col, scalar)
		}
		return fixedReaders[float64](h, col, scalar)
	}
}

func varReaders[T Numeric](h *BinaryTableHDU, col int) targetReaders[T] {
	return targetReaders[T]{
		masked: func() ([]T, []bool, error) {
			v, ok, err := ReadVarColumnMasked[T](h, col)
			if err != nil {
				return nil, nil, err
			}
			return v.values, ok.values, nil
		},
		unmasked: func() ([]T, error) {
			v, err := ReadVarColumn[T](h, col)
			if err != nil {
				return nil, err
			}
			return v.values, nil
		},
	}
}

func varReadersFor(h *BinaryTableHDU, col int) func(reflect.Kind) any {
	return func(k reflect.Kind) any {
		switch k {
		case reflect.Int8:
			return varReaders[int8](h, col)
		case reflect.Int16:
			return varReaders[int16](h, col)
		case reflect.Int32:
			return varReaders[int32](h, col)
		case reflect.Int64:
			return varReaders[int64](h, col)
		case reflect.Uint8:
			return varReaders[uint8](h, col)
		case reflect.Uint16:
			return varReaders[uint16](h, col)
		case reflect.Uint32:
			return varReaders[uint32](h, col)
		case reflect.Uint64:
			return varReaders[uint64](h, col)
		case reflect.Float32:
			return varReaders[float32](h, col)
		}
		return varReaders[float64](h, col)
	}
}

func asciiReaders[T Numeric](h *ASCIITableHDU, col int) targetReaders[T] {
	return targetReaders[T]{
		masked:   func() ([]T, []bool, error) { return ReadASCIIColumnMasked[T](h, col) },
		unmasked: func() ([]T, error) { return ReadASCIIColumn[T](h, col) },
	}
}

func asciiReadersFor(h *ASCIITableHDU, col int) func(reflect.Kind) any {
	return func(k reflect.Kind) any {
		switch k {
		case reflect.Int8:
			return asciiReaders[int8](h, col)
		case reflect.Int16:
			return asciiReaders[int16](h, col)
		case reflect.Int32:
			return asciiReaders[int32](h, col)
		case reflect.Int64:
			return asciiReaders[int64](h, col)
		case reflect.Uint8:
			return asciiReaders[uint8](h, col)
		case reflect.Uint16:
			return asciiReaders[uint16](h, col)
		case reflect.Uint32:
			return asciiReaders[uint32](h, col)
		case reflect.Uint64:
			return asciiReaders[uint64](h, col)
		case reflect.Float32:
			return asciiReaders[float32](h, col)
		}
		return asciiReaders[float64](h, col)
	}
}

// cardRat returns the exact value of a numeric card from cfitsio's value
// string (fits_parse_value), or nil when the card is absent.
func cardRat(t *testing.T, gh testutil.CfrefHDU, key string) *big.Rat {
	t.Helper()
	c := gh.Card(key)
	if c == nil {
		return nil
	}
	s := strings.ReplaceAll(strings.TrimSpace(string(c.Value)), "D", "E")
	r, ok := new(big.Rat).SetString(s)
	if !ok {
		t.Fatalf("%s = %q is not a number", key, c.Value)
	}
	return r
}

var two63 = new(big.Int).Lsh(big.NewInt(1), 63)

// colLiterals returns the exact TSCALn and TZEROn literals (1 and 0 when
// absent).
func colLiterals(t *testing.T, gh testutil.CfrefHDU, n int) (scale, zero *big.Rat) {
	t.Helper()
	scale = cardRat(t, gh, "TSCAL"+itoa(n))
	if scale == nil {
		scale = big.NewRat(1, 1)
	}
	zero = cardRat(t, gh, "TZERO"+itoa(n))
	if zero == nil {
		zero = new(big.Rat)
	}
	return scale, zero
}

// exactIntZero reports whether a column with this TSCAL/TZERO is on the
// exact integer path (TSCAL literal 1, TZERO literal an integer of any
// size), and the TZERO to add.
func exactIntZero(t *testing.T, gh testutil.CfrefHDU, gc testutil.CfrefColumn) (*big.Int, bool) {
	t.Helper()
	scale, zero := colLiterals(t, gh, gc.N)
	if scale.Cmp(big.NewRat(1, 1)) != 0 || !zero.IsInt() {
		return nil, false
	}
	return new(big.Int).Set(zero.Num()), true
}

func isIntTypecode(tc int) bool {
	switch tc {
	case cfByte, cfShort, cfLong, cfLongLong:
		return true
	}
	return false
}

// typecodeBytes is the size of one stored element of a cfitsio type code.
func typecodeBytes(tc int) int64 {
	switch tc {
	case cfByte, cfLogical, cfString:
		return 1
	case cfShort:
		return 2
	case cfLong, cfFloat:
		return 4
	case cfLongLong, cfDouble, cfComplex:
		return 8
	case cfDblComplex:
		return 16
	}
	return 0
}

// storedFloats decodes the stored E/D values of fixed column gc straight
// from the sample bytes, at the HDU data_start and the column layout
// cfitsio reports (each earlier column takes repeat elements of its type
// code, X columns ceil(repeat/8) bytes, P/Q descriptors 8/16 bytes).
func storedFloats(t *testing.T, path string, gh testutil.CfrefHDU, gc testutil.CfrefColumn) []float64 {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	var off int64
	for _, c := range gh.Table.Columns {
		if c.N == gc.N {
			break
		}
		switch {
		case c.Typecode < 0 && strings.Contains(string(c.TForm), "Q"):
			off += 16
		case c.Typecode < 0:
			off += 8
		case c.Typecode == cfBit:
			off += (c.Repeat + 7) / 8
		default:
			off += c.Repeat * typecodeBytes(c.Typecode)
		}
	}
	es := typecodeBytes(gc.Typecode)
	var out []float64
	for r := int64(0); r < gh.Table.NRows; r++ {
		base := gh.DataStart + r*gh.Table.RowLen + off
		for k := int64(0); k < gc.Repeat; k++ {
			p := b[base+k*es:]
			if gc.Typecode == cfFloat {
				out = append(out, float64(math.Float32frombits(binary.BigEndian.Uint32(p))))
			} else {
				out = append(out, math.Float64frombits(binary.BigEndian.Uint64(p)))
			}
		}
	}
	return out
}

// scaledRat is stored*scale + zero, exactly; nil for a non-finite stored
// value.
func scaledRat(stored *big.Rat, scale, zero *big.Rat) *big.Rat {
	if stored == nil {
		return nil
	}
	r := new(big.Rat).Mul(stored, scale)
	return r.Add(r, zero)
}

// binaryExact builds the exact cells of a numeric binary column (fixed or
// VLA) from the golden; ours/ourValid (our float64 read) is used only
// where cfitsio flushes a value, as described at the top of the file.
// It also checks cfitsio's own TDOUBLE read against the exact values and
// returns how many integer cells cfitsio does not read exactly.
func binaryExact(t *testing.T, where, path string, gh testutil.CfrefHDU, gc testutil.CfrefColumn, ours []float64, ourValid []bool) ([]exactCell, int) {
	t.Helper()
	tc := gc.Typecode
	if tc < 0 {
		tc = -tc
	}
	zero, exactPath := exactIntZero(t, gh, gc)
	exactPath = exactPath && isIntTypecode(tc)
	scaleR, zeroR := colLiterals(t, gh, gc.N)
	scaled := scaleR.Cmp(big.NewRat(1, 1)) != 0 || zeroR.Sign() != 0
	var stored []float64
	if !exactPath && scaled && !isIntTypecode(tc) {
		if gc.Typecode < 0 {
			t.Fatalf("%s: a scaled float VLA has no stored-value model", where)
		}
		stored = storedFloats(t, path, gh, gc)
	}
	var out []exactCell
	off := 0
	for r, cell := range gc.RowCells(t) {
		for k, p := range cell.Phys {
			i := len(out)
			if !p.Valid {
				if i < len(ours) && cfitsioFlushed(ours[i], ourValid[i], p) {
					out = append(out, exactCell{v: ours[i]})
				} else {
					out = append(out, exactCell{null: true})
				}
				continue
			}
			if !exactPath {
				v := float64(p.V)
				if i < len(ours) && cfitsioFlushed(ours[i], ourValid[i], p) {
					v = ours[i]
				}
				c := exactCell{v: v}
				if scaled {
					// Integer targets get the exact value, not v.
					var st *big.Rat
					switch {
					case isIntTypecode(tc):
						if len(cell.Raw) <= k || !cell.Raw[k].Valid {
							t.Fatalf("%s row %d elem %d: golden has no stored integer", where, r, k)
						}
						st = new(big.Rat).SetInt64(cell.Raw[k].V)
					case math.IsInf(stored[i], 0) || math.IsNaN(stored[i]):
					default:
						st = new(big.Rat).SetFloat64(stored[i])
					}
					c.r = scaledRat(st, scaleR, zeroR)
					c.noInt = c.r == nil
					if !isIntTypecode(tc) && zero != nil && c.r != nil {
						// E/D with TSCAL 1 and an integer TZERO: float
						// targets get stored + TZERO correctly rounded.
						// cfitsio adds TZERO's double, rounding twice
						// when no double holds TZERO.
						c.x = c.r
						c.v, _ = c.r.Float64()
						if c.v != float64(p.V) {
							off++
							if tol := 2 * math.Max(ulp(math.Abs(c.v)), ulp(math.Abs(stored[i]))); math.Abs(c.v-float64(p.V)) > tol {
								t.Errorf("%s row %d elem %d: cfitsio %v is far from stored+TZERO = %s; the exact model is wrong", where, r, k, p.V, c.r.RatString())
							}
						}
					}
				}
				out = append(out, c)
				continue
			}
			var stored *big.Int
			switch {
			case len(cell.Raw) > k && cell.Raw[k].Valid:
				stored = big.NewInt(cell.Raw[k].V)
			case zero.Sign() == 0 && len(cell.I64) > k && cell.I64[k].Valid:
				stored = big.NewInt(cell.I64[k].V)
			default:
				t.Fatalf("%s row %d elem %d: golden has no stored integer", where, r, k)
			}
			e := new(big.Int).Add(stored, zero)
			out = append(out, exactCell{i: e})
			// cfitsio adds TZERO in double precision: within a few ulp.
			ef, _ := new(big.Float).SetInt(e).Float64()
			if ef != float64(p.V) {
				off++
				sf, _ := new(big.Float).SetInt(stored).Float64()
				zf, _ := new(big.Float).SetInt(zero).Float64()
				tol := 4 * math.Max(ulp(math.Abs(ef)), math.Max(ulp(math.Abs(sf)), ulp(math.Abs(zf))))
				if math.Abs(ef-float64(p.V)) > tol {
					t.Errorf("%s row %d elem %d: cfitsio %v is far from stored+TZERO = %s; the exact model is wrong", where, r, k, p.V, e)
				}
			} else if (len(cell.I64) > k && cell.I64[k].Valid && (!e.IsInt64() || e.Int64() != cell.I64[k].V)) ||
				(len(cell.U64) > k && cell.U64[k].Valid && (!e.IsUint64() || e.Uint64() != cell.U64[k].V)) {
				off++
			}
		}
	}
	return out, off
}

func ulp(x float64) float64 {
	if x == 0 {
		return math.SmallestNonzeroFloat64
	}
	return math.Nextafter(x, math.Inf(1)) - x
}

// Samples where cfitsio's own reads are not exact (pinned so that a
// cfitsio change is noticed); everywhere else cfitsio must agree exactly.
var cfitsioInexactSamples = map[string]string{
	"tables/cf_bintable_kzero.fits":        "cfitsio adds an integer TZERO to K columns in double precision",
	"tables/astropy_ascii_i4-i20.fits":     "cfitsio reads ASCII I fields through a double: \"9223372036854775807\" becomes 2^63 and its TLONGLONG read wraps to -2^63",
	"tables/cf_ascii_exact.fits":           "cfitsio reads ASCII fields through a double, not correctly rounded for 20-digit fields, and wraps TLONGLONG reads beyond int64",
	"tables/cf_bintable_tzero_big.fits":    "cfitsio holds TZERO as a double: 9223372036854775807 acts as 2^63 (the uint64 convention) and 4611686018427387905 as 2^62",
	"tables/cf_bintable_scaled_exact.fits": "cfitsio holds TZERO as a double: 9223372036854775809 acts as 2^63 (uint64 XOR), and its integer reads of 2^63+n on a J column saturate",
	"tables/cf_ascii_scaled_exact.fits":    "cfitsio reads 19-digit ASCII fields through a double that is not correctly rounded (1152921504606846976 becomes 1152921504606846720)",
	"tables/cf_ascii_scaled_edge.fits":     "cfitsio holds TZERO as a double (2^53+1 acts as 2^53) and adds the field's double to it (rounding twice)",
	"tables/cf_ascii_tzero_round.fits":     "cfitsio holds TZERO as a double (2^53+1 acts as 2^53, 2^64+... rounds), reads 1.5E-900 as 0, and returns NaN for fields of more than 308 digits",
	"tables/cf_bintable_tzero_round.fits":  "cfitsio holds TZERO as a double (2^53+1 acts as 2^53, 2^64+2049 as 2^64+4096) and adds an integer TZERO to K columns in double precision",
}

var asciiFormRe = func(f string) (typ byte, frac int) {
	f = strings.TrimSpace(f)
	typ = f[0]
	if i := strings.IndexByte(f, '.'); i >= 0 {
		frac, _ = strconv.Atoi(f[i+1:])
	}
	return typ, frac
}

// asciiFieldRat reads an ASCII-table numeric field as an exact decimal
// (FITS 4.0 7.2.5): blanks are not significant, D or E starts the
// exponent, a mantissa without a decimal point gets the implied decimal
// of Fw.d / Ew.d / Dw.d, and a field with no digits is 0 (as cfitsio).
// neg reports a minus sign (for -0).
func asciiFieldRat(t *testing.T, field, tform string) (r *big.Rat, neg bool) {
	t.Helper()
	_, frac := asciiFormRe(tform)
	s := strings.ReplaceAll(field, " ", "")
	s = strings.ReplaceAll(s, "D", "E")
	mant, exp := s, ""
	if i := strings.IndexByte(s, 'E'); i >= 0 {
		mant, exp = s[:i], s[i+1:]
	}
	neg = strings.HasPrefix(mant, "-")
	if strings.Trim(mant, "+-.") == "" {
		return new(big.Rat), neg
	}
	r, ok := new(big.Rat).SetString(mant)
	if !ok {
		t.Fatalf("ASCII field %q: mantissa %q is not a decimal", field, mant)
	}
	e := 0
	if exp != "" {
		var err error
		if e, err = strconv.Atoi(exp); err != nil {
			t.Fatalf("ASCII field %q: exponent %q: %v", field, exp, err)
		}
	}
	if !strings.Contains(mant, ".") {
		e -= frac
	}
	p := new(big.Rat).SetInt(new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(abs(e))), nil))
	if e >= 0 {
		r.Mul(r, p)
	} else {
		r.Quo(r, p)
	}
	return r, neg
}

func abs(i int) int {
	if i < 0 {
		return -i
	}
	return i
}

// asciiExact builds the exact cells of a numeric ASCII column from the
// field text in cfitsio's raw rows. It also compares cfitsio's TDOUBLE
// read with the exact values (2 ulp) and returns how many cells cfitsio
// does not read exactly.
func asciiExact(t *testing.T, where string, gh testutil.CfrefHDU, gc testutil.CfrefColumn) ([]exactCell, int) {
	t.Helper()
	if gc.TBCol == nil {
		t.Fatalf("%s: no TBCOL in golden", where)
	}
	start := int(*gc.TBCol) - 1
	var tnull *string
	if gc.TNull != nil {
		s := string(*gc.TNull)
		tnull = &s
	}
	scale := float64(gc.TScal)
	scaleR, zeroR := colLiterals(t, gh, gc.N)
	zero, _ := exactIntZero(t, gh, gc)
	zf := float64(gc.TZero)
	var out []exactCell
	off := 0
	for r, raw := range gh.Table.RawRows {
		field := string(raw)[start : start+int(gc.Width)]
		if tnull != nil && len(*tnull) <= len(field) && field == *tnull+strings.Repeat(" ", len(field)-len(*tnull)) {
			out = append(out, exactCell{null: true})
			continue
		}
		v, neg := asciiFieldRat(t, field, string(gc.TForm))
		var c exactCell
		switch {
		case zero != nil:
			sum := new(big.Rat).Add(v, new(big.Rat).SetInt(zero))
			if sum.IsInt() {
				c.i = new(big.Int).Set(sum.Num())
				if sum.Sign() == 0 && neg && zero.Sign() == 0 {
					c = exactCell{v: math.Copysign(0, -1)}
				}
			} else {
				c.v, _ = sum.Float64()
				c.noInt = true
				c.x = sum
			}
		default:
			f, _ := v.Float64()
			if f == 0 && neg {
				f = math.Copysign(0, -1)
			}
			c.v = f*scale + zf
			c.r = scaledRat(v, scaleR, zeroR)
		}
		out = append(out, c)
		p := gc.RowCells(t)[r].Phys[0]
		if !p.Valid {
			continue // TNULL prefix divergence: checked in checkASCIITable
		}
		e := c.f64()
		if math.IsNaN(float64(p.V)) {
			// cfitsio accumulates the digits in a double, which overflows
			// past about 308 digits (Inf/Inf = NaN); the field is a number.
			if n := len(strings.Trim(strings.Map(func(r rune) rune {
				if r >= '0' && r <= '9' {
					return r
				}
				return -1
			}, field), "0")); n <= 308 {
				t.Errorf("%s row %d: field %q (%d significant digits): cfitsio reads NaN", where, r, field, n)
			}
			off++
			continue
		}
		if e == float64(p.V) {
			if len(gc.RowCells(t)[r].I64) == 1 && gc.RowCells(t)[r].I64[0].Valid {
				if b, ok := c.integral(); ok && (!b.IsInt64() || b.Int64() != gc.RowCells(t)[r].I64[0].V) {
					off++
				}
			}
			continue
		}
		off++
		// cfitsio adds the field's double to TZERO's double: within 2 ulp
		// of the largest of the result and the two operands (a sum that
		// cancels, -2^64 + (2^64+2051), keeps the operands' rounding).
		fv, _ := v.Float64()
		tol := 2 * math.Max(ulp(math.Abs(e)), math.Max(ulp(math.Abs(fv*scale)), ulp(math.Abs(zf))))
		if !(e == 0 && p.V == 0) && math.Abs(e-float64(p.V)) > tol {
			t.Errorf("%s row %d: field %q: cfitsio %v, exact decimal %v: more than 2 ulp apart, the test's field model is wrong", where, r, field, p.V, e)
		}
	}
	return out, off
}

// TestTargetRangesSample pins what cf_bintable_ranges.fits must exercise,
// so the generic target checks above cannot pass vacuously: for each
// column, which integer targets the golden values fit.
func TestTargetRangesSample(t *testing.T) {
	want := loadGoldenFile(t, testutil.TestdataPath(t, "tables", "golden", "cf_bintable_ranges.fits.json"))
	gh := want.HDUs[1]
	type exp struct{ col, fits, misses string }
	cases := []exp{
		{"K_I8_IN", "int8", ""}, {"K_I8_LO", "int16", "int8"}, {"K_I8_HI", "int16", "int8"},
		{"D_I16_LO", "int32", "int16"}, {"S_I16_HI", "int32", "int16"},
		{"S_I32_LO", "int64", "int32"}, {"D_U8_HI", "uint16", "uint8"}, {"D_U8_LO", "int8", "uint8"},
		{"S_U32_HI", "uint64", "uint32"}, {"K_U16_HI", "uint32", "uint16"},
		{"D_I64_IN", "int64", "uint64"}, {"D_I64_HI", "uint64", "int64"}, {"D_I64_LO", "", "int64"},
		{"D_U64_IN", "uint64", "int64"}, {"D_U64_HI", "", "uint64"}, {"S2_I64_HI", "uint64", "int64"},
		{"S2_U64_IN", "uint64", "int64"}, {"S2_U64_HI", "uint64", "int64"}, {"D_FRAC", "", "int64"},
	}
	for _, c := range cases {
		gc := gh.Table.Column(c.col)
		if gc == nil {
			t.Fatalf("no column %s", c.col)
		}
		cells, _ := binaryExact(t, c.col, testutil.TestdataPath(t, "tables", "cf_bintable_ranges.fits"), gh, *gc, nil, nil)
		for _, k := range []string{c.fits, c.misses} {
			if k == "" {
				continue
			}
			lo, hi := kindRange(map[string]reflect.Kind{"int8": reflect.Int8, "int16": reflect.Int16, "int32": reflect.Int32, "int64": reflect.Int64, "uint8": reflect.Uint8, "uint16": reflect.Uint16, "uint32": reflect.Uint32, "uint64": reflect.Uint64}[k])
			ok := true
			for _, e := range cells {
				b, isInt := e.integral()
				ok = ok && isInt && b.Cmp(lo) >= 0 && b.Cmp(hi) <= 0
			}
			if ok != (k == c.fits) {
				t.Errorf("%s: fits %s = %v, sample is not what gen_tables.c says (cells %v)", c.col, k, ok, cells)
			}
		}
	}
}
