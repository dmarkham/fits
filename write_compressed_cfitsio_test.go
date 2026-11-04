package fits

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/binary"
	"encoding/json"
	"fmt"
	"io"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/compress"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/crefgen"
	"github.com/dmarkham/fits/internal/testutil"
)

// Writer parity against libcfitsio's own tile-compression writer.
//
// compress/testdata/cref/golden/img holds, for every case, the file
// cfitsio 4.7.0 writes with fits_img_compress -- the code path fpack uses --
// under the parameters in <case>.json (algorithm, tile shape, quantize
// level, dither method, ZDITHER0, HCOMPRESS scale). The input pixels are
// not committed: internal/crefgen regenerates them from the generator
// parameters in <case>.json and they must hash to the sha256 cref.c
// recorded for the input it gave cfitsio. See compress/testdata/cref/cref.c.
//
// TestCompressedWriterMatchesCfitsio writes the same input with our
// writer and requires the result to be the same compressed image:
//
//   - the same structural compression keywords (ZBITPIX, ZNAXISn, ZTILEn,
//     ZCMPTYPE, ZNAMEi/ZVALi, ZQUANTIZ, ZDITHER0, ZBLANK, BZERO/BSCALE);
//   - the same table columns, row for row;
//   - RICE_1, HCOMPRESS_1 and PLIO_1 tiles byte-identical;
//   - GZIP_1/GZIP_2 tiles (and GZIP_COMPRESSED_DATA fallback tiles) that
//     inflate to identical bytes (DEFLATE output itself is not unique);
//   - identical ZSCALE/ZZERO per tile (float64 bit-exact).
//
// Identical tiles and keywords mean every reader (cfitsio, funpack,
// astropy) decodes our file exactly as it decodes cfitsio's. Cases where
// cfitsio refused the input ("status" != 0) require our writer to refuse
// it too. No external tool is needed at test time.

const imgGoldenDir = "compress/testdata/cref/golden/img"

// imgwGoldenDir holds images cfitsio compressed while they were written
// (fits_write_imgnull into a compressing HDU with null value -1e30)
// rather than through fits_img_compress: MANIFEST kind "imgw".
const imgwGoldenDir = "compress/testdata/cref/golden/imgw"

type imgGolden struct {
	Name             string  `json:"name"`
	Type             string  `json:"type"`
	Naxes            []int64 `json:"naxes"`
	Algo             string  `json:"algo"`
	Tile             []int64 `json:"tile"`
	QLevel           float64 `json:"qlevel"`
	QMethod          int     `json:"qmethod"`
	Seed             int     `json:"seed"`
	HScale           float64 `json:"hscale"`
	LosslessVerified bool    `json:"lossless_verified"`
	SrcHeader        bool    `json:"src_header"`
	Status           int     `json:"status"`
	Error            string  `json:"error"`
	Kind             int     `json:"kind"`      // cref.c generator kind
	GenSeed          uint64  `json:"gen_seed"`  // cref.c generator seed
	InSHA256         string  `json:"in_sha256"` // sha256 of the input cfitsio was given
	// Dec is cfitsio's decode of <case>.fits for lossy cases
	// (fits_read_img TDOUBLE; nulls flagged; float32 images rounded to
	// float32), an exact digest.
	Dec testutil.Values[testutil.NullFloat] `json:"dec"`
	// WriteImgnull marks an imgw case (see imgwGoldenDir).
	WriteImgnull bool `json:"write_imgnull"`
	// NullValue is the null value an imgw case gave fits_write_imgnull
	// (0: -1e30, with the NaN pixels passed as it). Our writer gets
	// Nullcheck with NullValue (float32) or NullValue64 (float64).
	NullValue float64 `json:"nullvalue"`

	dir string // golden directory the case was loaded from
}

// imgGoldenByName is one case of imgGoldenCases.
func imgGoldenByName(t *testing.T, name string) imgGolden {
	t.Helper()
	for _, g := range imgGoldenCases(t) {
		if g.Name == name {
			return g
		}
	}
	t.Fatalf("committed golden %s missing from MANIFEST (failure, not a skip)", name)
	return imgGolden{}
}

// imgGoldenCases reads the "img" entries of golden/MANIFEST. The count is
// pinned so a partial deletion fails rather than shrinking the test.
func imgGoldenCases(t *testing.T) []imgGolden {
	t.Helper()
	return manifestGoldens(t, "img", imgGoldenDir, 96)
}

// img64GoldenDir holds BITPIX -64 images compressed by fits_img_compress
// (fits_quantize_double): MANIFEST kind "img64".
const img64GoldenDir = "compress/testdata/cref/golden/img64"

// img64GoldenCases reads the "img64" entries of golden/MANIFEST.
func img64GoldenCases(t *testing.T) []imgGolden {
	t.Helper()
	return manifestGoldens(t, "img64", img64GoldenDir, 20)
}

// imgwGoldenCases reads the "imgw" entries of golden/MANIFEST.
func imgwGoldenCases(t *testing.T) []imgGolden {
	t.Helper()
	return manifestGoldens(t, "imgw", imgwGoldenDir, 5)
}

// img2GoldenDir holds more fits_img_compress cases (MANIFEST kind
// "img2"), kept out of golden/img so the set the reader tests pin is
// unchanged: SUBTRACTIVE_DITHER_2 with NaNs and zeros, a SIMPLE = F
// source, and three cfitsio quantizer overflows (knownCfitsioOverflow).
const img2GoldenDir = "compress/testdata/cref/golden/img2"

// img2GoldenCases reads the "img2" entries of golden/MANIFEST.
func img2GoldenCases(t *testing.T) []imgGolden {
	t.Helper()
	return manifestGoldens(t, "img2", img2GoldenDir, 5)
}

func manifestGoldens(t *testing.T, kind, dir string, want int) []imgGolden {
	t.Helper()
	mpath := filepath.Join("compress", "testdata", "cref", "golden", "MANIFEST")
	f, err := os.Open(mpath)
	if err != nil {
		t.Fatalf("committed golden missing (failure, not a skip): %v", err)
	}
	defer f.Close()
	var out []imgGolden
	sc := bufio.NewScanner(f)
	for sc.Scan() {
		fs := strings.Fields(sc.Text())
		if len(fs) != 2 || fs[0] != kind {
			continue
		}
		var g imgGolden
		b, err := os.ReadFile(filepath.Join(dir, fs[1]+".json"))
		if err != nil {
			t.Fatalf("committed golden missing (failure, not a skip): %v", err)
		}
		dec := json.NewDecoder(bytes.NewReader(b))
		dec.DisallowUnknownFields()
		if err := dec.Decode(&g); err != nil {
			t.Fatalf("parse %s.json: %v", fs[1], err)
		}
		if g.WriteImgnull != (kind == "imgw") {
			t.Fatalf("%s/%s.json: write_imgnull %v in MANIFEST kind %s", dir, fs[1], g.WriteImgnull, kind)
		}
		g.dir = dir
		out = append(out, g)
	}
	if len(out) != want {
		t.Fatalf("MANIFEST lists %d %s goldens, want %d", len(out), kind, want)
	}
	return out
}

func (g imgGolden) npix() int {
	n := 1
	for _, v := range g.Naxes {
		n *= int(v)
	}
	return n
}

func (g imgGolden) algorithm(t *testing.T) compress.Algorithm {
	a := compress.ParseAlgorithm(g.Algo)
	if a == compress.Unknown {
		t.Fatalf("%s: unknown algorithm %q", g.Name, g.Algo)
	}
	return a
}

// readInput regenerates the case's input pixels (little-endian, native
// type, the bytes cref.c gave cfitsio) and checks them against the
// recorded sha256.
func (g imgGolden) readInput(t *testing.T) []byte {
	t.Helper()
	n := int64(g.npix())
	var b []byte
	switch g.Type {
	case "f32":
		b = crefgen.LEFloat32(crefgen.ImgFloat(n, g.Naxes[0], g.Kind, g.GenSeed))
	case "f64":
		b = crefgen.LEFloat64(crefgen.ImgDouble(n, g.Naxes[0], g.Kind, g.GenSeed))
	default:
		b = crefgen.LEInts(crefgen.ImgInt(n, g.Naxes[0], g.Type, g.Kind, g.GenSeed), crefgen.TypeWidth(g.Type))
	}
	if err := crefgen.Check("img "+g.Name, b, g.InSHA256); err != nil {
		t.Fatal(err)
	}
	return b
}

// decTruth is cfitsio's decode of golden img/<case>.fits for a lossy case,
// element by element. The golden stores it as an exact digest, so the
// caller passes its own reading cand (float64, NaN = null; for a float32
// image the float32 values widened): it is verified against the digest
// exactly (every value and null flag, sha256 plus the sample) and, once
// equal, returned as the truth. A difference fails the test.
func (g imgGolden) decTruth(t *testing.T, where string, cand []float64) []float64 {
	t.Helper()
	if !g.Dec.Present() {
		t.Fatalf("%s: golden %s has no cfitsio decode (lossless case?)", where, g.Name)
	}
	got := g.Dec.Resolve(t, where+" vs cfitsio decode", func() []testutil.NullFloat {
		out := make([]testutil.NullFloat, len(cand))
		for i, v := range cand {
			if !math.IsNaN(v) {
				out[i] = testutil.NullFloat{V: v, Valid: true}
			}
		}
		return out
	})
	out := make([]float64, len(got))
	for i, v := range got {
		out[i] = v.V
		if !v.Valid {
			out[i] = math.NaN()
		}
	}
	return out
}

func leSlice[T any](b []byte, conv func([]byte) T, width int) []T {
	out := make([]T, len(b)/width)
	for i := range out {
		out[i] = conv(b[i*width:])
	}
	return out
}

func (g imgGolden) floatInput(t *testing.T) []float32 {
	return leSlice(g.readInput(t), func(b []byte) float32 { return math.Float32frombits(binary.LittleEndian.Uint32(b)) }, 4)
}

func (g imgGolden) doubleInput(t *testing.T) []float64 {
	return leSlice(g.readInput(t), func(b []byte) float64 { return math.Float64frombits(binary.LittleEndian.Uint64(b)) }, 8)
}

// writeOurs writes the golden's input with our writer, using the options
// equivalent to the cfitsio parameters, into dir/<case>.fits.
func (g imgGolden) writeOurs(t *testing.T, dir string) (string, error) {
	t.Helper()
	path := filepath.Join(dir, g.Name+".fits")
	f, err := Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := WriteImage(f, nil, []int64{}, []uint8{}); err != nil {
		t.Fatal(err)
	}
	var tile []int64
	if len(g.Tile) > 0 {
		tile = g.Tile
	}
	in := g.readInput(t)
	hdr := g.srcHeader(t)
	var werr error
	if g.Type == "f32" || g.Type == "f64" {
		opts := CompressFloatOptions{
			Algorithm: g.algorithm(t),
			TileShape: tile,
			QLevel:    float32(g.QLevel),
			ZDither0:  g.Seed,
		}
		if g.NullValue != 0 {
			opts.Nullcheck = true
			if g.Type == "f64" {
				opts.NullValue64 = g.NullValue
			} else {
				opts.NullValue = float32(g.NullValue)
			}
		}
		switch g.QMethod {
		case -1:
			opts.NoDither = true
		case 1:
			opts.DitherMethod = compress.SubtractiveDither1
		case 2:
			opts.DitherMethod = compress.SubtractiveDither2
		}
		if g.Type == "f64" {
			_, werr = AppendCompressedFloat64Image(f, hdr, g.Naxes, g.doubleInput(t), opts)
		} else {
			_, werr = AppendCompressedFloat32Image(f, hdr, g.Naxes, g.floatInput(t), opts)
		}
	} else {
		opts := CompressOptions{Algorithm: g.algorithm(t), TileShape: tile}
		if g.HScale < 0 {
			opts.HCompressScale = int(-g.HScale)
		} else if g.HScale > 0 {
			t.Fatalf("%s: noise-relative HCOMPRESS scale %g is not a writer option", g.Name, g.HScale)
		}
		switch g.Type {
		case "u8":
			_, werr = AppendCompressedImage(f, hdr, g.Naxes, in, opts)
		case "s8":
			_, werr = AppendCompressedImage(f, hdr, g.Naxes, leSlice(in, func(b []byte) int8 { return int8(b[0]) }, 1), opts)
		case "i16":
			_, werr = AppendCompressedImage(f, hdr, g.Naxes, leSlice(in, func(b []byte) int16 { return int16(binary.LittleEndian.Uint16(b)) }, 2), opts)
		case "u16":
			_, werr = AppendCompressedImage(f, hdr, g.Naxes, leSlice(in, binary.LittleEndian.Uint16, 2), opts)
		case "i32":
			_, werr = AppendCompressedImage(f, hdr, g.Naxes, leSlice(in, func(b []byte) int32 { return int32(binary.LittleEndian.Uint32(b)) }, 4), opts)
		case "u32":
			_, werr = AppendCompressedImage(f, hdr, g.Naxes, leSlice(in, binary.LittleEndian.Uint32, 4), opts)
		default:
			t.Fatalf("%s: type %q", g.Name, g.Type)
		}
	}
	if cerr := f.Close(); cerr != nil && werr == nil {
		t.Fatal(cerr)
	}
	return path, werr
}

// srcHeader is the primary header of <case>.src.fits -- the uncompressed
// image cfitsio compressed, written by cfitsio -- for "src_header" cases,
// else nil. Our writer gets it as the caller header.
func (g imgGolden) srcHeader(t *testing.T) *header.Header {
	t.Helper()
	if !g.SrcHeader {
		return nil
	}
	f, err := Open(goldenFITS(t, filepath.Join(g.dir, g.Name+".src.fits")))
	if err != nil {
		t.Fatalf("committed golden missing or unreadable (failure, not a skip): %v", err)
	}
	defer f.Close()
	h, err := f.HDU(0)
	if err != nil {
		t.Fatal(err)
	}
	return h.Header().Clone()
}

// tableLayoutKey reports the BINTABLE keywords the writer computes:
// XTENSION, BITPIX, NAXIS, NAXISn, PCOUNT, GCOUNT, TFIELDS, END and the
// per-column TFORMn ... TBCOLn.
func tableLayoutKey(k string) bool {
	switch k {
	case header.KeyXtension, header.KeyBitpix, header.KeyNaxis, header.KeyPcount,
		header.KeyGcount, header.KeyTfields, header.KeyEnd:
		return true
	}
	root, _, ok := indexedKey(k)
	if !ok {
		return false
	}
	switch root {
	case "NAXIS", "TFORM", "TTYPE", "TUNIT", "TDISP", "TDIM", "TSCAL", "TZERO", "TNULL", "TBCOL":
		return true
	}
	return false
}

// passThroughKeys lists every non-structural card of a compressed HDU
// header -- the caller's keywords as the writer copied or renamed them
// -- as "KEY=value" strings in order (COMMENT/HISTORY by text). EXTNAME
// is left out: cfitsio adds EXTNAME = 'COMPRESSED_IMAGE' when the source
// has none, our writer does not (a documented difference).
func passThroughKeys(h *header.Header) []string {
	var out []string
	for _, c := range h.Cards() {
		k := c.Key
		if tableLayoutKey(k) || k == "EXTNAME" || k == "ZIMAGE" || k == "ZCMPTYPE" ||
			k == "ZBITPIX" || k == "ZQUANTIZ" || k == "ZDITHER0" || k == "ZBLANK" ||
			k == "BZERO" || k == "BSCALE" || strings.HasPrefix(k, "ZNAXIS") ||
			strings.HasPrefix(k, "ZTILE") || strings.HasPrefix(k, "ZNAME") || strings.HasPrefix(k, "ZVAL") {
			continue
		}
		v := ""
		if c.Value != nil {
			v = strings.TrimSpace(fmt.Sprint(c.Value))
		}
		if k == "COMMENT" || k == "HISTORY" {
			v = strings.TrimSpace(c.Comment + v)
		}
		out = append(out, k+"="+v)
	}
	return out
}

// goldenFITS returns a readable path for a committed cfitsio-written FITS
// golden. The goldens are committed gzipped (<name>.fits.gz, gzip -n -9
// of the bytes cfitsio wrote; `make verify` compares the inflated bytes
// with a fresh cfitsio run): when path itself does not exist, path+".gz"
// is inflated into the test's temp directory and that copy's path is
// returned. Any other path is returned unchanged.
func goldenFITS(t *testing.T, path string) string {
	t.Helper()
	if _, err := os.Stat(path); err == nil {
		return path
	}
	gz, err := os.ReadFile(path + ".gz")
	if err != nil {
		t.Fatalf("committed golden %s(.gz) missing (failure, not a skip): %v", path, err)
	}
	b, err := gunzip(gz)
	if err != nil {
		t.Fatalf("%s.gz: %v", path, err)
	}
	out := filepath.Join(t.TempDir(), filepath.Base(path))
	if err := os.WriteFile(out, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return out
}

// openCompressed opens HDU 1 of path as a tile-compressed image. A
// committed golden stored as path+".gz" is read through goldenFITS.
func openCompressed(t *testing.T, path string) (*File, *CompressedImageHDU) {
	t.Helper()
	path = goldenFITS(t, path)
	f, err := Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	h, err := f.HDU(1)
	if err != nil {
		t.Fatalf("%s HDU 1: %v", path, err)
	}
	c, ok := h.(*CompressedImageHDU)
	if !ok {
		t.Fatalf("%s HDU 1 is %T, want *CompressedImageHDU", path, h)
	}
	return f, c
}

// compressionKeys are the keywords that define how the image is stored.
// A value is compared numerically when both sides parse as numbers.
func compressionKeys(h *header.Header, naxis int) map[string]string {
	keys := []string{"ZBITPIX", "ZNAXIS", "ZCMPTYPE", "ZQUANTIZ", "ZDITHER0", "ZBLANK", "BZERO", "BSCALE"}
	for i := 1; i <= naxis; i++ {
		keys = append(keys, fmt.Sprintf("ZNAXIS%d", i), fmt.Sprintf("ZTILE%d", i))
	}
	for i := 1; i <= 4; i++ {
		keys = append(keys, fmt.Sprintf("ZNAME%d", i), fmt.Sprintf("ZVAL%d", i))
	}
	out := map[string]string{}
	for _, k := range keys {
		c, ok := h.Get(k)
		if !ok {
			continue
		}
		switch v := c.Value.(type) {
		case string:
			out[k] = strings.TrimSpace(v)
		case bool:
			out[k] = fmt.Sprint(v)
		default:
			f, err := h.Float(k)
			if err != nil {
				out[k] = fmt.Sprint(v)
			} else {
				out[k] = fmt.Sprintf("%.17g", f)
			}
		}
	}
	return out
}

type tableCols struct {
	names []string
	data  map[string][][]byte  // VLA cells as raw bytes (int16 cells big-endian)
	f64   map[string][]float64 // ZSCALE / ZZERO
}

func readTableCols(t *testing.T, c *CompressedImageHDU) tableCols {
	t.Helper()
	cols, err := c.tbl.Columns()
	if err != nil {
		t.Fatal(err)
	}
	tc := tableCols{data: map[string][][]byte{}, f64: map[string][]float64{}}
	for _, col := range cols {
		name := strings.TrimSpace(col.Name)
		tc.names = append(tc.names, name)
		switch name {
		case "ZSCALE", "ZZERO":
			v, err := ReadColumn[float64](c.tbl, col.Index)
			if err != nil {
				t.Fatalf("read %s: %v", name, err)
			}
			tc.f64[name] = v
		default:
			if strings.Contains(col.TForm, "PI") || strings.Contains(col.TForm, "QI") {
				vc, err := ReadVarColumn[int16](c.tbl, col.Index)
				if err != nil {
					t.Fatalf("read %s: %v", name, err)
				}
				rows := make([][]byte, vc.Len())
				for r := range rows {
					cell := vc.At(r)
					b := make([]byte, 2*len(cell))
					for i, v := range cell {
						binary.BigEndian.PutUint16(b[2*i:], uint16(v))
					}
					rows[r] = b
				}
				tc.data[name] = rows
			} else {
				vc, err := ReadVarColumn[uint8](c.tbl, col.Index)
				if err != nil {
					t.Fatalf("read %s: %v", name, err)
				}
				rows := make([][]byte, vc.Len())
				for r := range rows {
					rows[r] = append([]byte(nil), vc.At(r)...)
				}
				tc.data[name] = rows
			}
		}
	}
	return tc
}

func gunzip(b []byte) ([]byte, error) {
	if len(b) == 0 {
		return nil, nil
	}
	r, err := gzip.NewReader(bytes.NewReader(b))
	if err != nil {
		return nil, err
	}
	defer r.Close()
	return io.ReadAll(r)
}

// knownRefusals are golden cases cfitsio writes but our writer refuses on
// purpose, with the reason (each is proven by its own test).
var knownRefusals = map[string]string{
	// TestPLIOUint16AstropyMisreadsCfitsio
	"plio_u16": "PLIO_1 uint16: astropy reads cfitsio's storage 32768 too high",
	// cfitsio's imcomp_init_table silently enlarges an HCOMPRESS_1 tile
	// whose last tile would be thinner than 4 pixels (8 -> 11 here); our
	// writer refuses rather than write a tile shape the caller did not
	// ask for (TestWriteCompressedRefusals/HCOMPRESS_edge_tile_of_3).
	"hcomp_i16_11x8_t8x8": "HCOMPRESS_1 3-pixel edge tile: cfitsio enlarges the tile, we refuse",
	// BLANK is not allowed with BITPIX < 0 (FITS 4.0 sect. 4.4.2.5);
	// cfitsio copies it into the compressed header and every reader then
	// nulls real pixels (TestFloatBLANKNullsRealPixels).
	"rice_f32_blank_hdr": "BLANK on a float32 image: cfitsio copies it and real pixels decode as null, we refuse",
}

// knownFloatNullValue are golden cases holding a real pixel equal to
// cfitsio's FLOATNULLVALUE (-9.11912e-36) in a GZIP_COMPRESSED_DATA tile
// that also has a NaN: cfitsio's writer stores it as a null (all-ones NaN
// bits, imcomp_float2nan), ours keeps the value (a documented
// difference). Value: the tile row (1-based) and the pixel's index in
// that tile. The tiles must be identical everywhere else.
var knownFloatNullValue = map[string][2]int{
	"rice_f32_floatnullvalue_t40x5": {1, 10},
	// float64: DOUBLENULLVALUE (-9.1191291391491e-36), the same difference
	"rice_f64_dnv_t40x5": {1, 10},
}

// knownCfitsioOverflow are golden cases where cfitsio 4.7.0's quantizer
// overflows and writes wrong pixels with status 0 (its decode, the golden,
// is far from the input). Our writer stores those tiles losslessly in
// GZIP_COMPRESSED_DATA instead (TestQuantizerOverflowStoredLosslessly).
var knownCfitsioOverflow = map[string]string{
	"rice_f32_zpt_overflow_abs":  "zeropt/delta = 1e22 overflows the int64 iqfactor",
	"rice_f32_center_overflow":   "the float sum minval+maxval rounds low and the maximum quantizes to 2^31",
	"rice_f64_zpt_overflow_abs":  "zeropt/delta = 1e22 overflows the int64 iqfactor",
	"rice_f32_null_sentinel_hit": "the centered minimum quantizes to NULL_VALUE and reads back as a null",
}

func TestCompressedWriterMatchesCfitsio(t *testing.T) {
	for _, g := range imgGoldenCases(t) {
		t.Run(g.Name, func(t *testing.T) { matchCfitsioGolden(t, g) })
	}
	for _, g := range img2GoldenCases(t) {
		t.Run("img2/"+g.Name, func(t *testing.T) {
			if _, ok := knownCfitsioOverflow[g.Name]; ok {
				t.Logf("cfitsio overflow case, see TestQuantizerOverflowStoredLosslessly")
				return
			}
			matchCfitsioGolden(t, g)
		})
	}
}

// TestReaderDecodesCfitsioImg2Goldens: our reader on the cfitsio-written
// golden/img2 files (TestReaderDecodesCfitsioWrittenGoldens covers
// golden/img): equal to cfitsio's decode of the same file -- the exact
// digest for lossy cases (also the two overflow cases, where cfitsio's
// pixels are wrong but are what the file holds), the cfitsio-verified
// input for lossless integer cases.
func TestReaderDecodesCfitsioImg2Goldens(t *testing.T) {
	for _, g := range img2GoldenCases(t) {
		t.Run(g.Name, func(t *testing.T) {
			f, h := openCompressed(t, filepath.Join(g.dir, g.Name+".fits"))
			defer f.Close()
			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatal(err)
			}
			if g.Type == "f32" {
				g.decTruth(t, "our reader on cfitsio's file", got)
				return
			}
			if !g.LosslessVerified || g.Type != "i16" {
				t.Fatalf("%s: unexpected lossless case type %s", g.Name, g.Type)
			}
			in := leSlice(g.readInput(t), func(b []byte) int16 { return int16(binary.LittleEndian.Uint16(b)) }, 2)
			for i, v := range in {
				if got[i] != float64(v) {
					t.Fatalf("pixel %d = %v, cfitsio-verified input %v", i, got[i], v)
				}
			}
		})
	}
}

// TestCompressedWriterFloat64MatchesCfitsio: AppendCompressedFloat64Image
// against fits_img_compress on BITPIX -64 images (golden/img64: every tile
// quantized by fits_quantize_double, or stored in GZIP_COMPRESSED_DATA):
// the same keywords, tiles, ZSCALE and ZZERO as for float32 (one
// documented difference, the real DOUBLENULLVALUE pixel of
// rice_f64_dnv_t40x5), and our reader's float64 decode of our file equal,
// value for value, to cfitsio's TDOUBLE decode of its own file.
// rice_f64_precise_1e8 has a 1e8 background with 0.01 noise, which only
// a double-precision quantizer reproduces.
func TestCompressedWriterFloat64MatchesCfitsio(t *testing.T) {
	for _, g := range img64GoldenCases(t) {
		t.Run(g.Name, func(t *testing.T) {
			if _, ok := knownCfitsioOverflow[g.Name]; ok {
				t.Logf("cfitsio overflow case, see TestQuantizerOverflowStoredLosslessly")
				return
			}
			matchCfitsioGolden(t, g)
			if t.Failed() {
				return
			}
			path, err := g.writeOurs(t, t.TempDir())
			if err != nil {
				t.Fatal(err)
			}
			ff, h := openCompressed(t, path)
			defer ff.Close()
			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatal(err)
			}
			if k, ok := knownFloatNullValue[g.Name]; ok {
				// cfitsio's file has the real DOUBLENULLVALUE pixel null;
				// ours keeps it, exactly (a lossless fallback tile).
				i := (k[0]-1)*int(g.Tile[0]) + k[1]
				if got[i] != crefgen.DoubleNullValue {
					t.Fatalf("pixel %d = %v, want the real DOUBLENULLVALUE %v", i, got[i], crefgen.DoubleNullValue)
				}
				got[i] = math.NaN() // the rest must equal cfitsio's decode
			}
			g.decTruth(t, "our float64 file read as float64", got)
		})
	}
}

// TestCompressedWriterFloatNullValueQuantizedMatchesCfitsio: a real pixel
// equal to FLOATNULLVALUE (-9.11912e-36) in QUANTIZED tiles that also hold
// a NaN. fits_img_compress cannot be the truth there (it uses
// FLOATNULLVALUE itself as the null), so cfitsio wrote these images with
// fits_write_imgnull and null value -1e30 (golden/imgw): it quantizes the
// pixel as a value. Our writer (NaN = null, the FLOATNULLVALUE pixel
// real) must produce the same tiles, ZSCALE and ZZERO, and cfitsio's
// decode of its file, which this also checks against ours, holds the
// pixel as a value.
func TestCompressedWriterFloatNullValueQuantizedMatchesCfitsio(t *testing.T) {
	ran := 0
	defer func() {
		if ran != 2 {
			t.Errorf("%d FLOATNULLVALUE imgw cases, want 2", ran)
		}
	}()
	for _, g := range imgwGoldenCases(t) {
		if g.NullValue != 0 {
			continue // TestCompressedWriterNullValueMatchesCfitsio
		}
		ran++
		t.Run(g.Name, func(t *testing.T) {
			matchCfitsioGolden(t, g)
			in := g.floatInput(t)
			path, err := g.writeOurs(t, t.TempDir())
			if err != nil {
				t.Fatal(err)
			}
			ff, h := openCompressed(t, path)
			defer ff.Close()
			got, err := ReadPixelsCompressed[float64](h)
			if err != nil {
				t.Fatal(err)
			}
			truth := g.decTruth(t, "our file read as float64", got)
			nfnv := 0
			for i, v := range in {
				if v != crefgen.FloatNullValue {
					continue
				}
				nfnv++
				if math.IsNaN(truth[i]) {
					t.Errorf("pixel %d: real FLOATNULLVALUE decoded as null by cfitsio", i)
				}
			}
			if nfnv == 0 {
				t.Fatal("input holds no FLOATNULLVALUE pixel")
			}
		})
	}
}

// matchCfitsioGolden writes golden g's input with our writer and requires
// the same compressed image as cfitsio's <case>.fits (see the comment at
// the top of this file).
func matchCfitsioGolden(t *testing.T, g imgGolden) {
	t.Helper()
	{
		{
			dir := t.TempDir()
			ours, werr := g.writeOurs(t, dir)
			if why, ok := knownRefusals[g.Name]; ok {
				if werr == nil {
					t.Fatalf("expected our writer to refuse this case (%s)", why)
				}
				t.Logf("refused on purpose: %v", werr)
				return
			}
			if g.Status != 0 {
				if werr == nil {
					t.Fatalf("cfitsio refuses this input (status %d: %s) but our writer accepted it", g.Status, g.Error)
				}
				t.Logf("refused like cfitsio: %v", werr)
				return
			}
			if werr != nil {
				t.Fatalf("our writer failed where cfitsio succeeded: %v", werr)
			}
			ff, fc := openCompressed(t, ours)
			defer ff.Close()
			gf, gc := openCompressed(t, filepath.Join(g.dir, g.Name+".fits"))
			defer gf.Close()

			wantK := compressionKeys(gc.Header(), len(g.Naxes))
			gotK := compressionKeys(fc.Header(), len(g.Naxes))
			for k, w := range wantK {
				if gotK[k] != w {
					t.Errorf("keyword %s: ours %q, cfitsio %q", k, gotK[k], w)
				}
			}
			for k, v := range gotK {
				if _, ok := wantK[k]; !ok {
					t.Errorf("keyword %s = %q written by us, absent from cfitsio's file", k, v)
				}
			}

			if g.SrcHeader {
				w, o := passThroughKeys(gc.Header()), passThroughKeys(fc.Header())
				if strings.Join(o, "\n") != strings.Join(w, "\n") {
					t.Errorf("caller keywords as written:\n  ours    %q\n  cfitsio %q", o, w)
				}
			}

			want := readTableCols(t, gc)
			got := readTableCols(t, fc)
			if strings.Join(got.names, ",") != strings.Join(want.names, ",") {
				t.Fatalf("columns: ours %v, cfitsio %v", got.names, want.names)
			}
			gzipAlgo := g.Algo == "GZIP_1" || g.Algo == "GZIP_2"
			for name, wrows := range want.data {
				grows := got.data[name]
				if len(grows) != len(wrows) {
					t.Fatalf("%s: %d rows, cfitsio %d", name, len(grows), len(wrows))
				}
				bad := 0
				for r := range wrows {
					w, gg := wrows[r], grows[r]
					if gzipAlgo || name == "GZIP_COMPRESSED_DATA" {
						var err error
						if w, err = gunzip(w); err != nil {
							t.Fatalf("%s row %d: cfitsio cell does not gunzip: %v", name, r+1, err)
						}
						if gg, err = gunzip(gg); err != nil {
							t.Fatalf("%s row %d: our cell does not gunzip: %v", name, r+1, err)
						}
					}
					if k, ok := knownFloatNullValue[g.Name]; ok && name == "GZIP_COMPRESSED_DATA" && r+1 == k[0] {
						wd := 4
						wantOurs := binary.BigEndian.AppendUint32(nil, math.Float32bits(crefgen.FloatNullValue))
						if g.Type == "f64" {
							wd = 8
							wantOurs = binary.BigEndian.AppendUint64(nil, math.Float64bits(crefgen.DoubleNullValue))
						}
						o := wd * k[1]
						allOnes := bytes.Repeat([]byte{0xFF}, wd)
						if len(w) < o+wd || len(gg) < o+wd ||
							!bytes.Equal(w[o:o+wd], allOnes) || !bytes.Equal(gg[o:o+wd], wantOurs) {
							t.Fatalf("%s row %d pixel %d: want cfitsio's null (all ones) and our FLOATNULLVALUE/DOUBLENULLVALUE bits (%x); got cfitsio %x, ours %x",
								name, r+1, k[1], wantOurs, w[o:min(o+wd, len(w))], gg[o:min(o+wd, len(gg))])
						}
						gg = append([]byte(nil), gg...)
						copy(gg[o:o+wd], w[o:o+wd]) // the one documented difference; the rest must match
					}
					if !bytes.Equal(gg, w) {
						if bad < 3 {
							t.Errorf("%s row %d differs from cfitsio: %s", name, r+1, diffBytes(gg, w))
						}
						bad++
					}
				}
				if bad > 0 {
					t.Errorf("%s: %d/%d tiles differ from cfitsio", name, bad, len(wrows))
				}
			}
			for name, w := range want.f64 {
				gv := got.f64[name]
				bad := 0
				for r := range w {
					if r >= len(gv) || math.Float64bits(gv[r]) != math.Float64bits(w[r]) {
						if bad < 3 {
							var o float64
							if r < len(gv) {
								o = gv[r]
							}
							t.Errorf("%s row %d: ours %.17g, cfitsio %.17g", name, r+1, o, w[r])
						}
						bad++
					}
				}
				if bad > 0 {
					t.Errorf("%s: %d/%d rows differ from cfitsio", name, bad, len(w))
				}
			}
		}
	}
}

func diffBytes(got, want []byte) string {
	n := len(got)
	if len(want) < n {
		n = len(want)
	}
	for i := 0; i < n; i++ {
		if got[i] != want[i] {
			return fmt.Sprintf("len %d vs %d, first difference at byte %d (%#02x vs %#02x)", len(got), len(want), i, got[i], want[i])
		}
	}
	return fmt.Sprintf("len %d vs %d, common prefix identical", len(got), len(want))
}

// TestHcompressDefaultTileMatchesCfitsio checks the default HCOMPRESS_1
// tile shape against the one cfitsio's imcomp_init_table picks, for every
// image height 4..2000 and for the first heights that reach the 26, 22,
// 18 and 17-row fallbacks (1681, 8401, 13443, 166323). The table
// golden/hcomp_default_tiles.txt was written by `cref encode` (cfitsio
// 4.7.0 fits_create_img on a compressing file; see cref.c
// run_hcomp_default_tiles). Heights 30 and 14 from cfitsio's candidate
// list are never chosen for NAXIS2 > 30 in this range; 30 is the
// whole-image boundary (naxis2 == 30 -> 30 rows).
func TestHcompressDefaultTileMatchesCfitsio(t *testing.T) {
	b, err := os.ReadFile(filepath.Join("compress", "testdata", "cref", "golden", "hcomp_default_tiles.txt"))
	if err != nil {
		t.Fatalf("committed golden missing (failure, not a skip): %v", err)
	}
	lines := strings.Split(strings.TrimSpace(string(b)), "\n")
	if len(lines) != 2001 {
		t.Fatalf("hcomp_default_tiles.txt has %d lines, want 2001", len(lines))
	}
	seen := map[int64]bool{}
	bad := 0
	for _, ln := range lines {
		var n1, n2, z1, z2 int64
		if _, err := fmt.Sscan(ln, &n1, &n2, &z1, &z2); err != nil {
			t.Fatalf("parse %q: %v", ln, err)
		}
		seen[z2] = true
		got, err := resolveTileShape("test", compress.HCOMPRESS1, []int64{n1, n2}, nil)
		if err != nil {
			bad++
			if bad < 5 {
				t.Errorf("%dx%d: our writer refuses cfitsio's default tiling %dx%d: %v", n1, n2, z1, z2, err)
			}
			continue
		}
		if got[0] != z1 || got[1] != z2 {
			bad++
			if bad < 5 {
				t.Errorf("%dx%d: default tile %v, cfitsio %dx%d", n1, n2, got, z1, z2)
			}
		}
	}
	if bad > 0 {
		t.Errorf("%d/%d image heights differ from cfitsio", bad, len(lines))
	}
	for _, h := range []int64{16, 17, 18, 20, 22, 24, 26, 28, 30} {
		if !seen[h] {
			t.Errorf("table never reaches tile height %d", h)
		}
	}
}
