package fits

import (
	"bytes"
	"context"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestCompressedHDUHCompressMetadata verifies that the HCOMPRESS metadata
// is accessible through the normal HDU surface (algorithm name, BITPIX,
// shape). Pixel correctness is verified bit-exact against cfitsio by
// TestCompressedFixtures/hcompress_* and TestCompressedCorpusRead.
func TestCompressedHDUHCompressMetadata(t *testing.T) {
	f, err := Open(filepath.Join("testdata", "comp_hcompress_i16.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	cimg, ok := h.(*CompressedImageHDU)
	if !ok {
		t.Fatalf("HDU(1) is %T, want *CompressedImageHDU", h)
	}
	if cimg.CompressionType() != "HCOMPRESS_1" {
		t.Fatalf("algorithm: %q", cimg.CompressionType())
	}
	if cimg.BITPIX() != 16 {
		t.Fatalf("BITPIX: %d", cimg.BITPIX())
	}
}

// TestCompressedHDUCompressedFlag verifies the HDU interface's new
// Compressed() and CompressionType() methods return the right values.
func TestCompressedHDUCompressedFlag(t *testing.T) {
	// Compressed case.
	f, _ := Open(filepath.Join("testdata", "comp_rice_i16.fits"))
	defer f.Close()
	h, _ := f.HDU(1)
	if !h.Compressed() {
		t.Error("compressed HDU.Compressed() should be true")
	}
	if h.CompressionType() != "RICE_1" {
		t.Errorf("CompressionType: %q", h.CompressionType())
	}
	// Uncompressed case.
	fr, _ := Open(filepath.Join("testdata", "comp_raw_rice_i16.fits"))
	defer fr.Close()
	hr, _ := fr.HDU(1)
	if hr.Compressed() {
		t.Error("plain ImageHDU.Compressed() should be false")
	}
	if hr.CompressionType() != "" {
		t.Errorf("CompressionType on plain image: %q", hr.CompressionType())
	}
}

// editedCompressed copies a committed compressed sample to a temp file,
// replacing the header card that starts with key (8 chars, padded) by
// card (padded to 80 bytes).
func editedCompressed(t *testing.T, sample, key, card string) *CompressedImageHDU {
	t.Helper()
	b, err := os.ReadFile(testutil.TestdataPath(t, sample))
	if err != nil {
		t.Fatal(err)
	}
	k := []byte(fmt.Sprintf("%-8s=", key))
	i := bytes.Index(b[2880:], k)
	if i < 0 || i%80 != 0 {
		t.Fatalf("%s: card %s not found in HDU 1", sample, key)
	}
	copy(b[2880+i:2880+i+80], fmt.Sprintf("%-80s", card))
	p := filepath.Join(t.TempDir(), "edited.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := Open(p)
	if err != nil {
		t.Fatal(err)
	}
	t.Cleanup(func() { f.Close() })
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	return h.(*CompressedImageHDU)
}

// Unproven or invalid layouts must fail loudly, never decode to guesses.
func TestCompressedReadRejects(t *testing.T) {
	cases := []struct {
		name, sample, key, card, want string
	}{
		// cfitsio refuses an unknown ZQUANTIZ (DATA_DECOMPRESSION_ERR).
		{"zquantiz", "comp_fz_f32_rice_sd1.fits", "ZQUANTIZ", "ZQUANTIZ= 'SUBTRACTIVE_DITHER_9'", "unknown ZQUANTIZ"},
		// RICE_1 tiles of an unquantized float image other than 4-byte
		// tiles (read: comp_fz_f32_rice_zqnone, comp_fz_f64_rice_zqnone):
		// cfitsio reads 1-/2-byte ones as floats past the end of the tile
		// (no defined value); refused at -32 and -64.
		{"rice_unquantized_f64_bp2", "comp_fz_f64_rice_zqnone.fits", "ZVAL2", "ZVAL2   =                    2", "RICE_1 BYTEPIX 2 tiles of an unquantized ZBITPIX -64 image are not supported"},
		{"rice_unquantized_bp2", "comp_fz_f32_rice_zqnone.fits", "ZVAL2", "ZVAL2   =                    2", "RICE_1 BYTEPIX 2 tiles of an unquantized ZBITPIX -32 image are not supported"},
		{"rice_unquantized_bp1", "comp_fz_f32_rice_zqnone.fits", "ZVAL2", "ZVAL2   =                    1", "RICE_1 BYTEPIX 1 tiles of an unquantized ZBITPIX -32 image are not supported"},
		// HCOMPRESS_1 with ZQUANTIZ 'NONE' (fpack refuses '-q 0' with it;
		// no reference reading): refused, not decoded as float bits.
		{"hcomp_unquantized", "comp_fz_f32_hcomp_q.fits", "ZQUANTIZ", "ZQUANTIZ= 'NONE'", "HCOMPRESS_1 tiles of an unquantized ZBITPIX -32 image are not supported"},
		// ZVAL2 gone from a file whose stream is 2-byte Rice: cfitsio then
		// decodes it as 4-byte Rice (BYTEPIX 4, as for the real
		// comp_fz_real_m13_rice) and gets values int16 cannot hold; a read
		// as the image type fails (status 412) and so does ours.
		{"bytepix_default", "comp_fz_u16_rice.fits", "ZVAL2", "COMMENT   ZVAL2 removed", "does not fit ZBITPIX 16"},
		// RICE_1 / HCOMPRESS_1 need ZVAL1 (cfitsio: required keyword).
		{"rice_zval1", "comp_fz_u16_rice.fits", "ZVAL1", "COMMENT   ZVAL1 removed", "ZVAL1"},
		{"hcomp_zval1", "comp_fz_i16_hcomp_s4.fits", "ZVAL1", "COMMENT   ZVAL1 removed", "ZVAL1"},
		// Only one of the ZSCALE/ZZERO columns: cfitsio's handling is
		// inconsistent (mixes in the keywords or fails); refuse.
		{"zscale_col_alone", "comp_fz_f32_rice_sd1.fits", "TTYPE3", "TTYPE3  = 'NOTZZERO'", "both ZSCALE and ZZERO"},
		// ZSCALE on an integer image: cfitsio would use it in place of
		// BSCALE; no writer makes it, refuse.
		{"zscale_int", "comp_fz_i16blank_rice.fits", "BLANK", "ZSCALE  =                   2.", "integer image"},
		// Lossless float with a ZSCALE keyword but no ZZERO.
		{"zscale_no_zzero", "comp_fz_f32_lossless_zscalekw.fits", "ZZERO", "COMMENT   ZZERO removed", "without ZZERO"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			h := editedCompressed(t, c.sample, c.key, c.card)
			_, err := ReadPixelsCompressed[float64](h)
			if err == nil || !strings.Contains(err.Error(), c.want) {
				t.Fatalf("err = %v, want one mentioning %q", err, c.want)
			}
			if _, err := ReadFloat32(h); err == nil || !strings.Contains(err.Error(), c.want) {
				t.Fatalf("ReadFloat32: err = %v, want one mentioning %q", err, c.want)
			}
		})
	}
}

// ReadPixelsCompressedContext honours cancellation before reading.
func TestReadPixelsCompressedContextCanceled(t *testing.T) {
	f, err := Open(testutil.TestdataPath(t, "comp_fz_u16_rice.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := ReadPixelsCompressedContext[uint16](ctx, h.(*CompressedImageHDU)); !errors.Is(err, context.Canceled) {
		t.Fatalf("err = %v, want context.Canceled", err)
	}
}

// The context is consulted before the read starts, before the HDU's
// metadata is even parsed: a canceled context on an HDU whose Z* keywords
// are invalid returns context.Canceled, not the metadata error (the
// between-tiles check alone would return the metadata error).
func TestReadPixelsCompressedContextCheckedBeforeRead(t *testing.T) {
	h := editedCompressed(t, "comp_fz_f32_rice_sd1.fits", "ZQUANTIZ", "ZQUANTIZ= 'SUBTRACTIVE_DITHER_9'")
	ctx, cancel := context.WithCancel(context.Background())
	cancel()
	if _, err := ReadPixelsCompressedContext[float32](ctx, h); !errors.Is(err, context.Canceled) {
		t.Fatalf("err = %v, want context.Canceled from the check before the read", err)
	}
	if _, err := ReadPixelsCompressedContext[float32](context.Background(), h); err == nil || !strings.Contains(err.Error(), "unknown ZQUANTIZ") {
		t.Fatalf("uncanceled read: err = %v, want the metadata error", err)
	}
}

// countdownCtx reports no error for its first n Err calls, then
// context.Canceled: a context canceled while the read is in progress.
type countdownCtx struct {
	context.Context
	n, calls int
}

func (c *countdownCtx) Err() error {
	c.calls++
	if c.calls > c.n {
		return context.Canceled
	}
	return nil
}

// ReadPixelsCompressedContext also checks the context between tiles: a
// cancellation after the first tile stops the read.
func TestReadPixelsCompressedContextCanceledBetweenTiles(t *testing.T) {
	f, err := Open(testutil.TestdataPath(t, "comp_fz_u16_rice.fits")) // 24 row tiles
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	// Call 1: before the read; call 2: before tile 0; call 3: before tile 1.
	ctx := &countdownCtx{Context: context.Background(), n: 2}
	if _, err := ReadPixelsCompressedContext[uint16](ctx, h.(*CompressedImageHDU)); !errors.Is(err, context.Canceled) {
		t.Fatalf("err = %v, want context.Canceled from the tile loop", err)
	}
	if ctx.calls != 3 {
		t.Fatalf("ctx.Err called %d times, want 3 (stop at tile 1)", ctx.calls)
	}
}

// TestCompressedHeaderParseErrorNoPanic: a tile-compressed HDU whose header
// does not parse gives errors, not panics. Open frames every header by its
// END card, so no committed file reaches this (the lazy parse in
// hduRecord.loadHeader fails only without END or on a length that is not
// a multiple of 80); the record is built here with no END card. Header()
// returns an empty header, CompressionType() "", the logical accessors
// zero values, and Validate, HeaderOf and every pixel reader the parse
// error (Header() and metadata() used to panic).
func TestCompressedHeaderParseErrorNoPanic(t *testing.T) {
	rec := &hduRecord{index: 1, kind: kindCompressed, xtension: "BINTABLE",
		rawHeader: bytes.Repeat([]byte("COMMENT   no END card"+strings.Repeat(" ", 59)), 36)}
	h := &CompressedImageHDU{rec: rec, tbl: &BinaryTableHDU{rec: rec}}
	_, herr := HeaderOf(h)
	if herr == nil {
		t.Fatalf("HeaderOf: no error for a header without END")
	}
	func() {
		defer func() {
			if r := recover(); r != nil {
				t.Fatalf("accessor panicked: %v", r)
			}
		}()
		if hdr := h.Header(); hdr == nil || len(hdr.Cards()) != 0 {
			t.Errorf("Header() = %v, want an empty header", hdr)
		}
		if ct := h.CompressionType(); ct != "" {
			t.Errorf("CompressionType() = %q, want \"\"", ct)
		}
		if h.BITPIX() != 0 || h.NAXIS() != 0 || h.Shape() != nil {
			t.Errorf("BITPIX/NAXIS/Shape = %d/%d/%v, want zero values", h.BITPIX(), h.NAXIS(), h.Shape())
		}
	}()
	if err := h.Validate(); err == nil || !strings.Contains(err.Error(), herr.Error()) {
		t.Errorf("Validate: %v, want the header parse error %v", err, herr)
	}
	if _, err := ReadPixelsCompressed[float64](h); err == nil || !strings.Contains(err.Error(), herr.Error()) {
		t.Errorf("ReadPixelsCompressed: err %v", err)
	}
	if _, _, err := ReadPixelsCompressedMasked[uint16](h); err == nil || !strings.Contains(err.Error(), herr.Error()) {
		t.Errorf("ReadPixelsCompressedMasked: err %v", err)
	}
	if _, err := ReadFloat32(h); err == nil || !strings.Contains(err.Error(), herr.Error()) {
		t.Errorf("ReadFloat32: err %v", err)
	}
}
