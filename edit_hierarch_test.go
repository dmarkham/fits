package fits

import (
	"bytes"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// u8Image is a primary BITPIX 8 4x2 image with pixels 1..8 and fill
// filler cards: with 29 the header (35 cards with END) keeps one block
// after an edit adds a card; with 30 that edit moves END to a second
// block. The file is input only: the truth is cfitsio's reading of what
// Flush writes.
func u8Image(t *testing.T, fill int) (string, []byte) {
	t.Helper()
	cards := []string{
		"SIMPLE  =                    T",
		"BITPIX  =                    8",
		"NAXIS   =                    2",
		"NAXIS1  =                    4",
		"NAXIS2  =                    2",
	}
	for i := range fill {
		cards = append(cards, fmt.Sprintf("FILL%04d=                    %d", i, i))
	}
	b := fitsBytes(struct {
		cards []string
		data  []byte
	}{cards, []byte{1, 2, 3, 4, 5, 6, 7, 8}})
	p := filepath.Join(t.TempDir(), "u8.fits")
	if err := os.WriteFile(p, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return p, b
}

// cfitsio (ffgknm) names a HIERARCH card by the text between "HIERARCH "
// and '=', so "HIERARCH BSCALE = 2" is the keyword BSCALE for every
// cfitsio routine: ffgphd applies it as scaling, ffgphd and ffbinit end
// the header at "HIERARCH END = 1", and fits_read_key finds it by the short
// name. astropy names it BSCALE (or END) too. This package's parser keys
// it "HIERARCH BSCALE" and ignores it, so writing such a card makes the
// readers disagree on the data. Flush must refuse every HIERARCH card
// whose path is itself a standard keyword name (at most 8 characters of
// A-Z 0-9 - _), leave the file byte-identical, and the refusal must be
// needed: cfitsio reads the card the encoder would have written as that
// short keyword. HIERARCH cards with a real path (blanks, or longer than 8
// characters) are still written, and cfitsio reads the file as we do.
//
// Truth: the cfref harness (cfitsio 4.7.0) on the flushed file, and on the
// header the encoder would have written.
func TestFlushRefusesHierarchAliasOfStandardKeyword(t *testing.T) {
	cases := []struct {
		key     string
		value   any
		refused bool
		fill    int // filler cards (default 29)
	}{
		{"HIERARCH BSCALE", 2.0, true, 0},
		{"HIERARCH BZERO", 100.0, true, 0},
		{"HIERARCH BLANK", int64(3), true, 0},
		{"HIERARCH END", int64(1), true, 0},
		// END moves to block 2: cfitsio would take the data unit to start
		// at 2880, inside the header.
		{"HIERARCH END", int64(1), true, 30},
		{"HIERARCH NAXIS1", int64(8), true, 0},
		{"HIERARCH EXTNAME", "SCI", true, 0},
		{"HIERARCH OBJECT", "m31", true, 0},
		{"HIERARCH TSCAL1", 2.0, true, 0},
		{"HIERARCH ESO DET BSCALE", 2.0, false, 0},
		{"HIERARCH LONGKEYWORD", "x", false, 0},
		{"HIERARCH BSCALE20X", 2.0, false, 0},
	}
	for _, tc := range cases {
		if tc.fill == 0 {
			tc.fill = 29
		}
		t.Run(fmt.Sprintf("%s/fill=%d", strings.ReplaceAll(tc.key, " ", "_"), tc.fill), func(t *testing.T) {
			path, orig := u8Image(t, tc.fill)
			f, err := OpenForEdit(path)
			if err != nil {
				t.Fatal(err)
			}
			h := hduHeader(t, f, 0)
			if err := h.Set(tc.key, tc.value, "edit"); err != nil {
				t.Fatalf("Set: %v", err)
			}
			encoded, encErr := header.Encode(h)
			ferr := f.Flush()
			f.Close()
			if !tc.refused {
				if ferr != nil {
					t.Fatalf("Flush: %v", ferr)
				}
				requireCfitsioReadsAsWe(t, path)
				return
			}
			if ferr == nil {
				// Show what the readers now disagree on.
				requireCfitsioReadsAsWe(t, path)
				t.Fatalf("Flush wrote %s; cfitsio names it %s", tc.key, strings.TrimPrefix(tc.key, "HIERARCH "))
			}
			if !strings.Contains(ferr.Error(), "HIERARCH") {
				t.Errorf("Flush error %q does not name the HIERARCH card", ferr)
			}
			if !bytes.Equal(mustRead(t, path), orig) {
				t.Fatal("refused Flush changed the file")
			}
			if encErr != nil {
				t.Fatalf("Encode: %v", encErr)
			}
			// The refusal is needed: cfitsio reads the card the encoder
			// would have written as the short keyword (or, for END, ends
			// the header there).
			short := strings.TrimPrefix(tc.key, "HIERARCH ")
			would := filepath.Join(t.TempDir(), "would.fits")
			if err := os.WriteFile(would, append(append([]byte(nil), encoded...), orig[2880:]...), 0o644); err != nil {
				t.Fatal(err)
			}
			cf := testutil.CfrefDump(t, would)
			if short == "END" {
				if n := cf.HDUs[0].NKeys; n >= len(h.Cards()) {
					t.Fatalf("cfitsio reads %d keywords past HIERARCH END (header has %d cards); the refusal is not needed", n, len(h.Cards()))
				}
				return
			}
			for _, c := range cf.HDUs[0].Cards {
				if strings.HasPrefix(string(c.Raw), tc.key+" ") {
					if string(c.Key) != short {
						t.Fatalf("cfitsio names the encoded card %q %q, not %s; the refusal is not needed", c.Raw, c.Key, short)
					}
					return
				}
			}
			t.Fatalf("cfitsio lists no encoded %s card", tc.key)
		})
	}
}

// requireCfitsioReadsAsWe fails unless cfitsio and this package agree on
// the primary HDU of path: where its data starts and every pixel value
// (cfitsio TDOUBLE with its scaling, ours ReadPixels[float64]).
func requireCfitsioReadsAsWe(t *testing.T, path string) {
	t.Helper()
	cf := testutil.CfrefDump(t, path)
	f, err := Open(path)
	if err != nil {
		t.Fatalf("Open of the flushed file: %v (cfitsio: data at %d)", err, cf.HDUs[0].DataStart)
	}
	defer f.Close()
	img, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	px, err := ReadPixels[float64](img)
	if err != nil {
		t.Fatal(err)
	}
	if got, want := recordOf(img).dataStart, cf.HDUs[0].DataStart; got != want {
		t.Errorf("data unit at %d, cfitsio %d", got, want)
	}
	if cf.HDUs[0].Image == nil {
		t.Fatalf("cfitsio reads no image (errors %v)", cf.HDUs[0].Errors)
	}
	phys := cf.HDUs[0].Image.Phys.Values(t)
	if len(phys) != len(px) {
		t.Fatalf("%d pixels, cfitsio %d", len(px), len(phys))
	}
	for i := range px {
		if !phys[i].Valid || math.Float64bits(phys[i].V) != math.Float64bits(px[i]) {
			t.Fatalf("pixels %v, cfitsio %v", px, phys)
		}
	}
}

// Keys that Set refuses (a blank inside the name, non-ASCII and control
// bytes) can still be in a header read from disk. A Flush of another edit
// to that header must write those cards byte-identical, where they were,
// and cfitsio must read every card of the header it read before exactly as
// before (raw bytes, name, value, comment), with the new card added. Truth:
// the cfref dump (cfitsio 4.7.0) before and after.
func TestFlushKeepsCardsSetWouldRefuse(t *testing.T) {
	odd := []string{
		"A B     =                    1 / blank inside the name",
		"K\xe9Y    =                    2 / non-ASCII byte",
		"K\x01Y    =                    3 / control byte",
		"\x7fDEL   = 'x'",
		"lower   =                    4",
	}
	cards := []string{
		"SIMPLE  =                    T",
		"BITPIX  =                    8",
		"NAXIS   =                    2",
		"NAXIS1  =                    4",
		"NAXIS2  =                    2",
	}
	cards = append(cards, odd...)
	b := fitsBytes(struct {
		cards []string
		data  []byte
	}{cards, []byte{1, 2, 3, 4, 5, 6, 7, 8}})
	path := filepath.Join(t.TempDir(), "odd.fits")
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	before := testutil.CfrefDump(t, path)
	f, err := OpenForEdit(path)
	if err != nil {
		t.Fatal(err)
	}
	h := hduHeader(t, f, 0)
	for _, k := range []string{"A B", "K\xe9Y", "K\x01Y"} {
		if err := h.Set(k, int64(9), ""); err == nil {
			t.Errorf("Set(%q) accepted a key FITS does not allow", k)
		}
	}
	if err := h.Set("NEWKEY", int64(7), "added"); err != nil {
		t.Fatal(err)
	}
	if err := f.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	f.Close()
	got := mustRead(t, path)
	for i := range len(cards) {
		if !bytes.Equal(got[i*80:(i+1)*80], b[i*80:(i+1)*80]) {
			t.Errorf("card %d changed: %q -> %q", i+1, b[i*80:(i+1)*80], got[i*80:(i+1)*80])
		}
	}
	after := testutil.CfrefDump(t, path)
	bc, ac := before.HDUs[0].Cards, after.HDUs[0].Cards
	if len(ac) != len(bc)+1 {
		for _, c := range ac {
			t.Logf("after: %q", c.Raw)
		}
		t.Fatalf("cfitsio reads %d cards after the Flush, %d before; want one more", len(ac), len(bc))
	}
	for i := range bc {
		x, y := bc[i], ac[i]
		if x.Raw != y.Raw || x.Key != y.Key || x.Value != y.Value || x.Comment != y.Comment || x.Type != y.Type {
			t.Errorf("cfitsio card %d before %+v, after %+v", i+1, x, y)
		}
	}
	if c := after.HDUs[0].Card("NEWKEY"); c == nil || c.Int == nil || *c.Int != 7 {
		t.Errorf("cfitsio reads NEWKEY as %+v", c)
	}
	requireCfitsioReadsAsWe(t, path)
}

// CopyHDUWithHeader runs the same per-card guard as Flush: a header with
// HIERARCH BSCALE (which cfitsio would apply as scaling, see
// TestFlushRefusesHierarchAliasOfStandardKeyword) is refused and the
// source is left byte-identical.
func TestEditFileCopyHDUWithHeaderRefusesHierarchAlias(t *testing.T) {
	path, orig := u8Image(t, 29)
	err := EditFile(path, func(in *File, out *Writer) error {
		h, err := in.HDU(0)
		if err != nil {
			return err
		}
		hdr := h.Header().Clone()
		if err := hdr.Set("HIERARCH BSCALE", 2.0, ""); err != nil {
			return err
		}
		return out.CopyHDUWithHeader(h, hdr)
	})
	if err == nil || !strings.Contains(err.Error(), "HIERARCH") {
		t.Fatalf("EditFile with a HIERARCH BSCALE header: %v, want the HIERARCH refusal", err)
	}
	if !bytes.Equal(mustRead(t, path), orig) {
		t.Fatal("refused EditFile changed the file")
	}
}
