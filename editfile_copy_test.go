package fits

import (
	"bytes"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// editCopyScenarios mirror SCENARIOS in testdata/cfcopy/Makefile: copy the
// listed HDUs of a committed sample, in order, through EditFile (skipping
// the rest). The truth is testdata/cfcopy/golden/<name>.fits, written by
// cfitsio's fits_copy_hdu with the same order.
var editCopyScenarios = []struct {
	name   string
	sample string
	order  []int
}{
	// IMAGE extension to position 0 (becomes the primary), primary to
	// position 1 (becomes an IMAGE extension), every other HDU kept.
	{"testprog_2_0_1-7", "testprog.std", []int{2, 0, 1, 3, 4, 5, 6, 7}},
	// Primary skipped: a binary table lands at position 0, so an empty
	// primary array must be put in front of it.
	{"testprog_skip0", "testprog.std", []int{1, 2, 3, 4, 5, 6, 7}},
	// Two IMAGE extensions only; the first becomes the primary.
	{"testprog_6_2", "testprog.std", []int{6, 2}},
	// The primary written twice; the second copy becomes an extension.
	{"testprog_0_3_0", "testprog.std", []int{0, 3, 0}},
	{"iter_image_0_0", "iter_image.fit", []int{0, 0}},
	// A tile-compressed image (a BINTABLE) at position 0.
	{"vari_skip0", "vari.fits", []int{1}},
	// cfcopy/reserve.fits (made by cfitsio with reserved header space, see
	// the Makefile): its IMAGE extension becomes the primary and its
	// primary an extension, both with reserved blank cards in the input;
	// the third HDU keeps its form and its reserved space.
	{"reserve_1_0_2", "cfcopy/reserve.fits", []int{1, 0, 2}},
	// cfcopy/comments.fits (made by cfitsio, see mkcomments.c): COMMENT
	// cards that do and do not match the 58- and 47-character prefixes by
	// which fits_copy_hdu drops the standard FITS comments when a primary
	// becomes an extension, and the reverse conversion with the standard
	// comments already present.
	{"comments_1_0", "cfcopy/comments.fits", []int{1, 0}},
	// An fpack-written tile-compressed image (uint16 RICE_1) at position 0
	// (an empty primary is put in front), and moved in front of its own
	// primary (which becomes an IMAGE extension).
	{"comp_skip0", "comp_fz_u16_rice.fits", []int{1}},
	{"comp_1_0", "comp_fz_u16_rice.fits", []int{1, 0}},
	// Data cubes (NAXIS 3, 4 and 6: cfcopy/cube*.fits, funpack -C of the
	// compressed cube samples, see the cfcopy Makefile): the primary copied to position 0 and to
	// position 1 (an IMAGE extension with NAXIS 3..6), and that extension
	// (from cfitsio's golden) moved back to position 0.
	{"cube3_0_0", "cfcopy/cube3.fits", []int{0, 0}},
	{"cube4_0_0", "cfcopy/cube4.fits", []int{0, 0}},
	{"cube6_0_0", "cfcopy/cube6.fits", []int{0, 0}},
	{"cube3_1_0", "cfcopy/golden/cube3_0_0.fits", []int{1, 0}},
	// A NAXIS = 0 primary to position 1 (an IMAGE extension with NAXIS 0),
	// and that extension (from cfitsio's golden) to position 0.
	{"naxis0_0_0", "vari.fits", []int{0, 0}},
	{"naxis0_1_0", "cfcopy/golden/naxis0_0_0.fits", []int{1, 0}},
}

func cfcopyGolden(t *testing.T, name string) string {
	t.Helper()
	p := testutil.TestdataPath(t, "cfcopy", "golden", name+".fits")
	testutil.RequireFile(t, p)
	return p
}

// hduParts splits a FITS file into per-HDU header cards (blank cards
// between the last keyword and END dropped: reserved space and cfitsio's
// blank fill are layout, checked separately by checkReservedSpace) and data
// bytes (including block padding).
func hduParts(t *testing.T, path string) (cards [][]string, data [][]byte) {
	t.Helper()
	f, err := Open(path)
	if err != nil {
		t.Fatalf("Open %s: %v", path, err)
	}
	defer f.Close()
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	for _, r := range f.hdus {
		cs, err := splitCards(r.rawHeader)
		if err != nil {
			t.Fatalf("%s HDU %d: %v", path, r.index, err)
		}
		for len(cs) > 0 && strings.TrimSpace(cs[len(cs)-1]) == "" {
			cs = cs[:len(cs)-1]
		}
		cards = append(cards, cs)
		data = append(data, raw[r.dataStart:r.paddedEnd])
	}
	return cards, data
}

func runEditCopy(t *testing.T, path string, order []int) error {
	t.Helper()
	return EditFile(path, func(in *File, out *Writer) error {
		used := map[int]bool{}
		for _, i := range order {
			used[i] = true
		}
		for i := range in.NumHDU() {
			if !used[i] {
				h, err := in.HDU(i)
				if err != nil {
					return err
				}
				if err := out.SkipHDU(h); err != nil {
					return err
				}
			}
		}
		for _, i := range order {
			h, err := in.HDU(i)
			if err != nil {
				return err
			}
			if err := out.CopyHDU(h); err != nil {
				return err
			}
		}
		return nil
	})
}

// TestEditFileCopyMatchesCfitsioCopyHDU proves that EditFile's CopyHDU and
// SkipHDU convert SIMPLE <-> XTENSION exactly like cfitsio when an HDU
// changes position: every header card and every data byte of every output
// HDU equals cfitsio's fits_copy_hdu output, cfitsio reads both files the
// same, and fitsverify finds no more in ours than in cfitsio's.
func TestEditFileCopyMatchesCfitsioCopyHDU(t *testing.T) {
	for _, sc := range editCopyScenarios {
		t.Run(sc.name, func(t *testing.T) {
			path, _ := copySample(t, sc.sample)
			if err := runEditCopy(t, path, sc.order); err != nil {
				t.Fatalf("EditFile: %v", err)
			}
			golden := cfcopyGolden(t, sc.name)
			gotCards, gotData := hduParts(t, path)
			wantCards, wantData := hduParts(t, golden)
			if len(gotCards) != len(wantCards) {
				t.Fatalf("%d HDUs, cfitsio wrote %d", len(gotCards), len(wantCards))
			}
			for i := range wantCards {
				if !reflect.DeepEqual(gotCards[i], wantCards[i]) {
					t.Errorf("HDU %d header differs from cfitsio's:\n got  %q\n want %q", i, firstCards(gotCards[i]), firstCards(wantCards[i]))
				}
				if !bytes.Equal(gotData[i], wantData[i]) {
					t.Errorf("HDU %d data differs from cfitsio's (len %d vs %d)", i, len(gotData[i]), len(wantData[i]))
				}
			}
			checkReservedSpace(t, path, testutil.TestdataPath(t, sc.sample), sc.order, wantCards)

			got := testutil.CfrefDump(t, path)
			want := testutil.CfrefDump(t, golden)
			compareDumpsIgnoringLayout(t, got, want)
			gv, wv := testutil.Fitsverify(t, path), testutil.Fitsverify(t, golden)
			if gv.Errors != wv.Errors || gv.Warnings != wv.Warnings {
				t.Errorf("fitsverify: ours %d warnings / %d errors, cfitsio's %d / %d\n%s", gv.Warnings, gv.Errors, wv.Warnings, wv.Errors, gv.Report)
			}
			checkCopiedCompressedPixels(t, path, sc.sample, sc.order)
		})
	}
}

// checkCopiedCompressedPixels: every tile-compressed HDU of the EditFile
// output decodes (ReadPixelsCompressed) to the physical values cfitsio
// reads from the same HDU of the source sample (its committed cfref
// golden, testdata/cfref/golden for top-level samples).
func checkCopiedCompressedPixels(t *testing.T, outPath, sample string, order []int) {
	t.Helper()
	if strings.Contains(sample, "/") {
		return // no cfref golden for scenario inputs under cfcopy/
	}
	src := testutil.LoadCfrefGolden(t, sample)
	out, err := Open(outPath)
	if err != nil {
		t.Fatal(err)
	}
	defer out.Close()
	shift := out.NumHDU() - len(order)
	for j, si := range order {
		h, err := out.HDU(j + shift)
		if err != nil {
			t.Fatal(err)
		}
		c, ok := h.(*CompressedImageHDU)
		if !ok {
			continue
		}
		img := src.HDUs[si].Image
		if img == nil || !img.Compressed || img.Phys.Len() == 0 {
			t.Fatalf("source HDU %d has no compressed-image golden", si)
		}
		got, err := ReadPixelsCompressed[float64](c)
		if err != nil {
			t.Fatalf("output HDU %d: %v", j+shift, err)
		}
		if len(got) != img.Phys.Len() {
			t.Fatalf("output HDU %d: %d pixels, cfitsio %d", j+shift, len(got), img.Phys.Len())
		}
		for k, w := range img.Phys.Values(t) {
			if !w.Valid || got[k] != w.V {
				t.Fatalf("output HDU %d pixel %d = %v, cfitsio %+v", j+shift, k, got[k], w)
			}
		}
	}
}

// checkReservedSpace checks what happens to blank cards before END. A
// converted header is rebuilt from its keywords: its untrimmed cards must
// be exactly cfitsio's keyword list (wantCards, which has no blank tail;
// cfitsio's ffghsp keyword count stops at the last non-blank card, so
// fits_copy_hdu never copies reserved space). An HDU that keeps its form
// must keep its header bytes, reserved space included.
func checkReservedSpace(t *testing.T, outPath, srcPath string, order []int, wantCards [][]string) {
	t.Helper()
	out, err := Open(outPath)
	if err != nil {
		t.Fatal(err)
	}
	defer out.Close()
	src, err := Open(srcPath)
	if err != nil {
		t.Fatal(err)
	}
	defer src.Close()
	shift := len(out.hdus) - len(order) // 1 when an empty primary was put in front
	for j, si := range order {
		oi := j + shift
		o, s := out.hdus[oi], src.hdus[si]
		converted := (oi == 0) != (si == 0)
		if !converted {
			if !bytes.Equal(o.rawHeader, s.rawHeader) {
				t.Errorf("HDU %d (input HDU %d) keeps its form but its header bytes changed", oi, si)
			}
			continue
		}
		cs, err := splitCards(o.rawHeader)
		if err != nil {
			t.Fatal(err)
		}
		if !reflect.DeepEqual(cs, wantCards[oi]) {
			t.Errorf("converted HDU %d (input HDU %d) has %d cards before END, cfitsio's keyword list has %d (reserved blank cards carried over?)", oi, si, len(cs), len(wantCards[oi]))
		}
	}
}

func firstCards(cs []string) []string {
	if len(cs) > 10 {
		return cs[:10]
	}
	return cs
}

// compareDumpsIgnoringLayout compares two cfitsio dumps HDU by HDU except
// for byte offsets (reserved header space may differ).
func compareDumpsIgnoringLayout(t *testing.T, got, want *testutil.CfrefFile) {
	t.Helper()
	if got.NHDU != want.NHDU {
		t.Fatalf("cfitsio reads %d HDUs, want %d", got.NHDU, want.NHDU)
	}
	for i := range want.HDUs {
		g, w := got.HDUs[i], want.HDUs[i]
		g.HeaderStart, g.DataStart, g.DataEnd = 0, 0, 0
		w.HeaderStart, w.DataStart, w.DataEnd = 0, 0, 0
		if !reflect.DeepEqual(g, w) {
			t.Errorf("HDU %d: cfitsio reads our output differently from the expected file (type %s/%s, %d/%d cards)", i, g.Type, w.Type, len(g.Cards), len(w.Cards))
		}
	}
}

// A new image written first through the Writer pushes the copied primary
// to position 1, where it must become an IMAGE extension: cfitsio must read
// it with the same cards and pixels as fits_copy_hdu's conversion (the
// third HDU of the testprog_0_3_0 golden is that same conversion).
func TestEditFileWriteImageThenCopyPrimary(t *testing.T) {
	path, _ := copySample(t, "testprog.std")
	pix := []int16{1, -2, 3, -4, 5, -6}
	err := EditFile(path, func(in *File, out *Writer) error {
		if _, err := WriteImage(out, nil, []int64{3, 2}, pix); err != nil {
			return err
		}
		p, err := in.Primary()
		if err != nil {
			return err
		}
		return out.CopyHDU(p)
	})
	if err != nil {
		t.Fatal(err)
	}
	gotCards, gotData := hduParts(t, path)
	wantCards, wantData := hduParts(t, cfcopyGolden(t, "testprog_0_3_0"))
	if len(gotCards) != 2 {
		t.Fatalf("%d HDUs, want 2", len(gotCards))
	}
	if !reflect.DeepEqual(gotCards[1], wantCards[2]) {
		t.Errorf("converted primary differs from cfitsio's:\n got  %q\n want %q", firstCards(gotCards[1]), firstCards(wantCards[2]))
	}
	if !bytes.Equal(gotData[1], wantData[2]) {
		t.Error("converted primary data differs from cfitsio's")
	}

	got := testutil.CfrefDump(t, path)
	want := testutil.CfrefDump(t, cfcopyGolden(t, "testprog_0_3_0"))
	if got.NHDU != 2 {
		t.Fatalf("cfitsio reads %d HDUs", got.NHDU)
	}
	img := got.HDUs[0].Image
	if img == nil || img.Bitpix != 16 || !reflect.DeepEqual(img.Naxes, []int64{3, 2}) {
		t.Fatalf("cfitsio HDU 0 image = %+v", img)
	}
	for i, v := range pix {
		if !img.I64.Values(t)[i].Valid || img.I64.Values(t)[i].V != int64(v) {
			t.Errorf("cfitsio pixel %d = %+v, want %d", i, img.I64.Values(t)[i], v)
		}
	}
	g, w := got.HDUs[1], want.HDUs[2]
	g.Index, g.HeaderStart, g.DataStart, g.DataEnd = 0, 0, 0, 0
	w.Index, w.HeaderStart, w.DataStart, w.DataEnd = 0, 0, 0, 0
	if !reflect.DeepEqual(g, w) {
		t.Error("cfitsio reads the converted primary differently from fits_copy_hdu's conversion")
	}
	if v := testutil.Fitsverify(t, path); v.Errors != 0 {
		t.Errorf("fitsverify: %d errors\n%s", v.Errors, v.Report)
	}
}

// New images written through the Writer after copied HDUs must be IMAGE
// extensions (XTENSION = 'IMAGE', PCOUNT, GCOUNT), not a second SIMPLE
// header. testprog.std HDU 0 is copied, a 3x2 int16 image written, HDU 2
// (an IMAGE extension) copied, and a 1-D float32 image written last.
// cfitsio must read four HDUs: the copied ones exactly as in the committed
// golden of testprog.std, the new ones as IMAGE HDUs with the pixels that
// were written; fitsverify must find no errors.
func TestEditFileWriteImageAfterCopy(t *testing.T) {
	path, _ := copySample(t, "testprog.std")
	pix16 := []int16{1, -2, 3, -4, 5, -32768}
	pix32 := []float32{0.5, -1.25, 3e10}
	err := EditFile(path, func(in *File, out *Writer) error {
		h0, err := in.HDU(0)
		if err != nil {
			return err
		}
		if err := out.CopyHDU(h0); err != nil {
			return err
		}
		if _, err := WriteImage(out, nil, []int64{3, 2}, pix16); err != nil {
			return err
		}
		h2, err := in.HDU(2)
		if err != nil {
			return err
		}
		if err := out.CopyHDU(h2); err != nil {
			return err
		}
		_, err = WriteImage(out, nil, []int64{3}, pix32)
		return err
	})
	if err != nil {
		t.Fatal(err)
	}
	got := testutil.CfrefDump(t, path)
	golden := testutil.LoadCfrefGolden(t, "testprog.std")
	if got.NHDU != 4 {
		t.Fatalf("cfitsio reads %d HDUs, want 4", got.NHDU)
	}
	for _, c := range []struct{ out, in int }{{0, 0}, {2, 2}} {
		g, w := got.HDUs[c.out], golden.HDUs[c.in]
		g.Index, g.HeaderStart, g.DataStart, g.DataEnd = 0, 0, 0, 0
		w.Index, w.HeaderStart, w.DataStart, w.DataEnd = 0, 0, 0, 0
		if !reflect.DeepEqual(g, w) {
			t.Errorf("output HDU %d: cfitsio reads it differently from testprog.std HDU %d", c.out, c.in)
		}
	}
	for _, i := range []int{1, 3} {
		h := got.HDUs[i]
		if h.Type != "IMAGE" || len(h.Cards) < 6 || !strings.HasPrefix(string(h.Cards[0].Raw), "XTENSION= 'IMAGE   '") {
			t.Fatalf("output HDU %d is not an IMAGE extension: type %s, first card %q", i, h.Type, h.Cards[0].Raw)
		}
		if p, g := h.Card("PCOUNT"), h.Card("GCOUNT"); p == nil || p.Int == nil || *p.Int != 0 || g == nil || g.Int == nil || *g.Int != 1 {
			t.Errorf("output HDU %d: PCOUNT/GCOUNT = %+v / %+v, want 0 / 1", i, p, g)
		}
		if h.Card("SIMPLE") != nil {
			t.Errorf("output HDU %d has a SIMPLE card", i)
		}
	}
	img := got.HDUs[1].Image
	if img == nil || img.Bitpix != 16 || !reflect.DeepEqual(img.Naxes, []int64{3, 2}) {
		t.Fatalf("cfitsio HDU 1 image = %+v", img)
	}
	for i, v := range pix16 {
		if !img.I64.Values(t)[i].Valid || img.I64.Values(t)[i].V != int64(v) {
			t.Errorf("cfitsio HDU 1 pixel %d = %+v, want %d", i, img.I64.Values(t)[i], v)
		}
	}
	img = got.HDUs[3].Image
	if img == nil || img.Bitpix != -32 || !reflect.DeepEqual(img.Naxes, []int64{3}) {
		t.Fatalf("cfitsio HDU 3 image = %+v", img)
	}
	for i, v := range pix32 {
		if !img.Phys.Values(t)[i].Valid || img.Phys.Values(t)[i].V != float64(v) {
			t.Errorf("cfitsio HDU 3 pixel %d = %+v, want %g", i, img.Phys.Values(t)[i], v)
		}
	}
	if v := testutil.Fitsverify(t, path); v.Errors != 0 {
		t.Errorf("fitsverify: %d errors\n%s", v.Errors, v.Report)
	}
}

// convertForPosition refuses a header whose NAXISn cards are out of order
// (FITS 4.0 sec. 4.4.1.1 fixes their order), rather than guessing which
// cards are the mandatory ones.
func TestConvertForPositionRefusesBadAxisOrder(t *testing.T) {
	for _, pos := range []int{0, 1} {
		first := "SIMPLE  =                    T"
		if pos == 0 {
			first = "XTENSION= 'IMAGE   '"
		}
		b := packHeader([]string{
			first, "BITPIX  =                    8", "NAXIS   =                    2",
			"NAXIS2  =                    1", "NAXIS1  =                    1",
		})
		if _, _, err := convertForPosition(b, pos); err == nil || !strings.Contains(err.Error(), "want NAXIS1") {
			t.Errorf("position %d: convertForPosition err = %v, want the NAXIS order refused", pos, err)
		}
	}
}

// A random-groups primary cannot become an extension: refused, source kept.
func TestEditFileRefusesRandomGroupsAsExtension(t *testing.T) {
	cards := []string{
		"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    2",
		"NAXIS1  =                    0", "NAXIS2  =                    1", "GROUPS  =                    T",
		"PCOUNT  =                    0", "GCOUNT  =                    1",
	}
	b := append(packHeader(cards), make([]byte, 2880)...)
	if _, _, err := convertForPosition(b[:2880], 1); err == nil || !strings.Contains(err.Error(), "random-groups") {
		t.Fatalf("convertForPosition(random groups, 1) err = %v", err)
	}
	path := filepath.Join(t.TempDir(), "rg.fits")
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	err := EditFile(path, func(in *File, out *Writer) error {
		if _, err := WriteImage(out, nil, []int64{1}, []uint8{0}); err != nil {
			return err
		}
		p, err := in.HDU(0)
		if err != nil {
			return err
		}
		return out.CopyHDU(p)
	})
	if err == nil {
		t.Fatal("EditFile wrote a random-groups primary as an extension")
	}
	if !bytes.Equal(mustRead(t, path), b) {
		t.Fatal("source modified")
	}
}
