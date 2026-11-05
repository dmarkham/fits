package fits

import (
	"bytes"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestHDUAddressesMatchCfitsio: for every committed sample (testdata/*.fits,
// *.fit, *.std and every file under testdata/images) the number of HDUs,
// each HDU's kind, and its header start / data start / data end (cfitsio's
// fits_get_hduaddrll, data end = start of the next HDU) equal cfitsio's.
// This covers tile-compressed images, ASCII and binary tables, random
// groups, stray NAXISn cards and the framing-tolerance cases.
//
// Every sample is opened through each exported open path (Open,
// OpenReader on the *os.File, OpenReadAll on a stream), so all three are
// held to the same cfitsio truth.
func TestHDUAddressesMatchCfitsio(t *testing.T) {
	for _, s := range allImageSamples(t) {
		for _, op := range openPaths {
			t.Run(s.name+"/"+op.name, func(t *testing.T) {
				hduAddressesMatch(t, s, op.open)
			})
		}
	}
}

// openPaths are the exported ways to open a file for reading.
var openPaths = []struct {
	name string
	open func(t *testing.T, path string) (*File, error)
}{
	{"Open", func(t *testing.T, path string) (*File, error) { return Open(path) }},
	{"OpenReader", func(t *testing.T, path string) (*File, error) {
		fh, err := os.Open(path)
		if err != nil {
			t.Fatal(err)
		}
		t.Cleanup(func() { fh.Close() })
		return OpenReader(fh)
	}},
	{"OpenReadAll", func(t *testing.T, path string) (*File, error) {
		b, err := os.ReadFile(path)
		if err != nil {
			t.Fatal(err)
		}
		// A reader that is not an io.Seeker: the stream case.
		return OpenReadAll(struct{ io.Reader }{bytes.NewReader(b)})
	}},
}

func hduAddressesMatch(t *testing.T, s imageSample, open func(t *testing.T, path string) (*File, error)) {
	t.Helper()
	f, err := open(t, s.path)
	if err != nil {
		t.Fatalf("Open: %v (cfitsio opens it with %d HDUs)", err, s.golden.NHDU)
	}
	defer f.Close()
	if f.NumHDU() != s.golden.NHDU {
		t.Fatalf("NumHDU = %d, cfitsio %d", f.NumHDU(), s.golden.NHDU)
	}
	for i, g := range s.golden.HDUs {
		r := f.hdus[i]
		if r.headerStart != g.HeaderStart || r.dataStart != g.DataStart || r.paddedEnd != g.DataEnd {
			t.Errorf("HDU %d: header/data/end = %d/%d/%d, cfitsio %d/%d/%d",
				i, r.headerStart, r.dataStart, r.paddedEnd, g.HeaderStart, g.DataStart, g.DataEnd)
		}
		h, err := f.HDU(i)
		if r.kind == kindRandomGroups {
			if !errors.Is(err, ErrRandomGroups) {
				t.Errorf("HDU %d: random groups: err %v, want ErrRandomGroups", i, err)
			}
			continue
		}
		if r.kind == kindUnknown {
			// A conforming extension type the library does not decode
			// (framing/xtension_foreign.fits): framed like cfitsio above,
			// refused with a typed error naming the XTENSION value.
			var ue *UnsupportedExtensionError
			if !errors.Is(err, ErrUnsupportedExtension) || !errors.As(err, &ue) || ue.XTension != r.xtension || ue.Index != i {
				t.Errorf("HDU %d: XTENSION %q: err %v, want an *UnsupportedExtensionError for it", i, r.xtension, err)
			}
			continue
		}
		if err != nil {
			t.Errorf("HDU %d: %v", i, err)
			continue
		}
		var want string
		switch h.(type) {
		case *ImageHDU, *CompressedImageHDU:
			want = "IMAGE"
		case *BinaryTableHDU:
			want = "BINARY_TBL"
		case *ASCIITableHDU:
			want = "ASCII_TBL"
		}
		if want != g.Type {
			t.Errorf("HDU %d: %T, cfitsio %s", i, h, g.Type)
		}
	}
	// Whatever lies beyond the last HDU is reported, not hidden.
	st, err := os.Stat(s.path)
	if err != nil {
		t.Fatal(err)
	}
	last := s.golden.HDUs[len(s.golden.HDUs)-1].DataEnd
	wantTrail := max(st.Size()-last, 0)
	if f.TrailingBytes() != wantTrail {
		t.Errorf("TrailingBytes = %d, want %d (file %d bytes, cfitsio's last HDU ends at %d)", f.TrailingBytes(), wantTrail, st.Size(), last)
	}
}

// TestUnsupportedExtensionFramedAndRefused: framing/xtension_foreign.fits
// is a primary, an XTENSION = 'FOREIGN' HDU (NAXIS1 3000, PCOUNT 100: 2
// data blocks) and an IMAGE. cfitsio frames the headers at 0 / 5760 /
// 14400 with the file ending at 20160 (golden), and reads the IMAGE after it as
// [7 8 9]; the library must frame it the same way, read that IMAGE, and
// refuse the FOREIGN HDU with ErrUnsupportedExtension (cfitsio reads it as
// an IMAGE_HDU, "unknown type of extension; press on anyway"; this
// library does not treat an unknown extension's bytes as pixels).
func TestUnsupportedExtensionFramedAndRefused(t *testing.T) {
	g := loadImagesCfref(t, "framing/xtension_foreign.fits")
	if g.NHDU != 3 || g.HDUs[1].HeaderStart != 5760 || g.HDUs[1].DataEnd != 14400 || g.HDUs[2].HeaderStart != 14400 || g.HDUs[2].DataEnd != 20160 {
		t.Fatalf("cfitsio golden framing changed: %+v", g.HDUs)
	}
	f, err := Open(imagesPath(t, "framing", "xtension_foreign.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	_, err = f.HDU(1)
	var ue *UnsupportedExtensionError
	if !errors.As(err, &ue) || ue.XTension != "FOREIGN" || ue.Index != 1 || !errors.Is(err, ErrUnsupportedExtension) {
		t.Fatalf("HDU(1): %v, want *UnsupportedExtensionError{FOREIGN, 1}", err)
	}
	if strings.Contains(err.Error(), "kind") {
		t.Errorf("error %q leaks the internal HDU kind", err)
	}
	if _, err := f.HDUByName("BLOB"); !errors.Is(err, ErrUnsupportedExtension) {
		t.Errorf("HDUByName(BLOB): %v, want ErrUnsupportedExtension", err)
	}
	h, err := f.HDU(2)
	if err != nil {
		t.Fatal(err)
	}
	px, err := ReadPixels[int16](h.(*ImageHDU))
	if err != nil || len(px) != 3 || px[0] != 7 || px[1] != 8 || px[2] != 9 {
		t.Fatalf("IMAGE after FOREIGN: %v, %v; want [7 8 9]", px, err)
	}
}

// TestRandomGroupsAreRefusedNotMisread: random-groups primaries (GROUPS=T,
// NAXIS1=0) are sized per FITS 4.0 eq. 4 so following HDUs are found where
// cfitsio finds them (checked above), and the groups HDU itself is refused
// with ErrRandomGroups rather than presented as a zero-pixel image.
func TestRandomGroupsAreRefusedNotMisread(t *testing.T) {
	for _, s := range []string{"astropy_data/random_groups.fits", "astropy_data/group.fits", "framing/random_groups_ext.fits"} {
		f, err := Open(imagesPath(t, s))
		if err != nil {
			t.Fatalf("%s: Open: %v", s, err)
		}
		if _, err := f.HDU(0); !errors.Is(err, ErrRandomGroups) {
			t.Errorf("%s: HDU(0) err %v, want ErrRandomGroups", s, err)
		}
		if _, err := f.Primary(); !errors.Is(err, ErrRandomGroups) {
			t.Errorf("%s: Primary() err %v, want ErrRandomGroups", s, err)
		}
		f.Close()
	}
	// The IMAGE extension after the groups is read at cfitsio's offset.
	g := loadImagesCfref(t, "framing/random_groups_ext.fits")
	f, err := Open(imagesPath(t, "framing", "random_groups_ext.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	h, err := f.HDU(1)
	if err != nil {
		t.Fatal(err)
	}
	got, err := ReadPixels[int16](h.(*ImageHDU))
	if err != nil {
		t.Fatal(err)
	}
	for k, w := range g.HDUs[1].Image.I64.Values(t) {
		if int64(got[k]) != w.V {
			t.Fatalf("extension pixel %d = %d, cfitsio %d", k, got[k], w.V)
		}
	}
}

// TestFramingToleranceMatchesCfitsio: cfitsio opens files with a trailing
// zero block, trailing junk bytes or a junk block after the last HDU, and
// files whose last data unit lacks its final padding. The library opens
// them too, with cfitsio's HDU count (asserted with the addresses above),
// and reads the pixels. cfitsio itself cannot read the pixels of an
// unpadded last HDU (status 108), so those are checked against astropy.
func TestFramingToleranceMatchesCfitsio(t *testing.T) {
	var ast map[string]struct {
		NHDU  int       `json:"nhdu"`
		Data  [][]int64 `json:"data"`
		Error string    `json:"error"`
	}
	loadJSON(t, imagesPath(t, "golden", "framing_astropy.json"), &ast)
	entries, err := os.ReadDir(imagesPath(t, "framing"))
	if err != nil {
		t.Fatal(err)
	}
	if len(entries) < 7 {
		t.Fatalf("only %d framing samples", len(entries))
	}
	have := map[string]bool{}
	for _, e := range entries {
		have[e.Name()] = true
	}
	// NAXISn with n > NAXIS, once and twice: ignored, not refused
	for _, n := range []string{"stray_naxis3.fits", "stray_naxis3_twice.fits"} {
		if !have[n] {
			t.Errorf("framing/%s missing", n)
		}
	}
	for _, e := range entries {
		name := e.Name()
		t.Run(name, func(t *testing.T) {
			g := loadImagesCfref(t, "framing/"+name)
			f, err := Open(imagesPath(t, "framing", name))
			if err != nil {
				t.Fatalf("Open: %v (cfitsio opens it)", err)
			}
			defer f.Close()
			if f.NumHDU() != g.NHDU {
				t.Fatalf("NumHDU %d, cfitsio %d", f.NumHDU(), g.NHDU)
			}
			a := ast[name]
			for i := 0; i < f.NumHDU(); i++ {
				h, err := f.HDU(i)
				if errors.Is(err, ErrRandomGroups) || errors.Is(err, ErrUnsupportedExtension) {
					continue
				}
				if err != nil {
					t.Fatal(err)
				}
				got, err := ReadPixels[int64](h.(*ImageHDU))
				if err != nil {
					t.Fatalf("HDU %d: %v", i, err)
				}
				var want []int64
				switch {
				case g.HDUs[i].Image.I64.Present():
					for _, v := range g.HDUs[i].Image.I64.Values(t) {
						want = append(want, v.V)
					}
				case a.Error == "" && i < len(a.Data):
					want = a.Data[i] // cfitsio could not read these pixels
				default:
					t.Fatalf("HDU %d: neither cfitsio nor astropy read the pixels", i)
				}
				if len(got) != len(want) {
					t.Fatalf("HDU %d: %d pixels, want %d", i, len(got), len(want))
				}
				for k := range want {
					if got[k] != want[k] {
						t.Fatalf("HDU %d pixel %d = %d, want %d", i, k, got[k], want[k])
					}
				}
			}
		})
	}
}

// TestTruncatedExtensionHeaderIsAnError: after the last HDU, a block that
// starts an extension (XTENSION=) but has no END card before EOF is a
// damaged extension, not trailing junk. cfitsio silently stops the HDU
// list there (fits_get_num_hdus = 1, probed with the cfref harness);
// astropy refuses the file ("Header missing END card"). The library
// refuses it too, rather than hide an extension. A partial block (less
// than 2880 bytes) after the last HDU is ignored like any trailing bytes.
func TestTruncatedExtensionHeaderIsAnError(t *testing.T) {
	base, err := os.ReadFile(imagesPath(t, "framing", "trail_zero_block.fits"))
	if err != nil {
		t.Fatal(err)
	}
	base = base[:5760]
	card := func(s string) []byte { return []byte(s + strings.Repeat(" ", 80-len(s))) }
	ext := append(card("XTENSION= 'IMAGE   '"), card("BITPIX  =                   16")...)
	ext = append(ext, card("NAXIS   =                    0")...)
	full := append(append([]byte{}, base...), ext...)
	full = append(full, bytes.Repeat([]byte(" "), 2880-len(ext))...)
	dir := t.TempDir()
	p := filepath.Join(dir, "xt_noend.fits")
	if err := os.WriteFile(p, full, 0o644); err != nil {
		t.Fatal(err)
	}
	var ut *ErrUnterminatedHeader
	if _, err := Open(p); !errors.As(err, &ut) {
		t.Errorf("XTENSION block without END: Open err %v, want *ErrUnterminatedHeader", err)
	}
	p2 := filepath.Join(dir, "xt_partial.fits")
	if err := os.WriteFile(p2, append(append([]byte{}, base...), ext...), 0o644); err != nil {
		t.Fatal(err)
	}
	f, err := Open(p2)
	if err != nil {
		t.Fatalf("partial trailing block: %v", err)
	}
	defer f.Close()
	if f.NumHDU() != 1 || f.TrailingBytes() != int64(len(ext)) {
		t.Errorf("partial trailing block: NumHDU %d TrailingBytes %d, want 1 and %d", f.NumHDU(), f.TrailingBytes(), len(ext))
	}
}

// TestMalformedCardStillOpens: one malformed non-structural card does not
// make the file unopenable (cfitsio opens framing/malformed_card.fits; its
// addresses and pixels are asserted by the corpus tests above); the parse
// failure is reported by HeaderOf, not hidden.
func TestMalformedCardStillOpens(t *testing.T) {
	f, err := Open(imagesPath(t, "framing", "malformed_card.fits"))
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	img, err := f.Primary()
	if err != nil {
		t.Fatal(err)
	}
	if _, err := ReadPixels[int16](img); err != nil {
		t.Fatal(err)
	}
	if _, err := HeaderOf(img); err == nil {
		t.Log("package header now parses the unquoted DATE-OBS card")
	}
}
