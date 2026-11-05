package fits

import (
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/checksum"
	"github.com/dmarkham/fits/internal/testutil"
)

// hduSums computes the data-unit and whole-HDU 32-bit ones'-complement sums
// of HDU i from the bytes on disk (zero padding beyond EOF, as FITS
// defines the pad), using the library's checksum package.
func hduSums(t *testing.T, path string, f *File, i int) (data, hdu uint32) {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatal(err)
	}
	r := f.hdus[i]
	pad := func(lo, hi int64) []byte {
		out := make([]byte, hi-lo)
		if lo < int64(len(b)) {
			copy(out, b[lo:min(hi, int64(len(b)))])
		}
		return out
	}
	data = checksum.Compute(pad(r.dataStart, r.paddedEnd))
	hdu = checksum.Update(data, pad(r.headerStart, r.dataStart))
	return data, hdu
}

// TestChecksumMatchesCfitsio: for every HDU of every sample (images,
// tables, compressed, files with valid, wrong and absent CHECKSUM/DATASUM,
// cfitsio- and astropy-written), the data sum and HDU sum equal cfitsio's
// fits_get_chksum, and ChecksumStatus equals fits_verify_chksum
// (1 correct, 0 absent, -1 wrong). VerifyChecksum fails exactly when
// cfitsio reports -1. Where fits_verify_chksum itself fails (a DATASUM or
// CHECKSUM card with no value, VALUE_UNDEFINED 204) ChecksumStatus returns
// ErrUndefinedChecksum.
func TestChecksumMatchesCfitsio(t *testing.T) {
	sawStatus := map[[2]int]bool{}
	sawUndefined := false
	samples := allImageSamples(t)
	// The DATASUM/CHECKSUM card-order, name-case and range samples must
	// stay in the corpus (gen.py gen_datasum).
	have := map[string]bool{}
	for _, s := range samples {
		have[s.name] = true
	}
	for _, n := range []string{"datasum_dup_after_checksum", "datasum_dup_before_checksum", "datasum_dup_no_checksum",
		"datasum_dup_wrap", "datasum_dup_undefined_after", "checksum_dup_empty_first", "datasum_lowercase_ok",
		"datasum_lowercase_bad", "datasum_2_32_empty", "datasum_plus_2_32", "datasum_neg_fraction_empty"} {
		if !have["checksum/"+n+".fits"] {
			t.Errorf("checksum/%s.fits is not in the corpus", n)
		}
	}
	for _, s := range samples {
		t.Run(s.name, func(t *testing.T) {
			f, err := Open(s.path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			name := s.name
			if !strings.Contains(name, "/") {
				name = "testdata/" + name
			}
			gck := loadChksumGolden(t, name)
			ck := &gck
			for i, g := range s.golden.HDUs {
				want := [2]int{g.Checksum.Datasum, g.Checksum.HDUSum}
				if ck != nil && ck.HDUs[i].Error == "" {
					d, h := hduSums(t, s.path, f, i)
					if d != ck.HDUs[i].Datasum || h != ck.HDUs[i].HDUSum {
						t.Errorf("HDU %d: data/HDU sum %d/%d, cfitsio %d/%d", i, d, h, ck.HDUs[i].Datasum, ck.HDUs[i].HDUSum)
					}
					want = [2]int{ck.HDUs[i].VerifyDatasum, ck.HDUs[i].VerifyHDUSum}
				}
				h, err := hduForChecksum(f, i)
				if errors.Is(err, ErrUnsupportedExtension) {
					continue // framing/xtension_foreign.fits: refused, see TestUnsupportedExtensionFramedAndRefused
				}
				if err != nil {
					t.Fatalf("HDU %d: %v", i, err)
				}
				ds, hs, err := ChecksumStatus(h)
				if ck != nil && ck.HDUs[i].VerifyError != "" {
					// cfitsio's fits_verify_chksum fails (a DATASUM or
					// CHECKSUM card with no value): so must we.
					if err == nil || !errors.Is(VerifyChecksum(h), ErrUndefinedChecksum) {
						t.Errorf("HDU %d: ChecksumStatus = %d/%d, %v; cfitsio fails: %s", i, ds, hs, err, ck.HDUs[i].VerifyError)
					}
					sawUndefined = true
					continue
				}
				if err != nil {
					t.Errorf("HDU %d: ChecksumStatus: %v", i, err)
					continue
				}
				if ds != want[0] || hs != want[1] {
					t.Errorf("HDU %d: ChecksumStatus = %d/%d, cfitsio fits_verify_chksum %d/%d", i, ds, hs, want[0], want[1])
				}
				sawStatus[want] = true
				verr := VerifyChecksum(h)
				if bad := want[0] < 0 || want[1] < 0; bad != (verr != nil) {
					t.Errorf("HDU %d: VerifyChecksum = %v, cfitsio status %v", i, verr, want)
				}
			}
		})
	}
	if !sawUndefined {
		t.Error("corpus has no HDU whose DATASUM/CHECKSUM has no value (cfitsio VALUE_UNDEFINED)")
	}
	for _, st := range [][2]int{{1, 1}, {0, 0}, {-1, -1}, {1, -1}, {0, 1}} {
		if !sawStatus[st] {
			t.Errorf("corpus has no HDU with cfitsio checksum status %v", st)
		}
	}
}

// hduForChecksum returns HDU i; random-groups HDUs (which HDU refuses) are
// wrapped as an image so their checksum can still be verified.
func hduForChecksum(f *File, i int) (HDU, error) {
	h, err := f.HDU(i)
	if errors.Is(err, ErrRandomGroups) {
		return &ImageHDU{rec: f.hdus[i]}, nil
	}
	return h, err
}

// TestWriteChecksumCfitsioVerifies: CHECKSUM/DATASUM written by the
// library verify in cfitsio (fits_verify_chksum 1/1, with cfitsio's own
// DATASUM value) and in fitsverify (no warnings, no errors), for Go-written
// images of every matrix case, a multi-HDU real file and an image+table
// file, after a header that must grow into a new block.
func TestWriteChecksumCfitsioVerifies(t *testing.T) {
	dir := t.TempDir()
	var files []string
	for _, c := range loadCases(t) {
		files = append(files, writeCaseFile(t, dir, c))
	}
	for _, src := range []string{"astropy_data/o4sp040b0_raw.fits", "astropy_data/checksum.fits", "astropy_data/arange.fits"} {
		dst := filepath.Join(dir, filepath.Base(src))
		copyFile(t, imagesPath(t, src), dst)
		files = append(files, dst)
	}
	for _, p := range files {
		e, err := OpenForEdit(p)
		if err != nil {
			t.Fatalf("%s: %v", p, err)
		}
		for i := 0; i < e.NumHDU(); i++ {
			h, err := e.HDU(i)
			if err != nil {
				t.Fatal(err)
			}
			if _, _, err := WriteChecksum(h); err != nil {
				t.Fatalf("%s HDU %d: WriteChecksum: %v", p, i, err)
			}
		}
		if err := e.Flush(); err != nil {
			t.Fatalf("%s: Flush: %v", p, err)
		}
		if err := e.Close(); err != nil {
			t.Fatal(err)
		}
		ck := imgrefChksum(t, p)
		f, err := Open(p)
		if err != nil {
			t.Fatal(err)
		}
		for i, g := range ck.HDUs {
			if g.VerifyDatasum != 1 || g.VerifyHDUSum != 1 {
				t.Errorf("%s HDU %d: cfitsio fits_verify_chksum %d/%d, want 1/1", filepath.Base(p), i, g.VerifyDatasum, g.VerifyHDUSum)
			}
			h, _ := f.HDU(i)
			if ds, hs, err := ChecksumStatus(h); err != nil || ds != 1 || hs != 1 {
				t.Errorf("%s HDU %d: ChecksumStatus %d/%d %v, want 1/1", filepath.Base(p), i, ds, hs, err)
			}
		}
		f.Close()
		if r := testutil.Fitsverify(t, p); !r.OK() {
			t.Errorf("%s: fitsverify %d warnings, %d errors:\n%s", filepath.Base(p), r.Warnings, r.Errors, r.Report)
		}
	}
}

// TestWriteChecksumRefusesCreateMode: in ModeCreate, Flush does not rewrite
// headers, so a checksum computed there would never reach the disk.
// WriteChecksum must say so instead of reporting success.
func TestWriteChecksumRefusesCreateMode(t *testing.T) {
	p := filepath.Join(t.TempDir(), "c.fits")
	f, err := Create(p)
	if err != nil {
		t.Fatal(err)
	}
	h, err := WriteImage(f, nil, []int64{3}, []int16{1, 2, 3})
	if err != nil {
		t.Fatal(err)
	}
	if _, _, err := WriteChecksum(h); err == nil {
		t.Error("WriteChecksum in ModeCreate succeeded although nothing would be persisted")
	}
	f.Close()
}

func copyFile(t *testing.T, src, dst string) {
	t.Helper()
	in, err := os.Open(src)
	if err != nil {
		t.Fatal(err)
	}
	defer in.Close()
	out, err := os.Create(dst)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := io.Copy(out, in); err != nil {
		t.Fatal(err)
	}
	if err := out.Close(); err != nil {
		t.Fatal(err)
	}
}
