package fits

import (
	"bytes"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/internal/testutil"
)

// headerPreserveSamples is every committed FITS sample whose headers the
// no-mutation tests below must reproduce byte for byte: all of testdata/
// plus the header package's CONTINUE / malformed-card / adversarial / cfitsio-written
// samples in header/testdata.
func headerPreserveSamples(t *testing.T) []string {
	t.Helper()
	var out []string
	for _, pat := range []string{"*.fits", "*.fit", "*.std"} {
		m, _ := filepath.Glob(testutil.TestdataPath(t, pat))
		out = append(out, m...)
	}
	for _, n := range []string{"continue.fits", "malformed.fits", "adversarial.fits", "cfwritten.fits"} {
		p := filepath.Join(testutil.RepoRoot(t), "header", "testdata", n)
		testutil.RequireFile(t, p)
		out = append(out, p)
	}
	sort.Strings(out)
	if len(out) < 10 {
		t.Fatalf("only %d samples found", len(out))
	}
	return out
}

func copyToTemp(t *testing.T, src string) string {
	t.Helper()
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	dst := filepath.Join(t.TempDir(), filepath.Base(src))
	if err := os.WriteFile(dst, b, 0o644); err != nil {
		t.Fatal(err)
	}
	return dst
}

// TestFlushWithoutMutationIsByteIdentical: OpenForEdit, load every HDU's
// header (which parses it, joining CONTINUE chains and decoding malformed
// cards), Flush and Close must leave the file byte-for-byte unchanged.
// Before the fix a CONTINUE'd string lost its continuation cards here.
func TestFlushWithoutMutationIsByteIdentical(t *testing.T) {
	for _, src := range headerPreserveSamples(t) {
		t.Run(filepath.Base(src), func(t *testing.T) {
			orig, _ := os.ReadFile(src)
			p := copyToTemp(t, src)
			f, err := OpenForEdit(p)
			if err != nil {
				t.Fatalf("OpenForEdit: %v", err)
			}
			for i := 0; i < f.NumHDU(); i++ {
				h, err := f.HDU(i)
				if err != nil {
					t.Fatalf("HDU %d: %v", i, err)
				}
				_ = h.Header()
			}
			if err := f.Flush(); err != nil {
				t.Fatalf("Flush: %v", err)
			}
			if err := f.Close(); err != nil {
				t.Fatalf("Close: %v", err)
			}
			got, _ := os.ReadFile(p)
			if !bytes.Equal(got, orig) {
				t.Fatalf("file changed by a no-op edit: %s", diffBlock(got, orig))
			}
		})
	}
}

// TestCopyHDUWithHeaderUnmodifiedIsByteIdentical: EditFile copying every
// HDU through CopyHDUWithHeader with its own parsed (unmodified) header
// reproduces the file exactly, tile-compressed HDUs included (hduRecordOf
// in checksum.go resolves a *CompressedImageHDU to its backing record).
func TestCopyHDUWithHeaderUnmodifiedIsByteIdentical(t *testing.T) {
	for _, src := range headerPreserveSamples(t) {
		t.Run(filepath.Base(src), func(t *testing.T) {
			orig, _ := os.ReadFile(src)
			in, err := Open(src)
			if err != nil {
				t.Fatal(err)
			}
			for i := 0; i < in.NumHDU(); i++ {
				h, err := in.HDU(i)
				if err != nil {
					t.Fatal(err)
				}
				rec := in.hdus[i]
				got, err := header.Encode(h.Header())
				if err != nil {
					t.Fatalf("HDU %d: Encode: %v", i, err)
				}
				if want := orig[rec.headerStart:rec.dataStart]; !bytes.Equal(got, want) {
					t.Errorf("HDU %d: encoded header differs from the file: %s", i, diffBlock(got, want))
				}
			}
			in.Close()

			p := copyToTemp(t, src)
			err = EditFile(p, func(in *File, out *Writer) error {
				for i := 0; i < in.NumHDU(); i++ {
					h, err := in.HDU(i)
					if err != nil {
						return err
					}
					if err := out.CopyHDUWithHeader(h, h.Header()); err != nil {
						return err
					}
				}
				return nil
			})
			if err != nil {
				t.Fatalf("EditFile: %v", err)
			}
			got, _ := os.ReadFile(p)
			if !bytes.Equal(got, orig) {
				t.Fatalf("CopyHDUWithHeader changed the file: %s", diffBlock(got, orig))
			}
		})
	}
}

// TestOpenToleratesWhatCfitsioOpens: the header samples with lone
// CONTINUE cards and malformed values open (cfitsio opens them: their
// goldens list every card and no error), Header() does not panic, and the
// card after the malformed ones is still decoded.
func TestOpenToleratesWhatCfitsioOpens(t *testing.T) {
	for _, n := range []string{"continue.fits", "malformed.fits", "adversarial.fits"} {
		t.Run(n, func(t *testing.T) {
			g := filepath.Join(testutil.RepoRoot(t), "header", "testdata", "golden", n+".json")
			testutil.RequireFile(t, g)
			gb, _ := os.ReadFile(g)
			cf, err := testutil.DecodeCfref(gb)
			if err != nil {
				t.Fatal(err)
			}
			if len(cf.HDUs) != 1 || len(cf.HDUs[0].Errors) != 0 {
				t.Fatalf("golden: cfitsio did not open %s cleanly: %+v", n, cf.HDUs)
			}
			f, err := Open(filepath.Join(testutil.RepoRoot(t), "header", "testdata", n))
			if err != nil {
				t.Fatalf("Open: %v (cfitsio opens this file)", err)
			}
			defer f.Close()
			if f.NumHDU() != cf.NHDU {
				t.Fatalf("NumHDU %d, cfitsio %d", f.NumHDU(), cf.NHDU)
			}
			h, err := f.HDU(0)
			if err != nil {
				t.Fatal(err)
			}
			var hdr interface{ String(string) (string, error) }
			func() {
				defer func() {
					if r := recover(); r != nil {
						t.Fatalf("Header() panicked: %v", r)
					}
				}()
				hdr = h.Header()
			}()
			want := cf.HDUs[0].Card("FINAL")
			got, err := hdr.String("FINAL")
			if want == nil || want.Str == nil || err != nil || got != string(*want.Str) {
				t.Fatalf("FINAL = %q, %v; cfitsio %v", got, err, want)
			}
		})
	}
}

func diffBlock(got, want []byte) string {
	if len(got) != len(want) {
		return "length " + strconv.Itoa(len(got)) + " bytes, original " + strconv.Itoa(len(want))
	}
	for i := range got {
		if got[i] != want[i] {
			c := i / 80 * 80
			return "first difference at byte " + strconv.Itoa(i) + "\n got " + strings.TrimRight(string(got[c:c+80]), " ") + "\nwant " + strings.TrimRight(string(want[c:c+80]), " ")
		}
	}
	return ""
}
