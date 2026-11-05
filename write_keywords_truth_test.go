package fits

import (
	"bytes"
	"encoding/hex"
	"errors"
	"path/filepath"
	"testing"
)

// keywordsOverwriteGolden is testdata/images/golden/keywords_overwrite.json:
// for each keywords/ sample, cfitsio's fits_write_img of the values over
// that sample's image, scaled with the BSCALE/BZERO cfitsio read from its
// header (imgref overwrite): its status and the data unit it wrote.
type keywordsOverwriteGolden map[string]struct {
	Datatype int      `json:"datatype"`
	Values   []string `json:"values"`
	Status   int      `json:"status"`
	Data     string   `json:"data"`
}

// TestWriteImageKeywordsMatchCfitsio: WriteImage given the header of each
// keywords/ sample that the scaling-card rule accepts as the caller header
// (D-exponent, '+', '.5', '2.', '1E1' values, BLANK '+2', the smallest
// normal BSCALE) stores the values exactly as cfitsio's fits_write_img
// stores them under that header (golden/keywords_overwrite.json), byte for
// byte. Headers the rule refuses (HIERARCH and lowercase names, hex
// values, '=' in odd columns, ...), duplicates and BSCALE 0 are refused
// (ErrInvalidScaling, TestWriteImageRefusesInvalidScaling), and where cfitsio reports
// NUM_OVERFLOW (412: a value outside the stored range, e.g. under a
// subnormal BSCALE) WriteImage returns ErrNumOverflow and writes nothing.
// (A writer that stores with the wrong scaling cannot pass by being read
// back by a reader with the same mistake: the bytes are compared with
// cfitsio's.)
func TestWriteImageKeywordsMatchCfitsio(t *testing.T) {
	var g keywordsOverwriteGolden
	loadJSON(t, imagesPath(t, "golden", "keywords_overwrite.json"), &g)
	samples := keywordSamples(t)
	refuse := writeRefusedKeywordSamples(t)
	dir := t.TempDir()
	written, overflow := 0, 0
	for _, name := range samples {
		c, ok := g[name]
		if !ok {
			t.Fatalf("%s: no keywords_overwrite golden", name)
		}
		if refuse[name] {
			continue
		}
		t.Run(name, func(t *testing.T) {
			hdr, err := HeaderOf(openImage(t, imagesPath(t, "keywords", name)))
			if err != nil {
				t.Fatalf("header does not parse: %v", err)
			}
			p := filepath.Join(dir, name)
			f, err := Create(p)
			if err != nil {
				t.Fatal(err)
			}
			var werr error
			switch c.Datatype {
			case 21:
				_, werr = WriteImage(f, hdr.Clone(), []int64{int64(len(c.Values))}, parseIntVals[int16](t, c.Values, 16))
			case 31:
				_, werr = WriteImage(f, hdr.Clone(), []int64{int64(len(c.Values))}, parseIntVals[int32](t, c.Values, 32))
			default:
				t.Fatalf("datatype %d", c.Datatype)
			}
			n := f.NumHDU()
			if err := f.Close(); err != nil && werr == nil {
				t.Fatal(err)
			}
			switch c.Status {
			case 412:
				overflow++
				if !errors.Is(werr, ErrNumOverflow) || n != 0 {
					t.Errorf("cfitsio reports NUM_OVERFLOW (412); WriteImage: err %v, %d HDUs; want ErrNumOverflow and none", werr, n)
				}
				return
			case 0:
			default:
				t.Fatalf("golden status %d", c.Status)
			}
			if werr != nil || n != 1 {
				t.Fatalf("WriteImage: err %v, %d HDUs (cfitsio writes, status 0)", werr, n)
			}
			written++
			_, data := dataUnit(t, p)
			want, err := hex.DecodeString(c.Data)
			if err != nil {
				t.Fatal(err)
			}
			if !bytes.Equal(data[:len(want)], want) {
				t.Errorf("data % x, cfitsio fits_write_img % x", data[:len(want)], want)
			}
		})
	}
	// bscale_min_normal_rounded (2^-1022): every nonzero value overflows
	if overflow != 1 || written < 6 {
		t.Errorf("%d overflow refusals (want 1), %d headers written (want >= 6)", overflow, written)
	}
}
