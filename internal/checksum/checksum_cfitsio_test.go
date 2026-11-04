package checksum

import (
	"bytes"
	"encoding/json"
	"os"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestEncodeMatchesCfitsio: Encode equals cfitsio's fits_encode_chksum
// (ffesum) for edge values, with and without complement, and Decode
// inverts it (testdata/images/golden/esum.json, made by imgref esum).
func TestEncodeMatchesCfitsio(t *testing.T) {
	p := testutil.TestdataPath(t, "images", "golden", "esum.json")
	testutil.RequireFile(t, p)
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	var vec []struct {
		Value      uint32 `json:"value"`
		Complement int    `json:"complement"`
		ASCII      string `json:"ascii"`
	}
	dec := json.NewDecoder(bytes.NewReader(b))
	dec.DisallowUnknownFields()
	if err := dec.Decode(&vec); err != nil {
		t.Fatal(err)
	}
	if len(vec) < 20 {
		t.Fatalf("only %d vectors", len(vec))
	}
	for _, v := range vec {
		got := Encode(v.Value, v.Complement == 1)
		if string(got[:]) != v.ASCII {
			t.Errorf("Encode(%#x, %v) = %q, cfitsio %q", v.Value, v.Complement == 1, got, v.ASCII)
		}
		if back := Decode(got, v.Complement == 1); back != v.Value && !(v.Value == 0 || v.Value == 0xFFFFFFFF) {
			// 0 and 0xFFFFFFFF are the two ones'-complement zeros.
			t.Errorf("Decode(Encode(%#x)) = %#x", v.Value, back)
		}
	}
}

// TestComputeLargeBufferMatchesCfitsio: summing a multi-block data unit in
// one call equals cfitsio's fits_get_chksum DATASUM. A 32-bit accumulator
// overflows after ~45 blocks; the real HST/STIS file and the 100x100
// uint16 image are several blocks each, and a synthetic 200-block buffer of
// 0xFF bytes is checked against its closed form.
func TestComputeLargeBufferMatchesCfitsio(t *testing.T) {
	for _, s := range []string{"astropy_data/fixed-1890.fits", "astropy_data/o4sp040b0_raw.fits", "checksum/stis_raw.fits"} {
		gp := testutil.TestdataPath(t, "images", "golden", goldenName(s)+".json")
		gb := testutil.MustReadGolden(t, gp)
		cf, err := testutil.DecodeCfref(gb)
		if err != nil {
			t.Fatal(err)
		}
		var ck struct {
			File    string `json:"file"`
			Cfitsio string `json:"cfitsio"`
			HDUs    []struct {
				Index         int    `json:"index"`
				Datasum       uint32 `json:"datasum"`
				HDUSum        uint32 `json:"hdusum"`
				VerifyDatasum int    `json:"verify_datasum"`
				VerifyHDUSum  int    `json:"verify_hdusum"`
				Error         string `json:"error"`
			} `json:"hdus"`
		}
		cb := testutil.MustReadGolden(t, testutil.TestdataPath(t, "images", "golden", goldenName(s)+".chksum.json"))
		if err := json.Unmarshal(cb, &ck); err != nil {
			t.Fatal(err)
		}
		raw, err := os.ReadFile(testutil.TestdataPath(t, "images", s))
		if err != nil {
			t.Fatal(err)
		}
		for i, h := range cf.HDUs {
			data := raw[h.DataStart:h.DataEnd]
			if got := Compute(data); got != ck.HDUs[i].Datasum {
				t.Errorf("%s HDU %d (%d blocks): Compute = %d, cfitsio DATASUM %d", s, i, len(data)/2880, got, ck.HDUs[i].Datasum)
			}
			whole := raw[h.HeaderStart:h.DataEnd]
			if got := Compute(whole); got != ck.HDUs[i].HDUSum {
				t.Errorf("%s HDU %d: Compute(header+data) = %d, cfitsio %d", s, i, got, ck.HDUs[i].HDUSum)
			}
		}
	}
	// 200 blocks of 0xFF: every 16-bit word is 0xFFFF (ones'-complement
	// negative zero), so the 32-bit ones'-complement sum is 0xFFFFFFFF.
	ff := bytes.Repeat([]byte{0xFF}, 200*2880)
	if got := Compute(ff); got != 0xFFFFFFFF {
		t.Errorf("Compute(200 blocks of 0xFF) = %#x, want 0xffffffff", got)
	}
	// 100 blocks of 0x00 0x01 words: hi half sums 720*100 ones, lo half the
	// same; 72000 < 0xFFFF so no carry: sum = 72000<<16 | 72000 folded.
	w := bytes.Repeat([]byte{0, 1}, 100*1440)
	got := Compute(w)
	// 72000 = 0x11940: fold once -> hi 0x1940 + carry(lo)=1, lo 0x1940 + carry(hi)=1.
	if want := uint32(0x1941<<16 | 0x1941); got != want {
		t.Errorf("Compute(100 blocks of 0x0001) = %#x, want %#x", got, want)
	}
}

func goldenName(s string) string {
	return string(bytes.ReplaceAll([]byte(s), []byte("/"), []byte("__")))
}
