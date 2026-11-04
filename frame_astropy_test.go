package fits

import (
	"bytes"
	"encoding/json"
	"math"
	"os"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// astropyF32 reads path's primary HDU with astropy and returns its shape
// (numpy order) and every pixel's float32 bit pattern in FITS order.
func astropyF32(t *testing.T, path string) (shape []int, bits []uint32) {
	t.Helper()
	script := `
import json, sys
import numpy as np
from astropy.io import fits
with fits.open(sys.argv[1]) as h:
    d = h[0].data
    assert d.dtype.kind == "f" and d.dtype.itemsize == 4, d.dtype
    print(json.dumps({"shape": list(d.shape), "bits": np.ascontiguousarray(d).view(">u4").ravel().astype(int).tolist()}))
`
	out := testutil.RunOK(t, testutil.ToolAstropy, testutil.Call{Args: []string{"-c", script, path}, Inputs: []string{path}}).Stdout
	var r struct {
		Shape []int    `json:"shape"`
		Bits  []uint32 `json:"bits"`
	}
	if err := json.Unmarshal(out, &r); err != nil {
		t.Fatalf("astropy output: %v\n%s", err, out)
	}
	return r.Shape, r.Bits
}

var awkwardFloats = []float32{
	0, float32(math.Copysign(0, -1)), 1, -1, 0.1, 1e-45, 1.1754942e-38,
	math.MaxFloat32, -math.MaxFloat32, float32(math.Inf(1)), float32(math.Inf(-1)),
	float32(math.NaN()), math.Float32frombits(0x7fc12345), 65535, 3.5,
	// Signaling NaNs: a float32 -> float64 -> float32 round trip quiets
	// them (sets bit 22), so they pin the verbatim float32 write path.
	math.Float32frombits(0x7f800001), math.Float32frombits(0x7f812345), math.Float32frombits(0xff800001),
}

// TestWriteMono_AstropyReadsVerbatim: WriteMono's "values are written
// verbatim (no scaling, no clamping)" and "NAXIS1=width, NAXIS2=height,
// pix[y*width+x]" claims, checked by astropy reading the file back bit for
// bit (NaN payloads, -0, ±Inf, denormals included) and by fitsverify.
func TestWriteMono_AstropyReadsVerbatim(t *testing.T) {
	const w, h = 6, 3
	pix := make([]float32, w*h)
	copy(pix, awkwardFloats)
	path := filepath.Join(t.TempDir(), "mono.fits")
	if err := WriteMono(path, pix, w, h); err != nil {
		t.Fatal(err)
	}
	if r := testutil.Fitsverify(t, path); !r.OK() {
		t.Errorf("fitsverify: %d warnings %d errors\n%s", r.Warnings, r.Errors, r.Report)
	}
	shape, bits := astropyF32(t, path)
	if len(shape) != 2 || shape[0] != h || shape[1] != w {
		t.Fatalf("astropy shape %v, want [%d %d] (NAXIS2, NAXIS1)", shape, h, w)
	}
	for i := range pix {
		if bits[i] != math.Float32bits(pix[i]) {
			t.Errorf("pixel %d: astropy 0x%08x, written 0x%08x", i, bits[i], math.Float32bits(pix[i]))
		}
	}
}

// TestWriteRGB_AstropyPlanes: WriteRGB writes NAXIS3=3 planes in R, G, B
// order, verbatim; astropy's data[k] must be plane k.
func TestWriteRGB_AstropyPlanes(t *testing.T) {
	const w, h = 3, 5
	n := w * h
	r, g, b := make([]float32, n), make([]float32, n), make([]float32, n)
	for i := 0; i < n; i++ {
		r[i] = awkwardFloats[i%len(awkwardFloats)]
		g[i] = float32(i) + 0.25
		b[i] = -float32(i) * 1e-3
	}
	path := filepath.Join(t.TempDir(), "rgb.fits")
	if err := WriteRGB(path, r, g, b, w, h); err != nil {
		t.Fatal(err)
	}
	if res := testutil.Fitsverify(t, path); !res.OK() {
		t.Errorf("fitsverify: %d warnings %d errors\n%s", res.Warnings, res.Errors, res.Report)
	}
	shape, bits := astropyF32(t, path)
	if len(shape) != 3 || shape[0] != 3 || shape[1] != h || shape[2] != w {
		t.Fatalf("astropy shape %v, want [3 %d %d]", shape, h, w)
	}
	for k, plane := range [][]float32{r, g, b} {
		for i := 0; i < n; i++ {
			if bits[k*n+i] != math.Float32bits(plane[i]) {
				t.Errorf("plane %d pixel %d: astropy 0x%08x, written 0x%08x", k, i, bits[k*n+i], math.Float32bits(plane[i]))
			}
		}
	}
}

// TestWriteMono_TruncatesExisting: "Any existing file at path is
// truncated" — a longer file at path must not leave trailing bytes.
func TestWriteMono_TruncatesExisting(t *testing.T) {
	path := filepath.Join(t.TempDir(), "old.fits")
	if err := os.WriteFile(path, bytes.Repeat([]byte("X"), 5*2880+17), 0o644); err != nil {
		t.Fatal(err)
	}
	pix := []float32{1, 2, 3, 4}
	if err := WriteMono(path, pix, 2, 2); err != nil {
		t.Fatal(err)
	}
	st, err := os.Stat(path)
	if err != nil {
		t.Fatal(err)
	}
	if st.Size() != 2*2880 {
		t.Fatalf("file is %d bytes after WriteMono, want 5760 (one header + one data block)", st.Size())
	}
	_, bits := astropyF32(t, path)
	for i := range pix {
		if bits[i] != math.Float32bits(pix[i]) {
			t.Errorf("pixel %d: 0x%08x", i, bits[i])
		}
	}
}
