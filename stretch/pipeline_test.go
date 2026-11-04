package stretch

import (
	"bytes"
	"errors"
	"fmt"
	"image/jpeg"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// TestDefaultPresets_Siril runs every default preset on the golden input
// and compares it with the siril-cli golden of the same name, whose
// recorded command must equal the preset's Command. Tolerance as for the
// individual operations: 0 ULP, except GHT (1-ULP glibc libm allowance,
// see requireGHT).
func TestDefaultPresets_Siril(t *testing.T) {
	presets := DefaultPresets()
	if len(presets) != 18 {
		t.Fatalf("%d presets, want 18", len(presets))
	}
	seen := map[string]bool{}
	for _, p := range presets {
		if seen[p.Name] {
			t.Fatalf("duplicate preset %q", p.Name)
		}
		seen[p.Name] = true
		t.Run(p.Name, func(t *testing.T) {
			c := loadCLI(t, p.Name)
			if p.Command != c.Command {
				t.Fatalf("preset Command %q, golden produced by %q", p.Command, c.Command)
			}
			if p.Label == "" {
				t.Error("empty Label")
			}
			requirePreset(t, p, "stretch_"+p.Name+".json")
		})
	}
}

// requirePreset applies p to the golden's input and compares the result
// with the golden: identical float32 bits, except at the libm-sensitive
// pixels of a GHT golden (see ghtLibm).
func requirePreset(t *testing.T, p Preset, rel string) {
	t.Helper()
	var c ghtCase
	readGolden(t, rel, &c)
	in, nx, ny := c.input(t)
	got := copyPlanes(in)
	p.ApplyFn(got[0], got[1], got[2], nx, ny)
	if !strings.HasPrefix(p.Command, "ght ") {
		c.check(t, p.Name, got)
		return
	}
	if c.Libm == nil {
		t.Fatalf("%s: GHT golden without a libm list", rel)
	}
	checkGHT(t, p.Name, c, in, got, ghtParamsOf(c.srefCase))
}

// TestDefaultPresets_SirilOracle runs every preset on the bright input
// (every channel median > 0.5: Siril's inverted autostretch branch) and on
// the mixed input (one bright channel of three: Siril's unlinked
// autostretch keeps the normal branch for all channels, which per-channel
// autostretch would not), against the verbatim Siril 1.4.0 oracle running
// the preset's command.
func TestDefaultPresets_SirilOracle(t *testing.T) {
	for _, tag := range []string{"bright", "mixed"} {
		for _, p := range DefaultPresets() {
			t.Run(p.Name+"_"+tag, func(t *testing.T) {
				c := loadSref(t, "preset_"+p.Name+"_"+tag)
				if c.Preset != p.Name || c.Command != p.Command {
					t.Fatalf("golden is for preset %q command %q; preset %q has %q", c.Preset, c.Command, p.Name, p.Command)
				}
				if c.Input.File != "inputs/"+tag+".f32" {
					t.Fatalf("golden input %s, want inputs/%s.f32", c.Input.File, tag)
				}
				requirePreset(t, p, "sref/preset_"+p.Name+"_"+tag+".json")
			})
		}
	}
}

func TestFindPreset(t *testing.T) {
	for _, p := range DefaultPresets() {
		q, err := FindPreset(p.Name)
		if err != nil || q.Name != p.Name || q.Command != p.Command {
			t.Errorf("FindPreset(%q) = %+v, %v", p.Name, q, err)
		}
	}
	if _, err := FindPreset("no_such_preset"); err == nil || !strings.Contains(err.Error(), "no_such_preset") {
		t.Errorf("FindPreset(unknown) error = %v", err)
	}
}

// TestToRGBImage: quantization floor(v*255+0.5) with clamping (the values
// below are worked out by hand), alpha 255, and orientation: data index
// y*nx+x lands at image point (x, y), row 0 on top.
func TestToRGBImage(t *testing.T) {
	nan, inf := float32(math.NaN()), float32(math.Inf(1))
	cases := []struct {
		v    float32
		want uint8
	}{
		{-1, 0}, {0, 0}, {nan, 0}, {float32(math.Inf(-1)), 0},
		{1, 255}, {2, 255}, {inf, 255},
		{0.5, 128},       // 127.5 + 0.5 = 128
		{1.0 / 255, 1},   // 1.0000001 + 0.5 -> 1
		{0.499 / 255, 0}, // 0.499 + 0.5 < 1
		{0.501 / 255, 1}, // 0.501 + 0.5 > 1
		{254.6 / 255, 255}, {254.4 / 255, 254},
		// float32 rounding boundaries (found by searching float32 values
		// next to (k+0.5)/255 with numpy float32 and exact rationals):
		// v*255 rounds up to exactly k+0.5 in float32 although the exact
		// product is just below it, so float32 gives k+1 where exact (or
		// float64) arithmetic gives k.
		{math.Float32frombits(989888640), 1},    // 0.0019607842: exact v*255 = 0.49999997, floor(+0.5) = 0
		{math.Float32frombits(1057030401), 129}, // 0.50392157: exact v*255 = 128.49999994, floor(+0.5) = 128
	}
	nx := len(cases)
	r := make([]float32, nx*2)
	g := make([]float32, nx*2)
	b := make([]float32, nx*2)
	for i, c := range cases {
		r[i] = c.v             // row 0
		g[nx+i] = c.v          // row 1
		b[i], b[nx+i] = 1, 0.5 // constant per row
	}
	img := ToRGBImage(r, g, b, nx, 2)
	if img.Bounds().Dx() != nx || img.Bounds().Dy() != 2 {
		t.Fatalf("bounds %v", img.Bounds())
	}
	for i, c := range cases {
		top, bot := img.RGBAAt(i, 0), img.RGBAAt(i, 1)
		if top.R != c.want || bot.G != c.want {
			t.Errorf("v=%v: R(row0)=%d G(row1)=%d, want %d", c.v, top.R, bot.G, c.want)
		}
		if top.G != 0 || bot.R != 0 || top.B != 255 || bot.B != 128 || top.A != 255 || bot.A != 255 {
			t.Errorf("pixel %d: top %+v bottom %+v (orientation or alpha wrong)", i, top, bot)
		}
	}
	defer func() {
		if recover() == nil {
			t.Error("short channel did not panic")
		}
	}()
	ToRGBImage(r[:3], g, b, nx, 2)
}

// TestEncodeJPEG: the output is a JPEG that an independent decoder
// (libjpeg-turbo, through opencv-python's cv2.imdecode) reads back at the
// right size and close to the source (lossy codec: the bound is PSNR >= 35
// dB at quality 95 on a smooth gradient; a broken encoder or swapped
// channels fall far below that). Go's own decoder must agree on the size.
func TestEncodeJPEG(t *testing.T) {
	nx, ny := 64, 48
	r := make([]float32, nx*ny)
	g := make([]float32, nx*ny)
	b := make([]float32, nx*ny)
	for y := 0; y < ny; y++ {
		for x := 0; x < nx; x++ {
			i := y*nx + x
			r[i] = float32(x) / float32(nx-1)
			g[i] = float32(y) / float32(ny-1)
			b[i] = 0.25
		}
	}
	img := ToRGBImage(r, g, b, nx, ny)
	var buf bytes.Buffer
	if err := EncodeJPEG(&buf, img, 95); err != nil {
		t.Fatal(err)
	}
	if d := buf.Bytes(); len(d) < 4 || d[0] != 0xFF || d[1] != 0xD8 || d[len(d)-2] != 0xFF || d[len(d)-1] != 0xD9 {
		t.Fatalf("not a JPEG stream (SOI/EOI markers missing)")
	}
	dec, err := jpeg.Decode(bytes.NewReader(buf.Bytes()))
	if err != nil || dec.Bounds() != img.Bounds() {
		t.Fatalf("image/jpeg decode: %v, bounds %v", err, dec.Bounds())
	}

	// writer errors are returned
	if err := EncodeJPEG(failWriter{}, img, 90); err == nil {
		t.Error("EncodeJPEG ignored a write error")
	}

	dir := t.TempDir()
	path := filepath.Join(dir, "x.jpg")
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatal(err)
	}
	src := make([]byte, 0, 3*nx*ny)
	for y := 0; y < ny; y++ {
		for x := 0; x < nx; x++ {
			c := img.RGBAAt(x, y)
			src = append(src, c.R, c.G, c.B)
		}
	}
	srcPath := filepath.Join(dir, "src.rgb")
	if err := os.WriteFile(srcPath, src, 0o644); err != nil {
		t.Fatal(err)
	}
	script := `import sys, cv2, numpy as np
im = cv2.imdecode(np.fromfile(sys.argv[1], np.uint8), cv2.IMREAD_COLOR)
src = np.fromfile(sys.argv[2], np.uint8).reshape(im.shape[0], im.shape[1], 3)
rgb = im[:, :, ::-1].astype(float)
mse = ((rgb - src) ** 2).mean()
print(im.shape[1], im.shape[0], 10 * np.log10(255.0 ** 2 / max(mse, 1e-12)))`
	out := testutil.RunOK(t, toolCv2, testutil.Call{Args: []string{"-c", script, path, srcPath}, Inputs: []string{path, srcPath}}).Stdout
	f := strings.Fields(string(out))
	if len(f) != 3 || f[0] != strconv.Itoa(nx) || f[1] != strconv.Itoa(ny) {
		t.Fatalf("cv2 decoded size %q, want %dx%d", out, nx, ny)
	}
	psnr, _ := strconv.ParseFloat(f[2], 64)
	t.Logf("cv2.imdecode: %sx%s, PSNR %.1f dB", f[0], f[1], psnr)
	if psnr < 35 {
		t.Errorf("PSNR %.1f dB < 35 dB", psnr)
	}
}

// toolCv2 is OpenCV's Python binding (opencv-python-headless in the
// FITS_ASTROPY_PYTHON venv) for tier 2; tier 1 replays its readings.
var toolCv2 = testutil.Tool{Name: "cv2", Python: true, Modules: []string{"cv2", "numpy"}}

type failWriter struct{}

func (failWriter) Write([]byte) (int, error) { return 0, errors.New("disk full") }

// jpegDQT returns the quantization tables of a JPEG stream in file order.
func jpegDQT(t *testing.T, b []byte) []dqtTable {
	t.Helper()
	var out []dqtTable
	for i := 2; i+4 <= len(b); {
		if b[i] != 0xFF {
			t.Fatalf("bad JPEG marker at %d", i)
		}
		m, n := b[i+1], int(b[i+2])<<8|int(b[i+3])
		if m == 0xDB {
			seg := b[i+4 : i+2+n]
			for k := 0; k < len(seg); {
				pq, tq := seg[k]>>4, int(seg[k]&15)
				tb := dqtTable{ID: tq, Bits: 8}
				for j := 0; j < 64; j++ {
					if pq != 0 {
						tb.Bits = 16
						tb.Q = append(tb.Q, int(seg[k+1+2*j])<<8|int(seg[k+2+2*j]))
					} else {
						tb.Q = append(tb.Q, int(seg[k+1+j]))
					}
				}
				out = append(out, tb)
				k += 1 + 64*int(1+pq)
			}
		}
		if m == 0xDA {
			break
		}
		i += 2 + n
	}
	return out
}

type dqtTable struct {
	ID   int   `json:"id"`
	Bits int   `json:"bits"`
	Q    []int `json:"q"`
}

// TestEncodeJPEG_QualityTables: EncodeJPEG uses its quality argument. At
// each quality the quantization tables it writes equal those libjpeg-turbo
// writes through cv2.imencode at that IMWRITE_JPEG_QUALITY (both scale the
// JPEG Annex K tables with the IJG rule), and out-of-range qualities are
// clamped to 1 and 100 as the doc says.
func TestEncodeJPEG_QualityTables(t *testing.T) {
	var g struct {
		Tables []struct {
			Quality int        `json:"quality"`
			DQT     []dqtTable `json:"dqt"`
		} `json:"tables"`
	}
	readGolden(t, "jpeg_dqt_cv2.json", &g)
	if len(g.Tables) < 5 {
		t.Fatalf("%d qualities", len(g.Tables))
	}
	nx, ny := 16, 16
	r, gg, b := make([]float32, nx*ny), make([]float32, nx*ny), make([]float32, nx*ny)
	for i := range r {
		r[i], gg[i], b[i] = float32(i%nx)/16, float32(i/nx)/16, 0.25
	}
	img := ToRGBImage(r, gg, b, nx, ny)
	enc := func(q int) []dqtTable {
		var buf bytes.Buffer
		if err := EncodeJPEG(&buf, img, q); err != nil {
			t.Fatal(err)
		}
		return jpegDQT(t, buf.Bytes())
	}
	byQ := map[int][]dqtTable{}
	for _, tq := range g.Tables {
		byQ[tq.Quality] = tq.DQT
		if got := enc(tq.Quality); fmt.Sprint(got) != fmt.Sprint(tq.DQT) {
			t.Errorf("quality %d: EncodeJPEG tables %v, libjpeg-turbo %v", tq.Quality, got, tq.DQT)
		}
	}
	for q, w := range map[int]int{0: 1, -20: 1, 101: 100, 1000: 100} {
		if got := enc(q); fmt.Sprint(got) != fmt.Sprint(byQ[w]) {
			t.Errorf("quality %d: tables %v, want libjpeg-turbo's at quality %d", q, got, w)
		}
	}
}
