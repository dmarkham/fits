package stretch

import (
	"fmt"
	"testing"
)

// CLAHE building blocks against opencv-python 5.0.0 (gen_clahe_cv2.py).
// All comparisons are exact: CLAHE is integer arithmetic plus float32
// interpolation rounded to integers, and the Lab conversions are
// deterministic float32 code.

func TestCLAHECore8bit_cv2(t *testing.T) {
	var g struct {
		Source       string    `json:"source"`
		Input        goldInput `json:"input"`
		Out          goldArray `json:"out"`
		Nx           int       `json:"nx"`
		Ny           int       `json:"ny"`
		ClipLimit    float64   `json:"clipLimit"`
		TileGridSize int       `json:"tileGridSize"`
	}
	readGolden(t, "clahe_mono_u8.json", &g)
	got := claheApply(g.Input.u16(t), g.Nx, g.Ny, g.TileGridSize, g.TileGridSize, g.ClipLimit, 256)
	checkU16(t, "8-bit CLAHE", g.Out, got)
}

func TestCLAHECore16bit_cv2(t *testing.T) {
	var g struct {
		Cases []struct {
			Nx        int       `json:"nx"`
			Ny        int       `json:"ny"`
			ClipLimit float64   `json:"clipLimit"`
			TilesX    int       `json:"tilesX"`
			TilesY    int       `json:"tilesY"`
			Input     goldInput `json:"input"`
			Out       goldArray `json:"out"`
		} `json:"cases"`
	}
	readGolden(t, "clahe_mono_u16.json", &g)
	if len(g.Cases) != 6 {
		t.Fatalf("%d cases, want 6", len(g.Cases))
	}
	for _, c := range g.Cases {
		got := claheApply(c.Input.u16(t), c.Nx, c.Ny, c.TilesX, c.TilesY, c.ClipLimit, 65536)
		checkU16(t, fmt.Sprintf("16-bit CLAHE %dx%d grid %dx%d clip %g", c.Nx, c.Ny, c.TilesX, c.TilesY, c.ClipLimit), c.Out, got)
	}
}

// convertRows runs a row converter over rows of `width` pixels (3 floats each).
func convertRows(in []float32, width int, row func(src, dst []float32)) []float32 {
	got := make([]float32, len(in))
	w := 3 * width
	for r := 0; r*w < len(in); r++ {
		row(in[r*w:(r+1)*w], got[r*w:(r+1)*w])
	}
	return got
}

type labCase struct {
	Width   int       `json:"width"`
	Grid    int       `json:"grid"`
	Vectors int       `json:"vectors"`
	Out     goldArray `json:"out"`
}

func loadLab(t *testing.T) (bgr, lab []float32, fwd, inv []labCase) {
	t.Helper()
	var g struct {
		BGRInput goldInput `json:"bgr_input"`
		LabInput goldInput `json:"lab_input"`
		Forward  []labCase `json:"forward"`
		Inverse  []labCase `json:"inverse"`
	}
	readGolden(t, "lab_cv2.json", &g)
	if len(g.Forward) == 0 || len(g.Inverse) == 0 {
		t.Fatal("lab_cv2.json has no vectors")
	}
	return g.BGRInput.flat(t), g.LabInput.flat(t), g.Forward, g.Inverse
}

// TestBGR2Lab_cv2 covers every node of OpenCV's 33^3 RGB->Lab LUT (which
// the Go port recomputes) and interpolated points, in rows of width 1, 13
// and 64 (scalar loop, SIMD body + tail).
func TestBGR2Lab_cv2(t *testing.T) {
	bgr, _, fwd, _ := loadLab(t)
	for _, v := range fwd {
		var in []float32
		if v.Grid > 0 {
			g := v.Grid
			in = make([]float32, 3*g*g*g)
			for p := 0; p < g; p++ {
				for q := 0; q < g; q++ {
					for r := 0; r < g; r++ {
						i := ((p*g+q)*g + r) * 3
						in[i], in[i+1], in[i+2] = float32(p)/32, float32(q)/32, float32(r)/32
					}
				}
			}
		} else {
			in = bgr[:3*v.Vectors]
		}
		checkF32(t, fmt.Sprintf("BGR2Lab width %d", v.Width), v.Out, convertRows(in, v.Width, bgr2labRow))
	}
}

func TestLab2BGR_cv2(t *testing.T) {
	_, lab, _, inv := loadLab(t)
	for _, v := range inv {
		checkF32(t, fmt.Sprintf("Lab2BGR width %d", v.Width), v.Out, convertRows(lab[:3*v.Vectors], v.Width, lab2bgrRow))
	}
}

// TestCLAHERGB_cv2Rebuild: the same Siril pipeline rebuilt from cv2 calls
// by gen_clahe_cv2.py (a second OpenCV build: the pip wheel).
func TestCLAHERGB_cv2Rebuild(t *testing.T) {
	var c srefCase
	readGolden(t, "stretch_clahe_cv2.json", &c)
	if c.Command != "clahe 3 8" {
		t.Fatalf("golden command %q", c.Command)
	}
	ch, nx, ny := c.input(t)
	CLAHERGB(ch[0], ch[1], ch[2], nx, ny, 8, 3.0)
	c.check(t, "clahe_cv2", ch)
}

func TestCLAHE_BadArgsPanic(t *testing.T) {
	for _, f := range []func(){
		func() { CLAHE(make([]float32, 4), 2, 2, 0, 2) },
		func() { CLAHE(make([]float32, 5), 2, 2, 2, 2) },
		func() { CLAHERGB(make([]float32, 4), make([]float32, 3), make([]float32, 4), 2, 2, 2, 2) },
	} {
		func() {
			defer func() {
				if recover() == nil {
					t.Error("expected panic")
				}
			}()
			f()
		}()
	}
}

// TestCLAHE_SirilOracle: Siril's `clahe` on float images
// (opencv/opencv.cpp:806 cvClahe_float, verbatim, linked against OpenCV
// 5.0.0). Includes mono, non-square tiles counts and image sizes not
// divisible by the tile grid (OpenCV pads with BORDER_REFLECT_101). In the
// small cases OpenCV's integer clip limit int(clip*tileTotal/65536) is
// below 1 and clamps to 1; the clahe_clip_* cases have tiles big enough for
// a limit of 3, and the generator checked that limits 2 and 4 change the
// output, so they pin the clip-limit arithmetic.
func TestCLAHE_SirilOracle(t *testing.T) {
	for _, name := range []string{
		"clahe_rgb_3_8", "clahe_rgb_2_4", "clahe_rgb_odd", "clahe_bright_rgb", "clahe_mono_3_8",
		"clahe_mono_odd", "clahe_clip_mono", "clahe_clip_rgb", "clahe_clip_rgb_pad",
	} {
		runSrefCase(t, name, func(c srefCase, ch [][]float32, nx, ny int) {
			if c.OpenCV != "5.0.0" {
				t.Fatalf("golden from OpenCV %q, the port follows 5.0.0", c.OpenCV)
			}
			var clip float64
			var tile int
			if n, err := fmt.Sscanf(c.Command, "clahe %g %d", &clip, &tile); n != 2 {
				t.Fatalf("command %q: %v", c.Command, err)
			}
			if len(ch) == 1 {
				CLAHE(ch[0], nx, ny, tile, clip)
				return
			}
			CLAHERGB(ch[0], ch[1], ch[2], nx, ny, tile, clip)
		})
	}
}

// TestCLAHECoreClip_cv2: OpenCV's 8- and 16-bit CLAHE where the integer
// clip limit is >= 2 (gen_clahe_cv2.py checked that cv2 at limit-1 and
// limit+1 gives a different output for every case). Exact.
func TestCLAHECoreClip_cv2(t *testing.T) {
	var g struct {
		Cases []struct {
			Bits      int       `json:"bits"`
			Nx        int       `json:"nx"`
			Ny        int       `json:"ny"`
			ClipLimit float64   `json:"clipLimit"`
			TilesX    int       `json:"tilesX"`
			TilesY    int       `json:"tilesY"`
			Limit     int       `json:"clip_limit_int"`
			Input     goldInput `json:"input"`
			Out       goldArray `json:"out"`
		} `json:"cases"`
	}
	readGolden(t, "clahe_clip_cv2.json", &g)
	if len(g.Cases) != 6 {
		t.Fatalf("%d cases, want 6", len(g.Cases))
	}
	for _, c := range g.Cases {
		label := fmt.Sprintf("%d-bit %dx%d grid %dx%d clip %g (limit %d)", c.Bits, c.Nx, c.Ny, c.TilesX, c.TilesY, c.ClipLimit, c.Limit)
		if c.Limit < 2 {
			t.Fatalf("%s: case does not exercise the clip limit", label)
		}
		hist := 256
		if c.Bits == 16 {
			hist = 65536
		}
		checkU16(t, label, c.Out, claheApply(c.Input.u16(t), c.Nx, c.Ny, c.TilesX, c.TilesY, c.ClipLimit, hist))
	}
}
