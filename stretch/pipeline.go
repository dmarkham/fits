package stretch

import (
	"fmt"
	"image"
	"image/color"
	"image/jpeg"
	"io"
)

// Preset is a named stretch configuration for an RGB image. Each of the
// DefaultPresets is one Siril command (Command) applied to a 3-channel
// float image.
type Preset struct {
	Name    string
	Label   string
	ApplyFn func(r, g, b []float32, nx, ny int)
	Command string // the siril-cli command the preset reproduces
}

// DefaultPresets returns the 18 standard presets (the stretch list of
// fitview, github.com/dmarkham/fitview). Each preset applies its Command
// the way siril-cli does on a 32-bit RGB image; the tests check every
// preset against a siril-cli golden for that command.
func DefaultPresets() []Preset {
	autostretch := func(sc, tb float32) func(r, g, b []float32, nx, ny int) {
		return func(r, g, b []float32, nx, ny int) {
			// the preset arguments are in range (TestDefaultPresets_Siril)
			_ = AutostretchUnlinked([][]float32{r, g, b}, sc, tb)
		}
	}
	ght := func(p GHTParams) func(r, g, b []float32, nx, ny int) {
		return func(r, g, b []float32, nx, ny int) {
			// the preset parameters are in range (TestDefaultPresets_Siril)
			_ = GHT(r, p)
			_ = GHT(g, p)
			_ = GHT(b, p)
		}
	}
	return []Preset{
		{"kstars1", "KStars 1 (Default)", autostretch(-2.8, 0.25), "autostretch -2.8 0.25"},
		{"kstars2", "KStars 2 (Bright BG)", autostretch(-1.5, 0.25), "autostretch -1.5 0.25"},
		{"kstars3", "KStars 3 (Much Brighter)", autostretch(-1.0, 0.125), "autostretch -1.0 0.125"},
		{"kstars4", "KStars 4 (Very Bright)", autostretch(-0.5, 0.125), "autostretch -0.5 0.125"},
		{"kstars5", "KStars 5 (Dark BG)", autostretch(-2.0, 0.125), "autostretch -2.0 0.125"},
		{"kstars6", "KStars 6 (Very Dark)", autostretch(-4.0, 0.125), "autostretch -4.0 0.125"},
		{"kstars7", "KStars 7 (Extreme Dark)", autostretch(-5.5, 0.25), "autostretch -5.5 0.25"},
		{"custom_dark", "Custom Dark", autostretch(-6.5, 0.03), "autostretch -6.5 0.03"},
		{"custom_extreme", "Custom Extreme", autostretch(-7.5, 0.02), "autostretch -7.5 0.02"},
		{"custom_ultra", "Custom Ultra Dark", autostretch(-9.0, 0.01), "autostretch -9.0 0.01"},
		{"custom_mask", "Custom Mask Only", autostretch(-12.0, 0.005), "autostretch -12.0 0.005"},
		{"auto_linked", "Auto Linked", func(r, g, b []float32, nx, ny int) {
			_ = AutostretchLinked([][]float32{r, g, b}, -2.8, 0.25)
		}, "autostretch -linked"},
		{"asinh_light", "Asinh Light", func(r, g, b []float32, nx, ny int) {
			_ = AsinhRGB(r, g, b, 1, 0) // in range (TestDefaultPresets_Siril)
		}, "asinh 1 0"},
		{"asinh_medium", "Asinh Medium", func(r, g, b []float32, nx, ny int) {
			_ = AsinhRGB(r, g, b, 50, 0) // in range (TestDefaultPresets_Siril)
		}, "asinh 50 0"},
		{"asinh_strong", "Asinh Strong", func(r, g, b []float32, nx, ny int) {
			_ = AsinhRGB(r, g, b, 500, 0) // in range (TestDefaultPresets_Siril)
		}, "asinh 500 0"},
		{"ght_mild", "GHT Mild", ght(GHTParams{D: 2.0, B: 0.0, SP: 0.5, HP: 1.0}), "ght -D=2.0 -B=0.0 -SP=0.5"},
		{"ght_strong", "GHT Strong", ght(GHTParams{D: 5.0, B: 2.0, SP: 0.3, HP: 1.0}), "ght -D=5.0 -B=2.0 -SP=0.3"},
		{"linear", "Linear", func(r, g, b []float32, nx, ny int) {
			_ = Linear(r, 0.05)
			_ = Linear(g, 0.05)
			_ = Linear(b, 0.05)
		}, "linstretch -BP=0.05"},
	}
}

// FindPreset returns the preset with the given name, or an error if
// not found.
func FindPreset(name string) (Preset, error) {
	for _, p := range DefaultPresets() {
		if p.Name == name {
			return p, nil
		}
	}
	return Preset{}, fmt.Errorf("stretch: unknown preset %q", name)
}

// ToRGBImage converts stretched float32 channels to an 8-bit image.RGBA
// (alpha 255). Pixel i = y*nx + x of the channels becomes image point
// (x, y): row 0 of the data is the top row of the image, with no flip.
// (FITS images are usually stored bottom row first, so a FITS image
// rendered this way appears upside down relative to DS9 or Siril unless
// the caller flips the rows.) Each value is clamped to [0, 1] and mapped
// to floor(v*255 + 0.5) in float32; NaN maps to 0.
//
// It panics if a channel holds fewer than nx*ny values.
func ToRGBImage(r, g, b []float32, nx, ny int) *image.RGBA {
	n := nx * ny
	if len(r) < n || len(g) < n || len(b) < n {
		panic(fmt.Sprintf("stretch: ToRGBImage %dx%d needs %d pixels per channel, got %d, %d, %d", nx, ny, n, len(r), len(g), len(b)))
	}
	img := image.NewRGBA(image.Rect(0, 0, nx, ny))
	for y := 0; y < ny; y++ {
		for x := 0; x < nx; x++ {
			i := y*nx + x
			img.SetRGBA(x, y, color.RGBA{
				R: quantize8(r[i]),
				G: quantize8(g[i]),
				B: quantize8(b[i]),
				A: 255,
			})
		}
	}
	return img
}

// EncodeJPEG writes img to w as a baseline JPEG with the given quality
// (1-100; out-of-range values are clamped by image/jpeg).
func EncodeJPEG(w io.Writer, img *image.RGBA, quality int) error {
	return jpeg.Encode(w, img, &jpeg.Options{Quality: quality})
}

func quantize8(v float32) uint8 {
	if !(v > 0) { // also NaN
		return 0
	}
	if v >= 1 {
		return 255
	}
	return uint8(float32(v*255) + 0.5)
}
