package stretch

import (
	"fmt"
	"image"
	"image/color"
	"image/jpeg"
	"io"
)

// Preset identifies a named stretch configuration matching the
// fitview presets. Each preset maps to a specific stretch function
// with fixed parameters.
type Preset struct {
	Name    string
	Label   string
	ApplyFn func(r, g, b []float32, nx, ny int)
}

// DefaultPresets returns the 18 standard presets matching fitview's
// siril.go preset list.
func DefaultPresets() []Preset {
	autostretchFn := func(sc, tb float32) func(r, g, b []float32, nx, ny int) {
		return func(r, g, b []float32, nx, ny int) {
			Autostretch(r, sc, tb)
			Autostretch(g, sc, tb)
			Autostretch(b, sc, tb)
		}
	}
	return []Preset{
		{"kstars1", "KStars 1 (Default)", autostretchFn(-2.8, 0.25)},
		{"kstars2", "KStars 2 (Bright BG)", autostretchFn(-1.5, 0.25)},
		{"kstars3", "KStars 3 (Much Brighter)", autostretchFn(-1.0, 0.125)},
		{"kstars4", "KStars 4 (Very Bright)", autostretchFn(-0.5, 0.125)},
		{"kstars5", "KStars 5 (Dark BG)", autostretchFn(-2.0, 0.125)},
		{"kstars6", "KStars 6 (Very Dark)", autostretchFn(-4.0, 0.125)},
		{"kstars7", "KStars 7 (Extreme Dark)", autostretchFn(-5.5, 0.25)},
		{"custom_dark", "Custom Dark", autostretchFn(-6.5, 0.03)},
		{"custom_extreme", "Custom Extreme", autostretchFn(-7.5, 0.02)},
		{"custom_ultra", "Custom Ultra Dark", autostretchFn(-9.0, 0.01)},
		{"custom_mask", "Custom Mask Only", autostretchFn(-12.0, 0.005)},
		{"auto_linked", "Auto Linked", func(r, g, b []float32, nx, ny int) {
			AutostretchLinked([][]float32{r, g, b}, -2.8, 0.25)
		}},
		{"asinh_light", "Asinh Light", func(r, g, b []float32, nx, ny int) {
			AsinhRGB(r, g, b, 1, 0)
		}},
		{"asinh_medium", "Asinh Medium", func(r, g, b []float32, nx, ny int) {
			AsinhRGB(r, g, b, 50, 0)
		}},
		{"asinh_strong", "Asinh Strong", func(r, g, b []float32, nx, ny int) {
			AsinhRGB(r, g, b, 500, 0)
		}},
		{"ght_mild", "GHT Mild", func(r, g, b []float32, nx, ny int) {
			p := GHTParams{D: 2.0, B: 0.0, SP: 0.5, HP: 1.0}
			GHT(r, p)
			GHT(g, p)
			GHT(b, p)
		}},
		{"ght_strong", "GHT Strong", func(r, g, b []float32, nx, ny int) {
			p := GHTParams{D: 5.0, B: 2.0, SP: 0.3, HP: 1.0}
			GHT(r, p)
			GHT(g, p)
			GHT(b, p)
		}},
		{"linear", "Linear", func(r, g, b []float32, nx, ny int) {
			Linear(r, 0.05)
			Linear(g, 0.05)
			Linear(b, 0.05)
		}},
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

// ToRGBImage converts stretched float32 channels to a standard Go
// image.RGBA. Each pixel is quantized from [0, 1] to [0, 255].
func ToRGBImage(r, g, b []float32, nx, ny int) *image.RGBA {
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

// EncodeJPEG writes the RGB image as a JPEG to w with the given quality.
func EncodeJPEG(w io.Writer, img *image.RGBA, quality int) error {
	return jpeg.Encode(w, img, &jpeg.Options{Quality: quality})
}

func quantize8(v float32) uint8 {
	if v <= 0 {
		return 0
	}
	if v >= 1 {
		return 255
	}
	return uint8(v*255 + 0.5)
}
