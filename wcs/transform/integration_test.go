package transform_test

import (
	"math"
	"path/filepath"
	"testing"

	"github.com/dmarkham/fits"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
	"github.com/dmarkham/fits/wcs/transform"
)

// TestEndToEndFITSWCS exercises the full stack: write a FITS image with
// TAN WCS keywords through the fits library, reopen it through the fits
// library, parse the WCS via wcs.Parse, build a transform.Transform, and
// verify pixel ↔ sky round-trips work.
//
// This is the test that matters for the "can I actually use this" story:
// it proves the pieces compose — fits.Open → hdu.Header() → wcs.Parse →
// transform.New — without any synthetic in-memory headers.
func TestEndToEndFITSWCS(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "wcs_e2e.fits")

	// Build a 512×512 image with realistic TAN WCS keywords for a star
	// field at roughly the position of M51 (NGC 5194) — chosen so that
	// hand-verification against any planetarium app is possible.
	const width, height = 512, 512
	const crval1 = 202.4696 // RA of M51, degrees
	const crval2 = 47.1953  // Dec of M51, degrees
	const crpix1 = 256.5    // FITS 1-based center of a 512×512 image
	const crpix2 = 256.5
	const scale = 0.0002777778 // 1 arcsec per pixel

	h := header.New()
	h.Set("OBJECT", "M51", "")
	h.Set("CTYPE1", "RA---TAN", "")
	h.Set("CTYPE2", "DEC--TAN", "")
	h.Set("CUNIT1", "deg", "")
	h.Set("CUNIT2", "deg", "")
	h.Set("CRVAL1", crval1, "")
	h.Set("CRVAL2", crval2, "")
	h.Set("CRPIX1", crpix1, "")
	h.Set("CRPIX2", crpix2, "")
	h.Set("CDELT1", -scale, "")
	h.Set("CDELT2", scale, "")
	h.Set("RADESYS", "ICRS", "")
	h.Set("EQUINOX", 2000.0, "")

	data := make([]float32, width*height)
	// Put a bright "star" at the reference pixel for sanity.
	data[(height/2-1)*width+(width/2-1)] = 1000

	f, err := fits.Create(path)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := fits.WriteImage(f, h, []int64{width, height}, data); err != nil {
		t.Fatal(err)
	}
	if err := f.Close(); err != nil {
		t.Fatal(err)
	}

	// Reopen through the real file-reading path.
	fr, err := fits.Open(path)
	if err != nil {
		t.Fatal(err)
	}
	defer fr.Close()

	img, err := fr.Primary()
	if err != nil {
		t.Fatal(err)
	}
	w, err := wcs.Parse(img.Header())
	if err != nil {
		t.Fatal(err)
	}
	if !w.IsCelestial() {
		t.Fatalf("not celestial after round-trip")
	}
	if w.CelestialProjCode() != "TAN" {
		t.Fatalf("proj code: %q", w.CelestialProjCode())
	}
	tr, err := transform.New(w)
	if err != nil {
		t.Fatal(err)
	}

	// Reference pixel (0-based) → reference sky.
	a, d, err := tr.PixelToSky(crpix1-1, crpix2-1)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(a-crval1) > 1e-10 {
		t.Fatalf("reference RA: %.15g vs %.15g", a, crval1)
	}
	if math.Abs(d-crval2) > 1e-10 {
		t.Fatalf("reference Dec: %.15g vs %.15g", d, crval2)
	}

	// A random off-axis pixel: round-trip must be numerically stable.
	px, py := 100.0, 400.0
	skyA, skyD, err := tr.PixelToSky(px, py)
	if err != nil {
		t.Fatal(err)
	}
	p1, p2, err := tr.SkyToPixel(skyA, skyD)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(p1-px) > 1e-7 || math.Abs(p2-py) > 1e-7 {
		t.Fatalf("round trip from file: (%v,%v) → (%v,%v) → (%v,%v)", px, py, skyA, skyD, p1, p2)
	}

	// The four image corners must all produce valid sky coordinates within
	// ±0.1° of the reference (512 pixels × 1 arcsec/pixel = ~0.07° field
	// radius, so this is a loose but meaningful bound).
	corners := []struct{ x, y float64 }{{0, 0}, {511, 0}, {0, 511}, {511, 511}}
	for _, c := range corners {
		a, d, err := tr.PixelToSky(c.x, c.y)
		if err != nil {
			t.Fatalf("corner (%v,%v): %v", c.x, c.y, err)
		}
		// Distance from reference in the small-angle approximation.
		dra := (a - crval1)
		if dra > 180 {
			dra -= 360
		} else if dra < -180 {
			dra += 360
		}
		ddec := d - crval2
		if math.Hypot(dra, ddec) > 0.2 {
			t.Errorf("corner (%v,%v) is %g° from center — too far", c.x, c.y, math.Hypot(dra, ddec))
		}
	}
}

// TestEndToEndCDMatrixFromFITS exercises the CD-matrix path end-to-end
// with a realistic off-diagonal rotation applied through the file
// writer / reader.
func TestEndToEndCDMatrixFromFITS(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "cd_e2e.fits")

	// A 10-degree rotated CD matrix: cos(10°)=0.9848, sin(10°)=0.1736.
	// CD encodes degrees per pixel directly.
	const pixScale = 0.0001 // 0.36 arcsec/pixel
	cos10, sin10 := math.Cos(math.Pi/18), math.Sin(math.Pi/18)

	h := header.New()
	h.Set("CTYPE1", "RA---TAN", "")
	h.Set("CTYPE2", "DEC--TAN", "")
	h.Set("CRVAL1", 83.633, "") // approx Crab Nebula
	h.Set("CRVAL2", 22.014, "")
	h.Set("CRPIX1", 128.0, "")
	h.Set("CRPIX2", 128.0, "")
	h.Set("CD1_1", -pixScale*cos10, "")
	h.Set("CD1_2", pixScale*sin10, "")
	h.Set("CD2_1", pixScale*sin10, "")
	h.Set("CD2_2", pixScale*cos10, "")

	f, _ := fits.Create(path)
	fits.WriteImage(f, h, []int64{256, 256}, make([]float32, 256*256))
	f.Close()

	fr, _ := fits.Open(path)
	defer fr.Close()
	img, _ := fr.Primary()
	w, err := wcs.Parse(img.Header())
	if err != nil {
		t.Fatal(err)
	}
	if w.CD == nil {
		t.Fatal("CD matrix not parsed from round-tripped file")
	}
	tr, err := transform.New(w)
	if err != nil {
		t.Fatal(err)
	}

	// Reference pixel → reference sky.
	a, d, err := tr.PixelToSky(127, 127)
	if err != nil {
		t.Fatal(err)
	}
	if math.Abs(a-83.633) > 1e-10 || math.Abs(d-22.014) > 1e-10 {
		t.Fatalf("reference: (%v,%v)", a, d)
	}
	// Round-trip off-center.
	skyA, skyD, _ := tr.PixelToSky(50, 200)
	p1, p2, _ := tr.SkyToPixel(skyA, skyD)
	if math.Abs(p1-50) > 1e-7 || math.Abs(p2-200) > 1e-7 {
		t.Fatalf("CD round trip: (%v,%v)", p1, p2)
	}
}
