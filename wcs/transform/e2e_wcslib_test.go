package transform

import (
	"bytes"
	"fmt"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/wcs"
)

// TestEndToEndFileAgainstWcslib: the whole user path -- header cards
// written into a real FITS file with fits.WriteImage, reopened with
// fits.Open, wcs.ParseAlt on the HDU header, New, PixelToSky/SkyToPixel --
// must reproduce wcslib's answers for a spread of golden cases (plain TAN,
// a real JWST SIP order-5 header, a real ZTF TPV header, a FREQ/RA/DEC
// cube, swapped axes, an alternate WCS). Structural keywords are written
// by WriteImage itself; the image is 2 pixels per axis because the pixel
// values play no part in the WCS.
func TestEndToEndFileAgainstWcslib(t *testing.T) {
	for _, name := range []string{
		"tan_cd_rot", "real_jwst_nircam_sip5", "real_ztf_tpv4",
		"axes_cube_freq_ra_dec", "axes_swapped", "alt_a_celestial", "lonpole_glon_ait_lp0",
	} {
		t.Run(name, func(t *testing.T) {
			g := loadHdrGolden(t, name)
			var buf bytes.Buffer
			naxis := 0
			for _, c := range g.Cards {
				buf.WriteString(c)
				if strings.HasPrefix(c, "NAXIS   =") {
					n, err := strconv.Atoi(strings.TrimSpace(strings.SplitN(c[10:], "/", 2)[0]))
					if err != nil {
						t.Fatalf("NAXIS card %q: %v", c, err)
					}
					naxis = n
				}
			}
			buf.WriteString(fmt.Sprintf("%-80s", "END"))
			for buf.Len()%2880 != 0 {
				buf.WriteByte(' ')
			}
			cards, _, err := header.ParseCards(buf.Bytes())
			if err != nil {
				t.Fatal(err)
			}
			if naxis == 0 {
				t.Fatalf("case %s has no NAXIS card", name)
			}
			shape := make([]int64, naxis)
			npix := 1
			for i := range shape {
				shape[i] = 2
				npix *= 2
			}
			path := filepath.Join(t.TempDir(), name+".fits")
			f, err := fits.Create(path)
			if err != nil {
				t.Fatal(err)
			}
			if _, err := fits.WriteImage(f, header.FromCards(cards), shape, make([]float32, npix)); err != nil {
				t.Fatal(err)
			}
			if err := f.Close(); err != nil {
				t.Fatal(err)
			}
			fr, err := fits.Open(path)
			if err != nil {
				t.Fatal(err)
			}
			defer fr.Close()
			img, err := fr.Primary()
			if err != nil {
				t.Fatal(err)
			}
			w, err := wcs.ParseAlt(img.Header(), g.Alt)
			if err != nil {
				t.Fatal(err)
			}
			tr, err := New(w)
			if err != nil {
				t.Fatal(err)
			}
			tolAS, tolPx := hdrDefaultTolArcsec, hdrDefaultTolPx
			if g.TolArcsec != nil {
				tolAS = *g.TolArcsec
			}
			if g.TolPx != nil {
				tolPx = *g.TolPx
			}
			p := hdrCheckP2S(g, tr, tolAS)
			s := hdrCheckS2P(g, tr, tolPx)
			for i, m := range append(p.mismatches, s.mismatches...) {
				if i == 10 {
					t.Error("...")
					break
				}
				t.Error(m)
			}
			t.Logf("via FITS file: p2s %d points max %.3g arcsec, s2p %d points max %.3g px", p.n, p.maxErr, s.n, s.maxErr)
		})
	}
}
