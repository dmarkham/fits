package healpix_test

import (
	"encoding/json"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/dmarkham/fits"
	"github.com/dmarkham/fits/header"
	"github.com/dmarkham/fits/healpix"
)

// These tests check the package's stated use: indexing HEALPix maps that
// are stored as FITS binary tables. See testdata/PROVENANCE.md.

func readCommitted(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("committed sample missing (this is a failure, not a skip): %v", err)
	}
	return b
}

func hpxScheme(t *testing.T, ordering string) healpix.Scheme {
	t.Helper()
	switch strings.TrimSpace(ordering) {
	case "NESTED":
		return healpix.Nested
	case "RING":
		return healpix.Ring
	}
	t.Fatalf("ORDERING %q is neither RING nor NESTED", ordering)
	return 0
}

// TestRealMapHeader_WMAP parses the real header of the WMAP 9-year ILC
// map (LAMBDA, both header blocks verbatim) and checks that its HEALPix
// keywords are self-consistent under this package: NSIDE gives NPIX and
// LASTPIX, NESTED is allowed for it, and the table row count equals the
// pixel count.
func TestRealMapHeader_WMAP(t *testing.T) {
	raw := readCommitted(t, filepath.Join("testdata", "maps", "wmap_ilc_9yr_v5.hdr"))
	if len(raw) != 2*2880 {
		t.Fatalf("header sample is %d bytes, want 5760 (primary + BINTABLE header blocks)", len(raw))
	}
	prim, _, err := header.ParseCards(raw[:2880])
	if err != nil {
		t.Fatalf("primary header: %v", err)
	}
	ext, _, err := header.ParseCards(raw[2880:])
	if err != nil {
		t.Fatalf("extension header: %v", err)
	}
	for name, h := range map[string]*header.Header{"primary": header.FromCards(prim), "BINTABLE": header.FromCards(ext)} {
		pixtype, err := h.String("PIXTYPE")
		if err != nil || strings.TrimSpace(pixtype) != "HEALPIX" {
			t.Fatalf("%s: PIXTYPE = %q, %v", name, pixtype, err)
		}
		ord, err := h.String("ORDERING")
		if err != nil {
			t.Fatalf("%s: ORDERING: %v", name, err)
		}
		scheme := hpxScheme(t, ord)
		nside, err := h.Int("NSIDE")
		if err != nil {
			t.Fatalf("%s: NSIDE: %v", name, err)
		}
		first, err1 := h.Int("FIRSTPIX")
		last, err2 := h.Int("LASTPIX")
		if err1 != nil || err2 != nil {
			t.Fatalf("%s: FIRSTPIX/LASTPIX: %v %v", name, err1, err2)
		}
		if nside != 512 || scheme != healpix.Nested {
			t.Fatalf("%s: NSIDE=%d ORDERING=%v, the published map is NSIDE 512 NESTED", name, nside, scheme)
		}
		if !healpix.NsideOK(int(nside), scheme) {
			t.Errorf("%s: NsideOK(%d, %v) = false", name, nside, scheme)
		}
		npix := healpix.Nside2Npix(int(nside))
		if first != 0 || last != npix-1 {
			t.Errorf("%s: FIRSTPIX..LASTPIX = %d..%d, want 0..%d", name, first, last, npix-1)
		}
		if got, ok := healpix.Npix2Nside(last + 1); !ok || int64(got) != nside {
			t.Errorf("%s: Npix2Nside(LASTPIX+1) = %d, %v", name, got, ok)
		}
		if v, err := h.Int("NPIX"); err == nil && v != npix {
			t.Errorf("%s: NPIX = %d, Nside2Npix = %d", name, v, npix)
		}
		if name == "BINTABLE" {
			rows, err := h.Int("NAXIS2")
			if err != nil || rows != npix {
				t.Errorf("BINTABLE NAXIS2 = %d (%v), want one row per pixel (%d)", rows, err, npix)
			}
		}
	}
}

// TestHealpyWrittenMap reads maps written by healpy.write_map (RING and
// NESTED, 1024-element vector rows) through the fits package, then checks
// every pixel's stored value against Pix2Vec and looks up seeded-random
// directions with Ang2Pix, comparing pixel and value with healpy.
func TestHealpyWrittenMap(t *testing.T) {
	var g struct {
		Nside int       `json:"nside"`
		Theta []float64 `json:"theta"`
		Phi   []float64 `json:"phi"`
		Nest  struct {
			Pix   []int64   `json:"pix"`
			Value []float64 `json:"value"`
		} `json:"nest"`
		Ring struct {
			Pix   []int64   `json:"pix"`
			Value []float64 `json:"value"`
		} `json:"ring"`
	}
	if err := json.Unmarshal(readCommitted(t, filepath.Join("testdata", "maps", "healpy_map_nside16.json")), &g); err != nil {
		t.Fatal(err)
	}
	for _, tc := range []struct {
		file string
		pix  []int64
		val  []float64
	}{
		{"healpy_map_nside16_nest.fits", g.Nest.Pix, g.Nest.Value},
		{"healpy_map_nside16_ring.fits", g.Ring.Pix, g.Ring.Value},
	} {
		t.Run(tc.file, func(t *testing.T) {
			path := filepath.Join("testdata", "maps", tc.file)
			readCommitted(t, path)
			f, err := fits.Open(path)
			if err != nil {
				t.Fatal(err)
			}
			defer f.Close()
			h, err := f.HDU(1)
			if err != nil {
				t.Fatal(err)
			}
			tbl, ok := h.(*fits.BinaryTableHDU)
			if !ok {
				t.Fatalf("HDU 1 is %T, want a binary table", h)
			}
			hdr := tbl.Header()
			ord, _ := hdr.String("ORDERING")
			scheme := hpxScheme(t, ord)
			nside64, err := hdr.Int("NSIDE")
			if err != nil {
				t.Fatal(err)
			}
			nside := int(nside64)
			if nside != g.Nside {
				t.Fatalf("NSIDE %d, golden %d", nside, g.Nside)
			}
			m, err := fits.ReadVectorColumn[float64](tbl, 1) // FITS columns are 1-based
			if err != nil {
				t.Fatal(err)
			}
			if int64(len(m)) != healpix.Nside2Npix(nside) {
				t.Fatalf("map has %d values, want %d", len(m), healpix.Nside2Npix(nside))
			}
			// Every stored value is the z of its pixel centre (healpy
			// pix2vec), so this checks file order == pixel order.
			for p := range m {
				_, _, z := healpix.Pix2Vec(nside, int64(p), scheme)
				if d := z - m[p]; d > 2*0x1p-52 || d < -2*0x1p-52 {
					t.Fatalf("pixel %d: stored z %.17g, Pix2Vec z %.17g", p, m[p], z)
				}
			}
			for i := range g.Theta {
				p := healpix.Ang2Pix(nside, g.Theta[i], g.Phi[i], scheme)
				if p != tc.pix[i] {
					t.Errorf("direction %d: Ang2Pix %d, healpy %d", i, p, tc.pix[i])
					continue
				}
				if m[p] != tc.val[i] {
					t.Errorf("direction %d: map[%d] = %.17g, healpy read_map %.17g", i, p, m[p], tc.val[i])
				}
			}
		})
	}
}
