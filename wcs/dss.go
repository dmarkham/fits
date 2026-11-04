package wcs

import (
	"fmt"
	"math"
	"regexp"
	"strconv"
	"strings"

	"github.com/dmarkham/fits/header"
)

// DSS holds the Digitized Sky Survey plate solution of a DSS image header
// (STScI GSSS: CNPIX1/2, PPO3/PPO6, XPIXELSZ/YPIXELSZ, PLTRAH/M/S,
// PLTDECSN/D/M/S, AMDX1..13, AMDY1..13), as read by wcslib 8.10 wcspih.
//
// When a header carries DSS keywords, wcslib does not use the header's own
// linear WCS: following Paper IV sect. 5.2 it rebuilds CRPIX, PC, CDELT,
// CRVAL and CTYPE (RA---TAN / DEC--TAN) of the primary description from
// the plate solution, erases CD, and applies the DSS plate polynomial as a
// sequent distortion (between the matrix and CDELT). Modern STScI cutouts
// carry both the plate solution and an approximate TAN WCS (CD, CROTA,
// PC00i00j); the approximate WCS is replaced. Parse does the same: the
// Header it returns for such a header holds the translated CTYPE, CRPIX,
// CRVAL, CDELT and PC (CD nil, CROTA not used) and the coefficients here,
// which wcs/transform applies.
type DSS struct {
	CNPix     [2]float64  // CNPIX1, CNPIX2: corner of the cutout on the plate, pixels
	PPO3      float64     // plate centre x, microns
	PPO6      float64     // plate centre y, microns
	PixelSize [2]float64  // XPIXELSZ, YPIXELSZ, microns
	RA        [3]float64  // PLTRAH, PLTRAM, PLTRAS
	DecSign   float64     // PLTDECSN: +1 or -1
	Dec       [3]float64  // PLTDECD, PLTDECM, PLTDECS
	AMDX      [13]float64 // AMDX1..AMDX13 (index 0 is AMDX1)
	AMDY      [13]float64 // AMDY1..AMDY13
}

// dssPlateKeys are the plate keywords wcslib reads into its DSS
// translation (wcspih.l); every one is required.
var dssPlateKeys = []string{
	"CNPIX1", "CNPIX2", "PPO3", "PPO6", "XPIXELSZ", "YPIXELSZ",
	"PLTRAH", "PLTRAM", "PLTRAS", "PLTDECSN", "PLTDECD", "PLTDECM", "PLTDECS",
}

var dssAMDRe = regexp.MustCompile(`^AMD([XY])([0-9]+)$`)

// HasDSSKeywords reports whether h carries any keyword that makes wcslib
// 8.10 apply the DSS translation: one of dssPlateKeys or an AMDXm / AMDYm
// coefficient. (PLATEID alone does not, as in wcslib.)
func HasDSSKeywords(h *header.Header) bool {
	for _, c := range h.Cards() {
		if dssAMDRe.MatchString(c.Key) {
			return true
		}
		for _, k := range dssPlateKeys {
			if c.Key == k {
				return true
			}
		}
	}
	return false
}

// parseDSS reads the DSS plate solution. It returns nil, nil when h has no
// DSS keyword. It is stricter than wcslib 8.10, and returns an error, when:
//
//   - any of the 13 plate keywords is missing or not numeric (wcslib then
//     reads an uninitialised value);
//   - PLTDECSN is not '+' or '-', or is '-' written without a trailing
//     blank inside the quotes (wcslib's sscanf then reads "-'" and takes
//     the declination as north);
//   - an AMDXm / AMDYm index is not 1..99 without a leading zero (wcslib
//     discards such a card) or a coefficient is not numeric;
//   - a NON-ZERO coefficient has m > 13 (wcslib refuses it in wcsset:
//     "Invalid parameter"; zero ones are skipped, as in wcslib);
//   - the plate scale sqrt(|A1 B1 - A2 B2|) is zero (wcslib refuses it).
func parseDSS(h *header.Header) (*DSS, error) {
	if !HasDSSKeywords(h) {
		return nil, nil
	}
	d := &DSS{}
	num := func(k string) (float64, error) {
		v, ok, err := getFloat(h, k)
		if err != nil {
			return 0, err
		}
		if !ok {
			return 0, fmt.Errorf("fits/wcs: DSS plate solution keyword %s is missing (wcslib 8.10 would read an uninitialised value)", k)
		}
		return v, nil
	}
	dst := []*float64{&d.CNPix[0], &d.CNPix[1], &d.PPO3, &d.PPO6, &d.PixelSize[0], &d.PixelSize[1],
		&d.RA[0], &d.RA[1], &d.RA[2], nil, &d.Dec[0], &d.Dec[1], &d.Dec[2]}
	for i, k := range dssPlateKeys {
		if dst[i] == nil {
			continue
		}
		v, err := num(k)
		if err != nil {
			return nil, err
		}
		*dst[i] = v
	}
	sn, ok, err := getString(h, "PLTDECSN")
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, fmt.Errorf("fits/wcs: DSS plate solution keyword PLTDECSN is missing")
	}
	switch sn {
	case "+":
		d.DecSign = 1
	case "-":
		d.DecSign = -1
		if c, _ := h.Get("PLTDECSN"); c.Raw != ([header.CardWidth]byte{}) && strings.Contains(string(c.Raw[:]), "'-'") {
			return nil, fmt.Errorf("fits/wcs: PLTDECSN = '-' without a trailing blank: wcslib 8.10 reads it as '+' (its sscanf takes \"-'\"); write '- '")
		}
	default:
		return nil, fmt.Errorf("fits/wcs: PLTDECSN = %q, want '+' or '-'", sn)
	}
	for _, c := range h.Cards() {
		m := dssAMDRe.FindStringSubmatch(c.Key)
		if m == nil {
			continue
		}
		if m[2][0] == '0' || len(m[2]) > 2 {
			return nil, fmt.Errorf("fits/wcs: %s is not a DSS coefficient keyword (index 1..99 without a leading zero; wcslib 8.10 ignores this card)", c.Key)
		}
		v, ok := floatValue(c.Value)
		if !ok {
			return nil, fmt.Errorf("fits/wcs: %s has a non-numeric value %v", c.Key, c.Value)
		}
		n, _ := strconv.Atoi(m[2])
		if n > 13 {
			if v != 0 {
				return nil, fmt.Errorf("fits/wcs: %s = %g: the DSS plate polynomial has 13 terms (wcslib 8.10 refuses it too)", c.Key, v)
			}
			continue
		}
		if m[1] == "X" {
			d.AMDX[n-1] = v
		} else {
			d.AMDY[n-1] = v
		}
	}
	if d.Scale() == 0 {
		return nil, fmt.Errorf("fits/wcs: DSS plate scale sqrt(|AMDX1*AMDY1 - AMDX2*AMDY2|) is zero")
	}
	return d, nil
}

// Scale is the DSS plate scale S = sqrt(|A1 B1 - A2 B2|) (arcsec per mm),
// Paper IV sect. 5.2.
func (d *DSS) Scale() float64 {
	a1, a2, b1, b2 := d.AMDX[0], d.AMDX[1], d.AMDY[0], d.AMDY[1]
	return math.Sqrt(math.Abs(float64(a1*b1) - float64(a2*b2)))
}

// applyDSS replaces the primary description's linear WCS with the one
// wcslib 8.10 wcspih_final derives from the plate solution (Paper IV
// sect. 5.2, same variable names).
func (w *Header) applyDSS(d *DSS) {
	xc, yc := d.PPO3/1000, d.PPO6/1000
	rx, ry := d.PixelSize[0]/1000, d.PixelSize[1]/1000
	a1, a2, a3 := d.AMDX[0], d.AMDX[1], d.AMDX[2]
	b1, b2, b3 := d.AMDY[0], d.AMDY[1], d.AMDY[2]
	s := d.Scale()
	det := float64(a1*b1) - float64(a2*b2)
	x0 := (float64(a2*b3) - float64(a3*b1)) / det
	y0 := (float64(a3*b2) - float64(a1*b3)) / det

	w.CRPix[0] = (xc-x0)/rx - (d.CNPix[0] - 0.5)
	w.CRPix[1] = (yc+y0)/ry - (d.CNPix[1] - 0.5)
	w.PC[0][0] = float64(a1*rx) / s
	w.PC[0][1] = float64(-a2*ry) / s
	w.PC[1][0] = float64(-b2*rx) / s
	w.PC[1][1] = float64(b1*ry) / s
	w.HasPC = true
	w.CD = nil
	w.CDelt[0] = -s / 3600
	w.CDelt[1] = s / 3600
	w.CRVal[0] = float64((d.RA[0] + (d.RA[1]+d.RA[2]/60)/60) * 15)
	w.CRVal[1] = d.Dec[0] + (d.Dec[1]+d.Dec[2]/60)/60
	if d.DecSign == -1 {
		w.CRVal[1] *= -1
	}
	w.CType[0], w.CType[1] = "RA---TAN", "DEC--TAN"
	w.DSS = d
}
