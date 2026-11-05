// fitsdump prints a human-readable summary of every HDU in a FITS file:
// its type, the image BITPIX/NAXIS/shape (the logical values for a
// tile-compressed image), the binary-table row count and columns (name,
// TFORM, unit) or the ASCII-table row count, and every header keyword with
// its decoded value and comment. Long string values continued with
// CONTINUE cards are printed as one keyword with the joined value; a
// keyword with no value is printed like a commentary card. Pixel and table
// values are not printed. An HDU, a compressed image's metadata or a
// table's columns that cannot be decoded are reported in place ("error:
// ..."), the rest of the file is still printed, and the exit status is 1
// with the file and the problems named on stderr. Every header of a file
// that fits.Open accepts parses (fits TestOpenedHeadersAlwaysParse), so a
// header itself is never an undecodable part; a file fits.Open refuses is
// reported as a whole.
//
// Usage:
//
//	fitsdump file.fits [file.fits ...]
package main

import (
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/dmarkham/fits"
)

func main() {
	os.Exit(run(os.Args[1:], os.Stdout, os.Stderr))
}

// run is main without the process exit, so tests can drive it. It returns
// the exit status: 0 on success, 1 if any file failed, 2 on bad usage.
func run(args []string, stdout, stderr io.Writer) int {
	if len(args) < 1 {
		fmt.Fprintln(stderr, "usage: fitsdump FILE [FILE ...]")
		return 2
	}
	exit := 0
	for _, path := range args {
		if err := dump(stdout, path); err != nil {
			fmt.Fprintf(stderr, "fitsdump: %s: %v\n", path, err)
			exit = 1
		}
	}
	return exit
}

// hduOf is f.HDU; tests replace it to reach the error report for an HDU
// that cannot be decoded.
var hduOf = (*fits.File).HDU

func dump(w io.Writer, path string) error {
	f, err := fits.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintf(w, "=== %s ===\n", path)
	fmt.Fprintf(w, "HDUs: %d\n\n", f.NumHDU())

	var problems []string
	for i := range f.NumHDU() {
		h, err := hduOf(f, i)
		if err != nil {
			fmt.Fprintf(w, "HDU %d: error: %v\n", i, err)
			problems = append(problems, fmt.Sprintf("HDU %d: %v", i, err))
			continue
		}
		fmt.Fprintf(w, "--- HDU %d: %s ---\n", i, h.Type())
		if h.Compressed() {
			fmt.Fprintf(w, "  Compressed: %s\n", h.CompressionType())
		}
		switch v := h.(type) {
		case *fits.ImageHDU:
			fmt.Fprintf(w, "  BITPIX = %d\n", v.BITPIX())
			fmt.Fprintf(w, "  NAXIS  = %d\n", v.NAXIS())
			if v.NAXIS() > 0 {
				fmt.Fprintf(w, "  Shape  = %v\n", v.Shape())
			}
		case *fits.CompressedImageHDU:
			if err := v.Validate(); err != nil {
				fmt.Fprintf(w, "  Compressed image: error: %v\n", err)
				problems = append(problems, fmt.Sprintf("HDU %d compressed image: %v", i, err))
				break
			}
			fmt.Fprintf(w, "  BITPIX = %d (logical, from ZBITPIX)\n", v.BITPIX())
			fmt.Fprintf(w, "  NAXIS  = %d\n", v.NAXIS())
			fmt.Fprintf(w, "  Shape  = %v\n", v.Shape())
		case *fits.BinaryTableHDU:
			cols, err := v.Columns()
			if err != nil {
				fmt.Fprintf(w, "  Columns: error: %v\n", err)
				problems = append(problems, fmt.Sprintf("HDU %d columns: %v", i, err))
			} else {
				fmt.Fprintf(w, "  Rows    = %d\n", v.NumRows())
				fmt.Fprintf(w, "  Columns = %d\n", len(cols))
				for _, c := range cols {
					fmt.Fprintf(w, "    %2d %-10s TFORM=%-8s unit=%s\n", c.Index, c.Name, c.TForm, c.Unit)
				}
			}
		case *fits.ASCIITableHDU:
			fmt.Fprintf(w, "  Rows = %d\n", v.NumRows())
		}
		fmt.Fprintln(w)
		// HeaderOf cannot fail for a file Open accepted (see the doc); the
		// error is still reported rather than ignored.
		hdr, err := fits.HeaderOf(h)
		if err != nil {
			fmt.Fprintf(w, "  Header: error: %v\n\n", err)
			problems = append(problems, fmt.Sprintf("HDU %d header: %v", i, err))
			continue
		}
		fmt.Fprintln(w, "  Header:")
		cards := hdr.Cards()
		// Blank cards between the last keyword and END are reserved space,
		// not keywords.
		for len(cards) > 0 && cards[len(cards)-1].Key == "" && strings.TrimSpace(cards[len(cards)-1].Comment) == "" {
			cards = cards[:len(cards)-1]
		}
		for _, c := range cards {
			switch c.Value.(type) {
			case nil:
				fmt.Fprintf(w, "    %-8s %s\n", c.Key, c.Comment)
			default:
				fmt.Fprintf(w, "    %-8s = %-20v / %s\n", c.Key, c.Value, c.Comment)
			}
		}
		fmt.Fprintln(w)
	}
	if len(problems) > 0 {
		return fmt.Errorf("%s", strings.Join(problems, "; "))
	}
	return nil
}
