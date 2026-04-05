// fitsdump prints a human-readable summary of every HDU in a FITS file,
// including the full header keyword list.
//
// Usage:
//
//	fitsdump file.fits [file.fits ...]
package main

import (
	"fmt"
	"os"

	"github.com/dmarkham/fits"
)

func main() {
	if len(os.Args) < 2 {
		fmt.Fprintln(os.Stderr, "usage: fitsdump FILE [FILE ...]")
		os.Exit(2)
	}
	exit := 0
	for _, path := range os.Args[1:] {
		if err := dump(path); err != nil {
			fmt.Fprintf(os.Stderr, "fitsdump: %s: %v\n", path, err)
			exit = 1
		}
	}
	os.Exit(exit)
}

func dump(path string) error {
	f, err := fits.Open(path)
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Printf("=== %s ===\n", path)
	fmt.Printf("HDUs: %d\n\n", f.NumHDU())

	for i := range f.NumHDU() {
		h, err := f.HDU(i)
		if err != nil {
			fmt.Printf("HDU %d: error: %v\n", i, err)
			continue
		}
		fmt.Printf("--- HDU %d: %s ---\n", i, h.Type())
		switch v := h.(type) {
		case *fits.ImageHDU:
			fmt.Printf("  BITPIX = %d\n", v.BITPIX())
			fmt.Printf("  NAXIS  = %d\n", v.NAXIS())
			if v.NAXIS() > 0 {
				fmt.Printf("  Shape  = %v\n", v.Shape())
			}
		case *fits.BinaryTableHDU:
			cols, err := v.Columns()
			if err == nil {
				fmt.Printf("  Rows    = %d\n", v.NumRows())
				fmt.Printf("  Columns = %d\n", len(cols))
				for _, c := range cols {
					fmt.Printf("    %2d %-10s TFORM=%-8s unit=%s\n", c.Index, c.Name, c.TForm, c.Unit)
				}
			}
		case *fits.ASCIITableHDU:
			fmt.Printf("  Rows = %d\n", v.NumRows())
		}
		fmt.Println()
		fmt.Println("  Header:")
		hdr := h.Header()
		for _, c := range hdr.Cards() {
			switch c.Value.(type) {
			case nil:
				fmt.Printf("    %-8s %s\n", c.Key, c.Comment)
			default:
				fmt.Printf("    %-8s = %-20v / %s\n", c.Key, c.Value, c.Comment)
			}
		}
		fmt.Println()
	}
	return nil
}
