// fitscopy reads an input FITS file through the fits library, walks every
// HDU, and writes them byte-for-byte to an output file via the library's
// public CopyTo API. When the library is working correctly and the input
// is well-formed, the output is byte-for-byte identical to the input.
//
// This is the primary regression tool for the byte-exact round-trip
// requirement in plan §"Done criteria for v1".
//
// Usage:
//
//	fitscopy input.fits output.fits
package main

import (
	"bytes"
	"fmt"
	"io"
	"os"

	"github.com/dmarkham/fits"
)

func main() {
	if len(os.Args) != 3 {
		fmt.Fprintln(os.Stderr, "usage: fitscopy INPUT OUTPUT")
		os.Exit(2)
	}
	in, out := os.Args[1], os.Args[2]
	if err := copyFile(in, out); err != nil {
		fmt.Fprintf(os.Stderr, "fitscopy: %v\n", err)
		os.Exit(1)
	}
	if in != out {
		if eq, err := filesEqual(in, out); err != nil {
			fmt.Fprintf(os.Stderr, "fitscopy: diff: %v\n", err)
			os.Exit(1)
		} else if !eq {
			fmt.Fprintln(os.Stderr, "fitscopy: WARNING: output differs from input byte-for-byte")
			os.Exit(3)
		}
	}
}

func copyFile(src, dst string) error {
	f, err := fits.Open(src)
	if err != nil {
		return err
	}
	defer f.Close()
	out, err := os.Create(dst)
	if err != nil {
		return err
	}
	if err := f.CopyTo(out); err != nil {
		out.Close()
		return err
	}
	return out.Close()
}

func filesEqual(a, b string) (bool, error) {
	ba, err := readAll(a)
	if err != nil {
		return false, err
	}
	bb, err := readAll(b)
	if err != nil {
		return false, err
	}
	return bytes.Equal(ba, bb), nil
}

func readAll(path string) ([]byte, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return io.ReadAll(f)
}
