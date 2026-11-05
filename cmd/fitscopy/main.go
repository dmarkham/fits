// fitscopy reads an input FITS file through the fits library, walks every
// HDU, and writes them to an output file via the library's public CopyTo
// API, then checks that the output is byte-for-byte identical to the input
// (it is not when the input has bytes after its last HDU or its last data
// unit lacks its final padding; see (*fits.File).CopyTo).
// It exercises HDU framing only (finding every header and data unit and
// passing its bytes through); it does not decode or check any value.
//
// Exit status: 0 when the copy is byte-identical, 1 on an error (including
// an input the library cannot open), 2 on bad usage or when INPUT and
// OUTPUT are the same file (which would be truncated before it is read),
// 3 when the copy differs from the input.
//
// Usage:
//
//	fitscopy input.fits output.fits
package main

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"

	"github.com/dmarkham/fits"
)

func main() {
	os.Exit(run(os.Args[1:], os.Stderr))
}

// run is main without the process exit, so tests can drive it.
func run(args []string, stderr io.Writer) int {
	if len(args) != 2 {
		fmt.Fprintln(stderr, "usage: fitscopy INPUT OUTPUT")
		return 2
	}
	in, out := args[0], args[1]
	if same, err := sameFile(in, out); err != nil {
		fmt.Fprintf(stderr, "fitscopy: %v\n", err)
		return 1
	} else if same {
		fmt.Fprintf(stderr, "fitscopy: INPUT and OUTPUT are the same file (%s); refusing to overwrite it\n", in)
		return 2
	}
	if err := copyFn(in, out); err != nil {
		fmt.Fprintf(stderr, "fitscopy: %v\n", err)
		return 1
	}
	eq, err := filesEqual(in, out)
	if err != nil {
		fmt.Fprintf(stderr, "fitscopy: diff: %v\n", err)
		return 1
	}
	if !eq {
		fmt.Fprintln(stderr, "fitscopy: WARNING: output differs from input byte-for-byte")
		return 3
	}
	return 0
}

// sameFile reports whether out names the existing file in (by inode, so
// different spellings and links are caught). A missing out is not the same.
func sameFile(in, out string) (bool, error) {
	si, err := os.Stat(in)
	if err != nil {
		return false, err
	}
	so, err := os.Stat(out)
	if errors.Is(err, os.ErrNotExist) {
		return false, nil
	}
	if err != nil {
		return false, err
	}
	return os.SameFile(si, so), nil
}

// copyFn is the copy step; a test replaces it to check that a copy that
// differs from its input is reported (no committed input makes CopyTo
// differ, and an unreadable one fails before the comparison).
var copyFn = copyFile

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
	ba, err := os.ReadFile(a)
	if err != nil {
		return false, err
	}
	bb, err := os.ReadFile(b)
	if err != nil {
		return false, err
	}
	return bytes.Equal(ba, bb), nil
}
