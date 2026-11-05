//go:build !unix

package fits

import "os"

// On platforms without Unix inodes, links and owners cannot be checked or
// carried over; editing is refused there (editSupported).
func checkReplaceable(string, os.FileInfo) error { return editSupported() }

func carryOver(string, os.FileInfo, *os.File, string) error { return editSupported() }
