package testutil

import (
	"bytes"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"io/fs"
	"os"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// Large text goldens are committed gzipped, as <name>.gz written with
// `gzip -9n` (no name, no timestamp), to keep the repository inside its
// size budget (TestRepositorySizeBudget). The check stays exact: the
// decompressed bytes are what the harness wrote, and every reproduce
// test and `make verify` compares the harness output with them byte for
// byte. (The compressed bytes depend on the zlib/gzip build, so they are
// not compared.)

// ReadGolden reads the committed golden path: path itself if it exists,
// otherwise path+".gz" decompressed. The error of a missing golden names
// both spellings.
func ReadGolden(path string) ([]byte, error) {
	b, err := os.ReadFile(path)
	if err == nil {
		return b, nil
	}
	if !errors.Is(err, fs.ErrNotExist) {
		return nil, err
	}
	z, zerr := os.ReadFile(path + ".gz")
	if zerr != nil {
		if errors.Is(zerr, fs.ErrNotExist) {
			return nil, fmt.Errorf("%w (nor %s.gz)", err, filepath.Base(path))
		}
		return nil, zerr
	}
	return Gunzip(z)
}

// Gunzip decompresses gzip data (the CRC-32 and length trailer are
// checked by compress/gzip).
func Gunzip(z []byte) ([]byte, error) {
	r, err := gzip.NewReader(bytes.NewReader(z))
	if err != nil {
		return nil, err
	}
	b, err := io.ReadAll(r)
	if err != nil {
		return nil, err
	}
	return b, r.Close()
}

// GoldenExists reports whether path or path+".gz" exists.
func GoldenExists(path string) bool {
	if _, err := os.Stat(path); err == nil {
		return true
	}
	_, err := os.Stat(path + ".gz")
	return err == nil
}

// GlobGoldens returns the goldens matching pattern (which must end in
// the uncompressed extension, e.g. "dir/*.json"), committed plain or
// gzipped, by their uncompressed name, sorted.
func GlobGoldens(pattern string) []string {
	seen := map[string]bool{}
	plain, _ := filepath.Glob(pattern)
	gz, _ := filepath.Glob(pattern + ".gz")
	for _, p := range plain {
		seen[p] = true
	}
	for _, p := range gz {
		seen[strings.TrimSuffix(p, ".gz")] = true
	}
	out := make([]string, 0, len(seen))
	for p := range seen {
		out = append(out, p)
	}
	sort.Strings(out)
	return out
}

// RequireGolden fails the test (never skips) unless the committed golden
// path or path+".gz" exists.
func RequireGolden(t testing.TB, path string) {
	t.Helper()
	if !GoldenExists(path) {
		t.Fatalf("committed golden missing (this is a failure, not a skip): %s (or %s.gz)", path, filepath.Base(path))
	}
}

// MustReadGolden is ReadGolden that fails the test on an error.
func MustReadGolden(t testing.TB, path string) []byte {
	t.Helper()
	b, err := ReadGolden(path)
	if err != nil {
		t.Fatalf("committed golden missing or unreadable (this is a failure, not a skip): %v", err)
	}
	return b
}
