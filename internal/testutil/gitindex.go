package testutil

import (
	"bufio"
	"bytes"
	"encoding/binary"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"sort"
	"strings"
)

// TrackedFiles lists the files git tracks in the checkout at root (slash
// separated, relative to root), read from the git index in pure Go so
// tier 1 needs no git program. ok is false when root is not a git
// checkout (a module download or an exported tree); the caller then
// walks the tree instead (WalkFiles).
//
// The index format is git's documented one (Documentation/gitformat-
// index.txt): versions 2, 3 and 4, SHA-1 or SHA-256 object names. Only
// stage-0 entries are listed (a merge in progress lists each path once).
func TrackedFiles(root string) (files []string, ok bool, err error) {
	gitDir, err := findGitDir(root)
	if err != nil || gitDir == "" {
		return nil, false, err
	}
	b, err := os.ReadFile(filepath.Join(gitDir, "index"))
	if err != nil {
		return nil, false, fmt.Errorf("git index: %w", err)
	}
	hashLen := 20
	if cfg, err := os.ReadFile(filepath.Join(commonDir(gitDir), "config")); err == nil && bytes.Contains(cfg, []byte("objectformat = sha256")) {
		hashLen = 32
	}
	files, err = parseGitIndex(b, hashLen)
	return files, err == nil, err
}

// findGitDir returns root/.git, or the directory a worktree's .git file
// points at ("gitdir: <path>"); "" when there is neither.
func findGitDir(root string) (string, error) {
	p := filepath.Join(root, ".git")
	st, err := os.Stat(p)
	if errors.Is(err, fs.ErrNotExist) {
		return "", nil
	}
	if err != nil {
		return "", err
	}
	if st.IsDir() {
		return p, nil
	}
	b, err := os.ReadFile(p)
	if err != nil {
		return "", err
	}
	line, _, _ := strings.Cut(string(b), "\n")
	dir, found := strings.CutPrefix(strings.TrimSpace(line), "gitdir:")
	if !found {
		return "", fmt.Errorf("%s: not a gitdir file", p)
	}
	dir = strings.TrimSpace(dir)
	if !filepath.IsAbs(dir) {
		dir = filepath.Join(root, dir)
	}
	return dir, nil
}

// commonDir is the repository's shared git directory (a worktree's
// "commondir" file names it relative to the worktree's git directory).
func commonDir(gitDir string) string {
	b, err := os.ReadFile(filepath.Join(gitDir, "commondir"))
	if err != nil {
		return gitDir
	}
	d := strings.TrimSpace(string(b))
	if !filepath.IsAbs(d) {
		d = filepath.Join(gitDir, d)
	}
	return d
}

func parseGitIndex(b []byte, hashLen int) ([]string, error) {
	if len(b) < 12 || string(b[:4]) != "DIRC" {
		return nil, errors.New("git index: bad signature")
	}
	version := binary.BigEndian.Uint32(b[4:8])
	if version < 2 || version > 4 {
		return nil, fmt.Errorf("git index: version %d not supported", version)
	}
	n := int(binary.BigEndian.Uint32(b[8:12]))
	off := 12
	var prev string
	files := make([]string, 0, n)
	for i := 0; i < n; i++ {
		start := off
		// ctime, mtime (8 bytes each), dev, ino, mode, uid, gid, size.
		off += 40 + hashLen
		if off+2 > len(b) {
			return nil, errors.New("git index: truncated entry")
		}
		flags := binary.BigEndian.Uint16(b[off:])
		off += 2
		if version >= 3 && flags&0x4000 != 0 {
			off += 2 // extended flags
		}
		var name string
		if version == 4 {
			// Prefix compression: a varint of bytes to drop from the
			// previous name, then the rest NUL-terminated.
			strip, k := gitVarint(b[off:])
			if k == 0 || int(strip) > len(prev) {
				return nil, errors.New("git index: bad v4 name prefix")
			}
			off += k
			end := bytes.IndexByte(b[off:], 0)
			if end < 0 {
				return nil, errors.New("git index: unterminated name")
			}
			name = prev[:len(prev)-int(strip)] + string(b[off:off+end])
			off += end + 1
		} else {
			end := bytes.IndexByte(b[off:], 0)
			if end < 0 {
				return nil, errors.New("git index: unterminated name")
			}
			name = string(b[off : off+end])
			// The entry is padded with 1..8 NULs to a multiple of 8.
			off = start + (off+end-start+8)&^7
		}
		prev = name
		if stage := (flags >> 12) & 3; stage == 0 {
			files = append(files, name)
		}
	}
	return files, nil
}

// gitVarint decodes git's offset varint (index v4 names).
func gitVarint(b []byte) (uint64, int) {
	if len(b) == 0 {
		return 0, 0
	}
	v := uint64(b[0] & 0x7f)
	k := 1
	for c := b[0]; c&0x80 != 0; k++ {
		if k >= len(b) {
			return 0, 0
		}
		c = b[k]
		v = (v+1)<<7 | uint64(c&0x7f)
	}
	return v, k
}

// WalkFiles lists every regular file under root (slash separated,
// relative), skipping the directories that are never part of the
// repository (.git, venv, node_modules) and anything matched by a
// top-level .gitignore pattern of the simple forms this repository uses
// ("/path", "dir/", "*.ext", "name").
func WalkFiles(root string) ([]string, error) {
	ignore := readGitignore(filepath.Join(root, ".gitignore"))
	var out []string
	err := filepath.WalkDir(root, func(p string, d fs.DirEntry, err error) error {
		if err != nil {
			return err
		}
		rel, _ := filepath.Rel(root, p)
		rel = filepath.ToSlash(rel)
		if rel == "." {
			return nil
		}
		if d.IsDir() {
			switch d.Name() {
			case ".git", "venv", "node_modules":
				return filepath.SkipDir
			}
			if ignore(rel, true) {
				return filepath.SkipDir
			}
			return nil
		}
		if d.Type().IsRegular() && !ignore(rel, false) {
			out = append(out, rel)
		}
		return nil
	})
	sort.Strings(out)
	return out, err
}

func readGitignore(path string) func(rel string, dir bool) bool {
	var pats []string
	if f, err := os.Open(path); err == nil {
		sc := bufio.NewScanner(f)
		for sc.Scan() {
			l := strings.TrimSpace(sc.Text())
			if l != "" && !strings.HasPrefix(l, "#") && !strings.HasPrefix(l, "!") {
				pats = append(pats, l)
			}
		}
		f.Close()
	}
	return func(rel string, dir bool) bool {
		base := rel[strings.LastIndexByte(rel, '/')+1:]
		for _, p := range pats {
			dirOnly := strings.HasSuffix(p, "/")
			p = strings.TrimSuffix(p, "/")
			if dirOnly && !dir {
				continue
			}
			var ok bool
			if strings.HasPrefix(p, "/") || strings.Contains(p, "/") {
				ok, _ = filepath.Match(strings.TrimPrefix(p, "/"), rel)
			} else {
				ok, _ = filepath.Match(p, base)
			}
			if ok {
				return true
			}
		}
		return false
	}
}
