package main

import (
	"bytes"
	"errors"
	"fmt"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits"
	"github.com/dmarkham/fits/internal/testutil"
)

// TestDumpMatchesCfitsio runs fitsdump on every committed sample and checks
// what it prints against the cfitsio golden for that sample
// (testdata/cfref/golden): the HDU count and kinds, image BITPIX / NAXIS /
// shape (logical values for tile-compressed images), table row and column
// counts and each column's name, TFORM and unit, and every header keyword
// in order with its value and comment.
func TestDumpMatchesCfitsio(t *testing.T) {
	var samples []string
	for _, pat := range []string{"*.fits", "*.fit", "*.std"} {
		m, _ := filepath.Glob(testutil.TestdataPath(t, pat))
		samples = append(samples, m...)
	}
	sort.Strings(samples)
	if len(samples) < 43 {
		t.Fatalf("found %d samples, want at least 43", len(samples))
	}
	for _, path := range samples {
		name := filepath.Base(path)
		t.Run(name, func(t *testing.T) {
			g := testutil.LoadCfrefGolden(t, name)
			var out, errb bytes.Buffer
			if code := run([]string{path}, &out, &errb); code != 0 {
				t.Fatalf("exit %d, stderr %q", code, errb.String())
			}
			if errb.Len() != 0 {
				t.Fatalf("stderr: %q", errb.String())
			}
			checkDump(t, path, out.String(), g)
		})
	}
}

func checkDump(t *testing.T, path, out string, g *testutil.CfrefFile) {
	t.Helper()
	pre := fmt.Sprintf("=== %s ===\nHDUs: %d\n\n", path, g.NHDU)
	if !strings.HasPrefix(out, pre) {
		t.Fatalf("output starts %q, want %q", firstLines(out, 2), pre)
	}
	out = strings.TrimPrefix(out, pre)
	for i := range g.HDUs {
		h := &g.HDUs[i]
		sep := fmt.Sprintf("--- HDU %d: ", i)
		if !strings.HasPrefix(out, sep) {
			t.Fatalf("HDU %d: section starts %q", i, firstLines(out, 1))
		}
		end := strings.Index(out, "\n--- HDU ")
		var sec string
		if end < 0 {
			sec, out = out, ""
		} else {
			sec, out = out[:end+1], out[end+1:]
		}
		want := expectedSection(t, i, h)
		if got := strings.TrimRight(sec, "\n"); got != strings.TrimRight(want, "\n") {
			gl, wl := strings.Split(got, "\n"), strings.Split(want, "\n")
			for j := 0; j < max(len(gl), len(wl)); j++ {
				var a, b string
				if j < len(gl) {
					a = gl[j]
				}
				if j < len(wl) {
					b = wl[j]
				}
				if a != b {
					t.Errorf("HDU %d line %d:\n got  %q\n want %q", i, j+1, a, b)
					break
				}
			}
		}
	}
	if out != "" {
		t.Errorf("output after the last HDU: %q", firstLines(out, 3))
	}
}

func firstLines(s string, n int) string {
	l := strings.SplitN(s, "\n", n+1)
	return strings.Join(l[:min(n, len(l))], "\n")
}

// expectedSection builds the text fitsdump must print for one HDU, from
// the cfitsio golden alone.
func expectedSection(t *testing.T, i int, h *testutil.CfrefHDU) string {
	t.Helper()
	var b strings.Builder
	kind := map[string]string{"IMAGE": "image", "BINARY_TBL": "binary_table", "ASCII_TBL": "ascii_table"}[h.Type]
	fmt.Fprintf(&b, "--- HDU %d: %s ---\n", i, kind)
	switch {
	case h.Image != nil && h.Image.Compressed:
		cmp := cardStr(t, h, "ZCMPTYPE")
		if cmp == "RICE_ONE" {
			cmp = "RICE_1" // cfitsio's alias; fitsdump names the algorithm
		}
		fmt.Fprintf(&b, "  Compressed: %s\n", cmp)
		fmt.Fprintf(&b, "  BITPIX = %d (logical, from ZBITPIX)\n", h.Image.Bitpix)
		fmt.Fprintf(&b, "  NAXIS  = %d\n", h.Image.Naxis)
		fmt.Fprintf(&b, "  Shape  = %v\n", h.Image.Naxes)
	case h.Image != nil:
		fmt.Fprintf(&b, "  BITPIX = %d\n", h.Image.Bitpix)
		fmt.Fprintf(&b, "  NAXIS  = %d\n", h.Image.Naxis)
		if h.Image.Naxis > 0 {
			fmt.Fprintf(&b, "  Shape  = %v\n", h.Image.Naxes)
		}
	case h.Type == "BINARY_TBL":
		fmt.Fprintf(&b, "  Rows    = %d\n", h.Table.NRows)
		fmt.Fprintf(&b, "  Columns = %d\n", h.Table.NCols)
		for _, c := range h.Table.Columns {
			fmt.Fprintf(&b, "    %2d %-10s TFORM=%-8s unit=%s\n", c.N, latin(c.TType), strings.TrimRight(string(c.TForm), " "), latin(c.TUnit))
		}
	case h.Type == "ASCII_TBL":
		fmt.Fprintf(&b, "  Rows = %d\n", h.Table.NRows)
	}
	b.WriteString("\n  Header:\n")
	cards := h.Cards
	for j := 0; j < len(cards); j++ {
		c := cards[j]
		key := string(c.Key)
		if strings.HasPrefix(string(c.Raw), "HIERARCH ") {
			key = "HIERARCH " + key // cfitsio strips the prefix; fitsdump keeps it
		}
		// The header parser keeps a comment exactly as cfitsio reads it,
		// leading blanks included.
		comment := string(c.Comment)
		if c.Type == "" {
			// Commentary card, or a keyword with no value: printed as the
			// keyword followed by the text after column 8 (cfitsio's
			// comment field for commentary cards is that text).
			text := ""
			if raw := string(c.Raw); len(raw) > 8 {
				text = raw[8:]
			}
			if c.Class != 0 && strings.Contains(string(c.Raw), "=") && len(c.Raw) > 9 && c.Raw[8] == '=' {
				text = comment // "KEY     =  / comment": no value
			}
			fmt.Fprintf(&b, "    %-8s %s\n", key, text)
			continue
		}
		var v string
		switch c.Type {
		case "C":
			v = string(*c.Str) // CONTINUE already joined by cfitsio
			// The joined card carries the comments of all its pieces. Only
			// a piece whose string ends in '&' continues (cfitsio ffgkls);
			// a CONTINUE card after any other string is a card of its own.
			for j+1 < len(cards) && string(cards[j+1].Key) == "CONTINUE" && continues(string(cards[j].Raw)) {
				j++
				if cc := strings.TrimSpace(contComment(string(cards[j].Raw))); cc != "" {
					if comment == "" {
						comment = cc
					} else {
						comment += " " + cc
					}
				}
			}
			v = strings.TrimRight(v, " ")
		case "L":
			v = strconv.FormatBool(*c.Logical)
		case "I":
			if c.Int == nil {
				// an integer cfitsio cannot hold in a long long (BZERO =
				// 9223372036854775808 of a uint64 image): its digits
				v = strings.TrimSpace(string(c.Value))
				break
			}
			v = strconv.FormatInt(*c.Int, 10)
		case "F":
			v = fmt.Sprint(float64(*c.Float))
		case "X":
			re, im := parseComplex(t, string(c.Value))
			v = fmt.Sprintf("(%g,%g)", re, im)
		default:
			t.Fatalf("card %s: unknown cfitsio type %q", key, c.Type)
		}
		fmt.Fprintf(&b, "    %-8s = %-20s / %s\n", key, v, comment)
	}
	b.WriteString("\n")
	return b.String()
}

func latin(p *testutil.Latin1) string {
	if p == nil {
		return ""
	}
	return strings.TrimRight(string(*p), " ")
}

func cardStr(t *testing.T, h *testutil.CfrefHDU, key string) string {
	t.Helper()
	c := h.Card(key)
	if c == nil || c.Str == nil {
		t.Fatalf("golden HDU has no %s string card", key)
	}
	return string(*c.Str)
}

// continues reports whether the string value of a card ("KEY     = 'x&'"
// or "CONTINUE  'x&'") ends in '&', i.e. the next CONTINUE card joins it.
func continues(raw string) bool {
	i := strings.IndexByte(raw, '\'')
	if i < 0 {
		return false
	}
	var v strings.Builder
	for i++; i < len(raw); i++ {
		if raw[i] == '\'' {
			if i+1 < len(raw) && raw[i+1] == '\'' {
				v.WriteByte('\'')
				i++
				continue
			}
			break
		}
		v.WriteByte(raw[i])
	}
	return strings.HasSuffix(strings.TrimRight(v.String(), " "), "&")
}

// contComment is the comment of a CONTINUE card ("CONTINUE  'text' / c").
func contComment(raw string) string {
	inStr := false
	for i := 10; i < len(raw); i++ {
		switch raw[i] {
		case '\'':
			inStr = !inStr
		case '/':
			if !inStr {
				return raw[i+1:]
			}
		}
	}
	return ""
}

func parseComplex(t *testing.T, s string) (float64, float64) {
	t.Helper()
	s = strings.Trim(strings.TrimSpace(s), "()")
	p := strings.Split(s, ",")
	if len(p) != 2 {
		t.Fatalf("bad complex value %q", s)
	}
	re, err1 := strconv.ParseFloat(strings.Replace(strings.TrimSpace(p[0]), "D", "E", 1), 64)
	im, err2 := strconv.ParseFloat(strings.Replace(strings.TrimSpace(p[1]), "D", "E", 1), 64)
	if err1 != nil || err2 != nil {
		t.Fatalf("bad complex value %q", s)
	}
	return re, im
}

// Usage errors and unreadable files give the documented exit status and a
// message on stderr; one bad file does not stop the others.
func TestDumpExitStatus(t *testing.T) {
	var out, errb bytes.Buffer
	if code := run(nil, &out, &errb); code != 2 || !strings.Contains(errb.String(), "usage") {
		t.Fatalf("no args: exit %d, stderr %q", code, errb.String())
	}
	out.Reset()
	errb.Reset()
	good := testutil.TestdataPath(t, "vari.fits")
	missing := filepath.Join(t.TempDir(), "missing.fits")
	if code := run([]string{missing, good}, &out, &errb); code != 1 {
		t.Fatalf("missing file: exit %d", code)
	}
	if !strings.Contains(errb.String(), "missing.fits") {
		t.Fatalf("stderr %q does not name the bad file", errb.String())
	}
	if !strings.Contains(out.String(), "=== "+good+" ===") {
		t.Fatal("the good file after the bad one was not dumped")
	}
}

// A binary table whose columns cannot be decoded (TFORM2 of iter_a.fit
// patched to the invalid code 'W') is reported, not skipped: the HDU's
// section says "Columns: error" with the reason, its header is still
// printed, and the exit status is 1 with the file named on stderr.
func TestDumpReportsColumnError(t *testing.T) {
	src := testutil.TestdataPath(t, "iter_a.fit")
	testutil.RequireFile(t, src)
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatal(err)
	}
	old := []byte("TFORM2  = 'E       '")
	i := bytes.Index(b, old)
	if i < 0 {
		t.Fatal("iter_a.fit has no TFORM2 = 'E' card")
	}
	copy(b[i:], "TFORM2  = 'W       '")
	path := filepath.Join(t.TempDir(), "badtform.fits")
	if err := os.WriteFile(path, b, 0o644); err != nil {
		t.Fatal(err)
	}
	var out, errb bytes.Buffer
	code := run([]string{path}, &out, &errb)
	s := out.String()
	sec := s[strings.Index(s, "--- HDU 1:"):]
	if !strings.Contains(sec, "  Columns: error: ") || !strings.Contains(sec, "W") {
		t.Fatalf("HDU 1 section does not report the column error:\n%s", sec)
	}
	if !strings.Contains(sec, "TFORM2   = W") {
		t.Fatalf("HDU 1 header not printed after the column error:\n%s", sec)
	}
	if code != 1 || !strings.Contains(errb.String(), "badtform.fits") {
		t.Fatalf("exit %d, stderr %q; want 1 and the file named", code, errb.String())
	}
}

// fitsFile packs header cards (80 columns each, END appended, blank fill
// to 2880 bytes) and zero data blocks into one FITS byte string per HDU.
func fitsFile(hdus ...[]string) []byte {
	var b bytes.Buffer
	for _, cards := range hdus {
		for _, c := range append(cards, "END") {
			b.WriteString(fmt.Sprintf("%-80s", c))
		}
		for b.Len()%2880 != 0 {
			b.WriteByte(' ')
		}
	}
	return b.Bytes()
}

// imageExt is a 2x1 BITPIX 8 IMAGE extension (one data block, zeros).
func imageExt() []byte {
	b := fitsFile([]string{
		"XTENSION= 'IMAGE   '", "BITPIX  =                    8", "NAXIS   =                    2",
		"NAXIS1  =                    2", "NAXIS2  =                    1", "PCOUNT  =                    0",
		"GCOUNT  =                    1", "EXTNAME = 'AFTER   '",
	})
	return append(b, make([]byte, 2880)...)
}

// A CONTINUE card after a string that does not end in '&' (a "lone"
// CONTINUE) is a header this package opens and parses like cfitsio: a
// commentary-like card of its own. fitsdump prints it as cfitsio reads it
// and exits 0, on a primary image and on a BINTABLE (whose Columns() call
// comes before the header is printed) and a tile-compressed image (whose
// CompressionType() does). Truth: the live cfref dump of each file.
func TestDumpLoneContinueMatchesCfitsio(t *testing.T) {
	lone := []string{"LONGSTR = 'abc'", "CONTINUE  'def'"}
	bintable := append([]string{
		"XTENSION= 'BINTABLE'", "BITPIX  =                    8", "NAXIS   =                    2",
		"NAXIS1  =                    4", "NAXIS2  =                    2", "PCOUNT  =                    0",
		"GCOUNT  =                    1", "TFIELDS =                    1", "TTYPE1  = 'J       '", "TFORM1  = '1J      '",
	}, lone...)
	primary := append([]string{
		"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0",
		"EXTEND  =                    T"}, lone...)
	tbl := append(fitsFile(bintable), make([]byte, 2880)...)
	comp, err := os.ReadFile(testutil.TestdataPath(t, "comp_rice_i16.fits"))
	if err != nil {
		t.Fatal(err)
	}
	// Put the lone CONTINUE into the compressed HDU's header after its last
	// keyword: it takes the END card's place and END moves one card on.
	compLone := append([]byte(nil), comp...)
	placed := false
	for off := 2880; off+160 <= len(compLone); off += 80 {
		if string(compLone[off:off+8]) == "END     " {
			if off%2880 == 2880-80 {
				break // END is the block's last card
			}
			copy(compLone[off:off+80], fmt.Sprintf("%-80s", "CONTINUE  'def'"))
			copy(compLone[off+80:off+160], fmt.Sprintf("%-80s", "END"))
			placed = true
			break
		}
	}
	if !placed {
		t.Fatal("no room in comp_rice_i16.fits HDU 1 for the lone CONTINUE")
	}
	for _, c := range []struct {
		name string
		b    []byte
	}{
		{"primary", append(fitsFile(primary), imageExt()...)},
		{"bintable", append(fitsFile([]string{"SIMPLE  =                    T", "BITPIX  =                    8", "NAXIS   =                    0", "EXTEND  =                    T"}), tbl...)},
		{"compressed", compLone},
	} {
		t.Run(c.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "lonecont.fits")
			if err := os.WriteFile(path, c.b, 0o644); err != nil {
				t.Fatal(err)
			}
			var out, errb bytes.Buffer
			if code := run([]string{path}, &out, &errb); code != 0 {
				t.Fatalf("exit %d, stderr %q; want 0\n%s", code, errb.String(), firstLines(out.String(), 60))
			}
			if !strings.Contains(out.String(), "CONTINUE") {
				t.Fatalf("lone CONTINUE card not printed:\n%s", out.String())
			}
			checkDump(t, path, out.String(), testutil.CfrefDump(t, path))
		})
	}
}

// An HDU that cannot be decoded (File.HDU fails; today only a random-groups
// primary does, and the committed samples have none, so the failure is
// injected for HDU 1 of testprog.std) is reported in place as "HDU 1:
// error", the HDUs after it are still printed, and the exit status is 1
// with the file and the HDU named on stderr.
func TestDumpReportsHDUError(t *testing.T) {
	path := testutil.TestdataPath(t, "testprog.std")
	testutil.RequireFile(t, path)
	injected := errors.New("injected HDU failure")
	hduOf = func(f *fits.File, i int) (fits.HDU, error) {
		if i == 1 {
			return nil, injected
		}
		return f.HDU(i)
	}
	defer func() { hduOf = (*fits.File).HDU }()
	var out, errb bytes.Buffer
	code := run([]string{path}, &out, &errb)
	s := out.String()
	if !strings.Contains(s, "HDU 1: error: injected HDU failure\n") {
		t.Fatalf("HDU 1 error not reported in place:\n%s", firstLines(s, 40))
	}
	g := testutil.LoadCfrefGolden(t, "testprog.std")
	for i := 2; i < g.NHDU; i++ {
		if !strings.Contains(s, fmt.Sprintf("--- HDU %d: ", i)) {
			t.Fatalf("HDU %d after the bad one was not printed", i)
		}
	}
	if code != 1 || !strings.Contains(errb.String(), "testprog.std") || !strings.Contains(errb.String(), "HDU 1: injected HDU failure") {
		t.Fatalf("exit %d, stderr %q; want 1 with the file and HDU 1 named", code, errb.String())
	}
}

// A compressed HDU whose Z* keywords cannot be read (ZNAXIS = -1; cfitsio
// refuses it with BAD_NAXIS) is reported in place, not printed as a zero
// shape (and does not panic, as it did before the ZNAXIS range check);
// the exit status is 1.
func TestDumpReportsCompressedMetadataError(t *testing.T) {
	path := filepath.Join("..", "..", "testdata", "comp_bad", "u16_rice_znaxism1.fits")
	testutil.RequireFile(t, path)
	var out, errb bytes.Buffer
	code := run([]string{path}, &out, &errb)
	s := out.String()
	if !strings.Contains(s, "  Compressed image: error: ") || !strings.Contains(s, "ZNAXIS -1") {
		t.Fatalf("compressed metadata error not reported in place:\n%s", firstLines(s, 40))
	}
	if code != 1 || !strings.Contains(errb.String(), "HDU 1 compressed image") {
		t.Fatalf("exit %d, stderr %q; want 1 with HDU 1 named", code, errb.String())
	}
}

// An extension of a type the library does not decode (XTENSION =
// 'FOREIGN', testdata/images/framing/xtension_foreign.fits) is reported in
// place with its XTENSION value, not an internal kind number; the IMAGE
// after it is still dumped, and the exit status is 1.
func TestDumpReportsUnsupportedExtension(t *testing.T) {
	path := filepath.Join("..", "..", "testdata", "images", "framing", "xtension_foreign.fits")
	testutil.RequireFile(t, path)
	var out, errb bytes.Buffer
	code := run([]string{path}, &out, &errb)
	s := out.String()
	if !strings.Contains(s, "HDU 1: error: fits: HDU 1: unsupported extension type XTENSION = 'FOREIGN'\n") {
		t.Fatalf("HDU 1 not reported as an unsupported FOREIGN extension:\n%s", firstLines(s, 40))
	}
	if strings.Contains(s, "kind") {
		t.Errorf("dump leaks an internal HDU kind:\n%s", firstLines(s, 40))
	}
	if !strings.Contains(s, "--- HDU 2: ") {
		t.Fatalf("the IMAGE after the FOREIGN HDU was not dumped:\n%s", firstLines(s, 60))
	}
	if code != 1 || !strings.Contains(errb.String(), "HDU 1: fits: HDU 1: unsupported extension type") {
		t.Fatalf("exit %d, stderr %q; want 1 with HDU 1 named", code, errb.String())
	}
}
