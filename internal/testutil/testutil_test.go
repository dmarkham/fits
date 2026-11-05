package testutil

import (
	"bytes"
	"crypto/sha256"
	"encoding/hex"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"runtime"
	"sort"
	"strings"
	"testing"
)

// recorder stands in for *testing.T so the gating helpers can be observed:
// Skipf/Fatalf record and stop the goroutine, like the real ones.
type recorder struct {
	testing.TB
	skipped, failed bool
	msg             string
}

func (r *recorder) Helper()             {}
func (r *recorder) Logf(string, ...any) {}
func (r *recorder) Skipf(f string, a ...any) {
	r.skipped, r.msg = true, fmt.Sprintf(f, a...)
	runtime.Goexit()
}
func (r *recorder) Fatalf(f string, a ...any) {
	r.failed, r.msg = true, fmt.Sprintf(f, a...)
	runtime.Goexit()
}

func observe(t *testing.T, fn func(testing.TB)) *recorder {
	r := &recorder{TB: t}
	done := make(chan struct{})
	go func() { defer close(done); fn(r) }()
	<-done
	return r
}

func TestDecodeCfrefValues(t *testing.T) {
	js := `{"file":"x.fits","cfitsio":"4.7.0","nhdu":1,"hdus":[{"index":0,"type":"IMAGE",
	"header_start":0,"data_start":2880,"data_end":5760,"checksum":{"datasum":1,"hdusum":-1},
	"nkeys":1,"cards":[{"n":1,"raw":"COMMENT \u00e9\u00ff","key":"COMMENT","class":10,"value":"","comment":""}],
	"image":{"bitpix":64,"equivtype":80,"compressed":false,"naxis":1,"naxes":[4],"npix":4,"bscale":1,"bzero":9223372036854775808,
	"phys":[null,"NaN","Inf","-Inf"],"u64":[18446744073709551615,null,0,1],"raw":[-9223372036854775808,null,0,9223372036854775807]},
	"errors":[]}]}`
	cf, err := DecodeCfref([]byte(js))
	if err != nil {
		t.Fatal(err)
	}
	h := cf.HDUs[0]
	if h.Checksum.Datasum != 1 || h.Checksum.HDUSum != -1 || h.DataEnd != 5760 {
		t.Fatalf("hdu fields: %+v", h)
	}
	if got := []byte(h.Cards[0].Raw); !bytes.Equal(got, []byte{'C', 'O', 'M', 'M', 'E', 'N', 'T', ' ', 0xE9, 0xFF}) {
		t.Fatalf("Latin1 bytes: % x", got)
	}
	im := h.Image
	phys, u64, raw := im.Phys.Values(t), im.U64.Values(t), im.Raw.Values(t)
	if phys[0].Valid || !math.IsNaN(phys[1].V) || !phys[1].Valid || !math.IsInf(phys[2].V, 1) || !math.IsInf(phys[3].V, -1) {
		t.Fatalf("phys: %+v", phys)
	}
	if u64[0].V != math.MaxUint64 || !u64[0].Valid || u64[1].Valid {
		t.Fatalf("u64: %+v", u64)
	}
	if raw[0].V != math.MinInt64 || raw[3].V != math.MaxInt64 || raw[1].Valid {
		t.Fatalf("raw: %+v", raw)
	}
	if float64(im.Bzero) != 9223372036854775808 {
		t.Fatalf("bzero %v", im.Bzero)
	}
	if _, err := DecodeCfref([]byte(`{"file":"x","cfitsio":"4.7.0","nhdu":0,"hdus":[],"surprise":1}`)); err == nil {
		t.Fatal("unknown field accepted; the decoder must be strict")
	}
	if _, err := DecodeCfref([]byte(`{"file":"x","cfitsio":"4.7.0","nhdu":2,"hdus":[]}`)); err == nil {
		t.Fatal("nhdu mismatch accepted")
	}
}

// samples lists every FITS file committed in testdata/.
func samples(t *testing.T) []string {
	var out []string
	for _, pat := range []string{"*.fits", "*.fit", "*.std"} {
		m, err := filepath.Glob(TestdataPath(t, pat))
		if err != nil {
			t.Fatal(err)
		}
		out = append(out, m...)
	}
	sort.Strings(out)
	if len(out) == 0 {
		t.Fatal("no FITS samples found in testdata/")
	}
	return out
}

// Every FITS sample in testdata/ must have a committed cfitsio golden and
// every golden must correspond to a sample: a sample without truth is a
// claim without proof.
func TestGoldensCoverEveryTestdataSample(t *testing.T) {
	have := map[string]bool{}
	for _, s := range samples(t) {
		name := filepath.Base(s)
		have[name] = true
		cf := LoadCfrefGolden(t, name)
		if cf.File != name {
			t.Errorf("%s: golden file field %q", name, cf.File)
		}
		if cf.Cfitsio != "4.7.0" {
			t.Errorf("%s: golden made with cfitsio %s, PROVENANCE.md pins 4.7.0", name, cf.Cfitsio)
		}
		if cf.NHDU == 0 {
			t.Errorf("%s: golden has no HDUs", name)
		}
	}
	gs, _ := filepath.Glob(filepath.Join(CfrefDir(t), "golden", "*.json"))
	for _, g := range gs {
		if n := strings.TrimSuffix(filepath.Base(g), ".json"); !have[n] {
			t.Errorf("golden %s has no sample testdata/%s", filepath.Base(g), n)
		}
	}
}

func TestCfrefDumpRejectsNonFITS(t *testing.T) {
	p := filepath.Join(t.TempDir(), "junk.fits")
	if err := os.WriteFile(p, bytes.Repeat([]byte("x"), 2880), 0o644); err != nil {
		t.Fatal(err)
	}
	if _, err := CfrefDumpErr(t, p); err == nil {
		t.Fatal("cfitsio accepted a non-FITS file")
	}
}

func TestParseFitsverify(t *testing.T) {
	r := parseFitsverify(t, "x", "...\n**** Verification found 3 warning(s) and 2 error(s). ****\n")
	if r.Warnings != 3 || r.Errors != 2 || r.OK() {
		t.Fatalf("%+v", r)
	}
	rec := observe(t, func(tb testing.TB) { parseFitsverify(tb, "x", "garbage") })
	if !rec.failed {
		t.Fatal("unparseable fitsverify output must fail")
	}
}

// fitsverify's own verdicts on two cfitsio-distributed samples (recorded
// in tier 1, live in tier 2).
func TestFitsverifySamples(t *testing.T) {
	r := Fitsverify(t, TestdataPath(t, "testprog.std"))
	if r.Errors != 0 || r.Warnings != 1 {
		t.Fatalf("fitsverify testprog.std: %d warnings %d errors, want 1/0:\n%s", r.Warnings, r.Errors, r.Report)
	}
	r = Fitsverify(t, TestdataPath(t, "comp_raw_rice_i16.fits"))
	if !r.OK() {
		t.Fatalf("fitsverify comp_raw_rice_i16.fits not clean:\n%s", r.Report)
	}
}

// selfCheck verifies a digest's internal consistency: the deterministic
// sample (cfref.c jdigest) covers the first and last elements and at
// least 64 indices, and each sampled element's null flag agrees with the
// null runs.
func (a Values[T]) selfCheck() error {
	if a.dig == nil {
		return nil
	}
	n := int(a.dig.N)
	if n <= 1024 {
		return fmt.Errorf("digest of %d values (DIGEST_MIN is 1024)", n)
	}
	s := a.Sampled()
	if _, ok := s[0]; !ok {
		return fmt.Errorf("sample lacks index 0")
	}
	if _, ok := s[n-1]; !ok {
		return fmt.Errorf("sample lacks index %d", n-1)
	}
	if len(s) < 64 {
		return fmt.Errorf("only %d samples", len(s))
	}
	for k, v := range s {
		if isNull(v) != a.IsNull(k) {
			return fmt.Errorf("sample %d null=%v, null runs say %v", k, isNull(v), a.IsNull(k))
		}
	}
	return nil
}

// Every digest in every committed cfref-format golden is self-consistent,
// and the large arrays really are digests (the size budget depends on it).
func TestCfrefDigestsConsistent(t *testing.T) {
	var files []string
	for _, pat := range []string{"cfref/golden/*.json", "images/golden/*.fits.json", "tables/golden/*.json", "comp_bad/golden/*.json"} {
		files = append(files, GlobGoldens(TestdataPath(t, pat))...)
	}
	files = append(files, GlobGoldens(filepath.Join(RepoRoot(t), "header", "testdata", "golden", "*.fits.json"))...)
	if len(files) < 350 { // 361 after the compressed read corpus deduped its funpack twins
		t.Fatalf("only %d cfref-format goldens found", len(files))
	}
	digests := 0
	for _, p := range files {
		b, err := ReadGolden(p)
		if err != nil {
			t.Fatal(err)
		}
		cf, err := DecodeCfref(b)
		if err != nil {
			t.Fatalf("%s: %v", p, err)
		}
		for _, h := range cf.HDUs {
			var errs []error
			check := func(digested bool, n int, err error) {
				if digested {
					digests++
				} else if n > 1024 {
					errs = append(errs, fmt.Errorf("explicit array of %d values", n))
				}
				if err != nil {
					errs = append(errs, err)
				}
			}
			if im := h.Image; im != nil {
				check(im.Phys.Digested(), im.Phys.Len(), im.Phys.selfCheck())
				check(im.I64.Digested(), im.I64.Len(), im.I64.selfCheck())
				check(im.U64.Digested(), im.U64.Len(), im.U64.selfCheck())
				check(im.Raw.Digested(), im.Raw.Len(), im.Raw.selfCheck())
				check(im.PhysF32.Digested(), im.PhysF32.Len(), im.PhysF32.selfCheck())
				if im.Phys.Digested() && im.Bitpix == -32 && !im.PhysF32.Present() {
					errs = append(errs, fmt.Errorf("BITPIX -32 digest without phys_f32"))
				}
			}
			if tb := h.Table; tb != nil {
				for i := range tb.Columns {
					c := &tb.Columns[i]
					if !c.Digested() {
						if c.Typecode > 0 && c.Flat().Len() > 1024 {
							errs = append(errs, fmt.Errorf("column %d: explicit rows with %d values", c.N, c.Flat().Len()))
						}
						continue
					}
					f := c.Flat()
					check(f.Str.Digested(), 0, f.Str.selfCheck())
					check(f.Log.Digested(), 0, f.Log.selfCheck())
					check(f.Bits.Digested(), 0, f.Bits.selfCheck())
					check(f.Cplx.Digested(), 0, f.Cplx.selfCheck())
					check(f.Phys.Digested(), 0, f.Phys.selfCheck())
					check(f.I64.Digested(), 0, f.I64.selfCheck())
					check(f.U64.Digested(), 0, f.U64.selfCheck())
					check(f.Raw.Digested(), 0, f.Raw.selfCheck())
				}
			}
			for _, e := range errs {
				t.Errorf("%s HDU %d: %v", filepath.Base(p), h.Index, e)
			}
		}
	}
	if digests < 50 {
		t.Errorf("only %d digested arrays in all goldens", digests)
	}
	t.Logf("%d goldens, %d digested arrays", len(files), digests)
}

// The canonical element encoding behind every digest is the one cfref.c
// specifies (flag byte, big-endian bits, one NaN, -0 kept; strings
// length-prefixed). That the C harness's own SHA-256 over these bytes
// equals Go's crypto/sha256 is shown by every digest the suite resolves
// (e.g. TestCfrefDigestCatchesOneChange on real M13).
func TestDigestCanonicalEncoding(t *testing.T) {
	// The canonical bytes of [1.5, null, NaN, -0] are fixed by the spec.
	v := ExplicitValues([]NullFloat{{V: 1.5, Valid: true}, {}, {V: math.NaN(), Valid: true}, {V: math.Copysign(0, -1), Valid: true}})
	want := []byte{
		1, 0x3f, 0xf8, 0, 0, 0, 0, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0,
		1, 0x7f, 0xf8, 0, 0, 0, 0, 0, 0,
		1, 0x80, 0, 0, 0, 0, 0, 0, 0,
	}
	var got []byte
	for _, e := range v.explicit {
		got = appendCanon(got, e)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("canonical f64 bytes % x, want % x", got, want)
	}
	s := ExplicitValues([]NullString{{V: "ab", Valid: true}, {}})
	got = nil
	for _, e := range s.explicit {
		got = appendCanon(got, e)
	}
	if want := []byte{1, 0, 0, 0, 2, 'a', 'b', 0, 0, 0, 0, 0}; !bytes.Equal(got, want) {
		t.Fatalf("canonical str bytes % x, want % x", got, want)
	}
}

// expandPooledCards rewrites every pooled card of a cfref golden (a line
// `    "<card text>"` inside a "cards" array) as the full card object the
// plain dump prints: {"n": <position>, "raw": <card text>, <the fields of
// pool/card-<sha256>.json>}.
func expandPooledCards(t *testing.T, golden []byte) []byte {
	t.Helper()
	lines := bytes.Split(golden, []byte("\n"))
	n, inCards := 0, false
	for i, ln := range lines {
		switch {
		case bytes.HasSuffix(ln, []byte(`"cards": [`)):
			inCards, n = true, 0
			continue
		case inCards && bytes.HasPrefix(ln, []byte("   ]")):
			inCards = false
			continue
		case !inCards || !bytes.HasPrefix(ln, []byte("    ")):
			continue
		}
		n++
		body := bytes.TrimPrefix(ln, []byte("    "))
		comma := bytes.HasSuffix(body, []byte(","))
		body = bytes.TrimSuffix(body, []byte(","))
		if body[0] != '"' {
			continue
		}
		var raw Latin1
		if err := raw.UnmarshalJSON(body); err != nil {
			t.Fatal(err)
		}
		sum := sha256.Sum256([]byte(raw))
		pb, err := ReadGolden(filepath.Join(cfrefPoolDir(), "card-"+hex.EncodeToString(sum[:])+".json"))
		if err != nil {
			t.Fatal(err)
		}
		out := fmt.Appendf(nil, "    {\"n\": %d, \"raw\": %s, %s", n, body, bytes.TrimSuffix(pb, []byte("\n"))[1:])
		if comma {
			out = append(out, ',')
		}
		lines[i] = out
	}
	return bytes.Join(lines, []byte("\n"))
}
