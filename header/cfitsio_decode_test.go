package header

import (
	"bytes"
	"errors"
	"math"
	"os"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// headerSamples are the samples committed in header/testdata (see
// gen_samples.py and hdrref write).
var headerSamples = []string{"malformed.fits", "continue.fits", "adversarial.fits", "cfwritten.fits"}

type goldenSample struct {
	name string
	data []byte
	cf   *testutil.CfrefFile
	ls   *longstrFile
}

func loadHeaderSample(t *testing.T, name string) goldenSample {
	t.Helper()
	p := headerTestdata(t, name)
	testutil.RequireFile(t, p)
	data, err := os.ReadFile(p)
	if err != nil {
		t.Fatal(err)
	}
	gp := headerTestdata(t, "golden", name+".json")
	testutil.RequireFile(t, gp)
	gb, err := os.ReadFile(gp)
	if err != nil {
		t.Fatal(err)
	}
	cf, err := testutil.DecodeCfref(gb)
	if err != nil {
		t.Fatalf("decode %s: %v", gp, err)
	}
	return goldenSample{name, data, cf, loadLongstrGolden(t, headerTestdata(t, "golden", name+".longstr.json"))}
}

// repoSamples lists every repo testdata/ sample that has a shared cfref
// golden, with this package's longstr golden for it.
func repoSamples(t *testing.T) []goldenSample {
	t.Helper()
	gl, err := filepath.Glob(filepath.Join(testutil.CfrefDir(t), "golden", "*.json"))
	if err != nil || len(gl) == 0 {
		t.Fatalf("no cfref goldens found (%v)", err)
	}
	sort.Strings(gl)
	var out []goldenSample
	for _, g := range gl {
		name := strings.TrimSuffix(filepath.Base(g), ".json")
		data, err := os.ReadFile(testutil.TestdataPath(t, name))
		if err != nil {
			t.Fatalf("sample for golden %s: %v", g, err)
		}
		cf := testutil.LoadCfrefGolden(t, name)
		ls := loadLongstrGolden(t, headerTestdata(t, "golden", "testdata", name+".longstr.json"))
		out = append(out, goldenSample{name, data, cf, ls})
	}
	return out
}

func (s goldenSample) headerBytes(t *testing.T, i int) []byte {
	t.Helper()
	h := s.cf.HDUs[i]
	if h.HeaderStart < 0 || h.DataStart > int64(len(s.data)) || h.HeaderStart >= h.DataStart {
		t.Fatalf("%s HDU %d: bad cfitsio addresses %d..%d (file %d bytes)", s.name, i, h.HeaderStart, h.DataStart, len(s.data))
	}
	return s.data[h.HeaderStart:h.DataStart]
}

func (s goldenSample) check(t *testing.T) {
	t.Helper()
	if len(s.ls.HDUs) != len(s.cf.HDUs) {
		t.Fatalf("%s: longstr golden has %d HDUs, cfref %d", s.name, len(s.ls.HDUs), len(s.cf.HDUs))
	}
	for i := range s.cf.HDUs {
		raw := s.headerBytes(t, i)
		cards, end, err := ParseCards(raw)
		if err != nil {
			t.Errorf("%s HDU %d: ParseCards: %v (cfitsio reads this header)", s.name, i, err)
			continue
		}
		if want := s.cf.HDUs[i].NKeys; end < want || !trailingBlanks(cards, want) {
			t.Errorf("%s HDU %d: END at card %d, cfitsio counts %d keys (only blank cards may come between)", s.name, i, end, want)
		}
		compareWithCfitsio(t, s.name+" HDU "+strconv.Itoa(i), cards, &s.cf.HDUs[i], &s.ls.HDUs[i])
		// Header.Float on integer cards is fits_read_key_dbl.
		h := FromCards(cards)
		seen := map[string]bool{}
		for _, c := range s.cf.HDUs[i].Cards {
			k := strings.ToUpper(string(c.Key))
			if seen[k] || c.Type != "I" || c.Int == nil || c.Float == nil || strings.Contains(k, " ") {
				seen[k] = true
				continue
			}
			seen[k] = true
			if f, err := h.Float(k); err != nil || math.Float64bits(f) != math.Float64bits(float64(*c.Float)) {
				t.Errorf("%s HDU %d: Float(%s) = %v, %v; cfitsio fits_read_key_dbl %v", s.name, i, k, f, err, float64(*c.Float))
			}
		}
	}
}

// TestDecodeMatchesCfitsioHeaderSamples: malformed cards, CONTINUE edge
// cases and a cfitsio-written header, each decoded exactly as cfitsio
// reads the same bytes (rules in compareWithCfitsio).
func TestDecodeMatchesCfitsioHeaderSamples(t *testing.T) {
	for _, name := range headerSamples {
		t.Run(name, func(t *testing.T) { loadHeaderSample(t, name).check(t) })
	}
}

// TestDecodeMatchesCfitsioRepoSamples runs the same comparison over every
// header of every FITS sample in the repository's testdata/.
func TestDecodeMatchesCfitsioRepoSamples(t *testing.T) {
	for _, s := range repoSamples(t) {
		t.Run(s.name, func(t *testing.T) { s.check(t) })
	}
}

// TestLongstrGoldensCoverRepoSamples: every repo sample with a cfref
// golden has a longstr golden here and vice versa.
func TestLongstrGoldensCoverRepoSamples(t *testing.T) {
	gl, _ := filepath.Glob(filepath.Join(testutil.CfrefDir(t), "golden", "*.json"))
	ll, _ := filepath.Glob(headerTestdata(t, "golden", "testdata", "*.longstr.json"))
	want := map[string]bool{}
	for _, g := range gl {
		want[strings.TrimSuffix(filepath.Base(g), ".json")] = true
	}
	have := map[string]bool{}
	for _, l := range ll {
		have[strings.TrimSuffix(filepath.Base(l), ".longstr.json")] = true
	}
	for k := range want {
		if !have[k] {
			t.Errorf("no header/testdata/golden/testdata/%s.longstr.json (run make -C header/testdata goldens)", k)
		}
	}
	for k := range have {
		if !want[k] {
			t.Errorf("longstr golden %s has no cfref golden / sample", k)
		}
	}
	if len(want) == 0 {
		t.Fatal("no cfref goldens")
	}
}

// TestReencodeUnmodifiedIsByteIdentical: parsing a header and encoding it
// again without any mutation gives back the exact bytes, for every header
// of every sample (CONTINUE chains, malformed cards and HIERARCH included).
func TestReencodeUnmodifiedIsByteIdentical(t *testing.T) {
	all := repoSamples(t)
	for _, n := range headerSamples {
		all = append(all, loadHeaderSample(t, n))
	}
	for _, s := range all {
		for i := range s.cf.HDUs {
			raw := s.headerBytes(t, i)
			cards, _, err := ParseCards(raw)
			if err != nil {
				t.Errorf("%s HDU %d: %v", s.name, i, err)
				continue
			}
			got, err := Encode(FromCards(cards))
			if err != nil {
				t.Errorf("%s HDU %d: Encode: %v", s.name, i, err)
				continue
			}
			if !bytes.Equal(got, raw) {
				t.Errorf("%s HDU %d: re-encoded header differs from the file\n%s", s.name, i, firstCardDiff(got, raw))
			}
		}
	}
}

func firstCardDiff(got, want []byte) string {
	if len(got) != len(want) {
		return "length " + strconv.Itoa(len(got)/2880) + " blocks, file " + strconv.Itoa(len(want)/2880) + " blocks; " + firstCardDiff(got[:min(len(got), len(want))], want[:min(len(got), len(want))])
	}
	for i := 0; i+CardWidth <= len(got); i += CardWidth {
		if !bytes.Equal(got[i:i+CardWidth], want[i:i+CardWidth]) {
			return "card " + strconv.Itoa(i/CardWidth) + ":\n got " + string(got[i:i+CardWidth]) + "\nfile " + string(want[i:i+CardWidth])
		}
	}
	return "(no card differs)"
}

// TestHeaderGettersMatchCfitsio: the typed Header getters give cfitsio's
// fits_read_key_dbl / TULONGLONG answers for out-of-int64 integers, and
// refuse (with ErrIntOverflow) where cfitsio's integer read overflows.
func TestHeaderGettersMatchCfitsio(t *testing.T) {
	s := loadHeaderSample(t, "malformed.fits")
	cards, _, err := ParseCards(s.headerBytes(t, 0))
	if err != nil {
		t.Fatal(err)
	}
	h := FromCards(cards)
	gh := &s.cf.HDUs[0]
	for _, key := range []string{"BIGINT", "BZEROU", "U64MAX", "NEGBIG"} {
		want := gh.Card(key)
		if want == nil || want.Float == nil || want.Int != nil {
			t.Fatalf("%s: golden lacks the overflow case: %+v", key, want)
		}
		f, err := h.Float(key)
		if err != nil || math.Float64bits(f) != math.Float64bits(float64(*want.Float)) {
			t.Errorf("Float(%s) = %v, %v; cfitsio fits_read_key_dbl %v", key, f, err, float64(*want.Float))
		}
		if _, err := h.Int(key); !errors.Is(err, ErrIntOverflow) {
			t.Errorf("Int(%s) err = %v, want ErrIntOverflow (cfitsio: %s)", key, err, want.ReadError)
		}
		b, err := h.BigInt(key)
		if err != nil || b.String() != strings.TrimPrefix(string(want.Value), "+") {
			t.Errorf("BigInt(%s) = %v, %v; card text %q", key, b, err, want.Value)
		}
		if ls := s.ls.HDUs[0].card(want.N); ls != nil && ls.U64 != nil && (b == nil || !b.IsUint64() || b.Uint64() != *ls.U64) {
			t.Errorf("BigInt(%s) = %v; cfitsio TULONGLONG %d", key, b, *ls.U64)
		}
	}
	// In-range integers: BigInt agrees with Int.
	if b, err := h.BigInt("I64MIN"); err != nil || !b.IsInt64() || b.Int64() != math.MinInt64 {
		t.Errorf("BigInt(I64MIN) = %v, %v", b, err)
	}
	// Unparsed values are kept, exposed, and refused by every typed getter.
	for _, key := range []string{"DATE-OBS", "X", "NANTOK", "FLTOVF"} {
		c, ok := h.Get(key)
		rv, isRaw := c.Value.(RawValue)
		if !ok || c.Type != TypeUnparsed || !isRaw || rv.Text != string(gh.Card(key).Value) {
			t.Errorf("Get(%s) = %+v, want TypeUnparsed with cfitsio's text %q", key, c, gh.Card(key).Value)
		}
		if _, err := h.Float(key); !errors.Is(err, ErrWrongType) {
			t.Errorf("Float(%s) err = %v, want ErrWrongType", key, err)
		}
		if _, err := h.String(key); !errors.Is(err, ErrWrongType) {
			t.Errorf("String(%s) err = %v, want ErrWrongType", key, err)
		}
	}
}

// TestLookupHistoryCommentsMatchCfitsio: case-insensitive lookup finds a
// lower-case keyword cfitsio reads by name, and History / Comments return
// exactly the HISTORY / COMMENT card texts cfitsio reports, in order, for
// the cfitsio-written sample and cfitsio's own test file.
func TestLookupHistoryCommentsMatchCfitsio(t *testing.T) {
	s := loadHeaderSample(t, "malformed.fits")
	cards, _, _ := ParseCards(s.headerBytes(t, 0))
	h := FromCards(cards)
	want := s.cf.HDUs[0].Card("lowkey")
	for _, name := range []string{"LOWKEY", "lowkey", "LowKey"} {
		if v, err := h.Int(name); err != nil || want == nil || want.Int == nil || v != *want.Int {
			t.Errorf("Int(%q) = %d, %v; cfitsio %v", name, v, err, want)
		}
	}
	if v, err := h.Int("HIERARCH ESO DBL SPACE"); err != nil || v != 7 {
		t.Errorf("HIERARCH lookup with collapsed blanks: %d, %v", v, err)
	}
	if _, ok := h.Get("COMMENT"); ok {
		t.Error("Get found a commentary card")
	}

	samples := []goldenSample{loadHeaderSample(t, "cfwritten.fits")}
	for _, r := range repoSamples(t) {
		if r.name == "testprog.std" {
			samples = append(samples, r)
		}
	}
	for _, s := range samples {
		for i := range s.cf.HDUs {
			cards, _, err := ParseCards(s.headerBytes(t, i))
			if err != nil {
				t.Fatal(err)
			}
			h := FromCards(cards)
			for _, kw := range []struct {
				key string
				got []string
			}{{"HISTORY", h.History()}, {"COMMENT", h.Comments()}} {
				var w []string
				for _, c := range s.cf.HDUs[i].Cards {
					if string(c.Key) == kw.key {
						w = append(w, string(c.Comment))
					}
				}
				if strings.Join(kw.got, "\n") != strings.Join(w, "\n") || len(kw.got) != len(w) {
					t.Errorf("%s HDU %d %s:\n got %q\nwant %q", s.name, i, kw.key, kw.got, w)
				}
			}
		}
	}
}

// TestCloneIsIndependent: Set/Add/Delete on a clone leave the original's
// cards and encoding unchanged.
func TestCloneIsIndependent(t *testing.T) {
	s := loadHeaderSample(t, "continue.fits")
	raw := s.headerBytes(t, 0)
	cards, _, _ := ParseCards(raw)
	h := FromCards(cards)
	c := h.Clone()
	if err := c.Set("LONGSTR", "changed", ""); err != nil {
		t.Fatal(err)
	}
	if err := c.Add("HISTORY", "added", ""); err != nil {
		t.Fatal(err)
	}
	if err := c.Delete("FINAL"); err != nil {
		t.Fatal(err)
	}
	got, err := Encode(h)
	if err != nil || !bytes.Equal(got, raw) {
		t.Fatalf("original changed by edits to its clone (err %v)", err)
	}
	if v, _ := c.String("LONGSTR"); v != "changed" {
		t.Fatalf("clone LONGSTR = %q", v)
	}
}

// TestBlockCountForMatchesEncode: BlockCountFor(n) is the number of
// 2880-byte blocks Encode writes for n single-card keywords plus END.
func TestBlockCountForMatchesEncode(t *testing.T) {
	for n := 0; n <= 80; n++ {
		h := New()
		for i := 0; i < n; i++ {
			if err := h.Add("K"+strconv.Itoa(i), int64(i), ""); err != nil {
				t.Fatal(err)
			}
		}
		b, err := Encode(h)
		if err != nil {
			t.Fatal(err)
		}
		if len(b)%2880 != 0 || len(b)/2880 != BlockCountFor(n) {
			t.Fatalf("n=%d: Encode wrote %d bytes, BlockCountFor = %d blocks", n, len(b), BlockCountFor(n))
		}
	}
}
