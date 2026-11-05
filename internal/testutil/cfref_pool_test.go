package testutil

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"os"
	"path/filepath"
	"regexp"
	"strings"
	"testing"
)

var poolRefRe = regexp.MustCompile(`"sha256": "([0-9a-f]{64})", "pool": "([a-z0-9]+)"`)

// TestCfrefPoolComplete: every file in testdata/cfref/golden/pool is
// referenced by a committed golden and every reference has its file, and
// each pool file loads (count and sha256 of its values equal its name).
func TestCfrefPoolComplete(t *testing.T) {
	gold := filepath.Join(cfrefPoolDir(), "..")
	goldens, _ := filepath.Glob(filepath.Join(gold, "*.json"))
	twins, _ := filepath.Glob(filepath.Join(gold, "twins", "*.json"))
	bad, _ := filepath.Glob(filepath.Join(gold, "..", "..", "comp_bad", "golden", "*.json"))
	goldens = append(append(goldens, twins...), bad...)
	if len(goldens) < 200 {
		t.Fatalf("found %d cfref goldens", len(goldens))
	}
	refs := map[string]bool{}
	cards := 0
	for _, g := range goldens {
		b, err := os.ReadFile(g)
		if err != nil {
			t.Fatal(err)
		}
		for _, m := range poolRefRe.FindAllStringSubmatch(string(b), -1) {
			refs[m[2]+"-"+m[1]+".json"] = true
		}
		// pooled cards: a card given as its text alone
		var doc struct {
			HDUs []struct {
				Cards []json.RawMessage `json:"cards"`
			} `json:"hdus"`
		}
		if err := json.Unmarshal(b, &doc); err != nil {
			t.Fatalf("%s: %v", g, err)
		}
		for _, h := range doc.HDUs {
			for _, c := range h.Cards {
				if len(c) == 0 || c[0] != '"' {
					continue
				}
				var raw Latin1
				if err := raw.UnmarshalJSON(c); err != nil {
					t.Fatalf("%s: %v", g, err)
				}
				sum := sha256.Sum256([]byte(raw))
				refs["card-"+hex.EncodeToString(sum[:])+".json"] = true
				cards++
			}
		}
		// and every golden decodes through the pool
		if _, err := DecodeCfref(b); err != nil {
			t.Errorf("%s: %v", g, err)
		}
	}
	if cards < 1000 {
		t.Errorf("only %d pooled cards referenced", cards)
	}
	files := GlobGoldens(filepath.Join(cfrefPoolDir(), "*.json"))
	have := map[string]bool{}
	for _, f := range files {
		have[filepath.Base(f)] = true
		if !refs[filepath.Base(f)] {
			t.Errorf("pool file %s is referenced by no golden", filepath.Base(f))
		}
	}
	for r := range refs {
		if !have[r] {
			t.Errorf("golden references missing pool file %s", r)
		}
	}
	if len(refs) < 100 {
		t.Fatalf("only %d pool references", len(refs))
	}
	// Each file's values hash to its name (loading checks it).
	for _, f := range files {
		kind, sha, _ := strings.Cut(strings.TrimSuffix(filepath.Base(f), ".json"), "-")
		if kind == "card" {
			continue // decoded through the goldens above
		}
		b, _ := ReadGolden(f)
		r := poolRef{N: int64(strings.Count(string(b), ",") + 1), SHA256: sha, Pool: kind}
		var err error
		switch kind {
		case "f64":
			_, err = loadPooled[NullFloat](cfrefPoolDir(), r)
		case "i64":
			_, err = loadPooled[NullInt64](cfrefPoolDir(), r)
		case "u64":
			_, err = loadPooled[NullUint64](cfrefPoolDir(), r)
		default:
			t.Errorf("pool file %s: unexpected kind", f)
		}
		if err != nil {
			t.Errorf("%v", err)
		}
	}
}

// TestCfrefPoolRefusesBadFiles: a pool file whose values do not hash to
// its name (one value changed), a wrong element count or kind, and a
// missing file are errors, so a pooled golden is as exact as an explicit
// one.
func TestCfrefPoolRefusesBadFiles(t *testing.T) {
	dir := t.TempDir()
	v := []NullInt64{{V: 1, Valid: true}, {}, {V: -7, Valid: true}}
	sha := digestOf(v)
	write := func(body string) {
		if err := os.WriteFile(filepath.Join(dir, "i64-"+sha+".json"), []byte(body), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	write("[1,null,-7]\n")
	if got, err := loadPooled[NullInt64](dir, poolRef{N: 3, SHA256: sha, Pool: "i64"}); err != nil || len(got) != 3 || got[2].V != -7 || got[1].Valid {
		t.Fatalf("good pool file: %v %v", got, err)
	}
	for name, c := range map[string]struct {
		body string
		ref  poolRef
		want string
	}{
		"value changed": {"[1,null,-8]\n", poolRef{N: 3, SHA256: sha, Pool: "i64"}, "hashes to"},
		"null changed":  {"[1,0,-7]\n", poolRef{N: 3, SHA256: sha, Pool: "i64"}, "hashes to"},
		"count":         {"[1,null,-7]\n", poolRef{N: 4, SHA256: sha, Pool: "i64"}, "holds 3 values"},
		"kind":          {"[1,null,-7]\n", poolRef{N: 3, SHA256: sha, Pool: "f64"}, "pooled f64 array"},
		"missing":       {"[1,null,-7]\n", poolRef{N: 3, SHA256: strings.Repeat("0", 64), Pool: "i64"}, "missing"},
	} {
		write(c.body)
		poolCache.Delete(filepath.Join(dir, "i64-"+sha+".json"))
		_, err := loadPooled[NullInt64](dir, c.ref)
		if err == nil || !strings.Contains(err.Error(), c.want) {
			t.Errorf("%s: err %v, want one mentioning %q", name, err, c.want)
		}
	}
}
