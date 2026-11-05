package transform

import (
	"bytes"
	"compress/zlib"
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"io"
	"math"
	"os"
	"path/filepath"
	"testing"
)

// The mpmath reference tables under testdata/wref/paper and
// testdata/wref/libm are stored as a small metadata <name>.json and a
// <name>.bin payload holding every double as its exact IEEE-754 bits
// (big-endian float64, zlib), written by testdata/wref/f64pack.py; see the
// format there. They used to be decimal JSON rows holding the same
// doubles (checked value for value, bit for bit, when the format changed).

// f64packMeta is the part of a table's metadata the reader checks; the
// generator-specific fields (row counts, source) are read by the caller.
type f64packMeta struct {
	Payload     string `json:"payload"`
	Values      int    `json:"values"`
	SHA256      string `json:"sha256"`
	ZlibVersion string `json:"zlib_version"`
}

// loadF64Pack reads metaPath into meta (which must embed or repeat the
// f64packMeta fields) and returns the payload's doubles, after checking
// the payload's length and sha256 against the metadata. A missing or
// corrupt table is a test failure.
func loadF64Pack(t *testing.T, metaPath string, meta any) []float64 {
	t.Helper()
	vals, _ := readF64Pack(t, metaPath, meta)
	return vals
}

func readF64Pack(t *testing.T, metaPath string, meta any) ([]float64, []byte) {
	t.Helper()
	data, err := os.ReadFile(metaPath)
	if err != nil {
		t.Fatalf("reference table missing: %v (regenerate with its generator in testdata/wref)", err)
	}
	if err := json.Unmarshal(data, meta); err != nil {
		t.Fatalf("%s: %v", metaPath, err)
	}
	var m f64packMeta
	if err := json.Unmarshal(data, &m); err != nil {
		t.Fatalf("%s: %v", metaPath, err)
	}
	z, err := os.ReadFile(filepath.Join(filepath.Dir(metaPath), m.Payload))
	if err != nil || m.Payload == "" {
		t.Fatalf("%s: payload %q missing: %v", metaPath, m.Payload, err)
	}
	zr, err := zlib.NewReader(bytes.NewReader(z))
	if err != nil {
		t.Fatalf("%s: payload: %v", metaPath, err)
	}
	raw, err := io.ReadAll(zr) // zlib checks its Adler-32 at EOF
	if err != nil {
		t.Fatalf("%s: payload: %v", metaPath, err)
	}
	if len(raw) != 8*m.Values {
		t.Fatalf("%s: payload is %d bytes, metadata says %d doubles", metaPath, len(raw), m.Values)
	}
	if sum := sha256.Sum256(raw); hex.EncodeToString(sum[:]) != m.SHA256 {
		t.Fatalf("%s: payload sha256 %x, metadata says %s", metaPath, sum, m.SHA256)
	}
	vals := make([]float64, m.Values)
	for i := range vals {
		vals[i] = math.Float64frombits(binary.BigEndian.Uint64(raw[8*i:]))
	}
	return vals, raw
}

// f64packTake cuts n*w doubles off the front of *vals.
func f64packTake(t *testing.T, vals *[]float64, n, w int) []float64 {
	t.Helper()
	if n < 0 || n*w > len(*vals) {
		t.Fatalf("reference table: %d rows of %d need %d doubles, %d left", n, w, n*w, len(*vals))
	}
	v := (*vals)[:n*w]
	*vals = (*vals)[n*w:]
	return v
}

// f64packRows4 and f64packRows3 cut n rows of 4 (3) doubles off *vals.
func f64packRows4(t *testing.T, vals *[]float64, n int) [][4]float64 {
	t.Helper()
	v := f64packTake(t, vals, n, 4)
	rows := make([][4]float64, n)
	for i := range rows {
		copy(rows[i][:], v[4*i:])
	}
	return rows
}

func f64packRows3(t *testing.T, vals *[]float64, n int) [][3]float64 {
	t.Helper()
	v := f64packTake(t, vals, n, 3)
	rows := make([][3]float64, n)
	for i := range rows {
		copy(rows[i][:], v[3*i:])
	}
	return rows
}

// f64packSame requires a regenerated table (genMeta, in a temp directory)
// to equal the committed one (committedMeta): the metadata byte for byte
// (with the zlib version normalised when the local zlib differs), the
// decoded payload byte for byte, and the compressed payload too when the
// local zlib is the one that wrote the committed table.
func f64packSame(t *testing.T, genMeta, committedMeta string) {
	t.Helper()
	var gm, cm f64packMeta
	_, graw := readF64Pack(t, genMeta, &gm)
	_, craw := readF64Pack(t, committedMeta, &cm)
	gj, _ := os.ReadFile(genMeta)
	cj, _ := os.ReadFile(committedMeta)
	sameZlib := gm.ZlibVersion == cm.ZlibVersion
	if !sameZlib {
		t.Logf("%s: local zlib %s, committed with %s: comparing the decoded payload only", committedMeta, gm.ZlibVersion, cm.ZlibVersion)
		gj = bytes.Replace(gj, []byte(`"zlib_version": "`+gm.ZlibVersion+`"`), []byte(`"zlib_version": "`+cm.ZlibVersion+`"`), 1)
	}
	if !bytes.Equal(gj, cj) {
		t.Errorf("%s: regenerated metadata differs from the committed one", committedMeta)
	}
	if !bytes.Equal(graw, craw) {
		t.Errorf("%s: regenerated payload differs from the committed one", committedMeta)
	}
	if sameZlib {
		gz, _ := os.ReadFile(filepath.Join(filepath.Dir(genMeta), gm.Payload))
		cz, _ := os.ReadFile(filepath.Join(filepath.Dir(committedMeta), cm.Payload))
		if !bytes.Equal(gz, cz) {
			t.Errorf("%s: same zlib %s but the compressed payload differs from the committed one", committedMeta, gm.ZlibVersion)
		}
	}
}
