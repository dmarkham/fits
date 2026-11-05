package stretch

import (
	"crypto/sha256"
	"encoding/binary"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"math"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/dmarkham/fits/internal/testutil"
)

// Golden loaders. Every expected value comes from an external truth source
// (siril-cli, the verbatim Siril 1.4.0 oracle in testdata/sref, or
// opencv-python); see testdata/PROVENANCE.md. A missing golden is a test
// failure, never a skip.
//
// Format (testdata/goldenfmt.py, which writes it): an array is its shape,
// the sha256 of the canonical encoding of every element (float32: big-endian
// IEEE bits, every NaN as 0x7fc00000, so -0 != +0; uint16: big-endian) and
// the canonical words in hex, either for every element or for a
// deterministic sample of indices. A hash mismatch fails; the sample makes
// the failure name pixels.

const goldenDir = "testdata/golden"

func readGolden(t *testing.T, rel string, v any) {
	t.Helper()
	path := filepath.Join(goldenDir, rel)
	data, err := testutil.ReadGolden(path)
	if err != nil {
		t.Fatalf("golden %s missing: %v (regenerate: see testdata/PROVENANCE.md)", path, err)
	}
	if err := json.Unmarshal(data, v); err != nil {
		t.Fatalf("golden %s: %v", path, err)
	}
}

// canonNaN32 is the one NaN the canonical encoding uses.
const canonNaN32 = 0x7fc00000

func canon32(v float32) uint32 {
	if v != v {
		return canonNaN32
	}
	return math.Float32bits(v)
}

// sha32 is the canonical digest of channel-major float32 data.
func sha32(chans ...[]float32) string {
	h := sha256.New()
	var b [4]byte
	for _, ch := range chans {
		for _, v := range ch {
			binary.BigEndian.PutUint32(b[:], canon32(v))
			h.Write(b[:])
		}
	}
	return hex.EncodeToString(h.Sum(nil))
}

func sha16(v []uint16) string {
	h := sha256.New()
	var b [2]byte
	for _, x := range v {
		binary.BigEndian.PutUint16(b[:], x)
		h.Write(b[:])
	}
	return hex.EncodeToString(h.Sum(nil))
}

// goldArray is one array in the golden format.
type goldArray struct {
	Shape  []int  `json:"shape"`
	SHA256 string `json:"sha256"`
	Idx    []int  `json:"idx"`
	Bits   string `json:"bits"`
}

func (a goldArray) size() int { return prod(a.Shape) }

func prod(s []int) int {
	n := 1
	for _, d := range s {
		n *= d
	}
	return n
}

// words decodes Bits into `width`-hex-digit words.
func (a goldArray) words(t *testing.T, width int) []uint32 {
	t.Helper()
	if len(a.Bits)%width != 0 {
		t.Fatalf("golden bits length %d is not a multiple of %d", len(a.Bits), width)
	}
	out := make([]uint32, len(a.Bits)/width)
	for i := range out {
		v, err := strconv.ParseUint(a.Bits[i*width:(i+1)*width], 16, 32)
		if err != nil {
			t.Fatalf("golden bits: %v", err)
		}
		out[i] = uint32(v)
	}
	n := len(a.Idx)
	if a.Idx == nil {
		n = a.size()
	}
	if len(out) != n {
		t.Fatalf("golden has %d words for %d values", len(out), n)
	}
	return out
}

// explicit32 returns every value of a small (explicit) float32 array.
func (a goldArray) explicit32(t *testing.T) []float32 {
	t.Helper()
	if a.Idx != nil {
		t.Fatalf("golden array of %d values is digested, not explicit", a.size())
	}
	w := a.words(t, 8)
	out := make([]float32, len(w))
	for i, u := range w {
		out[i] = math.Float32frombits(u)
	}
	if sha32(out) != a.SHA256 {
		t.Fatalf("golden array: explicit values do not match its own sha256")
	}
	return out
}

func (a goldArray) explicit16(t *testing.T) []uint16 {
	t.Helper()
	if a.Idx != nil {
		t.Fatalf("golden array of %d values is digested, not explicit", a.size())
	}
	w := a.words(t, 4)
	out := make([]uint16, len(w))
	for i, u := range w {
		out[i] = uint16(u)
	}
	return out
}

// pixelName formats flat index k of shape for a failure message.
func pixelName(shape []int, k int) string {
	if len(shape) == 3 {
		plane := shape[1] * shape[2]
		return fmt.Sprintf("%d (ch%d y%d x%d)", k, k/plane, k%plane/shape[2], k%shape[2])
	}
	return strconv.Itoa(k)
}

// checkF32 fails unless got (channel-major) is exactly the golden array:
// same element count, same canonical sha256, and the sampled values (all of
// them for an explicit array) bit for bit. NaN equals NaN; -0 != +0.
func checkF32(t *testing.T, label string, want goldArray, got ...[]float32) bool {
	t.Helper()
	var flat []float32
	for _, ch := range got {
		flat = append(flat, ch...)
	}
	if len(flat) != want.size() {
		t.Errorf("%s: %d values, golden shape %v has %d", label, len(flat), want.Shape, want.size())
		return false
	}
	words := want.words(t, 8)
	idx := want.Idx
	var bad []string
	nbad := 0
	for i, w := range words {
		k := i
		if idx != nil {
			k = idx[i]
		}
		if g := canon32(flat[k]); g != w {
			nbad++
			if len(bad) < 6 {
				bad = append(bad, fmt.Sprintf("pixel %s: got %v (%08x), want %v (%08x)",
					pixelName(want.Shape, k), flat[k], g, math.Float32frombits(w), w))
			}
		}
	}
	sum := sha32(flat)
	if sum == want.SHA256 && nbad == 0 {
		return true
	}
	what := "all values"
	if idx != nil {
		what = fmt.Sprintf("%d sampled values", len(idx))
	}
	t.Errorf("%s: sha256 %.16s, golden %.16s; %d of %s differ:\n  %s", label, sum, want.SHA256,
		nbad, what, strings.Join(bad, "\n  "))
	return false
}

func checkU16(t *testing.T, label string, want goldArray, got []uint16) {
	t.Helper()
	if len(got) != want.size() {
		t.Errorf("%s: %d values, golden shape %v has %d", label, len(got), want.Shape, want.size())
		return
	}
	words := want.words(t, 4)
	var bad []string
	nbad := 0
	for i, w := range words {
		k := i
		if want.Idx != nil {
			k = want.Idx[i]
		}
		if uint32(got[k]) != w {
			nbad++
			if len(bad) < 6 {
				bad = append(bad, fmt.Sprintf("value %s: got %d, want %d", pixelName(want.Shape, k), got[k], w))
			}
		}
	}
	if sum := sha16(got); sum != want.SHA256 || nbad > 0 {
		t.Errorf("%s: sha256 %.16s, golden %.16s; %d sampled values differ:\n  %s", label, sum, want.SHA256,
			nbad, strings.Join(bad, "\n  "))
	}
}

// goldInput is an input record (goldenfmt.py): inline values, a committed
// raw file (inputs/*.f32|*.u16, little-endian) transformed by Ops, or a
// constant. SHA256 is the digest of the final array.
type goldInput struct {
	goldArray
	File      string            `json:"file"`
	FileShape []int             `json:"file_shape"`
	Ops       []json.RawMessage `json:"ops"`
	Const     string            `json:"const"`
}

func readRaw(t *testing.T, rel string) []byte {
	t.Helper()
	b, err := os.ReadFile(filepath.Join(goldenDir, rel))
	if err != nil {
		t.Fatalf("golden input %s missing: %v", rel, err)
	}
	return b
}

// f32 returns the input as float32 channels (shape [nch, ny, nx]; other
// shapes come back as one flat channel, nx = ny = 0).
func (in goldInput) f32(t *testing.T) (ch [][]float32, nx, ny int) {
	t.Helper()
	var flat []float32
	shape := in.Shape
	switch {
	case in.File != "":
		if !strings.HasSuffix(in.File, ".f32") {
			t.Fatalf("input %s is not float32", in.File)
		}
		b := readRaw(t, in.File)
		if len(b) != 4*prod(in.FileShape) {
			t.Fatalf("input %s: %d bytes for shape %v", in.File, len(b), in.FileShape)
		}
		flat = make([]float32, len(b)/4)
		for i := range flat {
			flat[i] = math.Float32frombits(binary.LittleEndian.Uint32(b[4*i:]))
		}
		shape = in.FileShape
	case in.Const != "":
		v := hexF32(t, in.Const)
		flat = make([]float32, prod(in.Shape))
		for i := range flat {
			flat[i] = v
		}
	default:
		flat = in.explicit32(t)
	}
	if len(shape) != 3 {
		if len(in.Ops) > 0 {
			t.Fatalf("ops on a %d-D input", len(shape))
		}
		ch = [][]float32{flat}
	} else {
		n := shape[1] * shape[2]
		for c := 0; c < shape[0]; c++ {
			ch = append(ch, flat[c*n:(c+1)*n])
		}
		ny, nx = shape[1], shape[2]
		for _, raw := range in.Ops {
			ch, nx, ny = derive(t, raw, ch, nx, ny)
		}
	}
	if s := sha32(ch...); s != in.SHA256 {
		t.Fatalf("rebuilt input (%s %s) has sha256 %s, golden %s", in.File, opsString(in.Ops), s, in.SHA256)
	}
	return ch, nx, ny
}

// flat returns every input value in order (for vector inputs).
func (in goldInput) flat(t *testing.T) []float32 {
	t.Helper()
	ch, _, _ := in.f32(t)
	var out []float32
	for _, c := range ch {
		out = append(out, c...)
	}
	return out
}

func (in goldInput) u16(t *testing.T) []uint16 {
	t.Helper()
	var out []uint16
	if in.File != "" {
		b := readRaw(t, in.File)
		if !strings.HasSuffix(in.File, ".u16") || len(b) != 2*prod(in.FileShape) {
			t.Fatalf("input %s: not %v uint16", in.File, in.FileShape)
		}
		out = make([]uint16, len(b)/2)
		for i := range out {
			out[i] = binary.LittleEndian.Uint16(b[2*i:])
		}
	} else {
		out = in.explicit16(t)
	}
	if s := sha16(out); s != in.SHA256 {
		t.Fatalf("input %s has sha256 %s, golden %s", in.File, s, in.SHA256)
	}
	return out
}

func opsString(ops []json.RawMessage) string {
	s := make([]string, len(ops))
	for i, o := range ops {
		s[i] = string(o)
	}
	return strings.Join(s, " ")
}

func hexF32(t *testing.T, s string) float32 {
	t.Helper()
	u, err := strconv.ParseUint(s, 16, 32)
	if err != nil {
		t.Fatalf("bad float32 hex %q", s)
	}
	return math.Float32frombits(uint32(u))
}

// derive applies one input transformation (goldenfmt.py derive) to fresh
// copies of the channels.
func derive(t *testing.T, raw json.RawMessage, ch [][]float32, nx, ny int) ([][]float32, int, int) {
	t.Helper()
	var op []json.RawMessage
	if err := json.Unmarshal(raw, &op); err != nil || len(op) == 0 {
		t.Fatalf("bad input op %s", raw)
	}
	var kind string
	arg := func(i int, v any) {
		if i >= len(op) || json.Unmarshal(op[i], v) != nil {
			t.Fatalf("bad input op %s", raw)
		}
	}
	arg(0, &kind)
	out := copyPlanes(ch)
	switch kind {
	case "channels":
		var keep []int
		arg(1, &keep)
		sel := make([][]float32, len(keep))
		for i, c := range keep {
			sel[i] = out[c]
		}
		return sel, nx, ny
	case "fill":
		var c int
		var h string
		arg(1, &c)
		arg(2, &h)
		v := hexF32(t, h)
		for i := range out[c] {
			out[c][i] = v
		}
	case "nanpattern":
		nan := float32(math.NaN())
		for c := range out {
			for i := range out[c] {
				if (i+4*c)%11 == 0 {
					out[c][i] = nan
				}
			}
			out[c][100] = nan
		}
	case "crop":
		var y0, y1, x0, x1 int
		arg(1, &y0)
		arg(2, &y1)
		arg(3, &x0)
		arg(4, &x1)
		for c := range out {
			p := make([]float32, 0, (y1-y0)*(x1-x0))
			for y := y0; y < y1; y++ {
				p = append(p, ch[c][y*nx+x0:y*nx+x1]...)
			}
			out[c] = p
		}
		return out, x1 - x0, y1 - y0
	case "mul", "add":
		var h string
		arg(1, &h)
		v := hexF32(t, h)
		for c := range out {
			for i, x := range out[c] {
				if kind == "mul" {
					out[c][i] = x * v
				} else {
					out[c][i] = x + v
				}
			}
		}
	case "put":
		var c int
		var idx []int
		var vals []string
		arg(1, &c)
		arg(2, &idx)
		arg(3, &vals)
		for k := range out {
			if c >= 0 && k != c {
				continue
			}
			for j, i := range idx {
				out[k][i] = hexF32(t, vals[j])
			}
		}
	default:
		t.Fatalf("unknown input op %q", kind)
	}
	return out, nx, ny
}

// asParam is one parameter set Siril's autostretch computed:
// [retval, shadows, midtones, highlights], the floats as hex bits.
type asParam struct {
	Retval                        int
	Shadows, Midtones, Highlights float32
}

func (p *asParam) UnmarshalJSON(b []byte) error {
	var raw []json.RawMessage
	if err := json.Unmarshal(b, &raw); err != nil || len(raw) != 4 {
		return fmt.Errorf("parameter set %s: want [retval, s, m, h]", b)
	}
	if err := json.Unmarshal(raw[0], &p.Retval); err != nil {
		return err
	}
	for i, dst := range []*float32{&p.Shadows, &p.Midtones, &p.Highlights} {
		var h string
		if err := json.Unmarshal(raw[i+1], &h); err != nil {
			return err
		}
		u, err := strconv.ParseUint(h, 16, 32)
		if err != nil {
			return err
		}
		*dst = math.Float32frombits(uint32(u))
	}
	return nil
}

// same32 compares two float32 values by canonical bits (NaN == NaN, -0 != +0).
func same32(a, b float32) bool { return canon32(a) == canon32(b) }

// srefCase is one golden written by gen_siril_goldens.py from the
// verbatim Siril 1.4.0 oracle (the siril-cli preset goldens have the same
// layout).
type srefCase struct {
	Name string    `json:"name"`
	Op   string    `json:"op"`
	Args []float64 `json:"args"`
	// ArgsF64 holds every argument as float64 bits in hex when one of them
	// is NaN or ±Inf (JSON has neither; Args has null there).
	ArgsF64 []string  `json:"args_f64"`
	Command string    `json:"command"`
	Note    string    `json:"note"`
	Input   goldInput `json:"input"`
	Out     goldArray `json:"out"`
	Params  []asParam `json:"params"`
	// Preset names the DefaultPresets entry a preset_* golden is for.
	Preset string `json:"preset"`
	OpenCV string `json:"opencv"`
	// ClipLimitInt is OpenCV's integer clip limit of a clahe_clip_* case.
	ClipLimitInt int `json:"clip_limit_int"`
}

func loadSref(t *testing.T, name string) srefCase {
	t.Helper()
	var c srefCase
	readGolden(t, filepath.Join("sref", name+".json"), &c)
	if c.Name != name {
		t.Fatalf("golden sref/%s.json is named %q", name, c.Name)
	}
	if c.ArgsF64 != nil {
		if len(c.ArgsF64) != len(c.Args) {
			t.Fatalf("golden sref/%s.json: %d args_f64 for %d args", name, len(c.ArgsF64), len(c.Args))
		}
		for i, h := range c.ArgsF64 {
			u, err := strconv.ParseUint(h, 16, 64)
			if err != nil || len(h) != 16 {
				t.Fatalf("golden sref/%s.json: args_f64[%d] = %q", name, i, h)
			}
			c.Args[i] = math.Float64frombits(u)
		}
	}
	return c
}

// loadCLI returns the preset golden testdata/golden/stretch_<name>.json.
func loadCLI(t *testing.T, name string) srefCase {
	t.Helper()
	var c srefCase
	readGolden(t, "stretch_"+name+".json", &c)
	return c
}

// input returns fresh channels of the case's input.
func (c srefCase) input(t *testing.T) (ch [][]float32, nx, ny int) {
	t.Helper()
	ch, nx, ny = c.Input.f32(t)
	return copyPlanes(ch), nx, ny
}

// check compares channels with the case's output golden.
func (c srefCase) check(t *testing.T, label string, got [][]float32) bool {
	t.Helper()
	if len(c.Out.Shape) == 3 && len(got) != c.Out.Shape[0] {
		t.Errorf("%s: %d channels, golden %v", label, len(got), c.Out.Shape)
		return false
	}
	return checkF32(t, label, c.Out, got...)
}

func copyPlanes(src [][]float32) [][]float32 {
	out := make([][]float32, len(src))
	for i, p := range src {
		out[i] = append([]float32(nil), p...)
	}
	return out
}

// ulpDiff is the distance in float32 units in the last place.
func ulpDiff(a, b float32) int64 {
	if a == b || (a != a && b != b) {
		return 0
	}
	if a != a || b != b {
		return math.MaxInt64
	}
	ia, ib := int64(math.Float32bits(a)), int64(math.Float32bits(b))
	if ia < 0x80000000 == (ib < 0x80000000) {
		if ia > ib {
			return ia - ib
		}
		return ib - ia
	}
	// opposite signs: distance through zero
	return (ia & 0x7fffffff) + (ib & 0x7fffffff)
}
