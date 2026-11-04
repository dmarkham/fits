# healpix/testdata provenance

All expected values come from healpy, never from this package.

## Tools

| tool   | version | where |
|--------|---------|-------|
| Python | 3.11.11 (pyenv) | shared venv `fits-proof-audit/venv` |
| healpy | 1.20.1 (wraps healpix_cxx) | pip, in that venv |
| numpy  | 2.4.6  | pip |
| mpmath | 1.4.1  | pip (exact max_pixrad; the correctly rounded cos/sin/atan2 fed to the healpix_cxx harness and stored next to the polar, switch-over and 2/3 goldens; polar boundary solving) |
| healpix_cxx | 8.3.0 | source tree vendored by Siril, `~/git/dmarkham/Siril/subprojects/healpix_cxx` (Siril git `2cebd0d97615f96abeae10fa596d9088bc4e1834`, GPL-2.0, not copied here), compiled into `hpref/hpref` |
| g++ | 16.2.1 20260810 (system) | builds `hpref/hpref` with `-O2 -std=c++11 -ffp-contract=off` |
| glibc libm | 2.44 | what healpy (and the harness `a`/`v` queries) call for cos/sin/atan2 |

healpix_cxx sources used by `hpref` (`make -C hpref sha`):

```
942a9454e2e94496f3f94845a779760c22666fa4d935a05a8b9391f54a08f3b1  healpix_base.cc
2584388ffa39f11ac892ca5a24583ce6a069a5404029bc804ca53b6994515986  healpix_tables.cc
e14cc9eab5bf6a9664e5b94be486f7872d445b0b1378bb0d59be4ba231340e6f  error_handling.cc
394169810db3392f63fded182a090d3008e826040003105df83f8ff467cc4995  geom_utils.cc
aaaada283644a193c73159c44fd1afe948fd458da125e8189bf8a5b8d4905dbd  pointing.cc
3cc2c8b71bcc280537f29e789c7196739ff423683748a21c6cb5482202de2759  string_utils.cc
77df11e3f582831b33a50cf2bce03dc0ba5f17f9f77a36ae26d8f40d306721bf  healpix_base.h
46b7c2fe3fa3966b00ee9e9c0c8133e222f0e8ee761373b0dcb36673b0ae051a  math_utils.h
```

`gen_golden.py` and `maps/gen_map.py` exit if the versions differ.

## Reproduce

```sh
cd healpix/testdata
make -C hpref                               # healpix_cxx harness (gitignored binary)
<venv>/bin/python gen_golden.py golden      # every golden/*.json
cd maps && <venv>/bin/python gen_map.py     # healpy_map_nside16_* files (TestHealpixMapSamplesReproduce)
```

Both reproduce the committed files byte-for-byte (checked by regenerating
into a scratch directory and `cmp`-ing every file;
`TestHealpixGoldensReproduce`, `TestHealpixMapSamplesReproduce`). The
first part of `gen_golden.py` writes the original fixture set (ang2pix,
pix2ang, ordering, ordering_sampled, neighbors, utilities, xyf) with the
values of the healpy 1.19.0 generation (healpy 1.20.1 gives the same
values), now without indentation, `xyf` as columns.

## Encodings (exact; decoded by `golden_codec_test.go`)

To keep every file well under 1 MB without dropping any value:

- **float arrays** that a test needs element by element (inputs, values
  compared within a ULP budget, residuals) are `"f64:" + base64` of the
  big-endian IEEE float64 bytes (exact, ~10.7 characters per value).
- **integer result arrays** longer than 256 values (pixel indices, xyf,
  neighbours, ring <-> nest) are a digest `{"n", "sha256", "sample"}`:
  sha256 over the big-endian int64 bytes of every value, plus explicit
  `[index, value]` pairs at the first and last 4 indices, 32 evenly spaced
  ones and the first minimum and maximum. The test computes its full
  result and compares it with the digest: any single differing value
  fails (sha256), and the sampled values that differ are named. Shorter
  arrays stay explicit lists. Exception: the pixel arrays of
  `ang2pix_polar.json` and `ang2pix_twothird.json` are explicit lists
  whatever their length, because their tests accept a different pixel at
  a boundary-ambiguous direction (healpix/ambiguity_test.go) and so must
  know every expected value.

Before this encoding the golden directory was 7.9 MB (crmath.json 1.9 MB,
ang2pix_polar.json 1.2 MB); with it, 2.3 MB. crmath.json (the
correct-rounding truth for the package's own correctly rounded
cos/sin/atan2, 0.68 MB) was deleted with those functions in perf/trig.

## Golden files

| files | contents |
|-------|----------|
| `ang2pix_nside*.json`, `pix2ang_nside*.json`, `neighbors_nside*.json`, `ordering_nside*.json`, `ordering_sampled_nside*.json`, `utilities.json`, `xyf_nside*.json` | original set, Nside 1..8192 |
| `ang2pix_random.json` | 3000 seeded-random directions (1000 of them within 0.02 rad of a pole) x 21 Nside (pow2 1..8192, non-pow2 3,5,6,7,12,100,1000,3^10,3·2^27,2^29-1, and 2^20, 2^29) |
| `ang2pix_phiwrap.json` | 60 directions with phi < 0, = 2pi, > 2pi, +-100, +-1e6 x 6 Nside (1..2^29) |
| `vec2pix.json` | 1512 non-unit vectors (lengths 1e-3..1e6, axis and face-boundary vectors) x the same 21 Nside |
| `pixfuncs_nside*.json` | non-pow2 RING Nside incl. 3·2^27 and 2^29-1 (exhaustive when npix <= 1200) and Nside 2^20, 2^29 (ends, zone/face boundaries, random): pix2ang, pix2vec, pix2xyf, ordered neighbours, ring2nest/nest2ring |
| `pix2vec_nside*.json` | pix2vec for the original pixel samples |
| `neighbours_ordered_nside*.json` | get_all_neighbours (order SW,W,NW,N,NE,E,SE,S, -1 kept in place) for every pixel (Nside <= 16) or every face edge/corner pixel (64, 8192) |
| `ang2pix_polar.json` | 4000 directions on and within 2 ULP of polar-cap ring/pixel boundaries (solved in 50-digit mpmath; for each Nside and pole the 100 of 2000 candidate boundaries whose nearest double lies closest) x Nside 8192, 2^20, 2^29-1, 2^29; theta/phi and the unit vector (mpmath, rounded once), and the correctly rounded cos(theta), sin(theta), atan2(vy, vx) (`cos_cr`, `sin_cr`, `atan2_cr`). Expected (explicit lists): healpix_cxx loc2pix (`hpref`) with those correctly rounded cos/sin/atan2. The generator asserts that `hpref` calling glibc reproduces healpy's ang2pix/vec2pix on every input and lists any answer where healpy differs (only possible where glibc is not correctly rounded; 0 for the committed inputs, although glibc's cos/sin are not correctly rounded for 9 of them and atan2 for 3). |
| `branch_switch.json` | Ang2Pix/Vec2Pix at the cos -> sin switch-over points (theta < 0.01, theta > 3.14159-0.01, \|z\| > 0.99) x Nside 8192, 2^20, 2^29-1, 2^29: directions on (and 2 ULP around) polar-cap pixel boundaries with theta in [0.0098, 0.0102], [3.1304, 3.1318] or \|z\| in [0.9898, 0.9902] (both poles); the thetas 0.01, 3.14159-0.01 and pi-0.01 +-3 ULP with phi on a boundary; vectors with healpix nz exactly +-0.99. Only inputs where the two formulas give different pixels are kept (hpref `C`/`S`/`Q`/`P` forced-branch ops, used for selection only), 40 per group and Nside. Each record carries the correctly rounded cos/sin (`cos_cr`, `sin_cr`) or atan2 (`atan2_cr`) it was evaluated with. Expected: healpix_cxx (`hpref` `A`/`V`) with those correctly rounded cos/sin/atan2; the generator asserts `hpref` calling glibc equals healpy and drops any input where healpy differs from the expected value (0 dropped), so every expected value is healpy's own answer. |
| `vec2pix_extreme.json` | unit-scale vectors times 2^k (k = 0, +-600, 1000, 1019, -1000, -1070, and 501..599 / -599..-501 around the 2^+-511 overflow/underflow edge of the squares; exact) x Nside 4, 8192, 2^29; expected = healpy at unit scale; healpy's own answer for the scaled vector is recorded (wrong for 255 of 411) |
| `ang2pix_twothird.json` | Ang2Pix at the equatorial/polar split z = 2/3: theta = 0.8410686705679303 (the only double near acos(2/3) whose correctly rounded cos is the double nearest 2/3) and the doubles 1 and 2 ULP either side, and the 5 doubles around acos(-2/3) (none of which has cos == -2/3), each at phi = k pi/32 (k = 0..63), x every power-of-2 Nside 2^0..2^29, RING and NESTED; expected = healpy ang2pix, explicit lists (only thetas whose glibc cos is correctly rounded are kept: all 10); `cos_cr` is that correctly rounded cos per input. |
| `vec2pix_scaled_healpy.json` | healpy's own vec2pix on the `ang2pix_polar.json` vectors scaled by 2^k, every k = -490..-480 and -495, -500, at Nside 2^29-1 (RING) and 2^29 (RING, NESTED): the indices where its answer differs from its unscaled answer, with the scaled answer (none at -480 and -481; at Nside 2^29-1 RING 1, 2, 4, 10, 17, 28, 38, 39, 75 at -482..-490; 117 / 94 / 94 at -495; 219 / 262 / 262 at -500). `tiny`: healpy vec2pix of (1e-200, 0, 1e-200) and of (1, 0, 1) at Nside 1, 2, 16, 2^20, 2^29, RING and NESTED (RING gives 0 and NESTED another pixel at every Nside > 1). |
| `utilities_ext.json` | nside2npix, nside2pixarea, nside2resol, max_pixrad (+ mpmath exact) for 38 Nside; isnsideok for 137 values; isnpixok/npix2nside |

Known healpy defects recorded in the tests (Go deliberately differs):
`isnpixok(0)` is True; `isnpixok(12*2^58+12)` is True because healpy
rounds npix through a double; `isnpixok(12*(2^29+1)^2)` is True (no upper
Nside limit); `vec2pix` returns wrong pixels when x²+y²+z² underflows or
overflows (e.g. (2^-600, 0, 2^-600)).

## Map samples (`maps/`)

| file | source |
|------|--------|
| `wmap_ilc_9yr_v5.hdr` | REAL: the two header blocks (primary + BINTABLE, bytes 0..5759, verbatim) of the WMAP 9-year ILC map `https://lambda.gsfc.nasa.gov/data/map/dr5/dfp/ilc/wmap_ilc_9yr_v5.fits` (NASA LAMBDA, public domain; Last-Modified 2012-12-11, Content-Length 25174080), fetched with an HTTP range request. sha256 `faa8a0095726300eaeb5cf51dddafd74687179ab1ff8ff7be0bb6a13c745f067`. Data blocks not included. |
| `healpy_map_nside16_nest.fits`, `healpy_map_nside16_ring.fits`, `healpy_map_nside16.json` | generated by `gen_map.py` with healpy.write_map / read_map / ang2pix (TFORM 1024D vector rows, the healpy layout) |

## Storage

`gen_golden.py` writes a golden larger than 16384 bytes as `<name>.gz`
(Python `gzip.compress`, level 9, mtime 0) and every other one as plain
JSON. The Go tests read both through `testutil.ReadGolden`, and
`TestHealpixGoldensReproduce` compares the regenerated files with the
committed ones after decompression, byte for byte (the compressed bytes
depend on the zlib build). The JSON inside is unchanged.
