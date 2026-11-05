# wcslib reference goldens: provenance

Everything under `golden/` is produced by C programs in this directory
linked against **wcslib 8.10** (Arch package `wcslib 8.10-1`,
`/usr/include/wcslib`, `-lwcs`). No golden value is produced by the Go code
under test.

Generation environment for the goldens committed on branch
`proof/harness-wcslib`: gcc 16.2.1 (`-O2 -std=c99`), glibc 2.44, x86_64
Linux.

| Path | Generator | What it holds |
|---|---|---|
| `golden/<case>.json` + `golden/<case>.bin` | `wref.c` (`make wref-goldens`) | bare projections: `prjset` + `prjx2s` / `prjs2x` per projection code and PV set; `.json` is the metadata, `.bin` every input and wcslib output as exact IEEE-754 bits (big-endian, zlib level 9; layout in the header of `wref.c`) |
| `golden/prjset_status.json` | `wref.c` | `prjset` status of parameter sets wcslib refuses or accepts at a boundary |
| `golden/hdr/<name>.json.gz` (gzip -9n of the hdr2sky output) | `hdr2sky.c` (`make hdr-goldens`) from `hdrcases/<name>.case` | full header pipeline: `wcspih` -> `wcsset` -> `wcsp2s` / `wcss2p`, per-point status |

`make clean` removes only the binaries (`wref`, `hdr2sky`). Goldens are
committed truth and are overwritten in place by the `*-goldens` targets;
review a regeneration with `git diff golden/`.

Bare-projection golden format. Until round 2 of the proof audit each case
was one JSON file with a decimal row per point (61 cases, 5.4 MB, every
row repeating its input). They are now a small metadata `.json` and a
binary `.bin` per case (815 KB in all with the round-2 band probes,
1034 KB with the round-2b pole band probes, 1.15 MB with round 2c's
face-centre, seam and band probes and the azp_mu1_gamma45 case): the
same doubles bit for bit (checked value by value against the JSON
goldens when the format changed: 696580 values identical), inputs of the
round-trip rows not stored twice (the reader rebuilds them from the
forward rows), zlib on top. zlib 1.3.2 (Arch `zlib 1:1.3.2-3`) wrote the
committed payloads; its version is recorded in each `.json`.
`TestWrefGoldensReproduce` (and `make verify`) require the metadata and the
decoded payload byte for byte, and the compressed bytes too when the local
zlib is the recorded one (another zlib may deflate the same payload
differently). A corrupt payload fails zlib's Adler-32 check.

## Header-pipeline goldens (`hdr2sky`)

Each `hdrcases/<name>.case` is raw FITS header card text plus pixel and sky
points (format documented at the top of `hdr2sky.c`). The golden (format
`hdr2sky-2`) holds the wcslib version string, sha256 digests of the case
text, of the 80-column cards handed to wcspih and of the input
coordinates, what wcsset derived, and for every point wcslib's return
status, per-coordinate `stat` and output coordinates as exact big-endian
binary64 (base64-packed records), plus a readable `%.17g` sample of each
direction (`pickSample` in `hdr2sky.c`: first and last point, three
strided interior points, the first two points wcslib rejected, and the
accepted points holding the minimum and maximum of each output
coordinate), so a failure can be read without decoding. The case text, the cards and the
input points are not repeated: `hdr_golden_test.go` (package `transform`)
rebuilds them from the case file exactly as `hdr2sky` does, and the
digests prove the rebuild is bit-identical. It then parses the cards with
`header.ParseCards -> wcs.ParseAlt -> transform.New` and requires both
directions to agree, including which points are rejected. A missing
golden, a digest mismatch (a golden made from a different case file, or a
rebuild that differs from what wcslib was given), a sample that disagrees
with the packed records, an orphan golden, or a golden from a different
wcslib version is a test failure.

Size (round 2): the previous pretty-printed `%.17g` JSON, which also
embedded the case text and the cards, was 9.70 MB for the 192 cases; the
packed form was 1.53 MB (1.83 MB once the readable sample grew from the
first/last point to the `pickSample` set; 1.94 MB for the 208 cases of
round 3, largest file 48 kB) and carries the same values bit
for bit (checked
when it was introduced by decoding every new golden against the old one:
statuses, stats and every coordinate identical, inputs digest equal to the
digest of the old golden's input coordinates).

Reproducibility: `TestHdrGoldensReproduce` (and `make verify`) rebuilds
`hdr2sky` with `make` and requires byte-identical output for every case
(tier 2, `go test -tags interop`: a failure without a C compiler +
wcslib). Byte equality is required because the
coordinates are exact binary64; a failure that is only in the last bit or
two means a different compiler/libm/wcslib build, not a different answer,
and must be investigated rather than regenerated blindly.

Known caveat: SIP headers without AP/BP. wcslib 8.10 cannot read them, so
`hdr2sky` injects identity `AP_ORDER=1, AP_0_0=0, BP_ORDER=1, BP_0_0=0`
into wcslib's copy only (listed as `harness_injected` in the golden; the Go
side parses the original). wcslib inverts SIP by iterating on A/B with AP/BP
only as the starting guess, so the truth is the iterated inverse of A/B.
See the comment block at the top of `hdr2sky.c`.

### Cases

| Case | Source | License | Covers |
|---|---|---|---|
| `tan_cd_rot` | synthetic, hand-written (from the WCS audit probe `h/tan_cd_rot.hdr`) | n/a | TAN, CD matrix with 10 deg rotation, default LONPOLE; grid 256 px past every edge; sky points 90 and 180 deg from the tangent point (must be rejected) |
| `tan_sip3_apbp` | header cards of `astropy/wcs/tests/data/sip.fits` (astropy 8.0.2, file sha256 `eaf6be7d6a8748004106828326b5feb62043f6149a6d3a6cb20c2329c9f391d3`), primary HDU, copied card for card | BSD-3-Clause (astropy) | TAN-SIP, A/B order 3 with AP/BP order 3, PC + CDELT, explicit LONPOLE/LATPOLE |
| `car_allsky` | synthetic, hand-written | n/a | CAR all-sky 720x360, CRVAL2 = 0, no LONPOLE card (wcslib default), grid 20 px past every edge (off-map pixels rejected by wcslib) |

| `prj_*` (59 cases) | synthetic, written by `gen_prj_cases.py` (stdlib Python; re-running it reproduces the committed case files byte for byte) | n/a | every projection code at CRVAL2 = 0 and at CRVAL2 != theta_0 (most with non-default PVi_m: AZP mu/gamma, SZP, slant SIN, ZPN deg 3/5/12 (`prj_zpn_c0_deg12`: PV2_8, PV2_10, PV2_12), AIR theta_b, CYP Gall, CEA lambda, southern COP/COE/COO with eta, COD eta, BON theta_1 > 0 and < 0, HPX (3,2) and (6,1)), 1 deg/px, pixel grid spanning +-250 deg of intermediate coordinates (quad-cube -81..417) so every bounded projection is crossed, 13x13 sky grid; no LONPOLE card |
| `hiorder_*` (4 cases) | synthetic, written by `gen_highorder_cases.py` (stdlib Python; `TestHighOrderCaseGeneratorReproduce` re-runs it) | n/a | TAN-SIP order 7 and 9 with every A/B term 2 <= p+q <= order set (no AP/BP), TPV with every PVi_0..PVi_39 non-zero on both axes in CD and PC + CDELT form; coefficient formulas in each case file; each term moves the chip edge visibly; TOL_PX 1e-8 (measured <= 4.6e-10 px) |
| `pvfid_*` (3 cases) | synthetic, written by `gen_prj_cases.py` (`FIDUCIAL` list) | n/a | native fiducial point on the longitude axis: CAR with PV1_2 = theta_0 = 20, CAR with PV1_0 = 1, PV1_1 = 10, PV1_2 = 20, TAN with PV1_0 = 1, PV1_2 = 70 - agree with wcslib in both directions |
| `sample_1904-66_*` (28 cases) | header cards of Calabretta's Paper II example maps, `astropy/wcs/tests/data/maps/1904-66_<CODE>.hdr`, astropy 8.0.2, written by `gen_sample_cases.py` | BSD-3-Clause (astropy) | every Paper II projection incl. AIR, CSC and NCP-as-SIN, CRVAL2 = -90 with explicit LONPOLE/LATPOLE, real PV sets (e.g. ZPN degree 7) |
| `sample_msx_glon_car_lonpole0` | `astropy/visualization/wcsaxes/tests/data/msx_header` (astropy 8.0.2, sha256 `a9cccf6d...445b4c`), NAXIS1/2 = 150 added | BSD-3-Clause (astropy) | MSX survey GLON-CAR, LONPOLE = 0, LATPOLE = 90 |
| `sample_rosat_glon_ait_lonpole0` | `astropy/visualization/wcsaxes/tests/data/rosat_header` (astropy 8.0.2, sha256 `4b693fb7...6976c1`), NAXIS1/2 = 480/240 added | BSD-3-Clause (astropy) | ROSAT all-sky GLON-AIT, LONPOLE = 0; grid past the ellipse |
| `sample_hpxcvt_*` (5 cases) | `/usr/bin/HPXcvt` (wcslib 8.10) run by `gen_sample_cases.py` on an NSIDE = 8 HEALPix map `arange(768)` written by healpy 1.20.1 (RING galactic and NESTED equatorial; options none, `-q1`, `-xn`, `-xs`) | n/a (generated) | real HPXcvt HPX headers (PC rotated 45 deg, CDELT to 8 digits) and XPH headers (CRVAL2 = +-90, LONPOLE = 180); `hpxcvt_hpx_gal` grids every pixel centre |

`gen_sample_cases.py` must run with the venv interpreter (astropy 8.0.2,
healpy 1.20.1, numpy 2.4.6) and HPXcvt on PATH. Header cards are copied
verbatim except COMMENT, HISTORY, blank and DATE cards (no WCS content;
DATE and HPXcvt's HISTORY with the temp file name would make re-runs
differ). Both generators take an optional output root;
`TestPrjCaseGeneratorReproduce` and `TestSampleGeneratorsReproduce`
(tier 2: a failure without python3 / the venv + HPXcvt) re-run them into a temp directory and require
exactly the committed `prj_*`/`pvfid_*`/`sample_*` case files and both
`samples/` files back, byte for byte.

Every case, golden and sample the projection proof uses is also listed by
name in `prj_manifest_test.go`; `TestProjectionProofManifest` (no
external tools) fails on any listed file that is missing, any unlisted
file, a case without its wcslib golden, and any projection code without a
real header (`sample_1904-66_<code>`, or `sample_hpxcvt_xph_{n,s}_gal`
for XPH). The Calabretta and survey cases carry `TOL_PX 1e-7`: their sky
grid reaches points that land 1e5..3e7 px off the image, where the default
1e-9 px is below one ULP.

`samples/hpxcvt_hpx_gal_n8.fits` (sha256
`41040563a07e5ee644401e9ffb5039174b936a512bec08092de8c856d6b5654f`) is the
HPXcvt output image itself (40 x 40, each pixel = the RING index HPXcvt put
there); `samples/hpxcvt_hpx_gal_n8.healpy.json` (sha256
`6deb75cab475dcde11ae105adf8378b9e295b802462f1262a3f05c2e89ad796e`) is
healpy 1.20.1 `pix2ang(8, index)` for every filled pixel, the truth for
`TestHPXcvtPixelCentres`. HPXcvt writes the creation time into DATE and
the input file's (temp) path into a HISTORY card; `gen_sample_cases.py`
pins both (`PINNED_CARDS`: DATE = '2026-11-01T07:52:52', HISTORY
"Original input file: ring_gal.fits") so the image reproduces byte for
byte; no WCS card or pixel value is touched.

Observation, not a bug of this package: the XPH image HPXcvt writes with
`-xn` does NOT agree with healpy when read through its own header (via
astropy/wcslib, 736 of 768 pixel centres are 10-180 deg from `pix2ang` of
the stored index), while its HPX images agree to 2.6e-6 deg. The data
array is rotated by 180 deg relative to the header: with the array
reversed in both axes all 768 agree to 8.5e-7 deg (CDELT is written to
8 digits). wcslib's XPH itself agrees with the equations of the XPH
paper (see `paper/xph_paper.json` below), so the HPXcvt `-xn` header (or
its pixel layout) is the inconsistent part. The XPH sample cases are
therefore used as wcslib header truth only (our pipeline agrees with
wcslib on them), not as data-level truth; XPH's independent truth is the
paper.

`prjswap_*` (4 cases, written by `gen_prj_cases.py`, `SWAPPED` list):
latitude-first headers (`CTYPE1 = 'DEC--xxx'`, `CTYPE2 = 'RA---xxx'`,
projection parameters as PV1_m) for AZP mu/gamma, ZPN degree 3, CEA
lambda and a southern COP. wcslib accepts them, and so does
transform.New; both directions agree with wcslib.

#### Pipeline cases (WCS pipeline proof work, branch `proof/wcs-pipeline`)

Each targets one header rule; the comment block at the top of the case
file says what it pins and why. `TestHdrDerivedParameters` additionally
compares the WCS dimension, celestial axis indices and the derived
LONPOLE/LATPOLE (`cel_ref`) with wcslib for every case.

| Case | Source | License | Covers |
|---|---|---|---|
| `alt_a_cd_lonpole` | synthetic, hand-written | n/a | alternate WCS A, CDi_jA and LONPOLEA on CAR |
| `alt_a_celestial` | synthetic, hand-written | n/a | alternate WCS A with LONPOLEa, PVi_ma, WCSAXESa |
| `alt_a_no_wcsaxesa` | synthetic, hand-written | n/a | WCSAXES of the primary does not apply to alternate A |
| `alt_a_no_leak` | synthetic, hand-written | n/a | alternate WCS A whose answer depends on LONPOLEA/LATPOLEA (two-solution CAR); the primary carries different LONPOLE, LATPOLE, PCi_j, CDi_j and PVi_m that must not be read into A |
| `axes_cube_freq_ra_dec` | synthetic, hand-written | n/a | celestial axes at 2,3 behind a FREQ axis |
| `axes_index_above_naxis` | synthetic, hand-written | n/a | WCS axis count from the highest keyword index |
| `axes_swapped` | synthetic, hand-written | n/a | swapped celestial axes (DEC on axis 1) |
| `axes_swapped_pc_cdelt` | synthetic, hand-written | n/a | swapped celestial axes in PC + unequal CDELT form (CDELTi indexed by world axis) |
| `axes_wcsaxes3_naxis2` | synthetic, hand-written | n/a | WCSAXES > NAXIS |
| `crota2_tan` | synthetic, hand-written | n/a | CROTA2 on TAN |
| `crota2_with_cd` | synthetic, hand-written | n/a | CD matrix and CROTA2 = 25 together: CD wins, CROTA2 ignored |
| `crota2_with_pc` | synthetic, hand-written | n/a | PC1_2 and CROTA2 = 25 together: PC wins, CROTA2 ignored |
| `crota2_unequal_cdelt` | synthetic, hand-written | n/a | CROTA2 with unequal CDELT, CROTA1 ignored |
| `cunit_arcsec` | synthetic, hand-written | n/a | CUNIT arcsec with PC + CDELT |
| `cunit_arcmin` | synthetic, hand-written | n/a | CUNIT arcmin |
| `cunit_arcsec_cd` | synthetic, hand-written | n/a | CUNIT arcsec with CD |
| `cunit_mas` | synthetic, hand-written | n/a | CUNIT mas |
| `cunit_rad` | synthetic, hand-written | n/a | CUNIT rad |
| `frame_elon_tan` | synthetic, hand-written | n/a | ELON/ELAT |
| `frame_glon_car_lp0` | synthetic, hand-written | n/a | GLON-CAR, LONPOLE = 0 |
| `latpole_car_two_solutions` | synthetic, hand-written | n/a | LATPOLE choosing between two native-pole solutions |
| `latpole_coe_two_solutions` | synthetic, hand-written | n/a | conic with LONPOLE and LATPOLE both explicit |
| `lonpole_ait_gal_offset` | synthetic, hand-written | n/a | AIT with CRVAL2 = 20, default LONPOLE/LATPOLE |
| `lonpole_car_dec30_lp0` | synthetic, hand-written | n/a | LONPOLE = 0 on a cylindrical projection, CRVAL2 = 30 |
| `lonpole_car_decm30_default` | synthetic, hand-written | n/a | default LONPOLE on CAR with CRVAL2 < theta_0 |
| `lonpole_cea_lp90` | synthetic, hand-written | n/a | LONPOLE = 90 on CEA, pole latitude from LATPOLE only |
| `lonpole_glon_ait_lp0` | synthetic, hand-written | n/a | GLON-AIT all-sky with LONPOLE = 0 |
| `lonpole_mol_dec20_lp0` | synthetic, hand-written | n/a | MOL with LONPOLE = 0 |
| `lonpole_tan_lp150` | synthetic, hand-written | n/a | zenithal TAN with non-default LONPOLE |
| `lonpole_tan_dec40_lp0` | synthetic, hand-written | n/a | TAN at CRVAL2 = 40 with an explicit LONPOLE = 0, where wcslib's default is 180 (the only accepted case in which LONPOLE = 0 is not the default) |
| `lonpole_tan_southpole` | synthetic, hand-written | n/a | zenithal TAN at CRVAL2 = -89.5, default LONPOLE |
| `matrix_cd_partial` | synthetic, hand-written | n/a | partial CD matrix (missing off-diagonal CD1_2) |
| `matrix_cd_diag_missing` | synthetic, hand-written | n/a | CD matrix with a missing diagonal element (CD2_2 defaults to 0, not 1) |
| `matrix_pc_and_cd` | synthetic, hand-written | n/a | PC and CD both present, PC wins |
| `pvlon_car_offset` | synthetic, hand-written | n/a | PVi_0..2 on the longitude axis, fiducial offset applied |
| `pv1off_car_pv1_1_only` | synthetic, hand-written | n/a | PV1_0 = 1 with PV1_1 but no PV1_2: no fiducial offset (wcslib prjoff needs both) while PV1_1 still sets phi_0 for celset |
| `pv1off_car_pv1_2_only` | synthetic, hand-written | n/a | PV1_0 = 1 with PV1_2 but no PV1_1: no fiducial offset, PV1_2 still sets theta_0 |
| `pvlon_car_pv3_pv4` | synthetic, hand-written | n/a | PVi_3 / PVi_4 on the longitude axis override LONPOLE / LATPOLE |
| `pvlon_ait_phi0_offset` | synthetic, hand-written | n/a | PVi_0..2 offset on AIT, where phi_0 (PV1_1) changes the answer (CAR is translation-invariant in phi) |
| `pvlon_tan_phi0` | synthetic, hand-written | n/a | PVi_1 (phi_0) on TAN without offset |
| `pvlon_tan_theta0` | synthetic, hand-written | n/a | PVi_2 (theta_0) on the longitude axis without offset |
| `sip4_no_apbp` | synthetic, hand-written | n/a | SIP order 4 without AP/BP (Newton inverse) |
| `sip4_rough_apbp` | synthetic, hand-written | n/a | SIP AP/BP as starting guess only |
| `sip_pc_cdelt` | synthetic, hand-written | n/a | SIP with PC + CDELT |
| `sip_suffix_no_coeffs` | synthetic, hand-written | n/a | -SIP CTYPE without coefficients is undistorted |
| `tpv_explicit_pv1` | synthetic, hand-written | n/a | TPV with explicit PV1_1 = PV2_1 = 1 (registry default) |
| `tpv_pc_cdelt` | synthetic, hand-written | n/a | TPV with PC + CDELT, polynomial applied before CDELT (wcslib) |
| `tpv_radial_cd` | synthetic, hand-written | n/a | TPV radial terms r, r^3, r^5, r^7 on both axes (CD form) |
| `alias_ncp_c60` | synthetic, hand-written | n/a | AIPS RA---NCP/DEC--NCP at CRVAL2 = 60 (wcsset: SIN with PV2_2 = cot CRVAL2) |
| `alias_ncp_cm25_pv` | synthetic, hand-written | n/a | NCP at CRVAL2 = -25 with stray PV2_1/PV2_2, which wcsset overwrites |
| `alias_gls_c30` | synthetic, hand-written | n/a | AIPS RA---GLS (wcsset: SFL with offset, theta0 = CRVAL2) |
| `legacy_pc00i00j_rot` | synthetic, from the proof-audit probe `zz_pc00` | n/a | defunct PC00i00j rotation matrix + CDELT (wcspih WCSHDR_PC00i00j) |
| `legacy_cd00i00j` | synthetic, from probe `zz_cd00` | n/a | defunct CD00i00j matrix, no CDELT |
| `legacy_cd00i00j_cdelt` | synthetic, from probe `zz_crota_cd00` | n/a | defunct CD00i00j diagonal with CDELT present (CDELT ignored) |
| `legacy_projp_azp` | synthetic, from probe `zz_projp` | n/a | defunct PROJP1/PROJP2 on AZP (PV2_1 = mu, PV2_2 = gamma) |
| `legacy_projp_zpn` | synthetic, from probe `zz_projp_pv` | n/a | defunct PROJP1/PROJP3 on ZPN |
| `frame_slon_tan` | synthetic, from probe `zz_slon` | n/a | supergalactic SLON/SLAT with TAN |
| `frame_hlon_zea` | synthetic, from probe `zz_hlon` | n/a | helioecliptic HLON/HLAT with ZEA |
| `sip_on_sin` | synthetic, from probe `zz_sinsip` | n/a | `-SIP` on SIN (no AP/BP; hdr2sky injects for wcslib) |
| `sip_nosuffix_tan`, `sip_nosuffix_sin` | synthetic, hand-written (from the round-2 verifier probes) | n/a | SIP keywords on a plain `RA---TAN` / `RA---SIN` header, no `-SIP`: wcslib applies SIP (wcspih sets the SIP prior distortion from the A_p_q/B_p_q keywords whatever CTYPE says), ~8 arcsec at the corners |
| `tpv_sip_suffix` | synthetic, hand-written | n/a | `RA---TPV-SIP` (the suffix written out) with TPV and SIP keywords: wcslib applies SIP (prior) and TPV (sequent), and so do we; the same header without the suffix is `hdrcases_differ/tpv_sip_nosuffix` (below) |
| `tpv_sip_nosuffix_tpvonly` | synthetic, `hdrcases_differ/tpv_sip_nosuffix` without its SIP cards | n/a | the answer that `RA---TPV` header must give: its SIP set is ignored, as astropy 8.0.2 ignores it |
| `tpv_sip_suffix_atlas` | synthetic: the TPV cards of `real_atlas_k2_tpv` (real) with `RA---TPV-SIP` and a non-redundant SIP set (A/B order 3, no AP/BP; hdr2sky injects identity AP/BP) | n/a | wcslib 8.10 accepts `RA---TPV-SIP` and applies both; read back through the TPV solution its answer is the pixel plus the SIP vector (2.33-3.0 px) to 1.1e-10 px (`TestTPVSIPSuffixAppliesBoth`) |
| `sip_nosuffix_alt_a`, `sip_nosuffix_primary_of_alt` | synthetic, hand-written | n/a | the same header read as alternate A (undistorted: SIP has no alternate form, wcslib attaches it to the primary only) and as the primary (SIP applied) |
| `cunit_mixed_deg_arcsec`, `cunit_mixed_arcmin_deg` | synthetic, hand-written | n/a | a different CUNIT on each celestial axis (each axis converted with its own unit) |
| `pvlon_car_offset_phi0_only`, `pvlon_car_offset_theta0_only` | synthetic, hand-written | n/a | PV1_0 = 1 with only one of PV1_1 / PV1_2: no projection offset (wcslib prjoff), the given phi_0 / theta_0 still used for the native pole (celset) |
| `pvlon_car_theta0_tol` | synthetic, hand-written | n/a | PV1_2 = 90 + 5e-11, inside celset's 1e-10 tolerance: clamped to 90 |
| `celpole_car_eq_south`, `celpole_car_eq_latpole_tie` | synthetic, hand-written | n/a | CRVAL2 = 0, LONPOLE = 0: native-pole solutions +-90; LATPOLE = -90 and LATPOLE = 0 (an exact tie, which wcslib gives to latp2) both put the celestial south pole at the native pole (celset's z = 0, latp < 0 branch) |
| `sip_swapped_nosuffix`, `sip_swapped_suffix` | synthetic, from the round-2 verifier probes `zz_sip_swapped_*` | n/a | SIP with DEC--TAN on pixel axis 1, asymmetric A/B: A stays on pixel axis 1 (a swap mutant is 2.48 arcsec off) |
| `sip_cube_radec_freq` | synthetic, hand-written | n/a | SIP order 3 on a RA, DEC, FREQ cube (no AP/BP; hdr2sky injects) |
| `celpole_z_above_tol`, `celpole_z_within_tol`, `celpole_latp1_tol_out`, `celpole_latp2_tol_out`, `celpole_slz_within_tol` | synthetic, hand-written (round 2) | n/a | one case on each side of every celset tolerance (z = cos latp cos lat0 at 1.75e-10 and 1.6e-11; nearer native-pole latitude 1e-9 past +90 / -90 passed over; slz = 1 + 5e-11 clamped); the reject twins are below |
| `celpole_car_native_pole_sky` | synthetic, hand-written (round 2) | n/a | sky points exactly at and 1e-7 deg from both native poles of CAR, where the pixel follows wcslib sphs2x's longitude-origin convention and its rearranged formula |
| `dss_modern_cutout` | synthetic, from the round-2 verifier probe `zz_dss_probe` | n/a | DSS plate keywords over an approximate TAN CD WCS: wcslib replaces the approximate WCS (74 arcsec apart) |
| `dss_all_terms_deg5` | synthetic, hand-written | n/a | DSS plate solution with every AMDX/AMDY 1..13 non-zero (wcslib's degree-5 TPD form), southern sign, grid to +-5000 px (at +-6000 px wcslib's disx2p fails to converge at three corners: `dss_far_outside_plate`) |
| `dss_far_outside_plate` | synthetic, the plate solution of `dss_all_terms_deg5` | n/a | grid to +-6000 px: wcslib's p2s accepts every point but its s2p (disx2p) fails to converge (status 9) on three far corners; there `hdr_golden_test.go` (`hdrS2PBeyondWcslib`) requires our s2p to return the source pixel of wcslib's p2s (wcslib's forward direction is the truth), and every other point agrees with wcslib (`TOL_PX 3e-9`: wcslib's own round trip is 1.3e-9 px off there) |
| `celpole_z0_latpole95`, `celpole_z0_latpolem95` | synthetic, from the round-3 verifier probes `probecel/d_*`, `e_*` | n/a | CAR, CRVAL2 = 0, LONPOLE = 90: theta0 = 0 and \|phip - phi0\| = 90 put celset in its z = 0 branch (native-pole latitude from LATPOLE alone); LATPOLE = 95 / -95 are clamped to +90 / -90 by wcslib, which accepts both; the reject twin is `celpole_z0_lat30` |
| `lonpole_999_default` | synthetic, from the round-3 verifier probe `probecel/b_lonpole999` | n/a | LONPOLE = 999 is wcslib's "use the default" value (celset); read literally it would leave no valid native pole |

Cases whose sky -> pixel direction uses an iterative inverse at a fine
plate scale set `TOL_PX 1e-8` and say why in the case file: the s2p input
is wcslib's p2s sky position, whose angular round-off (~2e-10 arcsec) is
several 1e-9 px at 0.03-0.05 arcsec/px (wcslib's own round trip on those
cases is off by up to 4.8e-9 px).

#### Real archive headers (header cards only)

Generated by `realhdr/make_cases.py` (standard-library Python; downloads
each source, checks the pinned sha256, copies the HDU's cards verbatim,
removes only the cards listed in the case file, and replaces only the
cards listed there: `real_atlas_k2_tansip`'s two CTYPEs) and never edited
by hand (`real_atlas_k2_tpv_sip` and `real_ptf_tpv_plus_sip` are written
to `hdrcases_differ/`, below):
each case records the sha256 of its card block and
`TestRealHdrCasesIntact` checks it (this catches an accidental edit of
the cards, not a deliberate edit of cards and hash together; the pin
against the archive bytes is the source sha256 in the script). Re-running
the script (with network, or `--src DIR` holding the downloaded files)
reproduces the case files byte for byte; this was checked by hand for
every case (round 2: all 17 sources rerun, `--src` holding astropy 8.0.2's
test data and the downloads; `real_ptf_tpv_plus_sip` had carried four
hand-added TOL_PX lines and was regenerated), and is not run by a test,
because it needs the downloads. `gzip=True` sources (the DSS file) are
pinned by the sha256 of the compressed file. The
ROSAT and MSX sources are newline-separated header text (astropy test
data style), read with `text_header=True`. The ZTF science image is 37 MB and the JWST/HST products
are fetched with an HTTP Range request (first 300-400 kB, enough for the
headers); only header cards are committed.

| Case | Source URL | License | Source sha256 | What |
|---|---|---|---|---|
| `real_anet_sip3` | https://nova.astrometry.net/wcs_file/11126704 | nova.astrometry.net public job output (WCS solution only; no image data) | `87fb0ac3da96cfc424014caae15580a4fd8bd0e4bbcc74c13eaadc44b92a13b2, HDU 0` | astrometry.net (nova) job 11126704 wcs.fits: TAN-SIP A/B order 3 with AP/BP order 3, 2044x2044 image |
| `real_anet_sip3_wide` | https://nova.astrometry.net/wcs_file/11522654 | nova.astrometry.net public job output (WCS solution only; no image data) | `8c63705e10e3d7a59514a114970c30681faf04ce08bf86990f093dc1502b09ad, HDU 0` | astrometry.net (nova) job 11522654 wcs.fits: TAN-SIP order 3 with AP/BP order 3, 3920x1016 image at 2.5 arcsec/px (2.7 deg wide) |
| `real_anet_sip4` | https://nova.astrometry.net/wcs_file/11087109 | nova.astrometry.net public job output (WCS solution only; no image data) | `e148bc064534badece6b3d84e213c24d9b318a7964ace29f0ea3368856638671, HDU 0` | astrometry.net (nova) job 11087109 wcs.fits: TAN-SIP A/B order 4 with AP/BP order 4, 3192x2128 image |
| `real_hst_acs_wfc_sip4` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/j94f05bgq_flt.fits | astropy test data, BSD-3-Clause (header of a public HST/MAST product) | `900038e0d853828140a757e2656934cb268ff9f315c5c6f617de85a632ad526b, HDU 1` | HST ACS/WFC FLT j94f05bgq, SCI,1: TAN-SIP A/B order 4, no AP/BP, CD matrix; also carries an alternate WCS 'O' (OPUS) the primary parse must ignore |
| `real_hst_wfc3_uvis_sip4` | https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/ic3t33svq_flt.fits | public HST data from MAST (STScI); header cards only | `73994d4498b679d4c0140bbee70b3b59d87fc49004396ad181fbee79f478a380, HDU 1` | HST WFC3/UVIS FLT ic3t33svq, SCI,1: TAN-SIP A/B order 4, no AP/BP; lookup-table (CPDIS/D2IM) records removed |
| `real_jwst_nircam_sip4` | https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/jw02736001001_02105_00001_nrca1_cal.fits | public JWST data from MAST (STScI); header cards only | `6c4ddd65ab390ae070a7ddb4aed4fe2f26e7438c90b9bad99801ec72a0a3eafc, HDU 1` | JWST NIRCam NRCA1 cal jw02736001001_02105_00001 (ERO SMACS 0723), SCI: TAN-SIP A/B order 4 with approximate AP/BP order 4 |
| `real_jwst_nircam_sip5` | https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/jw01345001001_02201_00001_nrca1_cal.fits | public JWST data from MAST (STScI); header cards only | `a69f669dda4da9cb16c8afe9c98fc0edae1a909c78e43afe391c5ea16b89ae2f, HDU 1` | JWST NIRCam NRCA1 cal jw01345001001_02201_00001 (CEERS), SCI: TAN-SIP A/B order 5 with approximate AP/BP order 5 (the jwst pipeline's SIP approximation), PC+CDELT |
| `real_ptf_tpv4` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/tpvonly.hdr | astropy test data, BSD-3-Clause (header of a public PTF/IRSA product) | `ca3f753c704fadc7afb74b2c085afd4539c494271afbd212c93476a92025899d, HDU 0` | PTF/MOSAIC CCD 5 (2014-07-31): TPV degree 4 with CRPIX ~4000 px off the chip, CD matrix |
| `real_ptf_tpv_plus_sip` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/tpvonly.hdr plus the `(A\|B\|AP\|BP)_(ORDER\|DMAX\|p_q)` cards of https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/siponly.hdr | astropy test data, BSD-3-Clause (headers of a public PTF/IRSA product) | `ca3f753c704fadc7afb74b2c085afd4539c494271afbd212c93476a92025899d` + `ac786a124934a3cd00f11921c6f2647b470f05ef04a2f4033960b408635a34a9`, HDU 0 | the TPV header of real_ptf_tpv4 with the SIP order-4 A/B/AP/BP solution of the same exposure appended, no `-SIP` suffix (the TPV + redundant SIP form), in `hdrcases_differ/`: we and astropy 8.0.2 ignore the SIP set (the answer is `real_ptf_tpv4`'s), wcslib applies both, 72.7-133.4 arcsec from TPV alone over the chip grid |
| `real_atlas_k2_tpv_sip` | https://github.com/user-attachments/files/34297216/atlas_k2_tpv_sip_wcs.txt (attachment of Daryl Manning's comment on dmarkham/fits PR #4, 2026-11-03) | WCS cards only (no image data) of an ATLAS survey product, posted by Daryl Manning as a test sample | `3b537ea0fa4cb12a1f41539bef0abfc02ecfe022fbfd64534600c031cb8a7d4b` (134 80-column cards, newline-separated, no END) | ATLAS difference stamp 01a59742o0240o.diff (C/2017 K2, Mauna Loa unit, MJD 59742.348, 100x100 cutout): `RA---TPV` (no `-SIP`), CD matrix, CRPIX ~4000 px off the stamp, TPV PV1_0..22 / PV2_0..22 PLUS an order-5 A/B/AP/BP SIP set that is a copy of the TPV correction; the cards verbatim, in `hdrcases_differ/` (we and astropy 8.0.2 ignore the SIP set, wcslib applies both) |
| `real_atlas_k2_tpv` | the same file | as above | as above | the ATLAS header with its SIP cards removed (TPV alone): the answer of `real_atlas_k2_tpv_sip` |
| `real_atlas_k2_tansip` | the same file | as above | as above | the ATLAS header with its PV cards removed and CTYPE1/2 replaced by `RA---TAN-SIP` / `DEC--TAN-SIP` (the SIP set alone, CRPIX/CRVAL/CD unchanged) |
| `real_dss_plate_solution` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/dss.14.29.56-62.41.05.fits.gz | astropy test data, BSD-3-Clause (header of a public STScI DSS cutout) | `d9e8a9adced55200b418517d5fba5f61a403842aa5eb57df9281eb38beb8d68a` (gzip), HDU 0 | DSS plate solution (UK Schmidt S134 plate 0025) with the approximate TAN WCS (CD, PC00i00j, CROTA) that wcslib replaces; grid over the cutout and +-2000 px around it |
| `real_rosat_glon_ait_lp0` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/visualization/wcsaxes/tests/data/rosat_header | astropy test data, BSD-3-Clause (header of a public ROSAT survey map) | `4b693fb7d43dc7ee78f5031bda925de219715b3c00bc2ff1b300c4033c6976c1, HDU 0` | ROSAT all-sky survey map: GLON-AIT / GLAT-AIT, LONPOLE = 0, LATPOLE = 90, CDELT 0.675 deg; grid runs off the all-sky ellipse (off-ellipse pixels are rejected, as wcslib does) |
| `real_msx_glon_car_lp0` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/visualization/wcsaxes/tests/data/msx_header | astropy test data, BSD-3-Clause (header of a public MSX survey image) | `a9cccf6d7bf2365e24583bf6f10b7f7e75f4a731d35b9bade7ba6b6651445b4c, HDU 0` | MSX Galactic-plane image: GLON-CAR / GLAT-CAR, LONPOLE = 0, LATPOLE = 90 |
| `real_ps1_pc00i00j` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/defunct_keywords.hdr | astropy test data, BSD-3-Clause (header of a public Pan-STARRS1 product) | `30d1197ab65e77a2846f1117cf03143daa634488291a14e7ca6a3c9e435caf2b, HDU 0` (raw 80-column header, final END unpadded) | Pan-STARRS GPC1 skycell 0680.071 deep stack: RA---TAN with the defunct PC00i00j matrix (PC001001 = -1) and CDELT, CRPIX far off the image |
| `real_spitzer_irac_cryo_sip2` | https://raw.githubusercontent.com/astropy/astropy/v8.0.2/astropy/wcs/tests/data/irac_sip.hdr | astropy test data, BSD-3-Clause (header of a public Spitzer/IRSA product) | `e5ceeee30461dd50b2d674c8c923bf2e22ad3cb3f365b5cfee78a82afc7ddd2f, HDU 0` | Spitzer IRAC ch1 BCD (cryogenic mission, S8.9.0 pipeline): TAN-SIP A/B order 2 with AP/BP order 2, CD matrix |
| `real_spitzer_irac_warm_sip5` | https://irsa.ipac.caltech.edu/data/SPITZER/SHA/archive/proc/IRAC041600/r53849088/ch1/bcd/SPITZER_I1_53849088_0128_0000_1_cbcd.fits | public Spitzer Heritage Archive data (IRSA/IPAC); header cards only | `c171fe84e010d6a4c2e7c89dd777bdafcd13d7f63d21252c80b031121604d51d, HDU 0` | Spitzer IRAC ch1 CBCD (warm mission, AOR 53849088, 2016): TAN-SIP A/B order 5 with AP/BP order 5 |
| `real_ztf_tpv4` | https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/2020/0331/486227/ztf_20200331486227_000535_zr_c11_o_q3_sciimg.fits | public ZTF data (IRSA/IPAC); header cards only | `c516d7fe4292388352f366c02a91f5ae507056696e5b723e1095e2c07e52c89f, HDU 0` | ZTF sciimg field 535 zr c11 q3 (2020-03-31): SCAMP TPV, PV1_0..PV1_16 / PV2_0..PV2_16 (degree 4, no radial terms), CD matrix |

`real_hst_wfc3_uvis_sip4` has its CPDIS/D2IM lookup-table records removed
(the tables live in WCSDVARR/D2IMARR extensions); with them present
`transform.New` refuses the header (`TestNewRefusesUnimplementedDistortions`).

Spitzer IRAC: the warm-mission BCD found is order 5 (A/B and AP/BP), the
cryogenic one (astropy's `irac_sip.hdr`) order 2; no order-3 IRAC header
was located. ZTF science images use TPV (SCAMP), not SIP, so ZTF appears
as `real_ztf_tpv4`. No real TPV header with non-zero radial terms was
found (ZTF, PTF and the DECam CP header inspected write PVi_3/PVi_11 = 0
or omit them); radial terms are covered by the synthetic `tpv_radial_cd`.

TPV + redundant SIP (ATLAS): ATLAS pipeline headers carry a TPV solution
and, for SIP-only readers, a SIP copy of the same correction, with CTYPE
`RA---TPV` (no `-SIP`). Until PR #4 no such header was public and
`real_ptf_tpv_plus_sip` (the TPV and SIP solutions IRSA published for one
PTF exposure, merged the way astropy's own `test_tpv_ctype_tpv` merges
them) stood in for the form, read with both distortions as wcslib reads
it. Daryl Manning then posted the WCS cards of a real ATLAS stamp
(`real_atlas_k2_tpv_sip`), which shows that reading is wrong: wcslib on
its TPV part alone (`real_atlas_k2_tpv`) and on its SIP part alone
(`real_atlas_k2_tansip`) agree to 2.31e-7 arcsec over the stamp (pixel
offset (-12.28, 6.06) px at the corner (1, 1) from either), wcslib on the
header as written lands 24.45-26.84 arcsec (13.28-14.59 px) from that
answer at every point, and astropy 8.0.2, which drops the SIP set, gives
the TPV answer to 3.6e-10 arcsec, the positions that match the image
(Daryl Manning's report). A `RA---TPV` header without `-SIP` now ignores
its SIP keywords (`wcs.Header.SIPIgnored`); both real merged headers and
the synthetic `tpv_sip_nosuffix` are deliberate differences from wcslib
(`hdrcases_differ/`, `TestHdrDeliberateDifferences`,
`TestATLASTPVPlusRedundantSIP`, `TestATLASSIPCopiesTPV`).
`real_atlas_k2_tpv_sip`, `real_atlas_k2_tpv` and `real_atlas_k2_tansip`
were generated with `--src` holding the downloaded attachment (its sha256
checked against the pin, which the download reproduced on 2026-11-04);
`real_ptf_tpv_plus_sip` was regenerated into `hdrcases_differ/` from
astropy 8.0.2's `tpvonly.hdr` / `siponly.hdr` (its card block, and so its
`header cards sha256`, unchanged).

GLON with LONPOLE = 0 on real survey headers: ROSAT (AIT) and MSX (CAR).
Both have CRVAL2 = 0, where LONPOLE = 0 is also wcslib's default, so they
prove the real cards parse and agree, not that an explicit LONPOLE = 0 is
honoured. Neither do the synthetic `lonpole_car_dec30_lp0`,
`lonpole_glon_ait_lp0`, `lonpole_mol_dec20_lp0` and `frame_glon_car_lp0`:
all have CRVAL2 >= theta_0 = 0, where 0 is the default too (an earlier
version of this file called them non-default; a reader that ignores an
explicit LONPOLE = 0 passed all of them). For theta_0 = 0 (CAR, AIT, ...)
an explicit LONPOLE = 0 is the default whenever it is valid: with
phi_p = phi_0 the reference point needs cos(delta_p) = sin(delta_0), which
has no solution for delta_0 < 0. The rule is therefore pinned by
`lonpole_tan_dec40_lp0` (TAN, CRVAL2 = 40, default 180, wcslib accepts)
and by the GLON-CAR / GLON-AIT reject cases `lonpole_glon_car_bm30_lp0`
and `lonpole_glon_ait_bm20_lp0` (CRVAL2 < 0, wcslib refuses; a reader that
substitutes the default 180 would return coordinates).

#### Headers wcslib refuses (`hdrcases_reject/`)

Same case format; no golden. `TestHdrRejectCases` requires the built
`hdr2sky` to FAIL on each with the quoted wcslib error text (tier 1
replays hdr2sky's recorded reading of each case from
`testdata/recorded/wcs/transform`; tier 2 runs the built `hdr2sky` live
and compares) and `wcs.ParseAlt` / `transform.New` to refuse it too.
`TestHdrStricterThanWcslib` and `TestHdrLenientThanWcslib` read their
wcslib side the same way.

| Case | wcslib 8.10 error |
|---|---|
| `cunit_unknown` | Mismatched units type 'plane angle' |
| `cunit_upper` | Invalid symbol in EXPON context |
| `lonpole_car_dec30_lp180` | No valid solution for latp |
| `lonpole_glon_car_bm30_lp0` | Ill-conditioned coordinate transformation parameter (explicit LONPOLE = 0 at CRVAL2 = -30) |
| `lonpole_glon_ait_bm20_lp0` | Ill-conditioned coordinate transformation parameter (explicit LONPOLE = 0 at CRVAL2 = -20) |
| `proj_mismatch` | Inconsistent projection types |
| `two_lon_axes` | Inconsistent projection types |
| `pvlon_car_pv1_5` | Unrecognized coordinate transformation parameter |
| `pvlon_ait_theta0_95` | theta0 > 90 |
| `tpv_pv1_40` | Invalid parameter number (40) for TPV |
| `ncp_equator` | Invalid projection: NCP blows up on the equator |
| `tpv_pv2_40` | Invalid parameter number (40) for TPV (latitude axis) |
| `pvlon_axis2_pv5` | Unrecognized coordinate transformation parameter (PV2_5 on a longitude axis that is axis 2) |
| `pvlon_car_theta0_90p5` | theta0 > 90 (90.5 is beyond celset's 1e-10 tolerance) |
| `two_lat_axes` | Inconsistent projection types (two DEC axes) |
| `pvlon_car_theta0_tol_out` | theta0 > 90 (90 + 1.5e-10, just outside celset's 1e-10) |
| `celpole_slz_tol_out` | \|lat0\| <= 30.000 is required for these values of phip, phi0, and theta0 (slz = 1 + 2e-10) |
| `celpole_latp_tol_out` | No valid solution for latp (both solutions 1e-9 past +-90) |
| `dss_with_sip` | No distortion type, yet DPja keyvalues are present (DSS and SIP together; DQia with the SIP cards last) |
| `dss_zero_scale` | PCi_ja matrix is singular (DSS plate scale zero) |
| `celpole_z0_lat30` | lat0 == 0 is required for \|phip - phi0\| = 90 and theta0 == 0 (CAR, CRVAL2 = 30, LONPOLE = 90) |
| `sip_zero_coeff_deg10`, `sip_inverse_coeff_deg10` | Invalid powers (5, 5) for SIP on axis 1 / (9, 1) on axis 2: sipset refuses any SIP coefficient with p + q > 9, even 0.0, forward or inverse (round-3 verifier probe p21_zero_deg10) |
| `sip_ap_beyond_order` | Degree must not exceed 9, got -1 (BP has no coefficient card; we refuse the non-zero AP_2_0 above AP_ORDER = 1 instead; round-3 probe p40_ap_bp_beyond) |

#### Headers we refuse but wcslib accepts (`hdrcases_stricter/`)

Same case format; no golden. Deliberate policy (never silently drop or
guess at a declared distortion). `TestHdrStricterThanWcslib` requires
`wcs.ParseAlt` / `transform.New` to refuse each with the quoted
`# EXPECT-GO-ERROR:` text, and the built `hdr2sky` to ACCEPT the case as
written (no harness-injected cards, every pixel->sky point computed), which
proves the refusal is ours and not wcslib's. Each SIP case carries explicit
AP/BP cards: without them wcslib refuses every SIP header ("Degree must not
exceed 9, got -1", see the CAVEAT in `hdr2sky.c`), which says nothing about
A_ORDER/B_ORDER.

| Case | What wcslib 8.10 does | We refuse with |
|---|---|---|
| `sip_a_without_b` | reads the missing B polynomial as zero (p2s observed identical to the same header with `B_ORDER = 0`) | B_ORDER is missing |
| `sip_coeffs_without_orders` | applies `A_2_0` with no `A_ORDER` declared | A_ORDER is missing |
| `singular_cd` | wcsset succeeds and wcsp2s computes every point; wcss2p has no answer | linear transform is singular |
| `legacy_pc_conflict` | accepts PC1_1 and PC001001 with different values and keeps one | disagree |
| `legacy_projp_conflict` | accepts PROJP1 and PV2_1 with different values and keeps one | disagree |
| `sip_a_without_b_nosuffix` | as `sip_a_without_b`, on a plain `RA---TAN` header (wcslib applies SIP keywords without the suffix) | B_ORDER is missing |
| `sip_nosuffix_cube_axes23` | applies SIP (no axis mapping) to pixel axes 1 and 2 of a FREQ x RA x DEC cube, i.e. to FREQ and RA | SIP is defined on pixel axes 1,2 only |
| `tpv_cunit2_arcsec` | applies the TPV polynomial to a latitude coordinate in arcsec | TPV coefficients are defined in degrees |
| `couple_lng_row`, `couple_lat_row` | PC1_3 / PC2_3 mix the spectral pixel axis into RA / DEC and wcsp2s uses it | couples pixel axis 3 into celestial world axis 1 / 2 |
| `ps_numeric` | wcspih discards the numeric PS1_0 card and computes coordinates | PS1_0 has a non-string value |
| `sip_index_leading_zero` | the lexer discards `A_01_1` (single-digit indices only) and computes without that term | A_01_1 is not a SIP coefficient keyword |
| `dss_naxis3` | writes the 2x2 DSS matrix into pc[0..3] of a 3x3 matrix and computes | DSS plate solution on a 3-axis WCS |
| `dss_pltdecsn_unpadded` | reads `'-'` with sscanf as "-'" and puts the plate centre in the north | PLTDECSN = '-' without a trailing blank |
| `dss_amd_leading_zero` | discards `AMDX07` and computes without that term | AMDX07 is not a DSS coefficient keyword |
| `tnx_sip_nosuffix`, `tnx_sip_suffix` | an IRAF TNX header (WAT lngcor/latcor) with SIP keywords after the WAT cards, `RA---TNX` and `RA---TNX-SIP`: the lexer ends with distran = SIP, so wcslib applies SIP with TAN and ignores the WAT correction (SIP cards first: wcspih fails). From the round-3 verifier probes `probetnx/tnx_sip_*` | TNX header carries SIP keywords |
| `sip_beyond_order_nonzero` | A_3_0 = 1e-8 with A_ORDER = 2: sipset never reads A_ORDER (degree from the coefficients present) and applies A_3_0 | A_3_0 = 1e-08 has degree 3 > A_ORDER = 2 |
| `sip_order10`, `sip_neg_order`, `sip_order_float` | A_ORDER = 10, -1 and 2.0: never read, the header is computed from the coefficients present | A_ORDER = 10 is outside 0..9 / A_ORDER = -1 is outside 0..9 / A_ORDER = 2 is not an integer |
| `sip_coef_string` | A_0_2 = 'abc': the term is skipped (p2s records identical to the header without that card) | A_0_2 has a non-numeric value |
| `sip_mixed_suffix` | `RA---TAN-SIP` with `DEC--TAN`: applies SIP (the keywords alone make it the prior distortion) | inconsistent distortion codes |
| `sip_alt_suffix` | `CTYPE1A = 'RA---TAN-SIP'` read as description A: computed UNDISTORTED (p2s records identical to `hdrcases/sip_nosuffix_alt_a`), the declared -SIP silently not applied | SIP (A_ORDER/B_ORDER have no alternate-description form; CTYPE1A declares -SIP) |
| `sip_unknown_dist_code` | `RA---TAN-XYZ` with SIP keywords: the unknown code is ignored and SIP applied (p2s records identical to the `RA---TAN` header; round-3 probe p41_sip_unknown_dist) | unsupported distortion "XYZ" in CTYPE |
| `tpu_alias` | `RA---TPU`: wcslib's own alias for TPV applied as a prior distortion (wcs.c wcsset), computed (round-3 probe p36_tpu) | unknown projection code: "TPU" |

The SIP cases above come from the round-3 verifier probes
(`p01_beyond_order_nonzero`, `p02_order10`, `p38_neg_order`,
`p34_order_float`, `p35_coef_string`, `p05_mixed_suffix`) with explicit
identity AP/BP added.

#### Headers wcslib refuses that we accept (`hdrcases_lenient/`)

Same case format; no golden of their own. Each names the wcslib error
(`# EXPECT-WCSLIB-ERROR:`) and an `hdrcases` case (`# SAME-AS:`) whose
golden it must reproduce. `TestHdrLenientThanWcslib` builds the transform
from the case's own cards and requires both directions to agree with the
SAME-AS golden (rejections included, that golden's tolerances), and the
built `hdr2sky` to FAIL on the case as written with the quoted text
(interop as for `hdrcases_reject`). `NOINJECT` turns off the harness's
AP/BP injection for that case.

| Case | wcslib 8.10 error | Why we accept | Same as |
|---|---|---|---|
| `sip_no_apbp` | Degree must not exceed 9, got -1 (no AP/BP card: inverse degree -1) | AP/BP are only the starting guess of the inversion of A/B (in wcslib's disx2p too); most SIP headers (HST, astrometry.net without `-SIP` inverse) carry none | `sip4_no_apbp` |
| `sip_ap_one_axis` | Degree must not exceed 9, got -1 (AP_ORDER/BP_ORDER given, BP without a coefficient; round-3 probe p03_ap_only) | as above | `sip4_no_apbp` |
| `dss_plateid_last` | No distortion type, yet DQia keyvalues are present (the lexer's PLATEID rule resets distran; round-3 probe probedss/dss_plateid_last) | keyword order is free (FITS Standard 4.0 sect. 4.4.1); with PLATEID first, as written by the archive, wcslib applies the plate solution | `real_dss_plate_solution` |

#### Headers wcslib accepts that we read differently (`hdrcases_differ/`)

Same case format; no golden of their own (`make hdr-goldens` reads
`hdrcases/` only). A deliberate difference from wcslib with an
independent truth: each case names, with `# SAME-AS: <name>`, the
`hdrcases` case whose wcslib golden is the answer (the same header and
points without the cards we ignore), and pins how far wcslib's own reading
of the case as written is from it (`# WCSLIB-DIFFERS-ARCSEC: lo hi`) and
how closely astropy 8.0.2 must agree (`# ASTROPY-TOL-ARCSEC`,
`# ASTROPY-TOL-PX`). `TestHdrDeliberateDifferences` requires our transform
of the case to reproduce the SAME-AS golden in both directions (its
tolerances), astropy's reading of the case (`astropy_hdr.py`: its warning
that it removed the SIP set, `all_pix2world` / `all_world2pix` at origin 1)
to agree with the golden and with us, and hdr2sky to compute every point
of the case and land lo..hi arcsec from the golden at every point. The
astropy and hdr2sky readings are recorded tool readings
(`testdata/recorded/wcs/transform/TestHdrDeliberateDifferences.jsonl.gz`,
`TestATLASTPVPlusRedundantSIP.jsonl.gz`; tier 2 runs both live).

| Case | Same as | wcslib 8.10 on the case as written | astropy 8.0.2 / ours |
|---|---|---|---|
| `real_atlas_k2_tpv_sip` | `real_atlas_k2_tpv` | SIP and TPV both: 24.45-26.84 arcsec (13.28-14.59 px) off | SIP set dropped; vs the golden 3.6e-10 arcsec; ours vs astropy 3.7e-10 arcsec, 6.9e-11 px |
| `real_ptf_tpv_plus_sip` | `real_ptf_tpv4` | both: 72.7-133.4 arcsec off | dropped; 3.6e-10 arcsec; 3.4e-10 arcsec, 5.4e-10 px |
| `tpv_sip_nosuffix` | `tpv_sip_nosuffix_tpvonly` | both: 0-8.27 arcsec off (zero at CRPIX, a synthetic non-redundant SIP set) | dropped; 2.0e-10 arcsec; 2.1e-10 arcsec, 3.7e-10 px |

`astropy_hdr.py` reads a case as hdr2sky does (GRID and PIX points,
1-based pixels; the cards padded to 80 columns, END stops) and prints
every float as Python's shortest round-trip repr; it records astropy's
warnings and log messages (the log message saying astropy "is using the
SIP distortion coefficients" is printed while the WCS is built, before
the TPV rule removes them; the warning "Removed redundant SIP distortion
parameters because CTYPE explicitly specifies TPV distortions" and
`sip: false` are what the test checks).

Cases that are known to disagree today are listed, with their audit
reference, in `hdrKnownBad` in `hdr_golden_test.go` (and, for the
projection-proof cases, `prjHdrKnownBad` in `prj_hdr_test.go`, merged into
it); such an entry fails the test if the case starts agreeing. Both
maps are empty: every header case agrees with wcslib in both directions.

## Bare-projection goldens (`wref`)

`wref.c` defines 61 cases covering all 28 projection codes wcslib 8.10
implements (AZP SZP TAN STG SIN ARC ZPN ZEA AIR CYP CEA CAR MER SFL PAR
MOL AIT COP COE COD COO BON PCO TSC CSC QSC HPX XPH), each with defaults
and, where the projection has parameters, non-default PVi_m sets: AZP
mu/gamma (incl. |mu| > 1 and tilted), SZP mu/phi_c/theta_c, slant SIN and
NCP, ZPN degree 1, 2, 3 (with and without a stationary point), 5 and 29
(`zpn_pv29`: PV2_1 = 1, PV2_29 = -1e-12, the highest coefficient read), AIR
theta_b, CYP mu/lambda (incl. Gall), CEA lambda, COP/COE/COD/COO with
eta = 0 and eta != 0 (northern and southern cones), BON theta_1 > 0, < 0
and = 0, HPX (H, K) = (4, 3), (3, 2), (6, 1), (5, 5). Every case runs the
same fixed grids (forward: phi -180..180 step 15 plus +-270, theta
-90..90 step 15; inverse: a [-400, 400]^2 step-50 box, a [-180, 180]^2
step-15 box, and wcslib's own forward outputs), so every projection
boundary is crossed and wcslib's accept/reject decision is recorded per
point. On top of the grids, every case gets generic boundary probes
(`boundary_fwd`/`boundary_inv` in `wref.c`, the same code for every
case): on 80 inverse scan lines (28 + 28 off-grid horizontal/vertical
lines, the 9 + 9 exact lines x, y = 45 k, 24 rays from the origin) and 48
forward ones (24 parallels, 24 meridians), every accept/reject crossing
is bisected and the points 1e-7 deg either side of it are written, so a
boundary placed even 1e-6 deg too far out or in fails; crossings where
the boundary only grazes the line (a knife edge decided by last-bit
noise; found at PCO's x = +-180 edge) are skipped by a transversality
test, except where the accepted piece is itself an exact line (PAR's pole
lines). Every finite forward image of a pole gets points 1e-9 to 1e-3 deg
around it in 8 directions; the forward set gets theta = +-(90 - 10^-k),
k = 1..7, at six longitudes and phi = +-180 +- 1e-7 at four latitudes.
Before these probes, a CEA/ZEA/CSC/TSC/QSC/PAR/XPH boundary moved out by
1-25% (or by 1 deg) passed every test; now each such mutant fails.
Two branches the regular grids cannot reach get extra probe points
(`extra_fwd`/`extra_inv` in `wref.c`): points 0.0001-0.005 deg from all
six cube-face centres for TSC/CSC/QSC (QSC's small-angle formula for
1 - zeta < 1e-8, both directions), and points on MOL's pole line
|y| = sqrt(2) r0 at x = 0, 1e-13, 0.5, 2, 5, -3 and 1e-9 either side of
it (wcslib accepts only |x| <= tol there). AIR gets inverse points at
r = 1e6 .. 1e15 deg on two directions, past the r ~ 5e8 deg where airx2s's
30-step bracketing gives out (wcslib accepts 1e6 and 3e8 .. 1e10, rejects
1e8, 1e12 and 1e15; we must do the same).

Small-angle band probes (round 2; `band_fwd`/`band_inv`, every case):
inverse y = +-{5e-5, 9.9e-5, 1.01e-4, 2e-4, 1e-3, 5e-3, 2e-2, 5e-2, 9e-2}
deg at x = -170, -45, 10, 170, and forward theta = +-{5e-5, 9.9e-5,
1.01e-4, 2e-4, 5e-4, 9e-4, 3e-3, 9e-3, 3e-2} deg at phi = -170, -60, 30,
170. Before them PCO's small-angle switches (|y| < 1e-4 in pcox2s,
|theta| < 1e-4 in pcos2x) could be moved to 1e-3 or 1e-1 deg with every
test passing. `TestPCOSmallAngleBand` holds PCO's latitude-like output on
these rows to 4 ULP of wcslib (measured identical), which pins both
switches to 1%. Round 2b added 1e-4 itself and 1.005e-4 to both lists:
the `<=` and 1.01e-4 mutants of both switches had survived (no probe at or
between the thresholds); they now fail. The same test pins PCO's forward
y = -0.0 at theta = 0 bit for bit, which caught a real port fault: the Go
constant `-0.0` is +0, so the port returned +0 there (fixed with
`math.Copysign(0, -1)`).

Pole band probes (round 2b; `pole_fwd`, every case whose prjset theta0 is
90, i.e. the zenithal codes and XPH): forward theta = +-(90 - d) at
phi = -150.3, -41.7, 22.9, 113.1 for d = 1e-5 10^(k/10) deg, k = 0..40,
and d = T (1 +- 1e-3), T (1 +- 1e-2) for T = 1e-5 rad (SIN's forward
switch; the SIN/SZP inverse r^2 < 1e-10), 2e-4 rad (AIR's xi < 1e-4 rad,
both directions) and 1e-4 sqrt(2/3) rad (XPH's sigma < 1e-4), plus
theta = w[5] exactly for XPH; the round trips of these points straddle
the inverse switches. Before them the verifier moved the SIN/SZP inverse
switch to 1e-11 / 0.5e-10, AIR's tolerance to 1.2e-4 / 2e-4 and XPH's to
1.2e-4 .. 2e-4 and every test passed (results up to 5.4e-9 deg from
wcslib inside the annulus). `TestPoleSmallAngleBand` holds the
latitude-like outputs on these rows to 4 ULP of wcslib (measured: 0 ULP
for SIN, SZP, AIR and XPH, at most 2 ULP for STG/AZP); every one of those
thresholds moved by 2%, 20%, 2x or 10x either way, and XPH's forward `<=`,
now fails. Inverse inputs exactly ON the SIN/SZP r^2 = 1e-10 and AIR
r = w[5] thresholds (`tie_inv`: x stepped one ULP at a time until
wcslib's own expression, evaluated in the same order with its w[], ties)
kill the `<=` mutants of those three; the other three comparisons cannot
be tied (`TestSmallAngleSwitchTies`). The goldens grew from 815 KB to
1035 KB with these probes.

Round 2c (the verifier's leftovers of round 2b) added, all generic or
search-based in `wref.c`, the truth always wcslib's output at the probe:

- `band_fwd`/`band_inv`: 9.95e-5 and 9.99e-5 in both lists. PCO's two
  switches moved 0.5% and 0.9% DOWN had survived (nothing between 9.9e-5
  and 1e-4); every move of 0.1% or more either way now fails.
- `pole_fwd`: d = T (1 +- 5e-4) for all three T. SIN's forward switch
  moved exactly 0.1% up had survived (at T (1 + 1e-3) the computed t is
  not below 1.001e-5); it, and a 0.05% move, now fail.
- `centre_fwd` (TSC, CSC, QSC): forward points whose zeco = 1 - zeta, as
  qscs2x computes it (`qsc_zeco`), is 1e-8 (1 +- {5e-4, 1e-3, 1e-2,
  1e-1}), 0.5e-8 and 2e-8, on all six faces in eight directions (found
  by rescaling the angular offset; a probe on the wrong side aborts the
  run). QSC's zeco < 1e-8 switch moved to 0.5e-8 had survived every
  test; `TestQuadcubeFaceCentreBand` holds these rows to 4 ULP (measured
  0) and kills every move of 0.05% or more.
- `xph_seam_inv` (XPH): inverse inputs next to the north pole on the
  x = 0 seam whose raw native longitude (wcslib's xphx2s with its bounds
  check switched off) is past +-180 by at most 1e-12 (wcslib clamps and
  accepts) or by 1e-12 .. 1e-6 (wcslib rejects): x = +-k 1e-17 u at
  y = u = 1e-6 .. 0.1, the first k of each kind. A 2e7-point random
  search had found 47806 such inputs; none was in the goldens, so XPH's
  prjbchk tolerance widened to 1e-6 (or narrowed to 1e-14) survived.
- the case `azp_mu1_gamma45` (AZP mu = 1, gamma = 45): next to its
  horizon azpx2s's theta comes out a few 1e-14 .. 1e-9 deg below -90, so
  its prjbchk(1e-13) decides points (a 4e6-point search over eight AZP
  parameter sets found such inputs only for mu = 1 with a tilt). It
  exposed a real port fault: one ULP of Go's math.Atan2 in atan2d(1, s)
  moved theta from wcslib's -90 - 5.7e-14 (clamped, accepted) to
  -90 - 1.1e-13 (rejected) at (x, y) = (+-45, -165) and
  (408.3147032522826, -180). The proof audit then used a correctly
  rounded atan2 there; since perf/trig the port calls Go's math.Atan2
  again, and the two (+-45, -165) points are accepted as
  boundary-ambiguous decisions (our own decision flips under a 4-ULP
  libm perturbation, `libm_sens_test.go`); the third agrees.

The bare goldens are 1.15 MB in all (66 files of 2.3-37 KB per case).

`libm/wcstrig_guards.json` (`./wref trig libm`, also run by `make
wref-goldens` and checked by `TestWrefGoldensReproduce`) is wcslib 8.10's
own `asind`/`acosd` at +-1 moved in and out by 0 .. 1e-3, either side of
WCSTRIG_TOL = 1e-10, for `TestWcstrigGuards`.

`golden/prjset_status.json` (also from `wref gen`) records wcslib's
`prjset` status for 21 parameter sets, 18 of which wcslib refuses (AZP
mu = -1 or gamma = 90, SZP mu sin(theta_c) = -1, ZPN all-zero or
degree >= 2 with PV2_1 <= 0, AIR theta_b = -90, CYP lambda = 0 or
mu + lambda = 0, CEA lambda outside (0, 1], COP/COE/COD/COO/BON without
PV2_1, COP theta_a = 0, HPX H or K <= 0);
`TestSelectRejectsWhatPrjsetRejects` requires Select to refuse exactly
those.

Regenerate with `make wref-goldens`. `TestWrefGoldensReproduce` rebuilds
`wref`, regenerates into a temporary directory and requires every golden
to be byte-identical (and the case list to match the committed files);
tier 2 (`go test -tags interop`): a failure without a C compiler and
libwcs; `make verify` does the same from the shell. The
goldens committed on branch `proof2/wcs-projections` were generated with
gcc 16.2.1 (`-O2 -std=c99`), glibc 2.44, zlib 1.3.2, x86_64 Linux, wcslib
8.10-1, and reproduce byte for byte in that environment. The payload
holds exact IEEE-754 bits; metadata doubles are the shortest of
%.15g/%.16g/%.17g that reads back exactly.

The previous `wref.c` (commit 9d21799) had three case-definition bugs,
fixed here: AZP/SZP parameters one index low (`prj.pv[0]` is not PV2_1, so
"azp_tilted" was mu = 30, gamma = 0 and every SZP case had theta_c = 90);
a TSC inverse grid read 14 points from a 13-point array; and it covered 8
of the 28 codes with a forward grid that never reached phi > 120.

The Go side is checked against these goldens by `TestWrefProjections`
(budget 1e-11 deg; measured worst case 3.8e-12 deg over all 61 cases
including the boundary probes; inverse longitude is compared as arc
length dphi cos(theta); the latitude is compared as it is, with no
near-pole allowance (TSC's worst is 5.7e-14 deg) - see `wrefDiff`;
forward coordinates beyond 1000 deg, which occur only next to a forward
divergence, get a relative 1e-14 budget - see `wrefFwdScale`).
The probes next to poles and horizons showed that Go's math.Sin/Cos
(last bit) and math.Asin/Acos (up to 1024 and 2.5e7 ULP near +-1) were
magnified into the sky coordinates (STG 4e-5 deg at theta = -89.99, MOL
1.5e-7 deg next to its pole, COE 3e-9 deg). math.Asin/Acos are replaced
by cancellation-free identities (`libAsin`/`libAcos`, within 1 ULP,
`wcstrig.go`); the last-bit sin/cos/atan2 differences remain (the proof
audit's correctly rounded versions cost 3x to 17x in speed) and rows
where they are magnified are accepted within their own libm sensitivity,
counted (`libm_sens_test.go`). wcslib returns a NaN
or infinite coordinate with a zero status at a few points (CYP inverse
where asin's argument leaves [-1, 1], PCO inverse far off the map, COO at
the apex of a southern cone); the golden records null
coordinates and the test requires us to REJECT such points.

The wcslib 8.10 C source used to read the algorithms
(`wcslib-8.10.tar.bz2`, sha256
`b394d3306da9a52410f056136911496e599b5c48f1cff2f221563697077bbd8d`, from
https://www.atnf.csiro.au/computing/software/wcs/wcslib-releases/) is not
committed; the goldens come from the installed Arch `wcslib 8.10-1`
library.

## Reference tables in compact form (`paper/*.json` + `.bin`, `libm/trig.json` + `.bin`, `libm/atan2.json` + `.bin`)

The three mpmath tables below were decimal JSON rows until round 2b of the
proof audit (545 KB in all). They are now written by `f64pack.py`: a small
metadata `.json` (source, versions, row counts, column layout, the
payload's sha256 and the zlib version) and a `.bin` payload with every
double as its exact IEEE-754 bits (big-endian float64, zlib level 9), 108
KB in all. The conversion was checked value by value: the 15912 (HPX),
11840 (XPH) and 8886 (trig) doubles are bit-identical to the old JSON's.
The Go reader (`f64pack_test.go`) checks the length and the sha256; the
reproduce tests (`f64packSame`) require the metadata and the decoded
payload back byte for byte, and the compressed bytes too when the local
zlib is the recorded one (1.3.2, Arch `zlib 1:1.3.2-3`, also what the
venv's Python links).

## libm reference values (`libm/trig.json`)

`gen_libm.py` (mpmath 1.4.1, 200 bits; stdlib `random` with a fixed seed)
writes, for 2054 sin/cos arguments (degree angles as wcslib forms them,
arguments within 1e-16..1e-1 of 0, +-pi/2, +-pi, uniform ones) and 908
asin/acos arguments (uniform over [-1, 1] and within 1e-16..1e-1 of
+-1), the double nearest the exact value. `TestLibmAccuracy` requires
Go's math.Sin/Cos/Sincos within 2 ULP (2^-100 absolute next to their
zeros; measured 2 ULP) and `libAsin`/`libAcos` within 1 ULP (measured 1);
it logs math.Asin/Acos's error on the same arguments (973 and 2.5e7 ULP,
the reason they are not used). `TestLibmGeneratorReproduces` (interop)
re-runs the script and requires the table back (`f64packSame`). These
tables served as correct-rounding truth while the port carried its own
correctly rounded sin/cos/asin/acos/atan2 (removed in perf/trig).

## libm atan2 reference values (`libm/atan2.json`)

Also written by `gen_libm.py` (same seed, drawn after the sin/cos and
asin/acos arguments so `libm/trig.*` is unchanged): 2506 (y, x) pairs and
the double nearest atan2(y, x) (mpmath 1.4.1, 200 bits). The arguments
are (sin phi, cos phi) of degree angles within 1e-12..1e-1 deg of 0,
+-90 and +-180 (qscs2x's small-angle zeco) and over the circle, (1, s)
and (v, 1) for |s|, |v| in 1e-6..1e6 (atan2d(1, s) in azpx2s, atand),
random pairs of every sign with magnitudes 1e-3..1e3, and a few exact
ones. `TestAtan2Accuracy` requires Go's math.Atan2 (and math.Atan on the
(v, 1) rows) within 1 ULP of them (measured 1 ULP; math.Atan2 differs
from the correctly rounded value at 548 of them). glibc 2.44's atan2,
which wcslib calls, returns exactly those doubles for all 2506 (measured
with a C harness when the table was made).

## XPH from its defining paper (`paper/xph_paper.json` + `.bin`)

`gen_xph_paper.py` evaluates XPH straight from the equations of
Calabretta & Lowe 2013, PASA 30, e050, "Representing the 'butterfly'
projection in FITS - projection code XPH" (arXiv:1308.6647v1, PDF sha256
`fd5ff0cd630a713aaedd3029c08b67b6439c2ce28e0c1a62569e7d80c301a9e5`):
projection equations (2)-(9) and deprojection equations (10)-(17), in
mpmath 1.4.1 at 60 digits, rounded once. 1464 forward points (7.5 deg
grid; both poles approached as 90 - m 10^-k, m = 1, 3, k = 1..7;
theta = +-89.95 .. +-89.996 either side of sigma = 1e-4, where the
double-precision code switches to the small-angle form (19); the zone
boundary asin(2/3) +- 1e-7 deg) and 1496 inverse points (the images of
the forward points, south pole excluded since its image is the four map
corners, plus direct probes at sigma = m 10^-k, k = 1..10, next to both
poles). No wcslib and no Go code is involved. `TestXPHPaperEquations`
holds the port to it within 1e-11 deg, plus, where the paper's (9) is in
force and the port (like wcslib) evaluates it in double next to a pole,
the rounding bound of that formula (`xphSigmaRoundingDeg`; its zone ends
at the paper's own sigma < 1e-4, computed in the test, not taken from
the code under test). `TestXPHPolarPrecision` measures each polar zone
against the exact equations and requires the documented worst: forward
north 2.4e-11 deg (theta = 89.994, equation (9) in double, inside its
bound), inverse north and south ~1e-14 deg (the inverse uses (18) next
to both poles), forward SOUTH 4.08e-7 deg at theta = -89.9999997: wcslib's
xphs2x uses (19) only for theta > 89.9953 and keeps (9) next to the south
pole, where sin theta rounds to -1 and the point lands on the map corner.
The port is identical to wcslib there, so this is wcslib 8.10's precision
and is documented as such. `TestXPHPaperGeneratorReproduces` (interop)
re-runs the script and requires the table back (`f64packSame`).

## HPX from its defining paper (`paper/hpx_paper.json` + `.bin`)

`gen_hpx_paper.py` evaluates HPX for general (H, K) straight from
Calabretta & Roukema 2007, MNRAS 381, 865, "Mapping on the HEALPix grid"
(arXiv:0704.0626), section 6, in mpmath 1.4.1 at 60 digits, rounded once:
(H, K) = (4, 3) on a 7.5 deg grid and (5, 5) on a 15 deg grid, both
poles approached as 90 - m 10^-k (m = 1, 3, k = 1..8), the zone boundary
+-1e-7 deg, the inverse at every forward image except the poles (facet
apexes) and polar facet edges (knife edges: their double-rounded image
falls on either side, and wcslib rejects some), plus direct inverse
probes at sigma = m 10^-k (k = 1..10) next to the apexes. Odd K only (the
even-K southern offset is covered by the wcslib golden `hpx_h3_k2`). No
wcslib and no Go code is involved. `TestHPXPaperEquations` holds the
port to it within 5e-14 deg (`hpxPaperTolDeg`; 1e-11 deg until round 2c,
which let a forward sigma scaled by 1 + 1e-13 or an inverse sigma^2 / K
scaled by 1 + 1e-14 pass), plus, in the polar zones, the
double-precision rounding bound of sigma = sqrt(K (1 - |sin theta|))
(forward) and theta = asin(1 - sigma^2 / K) (inverse)
(`hpxPolarRoundingDeg`); measured worst 1.2e-14 deg outside the polar
zones, and in them at most 0.30 of the bound: forward 2.98e-7 deg (theta
= +-(90 - 3e-7), where sin theta rounds to 1 and the point lands on the
facet apex), inverse 4.67e-7 deg (sigma ~ 1e-8, where 1 - sigma^2 / K
rounds to 1). wcslib 8.10 gives the same values bit for bit (the port is
identical to it, `TestWrefProjections`), so these are wcslib's precision
next to the HPX poles, not a port fault; the test also fails if they grow
past the documented figures. `TestHPXPaperGeneratorReproduces` (interop)
re-runs the script and requires the table back (`f64packSame`).
