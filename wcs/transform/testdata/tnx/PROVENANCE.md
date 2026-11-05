# TNX truth: provenance

Truth source: **WCSTools 3.9.7** (11 August 2022, Jessica Mink, SAO),
Arch package `wcstools 3.9.7-1` (built from the upstream tarball under
`/var/lib/aurbuild/pkgbuild/wcstools`, installed binaries
`/usr/bin/xy2sky`, `/usr/bin/sky2xy`). Its `libwcs/tnxpos.c` is a C
translation of IRAF `mwcs/wftnx.x` + `wfgsurfit.x`, the reference TNX
implementation; it reads the WAT strings with `hget.c:hgetm`, which joins
`WATn_001`, `WATn_002`, ... verbatim with no blank trimming. wcslib is not
usable as TNX truth (multi-card WAT -> `WCSHDRERR_PARSER`, audit wcs.md).
No golden value is produced by the Go code under test.

| Path | Generator | What it holds |
|---|---|---|
| `headers/<name>.hdr` | real (see below) or `make_synthetic.py` | FITS header cards, one per line; trailing blanks after the last non-blank may be stripped (cards are padded back to 80 columns) |
| `refused/<name>.hdr` | `make_synthetic.py` | headers `NewTNX` must refuse (see "Refused headers" below) |
| `refused/<name>.json` | `gen_golden.py` | the WCSTools answer for that header, same format as `golden/` (an s2p pixel is `null` where sky2xy prints `nan`; a header WCSTools does not read as celestial records the xy2sky line as `wcstools_no_answer` with no points) |
| `golden/<name>.json` | `gen_golden.py` | `p2s`: `[x, y, ra, dec]` from `xy2sky -d -n 12` on a 7x7 grid + fractional + 10%-outside-the-edge pixels (1-based); `s2p`: `[ra, dec, x, y]` from `sky2xy -n 10` fed the exact printed p2s sky strings; `header_sha256` of the .hdr |

Regenerate (from this directory; stdlib python3 only):

    python3 make_synthetic.py      # rewrites headers/syn_*.hdr and refused/*.hdr
    python3 make_synthetic.py --check  # exit 1 if a syn_*.hdr or refused/*.hdr would change, is missing, or is extra
    python3 gen_golden.py          # rewrites golden/*.json and refused/*.json
    python3 gen_golden.py --check  # exit 1 if any golden would change
    python3 extract_kpno.py <downloaded .fits.fz>   # rewrites headers/kpno_kp4m_mosaic1_ccd1.hdr (needs funpack)

`TestTNXSyntheticHeadersReproduce` runs `make_synthetic.py --check`, then
deletes `headers/syn_*.hdr` and `refused/*.hdr` in a scratch copy and requires
`make_synthetic.py` to write them back byte for byte.
`TestTNXGoldensReproduce` runs `gen_golden.py --check` (tier 2, `go test
-tags interop`: a failure without python3/xy2sky/sky2xy).
`TestTNXWCSToolsGolden` needs only the committed files. It fails on a
missing/orphan golden, a header whose sha256 differs from the golden's,
and on any case of `tnxRequiredCases` (tnx_golden_test.go: every header
below by name, 56 points each) that is missing, thinned, or not listed
there, so deleting or thinning a case fails the default `go test`
without WCSTools.

## Headers

| Header | Source | License | Covers |
|---|---|---|---|
| `noao_ct4m_tnx1` | all 95 cards of the primary header of `https://fits.gsfc.nasa.gov/registry/tnx/tnx.fits` (FITS registry TNX sample; file sha256 `3e574e1cd2ddf6bd7b787eb8684aca2ff0823c8b31bfac474ef1cdf1002e10ab`, 2026-11-01), card for card; data not kept | NASA/GSFC FITS Support Office, public | real CTIO 4 m Mosaic II (NOAO) chip `im13`, MSCRED solution by F. Valdes: Chebyshev (type 1), 4x4, half cross terms, CD matrix, 5-card WAT per axis with chunk boundaries inside numbers (`0.305` / `7102603054315`, `2.327260049467782E-` / `5`) |
| `registry_tnxdoc_tnx3` | cards of figure 1 of the registry TNX document `https://fits.gsfc.nasa.gov/registry/tnx/tnxdoc.html` (copy of the IRAF document; page sha256 `cd555a2211c5bd619eb783343d6731a6e9d1d14d1d7d6c5d30d20b33e7613d2d`), verbatim, plus SIMPLE/BITPIX/NAXIS/NAXIS1=2048/NAXIS2=4096 (not in the figure; needed for a standalone file, they do not enter the WCS) | as above | real ct4m solution: simple polynomial (type 3, raw xi/eta, no normalisation), 4x4, half cross terms; chunk `WAT1_003` ends in a blank that separates two numbers (the header layer trims it: the concatWAT bug) |
| `kpno_kp4m_mosaic1_ccd1` | the `ccd1` extension header of the NOIRLab Astro Data Archive file md5 `42f68c5bef0524074d565f93b687c873` (KPNO 4 m Mosaic-1 InstCal, NGC4594, DATE-OBS 2006-02-26, `K4M06A_20060223-R-kp4m20060226T094108`; `.fits.fz` sha256 `fea97d3341a1bf3d8fa14594e7d83adaada9ab46132dde64ea031bce66a8d664`, 41 564 160 bytes, retrieved 2026-11-01 without login from `https://astroarchive.noirlab.edu/api/retrieve/42f68c5bef0524074d565f93b687c873/`), written by `extract_kpno.py` (cfitsio 4.7.0 `funpack -C -S`, the extension card for card except XTENSION -> SIMPLE and PCOUNT/GCOUNT/CHECKSUM/DATASUM dropped; `extract_kpno.py --check <file>` reproduces it byte for byte); data not kept | NOIRLab Astro Data Archive, public (served without authentication; NOIRLab data become public after the proprietary period) | real KPNO 4 m Mosaic-1 MSCRED solution: simple polynomial (type 3), 4x4, FULL cross terms (16 coefficients), CD matrix with CRPIX far off the chip (4198, 4113); 7-card WAT per axis, `WAT1_003` ends in a blank that separates two numbers, `WAT2_004` starts with one, and `WAT1_006/007`, `WAT2_006/007` come after `WAT2_005` and ~30 unrelated cards (interleaved) |
| `syn_leg_full` | `make_synthetic.py` | n/a | Legendre (type 2), 3x3 full cross terms, rotated CD; chunk ends in a blank and splits `-9.9e-0`/`5` |
| `syn_cheb_none_pc` | `make_synthetic.py` | n/a | Chebyshev, xorder 4 / yorder 3, no cross terms, PC + CDELT (no CD) |
| `syn_poly_half_pole` | `make_synthetic.py` | n/a | polynomial, xorder 3 / yorder 4 half cross terms, CRVAL2 = +90 with no LONPOLE (TNX default LONPOLE = 180, not Paper II's 0) |
| `syn_poly_full_watlongpole` | `make_synthetic.py` | n/a | polynomial 2x2 full cross terms, `longpole=170.` given only as a WAT attribute |
| `syn_leg_lonpole_card` | `make_synthetic.py` | n/a | `syn_leg_full` surfaces plus a `LONPOLE = 175` card and no WAT longpole: the card is used, not the TNX default 180 |
| `syn_poly_lonpole_card_and_wat` | `make_synthetic.py` | n/a | `syn_poly_full_watlongpole` plus a `LONPOLE = 175` card: the card beats the WAT `longpole=170.` (tnxinit reads the WAT value only when the card is absent) |
| `syn_poly_wat_precedence` | `make_synthetic.py` | n/a | `longpole` in both WAT1 (165) and WAT2 (195), and `lngcor` and `latcor` in both strings with different coefficients: WAT1 longpole wins, lngcor comes from WAT1, latcor from WAT2 (tnxinit order) |
| `syn_poly_cor_fallback` | `make_synthetic.py` | n/a | `lngcor` only in WAT2 and `latcor` only in WAT1: each is read from the other axis's string (tnxinit `igets(str1)` then `igets(str2)`) |
| `syn_poly_galactic`, `syn_poly_ecliptic`, `syn_poly_supergalactic`, `syn_poly_helioecliptic` | `make_synthetic.py` | n/a | the other celestial pairs: GLON/GLAT (xy2sky "galactic"), ELON/ELAT ("ecliptic"), SLON/SLAT and HLON/HLAT (no WCSTools system name; same TNX math in the native frame). xy2sky prints a system it cannot name in a fixed-width field that overflows at 12 decimals from 10 deg, so the SLON/HLON reference points are below 10 deg |
| `syn_tan_no_wat` | `make_synthetic.py` | n/a | `RA---TNX`/`DEC--TNX` with no WAT1/WAT2 cards: plain TAN with the TNX default LONPOLE 180 (WCSTools falls back to TAN) |
| `syn_poly_pv2_ignored` | `make_synthetic.py` | n/a | `syn_poly_full_watlongpole` plus `PV2_1`/`PV2_2`: latitude-axis PVi_m are ignored by WCSTools and by wcslib TAN, and accepted |
| `syn_poly_cunit_deg` | `make_synthetic.py` | n/a | `CUNIT1/2 = 'deg'`, the one CUNIT accepted besides blank |
| `syn_poly_lonpole_360` | `make_synthetic.py` | n/a | `LONPOLE = 360`, the largest value WCSTools reads (`tnxpos.c` `if (wcs->longpole > 360.0)` falls back to the WAT longpole; `wcsinit.c` sets 999 when the card is absent) |
| `syn_poly_lonpole_neg400` | `make_synthetic.py` | n/a | `LONPOLE = -400`: negative values are used as given by WCSTools |
| `syn_poly_wat_nondigit_keys` | `make_synthetic.py` | n/a | extra `WAT1_A` and `WAT2_1X` cards: not WAT cards (hgetm looks up numbered keys only), ignored |
| `syn_poly_projp1`, `syn_tan_no_wat_projp1` | `make_synthetic.py` | n/a | `PROJP1 = 0.5` with and without WAT cards: `wcs.Parse` reads PROJPn as the latitude-axis PVi_n (wcslib), which TAN ignores; WCSTools does not read PROJPn; accepted |
| `syn_poly_pc001001` | `make_synthetic.py` | n/a | CDELT + `PC001001..PC002002` (no PCi_j): WCSTools `wcsinit` reads PC001001 when PC1_1 is absent, wcslib reads PC00i00j as PCi_j |
| `syn_poly_pc_no_pc2_2` | `make_synthetic.py` | n/a | CDELT + PC1_1, PC1_2, PC2_1 (PC2_2 absent, 1 in both) |
| `syn_poly_cd_offdiag_only` | `make_synthetic.py` | n/a | only CD1_2 and CD2_1 (CD1_1, CD2_2 default 0 in both) |
| `syn_poly_cd_and_cdelt`, `syn_poly_cd_and_cdelt1` | `make_synthetic.py` | n/a | CD beside CDELT1/CDELT2, and beside CDELT1 alone: both readers use CD and ignore CDELT (the CDELT-pair refusal applies only without CD) |
| `syn_poly_cd_and_cd001001` | `make_synthetic.py` | n/a | CD00i00j beside CDi_j with equal values (WCSTools reads CDi_j, wcs.Parse both) |
| `syn_poly_no_scale` | `make_synthetic.py` | n/a | no CD, CDELT or PC on a 20 x 20 image: 1 deg/pixel in both (`wcsinit` no-scaling branch, wcslib CDELT default 1) |
| `syn_poly_lngcor_only_pv2` | `make_synthetic.py` | n/a | lngcor only (WAT2 has no latcor) plus PV2_1/PV2_2/PV2_4: tnxinit succeeds with one surface, so WCSTools ignores the PVi_m, as wcslib TAN does |
| `syn_leg_wat_1999` | `make_synthetic.py` | n/a | Legendre 9x9 full cross terms whose joined WAT1 string is exactly 1999 characters, the most WCSTools keeps (2000-byte buffer, `SZ_ATSTRING`), with `longpole=175.5` as its last token: a reader that cut one character would see `175.` (checked: at 2000 characters xy2sky gives exactly the `longpole=175.` answer) |

Rejections (no WCSTools equivalent: WCSTools never checks `wtype` and reads
past short coefficient lists) are pinned in `tnx_test.go`
(`TestTNXRejects`): wtype=zpx in WAT1 or WAT2, wrong coefficient count,
unknown type / cross code, fractional or zero order, empty domain,
unparsable number, missing quote or '=', non-default `ro` (including
57.3, within 1e-5 of 180/pi, and in WAT2), CROTA1 or CROTA2 without CD, a
PVi_0 fiducial offset or PV1_4 (refused as longitude-axis PVi_m), a mixed axis pair (RA/GLAT), and celestial axes in
DEC/RA order (WCSTools accepts that; `NewTNX` refuses it loudly,
inherited from `New`).

Refused headers (`refused/`, `TestTNXRefusedHeaders`): headers that
WCSTools reads one way and the wcslib conventions of `wcs.Parse` / `New`
another, or where WCSTools ignores a card. Each is a variant of
`syn_poly_full_watlongpole`; its `refused/<name>.json` is the WCSTools
answer, and before the refusal existed `NewTNX` returned a different
answer with no error (worst p2s disagreement over the 56 points):

| Case | Cards | WCSTools | NewTNX before | Disagreement |
|---|---|---|---|---|
| `pv1_fiducial` | `PV1_1 = 0`, `PV1_2 = 89` | ignores PVi_m | wcslib fiducial (phi0, theta0) | 3601 arcsec |
| `pv1_3_lonpole` | `PV1_3 = 150`, no WAT longpole | LONPOLE 180 | 150 | 253 arcsec |
| `cd_and_pc` | CD plus `CDELT` and `PC` | CD | PC (wcslib) | 69 arcsec |
| `lonpole_361`, `lonpole_400` | `LONPOLE` above 360 | WAT longpole 170 | the card | 975 / 887 arcsec |
| `cunit_arcmin` | `CUNIT1/2 = 'arcmin'`, `CDELT = 0.012` | CDELT read as degrees | arcmin scaled | 4.6e5 arcsec |
| `wat_duplicate_key` | a second `WAT1_002` | first card (`ksearch`) | first card | none (the second card's text is lost) |
| `cdelt1_only` | `CDELT1 = -0.0002`, no CDELT2, CD or PC | `wcsinit` sets CDELT2 = -CDELT1 | wcslib CDELT2 = 1 | 2.99e5 arcsec |
| `cdelt2_only` | `CDELT2 = 0.0002`, no CDELT1 | no-scaling branch: 1 deg/pixel, CDELT2 ignored (CDELT scaling read only when CDELT1 exists) | CDELT2 used | 5.22e5 arcsec |
| `pc_no_cdelt` | PC1_1..PC2_2, no CDELT or CD | no-scaling branch, PC ignored | PC with CDELT 1 | 3.22e5 arcsec |
| `pc_no_pc1_1` | CDELT + PC1_2, PC2_1 only | PC read only when PC1_1 (or PC001001) is present: CROTA path, identity | off-diagonal PC applied | 48.9 arcsec |
| `cd001001_only` | CDELT + CD001001..CD002002, no CDi_j | reads CDi_j only: CDELT | CD00i00j read as CD (wcslib) | 69.3 arcsec |
| `nowat_pv2`, `watnocor_pv2` | no WAT1/WAT2 cards, or WAT strings with no lngcor/latcor, plus `PV2_1 = 1`, `PV2_2 = 0.001`, `PV2_4 = 0.01` | tnxinit fails, the projection is demoted to TAN and the TPV block applies the PVi_m (PV1_m default 0: RA = 210.5 at every pixel; sky2xy prints `nan`) | plain TAN | 345 arcsec |
| `tnx_sip_suffix` | `RA---TNX-SIP`/`DEC--TNX-SIP` plus A_ORDER = 2, A_2_0, A_0_2, B_ORDER = 2, B_1_1 (no AP/BP) | `setdistcode` applies SIP to the pixel offsets before TNX (sky2xy without AP/BP inverts without SIP: 0.69 px from its own forward) | SIP dropped | 0.558 arcsec |
| `tnx_sip_nosuffix` | the same SIP cards with `RA---TNX` | SIP cards ignored | SIP not read (agreed, 2.0e-9 arcsec) | none; refused because WCSTools ignores the cards and wcslib 8.10 applies SIP with TAN and drops the WAT correction |
| `ctype_lowercase` | `CTYPE1/2 = 'RA---tnx'/'DEC--tnx'` | not TNX: LINEAR axes, xy2sky prints asterisks (`wcstools_no_answer`) | TNX answer | no truth to compare |

`NewTNX` refuses every PVi_m on the longitude axis (PV1_0 .. PV1_4 all
change the wcslib rotation; WCSTools reads none), latitude-axis PVi_m
when there is no correction surface, CD with PC, CD00i00j without its
CDi_j, (with no CD) one of CDELT1/CDELT2 alone or a PC matrix without
CDELT1 or without PC1_1/PC001001, SIP keywords with or without `-SIP`,
a CTYPE that is not upper case, LONPOLE above 360, a celestial CUNIT
other than blank or `deg`, and a repeated WATn_k key. Latitude-axis
PVi_m with a surface, PROJPn, LONPOLE = 360, negative LONPOLE and the
linear-step forms both readers agree on are accepted and proven by
goldens.

The `refused/` table rows above `cdelt1_only` came from the round-2
verifier probes; `cdelt1_only` .. `ctype_lowercase` from the round-3
verifier probes (`verify3-wcs-frames-tnx/probes`, `probes2`), each
reproduced here by `make_synthetic.py`.

What WCSTools silently drops is refused too (`TestConcatWATRefusesWhatWCSToolsDrops`):
WCSTools `hgetm` stops at the first missing `WATn_k` (a gap) and ignores
the cards after it, reads no key number above 499, reads only the first
card of a CONTINUE'd string, and cuts a joined string at 1999 characters.
`NewTNX` returns an error for each (gap, mixed `_001`/`_2` spellings, the
first card in two spellings, no first card, key 500, CONTINUE, 2000
characters) rather than a transform
that loses cards.

Real samples cover the CTIO 4 m Mosaic II (Chebyshev half cross,
polynomial half cross) and the KPNO 4 m Mosaic-1 (polynomial full
cross). No public real Legendre (type 2) TNX header was found; Legendre
and the LONPOLE / WAT precedence rules are proven on synthetic headers
only. Search record (2026-11-01): GitHub code search for TNX WAT strings
(`lngcor`/`latcor` with WAT keys, ~900 files fetched: wcslib
`TNX.keyrec`, Starlink AST `tnx.head` / `tnx-cheb.head`, LCOGT
neoexchange `example_kpno_mef_extnhdr`, IRAF forum posts) found only types
1 and 3 (the one type-2 string, in Starlink AST `testfitschanerrors.c`, is
a hand-made ZPX error case). A scan of 574 headers from the NOIRLab Astro
Data Archive header API (`/api/header/<md5>/`, raw / instcal / resampled
/ stacked products of 40 instruments including mosaic_1, mosaic_1_1,
mosaic_2, mosaic3, newfirm, 90prime, y4kcam, odi, whirc, kosmos) found
Chebyshev in raw Mosaic headers (MSCRED database solution), polynomial in
every pipeline InstCal Mosaic/NEWFIRM header, TAN/TPV elsewhere, and no
Legendre. Round 2 (2026-11-01) extended the search, again with no type 2
found: a by-epoch scan of 172 more raw headers (two per instrument per
half year, 1997-2024: mosaic, mosaic_1, mosaic_1_1, mosaic_2, newfirm,
kosmos, 90prime; Chebyshev or polynomial only); full clones of
iraf-community/iraf, iraf-mscred and pyraf, astropy, astropy-data,
photutils-datasets, ccdproc, stsci.tools and starlink-pyast, every file
and every gzip/bzip2 member grepped for `wtype=tnx` (only AST's
`tnx.head`, polynomial); web search for Legendre `lngcor`/`latcor`
strings (only the registry document and the Lowell copy of it, and a
CTIO ISPI page with a polynomial solution). Round 3 (2026-11-01), again
no type 2: the archive's instrument list (`/api/adv_search/cat_lists/`)
names arcon, ice, minimo/ice, mop/ice, pfccd and wfc, which the
earlier scans had not covered, but `adv_search/find` returns no files
for any of them; 15 high-level-science-product headers (hlsp_mosaic,
hlsp_mosaic_2, hlsp_ir_imager: NDWFS, the Local Group Survey and
ChaMPlane stacks, the user-reduced products most likely to carry a
hand-chosen ccmap function) are all resampled RA---TAN; GitHub code
search for `lngcor` (100 files) finds only TNX readers and writers
(IDLAstro, PHOTRED, SExtractor/SCAMP/SWarp, WCSTools copies), no header
with a Legendre string; web search repeats the registry and Lowell
copies.

The basis derivatives behind `TNX.Inverse`'s Newton Jacobian are checked
against the closed forms of T_n, P_n and x^n (n <= 5) and, for whole
surfaces of the real headers, against a central difference
(`tnx_deriv_test.go`).

Measured agreement (Go vs these goldens, `go test -v -run TestTNXWCSToolsGolden`):
p2s worst 2.5e-9 arcsec (the 12-decimal print rounding), budget 1e-8;
SkyToPixel of the golden sky vs its source pixel worst 7.4e-9 px, budget
1e-7; vs sky2xy worst 1.4e-5 px (syn_leg_wat_1999; WCSTools' own Newton tolerance, 2.8e-8
deg), budget 1e-4.
