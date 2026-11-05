#!/usr/bin/env python3
"""Generate the OpenCV (cv2) goldens used by the CLAHE tests.

Truth source: OpenCV itself (opencv-python-headless 5.0.0 in the project
venv). Siril's `clahe` command is a thin wrapper around OpenCV
(opencv/opencv.cpp:806 cvClahe_float), so the building blocks are pinned
against OpenCV directly and the whole Siril pipeline against the verbatim
Siril oracle (sref/sref_clahe, see gen_siril_goldens.py).

  golden/clahe_mono_u8.json    cv2 CLAHE on an 8-bit image (256 bins)
  golden/clahe_mono_u16.json   cv2 CLAHE on 16-bit images (65536 bins),
                               several sizes / grids incl. padded tiles
  golden/clahe_clip_cv2.json   cv2 CLAHE (8- and 16-bit) on images whose tiles
                               are large enough that the integer clip limit
                               int(clip*tileTotal/histSize) is >= 2, so the
                               clip limit decides the output (the cases above
                               all clamp it to 1)
  golden/lab_cv2.json          cv2 float BGR2Lab and Lab2BGR on test vectors
                               (all 33^3 LUT grid nodes, random, out of
                               range; rows of widths 1, 13 and 64 so both the
                               SIMD body and the scalar tail of OpenCV's
                               loops are covered)
  golden/nan_cv2.json          cv2 BGR2Lab / Lab2BGR with NaN components and
                               saturate_cast<int>(float) on NaN, +-inf and
                               out-of-range values
  golden/stretch_clahe_cv2.json  Siril's float CLAHE pipeline re-built from
                               cv2 calls on golden/inputs/synth.f32; must
                               equal the verbatim-Siril golden
                               golden/sref/clahe_rgb_3_8 (checked here:
                               independent OpenCV build)
  golden/jpeg_dqt_cv2.json     the quantization tables libjpeg-turbo (through
                               cv2.imencode) writes at each JPEG quality
  golden/inputs/cv2_*          the integer and vector inputs (little-endian)

Arrays use the compact exact format of goldenfmt.py (sha256 of every value
plus a deterministic sample).

Usage: ../../venv/bin/python gen_clahe_cv2.py [--out DIR]   (re-running
reproduces the committed files byte for byte, TestCLAHECv2GoldensReproduce).
"""

import json
import os
import sys

import cv2
import numpy as np

from goldenfmt import Inputs, arr16, arr32, dump, inline_input

HERE = os.path.dirname(os.path.abspath(__file__))
GOLD = os.path.join(HERE, "golden")
U16MAX = 65535.0


def convert_to(src, dst_dtype, alpha):
    """Mat::convertTo(dst, type, alpha) for 32F<->16U: OpenCV computes
    saturate_cast<dst>(src * (float)alpha + 0.f) in float32 (core
    convert_scale.simd.hpp cvt_32f) with round-half-even. Cross-checked
    against the C++ oracle by the final assertion of this script."""
    v = src.astype(np.float32) * np.float32(alpha)
    if dst_dtype == np.uint16:
        return np.clip(np.rint(v), 0, 65535).astype(np.uint16)
    return v.astype(np.float32)


def siril_clahe_float(planes, clip, tile):
    """Siril 1.4.0 cvClahe_float re-built from the same cv2 calls."""
    cl = cv2.createCLAHE()
    cl.setClipLimit(clip)
    cl.setTilesGridSize((tile, tile))
    if planes.shape[0] == 3:
        r, g, b = planes
        bgr = np.dstack([b, g, r]).astype(np.float32)          # fits_to_bgrbgr_float
        lab = cv2.cvtColor(bgr, cv2.COLOR_BGR2Lab)
        L = convert_to(lab[:, :, 0], np.uint16, U16MAX / 100.0)
        L = cl.apply(L)
        lab = lab.copy()
        lab[:, :, 0] = convert_to(L, np.float32, 100.0 / U16MAX)
        out = cv2.cvtColor(lab, cv2.COLOR_Lab2BGR)
        return np.stack([out[:, :, 2], out[:, :, 1], out[:, :, 0]])
    m = convert_to(planes[0], np.uint16, U16MAX)
    m = cl.apply(m)
    return convert_to(m, np.float32, 1.0 / U16MAX)[None]


def dqt_tables(jpeg):
    """The quantization tables (DQT segments) of a JPEG stream, in file
    order: table id, precision and the 64 values in zigzag order."""
    out, i = [], 2
    while i + 4 <= len(jpeg):
        if jpeg[i] != 0xFF:
            raise ValueError("bad JPEG marker")
        m, n = jpeg[i + 1], (jpeg[i + 2] << 8) | jpeg[i + 3]
        if m == 0xDB:
            seg, k = jpeg[i + 4:i + 2 + n], 0
            while k < len(seg):
                pq, tq = seg[k] >> 4, seg[k] & 15
                size = 128 if pq else 64
                v = seg[k + 1:k + 1 + size]
                if pq:
                    v = [(v[2 * j] << 8) | v[2 * j + 1] for j in range(64)]
                out.append({"id": tq, "bits": 16 if pq else 8, "q": [int(x) for x in v]})
                k += 1 + size
        if m == 0xDA:
            break
        i += 2 + n
    return out


def main():
    out = GOLD
    if "--out" in sys.argv:
        out = os.path.abspath(sys.argv[sys.argv.index("--out") + 1])
    os.makedirs(out, exist_ok=True)
    ins = Inputs(out)

    # --- 8-bit CLAHE ------------------------------------------------------
    rng = np.random.RandomState(3)
    nx, ny = 64, 48
    img8 = np.clip(np.linspace(10, 78, nx)[None, :] + rng.normal(0, 3, (ny, nx)), 0, 255).astype(np.uint8)
    res8 = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(8, 8)).apply(img8)
    ins.add("cv2_clahe_u8.u16", img8.astype(np.uint16))
    rin, _ = ins.ref("cv2_clahe_u8.u16")
    dump(os.path.join(out, "clahe_mono_u8.json"),
         {"source": f"cv2 {cv2.__version__} createCLAHE(3.0,(8,8)).apply(uint8)",
          "nx": nx, "ny": ny, "clipLimit": 3.0, "tileGridSize": 8, "input": rin, "out": arr16(res8)})

    # --- 16-bit CLAHE -----------------------------------------------------
    cases = []
    rng = np.random.RandomState(5)
    for k, (nx, ny, clip, tx, ty) in enumerate([(64, 48, 3.0, 8, 8), (37, 29, 2.5, 8, 8), (50, 41, 4.0, 5, 3),
                                                (16, 16, 1.0, 2, 2), (23, 17, 0.0, 4, 4), (31, 9, 40.0, 3, 2)]):
        base = np.linspace(2000, 30000, nx)[None, :] + 4000 * np.sin(np.arange(ny)[:, None] / 5.0)
        img = np.clip(base + rng.normal(0, 800, (ny, nx)), 0, 65535).astype(np.uint16)
        img[0, 0] = 0
        img[-1, -1] = 65535
        res = cv2.createCLAHE(clipLimit=clip, tileGridSize=(tx, ty)).apply(img)
        ins.add(f"cv2_clahe_u16_{k}.u16", img)
        rin, _ = ins.ref(f"cv2_clahe_u16_{k}.u16")
        cases.append({"nx": nx, "ny": ny, "clipLimit": clip, "tilesX": tx, "tilesY": ty,
                      "input": rin, "out": arr16(res)})
    dump(os.path.join(out, "clahe_mono_u16.json"),
         {"source": f"cv2 {cv2.__version__} createCLAHE(clip,(tx,ty)).apply(uint16)", "cases": cases})

    # --- clip limit in effect (8- and 16-bit) -----------------------------
    def tile_total(nx, ny, tx, ty):
        # clahe.cpp CLAHE_Impl::apply: pad to a multiple of the grid
        px = nx if nx % tx == 0 else nx + tx - nx % tx
        py = ny if ny % ty == 0 else ny + ty - ny % ty
        return (px // tx) * (py // ty)

    clip_cases = []
    rng = np.random.RandomState(13)
    for k, (bits, nx, ny, clip, tx, ty) in enumerate([(8, 128, 128, 3.0, 8, 8), (8, 160, 120, 2.5, 4, 4),
                                                      (8, 100, 90, 1.3, 3, 3), (16, 160, 128, 40.0, 2, 2),
                                                      (16, 128, 128, 20.0, 1, 1), (16, 150, 110, 300.0, 3, 3)]):
        hs = 256 if bits == 8 else 65536
        lo, span, sigma = (40.0, 30.0, 3.0) if bits == 8 else (12000.0, 900.0, 40.0)
        yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
        img = lo + span * (xx / nx + 0.5 * yy / ny) + rng.normal(0, sigma, (ny, nx))
        for sx, sy in rng.uniform(0, 1, (5, 2)) * [nx, ny]:
            img += (hs - 1 - lo) * 0.8 * np.exp(-((xx - sx) ** 2 + (yy - sy) ** 2) / 6.0)
        img = np.clip(np.rint(img), 0, hs - 1).astype(np.uint8 if bits == 8 else np.uint16)
        tt = tile_total(nx, ny, tx, ty)
        limit = int(clip * tt / hs)
        if limit < 2:
            sys.exit(f"clip case {bits}-bit {nx}x{ny}: effective clip limit {limit} < 2")
        res = cv2.createCLAHE(clipLimit=clip, tileGridSize=(tx, ty)).apply(img)
        for kk in (limit - 1, limit + 1):
            other = cv2.createCLAHE(clipLimit=(kk + 0.5) * hs / tt, tileGridSize=(tx, ty)).apply(img)
            if np.array_equal(other, res):
                sys.exit(f"clip case {bits}-bit {nx}x{ny}: output does not depend on the clip limit")
        ins.add(f"cv2_clahe_clip_{k}.u16", img.astype(np.uint16))
        rin, _ = ins.ref(f"cv2_clahe_clip_{k}.u16")
        clip_cases.append({"bits": bits, "nx": nx, "ny": ny, "clipLimit": clip, "tilesX": tx, "tilesY": ty,
                           "clip_limit_int": limit, "input": rin, "out": arr16(res)})
    dump(os.path.join(out, "clahe_clip_cv2.json"),
         {"source": f"cv2 {cv2.__version__} createCLAHE(clip,(tx,ty)).apply(uint8|uint16)",
          "check": "for each case cv2 at clip limits limit-1 and limit+1 gives a different output",
          "cases": clip_cases})

    # --- float Lab conversions -------------------------------------------
    rng = np.random.RandomState(9)
    g = np.arange(33, dtype=np.float32) / np.float32(32)
    grid = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(1, -1, 3)   # every LUT node
    rnd = rng.uniform(0, 1, (1500, 3)).astype(np.float32)
    edge = np.array([[0, 0, 0], [1, 1, 1], [-0.1, 0.5, 1.2], [1e-7, 2e-6, 0.004], [0.04045, 0.5, 0.999],
                     [np.float32(1) / 3, 0.25, 0.75], [2.0, -3.0, 0.5]], np.float32)
    bgr_vecs = np.concatenate([rnd, edge]).astype(np.float32)
    ins.add("cv2_lab_bgr.f32", bgr_vecs)
    bgr_ref, _ = ins.ref("cv2_lab_bgr.f32")
    # grid entry: input implied, node (p, q, r) = (B, G, R) = (p, q, r) / 32
    # at index (p*33 + q)*33 + r
    fwd = [{"width": grid.shape[1], "grid": 33, "out": arr32(cv2.cvtColor(grid, cv2.COLOR_BGR2Lab))}]
    for width in (1, 13, 64):
        n = (len(bgr_vecs) // width) * width
        img = bgr_vecs[:n].reshape(-1, width, 3)
        fwd.append({"width": width, "vectors": n, "out": arr32(cv2.cvtColor(img, cv2.COLOR_BGR2Lab))})
    lab_vecs = np.stack([rng.uniform(0, 100, 3000), rng.uniform(-110, 110, 3000), rng.uniform(-110, 110, 3000)], -1)
    lab_vecs = np.concatenate([lab_vecs, [[0, 0, 0], [100, 0, 0], [8, 0, 0], [7.99, 3, -3], [8.01, -3, 3],
                                          [50, 127, -128], [0.5, 100, 100], [99.9, -100, 90]]]).astype(np.float32)
    ins.add("cv2_lab_lab.f32", lab_vecs)
    lab_ref, _ = ins.ref("cv2_lab_lab.f32")
    inv = []
    for width in (1, 13, 64):
        n = (len(lab_vecs) // width) * width
        img = lab_vecs[:n].reshape(-1, width, 3)
        inv.append({"width": width, "vectors": n, "out": arr32(cv2.cvtColor(img, cv2.COLOR_Lab2BGR))})
    dump(os.path.join(out, "lab_cv2.json"),
         {"source": f"cv2 {cv2.__version__} cvtColor float32 COLOR_BGR2Lab / COLOR_Lab2BGR",
          "layout": "rows of `width` BGR (or Lab) triples taken from the first `vectors` input vectors",
          "bgr_input": bgr_ref, "lab_input": lab_ref, "forward": fwd, "inverse": inv})

    # --- NaN and out-of-range rounding ------------------------------------
    # Lab<->BGR with NaN components (OpenCV's SSE v_min/v_max return the
    # second operand for NaN; the scalar code clamps with std::min/max) and
    # saturate_cast<int>(float) = cvRound, which gives INT_MIN for NaN and
    # out-of-range values (cvtss2si / cvtps2dq "integer indefinite").
    # cv2.add(x, 0, dtype=CV_32S) converts float32 to int32 with that cast.
    nan = np.float32(np.nan)
    rng2 = np.random.RandomState(13)
    lab_nan = np.stack([rng2.uniform(0, 100, 96), rng2.uniform(-110, 110, 96), rng2.uniform(-110, 110, 96)], -1)
    bgr_nan = rng2.uniform(0, 1, (96, 3))
    for k in range(96):
        m = k % 8
        for c in range(3):
            if m & (1 << c):
                lab_nan[k, c] = nan
                bgr_nan[k, c] = nan
    lab_nan, bgr_nan = lab_nan.astype(np.float32), bgr_nan.astype(np.float32)
    nan_fwd, nan_inv = [], []
    for width in (1, 13, 64):
        n = (96 // width) * width
        b = bgr_nan[:n].reshape(-1, width, 3)
        l = lab_nan[:n].reshape(-1, width, 3)
        nan_fwd.append({"width": width, "input": inline_input(b), "out": arr32(cv2.cvtColor(b, cv2.COLOR_BGR2Lab))})
        nan_inv.append({"width": width, "input": inline_input(l), "out": arr32(cv2.cvtColor(l, cv2.COLOR_Lab2BGR))})
    rnd_in = np.array([np.nan, 3e9, -3e9, np.inf, -np.inf, 2.5, 3.5, -2.5, -3.5, 0.5, 1.5,
                       2147483520.0, -2147483648.0, 0.49999997, 65535.5, -0.5], np.float32)
    rnd_in = np.tile(rnd_in, 5)          # 80 values: SIMD body and scalar tail
    rnd_out = cv2.add(rnd_in.reshape(1, -1), 0.0, dtype=cv2.CV_32S).reshape(-1)
    dump(os.path.join(out, "nan_cv2.json"),
         {"source": f"cv2 {cv2.__version__}: cvtColor float32 BGR2Lab / Lab2BGR with NaN "
                    "components; cv2.add(x, 0, dtype=CV_32S) for saturate_cast<int>(float)",
          "forward": nan_fwd, "inverse": nan_inv,
          "round_in": inline_input(rnd_in), "round_out": [int(v) for v in rnd_out]})

    # --- Siril's float pipeline rebuilt from cv2 calls --------------------
    synth = np.fromfile(os.path.join(GOLD, "inputs", "synth.f32"), "<f4").reshape(3, 48, 64)
    sref = json.load(open(os.path.join(GOLD, "sref", "clahe_rgb_3_8.json")))
    res = siril_clahe_float(synth, 3.0, 8)
    rec = arr32(res)
    dump(os.path.join(out, "stretch_clahe_cv2.json"),
         {"name": "clahe", "command": "clahe 3 8",
          "note": "Siril 1.4.0 cvClahe_float rebuilt from cv2 calls (BGR2Lab float, L*655.35 -> "
                  "16U, CLAHE 65536 bins, back, Lab2BGR); clipLimit=3, tileGridSize=8",
          "source": f"cv2 {cv2.__version__}", "input": sref["input"], "out": rec})
    same = rec["sha256"] == sref["out"]["sha256"]
    print(f"cv2 rebuild vs verbatim Siril oracle: {'identical' if same else 'DIFFERENT'}")
    if not same:
        sys.exit("cv2 rebuild of the Siril pipeline disagrees with the oracle")

    # --- JPEG quantization tables (EncodeJPEG quality) --------------------
    # libjpeg-turbo scales the Annex K tables with the IJG quality rule
    # (jpeg_quality_scaling, baseline: clamped to 1..255); image/jpeg uses the
    # same rule, so EncodeJPEG at quality q must write these tables.
    yy, xx = np.mgrid[0:16, 0:16]
    bgr = np.dstack([np.full_like(xx, 64), yy * 16, xx * 16]).astype(np.uint8)
    tables = []
    for q in (1, 5, 10, 25, 49, 50, 51, 75, 90, 95, 100):
        ok, enc = cv2.imencode(".jpg", bgr, [cv2.IMWRITE_JPEG_QUALITY, q])
        if not ok:
            sys.exit(f"imencode quality {q} failed")
        tables.append({"quality": q, "dqt": dqt_tables(bytes(enc.reshape(-1)))})
    if len({json.dumps(t["dqt"]) for t in tables}) != len(tables):
        sys.exit("two qualities gave the same tables")
    dump(os.path.join(out, "jpeg_dqt_cv2.json"),
         {"source": f"cv2 {cv2.__version__} imencode .jpg IMWRITE_JPEG_QUALITY (libjpeg-turbo)",
          "note": "DQT tables in file order, values in zigzag order", "tables": tables})
    print("Done")


if __name__ == "__main__":
    main()
