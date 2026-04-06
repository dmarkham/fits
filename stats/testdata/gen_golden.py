#!/usr/bin/env python3
"""Generate golden test fixtures for fits/stats package.

Calls numpy/scipy/astropy to produce reference values for every stats
function. The Go test suite loads these and diffs against the pure-Go
implementation.

Usage:
    python gen_golden.py [output_dir]
"""

import json
import os
import sys

import numpy as np
from scipy import stats as scipy_stats

OUTDIR = sys.argv[1] if len(sys.argv) > 1 else "golden"
os.makedirs(OUTDIR, exist_ok=True)


def to_list(arr):
    """Convert numpy array to JSON-safe list."""
    return [float(x) if np.isfinite(x) else None for x in arr.flat]


def gen_test_arrays():
    """Generate a set of test arrays with various edge cases."""
    rng = np.random.RandomState(42)
    cases = {}

    # Standard float32 arrays
    cases["float32_random"] = rng.uniform(0, 1, 500).astype(np.float32)
    cases["float32_gradient"] = np.linspace(0.01, 0.3, 200).astype(np.float32)
    cases["float32_with_nan"] = np.array(
        [0.1, 0.2, float("nan"), 0.4, 0.5, float("nan"), 0.3, 0.15],
        dtype=np.float32,
    )
    cases["float32_with_zeros"] = np.array(
        [0.0, 0.1, 0.0, 0.2, 0.0, 0.3, 0.15, 0.25, 0.0, 0.05],
        dtype=np.float32,
    )
    cases["float32_all_same"] = np.full(100, 0.42, dtype=np.float32)
    cases["float32_single"] = np.array([3.14], dtype=np.float32)
    cases["float32_empty"] = np.array([], dtype=np.float32)
    cases["float32_all_nan"] = np.array(
        [float("nan")] * 5, dtype=np.float32
    )
    cases["float32_large"] = rng.normal(100, 15, 100000).astype(np.float32)
    cases["float32_negative"] = rng.uniform(-1, 1, 300).astype(np.float32)

    # Integer arrays
    cases["int16_random"] = rng.randint(0, 30000, 500).astype(np.int16)
    cases["uint8_random"] = rng.randint(0, 255, 500).astype(np.uint8)
    cases["int32_negative"] = rng.randint(-1000, 1000, 200).astype(np.int32)

    return cases


def compute_stats(data):
    """Compute all stats for a single array, NaN-aware."""
    result = {}
    result["dtype"] = str(data.dtype)
    result["len"] = len(data)
    result["input"] = to_list(data)

    # Filter NaN for computation
    if np.issubdtype(data.dtype, np.floating):
        clean = data[~np.isnan(data)]
    else:
        clean = data

    if len(clean) == 0:
        result["min"] = None
        result["max"] = None
        result["mean"] = None
        result["stdev"] = None
        result["median"] = None
        result["mad"] = None
        result["percentiles"] = {}
        return result

    result["min"] = float(clean.min())
    result["max"] = float(clean.max())
    result["mean"] = float(np.mean(clean.astype(np.float64)))
    result["stdev"] = float(np.std(clean.astype(np.float64), ddof=0))
    result["median"] = float(np.median(clean))
    # scipy MAD with scale=1 gives raw MAD (no 1.4826)
    result["mad"] = float(scipy_stats.median_abs_deviation(clean, scale=1))

    # Percentiles
    result["percentiles"] = {}
    for p in [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99]:
        result["percentiles"][str(p)] = float(np.percentile(clean, p * 100))

    return result


def compute_sigma_clip(data, k_low, k_high, max_iter, use_median):
    """Compute sigma-clipped stats matching astropy convention."""
    if np.issubdtype(data.dtype, np.floating):
        clean = data[~np.isnan(data)].copy()
    else:
        clean = data.astype(np.float64).copy()

    if len(clean) == 0:
        return {"mean": None, "stdev": None, "median": None, "ngood": 0}

    for _ in range(max_iter if max_iter > 0 else 100):
        if use_median:
            center = np.median(clean)
        else:
            center = np.mean(clean)
        stdev = np.std(clean, ddof=0)
        if stdev == 0:
            break
        mask = (clean >= center - k_low * stdev) & (
            clean <= center + k_high * stdev
        )
        if mask.sum() == len(clean):
            break
        clean = clean[mask]

    return {
        "mean": float(np.mean(clean)),
        "stdev": float(np.std(clean, ddof=0)),
        "median": float(np.median(clean)),
        "ngood": int(len(clean)),
    }


def compute_histogram(data, nbins):
    """Compute histogram matching numpy."""
    if np.issubdtype(data.dtype, np.floating):
        clean = data[~np.isnan(data)]
    else:
        clean = data.astype(np.float64)

    if len(clean) == 0:
        return {"counts": [], "edges": [], "nbins": nbins}

    counts, edges = np.histogram(clean, bins=nbins)
    return {
        "counts": counts.tolist(),
        "edges": edges.tolist(),
        "nbins": nbins,
    }


def main():
    print(f"Generating stats golden fixtures")
    cases = gen_test_arrays()

    all_golden = {}
    for name, data in cases.items():
        print(f"  {name}: len={len(data)} dtype={data.dtype}")
        entry = compute_stats(data)

        # Sigma clipping (only for arrays with enough data)
        if len(data) >= 10:
            entry["sigma_clip"] = {}
            for k in [2.0, 3.0, 5.0]:
                for center_name, use_med in [("mean", False), ("median", True)]:
                    key = f"k{k}_{center_name}"
                    entry["sigma_clip"][key] = compute_sigma_clip(
                        data, k, k, 10, use_med
                    )

        # Histogram
        if len(data) > 0:
            nbins = min(256, max(10, len(data) // 5))
            entry["histogram"] = compute_histogram(data, nbins)

        all_golden[name] = entry

    outpath = os.path.join(OUTDIR, "stats_golden.json")
    with open(outpath, "w") as f:
        json.dump(all_golden, f, indent=2)
    print(f"Written to {outpath}")


if __name__ == "__main__":
    main()
