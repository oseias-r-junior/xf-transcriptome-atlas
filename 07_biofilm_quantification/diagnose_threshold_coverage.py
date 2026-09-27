"""
diagnose_threshold_coverage.py — Visual QC: for each of the 22 background-
subtracted flask images (FIJI/masks/*.tif, the exact images the 70-145
threshold was applied to), overlay which pixels fall INSIDE the current
threshold window (counted as biofilm, shown red) versus pixels just below it
(30-70 gray levels, shown yellow) that are plausibly ring/foam texture being
missed ("gray-to-black" pixels the user noticed are not marked red).

This does not change the threshold or re-run the quantification -- it is
purely a per-image visual check so a human can decide, image by image,
whether 70 is cutting into real biofilm texture or correctly excluding
background/glass reflections.

Usage
-----
# --masks-dir defaults to 07_biofilm_quantification/FIJI/masks (see that
# folder's README.md -- these images are not distributed with this repo,
# you must place your own there first):
python diagnose_threshold_coverage.py \\
    --outdir     results/biofilm_quantification/threshold_qc \\
    --min-thr 70 --max-thr 145 --near-miss-low 30
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt


def find_ring_band(arr: np.ndarray, low=30, high=200, min_frac=0.15):
    """Locate the vertical row-range most likely to contain the ring band:
    the band of rows with the highest density of pixels in [low, high]
    (the ring's foam/bubble texture is bright relative to the very dark
    background and the very dark clear-liquid glass interior)."""
    h, w = arr.shape
    band_mask = (arr >= low) & (arr <= high)
    row_frac = band_mask.mean(axis=1)
    # smooth
    k = 15
    kernel = np.ones(k) / k
    smoothed = np.convolve(row_frac, kernel, mode="same")
    peak_row = int(np.argmax(smoothed))
    # expand outward while density stays above min_frac (or up to 20% of image height)
    max_half = int(h * 0.12)
    top = peak_row
    while top > 0 and smoothed[top] > min_frac and (peak_row - top) < max_half:
        top -= 1
    bottom = peak_row
    while bottom < h - 1 and smoothed[bottom] > min_frac and (bottom - peak_row) < max_half:
        bottom += 1
    margin = int(0.02 * h)
    top = max(0, top - margin)
    bottom = min(h, bottom + margin)
    return top, bottom


def make_overlay(arr: np.ndarray, min_thr: int, max_thr: int, near_miss_low: int):
    rgb = np.stack([arr, arr, arr], axis=-1).astype(np.uint8)
    counted = (arr >= min_thr) & (arr <= max_thr)
    near_miss = (arr >= near_miss_low) & (arr < min_thr)
    rgb[counted] = [220, 30, 30]     # red: currently counted as biofilm
    rgb[near_miss] = [235, 200, 0]   # yellow: just below threshold, possibly missed ring texture
    return rgb, counted, near_miss


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Per-image visual QC of the 70-145 biofilm threshold vs. near-miss pixels.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--masks-dir", type=Path,
                    default=Path("07_biofilm_quantification/FIJI/masks"),
                    help="Directory of background-subtracted 8-bit .tif images. "
                        "These are the authors' own manual measurements and are "
                        "NOT distributed with this repository -- see "
                        "07_biofilm_quantification/FIJI/README.md.")
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-thr", type=int, default=70)
    p.add_argument("--max-thr", type=int, default=145)
    p.add_argument("--near-miss-low", type=int, default=30,
                    help="Lower bound of the 'near miss' (yellow) band just below min-thr.")
    p.add_argument("--grid-cols", type=int, default=4)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    tif_paths = sorted(args.masks_dir.glob("*.tif"))
    print(f"[INFO] Found {len(tif_paths)} mask images in {args.masks_dir}")
    if not tif_paths:
        sys.exit(
            f"[ERROR] No .tif files found in {args.masks_dir}. These images "
            f"are the authors' own manual measurements and are not "
            f"distributed with this repository -- see "
            f"07_biofilm_quantification/FIJI/README.md for what to place "
            f"there, or pass --masks-dir pointing at your own copy."
        )

    summary_rows = []
    crops = []
    for tif_path in tif_paths:
        label = tif_path.stem
        arr = np.array(Image.open(tif_path))
        top, bottom = find_ring_band(arr)
        strip = arr[top:bottom, :]
        rgb, counted, near_miss = make_overlay(strip, args.min_thr, args.max_thr, args.near_miss_low)

        counted_px = int(counted.sum())
        near_miss_px = int(near_miss.sum())
        pct_near_miss = 100 * near_miss_px / max(counted_px, 1)
        summary_rows.append({
            "Label": label, "counted_px_in_strip": counted_px,
            "near_miss_px_in_strip": near_miss_px, "near_miss_pct_of_counted": round(pct_near_miss, 1),
            "strip_row_range": f"{top}-{bottom}",
        })

        fig, ax = plt.subplots(figsize=(6, 6 * strip.shape[0] / strip.shape[1] + 0.4))
        ax.imshow(rgb)
        ax.set_title(f"{label}\nred=counted ({args.min_thr}-{args.max_thr})  "
                      f"yellow=near-miss ({args.near_miss_low}-{args.min_thr})\n"
                      f"near-miss = {pct_near_miss:.0f}% of counted area", fontsize=8)
        ax.axis("off")
        plt.tight_layout()
        indiv_path = args.outdir / f"{label}_threshold_qc.png"
        fig.savefig(indiv_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        crops.append((label, rgb, pct_near_miss))
        print(f"[INFO] {label}: near-miss = {pct_near_miss:.1f}% of counted area -> {indiv_path.name}")

    # Composite contact sheet
    import math
    n = len(crops)
    cols = args.grid_cols
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4, rows * 2.6))
    axes = np.atleast_2d(axes)
    for i, (label, rgb, pct) in enumerate(crops):
        r, c = divmod(i, cols)
        ax = axes[r][c]
        ax.imshow(rgb)
        ax.set_title(f"{label}\nnear-miss {pct:.0f}%", fontsize=7)
        ax.axis("off")
    for j in range(n, rows * cols):
        r, c = divmod(j, cols)
        axes[r][c].axis("off")
    plt.tight_layout()
    contact_path = args.outdir / "threshold_qc_contact_sheet.png"
    fig.savefig(contact_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Contact sheet saved to {contact_path}")

    import pandas as pd
    summary = pd.DataFrame(summary_rows).sort_values("near_miss_pct_of_counted", ascending=False)
    summary_path = args.outdir / "threshold_qc_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[INFO] Summary saved to {summary_path}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
