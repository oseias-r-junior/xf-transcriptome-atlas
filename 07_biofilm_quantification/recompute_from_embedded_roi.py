"""
recompute_from_embedded_roi.py — Re-measure biofilm-ring area directly from
the FIJI/masks/*.tif files, using the EXACT per-photograph ROI that is
embedded inside each TIFF's ImageJ metadata (tag 50839, "IJMetadata"), and a
per-image Otsu threshold instead of one fixed global threshold.

Why this exists
-----------------
The original per-flask ROI selections were drawn interactively in Fiji and
were not exported to separate .roi files -- but ImageJ had already saved
them INSIDE each mask TIFF's metadata when the file was written (this is
the same mechanism as "Image > Overlay > To ROI Manager" / an active
selection saved with the image). This script decodes that embedded ROI
(standard ImageJ .roi binary format: "Iout" magic, header, then n x/n y
16-bit coordinates) and rasterizes it, so the exact same region the
original manual measurement used can be re-measured programmatically --
with no risk of drifting onto flask labels, volume-graduation marks, or
background, because we're reusing the real ROI polygon, not a guessed crop.

Sanity check: recomputing Area/Mean from the embedded ROI with the SAME
70-145 threshold as the original manual measurement reproduces
data/biofilm_ring_area_results.csv to within ~0.3% for 21/22 images (max
4.4% for one), confirming the decode is correct.

Why Otsu instead of a fixed 70-145 window
-------------------------------------------
A visual QC of all 22 photographs (diagnose_threshold_coverage.py) showed a
consistent "near-miss" halo of ring/foam texture just below the fixed
lower bound of 70 in nearly every image -- not a few outliers, a systematic
pattern. Computing Otsu's threshold separately within each photo's own ROI
pixel population confirms why: the optimal foreground/background split
varies from 37 to 74 across the 22 photos (mean 56.7, SD 8.9), i.e.
exposure/lighting was not perfectly identical between shots. Using each
photo's own Otsu threshold (rather than one global constant) adapts to
this while remaining a standard, parameter-free, non-arbitrary criterion.
No upper bound is applied (the original 145 cap was not found to exclude
any saturation artifacts -- per-ROI maxima across the dataset top out at
216/255 -- so it was dropped rather than kept as a second arbitrary
constant).

Usage
-----
# --masks-dir defaults to 07_biofilm_quantification/FIJI/masks (see that
# folder's README.md -- these images are not distributed with this repo,
# you must place your own there first):
python recompute_from_embedded_roi.py \\
    --output    data/biofilm_ring_area_results_otsu.csv

# Or point at a different copy of the masks:
python recompute_from_embedded_roi.py \\
    --masks-dir /path/to/your/masks \\
    --output    data/biofilm_ring_area_results_otsu.csv
"""
import argparse
import sys
from pathlib import Path
import struct
import glob

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw


def decode_embedded_roi(tif_path: Path):
    """Returns (PIL Image, list of (x,y) polygon vertices) or (im, None) if
    no ROI is embedded."""
    im = Image.open(tif_path)
    meta = im.tag_v2.get(50839) if hasattr(im, "tag_v2") else None
    if meta is None:
        return im, None
    idx = meta.find(b"Iout")
    if idx == -1:
        return im, None
    roi = meta[idx:]
    if len(roi) < 18:
        return im, None
    magic, version, rtype, _pad, top, left, bottom, right, n = struct.unpack(">4sh2b5h", roi[0:18])
    if n == 0 or len(roi) < 64 + 4 * n:
        return im, None
    xs = struct.unpack(">%dh" % n, roi[64:64 + 2 * n])
    ys = struct.unpack(">%dh" % n, roi[64 + 2 * n:64 + 4 * n])
    pts = [(left + x, top + y) for x, y in zip(xs, ys)]
    return im, pts


def rasterize_roi(im: Image.Image, pts) -> np.ndarray:
    mask_img = Image.new("L", im.size, 0)
    ImageDraw.Draw(mask_img).polygon(pts, outline=255, fill=255)
    return np.array(mask_img) > 0


def otsu_threshold(values: np.ndarray) -> int:
    hist, _ = np.histogram(values, bins=256, range=(0, 256))
    hist = hist.astype(float)
    total = hist.sum()
    sum_all = np.dot(np.arange(256), hist)
    sumB = wB = max_var = 0.0
    thresh = 0
    for t in range(256):
        wB += hist[t]
        if wB == 0:
            continue
        wF = total - wB
        if wF == 0:
            break
        sumB += t * hist[t]
        mB = sumB / wB
        mF = (sum_all - sumB) / wF
        var_between = wB * wF * (mB - mF) ** 2
        if var_between > max_var:
            max_var = var_between
            thresh = t
    return thresh


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Recompute biofilm-ring Area/Mean from embedded Fiji ROIs with per-image Otsu threshold.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--masks-dir", type=Path,
                    default=Path("07_biofilm_quantification/FIJI/masks"),
                    help="Directory with 8-bit grayscale flask-photo TIFFs, ROI "
                        "embedded in each file's metadata (see "
                        "07_biofilm_quantification/FIJI/README.md). These are "
                        "the authors' own manual measurements and are NOT "
                        "distributed with this repository -- place your own "
                        "here, or pass a different --masks-dir.")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--sanity-check-against", type=Path, default=None,
                    help="Optional original Results.csv (70-145 threshold) to print a "
                        "side-by-side reproduction check against.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tif_paths = sorted(args.masks_dir.glob("*.tif"))
    print(f"[INFO] Found {len(tif_paths)} images in {args.masks_dir}")
    if not tif_paths:
        sys.exit(
            f"[ERROR] No .tif files found in {args.masks_dir}. These images "
            f"are the authors' own manual measurements and are not "
            f"distributed with this repository -- see "
            f"07_biofilm_quantification/FIJI/README.md for what to place "
            f"there, or pass --masks-dir pointing at your own copy."
        )

    orig = None
    if args.sanity_check_against and args.sanity_check_against.exists():
        orig = pd.read_csv(args.sanity_check_against).set_index("Label")

    rows = []
    for tif_path in tif_paths:
        label = tif_path.stem
        im, pts = decode_embedded_roi(tif_path)
        if pts is None:
            print(f"[WARN] No embedded ROI found in {tif_path.name}, skipping.")
            continue
        arr = np.array(im)
        roi_mask = rasterize_roi(im, pts)
        roi_pixels = arr[roi_mask]

        otsu_t = otsu_threshold(roi_pixels)
        counted = roi_pixels >= otsu_t
        area = int(counted.sum())
        mean_gray = float(roi_pixels[counted].mean()) if area > 0 else float("nan")

        row = {"Label": label, "Area": area, "Mean": round(mean_gray, 3),
               "OtsuThr": otsu_t, "n_roi_px": int(roi_mask.sum())}

        if orig is not None and label in orig.index:
            orig_area = orig.loc[label, "Area"]
            check_counted = (roi_pixels >= 70) & (roi_pixels <= 145)
            reproduced_area = int(check_counted.sum())
            row["orig_area_70_145"] = orig_area
            row["reproduced_area_70_145"] = reproduced_area
            row["reproduction_pct_diff"] = round(100 * (reproduced_area - orig_area) / orig_area, 2)

        rows.append(row)
        print(f"[INFO] {label}: Otsu threshold={otsu_t}, Area={area} px, Mean={mean_gray:.2f}")

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"\n[INFO] Saved to {args.output}")
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
