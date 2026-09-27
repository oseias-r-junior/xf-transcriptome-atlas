"""
measure_biofilm_ring.py — Demonstration script for the quantitative image-
analysis pipeline used to measure Xylella fastidiosa biofilm-ring formation
on the glass wall of Erlenmeyer flasks (air-liquid interface).

This is a Python re-implementation of the Fiji/ImageJ workflow actually used
to generate the values in `data/biofilm_ring_area_results.csv` (exported
from ImageJ's Results window: Label, Area, Mean, MinThr, MaxThr columns).
It exists so the measurement is scriptable/reproducible outside Fiji and so
each step of the pipeline can be inspected and re-run on new photographs. It
does NOT re-analyse the flask photographs (those are not part of this
repository), it reproduces the exact processing chain that was applied to
them, and ships a `--demo` mode that runs the full pipeline on a synthetic
test image so the script can be verified without any input photo.

Pipeline (mirrors the ImageJ tutorial step for step)
-----------------------------------------------------
  1. Load the flask photograph and convert to 8-bit grayscale
     (Image > Type > 8-bit).
  2. Subtract uneven background caused by the curved glass
     (Process > Subtract Background, rolling-ball radius = 50 px by
     default). OpenCV has no native rolling-ball filter, so it is
     approximated here with a greyscale morphological opening using a
     disk-shaped structuring element of the same radius -- the standard,
     documented approximation to ImageJ's rolling-ball algorithm
     (background = grayscale opening; result = original - background).
  3. Restrict the measurement to a user-defined ROI covering only the
     air-liquid interface band on the glass wall (excludes the beads at
     the bottom of the flask and the empty headspace at the top) --
     equivalent to Analyze > Tools > ROI Manager > Add.
  4. Threshold the ROI with the same two-value (min, max) intensity window
     used for every image of a given trial (Image > Adjust > Threshold),
     default 70-145 as used for this dataset.
  5. Measure, with "Limit to threshold" behaviour: Area = number of
     thresholded pixels inside the ROI; Mean = mean gray value of only
     those thresholded pixels (Analyze > Set Measurements > Area, Mean
     gray value, Limit to threshold, Display label).
  6. Save a QC figure (original / grayscale / background-subtracted /
     thresholded-mask overlay) and append the measurement to a
     Results-style CSV with the same columns ImageJ exports
     (Label, Area, Mean, MinThr, MaxThr), so it is a drop-in match for
     `data/biofilm_ring_area_results.csv`.

Usage
-----
# Demo on a synthetic flask image (no input file needed):
python measure_biofilm_ring.py --demo --output-dir demo_output

# Real photograph, ROI in pixels (x, y, width, height) around the ring band:
python measure_biofilm_ring.py \\
    --image      photos/9a5c_PIM6_late_1.jpg \\
    --label      9a_pim6_late_1 \\
    --roi        180 220 260 140 \\
    --min-thr    70 --max-thr 145 \\
    --rolling-ball-radius 50 \\
    --output-dir figures/biofilm_qc \\
    --results-csv data/biofilm_ring_area_results.csv
"""
import argparse
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def to_8bit_grayscale(img_bgr: np.ndarray) -> np.ndarray:
    """Image > Type > 8-bit."""
    gray = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2GRAY)
    return gray.astype(np.uint8)


def subtract_background(gray: np.ndarray, rolling_ball_radius: int = 50) -> np.ndarray:
    """Process > Subtract Background (rolling ball), approximated with a
    greyscale morphological opening using a disk structuring element of
    the same radius -- background = opening(gray); result = gray - background,
    clipped to [0, 255]."""
    ksize = 2 * rolling_ball_radius + 1
    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (ksize, ksize))
    background = cv2.morphologyEx(gray, cv2.MORPH_OPEN, kernel)
    subtracted = cv2.subtract(gray, background)
    return subtracted


def apply_roi(img: np.ndarray, roi):
    """roi = (x, y, w, h) in pixels. Returns the cropped sub-image."""
    if roi is None:
        return img
    x, y, w, h = roi
    return img[y:y + h, x:x + w]


def threshold_and_measure(gray_roi: np.ndarray, min_thr: int, max_thr: int):
    """Image > Adjust > Threshold + Analyze > Measure with 'Limit to
    threshold'. Returns (binary_mask, area_px, mean_gray_in_mask)."""
    mask = cv2.inRange(gray_roi, min_thr, max_thr)
    area_px = int(np.count_nonzero(mask))
    if area_px > 0:
        mean_gray = float(gray_roi[mask > 0].mean())
    else:
        mean_gray = float("nan")
    return mask, area_px, mean_gray


def make_qc_figure(img_bgr, gray, subtracted, roi, mask_full, label, area_px,
                    mean_gray, min_thr, max_thr, output_path: Path):
    """Four-panel QC figure: original / grayscale / background-subtracted /
    thresholded mask overlay (red), mirroring what you'd see in Fiji."""
    img_rgb = cv2.cvtColor(img_bgr, cv2.COLOR_BGR2RGB)

    overlay = cv2.cvtColor(subtracted, cv2.COLOR_GRAY2RGB).astype(np.uint8)
    if roi is not None:
        x, y, w, h = roi
        mask_canvas = np.zeros(subtracted.shape, dtype=np.uint8)
        mask_canvas[y:y + h, x:x + w] = mask_full
        cv2.rectangle(overlay, (x, y), (x + w, y + h), (0, 255, 0), 3)
    else:
        mask_canvas = mask_full
    red_layer = overlay.copy()
    red_layer[mask_canvas > 0] = [255, 0, 0]
    overlay = cv2.addWeighted(overlay, 0.5, red_layer, 0.5, 0)

    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
    axes[0].imshow(img_rgb); axes[0].set_title("1. Original (8-bit color)")
    axes[1].imshow(gray, cmap="gray"); axes[1].set_title("2. 8-bit grayscale")
    axes[2].imshow(subtracted, cmap="gray"); axes[2].set_title("3. Background-subtracted\n(ROI in green)")
    axes[3].imshow(overlay); axes[3].set_title(
        f"4. Threshold {min_thr}-{max_thr} (red)\nArea = {area_px} px, Mean = {mean_gray:.1f}"
    )
    for ax in axes:
        ax.axis("off")
    fig.suptitle(f"Biofilm ring quantification — {label}", fontsize=13)
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] QC figure saved to {output_path}")


def make_demo_image(width=500, height=350, seed=0) -> np.ndarray:
    """Synthetic Erlenmeyer-flask-wall photo for --demo: a glass-coloured
    background with an uneven lighting gradient (mimicking curved-glass
    reflections) and a brighter, textured horizontal band representing the
    biofilm ring at the air-liquid interface."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:height, 0:width]
    # Uneven background: radial-ish gradient from a bright reflection spot
    gradient = 120 + 40 * np.exp(-((xx - width * 0.25) ** 2 + (yy - height * 0.3) ** 2) / (2 * 120 ** 2))
    gradient += 20 * (yy / height)
    background = gradient + rng.normal(0, 4, size=(height, width))

    ring_center = height * 0.55
    ring_half_width = height * 0.06
    ring_mask = np.abs(yy - ring_center) < ring_half_width
    ring_signal = np.zeros((height, width))
    ring_signal[ring_mask] = 55 + rng.normal(0, 12, size=ring_mask.sum())

    gray = np.clip(background + ring_signal, 0, 255).astype(np.uint8)
    img_bgr = cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR)
    return img_bgr


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Quantify the X. fastidiosa biofilm ring on an Erlenmeyer flask wall "
                    "from a photograph, reproducing the Fiji/ImageJ workflow.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--image", type=Path, help="Input flask photograph (any format OpenCV can read).")
    p.add_argument("--demo", action="store_true",
                    help="Run the full pipeline on a synthetic demo image instead of --image.")
    p.add_argument("--label", type=str, default=None,
                    help="Sample label to record (e.g. '9a_pim6_late_1'). Defaults to the image filename stem.")
    p.add_argument("--roi", type=int, nargs=4, metavar=("X", "Y", "W", "H"), default=None,
                    help="Rectangular ROI (pixels) restricting the measurement to the ring band "
                        "on the flask wall. If omitted, the whole image is used.")
    p.add_argument("--min-thr", type=int, default=70, help="Lower threshold bound.")
    p.add_argument("--max-thr", type=int, default=145, help="Upper threshold bound.")
    p.add_argument("--rolling-ball-radius", type=int, default=50,
                    help="Rolling-ball radius (px) for background subtraction.")
    p.add_argument("--output-dir", type=Path, default=Path("figures/biofilm_qc"),
                    help="Directory for the QC figure.")
    p.add_argument("--results-csv", type=Path, default=None,
                    help="If given, append/update this row in a Results-style CSV "
                        "(Label, Area, Mean, MinThr, MaxThr) matching the ImageJ export format.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.demo:
        img_bgr = make_demo_image()
        label = args.label or "demo_flask"
        roi = args.roi or (40, 150, 420, 80)  # covers the synthetic ring band
    else:
        if args.image is None:
            raise SystemExit("Provide --image <path> or use --demo.")
        img_bgr = cv2.imread(str(args.image))
        if img_bgr is None:
            raise SystemExit(f"Could not read image: {args.image}")
        label = args.label or args.image.stem
        roi = tuple(args.roi) if args.roi else None

    gray = to_8bit_grayscale(img_bgr)
    subtracted = subtract_background(gray, args.rolling_ball_radius)
    roi_img = apply_roi(subtracted, roi)
    mask, area_px, mean_gray = threshold_and_measure(roi_img, args.min_thr, args.max_thr)

    print(f"[INFO] {label}: Area = {area_px} px, Mean (in-threshold) = {mean_gray:.3f}, "
          f"MinThr = {args.min_thr}, MaxThr = {args.max_thr}")

    qc_path = args.output_dir / f"{label}_qc.png"
    make_qc_figure(img_bgr, gray, subtracted, roi, mask, label, area_px, mean_gray,
                    args.min_thr, args.max_thr, qc_path)

    if args.results_csv:
        row = pd.DataFrame([{
            "Label": label, "Area": area_px, "Mean": round(mean_gray, 3),
            "MinThr": args.min_thr, "MaxThr": args.max_thr,
        }])
        if args.results_csv.exists():
            existing = pd.read_csv(args.results_csv)
            existing = existing[existing["Label"] != label]
            out = pd.concat([existing, row], ignore_index=True)
        else:
            out = row
        args.results_csv.parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.results_csv, index=False)
        print(f"[INFO] Result appended to {args.results_csv}")


if __name__ == "__main__":
    main()
