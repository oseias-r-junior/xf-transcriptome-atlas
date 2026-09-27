"""
normalize_by_flask_width.py — Normalize ImageJ/Fiji biofilm-ring pixel areas
(data/biofilm_ring_area_results.csv) by the apparent flask width measured in
each photograph (data/flask_width_measurements.csv), to remove the
photographic zoom/camera-distance confound from raw pixel-area comparisons.

Why this is needed
-------------------
All flasks are physically identical 250 mL Erlenmeyers, but photos were taken
at slightly different camera distances/zoom levels, so the same physical ring
projects to a different pixel area in different photos. To correct for this,
a straight-line measurement of the flask's outer width (same row height as
the ring band, "Analyze > Measure" with only Length checked) was taken for
every photo in Fiji. If a photo is scaled by some factor k relative to
another (pure zoom/distance difference, same physical flask), length scales
by k and area scales by k^2. Dividing each flask's ring area by the SQUARE of
its own flask-width measurement therefore yields a dimensionless,
scale-invariant quantity that is comparable across photos regardless of
camera distance:

    normalized_area = Area_px / (flask_width_px)^2

This is multiplied by a constant (1e6 by default) purely for readability
(so values aren't tiny decimals); the constant cancels out of all pairwise
comparisons and does not change any p-value.

Usage
-----
python normalize_by_flask_width.py \\
    --results      data/biofilm_ring_area_results.csv \\
    --flask-width  data/flask_width_measurements.csv \\
    --outdir       results/biofilm_quantification \\
    --scale-factor 1e6
"""
import argparse
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent.parent / "06_stats"))
sys.path.insert(0, str(Path(__file__).parent))
from pure_stats import t_sf_twotailed, benjamini_hochberg
from analyze_biofilm_areas import (
    label_strain, label_medium, label_phase, welch_t_test, STRAIN_MEDIUM_COLORS,
)


# Label spelling in the Fiji flask-width export doesn't always exactly match
# the Results.csv Label for the same flask (e.g. a missing underscore).
# Map explicit known aliases here rather than guessing with fuzzy matching.
LABEL_ALIASES = {
    "tem1_early_pwg3": "tem1_early_pwg_3",
}


def load_flask_widths(path: Path) -> pd.DataFrame:
    """Parses either a proper CSV (Label, Length columns) or a raw
    Fiji Results-window text export (Label like 'foo.jpg:flask_measruement',
    tab-separated, with a leading unnamed index column)."""
    text = path.read_text(encoding="utf-8")
    if "\t" in text.splitlines()[0]:
        rows = []
        for line in text.splitlines()[1:]:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\t+", line)
            if len(parts) < 6:
                continue
            _, raw_label, _area, _mean, _angle, length = parts[:6]
            label = raw_label.split(".jpg")[0].split(":")[0].strip()
            label = LABEL_ALIASES.get(label, label)
            rows.append({"Label": label, "flask_width_px": float(length)})
        return pd.DataFrame(rows)
    else:
        df = pd.read_csv(path)
        df.columns = [c.strip() for c in df.columns]
        df["Label"] = df["Label"].apply(lambda s: LABEL_ALIASES.get(str(s).strip(), str(s).strip()))
        return df[["Label", "flask_width_px" if "flask_width_px" in df.columns else "Length"]].rename(
            columns={"Length": "flask_width_px"}
        )


def pairwise_tests(df: pd.DataFrame, value_col: str) -> pd.DataFrame:
    rows = []
    for medium in ["PIM6", "PWG"]:
        for phase in ["early", "late"]:
            a = df[(df.strain == "9a5c") & (df.medium == medium) & (df.phase == phase)][value_col]
            b = df[(df.strain == "Temecula1") & (df.medium == medium) & (df.phase == phase)][value_col]
            if len(a) >= 2 and len(b) >= 2:
                t, dfree, p = welch_t_test(a, b)
                rows.append({
                    "contrast": "strain (9a5c vs Temecula1)", "medium": medium, "phase": phase,
                    "n_9a5c": len(a), "n_Temecula1": len(b),
                    "mean_9a5c": a.mean(), "mean_Temecula1": b.mean(),
                    "t": t, "df": dfree, "p_raw": p,
                })
    for strain in ["9a5c", "Temecula1"]:
        for medium in ["PIM6", "PWG"]:
            a = df[(df.strain == strain) & (df.medium == medium) & (df.phase == "early")][value_col]
            b = df[(df.strain == strain) & (df.medium == medium) & (df.phase == "late")][value_col]
            if len(a) >= 2 and len(b) >= 2:
                t, dfree, p = welch_t_test(a, b)
                rows.append({
                    "contrast": f"phase (early vs late), {strain}", "medium": medium, "phase": "-",
                    "n_early": len(a), "n_late": len(b),
                    "mean_early": a.mean(), "mean_late": b.mean(),
                    "t": t, "df": dfree, "p_raw": p,
                })
    res = pd.DataFrame(rows)
    res["p_BH"] = benjamini_hochberg(res["p_raw"].values)
    res["significant_BH_0.05"] = res["p_BH"] < 0.05
    return res


def plot_bars(desc: pd.DataFrame, value_label: str, output: Path):
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)
    phase_order = ["early", "late"]
    x = np.arange(len(phase_order))
    width = 0.18
    for ax, medium in zip(axes, ["PIM6", "PWG"]):
        offset = -1.5
        for strain in ["9a5c", "Temecula1"]:
            sub = desc[(desc.strain == strain) & (desc.medium == medium)].set_index("phase").reindex(phase_order)
            color = STRAIN_MEDIUM_COLORS[(strain, medium)]
            ax.bar(x + offset * width, sub["mean"], width, yerr=sub["se"],
                   color=color, capsize=4, label=strain, edgecolor="black", linewidth=0.5)
            offset += 1
        ax.set_xticks(x); ax.set_xticklabels(["Early", "Late"])
        ax.set_title(f"Medium: {medium}"); ax.set_xlabel("Growth phase")
        ax.legend(fontsize=9)
    axes[0].set_ylabel(value_label)
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Bar chart saved to {output}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Normalize biofilm ring area by flask width to remove photo zoom/distance confound.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--results", required=True, type=Path)
    p.add_argument("--flask-width", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--scale-factor", type=float, default=1e6,
                    help="Multiplier applied to Area/width^2 purely for readable numbers "
                        "(cancels out of every comparison/p-value).")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    results = pd.read_csv(args.results)
    results.columns = [c.strip() for c in results.columns]
    widths = load_flask_widths(args.flask_width)

    merged = results.merge(widths, on="Label", how="left", validate="one_to_one")
    missing = merged[merged["flask_width_px"].isna()]
    if len(missing):
        print(f"[WARN] {len(missing)} rows have no matching flask-width measurement:")
        print(missing[["Label"]].to_string(index=False))
    merged = merged.dropna(subset=["flask_width_px"])

    merged["strain"] = merged["Label"].map(label_strain)
    merged["medium"] = merged["Label"].map(label_medium)
    merged["phase"] = merged["Label"].map(label_phase)
    merged["normalized_area"] = args.scale_factor * merged["Area"] / (merged["flask_width_px"] ** 2)

    merged_out = args.outdir / "biofilm_area_normalized.csv"
    merged.round(4).to_csv(merged_out, index=False)
    print(f"[INFO] Per-flask normalized data saved to {merged_out}")
    print(merged[["Label", "Area", "flask_width_px", "normalized_area"]].round(2).to_string(index=False))

    desc = (
        merged.groupby(["strain", "medium", "phase"])["normalized_area"]
        .agg(mean="mean", sd="std", n="count").reset_index()
    )
    desc["se"] = desc["sd"] / np.sqrt(desc["n"])
    desc_out = args.outdir / "biofilm_area_normalized_descriptive_stats.csv"
    desc.round(4).to_csv(desc_out, index=False)
    print("\n[INFO] Descriptive stats (normalized):")
    print(desc.round(3).to_string(index=False))

    tests = pairwise_tests(merged, "normalized_area")
    tests_out = args.outdir / "biofilm_area_normalized_pairwise_tests.csv"
    tests.round(4).to_csv(tests_out, index=False)
    print("\n[INFO] Pairwise Welch t-tests on NORMALIZED area (BH-corrected across all tests):")
    print(tests.round(4).to_string(index=False))

    plot_bars(desc, f"Normalized ring area (Area / flask width² × {args.scale_factor:.0e})",
               args.outdir / "fig_biofilm_ring_area_normalized.tiff")


if __name__ == "__main__":
    main()
