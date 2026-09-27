"""
analyze_biofilm_areas.py — Descriptive statistics and pairwise comparisons
for the ImageJ/Fiji biofilm-ring area measurements
(data/biofilm_ring_area_results.csv), split by strain (9a5c / Temecula1),
medium (PIM6 / PWG) and growth phase (early / late).

Mirrors the dependency-light convention used throughout this repository
(06_stats/pure_stats.py): everything here runs on NumPy/pandas/matplotlib
only, no scipy/statsmodels required.

Because each strain x medium x phase cell has only 2-5 flasks, this script
does NOT run a single omnibus model over all 8 groups (too few df for a
reliable 3-way interaction with cells this small). Instead it runs targeted,
pre-specified Welch two-sample t-tests for exactly the two contrasts the
manuscript's Results/Discussion make claims about:
  1. strain effect (9a5c vs Temecula1), computed separately within each
     medium x phase cell (so it isn't confounded by medium or phase);
  2. phase effect (early vs late), computed separately within each
     strain x medium cell (so it isn't confounded by strain or medium).
All 8 resulting p-values are then jointly Benjamini-Hochberg-corrected
(multiple-comparison control across every test run), matching the FDR
convention already used for DESeq2 and GO enrichment in this pipeline.

Usage
-----
python analyze_biofilm_areas.py \\
    --results  data/biofilm_ring_area_results.csv \\
    --outdir   results/biofilm_quantification
"""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent.parent / "06_stats"))
from pure_stats import t_sf_twotailed, benjamini_hochberg


STRAIN_MEDIUM_COLORS = {
    ("9a5c", "PIM6"):      "#FF8C00",
    ("9a5c", "PWG"):       "#FF0000",
    ("Temecula1", "PIM6"): "#ADD8E6",
    ("Temecula1", "PWG"):  "#00008B",
}


def label_strain(label: str) -> str:
    return "9a5c" if label.lower().startswith("9a") else "Temecula1"


def label_medium(label: str) -> str:
    return "PIM6" if "pim6" in label.lower() else "PWG"


def label_phase(label: str) -> str:
    return "early" if "early" in label.lower() else "late"


def welch_t_test(a, b):
    """Welch's two-sample t-test (unequal variances), matching
    scipy.stats.ttest_ind(equal_var=False)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    na, nb = len(a), len(b)
    va, vb = a.var(ddof=1), b.var(ddof=1)
    se = np.sqrt(va / na + vb / nb)
    t = (a.mean() - b.mean()) / se
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    p = t_sf_twotailed(t, df)
    return t, df, p


def load_and_annotate(results_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(results_csv)
    df.columns = [c.strip() for c in df.columns]
    df["strain"] = df["Label"].map(label_strain)
    df["medium"] = df["Label"].map(label_medium)
    df["phase"] = df["Label"].map(label_phase)
    return df


def descriptive_table(df: pd.DataFrame) -> pd.DataFrame:
    g = (
        df.groupby(["strain", "medium", "phase"])["Area"]
        .agg(mean_area="mean", sd_area="std", n="count")
        .reset_index()
    )
    g["se_area"] = g["sd_area"] / np.sqrt(g["n"])
    return g


def pairwise_tests(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for medium in ["PIM6", "PWG"]:
        for phase in ["early", "late"]:
            a = df[(df.strain == "9a5c") & (df.medium == medium) & (df.phase == phase)]["Area"]
            b = df[(df.strain == "Temecula1") & (df.medium == medium) & (df.phase == phase)]["Area"]
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
            a = df[(df.strain == strain) & (df.medium == medium) & (df.phase == "early")]["Area"]
            b = df[(df.strain == strain) & (df.medium == medium) & (df.phase == "late")]["Area"]
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


def plot_bars(desc: pd.DataFrame, output: Path):
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
            ax.bar(x + offset * width, sub["mean_area"], width, yerr=sub["se_area"],
                   color=color, capsize=4, label=f"{strain}", edgecolor="black", linewidth=0.5)
            offset += 1
        ax.set_xticks(x)
        ax.set_xticklabels(["Early", "Late"])
        ax.set_title(f"Medium: {medium}")
        ax.set_xlabel("Growth phase")
        ax.legend(fontsize=9)
    axes[0].set_ylabel("Biofilm ring area (px, thresholded 70-145)")
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Bar chart saved to {output}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Descriptive stats + pairwise Welch t-tests for biofilm-ring ImageJ areas.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--results", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    df = load_and_annotate(args.results)
    args.outdir.mkdir(parents=True, exist_ok=True)

    desc = descriptive_table(df)
    desc_out = args.outdir / "biofilm_area_descriptive_stats.csv"
    desc.round(2).to_csv(desc_out, index=False)
    print("[INFO] Descriptive stats:")
    print(desc.round(1).to_string(index=False))
    print(f"[INFO] saved to {desc_out}")

    tests = pairwise_tests(df)
    tests_out = args.outdir / "biofilm_area_pairwise_tests.csv"
    tests.round(4).to_csv(tests_out, index=False)
    print("\n[INFO] Pairwise Welch t-tests (BH-corrected across all 8 tests):")
    print(tests.round(4).to_string(index=False))
    print(f"[INFO] saved to {tests_out}")

    plot_bars(desc, args.outdir / "fig_biofilm_ring_area.tiff")


if __name__ == "__main__":
    main()
