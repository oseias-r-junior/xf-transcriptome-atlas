"""
fig_upset.py — UpSet plot of gene-expression intersections across conditions
(Figure 2C).

Each condition is binarised (gene expressed if TPM ≥ --tpm-min, default 0).
Vertical bars show intersection size; horizontal bars show total expressed
genes per condition.  Bars are coloured by strain membership:
  single-strain intersection → strain colour
  mixed intersection         → grey

Usage
-----
python fig_upset.py \\
    --tpm          data/tpm_expression.csv \\
    --metadata     data/sample_info.csv \\
    --condition-col condition \\
    --tpm-min      0 \\
    --sort-by      degree \\
    --output       figures/fig_2C_upset.tiff

Arguments
---------
--tpm-min       Minimum TPM to count a gene as expressed (default: 0,
                i.e., TPM > 0).
--strain-col    Metadata column that identifies the strain (default: strain).
--condition-col Metadata column used to define the sets (default: condition).
--sort-by       UpSet sort order: 'degree' or 'cardinality' (default: degree).

Input format
------------
TPM matrix (CSV, genes × samples):
    gene_id     s1   s2   ...
    Xf9a_00001  45   48   ...

Sample metadata (CSV, samples × traits):
    sample_id        condition      strain
    9a5c_PIM6_1d_1   9a5c_PIM6_1d   9a5c
    ...

The script averages replicates within each condition before binarising.
"""

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

try:
    from upsetplot import from_contents, UpSet
except ImportError:
    sys.exit("[ERROR] upsetplot is not installed. Run: pip install upsetplot")


# ---------------------------------------------------------------------------
# Colour palette
# ---------------------------------------------------------------------------

# Default two-strain palette; extended automatically for >2 strains
_DEFAULT_STRAIN_COLORS = [
    "#3399FF",   # 9a5c  (blue)
    "#FF6666",   # Temecula1 (red)
    "#59a86f",   # extra strains
    "#f58231",
    "#911eb4",
]
MIXED_COLOR = "#888888"


def build_strain_palette(strains: list[str]) -> dict[str, str]:
    return {s: _DEFAULT_STRAIN_COLORS[i % len(_DEFAULT_STRAIN_COLORS)]
            for i, s in enumerate(sorted(strains))}


# ---------------------------------------------------------------------------
# Data loading and binarisation
# ---------------------------------------------------------------------------

def load_data(tpm_path: Path, meta_path: Path) -> tuple:
    sep_t = "\t" if tpm_path.suffix in {".tsv", ".txt"} else ","
    sep_m = "\t" if meta_path.suffix in {".tsv", ".txt"} else ","
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sep_t)
    meta = pd.read_csv(meta_path, index_col=0, sep=sep_m)
    shared = tpm.columns.intersection(meta.index)
    return tpm[shared], meta.loc[shared]


def make_binary_per_condition(
    tpm: pd.DataFrame,
    meta: pd.DataFrame,
    condition_col: str,
    tpm_min: float,
) -> pd.DataFrame:
    """Return bool DataFrame (genes × conditions); True = expressed in condition."""
    result: dict[str, pd.Series] = {}
    for cond, rows in meta.groupby(condition_col):
        sids = rows.index.intersection(tpm.columns)
        mean_tpm = tpm[sids].mean(axis=1)
        result[cond] = mean_tpm > tpm_min
    return pd.DataFrame(result)


# ---------------------------------------------------------------------------
# Intersection colour logic
# ---------------------------------------------------------------------------

def classify_intersection(
    active_conditions: list[str],
    condition_to_strain: dict[str, str],
    strain_colors: dict[str, str],
) -> str:
    strains = {condition_to_strain.get(c, "unknown") for c in active_conditions}
    if len(strains) == 1:
        return strain_colors.get(strains.pop(), MIXED_COLOR)
    return MIXED_COLOR


def sorted_intersection_colors(
    upset_data,
    condition_to_strain: dict[str, str],
    strain_colors: dict[str, str],
    sort_by: str,
) -> list[str]:
    """
    Return a list of colours in the same left-to-right order that UpSet will
    draw the vertical bars.
    """
    cat_names = list(upset_data.index.names)
    df = upset_data.reset_index()
    df["_degree"] = df[cat_names].astype(int).sum(axis=1)
    count_col = [c for c in df.columns if c not in cat_names + ["_degree"]][0]

    if sort_by == "degree":
        df_sorted = df.sort_values(["_degree", count_col], ascending=[True, False])
    else:  # cardinality
        df_sorted = df.sort_values(count_col, ascending=False)

    colors = []
    for _, row in df_sorted.iterrows():
        active = [c for c in cat_names if row[c]]
        colors.append(classify_intersection(active, condition_to_strain, strain_colors))
    return colors


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_upset(
    bool_df: pd.DataFrame,
    condition_to_strain: dict[str, str],
    strain_colors: dict[str, str],
    sort_by: str,
    output: Path,
):
    # Build upsetplot data structure
    data_for_plot = {
        cond: set(bool_df.index[bool_df[cond]])
        for cond in bool_df.columns
    }
    upset_data = from_contents(data_for_plot)

    n_intersections = len(upset_data)
    n_conditions    = len(bool_df.columns)
    element_size    = max(25, min(45, 600 // max(n_intersections, 1)))

    upset = UpSet(
        upset_data,
        sort_by=sort_by,
        show_counts=True,
        element_size=element_size,
        totals_plot_elements=3,
    )

    fig = plt.figure(figsize=(max(12, 0.35 * n_intersections + 4), 6))
    axes_dict = upset.plot(fig)

    # --- Colour vertical (intersection) bars ---------------------------------
    inter_ax = axes_dict["intersections"]
    bar_colors = sorted_intersection_colors(
        upset_data, condition_to_strain, strain_colors, sort_by
    )
    for bar, color in zip(inter_ax.patches, bar_colors):
        bar.set_facecolor(color)
        bar.set_edgecolor("white")
        bar.set_linewidth(0.5)

    inter_ax.set_ylabel("Intersection size", fontsize=10)
    inter_ax.spines["top"].set_visible(False)
    inter_ax.spines["right"].set_visible(False)

    # --- Colour horizontal (totals) bars -------------------------------------
    totals_ax = axes_dict["totals"]
    # Condition order: read from matrix y-axis labels (bottom → top in figure)
    matrix_ax = axes_dict["matrix"]
    condition_order = [t.get_text() for t in matrix_ax.get_yticklabels()]

    # Totals patches are in same order as condition_order (bottom to top → reversed for patches)
    total_patches = totals_ax.patches
    if len(total_patches) == len(condition_order):
        for bar, cond in zip(total_patches, condition_order):
            strain = condition_to_strain.get(cond, "unknown")
            bar.set_facecolor(strain_colors.get(strain, MIXED_COLOR))
            bar.set_edgecolor("white")
            bar.set_linewidth(0.5)

    totals_ax.set_xlabel("Set size", fontsize=10)

    # --- Legend --------------------------------------------------------------
    handles = [
        mpatches.Patch(facecolor=color, label=strain)
        for strain, color in sorted(strain_colors.items())
    ]
    handles.append(mpatches.Patch(facecolor=MIXED_COLOR, label="mixed"))
    inter_ax.legend(
        handles=handles,
        title="Strain",
        fontsize=8,
        title_fontsize=8,
        loc="upper right",
        frameon=True,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] UpSet plot saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="UpSet plot of expressed-gene intersections across conditions.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm",           required=True, type=Path)
    p.add_argument("--metadata",      required=True, type=Path)
    p.add_argument("--condition-col", default="condition",
                   help="Metadata column defining the sets (x-axis groups).")
    p.add_argument("--strain-col",    default="strain",
                   help="Metadata column for strain identity (used for colouring).")
    p.add_argument("--tpm-min",       type=float, default=0.0,
                   help="Threshold: gene expressed if mean TPM > this value.")
    p.add_argument("--sort-by",       choices=["degree", "cardinality"],
                   default="degree")
    p.add_argument("--output",        required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm, meta = load_data(args.tpm, args.metadata)
    print(f"[INFO] {tpm.shape[0]} genes × {tpm.shape[1]} samples")

    if args.condition_col not in meta.columns:
        sys.exit(f"[ERROR] --condition-col '{args.condition_col}' not in metadata.")
    if args.strain_col not in meta.columns:
        sys.exit(f"[ERROR] --strain-col '{args.strain_col}' not in metadata.")

    # Condition → strain mapping (majority vote per condition)
    condition_to_strain: dict[str, str] = (
        meta.groupby(args.condition_col)[args.strain_col]
        .agg(lambda x: x.mode().iloc[0])
        .to_dict()
    )
    unique_strains = sorted(set(condition_to_strain.values()))
    strain_colors  = build_strain_palette(unique_strains)
    print(f"[INFO] {len(condition_to_strain)} conditions | "
          f"strains: {', '.join(unique_strains)}")

    bool_df = make_binary_per_condition(
        tpm, meta, args.condition_col, args.tpm_min
    )
    n_expressed = {c: bool_df[c].sum() for c in bool_df.columns}
    for c, n in n_expressed.items():
        print(f"  {c}: {n} expressed genes")

    plot_upset(bool_df, condition_to_strain, strain_colors, args.sort_by, args.output)


if __name__ == "__main__":
    main()
