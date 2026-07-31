"""
fig_pearson_heatmap.py — Hierarchically clustered Pearson correlation heatmap
of condition-level expression profiles (Figure S1).

Two input modes:
  (a) Pre-computed correlation matrix (CSV) — use when the values come from a
      supplementary table or were computed externally.
  (b) TPM matrix + sample metadata — the script averages replicates per
      condition, log2-transforms, and computes Pearson r between conditions.

Usage
-----
# From a pre-computed correlation CSV (recommended for reproducibility with
# the values published in the supplementary material):
python fig_pearson_heatmap.py \\
    --correlation-csv data/pearson_correlation.csv \\
    --output          figures/fig_S1_pearson_heatmap.tiff

# Compute directly from the TPM matrix:
python fig_pearson_heatmap.py \\
    --tpm       data/tpm_expression.csv \\
    --metadata  data/sample_info.csv \\
    --output    figures/fig_S1_pearson_heatmap.tiff

Input format for --correlation-csv (CSV with header and row index):
    ,9a-PIM6-1d,9a-PIM6-3d,...
    9a-PIM6-1d,1.0000,0.6911,...
    ...

Strain labels are auto-detected from row/column names using --strain-prefixes
(default: "9a" → 9a5c, "Tem1" → Temecula1) for colour-coding the annotation
strip.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ---------------------------------------------------------------------------
# Strain detection for annotation strip
# ---------------------------------------------------------------------------

DEFAULT_STRAIN_PREFIXES = {
    "9a":   "#3399FF",   # blue
    "Tem1": "#FF6666",   # red
}


def detect_strain_color(label: str, prefix_colors: dict) -> str:
    for prefix, color in prefix_colors.items():
        if label.startswith(prefix):
            return color
    return "#aaaaaa"


def make_row_colors(labels: list, prefix_colors: dict) -> pd.Series:
    return pd.Series(
        [detect_strain_color(lbl, prefix_colors) for lbl in labels],
        index=labels,
        name="strain",
    )


# ---------------------------------------------------------------------------
# Correlation matrix construction
# ---------------------------------------------------------------------------

def load_from_csv(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, index_col=0, sep=sep)


def compute_from_tpm(tpm_path: Path, meta_path: Path, condition_col: str) -> pd.DataFrame:
    """Average replicates → log2(TPM+1) → Pearson r between conditions."""
    sep_t = "\t" if tpm_path.suffix in {".tsv", ".txt"} else ","
    sep_m = "\t" if meta_path.suffix in {".tsv", ".txt"} else ","
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sep_t)
    meta = pd.read_csv(meta_path, index_col=0, sep=sep_m)
    shared = tpm.columns.intersection(meta.index)
    tpm, meta = tpm[shared], meta.loc[shared]

    expr_log = np.log2(tpm + 1)
    # Condition-level means
    cond_means = {}
    for cond, rows in meta.groupby(condition_col):
        sids = rows.index.intersection(expr_log.columns)
        cond_means[cond] = expr_log[sids].mean(axis=1)
    df_cond = pd.DataFrame(cond_means)   # genes × conditions

    corr = df_cond.corr(method="pearson")
    print(f"[INFO] Computed Pearson r for {len(corr)} conditions")
    return corr


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_clustermap(corr: pd.DataFrame, prefix_colors: dict, output: Path):
    labels = list(corr.index)
    row_colors = make_row_colors(labels, prefix_colors)

    vmin = corr.values[~np.eye(len(corr), dtype=bool)].min()
    vmax = 1.0

    g = sns.clustermap(
        corr,
        cmap="Oranges",
        vmin=vmin,
        vmax=vmax,
        figsize=(9, 8),
        annot=False,
        linewidths=0.0,
        row_colors=row_colors,
        col_colors=row_colors,
        dendrogram_ratio=(0.15, 0.15),
        metric="euclidean",
        method="average",
        cbar_kws={"shrink": 0.8},
    )

    # Fix x-axis labels (reordered by dendrogram)
    g.ax_heatmap.set_xticklabels(
        [labels[i] for i in g.dendrogram_col.reordered_ind],
        rotation=25, ha="right", fontsize=9,
    )
    g.ax_heatmap.set_yticklabels(
        [labels[i] for i in g.dendrogram_row.reordered_ind],
        rotation=0, fontsize=9,
    )
    g.ax_heatmap.tick_params(axis="both", length=0)

    # Reposition colorbar
    g.fig.canvas.draw()
    cbar_ax = g.ax_cbar
    cbar_ax.set_position([1.02, 0.65, 0.03, 0.20])
    cbar_ax.set_title("Pearson's r", fontsize=9, pad=8)

    # Strain legend
    from matplotlib.patches import Patch
    unique_colors = {detect_strain_color(l, prefix_colors) for l in labels}
    strain_labels = {v: k for k, v in prefix_colors.items() if v in unique_colors}
    handles = [Patch(facecolor=c, label=strain_labels.get(c, c))
               for c in sorted(unique_colors)]
    g.ax_col_dendrogram.legend(
        handles=handles, title="Strain",
        fontsize=8, title_fontsize=8,
        loc="upper left", frameon=True,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    g.fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"[INFO] Pearson heatmap saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Hierarchically clustered Pearson correlation heatmap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--correlation-csv", type=Path,
                     help="Pre-computed condition × condition Pearson r matrix.")
    grp.add_argument("--tpm", type=Path,
                     help="TPM matrix (genes × samples); Pearson r computed on-the-fly.")
    p.add_argument("--metadata", type=Path, default=None,
                   help="Required when using --tpm.")
    p.add_argument("--condition-col", default="condition",
                   help="Metadata column with condition labels.")
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if args.correlation_csv:
        corr = load_from_csv(args.correlation_csv)
        print(f"[INFO] Loaded pre-computed correlation matrix ({corr.shape[0]} × {corr.shape[1]})")
    else:
        if args.metadata is None:
            sys.exit("[ERROR] --metadata is required when using --tpm.")
        corr = compute_from_tpm(args.tpm, args.metadata, args.condition_col)

    plot_clustermap(corr, DEFAULT_STRAIN_PREFIXES, args.output)


if __name__ == "__main__":
    main()
