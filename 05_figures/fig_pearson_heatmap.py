"""
fig_pearson_heatmap.py — Hierarchically clustered Pearson correlation heatmap
of condition-level expression profiles (Figure S3).

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
    --output          figures/fig_S3_pearson_heatmap.tiff

# Compute directly from the TPM matrix:
python fig_pearson_heatmap.py \\
    --tpm       data/tpm_expression.csv \\
    --metadata  data/sample_info.csv \\
    --output    figures/fig_S3_pearson_heatmap.tiff

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

# Ordered (prefix, canonical label, color) triples -- checked in order, so
# longer/more specific prefixes must come first. Two naming conventions show
# up for condition/sample labels across this repo: the short one used by
# data/sample_info.csv's own "condition" column ("9a_...", "Tem1_...") and
# the strain column's own full values used when fig_pearson_heatmap.py has
# to derive "condition" itself from strain/medium/timepoint (see
# compute_from_tpm), which are "9a5c"/"Temecula1" written out in full.
# "Temecula1" does NOT start with "Tem1" (5th char is "e", not "1"), so a
# naive {"9a": ..., "Tem1": ...} prefix dict silently fails to match
# "Temecula1_..." labels and falls through to the grey "no match" colour --
# that was the "#aaaaaa" and missing-legend-entry bug. Longer prefixes first
# so "9a5c"/"Temecula1" are tried before the shorter "9a"/"Tem1".
DEFAULT_STRAIN_PREFIXES = [
    ("9a5c",      "9a5c",      "#FF6666"),   # red
    ("9a",        "9a5c",      "#FF6666"),
    ("Temecula1", "Temecula1", "#3399FF"),   # blue
    ("Tem1",      "Temecula1", "#3399FF"),
]


def detect_strain(label: str, prefixes: list) -> tuple[str, str]:
    """Return (canonical_label, color) for the first matching prefix."""
    for prefix, canon, color in prefixes:
        if label.startswith(prefix):
            return canon, color
    return label, "#aaaaaa"


def detect_strain_color(label: str, prefixes: list) -> str:
    return detect_strain(label, prefixes)[1]


def make_row_colors(labels: list, prefixes: list) -> pd.Series:
    return pd.Series(
        [detect_strain_color(lbl, prefixes) for lbl in labels],
        index=labels,
        name="strain",
    )


# ---------------------------------------------------------------------------
# Correlation matrix construction
# ---------------------------------------------------------------------------

def sniff_sep(path: Path) -> str:
    """Detect the field separator from the file's own header line.

    File extension alone is unreliable: compute_tpm.py writes tab-separated
    data into a file named tpm_expression.csv. See the same helper (and the
    bug it fixes) in 04_wgcna/run_wgcna.py.
    """
    with open(path, "r", newline="") as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ("\t", ",", ";")}
    best_sep, best_count = max(counts.items(), key=lambda kv: kv[1])
    return best_sep if best_count > 0 else ","


def load_from_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0, sep=sniff_sep(path))


def compute_from_tpm(tpm_path: Path, meta_path: Path, condition_col: str) -> pd.DataFrame:
    """Average replicates → log2(TPM+1) → Pearson r between conditions."""
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sniff_sep(tpm_path))
    meta = pd.read_csv(meta_path, index_col=0, sep=sniff_sep(meta_path))
    meta.columns = meta.columns.str.strip()
    for c in meta.select_dtypes(include="object").columns:
        meta[c] = meta[c].str.strip()

    if condition_col not in meta.columns:
        # data/sample_info_original.csv (WGCNA_paper_pierry_feitosa.ipynb's
        # own metadata) has no "condition" column, only strain/medium/
        # timepoint(early|late)/replicate. strain+medium+timepoint together
        # still uniquely identify each of the 8 conditions (each
        # strain x medium pair only has two actual timepoints, one "early"
        # and one "late"), so derive "condition" from those instead of
        # failing outright.
        needed = {"strain", "medium", "timepoint"}
        if not needed.issubset(meta.columns):
            sys.exit(
                f"[ERROR] --condition-col {condition_col!r} not found in "
                f"{meta_path.name}, and no strain/medium/timepoint columns "
                f"to derive it from either. Available columns: "
                f"{list(meta.columns)}"
            )
        meta[condition_col] = (
            meta["strain"].astype(str) + "_" +
            meta["medium"].astype(str) + "_" +
            meta["timepoint"].astype(str)
        )
        print(f"[INFO] {condition_col!r} not in {meta_path.name}; derived "
              f"it from strain+medium+timepoint instead "
              f"(e.g. {meta[condition_col].iloc[0]!r}).")

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

def plot_clustermap(corr: pd.DataFrame, prefixes: list, output: Path):
    labels = list(corr.index)
    row_colors = make_row_colors(labels, prefixes)

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

    # row_colors and col_colors are the SAME Series here (both named
    # "strain"), so each gets its own single auto-generated tick label
    # reading "strain" -- but on DIFFERENT axes:
    #   - ax_row_colors (narrow strip to the LEFT of the heatmap) labels
    #     itself with an X-TICK, at the bottom -- i.e. at the same height
    #     as ax_heatmap's own x-tick labels, right next to the leftmost
    #     one, which is what was actually colliding/showing vertical.
    #   - ax_col_colors (narrow strip ABOVE the heatmap) labels itself with
    #     a Y-TICK, on the left -- unrelated to the heatmap's x-axis, and
    #     was not the one causing the overlap; leave it as seaborn drew it.
    if g.ax_row_colors is not None:
        for tick in g.ax_row_colors.get_xticklabels():
            tick.set_rotation(25)
            tick.set_ha("right")
            tick.set_fontsize(9)

    # Reposition colorbar
    g.fig.canvas.draw()
    cbar_ax = g.ax_cbar
    cbar_ax.set_position([1.02, 0.65, 0.03, 0.20])
    cbar_ax.set_title("Pearson's r", fontsize=9, pad=8)

    # Strain legend -- canonical label per colour (detect_strain's 2nd
    # element), not the raw matched prefix, so "Temecula1" shows in full
    # rather than a partial/failed-match artifact.
    from matplotlib.patches import Patch
    color_to_label: dict[str, str] = {}
    for lbl in labels:
        canon, color = detect_strain(lbl, prefixes)
        color_to_label.setdefault(color, canon)
    handles = [Patch(facecolor=c, label=color_to_label[c])
               for c in sorted(color_to_label)]
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
