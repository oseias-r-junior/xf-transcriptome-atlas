"""
fig_tpm_bubble.py — TPM expression distribution bubble plot (Figure 2B).

For each condition, the script plots a bubble whose:
  - x-axis  = condition label (grouped by strain/medium/time)
  - y-axis  = mean log2(TPM+1) across replicates
  - size    = percentage of expressed genes (TPM ≥ 1 cutoff)
  - colour  = strain identity

Usage
-----
python fig_tpm_bubble.py \\
    --tpm       data/tpm_expression.csv \\
    --metadata  data/sample_info.csv \\
    --tpm-min   1.0 \\
    --group-col condition \\
    --colour-col strain \\
    --output    figures/fig_tpm_bubble.tiff
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


PALETTE = [
    "#e6194B", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#469990", "#9A6324",
]


def load_data(tpm_path: Path, meta_path: Path) -> tuple:
    sep_t = "\t" if tpm_path.suffix in {".tsv", ".txt"} else ","
    sep_m = "\t" if meta_path.suffix in {".tsv", ".txt"} else ","
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sep_t)
    meta = pd.read_csv(meta_path, index_col=0, sep=sep_m)
    shared = tpm.columns.intersection(meta.index)
    return tpm[shared], meta.loc[shared]


def compute_stats(
    tpm: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    tpm_min: float,
) -> pd.DataFrame:
    """Per-condition mean log2(TPM+1) and % expressed genes."""
    expr_log = np.log2(tpm + 1)
    rows = []
    for cond, samples in meta.groupby(group_col):
        sids = samples.index.intersection(tpm.columns)
        sub  = tpm[sids]
        sub_log = expr_log[sids]
        mean_log2 = sub_log.mean().mean()
        pct_expr  = ((sub >= tpm_min).sum(axis=0).mean() / len(sub)) * 100
        rows.append({group_col: cond,
                     "mean_log2_tpm": round(mean_log2, 3),
                     "pct_expressed":  round(pct_expr, 2)})
    return pd.DataFrame(rows)


def plot(
    stats: pd.DataFrame,
    meta: pd.DataFrame,
    group_col: str,
    colour_col: str | None,
    output: Path,
):
    if colour_col and colour_col in meta.columns:
        # Map condition → colour group via metadata
        cond_to_color_group = (
            meta.groupby(group_col)[colour_col]
            .agg(lambda x: x.mode().iloc[0])
        )
        unique_groups = sorted(cond_to_color_group.unique())
        color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(unique_groups)}
        colors = [color_map[cond_to_color_group.get(c, unique_groups[0])]
                  for c in stats[group_col]]
    else:
        colors = PALETTE[0]
        unique_groups = None

    fig, ax = plt.subplots(figsize=(max(8, 0.7 * len(stats)), 5))

    x = range(len(stats))
    sizes = stats["pct_expressed"] * 4   # scale for visibility

    ax.scatter(
        x, stats["mean_log2_tpm"],
        s=sizes,
        c=colors,
        alpha=0.85,
        edgecolors="white",
        linewidths=0.5,
        zorder=3,
    )

    ax.set_xticks(list(x))
    ax.set_xticklabels(stats[group_col], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Mean log₂(TPM+1)", fontsize=10)
    ax.set_xlabel("Condition", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.yaxis.grid(True, linestyle="--", alpha=0.5, zorder=0)

    # Size legend
    for pct in [25, 50, 75, 100]:
        ax.scatter([], [], s=pct * 4, c="grey", alpha=0.6,
                   label=f"{pct}%")
    leg1 = ax.legend(title="% expressed genes\n(TPM ≥ 1)",
                     fontsize=8, title_fontsize=8,
                     loc="upper right", frameon=True)

    # Colour legend
    if unique_groups:
        handles = [mpatches.Patch(color=color_map[g], label=g)
                   for g in unique_groups]
        ax.legend(handles=handles, title=colour_col, fontsize=8, title_fontsize=8,
                  loc="upper left", frameon=True)
        ax.add_artist(leg1)

    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Bubble plot saved to {output}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="TPM distribution bubble plot.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm",        required=True, type=Path)
    p.add_argument("--metadata",   required=True, type=Path)
    p.add_argument("--tpm-min",    type=float, default=1.0,
                   help="TPM cutoff for 'expressed' genes.")
    p.add_argument("--group-col",  default="condition",
                   help="Metadata column defining conditions (x-axis).")
    p.add_argument("--colour-col", default=None,
                   help="Metadata column for bubble colour grouping (e.g. strain).")
    p.add_argument("--output",     required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm, meta = load_data(args.tpm, args.metadata)
    print(f"[INFO] {tpm.shape[0]} genes × {tpm.shape[1]} samples, "
          f"{meta[args.group_col].nunique()} conditions")

    stats = compute_stats(tpm, meta, args.group_col, args.tpm_min)
    plot(stats, meta, args.group_col, args.colour_col, args.output)


if __name__ == "__main__":
    main()
