"""
fig_virulence_trends.py — TPM expression trend lines for virulence-associated
genes across growth phases (Figure 4D).

For each gene the script plots mean log2(TPM+1) across conditions ordered
by a user-defined progression (e.g. early → mobile → sessile timepoints).
Genes are grouped by phase and drawn as overlapping semi-transparent lines,
with a bold mean trend per group.

Usage
-----
python fig_virulence_trends.py \\
    --tpm             data/tpm_expression.csv \\
    --metadata        data/sample_info.csv \\
    --virulence-table data/virulence_table.csv \\
    --condition-order 9a5c_PIM6_1d 9a5c_PIM6_3d 9a5c_XFM2_1d 9a5c_XFM2_3d \\
    --output          figures/fig_4D_virulence_trends.tiff

Arguments
---------
--condition-order   Ordered list of conditions for the x-axis. Must match
                    the condition column in --metadata. If omitted, conditions
                    are sorted alphabetically.
--condition-col     Column in metadata that holds condition labels (default: condition).
--phase-col         Column in virulence table with phase codes (default: phase).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PHASE_COLORS = {
    0: "#888888",  # unspecific
    1: "#4d9ad1",  # mobile
    2: "#e07b54",  # sessile
    3: "#59a86f",  # early
    4: "#c457b5",  # late
}
PHASE_LABELS = {0: "unspecific", 1: "mobile", 2: "sessile",
                3: "early", 4: "late"}


def load_data(tpm_path, meta_path, vir_path):
    sep_t = "\t" if tpm_path.suffix in {".tsv", ".txt"} else ","
    sep_m = "\t" if meta_path.suffix in {".tsv", ".txt"} else ","
    sep_v = "\t" if vir_path.suffix in {".tsv", ".txt"} else ","
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sep_t)
    meta = pd.read_csv(meta_path, index_col=0, sep=sep_m)
    vir  = pd.read_csv(vir_path, sep=sep_v)
    return tpm, meta, vir


def compute_mean_log2(
    tpm: pd.DataFrame,
    meta: pd.DataFrame,
    condition_col: str,
    conditions: list[str],
) -> pd.DataFrame:
    """Return DataFrame (genes × conditions) of mean log2(TPM+1)."""
    expr_log = np.log2(tpm + 1)
    shared   = tpm.columns.intersection(meta.index)
    expr_log = expr_log[shared]
    meta     = meta.loc[shared]

    cols = {}
    for cond in conditions:
        sids = meta[meta[condition_col] == cond].index.intersection(expr_log.columns)
        if len(sids) == 0:
            cols[cond] = np.nan
        else:
            cols[cond] = expr_log[sids].mean(axis=1)
    return pd.DataFrame(cols)   # genes × conditions


def plot_trends(
    mean_df: pd.DataFrame,
    virulence: pd.DataFrame,
    conditions: list[str],
    phase_col: str,
    output: Path,
):
    id_col = next(c for c in virulence.columns if "gene_id" in c.lower())
    phase_map = dict(zip(virulence[id_col], virulence[phase_col]))

    phases = sorted(set(phase_map.values()))
    n_phases = len(phases)

    fig, axes = plt.subplots(
        1, n_phases,
        figsize=(4 * n_phases, 4),
        sharey=True,
    )
    if n_phases == 1:
        axes = [axes]

    x = np.arange(len(conditions))

    for ax, phase in zip(axes, phases):
        genes_in_phase = [
            g for g in virulence[id_col]
            if phase_map.get(g) == phase and g in mean_df.index
        ]
        color = PHASE_COLORS.get(phase, "#aaaaaa")

        for gene in genes_in_phase:
            y = mean_df.loc[gene, conditions].values.astype(float)
            ax.plot(x, y, color=color, alpha=0.25, linewidth=0.8)

        if genes_in_phase:
            mean_trend = mean_df.loc[genes_in_phase, conditions].mean(axis=0).values
            ax.plot(x, mean_trend, color=color, linewidth=2.5,
                    label=f"n={len(genes_in_phase)}")

        ax.set_title(f"{PHASE_LABELS.get(phase, phase)}\n(n={len(genes_in_phase)})",
                     fontsize=10, color=color, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(conditions, rotation=45, ha="right", fontsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    axes[0].set_ylabel("Mean log₂(TPM+1)", fontsize=10)
    fig.suptitle("Virulence gene expression trends", fontsize=11, y=1.02)

    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Trend plot saved to {output}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Virulence gene expression trend lines.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm",             required=True, type=Path)
    p.add_argument("--metadata",        required=True, type=Path)
    p.add_argument("--virulence-table", required=True, type=Path)
    p.add_argument("--condition-order", nargs="*", default=None,
                   help="Ordered condition labels for x-axis.")
    p.add_argument("--condition-col",   default="condition")
    p.add_argument("--phase-col",       default="phase")
    p.add_argument("--output",          required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm, meta, vir = load_data(args.tpm, args.metadata, args.virulence_table)

    if args.condition_order:
        conditions = args.condition_order
    else:
        conditions = sorted(meta[args.condition_col].unique())

    print(f"[INFO] {len(vir)} virulence genes; {len(conditions)} conditions")
    mean_df = compute_mean_log2(tpm, meta, args.condition_col, conditions)
    plot_trends(mean_df, vir, conditions, args.phase_col, args.output)


if __name__ == "__main__":
    main()
