"""
fig_virulence_clustermap.py — Log2FC heatmap for virulence-associated genes
(Figure 4A).

Columns (comparisons) are ordered by phase: unspecific (0) → mobile (1) →
sessile (2), with vertical separator lines between phases. Column clustering
is disabled (col_cluster=False) because the phase-based grouping provides the
biological axis of interest.

Rows (comparisons) are hierarchically clustered.

Usage
-----
python fig_virulence_clustermap.py \\
    --virulence-table data/virulence_table.csv \\
    --deseq-dir       results/DESeq2_results \\
    --output          figures/fig_4A_virulence_clustermap.tiff

Input: virulence_table.csv
    Columns: gene_id, gene_name, function, phase, source (and any others)
    phase: 0=unspecific, 1=mobile, 2=sessile, 3=early, 4=late
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------------------------
# Palette / style constants
# ---------------------------------------------------------------------------

PHASE_COLORS = {
    0: "#888888",   # unspecific (grey)
    1: "#4d9ad1",   # mobile (blue)
    2: "#e07b54",   # sessile (orange)
    3: "#59a86f",   # early (green)
    4: "#c457b5",   # late (purple)
}
PHASE_LABELS = {0: "unspecific", 1: "mobile", 2: "sessile",
                3: "early", 4: "late"}

FUNCTION_COLORS = {
    "adhesion":         "#e6194B",
    "regulation":       "#3cb44b",
    "motility":         "#4363d8",
    "stress response":  "#f58231",
    "secretion":        "#911eb4",
    "metabolism":       "#42d4f4",
    "membrane":         "#f032e6",
    "others":           "#aaaaaa",
}

CMAP = LinearSegmentedColormap.from_list(
    "diverging_lfc", ["#2c7bb6", "#f7f7f7", "#d7191c"]
)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_virulence_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def build_lfc_matrix(
    virulence: pd.DataFrame,
    deseq_dir: Path,
) -> pd.DataFrame:
    """Collect log2FC values for each virulence gene across all comparisons."""
    id_col   = next(c for c in virulence.columns if "gene_id" in c.lower())
    gene_ids = virulence[id_col].tolist()

    lfc_dict: dict[str, dict] = {}

    for comp_dir in sorted(deseq_dir.iterdir()):
        if not comp_dir.is_dir() or "_vs_" not in comp_dir.name:
            continue
        # Try annotated DEG file first, then full results
        cands = sorted(comp_dir.glob("DESeq2_results_full.csv"))
        if not cands:
            continue
        df = pd.read_csv(cands[0])
        gcol = next(
            (c for c in df.columns if c.lower() in {"gene_id", "geneid"}),
            df.columns[0],
        )
        df = df.set_index(gcol)
        col_vals: dict[str, float] = {}
        for gid in gene_ids:
            if gid in df.index and "log2FoldChange" in df.columns:
                col_vals[gid] = df.loc[gid, "log2FoldChange"]
            else:
                col_vals[gid] = np.nan
        lfc_dict[comp_dir.name] = col_vals

    if not lfc_dict:
        raise RuntimeError(
            f"[ERROR] No DESeq2 result directories found under {deseq_dir}."
        )

    lfc_df = pd.DataFrame(lfc_dict).T   # comparisons × genes
    lfc_df.columns = gene_ids
    return lfc_df


def order_columns_by_phase(
    lfc_df: pd.DataFrame,
    virulence: pd.DataFrame,
) -> tuple:
    """Return (reordered lfc_df, row_colors_series, phase_boundaries)."""
    id_col   = next(c for c in virulence.columns if "gene_id" in c.lower())
    phase_map = dict(zip(virulence[id_col], virulence.get("phase", 0)))
    func_map  = dict(zip(virulence[id_col], virulence.get("function", "others")))

    # Sort genes by phase, then by function within phase
    gene_order = (
        pd.DataFrame({
            "gene_id": lfc_df.columns,
            "phase":    [phase_map.get(g, 0) for g in lfc_df.columns],
            "function": [func_map.get(g, "others") for g in lfc_df.columns],
        })
        .sort_values(["phase", "function"])
        ["gene_id"]
        .tolist()
    )
    lfc_ordered = lfc_df[gene_order]

    # Phase separator positions (0-based column index where phase changes)
    phases_ordered = [phase_map.get(g, 0) for g in gene_order]
    boundaries: list[int] = []
    for i in range(1, len(phases_ordered)):
        if phases_ordered[i] != phases_ordered[i - 1]:
            boundaries.append(i)

    # Row color strips: phase
    phase_colors_series = pd.Series(
        [PHASE_COLORS.get(phase_map.get(g, 0), "#cccccc") for g in gene_order],
        index=gene_order,
        name="phase",
    )
    func_colors_series = pd.Series(
        [FUNCTION_COLORS.get(func_map.get(g, "others"), "#aaaaaa") for g in gene_order],
        index=gene_order,
        name="function",
    )
    col_colors = pd.DataFrame({"phase": phase_colors_series,
                               "function": func_colors_series})
    return lfc_ordered, col_colors, boundaries


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_clustermap(
    lfc_df: pd.DataFrame,
    col_colors: pd.DataFrame,
    boundaries: list,
    output: Path,
):
    ncols = lfc_df.shape[1]
    nrows = lfc_df.shape[0]
    vmax  = np.nanpercentile(lfc_df.abs().values, 95)
    vmin  = -vmax

    fig_w = max(22, 0.42 * ncols + 10)
    fig_h = max(8,  0.35 * nrows + 4)

    g = sns.clustermap(
        lfc_df,
        cmap=CMAP,
        vmin=vmin,
        vmax=vmax,
        center=0,
        col_colors=col_colors,
        row_cluster=True,
        col_cluster=False,
        yticklabels=True,
        xticklabels=True,
        figsize=(fig_w, fig_h),
        dendrogram_ratio=(0.10, 0.02),
        cbar_pos=None,
        linewidths=0,
        na_col="#eeeeee",
    )

    # Phase separator lines
    ax = g.ax_heatmap
    for b in boundaries:
        ax.axvline(b, color="white", linewidth=2)

    ax.set_xticklabels(
        [t.get_text() for t in ax.get_xticklabels()],
        rotation=45, ha="right", fontsize=6,
    )
    ax.set_yticklabels(
        [t.get_text() for t in ax.get_yticklabels()],
        rotation=0, fontsize=7,
    )
    ax.set_xlabel("Virulence-associated gene", fontsize=9)
    ax.set_ylabel("Comparison", fontsize=9)

    # Colorbar
    g.fig.canvas.draw()
    hm_pos = ax.get_position()
    cbar_ax = g.fig.add_axes([hm_pos.x1 + 0.01, hm_pos.y0 + 0.3,
                               0.012, hm_pos.height * 0.4])
    sm = plt.cm.ScalarMappable(cmap=CMAP,
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cb = g.fig.colorbar(sm, cax=cbar_ax)
    cb.set_label("log₂FC", fontsize=8)

    # Phase legend
    phase_handles = [
        mpatches.Patch(color=v, label=PHASE_LABELS[k])
        for k, v in PHASE_COLORS.items()
        if k in {0, 1, 2}
    ]
    g.ax_col_dendrogram.legend(
        handles=phase_handles, title="Phase",
        fontsize=7, title_fontsize=7,
        loc="upper left", bbox_to_anchor=(0, 1),
        frameon=True, ncol=len(phase_handles),
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    g.fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(g.fig)
    print(f"[INFO] Virulence clustermap saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Log2FC heatmap for virulence-associated genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--virulence-table", required=True, type=Path,
                   help="CSV table of virulence gene annotations "
                        "(gene_id, gene_name, function, phase).")
    p.add_argument("--deseq-dir",       required=True, type=Path,
                   help="Root directory of DESeq2 results (run_deseq2.py output).")
    p.add_argument("--output",          required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    virulence = load_virulence_table(args.virulence_table)
    print(f"[INFO] {len(virulence)} virulence genes loaded")

    lfc_df = build_lfc_matrix(virulence, args.deseq_dir)
    print(f"[INFO] LFC matrix: {lfc_df.shape[0]} comparisons × {lfc_df.shape[1]} genes")

    lfc_ordered, col_colors, boundaries = order_columns_by_phase(lfc_df, virulence)
    plot_clustermap(lfc_ordered, col_colors, boundaries, args.output)


if __name__ == "__main__":
    main()
