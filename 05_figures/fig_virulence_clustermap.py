"""
fig_virulence_clustermap.py — Log2FC heatmap for virulence-associated genes
(Figure 4A).

Design:
  COLUMNS  = 96 virulence genes, ordered by phase (unspecific → mobile →
             sessile) then by functional group. No column clustering.
  ROWS     = all pairwise comparisons (typically 12), ordered by comparison
             group (A: within-strain/temporal → B: within-strain/medium →
             C: cross-strain). No row clustering.
  COLOUR   = log₂FC where significant (padj ≤ α and |LFC| ≥ cutoff);
             non-significant cells appear in grey (NaN → na_col).
  SEPARATORS = dashed lines between phase groups (columns) and comparison
             groups (rows).

Usage
-----
python fig_virulence_clustermap.py \\
    --virulence-table  data/virulence_table.csv \\
    --deseq-dir        results/DESeq2_results \\
    --comparison-groups data/comparison_groups.tsv \\
    --output           figures/fig_4A_virulence_clustermap.tiff

--comparison-groups (TSV, no header):
    comparison_name<TAB>group_label
    9a5c_PIM6_1d_vs_9a5c_PIM6_3d<TAB>A
    9a5c_PWG_3d_vs_9a5c_PWG_10d<TAB>A
    ...
If omitted, comparisons are auto-detected from --deseq-dir sub-folders
and assigned to groups A / B / C based on strain patterns.

--alpha / --lfc     Significance thresholds (default: 0.05 / 1.0).
                    Only significant LFC values are shown; others → grey.
"""

import argparse
import re
from pathlib import Path

import matplotlib.cm as mcm
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap


# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------

PHASE_COLORS = {
    0: "#aaaaaa",   # unspecific
    1: "#4d9ad1",   # mobile
    2: "#e07b54",   # sessile
    3: "#59a86f",   # early
    4: "#c457b5",   # late
}
PHASE_LABELS = {0: "unspecific", 1: "mobile", 2: "sessile",
                3: "early", 4: "late"}

GROUP_COLORS = {
    "A": "#4292c6",  # within-strain temporal  (blue)
    "B": "#e6550d",  # within-strain medium     (orange)
    "C": "#31a354",  # cross-strain             (green)
}
GROUP_LABELS = {
    "A": "A – within-strain / temporal",
    "B": "B – within-strain / medium",
    "C": "C – cross-strain",
}

CMAP = LinearSegmentedColormap.from_list(
    "div_lfc", ["#2c7bb6", "#f7f7f7", "#d7191c"]
)
NA_COLOR = "#d0d0d0"


# ---------------------------------------------------------------------------
# Virulence table
# ---------------------------------------------------------------------------

def load_virulence_table(path: Path) -> pd.DataFrame:
    try:
        df = pd.read_csv(path, sep=None, engine="python")
    except UnicodeDecodeError:
        df = pd.read_csv(path, sep=None, engine="python", encoding="latin-1")

    col_map = {}
    for c in df.columns:
        lc = c.lower().strip()
        if "gene_id" in lc or ("gene" in lc and ("id" in lc or "pd" in lc)):
            col_map[c] = "gene_id"
        elif "function" in lc:
            col_map[c] = "function"
        elif "phase" in lc:
            col_map[c] = "phase"
    df = df.rename(columns=col_map)
    for req in ("gene_id", "function", "phase"):
        if req not in df.columns:
            df[req] = "Unknown" if req != "phase" else 0
    df["gene_id"]  = df["gene_id"].astype(str).str.strip().str.upper()
    df["phase"]    = pd.to_numeric(df["phase"], errors="coerce").fillna(0).astype(int)
    df["function"] = df["function"].astype(str).str.strip()
    return df


# ---------------------------------------------------------------------------
# Comparison groups
# ---------------------------------------------------------------------------

def auto_assign_group(comp_name: str) -> str:
    """Heuristic A/B/C assignment based on strain/condition patterns."""
    parts = comp_name.split("_vs_")
    if len(parts) != 2:
        return "?"
    c1, c2 = parts
    strain1 = c1.split("_")[0]
    strain2 = c2.split("_")[0]
    if strain1 != strain2:
        return "C"
    # Same strain: check if medium changes
    med1 = next((p for p in c1.split("_") if p in ("PIM6", "PWG")), "")
    med2 = next((p for p in c2.split("_") if p in ("PIM6", "PWG")), "")
    return "A" if med1 == med2 else "B"


def load_comparison_groups(
    deseq_dir: Path,
    groups_file: Path | None,
) -> dict[str, list[str]]:
    """Return ordered dict {group: [comp_name, ...]}."""
    if groups_file and groups_file.exists():
        df = pd.read_csv(groups_file, sep="\t", header=None,
                         names=["comp", "group"])
        groups: dict[str, list] = {}
        for _, row in df.iterrows():
            groups.setdefault(row["group"], []).append(row["comp"])
        return groups

    # Auto-detect from directory
    comps = sorted(
        d.name for d in deseq_dir.iterdir()
        if d.is_dir() and "_vs_" in d.name
    )
    groups: dict[str, list] = {"A": [], "B": [], "C": [], "?": []}
    for c in comps:
        groups[auto_assign_group(c)].append(c)
    return {k: v for k, v in groups.items() if v}


# ---------------------------------------------------------------------------
# LFC matrix (significant values only)
# ---------------------------------------------------------------------------

def build_lfc_matrix(
    virulence: pd.DataFrame,
    deseq_dir: Path,
    groups: dict[str, list],
    alpha: float,
    lfc_cutoff: float,
) -> tuple:
    gene_ids  = virulence["gene_id"].tolist()
    all_comps = [c for grp in groups.values() for c in grp]

    rows = {}
    for comp in all_comps:
        full_csv = deseq_dir / comp / "DESeq2_results_full.csv"
        row = dict.fromkeys(gene_ids, np.nan)
        if not full_csv.exists():
            print(f"  [INFO] {comp}: no results → all-grey row")
            rows[comp] = row
            continue
        try:
            res = pd.read_csv(full_csv)
        except Exception as e:
            print(f"  [WARN] {comp}: {e}")
            rows[comp] = row
            continue

        gcol   = next((c for c in res.columns if "gene_id" in c.lower()), res.columns[0])
        lfc_c  = next((c for c in res.columns if "log2foldchange" in c.lower()), None)
        padj_c = next((c for c in res.columns if "padj" in c.lower()), None)
        if not lfc_c:
            rows[comp] = row
            continue

        res = res.set_index(gcol)
        for gid in gene_ids:
            if gid not in res.index:
                continue
            lfc = res.loc[gid, lfc_c]
            if padj_c:
                padj = res.loc[gid, padj_c]
                if pd.notna(padj) and padj <= alpha and abs(lfc) >= lfc_cutoff:
                    row[gid] = float(lfc)
            else:
                row[gid] = float(lfc)
        rows[comp] = row

    matrix = pd.DataFrame(rows, index=gene_ids).T
    matrix.index = all_comps
    n_sig = matrix.notna().sum().sum()
    print(f"[INFO] LFC matrix: {matrix.shape[0]} × {matrix.shape[1]}; "
          f"{n_sig} significant cells ({100*n_sig/matrix.size:.0f}%)")
    return matrix


# ---------------------------------------------------------------------------
# Column ordering and colour strips
# ---------------------------------------------------------------------------

def order_columns(
    matrix: pd.DataFrame,
    virulence: pd.DataFrame,
) -> tuple:
    col_meta = virulence.set_index("gene_id").reindex(matrix.columns)
    col_meta["function"] = col_meta["function"].fillna("Unknown")
    col_meta["phase"]    = col_meta["phase"].fillna(0).astype(int)

    gene_order = (
        col_meta.sort_values(["phase", "function"]).index.tolist()
    )
    gene_order = [g for g in gene_order if g in matrix.columns]
    matrix     = matrix[gene_order]
    col_meta   = col_meta.loc[gene_order]

    # Phase boundaries
    phases     = col_meta["phase"].tolist()
    boundaries = [i for i in range(1, len(phases)) if phases[i] != phases[i - 1]]

    # Function colour strip
    unique_funcs = sorted(col_meta["function"].unique())
    func_pal     = dict(zip(unique_funcs,
                             sns.color_palette("tab20", n_colors=max(2, len(unique_funcs)))))
    func_colors  = [func_pal[col_meta.loc[g, "function"]] for g in gene_order]
    phase_colors = [PHASE_COLORS.get(col_meta.loc[g, "phase"], "#cccccc") for g in gene_order]

    col_colors_df = pd.DataFrame(
        {"Function": func_colors, "Phase": phase_colors},
        index=gene_order,
    )
    return matrix, col_meta, col_colors_df, boundaries, func_pal


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_heatmap(
    matrix: pd.DataFrame,
    col_colors_df: pd.DataFrame,
    groups: dict,
    phase_col_boundaries: list,
    func_pal: dict,
    output: Path,
):
    # Row colour strip
    group_of: dict[str, str] = {
        c: grp for grp, comps in groups.items() for c in comps
    }
    row_colors = pd.Series(
        [GROUP_COLORS.get(group_of.get(c, "?"), "#cccccc") for c in matrix.index],
        index=matrix.index,
        name="Group",
    )

    vmax = np.nanpercentile(matrix.abs().values, 95)
    if vmax == 0 or np.isnan(vmax):
        vmax = 2.0

    ncols = matrix.shape[1]
    nrows = matrix.shape[0]
    fig_w = max(22, 0.38 * ncols + 10)
    fig_h = max(6,  0.50 * nrows + 3)

    g = sns.clustermap(
        matrix,
        cmap=CMAP,
        vmin=-vmax, vmax=vmax, center=0,
        col_colors=col_colors_df,
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=False,
        xticklabels=True, yticklabels=True,
        dendrogram_ratio=(0.01, 0.01),
        cbar_pos=None,
        linewidths=0,
        figsize=(fig_w, fig_h),
        na_col=NA_COLOR,
    )
    ax = g.ax_heatmap

    # Group separator lines (rows)
    n_in_group = [len(v) for v in groups.values()]
    cumulative  = [sum(n_in_group[:i]) for i in range(1, len(n_in_group))]
    for sep_y in cumulative:
        ax.axhline(y=sep_y, color="white",   linewidth=2.5, zorder=5)
        ax.axhline(y=sep_y, color="#444444", linewidth=0.9,
                   linestyle="--", zorder=6)

    # Phase separator lines (columns)
    for sep_x in phase_col_boundaries:
        ax.axvline(x=sep_x, color="white",   linewidth=2.5, zorder=5)
        ax.axvline(x=sep_x, color="#444444", linewidth=0.8,
                   linestyle=":",  zorder=6)

    # Labels
    ax.set_xticklabels(
        [t.get_text() for t in ax.get_xticklabels()],
        rotation=45, ha="right", fontsize=6,
    )
    ax.set_yticklabels(
        [t.get_text() for t in ax.get_yticklabels()],
        rotation=0, fontsize=8.5,
    )
    ax.tick_params(axis="both", length=0)
    ax.set_xlabel("Virulence-associated gene", fontsize=10)
    ax.set_ylabel("Comparison", fontsize=10)

    # Colorbar
    g.fig.canvas.draw()
    hm_pos = ax.get_position()
    cbar_ax = g.fig.add_axes(
        [hm_pos.x1 + 0.015, hm_pos.y0 + 0.35,
         0.012, hm_pos.height * 0.35]
    )
    sm = mcm.ScalarMappable(
        cmap=CMAP,
        norm=mcolors.Normalize(vmin=-vmax, vmax=vmax),
    )
    sm.set_array([])
    cb = g.fig.colorbar(sm, cax=cbar_ax)
    cb.set_label("log₂FC (significant only)", fontsize=8, rotation=270, labelpad=14)
    cb.ax.tick_params(labelsize=7)

    # Legends
    group_handles = [
        mpatches.Patch(facecolor=GROUP_COLORS[grp], label=GROUP_LABELS[grp])
        for grp in groups if grp in GROUP_COLORS
    ]
    phase_handles = [
        mpatches.Patch(facecolor=PHASE_COLORS[k], label=PHASE_LABELS.get(k, f"Phase {k}"))
        for k in sorted(PHASE_COLORS)
        if k in [0, 1, 2]
    ]
    na_handle = mpatches.Patch(facecolor=NA_COLOR, label="Not significant / no data")

    legend = g.fig.legend(
        handles=group_handles + [mpatches.Patch(visible=False)]
                + phase_handles + [mpatches.Patch(visible=False)] + [na_handle],
        loc="lower left",
        bbox_to_anchor=(hm_pos.x1 + 0.01, hm_pos.y0),
        bbox_transform=g.fig.transFigure,
        fontsize=7.5, frameon=True, framealpha=0.9,
    )
    func_handles = [mpatches.Patch(facecolor=c, label=k)
                    for k, c in sorted(func_pal.items())]
    func_legend = g.fig.legend(
        handles=func_handles, title="Function",
        loc="upper left",
        bbox_to_anchor=(hm_pos.x1 + 0.01, hm_pos.y1),
        bbox_transform=g.fig.transFigure,
        fontsize=7, title_fontsize=8, ncol=2,
        frameon=True, framealpha=0.9,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    g.fig.savefig(output, dpi=300, bbox_inches="tight",
                  bbox_extra_artists=(legend, func_legend))
    plt.close(g.fig)
    print(f"[INFO] Heatmap saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Log2FC heatmap for virulence-associated genes (Fig 4A).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--virulence-table",    required=True, type=Path)
    p.add_argument("--deseq-dir",          required=True, type=Path)
    p.add_argument("--comparison-groups",  type=Path, default=None,
                   help="TSV: comparison_name<TAB>group. Auto-detected if omitted.")
    p.add_argument("--alpha",  type=float, default=0.05)
    p.add_argument("--lfc",    type=float, default=1.0)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    virulence = load_virulence_table(args.virulence_table)
    print(f"[INFO] {len(virulence)} virulence genes loaded")

    groups = load_comparison_groups(args.deseq_dir, args.comparison_groups)
    for grp, comps in groups.items():
        print(f"  Group {grp}: {len(comps)} comparisons")

    matrix = build_lfc_matrix(virulence, args.deseq_dir, groups, args.alpha, args.lfc)
    matrix, _, col_colors_df, phase_boundaries, func_pal = \
        order_columns(matrix, virulence)

    plot_heatmap(matrix, col_colors_df, groups, phase_boundaries, func_pal, args.output)


if __name__ == "__main__":
    main()
