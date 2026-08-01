"""
fig_virulence_clustermap.py — Log2FC heatmap for virulence-associated genes
(Figure 4A).

Design:
  COLUMNS  = virulence genes found in at least one significant comparison,
             ordered by phase (0 → 1 → 2) then by within-phase hierarchical
             clustering (average/euclidean). No column clustering.
  ROWS     = all pairwise comparisons (typically 12), ordered by comparison
             group (A: within-strain/temporal → B: within-strain/medium →
             C: cross-strain). No row clustering.
  COLOUR   = log₂FC where significant (padj ≤ α and |LFC| ≥ cutoff);
             non-significant / absent cells appear in grey.
  SEPARATORS = dotted vertical lines between phases, dashed horizontal lines
             between comparison groups.

ID resolution:
  DESeq2 results use IMG locus tags (XF9A_* for 9a5c, XFTem_* for
  Temecula1).  The virulence table uses NCBI old locus tags (PD####).
  --gene-dict provides the RBH BLASTP table that bridges the two systems.

Usage
-----
python fig_virulence_clustermap.py \\
    --virulence-table  data/virulence_table.csv \\
    --gene-dict        data/gene_dictionary.csv \\
    --deseq-dir        results/DESeq2_results \\
    --comparison-groups data/comparison_groups.tsv \\
    --output           figures/fig_4A_virulence_clustermap.tiff

--comparison-groups (TSV, no header):
    comparison_name<TAB>group_label
    9a5c_PIM6_1d_vs_9a5c_PIM6_3d<TAB>A
    ...
If omitted, comparisons are auto-detected from --deseq-dir sub-folders
and assigned to A / B / C based on strain/medium patterns.

--alpha / --lfc     Significance thresholds (default: 0.05 / 1.0).
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
from scipy.cluster.hierarchy import leaves_list, linkage


# ---------------------------------------------------------------------------
# Visual constants
# ---------------------------------------------------------------------------

PHASE_COLORS = {0: "#aaaaaa", 1: "#4d9ad1", 2: "#e07b54"}
PHASE_LABELS = {0: "unspecific", 1: "mobile", 2: "sessile"}

GROUP_COLORS = {
    "A": "#4292c6",   # within-strain temporal
    "B": "#e6550d",   # within-strain medium
    "C": "#31a354",   # cross-strain
}
GROUP_LABELS = {
    "A": "A – within-strain / temporal",
    "B": "B – within-strain / medium",
    "C": "C – cross-strain",
}

CMAP     = LinearSegmentedColormap.from_list("div_lfc", ["#2c7bb6", "#f7f7f7", "#b2182b"])
NA_COLOR = "#d0d0d0"


# ---------------------------------------------------------------------------
# Virulence table
# ---------------------------------------------------------------------------

def load_virulence_table(path: Path) -> tuple[pd.DataFrame, str]:
    """Load virulence table; return (df, id_col)."""
    try:
        df = pd.read_csv(path, sep=None, engine="python", encoding="utf-8")
    except UnicodeDecodeError:
        df = pd.read_csv(path, sep=None, engine="python", encoding="latin-1")

    col_map: dict[str, str] = {}
    for c in df.columns:
        lc = c.lower().strip()
        if "tem" in lc and ("gene" in lc or "id" in lc):
            col_map[c] = "tem_gene_id"
        elif "function" in lc:
            col_map[c] = "Function"
        elif "phase" in lc:
            col_map[c] = "Phase"
        elif "gene" in lc and "id" in lc:
            col_map.setdefault(c, "gene_id")
    df = df.rename(columns=col_map)

    id_col = "tem_gene_id" if "tem_gene_id" in df.columns else "gene_id"
    if id_col not in df.columns:
        df[id_col] = df.iloc[:, 0].astype(str)

    df[id_col]      = df[id_col].astype(str).str.strip().str.upper()
    df["Function"]  = df.get("Function",
                              pd.Series("Unknown", index=df.index)).fillna("Unknown")
    df["Phase"]     = pd.to_numeric(
        df.get("Phase", pd.Series(0, index=df.index)), errors="coerce"
    ).fillna(0).astype(int)
    return df, id_col


# ---------------------------------------------------------------------------
# Gene dictionary (IMG locus tags → PD#### NCBI old locus tags)
# ---------------------------------------------------------------------------

def _extract_pd(s: object) -> str | None:
    """Extract PD#### (no underscore, 4 digits) from an old_locus_tags string."""
    if pd.isna(s):
        return None
    s = str(s)
    m = re.search(r"\bPD(\d{4})\b", s)
    if m:
        return f"PD{m.group(1)}"
    m = re.search(r"\bPD_(\d{4})\b", s)
    if m:
        return f"PD{m.group(1)}"
    return None


def build_img_to_pd(gene_dict_path: Path, valid_pd_ids: set[str]) -> dict[str, str]:
    """
    Parse gene_dictionary.csv and build a mapping:
        IMG locus tag (uppercase) → PD#### NCBI old locus tag
    Only entries whose PD#### is in valid_pd_ids (virulence gene set) are kept.
    """
    try:
        gd = pd.read_csv(gene_dict_path, sep="\t", dtype=str)
    except Exception:
        gd = pd.read_csv(gene_dict_path, sep=None, engine="python", dtype=str)

    # Locate the Temecula1 old_locus_tags column
    tem_old_col = next(
        (c for c in gd.columns if "old" in c.lower() and "tem" in c.lower()),
        next((c for c in gd.columns if "old_locus" in c.lower()), None),
    )
    if tem_old_col is None:
        raise ValueError(
            f"Could not find Temecula1 old_locus_tags column in {gene_dict_path}.\n"
            f"Available columns: {gd.columns.tolist()}"
        )

    mapping: dict[str, str] = {}
    for _, row in gd.iterrows():
        pd_id = _extract_pd(row.get(tem_old_col))
        if pd_id is None or pd_id not in valid_pd_ids:
            continue
        # 9a5c IMG IDs → PD####
        for col in ("9a5c_IMG_ID", "9a5c_locus_tag", "9a5c_old_locus_tags"):
            val = str(row.get(col, "")).strip().upper()
            if val and val != "NAN":
                mapping[val] = pd_id
        # Temecula1 IMG IDs → PD####
        for col in ("Temecula1_IMG_ID", "Temecula1_locus_tag"):
            val = str(row.get(col, "")).strip().upper()
            if val and val != "NAN":
                mapping[val] = pd_id
        # Direct PD#### identity (in case results already use NCBI IDs)
        mapping[pd_id] = pd_id

    print(f"[INFO] Gene dictionary: {len(gd)} RBH pairs → {len(mapping)} IMG→PD entries")
    return mapping


# ---------------------------------------------------------------------------
# Comparison groups
# ---------------------------------------------------------------------------

def auto_assign_group(comp_name: str) -> str:
    parts = comp_name.split("_vs_")
    if len(parts) != 2:
        return "?"
    c1, c2 = parts
    if c1.split("_")[0] != c2.split("_")[0]:
        return "C"
    m1 = next((p for p in c1.split("_") if p in ("PIM6", "PWG")), "")
    m2 = next((p for p in c2.split("_") if p in ("PIM6", "PWG")), "")
    return "A" if m1 == m2 else "B"


def load_comparison_groups(
    deseq_dir: Path,
    groups_file: Path | None,
) -> dict[str, list[str]]:
    if groups_file and groups_file.exists():
        df = pd.read_csv(groups_file, sep="\t", header=None,
                         names=["comp", "group"])
        groups: dict[str, list] = {}
        for _, row in df.iterrows():
            groups.setdefault(row["group"], []).append(row["comp"])
        return groups

    comps = sorted(
        d.name for d in deseq_dir.iterdir()
        if d.is_dir() and "_vs_" in d.name
    )
    groups = {"A": [], "B": [], "C": [], "?": []}
    for c in comps:
        groups[auto_assign_group(c)].append(c)
    return {k: v for k, v in groups.items() if v}


# ---------------------------------------------------------------------------
# LFC matrix  (inside-out approach + gene dictionary mapping)
# ---------------------------------------------------------------------------

def build_lfc_matrix(
    virulence: pd.DataFrame,
    id_col: str,
    img_to_pd: dict[str, str],
    deseq_dir: Path,
    groups: dict[str, list],
    alpha: float,
    lfc_cutoff: float,
) -> pd.DataFrame:
    """
    1. Load significant DEGs from each comparison (inside-out).
    2. Map IMG locus tags → PD#### via gene dictionary.
    3. Drop genes with no significant value in any comparison.
    """
    all_comps    = [c for grp in groups.values() for c in grp]
    vir_ids_list = virulence[id_col].tolist()

    # --- Step 1: collect significant genes per comparison ---
    comp_sig: dict[str, dict] = {}
    for comp in all_comps:
        full_csv = deseq_dir / comp / "DESeq2_results_full.csv"
        if not full_csv.exists():
            print(f"  [INFO] {comp}: file not found → all-grey row")
            comp_sig[comp] = {}
            continue
        try:
            df = pd.read_csv(full_csv, index_col=0)
        except Exception as e:
            print(f"  [WARN] {comp}: {e}")
            comp_sig[comp] = {}
            continue

        df.index = df.index.astype(str).str.strip().str.upper()
        lfc_c  = next((c for c in df.columns if "log2foldchange" in c.lower()), None)
        padj_c = next((c for c in df.columns if "padj" in c.lower()), None)
        if lfc_c is None:
            print(f"  [WARN] {comp}: no log2FoldChange column → {df.columns.tolist()[:6]}")
            comp_sig[comp] = {}
            continue

        if padj_c:
            mask = df[padj_c].notna() & (df[padj_c] <= alpha) & (df[lfc_c].abs() >= lfc_cutoff)
        else:
            mask = df[lfc_c].abs() >= lfc_cutoff

        comp_sig[comp] = df.loc[mask, lfc_c].to_dict()
        print(f"  {comp[:55]:<55}  {len(comp_sig[comp]):>4} sig. genes")

    # --- Step 2: build matrix with IMG→PD mapping ---
    rows: dict[str, dict] = {}
    for comp in all_comps:
        row: dict[str, float] = dict.fromkeys(vir_ids_list, float("nan"))
        for img_id, lfc in comp_sig[comp].items():
            pd_id = img_to_pd.get(img_id)
            if pd_id is not None and pd_id in row:
                row[pd_id] = float(lfc)
        rows[comp] = row

    matrix = pd.DataFrame(rows, index=vir_ids_list).T
    matrix.index = pd.Index(all_comps, name="comparison")

    # --- Step 3: drop genes absent from all comparisons ---
    before = matrix.shape[1]
    matrix = matrix.dropna(axis=1, how="all")
    after  = matrix.shape[1]
    n_sig  = int(matrix.notna().sum().sum())
    print(f"\n[INFO] Matrix: {matrix.shape[0]} × {after}  "
          f"(removed {before - after} all-grey genes; {n_sig} sig. cells)")

    if after == 0:
        all_ids = {k for d in comp_sig.values() for k in d}
        print("[ERROR] No virulence genes matched. Diagnostics:")
        print(f"  DESeq2 sample IDs : {sorted(all_ids)[:5]}")
        print(f"  Dict key sample   : {list(img_to_pd.keys())[:5]}")
        print(f"  Virulence ID sample: {vir_ids_list[:5]}")
        raise RuntimeError("No virulence genes found in DESeq2 results after mapping.")

    return matrix


# ---------------------------------------------------------------------------
# Column ordering (phase order + within-phase hierarchical clustering)
# ---------------------------------------------------------------------------

def order_columns(
    matrix: pd.DataFrame,
    virulence: pd.DataFrame,
    id_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list, dict]:
    col_meta = virulence.set_index(id_col).reindex(matrix.columns)
    col_meta["Function"] = col_meta["Function"].fillna("Unknown")
    col_meta["Phase"]    = col_meta["Phase"].fillna(0).astype(int)

    data_clust = matrix.fillna(0)
    final_order: list[str] = []
    boundaries: list[int]  = []
    cursor = 0
    prev_ph = None

    for ph in sorted(col_meta["Phase"].unique()):
        cols_ph = col_meta[col_meta["Phase"] == ph].index.tolist()
        if not cols_ph:
            continue
        if prev_ph is not None:
            boundaries.append(cursor)
        prev_ph = ph
        if len(cols_ph) > 1:
            Z   = linkage(data_clust[cols_ph].T, method="average", metric="euclidean")
            ord_ = [cols_ph[i] for i in leaves_list(Z)]
        else:
            ord_ = cols_ph
        final_order.extend(ord_)
        cursor += len(cols_ph)

    matrix   = matrix[final_order]
    col_meta = col_meta.loc[final_order]

    unique_funcs = sorted(col_meta["Function"].unique())
    func_pal     = dict(zip(unique_funcs,
                             sns.color_palette("tab20", n_colors=max(2, len(unique_funcs)))))
    func_colors  = [func_pal[col_meta.loc[g, "Function"]] for g in final_order]
    phase_colors = [PHASE_COLORS.get(col_meta.loc[g, "Phase"], "#cccccc") for g in final_order]

    col_colors_df = pd.DataFrame(
        {"Function": func_colors, "Phase": phase_colors},
        index=final_order,
    )
    return matrix, col_meta, col_colors_df, boundaries, func_pal


# ---------------------------------------------------------------------------
# Row label helper
# ---------------------------------------------------------------------------

def clean_label(s: str) -> str:
    return s.replace("_vs_", " ¶ ").replace("_", " ").replace(" ¶ ", " vs ")


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
) -> None:
    CMAP.set_bad(NA_COLOR)

    # Row colour strip and pretty row labels
    group_of: dict[str, str] = {c: grp for grp, comps in groups.items() for c in comps}
    row_colors = pd.Series(
        [GROUP_COLORS.get(group_of.get(c, "?"), "#cccccc") for c in matrix.index],
        index=matrix.index, name="Group",
    )
    matrix.index         = [clean_label(c) for c in matrix.index]
    row_colors.index     = matrix.index

    vals = matrix.stack(dropna=True).abs()
    vmax = float(vals.max()) if len(vals) > 0 else 2.0

    n_rows = matrix.shape[0]
    n_cols = matrix.shape[1]
    fig_h  = max(7.0, 0.50 * n_rows + 5.5)
    fig_w  = max(12,  min(40, 0.23 * n_cols + 8))
    leg_h  = 3.8
    bot_fr = leg_h / fig_h

    sns.set(font_scale=0.9)
    g = sns.clustermap(
        matrix,
        cmap=CMAP,
        mask=matrix.isna(),          # masked cells → transparent → facecolor shows
        col_colors=col_colors_df,
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=False,
        xticklabels=True, yticklabels=True,
        dendrogram_ratio=(0.005, 0.005),
        cbar_pos=None,
        linewidths=0,
        figsize=(fig_w, fig_h),
    )
    ax = g.ax_heatmap
    ax.set_facecolor(NA_COLOR)       # grey background for masked (NaN) cells

    # Group separators (horizontal)
    n_in_group = [len(v) for v in groups.values()]
    for sep_y in [sum(n_in_group[:i]) for i in range(1, len(n_in_group))]:
        ax.axhline(y=sep_y, color="white",   linewidth=0.5,  zorder=5)
        ax.axhline(y=sep_y, color="#333333", linewidth=0.35, linestyle="--", zorder=6)

    # Phase separators (vertical)
    for sep_x in phase_col_boundaries:
        ax.axvline(x=sep_x, color="white",   linewidth=0.5,  zorder=5)
        ax.axvline(x=sep_x, color="#555555", linewidth=0.75, linestyle=":", zorder=6)

    # Tick labels (after canvas.draw)
    g.fig.canvas.draw()
    if g.ax_col_colors is not None:
        g.ax_col_colors.set_xticklabels([])
        g.ax_col_colors.tick_params(axis="x", length=0)
    ax.set_xticklabels(
        [t.get_text() for t in ax.get_xticklabels()],
        rotation=35, ha="right", fontsize=8,
    )
    ax.set_yticklabels(
        [t.get_text() for t in ax.get_yticklabels()],
        rotation=0, fontsize=9,
    )
    ax.tick_params(axis="both", length=0)
    if g.ax_row_colors is not None:
        for lbl in g.ax_row_colors.get_xticklabels():
            lbl.set_rotation(35)
            lbl.set_ha("right")
            lbl.set_fontsize(7)

    # Legend space at bottom
    g.fig.subplots_adjust(bottom=bot_fr + 0.03)

    # Horizontal colorbar
    cbar_y  = (bot_fr + 0.03) * 0.62
    cbar_ax = g.fig.add_axes([0.28, cbar_y, 0.45, 0.025])
    sm = mcm.ScalarMappable(cmap=CMAP, norm=mcolors.Normalize(vmin=-vmax, vmax=vmax))
    sm.set_array([])
    cb = g.fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cb.ax.tick_params(labelsize=8, length=0)
    cb.set_label("log₂FC (significant only)", fontsize=9)

    # Legends
    func_handles  = [mpatches.Patch(facecolor=func_pal[k], label=k)
                     for k in sorted(func_pal)]
    phase_handles = [mpatches.Patch(facecolor=PHASE_COLORS[k], label=f"Phase {k}")
                     for k in [0, 1, 2] if k in PHASE_COLORS]
    grp_handles   = [mpatches.Patch(facecolor=GROUP_COLORS[k], label=GROUP_LABELS[k])
                     for k in ["A", "B", "C"] if k in groups]
    na_handle     = mpatches.Patch(facecolor=NA_COLOR, label="Not significant / no data")

    func_leg = g.fig.legend(
        handles=func_handles, title="Function",
        loc="lower left", bbox_to_anchor=(0.07, 0.02),
        bbox_transform=g.fig.transFigure,
        fontsize="small", frameon=True,
    )
    right_leg = g.fig.legend(
        handles=grp_handles + [mpatches.Patch(visible=False)]
                + phase_handles + [mpatches.Patch(visible=False)] + [na_handle],
        title="Group / Phase",
        loc="lower right", bbox_to_anchor=(0.98, 0.02),
        bbox_transform=g.fig.transFigure,
        fontsize="small", frameon=True,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    g.fig.savefig(output, dpi=300, bbox_inches="tight",
                  bbox_extra_artists=(func_leg, right_leg))
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
    p.add_argument("--virulence-table",   required=True, type=Path,
                   help="CSV/TSV with gene_id (NCBI PD####), Function, Phase.")
    p.add_argument("--gene-dict",         required=True, type=Path,
                   help="gene_dictionary.csv (RBH BLASTP table mapping IMG↔NCBI IDs).")
    p.add_argument("--deseq-dir",         required=True, type=Path,
                   help="Directory with one sub-folder per comparison containing "
                        "DESeq2_results_full.csv.")
    p.add_argument("--comparison-groups", type=Path, default=None,
                   help="TSV (no header): comparison_name<TAB>group. "
                        "Auto-detected if omitted.")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--lfc",   type=float, default=1.0)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    virulence, id_col = load_virulence_table(args.virulence_table)
    print(f"[INFO] {len(virulence)} virulence genes loaded (ID column: '{id_col}')")

    valid_pd_ids = set(virulence[id_col])
    img_to_pd    = build_img_to_pd(args.gene_dict, valid_pd_ids)

    groups = load_comparison_groups(args.deseq_dir, args.comparison_groups)
    for grp, comps in groups.items():
        print(f"  Group {grp}: {len(comps)} comparisons")

    matrix = build_lfc_matrix(
        virulence, id_col, img_to_pd,
        args.deseq_dir, groups, args.alpha, args.lfc,
    )
    matrix, _, col_colors_df, phase_boundaries, func_pal = \
        order_columns(matrix, virulence, id_col)

    plot_heatmap(matrix, col_colors_df, groups, phase_boundaries, func_pal, args.output)


if __name__ == "__main__":
    main()
