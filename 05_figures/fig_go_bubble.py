"""
fig_go_bubble.py — Combined GO enrichment bubble plot for DESeq2 pairwise
comparisons AND WGCNA co-expression modules (Figure 4C).

Each bubble encodes:
  x-axis  → comparison label or WGCNA module name
  y-axis  → GO term description
  size    → −log10(FDR-adjusted p-value)
  colour  → ontology category (BP / MF / CC)  OR  enrichment direction
             (up-regulated vs. down-regulated / module enriched)

The script reads GO enrichment CSV files produced by run_go_enrichment.py
(for DESeq2 comparisons) and, optionally, the same script run on WGCNA
module member gene lists (see note below).

Usage
-----
python fig_go_bubble.py \\
    --deseq-dir      results/DESeq2_results \\
    --wgcna-dir      results/WGCNA_go \\
    --top-n-terms    30 \\
    --alpha          0.05 \\
    --output         figures/fig_4C_go_bubble.tiff

Arguments
---------
--deseq-dir    Root directory from run_go_enrichment.py (DESeq2 comparisons).
               Each sub-folder named <c1>_vs_<c2> must contain
               GO_enrichment_up_in_<label>.csv files.
--wgcna-dir    Optional. Directory containing per-module GO enrichment files
               named GO_enrichment_module_<name>.csv (same column format as
               DESeq2 GO output).  Produced by running run_go_enrichment.py
               with module member gene lists (see note).
--top-n-terms  Select the N most significant GO terms across all sources
               (default: 30).
--alpha        FDR threshold for inclusion (default: 0.05).

Note on WGCNA module GO enrichment
-----------------------------------
run_go_enrichment.py expects DESeq2 result folders, but the same Fisher's
exact test can be applied to WGCNA module gene sets by creating a synthetic
folder structure:

    results/WGCNA_go/
        module_turquoise/
            GO_enrichment_up_in_turquoise.csv   ← members of turquoise module
        module_blue/
            GO_enrichment_up_in_blue.csv
        ...

Each CSV should be produced by running run_go_enrichment.py where the "DEGs"
are the module member genes (all treated as up-regulated for enrichment
purposes).  A helper script or a simple wrapper around run_go_enrichment.py
with a gene-list-per-module input is the recommended approach.
"""

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Category colours
# ---------------------------------------------------------------------------

ONTOLOGY_COLORS = {
    "BP": "#e95050",  # Biological Process (red)
    "MF": "#4e904a",  # Molecular Function (green)
    "CC": "#3367e1",  # Cellular Component (blue)
    "?":  "#aaaaaa",  # unknown
}

DIRECTION_COLORS = {
    "up":     "#d7191c",
    "down":   "#2c7bb6",
    "module": "#59a86f",
}

# Simple keyword-based ontology guesser (fallback when no category info)
_BP_KEYWORDS = {"process", "pathway", "response", "regulation", "biosynthesis",
                "catabolism", "signaling", "metabolism", "cycle"}
_MF_KEYWORDS = {"binding", "activity", "catalytic", "transporter", "kinase",
                "transferase", "hydrolase", "reductase", "synthase"}
_CC_KEYWORDS = {"membrane", "cytoplasm", "ribosome", "periplasm", "envelope",
                "chromosome", "flagellum", "pilus"}


def guess_ontology(term: str) -> str:
    t = term.lower()
    if any(k in t for k in _CC_KEYWORDS):
        return "CC"
    if any(k in t for k in _MF_KEYWORDS):
        return "MF"
    if any(k in t for k in _BP_KEYWORDS):
        return "BP"
    return "?"


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_enrichment_files(
    root: Path,
    source_type: str,   # "deseq2" or "wgcna"
    alpha: float,
) -> pd.DataFrame:
    """
    Parse GO enrichment CSVs from a directory tree.
    Returns long-format DataFrame: term, p_adj, source_label, direction, ontology.
    """
    rows = []

    for sub in sorted(root.iterdir()):
        if not sub.is_dir():
            continue

        for csv_file in sorted(sub.glob("GO_enrichment_*.csv")):
            try:
                df = pd.read_csv(csv_file)
            except Exception as e:
                print(f"[WARN] Cannot read {csv_file}: {e}")
                continue

            # Identify required columns flexibly
            padj_col = next(
                (c for c in df.columns if "p_adj" in c.lower()), None
            )
            term_col = next(
                (c for c in df.columns if "go_term" in c.lower()
                 or "go_description" in c.lower()
                 or "term" in c.lower()), None
            )
            if padj_col is None or term_col is None:
                continue

            df = df[df[padj_col] <= alpha].copy()
            if df.empty:
                continue

            # Extract direction / source label from filename
            fname = csv_file.stem  # e.g. "GO_enrichment_up_in_9a5c_PIM6_1d"
            if source_type == "deseq2":
                # Direction: up_in or down_in
                if "up_in" in fname:
                    direction = "up"
                    label_part = re.sub(r"GO_enrichment_up_in_", "", fname)
                else:
                    direction = "down"
                    label_part = re.sub(r"GO_enrichment_(down|up)_in_", "", fname)
                source_label = sub.name   # <c1>_vs_<c2>
            else:
                direction = "module"
                label_part = re.sub(r"GO_enrichment_module_", "", fname)
                source_label = label_part   # module name

            for _, row in df.iterrows():
                term = str(row[term_col])
                rows.append(
                    {
                        "term":         term,
                        "p_adj":        float(row[padj_col]),
                        "source_label": source_label,
                        "direction":    direction,
                        "ontology":     guess_ontology(term),
                        "source_type":  source_type,
                    }
                )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Term selection
# ---------------------------------------------------------------------------

def select_top_terms(df: pd.DataFrame, top_n: int) -> list:
    """Pick terms with the lowest minimum p_adj across all sources."""
    best = df.groupby("term")["p_adj"].min().nsmallest(top_n)
    return best.index.tolist()


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_bubble(
    df: pd.DataFrame,
    top_terms: list,
    color_by: str,   # "ontology" or "direction"
    output: Path,
):
    df_plot = df[df["term"].isin(top_terms)].copy()
    df_plot["-log10_padj"] = -np.log10(df_plot["p_adj"].clip(lower=1e-300))

    # Order axes
    term_order = (
        df_plot.groupby("term")["-log10_padj"].max().sort_values(ascending=False).index.tolist()
    )
    src_order = sorted(
        df_plot["source_label"].unique(),
        key=lambda s: (0 if "vs" in s else 1, s),
    )

    term_idx = {t: i for i, t in enumerate(term_order)}
    src_idx  = {s: i for i, s in enumerate(src_order)}

    color_map = ONTOLOGY_COLORS if color_by == "ontology" else DIRECTION_COLORS
    color_key = "ontology" if color_by == "ontology" else "direction"

    fig_w = max(10, 0.6 * len(src_order) + 4)
    fig_h = max(8,  0.3 * len(term_order) + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    for _, row in df_plot.iterrows():
        x = src_idx[row["source_label"]]
        y = term_idx[row["term"]]
        s = row["-log10_padj"] * 30   # scale factor
        c = color_map.get(row[color_key], "#aaaaaa")
        ax.scatter(x, y, s=s, color=c, alpha=0.8,
                   edgecolors="white", linewidths=0.4, zorder=3)

    ax.set_xticks(range(len(src_order)))
    ax.set_xticklabels(src_order, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(term_order)))
    ax.set_yticklabels(term_order, fontsize=8)
    ax.set_xlim(-0.6, len(src_order) - 0.4)
    ax.set_ylim(-0.6, len(term_order) - 0.4)
    ax.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
    ax.yaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Size legend
    ref_vals = [1, 2, 5]
    for v in ref_vals:
        ax.scatter([], [], s=v * 30, c="grey", alpha=0.7,
                   label=f"−log₁₀(FDR) = {v}")
    leg1 = ax.legend(title="Significance", fontsize=7, title_fontsize=7,
                     loc="lower right", frameon=True)

    # Colour legend
    handles = [mpatches.Patch(facecolor=v, label=k)
               for k, v in color_map.items()
               if k in df_plot[color_key].values]
    ax.legend(handles=handles, title=color_by.capitalize(),
              fontsize=7, title_fontsize=7,
              loc="upper right", frameon=True)
    ax.add_artist(leg1)

    ax.set_xlabel("Comparison / Module", fontsize=10)
    ax.set_ylabel("GO term", fontsize=10)
    ax.set_title("GO enrichment — DESeq2 comparisons & WGCNA modules", fontsize=11)

    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] GO bubble plot saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Combined GO enrichment bubble plot (DESeq2 + WGCNA).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--deseq-dir",    required=True, type=Path,
                   help="Root directory of run_go_enrichment.py output.")
    p.add_argument("--wgcna-dir",    type=Path, default=None,
                   help="Optional: directory with module GO enrichment CSVs.")
    p.add_argument("--top-n-terms",  type=int, default=30,
                   help="Number of most-significant GO terms to display.")
    p.add_argument("--alpha",        type=float, default=0.05)
    p.add_argument("--color-by",     choices=["ontology", "direction"],
                   default="ontology",
                   help="Bubble colour: 'ontology' (BP/MF/CC) or 'direction' (up/down/module).")
    p.add_argument("--output",       required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    frames = []

    if args.deseq_dir.exists():
        df_deseq = load_enrichment_files(args.deseq_dir, "deseq2", args.alpha)
        print(f"[INFO] DESeq2: {len(df_deseq)} significant GO associations")
        frames.append(df_deseq)
    else:
        print(f"[WARN] --deseq-dir not found: {args.deseq_dir}")

    if args.wgcna_dir and args.wgcna_dir.exists():
        df_wgcna = load_enrichment_files(args.wgcna_dir, "wgcna", args.alpha)
        print(f"[INFO] WGCNA: {len(df_wgcna)} significant GO associations")
        frames.append(df_wgcna)
    elif args.wgcna_dir:
        print(f"[WARN] --wgcna-dir not found: {args.wgcna_dir}")

    if not frames:
        sys.exit("[ERROR] No enrichment data could be loaded.")

    df_all = pd.concat(frames, ignore_index=True)
    top_terms = select_top_terms(df_all, args.top_n_terms)
    print(f"[INFO] Selected {len(top_terms)} GO terms for plotting")

    plot_bubble(df_all, top_terms, args.color_by, args.output)


if __name__ == "__main__":
    main()
