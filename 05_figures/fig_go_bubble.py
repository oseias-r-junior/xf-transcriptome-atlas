"""
fig_go_bubble.py — Combined GO enrichment bubble plot: DESeq2 pairwise
comparisons (left panel) vs. WGCNA co-expression modules (right panel),
both collapsed to GO-slim "Level-2" terms (Figure 4C).

This reproduces the exact analysis and plot submitted to MDPI Pathogens
(see figure_2_3_4_5_composites_mdpi_pathogens.pptx, slide 4, panel C),
whose data pipeline traces back to the "map_to_2nd_level" /
"build_signed_column" cells of GO_enrichment_pierry_feitosa_paper_v2.ipynb:

  1. Each specific GO term significant in a DESeq2 "up in <condition>"
     enrichment file (p_adj < alpha) is walked up its GO ancestor chain
     (cached in data/ancestor_cache.json, originally fetched from QuickGO)
     until a GO-slim Level-2 term is found (data/dictionary_2_level.csv).
     If both directions of a comparison hit the same Level-2 term, the
     direction with the larger |-log10(p)| wins ("signed" score).
  2. For each WGCNA module's GO enrichment result
     (data/wgcna_enrichment/enrichment_<module>.csv), the Level2_ID column
     (already a real Level-2 collapse computed when that file was
     generated -- verified to differ from GO_ID for most rows, so it is
     NOT a no-op copy) is mapped through go_dict directly to get the same
     term labels the DESeq2 panel uses. No ancestor-chain climbing is
     applied to the WGCNA panel; climbing an already-collapsed ID through
     the same dictionary would silently drop valid terms and was the
     source of a mismatch with the originally submitted figure.
  3. Bubble size = -log10(FDR), capped at MAX_BUBBLE=300 pt^2 and scaled
     relative to the largest -log10(FDR) within each panel separately
     (fixes bubbles overflowing their column, which happened when size was
     scaled linearly with no cap); colour = direction (up in cond. 1 / up
     in cond. 2 / enriched in WGCNA module); a left-hand colour strip marks
     each Level-2 term's GO ontology (Biological Process / Molecular
     Function / Cellular Component).

Usage
-----
python fig_go_bubble.py \\
    --deseq-dir        results/DESeq2_results \\
    --wgcna-dir        data/wgcna_enrichment \\
    --dictionary       data/dictionary_2_level.csv \\
    --ancestor-cache   data/ancestor_cache.json \\
    --alpha            0.05 \\
    --output           figures/fig_4C_go_bubble.tiff

(results/DESeq2_results/<comparison>/GO_enrichment_up_in_<cond>.csv is
where run_go_enrichment.py actually writes its output in this repo --
"data/DESeq2_results_src" does not exist. There is also
data/deseq2_go_results/, an older copy keyed by old-locus-tag gene IDs
with a "_old_locus_tags" filename suffix; both have compatible go_term/
p_adj columns, but results/DESeq2_results is the canonical, currently
regenerated output and is what --deseq-dir should point to.)
"""
import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Colours (sampled from the original submitted figure's legend swatches;
# GO-category colours match GO_enrichment_pierry_feitosa_paper_v2.ipynb's
# `category_colors` dict exactly)
# ---------------------------------------------------------------------------

DIRECTION_COLORS = {
    "up1":    "#F768A1",   # Up in cond. 1 (DESeq2)
    "up2":    "#7B3294",   # Up in cond. 2 (DESeq2)
    "wgcna":  "#5B7FBD",   # Enriched (WGCNA)
}
DIRECTION_LABELS = {
    "up1":   "Up in cond. 1 (DESeq2)",
    "up2":   "Up in cond. 2 (DESeq2)",
    "wgcna": "Enriched (WGCNA)",
}
CATEGORY_COLORS = {"P": "#E95050", "F": "#4E904A", "C": "#3367E1"}
CATEGORY_LABELS = {"P": "Biological Process", "F": "Molecular Function", "C": "Cellular Component"}
CATEGORY_ORDER = ["C", "F", "P"]   # Cellular Component (top) -> Molecular Function -> Biological Process (bottom)


# ---------------------------------------------------------------------------
# Level-2 GO dictionary + ancestor-based collapsing
# ---------------------------------------------------------------------------

def load_level2_dictionary(path: Path):
    df = pd.read_csv(path, header=None,
                      names=["level", "category_name", "GO_id", "term", "relation"])
    df = df[df["relation"] == "is_a"]
    go_to_term = dict(zip(df["GO_id"], df["term"]))
    cat_map = {"biological_process": "P", "molecular_function": "F", "cellular_component": "C"}
    go_to_category = {go: cat_map.get(cat, "NA") for go, cat in zip(df["GO_id"], df["category_name"])}
    return go_to_term, go_to_category


def map_to_2nd_level(go_pvals: dict, go_dict: dict, go_category: dict, ancestors: dict):
    """For each (GO id -> p_adj), climb the ancestor chain (self included) and
    keep the best (min) p-value per Level-2 term reached. Terms with no path
    to a Level-2 ancestor in `ancestors` are skipped (documented limitation:
    ancestor_cache.json only covers GO ids actually queried when this
    analysis was originally run)."""
    valid = set(go_dict.keys())
    term_pvals: dict[str, list] = {}
    term_category: dict[str, str] = {}

    for go, pval in go_pvals.items():
        candidates = [go] + ancestors.get(go, [])
        for anc in candidates:
            if anc in valid:
                term = go_dict[anc]
                term_pvals.setdefault(term, []).append(pval)
                cat = go_category.get(anc, "NA")
                if cat != "NA":
                    term_category[term] = cat

    term_best = {term: min(pvals) for term, pvals in term_pvals.items()}
    return term_best, term_category


def load_go_enrichment_csv(path: Path, alpha: float, id_col: str, padj_col: str):
    df = pd.read_csv(path)
    df = df[df[padj_col] < alpha]
    return dict(zip(df[id_col], df[padj_col]))


# ---------------------------------------------------------------------------
# DESeq2 panel
# ---------------------------------------------------------------------------

def build_deseq2_rows(deseq_dir: Path, go_dict, go_category, ancestors, alpha: float):
    rows = []
    for comp_dir in sorted(p for p in deseq_dir.iterdir() if p.is_dir()):
        files = sorted(comp_dir.glob("GO_enrichment_up_in_*.csv"))
        if len(files) != 2:
            continue
        try:
            cond_a, cond_b = comp_dir.name.split("_vs_")
        except ValueError:
            continue
        file_a = next((f for f in files if cond_a in f.name), None)
        file_b = next((f for f in files if cond_b in f.name), None)
        if file_a is None or file_b is None:
            continue

        pvals_a = load_go_enrichment_csv(file_a, alpha, "go_term", "p_adj")
        pvals_b = load_go_enrichment_csv(file_b, alpha, "go_term", "p_adj")

        terms_a, cat_a = map_to_2nd_level(pvals_a, go_dict, go_category, ancestors)
        terms_b, cat_b = map_to_2nd_level(pvals_b, go_dict, go_category, ancestors)

        term_scores: dict[str, float] = {}
        term_cat: dict[str, str] = {}
        for term, p in terms_a.items():
            term_scores[term] = -np.log10(p)
            term_cat[term] = cat_a.get(term, "NA")
        for term, p in terms_b.items():
            score = -(-np.log10(p))
            if term in term_scores:
                if abs(score) > abs(term_scores[term]):
                    term_scores[term] = score
            else:
                term_scores[term] = score
            term_cat.setdefault(term, cat_b.get(term, "NA"))

        for term, score in term_scores.items():
            rows.append({
                "term": term,
                "source_label": comp_dir.name,
                "panel": "deseq2",
                "direction": "up1" if score > 0 else "up2",
                "neglog10_fdr": abs(score),
                "category": term_cat.get(term, "NA"),
            })
    return rows


# ---------------------------------------------------------------------------
# WGCNA panel
# ---------------------------------------------------------------------------

def build_wgcna_rows(wgcna_dir: Path, go_dict, go_category, alpha: float):
    """Unlike the DESeq2 panel, the WGCNA enrichment files already carry a
    Level2_ID column that is the Level-2 collapse of GO_ID computed when the
    enrichment analysis itself was originally run (verified: Level2_ID !=
    GO_ID for most rows, so it is a real collapse, not a copy). The
    original submitted figure's notebook cell maps WGCNA terms via
    go_dict[row["Level2_ID"]] directly -- no ancestor-chain climbing needed
    or wanted here, since that would re-collapse an already-collapsed ID
    through the same dictionary and silently drop rows whose Level2_ID
    (while a valid Level-2 term) isn't itself a key of `ancestors`."""
    rows = []
    for csv_file in sorted(wgcna_dir.glob("enrichment_*.csv")):
        module = csv_file.stem.replace("enrichment_", "")
        df = pd.read_csv(csv_file)
        if not {"Level2_ID", "P_adj_BH"}.issubset(df.columns):
            print(f"[WARN] {csv_file.name}: missing Level2_ID/P_adj_BH columns, skipping")
            continue
        df_sig = df[df["P_adj_BH"] < alpha].dropna(subset=["Level2_ID"]).copy()
        df_sig["term"] = df_sig["Level2_ID"].map(go_dict)
        df_sig = df_sig.dropna(subset=["term"])
        if df_sig.empty:
            continue
        agg = df_sig.groupby("term")["P_adj_BH"].min()
        lvl2_by_term = df_sig.drop_duplicates("term").set_index("term")["Level2_ID"]
        for term, p in agg.items():
            cat = go_category.get(lvl2_by_term[term], "NA")
            rows.append({
                "term": term,
                "source_label": module,
                "panel": "wgcna",
                "direction": "wgcna",
                "neglog10_fdr": -np.log10(p + 1e-300),
                "category": cat,
            })
    return rows


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def select_top_terms(df: pd.DataFrame, top_n: int) -> list:
    """Keep the top_n Level-2 terms by best (max) significance reached in
    either panel -- matches the original figure's term shortlist size."""
    best = df.groupby("term")["neglog10_fdr"].max().sort_values(ascending=False)
    return best.head(top_n).index.tolist()


def plot_dual_panel(df: pd.DataFrame, output: Path, top_n_terms: int | None):
    df = df[df["category"].isin(CATEGORY_COLORS)].copy()
    df["term"] = df["term"].str.replace("_", " ", regex=False)

    if top_n_terms:
        keep = set(select_top_terms(df, top_n_terms))
        df = df[df["term"].isin(keep)].copy()

    # Term order: grouped by category (CC -> MF -> BP), each group ordered by
    # descending max significance (most-significant terms first within a block).
    term_cat = df.groupby("term")["category"].first()
    term_best = df.groupby("term")["neglog10_fdr"].max()
    term_order = []
    for cat in CATEGORY_ORDER:
        cat_terms = [t for t in term_cat.index if term_cat[t] == cat]
        cat_terms.sort(key=lambda t: -term_best[t])
        term_order.extend(cat_terms)
    term_idx = {t: i for i, t in enumerate(term_order)}

    deseq_comparisons = sorted(df.loc[df["panel"] == "deseq2", "source_label"].unique())
    wgcna_modules = sorted(df.loc[df["panel"] == "wgcna", "source_label"].unique())

    n_terms = len(term_order)
    fig_h = max(6, 0.45 * n_terms + 1.5)
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(max(14, 0.55 * len(deseq_comparisons) + 6), fig_h),
        gridspec_kw={"width_ratios": [max(len(deseq_comparisons), 1), max(len(wgcna_modules), 1) * 0.6]},
    )

    # Cap bubble AREA at MAX_BUBBLE pt^2, scaled per panel relative to that
    # panel's own largest -log10(FDR). This is the fix for bubbles
    # overflowing their column: the previous version scaled size linearly
    # with no cap (s = neglog10_fdr * 40), so a very small p-value (large
    # -log10) produced a marker far bigger than the plot area, especially
    # in the narrower WGCNA panel. Matches the original submitted figure's
    # notebook cell (MAX_BUBBLE = 300, size = min((|val|/max_val)*MAX_BUBBLE, MAX_BUBBLE)).
    MAX_BUBBLE = 300

    def draw_panel(ax, sub_df, x_labels, title, xlabel):
        x_idx = {c: i for i, c in enumerate(x_labels)}
        max_val = sub_df["neglog10_fdr"].abs().max() if len(sub_df) else 1
        max_val = max_val or 1
        for _, row in sub_df.iterrows():
            x = x_idx[row["source_label"]]
            y = term_idx[row["term"]]
            s = min((abs(row["neglog10_fdr"]) / max_val) * MAX_BUBBLE, MAX_BUBBLE)
            c = DIRECTION_COLORS[row["direction"]]
            ax.scatter(x, y, s=s, color=c, alpha=0.85, edgecolors="white",
                       linewidths=0.5, zorder=3)
        ax.set_xticks(range(len(x_labels)))
        ax.set_xticklabels(x_labels, rotation=45, ha="right", fontsize=8)
        ax.set_xlim(-0.6, max(len(x_labels), 1) - 0.4)
        ax.set_ylim(-0.6, n_terms - 0.4)
        ax.xaxis.grid(True, linestyle="--", alpha=0.3, zorder=0)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(xlabel, fontsize=10)

    draw_panel(ax1, df[df["panel"] == "deseq2"], deseq_comparisons,
               "DESeq2 GO enrichment", "DESeq2 comparisons")
    draw_panel(ax2, df[df["panel"] == "wgcna"], wgcna_modules,
               "WGCNA GO enrichment", "WGCNA modules")

    ax1.set_yticks(range(n_terms))
    ax1.set_yticklabels(term_order, fontsize=9)
    ax2.set_yticks(range(n_terms))
    ax2.set_yticklabels([])

    # Category colour strip (left of ax1's y-tick labels)
    for i, term in enumerate(term_order):
        ax1.add_patch(mpatches.Rectangle(
            (-0.6, i - 0.4), 0.06, 0.8,
            transform=ax1.get_yaxis_transform(), clip_on=False,
            facecolor=CATEGORY_COLORS[term_cat[term]], edgecolor="none",
        ))

    fig.suptitle("GO enrichment — DESeq2 vs WGCNA (Level-2 terms)", fontsize=13)

    # Legends (on a dedicated area to the right, mirroring the original figure)
    direction_handles = [mpatches.Patch(facecolor=DIRECTION_COLORS[k], label=DIRECTION_LABELS[k])
                          for k in ["up1", "up2", "wgcna"]]
    category_handles = [mpatches.Patch(facecolor=CATEGORY_COLORS[k], label=CATEGORY_LABELS[k])
                         for k in CATEGORY_ORDER]
    # Size-legend handles: Line2D markers, not plt.scatter([], []). An
    # empty PathCollection (what plt.scatter([], []) returns) crashes during
    # fig.savefig(..., bbox_inches="tight") in this matplotlib version --
    # bbox_inches="tight" calls get_extents() on every legend handle
    # (including the invisible "blank" spacer handles), and an empty
    # PathCollection's Path.get_extents() does
    # `np.concatenate(xys)` on an empty list, raising
    # "ValueError: need at least one array to concatenate". Line2D handles
    # with no data don't hit that code path and render identically for
    # legend purposes. markersize is a diameter (points), unlike scatter's
    # area-based `s` (points^2), hence the sqrt conversion below.
    deseq_max_val = df.loc[df["panel"] == "deseq2", "neglog10_fdr"].abs().max() or 1
    size_handles = [
        Line2D([], [], marker="o", linestyle="",
               markerfacecolor="grey", markeredgecolor="white", alpha=0.7,
               markersize=min((v / deseq_max_val) * MAX_BUBBLE, MAX_BUBBLE) ** 0.5,
               label=f"−log₁₀(FDR) = {v}")
        for v in [1, 2, 4]
    ]

    blank = Line2D([], [], linestyle="")
    all_handles = ([blank] + direction_handles +
                   [blank, blank] + category_handles +
                   [blank, blank] + size_handles)
    all_labels = (["— Direction —"] + [h.get_label() for h in direction_handles] +
                  ["", "— GO category —"] + [h.get_label() for h in category_handles] +
                  ["", "— Significance —"] + [h.get_label() for h in size_handles])

    fig.legend(
        handles=all_handles, labels=all_labels,
        loc="upper left", bbox_to_anchor=(0.985, 0.95), frameon=True, fontsize=9,
        handletextpad=0.8, labelspacing=0.6,
    )

    plt.tight_layout(rect=[0, 0, 0.98, 0.96])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Figure 4C saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Figure 4C: dual-panel GO enrichment bubble plot (DESeq2 vs WGCNA, Level-2 terms).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--deseq-dir", required=True, type=Path,
                   help="Directory with one <cond_A>_vs_<cond_B>/ sub-folder per comparison, "
                        "each containing GO_enrichment_up_in_<cond>*.csv (run_go_enrichment.py output).")
    p.add_argument("--wgcna-dir", required=True, type=Path,
                   help="Directory with enrichment_<module>.csv files (GO_ID, P_adj_BH columns).")
    p.add_argument("--dictionary", required=True, type=Path, help="dictionary_2_level.csv")
    p.add_argument("--ancestor-cache", required=True, type=Path,
                   help="JSON cache of GO id -> list of ancestor GO ids (from QuickGO).")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--top-n-terms", type=int, default=15,
                   help="Number of most-significant Level-2 GO terms to display "
                        "(0 = show all). Default 15 matches the original figure.")
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    go_dict, go_category = load_level2_dictionary(args.dictionary)
    ancestors = json.loads(args.ancestor_cache.read_text())

    rows = []
    if args.deseq_dir.exists():
        rows.extend(build_deseq2_rows(args.deseq_dir, go_dict, go_category, ancestors, args.alpha))
    else:
        print(f"[WARN] --deseq-dir not found: {args.deseq_dir}")

    if args.wgcna_dir.exists():
        rows.extend(build_wgcna_rows(args.wgcna_dir, go_dict, go_category, args.alpha))
    else:
        print(f"[WARN] --wgcna-dir not found: {args.wgcna_dir}")

    if not rows:
        sys.exit("[ERROR] No enrichment data could be loaded.")

    df = pd.DataFrame(rows)
    csv_out = args.output.with_suffix(".csv")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(csv_out, index=False)
    print(f"[INFO] {df['term'].nunique()} Level-2 terms, "
          f"{df.loc[df.panel == 'deseq2', 'source_label'].nunique()} DESeq2 comparisons, "
          f"{df.loc[df.panel == 'wgcna', 'source_label'].nunique()} WGCNA modules "
          f"-> {csv_out}")

    plot_dual_panel(df, args.output, args.top_n_terms or None)


if __name__ == "__main__":
    main()
