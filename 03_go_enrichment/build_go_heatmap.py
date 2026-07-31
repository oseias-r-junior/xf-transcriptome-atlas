"""
build_go_heatmap.py — Integrate GO enrichment results across all DESeq2
comparisons into a signed heatmap (Figure 5 / Supplementary).

Each comparison contributes a signed score column:
    +score = GO term enriched in up-regulated DEGs (c1 > c2)
    −score = GO term enriched in down-regulated DEGs (c2 > c1)
    score  = −log10(p_adj)

GO terms are mapped to their second-level ontology ancestors via the QuickGO
REST API and grouped by ontology category (Biological Process / Molecular
Function / Cellular Component).  Within each category rows are ordered by
hierarchical clustering.

Usage
-----
python build_go_heatmap.py \\
    --deseq-dir results/DESeq2_results \\
    --go-dict   data/dictionary_2_level.csv \\
    --cache     results/go_ancestor_cache.json \\
    --output    figures/go_heatmap.tiff

GO level-2 dictionary (--go-dict): CSV produced by the goatools-based
pipeline.  Four-column, no header:
    level, category (biological_process|...), GO_id, term_label, relation

QuickGO API is used to fetch ancestors on first run and results are cached.
"""

import argparse
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage, leaves_list


# ---------------------------------------------------------------------------
# QuickGO ancestor cache
# ---------------------------------------------------------------------------

def load_cache(path: Path) -> dict:
    if path and path.exists():
        with open(path) as f:
            return json.load(f)
    return {}


def save_cache(cache: dict, path: Path):
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(cache, f)


def get_ancestors(go_id: str, cache: dict, retries: int = 3) -> list[str]:
    if go_id in cache:
        return cache[go_id]
    number = go_id.replace("GO:", "")
    url = (
        f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO%3A"
        f"{number}/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates"
    )
    for attempt in range(retries):
        try:
            r = requests.get(url, headers={"Accept": "application/json"}, timeout=20)
            r.raise_for_status()
            data = r.json()
            ancestors = []
            for res in data.get("results", []):
                ancestors.extend(res.get("ancestors", []))
            cache[go_id] = ancestors
            return ancestors
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
            else:
                print(f"[WARN] QuickGO failed for {go_id}: {e}", file=sys.stderr)
                cache[go_id] = []
                return []


# ---------------------------------------------------------------------------
# GO level-2 dictionary
# ---------------------------------------------------------------------------

def load_go_dict(path: Path) -> tuple[dict, dict]:
    """Return (go_id → term_label, go_id → category_code P/F/C)."""
    df = pd.read_csv(path, header=None,
                     names=["level", "category_name", "GO_id", "term", "relation"])
    df = df[df["relation"] == "is_a"]
    cat_map = {
        "biological_process": "P",
        "molecular_function": "F",
        "cellular_component": "C",
    }
    go_dict = dict(zip(df["GO_id"], df["term"]))
    go_cat  = {go: cat_map.get(cat, "NA")
               for go, cat in zip(df["GO_id"], df["category_name"])}
    return go_dict, go_cat


# ---------------------------------------------------------------------------
# Build signed column for one comparison
# ---------------------------------------------------------------------------

def map_to_level2(
    go_pvals: dict[str, float],
    go_dict: dict,
    go_cat: dict,
    cache: dict,
) -> tuple[dict, dict]:
    """Map raw GO p-values to level-2 ancestors, keeping the best (min) p."""
    term_pvals: dict[str, list] = defaultdict(list)
    term_cat: dict[str, str] = {}
    valid = set(go_dict.keys())

    for go, pval in go_pvals.items():
        for anc in [go] + get_ancestors(go, cache):
            if anc in valid:
                term = go_dict[anc]
                term_pvals[term].append(pval)
                cat = go_cat.get(anc, "NA")
                if cat != "NA":
                    term_cat[term] = cat

    return {t: min(ps) for t, ps in term_pvals.items()}, term_cat


def build_signed_column(
    go_up: dict, go_down: dict, go_dict: dict, go_cat: dict, cache: dict
) -> tuple[dict, dict]:
    up_terms,   up_cat   = map_to_level2(go_up,   go_dict, go_cat, cache)
    down_terms, down_cat = map_to_level2(go_down, go_dict, go_cat, cache)

    scores: dict[str, float] = {}
    cats:   dict[str, str]   = {}

    for term, p in up_terms.items():
        scores[term] = -np.log10(p + 1e-300)
        cats[term]   = up_cat.get(term, "NA")

    for term, p in down_terms.items():
        score = -np.log10(p + 1e-300) * -1
        if term in scores:
            if abs(score) > abs(scores[term]):
                scores[term] = score
        else:
            scores[term] = score
        cats.setdefault(term, down_cat.get(term, "NA"))

    return scores, cats


# ---------------------------------------------------------------------------
# Matrix construction
# ---------------------------------------------------------------------------

def build_matrix(deseq_dir: Path, go_dict: dict, go_cat: dict, cache: dict
                 ) -> tuple[pd.DataFrame, dict]:
    all_columns: dict[str, dict] = {}
    global_cat: dict[str, str] = {}

    for comp_dir in sorted(deseq_dir.iterdir()):
        if not comp_dir.is_dir() or "_vs_" not in comp_dir.name:
            continue
        files = [f for f in comp_dir.iterdir() if "GO_enrichment_up_in" in f.name]
        c1, c2 = comp_dir.name.split("_vs_", 1)
        file_c1 = [f for f in files if c1 in f.name]
        file_c2 = [f for f in files if c2 in f.name]
        if not file_c1 or not file_c2:
            continue

        def load(p):
            df = pd.read_csv(p)
            col = next((c for c in df.columns if "p_adj" in c.lower()), None)
            if col is None:
                return {}
            df = df[df[col] < 0.05]
            return dict(zip(df["go_term"], df[col]))

        go_up   = load(file_c1[0])
        go_down = load(file_c2[0])
        col, cat = build_signed_column(go_up, go_down, go_dict, go_cat, cache)
        all_columns[comp_dir.name] = col
        global_cat.update({k: v for k, v in cat.items() if v != "NA"})

    df = pd.DataFrame(all_columns).fillna(0)
    return df, global_cat


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

CAT_COLORS = {"P": "#e95050", "F": "#4e904a", "C": "#3367e1"}

CMAP = LinearSegmentedColormap.from_list(
    "signed_go", ["#2c7bb6", "#f7f7f7", "#b2182b"]
)


def cluster_rows(df: pd.DataFrame) -> pd.DataFrame:
    if len(df) <= 1:
        return df
    Z = linkage(df.values, method="average", metric="euclidean")
    return df.iloc[leaves_list(Z)]


def plot_heatmap(df: pd.DataFrame, global_cat: dict, output: Path):
    df_bp = cluster_rows(df.loc[[t for t in df.index if global_cat.get(t) == "P"]])
    df_mf = cluster_rows(df.loc[[t for t in df.index if global_cat.get(t) == "F"]])
    df_cc = cluster_rows(df.loc[[t for t in df.index if global_cat.get(t) == "C"]])
    df_ordered = pd.concat([df_bp, df_mf, df_cc])

    if df_ordered.empty:
        print("[WARN] No terms to plot after filtering.", file=sys.stderr)
        return

    row_colors = pd.Series(
        [CAT_COLORS.get(global_cat.get(t, "NA"), "#cccccc") for t in df_ordered.index],
        index=df_ordered.index,
    )

    vmax = df_ordered.abs().max().max()
    vmin = -vmax

    fig_h = max(6, 0.3 * len(df_ordered))
    fig_w = max(10, 0.5 * len(df_ordered.columns))

    sns.set(font_scale=0.85)
    g = sns.clustermap(
        df_ordered,
        cmap=CMAP,
        center=0,
        vmin=vmin,
        vmax=vmax,
        row_colors=row_colors,
        row_cluster=False,
        col_cluster=True,
        linewidths=0.0,
        figsize=(fig_w, fig_h),
        xticklabels=True,
        yticklabels=True,
        cbar_pos=(0.02, 0.82, 0.015, 0.12),
    )

    g.ax_heatmap.set_xticklabels(
        [t.get_text() for t in g.ax_heatmap.get_xticklabels()],
        rotation=35, ha="right", fontsize=8,
    )
    g.ax_heatmap.set_yticklabels(
        [t.get_text() for t in g.ax_heatmap.get_yticklabels()],
        rotation=0, fontsize=8,
    )
    g.cax.set_ylabel("±log₁₀(FDR)", fontsize=8, rotation=270, labelpad=10)

    handles = [Patch(facecolor=v, label=k)
               for k, v in {"Biological Process (P)": "#e95050",
                            "Molecular Function (F)": "#4e904a",
                            "Cellular Component (C)": "#3367e1"}.items()]
    g.ax_row_dendrogram.legend(handles=handles, loc="center",
                               fontsize="small", frameon=True)

    output.parent.mkdir(parents=True, exist_ok=True)
    g.fig.savefig(output, dpi=300, bbox_inches="tight")
    print(f"[INFO] GO heatmap saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build integrative GO enrichment heatmap.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--deseq-dir", required=True, type=Path)
    p.add_argument("--go-dict",   required=True, type=Path,
                   help="Level-2 GO dictionary CSV (dictionary_2_level.csv).")
    p.add_argument("--cache",     type=Path, default=None,
                   help="JSON file for QuickGO ancestor cache (read+write).")
    p.add_argument("--output",    required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cache = load_cache(args.cache)
    go_dict, go_cat = load_go_dict(args.go_dict)

    print(f"[INFO] Level-2 GO dictionary: {len(go_dict)} terms")
    print(f"[INFO] Ancestor cache: {len(cache)} entries")

    df, global_cat = build_matrix(args.deseq_dir, go_dict, go_cat, cache)
    save_cache(cache, args.cache)

    print(f"[INFO] Matrix: {df.shape[0]} GO terms × {df.shape[1]} comparisons")
    plot_heatmap(df, global_cat, args.output)

    # Save matrix
    matrix_out = args.output.with_suffix(".tsv")
    df.to_csv(matrix_out, sep="\t")
    print(f"[INFO] Matrix saved to {matrix_out}")


if __name__ == "__main__":
    main()
