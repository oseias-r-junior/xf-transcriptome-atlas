"""
run_go_enrichment.py — GO term enrichment analysis for DESeq2 output.

For each comparison folder produced by run_deseq2.py the script performs:
  1. Splits DEGs into up-regulated (in c1) and down-regulated (in c1)
  2. Runs Fisher's exact test against a background gene set
  3. Applies Benjamini-Hochberg FDR correction
  4. Writes one CSV per direction per comparison

Background strategy:
  - Within-strain comparisons: all genes detected in that strain
  - Cross-strain comparisons:  union of orthologous gene pairs (combined background)

This follows the approach in Feitosa-Junior et al. (2025), Methods section
"GO enrichment and clustering analyses".

Usage
-----
python run_go_enrichment.py \\
    --deseq-dir   results/DESeq2_results \\
    --go-annot    results/go_annotations.tsv \\
    --dictionary  data/gene_dictionary.tsv \\
    --alpha       0.05 \\
    --lfc         1.0

Output (per comparison, written inside --deseq-dir/<comp>/):
    GO_enrichment_up_in_<c1>.csv
    GO_enrichment_up_in_<c2>.csv

Each CSV contains:
    go_term, go_description, deg_count, bg_count,
    deg_pct, bg_pct, fold_enrichment, p_value, p_adj
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_go_annotations(path: Path) -> dict[str, set[str]]:
    """Return gene_id → set of GO terms."""
    df = pd.read_csv(path, sep="\t")
    go_map: dict[str, set] = {}
    for _, row in df.iterrows():
        gid = str(row["gene_id"]).strip()
        go  = str(row["go_term"]).strip()
        go_map.setdefault(gid, set()).add(go)
    return go_map


def load_dictionary(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def unify_id_set(gene_ids: list[str], dictionary: pd.DataFrame) -> set[str]:
    """Expand a list of gene IDs to include all synonyms from the dictionary."""
    img_cols = [c for c in dictionary.columns if "IMG_ID" in c]
    out: set[str] = set(gene_ids)
    for col in img_cols:
        mask = dictionary[col].isin(out)
        for other_col in img_cols:
            out.update(dictionary.loc[mask, other_col].dropna().astype(str))
    return out


def enrichment_test(
    query_genes: set[str],
    background_genes: set[str],
    go_map: dict[str, set[str]],
    alpha: float = 0.05,
) -> pd.DataFrame:
    """Fisher's exact test for GO term enrichment."""
    # Collect all GO terms represented in the background
    all_terms: set[str] = set()
    for g in background_genes:
        all_terms.update(go_map.get(g, set()))

    if not all_terms or not query_genes:
        return pd.DataFrame()

    n_query = len(query_genes)
    n_bg    = len(background_genes)

    rows = []
    for term in all_terms:
        q_with = len([g for g in query_genes    if term in go_map.get(g, set())])
        b_with = len([g for g in background_genes if term in go_map.get(g, set())])

        if q_with == 0:
            continue

        # 2×2 contingency table
        table = [
            [q_with,            n_query - q_with],
            [b_with - q_with,   n_bg - n_query - (b_with - q_with)],
        ]
        _, p = fisher_exact(table, alternative="greater")

        rows.append(
            {
                "go_term":         term,
                "deg_count":       q_with,
                "bg_count":        b_with,
                "deg_pct":         100 * q_with / n_query,
                "bg_pct":          100 * b_with / n_bg,
                "fold_enrichment": (q_with / n_query) / (b_with / n_bg + 1e-9),
                "p_value":         p,
            }
        )

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows).sort_values("p_value")
    _, df["p_adj"], _, _ = multipletests(df["p_value"], method="fdr_bh")
    return df[df["p_adj"] <= alpha].reset_index(drop=True)


# ---------------------------------------------------------------------------
# Scan comparison folders
# ---------------------------------------------------------------------------

def find_comparisons(deseq_dir: Path, lfc: float, alpha: float):
    lfc_str = str(lfc).replace(".", "p")
    for comp_dir in sorted(deseq_dir.iterdir()):
        if not comp_dir.is_dir() or "_vs_" not in comp_dir.name:
            continue
        # Find significant DEG file
        cands = list(comp_dir.glob(f"DEGs_annotated_p{alpha}_LFC{lfc_str}.csv"))
        if not cands:
            cands = list(comp_dir.glob("DEGs_annotated_*.csv"))
        if not cands:
            cands = list(comp_dir.glob("DESeq2_results_full.csv"))
        if cands:
            yield comp_dir.name, comp_dir, cands[0]


def is_cross_strain(comp_name: str) -> bool:
    parts = comp_name.split("_vs_")
    if len(parts) != 2:
        return False
    strain_a = parts[0].split("_")[0]
    strain_b = parts[1].split("_")[0]
    return strain_a != strain_b


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="GO enrichment for DESeq2 output directories.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--deseq-dir",  required=True, type=Path,
                   help="Root directory produced by run_deseq2.py.")
    p.add_argument("--go-annot",   required=True, type=Path,
                   help="GO annotation TSV from parse_go_from_genbank.py.")
    p.add_argument("--dictionary", required=True, type=Path,
                   help="Gene dictionary TSV from build_gene_dictionary.py.")
    p.add_argument("--alpha",      type=float, default=0.05)
    p.add_argument("--lfc",        type=float, default=1.0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    go_map    = load_go_annotations(args.go_annot)
    dictionary = load_dictionary(args.dictionary)

    all_ids_9a  = set(dictionary["9a5c_IMG_ID"].dropna().astype(str))
    all_ids_tem = set(dictionary["Temecula1_IMG_ID"].dropna().astype(str))
    combined_bg = all_ids_9a | all_ids_tem

    comparisons = list(find_comparisons(args.deseq_dir, args.lfc, args.alpha))
    print(f"[INFO] Found {len(comparisons)} comparison folders")

    for comp_name, comp_dir, deg_file in comparisons:
        print(f"[INFO] Processing: {comp_name}")
        df = pd.read_csv(deg_file)

        if "log2FoldChange" not in df.columns:
            print(f"  [SKIP] No log2FoldChange column in {deg_file.name}")
            continue

        if "padj" in df.columns:
            sig = df[df["padj"].notna() & (df["padj"] <= args.alpha)
                     & (df["log2FoldChange"].abs() >= args.lfc)]
        else:
            sig = df  # assume file is already filtered

        gene_col = "gene_id" if "gene_id" in df.columns else df.columns[0]
        up_genes   = set(sig.loc[sig["log2FoldChange"] > 0, gene_col].astype(str))
        down_genes = set(sig.loc[sig["log2FoldChange"] < 0, gene_col].astype(str))

        bg = combined_bg if is_cross_strain(comp_name) else (
            all_ids_9a if "9a5c" in comp_name else all_ids_tem
        )

        c1, c2 = comp_name.split("_vs_", 1)

        for direction, genes, label in [
            ("up_in_c1", up_genes,   c1),
            ("up_in_c2", down_genes, c2),
        ]:
            result = enrichment_test(genes, bg, go_map, alpha=args.alpha)
            if result.empty:
                print(f"  No significant GO terms for {direction}")
                continue
            out = comp_dir / f"GO_enrichment_up_in_{label}.csv"
            result.to_csv(out, index=False)
            print(f"  [{label}] {len(result)} significant terms → {out.name}")


if __name__ == "__main__":
    main()
