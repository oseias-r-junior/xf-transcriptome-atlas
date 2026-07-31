"""
run_deseq2.py — Differential expression analysis with PyDESeq2 for all
pairwise comparisons defined in a metadata file.

For each comparison the script:
  1. Runs PyDESeq2 (Love et al. 2014; Muzellec et al. 2023)
  2. Filters results at the specified FDR and |log2FC| thresholds
  3. Annotates DEGs with old locus tags using the cross-strain gene dictionary
  4. Writes per-comparison result tables and a formatted summary CSV

Usage
-----
python run_deseq2.py \\
    --counts   data/raw_counts_combined.tsv \\
    --metadata data/sample_info.tsv \\
    --dictionary data/gene_dictionary.tsv \\
    --comparisons data/comparisons.tsv \\
    --outdir   results/DESeq2_results \\
    --alpha    0.05 \\
    --lfc      1.0

Input formats
-------------
Counts matrix (genes × samples, TSV):
    gene_id       sample_A_1   sample_A_2   sample_B_1   ...
    Xf9a_00001    245          312          98           ...

Sample metadata (TSV):
    sample_id    condition        strain   medium  timepoint
    sample_A_1   9a5c_PIM6_1d     9a5c     PIM6    1d
    ...

Comparisons file (TSV, two columns, no header):
    9a5c_PIM6_1d    9a5c_PIM6_3d
    9a5c_PIM6_1d    Temecula1_PIM6_1d
    ...

Gene dictionary (from build_gene_dictionary.py): used to add old locus tags
to DEG tables for cross-strain interpretation.

Output (per comparison under --outdir/<c1>_vs_<c2>/):
    DESeq2_results_full.csv       — complete results table
    DEGs_annotated_*.csv          — significant DEGs with annotations
    comparisons_summary_formatted.csv — one row per comparison (root outdir)
"""

import argparse
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Dictionary helpers
# ---------------------------------------------------------------------------

def load_dictionary(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _choose_old_tag(raw: str) -> str:
    """From a space-separated list of old locus tags pick the canonical form
    (no underscore when available, e.g. 'XF0001' over 'XF_0001')."""
    if not isinstance(raw, str) or not raw.strip():
        return ""
    parts = re.split(r"[\s,;|]+", raw.strip())
    for p in parts:
        if p and "_" not in p:
            return p
    return parts[0] if parts else ""


def build_id_maps(dictionary: pd.DataFrame):
    """Return two dicts: img_id → tem_old_tag and img_id → 9a5c_old_tag."""
    img_to_tem, img_to_9a = {}, {}
    for _, row in dictionary.iterrows():
        tem_tag = _choose_old_tag(str(row.get("Temecula1_old_locus_tags", "")))
        x9a_tag = _choose_old_tag(str(row.get("9a5c_old_locus_tags", "")))
        for col in ["9a5c_IMG_ID", "Temecula1_IMG_ID"]:
            v = str(row.get(col, "")).strip()
            if v:
                if tem_tag:
                    img_to_tem[v] = tem_tag
                if x9a_tag:
                    img_to_9a[v] = x9a_tag
    return img_to_tem, img_to_9a


# ---------------------------------------------------------------------------
# PyDESeq2 wrapper
# ---------------------------------------------------------------------------

def run_deseq2_comparison(
    counts: pd.DataFrame,
    meta: pd.DataFrame,
    condition_col: str,
    c1: str,
    c2: str,
    alpha: float,
    lfc_cutoff: float,
) -> pd.DataFrame:
    """Run PyDESeq2 for samples in conditions c1 and c2."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    mask = meta[condition_col].isin([c1, c2])
    sub_meta = meta.loc[mask].copy()
    sub_meta["_cond"] = sub_meta[condition_col].astype(str)

    sub_counts = counts[sub_meta.index].T.astype(int)   # samples × genes

    dds = DeseqDataSet(
        counts=sub_counts,
        metadata=sub_meta,
        design_factors="_cond",
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()

    stat = DeseqStats(dds, contrast=["_cond", c1, c2], quiet=True)
    stat.summary()
    res = stat.results_df.copy()
    res.index.name = "gene_id"
    res = res.reset_index()
    return res


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def parse_comparisons(path: Path) -> list[tuple[str, str]]:
    df = pd.read_csv(path, sep="\t", header=None)
    return [(str(r[0]).strip(), str(r[1]).strip()) for _, r in df.iterrows()]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run PyDESeq2 for all pairwise comparisons.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--counts",      required=True, type=Path)
    p.add_argument("--metadata",    required=True, type=Path)
    p.add_argument("--dictionary",  required=True, type=Path,
                   help="Gene dictionary (build_gene_dictionary.py output).")
    p.add_argument("--comparisons", required=True, type=Path,
                   help="Two-column TSV: condition_1, condition_2.")
    p.add_argument("--condition-col", default="condition",
                   help="Column in metadata that holds condition labels.")
    p.add_argument("--outdir",      required=True, type=Path)
    p.add_argument("--alpha",       type=float, default=0.05,
                   help="FDR threshold.")
    p.add_argument("--lfc",         type=float, default=1.0,
                   help="|log2FC| threshold.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # -- Load inputs ---------------------------------------------------------
    sep = "\t" if args.counts.suffix in {".tsv", ".txt"} else ","
    counts = pd.read_csv(args.counts, sep=sep, index_col=0)
    meta   = pd.read_csv(args.metadata, sep="\t", index_col=0)
    comparisons = parse_comparisons(args.comparisons)
    dictionary  = load_dictionary(args.dictionary)
    img_to_tem, img_to_9a = build_id_maps(dictionary)

    print(f"[INFO] {counts.shape[0]} genes × {counts.shape[1]} samples")
    print(f"[INFO] {len(comparisons)} comparisons to run")

    summary_rows = []

    for c1, c2 in comparisons:
        comp_name = f"{c1}_vs_{c2}"
        comp_dir  = args.outdir / comp_name
        comp_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Running: {comp_name}")

        try:
            res = run_deseq2_comparison(
                counts, meta, args.condition_col, c1, c2, args.alpha, args.lfc
            )
        except Exception as e:
            print(f"[WARN] {comp_name} failed: {e}", file=sys.stderr)
            continue

        # Annotate with old locus tags
        res["tem_old_tag"] = res["gene_id"].map(lambda g: img_to_tem.get(str(g), ""))
        res["9a5c_old_tag"] = res["gene_id"].map(lambda g: img_to_9a.get(str(g), ""))

        # Save full results
        full_path = comp_dir / "DESeq2_results_full.csv"
        res.to_csv(full_path, index=False)

        # Significant DEGs
        sig = res[
            res["padj"].notna()
            & (res["padj"] <= args.alpha)
            & (res["log2FoldChange"].abs() >= args.lfc)
        ].copy()

        lfc_str = str(args.lfc).replace(".", "p")
        sig_path = comp_dir / f"DEGs_annotated_p{args.alpha}_LFC{lfc_str}.csv"
        sig.to_csv(sig_path, index=False)

        n_up   = (sig["log2FoldChange"] > 0).sum()
        n_down = (sig["log2FoldChange"] < 0).sum()
        print(f"  → {len(sig)} DEGs  (↑{n_up}  ↓{n_down})")

        summary_rows.append(
            {
                "Condition 1 (c1)": c1,
                "Condition 2 (c2)": c2,
                "comparison":       comp_name,
                "n_tested":         len(res),
                "n_DEG":            len(sig),
                "n_up_in_c1":       n_up,
                "n_down_in_c1":     n_down,
            }
        )

    # -- Write summary -------------------------------------------------------
    summary_df = pd.DataFrame(summary_rows)
    summary_path = args.outdir / "comparisons_summary_formatted.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[INFO] Summary saved to {summary_path}")


if __name__ == "__main__":
    main()
