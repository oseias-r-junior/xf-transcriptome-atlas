"""
run_deseq2.py — Differential expression analysis with PyDESeq2 for all
pairwise comparisons defined in a metadata file.

This reproduces Cells 3-4 and 7-8 of the original analysis notebook
(deseq2_pierry_feitosa_paper.ipynb) as faithfully as possible, so that
re-running this script on the raw inputs reproduces the DEG counts actually
submitted to MDPI Pathogens. Two points earlier versions of this script got
wrong are fixed here:

1. Gene universe per comparison (ortholog handling)
   --------------------------------------------------
   `--counts` (data/raw_counts_combined.tsv) is a strain-native, OUTER-JOIN
   matrix (see 01_preprocessing/build_raw_counts_combined.py): a 9a5c gene's
   columns are 0 for every Temecula1 sample, and vice-versa, because that
   gene ID simply never appears in a Temecula1 count file. No ortholog
   merge has happened yet.
     - WITHIN-strain comparisons (both conditions map to the same `strain`
       in --metadata): samples are subset, then rows that sum to 0 across
       *those* columns are dropped. Because the other strain's genes are
       structurally all-zero there, this alone restores that strain's full
       native gene universe -- no ortholog-availability restriction, no
       explicit filtering needed.
     - CROSS-strain comparisons (conditions map to different `strain`
       values): genes are mapped Temecula1 IMG ID -> 9a5c IMG ID via
       --dictionary (summed if many-to-one, matching the original
       safe_map_gene_ids()), then intersected on the common (now-9a5c-IMG)
       gene_id before the two strains' sample columns are concatenated.
   This mirrors run_go_enrichment.py's documented "Background strategy"
   (within-strain: all genes detected in that strain; cross-strain: union
   of orthologous gene pairs) -- that distinction was already the repo's
   stated convention for GO backgrounds, it just hadn't been applied to the
   DE step itself until now. Cross/within status is read from --metadata's
   `strain` column (not string-matched on condition names), since condition
   labels use inconsistent strain abbreviations across files (e.g. "Tem1"
   vs "Temecula1").

2. The statistical test itself
   -----------------------------
   The original notebook does NOT run a plain Wald test against LFC = 0
   and filter afterwards. It runs PyDESeq2's TREAT-style test against the
   LFC threshold itself:
       DeseqStats(dds, contrast=["group","group1","group2"], alpha=ALPHA,
                  lfc_null=LFC_CUTOFF, alt_hypothesis="greaterAbs",
                  cooks_filter=True, independent_filter=True)
   with `design="group"` (group2, group1 as levels, in that order) rather
   than a "_cond" factor with the raw condition strings. DEGs are then
   `padj <= alpha & abs(log2FoldChange) > lfc_cutoff` -- confirmed by
   reading Cell 8 of deseq2_pierry_feitosa_paper.ipynb directly (the
   `deg_mask` there uses `res_df["padj"] <= ALPHA`, non-strict). A
   previous version of this script switched this to a strict `padj <
   alpha` based on an unverified claim that separate "original R scripts"
   used strict `<`; no such R scripts exist in this repo or were located,
   and the actual Python notebook -- which this script's own docstring
   above states it reproduces "as faithfully as possible" -- is the
   ground truth and uses `<=`. That earlier "fix" was itself a regression:
   it silently dropped genes with padj exactly at the 0.05 boundary,
   which explains small (~1-gene) DEG-count deficits seen against the
   submitted Table 2 (e.g. 86 vs 87 for 9a5c_PIM6_1d_vs_9a5c_PIM6_3d).
   Reverted to `<=` here to match Cell 8 exactly. The LFC side of the
   filter (`abs(log2FoldChange) > lfc_cutoff`, strict) does match Cell 8
   as well (`res_df["log2FoldChange"].abs() > LFC_CUTOFF`) and is
   unchanged.

Usage
-----
python run_deseq2.py \\
    --counts      data/raw_counts_combined.tsv \\
    --metadata    data/sample_info.csv \\
    --dictionary  data/gene_dictionary.tsv \\
    --comparisons data/comparisons.tsv \\
    --outdir      results/DESeq2_results \\
    --alpha       0.05 \\
    --lfc         1.0

Input formats
-------------
Counts matrix (genes x samples, TSV, strain-native IDs, outer-joined,
zero-filled -- see build_raw_counts_combined.py):
    gene_id       sample_A_1   sample_A_2   sample_B_1   ...
    XF9a_00001    245          312          0            ...

Sample metadata (CSV or TSV, must include `condition` and `strain`):
    sample_id    condition        strain   medium  timepoint
    sample_A_1   9a5c_PIM6_1d     9a5c     PIM6    1d
    ...

Comparisons file (TSV, two columns, no header):
    9a5c_PIM6_1d    9a5c_PIM6_3d
    9a5c_PIM6_1d    Temecula1_PIM6_1d
    ...

Gene dictionary (from build_gene_dictionary.py): used both (a) to map
Temecula1 IMG ID -> 9a5c IMG ID for cross-strain merging, and (b) to
annotate DEG tables with old locus tags for cross-strain interpretation.

Output (per comparison under --outdir/<c1>_vs_<c2>/):
    DESeq2_results_full.csv            — complete results table
    DEGs_annotated_p<alpha>_LFC<lfc>.csv — significant DEGs with annotations
    deseq2_normalized_counts.csv       — size-factor-normalized counts
    size_factors.csv
    vst_counts.csv
    comparisons_summary_formatted.csv  — one row per comparison (root outdir)
"""

import argparse
import re
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Dictionary helpers
# ---------------------------------------------------------------------------

def load_dictionary(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, sep=sep)


def _choose_old_tag(raw: str) -> str:
    """From a space/comma/semicolon/pipe-separated list of old locus tags,
    pick the canonical form (no underscore when available)."""
    if not isinstance(raw, str) or not raw.strip():
        return ""
    parts = re.split(r"[\s,;|]+", raw.strip())
    for p in parts:
        if p and "_" not in p:
            return p
    return parts[0] if parts else ""


def build_id_maps(dictionary: pd.DataFrame):
    """Return (tem_to_9a, img_to_tem, img_to_9a):
    - tem_to_9a: Temecula1_IMG_ID -> 9a5c_IMG_ID (for cross-strain merging)
    - img_to_tem / img_to_9a: either strain's IMG_ID -> old locus tag
      (for annotating DEG tables)."""
    tem_to_9a = dict(zip(dictionary.get("Temecula1_IMG_ID", []),
                          dictionary.get("9a5c_IMG_ID", [])))

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
    return tem_to_9a, img_to_tem, img_to_9a


def safe_map_gene_ids(counts_sub: pd.DataFrame, mapper: dict) -> pd.DataFrame:
    """Map gene IDs in counts_sub (index = gene IDs) via mapper (old -> new),
    summing rows that collapse many-to-one. Unmapped genes are dropped.
    Identical logic to Cell 7's safe_map_gene_ids()."""
    mapped = counts_sub.reset_index()
    gene_col = mapped.columns[0]
    mapped = mapped.rename(columns={gene_col: "original_gene_id"})
    mapped["gene_id"] = mapped["original_gene_id"].map(mapper)
    mapped = mapped.dropna(subset=["gene_id"]).drop(columns=["original_gene_id"])
    mapped = mapped.set_index("gene_id").groupby("gene_id").sum()
    return mapped


# ---------------------------------------------------------------------------
# PyDESeq2 wrapper
# ---------------------------------------------------------------------------

def run_deseq2_comparison(
    sub_counts: pd.DataFrame,   # genes x samples, already strain-resolved
    sub_meta: pd.DataFrame,     # samples x traits, includes `condition`
    c1: str,
    c2: str,
    condition_col: str,
    alpha: float,
    lfc_cutoff: float,
):
    """Run PyDESeq2 exactly as Cell 8 does: design='group' with levels
    [group2, group1] (group2 = c2 = reference), TREAT-style test against
    lfc_null via alt_hypothesis='greaterAbs'."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    from pydeseq2.preprocessing import deseq2_norm

    sub_meta = sub_meta.copy()
    sub_meta["group"] = np.where(sub_meta[condition_col] == c1, "group1", "group2")
    sub_meta["group"] = pd.Categorical(sub_meta["group"], categories=["group2", "group1"])

    dds = DeseqDataSet(
        counts=sub_counts.T,      # samples x genes
        metadata=sub_meta,
        design="group",
        refit_cooks=True,
    )
    dds.deseq2()

    stat = DeseqStats(
        dds,
        contrast=["group", "group1", "group2"],   # c1 vs c2
        alpha=alpha,
        lfc_null=lfc_cutoff,
        alt_hypothesis="greaterAbs",
        cooks_filter=True,
        independent_filter=True,
    )
    stat.summary()

    res = stat.results_df.copy()
    if len(res) != sub_counts.shape[0]:
        raise RuntimeError("Mismatch between result rows and gene IDs.")
    res.index = sub_counts.index
    res.index.name = "gene_id"
    res = res.reset_index()

    norm_df = sizef_df = vst_df = None
    try:
        norm_counts, size_factors = deseq2_norm(sub_counts.T)
        norm_df = pd.DataFrame(norm_counts, index=sub_counts.columns, columns=sub_counts.index)
        sizef_df = pd.DataFrame(size_factors, index=sub_counts.columns, columns=["size_factor"])
    except Exception:
        pass
    try:
        dds.vst(use_design=False, fit_type=None)
        vst_df = pd.DataFrame(dds.layers["vst_counts"],
                               index=sub_counts.columns, columns=sub_counts.index)
    except Exception:
        pass

    return res, norm_df, sizef_df, vst_df


# ---------------------------------------------------------------------------
# Per-comparison gene-universe resolution (within-strain vs cross-strain)
# ---------------------------------------------------------------------------

def resolve_comparison_counts(
    counts_df: pd.DataFrame,
    sample_df: pd.DataFrame,
    c1: str,
    c2: str,
    condition_col: str,
    tem_to_9a: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
    """Subset counts_df/sample_df to samples in {c1, c2}. If c1/c2 belong to
    different `strain` values, merge Temecula1 IMG IDs -> 9a5c IMG IDs
    (summed, many-to-one) and intersect on common genes before concatenating
    the two strains' sample columns -- exactly Cell 8's is_strain_comparison
    branch. If c1/c2 are the same strain, just drop rows that are all-zero
    across the selected columns (restores that strain's full native gene
    set, since the other strain's genes are structurally zero here)."""
    mask = sample_df[condition_col].isin([c1, c2])
    sub_samples = sample_df.loc[mask].copy()
    if sub_samples.shape[0] == 0:
        raise RuntimeError(f"No samples found for {c1!r} or {c2!r}")

    sub_counts = counts_df.loc[:, sub_samples.index]
    sub_counts = sub_counts.loc[sub_counts.sum(axis=1) > 0, :]
    if sub_counts.shape[0] == 0:
        raise RuntimeError("No genes left after filtering zero-sum across selected samples")

    strains_present = sub_samples["strain"].unique().tolist()
    is_cross_strain = len(strains_present) > 1

    if is_cross_strain:
        tem_samples = sub_samples[sub_samples["strain"] == "Temecula1"].index.tolist()
        a9_samples  = sub_samples[sub_samples["strain"] == "9a5c"].index.tolist()
        if not tem_samples or not a9_samples:
            raise RuntimeError(
                f"Cross-strain comparison {c1!r} vs {c2!r} but strain values "
                f"present are {strains_present!r} (expected both '9a5c' and "
                f"'Temecula1')."
            )

        counts_tem = sub_counts.loc[:, tem_samples]
        counts_a9  = sub_counts.loc[:, a9_samples]

        counts_tem_mapped = safe_map_gene_ids(counts_tem, tem_to_9a)
        common = counts_a9.index.intersection(counts_tem_mapped.index)
        if len(common) == 0:
            raise RuntimeError("No common genes after Temecula1->9a5c ortholog mapping")

        sub_counts = pd.concat([counts_a9.loc[common], counts_tem_mapped.loc[common]], axis=1)
        sub_samples = sub_samples.loc[sub_counts.columns]
        sub_counts = sub_counts.loc[sub_counts.sum(axis=1) > 0, :]
        if sub_counts.shape[0] == 0:
            raise RuntimeError("No genes remain after ortholog mapping + filtering")

    return sub_counts, sub_samples, is_cross_strain


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

def parse_comparisons(path: Path) -> list[tuple[str, str]]:
    df = pd.read_csv(path, sep="\t", header=None)
    return [(str(r[0]).strip(), str(r[1]).strip()) for _, r in df.iterrows()]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run PyDESeq2 for all pairwise comparisons, reproducing "
                     "the original notebook's cross-strain ortholog handling "
                     "and TREAT-style (greaterAbs / lfc_null) test.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--counts",      required=True, type=Path,
                   help="Strain-native, outer-joined counts matrix "
                        "(build_raw_counts_combined.py output).")
    p.add_argument("--metadata",    required=True, type=Path,
                   help="sample_info.csv/tsv; must include `condition` and `strain`.")
    p.add_argument("--dictionary",  required=True, type=Path,
                   help="Gene dictionary (build_gene_dictionary.py output).")
    p.add_argument("--comparisons", required=True, type=Path,
                   help="Two-column TSV: condition_1, condition_2.")
    p.add_argument("--condition-col", default="condition",
                   help="Column in metadata that holds condition labels.")
    p.add_argument("--sample-id-col", default=None,
                   help="Column to use as sample index if metadata isn't "
                        "already indexed by sample_id (default: first column).")
    p.add_argument("--outdir",      required=True, type=Path)
    p.add_argument("--alpha",       type=float, default=0.05,
                   help="FDR threshold (passed to DeseqStats as `alpha`).")
    p.add_argument("--lfc",         type=float, default=1.0,
                   help="|log2FC| threshold (passed to DeseqStats as `lfc_null`; "
                        "final DEG filter uses strict `>`, matching the "
                        "original notebook).")
    p.add_argument("--min-replicates", type=int, default=2,
                   help="Minimum samples required per condition.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # -- Load inputs ----------------------------------------------------
    counts_sep = "\t" if args.counts.suffix in {".tsv", ".txt"} else ","
    counts_df = pd.read_csv(args.counts, sep=counts_sep, index_col=0)

    meta_sep = "\t" if args.metadata.suffix in {".tsv", ".txt"} else ","
    if args.sample_id_col:
        sample_df = pd.read_csv(args.metadata, sep=meta_sep).set_index(args.sample_id_col)
    else:
        sample_df = pd.read_csv(args.metadata, sep=meta_sep, index_col=0)

    if "strain" not in sample_df.columns:
        sys.exit("[ERROR] --metadata must include a `strain` column "
                  "(e.g. '9a5c' / 'Temecula1') to resolve within- vs "
                  "cross-strain gene universes.")

    comparisons = parse_comparisons(args.comparisons)
    dictionary = load_dictionary(args.dictionary)
    tem_to_9a, img_to_tem, img_to_9a = build_id_maps(dictionary)

    print(f"[INFO] {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
    print(f"[INFO] {len(comparisons)} comparisons to run")

    summary_rows = []

    for c1, c2 in comparisons:
        comp_name = f"{c1}_vs_{c2}".replace("/", "_")
        comp_dir = args.outdir / comp_name
        comp_dir.mkdir(parents=True, exist_ok=True)
        print(f"[INFO] Running: {comp_name}")

        n_genes_tested = n_DEG = n_up = n_down = 0
        error_msg = None

        try:
            group_counts = sample_df.loc[
                sample_df[args.condition_col].isin([c1, c2]), args.condition_col
            ].value_counts().to_dict()
            if group_counts.get(c1, 0) < args.min_replicates or group_counts.get(c2, 0) < args.min_replicates:
                raise RuntimeError(f"Not enough replicates per group: {group_counts}")

            sub_counts, sub_samples, is_cross_strain = resolve_comparison_counts(
                counts_df, sample_df, c1, c2, args.condition_col, tem_to_9a,
            )
            n_genes_tested = sub_counts.shape[0]
            print(f"  [{'cross-strain, ortholog-merged' if is_cross_strain else 'within-strain, native genes'}] "
                  f"{n_genes_tested} genes, {sub_counts.shape[1]} samples")

            res, norm_df, sizef_df, vst_df = run_deseq2_comparison(
                sub_counts, sub_samples, c1, c2, args.condition_col, args.alpha, args.lfc,
            )

            res.to_csv(comp_dir / "DESeq2_results_full.csv", index=False)

            sig = res[
                res["padj"].notna()
                & (res["padj"] <= args.alpha)   # non-strict <=, matches Cell 8's
                                                 # `res_df["padj"] <= ALPHA` exactly
                & (res["log2FoldChange"].abs() > args.lfc)   # strict >, matches `res_df["log2FoldChange"].abs() > LFC_CUTOFF`
            ].copy()

            sig["tem_old_tag"] = sig["gene_id"].map(lambda g: img_to_tem.get(str(g), ""))
            sig["9a5c_old_tag"] = sig["gene_id"].map(lambda g: img_to_9a.get(str(g), ""))

            lfc_str = str(args.lfc).replace(".", "p")
            sig_path = comp_dir / f"DEGs_annotated_p{args.alpha}_LFC{lfc_str}.csv"
            sig.to_csv(sig_path, index=False)

            if norm_df is not None:
                norm_df.to_csv(comp_dir / "deseq2_normalized_counts.csv")
                sizef_df.to_csv(comp_dir / "size_factors.csv")
            if vst_df is not None:
                vst_df.to_csv(comp_dir / "vst_counts.csv")

            n_DEG = len(sig)
            n_up = int((sig["log2FoldChange"] > 0).sum())
            n_down = int((sig["log2FoldChange"] < 0).sum())
            print(f"  -> {n_DEG} DEGs (up {n_up} / down {n_down})")

        except Exception as e:
            error_msg = str(e)
            (comp_dir / "error.txt").write_text(traceback.format_exc())
            print(f"[WARN] {comp_name} failed: {error_msg}", file=sys.stderr)

        summary_rows.append({
            "Condition 1 (c1)": c1,
            "Condition 2 (c2)": c2,
            "comparison":       comp_name,
            "n_tested":         n_genes_tested,
            "n_DEG":            n_DEG,
            "n_up_in_c1":       n_up,
            "n_down_in_c1":     n_down,
            "error":            error_msg,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_path = args.outdir / "comparisons_summary_formatted.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"[INFO] Summary saved to {summary_path}")


if __name__ == "__main__":
    main()
