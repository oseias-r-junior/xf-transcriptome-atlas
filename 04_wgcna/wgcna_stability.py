"""
wgcna_stability.py — Leave-one-out module stability analysis.

With only 24 transcriptomes and 640 selected genes, WGCNA module
assignments could in principle be driven by one or a few influential
samples. This script provides direct evidence for (or against) that by
checking whether the module structure holds up when any single sample is
removed.

Method
------
For each of the 24 samples, remove it and re-run the FULL WGCNA pipeline
(same top-640-variable-gene selection, same power/networkType/TOMType/
minModuleSize/MEDissThres/deepSplit as the main analysis) on the remaining
23 samples. After each sub-run, match the resulting modules back to the
three reference modules (dimgrey / darkgrey / silver) by maximum gene-set
overlap (Jaccard index), then report:
  - number of modules recovered per sub-run
  - Jaccard overlap of each recovered module against its matched reference
  - fraction of reference hub genes retained in the matched sub-run module

REQUIRES PyWGCNA — NOT reproducible in a plain NumPy/pandas environment
-------------------------------------------------------------------------
Unlike permanova_factors.py and module_trait_heatmap.py in this repo (which
are pure NumPy/pandas so they run anywhere), this script calls WGCNA's
dynamic tree cut algorithm (dynamicTreeCut / cutreeHybrid), which is a
nontrivial hybrid dendrogram+dissimilarity procedure. Reimplementing it by
hand risks silently producing DIFFERENT module boundaries than PyWGCNA
would — which would make a reported "stability" number scientifically
meaningless. Run this script in the same PyWGCNA environment used for the
main analysis (see environment.yml) rather than approximating it.

IMPORTANT: --tpm must be an ortholog-merged matrix, not the strain-native
data/tpm_expression.csv -- see the long comment in 04_wgcna/run_wgcna.py's
own docstring for why, and data/README.md's "WGCNA input matrix" section
for which of data/tpm_expression_original.csv (the file that actually
produced the manuscript's WGCNA result -- use this one to reproduce/extend
it) vs data/tpm_expression_wgcna.csv (rebuilt from raw data, ~3% fewer
genes) to use.

Usage
-----
python wgcna_stability.py \\
    --tpm       data/tpm_expression_original.csv \\
    --metadata  data/sample_info_original.csv \\
    --top-n     640 \\
    --power     <same beta as main run, e.g. from soft_threshold_plot.tiff> \\
    --reference-assignments results/WGCNA/module_assignments.csv \\
    --outdir    results/WGCNA_stability
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def sniff_sep(path: Path) -> str:
    """Detect the field separator from the file's own header line.

    File extension alone is unreliable: compute_tpm.py writes tab-separated
    data into a file named tpm_expression.csv. Reading it with sep=","
    (an extension-based guess) silently collapses every row into a single
    column, which then shows up downstream as "N genes x 0 samples (after
    alignment)" -- see the same bug and fix in 04_wgcna/run_wgcna.py.
    """
    with open(path, "r", newline="") as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ("\t", ",", ";")}
    best_sep, best_count = max(counts.items(), key=lambda kv: kv[1])
    return best_sep if best_count > 0 else ","


def jaccard(set_a: set, set_b: set) -> float:
    if not set_a and not set_b:
        return 1.0
    inter = len(set_a & set_b)
    union = len(set_a | set_b)
    return inter / union if union else 0.0


def run_one_wgcna(expr_for_wgcna: pd.DataFrame, meta: pd.DataFrame, power: int | None):
    """Run PyWGCNA on one (samples x genes) matrix; return module assignments (Series)."""
    import PyWGCNA as pw

    geneInfo = pd.DataFrame(index=expr_for_wgcna.columns)
    geneInfo["gene_id"] = geneInfo.index

    w = pw.WGCNA(name="Xf_WGCNA_LOO", geneExp=expr_for_wgcna,
                 sampleInfo=meta.loc[expr_for_wgcna.index], geneInfo=geneInfo, level=1)
    w.networkType = "signed"
    w.TOMType = "signed"
    w.minModuleSize = 30
    w.MEDissThres = 0.25

    if power is None:
        power, _ = w.pickSoftThreshold(data=expr_for_wgcna)
    w.power = power

    w.findModules(kwargs_function={"cutreeHybrid": {"deepSplit": 2, "pamRespectsDendro": False}})

    color_col = next((c for c in w.datExpr.var.columns
                       if "color" in c.lower() or "module" in c.lower()),
                      w.datExpr.var.columns[0])
    return w.datExpr.var[color_col].rename("module")


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Leave-one-out WGCNA module stability analysis.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--top-n", type=int, default=640)
    p.add_argument("--min-var", type=float, default=0.1)
    p.add_argument("--power", type=int, default=None,
                    help="Fix the soft-threshold power to match the main run "
                         "(recommended, so LOO variation isn't confounded by "
                         "power re-estimation each time).")
    p.add_argument("--reference-assignments", required=True, type=Path,
                    help="module_assignments.csv from the main (all-24-sample) run.")
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(args.tpm, index_col=0, sep=sniff_sep(args.tpm))
    meta = pd.read_csv(args.metadata, index_col=0, sep=sniff_sep(args.metadata))
    shared = expr.columns.intersection(meta.index)
    expr, meta = expr[shared], meta.loc[shared]
    samples = list(shared)
    print(f"[INFO] {len(samples)} samples total")

    ref = pd.read_csv(args.reference_assignments)
    ref_col = "module" if "module" in ref.columns else ref.columns[-1]
    gene_col = "gene" if "gene" in ref.columns else ref.columns[0]
    ref_sets = {m: set(g[gene_col]) for m, g in ref.groupby(ref_col)}
    print(f"[INFO] Reference modules: { {k: len(v) for k, v in ref_sets.items()} }")

    all_rows = []
    for held_out in samples:
        sub_samples = [s for s in samples if s != held_out]
        expr_sub = expr[sub_samples]

        keep = (expr_sub >= 1.0).sum(axis=1) >= 3
        expr_f = expr_sub.loc[keep]
        expr_log = np.log2(expr_f + 1)
        expr_log = expr_log[expr_log.var(axis=1) > args.min_var]
        top_genes = expr_log.var(axis=1).nlargest(args.top_n).index
        expr_for_wgcna = expr_log.loc[top_genes].T  # samples x genes

        try:
            assign = run_one_wgcna(expr_for_wgcna, meta, args.power)
        except SystemExit as e:
            print(f"[WARN] held_out={held_out}: WGCNA collapsed to one module ({e}); skipping.")
            continue

        sub_sets = {m: set(g.index) for m, g in
                    [(mod, assign[assign == mod]) for mod in assign.unique()]}

        for ref_mod, ref_genes in ref_sets.items():
            best_match, best_j = None, -1.0
            for sub_mod, sub_genes in sub_sets.items():
                j = jaccard(ref_genes, sub_genes)
                if j > best_j:
                    best_j, best_match = j, sub_mod
            all_rows.append({
                "held_out_sample": held_out,
                "reference_module": ref_mod,
                "n_ref_genes": len(ref_genes),
                "best_matching_submodule": best_match,
                "jaccard": best_j,
                "n_submodules_total": len(sub_sets),
            })

    result = pd.DataFrame(all_rows)
    result.to_csv(args.outdir / "loo_stability_results.csv", index=False)

    summary = result.groupby("reference_module")["jaccard"].agg(["mean", "std", "min", "max", "count"])
    summary.to_csv(args.outdir / "loo_stability_summary.csv")

    print(f"\n[INFO] Per-module stability (Jaccard overlap vs. reference, across {len(samples)} LOO runs):")
    print(summary.round(3).to_string())
    print(f"\n[INFO] Full results: {args.outdir / 'loo_stability_results.csv'}")


if __name__ == "__main__":
    main()
