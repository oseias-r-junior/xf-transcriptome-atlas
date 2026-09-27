"""
wgcna_sensitivity.py — Parameter sensitivity sweep for the WGCNA module
detection thresholds.

Checks whether the reported 3-module structure (dimgrey/darkgrey/silver)
is an artifact of the specific soft-power (beta) and deepSplit values
chosen for the main run, or is robust to reasonable alternative choices.

Sweep
-----
  power (beta):  main_power - 1, main_power, main_power + 1
  deepSplit:     1, 2, 3
  edge threshold (|r| for build_coexpression_network.py, not run here):
                 0.80, 0.85 (main), 0.90  -- see note below

For each (power, deepSplit) combination this script re-runs module
detection on the full 24-sample / top-640-gene set and reports:
  - number of modules detected
  - module sizes
  - Jaccard overlap of each new module against the best-matching reference
    module (dimgrey / darkgrey / silver), so you can see whether the same
    three biological modules re-emerge under each parameter combination.

The edge-threshold part of the sweep (0.80 vs 0.90 |r| cutoff) does NOT
require re-running WGCNA -- it only changes which edges are kept in
build_coexpression_network.py for the Cytoscape/Gephi export. Re-run that
script directly with --min-cor 0.80 and --min-cor 0.90 and compare the
resulting edge/node counts; no new script is needed for that part.

REQUIRES PyWGCNA — see the note in wgcna_stability.py for why this is not
reimplemented in pure NumPy.

IMPORTANT: --tpm must be an ortholog-merged matrix, not the strain-native
data/tpm_expression.csv -- see the long comment in 04_wgcna/run_wgcna.py's
own docstring for why, and data/README.md's "WGCNA input matrix" section
for which of data/tpm_expression_original.csv (the file that actually
produced the manuscript's WGCNA result -- use this one to reproduce/extend
it) vs data/tpm_expression_wgcna.csv (rebuilt from raw data, ~3% fewer
genes) to use.

Usage
-----
python wgcna_sensitivity.py \\
    --tpm        data/tpm_expression_original.csv \\
    --metadata   data/sample_info_original.csv \\
    --top-n      640 \\
    --main-power <beta used in the main run> \\
    --reference-assignments results/WGCNA/module_assignments.csv \\
    --outdir     results/WGCNA_sensitivity
"""
import argparse
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


def run_one_wgcna(expr_for_wgcna: pd.DataFrame, meta: pd.DataFrame, power: int, deep_split: int):
    import PyWGCNA as pw

    geneInfo = pd.DataFrame(index=expr_for_wgcna.columns)
    geneInfo["gene_id"] = geneInfo.index

    w = pw.WGCNA(name=f"Xf_WGCNA_p{power}_ds{deep_split}", geneExp=expr_for_wgcna,
                 sampleInfo=meta.loc[expr_for_wgcna.index], geneInfo=geneInfo, level=1)
    w.networkType = "signed"
    w.TOMType = "signed"
    w.minModuleSize = 30
    w.MEDissThres = 0.25
    w.power = power

    w.findModules(kwargs_function={"cutreeHybrid": {"deepSplit": deep_split, "pamRespectsDendro": False}})

    color_col = next((c for c in w.datExpr.var.columns
                       if "color" in c.lower() or "module" in c.lower()),
                      w.datExpr.var.columns[0])
    return w.datExpr.var[color_col].rename("module")


def main(argv=None):
    p = argparse.ArgumentParser(
        description="WGCNA parameter sensitivity sweep (power x deepSplit).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm", required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--top-n", type=int, default=640)
    p.add_argument("--min-var", type=float, default=0.1)
    p.add_argument("--main-power", required=True, type=int,
                    help="Soft-threshold power (beta) used in the main run; sweep tests beta-1, beta, beta+1.")
    p.add_argument("--deep-splits", nargs="*", type=int, default=[1, 2, 3])
    p.add_argument("--reference-assignments", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    args = p.parse_args(argv)

    args.outdir.mkdir(parents=True, exist_ok=True)

    expr = pd.read_csv(args.tpm, index_col=0, sep=sniff_sep(args.tpm))
    meta = pd.read_csv(args.metadata, index_col=0, sep=sniff_sep(args.metadata))
    shared = expr.columns.intersection(meta.index)
    expr, meta = expr[shared], meta.loc[shared]

    keep = (expr >= 1.0).sum(axis=1) >= 3
    expr_f = expr.loc[keep]
    expr_log = np.log2(expr_f + 1)
    expr_log = expr_log[expr_log.var(axis=1) > args.min_var]
    top_genes = expr_log.var(axis=1).nlargest(args.top_n).index
    expr_for_wgcna = expr_log.loc[top_genes].T
    print(f"[INFO] {expr_for_wgcna.shape[1]} genes x {expr_for_wgcna.shape[0]} samples for sweep")

    ref = pd.read_csv(args.reference_assignments)
    ref_col = "module" if "module" in ref.columns else ref.columns[-1]
    gene_col = "gene" if "gene" in ref.columns else ref.columns[0]
    ref_sets = {m: set(g[gene_col]) for m, g in ref.groupby(ref_col)}

    powers = sorted({args.main_power - 1, args.main_power, args.main_power + 1})
    powers = [p_ for p_ in powers if p_ >= 1]

    rows = []
    for power in powers:
        for ds in args.deep_splits:
            try:
                assign = run_one_wgcna(expr_for_wgcna, meta, power, ds)
            except SystemExit as e:
                print(f"[WARN] power={power}, deepSplit={ds}: collapsed to one module ({e}); skipping.")
                continue
            sub_sets = {m: set(assign[assign == m].index) for m in assign.unique()}
            for ref_mod, ref_genes in ref_sets.items():
                best_match, best_j = None, -1.0
                for sub_mod, sub_genes in sub_sets.items():
                    j = jaccard(ref_genes, sub_genes)
                    if j > best_j:
                        best_j, best_match = j, sub_mod
                rows.append({
                    "power": power, "deepSplit": ds,
                    "reference_module": ref_mod, "n_ref_genes": len(ref_genes),
                    "best_matching_submodule": best_match, "jaccard": best_j,
                    "n_modules_detected": len(sub_sets),
                    "is_main_params": (power == args.main_power and ds == 2),
                })

    result = pd.DataFrame(rows)
    result.to_csv(args.outdir / "sensitivity_results.csv", index=False)
    print("\n[INFO] Sensitivity sweep results:")
    print(result.round(3).to_string(index=False))
    print(f"\n[INFO] Full results: {args.outdir / 'sensitivity_results.csv'}")
    print("\n[NOTE] To check sensitivity to the network edge threshold "
          "(0.80 vs 0.90 |r|), re-run build_coexpression_network.py with "
          "--min-cor 0.80 and --min-cor 0.90 directly -- no WGCNA re-run "
          "needed for that part.")


if __name__ == "__main__":
    main()
