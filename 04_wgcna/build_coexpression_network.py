"""
build_coexpression_network.py — Build and export a co-expression network from
WGCNA module assignments for visualization in Cytoscape or Gephi.

This script reads the outputs of run_wgcna.py and constructs gene-gene edges
based on their pairwise co-expression (Pearson correlation on the
log2(TPM+1) matrix), retaining only edges that exceed --min-cor and connect
genes within the same module.  An optional hub-gene subnetwork mode exports
only the top-kME genes per module.

Usage
-----
python build_coexpression_network.py \\
    --expr          results/WGCNA/expr_log.csv \\
    --assignments   results/WGCNA/module_assignments.csv \\
    --kme           results/WGCNA/kME.csv \\
    --min-cor       0.80 \\
    --hub-n         20 \\
    --outdir        results/network

Outputs
-------
results/network/
    edges.csv         gene_a, gene_b, correlation, module
    nodes.csv         gene_id, module, kME_to_module
    edges_hub.csv     edges restricted to hub genes (--hub-n per module)
    nodes_hub.csv     node table for hub subnetwork
"""

import argparse
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd


def load_inputs(
    expr_path: Path,
    assignments_path: Path,
    kme_path: Path,
) -> tuple:
    expr = pd.read_csv(expr_path, index_col=0)         # genes × samples
    asgn = pd.read_csv(assignments_path)
    asgn = asgn.set_index(asgn.columns[0])["module"]  # gene → module
    kme  = pd.read_csv(kme_path, index_col=0)          # genes × modules
    return expr, asgn, kme


def kme_to_module(kme: pd.DataFrame, asgn: pd.Series) -> pd.Series:
    """For each gene, extract kME to its own module."""
    result = {}
    for gene, mod in asgn.items():
        col = f"ME{mod}" if f"ME{mod}" in kme.columns else mod
        if col in kme.columns and gene in kme.index:
            result[gene] = kme.loc[gene, col]
        else:
            result[gene] = np.nan
    return pd.Series(result, name="kME_to_module")


def build_edges(
    expr: pd.DataFrame,
    asgn: pd.Series,
    min_cor: float,
) -> pd.DataFrame:
    """Compute within-module pairwise correlations; keep edges ≥ min_cor."""
    rows = []
    for mod in asgn.unique():
        genes = asgn[asgn == mod].index.intersection(expr.index)
        if len(genes) < 2:
            continue
        X = expr.loc[genes].values
        corr = np.corrcoef(X)
        for i, j in combinations(range(len(genes)), 2):
            r = corr[i, j]
            if abs(r) >= min_cor:
                rows.append(
                    {
                        "gene_a":      genes[i],
                        "gene_b":      genes[j],
                        "correlation": round(r, 4),
                        "module":      mod,
                    }
                )
    return pd.DataFrame(rows)


def hub_subnetwork(
    edges: pd.DataFrame,
    asgn: pd.Series,
    kme: pd.DataFrame,
    hub_n: int,
) -> tuple:
    hub_genes: set = set()
    for mod in asgn.unique():
        genes = asgn[asgn == mod].index
        col   = f"ME{mod}" if f"ME{mod}" in kme.columns else mod
        if col in kme.columns:
            top = kme.loc[genes.intersection(kme.index), col].nlargest(hub_n).index
            hub_genes.update(top)
    edges_hub = edges[
        edges["gene_a"].isin(hub_genes) & edges["gene_b"].isin(hub_genes)
    ].copy()
    return edges_hub, hub_genes


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build co-expression network from WGCNA outputs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--expr",        required=True, type=Path)
    p.add_argument("--assignments", required=True, type=Path)
    p.add_argument("--kme",         required=True, type=Path)
    p.add_argument("--min-cor",     type=float, default=0.80,
                   help="Minimum |Pearson r| to include an edge.")
    p.add_argument("--hub-n",       type=int,   default=20,
                   help="Top-N hub genes per module for hub subnetwork.")
    p.add_argument("--outdir",      required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    expr, asgn, kme = load_inputs(args.expr, args.assignments, args.kme)
    print(f"[INFO] {expr.shape[0]} genes, {asgn.nunique()} modules")

    # Full network
    edges = build_edges(expr, asgn, args.min_cor)
    print(f"[INFO] {len(edges)} edges (|r| ≥ {args.min_cor})")
    edges.to_csv(args.outdir / "edges.csv", index=False)

    kme_series = kme_to_module(kme, asgn)
    nodes = pd.DataFrame({"gene_id": asgn.index, "module": asgn.values,
                          "kME_to_module": kme_series.reindex(asgn.index).values})
    nodes.to_csv(args.outdir / "nodes.csv", index=False)

    # Hub subnetwork
    edges_hub, hub_genes = hub_subnetwork(edges, asgn, kme, args.hub_n)
    edges_hub.to_csv(args.outdir / "edges_hub.csv", index=False)
    nodes_hub = nodes[nodes["gene_id"].isin(hub_genes)]
    nodes_hub.to_csv(args.outdir / "nodes_hub.csv", index=False)

    print(f"[INFO] Hub subnetwork: {len(hub_genes)} genes, {len(edges_hub)} edges")
    print(f"[INFO] Network files saved to {args.outdir}")


if __name__ == "__main__":
    main()
