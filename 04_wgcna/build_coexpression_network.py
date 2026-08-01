"""
build_coexpression_network.py — Build and export a co-expression network from
WGCNA module assignments for visualization in Cytoscape or Gephi.

This script reads the outputs of run_wgcna.py and constructs gene-gene edges
based on their pairwise co-expression (Pearson correlation on the
log2(TPM+1) matrix), retaining only edges that exceed --min-cor and connect
genes within the same module.  An optional hub-gene subnetwork mode exports
only the top-kME genes per module.

When --gene-dict and --virulence-table are supplied, the node table is enriched
with ``pd_label`` (NCBI PD#### identifier or IMG ID fallback),
``is_virulence`` (boolean), and ``virulence_phase`` ("mobile" / "sessile").
Gene IDs in WGCNA output use the concatenated IMG format produced by
run_wgcna.py (e.g. "XF9a_01234|XFTem_05678"); the gene dictionary maps
XFTem_* → PD#### via Temecula1_old_locus_tags.

Usage
-----
python build_coexpression_network.py \\
    --expr            results/WGCNA/expr_log.csv \\
    --assignments     results/WGCNA/module_assignments.csv \\
    --kme             results/WGCNA/kME.csv \\
    --min-cor         0.80 \\
    --hub-n           20 \\
    --gene-dict       data/gene_dictionary.csv \\
    --virulence-table data/virulence_table.csv \\
    --outdir          results/network

Outputs
-------
results/network/
    edges.csv         gene_a, gene_b, correlation, module
    nodes.csv         gene_id, module, kME_to_module[, pd_label, is_virulence, virulence_phase]
    edges_hub.csv     edges restricted to hub genes (--hub-n per module)
    nodes_hub.csv     node table for hub subnetwork
"""

import argparse
import re
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


# ---------------------------------------------------------------------------
# Gene-dictionary helpers (optional node enrichment)
# ---------------------------------------------------------------------------

def _split_tags(raw: str) -> list[str]:
    """Split a tab/space/comma/semicolon-delimited tag string."""
    return [t.strip().strip('"') for t in re.split(r"[\t\s;,]+", str(raw)) if t.strip()]


def build_xftem_to_pd(gene_dict_path: Path) -> dict[str, str]:
    """Return {XFTem_IMG_ID: PD####} from the gene dictionary TSV/CSV."""
    dct = pd.read_csv(gene_dict_path, sep=None, engine="python")
    # Identify the Temecula1 IMG-ID column
    img_col = next(
        (c for c in dct.columns if "temecula" in c.lower() and "img" in c.lower()),
        None,
    )
    # Identify the old-locus-tag column holding PD#### identifiers
    old_col = next(
        (c for c in dct.columns if "temecula" in c.lower() and "old_locus" in c.lower()),
        None,
    )
    if img_col is None or old_col is None:
        return {}

    mapping: dict[str, str] = {}
    for _, row in dct.iterrows():
        img_id = str(row.get(img_col, "")).strip()
        if not img_id or img_id == "nan":
            continue
        for tag in _split_tags(row.get(old_col, "")):
            t = tag.upper()
            m = re.match(r"PD_?(\d+)", t)
            if m:
                mapping[img_id] = f"PD{m.group(1)}"
                break
    return mapping


def pd_label(gene_concat: str, xftem_to_pd: dict[str, str]) -> str:
    """Derive a human-readable label for a WGCNA gene (concatenated IMG ID)."""
    mT = re.search(r"(XFTem_\d+)", gene_concat)
    if mT:
        tag = xftem_to_pd.get(mT.group(1))
        if tag:
            return tag
    m9 = re.search(r"(XF9a_\d+)", gene_concat, re.IGNORECASE)
    return m9.group(1) if m9 else gene_concat


def load_virulence_phases(virulence_path: Path) -> dict[str, int]:
    """Return {PD####: phase_int} for phase 1 and 2 genes."""
    vir = pd.read_csv(virulence_path, sep=None, engine="python", encoding="latin-1")
    phase_col = next((c for c in vir.columns if "phase" in c.lower()), None)
    id_col    = next((c for c in vir.columns if "gene id" in c.lower()), None)
    if phase_col is None or id_col is None:
        return {}
    result: dict[str, int] = {}
    for _, row in vir.iterrows():
        try:
            ph = int(row[phase_col])
        except (ValueError, TypeError):
            continue
        if ph in (1, 2):
            result[str(row[id_col]).strip().upper()] = ph
    return result


def enrich_nodes(
    nodes: pd.DataFrame,
    gene_dict_path: Path,
    virulence_path: Path | None,
) -> pd.DataFrame:
    """Add pd_label, is_virulence, virulence_phase columns to the node table."""
    xftem_to_pd = build_xftem_to_pd(gene_dict_path)
    vir_phase: dict[str, int] = {}
    if virulence_path is not None:
        vir_phase = load_virulence_phases(virulence_path)

    phase_label = {1: "mobile", 2: "sessile"}

    labels, is_vir, phases = [], [], []
    for gid in nodes["gene_id"]:
        lbl = pd_label(str(gid), xftem_to_pd)
        labels.append(lbl)
        ph = vir_phase.get(lbl.upper()) if vir_phase else None
        is_vir.append(ph is not None)
        phases.append(phase_label.get(ph, "") if ph is not None else "")

    nodes = nodes.copy()
    nodes["pd_label"]       = labels
    nodes["is_virulence"]   = is_vir
    nodes["virulence_phase"] = phases
    return nodes


# ---------------------------------------------------------------------------
# Hub subnetwork
# ---------------------------------------------------------------------------

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
    p.add_argument(
        "--gene-dict", type=Path, default=None,
        help="gene_dictionary.csv (TSV) from build_gene_dictionary.py.  "
             "When supplied, adds pd_label, is_virulence, and virulence_phase "
             "columns to nodes.csv.",
    )
    p.add_argument(
        "--virulence-table", type=Path, default=None,
        help="Virulence gene table (CSV/TSV) with 'Gene ID' and 'Phase' columns.  "
             "Requires --gene-dict.  Adds virulence phase annotation to nodes.",
    )
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
    if args.gene_dict is not None:
        vir_path = args.virulence_table if args.virulence_table is not None else None
        nodes = enrich_nodes(nodes, args.gene_dict, vir_path)
        print(f"[INFO] Node table enriched with pd_label / virulence annotations")

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
