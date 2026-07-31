"""
run_wgcna.py — Weighted Gene Co-expression Network Analysis (WGCNA) using
PyWGCNA (Morabito et al. 2023; Rezaie et al. 2023).

Pipeline:
  1. Load TPM expression matrix and sample metadata
  2. Filter: genes with TPM ≥ 1 in ≥ 3 samples; variance > --min-var
  3. Log2(TPM+1) transform; select top --top-n most variable genes
  4. Estimate soft-thresholding power β (scale-free topology criterion R² ≥ 0.8)
  5. Detect co-expression modules (deepSplit=2, average linkage)
  6. Compute module eigengenes (MEs) and correlate with phenotypic traits
  7. Identify hub genes per module (top kME)
  8. Export: module assignments, kME table, ME matrix, trait correlations

Usage
-----
python run_wgcna.py \\
    --tpm       data/tpm_expression.csv \\
    --metadata  data/sample_info.csv \\
    --top-n     640 \\
    --min-var   0.1 \\
    --outdir    results/WGCNA

Input formats
-------------
TPM matrix (CSV, genes × samples):
    gene_id       9a5c_PIM6_1d_1   9a5c_PIM6_1d_2   ...
    Xf9a_00001    45.3             48.1             ...

Sample metadata (CSV, samples × traits):
    sample_id         strain   medium   timepoint   is_mobile   is_sessile
    9a5c_PIM6_1d_1    9a5c     PIM6     1d          1           0
    ...
Trait columns should be binary (0/1) for module-trait correlation.

Outputs
-------
results/WGCNA/
    module_assignments.csv   gene, module
    kME.csv                  kME per gene per module
    module_eigengenes.csv    ME per sample per module
    trait_correlations.csv   Pearson r and p-value for ME × trait
    hub_genes.csv            top-kME genes per module
    expr_log.csv             log2(TPM+1) matrix used for WGCNA
    sample_info_used.csv     metadata aligned to expression matrix
    soft_threshold_plot.tiff scale-free fit vs power
    module_trait_heatmap.tiff
"""

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Expression loading and preprocessing
# ---------------------------------------------------------------------------

def load_expression(path: Path) -> pd.DataFrame:
    """Load TPM matrix (genes × samples)."""
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, index_col=0, sep=sep)


def preprocess(
    expr: pd.DataFrame,
    top_n: int,
    min_var: float,
    min_tpm: float = 1.0,
    min_samples: int = 3,
) -> pd.DataFrame:
    """Filter and log-transform; return log2(TPM+1) for top_n variable genes."""
    # Filter low-expressed
    keep = (expr >= min_tpm).sum(axis=1) >= min_samples
    expr_f = expr.loc[keep]
    print(f"[INFO] After expression filter: {expr_f.shape[0]} genes")

    # Log transform
    expr_log = np.log2(expr_f + 1)

    # Variance filter
    expr_log = expr_log[expr_log.var(axis=1) > min_var]
    print(f"[INFO] After variance filter (>{min_var}): {expr_log.shape[0]} genes")

    # Top-N most variable
    top_genes = expr_log.var(axis=1).nlargest(top_n).index
    expr_sel  = expr_log.loc[top_genes]
    print(f"[INFO] Selecting top {len(expr_sel)} variable genes for WGCNA")
    return expr_sel


# ---------------------------------------------------------------------------
# Module eigengenes via PCA (robust fallback)
# ---------------------------------------------------------------------------

def compute_module_eigengenes(
    expr_log: pd.DataFrame,   # genes × samples
    assignments: pd.Series,   # gene → module
) -> pd.DataFrame:
    """Compute module eigengene as first PC of intra-module expression."""
    mes: dict[str, np.ndarray] = {}
    for mod in assignments.unique():
        genes = assignments[assignments == mod].index
        sub = expr_log.loc[genes].T.values   # samples × genes
        if sub.shape[1] < 2:
            mes[mod] = sub[:, 0]
            continue
        pca = PCA(n_components=1)
        me  = pca.fit_transform(sub).squeeze()
        # Convention: ME correlated with mean expression
        mean_exp = sub.mean(axis=1)
        if np.corrcoef(me, mean_exp)[0, 1] < 0:
            me = -me
        mes[mod] = me

    return pd.DataFrame(mes, index=expr_log.columns)


# ---------------------------------------------------------------------------
# Trait correlation
# ---------------------------------------------------------------------------

def correlate_traits(
    MEs: pd.DataFrame,        # samples × modules
    meta: pd.DataFrame,       # samples × traits
    trait_cols: list[str],
) -> pd.DataFrame:
    """Pearson r between each ME and each trait, with BH FDR correction."""
    rows = []
    for mod in MEs.columns:
        for trait in trait_cols:
            shared = MEs.index.intersection(meta.index)
            me_vals  = MEs.loc[shared, mod].values
            tr_vals  = meta.loc[shared, trait].values.astype(float)
            mask = ~np.isnan(me_vals) & ~np.isnan(tr_vals)
            if mask.sum() < 4:
                continue
            r, p = pearsonr(me_vals[mask], tr_vals[mask])
            rows.append({"module": mod, "trait": trait, "r": r, "p_value": p})

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    _, df["p_adj"], _, _ = multipletests(df["p_value"], method="fdr_bh")
    return df


def plot_trait_heatmap(df_corr: pd.DataFrame, output: Path):
    pivot_r = df_corr.pivot(index="module", columns="trait", values="r")
    pivot_p = df_corr.pivot(index="module", columns="trait", values="p_adj")

    annot = pivot_r.applymap(lambda v: f"{v:.2f}") + "\n" + \
            pivot_p.applymap(lambda p: ("*" if p < 0.05 else ""))

    fig, ax = plt.subplots(figsize=(max(6, len(pivot_r.columns) * 0.8),
                                    max(4, len(pivot_r) * 0.5)))
    sns.heatmap(
        pivot_r,
        annot=annot, fmt="",
        cmap="RdBu_r", center=0, vmin=-1, vmax=1,
        linewidths=0.5, ax=ax,
    )
    ax.set_title("Module–trait correlation (r; * p_adj < 0.05)")
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Trait heatmap saved to {output}")


# ---------------------------------------------------------------------------
# PyWGCNA wrapper
# ---------------------------------------------------------------------------

def run_pywgcna(
    expr_for_wgcna: pd.DataFrame,   # samples × genes
    meta: pd.DataFrame,
    outdir: Path,
    power: int | None = None,
) -> tuple:
    try:
        import PyWGCNA as pw
    except ImportError:
        sys.exit("[ERROR] PyWGCNA is not installed. Run: pip install PyWGCNA")

    geneInfo = pd.DataFrame(index=expr_for_wgcna.columns)
    geneInfo["gene_id"] = geneInfo.index

    w = pw.WGCNA(
        name="Xf_WGCNA",
        geneExp=expr_for_wgcna,
        sampleInfo=meta.loc[expr_for_wgcna.index],
        geneInfo=geneInfo,
        level=1,
    )
    w.networkType  = "signed"
    w.TOMType      = "signed"
    w.minModuleSize = 30
    w.MEDissThres  = 0.25

    # Soft-threshold selection
    estimated_power, sft_df = w.pickSoftThreshold(data=expr_for_wgcna)
    w.power = power or estimated_power
    print(f"[INFO] Using soft-threshold power β = {w.power}")

    # Plot scale-free fit
    fig, ax = plt.subplots()
    ax.plot(sft_df["Power"], sft_df["SFT.R.sq"], "o-")
    ax.axhline(0.8, color="red", linestyle="--", label="R²=0.8 threshold")
    ax.axvline(w.power, color="orange", linestyle="--", label=f"β={w.power}")
    ax.set_xlabel("Soft-threshold power (β)")
    ax.set_ylabel("Scale-free fit R²")
    ax.legend()
    sft_path = outdir / "soft_threshold_plot.tiff"
    fig.savefig(sft_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Module detection
    w.findModules(kwargs_function={
        "cutreeHybrid": {"deepSplit": 2, "pamRespectsDendro": False}
    })

    # Extract module assignments
    color_col = next(
        (c for c in w.datExpr.var.columns
         if "color" in c.lower() or "module" in c.lower()),
        w.datExpr.var.columns[0],
    )
    assignments = w.datExpr.var[color_col].rename("module")
    return w, assignments


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run WGCNA on TPM expression data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm",      required=True, type=Path)
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--top-n",    type=int,   default=640)
    p.add_argument("--min-var",  type=float, default=0.1)
    p.add_argument("--power",    type=int,   default=None,
                   help="Force soft-threshold power (default: auto-detect).")
    p.add_argument("--trait-cols", nargs="*", default=None,
                   help="Metadata columns for trait correlation (default: all numeric).")
    p.add_argument("--outdir",   required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # -- Load ----------------------------------------------------------------
    expr = load_expression(args.tpm)
    meta_sep = "\t" if args.metadata.suffix in {".tsv", ".txt"} else ","
    meta = pd.read_csv(args.metadata, index_col=0, sep=meta_sep)

    # Align samples
    shared = expr.columns.intersection(meta.index)
    expr = expr[shared]
    meta = meta.loc[shared]
    print(f"[INFO] {expr.shape[0]} genes × {len(shared)} samples (after alignment)")

    # -- Preprocess ----------------------------------------------------------
    expr_sel = preprocess(expr, args.top_n, args.min_var)
    expr_for_wgcna = expr_sel.T   # samples × genes

    # Save preprocessed matrix
    expr_sel.to_csv(args.outdir / "expr_log.csv")

    # -- WGCNA ---------------------------------------------------------------
    w, assignments = run_pywgcna(expr_for_wgcna, meta, args.outdir, args.power)

    # -- MEs and hub genes ---------------------------------------------------
    MEs = compute_module_eigengenes(expr_sel, assignments)
    MEs.to_csv(args.outdir / "module_eigengenes.csv")

    # kME (correlation of each gene with each ME)
    kME = pd.DataFrame(
        np.corrcoef(expr_sel.values, MEs.T.values)[: len(expr_sel), len(expr_sel):],
        index=expr_sel.index,
        columns=MEs.columns,
    )
    kME.to_csv(args.outdir / "kME.csv")

    # Hub genes: top-10 per module by kME
    hub_rows = []
    for mod in assignments.unique():
        genes = assignments[assignments == mod].index
        if f"ME{mod}" not in kME.columns:
            continue
        top = kME.loc[genes, f"ME{mod}"].nlargest(10)
        for g, k in top.items():
            hub_rows.append({"gene": g, "module": mod, "kME": k})
    pd.DataFrame(hub_rows).to_csv(args.outdir / "hub_genes.csv", index=False)

    # Module assignments
    assignments.reset_index().rename(columns={"index": "gene"}) \
        .to_csv(args.outdir / "module_assignments.csv", index=False)

    # -- Trait correlation ---------------------------------------------------
    trait_cols = args.trait_cols or [
        c for c in meta.columns if pd.api.types.is_numeric_dtype(meta[c])
    ]
    if trait_cols:
        corr_df = correlate_traits(MEs, meta, trait_cols)
        corr_df.to_csv(args.outdir / "trait_correlations.csv", index=False)
        if not corr_df.empty:
            plot_trait_heatmap(corr_df, args.outdir / "module_trait_heatmap.tiff")

    meta.to_csv(args.outdir / "sample_info_used.csv")
    print(f"[INFO] WGCNA outputs saved to {args.outdir}")


if __name__ == "__main__":
    main()
