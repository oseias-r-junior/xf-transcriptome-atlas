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

IMPORTANT -- which TPM matrix to pass in
-----------------------------------------
--tpm must be an ORTHOLOG-MERGED matrix, gene_id = "<9a5c_IMG_ID>_<Temecula1_IMG_ID>"
(e.g. "XF9a_00002_XFTem_00002"), NOT the strain-native data/tpm_expression.csv
produced directly by compute_tpm.py. That strain-native file is genes x
samples where "genes" is the union of each strain's own native gene IDs:
every row is only ever populated for ONE strain's 12 columns (the other
strain's 12 columns are hard structural zeros, not biology). Feeding it
straight into WGCNA across all 24 samples confounds strain identity with
co-expression.

Two ortholog-merged files exist in data/ -- see data/README.md's "WGCNA
input matrix" section for the full comparison, but in short:

  data/tpm_expression_original.csv  (1692 genes / 1685 after filtering)
      The exact file from WGCNA_paper_pierry_feitosa.ipynb. THIS is what
      produced the published WGCNA result (3 modules: dimgrey/darkgrey/
      silver) -- use this one to reproduce or extend that result (module
      stability, power/deepSplit sensitivity, module-trait correlation).
      Treat it as a frozen, canonical input; don't try to regenerate it.

  data/tpm_expression_wgcna.csv  (1639 genes)
      Built from raw data by 01_preprocessing/build_wgcna_tpm_matrix.py.
      Fully reproducible from the raw CLC .xlsx exports through this repo's
      own scripts, but differs from the original by 53 gene pairs (~3.1%)
      excluded for a legitimate reason (ambiguous protein_id -- see that
      script's docstring). Confirmed (2026-09) that this 53-gene gap is the
      full explanation for a module-count difference between the two files
      (same PyWGCNA version, same power/cutHeight auto-detection either
      way) -- not a code or environment bug. Use this one only to
      demonstrate/validate end-to-end reproducibility from raw data.

    python 04_wgcna/run_wgcna.py \\
        --tpm       data/tpm_expression_original.csv \\
        --metadata  data/sample_info_original.csv \\
        --top-n     640 \\
        --outdir    results/WGCNA
    # (omit --power/--cut-height to auto-detect, as the original analysis
    # effectively did -- see run_pywgcna()'s own comment on why)

Input formats
-------------
TPM matrix (CSV, genes × samples, ortholog-pair-keyed -- see above):
    gene_id                        9a5c_PIM6_1d_rep1   ...  Tem1_PWG_7d_rep3
    XF9a_00002_XFTem_00002         45.3                ...  12.1

Sample metadata (CSV, samples × traits) -- two schemas are both accepted,
in either the old or new sample-naming convention (see normalize_sample_name /
load_metadata below, which reconcile whichever combination is passed in):

  new/rebuilt schema (data/sample_info.csv):
    sample_id           strain   medium   timepoint   is_mobile   is_sessile
    9a5c_PIM6_1d_rep1    9a5c     PIM6     1d          1           0
    ...

  original notebook schema (data/sample_info_original.csv, from
  WGCNA_paper_pierry_feitosa.ipynb):
    sampleID        strain   medium   timepoint   replicate
    9a_PIM6_1d_1     9a5c     PIM6     early       1
    ...

Trait handling for module-trait correlation does NOT rely on binary
is_mobile/is_sessile/is_early/is_late columns being present (those names
mean something different elsewhere in this repo -- a virulence GENE's
biofilm phase, not a sample condition). Instead build_trait_matrix() below
one-hot-encodes every categorical metadata column (strain, medium,
timepoint, ...) and passes numeric columns (e.g. "replicate") through
as-is, exactly matching WGCNA_paper_pierry_feitosa.ipynb's own Cell 9.

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
import re
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

def sniff_sep(path: Path) -> str:
    """Detect the field separator from the file's own header line.

    File extension alone is NOT reliable here: compute_tpm.py (and some
    manually-exported files) write tab-separated data into a file named
    '*.csv'. Reading such a file with sep="," (as a naive extension-based
    check would do) silently collapses every row into a single column,
    which downstream shows up as "N genes x 0 samples (after alignment)"
    -- an easy trap to fall into and hard to diagnose from that message
    alone. Sniff the header line instead: whichever of tab/comma/semicolon
    appears more often wins; default to comma if neither appears.
    """
    with open(path, "r", newline="") as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ("\t", ",", ";")}
    best_sep, best_count = max(counts.items(), key=lambda kv: kv[1])
    return best_sep if best_count > 0 else ","


def load_expression(path: Path) -> pd.DataFrame:
    """Load TPM matrix (genes × samples)."""
    sep = sniff_sep(path)
    return pd.read_csv(path, index_col=0, sep=sep)


# ---------------------------------------------------------------------------
# Sample-name normalisation
# ---------------------------------------------------------------------------

# Two sample-naming conventions show up across this repo's history:
#   - old/original (e.g. WGCNA_paper_pierry_feitosa.ipynb, its tpm_expression.csv
#     and sample_info.csv): "9a_PIM6_1d_1", "Tem1_PWG_7d_3"
#   - new/rebuilt (build_matrices_from_fpkm_xlsx.py's own convention, used by
#     data/tpm_expression.csv and data/sample_info.csv): "9a5c_PIM6_1d_rep1",
#     "Tem1_PWG_7d_rep3"
# The only differences are the 9a5c strain prefix ("9a" vs "9a5c") and the
# replicate suffix ("_1" vs "_rep1") -- Temecula1's own prefix ("Tem1") is
# already the same in both. Whichever --tpm/--metadata pair is passed in,
# normalise both sides to the SAME canonical form before aligning, so a
# mismatch in convention (e.g. an old-style sample_info.csv against the
# new-style tpm_expression_wgcna.csv) doesn't silently produce a 0-sample
# alignment (or worse, a partial one) the way a plain string match would.
_SAMPLE_RE = re.compile(
    r"^(?P<strain>9a5c|9a|Tem1|Temecula1)_(?P<medium>[A-Za-z0-9]+)_"
    r"(?P<timepoint>\d+[dh])_(?:rep)?(?P<rep>\d+)$"
)
_STRAIN_CANON = {"9a5c": "9a5c", "9a": "9a5c", "Tem1": "Tem1", "Temecula1": "Tem1"}


def normalize_sample_name(name: str) -> str:
    """Map either naming convention to one canonical form.

    '9a_PIM6_1d_1' and '9a5c_PIM6_1d_rep1' both -> '9a5c_PIM6_1d_rep1'.
    Names that don't match the expected pattern are returned unchanged
    (stripped of whitespace) so this degrades gracefully instead of hiding
    a genuinely different naming scheme.
    """
    name = str(name).strip()
    m = _SAMPLE_RE.match(name)
    if not m:
        return name
    strain = _STRAIN_CANON.get(m.group("strain"), m.group("strain"))
    return f"{strain}_{m.group('medium')}_{m.group('timepoint')}_rep{int(m.group('rep'))}"


def load_metadata(path: Path) -> pd.DataFrame:
    """Load sample metadata, tolerant of both naming/schema conventions.

    Handles:
      - either naming convention in the sample index (normalised via
        normalize_sample_name so it lines up with load_expression()'s
        columns regardless of which convention --tpm used)
      - stray leading/trailing whitespace in column names (the original
        notebook's own sample_info.csv has literal ", 9a5c, PIM6, early, 1"
        rows -- a plain read_csv leaves every string value with a leading
        space, e.g. " 9a5c" != "9a5c"). Column *values* are intentionally
        NOT stripped here: build_trait_matrix() one-hot-encodes them with
        pd.get_dummies() exactly as the original notebook's own Cell 9 does
        (WGCNA_paper_pierry_feitosa.ipynb), and that notebook never strips
        the values either -- so the resulting dummy column names carry the
        same leading space (e.g. "strain_ 9a5c") there as here. Stripping
        would produce prettier column names but no longer match a
        byte-for-byte diff against the original notebook's own output.

    Do NOT invent trait columns here (no is_mobile/is_sessile/is_early/
    is_late derivation). Those names are used elsewhere in this repo for a
    completely different thing -- a virulence GENE's biofilm phase
    (mobile/sessile), not a sample condition -- and inventing a
    medium-based analogue for samples was wrong. Trait handling for
    module-trait correlation happens in build_trait_matrix() instead,
    which mirrors the original notebook's own dummy-encoding logic exactly.
    """
    meta = pd.read_csv(path, index_col=0, sep=sniff_sep(path))
    meta.columns = meta.columns.str.strip()
    meta.index = meta.index.map(normalize_sample_name)
    return meta


def build_trait_matrix(meta: pd.DataFrame) -> pd.DataFrame:
    """Reproduce WGCNA_paper_pierry_feitosa.ipynb's own trait-matrix logic
    (Cell 9 / "WGCNA module-trait heatmap (REFATORADA + FIXED)") exactly:

        traits = meta.copy()
        for c in traits.columns:
            if traits[c].dtype == object or str(traits[c].dtype).startswith("category"):
                d = pd.get_dummies(traits[c], prefix=c)
                traits = traits.drop(columns=c).join(d)

    i.e. every categorical/string column (strain, medium, timepoint, ...) is
    one-hot encoded into one 0/1 column per observed category; any column
    that is already numeric (e.g. "replicate", or a pre-existing binary
    trait column) is kept as-is, untouched. No column is excluded -- the
    original notebook correlates "replicate" against module eigengenes too
    (weak/non-significant r values in practice, but it's part of the
    original output and this reproduces it faithfully).
    """
    traits = meta.copy()
    for c in traits.columns:
        if traits[c].dtype == object or str(traits[c].dtype).startswith("category"):
            d = pd.get_dummies(traits[c], prefix=c)
            traits = traits.drop(columns=c).join(d)
    return traits


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
    cut_height: float | None = None,
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

    # Module detection.
    # cutreeHybrid's own "cutHeight" (dynamicTreeCut) defaults to None, which
    # PyWGCNA/dynamicTreeCut resolves at runtime to 99% of the (truncated)
    # dendrogram height range -- that's the ".. cutHeight not given, setting
    # it to X ===> 99% of the (truncated) height range in dendro." line you
    # see logged, and X varies run to run with the dendrogram itself (gene
    # set, power, etc.). It is NOT a parameter PyWGCNA lets you read back out
    # after the fact and is not stored on `w` -- the only way to pin it to a
    # specific value (e.g. to match a past run's own auto-picked 0.436) is to
    # pass that value in explicitly here, via --cut-height.
    cutree_kwargs = {"deepSplit": 2, "pamRespectsDendro": False}
    if cut_height is not None:
        cutree_kwargs["cutHeight"] = cut_height
        print(f"[INFO] Forcing cutreeHybrid cutHeight = {cut_height}")
    w.findModules(kwargs_function={"cutreeHybrid": cutree_kwargs})

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
    p.add_argument("--cut-height", type=float, default=None,
                   help="Force cutreeHybrid's cutHeight (dynamicTreeCut dendrogram "
                        "cut height). Default: PyWGCNA auto-picks 99%% of the "
                        "truncated dendrogram height range at runtime, which "
                        "varies with the gene set/power -- pass this to pin it "
                        "to a specific past run's value (e.g. 0.436).")
    p.add_argument("--trait-cols", nargs="*", default=None,
                   help="Column names in the DUMMY-ENCODED trait matrix "
                        "(build_trait_matrix output, e.g. 'strain_9a5c', "
                        "'medium_PIM6') to correlate against module "
                        "eigengenes. Default: every numeric/dummy column.")
    p.add_argument("--outdir",   required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)

    # -- Load ----------------------------------------------------------------
    expr = load_expression(args.tpm)
    expr.columns = expr.columns.map(normalize_sample_name)
    meta = load_metadata(args.metadata)

    # Align samples (both sides already normalised to the same naming
    # convention by normalize_sample_name, so this works regardless of
    # whether --tpm/--metadata individually use the old or new convention)
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
    w, assignments = run_pywgcna(
        expr_for_wgcna, meta, args.outdir, args.power, args.cut_height
    )

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
    # kME's columns come straight from MEs.columns (compute_module_eigengenes()
    # keys its dict by the bare module name), so there is no "ME" prefix to
    # match here -- f"ME{mod}" never matched, silently skipping every module
    # and leaving hub_genes.csv empty.
    hub_rows = []
    for mod in assignments.unique():
        genes = assignments[assignments == mod].index
        if mod not in kME.columns:
            continue
        top = kME.loc[genes, mod].nlargest(10)
        for g, k in top.items():
            hub_rows.append({"gene": g, "module": mod, "kME": k})
    pd.DataFrame(hub_rows).to_csv(args.outdir / "hub_genes.csv", index=False)

    # Module assignments
    assignments.reset_index().rename(columns={"index": "gene"}) \
        .to_csv(args.outdir / "module_assignments.csv", index=False)

    # -- Trait correlation ---------------------------------------------------
    # Dummy-encode categorical columns (strain, medium, timepoint, ...) into
    # one 0/1 column per category, exactly as the original notebook's own
    # Cell 9 does (see build_trait_matrix's docstring) -- rather than only
    # picking out columns that already happen to be numeric/binary.
    traits = build_trait_matrix(meta)
    trait_cols = args.trait_cols or [
        c for c in traits.columns if pd.api.types.is_numeric_dtype(traits[c])
    ]
    if trait_cols:
        corr_df = correlate_traits(MEs, traits, trait_cols)
        corr_df.to_csv(args.outdir / "trait_correlations.csv", index=False)
        if not corr_df.empty:
            plot_trait_heatmap(corr_df, args.outdir / "module_trait_heatmap.tiff")

    meta.to_csv(args.outdir / "sample_info_used.csv")
    print(f"[INFO] WGCNA outputs saved to {args.outdir}")


if __name__ == "__main__":
    main()
