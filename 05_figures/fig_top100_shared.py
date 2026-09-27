"""
fig_top100_shared.py — Mean TPM trend (early vs. late growth phase) of the
top-N ANOVA-selected shared genes between strains 9a5c and Temecula1
(Figure 3A).

This reproduces the exact analysis and plot from the original notebook cell
(PCA_PCoA_heatmap_Pierry_Transcriptome.ipynb, "Figure 3A pipeline (Option A)"
/ `process_data()` + `plot_expression_trends()`), whose output file was
`gene_expression_trends_top100_shared.tiff` — the figure actually submitted
to MDPI Pathogens (see figure_2_3_4_5_composites_mdpi_pathogens.pptx, slide 3,
panel A).

IMPORTANT — this is NOT a simple "top TPM" ranking. Gene selection
("Option A") is significance-first:
  1. Collapse 9a5c/Temecula1 orthologs to a single unified_gene_name (RBH
     pairs from the gene dictionary). `--tpm` is expected to already be
     indexed by this unified id (gene_id column = "<9a5c_IMG_ID>_<Temecula1_IMG_ID>"),
     matching data/tpm_expression.csv as produced by 01_preprocessing.
  2. For each gene, run a one-way ANOVA of TPM ~ (strain x medium x phase)
     group (8 groups, 3 replicates each = 24 samples) and BH-FDR-correct the
     p-values across all genes.
  3. Select genes as: all genes with FDR < alpha, sorted by mean TPM
     (descending); if fewer than --top-n, fill the remaining slots with the
     highest-mean-TPM non-significant genes until --top-n is reached.
  4. Aggregate mean +/- SE TPM per strain x medium x phase group across the
     selected gene set, then plot.

The ANOVA F-test is implemented in 06_stats/pure_stats.py (f_oneway_manual)
so this script has no scipy/statsmodels dependency, matching the
dependency-light convention already used by permanova_factors.py and
module_trait_heatmap.py in this repo.

Usage
-----
python fig_top100_shared.py \\
    --tpm      data/tpm_expression.csv \\
    --top-n    100 \\
    --alpha    0.05 \\
    --output   figures/fig_3A_top100_shared.tiff
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent.parent / "06_stats"))
from pure_stats import f_oneway_manual, benjamini_hochberg


# ---------------------------------------------------------------------------
# Sample -> strain / medium / phase parsing (column-name convention shared
# with fig_virulence_trends.py: "<strain prefix>_<medium>_<timepoint>_<rep>")
# ---------------------------------------------------------------------------

def default_strain(s: str) -> str:
    return "9a5c" if s.lower().startswith("9a") else "Temecula1"


def default_medium(s: str) -> str:
    return "PIM6" if "PIM6" in s else "PWG"


def default_phase(s: str) -> str:
    return "early" if ("PIM6_1d" in s or "PWG_3d" in s) else "late"


# Plot colours / jitter, exactly as in the original notebook cell
STRAIN_MEDIUM_COLORS = {
    ("9a5c", "PIM6"):      "#FF8C00",   # orange
    ("9a5c", "PWG"):       "#FF0000",   # red
    ("Temecula1", "PIM6"): "#ADD8E6",   # light blue
    ("Temecula1", "PWG"):  "#00008B",   # dark blue
}
MEDIUM_JITTER = {"PIM6": 0.1, "PWG": 0.1075}


# ---------------------------------------------------------------------------
# Data pipeline
# ---------------------------------------------------------------------------

def load_tpm(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, index_col=0, sep=sep)


def build_long_table(tpm: pd.DataFrame) -> pd.DataFrame:
    samples = tpm.columns.tolist()
    meta = pd.DataFrame({
        "sample": samples,
        "strain": [default_strain(s) for s in samples],
        "medium": [default_medium(s) for s in samples],
        "phase": [default_phase(s) for s in samples],
    })
    meta["strain_condition"] = meta["strain"] + " | " + meta["medium"] + " - " + meta["phase"]

    long_rows = tpm.reset_index().melt(id_vars=tpm.index.name or "index",
                                        var_name="sample", value_name="tpm")
    long_rows = long_rows.rename(columns={tpm.index.name or "index": "gene_id"})
    long_rows = long_rows.merge(meta, on="sample")
    return long_rows


def select_genes_option_a(long_df: pd.DataFrame, top_n: int, alpha: float) -> pd.DataFrame:
    """Reproduces process_data()'s Option-A gene selection: ANOVA + BH-FDR,
    significant genes first (by mean TPM desc), filled with highest-mean-TPM
    non-significant genes if fewer than top_n are significant."""
    stats_rows = []
    for gene, gdf in long_df.groupby("gene_id"):
        mean_tpm = gdf["tpm"].mean()
        groups = [g["tpm"].values for _, g in gdf.groupby("strain_condition")]
        if len(groups) > 1 and gdf["tpm"].nunique() > 1:
            f_stat, pval = f_oneway_manual(*groups)
        else:
            pval = np.nan
        stats_rows.append({"gene_id": gene, "mean_tpm": mean_tpm, "pval_anova": pval})

    stats_df = pd.DataFrame(stats_rows)

    n_genes_total = stats_df["gene_id"].nunique()
    top_n = min(top_n, n_genes_total)

    mask_valid = stats_df["pval_anova"].notna()
    stats_df["pval_fdr"] = np.nan
    if mask_valid.any():
        stats_df.loc[mask_valid, "pval_fdr"] = benjamini_hochberg(
            stats_df.loc[mask_valid, "pval_anova"].values
        )

    sig = stats_df[stats_df["pval_fdr"].notna() & (stats_df["pval_fdr"] < alpha)].copy()
    sig = sig.sort_values("mean_tpm", ascending=False)
    selected = list(sig["gene_id"])

    if len(selected) < top_n:
        non_sig = stats_df[~stats_df["gene_id"].isin(selected)].copy()
        non_sig = non_sig.sort_values("mean_tpm", ascending=False)
        for g in non_sig["gene_id"]:
            selected.append(g)
            if len(selected) == top_n:
                break
    else:
        selected = selected[:top_n]

    stats_df["selected"] = stats_df["gene_id"].isin(selected)
    print(f"[INFO] {n_genes_total} shared (unified) genes tested; "
          f"{int((stats_df['pval_fdr'] < alpha).sum())} FDR-significant; "
          f"{len(selected)} genes selected for the top-{top_n} set "
          f"({(stats_df.loc[stats_df['gene_id'].isin(selected), 'pval_fdr'] < alpha).sum()} "
          f"of those are FDR-significant).")
    return long_df[long_df["gene_id"].isin(selected)].copy()


def aggregate(top_df: pd.DataFrame) -> pd.DataFrame:
    agg = (
        top_df.groupby(["strain", "medium", "phase"])
        .agg(avg_tpm=("tpm", "mean"), std_tpm=("tpm", "std"), n=("tpm", "size"))
        .reset_index()
    )
    agg["std_error"] = agg["std_tpm"] / np.sqrt(agg["n"])
    return agg


# ---------------------------------------------------------------------------
# Plot — exact style of the original submitted figure
# ---------------------------------------------------------------------------

def plot_expression_trends(agg: pd.DataFrame, top_n: int, output: Path):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 8))

    for (st, md), col in STRAIN_MEDIUM_COLORS.items():
        sub = agg[(agg["strain"] == st) & (agg["medium"] == md)].sort_values("phase")
        if sub.empty:
            continue
        jitter = MEDIUM_JITTER[md]
        x = [0 + jitter if p == "early" else 1 + jitter for p in sub["phase"]]
        y = sub["avg_tpm"].values
        y_err = sub["std_error"].values

        ax.errorbar(x, y, yerr=y_err, fmt="o", color=col,
                    markersize=8, capsize=5, label=f"{st} - {md}")
        ax.plot(x, y, color=col, linestyle="-", alpha=0.7, linewidth=2)
        ax.fill_between(x, y - y_err, y + y_err, color=col, alpha=0.1)

    ax.set_xticks([0.1, 1.1])
    ax.set_xticklabels(["early", "late"])
    ax.set_xlabel("(exponential phase) growth")
    ax.set_ylabel(f"Mean TPM (top {top_n} shared genes)")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(fontsize=14)
    plt.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Figure 3A saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Figure 3A: mean TPM trend (early/late) of top-N ANOVA-selected shared genes.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm", required=True, type=Path,
                   help="TPM matrix indexed by unified gene_id (9a5c_id_Temecula1_id), "
                        "24 sample columns (data/tpm_expression.csv).")
    p.add_argument("--top-n", type=int, default=100)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm = load_tpm(args.tpm)
    long_df = build_long_table(tpm)

    top_df = select_genes_option_a(long_df, args.top_n, args.alpha)
    csv_out = args.output.with_suffix(".csv")
    agg = aggregate(top_df)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    agg.to_csv(csv_out, index=False)
    print(f"[INFO] Aggregated data saved to {csv_out}")
    print(agg.round(2).to_string(index=False))

    plot_expression_trends(agg, args.top_n, args.output)


if __name__ == "__main__":
    main()
