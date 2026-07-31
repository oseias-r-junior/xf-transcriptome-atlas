"""
fig_top100_shared.py — Mean transcript abundance of the top 100 shared genes
between strains 9a5c and Temecula1 (Figure 3A).

"Shared genes" are reciprocal best-hit (RBH) ortholog pairs identified by
build_gene_dictionary.py. For each pair the script uses the mean TPM of the
9a5c gene (or optionally both strains combined) averaged across all conditions.
The top 100 by mean TPM are plotted as a horizontal bar chart coloured by
functional category (if a --virulence-table with a 'function' column is
provided; otherwise coloured by strain).

Usage
-----
python fig_top100_shared.py \\
    --tpm          data/tpm_expression.csv \\
    --metadata     data/sample_info.csv \\
    --dictionary   data/gene_dictionary.tsv \\
    --top-n        100 \\
    --output       figures/fig_3A_top100_shared.tiff

Optional:
    --virulence-table data/virulence_table.csv   # adds functional colour coding
    --strain-col      strain                     # metadata column (default: strain)
    --gene-labels     9a5c_old_locus_tags        # dictionary column for gene labels
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd


FUNCTION_COLORS = {
    "adhesion":        "#e6194B",
    "regulation":      "#3cb44b",
    "motility":        "#4363d8",
    "stress response": "#f58231",
    "secretion":       "#911eb4",
    "metabolism":      "#42d4f4",
    "membrane":        "#f032e6",
    "others":          "#cccccc",
}
STRAIN_COLORS = {"9a5c": "#3399FF", "temecula1": "#FF6666"}
DEFAULT_COLOR = "#888888"


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------

def load_tpm(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    return pd.read_csv(path, index_col=0, sep=sep)


def condition_mean_tpm(
    tpm: pd.DataFrame,
    meta: pd.DataFrame,
    condition_col: str,
    strain_col: str,
    strain: str | None = None,
) -> pd.Series:
    """Mean TPM per gene across all conditions (optionally filtered to one strain)."""
    shared = tpm.columns.intersection(meta.index)
    tpm_s, meta_s = tpm[shared], meta.loc[shared]
    if strain:
        sids = meta_s[meta_s[strain_col].str.lower() == strain.lower()].index
        tpm_s = tpm_s[sids]
    return tpm_s.mean(axis=1)


def load_dictionary(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def shared_gene_mean_tpm(
    tpm: pd.DataFrame,
    meta: pd.DataFrame,
    dictionary: pd.DataFrame,
    strain_col: str,
    gene_label_col: str,
) -> pd.DataFrame:
    """
    For each RBH pair compute the mean TPM of the 9a5c representative gene
    across all 9a5c conditions.  Returns a DataFrame with columns:
      9a5c_id, label, mean_tpm
    """
    mean_9a = condition_mean_tpm(tpm, meta, "condition", strain_col, strain="9a5c")

    id_col_9a = "9a5c_IMG_ID"
    if id_col_9a not in dictionary.columns:
        raise ValueError(f"[ERROR] Column '{id_col_9a}' not found in dictionary.")

    rows = []
    for _, row in dictionary.iterrows():
        gid = str(row[id_col_9a]).strip()
        if gid not in mean_9a.index:
            continue
        # Preferred display label: old locus tag
        if gene_label_col and gene_label_col in dictionary.columns:
            raw_tag = str(row.get(gene_label_col, ""))
            label = raw_tag.split()[0] if raw_tag and raw_tag != "nan" else gid
        else:
            label = gid
        rows.append({"gene_id": gid, "label": label, "mean_tpm": mean_9a[gid]})

    return pd.DataFrame(rows).drop_duplicates("gene_id").sort_values(
        "mean_tpm", ascending=False
    )


# ---------------------------------------------------------------------------
# Colour mapping
# ---------------------------------------------------------------------------

def add_function_color(
    df: pd.DataFrame,
    vir_path: Path | None,
) -> pd.DataFrame:
    df = df.copy()
    if vir_path is None:
        df["color"] = DEFAULT_COLOR
        return df

    sep = "\t" if vir_path.suffix in {".tsv", ".txt"} else ","
    vir = pd.read_csv(vir_path, sep=sep)
    id_col = next(c for c in vir.columns if "gene_id" in c.lower())
    func_col = next((c for c in vir.columns if "function" in c.lower()), None)
    if func_col is None:
        df["color"] = DEFAULT_COLOR
        return df

    func_map = dict(zip(vir[id_col], vir[func_col].str.lower()))
    df["color"] = df["gene_id"].map(
        lambda g: FUNCTION_COLORS.get(func_map.get(g, "others"), "#cccccc")
    )
    return df


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

def plot_bar(df: pd.DataFrame, top_n: int, output: Path, use_function_colors: bool):
    top = df.head(top_n).copy()
    top = top.iloc[::-1]   # highest TPM at top

    fig_h = max(8, 0.22 * top_n)
    fig, ax = plt.subplots(figsize=(7, fig_h))

    ax.barh(
        range(len(top)),
        top["mean_tpm"],
        color=top["color"] if "color" in top.columns else DEFAULT_COLOR,
        edgecolor="white",
        linewidth=0.3,
        height=0.7,
    )
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["label"], fontsize=6.5)
    ax.set_xlabel("Mean TPM (9a5c, all conditions)", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.grid(True, linestyle="--", alpha=0.4, zorder=0)

    if use_function_colors:
        handles = [
            mpatches.Patch(facecolor=v, label=k)
            for k, v in FUNCTION_COLORS.items()
            if v in top["color"].values
        ]
        if handles:
            ax.legend(handles=handles, title="Function",
                      fontsize=7, title_fontsize=7,
                      loc="lower right", frameon=True)

    ax.set_title(f"Top {top_n} shared genes (9a5c × Temecula1 RBH pairs)", fontsize=10)
    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Top-{top_n} shared genes plot saved to {output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Mean TPM bar chart of top shared genes between strains.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm",             required=True, type=Path)
    p.add_argument("--metadata",        required=True, type=Path)
    p.add_argument("--dictionary",      required=True, type=Path,
                   help="Gene dictionary (build_gene_dictionary.py output).")
    p.add_argument("--top-n",           type=int, default=100)
    p.add_argument("--strain-col",      default="strain")
    p.add_argument("--condition-col",   default="condition")
    p.add_argument("--gene-labels",     default="9a5c_old_locus_tags",
                   help="Dictionary column to use as gene labels in the plot.")
    p.add_argument("--virulence-table", type=Path, default=None,
                   help="Optional: adds function-based colour coding.")
    p.add_argument("--output",          required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    sep_m = "\t" if args.metadata.suffix in {".tsv", ".txt"} else ","
    meta  = pd.read_csv(args.metadata, index_col=0, sep=sep_m)
    tpm   = load_tpm(args.tpm)
    dic   = load_dictionary(args.dictionary)

    print(f"[INFO] {len(dic)} RBH pairs in dictionary")

    df = shared_gene_mean_tpm(tpm, meta, dic, args.strain_col, args.gene_labels)
    print(f"[INFO] {len(df)} shared genes with TPM data; selecting top {args.top_n}")

    df = add_function_color(df, args.virulence_table)
    plot_bar(df, args.top_n, args.output,
             use_function_colors=args.virulence_table is not None)


if __name__ == "__main__":
    main()
