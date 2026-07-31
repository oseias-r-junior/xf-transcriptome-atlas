"""
fig_pcoa.py — Principal Coordinates Analysis (PCoA) of TPM expression profiles
with PERMANOVA and PERMDISP significance tests (Figure 2A).

Distance metric: Bray-Curtis dissimilarity on log2(TPM+1)-transformed data.

Uses scikit-bio for PCoA + PERMANOVA + PERMDISP.

Usage
-----
python fig_pcoa.py \\
    --tpm        data/tpm_expression.csv \\
    --metadata   data/sample_info.csv \\
    --group-col  condition \\
    --n-perms    999 \\
    --output     figures/fig_pcoa.tiff

Arguments
---------
--group-col   Column in metadata used to colour/group samples (default: condition).
--shape-col   Optional second metadata column for point shapes (e.g., strain).
--n-perms     Number of PERMANOVA/PERMDISP permutations (default: 999).
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

# scikit-bio imports
try:
    from skbio import DistanceMatrix
    from skbio.stats.ordination import pcoa
    from skbio.stats.distance import permanova, permdisp
except ImportError:
    raise SystemExit(
        "[ERROR] scikit-bio is required. Install with:\n"
        "  pip install scikit-bio"
    )


PALETTE = [
    "#e6194B", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
]


def load_data(tpm_path: Path, meta_path: Path) -> tuple:
    sep_t = "\t" if tpm_path.suffix in {".tsv", ".txt"} else ","
    sep_m = "\t" if meta_path.suffix in {".tsv", ".txt"} else ","
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sep_t)
    meta = pd.read_csv(meta_path, index_col=0, sep=sep_m)
    shared = tpm.columns.intersection(meta.index)
    return tpm[shared], meta.loc[shared]


def log_transform(tpm: pd.DataFrame) -> pd.DataFrame:
    return np.log2(tpm + 1)


def bray_curtis_dm(expr_log: pd.DataFrame) -> DistanceMatrix:
    """Bray-Curtis dissimilarity; samples are rows."""
    mat = expr_log.T.values
    bc  = pdist(mat, metric="braycurtis")
    ids = list(expr_log.columns)
    return DistanceMatrix(squareform(bc), ids=ids)


def run_ordination(dm: DistanceMatrix) -> tuple:
    result = pcoa(dm)
    return result


def run_tests(dm: DistanceMatrix, grouping: pd.Series, n_perms: int) -> dict:
    g = grouping.loc[dm.ids]   # align order
    perm_result = permanova(dm, g, permutations=n_perms)
    disp_result = permdisp(dm, g, permutations=n_perms)
    return {
        "PERMANOVA": {
            "pseudo-F": round(perm_result["test statistic"], 3),
            "p":        round(perm_result["p-value"], 4),
        },
        "PERMDISP": {
            "F":  round(disp_result["test statistic"], 3),
            "p":  round(disp_result["p-value"], 4),
        },
    }


def plot(
    pcoa_result,
    meta: pd.DataFrame,
    group_col: str,
    shape_col: str | None,
    stats: dict,
    output: Path,
):
    coords = pcoa_result.samples[["PC1", "PC2"]]
    prop   = pcoa_result.proportion_explained

    groups = meta[group_col].astype(str)
    unique_groups = sorted(groups.unique())
    color_map = {g: PALETTE[i % len(PALETTE)] for i, g in enumerate(unique_groups)}

    marker_map: dict = {}
    if shape_col and shape_col in meta.columns:
        shapes = meta[shape_col].astype(str)
        unique_shapes = sorted(shapes.unique())
        _markers = ["o", "s", "D", "^", "v", "<", ">", "P"]
        marker_map = {s: _markers[i % len(_markers)] for i, s in enumerate(unique_shapes)}

    fig, ax = plt.subplots(figsize=(7, 6))

    for sample in coords.index:
        x, y = coords.loc[sample, "PC1"], coords.loc[sample, "PC2"]
        c = color_map.get(groups.loc[sample], "#888888")
        m = marker_map.get(shapes.loc[sample], "o") if marker_map else "o"
        ax.scatter(x, y, color=c, marker=m, s=80, edgecolors="white", linewidths=0.5, zorder=3)

    ax.set_xlabel(f"PC1 ({prop['PC1']*100:.1f}%)", fontsize=11)
    ax.set_ylabel(f"PC2 ({prop['PC2']*100:.1f}%)", fontsize=11)
    ax.axhline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#cccccc", linewidth=0.8, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Stat annotation
    perm_str = (
        f"PERMANOVA: F={stats['PERMANOVA']['pseudo-F']}, p={stats['PERMANOVA']['p']}\n"
        f"PERMDISP:  F={stats['PERMDISP']['F']}, p={stats['PERMDISP']['p']}"
    )
    ax.text(0.02, 0.98, perm_str, transform=ax.transAxes,
            va="top", ha="left", fontsize=8,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))

    # Legend — groups (colour)
    handles = [mpatches.Patch(color=color_map[g], label=g) for g in unique_groups]
    if marker_map:
        from matplotlib.lines import Line2D
        shape_handles = [
            Line2D([0], [0], marker=marker_map[s], color="w",
                   markerfacecolor="grey", markersize=9, label=s)
            for s in sorted(marker_map)
        ]
        handles += shape_handles
    ax.legend(handles=handles, fontsize=8, frameon=True,
              loc="lower right", borderaxespad=0.5)

    plt.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] PCoA plot saved to {output}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="PCoA + PERMANOVA + PERMDISP for expression data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm",       required=True, type=Path)
    p.add_argument("--metadata",  required=True, type=Path)
    p.add_argument("--group-col", default="condition",
                   help="Metadata column for grouping (colour).")
    p.add_argument("--shape-col", default=None,
                   help="Optional metadata column for point shape.")
    p.add_argument("--n-perms",   type=int, default=999)
    p.add_argument("--output",    required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm, meta = load_data(args.tpm, args.metadata)
    print(f"[INFO] {tpm.shape[0]} genes × {tpm.shape[1]} samples")

    expr_log = log_transform(tpm)
    dm       = bray_curtis_dm(expr_log)
    result   = run_ordination(dm)

    if args.group_col not in meta.columns:
        raise SystemExit(f"[ERROR] --group-col '{args.group_col}' not in metadata.")

    stats = run_tests(dm, meta[args.group_col].astype(str), args.n_perms)
    print(f"[INFO] PERMANOVA pseudo-F={stats['PERMANOVA']['pseudo-F']}  "
          f"p={stats['PERMANOVA']['p']}")
    print(f"[INFO] PERMDISP   F={stats['PERMDISP']['F']}  "
          f"p={stats['PERMDISP']['p']}")

    plot(result, meta, args.group_col, args.shape_col, stats, args.output)


if __name__ == "__main__":
    main()
