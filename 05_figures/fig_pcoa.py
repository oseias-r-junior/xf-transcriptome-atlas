"""
fig_pcoa.py — Principal Coordinates Analysis (PCoA) of TPM expression profiles
with PERMANOVA and PERMDISP significance tests (Figure 2 / Supplementary).

Reproduces the original notebook workflow
(PCA_PCoA_heatmap_Pierry_Transcriptome.ipynb), not a from-scratch PCoA design:

  * Distance metric: EUCLIDEAN distance on the raw TPM matrix (scipy.spatial.
    distance.pdist's default metric; samples are rows). No log-transform is
    applied before the distance computation -- the notebook calls
    ``pdist(transcriptome_pcoa_t)`` directly on the untransformed, transposed
    TPM table (cells 5-7 of the notebook). This intentionally differs from
    the WGCNA pipeline's log2(TPM+1) convention; do not "fix" it to match.
  * Ordination: skbio.stats.ordination.pcoa() on the resulting DistanceMatrix
    -- same call the notebook makes (cell 11). Emperor (cell 12-13) is purely
    a visualization layer on top of this same pcoa() result; it does not
    change the ordination itself. This script reproduces Emperor's output by
    (a) keeping the identical pcoa() computation, and (b) optionally writing
    Emperor's own standalone interactive HTML via --emperor-html, using the
    same ``Emperor(pcoa_result, metadata_df, remote='.')`` /
    ``make_emperor(standalone=True)`` call as the notebook. The static
    publication figure (matplotlib scatter, --output) is generated from the
    same coordinates for the manuscript's print figure.
  * PERMANOVA: skbio.stats.distance.permanova(dm, grouping), default
    permutations=999 (notebook cell 15 -- no explicit `permutations` kwarg).
  * PERMDISP: skbio.stats.distance.permdisp(dm, grouping, permutations=9999)
    (notebook cell 16 -- note this is a DIFFERENT permutation count than
    PERMANOVA's; --n-perms-permanova / --n-perms-permdisp are separate flags
    for this reason, defaulting to 999 / 9999 respectively).
  * Grouping: a plain per-sample condition label (e.g. "9a_PIM6_1d"), not
    including replicate number -- matches --group-col.

Usage
-----
python fig_pcoa.py \\
    --tpm                 data/tpm_expression.csv \\
    --metadata            data/sample_info.csv \\
    --group-col           condition \\
    --n-perms-permanova   999 \\
    --n-perms-permdisp    9999 \\
    --output              figures/fig_pcoa.tiff \\
    --emperor-html        figures/emperor_pcoa

Arguments
---------
--group-col   Column in metadata used to colour/group samples (default: condition).
--shape-col   Optional second metadata column for point shapes (e.g., strain).
--emperor-html  Optional output directory; when given, also writes Emperor's
                standalone interactive HTML (index.html + support files),
                exactly mirroring notebook cells 12-13.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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


def sniff_sep(path: Path) -> str:
    """Detect the field separator from the file's own header line.

    File extension alone is unreliable (see the same helper/bug in
    04_wgcna/run_wgcna.py and 05_figures/fig_pearson_heatmap.py).
    """
    with open(path, "r", newline="") as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ("\t", ",", ";")}
    best_sep, best_count = max(counts.items(), key=lambda kv: kv[1])
    return best_sep if best_count > 0 else ","


def load_data(tpm_path: Path, meta_path: Path, condition_col: str = "condition") -> tuple:
    tpm  = pd.read_csv(tpm_path, index_col=0, sep=sniff_sep(tpm_path))
    meta = pd.read_csv(meta_path, index_col=0, sep=sniff_sep(meta_path))
    meta.columns = [c.strip() for c in meta.columns]
    for c in meta.select_dtypes(include="object").columns:
        meta[c] = meta[c].str.strip()

    if condition_col not in meta.columns:
        # data/sample_info_original.csv has no "condition" column, only
        # strain/medium/timepoint/replicate -- same situation already
        # handled in fig_pearson_heatmap.py. strain+medium+timepoint
        # together uniquely identify each of the 8 conditions.
        needed = {"strain", "medium", "timepoint"}
        if needed.issubset(meta.columns):
            meta[condition_col] = (
                meta["strain"].astype(str) + "_" +
                meta["medium"].astype(str) + "_" +
                meta["timepoint"].astype(str)
            )
            print(f"[INFO] {condition_col!r} not in {meta_path.name}; derived "
                  f"it from strain+medium+timepoint instead "
                  f"(e.g. {meta[condition_col].iloc[0]!r}).")

    shared = tpm.columns.intersection(meta.index)
    return tpm[shared], meta.loc[shared]


def euclidean_dm(tpm: pd.DataFrame) -> DistanceMatrix:
    """Euclidean distance on raw TPM values; samples are rows.

    Mirrors notebook cells 6-8 exactly: transpose (genes x samples ->
    samples x genes), then pdist with scipy's default metric (euclidean),
    with NO log-transform applied beforehand.
    """
    mat = tpm.T.values
    euc = pdist(mat)  # default metric="euclidean"
    ids = list(tpm.columns)
    return DistanceMatrix(squareform(euc), ids=ids)


def run_ordination(dm: DistanceMatrix) -> tuple:
    result = pcoa(dm)
    return result


def write_emperor_html(pcoa_result, meta: pd.DataFrame, outdir: Path):
    """Write Emperor's standalone interactive HTML, mirroring notebook
    cells 12-13 (Emperor(...).make_emperor(standalone=True))."""
    try:
        from emperor import Emperor
    except ImportError:
        print("[WARN] emperor package not installed; skipping --emperor-html "
              "output (pip install emperor). The static figure and stats "
              "are unaffected -- Emperor is a visualization layer only.")
        return
    outdir.mkdir(parents=True, exist_ok=True)
    emp = Emperor(pcoa_result, meta, remote=".")
    with open(outdir / "index.html", "w") as f:
        f.write(emp.make_emperor(standalone=True))
    emp.copy_support_files(str(outdir))
    print(f"[INFO] Emperor interactive PCoA saved to {outdir / 'index.html'}")


def run_tests(
    dm: DistanceMatrix,
    grouping: pd.Series,
    n_perms_permanova: int,
    n_perms_permdisp: int,
) -> dict:
    # dm.ids is a tuple, not a list -- pandas .loc[a_tuple] is parsed as a
    # multi-axis indexer ("Too many indexers"), not as a list of row labels.
    # Wrap in list() to force label-based fancy indexing.
    g = grouping.loc[list(dm.ids)]   # align order
    perm_result = permanova(dm, g, permutations=n_perms_permanova)
    disp_result = permdisp(dm, g, permutations=n_perms_permdisp)
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
    p.add_argument("--n-perms-permanova", type=int, default=999,
                   help="PERMANOVA permutations (notebook default: 999).")
    p.add_argument("--n-perms-permdisp", type=int, default=9999,
                   help="PERMDISP permutations (notebook default: 9999 -- "
                        "note this differs from PERMANOVA's).")
    p.add_argument("--output",    required=True, type=Path)
    p.add_argument("--emperor-html", type=Path, default=None,
                   help="Optional output directory for Emperor's standalone "
                        "interactive HTML (mirrors notebook cells 12-13). "
                        "Requires the 'emperor' package.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm, meta = load_data(args.tpm, args.metadata, args.group_col)
    print(f"[INFO] {tpm.shape[0]} genes × {tpm.shape[1]} samples")

    dm     = euclidean_dm(tpm)
    result = run_ordination(dm)

    if args.group_col not in meta.columns:
        raise SystemExit(f"[ERROR] --group-col '{args.group_col}' not in metadata.")

    stats = run_tests(
        dm, meta[args.group_col].astype(str),
        args.n_perms_permanova, args.n_perms_permdisp,
    )
    print(f"[INFO] PERMANOVA pseudo-F={stats['PERMANOVA']['pseudo-F']}  "
          f"p={stats['PERMANOVA']['p']}  (permutations={args.n_perms_permanova})")
    print(f"[INFO] PERMDISP   F={stats['PERMDISP']['F']}  "
          f"p={stats['PERMDISP']['p']}  (permutations={args.n_perms_permdisp})")

    plot(result, meta, args.group_col, args.shape_col, stats, args.output)

    if args.emperor_html is not None:
        write_emperor_html(result, meta, args.emperor_html)


if __name__ == "__main__":
    main()
