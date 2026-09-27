"""
module_trait_heatmap.py — Module-trait correlation heatmap with BH-adjusted
p-values.

Reporting adjusted p-values alongside the correlation coefficients (rather
than the coefficients alone) avoids overstating the confidence of module-
trait associations tested across multiple modules and traits at once.

Colour scheme: green (#385c42) -> white (#f7f7f7) -> orange (#c77a30),
matching Cell 9 of the original WGCNA notebook and the figure submitted to
MDPI Pathogens. The BH-adjusted-p significance annotations added to the
cell text do not change this palette.

This is a pure NumPy/pandas reimplementation of run_wgcna.py's
correlate_traits() + plot_trait_heatmap() (which use scipy.stats.pearsonr
and statsmodels multipletests) — kept dependency-light so it can be re-run
directly from module_eigengenes.csv + sample_info.csv without needing
scipy/statsmodels installed. Produces numerically identical results
(validated against scipy.stats.pearsonr / statsmodels fdr_bh to <1e-10).

Usage
-----
python module_trait_heatmap.py \\
    --eigengenes  results/WGCNA/module_eigengenes.csv \\
    --metadata    data/sample_info.csv \\
    --traits      strain medium timepoint \\
    --output      figures/fig_module_trait_heatmap.tiff
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).parent.parent / "06_stats"))
from pure_stats import pearsonr_manual, benjamini_hochberg


# Same old-vs-new sample-naming reconciliation as 04_wgcna/run_wgcna.py
# (see that module's own comment for the full rationale). Needed here too:
# run_wgcna.py normalises sample names to the canonical "9a5c_..._rep#" form
# before saving module_eigengenes.csv, so that file's index is ALREADY
# canonical even when it was built from an old-style --tpm/--metadata pair
# (e.g. data/tpm_expression_original.csv + data/sample_info_original.csv,
# whose own sample names are "9a_PIM6_1d_1"-style). Loading
# sample_info_original.csv fresh here without the same normalisation left
# its index in the old form -- zero overlap with the already-canonical
# module_eigengenes.csv index, hence "No overlapping sample IDs".
_SAMPLE_RE = re.compile(
    r"^(?P<strain>9a5c|9a|Tem1|Temecula1)_(?P<medium>[A-Za-z0-9]+)_"
    r"(?P<timepoint>\d+[dh])_(?:rep)?(?P<rep>\d+)$"
)
_STRAIN_CANON = {"9a5c": "9a5c", "9a": "9a5c", "Tem1": "Tem1", "Temecula1": "Tem1"}


def normalize_sample_name(name: str) -> str:
    name = str(name).strip()
    m = _SAMPLE_RE.match(name)
    if not m:
        return name
    strain = _STRAIN_CANON.get(m.group("strain"), m.group("strain"))
    return f"{strain}_{m.group('medium')}_{m.group('timepoint')}_rep{int(m.group('rep'))}"


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

# Same green -> white -> orange diverging colormap as the original submitted
# figure (Cell 9 of the WGCNA notebook); adding the BH-adjusted-p
# annotations only changed the cell text/statistics, not the color scheme.
WGCNA_CMAP = LinearSegmentedColormap.from_list(
    "green_white_orange", ["#385c42", "#f7f7f7", "#c77a30"], N=256,
)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Module-trait correlation heatmap with BH-adjusted p-values.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--eigengenes", required=True, type=Path,
                    help="module_eigengenes.csv (samples x modules) from run_wgcna.py.")
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--traits", nargs="*", default=None,
                    help="Metadata columns to correlate against (default: all non-numeric "
                         "columns are one-hot encoded, all numeric columns used directly).")
    p.add_argument("--drop-cols", nargs="*", default=["replicate"],
                    help="Metadata columns to exclude (e.g. sample-level nuisance vars).")
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args(argv)

    ME = pd.read_csv(args.eigengenes, index_col=0)
    ME.index = ME.index.map(normalize_sample_name)
    meta = pd.read_csv(args.metadata, index_col=0, sep=sniff_sep(args.metadata))
    meta.columns = [c.strip() for c in meta.columns]
    meta.index = meta.index.map(normalize_sample_name)
    for c in meta.select_dtypes(include="object").columns:
        meta[c] = meta[c].str.strip()

    common = ME.index.intersection(meta.index)
    if len(common) == 0:
        raise SystemExit("[ERROR] No overlapping sample IDs between eigengenes and metadata.")
    ME, meta = ME.loc[common], meta.loc[common]

    meta = meta.drop(columns=[c for c in args.drop_cols if c in meta.columns])
    if args.traits:
        meta = meta[[c for c in args.traits if c in meta.columns]]

    traits = meta.copy()
    for c in traits.columns:
        if traits[c].dtype == object:
            d = pd.get_dummies(traits[c], prefix=c)
            traits = traits.drop(columns=c).join(d)
    traits = traits.loc[ME.index]

    rows = []
    for m in ME.columns:
        for t in traits.columns:
            r, pv = pearsonr_manual(ME[m].values.astype(float), traits[t].values.astype(float))
            rows.append({"module": m, "trait": t, "r": r, "p": pv})
    res = pd.DataFrame(rows)
    res["p_BH"] = benjamini_hochberg(res["p"].values)

    csv_out = args.output.with_suffix(".csv")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    res.to_csv(csv_out, index=False)
    print(f"[INFO] Module-trait correlations (with BH-adjusted p) saved to {csv_out}")
    print(res.round(4).to_string(index=False))

    pivot_r = res.pivot(index="module", columns="trait", values="r")
    pivot_p = res.pivot(index="module", columns="trait", values="p_BH")

    fig, ax = plt.subplots(figsize=(max(5, len(pivot_r.columns) * 1.4),
                                     max(3.5, len(pivot_r) * 1.3)))
    im = ax.imshow(pivot_r.values, cmap=WGCNA_CMAP, vmin=-1, vmax=1, aspect="auto")
    for i in range(pivot_r.shape[0]):
        for j in range(pivot_r.shape[1]):
            r = pivot_r.values[i, j]
            pv = pivot_p.values[i, j]
            sig = "***" if pv < 0.001 else ("**" if pv < 0.01 else ("*" if pv < 0.05 else "ns"))
            ax.text(j, i, f"{r:.2f}\n(FDR {sig})", ha="center", va="center",
                     fontsize=8.5, color="white" if abs(r) > 0.6 else "black")
    ax.set_xticks(range(len(pivot_r.columns)))
    ax.set_xticklabels(pivot_r.columns, fontsize=9, rotation=20, ha="right")
    ax.set_yticks(range(len(pivot_r.index)))
    ax.set_yticklabels(pivot_r.index, fontsize=9)
    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Pearson r (module eigengene vs. trait)", fontsize=9)
    ax.set_title("Module–trait correlations\n"
                  "(BH-adjusted P: *** <0.001, ** <0.01, * <0.05, ns ≥ 0.05)", fontsize=10)
    plt.tight_layout()
    fig.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Heatmap saved to {args.output}")


if __name__ == "__main__":
    main()
