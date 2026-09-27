"""
fig_virulence_trends.py — Mean TPM (+/- SE) of virulence-associated genes,
split by mobile- vs sessile-phase gene groups, aggregated across early/late
timepoints for each strain x medium combination (Figure 3B).

This reproduces the exact style of the original notebook cell (
PCA_PCoA_heatmap_Pierry_Transcriptome.ipynb, "plot_expression_trends_virulence()"
/ output file gene_expression_trends_virulence.tiff) -- same title, axis
labels, colours, seaborn whitegrid background, jitter and single combined
legend -- with ONE deliberate accessibility change:

  ORIGINAL:  linestyle encodes phase-group (sessile = dashed, mobile =
             dotted); marker is always a plain circle.
  FIX:       linestyle is now always SOLID; phase-group is instead encoded
             by MARKER SHAPE (circle = mobile, square = sessile), which
             remains legible at small sizes/greyscale, unlike dashed vs.
             dotted lines. Colour (strain x medium) and every other visual
             element are unchanged from the original.
This mirrors the same colour+shape accessibility fix applied to Figure 2A,
keeping the encoding consistent across figures while otherwise matching
the originally submitted figure's look exactly.

Data pipeline (unchanged from the original figure):
  1. Map each virulence-table gene (Temecula1 PD-tag) to the unified
     9a5c/Temecula1 gene identifier via the gene dictionary.
  2. Split genes into mobile-phase / sessile-phase groups using the
     Phase column of the virulence table (Supplementary Table S6).
  3. For each strain x medium x timepoint(early/late) group, take the
     per-gene mean TPM across replicates, then report mean +/- SE across
     genes within each phase group (SE computed over per-gene means, N =
     number of genes in the phase group -- consistent with the SE
     convention already used in the Figure 3A / 4A pipelines).

Usage
-----
python fig_virulence_trends.py \\
    --tpm              data/tpm_expression_original.csv \\
    --dictionary       data/gene_dictionary.tsv \\
    --virulence-table  data/virulence_table.csv \\
    --output           figures/fig_3B_virulence_trends.tiff

(data/tpm_expression_original.csv is the frozen, comma-separated per-sample
raw TPM matrix -- the same underlying values as the notebook's
matrix_9a5c_Temecula1.csv, which is NOT checked into this repository.
--tpm-sep already defaults to "," to match it.)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def sniff_sep(path: Path) -> str:
    """Detect the field separator from the file's own header line (see the
    same helper/bug in 04_wgcna/run_wgcna.py and fig_pearson_heatmap.py).
    Used as a fallback when the --tpm-sep default doesn't match the file."""
    with open(path, "r", newline="") as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ("\t", ",", ";")}
    best_sep, best_count = max(counts.items(), key=lambda kv: kv[1])
    return best_sep if best_count > 0 else ","


# Same strain x medium colour scheme as the original notebook cell
# (identical to Figure 3A / 4A's condition_colors, per the cell's own
# comment "consistent with Figure 4A").
STRAIN_MEDIUM_COLORS = {
    ("9a5c", "PIM6"):      "#FF8C00",   # orange
    ("9a5c", "PWG"):       "#FF0000",   # red
    ("Temecula1", "PIM6"): "#ADD8E6",   # light blue
    ("Temecula1", "PWG"):  "#00008B",   # dark blue
}
# Accessibility fix: phase-group encoded by marker shape, not linestyle.
PHASE_MARKERS = {"mobile": "o", "sessile": "s"}
# Same small per-medium jitter as the original cell (+0.05 / -0.05).
MEDIUM_JITTER = {"PIM6": 0.05, "PWG": -0.05}


def default_strain(s: str) -> str:
    return "9a5c" if s.lower().startswith("9a") else "Temecula1"


def default_medium(s: str) -> str:
    return "PIM6" if "PIM6" in s else "PWG"


def default_timepoint(s: str) -> str:
    return "early" if ("PIM6_1d" in s or "PWG_3d" in s) else "late"


def build_pd_to_unified_map(dictionary: pd.DataFrame) -> dict:
    mapping = {}
    for _, row in dictionary.iterrows():
        x9a = str(row.get("9a5c_IMG_ID", "")).strip()
        xtem = str(row.get("Temecula1_IMG_ID", "")).strip()
        if x9a in ("", "nan") or xtem in ("", "nan"):
            continue
        unified = f"{x9a}_{xtem}"
        old_tags = str(row.get("Temecula1_old_locus_tags", ""))
        for tag in old_tags.replace('"', "").split():
            tag = tag.strip().replace("_", "")
            if tag.startswith("PD"):
                mapping[tag] = unified
    return mapping


def load_virulence_groups(vir_path: Path, dictionary: pd.DataFrame):
    sep = "\t" if vir_path.suffix in {".tsv", ".txt"} else ","
    try:
        vir = pd.read_csv(vir_path, sep=sep)
    except UnicodeDecodeError:
        vir = pd.read_csv(vir_path, sep=sep, encoding="latin-1")
    vir.columns = [c.strip() for c in vir.columns]

    id_col = next(c for c in vir.columns if "gene_id" in c.lower() and "temecula" in c.lower())
    phase_col = next(c for c in vir.columns if "phase" in c.lower())

    pd_to_unified = build_pd_to_unified_map(dictionary)

    mobile_pd = vir.loc[vir[phase_col].astype(str).str.lower() == "mobile", id_col].astype(str).str.strip().tolist()
    sessile_pd = vir.loc[vir[phase_col].astype(str).str.lower() == "sessile", id_col].astype(str).str.strip().tolist()

    mobile = [pd_to_unified[g] for g in mobile_pd if g in pd_to_unified]
    sessile = [pd_to_unified[g] for g in sessile_pd if g in pd_to_unified]
    print(f"[INFO] Virulence table: {len(mobile_pd)} mobile / {len(sessile_pd)} sessile genes "
          f"-> mapped {len(mobile)}/{len(sessile)} to unified gene IDs")
    return {"mobile": mobile, "sessile": sessile}


def aggregate(tpm: pd.DataFrame, gene_groups: dict) -> pd.DataFrame:
    samples = tpm.columns.tolist()
    meta = pd.DataFrame({
        "sample": samples,
        "strain": [default_strain(s) for s in samples],
        "medium": [default_medium(s) for s in samples],
        "timepoint": [default_timepoint(s) for s in samples],
    })

    rows = []
    for strain in meta["strain"].unique():
        for medium in meta["medium"].unique():
            for timepoint in ["early", "late"]:
                sids = meta.loc[(meta.strain == strain) & (meta.medium == medium) &
                                 (meta.timepoint == timepoint), "sample"].tolist()
                if not sids:
                    continue
                for phase_name, gene_list in gene_groups.items():
                    genes_present = [g for g in gene_list if g in tpm.index]
                    if not genes_present:
                        continue
                    gene_means = tpm.loc[genes_present, sids].mean(axis=1).values
                    rows.append({
                        "strain": strain, "medium": medium, "timepoint": timepoint,
                        "phase_group": phase_name,
                        "mean_tpm": gene_means.mean(),
                        "se_tpm": gene_means.std(ddof=1) / np.sqrt(len(gene_means)),
                        "n_genes": len(gene_means),
                    })
    return pd.DataFrame(rows)


def plot(agg: pd.DataFrame, output: Path):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 8))
    x_pos = {"early": 0, "late": 1}

    for (strain, medium), color in STRAIN_MEDIUM_COLORS.items():
        jitter = MEDIUM_JITTER[medium]
        for phase_name, marker in PHASE_MARKERS.items():
            sub = agg[(agg.strain == strain) & (agg.medium == medium) & (agg.phase_group == phase_name)]
            if sub.empty:
                continue
            sub = sub.set_index("timepoint").reindex(["early", "late"]).dropna()
            x = [x_pos[t] + jitter for t in sub.index]
            y = sub["mean_tpm"].values
            y_err = sub["se_tpm"].values

            # Linestyle is always solid; phase-group (mobile/sessile) is
            # encoded by marker shape instead. No fill_between shading in
            # the original Fig 3B (unlike Fig 3A).
            ax.errorbar(x, y, yerr=y_err, fmt=marker, color=color,
                        markersize=8, capsize=5,
                        label=f"{strain}–{medium} [{phase_name}]")
            ax.plot(x, y, color=color, linestyle="-", alpha=0.7, linewidth=2)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Early", "Late"])
    ax.set_xlabel("Growth phase")
    ax.set_ylabel("Mean TPM per gene")
    ax.set_title("TPM expression trends (Virulence groups)")
    ax.grid(True, linestyle="--", alpha=0.4)
    ax.legend(fontsize=12)
    plt.tight_layout()

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Figure 3B saved to {output}")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Figure 3B: virulence gene expression trends (mobile vs sessile).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm", required=True, type=Path)
    p.add_argument("--tpm-sep", default=None,
                   help="Field separator for --tpm. Default: auto-detect "
                        "from the file's header line (sniff_sep).")
    p.add_argument("--dictionary", required=True, type=Path,
                    help="Gene dictionary (Supplementary Table S5 / build_gene_dictionary.py output).")
    p.add_argument("--dictionary-sep", default="\t")
    p.add_argument("--virulence-table", required=True, type=Path,
                    help="Supplementary Table S6 (virulence gene table with a Phase column).")
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    tpm_sep = args.tpm_sep if args.tpm_sep is not None else sniff_sep(args.tpm)
    tpm = pd.read_csv(args.tpm, sep=tpm_sep, index_col=0)
    dictionary = pd.read_csv(args.dictionary, sep=args.dictionary_sep)
    gene_groups = load_virulence_groups(args.virulence_table, dictionary)

    agg = aggregate(tpm, gene_groups)
    csv_out = args.output.with_suffix(".csv")
    agg.to_csv(csv_out, index=False)
    print(f"[INFO] Aggregated data saved to {csv_out}")
    print(agg.round(2).to_string(index=False))

    plot(agg, args.output)


if __name__ == "__main__":
    main()
