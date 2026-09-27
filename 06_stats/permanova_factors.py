"""
permanova_factors.py — Factor-specific PERMANOVA (strain / medium / growth
phase).

The global 8-group PERMANOVA/PERMDISP run in fig_pcoa.py establishes that
the eight strain x medium x phase groups differ significantly overall, but
does not by itself say which factor(s) drive that separation, nor whether
one factor's effect is a confound of another. This script runs THREE
complementary tests per factor, so the confound each addresses is
explicit:

  1. Global one-way PERMANOVA per factor (all 24 samples, ignoring the
     other factors) — matches a naive "does this factor separate groups"
     test, but can be inflated/masked by the other factors' variance.
  2. Stratified PERMANOVA (permutations restricted within strain) — tests
     medium/phase while controlling for the strain confound, analogous to
     vegan::adonis2(..., strata = strain).
  3. Within-strain subset PERMANOVA — an even more conservative check:
     re-runs the medium/phase test independently inside each strain's
     12-sample subset.

Distance metric
----------------
The original analysis notebook (PCA_PCoA_heatmap_Pierry_Transcriptome.ipynb)
computes pdist() on RAW (non-log-transformed) TPM with the scipy default
metric, i.e. Euclidean distance on raw TPM -- this is what the published
Figure 2A / Results text ("PCoA on Euclidean distances... PC1 52.9%")
reflects, and what fig_pcoa.py and this script both use by default
(--log-transform is OFF unless explicitly passed). Log2(TPM+1) or
Bray-Curtis are available via --log-transform/--metric for sensitivity
comparison, but are not the metric used for the reported results.

Usage
-----
python permanova_factors.py \\
    --tpm       data/tpm_expression_original.csv --tpm-sep "," \\
    --output    results/permanova/permanova_factors.csv \\
    --n-perms   999

(data/tpm_expression_original.csv is the same per-sample raw TPM matrix as
the notebook's matrix_9a5c_Temecula1.csv, just comma- instead of
semicolon-separated -- values are identical. Default flags now already
match the notebook: euclidean metric, no log-transform.)

Input
-----
--tpm : genes x samples TPM matrix (CSV). Sample names must encode strain
        (starts with "9a" / "Tem" or contains "9a5c"/"Temecula1"), medium
        ("PIM6"/"PWG") and timepoint in a way --strain-fn/--medium-fn/
        --phase-fn (or the defaults below) can parse. Defaults match the
        sample-naming convention already used throughout this project
        (e.g. "9a_PIM6_1d_1", "Tem1_PWG_7d_3").
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from pure_stats import permanova, permanova_stratified, euclidean_distance_matrix, bray_curtis_distance_matrix


def default_strain(s: str) -> str:
    return "9a5c" if s.lower().startswith("9a") else "Temecula1"


def default_medium(s: str) -> str:
    return "PIM6" if "PIM6" in s else "PWG"


def default_phase(s: str) -> str:
    # Matches the sample_info.csv / meta.csv convention used in this project:
    # PIM6 early=1d, late=3d; PWG early=3d, late=7d/10d.
    return "early" if ("PIM6_1d" in s or "PWG_3d" in s) else "late"


def build_distance_matrix(expr_samples_x_genes: np.ndarray, metric: str) -> np.ndarray:
    if metric == "euclidean":
        return euclidean_distance_matrix(expr_samples_x_genes)
    elif metric == "braycurtis":
        return bray_curtis_distance_matrix(expr_samples_x_genes)
    else:
        raise ValueError(f"Unknown metric: {metric}")


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Factor-specific PERMANOVA (strain, medium, phase).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm", required=True, type=Path)
    p.add_argument("--tpm-sep", default=",", help="Field separator for --tpm (',' or ';').")
    p.add_argument("--metric", choices=["euclidean", "braycurtis"], default="euclidean")
    p.add_argument("--log-transform", action="store_true",
                    help="Apply log2(TPM+1) before the distance calc. OFF by default, "
                         "matching the original notebook and fig_pcoa.py (raw TPM, "
                         "no log-transform).")
    p.add_argument("--n-perms", type=int, default=999)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output", required=True, type=Path)
    args = p.parse_args(argv)

    tpm = pd.read_csv(args.tpm, sep=args.tpm_sep, index_col=0)
    print(f"[INFO] TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")

    expr = np.log2(tpm + 1) if args.log_transform else tpm
    X = expr.T  # samples x genes
    samples = X.index.tolist()

    strain = np.array([default_strain(s) for s in samples])
    medium = np.array([default_medium(s) for s in samples])
    phase = np.array([default_phase(s) for s in samples])

    print(f"[INFO] strain groups: {dict(zip(*np.unique(strain, return_counts=True)))}")
    print(f"[INFO] medium groups: {dict(zip(*np.unique(medium, return_counts=True)))}")
    print(f"[INFO] phase groups:  {dict(zip(*np.unique(phase, return_counts=True)))}")

    D = build_distance_matrix(X.values, args.metric)
    print(f"[INFO] Distance matrix built (metric={args.metric}, log-transform={args.log_transform})")

    rows = []

    # 1. Global one-way per factor
    for name, labels in [("strain", strain), ("medium", medium), ("phase", phase)]:
        res = permanova(D, labels, permutations=args.n_perms, seed=args.seed)
        rows.append({"test": "global_one_way", "factor": name,
                     "F": res["test statistic"], "R2": res["R2"], "p": res["p-value"],
                     "n_perms": args.n_perms})

    # 2. Global 8-group sanity check (strain x medium x phase)
    grouping_8 = np.array([f"{a}_{b}_{c}" for a, b, c in zip(strain, medium, phase)])
    res8 = permanova(D, grouping_8, permutations=args.n_perms, seed=args.seed)
    rows.append({"test": "global_8group_sanity_check", "factor": "strain*medium*phase",
                 "F": res8["test statistic"], "R2": res8["R2"], "p": res8["p-value"],
                 "n_perms": args.n_perms})

    # 3. Stratified (permute within strain) for medium and phase
    for name, labels in [("medium", medium), ("phase", phase)]:
        res = permanova_stratified(D, labels, strata=strain, permutations=args.n_perms, seed=args.seed)
        rows.append({"test": "stratified_by_strain", "factor": name,
                     "F": res["test statistic"], "R2": res["R2"], "p": res["p-value"],
                     "n_perms": args.n_perms})

    # 4. Within-strain independent subsets for medium and phase
    for s in np.unique(strain):
        idx = np.where(strain == s)[0]
        Dsub = D[np.ix_(idx, idx)]
        for name, labels in [("medium", medium[idx]), ("phase", phase[idx])]:
            res = permanova(Dsub, labels, permutations=args.n_perms, seed=args.seed)
            rows.append({"test": f"within_strain_{s}", "factor": name,
                         "F": res["test statistic"], "R2": res["R2"], "p": res["p-value"],
                         "n_perms": args.n_perms})

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"\n[INFO] Results written to {args.output}\n")
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
