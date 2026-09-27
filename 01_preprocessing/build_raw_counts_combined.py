"""
DEPRECATED / UNUSED for this project -- kept for reference only.

This script assumed the raw per-sample data was a flat directory of .txt
count files. That turned out not to match how the actual raw data is
organized: it is a nested tree of CLC Genomics Workbench .xlsx "Gene
Expression" track exports (FPKM/<strain>/<medium>/*.xlsx, +pXF51/ for the
9a5c plasmid). Use `build_matrices_from_fpkm_xlsx.py` instead -- see its
docstring for the actual raw-data layout and why TPM/RPKM must be
recomputed (via compute_tpm.py) rather than copied from the .xlsx columns.

Everything below this notice describes the now-superseded .txt-based
approach; it is left in place only in case a future raw-data batch really
does arrive as flat per-sample .txt files.

---

build_raw_counts_combined.py — Build data/raw_counts_combined.tsv from the
per-sample raw read-count .txt files (CLC Genomics Workbench exports), one
file per sample, each already quantified against its OWN strain's reference
genome (9a5c or Temecula1 IMG locus tags).

This reproduces Cell 3 ("load all count .txt files and build counts_df +
sample_df") of the original analysis notebook
(deseq2_pierry_feitosa_paper.ipynb) EXACTLY for the counts-matrix part:

    counts_df = pd.concat(count_dfs, axis=1, join="outer").fillna(0).astype(int)
    counts_df = counts_df.loc[counts_df.sum(axis=1) > 0]

No ortholog merging happens here. Each row is a strain-native gene_id (9a5c
IMG ID or Temecula1 IMG ID); a 9a5c gene's columns for Temecula1 samples are
0 (and vice-versa) simply because that gene ID was never in a Temecula1
sample's own count file — this is an outer join, not an intersection. The
ortholog merge (Temecula1 IMG -> 9a5c IMG, summed if many-to-one) is applied
later, per-comparison, ONLY for cross-strain comparisons, inside
run_deseq2.py -- exactly as it was in Cell 8 of the original notebook. Doing
the merge there instead of here is what keeps within-strain comparisons
(e.g. 9a5c_PIM6_1d vs 9a5c_PIM6_3d) tested against that strain's FULL native
gene set, with no ortholog-availability restriction -- see the "Background
strategy" note in 03_go_enrichment/run_go_enrichment.py, which documents the
same within-strain / cross-strain distinction for GO background sets.

Column naming
-------------
Each output column is named after the input filename's stem (extension
stripped) by default. This MUST end up matching the `sample_id` values used
in data/sample_info.csv, since run_deseq2.py indexes samples by
`sample_info.csv`'s sample_id. If your raw filenames don't already match
those IDs one-to-one, pass --rename-map (two-column TSV, no header:
<filename_stem><TAB><sample_id>) to rename columns during the build.

Usage
-----
python build_raw_counts_combined.py \\
    --counts-dir  "DeSeq2/" \\
    --output      data/raw_counts_combined.tsv

# With filename -> sample_id renaming:
python build_raw_counts_combined.py \\
    --counts-dir  "DeSeq2/" \\
    --rename-map  data/raw_counts_rename_map.tsv \\
    --output      data/raw_counts_combined.tsv

Input format (each .txt file, whitespace/tab-delimited, no header assumed;
only the first two columns are used):
    gene_id       count
    Xf9a_00001    245
    Xf9a_00002    0
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_one_count_file(path: Path, sample_name: str) -> pd.DataFrame:
    """Robustly read a 2-column (gene_id, count) whitespace- or tab-delimited
    file, exactly as Cell 3 does (whitespace-delim first, tab fallback)."""
    try:
        df = pd.read_csv(
            path, sep=r"\s+", header=None,
            usecols=[0, 1], names=["gene_id", sample_name],
        )
    except Exception:
        tmp = pd.read_csv(path, sep="\t", header=None, engine="python")
        tmp = tmp.iloc[:, :2]
        tmp.columns = ["gene_id", sample_name]
        df = tmp
    df = df.dropna(subset=["gene_id"]).set_index("gene_id")
    return df


def load_rename_map(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", header=None, names=["stem", "sample_id"])
    return dict(zip(df["stem"].astype(str).str.strip(),
                     df["sample_id"].astype(str).str.strip()))


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Concatenate per-sample raw count .txt files into "
                     "data/raw_counts_combined.tsv (Cell 3 of the original "
                     "notebook; no ortholog merge at this stage).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--counts-dir", required=True, type=Path,
                   help="Directory containing per-sample raw-count .txt files.")
    p.add_argument("--pattern", default="*.txt",
                   help="Glob pattern (relative to --counts-dir) selecting "
                        "which files to include.")
    p.add_argument("--rename-map", type=Path, default=None,
                   help="Optional TSV (no header): <filename_stem>\\t<sample_id>. "
                        "Renames output columns from filename stems to the "
                        "sample_id values used in data/sample_info.csv. "
                        "Files whose stem is not in the map are kept with "
                        "their original stem as the column name (with a "
                        "warning), NOT dropped.")
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--keep-all-zero-genes", action="store_true",
                   help="By default, genes with 0 total counts across every "
                        "sample are dropped (matches the original notebook's "
                        "counts_df.loc[counts_df.sum(axis=1) > 0]). Pass this "
                        "flag to keep them.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    files = sorted(args.counts_dir.glob(args.pattern))
    if not files:
        sys.exit(f"[ERROR] No files matching {args.pattern!r} found in {args.counts_dir}")
    print(f"[INFO] Found {len(files)} count files in {args.counts_dir}")

    rename_map = load_rename_map(args.rename_map) if args.rename_map else {}

    count_dfs = []
    used_names = []
    for path in files:
        stem = path.stem
        sample_name = rename_map.get(stem, stem)
        if rename_map and stem not in rename_map:
            print(f"[WARN] {path.name}: stem {stem!r} not found in --rename-map; "
                  f"keeping column name as-is.")
        count_dfs.append(load_one_count_file(path, sample_name))
        used_names.append(sample_name)

    dupes = {n for n in used_names if used_names.count(n) > 1}
    if dupes:
        sys.exit(f"[ERROR] Duplicate output column names after renaming: {sorted(dupes)}. "
                  f"Fix --rename-map or filenames.")

    # Outer join across all samples, zero-fill missing (strain-exclusive genes),
    # exactly as Cell 3: counts_df = pd.concat(..., axis=1, join="outer").fillna(0).astype(int)
    counts_df = pd.concat(count_dfs, axis=1, join="outer").fillna(0).astype(int)
    print(f"[INFO] Combined matrix (before zero-row filter): "
          f"{counts_df.shape[0]} genes x {counts_df.shape[1]} samples")

    if not args.keep_all_zero_genes:
        before = counts_df.shape[0]
        counts_df = counts_df.loc[counts_df.sum(axis=1) > 0]
        print(f"[INFO] Dropped {before - counts_df.shape[0]} all-zero genes "
              f"(matches original notebook behaviour).")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    counts_df.to_csv(args.output, sep="\t")
    print(f"[INFO] Saved {counts_df.shape[0]} genes x {counts_df.shape[1]} samples "
          f"to {args.output}")
    print("[INFO] Columns:", ", ".join(counts_df.columns))
    print("[INFO] Reminder: column names must match sample_info.csv's sample_id "
          "values for run_deseq2.py to find them.")


if __name__ == "__main__":
    main()
