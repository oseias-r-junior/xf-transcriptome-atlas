"""
compute_tpm.py — Convert raw read counts or FPKM tables to TPM.

TPM (Transcripts Per Million) normalises for both sequencing depth and
gene length, allowing within- and cross-sample comparisons.

Formula:
    RPK_i  = counts_i / gene_length_kb_i
    TPM_i  = RPK_i / sum(RPK) * 1e6

If the input is already in FPKM/RPKM, gene-length normalisation has already
been applied; set --from-fpkm to skip the length step and only rescale to
per-million.

Usage examples
--------------
# From raw counts (requires --lengths):
python compute_tpm.py \\
    --counts data/raw_counts_9a5c.tsv \\
    --lengths data/gene_lengths.tsv \\
    --output results/tpm_9a5c.tsv

# From FPKM table (no length file needed):
python compute_tpm.py \\
    --counts data/fpkm_temecula1.tsv \\
    --from-fpkm \\
    --output results/tpm_temecula1.tsv

Input format
------------
TSV/CSV with genes as rows and samples as columns.
First column must be the gene identifier (IMG ID or locus tag).

    gene_id         sample_1   sample_2   ...
    Xf9a_00001      245        312        ...
    Xf9a_00002      0          18         ...

Gene-length file (--lengths): two-column TSV, no header:
    Xf9a_00001   1245
    Xf9a_00002   987
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def fpkm_to_tpm(fpkm: pd.DataFrame) -> pd.DataFrame:
    """Rescale an FPKM matrix to TPM (column-wise)."""
    col_sums = fpkm.sum(axis=0)
    return fpkm.div(col_sums, axis=1) * 1e6


def counts_to_tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """Convert raw counts to TPM given gene lengths (in bp)."""
    lengths_kb = lengths.reindex(counts.index)
    if lengths_kb.isna().any():
        missing = lengths_kb[lengths_kb.isna()].index.tolist()
        print(
            f"[WARN] {len(missing)} genes have no length information and will "
            f"receive TPM = 0. First 5: {missing[:5]}",
            file=sys.stderr,
        )
        lengths_kb = lengths_kb.fillna(np.nan)

    lengths_kb = lengths_kb / 1000.0          # bp → kb
    rpk = counts.div(lengths_kb, axis=0)       # reads per kilobase
    col_sums = rpk.sum(axis=0)
    tpm = rpk.div(col_sums, axis=1) * 1e6
    return tpm.fillna(0.0)


# ---------------------------------------------------------------------------
# Coverage helper (reported in Supplementary Table S1)
# ---------------------------------------------------------------------------

def compute_coverage(read_length: int, n_reads: pd.Series, genome_size: int) -> pd.Series:
    """
    Coverage = (L × N) / G
        L = read length (bp)
        N = number of mapped reads per sample
        G = genome size (bp)
    """
    return (read_length * n_reads) / genome_size


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def read_table(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    # Force numeric (some CLC exports use commas as decimal separator)
    for col in df.columns:
        df[col] = (
            df[col]
            .astype(str)
            .str.replace(",", ".", regex=False)
            .pipe(pd.to_numeric, errors="coerce")
            .fillna(0.0)
        )
    return df


def read_lengths(path: Path) -> pd.Series:
    df = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "length"])
    return df.set_index("gene_id")["length"]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Convert raw counts or FPKM to TPM.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--counts",
        required=True,
        type=Path,
        metavar="FILE",
        help="Count/FPKM matrix (TSV or CSV; genes × samples).",
    )
    p.add_argument(
        "--lengths",
        type=Path,
        default=None,
        metavar="FILE",
        help="Two-column TSV with gene_id and length in bp. Required unless --from-fpkm.",
    )
    p.add_argument(
        "--from-fpkm",
        action="store_true",
        help="Input is already in FPKM/RPKM — skip length normalisation.",
    )
    p.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="FILE",
        help="Output TPM matrix (TSV).",
    )
    p.add_argument(
        "--coverage",
        nargs=3,
        metavar=("READ_LEN", "MAPPED_READS_CSV", "GENOME_SIZE"),
        default=None,
        help=(
            "Optionally compute per-sample coverage. Provide: read length (int), "
            "path to a CSV with sample → mapped_reads, and genome size in bp."
        ),
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # -- Load counts ---------------------------------------------------------
    print(f"[INFO] Reading input: {args.counts}")
    mat = read_table(args.counts)
    print(f"[INFO] Shape: {mat.shape[0]} genes × {mat.shape[1]} samples")

    # -- Convert to TPM ------------------------------------------------------
    if args.from_fpkm:
        print("[INFO] Mode: FPKM → TPM")
        tpm = fpkm_to_tpm(mat)
    else:
        if args.lengths is None:
            sys.exit("[ERROR] --lengths is required when not using --from-fpkm.")
        print(f"[INFO] Mode: counts → TPM  (lengths: {args.lengths})")
        lengths = read_lengths(args.lengths)
        tpm = counts_to_tpm(mat, lengths)

    # -- Optional coverage ---------------------------------------------------
    if args.coverage is not None:
        read_len, mapped_path, genome_size = args.coverage
        read_len = int(read_len)
        genome_size = int(genome_size)
        mapped = pd.read_csv(mapped_path, index_col=0).squeeze("columns")
        cov = compute_coverage(read_len, mapped, genome_size)
        cov_path = args.output.with_suffix("").parent / "coverage_stats.tsv"
        cov.to_csv(cov_path, sep="\t", header=["coverage"])
        print(f"[INFO] Coverage stats saved to {cov_path}")

    # -- Save ----------------------------------------------------------------
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tpm.to_csv(args.output, sep="\t")
    print(f"[INFO] TPM matrix saved to {args.output}")


if __name__ == "__main__":
    main()
