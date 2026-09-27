"""
build_matrices_from_fpkm_xlsx.py — Ingest the per-sample CLC Genomics
Workbench "Gene Expression" (GE) track exports (.xlsx) into
data/raw_counts_combined.tsv (+ data/gene_lengths.tsv), ready for
compute_tpm.py.

Background
----------
The raw per-sample data for this project is NOT a flat directory of .txt
count files (an earlier assumption in this repo, since corrected). It is a
directory tree of CLC-exported .xlsx files, one per sample, e.g.:

    FPKM/<strain>/<medium>/<strain>-<medium>-<timepoint>-<replicate>_..._GE.xlsx
    FPKM/<strain>/<medium>/pXF51/<...>_GE_pXF51.xlsx   (9a5c only: plasmid genes)

Each chromosome-level file has two sheets; the one actually carrying
computed TPM values (and a "Soma de FPKM" footer row) is named **"Bruto"**
— use that one, not the other (which lacks the TPM column). Each pXF51
file has a single sheet with the same layout as "Bruto". Columns used:
    Name, Gene Product Name, Gene length, RPKM, Unique gene reads,
    Total gene reads, TPM

Why this script only outputs counts + gene lengths (not TPM directly)
-----------------------------------------------------------------------
The chromosome and pXF51 sheets were exported as two SEPARATE CLC
expression tracks, each normalised (RPKM/TPM) against its own total
(e.g. "Soma de FPKM" = 1,046,570 for one chromosome sample vs 915,326 for
its pXF51 counterpart in the same sample). Concatenating their
already-computed TPM/RPKM columns directly would silently combine two
values normalised against two different denominators -- not a valid
genome-wide TPM. Raw read counts have no such issue (they're just counts),
so this script unions chromosome + pXF51 raw counts and gene lengths per
sample, and leaves TPM computation to `compute_tpm.py`, which recomputes
it once, consistently, over the full combined gene set per sample.

Filename parsing
-----------------
Sample identity is parsed from the filename, tolerant of the inconsistent
separators/parentheses CLC produces (e.g. both
"9a5c-PIM6-1d-3_..._paired_GE.xlsx" and
"9a5c-PIM6-1h-2_... (paired) (GE).xlsx"). Files with no <timepoint>-<rep>
suffix (e.g. "9a5c-PIM6_S2_..._GE.xlsx") are per-medium CLC summary
exports, not single-sample data -- they are skipped, with a warning.

A number of files use "<N>h" (hour) instead of "<N>d" (day) in the
timepoint field (e.g. "9a5c-PIM6-1h-2"), which does not match any
condition in sample_info.csv/comparisons.tsv (only 1d/3d/10d/7d exist).
Per project decision, every "<N>h" is treated as a typo for "<N>d" and
normalised accordingly; every such normalisation is logged so it can be
audited against the lab notebook.

Usage
-----
python build_matrices_from_fpkm_xlsx.py \\
    --fpkm-dir          FPKM/ \\
    --counts-column     "Unique gene reads" \\
    --counts-output     data/raw_counts_combined.tsv \\
    --lengths-output     data/gene_lengths.tsv \\
    --manifest-output   data/fpkm_ingest_manifest.csv

Then, to get TPM:
python compute_tpm.py \\
    --counts data/raw_counts_combined.tsv \\
    --lengths data/gene_lengths.tsv \\
    --output data/tpm_expression.csv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

try:
    import openpyxl
except ImportError:
    sys.exit("[ERROR] openpyxl is required (`pip install openpyxl`).")


# ---------------------------------------------------------------------------
# Filename parsing
# ---------------------------------------------------------------------------

FNAME_RE = re.compile(
    r"(?i)^(?P<strain>9a5c|tem1?|temecula1)"
    r"[-_](?P<medium>PIM6|PWG)"
    r"[-_](?P<timepoint>\d+)(?P<tunit>[dh])"
    r"[-_](?P<replicate>\d+)"
)

STRAIN_CANON = {  # raw filename prefix (lowered) -> (condition_prefix, sample_info strain value)
    "9a5c": ("9a5c", "9a5c"),
    "tem1": ("Tem1", "Temecula1"),
    "tem": ("Tem1", "Temecula1"),
    "temecula1": ("Tem1", "Temecula1"),
}


def parse_filename(path: Path) -> dict | None:
    """Returns None (and the caller should log a skip) for files that don't
    match the <strain>-<medium>-<timepoint>-<replicate> pattern -- these are
    CLC per-medium summary exports, not per-sample data."""
    stem = path.stem
    if path.name.startswith("~$"):
        return None  # Excel lock file
    m = FNAME_RE.match(stem)
    if not m:
        return None

    g = m.groupdict()
    cond_prefix, strain = STRAIN_CANON[g["strain"].lower()]
    medium = g["medium"].upper()
    timepoint = int(g["timepoint"])
    tunit = g["tunit"].lower()
    replicate = int(g["replicate"])

    normalized_from_hour = False
    if tunit == "h":
        # Project decision: every "<N>h" file found is a data-entry typo
        # for "<N>d" (no condition in sample_info.csv/comparisons.tsv uses
        # an hour unit). Treated as such here; flagged in the manifest.
        normalized_from_hour = True
        tunit = "d"

    condition = f"{cond_prefix}_{medium}_{timepoint}{tunit}"
    sample_id = f"{condition}_rep{replicate}"

    is_plasmid = "pxf51" in stem.lower() or "pxf51" in path.parent.name.lower()

    return dict(
        sample_id=sample_id, condition=condition, strain=strain, medium=medium,
        timepoint=f"{timepoint}{tunit}", replicate=replicate,
        is_plasmid=is_plasmid, normalized_from_hour=normalized_from_hour,
        source_file=str(path),
    )


# ---------------------------------------------------------------------------
# XLSX reading
# ---------------------------------------------------------------------------

def read_ge_sheet(path: Path, is_plasmid: bool) -> pd.DataFrame:
    """Read the sheet carrying real TPM values: 'Bruto' for chromosome
    files, the (only) sheet present for pXF51 files. Drops the trailing
    blank row and the 'Soma de FPKM' footer row."""
    wb = openpyxl.load_workbook(str(path), read_only=True, data_only=True)
    if is_plasmid:
        ws = wb[wb.sheetnames[0]]
    else:
        if "Bruto" not in wb.sheetnames:
            raise RuntimeError(f"{path.name}: no 'Bruto' sheet found (sheets: {wb.sheetnames})")
        ws = wb["Bruto"]

    rows = list(ws.iter_rows(values_only=True))
    header = [str(h).strip() if h is not None else "" for h in rows[0]]
    df = pd.DataFrame(rows[1:], columns=header)
    df = df[df["Name"].notna()]
    wb.close()
    return df


# ---------------------------------------------------------------------------
# Main ingest
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build raw_counts_combined.tsv + gene_lengths.tsv from "
                     "the CLC .xlsx FPKM export tree.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--fpkm-dir", required=True, type=Path,
                   help="Root of the FPKM/<strain>/<medium>/*.xlsx (+pXF51/) tree.")
    p.add_argument("--counts-column", default="Unique gene reads",
                   choices=["Unique gene reads", "Total gene reads"],
                   help="Which column to treat as the raw count for DESeq2.")
    p.add_argument("--include-plasmid", dest="include_plasmid",
                   action="store_true", default=True,
                   help="Union pXF51 gene rows into the 9a5c samples that have them.")
    p.add_argument("--exclude-plasmid", dest="include_plasmid", action="store_false")
    p.add_argument("--counts-output", required=True, type=Path)
    p.add_argument("--lengths-output", required=True, type=Path)
    p.add_argument("--manifest-output", type=Path, default=None,
                   help="Optional CSV logging every file's parsed sample_id/"
                        "condition/strain/medium/timepoint/replicate, plus "
                        "skip/normalisation notes -- cross-check this "
                        "against data/sample_info.csv.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    all_files = sorted(args.fpkm_dir.rglob("*.xlsx"))
    print(f"[INFO] Found {len(all_files)} .xlsx files under {args.fpkm_dir}")

    manifest_rows = []
    per_sample_counts: dict[str, list[pd.DataFrame]] = {}
    per_sample_lengths: dict[str, pd.DataFrame] = {}

    for path in all_files:
        meta = parse_filename(path)
        if meta is None:
            if not path.name.startswith("~$"):
                print(f"[SKIP] {path.relative_to(args.fpkm_dir)}: doesn't match "
                      f"<strain>-<medium>-<timepoint>-<replicate> (likely a "
                      f"per-medium CLC summary export, not per-sample data)")
                manifest_rows.append({"source_file": str(path), "status": "SKIPPED (no filename match)"})
            continue

        if meta["is_plasmid"] and not args.include_plasmid:
            manifest_rows.append({**meta, "status": "SKIPPED (--exclude-plasmid)"})
            continue

        try:
            df = read_ge_sheet(path, meta["is_plasmid"])
        except Exception as e:
            print(f"[WARN] {path.name}: failed to read ({e}); skipping.", file=sys.stderr)
            manifest_rows.append({**meta, "status": f"ERROR: {e}"})
            continue

        counts = df[["Name", args.counts_column]].rename(
            columns={"Name": "gene_id", args.counts_column: meta["sample_id"]}
        ).set_index("gene_id")
        lengths = df[["Name", "Gene length"]].rename(
            columns={"Name": "gene_id", "Gene length": "length"}
        ).set_index("gene_id")["length"]

        per_sample_counts.setdefault(meta["sample_id"], []).append(counts)
        if meta["sample_id"] not in per_sample_lengths:
            per_sample_lengths[meta["sample_id"]] = lengths
        else:
            per_sample_lengths[meta["sample_id"]] = pd.concat(
                [per_sample_lengths[meta["sample_id"]], lengths]
            ).groupby(level=0).first()

        note = " [1h->1d]" if meta["normalized_from_hour"] else ""
        print(f"[OK]   {path.relative_to(args.fpkm_dir)} -> sample_id={meta['sample_id']} "
              f"({'plasmid' if meta['is_plasmid'] else 'chromosome'}, {len(df)} genes){note}")
        manifest_rows.append({**meta, "status": "OK", "n_genes": len(df)})

    if not per_sample_counts:
        sys.exit("[ERROR] No usable samples parsed -- check --fpkm-dir and filename patterns.")

    # -- Union chromosome + plasmid rows per sample, then pivot to a wide matrix
    sample_series = {}
    for sample_id, dfs in per_sample_counts.items():
        merged = pd.concat(dfs)
        dupes = merged.index[merged.index.duplicated()]
        if len(dupes) > 0:
            print(f"[WARN] {sample_id}: {len(dupes)} gene_id(s) appear in both "
                  f"chromosome and plasmid exports (unexpected overlap): "
                  f"{list(dupes[:5])}{'...' if len(dupes) > 5 else ''}. "
                  f"Keeping the first occurrence.")
            merged = merged[~merged.index.duplicated(keep="first")]
        sample_series[sample_id] = merged[sample_id]

    counts_df = pd.DataFrame(sample_series).fillna(0).astype(int)
    counts_df.index.name = "gene_id"

    lengths_combined = pd.concat(per_sample_lengths.values())
    lengths_combined = lengths_combined.groupby(level=0).first()
    lengths_combined = lengths_combined.reindex(counts_df.index)
    missing_len = lengths_combined.isna().sum()
    if missing_len:
        print(f"[WARN] {missing_len} genes have no length recovered; "
              f"they will get TPM=0 downstream in compute_tpm.py.")

    args.counts_output.parent.mkdir(parents=True, exist_ok=True)
    counts_df.to_csv(args.counts_output, sep="\t")
    print(f"[INFO] Saved {counts_df.shape[0]} genes x {counts_df.shape[1]} samples "
          f"to {args.counts_output}")

    args.lengths_output.parent.mkdir(parents=True, exist_ok=True)
    lengths_combined.dropna().astype(int).to_csv(
        args.lengths_output, sep="\t", header=False,
    )
    print(f"[INFO] Saved gene lengths to {args.lengths_output}")

    if args.manifest_output:
        pd.DataFrame(manifest_rows).to_csv(args.manifest_output, index=False)
        print(f"[INFO] Manifest saved to {args.manifest_output} -- "
              f"cross-check sample counts per condition against data/sample_info.csv.")

    print(f"\n[INFO] Samples per condition (verify replicate counts against "
          f"data/sample_info.csv):")
    cond_counts = pd.Series([m["condition"] for m in manifest_rows if m.get("status") == "OK"
                              and not m.get("is_plasmid")]).value_counts().sort_index()
    print(cond_counts.to_string())

    print(f"\n[INFO] Next step:\n"
          f"  python compute_tpm.py --counts {args.counts_output} "
          f"--lengths {args.lengths_output} --output data/tpm_expression.csv")


if __name__ == "__main__":
    main()
