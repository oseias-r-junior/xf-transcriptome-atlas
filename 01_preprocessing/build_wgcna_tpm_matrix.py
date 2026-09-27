"""
build_wgcna_tpm_matrix.py — Build the cross-strain, ortholog-merged TPM matrix
that 04_wgcna/run_wgcna.py (and the other 04_wgcna/* scripts) actually need.

Why this script exists
-----------------------
data/tpm_expression.csv, as produced by build_matrices_from_fpkm_xlsx.py +
compute_tpm.py, is genes x samples where "genes" is the UNION of each
strain's own native gene IDs (XF9a_#####, XFTem_#####, pXF51_#####, rRNA
names like "16S_1", ...). Every gene row is therefore only ever populated
for ONE strain's 12 columns -- the other strain's 12 columns are hard
structural zeros, not biological zeros. Running WGCNA on that matrix across
all 24 samples at once is wrong: co-expression would partly reflect "which
strain has data for this gene" rather than real biology, and the original
analysis (see WGCNA_paper_pierry_feitosa.ipynb) never did this.

The original notebook instead reads an already-merged tpm_expression.csv
whose gene_id is "<9a5c_IMG_ID>_<Temecula1_IMG_ID>" (e.g.
"XF9a_00002_XFTem_00002") -- one row per RBH ortholog PAIR, with real TPM
values for both strains on every row. That pairing is exactly what
data/gene_dictionary.tsv already computes (build_gene_dictionary.py), so
this script just needs to re-key data/tpm_expression.csv through it:

    for each ortholog pair (9a5c_IMG_ID, Temecula1_IMG_ID) in gene_dictionary.tsv:
        new_gene_id = f"{9a5c_IMG_ID}_{Temecula1_IMG_ID}"
        9a5c columns   <- tpm_expression.csv.loc[9a5c_IMG_ID,      9a5c_sample_cols]
        Temecula1 cols <- tpm_expression.csv.loc[Temecula1_IMG_ID, Temecula1_sample_cols]

giving a matrix with one row per ortholog pair (~1600-1700 genes, matching
the original notebook's scale) where every value is a real measurement,
not a structural zero.

Usage
-----
python 01_preprocessing/build_wgcna_tpm_matrix.py \\
    --tpm        data/tpm_expression.csv \\
    --dictionary data/gene_dictionary.tsv \\
    --output     data/tpm_expression_wgcna.csv

Output: CSV (comma-separated, to match the original notebook's own
tpm_expression.csv exactly), genes(pairs) x 24 samples, gene_id =
"<9a5c_IMG_ID>_<Temecula1_IMG_ID>".
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def sniff_sep(path: Path) -> str:
    """Detect the field separator from the file's own header line.

    File extension alone is unreliable: compute_tpm.py writes tab-separated
    data into a file named tpm_expression.csv. See the same helper (and the
    bug it fixes) in 04_wgcna/run_wgcna.py.
    """
    with open(path, "r", newline="") as f:
        header = f.readline()
    counts = {sep: header.count(sep) for sep in ("\t", ",", ";")}
    best_sep, best_count = max(counts.items(), key=lambda kv: kv[1])
    return best_sep if best_count > 0 else ","


def build_merged_matrix(
    tpm: pd.DataFrame,
    gene_dict: pd.DataFrame,
    strain_9a5c_prefix: str = "9a5c",
    strain_tem_prefix: str = "Tem1",
) -> pd.DataFrame:
    cols_9a5c = [c for c in tpm.columns if c.startswith(strain_9a5c_prefix)]
    cols_tem  = [c for c in tpm.columns if c.startswith(strain_tem_prefix)]
    if not cols_9a5c or not cols_tem:
        sys.exit(
            f"[ERROR] Could not find sample columns for both strains "
            f"(9a5c prefix={strain_9a5c_prefix!r}: {len(cols_9a5c)} cols, "
            f"Temecula1 prefix={strain_tem_prefix!r}: {len(cols_tem)} cols). "
            f"Check --strain-9a5c-prefix/--strain-tem-prefix against the "
            f"actual column names in --tpm."
        )

    rows = []
    gene_ids = []
    n_missing_9a5c = n_missing_tem = 0

    for _, r in gene_dict.iterrows():
        id_9a5c = str(r["9a5c_IMG_ID"]).strip()
        id_tem  = str(r["Temecula1_IMG_ID"]).strip()
        if not id_9a5c or id_9a5c.lower() == "nan" or not id_tem or id_tem.lower() == "nan":
            continue
        if id_9a5c not in tpm.index:
            n_missing_9a5c += 1
            continue
        if id_tem not in tpm.index:
            n_missing_tem += 1
            continue

        new_id = f"{id_9a5c}_{id_tem}"
        vals_9a5c = tpm.loc[id_9a5c, cols_9a5c]
        vals_tem  = tpm.loc[id_tem, cols_tem]
        # A gene_id can appear on >1 dictionary row (paralogs sharing an
        # ortholog partner via RBH quirks); loc[] on a duplicated index
        # would return a DataFrame instead of a Series. Guard against that
        # by taking the first match rather than silently erroring below.
        if isinstance(vals_9a5c, pd.DataFrame):
            vals_9a5c = vals_9a5c.iloc[0]
        if isinstance(vals_tem, pd.DataFrame):
            vals_tem = vals_tem.iloc[0]

        row = pd.concat([vals_9a5c, vals_tem])
        rows.append(row)
        gene_ids.append(new_id)

    if n_missing_9a5c or n_missing_tem:
        print(
            f"[WARN] {n_missing_9a5c} dictionary rows had a 9a5c_IMG_ID not "
            f"found in --tpm; {n_missing_tem} had a Temecula1_IMG_ID not "
            f"found in --tpm. These pairs were skipped.",
            file=sys.stderr,
        )

    merged = pd.DataFrame(rows, index=pd.Index(gene_ids, name="gene_id"))
    # Restore the original sample column order (9a5c cols then Tem1 cols,
    # each in the order they appeared in --tpm).
    merged = merged[cols_9a5c + cols_tem]

    n_dupe = merged.index.duplicated().sum()
    if n_dupe:
        print(
            f"[WARN] {n_dupe} duplicate ortholog-pair gene_ids after merge "
            f"(one 9a5c/Temecula1 gene resolved via >1 dictionary row); kept "
            f"all as separate rows -- inspect data/gene_dictionary.tsv if "
            f"this count is unexpectedly large.",
            file=sys.stderr,
        )

    return merged


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build an ortholog-merged, cross-strain TPM matrix for WGCNA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--tpm", required=True, type=Path,
                   help="Strain-native TPM matrix (data/tpm_expression.csv from compute_tpm.py).")
    p.add_argument("--dictionary", required=True, type=Path,
                   help="data/gene_dictionary.tsv from build_gene_dictionary.py.")
    p.add_argument("--strain-9a5c-prefix", default="9a5c")
    p.add_argument("--strain-tem-prefix", default="Tem1")
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    tpm = pd.read_csv(args.tpm, index_col=0, sep=sniff_sep(args.tpm))
    gene_dict = pd.read_csv(args.dictionary, sep="\t")

    print(f"[INFO] Loaded TPM matrix: {tpm.shape[0]} genes x {tpm.shape[1]} samples")
    print(f"[INFO] Loaded gene dictionary: {len(gene_dict)} ortholog pairs")

    merged = build_merged_matrix(
        tpm, gene_dict, args.strain_9a5c_prefix, args.strain_tem_prefix
    )
    print(f"[INFO] Merged matrix: {merged.shape[0]} ortholog-pair genes x {merged.shape[1]} samples")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output)
    print(f"[INFO] WGCNA-ready TPM matrix saved to {args.output}")


if __name__ == "__main__":
    main()
