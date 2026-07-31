"""
build_gene_dictionary.py — Construct a cross-strain gene correspondence table
for Xylella fastidiosa strains 9a5c and Temecula1 using reciprocal best-hit
(RBH) BLASTP results and GenBank/IMG annotations.

The output dictionary is used by all downstream scripts to reconcile gene
identifiers across annotation systems (IMG ID, old locus tag, NCBI locus tag).

Usage
-----
python build_gene_dictionary.py \\
    --blastp-fwd data/blast_9a5c_vs_temecula1.tsv \\
    --blastp-rev data/blast_temecula1_vs_9a5c.tsv \\
    --gbk-9a5c  data/9a5c.gbff \\
    --gbk-tem   data/temecula1.gbff \\
    --output    results/gene_dictionary.tsv

BLASTP input format (standard tabular, -outfmt 6):
    qseqid  sseqid  pident  length  mismatch  gapopen  qstart  qend  sstart  send  evalue  bitscore

GenBank files: NCBI .gbff or IMG-exported .gbff for each strain.

Output columns
--------------
9a5c_IMG_ID         : IMG gene identifier for 9a5c
9a5c_old_locus_tags : space-separated old locus tags (e.g. XF0001 XF_0001)
9a5c_locus_tag      : NCBI locus tag
9a5c_product        : protein product annotation
Temecula1_IMG_ID    : IMG gene identifier for Temecula1
Temecula1_old_locus_tags : space-separated old locus tags (e.g. PD0001 PD_0001)
Temecula1_locus_tag : NCBI locus tag
Temecula1_product   : protein product annotation
identity_pct        : BLASTP % identity of the RBH pair
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# GenBank parsing
# ---------------------------------------------------------------------------

def parse_genbank(gbk_path: Path) -> pd.DataFrame:
    """
    Extract per-gene annotation from a GenBank file.

    Returns a DataFrame indexed by protein_id (IMG ID) with columns:
        locus_tag, old_locus_tags, product
    """
    records = []
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            protein_id = q.get("protein_id", [None])[0]
            locus_tag  = q.get("locus_tag",  [None])[0]
            old_tags   = " ".join(q.get("old_locus_tag", []))
            product    = q.get("product", ["hypothetical protein"])[0]
            if protein_id is None and locus_tag is None:
                continue
            records.append(
                dict(
                    img_id=protein_id or locus_tag,
                    locus_tag=locus_tag,
                    old_locus_tags=old_tags,
                    product=product,
                )
            )
    df = pd.DataFrame(records).drop_duplicates(subset="img_id").set_index("img_id")
    print(f"[INFO] {gbk_path.name}: {len(df)} CDS features parsed")
    return df


# ---------------------------------------------------------------------------
# RBH identification
# ---------------------------------------------------------------------------

BLAST_COLS = [
    "qseqid", "sseqid", "pident", "length",
    "mismatch", "gapopen", "qstart", "qend",
    "sstart", "send", "evalue", "bitscore",
]

def load_blast(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", names=BLAST_COLS, comment="#")
    # Best hit per query (highest bitscore)
    best = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid")
    return best.set_index("qseqid")[["sseqid", "pident", "evalue"]]


def reciprocal_best_hits(
    fwd: pd.DataFrame,   # 9a5c → Tem1
    rev: pd.DataFrame,   # Tem1 → 9a5c
    evalue_cutoff: float = 1e-5,
    identity_cutoff: float = 30.0,
) -> pd.DataFrame:
    """Return gene pairs that are each other's best BLASTP hit."""
    fwd_f = fwd[(fwd["evalue"] <= evalue_cutoff) & (fwd["pident"] >= identity_cutoff)]
    rev_f = rev[(rev["evalue"] <= evalue_cutoff) & (rev["pident"] >= identity_cutoff)]

    rbh_rows = []
    for q9a, row in fwd_f.iterrows():
        tem_hit = row["sseqid"]
        if tem_hit in rev_f.index and rev_f.loc[tem_hit, "sseqid"] == q9a:
            rbh_rows.append(
                dict(
                    img_9a5c=q9a,
                    img_tem=tem_hit,
                    identity_pct=row["pident"],
                )
            )

    df_rbh = pd.DataFrame(rbh_rows)
    print(f"[INFO] {len(df_rbh)} reciprocal best-hit pairs found")
    return df_rbh


# ---------------------------------------------------------------------------
# Main builder
# ---------------------------------------------------------------------------

def build_dictionary(
    rbh: pd.DataFrame,
    ann9a: pd.DataFrame,
    anntm: pd.DataFrame,
) -> pd.DataFrame:
    rows = []
    for _, pair in rbh.iterrows():
        r9 = ann9a.loc[pair["img_9a5c"]] if pair["img_9a5c"] in ann9a.index else pd.Series()
        rt = anntm.loc[pair["img_tem"]]  if pair["img_tem"]  in anntm.index  else pd.Series()
        rows.append(
            {
                "9a5c_IMG_ID":              pair["img_9a5c"],
                "9a5c_old_locus_tags":      r9.get("old_locus_tags", ""),
                "9a5c_locus_tag":           r9.get("locus_tag", ""),
                "9a5c_product":             r9.get("product", ""),
                "Temecula1_IMG_ID":         pair["img_tem"],
                "Temecula1_old_locus_tags": rt.get("old_locus_tags", ""),
                "Temecula1_locus_tag":      rt.get("locus_tag", ""),
                "Temecula1_product":        rt.get("product", ""),
                "identity_pct":             pair["identity_pct"],
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build cross-strain gene correspondence dictionary via RBH BLASTP.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--blastp-fwd", required=True, type=Path,
                   help="BLASTP results: 9a5c (query) vs Temecula1 (subject), -outfmt 6.")
    p.add_argument("--blastp-rev", required=True, type=Path,
                   help="BLASTP results: Temecula1 (query) vs 9a5c (subject), -outfmt 6.")
    p.add_argument("--gbk-9a5c",  required=True, type=Path,
                   help="GenBank file for strain 9a5c.")
    p.add_argument("--gbk-tem",   required=True, type=Path,
                   help="GenBank file for strain Temecula1.")
    p.add_argument("--output",    required=True, type=Path,
                   help="Output TSV file (gene_dictionary.tsv).")
    p.add_argument("--evalue",    type=float, default=1e-5,
                   help="E-value cutoff for BLASTP hits.")
    p.add_argument("--identity",  type=float, default=30.0,
                   help="Minimum %% identity for BLASTP hits.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    ann9a = parse_genbank(args.gbk_9a5c)
    anntm = parse_genbank(args.gbk_tem)

    fwd = load_blast(args.blastp_fwd)
    rev = load_blast(args.blastp_rev)

    rbh = reciprocal_best_hits(fwd, rev,
                                evalue_cutoff=args.evalue,
                                identity_cutoff=args.identity)

    dictionary = build_dictionary(rbh, ann9a, anntm)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    dictionary.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Gene dictionary saved to {args.output}  ({len(dictionary)} pairs)")


if __name__ == "__main__":
    main()
