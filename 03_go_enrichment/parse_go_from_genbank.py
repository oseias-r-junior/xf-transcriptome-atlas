"""
parse_go_from_genbank.py — Extract GO term annotations from NCBI GenBank files.

GO terms are read from the 'db_xref' qualifier (format: GeneID:XXXX) and
matched against GO annotations available in the GenBank feature qualifiers.
For Xylella fastidiosa, GO terms are embedded in the 'note' or 'db_xref' fields
of CDS features in GenBank/RefSeq records.

Usage
-----
python parse_go_from_genbank.py \\
    --gbk  data/9a5c.gbff data/temecula1.gbff \\
    --output results/go_annotations.tsv

Output (TSV):
    gene_id   go_term   go_description   strain
    Xf9a_00001   GO:0003677   DNA binding   9a5c
    ...

If a gene has multiple GO terms it appears on multiple rows.
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


GO_RE = re.compile(r"GO:(\d{7})")


def parse_go_from_gbk(path: Path, strain_label: str | None = None) -> pd.DataFrame:
    """
    Parse all CDS features in a GenBank file and return a long-format
    DataFrame of (gene_id, go_term) pairs.

    GO terms are extracted from:
      - 'db_xref' qualifiers containing 'GO:' prefixes
      - 'note' qualifiers containing GO term IDs
    """
    records = []
    for rec in SeqIO.parse(str(path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            protein_id = q.get("protein_id", [None])[0]
            locus_tag  = q.get("locus_tag",  [None])[0]
            gene_id = protein_id or locus_tag
            if gene_id is None:
                continue

            # Collect raw text to mine GO IDs from
            text_pool = []
            for field in ("db_xref", "note", "function"):
                text_pool.extend(q.get(field, []))

            go_terms = set()
            for text in text_pool:
                for match in GO_RE.finditer(text):
                    go_terms.add(f"GO:{match.group(1)}")

            for go in go_terms:
                records.append(
                    {
                        "gene_id":  gene_id,
                        "go_term":  go,
                        "strain":   strain_label or path.stem,
                    }
                )

    df = pd.DataFrame(records)
    if df.empty:
        print(
            f"[WARN] No GO annotations found in {path.name}. "
            "Check that the GenBank file contains db_xref or note fields "
            "with GO IDs (e.g. 'GO:0003677').",
            file=sys.stderr,
        )
    else:
        print(f"[INFO] {path.name}: {len(df)} gene–GO pairs "
              f"({df['gene_id'].nunique()} genes, {df['go_term'].nunique()} unique terms)")
    return df


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Extract GO annotations from GenBank files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "--gbk",
        nargs="+",
        required=True,
        type=Path,
        metavar="FILE",
        help="One or more GenBank files (.gbff / .gb).",
    )
    p.add_argument(
        "--strain-labels",
        nargs="*",
        default=None,
        metavar="LABEL",
        help=(
            "Optional labels for each --gbk file (same order). "
            "Defaults to file stem."
        ),
    )
    p.add_argument(
        "--output",
        required=True,
        type=Path,
        metavar="FILE",
        help="Output TSV file (gene_id × go_term, long format).",
    )
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    labels = args.strain_labels or [None] * len(args.gbk)
    if len(labels) != len(args.gbk):
        sys.exit("[ERROR] --strain-labels must match the number of --gbk files.")

    dfs = []
    for path, label in zip(args.gbk, labels):
        dfs.append(parse_go_from_gbk(path, label))

    combined = pd.concat(dfs, ignore_index=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] GO annotation table saved to {args.output} ({len(combined)} rows)")


if __name__ == "__main__":
    main()
