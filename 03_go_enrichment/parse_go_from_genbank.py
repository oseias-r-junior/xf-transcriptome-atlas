"""
parse_go_from_genbank.py — Extract GO term annotations from GenBank files,
keyed by IMG ID (matching the gene_id convention used everywhere else in
this repo's DESeq2/GO/WGCNA pipeline).

Two things had to be fixed here relative to a naive "grep GO: out of the
GenBank file" approach:

1. GO terms in these files are NOT in the plain NCBI RefSeq style (a bare
   'GO:0003677' inside 'db_xref'/'note'). They are IMG-style qualifiers with
   the ID and a human-readable description packed into one string, one
   qualifier per GO aspect:

       /GO_function="GO:0003697 - single-stranded DNA binding [Evidence IEA]"
       /GO_function="GO:0005524 - ATP binding [Evidence IEA]"
       /GO_process="GO:0006260 - DNA replication [Evidence IEA]"
       /GO_component="GO:0005737 - cytoplasm [Evidence IEA]"

   A CDS can carry several qualifiers of the same kind. 'db_xref'/'note'/
   'function' are still scanned too (bare 'GO:0003677' style), in case a
   future GenBank export uses the plain RefSeq convention instead.

2. GenBank CDS features only carry old_locus_tag / locus_tag / protein_id --
   never the IMG ID (XF9a_##### / XFTem_#####) that run_go_enrichment.py's
   DEG tables and gene universe are keyed by (they come from
   raw_counts_combined.tsv, which is IMG-ID-indexed). So this script also
   takes the same --annot-9a5c/--annot-temecula1 "annotation comparator"
   tables used by build_gene_dictionary.py, and translates
   old_locus_tag -> IMG_ID (matching both the "XF0677" and "XF_0677"
   underscore conventions) before writing gene_id. A GO annotation table
   keyed by locus_tag or protein_id would silently produce ZERO matches in
   run_go_enrichment.py's enrichment_test(), since query/background gene
   sets there are always IMG IDs.

Usage
-----
python parse_go_from_genbank.py \\
    --gbk-9a5c         data/9a5c.gbff \\
    --gbk-temecula1    data/Temecula1.gbff \\
    --annot-9a5c       data/annot_comprator_9a5c.csv \\
    --annot-temecula1  data/annot_comprator_Temecula1.csv \\
    --output           results/go_annotations.tsv

Output (TSV):
    gene_id      go_term      go_description                go_category  strain
    XF9a_00003   GO:0003697   single-stranded DNA binding    function     9a5c
    XF9a_00003   GO:0006260   DNA replication                process      9a5c
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
# "GO:0003697 - single-stranded DNA binding [Evidence IEA]" -> description
GO_DESC_RE = re.compile(r"GO:\d{7}\s*-\s*([^\[]+?)\s*(?:\[|$)")

# IMG-style qualifiers, one GO aspect each.
GO_ASPECT_QUALIFIERS = {
    "GO_function":  "function",
    "GO_process":   "process",
    "GO_component": "component",
}


# ---------------------------------------------------------------------------
# old_locus_tag normalisation (shared convention with build_gene_dictionary.py)
# ---------------------------------------------------------------------------

_TAG_RE = re.compile(r"^([A-Za-z]+)_?(\d+)$")


def tag_variants(tag: str) -> set[str]:
    """'XF_0677' or 'XF0677' -> {'XF0677', 'XF_0677'}."""
    tag = (tag or "").strip()
    if not tag:
        return set()
    m = _TAG_RE.match(tag)
    if not m:
        return {tag}
    prefix, num = m.groups()
    return {f"{prefix}{num}", f"{prefix}_{num}"}


# ---------------------------------------------------------------------------
# GenBank parsing: old_locus_tag -> [(go_term, go_description, go_category), ...]
# ---------------------------------------------------------------------------

def parse_go_by_old_locus_tag(path: Path) -> dict[str, list[tuple[str, str, str]]]:
    lookup: dict[str, list[tuple[str, str, str]]] = {}
    n_cds = n_with_go = 0
    for rec in SeqIO.parse(str(path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            old_tags = q.get("old_locus_tag", [])
            if not old_tags:
                continue
            n_cds += 1

            terms: dict[str, tuple[str, str]] = {}  # go_term -> (desc, category)

            # IMG-style: /GO_function=, /GO_process=, /GO_component=
            for qualifier_name, category in GO_ASPECT_QUALIFIERS.items():
                for text in q.get(qualifier_name, []):
                    m = GO_RE.search(text)
                    if not m:
                        continue
                    go_term = f"GO:{m.group(1)}"
                    dm = GO_DESC_RE.search(text)
                    desc = dm.group(1).strip() if dm else ""
                    terms[go_term] = (desc, category)

            # Fallback: plain NCBI style (bare 'GO:xxxxxxx' in db_xref/note/function)
            for field in ("db_xref", "note", "function"):
                for text in q.get(field, []):
                    for match in GO_RE.finditer(text):
                        go_term = f"GO:{match.group(1)}"
                        terms.setdefault(go_term, ("", ""))

            if not terms:
                continue
            n_with_go += 1
            term_list = [(t, d, c) for t, (d, c) in terms.items()]
            for raw_tag in old_tags:
                for variant in tag_variants(raw_tag):
                    lookup[variant] = term_list

    print(f"[INFO] {path.name}: {n_cds} CDS with an old_locus_tag, "
          f"{n_with_go} carry >=1 GO term "
          f"({len(lookup)} lookup keys after underscore-variant expansion)")
    return lookup


# ---------------------------------------------------------------------------
# old_locus_tag -> IMG_ID via the annotation-comparator table
# ---------------------------------------------------------------------------

def load_old_tag_to_img_id(annot_path: Path) -> dict[str, str]:
    sep = "\t" if annot_path.suffix in {".tsv", ".txt"} else None
    df = pd.read_csv(annot_path, sep=sep, engine="python" if sep is None else "c")
    required = ["annot_1.1", "annot_3.1"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {annot_path.name} missing expected columns: {missing}")

    mapping: dict[str, str] = {}
    for _, r in df.iterrows():
        old_tag = str(r["annot_1.1"]).strip()
        img_id = str(r["annot_3.1"]).strip()
        if not old_tag or old_tag.lower() == "nan" or not img_id or img_id.lower() == "nan":
            continue
        for variant in tag_variants(old_tag):
            mapping[variant] = img_id
    return mapping


# ---------------------------------------------------------------------------
# Combine: GBFF GO terms (by old_locus_tag) + annot_comprator (old_locus_tag -> IMG_ID)
# ---------------------------------------------------------------------------

def build_strain_go_table(
    gbk_path: Path, annot_path: Path, strain_label: str,
) -> pd.DataFrame:
    go_lookup = parse_go_by_old_locus_tag(gbk_path)
    tag_to_img = load_old_tag_to_img_id(annot_path)

    records = []
    matched_genes = set()
    for old_tag_variant, img_id in tag_to_img.items():
        if old_tag_variant not in go_lookup:
            continue
        matched_genes.add(img_id)
        for go_term, desc, category in go_lookup[old_tag_variant]:
            records.append({
                "gene_id": img_id, "go_term": go_term,
                "go_description": desc, "go_category": category,
                "strain": strain_label,
            })

    df = pd.DataFrame(records).drop_duplicates(subset=["gene_id", "go_term"])
    if df.empty:
        print(
            f"[WARN] No GO annotations resolved to an IMG ID for {strain_label} "
            f"({gbk_path.name} x {annot_path.name}). Checked /GO_function, "
            f"/GO_process, /GO_component qualifiers (IMG style) and "
            f"db_xref/note/function (plain 'GO:0003677' style); joined via "
            f"old_locus_tag against annot_1.1.",
            file=sys.stderr,
        )
    else:
        print(f"[INFO] {strain_label}: {len(df)} gene-GO pairs "
              f"({len(matched_genes)} genes with an IMG ID and >=1 GO term, "
              f"{df['go_term'].nunique()} unique terms)")
    return df


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Extract GO annotations from GenBank files, keyed by IMG ID.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--gbk-9a5c", required=True, type=Path)
    p.add_argument("--gbk-temecula1", required=True, type=Path)
    p.add_argument("--annot-9a5c", required=True, type=Path,
                   help="annot_comprator_9a5c.csv (old_locus_tag <-> IMG_ID per gene).")
    p.add_argument("--annot-temecula1", required=True, type=Path,
                   help="annot_comprator_Temecula1.csv (old_locus_tag <-> IMG_ID per gene).")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV file (gene_id x go_term, long format).")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    df_9a5c = build_strain_go_table(args.gbk_9a5c, args.annot_9a5c, "9a5c")
    df_tem  = build_strain_go_table(args.gbk_temecula1, args.annot_temecula1, "Temecula1")

    combined = pd.concat([df_9a5c, df_tem], ignore_index=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] GO annotation table saved to {args.output} ({len(combined)} rows)")


if __name__ == "__main__":
    main()
