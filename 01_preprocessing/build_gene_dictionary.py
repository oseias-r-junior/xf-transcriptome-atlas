"""
build_gene_dictionary.py — Construct the cross-strain gene correspondence
table (gene_dictionary.tsv / Supplementary Table S4) for Xylella fastidiosa
strains 9a5c and Temecula1.

This does NOT run BLASTP itself. Reciprocal-best-hit (RBH) identification
between the two strains' proteomes is a separate, already-completed step —
supply its result via --rbh. What this script does is join three
independent per-gene annotation sources into one cross-strain table:

  1. --rbh                 : precomputed RBH pairs, protein_ID vs protein_ID
                              (e.g. "WP_004083467.1,WP_004083467.1" — one row
                              per orthologous gene pair). This is the ONLY
                              cross-strain link; everything else below is
                              purely per-strain.
  2. --annot-9a5c/--annot-temecula1 : per-strain "annotation comparator"
                              tables that already link, for each gene, three
                              parallel ID systems used across the project:
                                annot_1.* = old locus tag (XF_####/PD_####)
                                            + product + length
                                annot_2.* = secondary NCBI-style numbering
                                            (kept as a reference column;
                                            not used elsewhere downstream)
                                annot_3.* = the IMG ID (XF9a_#####/XFTem_#####)
                                            used throughout this repository
                                            + IMG product + IMG length
                              These tables do NOT carry protein_id, so they
                              cannot be joined directly to --rbh.
  3. --gbk-9a5c/--gbk-tem   : the NCBI GenBank flat files, parsed only to
                              recover, per old_locus_tag, the protein_id and
                              GenBank product/gene/locus_tag -- this is the
                              bridge that lets the RBH pairs (protein_id)
                              reach the IMG IDs (from step 2).

Join logic, per strain:
    gbff CDS features (old_locus_tag -> protein_id, locus_tag, product, gene)
        left-joined onto
    annot_comprator_<strain>.csv (old_locus_tag -> IMG_ID, IMG product/length)
        on old_locus_tag (matched allowing both "XF0677" and "XF_0677" forms,
        since the two annotation sources don't always agree on the
        underscore convention).

The two resulting per-strain tables (now both carrying protein_id) are then
merged on the RBH pairs to produce one row per orthologous gene pair.

Usage
-----
python build_gene_dictionary.py \\
    --rbh              data/9a5c_vs_Temecula1_rbh.csv \\
    --annot-9a5c       data/annot_comprator_9a5c.csv \\
    --annot-temecula1  data/annot_comprator_Temecula1.csv \\
    --gbk-9a5c         data/9a5c.gbff \\
    --gbk-tem          data/temecula1.gbff \\
    --output           data/gene_dictionary.tsv

Input formats
-------------
--rbh (CSV, header): 9a5c_protein_ID,Temecula1_protein_ID
    WP_004083467.1,WP_004083467.1

--annot-9a5c / --annot-temecula1 (TSV, header): one row per gene, with three
parallel annotation blocks (see docstring above):
    gene    annot_1.1  annot_1.2  annot_1.3  annot_1.4  annot_2.1  ...  annot_3.4
    gene_1  XF_0001    CDS        <product>  1320       9a5c_00001 ...  1332

--gbk-9a5c / --gbk-tem: standard NCBI GenBank flat files (.gbff).

Output columns (matches Supplementary Table S4 / the dictionary format
already consumed by every downstream script in this repo):
    9a5c_IMG_ID, 9a5c_protein_ID, 9a5c_locus_tag, 9a5c_old_locus_tags,
    9a5c_gene, 9a5c_product_GBFF, 9a5c_product_IMG, 9a5c_IMG_length,
    Temecula1_protein_ID, Temecula1_locus_tag, Temecula1_old_locus_tags,
    Temecula1_gene, Temecula1_product_GBFF, Temecula1_IMG_ID,
    Temecula1_product_IMG, Temecula1_IMG_length
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO


# ---------------------------------------------------------------------------
# Old-locus-tag normalisation (the two annotation sources don't always agree
# on "XF0677" vs "XF_0677")
# ---------------------------------------------------------------------------

_TAG_RE = re.compile(r"^([A-Za-z]+)_?(\d+)$")


def tag_variants(tag: str) -> set[str]:
    """Return both the underscored and non-underscored form of a locus tag,
    e.g. 'XF_0677' or 'XF0677' -> {'XF0677', 'XF_0677'}."""
    tag = (tag or "").strip()
    if not tag:
        return set()
    m = _TAG_RE.match(tag)
    if not m:
        return {tag}
    prefix, num = m.groups()
    return {f"{prefix}{num}", f"{prefix}_{num}"}


# ---------------------------------------------------------------------------
# 1. GenBank parsing: old_locus_tag -> protein_id, locus_tag, product, gene
# ---------------------------------------------------------------------------

def parse_genbank_by_old_locus_tag(gbk_path: Path) -> dict[str, dict]:
    """Returns {old_locus_tag_variant: {protein_id, locus_tag, product, gene}},
    with every CDS keyed under BOTH underscore variants of each of its
    old_locus_tag qualifier values, so lookups from either annotation
    convention succeed."""
    lookup: dict[str, dict] = {}
    n_cds = 0
    for rec in SeqIO.parse(str(gbk_path), "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            q = feat.qualifiers
            protein_id = q.get("protein_id", [None])[0]
            locus_tag  = q.get("locus_tag", [None])[0]
            old_tags   = q.get("old_locus_tag", [])
            product    = q.get("product", ["hypothetical protein"])[0]
            gene       = q.get("gene", [""])[0]
            if not old_tags:
                continue
            n_cds += 1
            record = dict(protein_id=protein_id, locus_tag=locus_tag,
                          product=product, gene=gene)
            for raw_tag in old_tags:
                for variant in tag_variants(raw_tag):
                    lookup[variant] = record
    print(f"[INFO] {gbk_path.name}: {n_cds} CDS features with an old_locus_tag "
          f"({len(lookup)} lookup keys after underscore-variant expansion)")
    return lookup


# ---------------------------------------------------------------------------
# 2. Annotation-comparator table -> per-strain table with protein_id attached
# ---------------------------------------------------------------------------

def load_annot_comparator(path: Path) -> pd.DataFrame:
    sep = "\t" if path.suffix in {".tsv", ".txt"} else None
    df = pd.read_csv(path, sep=sep, engine="python" if sep is None else "c")
    required = ["annot_1.1", "annot_1.3", "annot_1.4",
                "annot_3.1", "annot_3.3", "annot_3.4"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        sys.exit(f"[ERROR] {path.name} missing expected columns: {missing}")
    return df


def build_strain_table(annot_path: Path, gbk_lookup: dict[str, dict], strain_label: str) -> pd.DataFrame:
    annot = load_annot_comparator(annot_path)

    rows = []
    n_matched = 0
    for _, r in annot.iterrows():
        old_tag = str(r["annot_1.1"]).strip()
        img_id  = str(r["annot_3.1"]).strip()
        if not img_id or img_id.lower() == "nan":
            continue  # no IMG ID for this gene, unusable downstream

        gbk_rec = None
        for variant in tag_variants(old_tag):
            if variant in gbk_lookup:
                gbk_rec = gbk_lookup[variant]
                break
        if gbk_rec is not None:
            n_matched += 1

        rows.append({
            "IMG_ID":         img_id,
            "protein_ID":     (gbk_rec or {}).get("protein_id") or "",
            "locus_tag":      (gbk_rec or {}).get("locus_tag") or "",
            "old_locus_tags": " ".join(sorted(tag_variants(old_tag))) if old_tag and old_tag.lower() != "nan" else "",
            "gene":           (gbk_rec or {}).get("gene") or "",
            "product_GBFF":   (gbk_rec or {}).get("product") or "",
            "product_IMG":    str(r["annot_3.3"]) if pd.notna(r["annot_3.3"]) else "",
            "IMG_length":     r["annot_3.4"] if pd.notna(r["annot_3.4"]) else "",
        })

    table = pd.DataFrame(rows).drop_duplicates(subset="IMG_ID")

    # A handful of RefSeq protein_id values are assigned to more than one
    # genomic locus within the SAME strain (e.g. two CDS with byte-identical
    # translations, such as duplicated/paralogous genes). protein_ID is the
    # only key available to join to the RBH pairs, so a locus sharing its
    # protein_id with another locus cannot be disambiguated here -- joining
    # on it naively would cross-product every (9a5c locus, Temecula1 locus)
    # combination that happens to share that protein_id, silently injecting
    # spurious ortholog pairs. Such loci are excluded from the RBH join
    # (protein_ID blanked, same treatment as "no GenBank match") rather than
    # guessed at; they still appear in `table` un-joined.
    dup_mask = (table["protein_ID"] != "") & table["protein_ID"].duplicated(keep=False)
    n_ambiguous = dup_mask.sum()
    if n_ambiguous:
        print(f"[WARN] {strain_label}: {n_ambiguous} genes share a protein_id with "
              f"another locus in the same strain ({table.loc[dup_mask, 'protein_ID'].nunique()} "
              f"distinct protein_ids affected) -- excluded from the RBH join to avoid "
              f"cross-product false pairs; kept in the per-strain table with protein_ID blanked.")
        table.loc[dup_mask, "protein_ID"] = ""

    print(f"[INFO] {strain_label}: {len(table)} genes with an IMG ID; "
          f"{n_matched}/{len(rows)} matched to a GenBank protein_id via old_locus_tag "
          f"({len(rows) - n_matched} have no protein_id and therefore cannot enter "
          f"the RBH join below; of the matched ones, {n_ambiguous} were ambiguous "
          f"multi-locus protein_ids and were also excluded from the join)")
    return table


# ---------------------------------------------------------------------------
# 3. Join via precomputed RBH pairs
# ---------------------------------------------------------------------------

def load_rbh(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    col_9a = next((cols[c] for c in cols if "9a5c" in c), None)
    col_tem = next((cols[c] for c in cols if "tem" in c), None)
    if col_9a is None or col_tem is None:
        sys.exit(f"[ERROR] Could not identify 9a5c/Temecula1 protein_ID columns in {path}. "
                  f"Found columns: {list(df.columns)}")
    df = df.rename(columns={col_9a: "9a5c_protein_ID", col_tem: "Temecula1_protein_ID"})
    df = df.dropna(subset=["9a5c_protein_ID", "Temecula1_protein_ID"])
    df["9a5c_protein_ID"] = df["9a5c_protein_ID"].astype(str).str.strip()
    df["Temecula1_protein_ID"] = df["Temecula1_protein_ID"].astype(str).str.strip()
    print(f"[INFO] {len(df)} precomputed RBH pairs loaded from {path.name}")
    return df[["9a5c_protein_ID", "Temecula1_protein_ID"]]


def build_dictionary(rbh: pd.DataFrame, tbl_9a5c: pd.DataFrame, tbl_tem: pd.DataFrame) -> pd.DataFrame:
    tbl_9a5c_idx = tbl_9a5c[tbl_9a5c["protein_ID"] != ""].set_index("protein_ID")
    tbl_tem_idx  = tbl_tem[tbl_tem["protein_ID"] != ""].set_index("protein_ID")

    merged = rbh.merge(
        tbl_9a5c_idx, left_on="9a5c_protein_ID", right_index=True, how="inner",
    ).merge(
        tbl_tem_idx, left_on="Temecula1_protein_ID", right_index=True, how="inner",
        suffixes=("_9a5c", "_Temecula1"),
    )
    print(f"[INFO] {len(merged)}/{len(rbh)} RBH pairs resolved to a gene on "
          f"both sides (protein_id present + matched via old_locus_tag "
          f"in both annotation tables)")

    out = pd.DataFrame({
        "9a5c_IMG_ID":              merged["IMG_ID_9a5c"],
        "9a5c_protein_ID":          merged["9a5c_protein_ID"],
        "9a5c_locus_tag":           merged["locus_tag_9a5c"],
        "9a5c_old_locus_tags":      merged["old_locus_tags_9a5c"],
        "9a5c_gene":                merged["gene_9a5c"],
        "9a5c_product_GBFF":        merged["product_GBFF_9a5c"],
        "9a5c_product_IMG":         merged["product_IMG_9a5c"],
        "9a5c_IMG_length":          merged["IMG_length_9a5c"],
        "Temecula1_protein_ID":     merged["Temecula1_protein_ID"],
        "Temecula1_locus_tag":      merged["locus_tag_Temecula1"],
        "Temecula1_old_locus_tags": merged["old_locus_tags_Temecula1"],
        "Temecula1_gene":           merged["gene_Temecula1"],
        "Temecula1_product_GBFF":   merged["product_GBFF_Temecula1"],
        "Temecula1_IMG_ID":         merged["IMG_ID_Temecula1"],
        "Temecula1_product_IMG":    merged["product_IMG_Temecula1"],
        "Temecula1_IMG_length":     merged["IMG_length_Temecula1"],
    })
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Join precomputed RBH pairs with per-strain GenBank/IMG "
                     "annotations into the cross-strain gene dictionary.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--rbh", required=True, type=Path,
                   help="Precomputed reciprocal-best-hit pairs (protein_ID vs protein_ID).")
    p.add_argument("--annot-9a5c", required=True, type=Path,
                   help="annot_comprator_9a5c.csv (old_locus_tag <-> IMG_ID per gene).")
    p.add_argument("--annot-temecula1", required=True, type=Path,
                   help="annot_comprator_Temecula1.csv (old_locus_tag <-> IMG_ID per gene).")
    p.add_argument("--gbk-9a5c", required=True, type=Path,
                   help="GenBank file for strain 9a5c (source of protein_id).")
    p.add_argument("--gbk-tem", required=True, type=Path,
                   help="GenBank file for strain Temecula1 (source of protein_id).")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV file (gene_dictionary.tsv).")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    gbk_9a5c = parse_genbank_by_old_locus_tag(args.gbk_9a5c)
    gbk_tem  = parse_genbank_by_old_locus_tag(args.gbk_tem)

    tbl_9a5c = build_strain_table(args.annot_9a5c, gbk_9a5c, "9a5c")
    tbl_tem  = build_strain_table(args.annot_temecula1, gbk_tem, "Temecula1")

    rbh = load_rbh(args.rbh)
    dictionary = build_dictionary(rbh, tbl_9a5c, tbl_tem)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    dictionary.to_csv(args.output, sep="\t", index=False)
    print(f"[INFO] Gene dictionary saved to {args.output} ({len(dictionary)} ortholog pairs)")


if __name__ == "__main__":
    main()
