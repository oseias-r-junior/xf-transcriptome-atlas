# Data directory

This directory holds the input data files required to run the pipeline.
**Raw data are not included in this repository** to respect data-sharing
agreements and file-size constraints. The table below describes the expected
format and origin of each file.

Place the corresponding files at the paths indicated before running any script.

---

## Required files

### Expression data

| File | Format | Description |
|------|--------|-------------|
| `tpm_expression.csv` | CSV, genes × samples | TPM values computed by `01_preprocessing/compute_tpm.py`. Gene IDs from CLC Genomics Workbench IMG locus tags. |
| `raw_counts_combined.tsv` | TSV, genes × samples | Raw read counts (integer) used as input for DESeq2. Built by `01_preprocessing/build_raw_counts_combined.py` from the per-sample raw-count `.txt` files (CLC exports, one per sample, each already quantified against its own strain's genome). This is a **strain-native outer join** — no ortholog merging happens at this stage: a 9a5c gene's columns are 0 for every Temecula1 sample and vice-versa, since that gene ID never appears in the other strain's own count files. The ortholog merge (via `gene_dictionary.tsv`) is applied later, per comparison, only for cross-strain comparisons, inside `02_differential_expression/run_deseq2.py` — see that script's docstring for why (within-strain comparisons must keep that strain's full native gene universe, matching `03_go_enrichment/run_go_enrichment.py`'s documented background-set convention). |

Example header for `tpm_expression.csv`:
```
gene_id,9a5c_PIM6_1d_1,9a5c_PIM6_1d_2,9a5c_PIM6_1d_3,...
Xf9a_00001,45.3,48.1,41.7,...
```

### WGCNA input matrix -- do NOT use `tpm_expression.csv` directly

`tpm_expression.csv` above is **strain-native**: "genes" is the union of
each strain's own native gene IDs, so every row only has real values for
ONE strain's 12 columns (the other strain's 12 are hard structural zeros,
not biology). Feeding it straight into `04_wgcna/run_wgcna.py` across all
24 samples at once confounds strain identity with co-expression and does
not match the original analysis. WGCNA needs a matrix keyed by a shared
cross-strain ortholog identity instead. Two such files exist, and which one
you use matters:

| File | Genes | Source | Use for |
|------|-------|--------|---------|
| `tpm_expression_original.csv` | 1692 (1685 after the standard TPM≥1-in-≥3-samples + var>0.1 filter) | The exact ortholog-merged matrix from `WGCNA_paper_pierry_feitosa.ipynb`, gene_id = `<9a5c_IMG_ID>_<Temecula1_IMG_ID>`. **This is the file that produced the WGCNA figure/results already in the manuscript.** | Reproducing/consolidating the submitted WGCNA result (3 modules: dimgrey/darkgrey/silver) and the downstream robustness analyses that characterize it (leave-one-out stability, power/deepSplit sensitivity, module-trait correlation). Treat as a frozen, canonical input -- do not regenerate it from raw data. |
| `tpm_expression_wgcna.csv` | 1639 | Built from raw data by `01_preprocessing/build_wgcna_tpm_matrix.py --tpm tpm_expression.csv --dictionary gene_dictionary.tsv`, i.e. fully reproducible from the raw CLC `.xlsx` exports through this repo's own scripts. | Demonstrating end-to-end reproducibility from raw data in this GitHub repo. Differs from `tpm_expression_original.csv` by 53 gene pairs (~3.1%) because `build_gene_dictionary.py` excludes ortholog pairs whose RefSeq `protein_id` is shared by more than one locus in the same strain (ambiguous 1:1 ortholog assignment) -- a deliberate correctness fix, not a bug. Confirmed (2026-09) that running WGCNA on `tpm_expression_original.csv` with auto-detected power/cutHeight reproduces the original 3-module result exactly (same PyWGCNA 2.2.1 in both environments), so this residual 53-gene gap is the full and only explanation for the module-count difference between the two files -- not a code or environment issue. |

Both are ortholog-pair-keyed (`gene_id` = `<9a5c_IMG_ID>_<Temecula1_IMG_ID>`,
e.g. `XF9a_00002_XFTem_00002`) and both align against `sample_info.csv`'s
sample naming after `04_wgcna/run_wgcna.py`'s own `normalize_sample_name()`.

### Sample metadata

| File | Format | Description |
|------|--------|-------------|
| `sample_info.csv` | CSV, samples × traits | One row per sample. Required columns: `sample_id`, `condition`, `strain`, `medium`, `timepoint`. Binary trait columns (`is_mobile`, `is_sessile`, `is_early`, `is_late`) are used for WGCNA module–trait correlations. |

Example:
```
sample_id,condition,strain,medium,timepoint,is_mobile,is_sessile
9a5c_PIM6_1d_1,9a5c_PIM6_1d,9a5c,PIM6,1d,0,0
```

### Reference genome files

| File | Format | Description |
|------|--------|-------------|
| `9a5c.gbff` | GenBank flat file | *X. fastidiosa* 9a5c annotated genome (RefSeq/NCBI). |
| `Temecula1.gbff` | GenBank flat file | *X. fastidiosa* Temecula1 annotated genome. |

Both files are available from NCBI under their respective accession numbers.

### Cross-strain gene dictionary

`gene_dictionary.tsv` (Supplementary Table S4) is built by
`01_preprocessing/build_gene_dictionary.py` from three independent
per-strain/cross-strain sources — it does **not** run BLASTP itself, and it
does **not** treat a GenBank `protein_id` as the IMG ID (an earlier version
of this script conflated the two, which is wrong: they are different ID
systems for the same gene). The three required inputs:

| File | Format | Origin |
|------|--------|--------|
| `9a5c_vs_Temecula1_rbh.csv` | CSV, 2 columns (`9a5c_protein_ID`, `Temecula1_protein_ID`) | Precomputed reciprocal-best-hit (RBH) pairs from BLASTP of the two strains' proteomes (protein sequence vs protein sequence), run separately/upstream of this repo. This is the **only** cross-strain link the dictionary uses. **Bundled in this repository.** |
| `annot_comprator_9a5c.csv` / `annot_comprator_Temecula1.csv` | TSV (despite the `.csv` extension), one row per gene | Per-strain table linking, for each gene, its old locus tag + product/length (`annot_1.*`), a secondary NCBI-style numbering kept only for reference (`annot_2.*`), and the **IMG ID** used throughout this repository (`XF9a_#####` / `XFTem_#####`) + IMG product/length (`annot_3.*`). Pre-existing colleague-computed NCBI↔IMG cross-references; does not carry `protein_id`. **Bundled in this repository.** |
| `9a5c.gbff` / `Temecula1.gbff` | GenBank flat file | See "Reference genome files" above. Parsed only to recover, per `old_locus_tag`, the `protein_id` that bridges the RBH pairs to the IMG IDs in the annotation-comparator tables. |

Join logic: for each strain, `old_locus_tag` (from the GenBank file, matched
against `annot_comprator_<strain>.csv`'s `annot_1.1` allowing both
`XF0677`/`XF_0677` underscore conventions) attaches a `protein_id` to every
IMG ID; the two strains' resulting tables are then merged on the RBH pairs.
Genes with no IMG ID in the annotation-comparator table, or no `protein_id`
match in the GenBank file, are dropped from that step (reported in the
script's console output) rather than silently mis-joined.

Required columns used by downstream scripts: `9a5c_IMG_ID`,
`Temecula1_IMG_ID`, `Temecula1_old_locus_tags`. Extra columns (product
names, protein IDs, IMG lengths, gene symbols) are carried through for
annotation but ignored by most scripts.

### Sample metadata and comparisons

These two files are **included in the repository** under `data/` and do not need to be re-created. They describe the experimental design and must match the column names of the TPM/counts matrices exported from CLC Genomics Workbench.

| File | Format | Description |
|------|--------|-------------|
| `sample_info.csv` | CSV, samples × traits | One row per sample. Required columns: `sample_id`, `condition`, `strain`, `medium`, `timepoint`. Binary trait columns (`is_mobile`, `is_sessile`, `is_early`, `is_late`) are used for WGCNA module–trait correlations. |
| `comparisons.tsv` | TSV, two columns, no header | Each row is a pairwise comparison: `condition_1<TAB>condition_2`. Condition names must match the `condition` column in `sample_info.csv`. |

### Virulence gene table

| File | Format | Origin |
|------|--------|--------|
| `virulence_table.csv` | CSV | **Supplementary Table S5** of Feitosa-Junior et al. (2025). Export Table S5 from the supplementary XLSX and save as `virulence_table.csv`. Column names are auto-detected by scripts (case-insensitive): any column containing `gene_id` serves as the gene identifier (Temecula1 PD#### format), `Phase` as virulence phase, `Function` as functional group. Phase codes: 0=unspecific, 1=mobile, 2=sessile, 3=early, 4=late. |

### GO term resources

| File | Format | Origin |
|------|--------|--------|
| `dictionary_2_level.csv` | CSV, no header | Level-2 GO term dictionary. Columns: `level`, `category_name`, `GO_id`, `term`, `relation`. Used by `03_go_enrichment/build_go_heatmap.py` and `05_figures/fig_go_bubble.py`. **Bundled in this repository.** Originally generated by [oseias-r-junior/Gene_Ontology_2nd_Level](https://github.com/oseias-r-junior/Gene_Ontology_2nd_Level); the version included here is the one used in the analyses reported in Feitosa-Junior et al. (2025). |
| `ancestor_cache.json` | JSON | GO id → list of ancestor GO ids, as originally fetched from QuickGO when the Figure 4C analysis was run. Used by `fig_go_bubble.py` to collapse specific GO terms to their Level-2 ancestor without needing internet access. **Bundled in this repository** (339 entries — covers the GO ids that were actually significant in the original run; a GO id not in the cache is only matched against the dictionary as itself, not via its ancestors, which is a known limitation for any newly-computed enrichment result outside the original analysis). |

### Fig 4C inputs (DESeq2 + WGCNA GO enrichment, Level-2 terms)

| Path | Format | Origin |
|------|--------|--------|
| `deseq2_go_results/<c1>_vs_<c2>/GO_enrichment_up_in_<cond>_old_locus_tags.csv` | CSV | Per-comparison, per-direction GO enrichment (`go_term`, `go_description`, `p_adj`, …), one sub-folder per pairwise comparison. This is a **lean, CSV-only** copy of the relevant files from `DeSeq2/DESeq2_results/<comparison>/` in the source analysis folder — the source folders also contain large `.tiff` renders and count matrices that `fig_go_bubble.py` does not need. |
| `wgcna_enrichment/enrichment_<module>.csv` | CSV | Per-module GO enrichment (`GO_ID`, `P_adj_BH`, …) for the `dimgrey`/`darkgrey`/`silver` WGCNA modules. Note: this file's own `Level2_ID`/`Level2_Label` columns are a no-op copy of `GO_ID`/`GO_Label` in the source data (not actually collapsed) — `fig_go_bubble.py` re-derives the true Level-2 term via `dictionary_2_level.csv` + `ancestor_cache.json` instead of trusting those columns. |

### Biofilm ring quantification input

| File | Format | Origin |
|------|--------|--------|
| `biofilm_ring_area_results.csv` | CSV | ImageJ/Fiji `Results` window export (Label, Area, Mean, MinThr, MaxThr) for the biofilm-ring densitometry described in Methods 2.X. One row per flask photograph; `Label` encodes strain (`9a`/`tem1`), medium (`pim6`/`pwg`) and phase (`early`/`late`), e.g. `9a_pim6_late_1`. Used by `07_biofilm_quantification/analyze_biofilm_areas.py`. |
| `flask_width_measurements.csv` | CSV (Label, flask_width_px) | Fiji straight-line length measurement of each flask's outer width at the same row height as the ring band (Analyze > Measure, Length only), one row per photograph. Used by `07_biofilm_quantification/normalize_by_flask_width.py` to normalize ring area by (flask width)² and remove the camera-distance/zoom confound between photographs before statistical comparison — see Methods 2.X. `data/flask_width_measurements.txt` is kept alongside it as the raw Fiji Results-window export this CSV was parsed from. |

### Fig 5 network input

| File | Format | Origin |
|------|--------|--------|
| `expr_log.csv` (passed via `--wgcna-dir`, alongside `module_assignments.csv` and `kME.csv`) | CSV, samples × genes | log-transformed expression matrix used to compute the pairwise gene correlations for network edges. In the source analysis folder this is `WGCNA/expr_mat.csv` (same content, samples-as-rows) — copy/rename it to `expr_log.csv` in whichever directory you pass as `--wgcna-dir`. |

---

## Directory layout expected by scripts

```
data/
├── tpm_expression.csv
├── raw_counts_combined.tsv      # build_raw_counts_combined.py output
├── sample_info.csv
├── comparisons.tsv
├── gene_dictionary.tsv          # build_gene_dictionary.py output
├── 9a5c_vs_Temecula1_rbh.csv    # build_gene_dictionary.py input (bundled)
├── annot_comprator_9a5c.csv     # build_gene_dictionary.py input (bundled)
├── annot_comprator_Temecula1.csv # build_gene_dictionary.py input (bundled)
├── virulence_table.csv
├── dictionary_2_level.csv
├── ancestor_cache.json
├── deseq2_go_results/
│   └── <c1>_vs_<c2>/GO_enrichment_up_in_<cond>_old_locus_tags.csv
├── wgcna_enrichment/
│   └── enrichment_<module>.csv
├── 9a5c.gbff
└── Temecula1.gbff
```

All scripts accept `--tpm`, `--metadata`, `--dictionary` and similar arguments
so paths are fully configurable. The paths above are the defaults used in
`README.md` usage examples.
