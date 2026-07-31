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
| `raw_counts_combined.tsv` | TSV, genes × samples | Raw read counts (integer) used as input for DESeq2. |

Example header for `tpm_expression.csv`:
```
gene_id,9a5c_PIM6_1d_1,9a5c_PIM6_1d_2,9a5c_PIM6_1d_3,...
Xf9a_00001,45.3,48.1,41.7,...
```

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
| `temecula1.gbff` | GenBank flat file | *X. fastidiosa* Temecula1 annotated genome. |

Both files are available from NCBI under their respective accession numbers.

### Cross-strain gene dictionary

| File | Format | Description |
|------|--------|-------------|
| `gene_dictionary.tsv` | TSV | Output of `01_preprocessing/build_gene_dictionary.py`. Columns: `9a5c_IMG_ID`, `9a5c_old_locus_tags`, `Temecula1_IMG_ID`, `Temecula1_old_locus_tags`, `identity_pct`. |

### Comparisons list

| File | Format | Description |
|------|--------|-------------|
| `comparisons.tsv` | TSV, two columns, no header | Each row is a pairwise comparison: `condition_1<TAB>condition_2`. Condition names must match the `condition` column in `sample_info.csv`. |

### Virulence gene table

| File | Format | Description |
|------|--------|-------------|
| `virulence_table.csv` | CSV | One row per virulence-associated gene. Required columns: `gene_id`, `gene_name`, `function`, `phase`, `source`. Phase codes: 0=unspecific, 1=mobile, 2=sessile, 3=early, 4=late. |

### GO term resources

| File | Format | Description |
|------|--------|-------------|
| `dictionary_2_level.csv` | CSV, no header | Level-2 GO term dictionary. Columns: `level`, `category_name`, `GO_id`, `term`, `relation`. Used by `03_go_enrichment/build_go_heatmap.py`. |

---

## Directory layout expected by scripts

```
data/
├── tpm_expression.csv
├── raw_counts_combined.tsv
├── sample_info.csv
├── comparisons.tsv
├── gene_dictionary.tsv
├── virulence_table.csv
├── dictionary_2_level.csv
├── 9a5c.gbff
└── temecula1.gbff
```

All scripts accept `--tpm`, `--metadata`, `--dictionary` and similar arguments
so paths are fully configurable. The paths above are the defaults used in
`README.md` usage examples.
