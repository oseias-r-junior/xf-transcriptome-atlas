# Data Access

## SRA (public after 2026-10-01)

All 24 raw RNA-Seq libraries are deposited at:  
**BioProject: [PRJNA1344811](https://www.ncbi.nlm.nih.gov/sra/PRJNA1344811)**

Download with SRA Toolkit:
```bash
# Install SRA toolkit, then:
prefetch PRJNA1344811
fasterq-dump --outdir fastq/ SRR46283131  # repeat for each run
```

## Local data structure (for notebook execution)

Set `DATA_DIR` in each notebook to the root of your local data folder:

```
DATA_DIR/
├── Arquivos brutos/
│   ├── TPM/              # TPM tables per sample (*.csv)
│   ├── FPKM/             # FPKM tables per sample
│   └── genomes_9a5c_temecula1/   # Reference GBK/FASTA files
│       ├── 9a5c.gbff
│       ├── Temecula1.gbff
│       └── 9a5c_vs_Temecula1_rbh.csv
├── DeSeq2/
│   ├── *.txt             # Raw count tables (24 files: strain_medium_day_replicate)
│   ├── gene_dictionary.csv
│   ├── virulence_table.csv
│   └── DESeq2_results/   # Output tables from DESeq2
└── Supplementary tables/
    └── Pierry_Feitosa_et_al_2025_supplementary tables.xlsx
```

## Reference genomes

| Strain | GenBank accession | Description |
|--------|-------------------|-------------|
| 9a5c   | AE003849 + AE003850 | *X. fastidiosa* 9a5c chromosome + pXF51 plasmid |
| Temecula1 | AE009442 | *X. fastidiosa* Temecula1 |
