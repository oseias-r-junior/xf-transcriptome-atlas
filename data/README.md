# xf-transcriptome-atlas

Reproducible bioinformatics pipeline for **Feitosa-Junior et al. (2025)**:  
*"Transcriptome atlas of* Xylella fastidiosa *across biofilm formation stages reveals a common core of virulence-associated genes."*

All code is organized as self-contained Python scripts with `argparse` interfaces so that any step can be re-run independently.

---

## Pipeline overview

```
RAW DATA
   │
   ▼
01_preprocessing/
   ├─ compute_tpm.py            Raw counts / FPKM  →  TPM matrix
   └─ build_gene_dictionary.py  RBH BLASTP  →  cross-strain gene dictionary
   │
   ▼
02_differential_expression/
   └─ run_deseq2.py             PyDESeq2  →  per-comparison DEG tables
   │
   ▼
03_go_enrichment/
   ├─ parse_go_from_genbank.py  GenBank  →  GO annotation table
   ├─ run_go_enrichment.py      Fisher + BH FDR  →  enriched GO terms
   └─ build_go_heatmap.py       QuickGO API + signed score  →  GO heatmap
   │
   ▼
04_wgcna/
   ├─ run_wgcna.py              PyWGCNA  →  modules, MEs, kME, trait correlations
   └─ build_coexpression_network.py  →  edge/node tables for Cytoscape / Gephi
   │
   ▼
05_figures/
   ├─ fig_pearson_heatmap.py    Pearson r clustermap between conditions (Fig S1)
   ├─ fig_pcoa.py               Bray-Curtis PCoA + PERMANOVA/PERMDISP (Fig 2A)
   ├─ fig_tpm_bubble.py         TPM distribution bubble plot (Fig 2B)
   ├─ fig_upset.py              UpSet plot of expressed-gene intersections (Fig 2C)
   ├─ fig_top100_shared.py      Mean TPM of top 100 RBH shared genes (Fig 3A)
   ├─ fig_virulence_trends.py   TPM trend lines per virulence phase (Fig 3B)
   ├─ fig_virulence_clustermap.py  Log2FC heatmap, phase-ordered (Fig 4A)
   └─ fig_go_bubble.py          Combined GO bubble — DESeq2 + WGCNA modules (Fig 4C)
```

---

## Repository structure

```
xf-transcriptome-atlas/
├── 01_preprocessing/
│   ├── compute_tpm.py
│   └── build_gene_dictionary.py
├── 02_differential_expression/
│   └── run_deseq2.py
├── 03_go_enrichment/
│   ├── parse_go_from_genbank.py
│   ├── run_go_enrichment.py
│   └── build_go_heatmap.py
├── 04_wgcna/
│   ├── run_wgcna.py
│   └── build_coexpression_network.py
├── 05_figures/
│   ├── fig_pearson_heatmap.py
│   ├── fig_pcoa.py
│   ├── fig_tpm_bubble.py
│   ├── fig_upset.py
│   ├── fig_top100_shared.py
│   ├── fig_virulence_trends.py
│   ├── fig_virulence_clustermap.py
│   └── fig_go_bubble.py
├── data/
│   └── README.md               ← expected file formats and origins
├── environment.yml
└── README.md
```

---

## Installation

```bash
conda env create -f environment.yml
conda activate xf-transcriptome
```

---

## Usage

Every script accepts `--help` for full argument documentation.  
Below is the recommended execution order.

### 1 · Preprocessing

```bash
# Convert CLC Workbench FPKM export to TPM
python 01_preprocessing/compute_tpm.py \
    --from-fpkm data/fpkm_matrix.tsv \
    --output    data/tpm_expression.csv

# Build cross-strain gene dictionary via RBH BLASTP
python 01_preprocessing/build_gene_dictionary.py \
    --blastp-fwd data/blast_9a5c_vs_tem.tsv \
    --blastp-rev data/blast_tem_vs_9a5c.tsv \
    --gbk-9a5c   data/9a5c.gbff \
    --gbk-tem    data/temecula1.gbff \
    --output     data/gene_dictionary.tsv
```

### 2 · Differential expression

```bash
python 02_differential_expression/run_deseq2.py \
    --counts      data/raw_counts_combined.tsv \
    --metadata    data/sample_info.csv \
    --dictionary  data/gene_dictionary.tsv \
    --comparisons data/comparisons.tsv \
    --outdir      results/DESeq2_results \
    --alpha 0.05 --lfc 1.0
```

### 3 · GO enrichment

```bash
# Extract GO terms from GenBank files
python 03_go_enrichment/parse_go_from_genbank.py \
    --gbk data/9a5c.gbff data/temecula1.gbff \
    --strain-labels 9a5c temecula1 \
    --output results/go_annotations.tsv

# Run Fisher's exact test + BH FDR
python 03_go_enrichment/run_go_enrichment.py \
    --deseq-dir  results/DESeq2_results \
    --go-annot   results/go_annotations.tsv \
    --dictionary data/gene_dictionary.tsv

# Build integrative GO heatmap (Fig 5)
python 03_go_enrichment/build_go_heatmap.py \
    --deseq-dir results/DESeq2_results \
    --go-dict   data/dictionary_2_level.csv \
    --cache     results/go_ancestor_cache.json \
    --output    figures/fig_go_heatmap.tiff
```

### 4 · WGCNA co-expression network

```bash
python 04_wgcna/run_wgcna.py \
    --tpm      data/tpm_expression.csv \
    --metadata data/sample_info.csv \
    --top-n    640 \
    --outdir   results/WGCNA

python 04_wgcna/build_coexpression_network.py \
    --expr        results/WGCNA/expr_log.csv \
    --assignments results/WGCNA/module_assignments.csv \
    --kme         results/WGCNA/kME.csv \
    --min-cor     0.80 \
    --outdir      results/network
```

### 5 · Figures

```bash
# Fig S1 — Pearson correlation clustermap
python 05_figures/fig_pearson_heatmap.py \
    --correlation-csv data/pearson_correlation.csv \
    --output          figures/fig_S1_pearson_heatmap.tiff
# (or compute from TPM: --tpm data/tpm_expression.csv --metadata data/sample_info.csv)

# Fig 3A — Top 100 shared genes (mean TPM)
python 05_figures/fig_top100_shared.py \
    --tpm        data/tpm_expression.csv \
    --metadata   data/sample_info.csv \
    --dictionary data/gene_dictionary.tsv \
    --top-n      100 \
    --output     figures/fig_3A_top100_shared.tiff

# Fig 2C — UpSet plot of expressed-gene intersections
python 05_figures/fig_upset.py \
    --tpm          data/tpm_expression.csv \
    --metadata     data/sample_info.csv \
    --condition-col condition \
    --strain-col   strain \
    --tpm-min      0 \
    --sort-by      degree \
    --output       figures/fig_2C_upset.tiff

# Fig 2 — PCoA
python 05_figures/fig_pcoa.py \
    --tpm      data/tpm_expression.csv \
    --metadata data/sample_info.csv \
    --group-col condition \
    --output   figures/fig_pcoa.tiff

# Fig 1 — TPM bubble plot
python 05_figures/fig_tpm_bubble.py \
    --tpm        data/tpm_expression.csv \
    --metadata   data/sample_info.csv \
    --colour-col strain \
    --output     figures/fig_tpm_bubble.tiff

# Fig 4A — Virulence gene log2FC clustermap (phase-ordered)
python 05_figures/fig_virulence_clustermap.py \
    --virulence-table data/virulence_table.csv \
    --deseq-dir       results/DESeq2_results \
    --output          figures/fig_4A_virulence_clustermap.tiff

# Fig 4C — Combined GO enrichment bubble (DESeq2 + WGCNA modules)
python 05_figures/fig_go_bubble.py \
    --deseq-dir    results/DESeq2_results \
    --wgcna-dir    results/WGCNA_go \
    --top-n-terms  30 \
    --color-by     ontology \
    --output       figures/fig_4C_go_bubble.tiff

# Fig 4D — Virulence gene expression trends
python 05_figures/fig_virulence_trends.py \
    --tpm             data/tpm_expression.csv \
    --metadata        data/sample_info.csv \
    --virulence-table data/virulence_table.csv \
    --condition-order 9a5c_PIM6_1d 9a5c_PIM6_3d 9a5c_XFM2_1d 9a5c_XFM2_3d \
    --output          figures/fig_4D_virulence_trends.tiff
```

---

## Input data formats

See [`data/README.md`](data/README.md) for complete format specifications of
all required input files.

---

## Methods summary

**Expression quantification.** Raw reads were mapped to the *X. fastidiosa*
9a5c (CVC strain) and Temecula1 reference genomes using CLC Genomics Workbench.
FPKM values were converted to TPM via `compute_tpm.py`.

**Cross-strain gene dictionary.** Orthologous gene pairs between 9a5c and
Temecula1 were identified by reciprocal best-hit (RBH) BLASTP (E-value ≤ 1×10⁻⁵,
identity ≥ 30%).

**Differential expression.** PyDESeq2 (Love et al. 2014; Muzellec et al. 2023)
was applied to all pairwise comparisons of growth conditions. DEGs were called
at FDR ≤ 0.05 and |log₂FC| ≥ 1.0.

**GO enrichment.** Fisher's exact test with Benjamini-Hochberg FDR correction
was performed separately for up- and down-regulated DEGs. The background set
was all detected genes for within-strain comparisons and the union of orthologous
pairs for cross-strain comparisons. GO terms were mapped to level-2 ancestors
via the QuickGO REST API.

**Co-expression network.** The top 640 most variable genes (log₂(TPM+1),
variance > 0.1, TPM ≥ 1 in ≥ 3 samples) were analysed with PyWGCNA
(signed network; soft-threshold power β selected for scale-free R² ≥ 0.8;
deepSplit = 2; MEDissThres = 0.25).

**Ordination.** Bray-Curtis dissimilarity on log₂(TPM+1) profiles was
visualized by PCoA (scikit-bio). Group differences were tested with PERMANOVA
and homogeneity of dispersions with PERMDISP (999 permutations).

---

## Citation

Feitosa-Junior et al. (2025) *eLife* — DOI: [TBD upon acceptance]

---

## License

MIT License. See `LICENSE` for details.
