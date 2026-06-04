# *Xylella fastidiosa* Transcriptome Atlas

[![bioRxiv](https://img.shields.io/badge/bioRxiv-2025.08.15.669762-red)](https://doi.org/10.1101/2025.08.15.669762)
[![SRA](https://img.shields.io/badge/SRA-PRJNA1344811-blue)](https://www.ncbi.nlm.nih.gov/sra/PRJNA1344811)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey)](https://creativecommons.org/licenses/by/4.0/)

Reproducibility repository for:

> **Transcriptome profiling reveals differential expression of virulence genes in *Xylella fastidiosa* under nutrient-rich and xylem-like conditions**  
> Paulo M. Pierry‡, Oseias R. Feitosa-Junior‡\*, Joaquim Martins-Junior, Deibs Barbosa, Aline M. da Silva†, Paulo A. Zaini\*  
> *bioRxiv* (2025). doi: [10.1101/2025.08.15.669762](https://doi.org/10.1101/2025.08.15.669762)  
> ‡ Equal contribution · † In memoriam

---

## Overview

This repository contains all Python/Jupyter analysis code used to generate the figures and statistical results of the paper. We profiled transcriptomes of *X. fastidiosa* strains **9a5c** and **Temecula1** across two culture media (PWG and PIM6) at early and late exponential growth phases — 24 libraries in total.

## Repository structure

```
xf-transcriptome-atlas/
├── notebooks/
│   ├── 01_tpm_processing_bubble_plot.ipynb    # Fig. 3B — TPM distribution bubble plots
│   ├── 02_pcoa_correlation_heatmap_upset.ipynb # Fig. 2, 3A, 3C — PCoA, heatmap, UpSet
│   ├── 03_deseq2_differential_expression.ipynb # Table 2, Fig. 4B, 6 — DESeq2 + GO
│   ├── 04_go_enrichment.ipynb                 # Fig. 5 — GO enrichment (BayGO)
│   ├── 05_all_figures_integrated.ipynb        # Master notebook: all publication figures
│   └── 06_wgcna_analysis.ipynb               # Supplementary WGCNA analysis
├── data/
│   └── README.md                             # Data access instructions (SRA + local)
├── figures/                                  # Rendered publication figures
└── environment.yml                           # Conda environment
```

## Raw data

Raw RNA-Seq reads are deposited at NCBI SRA under **BioProject [PRJNA1344811](https://www.ncbi.nlm.nih.gov/sra/PRJNA1344811)** (24 paired-end libraries, 2×250 bp MiSeq).

> **Note:** Data become publicly available on **2026-10-01**, coinciding with journal submission.

| Run | Condition |
|-----|-----------|
| SRR46283131–SRR46283154 | 9a5c and Temecula1 × PIM6/PWG × early/late (3 replicates each) |

## Setup

```bash
# Clone
git clone https://github.com/oseias-r-junior/xf-transcriptome-atlas.git
cd xf-transcriptome-atlas

# Create environment
conda env create -f environment.yml
conda activate xf-transcriptome

# Launch
jupyter lab
```

## Reproducing figures

Each notebook is self-contained and numbered in analysis order. Set `DATA_DIR` at the top of each notebook to your local data path (see `data/README.md`).

## Citation

```bibtex
@article{pierry2025xylella,
  title   = {Transcriptome profiling reveals differential expression of virulence genes
             in \textit{Xylella fastidiosa} under nutrient-rich and xylem-like conditions},
  author  = {Pierry, Paulo M. and Feitosa-Junior, Oseias R. and Martins-Junior, Joaquim
             and Barbosa, Deibs and da Silva, Aline M. and Zaini, Paulo A.},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.08.15.669762}
}
```

## License

Code: [MIT](LICENSE-CODE)  
Text and figures: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
