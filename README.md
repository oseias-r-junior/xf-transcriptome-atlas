# xf-transcriptome-atlas

Reproducible bioinformatics pipeline for:

**"Transcriptome profiling reveals differential expression of virulence genes in *Xylella fastidiosa* under nutrient-rich and xylem-like conditions"**

Paulo M. Pierry<sup>‡</sup>, Oseias R. Feitosa-Junior<sup>‡*</sup>, Joaquim Martins-Junior, Deibs Barbosa, Aline M. da Silva<sup>†</sup>, Paulo A. Zaini

<sup>‡</sup> These authors contributed equally. <sup>†</sup> in memorian. <sup>*</sup> Corresponding author. 

Submitted to *Pathogens* (MDPI).

This repository takes the two *X. fastidiosa* strains' raw per-sample expression data all the way through to every figure and statistical result in the manuscript: read counts → TPM normalization → differential expression → GO enrichment → co-expression network (WGCNA) → ordination/virulence-gene figures → biofilm-phenotype quantification. Every stage is a standalone Python script with an `argparse` interface, so any single step can be re-run, inspected, or extended on its own without re-running the whole pipeline.

---

## Environment

```
conda env create -f xf-transcriptome-atlas.yml
conda activate xf-transcriptome
```

---

## Pipeline overview

```
RAW DATA (CLC Genomics Workbench "Gene Expression" .xlsx exports, one per sample)
   │
   ▼
01_preprocessing/
   ├─ build_matrices_from_fpkm_xlsx.py   FPKM/<strain>/<medium>/*.xlsx [+pXF51/]  →  raw_counts_combined.tsv + gene_lengths.tsv
   │                                     (strain-native outer join; no ortholog merge yet)
   ├─ compute_tpm.py                     raw_counts_combined.tsv + gene_lengths.tsv  →  tpm_expression.csv
   ├─ build_gene_dictionary.py           RBH pairs + per-strain NCBI/IMG annotations  →  gene_dictionary.tsv (Supplementary Table S5)
   └─ build_wgcna_tpm_matrix.py          tpm_expression.csv + gene_dictionary.tsv  →  tpm_expression_wgcna.csv (ortholog-merged, for WGCNA)

   (build_raw_counts_combined.py is superseded: an earlier version of this
   pipeline assumed a flat directory of per-sample .txt count files, which is
   not how the raw data is actually organized. Kept only for reference.)
   │
   ▼
02_differential_expression/
   └─ run_deseq2.py     raw_counts_combined.tsv + sample_info.csv + comparisons.tsv  →  results/DESeq2_results/<c1>_vs_<c2>/
   │
   ▼
03_go_enrichment/
   ├─ parse_go_from_genbank.py   9a5c.gbff + Temecula1.gbff + annotation tables  →  go_annotations.tsv
   ├─ run_go_enrichment.py       DESeq2 output + go_annotations.tsv  →  per-comparison GO enrichment CSVs
   └─ build_go_heatmap.py        GO enrichment across all comparisons  →  figures/go_heatmap.tiff
   │
   ▼
04_wgcna/
   ├─ run_wgcna.py                    tpm_expression_wgcna.csv + sample_info.csv  →  results/WGCNA/ (modules, kME, hub genes, trait correlations)
   ├─ module_trait_heatmap.py         module_eigengenes.csv + sample_info.csv  →  figures/fig_module_trait_heatmap.tiff (BH-adjusted p-values)
   ├─ build_coexpression_network.py   results/WGCNA/  →  results/network/ (Cytoscape/Gephi edge & node tables)
   ├─ wgcna_stability.py              leave-one-out module-stability check
   └─ wgcna_sensitivity.py            power/deepSplit parameter sweep
   │
   ▼
05_figures/
   ├─ fig_pcoa.py                  Figure 2A — PCoA + PERMANOVA/PERMDISP (Euclidean distance, raw TPM)
   ├─ fig_tpm_bubble.py            Figure 2B — TPM expression distribution bubble plot
   ├─ fig_upset.py                 Figure 2C — UpSet plot of expressed-gene intersections
   ├─ fig_top100_shared.py         Figure 3A — top-100 shared-gene expression trends
   ├─ fig_virulence_trends.py      Figure 3B — virulence gene trends (mobile vs. sessile)
   ├─ fig_virulence_clustermap.py  Figure 4A — virulence gene log2FC heatmap
   ├─ fig_go_bubble.py             Figure 4C — combined DESeq2/WGCNA GO-term bubble plot
   ├─ fig_network.py               Figure 5 — co-expression network (WGCNA hubs, DEGs, virulence genes)
   └─ fig_pearson_heatmap.py       Figure S3 — Pearson correlation clustermap
   │
   ▼
06_stats/
   ├─ pure_stats.py            dependency-light NumPy/pandas reimplementation of PERMANOVA, Pearson-r, BH-FDR
   └─ permanova_factors.py     factor-specific PERMANOVA (strain / medium / growth phase), global + stratified + within-strain
   │
   ▼
07_biofilm_quantification/
   ├─ recompute_from_embedded_roi.py    FIJI/masks/*.tif (embedded ImageJ ROI)  →  biofilm_ring_area_results_otsu.csv (per-image Otsu threshold)
   ├─ diagnose_threshold_coverage.py    visual QC contact sheet  →  figures/fig_S1_biofilm_threshold_qc_contact_sheet.png
   ├─ normalize_by_flask_width.py       ring area ÷ (flask width)²  →  figures/fig_S2_biofilm_ring_area_normalized.tiff
   ├─ analyze_biofilm_areas.py          Welch t-tests + BH correction, strain/medium/phase contrasts
   └─ measure_biofilm_ring.py           single-image demonstration/verification of the Fiji processing chain (--demo mode needs no input photo)
```

---

## Repository structure

```
xf-transcriptome-atlas/
├── 01_preprocessing/
│   ├── build_matrices_from_fpkm_xlsx.py
│   ├── build_raw_counts_combined.py   (superseded, kept for reference)
│   ├── build_gene_dictionary.py
│   ├── build_wgcna_tpm_matrix.py
│   └── compute_tpm.py
├── 02_differential_expression/
│   └── run_deseq2.py
├── 03_go_enrichment/
│   ├── parse_go_from_genbank.py
│   ├── run_go_enrichment.py
│   └── build_go_heatmap.py
├── 04_wgcna/
│   ├── run_wgcna.py
│   ├── module_trait_heatmap.py
│   ├── build_coexpression_network.py
│   ├── wgcna_stability.py
│   └── wgcna_sensitivity.py
├── 05_figures/
│   ├── fig_tpm_bubble.py
│   ├── fig_pcoa.py
│   ├── fig_upset.py
│   ├── fig_top100_shared.py
│   ├── fig_virulence_trends.py
│   ├── fig_virulence_clustermap.py
│   ├── fig_go_bubble.py
│   ├── fig_network.py
│   └── fig_pearson_heatmap.py
├── 06_stats/
│   ├── pure_stats.py
│   └── permanova_factors.py
├── 07_biofilm_quantification/
│   ├── recompute_from_embedded_roi.py
│   ├── diagnose_threshold_coverage.py
│   ├── normalize_by_flask_width.py
│   ├── analyze_biofilm_areas.py
│   ├── measure_biofilm_ring.py
│   └── FIJI/
│       ├── README.md          (explains what belongs in masks/)
│       └── masks/             (empty placeholder — see below)
├── data/                      (inputs; see data/README.md for the full file-by-file description)
├── figures/                   (generated by 05_figures/*.py; not distributed, see .gitignore)
├── results/                   (generated intermediate outputs of every stage; not distributed)
├── photos/                    (Fiji-converted 8-bit flask photographs for the biofilm pipeline; authors' own bench images, not distributed -- see 07_biofilm_quantification/FIJI/README.md)
└── xf-transcriptome-atlas.yml (conda environment)
```

---

## 1 · Preprocessing

```bash
python 01_preprocessing/build_matrices_from_fpkm_xlsx.py \
    --fpkm-dir          FPKM/ \
    --counts-column     "Unique gene reads" \
    --counts-output     data/raw_counts_combined.tsv \
    --lengths-output    data/gene_lengths.tsv \
    --manifest-output   data/fpkm_ingest_manifest.csv

python 01_preprocessing/compute_tpm.py \
    --counts  data/raw_counts_combined.tsv \
    --lengths data/gene_lengths.tsv \
    --output  data/tpm_expression.csv

python 01_preprocessing/build_gene_dictionary.py \
    --rbh              data/9a5c_vs_Temecula1_rbh.csv \
    --annot-9a5c        data/annot_comprator_9a5c.csv \
    --annot-temecula1   data/annot_comprator_Temecula1.csv \
    --gbk-9a5c          data/9a5c.gbff \
    --gbk-tem           data/Temecula1.gbff \
    --output            data/gene_dictionary.tsv
```

`compute_tpm.py` is a fully vectorized pandas implementation (no per-sample loop, so there is nothing to parallelize); it also has a `--from-fpkm` mode for tables that are already length-normalized, in which case `--lengths` is not needed.

The WGCNA-ready, ortholog-merged matrix (`04_wgcna/run_wgcna.py` and its companion scripts) is built separately — see "4 · WGCNA" below for why `tpm_expression.csv` cannot be used directly for that step.

---

## 2 · Differential expression (DESeq2)

```bash
python 02_differential_expression/run_deseq2.py \
    --counts    data/raw_counts_combined.tsv \
    --metadata  data/sample_info.csv \
    --comparisons data/comparisons.tsv \
    --dictionary  data/gene_dictionary.tsv \
    --alpha     0.05 \
    --lfc       1.0 \
    --outdir    results/DESeq2_results
```

Within-strain comparisons use that strain's full native gene universe; cross-strain comparisons use the ortholog union from `gene_dictionary.tsv`. Significance is `padj <= alpha` (non-strict), matching the analysis notebook this script reproduces.

---

## 3 · GO enrichment

```bash
python 03_go_enrichment/parse_go_from_genbank.py \
    --gbk-9a5c         data/9a5c.gbff \
    --gbk-temecula1    data/Temecula1.gbff \
    --annot-9a5c       data/annot_comprator_9a5c.csv \
    --annot-temecula1  data/annot_comprator_Temecula1.csv \
    --output           results/go_annotations.tsv

python 03_go_enrichment/run_go_enrichment.py \
    --deseq-dir   results/DESeq2_results \
    --go-annot    results/go_annotations.tsv \
    --dictionary  data/gene_dictionary.tsv \
    --alpha       0.05 \
    --lfc         1.0

python 03_go_enrichment/build_go_heatmap.py \
    --deseq-dir results/DESeq2_results \
    --go-dict   data/dictionary_2_level.csv \
    --cache     results/go_ancestor_cache.json \
    --output    figures/go_heatmap.tiff
```

---

## 4 · WGCNA co-expression network

`04_wgcna/run_wgcna.py` needs a cross-strain, **ortholog-merged** TPM matrix (`gene_id = "<9a5c_IMG_ID>_<Temecula1_IMG_ID>"`), not the strain-native `data/tpm_expression.csv` (which is a union of each strain's own gene IDs and would confound co-expression with "which strain has data for this gene"). Build it first from the shipped `tpm_expression.csv`:

```bash
python 01_preprocessing/build_wgcna_tpm_matrix.py \
    --tpm        data/tpm_expression.csv \
    --dictionary data/gene_dictionary.tsv \
    --output     data/tpm_expression_wgcna.csv
```

This reproduces the 3-module WGCNA result (dimgrey/darkgrey/silver) to within a small residual (~3% fewer ortholog pairs than the exact matrix behind the submitted manuscript figure -- see `data/README.md`'s "WGCNA input matrix" section for why). Then:

```bash
python 04_wgcna/run_wgcna.py \
    --tpm       data/tpm_expression_wgcna.csv \
    --metadata  data/sample_info.csv \
    --top-n     640 \
    --outdir    results/WGCNA

python 04_wgcna/module_trait_heatmap.py \
    --eigengenes  results/WGCNA/module_eigengenes.csv \
    --metadata    data/sample_info.csv \
    --traits      strain medium timepoint \
    --output      figures/fig_module_trait_heatmap.tiff

python 04_wgcna/build_coexpression_network.py \
    --expr            results/WGCNA/expr_log.csv \
    --assignments     results/WGCNA/module_assignments.csv \
    --kme             results/WGCNA/kME.csv \
    --min-cor         0.80 \
    --hub-n           20 \
    --gene-dict       data/gene_dictionary.tsv \
    --virulence-table data/virulence_table.csv \
    --outdir          results/network
```

Robustness checks on the module structure:

```bash
python 04_wgcna/wgcna_stability.py \
    --tpm       data/tpm_expression_wgcna.csv \
    --metadata  data/sample_info.csv \
    --top-n     640 \
    --power     <beta from soft_threshold_plot.tiff> \
    --reference-assignments results/WGCNA/module_assignments.csv \
    --outdir    results/WGCNA_stability

python 04_wgcna/wgcna_sensitivity.py \
    --tpm        data/tpm_expression_wgcna.csv \
    --metadata   data/sample_info.csv \
    --top-n      640 \
    --main-power <same beta> \
    --reference-assignments results/WGCNA/module_assignments.csv \
    --outdir     results/WGCNA_sensitivity
```

`wgcna_stability.py` (leave-one-out: does the 3-module structure survive removing any single sample?) and `wgcna_sensitivity.py` (does it survive a power ± 1 / deepSplit 1-3 sweep?) both require PyWGCNA — they are not reimplemented in pure NumPy, since re-deriving `cutreeHybrid`'s dendrogram-cut behaviour by hand risks silently producing different module boundaries than PyWGCNA would. These two checks are kept as optional, on-demand analyses -- runnable from this repository but not tied to a numbered manuscript figure or table.

---

## 5 · Figures

```bash
python 05_figures/fig_tpm_bubble.py \
    --tpm data/tpm_expression.csv --metadata data/sample_info.csv \
    --tpm-min 1.0 --group-col condition --colour-col strain \
    --output figures/fig_tpm_bubble.tiff

python 05_figures/fig_pcoa.py \
    --tpm data/tpm_expression.csv --metadata data/sample_info.csv \
    --group-col condition --n-perms-permanova 999 --n-perms-permdisp 9999 \
    --output figures/fig_pcoa.tiff --emperor-html figures/emperor_pcoa

python 05_figures/fig_upset.py \
    --tpm data/tpm_expression.csv --metadata data/sample_info.csv \
    --condition-col condition --tpm-min 0 --sort-by degree \
    --output figures/fig_2C_upset.tiff

python 05_figures/fig_top100_shared.py \
    --tpm data/tpm_expression.csv --top-n 100 --alpha 0.05 \
    --output figures/fig_3A_top100_shared.tiff

python 05_figures/fig_virulence_trends.py \
    --tpm data/tpm_expression_wgcna.csv \
    --dictionary data/gene_dictionary.tsv \
    --virulence-table data/virulence_table.csv \
    --output figures/fig_3B_virulence_trends.tiff

python 05_figures/fig_virulence_clustermap.py \
    --virulence-table data/virulence_table.csv \
    --gene-dict data/gene_dictionary.tsv \
    --deseq-dir results/DESeq2_results \
    --output figures/fig_4A_virulence_clustermap.tiff

python 05_figures/fig_go_bubble.py \
    --deseq-dir      data/deseq2_go_results \
    --wgcna-dir      data/wgcna_enrichment \
    --dictionary     data/dictionary_2_level.csv \
    --ancestor-cache data/ancestor_cache.json \
    --output         figures/fig_4C_go_bubble.tiff

python 05_figures/fig_network.py \
    --wgcna-dir results/WGCNA \
    --deseq-dir results/DESeq2_results \
    --dictionary data/gene_dictionary.tsv \
    --virulence-table data/virulence_table.csv \
    --output figures/fig_5_network.tiff

python 05_figures/fig_pearson_heatmap.py \
    --tpm data/tpm_expression.csv --metadata data/sample_info.csv \
    --output figures/fig_S3_pearson_heatmap.tiff
```

The PCoA (Figure 2A) and the factor-specific PERMANOVA below both use Euclidean distance on the raw (non-log-transformed) TPM matrix, matching the original analysis notebook's `pdist()` call exactly. Log2(TPM+1)/Bray-Curtis alternatives exist as flags (`--log-transform`, `--metric`) purely for sensitivity comparison; they are off by default and are not the metric behind the reported results.

---

## 6 · Statistics

```bash
python 06_stats/permanova_factors.py \
    --tpm    data/tpm_expression_wgcna.csv --tpm-sep "," \
    --output results/permanova/permanova_factors.csv \
    --n-perms 999
```

`permanova_factors.py` runs three complementary tests per factor (strain / medium / growth phase) so that each addresses a specific confound: a global one-way PERMANOVA per factor, a stratified PERMANOVA that permutes only within strain (controlling for the strain confound when testing medium or phase), and an independent PERMANOVA re-run inside each strain's own 12-sample subset. Together these establish which factor(s) actually drive the separation seen in the global 8-group PERMANOVA/PERMDISP test in `fig_pcoa.py`, and whether one factor's apparent effect is a confound of another.

`06_stats/pure_stats.py` is a dependency-light NumPy/pandas reimplementation of the PERMANOVA, Pearson-r, and Benjamini-Hochberg routines used above and in `module_trait_heatmap.py`, validated against known reference values in its own `__main__` block (`python 06_stats/pure_stats.py`). It exists so these analyses can be reproduced or audited in any Python environment, without requiring scikit-bio/statsmodels.

---

## 7 · Biofilm ring quantification

*X. fastidiosa* forms a visible biofilm ring at the air-liquid interface of static flask cultures. This was quantified from photographs of each flask, processed and measured in Fiji/ImageJ, and re-derived here programmatically for a fully scriptable, auditable version of that measurement.

The 22 flask photographs (`photos/*.tif`) are Fiji-converted 8-bit grayscale intermediates, not raw camera images and not background-subtracted/thresholded outputs — those two processing steps are applied by the scripts below, which is also why this stage isn't a single batch/threaded pass over all images: the original workflow required a human to visually confirm the threshold window per photograph (the per-image threshold variation captured in `data/biofilm_ring_area_results.csv` reflects that).

`measure_biofilm_ring.py` is a single-image demonstration of the full Fiji processing chain (background subtraction → threshold → area measurement), including a `--demo` mode that runs on a synthetic test image so the pipeline can be verified with no input photo at all. The batch-capable path for re-measuring all 22 images at once is:

```bash
# 1. Re-measure ring area from each mask's embedded ImageJ ROI, per-image Otsu threshold
python 07_biofilm_quantification/recompute_from_embedded_roi.py \
    --output data/biofilm_ring_area_results_otsu.csv

# 2. Visual QC: confirm the threshold isn't cutting into real biofilm texture
python 07_biofilm_quantification/diagnose_threshold_coverage.py

# 3. Normalize by flask width (removes camera-distance/zoom confound)
python 07_biofilm_quantification/normalize_by_flask_width.py \
    --results     data/biofilm_ring_area_results.csv \
    --flask-width data/flask_width_measurements.csv \
    --outdir      results/biofilm_quantification \
    --scale-factor 1e6

# 4. Statistics: Welch t-tests (strain / medium / phase), BH-corrected
python 07_biofilm_quantification/analyze_biofilm_areas.py \
    --results data/biofilm_ring_area_results.csv \
    --outdir  results/biofilm_quantification
```

Steps 1-2 read from `07_biofilm_quantification/FIJI/masks/` by default. **These 22 background-subtracted flask TIFFs are the authors' own manual Fiji measurements and are not distributed with this repository** — `FIJI/masks/` ships as an empty placeholder (see `FIJI/README.md` for the expected file format) so the folder's role is unambiguous; place your own copies there, or pass `--masks-dir /path/to/your/masks`, before running steps 1-2.

- Figure S1 — `figures/fig_S1_biofilm_threshold_qc_contact_sheet.png` (step 2 output)
- Figure S2 — `figures/fig_S2_biofilm_ring_area_normalized.tiff` (step 3 output)

---

## Sample metadata schema

`data/sample_info.csv` — one row per sample, standardized schema used by every script's `--metadata` argument:

```csv
sample_id,strain,medium,timepoint,replicate,condition
9a5c_PIM6_1d_rep1,9a5c,PIM6,1d,1,9a5c_PIM6_1d
```

`sample_id` matches the sample columns in `data/raw_counts_combined.tsv`/`data/tpm_expression.csv` exactly; `condition` matches the comparison labels in `data/comparisons.tsv`. See `data/README.md` for the full column reference, including the binary trait columns (`is_mobile`, `is_sessile`, `is_early`, `is_late`) used for WGCNA module-trait correlations.

---

## Outputs are not distributed here

`figures/`, `photos/`, and `results/` are excluded from version control (`.gitignore`). `figures/` and `results/` hold nothing but generated outputs — every file in them is reproduced exactly by running the corresponding script in the sections above (several of the `.tiff` figures also exceed GitHub's 100 MB per-file limit on their own). `photos/` holds the 22 raw flask photographs behind the biofilm quantification; those are the authors' own bench images and aren't redistributed, so `07_biofilm_quantification/*.py` needs your own copies placed in `07_biofilm_quantification/FIJI/masks/` (see that folder's `README.md`) to reproduce that stage end-to-end. `data/` is likewise shipped as format-only examples (header + a couple of rows) rather than the real dataset -- see `data/README.md` for what each file needs replacing with and how to regenerate it.

---

## Citation

Pierry, P.M., Feitosa-Junior, O.R., Martins-Junior, J., Barbosa, D., da Silva, A.M., Zaini, P.A. "Transcriptome profiling reveals differential expression of virulence genes in *Xylella fastidiosa* under nutrient-rich and xylem-like conditions." *Pathogens* (MDPI), submitted.

## License

See `LICENSE` (if present) or contact the corresponding author.
