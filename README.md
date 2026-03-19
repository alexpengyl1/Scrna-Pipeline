# scRNA-seq Analysis Pipeline in R

This repository contains a collection of R scripts for common single-cell RNA-seq analyses, including Seurat-based preprocessing and integration, cell type annotation, differential expression, pseudobulk analysis, pathway scoring, trajectory inference, cell-cell communication, regulatory network inference, CNV inference, and RNA velocity.

The repository has been renamed and documented to improve readability, GitHub usability, and FAIR alignment.

## Repository scope

The scripts in this repository are examples and analysis templates collected across multiple single-cell projects. They are most useful as modular workflow references rather than as a single end-to-end pipeline that can be run unchanged on any dataset.

## FAIR-oriented design choices

This repository has been improved to better support FAIR principles:

- `Findable`: script names are descriptive, ordered, and standardized.
- `Accessible`: files are stored in plain-text R or notebook formats suitable for Git-based sharing.
- `Interoperable`: common R/Bioconductor tools are used, and scripts rely on broadly used data structures such as Seurat objects and sparse count matrices.
- `Reusable`: this README documents the purpose, expected inputs, outputs, dependencies, and limitations of the scripts.

## Current limitations

Several scripts still contain project-specific assumptions that should be adapted before reuse:

- hard-coded `setwd()` paths pointing to local Windows or macOS directories
- hard-coded object names such as `scRNA1`, `scRNA_harmony`, `ifnb`, or `sc.combined`
- package installation commands embedded inside analysis scripts
- assumptions about metadata column names such as `orig.ident`, `celltype`, `seurat_clusters`, and custom cluster labels
- references to local reference files such as `ref_Human_all.RData`, `ref_Mouse_all.RData`, `.rds` resources, or GMT files not included in this repository

For FAIR reuse, these paths and object-specific assumptions should ideally be parameterized in a future cleanup pass.

## Suggested repository structure

The repository now stores analysis scripts in a dedicated `scripts/` directory. A future improvement would be to further separate reference data, example outputs, and environment files into dedicated folders such as:

```text
Scrna-Pipeline/
├── README.md
├── scripts/
├── references/
├── examples/
├── env/
└── docs/
```

## Script index

### Core preprocessing and integration

- `scripts/01_seurat_integration_basics.R`
  Basic Seurat workflow for reading 10x data, QC, integration, dimensionality reduction, clustering, and visualization.

- `scripts/02_qc_and_doublet_removal.R`
  Single-sample and multi-sample doublet detection workflows using DoubletFinder, including object slimming with `DietSeurat()`.

- `scripts/03_integration_and_singler_annotation.R`
  Multi-sample Seurat integration followed by SingleR-based cell type annotation.

- `scripts/17_seurat_v4_v5_conversion.R`
  Minimal examples showing how to convert Seurat assays between v3/v4 and v5 assay classes.

### Exploratory analysis and differential expression

- `scripts/04_3d_pca_plot.R`
  3D PCA visualization for averaged expression profiles across samples and cell types.

- `scripts/05_differential_expression_analysis.R`
  Differential expression analysis between selected groups using Seurat `FindMarkers()` and annotated scatter plots.

- `scripts/06_pseudobulk_analysis.R`
  Pseudobulk differential expression using `AggregateExpression()` and `DESeq2`, with comparison against single-cell DE results.

- `scripts/07_multisample_volcano_plot.R`
  Multi-sample differential expression and volcano-style comparison workflow.

- `scripts/08_enrichment_and_signature_scoring.R`
  AUCell and msigdbr-based pathway or signature scoring.

### Trajectory and dynamic state analysis

- `scripts/09_monocle2_trajectory_analysis.R`
  Monocle2 workflow for cell ordering and pseudotime trajectory analysis.

- `scripts/10_monocle3_basic_analysis.R`
  Introductory Monocle3 workflow covering preprocessing, reduction, and visualization.

- `scripts/11_monocle3_pseudotime_analysis.R`
  A more advanced Monocle3 workflow including graph learning, pseudotime ordering, and reuse of Seurat UMAP embeddings.

- `scripts/16_velocyto_analysis.R`
  RNA velocity analysis using `velocyto.R`, Seurat, and SeuratWrappers.

### Cell-cell communication and regulatory analysis

- `scripts/12_cellchat_analysis.R`
  CellChat-based inference of intercellular communication.

- `scripts/13_nichenet_analysis.R`
  NicheNet workflow for ligand-target prioritization and sender-receiver analysis.

- `scripts/14_scenic_analysis.R`
  SCENIC workflow for transcription factor regulatory network inference.

- `scripts/15_infercnv_analysis.R`
  inferCNV workflow for copy number variation inference in single-cell data.

## Software dependencies

Most scripts rely on the following R ecosystem packages:

- `Seurat`
- `tidyverse`
- `dplyr`
- `patchwork`
- `ggplot2`
- `SingleR`
- `harmony`

Specialized workflows additionally use packages such as:

- `DoubletFinder`
- `AUCell`
- `clusterProfiler`
- `msigdbr`
- `monocle`
- `monocle3`
- `CellChat`
- `nichenetr`
- `SCENIC`
- `infercnv`
- `velocyto.R`
- `SeuratWrappers`

## Expected inputs

Depending on the script, expected inputs include one or more of the following:

- 10x Genomics matrix directories readable by `Read10X()`
- serialized R objects such as `.Rdata` or `.rdata`
- Seurat objects already containing normalized data, metadata, embeddings, and cluster labels
- reference annotation files for SingleR
- ligand-receptor or pathway reference resources such as `.rds` or `.gmt` files
- gene position annotation files for inferCNV

## Expected outputs

Common outputs include:

- QC plots and dimensionality reduction plots
- cluster or cell type annotations stored in Seurat metadata
- differential expression tables
- pseudobulk comparisons
- enrichment or AUCell scores
- trajectory and pseudotime plots
- cell-cell communication summaries
- CNV, regulon, or velocity outputs

The scripts currently save outputs inconsistently. A future FAIR improvement would be to standardize output directories such as `results/figures/` and `results/tables/`.

## Reuse recommendations

Before reusing any script:

1. replace `setwd()` calls with project-relative paths
2. move package installation commands to a setup script or environment file
3. document input files and expected metadata columns for each workflow
4. standardize object names across scripts
5. pin package versions using `renv` or another environment manager
6. add small example datasets or mock inputs where licensing permits

## Recommended next improvements

To make this repository more FAIR and publication-ready, the most valuable next steps would be:

1. add an environment lock file such as `renv.lock`
2. add a `LICENSE`
3. add per-script header metadata with author, purpose, required inputs, and outputs
4. remove embedded installation commands from analysis scripts
5. create a minimal reproducible example for at least one workflow
6. separate reusable utility functions from project-specific analysis notebooks

## Citation and attribution

If you use this repository in teaching, collaboration, or downstream analysis, please cite the underlying tools and databases used by each workflow, including Seurat, SingleR, Monocle, CellChat, NicheNet, SCENIC, inferCNV, and velocyto where relevant.

## Original notes preserved from the previous README

- For `FindAllMarkers`, consider using Presto for faster Wilcoxon and auROC calculations: https://github.com/immunogenomics/presto
- For doublet removal, consider using `scDblFinder` as an alternative to DoubletFinder: https://biostatsquid.com/scdblfinder-tutorial/
