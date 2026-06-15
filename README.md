# Analysis for Spatiotemporal transcriptomics of water lily (Nymphaea colorata)

[![R 99.2%](https://img.shields.io/badge/R-99.2%25-276DC3?logo=r&logoColor=white)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the complete analysis pipeline for the **spatial single‑cell transcriptomics** of the tropical water lily (*Nymphaea colorata*). The project aims to dissect the cellular heterogeneity, developmental trajectories, and molecular mechanisms underlying the formation of complex floral organs (sepals, petals, stamens, and carpels).

## Objectives

- **Cellular heterogeneity** – Identify and characterise distinct cell types within each floral organ using clustering analysis.
- **Developmental trajectories** – Reconstruct the differentiation paths from meristematic cells to mature floral organs.
- **Molecular mechanisms** – Investigate the expression dynamics of key transcription factors (e.g., MADS‑box genes) and their role in floral organ identity determination.

## Key features

- Single‑cell resolution combined with spatial positional information (barcodes coordinates).
- Trajectory inference using Monocle3 to reconstruct organ‑specific differentiation paths.
- Focus on MADS‑box tetramer optimisation and cross‑organ expression correlation.
- Fully reproducible R‑based workflow with shell scripts for clustering.

## Repository structure
st_water_lily/
├── data_for_analysis/ # Intermediate data files (spatial barcodes, etc.)
│ ├── L7_petal_barcodes_pos.tsv
│ ├── L7_sepal_barcodes_pos.tsv
│ ├── L7_stamen1-4_barcodes_pos.tsv
│ └── ...
├── MADS-box_tetramer.R # MADS‑box tetramer optimisation model
├── S4_clustering.R # Clustering analysis for carpel (S4)
├── cor_nym_vs_pha.R # Cross‑species MADS‑box co‑expression analysis
├── cor_ot_vs_it_vs_st.R # Expression correlation across different organs
├── meristem_cell_monocel.R # Monocle3 trajectory from meristem to organs
├── ot2st_trajectory.R # Trajectory from ovary to stamen
├── steel_clustering.sh # Shell script calling the STEEL clustering tool
└── README.md # This file


## Data sources

- **Expression matrix** – All spatial transcriptomic data are publicly available at:  
  [http://osf.io/m68cn/overview](http://osf.io/m68cn/overview)
- **Intermediate analysis files** – Pre‑processed barcode position tables and other intermediate files are stored in the [`data_for_analysis/`](data_for_analysis) folder of this repository.

## Getting started

### Prerequisites

- R (≥ 4.2) with the following core packages:
  ```r
  install.packages(c("Seurat", "monocle", "dplyr", "ggplot2", "ggsci", "clustree", "reshape2", "pheatmap"))
- STEEL (Spatial Transcriptome based cEll typE cLustering)

STEEL is an unsupervised manifold learning algorithm designed for spatial transcriptome data analysis.

- **Project homepage**: [http://steel-st.sourceforge.io](http://steel-st.sourceforge.io)
- **Download**: Get the latest version (source code and precompiled binaries for Linux/macOS) from its [SourceForge page](https://sourceforge.net/projects/steel-st/).
- **Installation**:
  - Compile from source:
    ```bash
    g++ src/STEEL.cpp -o steel -O3

### Analysis workflow

1. Initial clustering – STEEL algorithm applied to spatial barcodes (steel_clustering.sh).
2. Cell type annotation – Based on known marker genes and spatial coordinates.
3. Trajectory inference – Monocle3 constructs developmental paths (see meristem_cell_monocel.R and ot2st_trajectory.R).
4. Cross‑organ comparison – Correlation of gene expression between ovary, stamen, and pistil (cor_ot_vs_it_vs_st.R).
5. Mechanistic investigation – MADS‑box expression patterns and tetramer optimisation (MADS-box_tetramer.R).
