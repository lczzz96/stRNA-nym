# Analysis for Spatiotemporal transcriptomics of water lily (Nymphaea colorata)

This repository contains the complete analysis pipeline for the **spatial single‑cell transcriptomics** of the tropical water lily (*Nymphaea colorata*). The project aims to dissect the cellular heterogeneity, developmental trajectories, and molecular mechanisms underlying the formation of complex floral organs (sepals, petals, stamens, and carpels).

## Key features

- Single‑cell resolution combined with spatial positional information (barcodes coordinates).
- Trajectory inference using Monocle3 to reconstruct organ‑specific differentiation paths.
- Focus on MADS‑box tetramer optimisation and cross‑organ expression correlation.
- Fully reproducible R‑based workflow with shell scripts for clustering.

## Repository structure
```
st_water_lily/
├── data_for_analysis/ # Intermediate data files
│ ├── L7_petal_barcodes_pos.tsv
│ ├── L7_sepal_barcodes_pos.tsv
│ ├── L7_stamen1-4_barcodes_pos.tsv
│ └── ...
├── MADS-box_tetramer.R # MADS‑box tetramer optimisation model
├── S4_clustering.R # Clustering analysis for Stage4 samples
├── cor_nym_vs_pha.R # Cross‑species co‑expression analysis
├── cor_ot_vs_it_vs_st.R # Expression correlation across different organs
├── meristem_cell_trajectory.R # Monocle trajectory from meristem to organs
├── ot2st_trajectory.R # Trajectory from outer tepal to stamen
├── steel_clustering.sh # Cell clustering by STEEL
└── README.md
```

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

1. Initial clustering – STEEL and Seurat algorithm applied to spatial barcodes (steel_clustering.sh and S4_clustering.R).
2. Trajectory inference – Monocle2 constructs developmental trajectory (see meristem_cell_monocel.R and ot2st_trajectory.R).
3. Cross‑organ comparison – Correlation of gene expression between ovary, stamen, and pistil (cor_ot_vs_it_vs_st.R).
4. Cross‑species comparison - Correlation of organ co-expression genes between N.colorata and P.aphrodite (cor_nym_vs_pha.R).
5. Mechanistic investigation – MADS‑box expression patterns and tetramer optimisation (MADS-box_tetramer.R).
