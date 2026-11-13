# scSimEval 0.99.10

* Updated the data input option which align with Seurat and Single cell experiment objects. 


# scSimEval 0.99.9

* **Category VII Trajectory Standardization**: Upgraded Category VII to the `dynverse` trajectory inference protocol (`him`, `edge_flip`, `geodesic_correlation`, `F1_branches`, and `wcor_features`).
* **Category VIII Scalability Normalization**: Standardized computational efficiency to core-normalized metrics (`core_seconds` and `memory_per_core_mib`).
* **S4 Object Standardization**: Standardized input handling to require `Seurat` or `SingleCellExperiment` objects with precomputed dimensional reductions.
* **Category VI Integration Metrics**: Standardized multiomics integration to 5 core metrics (`batch_asw`, `cell_type_asw`, `graph_connectivity`, `ari`, and `nmi`).
* **Category V Pseudobulk DE/DA**: Structured Category V into scRNA-seq DE and scATAC-seq DA sub-categories using biological replicate pseudobulk testing via `edgeR-QLF`.
* **Multi-Dataset Benchmarking**: Added `combine_benchmark_results()` to aggregate and average performance across multiple experimental datasets.
* **Interactive Diagnostic Logs**: Added real-time log streaming and diagnostic console in the Shiny application.

# scSimEval 0.99.8

* **Bubble Matrix Visual Refinements**: High-contrast bold circular outlines highlight best-performing simulators.
* **Modality Selectors**: Independent assay modality filters (scRNA-seq, scATAC-seq, Combined) across visualization tabs.
* **Export Enhancements**: Standalone reproducible R script generators for all diagnostic plots.

# scSimEval 0.99.7

* **Input Capacity Guidance**: Context-aware guidance detailing available metrics based on uploaded matrix, cell-type, and batch annotations.
* **Reproducible Script Downloads**: 1-click generation of self-contained ggplot2 R scripts for each visualization panel.

# scSimEval 0.99.6

* **Cross-Dataset Averaging**: Seamless toggling between cross-dataset averages and individual datasets across all visualization panels.
* **Memory Standardization**: Standardized computational memory reporting to binary MiB RAM ($1\text{ MiB} = 1,048,576\text{ bytes}$).

# scSimEval 0.99.5

* **Multi-Dataset Benchmark Suite**: Added `benchmark_simulation_suite()` for high-throughput multi-dataset evaluations.
* **Redesigned Data Hub**: Four streamlined workflow modes: Example Datasets, Single-Cell Upload, Multiomics Upload, and Saved Results.

# scSimEval 0.99.4

* **Large Matrix Uploads**: Increased Shiny upload limit to 10 GB (`shiny.maxRequestSize`).
* **Organism-Agnostic Metrics**: Removed species-specific priors to ensure universal applicability to plant and animal datasets.

# scSimEval 0.99.3

* **Paired & Unpaired Multiomics**: Dedicated evaluation workflows for paired (co-assay) and unpaired single-cell multiomics.
* **Cell Embeddings**: Added UMAP, t-SNE, and PCA visualization panels with quantitative embedding quality metrics.
* **Complete Benchmark Archive**: Multi-format ZIP export containing CSV tables, Excel workbooks, RDS objects, vector PDFs, and 600 DPI figures.

# scSimEval 0.99.2

* **Multiomics Methodology**: Documented metric applicability across paired, unpaired, and mosaic multiomics experimental designs.
* **Package Check Fixes**: Corrected documentation examples and vignette chunk settings for Bioconductor compliance.

# scSimEval 0.99.1

* **Documentation & Build Fixes**: Cleaned up documentation files and build ignore paths for R-universe / Bioconductor staging.

# scSimEval 0.99.0

* **Initial Submission**: Initial submission of `scSimEval` to Bioconductor for benchmarking single-cell and multiomics simulations across 8 biological and computational categories.
