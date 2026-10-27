# scSimEval 0.99.9

## Major Improvements

* **Interactive Step-by-Step Zoom Controls (Shiny App & R Package):**
  - Integrated an interactive step-by-step Zoom Toolbar on **every plot figure and architectural diagram** in `scSimEval Studio` (`launch_scSimEval_app()`):
    - Comparative Bubble Matrix (`plot_bubble_matrix`)
    - Radar & Evaluation Summary (`plot_eval_summary`)
    - Distribution QC Curves (`plot_dist_qc`)
    - Scalability Benchmarking (`plot_scale_bench`)
    - Metric Boxplots & Barplots (`plot_metric_box`)
    - Performance Heatmap (`plot_metric_heat`)
    - PCA Ordination (`plot_metric_pca_out`)
    - MDS Metric Space (`plot_metric_mds_out`)
    - Cell Embeddings UMAP/t-SNE/PCA (`plot_cell_embeddings`)
    - Benchmarking Workflow Diagram & Score Normalization Diagram (`workflow_diagram`, `score_normalization_workflow`)
  - Provides instant, client-side step zoom (+15% step in, -15% step out, 50% to 300% range) with a live magnification percentage badge (`100%`, `115%`, `130%`, etc.), 1-click reset, and a smooth scrollable/pannable viewport container.
  - Added and exported `zoom_plot(p, zoom = 1.0, step = 0, xlim = NULL, ylim = NULL)` in the R package to allow programmatic step-by-step zooming on any ggplot2 benchmark figure without altering underlying data or scaling factors.

* **Multi-Dataset Benchmarking & "Add to Current Benchmark Results" in Data Mode 4:**
  - Enabled the **"Add to current benchmark results"** workflow for Data Mode 4 (**"Upload Saved Results"**) in both the sidebar and the Data Hub main panel. Users can now sequentially upload multiple RDS or CSV benchmark datasets (e.g., Rice Bud, Soybean Root, PBMC) into a single unified benchmarking session.
  - Resolved a critical deduplication issue in single-cell (Mode 2) and multiomics (Mode 3) evaluation pipelines where shared simulator names across distinct datasets were inadvertently overwritten.
  - When the same simulator (e.g., `scDesign3`, `scMultiSim`, or `MOSim`) is evaluated across multiple datasets, `scSimEval` preserves each dataset's distinct records and automatically calculates the simulator's cross-dataset average performance across the 48 Gold-Standard metrics, updates overall fidelity scores and rankings in the leaderboard, and renders the combined comparative bubble matrix.
  - Added and exported `combine_benchmark_results(..., list_of_benchmarks, dataset_names, filter_gold_standard, recompute_normalized)` and its alias `append_benchmark_results()` in the R package to consolidate, standardize, and statistically aggregate multiple benchmark datasets or session files into a unified `scSimEval_multidataset_benchmark` object.

* **Streamlined 48 Gold-Standard Curated Evaluation Metrics Architecture (Option A):**
  - Following systematic sensitivity and collinearity audits, transitioned the evaluation engine to a streamlined, publication-grade **48 Gold-Standard Metrics Suite** (evaluating **39 metrics for unimodal scRNA-seq / scATAC-seq** and **48 metrics for paired single-cell multiomics co-assays** across the 8 canonical categories).
  - **Eliminated Collinear Redundancies with Empirical Proof:**
    - Pruned `ECDF_DiffArea` due to mathematical identity with 1D Wasserstein distance by Kantorovich-Rubinstein duality ($r = 1.000$).
    - Pruned redundant summary moment error scores (`MAE`, `RMSE`) in favor of outlier-resistant `MAD` and global `Wasserstein`.
    - Pruned collinear clustering partition metrics (`AMI`, `v_measure`, `fmi`, `homogeneity`, `completeness`) in favor of `ari` and `nmi`.
    - Pruned redundant `shannon_entropy` in favor of hypothesis-tested `cms` (CellMixS) and `isi` (LISI).
    - Pruned raw accuracy measures in favor of `predictive_macro_F1`, `cross_batch_F1`, and `cross_modal_F1` to eliminate majority-class distortion on rare single-cell populations.
  - Added new exported function `filter_gold_standard_metrics(df)` which dynamically filters any benchmark dataframe or session object to the 48 Gold-Standard metrics, eliminating collinear duplicates and ensuring clean backwards compatibility when loading historical RDS files.
  - Updated demo benchmark database (`demo_benchmark_data.rds`) to 48 Gold-Standard metrics (1,569 rows across 6 simulators on 2 datasets).
  - Synchronized `R/10_visualizations.R` and `inst/shiny/scSimEvalApp/10_visualizations.R` with comprehensive display labels, directional polarity rules, and category mappings.

* **In-App Live Running Logs & Diagnostics Suite:**
  - Added a dedicated **"Logs"** tab immediately after **"Help & Getting Started"** in the interactive Shiny application (`launch_scSimEval_app()`).
  - Implemented real-time log streaming that captures all execution events, data ingestion and validation steps, metric calculation progress across all 48 gold-standard benchmark metrics, visualization rendering, and file export actions.
  - Automatically captures and formats runtime `ERRORS`, `WARNINGS`, `NOTES`, and `INFO` events with timestamps, components, severity badges, and error tracebacks to simplify interactive debugging.
  - Provided interactive control tools including severity-level filtering (`All`, `ERROR`, `WARNING`, `NOTE/INFO`), real-time search filtering, auto-scroll toggle, manual log refresh, instant log clearing, and full text `.txt` download export handler for reproducibility.
  - Initialized automatic system diagnostics on startup (R version, architecture, OS platform, memory/CPU cores detection, cache status, and mounted asset paths).

* **Visual & Documentation Harmonization:**
  - Standardized Figure 2 (`score_normalization_workflow.png`) across pkgdown articles, vignette references, and Shiny Studio help pages to feature high-contrast bold circular outlines (`shape = 21`, `stroke = 1.6`) for best-performing methods, eliminating legacy square glyphs.
  - Synchronized documentation, vignette captions, and schema summaries to reflect the updated circular visual encodings.

# scSimEval 0.99.8

## Major Improvements

* **Visual Refinements & Usability Streamlining:**
  - **Comparative Bubble Matrix:** Replaced square glyphs with high-contrast bold circular outlines (`shape = 21`, `stroke = 1.6`) highlighting the best-performing simulator for each evaluation metric, eliminating visual clutter.
  - **Cell Embeddings (UMAP, t-SNE, PCA):** Fixed facet title for empirical reference to cleanly show as `"Reference"` instead of `"Reference (Reference)"`.
  - **Per-Tab Assay Modality Selectors:** Added independent Modality dropdowns (scRNA-seq, scATAC-seq, Combined) across Metric Plots, Metric Heatmap, PCA Ordination, MDS Metric Space, and Cell Embeddings.
  - **Metric Plots Clean Wordings:** Direction labels (`Higher is better` / `Lower is better`) are now displayed exclusively for raw scores, while normalized [0, 1] barplots display clean titles and standard `"Standardized Fidelity Score"` y-axis labels.
  - **Unpaired Scalability Benchmarking:** Added modality breakdown options to plot scRNA-seq, scATAC-seq, or combined runtime and peak RAM for unpaired multiomics datasets.
  - **Metric Heatmap Normalized Scores:** Ensured heatmap cell labels display normalized [0, 1] fidelity scores when relative scaling is active.
  - **Cell Embeddings (Side-by-Side Comparison):** Resolved side-by-side layout visualization where only the reference panel was displayed; enabled robust dual-facet rendering (Reference vs. Single Simulator) with dynamic simulator selection and automatic fallback.
  - **Hardware & Memory Guidelines:** Streamlined the Help page by replacing the bulky table with concise 3–4 line general guidelines for CPU and memory sizing across exploratory, large-scale multiomics, and HPC atlas-scale benchmarks.
  - **Benchmark Export Streamlining:** Streamlined Benchmark Results Export to provide the comprehensive `.rds` session object for full analytical reproducibility, and added a dedicated download handler for the Master Benchmark Scores table.
* **Data Transformation Audit & Empirical Ground-Truth Protocol:**
  - Audited and updated metric calculations to ensure proper single-cell data transformations (CPM normalization, log2(CPM+1), library-size scaling, and low-dimensional PCA projections) across all 62 evaluation measures according to benchmark literature.
  - Refined documentation and Shiny Studio descriptions to consistently reflect empirical reference-based evaluation across all 62 metrics without artificial zero-ground-truth claims.

# scSimEval 0.99.7

## Major Improvements

* **Comprehensive Metric Measurement Capacity Guidance Across All Interfaces:**
  - Embedded rigorous, multi-layered guidance documenting the exact quantitative metric measurement capacity across single-cell and multiomics experiment modes based on uploaded inputs:
    - **Single-cell count matrix alone:** Evaluates **~23 baseline metrics** spanning Category I (Distributional Properties: KS, Wasserstein, MAD, MAE, RMSE, OV, Bhattacharyya, ECDF diff, Runs test, 2D Fasano-Franceschini, 2D Peacock, KDE test, MMD, Fréchet single-cell distance), Category II (Correlations & Zeroes: BCV discrepancy, dropout curve midpoint/slope, mean-variance $R^2$, cell/gene Pearson correlations), Category VI (Trajectory Dynamics), and Category VIII (Scalability).
    - **+ Cell Types (Optional):** Adds **+12–15 metrics** unlocking Category III (Cellular Structure & Clustering / ARI, NMI, AMI, Silhouette width, Dunn index) and Category V (Biological Signal Fidelity / SimBench DE, DP, DV, top DEG Jaccard overlap, log2FC correlation).
    - **+ Batch Annotations (Optional):** Adds **+7 metrics** unlocking Category IV (Batch Mixing & Technical Confounders / Batch Silhouette, Shannon Entropy, PCR $R^2$, CMS CellMixS, ISI, Seurat mixing metric, LDF difference).
    - **Full Single-Cell Metadata:** Unlocks up to **~55 metrics** across Categories I–VI and VIII.
    - **Paired Multiomics (RNA + ATAC alone):** Evaluates **~48 metrics** across unimodal layers plus Category VII Cross-Modal Coupling (FOSCTTM cell alignment score, Match@1 nearest neighbor accuracy, cross-modal generative prediction accuracy/F1, $R_V$ matrix correlation, peak-to-gene regulatory linkage), scaling to the **full suite of 62 metrics** with metadata.
    - **Unpaired Multiomics:** Evaluates **~38–42 metrics** baseline, with cell-type labels unlocking cluster-matched label transfer and cell-type regulatory linkage (up to ~58 metrics).

* **End-to-End Analytical Workflow R Script Generators:**
  - Upgraded the standalone R script downloaders across all 8 diagnostic visualization panels and the Comparative Bubble Matrix. Rather than outputting pre-computed static plot scripts, generated scripts now provide the complete logical workflow of the analysis: reading uploaded or built-in count matrices, executing metric calculations, standardizing benchmark scores, and rendering publication-ready ggplot2 figures with custom styling and export parameters.

* **Inbuilt and Standalone Shiny Studio Ecosystem Alignment:**
  - Fully synchronized the package's inbuilt Shiny application (`launch_scSimEval_app()`) and the standalone deployment repository (`scSimEval_shiny`) to complete feature parity, including the four-mode Data Hub, confirmation buttons, app-wide loading indicators, and HPC batch execution script generators.
  - Retained the standalone application's Google Analytics 4 pageview counter with persistent container-lifetime session tracking initialized to 0 views.

# scSimEval 0.99.6

## Major Improvements

* **Multi-Dataset Options & Standalone Reproducible R Scripts Across All Visualizations:**
  - Extended multi-dataset selection options to all 8 diagnostic visualization panels in the Visualizations tab (Evaluation Summary, Distribution QC, Scalability Benchmark, Metric Boxplots, Metric Heatmap, PCA Ordination, MDS Metric Space, and Cell Embeddings). Users can seamlessly toggle between **All Datasets Combined (Cross-Dataset Average)** and individual datasets (`Dataset 1`, `Dataset 2`, ..., `Dataset N`) with live active dataset badge synchronization.
  - Added dedicated **Download Standalone R Script** (`download_*_r_script`) buttons to all 8 visualization subpanels and the Comparative Bubble Matrix. These generate 100% self-contained, reproducible R scripts with embedded benchmark data, user choices, and publication-ready ggplot2 code that can be executed independently outside of Shiny.

* **Standardization to MiB for Peak RAM & Dual-Value Reporting in Master Benchmark Scores Table:**
  - Standardized computational scalability memory measurement to **MiB RAM** ($1\text{ MiB} = 1,048,576\text{ bytes}$) across all scalability plots, selectors, table wordings, confirmation dialogs, and documentation, resolving ambiguity with decimal megabytes (MB).
  - Enhanced the Master Benchmark Scores Table to explicitly report both unnormalized raw measurements (`Score_Raw` in native statistical and physical units: KS distance, runtime in seconds, peak RAM in MiB, Chi-Square values, etc.) alongside direction-aware normalized scores (`Score` in $[0, 1]$ where 1.0 represents best performance), equipped with per-column dropdown filtering.

* **Modernized Help & Getting Started Tab and Ecosystem-Wide Documentation Alignment:**
  - Fully restructured the in-app **Help & Getting Started** tab to reflect the updated **8-tab architecture** (Home, Data Hub, Leaderboard & Simulator Rankings, Comparative Bubble Matrix, Visualizations, Benchmark Results & Export, Help & Getting Started, and Contact).
  - Synchronized all documentation across the R package vignettes, manuals, README, and the GitHub Pages site (https://kabilanbio.github.io/scSimEval), explicitly clarifying the distinction between the lightweight command-line testing datasets (`example_scrna`, `example_scatac`, `example_multiomics`) and the rich pre-computed multi-dataset benchmark suite bundled in the Shiny application.


# scSimEval 0.99.5

## Major Improvements

* **Multi-Dataset Benchmarking Study Support (e.g., 16 Datasets × 6 Simulators):**
  - Added dedicated support for multi-dataset benchmark studies (e.g., evaluating 6 single-cell simulators across 16 diverse experimental datasets or cell types, yielding 5,952 metric observations across all 62 canonical measures).
  - Implemented `benchmark_simulation_suite()`: high-performance R function engineered for beginner users and high-throughput HPC clusters with automatic sparse matrix coercion (`dgCMatrix`), proactive garbage collection (`gc()`) per dataset to prevent RAM accumulation, and multi-threading.
  - Automatically calculates both a consolidated granular benchmark table (`granular_benchmark_table`) and a cross-dataset averaged summary table (`benchmark_summary_table`).

* **Cross-Dataset Averaged Comparative Bubble Matrix:**
  - In `plot_benchmark_bubble_matrix()` and Shiny Studio Tab 4 (Comparative Bubble Matrix), multi-dataset evaluations automatically average fidelity scores across all evaluated datasets for each simulator.
  - Keeps the flagship bubble plot uncluttered by displaying strictly simulators on the axis (no dataset names on axes), with a dynamic subtitle indicating the number of datasets averaged.
  - Extended cross-dataset averaging to heatmaps, PCA, and MDS ordination plots.

* **Streamlined Data Hub Ingestion (Four Canonical Modes):**
  - Redesigned the Data Hub tab to offer strictly four clear, uncluttered dataset options:
    1. **Example dataset analysis**: Instant 1-click loading of pre-computed benchmark results across 6 simulators (Splatter, scDesign3, SCRIP, SymSim, dyngen, simATAC) evaluated on all 62 canonical measures.
    2. **Single-Cell (scRNA or scATAC)**: Evaluates 1 or multiple original reference datasets (e.g., up to 50 datasets) against multiple simulators with per-simulator simulated dataset uploads, tracking runtime and peak RAM.
    3. **Multiomics (scRNA + scATAC)**: Evaluates 1 or multiple original reference multiomics datasets (RNA + ATAC, paired or unpaired) against multiple simulators with per-simulator simulated RNA and ATAC dataset uploads.
    4. **Upload Saved Results (.rds / .csv)**: Instant loading of previously computed benchmark files (.rds or .csv) from scSimEval or HPC batch runs, supporting both single-dataset and multi-dataset benchmarks.
  - Interactive Dataset Preview in Subtab 3 allowing toggling between overall cross-dataset averages and granular dataset tables, with direct CSV export.
  - Direct matrix inspection in Subtab 2 supporting browsing through all loaded reference and simulated count matrices.

* **Minimalist & Streamlined Shiny Studio Interface:**
  - Modernized the web interface across all tabs with a clean, clutter-free scientific layout, removing redundant controls in the Home, Data Hub, and Bubble Matrix panels.
  - Simplified comparative metric filtering into two intuitive modes: **All 62 Measures** and **Custom Selection** (with direct category dropdown and measure toggles).
  - Streamlined Data Hub ingestion controls and enhanced navigation linking directly to the methodology guide.

* **Refined Multiomics Scope & Natural Language:**
  - Clearly unified application terminology to emphasize simulation evaluation across single-cell RNA-seq (scRNA-seq), single-cell ATAC-seq (scATAC-seq), and single-cell multiomics (scRNA-seq + scATAC-seq) against empirical references across 62 evaluation measures.
  - Replaced complex technical jargon with direct, human-friendly scientific explanations across all tabs.

* **Methodology Guide & Normalization Pipeline:**
  - Integrated the comprehensive **Score Normalization Pipeline** workflow figure (Figure 2) and structured methodology cards detailing directional transformation, min-max feature scaling, outlier winsorization, and composite rankings.

* **Package & Visualization Synchronization:**
  - Fully synchronized standalone Shiny Studio visualization components with core package routines, ensuring complete integration of low-dimensional cell embedding algorithms (`compute_dataset_embeddings()`, `plot_dataset_embeddings()`, `compute_embedding_quality_metrics()`).

* **Example Metadata Inputs & Data Format Documentation:**
  - Added export of example metadata datasets `example_cell_types` and `example_batch_info` matching the package's benchmark reference matrices (`example_scrna`).
  - Added interactive metadata format modals and direct sample CSV template downloads (`example_cell_types.csv`, `example_batch_annotations.csv`) in Data Hub (Modes 2 & 3) and Help & Getting Started (Section 5) of scSimEval Studio.
  - Enhanced in-app parser `read_uploaded_labels()` to automatically parse both 1-column vectors and standard 2-column tables (`Cell_ID`, `Cell_Type` or `Batch`) exported from Seurat and SingleCellExperiment.
  - Added dedicated documentation and code recipes in vignettes and help panels detailing exact count matrix, cell-type, and batch input specifications.


# scSimEval 0.99.4

## Bug Fixes & Improvements

* **`compute_method_leaderboard()` & `.ingest_bubble_data()` Robustness:**
  - Added automatic fallback to `Method = "Simulation"` and `Category = "Uncategorized"` when evaluating single simulation accuracy results (`evaluate_simulation_accuracy()`) or pre-formed benchmark data frames lacking explicit method labels.
  - Resolved `R CMD check` example failure in `man/compute_method_leaderboard.Rd`.

* **Expanded Upload Size Limit (10 GB):**
  - Increased `shiny.maxRequestSize` from 500 MB to 10 GB (`10 * 1024^3` bytes) in `inst/shiny/scSimEvalApp/app.R` to accommodate large-scale single-cell and single-nucleus multiomics count matrices.
  - Added configurable `maxRequestSize` argument to `launch_scSimEval_app()` in `R/11_shiny_app.R` (defaulting to 10 GB).

* **Organism-Agnostic Metric Harmonization:**
  - Removed `calc_cell_cycle_phase_fidelity()` which was reliant on species-specific cell-cycle phase definitions and marker genes (primarily human and mouse). All 62 core evaluation metrics across the 8 categories are verified to operate purely on mathematical, distributional, and machine-learning properties of single-cell matrices without species or organism-specific priors, guaranteeing universal applicability to plant (e.g. soybean, rice, Arabidopsis) and animal datasets alike.

* **Dependency & Package Best Practices:**
  - Added `irlba` to `Suggests:` in `DESCRIPTION` to properly declare conditional SVD acceleration in `compute_dataset_embeddings()`.
  - Replaced direct `set.seed()` invocation with RNG-safe execution that preserves and restores `.Random.seed` on function exit in accordance with Bioconductor guidelines.
  - Removed global `eval = FALSE` in `vignettes/shiny-app.Rmd` to satisfy Bioconductor vignette check requirements.


# scSimEval 0.99.3

## New Features & Enhancements

* **Paired, Unpaired & Mosaic Multiomics Framework:**
  - Added explicit `pairing = c("paired", "unpaired")` parameter (with `is_paired` alias) to `evaluate_multiomics_accuracy()` and `evaluate_multiple_datasets()`.
  - In **Paired mode** (default; e.g. 10x Chromium Multiome, SHARE-seq), all 62 measures across all 8 canonical categories are evaluated, including cell-level pairing metrics (FOSCTTM, Match@1, in silico cross-modal generation, direct peak-to-gene linkage).
  - In **Unpaired mode** (separate cells from the same tissue/condition), cell-level pairing metrics that require 1-to-1 matching cell barcodes are safely omitted with clear informational logging, while population-level cross-modal metrics and unimodal evaluations are fully computed and standardized.
  - Added comprehensive methodological guidance for decomposing mosaic multiomics datasets into paired and unpaired evaluation blocks.
  - Integrated dynamic **Multiomics Dataset Type** selector (radio buttons) in Shiny Studio Tab 2 (Data Hub: Mode 3) with live contextual warning/info banners.

* **Sub-panel 8: Cell Embeddings (UMAP, t-SNE, PCA) & Quality Metrics:**
  - Implemented high-level dimensionality reduction and visualization functions: `compute_dataset_embeddings()`, `plot_dataset_embeddings()`, and `compute_embedding_quality_metrics()`.
  - Interactive layout options: Multi-Simulator Faceted Grid and Direct 1-to-1 Comparison (Reference vs. Selected Simulator).
  - Flexible cell coloring: Cell Type labels, Unsupervised cluster recovery, Library size gradient, and Detected features gradient.
  - Integrated interactive summary table of quantitative embedding fidelity metrics (Mean Silhouette width, Adjusted Rand Index cluster concordance, and library size discrepancy %).

* **Comprehensive All-in-One Benchmark Archive (.zip) & Export Suite:**
  - Overhauled the complete benchmark bundle (`download_complete_zip`) in Shiny Studio Tab 5 ("Download Results") to deliver all quantitative results, spreadsheets, reports, and figures in an organized hierarchy.
  - **All Metrics in Multiple Formats:** Bundles complete evaluation results in Excel (`.xlsx`), CSV (`.csv`), tab-delimited text (`.txt`), and native R object (`.rds`) formats in both root and `metrics/` subfolders.
  - **Rankings & Metadata Tables:** Automatically generates and packages `method_rankings_leaderboard` (`.csv` & `.txt`) and `dataset_properties_summary` (`.csv` & `.txt`).
  - **Multi-Page Compiled PDF Report:** Compiles all 9 diagnostic and comparative figures into a single publication-quality vector PDF (`scSimEval_all_plots_report.pdf`).
  - **All Kinds of Figures (Grouped & Individual):**
    - `figures/grouped/`: Both 600 DPI publication-grade JPEGs and vectorized PDFs for comparative bubble matrix, overall evaluation summary, scalability benchmark, metric boxplots by category, performance heatmap, PCA simulator ordination, MDS metric space, comparative distribution QC, and cell embeddings comparison grid.
    - `figures/individual/`: Individual publication-grade barplots for every evaluated metric (`metric_bar_<MetricName>.jpeg` & `.pdf`) and simulator-specific 1-to-1 comparison cell embeddings against empirical reference (`cell_embeddings_compare_<SimulatorName>.jpeg` & `.pdf`).
  - Added dedicated **Download TXT Table (.txt)** action button to Card 2 in Tab 5 for one-click tab-delimited exports.

* **Documentation, Search & Site Integration:**
  - Enabled client-side Fuse.js full-text search across the GitHub Pages documentation site, allowing instant search of all functions, parameters, vignettes, and metric descriptions.
  - Completely rewritten `vignettes/shiny-app.Rmd` documenting all 6 navigation tabs, 8 visualization panels, UI controls, export options, and an R-API cross-reference table.
  - Integrated official documentation links and quick-launch actions into Tab 6 ("Help & Getting Started").

* **Dedicated Individual Figure Export (Shiny Tab 5 - Card 4):**
  - Added on-demand figure export suite in Tab 5 ("Download Results") allowing researchers to download any single diagnostic or comparative visualization in either publication-grade vector PDF (`.pdf`) or ultra-high-resolution 600 DPI JPEG (`.jpeg`).
  - Covers all 14 figure types: Comparative Bubble Matrix, Overall Evaluation Summary, Scalability Benchmark, Metric Boxplots, Performance Heatmap, PCA Ordination, MDS Metric Space, Comparative Distribution QC, UMAP Grid, t-SNE Grid, PCA Grid, and 1-to-1 simulator vs reference comparisons.
  - Dynamically populated simulator selector for custom 1-to-1 side-by-side comparison figure downloads.

* **Expanded Embeddings in Complete Benchmark Archive (`download_complete_zip`):**
  - Upgraded the comprehensive `.zip` bundle to systematically generate and organize cell embedding figures for all three supported dimensionality reduction techniques (UMAP, t-SNE, and PCA) in `figures/grouped/` (`09_cell_embeddings_umap_grid`, `10_cell_embeddings_tsne_grid`, `11_cell_embeddings_pca_grid`).
  - Systematically exports 1-to-1 comparison plots for every simulator against the empirical reference across UMAP, t-SNE, and PCA to `figures/individual/`.

* **Documentation & Vignette Synchronization:**
  - Fully updated `vignettes/shiny-app.Rmd` to document the new Card 4 export controls and the expanded zip bundle directory structure.

## Bug Fixes

* **`plot_metric_pca()` and `plot_metric_mds()` — `category = "all"` support:**
  - Fixed category filter edge case where passing `category = "all"` or `NULL` now correctly retains all metrics across all categories without filtering out rows.
* **Self-Contained Shiny Studio Runtime:**
  - Embedded runtime definitions directly within `inst/shiny/scSimEvalApp/app.R` to ensure seamless execution across any R environment without requiring package re-installation.


# scSimEval 0.99.2

## Bug Fixes
* Fixed the Bioconductor staging build failure (`Error in x[clustering == i, ]
  : (subscript) logical subscript too long`): corrected all roxygen `@examples`
  blocks in the clustering/batch/DE/trajectory/multiomics metric functions so
  cluster labels align with matrix dimensions and all required arguments are
  supplied. All examples now pass `R CMD check`.
* Removed the example from the internal helper `calc_expected_mi`, which errored
  in an installed-namespace context ("could not find function").
* Cleaned up Rd generation: `\figure` widths now declared in pixels, explicit
  `@title` tags for distribution metrics, ASCII-only documentation text, and
  regenerated `man/` pages.
* Labeled all vignette code chunks and added `sessionInfo()` blocks to every
  vignette as required by BiocCheck.

## Documentation

* **Paired / Unpaired / Mosaic Multiomics:** Added comprehensive, scientifically
  accurate documentation throughout the package clarifying which evaluation metrics
  are valid for each multiomics experimental design:
  - **Paired** (10x Multiome, SHARE-seq, SNARE-seq): Full 62-metric evaluation
    including all Category 7 cross-modal coupling metrics (FOSCTTM, Match@1,
    cross-modal generation fidelity, peak-to-gene linkage).
  - **Unpaired** (independent scRNA-seq + scATAC-seq from separate cells): Full
    unimodal evaluation (Categories 1-6, 8) plus population-level cross-modal
    metrics (co-expression module fidelity, network Jaccard, cross-modal label
    transfer). Cell-pairing metrics (FOSCTTM, cross-modal generation) are not
    applicable and excluded.
  - **Mosaic** (partial co-measurement, e.g. DOGMA-seq): Per-modality unimodal
    evaluation on the full matrix; paired coupling metrics applied only to the
    co-assayed cell subset.
* Updated `DESCRIPTION`, `README.md`, `vignettes/demo-multiomics.Rmd`, and
  `vignettes/scSimEval-workflow.Rmd` to reflect these distinctions.

# scSimEval 0.99.1

## Bug Fixes
* Fixed R-universe / Bioconductor build failure: removed `^doc$` and `^inst/doc$`
  from `.Rbuildignore` so that pre-rendered vignette HTML files are correctly
  included in the source tarball and `R CMD INSTALL --html` succeeds.
* Cleaned up `docs/` folder by removing scratch files (`git_help.Rmd`),
  redundant source Markdown files, and the empty `tutorials/` directory.

# scSimEval 0.99.0

## New Features
* Initial submission to Bioconductor.
* Unified benchmarking and fidelity assessment for single-cell and multiomics simulation algorithms (scRNA-seq, scATAC-seq, and paired multiomics).
* 62 quantitative fidelity metrics spanning 8 biological and computational categories:
  - (I) Distributional Properties
  - (II) Correlations & Zero-Inflation
  - (III) Cellular Structure & Concordance
  - (IV) Batch Effects & Confounder Mixing
  - (V) Biological Signal & Downstream Fidelity
  - (VI) Trajectory & Lineage Dynamics
  - (VII) Cross-Modal Coupling & Modularity
  - (VIII) Computational Scalability
* Comprehensive visualization suite with publication-ready comparative bubble matrices, individual metric barplots, heatmaps, PCA, and MDS projections.
* Interactive Shiny Studio (`launch_scSimEval_app()`) for real-time benchmark exploration, custom data upload, and 600 DPI publication exports.
