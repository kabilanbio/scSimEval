# Demonstration 4: Analyzing Pre-Saved Benchmark Results (.rds)

## Introduction

In large-scale benchmarking studies involving multiple simulation
methods, diverse biological datasets, or multi-core high-performance
computing (HPC) clusters, running the full 56 curated evaluation metric
suite can take considerable time. Once an evaluation is complete,
researchers typically serialize the output using
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) for long-term
storage, publication sharing, and reproducible downstream analysis.

**`scSimEval`** is architected to allow instant reloading, re-analysis,
leaderboard ranking, and publication-grade visualization of any
pre-saved `.rds` benchmark file — **without re-computing a single
metric**.

This tutorial demonstrates how to load, inspect, rank, and visualize
pre-saved benchmark results in R console sessions and prepare them for
interactive exploration.

------------------------------------------------------------------------

## 1. Loading a Pre-Saved Benchmark Object

`scSimEval` bundles a pre-computed multi-dataset benchmark archive in
its package directory. You can load it directly via base R’s
[`readRDS()`](https://rdrr.io/r/base/readRDS.html):

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`scSimEval`](https://kabilanbio.github.io/scSimEval/)`)`\
\
`# Locate the bundled benchmark file or supply your own path`\
`rds_path`` ``<-`` `[`system.file`](https://rdrr.io/r/base/system.file.html)`(``"shiny"``, ``"scSimEvalApp"``, ``"data"``, ``"demo_multidataset_benchmark.rds"``, package ``=`` ``"scSimEval"``)`\
`if`` ``(``rds_path`` ``==`` ``""``)`` ``{`\
`  ``rds_path`` ``<-`` `[`file.path`](https://rdrr.io/r/base/file.path.html)`(``".."``, ``"inst"``, ``"shiny"``, ``"scSimEvalApp"``, ``"data"``, ``"demo_multidataset_benchmark.rds"``)`\
`}`\
\
`# Ingest the serialized benchmark object`\
`bench_data`` ``<-`` `[`readRDS`](https://rdrr.io/r/base/readRDS.html)`(``rds_path``)`\
\
`# Verify class and contents`\
[`class`](https://rdrr.io/r/base/class.html)`(``bench_data``)`\
`#> [1] "scSimEval_multidataset_benchmark" "scSimEval_consolidated"          `\
`#> [3] "list"`\
[`names`](https://rdrr.io/r/base/names.html)`(``bench_data``)`\
`#>  [1] "benchmark_df"             "granular_benchmark_table"`\
`#>  [3] "benchmark_summary_table"  "cross_dataset_stats"     `\
`#>  [5] "dataset_summary"          "dataset_summary_df"      `\
`#>  [7] "methods"                  "datasets"                `\
`#>  [9] "is_multidataset"          "n_datasets"              `\
`#> [11] "n_simulators"`

The loaded object contains: \* `benchmark_df`: The complete tidy
evaluations table across all simulators and metrics. \*
`benchmark_summary_table`: Formatted scores ready for reporting. \*
`methods`: Character vector of evaluated simulators (e.g., `scDesign3`,
`SCRIP`, `dyngen`, `Splatter`, `SymSim`, `simATAC`). \* `datasets`:
Evaluated datasets (e.g., `Rice Bud (Paired)`,
`Soybean Root (Unpaired)`).

------------------------------------------------------------------------

## 2. Inspecting the Benchmark Structure

The core of every `scSimEval` benchmark object is a tidy data frame
following standard tidyverse conventions:

\
`# Preview top rows of the benchmark scores`\
[`head`](https://rdrr.io/r/utils/head.html)`(``bench_data``$``benchmark_df``[``, `[`c`](https://rdrr.io/r/base/c.html)`(``"Dataset"``, ``"Method"``, ``"Category"``, ``"Metric"``, ``"Score_Raw"``, ``"Score"``)``]``, ``8``)`\
`#>             Dataset     Method                      Category      Metric`\
`#> 1 Rice Bud (Paired)  scDesign3 (I) Distributional Properties Wasserstein`\
`#> 2 Rice Bud (Paired) scMultiSim (I) Distributional Properties Wasserstein`\
`#> 3 Rice Bud (Paired)    scMoMtF (I) Distributional Properties Wasserstein`\
`#> 4 Rice Bud (Paired)    Matilda (I) Distributional Properties Wasserstein`\
`#> 5 Rice Bud (Paired)      MOSim (I) Distributional Properties Wasserstein`\
`#> 6 Rice Bud (Paired)    scCross (I) Distributional Properties Wasserstein`\
`#> 7 Rice Bud (Paired)  scDesign3 (I) Distributional Properties          KS`\
`#> 8 Rice Bud (Paired) scMultiSim (I) Distributional Properties          KS`\
`#>   Score_Raw  Score`\
`#> 1    0.0842 1.0000`\
`#> 2    0.1014 0.8289`\
`#> 3    0.1411 0.4334`\
`#> 4    0.1593 0.2522`\
`#> 5    0.1847 0.0500`\
`#> 6    0.1299 0.5452`\
`#> 7    0.0541 1.0000`\
`#> 8    0.0818 0.7539`

Each row records: \* **`Dataset`**: Experimental source or condition. \*
**`Method`**: Simulation algorithm evaluated. \* **`Category`**: One of
the 8 canonical evaluation categories (e.g.,
`(I) Distributional Properties`,
`(V) Biological Signal & Differential Fidelity`). \* **`Metric`**: Name
of the curated evaluation metric. \* **`Score_Raw`**: Unnormalized
empirical value in physical units (seconds, peak RAM in MiB, correlation
coefficient, Kolmogorov-Smirnov distance). \* **`Score`**: Standardized
direction-aware fidelity score normalized to $`[0, 1]`$, where $`1.0`$
represents top fidelity.

------------------------------------------------------------------------

## 3. Generating Simulator Leaderboards

Use
[`compute_method_leaderboard()`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md)
to compute composite rankings across all 56 curated evaluation metrics:

\
`# Compute standard overall composite leaderboard`\
`leaderboard`` ``<-`` `[`compute_method_leaderboard`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md)`(``bench_data``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``leaderboard``)`\
`#>   Overall_Rank     Method Average_Fidelity Fidelity_Score`\
`#> 4            1  scDesign3            85.3%         0.8529`\
`#> 6            2 scMultiSim            83.9%         0.8387`\
`#> 5            3    scMoMtF            83.1%         0.8310`\
`#> 3            4    scCross            80.2%         0.8023`\
`#> 2            5      MOSim            77.0%         0.7699`\
`#> 1            6    Matilda            76.6%         0.7658`

#### 3.1 Category-Specific Leaderboards

In addition to the global composite leaderboard, researchers often
evaluate simulator rankings within specific biological domains. For
example, to assess which method excels at capturing **Category I
(Distributional Properties)**:

\
`# Extract Category I evaluation rows`\
`cat1_df`` ``<-`` `[`subset`](https://rdrr.io/r/base/subset.html)`(``bench_data``$``benchmark_df``, ``Category`` ``==`` ``"(I) Distributional Properties"``)`\
\
`# Compute Category I leaderboard`\
`cat1_leaderboard`` ``<-`` `[`compute_method_leaderboard`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md)`(``cat1_df``)`\
[`print`](https://rdrr.io/r/base/print.html)`(``cat1_leaderboard``)`\
`#>   Overall_Rank     Method Average_Fidelity Fidelity_Score`\
`#> 4            1  scDesign3            94.3%         0.9430`\
`#> 6            2 scMultiSim            89.6%         0.8956`\
`#> 5            3    scMoMtF            86.4%         0.8636`\
`#> 3            4    scCross            84.2%         0.8422`\
`#> 1            5    Matilda            81.3%         0.8125`\
`#> 2            6      MOSim            77.0%         0.7701`

------------------------------------------------------------------------

## 4. Generating Publication-Quality Figures

All visualization functions in `scSimEval` accept the loaded `.rds`
object directly:

#### 4.1 Comparative Bubble Matrix

Generate the flagship multi-dimensional bubble matrix mapping every
simulator across all curated evaluation measures:

\
`# Multi-simulator comparative bubble matrix`\
`p_bubble`` ``<-`` `[`plot_benchmark_bubble_matrix`](https://kabilanbio.github.io/scSimEval/reference/plot_benchmark_bubble_matrix.md)`(`\
`  data       ``=`` ``bench_data``,`\
`  title      ``=`` ``"Cross-Dataset Simulation Fidelity Benchmarking"``,`\
`  theme_mode ``=`` ``"publication"`\
`)`\
[`print`](https://rdrr.io/r/base/print.html)`(``p_bubble``)`

#### 4.2 Category Summary Bars

Generate rank-ordered horizontal bar charts comparing simulator
performance across each category:

\
`# Category-wise summary horizontal bar chart`\
`p_summary`` ``<-`` `[`plot_evaluation_summary`](https://kabilanbio.github.io/scSimEval/reference/plot_evaluation_summary.md)`(`\
`  data  ``=`` ``bench_data``,`\
`  title ``=`` ``"Overall Simulation Fidelity by Category"`\
`)`\
[`print`](https://rdrr.io/r/base/print.html)`(``p_summary``)`

#### 4.3 Interactive Figure Zooming with `zoom_plot()`

When inspecting dense figures with many metrics or simulators, use
[`zoom_plot()`](https://kabilanbio.github.io/scSimEval/reference/zoom_plot.md)
to focus on specific coordinate windows:

\
`# Magnify top-left corner of the bubble matrix (first 10 metrics, top 3 simulators)`\
[`zoom_plot`](https://kabilanbio.github.io/scSimEval/reference/zoom_plot.md)`(``p_bubble``, x_range ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``1``, ``10``)``, y_range ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``4``, ``6``)``)`

------------------------------------------------------------------------

## 5. Merging Multiple Saved Benchmark Runs

When benchmark evaluations are executed across different computing nodes
or datasets at different times, use
[`combine_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
or
[`append_benchmark_results()`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)
to merge individual `.rds` objects into a single consolidated benchmark:

\
`# Load two independently computed benchmark runs`\
`run1`` ``<-`` `[`readRDS`](https://rdrr.io/r/base/readRDS.html)`(``"benchmark_dataset_A.rds"``)`\
`run2`` ``<-`` `[`readRDS`](https://rdrr.io/r/base/readRDS.html)`(``"benchmark_dataset_B.rds"``)`\
\
`# Merge into a single cross-dataset benchmark object`\
`merged_benchmark`` ``<-`` `[`combine_benchmark_results`](https://kabilanbio.github.io/scSimEval/reference/combine_benchmark_results.md)`(`\
`  ``run1``, ``run2``,`\
`  dataset_names ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``"Dataset_A"``, ``"Dataset_B"``)`\
`)`\
\
`# Compute unified cross-dataset leaderboard`\
[`compute_method_leaderboard`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md)`(``merged_benchmark``)`

------------------------------------------------------------------------

## 6. Instant Exploration in Shiny Studio

Any `.rds` benchmark file saved from
[`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md),
[`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md),
or
[`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md)
can also be loaded into the **scSimEval Shiny Studio** with zero
latency:

1.  Launch the app in R:

    \
    [`launch_scSimEval_app`](https://kabilanbio.github.io/scSimEval/reference/launch_scSimEval_app.md)`(``)`

2.  Navigate to **Data Hub** $`\rightarrow`$**Mode 4: Upload
    Pre-Computed Benchmark Table (.rds / .csv)**.

3.  Select your `.rds` file.

4.  The entire Shiny interface — including the 56-metric Bubble Matrix,
    Fidelity Leaderboards, and all 8 Visualization Panels — populates
    instantaneously without re-running any computations.

------------------------------------------------------------------------

## 7. Summary

Analyzing pre-saved `.rds` benchmark files allows you to: \* **Avoid
redundant computation**: Calculate metrics once on an HPC cluster and
re-analyze instantly on a local workstation. \* **Explore custom
biological priorities**: Re-weight categories dynamically using
[`compute_method_leaderboard()`](https://kabilanbio.github.io/scSimEval/reference/compute_method_leaderboard.md).
\* **Export publication figures**: Render high-resolution (600 DPI)
bubble matrices, heatmaps, and summary bars on demand. \* **Seamlessly
bridge CLI and GUI**: Move effortlessly between R console workflows and
the interactive Shiny Studio interface.

\
[`sessionInfo`](https://rdrr.io/r/utils/sessionInfo.html)`(``)`\
`#> R version 4.6.1 (2026-06-24 ucrt)`\
`#> Platform: x86_64-w64-mingw32/x64`\
`#> Running under: Windows 11 x64 (build 26300)`\
`#> `\
`#> Matrix products: default`\
`#>   LAPACK version 3.12.1`\
`#> `\
`#> locale:`\
`#> [1] LC_COLLATE=English_India.utf8  LC_CTYPE=English_India.utf8   `\
`#> [3] LC_MONETARY=English_India.utf8 LC_NUMERIC=C                  `\
`#> [5] LC_TIME=English_India.utf8    `\
`#> `\
`#> time zone: Asia/Calcutta`\
`#> tzcode source: internal`\
`#> `\
`#> attached base packages:`\
`#> [1] stats     graphics  grDevices utils     datasets  methods   base     `\
`#> `\
`#> other attached packages:`\
`#> [1] scSimEval_0.99.10`\
`#> `\
`#> loaded via a namespace (and not attached):`\
`#>   [1] RColorBrewer_1.1-3          jsonlite_2.0.0             `\
`#>   [3] magrittr_2.0.5              spatstat.utils_3.2-5       `\
`#>   [5] farver_2.1.2                rmarkdown_2.32             `\
`#>   [7] fs_2.1.0                    ragg_1.5.2                 `\
`#>   [9] vctrs_0.7.3                 ROCR_1.0-12                `\
`#>  [11] spatstat.explore_3.8-3      htmltools_0.5.9            `\
`#>  [13] S4Arrays_1.13.0             SparseArray_1.13.2         `\
`#>  [15] sass_0.4.10                 sctransform_0.4.3          `\
`#>  [17] parallelly_1.48.0           KernSmooth_2.23-26         `\
`#>  [19] bslib_0.12.0                htmlwidgets_1.6.4          `\
`#>  [21] desc_1.4.3                  ica_1.0-3                  `\
`#>  [23] plyr_1.8.9                  plotly_4.12.1              `\
`#>  [25] zoo_1.9-1                   cachem_1.1.0               `\
`#>  [27] igraph_2.3.3                mime_0.13                  `\
`#>  [29] lifecycle_1.0.5             pkgconfig_2.0.3            `\
`#>  [31] Matrix_1.7-6                R6_2.6.1                   `\
`#>  [33] fastmap_1.2.0               MatrixGenerics_1.25.0      `\
`#>  [35] fitdistrplus_1.2-6          future_1.76.0              `\
`#>  [37] shiny_1.14.0                digest_0.6.39              `\
`#>  [39] patchwork_1.3.2             S4Vectors_0.51.10          `\
`#>  [41] Seurat_5.5.1                tensor_1.5.1               `\
`#>  [43] RSpectra_0.16-2             irlba_2.3.7                `\
`#>  [45] textshaping_1.0.5           GenomicRanges_1.65.4       `\
`#>  [47] progressr_1.0.0             spatstat.sparse_3.2-0      `\
`#>  [49] polyclip_1.10-7             httr_1.4.9                 `\
`#>  [51] abind_1.4-8                 compiler_4.6.1             `\
`#>  [53] S7_0.2.2                    BiocParallel_1.47.0        `\
`#>  [55] fastDummies_1.7.6           MASS_7.3-66                `\
`#>  [57] DelayedArray_0.39.8         ggsci_5.2.0                `\
`#>  [59] tools_4.6.1                 lmtest_0.9-40              `\
`#>  [61] otel_0.2.0                  httpuv_1.6.17              `\
`#>  [63] future.apply_1.20.2         goftest_1.2-3              `\
`#>  [65] glue_1.8.1                  nlme_3.1-171               `\
`#>  [67] promises_1.5.0              grid_4.6.1                 `\
`#>  [69] Rtsne_0.17                  cluster_2.1.8.3            `\
`#>  [71] reshape2_1.4.5              generics_0.1.4             `\
`#>  [73] gtable_0.3.6                spatstat.data_3.1-9        `\
`#>  [75] tidyr_1.3.2                 data.table_1.18.6.1        `\
`#>  [77] sp_2.2-3                    XVector_0.52.0             `\
`#>  [79] BiocGenerics_0.59.12        spatstat.geom_3.8-3        `\
`#>  [81] RcppAnnoy_0.0.23            ggrepel_0.9.8              `\
`#>  [83] RANN_2.6.3                  pillar_1.11.1              `\
`#>  [85] stringr_1.6.0               spam_2.11-4                `\
`#>  [87] RcppHNSW_0.7.0              later_1.4.8                `\
`#>  [89] splines_4.6.1               dplyr_1.2.1                `\
`#>  [91] lattice_0.22-9              deldir_2.0-4               `\
`#>  [93] survival_3.8-12             tidyselect_1.2.1           `\
`#>  [95] SingleCellExperiment_1.35.2 CDI_1.11.0                 `\
`#>  [97] miniUI_0.1.2                pbapply_1.7-5              `\
`#>  [99] knitr_1.52                  gridExtra_2.3.1            `\
`#> [101] IRanges_2.47.5              Seqinfo_1.3.2              `\
`#> [103] SummarizedExperiment_1.43.0 scattermore_1.2            `\
`#> [105] stats4_4.6.1                xfun_0.60                  `\
`#> [107] Biobase_2.73.2              matrixStats_1.5.0          `\
`#> [109] stringi_1.8.9               yaml_2.3.12                `\
`#> [111] evaluate_1.0.5              codetools_0.2-20           `\
`#> [113] tibble_3.3.1                cli_3.6.6                  `\
`#> [115] uwot_0.2.5                  xtable_1.8-8               `\
`#> [117] reticulate_1.47.0           systemfonts_1.3.2          `\
`#> [119] jquerylib_0.1.4             dichromat_2.0-1            `\
`#> [121] Rcpp_1.1.2                  spatstat.random_3.5-2      `\
`#> [123] globals_0.19.1              png_0.1-9                  `\
`#> [125] spatstat.univar_3.2-0       parallel_4.6.1             `\
`#> [127] pkgdown_2.2.1               ggplot2_4.0.3              `\
`#> [129] dotCall64_1.2               listenv_1.0.0              `\
`#> [131] viridisLite_0.4.3           scales_1.4.0               `\
`#> [133] ggridges_0.5.7              SeuratObject_5.4.0         `\
`#> [135] purrr_1.2.2                 rlang_1.3.0                `\
`#> [137] cowplot_1.2.0`
