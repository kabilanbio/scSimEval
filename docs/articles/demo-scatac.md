# Demonstration 2: scATAC-seq Simulation Evaluation

## Introduction

Single-cell assay for transposase-accessible chromatin sequencing
(scATAC-seq) measures chromatin accessibility at single-cell resolution.
Unlike scRNA-seq expression counts, scATAC-seq data have unique
properties: \* **High Sparsity:** Over 90% to 98% of values in
peak-by-cell matrices are zeroes because each locus has at most 2 copies
per cell. \* **Near-Binary Nature:** Values primarily represent open
(accessible) vs. closed chromatin states. \* **Peak Co-Accessibility:**
Genomic loci regulated together exhibit correlated accessibility across
cells.

This tutorial shows how to evaluate synthetic scATAC-seq data using
**`scSimEval`**.

------------------------------------------------------------------------

## 1. Loading the Packaged scATAC-seq Dataset

`scSimEval` includes an empirical reference and simulated scATAC-seq
peak accessibility dataset (`example_scatac`):

\
[`library`](https://rdrr.io/r/base/library.html)`(`[`scSimEval`](https://kabilanbio.github.io/scSimEval/)`)`\
\
`# Load packaged scATAC-seq example data`\
[`data`](https://rdrr.io/r/utils/data.html)`(``example_scatac``)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Reference scATAC-seq dimensions:"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``example_scatac``$``ref``)``, ``"\n"``)`\
`#> Reference scATAC-seq dimensions: 60 80`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Simulated scATAC-seq dimensions:"``, `[`dim`](https://rdrr.io/r/base/dim.html)`(``example_scatac``$``sim``)``, ``"\n"``)`\
`#> Simulated scATAC-seq dimensions: 60 80`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Cell types:"``, `[`levels`](https://rdrr.io/r/base/levels.html)`(``example_scatac``$``cell_types``)``, ``"\n"``)`\
`#> Cell types: TypeA TypeB`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Batches:"``, `[`levels`](https://rdrr.io/r/base/levels.html)`(``example_scatac``$``batch_info``)``, ``"\n"``)`\
`#> Batches: Batch1 Batch2`

The matrices contain integer insertion counts across genomic peaks
($`60\text{ peaks} \times 80\text{ cells}`$).

------------------------------------------------------------------------

## 2. Epigenomic Distribution & Summary Properties

We evaluate whether simulated peak insertion frequencies and cell
library sizes match real distributions using
[`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md):

\
`atac_res`` ``<-`` `[`evaluate_simulation_accuracy`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)`(`\
`  ref_data          ``=`` ``example_scatac``$``ref``,`\
`  sim_data          ``=`` ``example_scatac``$``sim``,`\
`  compute_bivariate ``=`` ``FALSE``,`\
`  verbose           ``=`` ``FALSE`\
`)`\
`#> calcNormFactors has been renamed to normLibSizes`\
`#> calcNormFactors has been renamed to normLibSizes`\
\
`# Inspect peak accessibility and library size statistical distances`\
[`head`](https://rdrr.io/r/utils/head.html)`(``atac_res``$``metrics_summary_table``, ``8``)`\
`#>                    Category     Property                      Metric`\
`#> 1 Distributional Properties library_size                 Wasserstein`\
`#> 2 Distributional Properties library_size                          KS`\
`#> 3 Distributional Properties library_size                         MAD`\
`#> 4 Distributional Properties library_size               Bhattacharyya`\
`#> 5 Distributional Properties library_size              Runs_Statistic`\
`#> 6 Distributional Properties library_size                 NN_Mismatch`\
`#> 7 Distributional Properties library_size Between_Dataset_Silh_Global`\
`#> 8 Distributional Properties library_size                         MAE`\
`#>           Value`\
`#> 1  2.3375000000`\
`#> 2  0.1750000000`\
`#> 3  2.0000000000`\
`#> 4  0.0001055641`\
`#> 5 -6.0273112590`\
`#> 6  0.1937500000`\
`#> 7  0.0218867260`\
`#> 8  2.3375000000`

------------------------------------------------------------------------

## 3. Epigenomic Cluster & Cell-Type Concordance

In scATAC-seq, identifying cell types requires clustering cells based on
chromatin accessibility profiles. We test whether the simulator
preserves genuine biological separation using
[`evaluate_clustering_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md):

\
`atac_clust`` ``<-`` `[`evaluate_clustering_metrics`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md)`(`\
`  data       ``=`` ``example_scatac``$``sim``,`\
`  cell_types ``=`` ``example_scatac``$``cell_types``,`\
`  ref_data   ``=`` ``example_scatac``$``ref`\
`)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Average Silhouette Width (Simulated):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``silhouette``, ``4``)``, ``"\n"``)`\
`#> Average Silhouette Width (Simulated): 0.0415`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Neighborhood Purity:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``neighborhood_purity``, ``4``)``, ``"\n"``)`\
`#> Neighborhood Purity: 0.7275`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Adjusted Rand Index (ARI):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``ARI``, ``4``)``, ``"\n"``)`\
`#> Adjusted Rand Index (ARI): 0.5955`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Normalized Mutual Information (NMI):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_clust``$``NMI``, ``4``)``, ``"\n"``)`\
`#> Normalized Mutual Information (NMI): 0.4934`

------------------------------------------------------------------------

## 4. Epigenomic Batch Integration & Confounder Handling

When simulating multi-sample experiments, batch variations can confound
biological signals. We evaluate batch mixing metrics using
[`evaluate_batch_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md):

\
`atac_batch`` ``<-`` `[`evaluate_batch_metrics`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md)`(`\
`  data       ``=`` ``example_scatac``$``sim``,`\
`  batch_info ``=`` ``example_scatac``$``batch_info``,`\
`  cell_types ``=`` ``example_scatac``$``cell_types``,`\
`  verbose    ``=`` ``FALSE`\
`)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Batch Silhouette:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_batch``$``batch_silhouette``, ``4``)``, ``"\n"``)`\
`#> Batch Silhouette: 9e-04`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Principal Component Regression (PCR R2):"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_batch``$``pcr_r2``, ``4``)``, ``"\n"``)`\
`#> Principal Component Regression (PCR R2): 0.0138`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Cross-Batch Transfer F1:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``atac_batch``$``cross_batch_F1``, ``4``)``, ``"\n"``)`\
`#> Cross-Batch Transfer F1: 0.7347`

------------------------------------------------------------------------

## 5. Chromatin Accessibility-Sparsity Curve Dynamics (Category 2)

In real scATAC-seq experiments, peak detection proportion follows a
quadratic relationship with mean peak accessibility. We evaluate how
well the simulated peak matrix reproduces this accessibility-sparsity
curve using
[`evaluate_accessibility_sparsity_curve()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_accessibility_sparsity_curve.md):

\
`# Evaluate quadratic accessibility-sparsity curve fidelity`\
`sparsity_res`` ``<-`` `[`evaluate_accessibility_sparsity_curve`](https://kabilanbio.github.io/scSimEval/reference/evaluate_accessibility_sparsity_curve.md)`(`\
`  ref_data ``=`` ``example_scatac``$``ref``,`\
`  sim_data ``=`` ``example_scatac``$``sim`\
`)`\
\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Accessibility-Sparsity Curve RMSE:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``sparsity_res``$``curve_rmse``, ``4``)``, ``"\n"``)`\
`#> Accessibility-Sparsity Curve RMSE: 0.0123`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Quadratic Fit R2 Discrepancy:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``sparsity_res``$``r2_discrepancy``, ``4``)``, ``"\n"``)`\
`#> Quadratic Fit R2 Discrepancy: 9e-04`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Reference Curve R2:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``sparsity_res``$``ref_r_squared``, ``4``)``, ``"\n"``)`\
`#> Reference Curve R2: 0.9693`\
[`cat`](https://rdrr.io/r/base/cat.html)`(``"Simulated Curve R2:"``, `[`round`](https://rdrr.io/r/base/Round.html)`(``sparsity_res``$``sim_r_squared``, ``4``)``, ``"\n"``)`\
`#> Simulated Curve R2: 0.9684`

------------------------------------------------------------------------

## 6. Master Unimodal scATAC-seq Pipeline

Execute the comprehensive evaluation workflow in a single coordinated
command:

\
`master_atac`` ``<-`` `[`evaluate_simulation_accuracy`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md)`(`\
`  ref_data          ``=`` ``example_scatac``$``ref``,`\
`  sim_data          ``=`` ``example_scatac``$``sim``,`\
`  elapsed_time      ``=`` ``18.2``,`\
`  memory_mb         ``=`` ``512.4``,`\
`  compute_bivariate ``=`` ``FALSE``,`\
`  verbose           ``=`` ``FALSE`\
`)`\
`#> calcNormFactors has been renamed to normLibSizes`\
`#> calcNormFactors has been renamed to normLibSizes`\
\
`# Inspect tidy summary table`\
[`head`](https://rdrr.io/r/utils/head.html)`(``master_atac``$``metrics_summary_table``)`\
`#>                    Category     Property         Metric         Value`\
`#> 1 Distributional Properties library_size    Wasserstein  2.3375000000`\
`#> 2 Distributional Properties library_size             KS  0.1750000000`\
`#> 3 Distributional Properties library_size            MAD  2.0000000000`\
`#> 4 Distributional Properties library_size  Bhattacharyya  0.0001055641`\
`#> 5 Distributional Properties library_size Runs_Statistic -6.0273112590`\
`#> 6 Distributional Properties library_size    NN_Mismatch  0.1937500000`

------------------------------------------------------------------------

## 7. Conclusion

`scSimEval` effectively benchmarks scATAC-seq simulations, validating
peak insertion distributions, extreme sparsity, cell type clustering,
batch mixing, and chromatin peak co-accessibility across the curated
evaluation measures suite.

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
`#>   [1] RcppAnnoy_0.0.23            CellMixS_1.29.0            `\
`#>   [3] splines_4.6.1               later_1.4.8                `\
`#>   [5] filelock_1.0.3              tibble_3.3.1               `\
`#>   [7] polyclip_1.10-7             fastDummies_1.7.6          `\
`#>   [9] lifecycle_1.0.5             httr2_1.3.0                `\
`#>  [11] Rdpack_2.6.6                edgeR_4.99.6               `\
`#>  [13] globals_0.19.1              lattice_0.22-9             `\
`#>  [15] MASS_7.3-66                 magrittr_2.0.5             `\
`#>  [17] limma_3.99.0                plotly_4.12.1              `\
`#>  [19] sass_0.4.10                 rmarkdown_2.32             `\
`#>  [21] jquerylib_0.1.4             yaml_2.3.12                `\
`#>  [23] httpuv_1.6.17               otel_0.2.0                 `\
`#>  [25] Seurat_5.5.1                sctransform_0.4.3          `\
`#>  [27] spam_2.11-4                 sp_2.2-3                   `\
`#>  [29] spatstat.sparse_3.2-0       reticulate_1.47.0          `\
`#>  [31] minqa_1.2.8                 cowplot_1.2.0              `\
`#>  [33] pbapply_1.7-5               DBI_1.3.0                  `\
`#>  [35] RColorBrewer_1.1-3          ade4_1.7-24                `\
`#>  [37] abind_1.4-8                 Rtsne_0.17                 `\
`#>  [39] GenomicRanges_1.65.4        purrr_1.2.2                `\
`#>  [41] BiocGenerics_0.59.12        IRanges_2.47.5             `\
`#>  [43] S4Vectors_0.51.10           ggrepel_0.9.8              `\
`#>  [45] irlba_2.3.7                 listenv_1.0.0              `\
`#>  [47] spatstat.utils_3.2-5        arm_1.15-3                 `\
`#>  [49] goftest_1.2-3               RSpectra_0.16-2            `\
`#>  [51] spatstat.random_3.5-2       fitdistrplus_1.2-6         `\
`#>  [53] parallelly_1.48.0           provenance_4.4             `\
`#>  [55] pkgdown_2.2.1               codetools_0.2-20           `\
`#>  [57] DelayedArray_0.39.8         scuttle_1.23.2             `\
`#>  [59] SuppDists_1.1-9.9           tidyselect_1.2.1           `\
`#>  [61] farver_2.1.2                viridis_0.6.5              `\
`#>  [63] ScaledMatrix_1.21.0         lme4_2.0-6                 `\
`#>  [65] matrixStats_1.5.0           stats4_4.6.1               `\
`#>  [67] BiocFileCache_3.3.0         spatstat.explore_3.8-3     `\
`#>  [69] Seqinfo_1.3.2               jsonlite_2.0.0             `\
`#>  [71] BiocNeighbors_2.7.3         e1071_1.7-17               `\
`#>  [73] progressr_1.0.0             scater_1.41.2              `\
`#>  [75] ggridges_0.5.7              survival_3.8-12            `\
`#>  [77] systemfonts_1.3.2           tools_4.6.1                `\
`#>  [79] ragg_1.5.2                  ica_1.0-3                  `\
`#>  [81] Rcpp_1.1.2                  glue_1.8.1                 `\
`#>  [83] gridExtra_2.3.1             SparseArray_1.13.2         `\
`#>  [85] xfun_0.60                   kSamples_1.2-12            `\
`#>  [87] MatrixGenerics_1.25.0       dplyr_1.2.1                `\
`#>  [89] withr_3.0.3                 fastmap_1.2.0              `\
`#>  [91] bluster_1.23.1              boot_1.3-32                `\
`#>  [93] rsvd_1.0.5                  digest_0.6.39              `\
`#>  [95] R6_2.6.1                    mime_0.13                  `\
`#>  [97] textshaping_1.0.5           scattermore_1.2            `\
`#>  [99] tensor_1.5.1                dichromat_2.0-1            `\
`#> [101] spatstat.data_3.1-9         RSQLite_3.53.3             `\
`#> [103] RhpcBLASctl_0.23-42         tidyr_1.3.2                `\
`#> [105] generics_0.1.4              ggsci_5.2.0                `\
`#> [107] data.table_1.18.6.1         class_7.3-24               `\
`#> [109] httr_1.4.9                  htmlwidgets_1.6.4          `\
`#> [111] S4Arrays_1.13.0             uwot_0.2.5                 `\
`#> [113] pkgconfig_2.0.3             gtable_0.3.6               `\
`#> [115] blob_1.3.0                  lmtest_0.9-40              `\
`#> [117] S7_0.2.2                    MLmetrics_1.1.3            `\
`#> [119] SingleCellExperiment_1.35.2 XVector_0.52.0             `\
`#> [121] htmltools_0.5.9             clusterSim_0.51-6          `\
`#> [123] dotCall64_1.2               SeuratObject_5.4.0         `\
`#> [125] scales_1.4.0                Biobase_2.73.2             `\
`#> [127] png_0.1-9                   reformulas_0.4.4           `\
`#> [129] spatstat.univar_3.2-0       knitr_1.52                 `\
`#> [131] reshape2_1.4.5              nloptr_2.2.1               `\
`#> [133] coda_0.19-4.1               curl_8.0.0                 `\
`#> [135] nlme_3.1-171                proxy_0.4-29               `\
`#> [137] cachem_1.1.0                zoo_1.9-1                  `\
`#> [139] stringr_1.6.0               KernSmooth_2.23-26         `\
`#> [141] vipor_0.4.7                 parallel_4.6.1             `\
`#> [143] miniUI_0.1.2                desc_1.4.3                 `\
`#> [145] CDI_1.11.0                  pillar_1.11.1              `\
`#> [147] grid_4.6.1                  vctrs_0.7.3                `\
`#> [149] RANN_2.6.3                  promises_1.5.0             `\
`#> [151] BiocSingular_1.29.1         beachmat_2.29.3            `\
`#> [153] dbplyr_2.6.0                xtable_1.8-8               `\
`#> [155] cluster_2.1.8.3             beeswarm_0.4.0             `\
`#> [157] evaluate_1.0.5              cli_3.6.6                  `\
`#> [159] locfit_1.5-9.12             compiler_4.6.1             `\
`#> [161] rlang_1.3.0                 future.apply_1.20.2        `\
`#> [163] mclust_6.1.3                ggbeeswarm_0.7.3           `\
`#> [165] plyr_1.8.9                  fs_2.1.0                   `\
`#> [167] stringi_1.8.9               viridisLite_0.4.3          `\
`#> [169] deldir_2.0-4                BiocParallel_1.47.0        `\
`#> [171] philentropy_0.10.0          IsoplotR_7.0               `\
`#> [173] waddR_1.27.0                spatstat.geom_3.8-3        `\
`#> [175] Matrix_1.7-6                RcppHNSW_0.7.0             `\
`#> [177] patchwork_1.3.2             bit64_4.8.6                `\
`#> [179] future_1.76.0               eva_0.2.7                  `\
`#> [181] ggplot2_4.0.3               statmod_1.5.2              `\
`#> [183] shiny_1.14.0                SummarizedExperiment_1.43.0`\
`#> [185] rbibutils_2.4.1             ROCR_1.0-12                `\
`#> [187] igraph_2.3.3                memoise_2.0.1              `\
`#> [189] RcppParallel_6.2.1          bslib_0.12.0               `\
`#> [191] bit_4.6.0                   randtests_1.0.2`
