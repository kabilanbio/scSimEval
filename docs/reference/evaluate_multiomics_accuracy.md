# Master Multiomics Simulation Benchmarking Pipeline

Orchestrates comprehensive benchmarking across all 8 evaluation
categories using strictly the user's provided inputs:

- `ref_multi`: List of reference omics objects/matrices (e.g.
  `list(RNA = ..., ATAC = ...)`).

- `sim_multi`: List of simulated omics objects/matrices (e.g.
  `list(RNA = ..., ATAC = ...)`).

- `cell_types`: Vector/factor of cell-type identities.

- `batch_info`: Vector/factor of batch or biological replicate
  annotations.

- Computational resource usage: `elapsed_time`, `memory_mb`.

- Trajectory features inferred directly from scRNA-seq counts.

## Usage

``` r
evaluate_multiomics_accuracy(
  ref_multi,
  sim_multi,
  cell_types = NULL,
  batch_info = NULL,
  pairing = c("paired", "unpaired"),
  is_paired = NULL,
  memory_mb = NULL,
  elapsed_time = NULL,
  feature_pairs = NULL,
  compute_bivariate = FALSE,
  threads = 1,
  verbose = TRUE,
  time_unit = c("seconds", "minutes"),
  memory_unit = c("MiB", "MB"),
  sim_cores = NULL,
  elapsed_time_rna = NULL,
  elapsed_time_atac = NULL,
  memory_mb_rna = NULL,
  memory_mb_atac = NULL,
  time_unit_rna = NULL,
  time_unit_atac = NULL,
  memory_unit_rna = NULL,
  memory_unit_atac = NULL,
  sim_cores_rna = NULL,
  sim_cores_atac = NULL,
  trajectory_num_trees = 10000L
)
```

## Arguments

- ref_multi:

  Named list of reference matrices for Modality 1 and 2.

- sim_multi:

  Named list of simulated matrices for Modality 1 and 2.

- cell_types:

  Optional factor or vector of cell type labels for cells.

- batch_info:

  Optional factor or vector of batch labels for cells.

- pairing:

  Character string specifying whether the multiomics dataset is
  `"paired"` (simultaneous co-assay from the same individual cells,
  e.g., 10x Chromium Multiome, SHARE-seq, SNARE-seq) or `"unpaired"`
  (independent profiling of modalities from separate cells of the same
  biological tissue). Default is `"paired"`. In `"unpaired"` mode,
  cell-level pairing metrics (FOSCTTM, Match@1, cross-modal generation,
  and direct peak-to-gene coupling) are skipped, while population-level
  cross-modal metrics (cross-modal label transfer, co-expression module
  fidelity, peak co-accessibility, and accessibility profile
  concordance) are computed alongside all unimodal metrics.

- is_paired:

  Optional logical alias for `pairing`. If provided, `TRUE` corresponds
  to `"paired"` and `FALSE` to `"unpaired"`.

- memory_mb:

  Optional numeric value of peak memory allocation (in `memory_unit`).

- elapsed_time:

  Optional numeric value of wall-clock elapsed time (in `time_unit`).

- feature_pairs:

  Optional 2-column data.frame of linked feature pairs.

- compute_bivariate:

  Logical, whether to compute 2D bivariate tests. Default FALSE for
  speed.

- threads:

  Number of CPU threads. Default is 1.

- verbose:

  Logical, whether to print execution progress. Default is TRUE.

- time_unit:

  Character string specifying the unit of `elapsed_time`: `"seconds"`
  (default, or `"s"`) or `"minutes"` (or `"min"`).

- memory_unit:

  Character string specifying the unit of `memory_mb`: `"MiB"` (default,
  binary mebibytes) or `"MB"` (decimal megabytes, converted to MiB).

- sim_cores:

  Optional positive integer specifying the number of CPU cores used by
  the simulator to generate `sim_multi`. Default is NULL.

- elapsed_time_rna, elapsed_time_atac:

  Optional separate wall-clock runtimes for RNA and ATAC layers when
  modalities are simulated independently (used if combined
  `elapsed_time` is NULL).

- memory_mb_rna, memory_mb_atac:

  Optional separate peak RAM allocations for RNA and ATAC layers when
  modalities are simulated independently (used if combined `memory_mb`
  is NULL).

- time_unit_rna, time_unit_atac:

  Optional time units (`"seconds"` or `"minutes"`) for
  `elapsed_time_rna` and `elapsed_time_atac`; defaults to `time_unit`.

- memory_unit_rna, memory_unit_atac:

  Optional memory units (`"MiB"` or `"MB"`) for `memory_mb_rna` and
  `memory_mb_atac`; defaults to `memory_unit`.

- sim_cores_rna, sim_cores_atac:

  Optional positive integers specifying simulation CPU cores for
  separate RNA and ATAC runs (used if combined `sim_cores` is NULL).

- trajectory_num_trees:

  Number of Random Forest trees for Category VII `wcor_features`
  (default `10000L`, per Saelens et al. 2019).

## Value

A list containing detailed results and a unified tidy master summary
table.

## Details

![scSimEval Unified Benchmarking Workflow
Diagram](figures/workflow_diagram.png)

**Multiomics Metric Measurement Capacity Across the 56 Evaluation
Measures:**

- **Paired or Unpaired Count Data Alone:** Evaluates unimodal Category I
  (10 univariate + 4 bivariate metrics when `compute_bivariate = TRUE`),
  Category II (5 sparsity/zero-inflation metrics across RNA & ATAC),
  Category III manifold metrics (MMD, Generative Precision/Recall), and
  Category VII trajectory metrics inferred from scRNA-seq via the
  HVG-filtered dynverse workflow (`him`, `edge_flip`,
  `geodesic_correlation`, `F1_branches`, `wcor_features`).

- **+ Cell Types:** Unlocks Category III supervised clustering
  concordance (Silhouette ASW, CDI, Neighborhood Purity, ARI, NMI) and
  Category VI (Multiomics Integration & Alignment: Batch ASW, Cell-Type
  ASW, Graph Connectivity, Integration ARI, Integration NMI).

- **+ Batch / Biological Replicate Annotations:** Unlocks Category IV
  batch mixing across both layers and Category V (12 differential
  fidelity metrics: null P-value uniformity, accuracy, precision,
  recall, F1-score, ROC-AUC across both scRNA-seq DE and scATAC-seq DA).

- **+ Runtime & Memory:** Logs Category VIII (`core_seconds`,
  `memory_per_core_mib`) to complete the full 56-measure evaluation
  suite across all 8 canonical categories.

## Examples

``` r
data(example_multiomics, package = "scSimEval")
evaluate_multiomics_accuracy(example_multiomics$ref_multi, example_multiomics$sim_multi)
#> === [Step 1] Evaluating rna Layer Accuracy ===
#> [1/5] Extracting cell-level properties...
#> [2/5] Extracting feature-level properties...
#> [3/5] Computing univariate accuracy metrics...
#> [5/5] Computing zero-probability and manifold distances...
#> Evaluating batch integration & confounder mixing...
#> Evaluating scRNA-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> Unimodal accuracy evaluation complete.
#> === [Step 2] Evaluating atac Layer Accuracy ===
#> [1/5] Extracting cell-level properties...
#> [2/5] Extracting feature-level properties...
#> [3/5] Computing univariate accuracy metrics...
#> [5/5] Computing zero-probability and manifold distances...
#> Evaluating batch integration & confounder mixing...
#> Evaluating scATAC-seq biological signal fidelity (Libra Pseudobulk edgeR-QLF)...
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> [6/7] Category VII (Trajectory & Lineage Topology) requires scRNA-seq data: skipping for scATAC-seq.
#> Unimodal accuracy evaluation complete.
#> === [Step 3] Evaluating Cellular Structure & Clustering Preservation (rna & atac) ===
#> === [Step 4] Evaluating Batch Integration & Technical Confounder Mixing ===
#> === [Step 5] Evaluating rna Biological Signal Fidelity (Libra Pseudobulk edgeR-QLF) ===
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> === [Step 5] Evaluating atac Biological Signal Fidelity (Libra Pseudobulk edgeR-QLF) ===
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
#> === [Step 6] Evaluating Trajectory & Lineage Topology Inferred from rna (HVG + dynverse) ===
#> === [Step 7] Evaluating Multiomics Integration & Alignment (Category VI) ===
#> Using precomputed dimensional reduction from sim_multi for Category VI.
#> Master multiomics benchmarking complete.
#> $benchmark_summary_table
#>                                         Category
#> 1                      Distributional Properties
#> 2                      Distributional Properties
#> 3                      Distributional Properties
#> 4                      Distributional Properties
#> 5                      Distributional Properties
#> 6                      Distributional Properties
#> 7                      Distributional Properties
#> 8                      Distributional Properties
#> 9                      Distributional Properties
#> 10                     Distributional Properties
#> 11                     Distributional Properties
#> 12                     Distributional Properties
#> 13                     Distributional Properties
#> 14                     Distributional Properties
#> 15                     Distributional Properties
#> 16                     Distributional Properties
#> 17                     Distributional Properties
#> 18                     Distributional Properties
#> 19                     Distributional Properties
#> 20                     Distributional Properties
#> 21                     Distributional Properties
#> 22                     Distributional Properties
#> 23                     Distributional Properties
#> 24                     Distributional Properties
#> 25                     Distributional Properties
#> 26                     Distributional Properties
#> 27                     Distributional Properties
#> 28                     Distributional Properties
#> 29                     Distributional Properties
#> 30                     Distributional Properties
#> 31                     Distributional Properties
#> 32                     Distributional Properties
#> 33                     Distributional Properties
#> 34                     Distributional Properties
#> 35                     Distributional Properties
#> 36                     Distributional Properties
#> 37                     Distributional Properties
#> 38                     Distributional Properties
#> 39                     Distributional Properties
#> 40                     Distributional Properties
#> 41                     Distributional Properties
#> 42                     Distributional Properties
#> 43                     Distributional Properties
#> 44                     Distributional Properties
#> 45                     Distributional Properties
#> 46                     Distributional Properties
#> 47                     Distributional Properties
#> 48                     Distributional Properties
#> 49                     Distributional Properties
#> 50                     Distributional Properties
#> 51                     Distributional Properties
#> 52                     Distributional Properties
#> 53                     Distributional Properties
#> 54                     Distributional Properties
#> 55                     Distributional Properties
#> 56                     Distributional Properties
#> 57                     Distributional Properties
#> 58                     Distributional Properties
#> 59                     Distributional Properties
#> 60                     Distributional Properties
#> 61                     Distributional Properties
#> 62                     Distributional Properties
#> 63                     Distributional Properties
#> 64                     Distributional Properties
#> 65                     Distributional Properties
#> 66                     Distributional Properties
#> 67                     Distributional Properties
#> 68                     Distributional Properties
#> 69                     Distributional Properties
#> 70                     Distributional Properties
#> 71                     Distributional Properties
#> 72                     Distributional Properties
#> 73                     Distributional Properties
#> 74                     Distributional Properties
#> 75                     Distributional Properties
#> 76                     Distributional Properties
#> 77                     Distributional Properties
#> 78                     Distributional Properties
#> 79                     Distributional Properties
#> 80                     Distributional Properties
#> 81                     Distributional Properties
#> 82                     Distributional Properties
#> 83                     Distributional Properties
#> 84                     Distributional Properties
#> 85                     Distributional Properties
#> 86                     Distributional Properties
#> 87                     Distributional Properties
#> 88                     Distributional Properties
#> 89                     Distributional Properties
#> 90                     Distributional Properties
#> 91                     Distributional Properties
#> 92                     Distributional Properties
#> 93                     Distributional Properties
#> 94                     Distributional Properties
#> 95                     Distributional Properties
#> 96                     Distributional Properties
#> 97                     Distributional Properties
#> 98                     Distributional Properties
#> 99                     Distributional Properties
#> 100                    Distributional Properties
#> 101                    Distributional Properties
#> 102                    Distributional Properties
#> 103                    Distributional Properties
#> 104                    Distributional Properties
#> 105                    Distributional Properties
#> 106                    Distributional Properties
#> 107                    Distributional Properties
#> 108                    Distributional Properties
#> 109                    Distributional Properties
#> 110                    Distributional Properties
#> 111                    Distributional Properties
#> 112                    Distributional Properties
#> 113                    Distributional Properties
#> 114                    Distributional Properties
#> 115                    Distributional Properties
#> 116                    Distributional Properties
#> 117                    Distributional Properties
#> 118                    Distributional Properties
#> 119                    Distributional Properties
#> 120                    Distributional Properties
#> 121                    Distributional Properties
#> 122                    Distributional Properties
#> 123                    Distributional Properties
#> 124                    Distributional Properties
#> 125                    Distributional Properties
#> 126                    Distributional Properties
#> 127                    Distributional Properties
#> 128                    Distributional Properties
#> 129                    Distributional Properties
#> 130                    Distributional Properties
#> 131                    Distributional Properties
#> 132                    Distributional Properties
#> 133                    Distributional Properties
#> 134                    Distributional Properties
#> 135                    Distributional Properties
#> 136                    Distributional Properties
#> 137                    Distributional Properties
#> 138                    Distributional Properties
#> 139                    Distributional Properties
#> 140                    Distributional Properties
#> 141           Sparsity & Zero-Inflation Dynamics
#> 142           Sparsity & Zero-Inflation Dynamics
#> 143           Sparsity & Zero-Inflation Dynamics
#> 144           Sparsity & Zero-Inflation Dynamics
#> 145   Cellular Manifold & Clustering Concordance
#> 146   Cellular Manifold & Clustering Concordance
#> 147   Cellular Manifold & Clustering Concordance
#> 151  Technical Variation & Batch Mixing Dynamics
#> 152  Technical Variation & Batch Mixing Dynamics
#> 153  Technical Variation & Batch Mixing Dynamics
#> 154  Technical Variation & Batch Mixing Dynamics
#> 155  Technical Variation & Batch Mixing Dynamics
#> 156    Biological Signal & Differential Fidelity
#> 157    Biological Signal & Differential Fidelity
#> 158    Biological Signal & Differential Fidelity
#> 159    Biological Signal & Differential Fidelity
#> 160    Biological Signal & Differential Fidelity
#> 161    Biological Signal & Differential Fidelity
#> 162                    Distributional Properties
#> 210                    Distributional Properties
#> 310                    Distributional Properties
#> 410                    Distributional Properties
#> 510                    Distributional Properties
#> 610                    Distributional Properties
#> 710                    Distributional Properties
#> 810                    Distributional Properties
#> 910                    Distributional Properties
#> 1010                   Distributional Properties
#> 1110                   Distributional Properties
#> 1210                   Distributional Properties
#> 1310                   Distributional Properties
#> 1410                   Distributional Properties
#> 1510                   Distributional Properties
#> 163                    Distributional Properties
#> 171                    Distributional Properties
#> 181                    Distributional Properties
#> 191                    Distributional Properties
#> 201                    Distributional Properties
#> 211                    Distributional Properties
#> 221                    Distributional Properties
#> 231                    Distributional Properties
#> 241                    Distributional Properties
#> 251                    Distributional Properties
#> 261                    Distributional Properties
#> 271                    Distributional Properties
#> 281                    Distributional Properties
#> 291                    Distributional Properties
#> 301                    Distributional Properties
#> 311                    Distributional Properties
#> 321                    Distributional Properties
#> 331                    Distributional Properties
#> 341                    Distributional Properties
#> 351                    Distributional Properties
#> 361                    Distributional Properties
#> 371                    Distributional Properties
#> 381                    Distributional Properties
#> 391                    Distributional Properties
#> 401                    Distributional Properties
#> 411                    Distributional Properties
#> 421                    Distributional Properties
#> 431                    Distributional Properties
#> 441                    Distributional Properties
#> 451                    Distributional Properties
#> 461                    Distributional Properties
#> 471                    Distributional Properties
#> 481                    Distributional Properties
#> 491                    Distributional Properties
#> 501                    Distributional Properties
#> 511                    Distributional Properties
#> 521                    Distributional Properties
#> 531                    Distributional Properties
#> 541                    Distributional Properties
#> 551                    Distributional Properties
#> 561                    Distributional Properties
#> 571                    Distributional Properties
#> 581                    Distributional Properties
#> 591                    Distributional Properties
#> 601                    Distributional Properties
#> 611                    Distributional Properties
#> 621                    Distributional Properties
#> 631                    Distributional Properties
#> 641                    Distributional Properties
#> 651                    Distributional Properties
#> 661                    Distributional Properties
#> 671                    Distributional Properties
#> 681                    Distributional Properties
#> 691                    Distributional Properties
#> 701                    Distributional Properties
#> 711                    Distributional Properties
#> 721                    Distributional Properties
#> 731                    Distributional Properties
#> 741                    Distributional Properties
#> 751                    Distributional Properties
#> 761                    Distributional Properties
#> 771                    Distributional Properties
#> 781                    Distributional Properties
#> 791                    Distributional Properties
#> 801                    Distributional Properties
#> 811                    Distributional Properties
#> 821                    Distributional Properties
#> 831                    Distributional Properties
#> 841                    Distributional Properties
#> 851                    Distributional Properties
#> 861                    Distributional Properties
#> 871                    Distributional Properties
#> 881                    Distributional Properties
#> 891                    Distributional Properties
#> 901                    Distributional Properties
#> 911                    Distributional Properties
#> 921                    Distributional Properties
#> 931                    Distributional Properties
#> 941                    Distributional Properties
#> 951                    Distributional Properties
#> 961                    Distributional Properties
#> 971                    Distributional Properties
#> 981                    Distributional Properties
#> 991                    Distributional Properties
#> 1001                   Distributional Properties
#> 1011                   Distributional Properties
#> 1021                   Distributional Properties
#> 1031                   Distributional Properties
#> 1041                   Distributional Properties
#> 1051                   Distributional Properties
#> 1061                   Distributional Properties
#> 1071                   Distributional Properties
#> 1081                   Distributional Properties
#> 1091                   Distributional Properties
#> 1101                   Distributional Properties
#> 1111                   Distributional Properties
#> 1121                   Distributional Properties
#> 1131                   Distributional Properties
#> 1141                   Distributional Properties
#> 1151                   Distributional Properties
#> 1161                   Distributional Properties
#> 1171                   Distributional Properties
#> 1181                   Distributional Properties
#> 1191                   Distributional Properties
#> 1201                   Distributional Properties
#> 1211                   Distributional Properties
#> 1221                   Distributional Properties
#> 1231                   Distributional Properties
#> 1241                   Distributional Properties
#> 1251                   Distributional Properties
#> 1261                   Distributional Properties
#> 1271                   Distributional Properties
#> 1281                   Distributional Properties
#> 1291                   Distributional Properties
#> 1301                   Distributional Properties
#> 1311          Sparsity & Zero-Inflation Dynamics
#> 1321          Sparsity & Zero-Inflation Dynamics
#> 1331          Sparsity & Zero-Inflation Dynamics
#> 1341  Cellular Manifold & Clustering Concordance
#> 1351  Cellular Manifold & Clustering Concordance
#> 1361  Cellular Manifold & Clustering Concordance
#> 1401 Technical Variation & Batch Mixing Dynamics
#> 1411 Technical Variation & Batch Mixing Dynamics
#> 1421 Technical Variation & Batch Mixing Dynamics
#> 1431 Technical Variation & Batch Mixing Dynamics
#> 1441 Technical Variation & Batch Mixing Dynamics
#> 1451   Biological Signal & Differential Fidelity
#> 1461   Biological Signal & Differential Fidelity
#> 1471   Biological Signal & Differential Fidelity
#> 148    Biological Signal & Differential Fidelity
#> 149    Biological Signal & Differential Fidelity
#> 150    Biological Signal & Differential Fidelity
#> 164   Cellular Manifold & Clustering Concordance
#> 165   Cellular Manifold & Clustering Concordance
#> 166   Cellular Manifold & Clustering Concordance
#> 167   Cellular Manifold & Clustering Concordance
#> 168   Cellular Manifold & Clustering Concordance
#> 169   Cellular Manifold & Clustering Concordance
#> 170  Technical Variation & Batch Mixing Dynamics
#> 172  Technical Variation & Batch Mixing Dynamics
#> 173  Technical Variation & Batch Mixing Dynamics
#> 174  Technical Variation & Batch Mixing Dynamics
#> 175  Technical Variation & Batch Mixing Dynamics
#> 176  Technical Variation & Batch Mixing Dynamics
#> 177  Technical Variation & Batch Mixing Dynamics
#> 178  Technical Variation & Batch Mixing Dynamics
#> 179  Technical Variation & Batch Mixing Dynamics
#> 180  Technical Variation & Batch Mixing Dynamics
#> 182    Biological Signal & Differential Fidelity
#> 183    Biological Signal & Differential Fidelity
#> 184    Biological Signal & Differential Fidelity
#> 185    Biological Signal & Differential Fidelity
#> 186    Biological Signal & Differential Fidelity
#> 187    Biological Signal & Differential Fidelity
#> 188    Biological Signal & Differential Fidelity
#> 189    Biological Signal & Differential Fidelity
#> 190    Biological Signal & Differential Fidelity
#> 192    Biological Signal & Differential Fidelity
#> 193    Biological Signal & Differential Fidelity
#> 194    Biological Signal & Differential Fidelity
#> 195          (VII) Trajectory & Lineage Topology
#> 196          (VII) Trajectory & Lineage Topology
#> 197          (VII) Trajectory & Lineage Topology
#> 198          (VII) Trajectory & Lineage Topology
#> 199          (VII) Trajectory & Lineage Topology
#> 1100     (VI) Multiomics Integration & Alignment
#> 1102     (VI) Multiomics Integration & Alignment
#> 1103     (VI) Multiomics Integration & Alignment
#> 1104     (VI) Multiomics Integration & Alignment
#> 1105     (VI) Multiomics Integration & Alignment
#>                               Property                         Metric
#> 1                         library_size                    Wasserstein
#> 2                         library_size                             KS
#> 3                         library_size                            MAD
#> 4                         library_size                  Bhattacharyya
#> 5                         library_size                 Runs_Statistic
#> 6                         library_size                    NN_Mismatch
#> 7                         library_size    Between_Dataset_Silh_Global
#> 8                         library_size                            MAE
#> 9                         library_size                           RMSE
#> 10                        library_size                        KDE_ISE
#> 11                    log_library_size                    Wasserstein
#> 12                    log_library_size                             KS
#> 13                    log_library_size                            MAD
#> 14                    log_library_size                  Bhattacharyya
#> 15                    log_library_size                 Runs_Statistic
#> 16                    log_library_size                    NN_Mismatch
#> 17                    log_library_size    Between_Dataset_Silh_Global
#> 18                    log_library_size                            MAE
#> 19                    log_library_size                           RMSE
#> 20                    log_library_size                        KDE_ISE
#> 21                          tmm_factor                    Wasserstein
#> 22                          tmm_factor                             KS
#> 23                          tmm_factor                            MAD
#> 24                          tmm_factor                  Bhattacharyya
#> 25                          tmm_factor                 Runs_Statistic
#> 26                          tmm_factor                    NN_Mismatch
#> 27                          tmm_factor    Between_Dataset_Silh_Global
#> 28                          tmm_factor                            MAE
#> 29                          tmm_factor                           RMSE
#> 30                          tmm_factor                        KDE_ISE
#> 31                  zero_fraction_cell                    Wasserstein
#> 32                  zero_fraction_cell                             KS
#> 33                  zero_fraction_cell                            MAD
#> 34                  zero_fraction_cell                  Bhattacharyya
#> 35                  zero_fraction_cell                 Runs_Statistic
#> 36                  zero_fraction_cell                    NN_Mismatch
#> 37                  zero_fraction_cell    Between_Dataset_Silh_Global
#> 38                  zero_fraction_cell                            MAE
#> 39                  zero_fraction_cell                           RMSE
#> 40                  zero_fraction_cell                        KDE_ISE
#> 41                 detection_freq_cell                    Wasserstein
#> 42                 detection_freq_cell                             KS
#> 43                 detection_freq_cell                            MAD
#> 44                 detection_freq_cell                  Bhattacharyya
#> 45                 detection_freq_cell                 Runs_Statistic
#> 46                 detection_freq_cell                    NN_Mismatch
#> 47                 detection_freq_cell    Between_Dataset_Silh_Global
#> 48                 detection_freq_cell                            MAE
#> 49                 detection_freq_cell                           RMSE
#> 50                 detection_freq_cell                        KDE_ISE
#> 51                            cell_cor                    Wasserstein
#> 52                            cell_cor                             KS
#> 53                            cell_cor                            MAD
#> 54                            cell_cor                  Bhattacharyya
#> 55                            cell_cor                 Runs_Statistic
#> 56                            cell_cor                    NN_Mismatch
#> 57                            cell_cor    Between_Dataset_Silh_Global
#> 58                            cell_cor                            MAE
#> 59                            cell_cor                           RMSE
#> 60                            cell_cor                        KDE_ISE
#> 61                            pca_dist                    Wasserstein
#> 62                            pca_dist                             KS
#> 63                            pca_dist                            MAD
#> 64                            pca_dist                  Bhattacharyya
#> 65                            pca_dist                 Runs_Statistic
#> 66                            pca_dist                    NN_Mismatch
#> 67                            pca_dist    Between_Dataset_Silh_Global
#> 68                            pca_dist                            MAE
#> 69                            pca_dist                           RMSE
#> 70                            pca_dist                        KDE_ISE
#> 71                     mean_expression                    Wasserstein
#> 72                     mean_expression                             KS
#> 73                     mean_expression                            MAD
#> 74                     mean_expression                  Bhattacharyya
#> 75                     mean_expression                 Runs_Statistic
#> 76                     mean_expression                    NN_Mismatch
#> 77                     mean_expression    Between_Dataset_Silh_Global
#> 78                     mean_expression                            MAE
#> 79                     mean_expression                           RMSE
#> 80                     mean_expression                        KDE_ISE
#> 81                            variance                    Wasserstein
#> 82                            variance                             KS
#> 83                            variance                            MAD
#> 84                            variance                  Bhattacharyya
#> 85                            variance                 Runs_Statistic
#> 86                            variance                    NN_Mismatch
#> 87                            variance    Between_Dataset_Silh_Global
#> 88                            variance                            MAE
#> 89                            variance                           RMSE
#> 90                            variance                        KDE_ISE
#> 91                     scaled_variance                    Wasserstein
#> 92                     scaled_variance                             KS
#> 93                     scaled_variance                            MAD
#> 94                     scaled_variance                  Bhattacharyya
#> 95                     scaled_variance                 Runs_Statistic
#> 96                     scaled_variance                    NN_Mismatch
#> 97                     scaled_variance    Between_Dataset_Silh_Global
#> 98                     scaled_variance                            MAE
#> 99                     scaled_variance                           RMSE
#> 100                    scaled_variance                        KDE_ISE
#> 101                                 cv                    Wasserstein
#> 102                                 cv                             KS
#> 103                                 cv                            MAD
#> 104                                 cv                  Bhattacharyya
#> 105                                 cv                 Runs_Statistic
#> 106                                 cv                    NN_Mismatch
#> 107                                 cv    Between_Dataset_Silh_Global
#> 108                                 cv                            MAE
#> 109                                 cv                           RMSE
#> 110                                 cv                        KDE_ISE
#> 111              zero_fraction_feature                    Wasserstein
#> 112              zero_fraction_feature                             KS
#> 113              zero_fraction_feature                            MAD
#> 114              zero_fraction_feature                  Bhattacharyya
#> 115              zero_fraction_feature                 Runs_Statistic
#> 116              zero_fraction_feature                    NN_Mismatch
#> 117              zero_fraction_feature    Between_Dataset_Silh_Global
#> 118              zero_fraction_feature                            MAE
#> 119              zero_fraction_feature                           RMSE
#> 120              zero_fraction_feature                        KDE_ISE
#> 121                         dispersion                    Wasserstein
#> 122                         dispersion                             KS
#> 123                         dispersion                            MAD
#> 124                         dispersion                  Bhattacharyya
#> 125                         dispersion                 Runs_Statistic
#> 126                         dispersion                    NN_Mismatch
#> 127                         dispersion    Between_Dataset_Silh_Global
#> 128                         dispersion                            MAE
#> 129                         dispersion                           RMSE
#> 130                         dispersion                        KDE_ISE
#> 131                           gene_cor                    Wasserstein
#> 132                           gene_cor                             KS
#> 133                           gene_cor                            MAD
#> 134                           gene_cor                  Bhattacharyya
#> 135                           gene_cor                 Runs_Statistic
#> 136                           gene_cor                    NN_Mismatch
#> 137                           gene_cor    Between_Dataset_Silh_Global
#> 138                           gene_cor                            MAE
#> 139                           gene_cor                           RMSE
#> 140                           gene_cor                        KDE_ISE
#> 141                      dropout_curve                    slope_error
#> 142                      dropout_curve                 midpoint_error
#> 143                      dropout_curve                 r2_discrepancy
#> 144                excess_zero_weights excess_zero_weight_discrepancy
#> 145                           manifold                            MMD
#> 146                           manifold           Generative_Precision
#> 147                           manifold              Generative_Recall
#> 151                       batch_mixing               batch_silhouette
#> 152                       batch_mixing                         pcr_r2
#> 153                       batch_mixing                            cms
#> 154                       batch_mixing                            isi
#> 155                       batch_mixing                 cross_batch_F1
#> 156      scRNA_differential_expression                pval_dist_score
#> 157      scRNA_differential_expression                       accuracy
#> 158      scRNA_differential_expression                      precision
#> 159      scRNA_differential_expression                         recall
#> 160      scRNA_differential_expression                       f1_score
#> 161      scRNA_differential_expression                            auc
#> 162                       library_size                    Wasserstein
#> 210                       library_size                             KS
#> 310                       library_size                            MAD
#> 410                       library_size                  Bhattacharyya
#> 510                       library_size                 Runs_Statistic
#> 610                       library_size                    NN_Mismatch
#> 710                       library_size    Between_Dataset_Silh_Global
#> 810                       library_size                            MAE
#> 910                       library_size                           RMSE
#> 1010                      library_size                        KDE_ISE
#> 1110                  log_library_size                    Wasserstein
#> 1210                  log_library_size                             KS
#> 1310                  log_library_size                            MAD
#> 1410                  log_library_size                  Bhattacharyya
#> 1510                  log_library_size                 Runs_Statistic
#> 163                   log_library_size                    NN_Mismatch
#> 171                   log_library_size    Between_Dataset_Silh_Global
#> 181                   log_library_size                            MAE
#> 191                   log_library_size                           RMSE
#> 201                   log_library_size                        KDE_ISE
#> 211                 zero_fraction_cell                    Wasserstein
#> 221                 zero_fraction_cell                             KS
#> 231                 zero_fraction_cell                            MAD
#> 241                 zero_fraction_cell                  Bhattacharyya
#> 251                 zero_fraction_cell                 Runs_Statistic
#> 261                 zero_fraction_cell                    NN_Mismatch
#> 271                 zero_fraction_cell    Between_Dataset_Silh_Global
#> 281                 zero_fraction_cell                            MAE
#> 291                 zero_fraction_cell                           RMSE
#> 301                 zero_fraction_cell                        KDE_ISE
#> 311                detection_freq_cell                    Wasserstein
#> 321                detection_freq_cell                             KS
#> 331                detection_freq_cell                            MAD
#> 341                detection_freq_cell                  Bhattacharyya
#> 351                detection_freq_cell                 Runs_Statistic
#> 361                detection_freq_cell                    NN_Mismatch
#> 371                detection_freq_cell    Between_Dataset_Silh_Global
#> 381                detection_freq_cell                            MAE
#> 391                detection_freq_cell                           RMSE
#> 401                detection_freq_cell                        KDE_ISE
#> 411                           cell_cor                    Wasserstein
#> 421                           cell_cor                             KS
#> 431                           cell_cor                            MAD
#> 441                           cell_cor                  Bhattacharyya
#> 451                           cell_cor                 Runs_Statistic
#> 461                           cell_cor                    NN_Mismatch
#> 471                           cell_cor    Between_Dataset_Silh_Global
#> 481                           cell_cor                            MAE
#> 491                           cell_cor                           RMSE
#> 501                           cell_cor                        KDE_ISE
#> 511                           pca_dist                    Wasserstein
#> 521                           pca_dist                             KS
#> 531                           pca_dist                            MAD
#> 541                           pca_dist                  Bhattacharyya
#> 551                           pca_dist                 Runs_Statistic
#> 561                           pca_dist                    NN_Mismatch
#> 571                           pca_dist    Between_Dataset_Silh_Global
#> 581                           pca_dist                            MAE
#> 591                           pca_dist                           RMSE
#> 601                           pca_dist                        KDE_ISE
#> 611                    mean_expression                    Wasserstein
#> 621                    mean_expression                             KS
#> 631                    mean_expression                            MAD
#> 641                    mean_expression                  Bhattacharyya
#> 651                    mean_expression                 Runs_Statistic
#> 661                    mean_expression                    NN_Mismatch
#> 671                    mean_expression    Between_Dataset_Silh_Global
#> 681                    mean_expression                            MAE
#> 691                    mean_expression                           RMSE
#> 701                    mean_expression                        KDE_ISE
#> 711                           variance                    Wasserstein
#> 721                           variance                             KS
#> 731                           variance                            MAD
#> 741                           variance                  Bhattacharyya
#> 751                           variance                 Runs_Statistic
#> 761                           variance                    NN_Mismatch
#> 771                           variance    Between_Dataset_Silh_Global
#> 781                           variance                            MAE
#> 791                           variance                           RMSE
#> 801                           variance                        KDE_ISE
#> 811                    scaled_variance                    Wasserstein
#> 821                    scaled_variance                             KS
#> 831                    scaled_variance                            MAD
#> 841                    scaled_variance                  Bhattacharyya
#> 851                    scaled_variance                 Runs_Statistic
#> 861                    scaled_variance                    NN_Mismatch
#> 871                    scaled_variance    Between_Dataset_Silh_Global
#> 881                    scaled_variance                            MAE
#> 891                    scaled_variance                           RMSE
#> 901                    scaled_variance                        KDE_ISE
#> 911                                 cv                    Wasserstein
#> 921                                 cv                             KS
#> 931                                 cv                            MAD
#> 941                                 cv                  Bhattacharyya
#> 951                                 cv                 Runs_Statistic
#> 961                                 cv                    NN_Mismatch
#> 971                                 cv    Between_Dataset_Silh_Global
#> 981                                 cv                            MAE
#> 991                                 cv                           RMSE
#> 1001                                cv                        KDE_ISE
#> 1011             zero_fraction_feature                    Wasserstein
#> 1021             zero_fraction_feature                             KS
#> 1031             zero_fraction_feature                            MAD
#> 1041             zero_fraction_feature                  Bhattacharyya
#> 1051             zero_fraction_feature                 Runs_Statistic
#> 1061             zero_fraction_feature                    NN_Mismatch
#> 1071             zero_fraction_feature    Between_Dataset_Silh_Global
#> 1081             zero_fraction_feature                            MAE
#> 1091             zero_fraction_feature                           RMSE
#> 1101             zero_fraction_feature                        KDE_ISE
#> 1111                        dispersion                    Wasserstein
#> 1121                        dispersion                             KS
#> 1131                        dispersion                            MAD
#> 1141                        dispersion                  Bhattacharyya
#> 1151                        dispersion                 Runs_Statistic
#> 1161                        dispersion                    NN_Mismatch
#> 1171                        dispersion    Between_Dataset_Silh_Global
#> 1181                        dispersion                            MAE
#> 1191                        dispersion                           RMSE
#> 1201                        dispersion                        KDE_ISE
#> 1211                          gene_cor                    Wasserstein
#> 1221                          gene_cor                             KS
#> 1231                          gene_cor                            MAD
#> 1241                          gene_cor                  Bhattacharyya
#> 1251                          gene_cor                 Runs_Statistic
#> 1261                          gene_cor                    NN_Mismatch
#> 1271                          gene_cor    Between_Dataset_Silh_Global
#> 1281                          gene_cor                            MAE
#> 1291                          gene_cor                           RMSE
#> 1301                          gene_cor                        KDE_ISE
#> 1311      accessibility_sparsity_curve                     curve_rmse
#> 1321      accessibility_sparsity_curve                 r2_discrepancy
#> 1331               excess_zero_weights excess_zero_weight_discrepancy
#> 1341                          manifold                            MMD
#> 1351                          manifold           Generative_Precision
#> 1361                          manifold              Generative_Recall
#> 1401                      batch_mixing               batch_silhouette
#> 1411                      batch_mixing                         pcr_r2
#> 1421                      batch_mixing                            cms
#> 1431                      batch_mixing                            isi
#> 1441                      batch_mixing                 cross_batch_F1
#> 1451 scATAC_differential_accessibility                pval_dist_score
#> 1461 scATAC_differential_accessibility                       accuracy
#> 1471 scATAC_differential_accessibility                      precision
#> 148  scATAC_differential_accessibility                         recall
#> 149  scATAC_differential_accessibility                       f1_score
#> 150  scATAC_differential_accessibility                            auc
#> 164                         clustering                 silhouette_sim
#> 165                         clustering                 silhouette_sim
#> 166                         clustering                            cdi
#> 167                         clustering                            cdi
#> 168                         clustering            neighborhood_purity
#> 169                         clustering            neighborhood_purity
#> 170                       batch_mixing               batch_silhouette
#> 172                       batch_mixing                         pcr_r2
#> 173                       batch_mixing                            cms
#> 174                       batch_mixing                            isi
#> 175                       batch_mixing                 cross_batch_F1
#> 176                       batch_mixing               batch_silhouette
#> 177                       batch_mixing                         pcr_r2
#> 178                       batch_mixing                            cms
#> 179                       batch_mixing                            isi
#> 180                       batch_mixing                 cross_batch_F1
#> 182      scRNA_differential_expression                pval_dist_score
#> 183      scRNA_differential_expression                       accuracy
#> 184      scRNA_differential_expression                      precision
#> 185      scRNA_differential_expression                         recall
#> 186      scRNA_differential_expression                       f1_score
#> 187      scRNA_differential_expression                            auc
#> 188  scATAC_differential_accessibility                pval_dist_score
#> 189  scATAC_differential_accessibility                       accuracy
#> 190  scATAC_differential_accessibility                      precision
#> 192  scATAC_differential_accessibility                         recall
#> 193  scATAC_differential_accessibility                       f1_score
#> 194  scATAC_differential_accessibility                            auc
#> 195                trajectory_topology                            him
#> 196                trajectory_topology                      edge_flip
#> 197                  geodesic_ordering           geodesic_correlation
#> 198                 branch_concordance                    F1_branches
#> 199                   dynamic_features                  wcor_features
#> 1100            multiomics_integration                      batch_asw
#> 1102            multiomics_integration                  cell_type_asw
#> 1103            multiomics_integration             graph_connectivity
#> 1104            multiomics_integration                            ari
#> 1105            multiomics_integration                            nmi
#>              Value            Modality
#> 1     1.506250e+01                 rna
#> 2     2.125000e-01                 rna
#> 3     1.500000e+01                 rna
#> 4     1.188839e-04                 rna
#> 5    -1.744748e+00                 rna
#> 6     1.750000e-01                 rna
#> 7     2.247078e-02                 rna
#> 8     1.506250e+01                 rna
#> 9     1.571504e+01                 rna
#> 10    7.980370e-04                 rna
#> 11    6.194965e-02                 rna
#> 12    2.125000e-01                 rna
#> 13    6.363671e-02                 rna
#> 14    3.834591e-06                 rna
#> 15   -1.744748e+00                 rna
#> 16    1.875000e-01                 rna
#> 17    2.176201e-02                 rna
#> 18    6.194965e-02                 rna
#> 19    6.566950e-02                 rna
#> 20    1.962832e-01                 rna
#> 21    3.268387e-02                 rna
#> 22    1.018380e-01                 rna
#> 23    2.555864e-02                 rna
#> 24    2.186586e-04                 rna
#> 25   -9.516807e-01                 rna
#> 26    5.000000e-02                 rna
#> 27   -6.195226e-04                 rna
#> 28    3.268387e-02                 rna
#> 29    4.613293e-02                 rna
#> 30    7.601569e-02                 rna
#> 31    1.562500e-02                 rna
#> 32    1.250000e-01                 rna
#> 33    1.666667e-02                 rna
#> 34    8.854773e-04                 rna
#> 35   -9.358194e+00                 rna
#> 36    6.875000e-01                 rna
#> 37    7.524680e-03                 rna
#> 38    1.562500e-02                 rna
#> 39    1.980811e-02                 rna
#> 40    3.253205e-01                 rna
#> 41    1.562500e-02                 rna
#> 42    1.250000e-01                 rna
#> 43    1.666667e-02                 rna
#> 44    2.679272e-05                 rna
#> 45   -9.040967e+00                 rna
#> 46    6.875000e-01                 rna
#> 47    7.524680e-03                 rna
#> 48    1.562500e-02                 rna
#> 49    1.980811e-02                 rna
#> 50    3.253205e-01                 rna
#> 51    5.728994e-03                 rna
#> 52    3.323906e-02                 rna
#> 53    5.684572e-03                 rna
#> 54    2.566053e-05                 rna
#> 55   -2.188896e+00                 rna
#> 56    1.666667e-03                 rna
#> 57    1.423443e-03                 rna
#> 58    5.728994e-03                 rna
#> 59    6.328019e-03                 rna
#> 60    6.120211e-03                 rna
#> 61    1.594340e+00                 rna
#> 62    8.747786e-02                 rna
#> 63    1.428096e+00                 rna
#> 64    3.105370e-05                 rna
#> 65    2.515972e-02                 rna
#> 66    1.500000e-02                 rna
#> 67    9.907338e-03                 rna
#> 68    1.594340e+00                 rna
#> 69    1.676818e+00                 rna
#> 70    8.585380e-04                 rna
#> 71    2.473408e-01                 rna
#> 72    2.055826e-01                 rna
#> 73    2.571903e-01                 rna
#> 74    3.690266e-05                 rna
#> 75    3.666923e-01                 rna
#> 76    5.000000e-02                 rna
#> 77   -9.070018e-05                 rna
#> 78    2.473408e-01                 rna
#> 79    2.887192e-01                 rna
#> 80    1.208325e-02                 rna
#> 81    2.036911e+00                 rna
#> 82    1.533879e-01                 rna
#> 83    1.969816e+00                 rna
#> 84    3.196235e-04                 rna
#> 85   -3.666923e-01                 rna
#> 86    4.166667e-02                 rna
#> 87    9.752197e-04                 rna
#> 88    2.036911e+00                 rna
#> 89    2.271370e+00                 rna
#> 90    1.659081e-03                 rna
#> 91    1.119261e-01                 rna
#> 92    6.742514e-02                 rna
#> 93    7.339222e-02                 rna
#> 94    2.010495e-03                 rna
#> 95   -1.283423e+00                 rna
#> 96    1.666667e-02                 rna
#> 97   -1.234822e-02                 rna
#> 98    1.119261e-01                 rna
#> 99    1.504922e-01                 rna
#> 100   6.265612e-03                 rna
#> 101   7.432972e+00                 rna
#> 102   2.906651e-01                 rna
#> 103   7.484040e+00                 rna
#> 104   7.165290e-05                 rna
#> 105  -1.100077e+00                 rna
#> 106   4.166667e-02                 rna
#> 107   5.142707e-02                 rna
#> 108   7.432972e+00                 rna
#> 109   7.844679e+00                 rna
#> 110   4.198168e-03                 rna
#> 111   1.687500e-02                 rna
#> 112   1.333333e-01                 rna
#> 113   1.250000e-02                 rna
#> 114   6.100311e-04                 rna
#> 115  -4.400308e+00                 rna
#> 116   9.166667e-02                 rna
#> 117  -1.887003e-03                 rna
#> 118   1.687500e-02                 rna
#> 119   1.916214e-02                 rna
#> 120   1.732278e-01                 rna
#> 121   4.456931e-01                 rna
#> 122   1.741354e-01                 rna
#> 123   4.388661e-01                 rna
#> 124   5.586790e-04                 rna
#> 125   1.466769e+00                 rna
#> 126   2.500000e-02                 rna
#> 127   1.024133e-02                 rna
#> 128   4.456931e-01                 rna
#> 129   4.890173e-01                 rna
#> 130   1.955467e-02                 rna
#> 131   3.299875e-03                 rna
#> 132   2.672058e-02                 rna
#> 133   2.118671e-03                 rna
#> 134   6.113770e-05                 rna
#> 135  -2.050782e+00                 rna
#> 136   3.000000e-02                 rna
#> 137  -2.482463e-06                 rna
#> 138   3.299875e-03                 rna
#> 139   5.067781e-03                 rna
#> 140   1.403010e-02                 rna
#> 141   2.560909e-01                 rna
#> 142   2.624144e-01                 rna
#> 143   1.646825e-01                 rna
#> 144   3.712431e-02                 rna
#> 145   1.103894e-01                 rna
#> 146   7.625000e-01                 rna
#> 147   9.625000e-01                 rna
#> 151  -2.710457e-03                 rna
#> 152   9.476487e-03                 rna
#> 153   0.000000e+00                 rna
#> 154   1.938397e+00                 rna
#> 155   5.126418e-01                 rna
#> 156   8.719073e-01                 rna
#> 157   9.500000e-01                 rna
#> 158   0.000000e+00                 rna
#> 159   0.000000e+00                 rna
#> 160   0.000000e+00                 rna
#> 161   6.465517e-01                 rna
#> 162   2.337500e+00                atac
#> 210   1.750000e-01                atac
#> 310   2.000000e+00                atac
#> 410   1.055641e-04                atac
#> 510  -6.027311e+00                atac
#> 610   1.937500e-01                atac
#> 710   2.188673e-02                atac
#> 810   2.337500e+00                atac
#> 910   2.487469e+00                atac
#> 1010  4.478217e-03                atac
#> 1110  7.991730e-02                atac
#> 1210  1.750000e-01                atac
#> 1310  8.223810e-02                atac
#> 1410  1.155313e-05                atac
#> 1510 -6.027311e+00                atac
#> 163   1.937500e-01                atac
#> 171   1.990590e-02                atac
#> 181   7.991730e-02                atac
#> 191   8.389317e-02                atac
#> 201   1.268716e-01                atac
#> 211   2.750000e-02                atac
#> 221   2.125000e-01                atac
#> 231   3.333333e-02                atac
#> 241   6.907067e-05                atac
#> 251  -7.296219e+00                atac
#> 261   2.812500e-01                atac
#> 271   2.566981e-02                atac
#> 281   2.750000e-02                atac
#> 291   3.016160e-02                atac
#> 301   4.176638e-01                atac
#> 311   2.750000e-02                atac
#> 321   2.125000e-01                atac
#> 331   3.333333e-02                atac
#> 341   1.795082e-04                atac
#> 351  -6.978992e+00                atac
#> 361   2.812500e-01                atac
#> 371   2.566981e-02                atac
#> 381   2.750000e-02                atac
#> 391   3.016160e-02                atac
#> 401   4.176638e-01                atac
#> 411   3.342856e-03                atac
#> 421   2.056815e-02                atac
#> 431   2.270508e-03                atac
#> 441   5.400110e-05                atac
#> 451  -2.515972e-01                atac
#> 461   4.333333e-02                atac
#> 471  -2.274271e-03                atac
#> 481   3.342856e-03                atac
#> 491   6.619754e-03                atac
#> 501   4.263335e-03                atac
#> 511   1.502554e-01                atac
#> 521   2.948321e-02                atac
#> 531   1.327661e-01                atac
#> 541   2.878871e-06                atac
#> 551   7.547916e-01                atac
#> 561   2.666667e-02                atac
#> 571   6.595987e-04                atac
#> 581   1.502554e-01                atac
#> 591   1.883752e-01                atac
#> 601   3.459236e-04                atac
#> 611   1.501367e-01                atac
#> 621   1.879472e-01                atac
#> 631   1.536693e-01                atac
#> 641   2.959894e-04                atac
#> 651  -1.833462e+00                atac
#> 661   2.500000e-02                atac
#> 671   3.291150e-03                atac
#> 681   1.501367e-01                atac
#> 691   1.677762e-01                atac
#> 701   2.809910e-02                atac
#> 711   2.170317e-01                atac
#> 721   2.753225e-01                atac
#> 731   2.419452e-01                atac
#> 741   5.389967e-05                atac
#> 751  -2.566846e+00                atac
#> 761   1.166667e-01                atac
#> 771   3.168098e-02                atac
#> 781   2.170317e-01                atac
#> 791   2.541541e-01                atac
#> 801   1.429083e-01                atac
#> 811   1.993187e-01                atac
#> 821   1.174387e-01                atac
#> 831   1.419895e-01                atac
#> 841   2.917823e-03                atac
#> 851  -9.167309e-01                atac
#> 861   8.333333e-02                atac
#> 871  -3.112124e-03                atac
#> 881   1.993187e-01                atac
#> 891   2.711405e-01                atac
#> 901   3.669152e-02                atac
#> 911   8.864353e+00                atac
#> 921   1.763639e-01                atac
#> 931   8.752987e+00                atac
#> 941   2.185211e-04                atac
#> 951  -1.833462e+00                atac
#> 961   3.333333e-02                atac
#> 971   6.571696e-03                atac
#> 981   8.864353e+00                atac
#> 991   1.046915e+01                atac
#> 1001  7.968311e-04                atac
#> 1011  2.958333e-02                atac
#> 1021  2.000000e-01                atac
#> 1031  2.500000e-02                atac
#> 1041  1.230865e-04                atac
#> 1051 -1.650116e+00                atac
#> 1061  5.000000e-02                atac
#> 1071  7.445303e-03                atac
#> 1081  2.958333e-02                atac
#> 1091  3.307189e-02                atac
#> 1101  1.835132e-01                atac
#> 1111  2.544810e-02                atac
#> 1121  1.136364e-01                atac
#> 1131  2.119535e-02                atac
#> 1141  1.193835e-04                atac
#> 1151 -5.500385e-01                atac
#> 1161  8.333333e-03                atac
#> 1171 -7.584339e-04                atac
#> 1181  2.544810e-02                atac
#> 1191  2.978845e-02                atac
#> 1201  1.190016e-01                atac
#> 1211  2.718776e-03                atac
#> 1221  1.759557e-02                atac
#> 1231  1.838118e-03                atac
#> 1241  8.878612e-05                atac
#> 1251 -6.723876e-02                atac
#> 1261  4.666667e-02                atac
#> 1271 -2.426294e-04                atac
#> 1281  2.718776e-03                atac
#> 1291  5.260476e-03                atac
#> 1301  6.132365e-03                atac
#> 1311  1.233884e-02                atac
#> 1321  8.535249e-04                atac
#> 1331  8.480348e-04                atac
#> 1341  1.096098e-01                atac
#> 1351  1.000000e+00                atac
#> 1361  9.500000e-01                atac
#> 1401  8.780031e-04                atac
#> 1411  1.375709e-02                atac
#> 1421  0.000000e+00                atac
#> 1431  1.931492e+00                atac
#> 1441  7.346700e-01                atac
#> 1451  8.563379e-01                atac
#> 1461  8.166667e-01                atac
#> 1471  5.000000e-01                atac
#> 148   9.090909e-02                atac
#> 149   1.538462e-01                atac
#> 150   8.367347e-01                atac
#> 164   3.402312e-03                 rna
#> 165   4.149183e-02                atac
#> 166   4.984432e+00                 rna
#> 167   1.687605e+00                atac
#> 168   4.775000e-01                 rna
#> 169   7.275000e-01                atac
#> 170  -2.710457e-03                 rna
#> 172   9.476487e-03                 rna
#> 173   0.000000e+00                 rna
#> 174   1.938397e+00                 rna
#> 175   5.126418e-01                 rna
#> 176   8.780031e-04                atac
#> 177   1.375709e-02                atac
#> 178   0.000000e+00                atac
#> 179   1.931492e+00                atac
#> 180   7.346700e-01                atac
#> 182   8.719073e-01                 rna
#> 183   9.500000e-01                 rna
#> 184   0.000000e+00                 rna
#> 185   0.000000e+00                 rna
#> 186   0.000000e+00                 rna
#> 187   6.465517e-01                 rna
#> 188   8.563379e-01                atac
#> 189   8.166667e-01                atac
#> 190   5.000000e-01                atac
#> 192   9.090909e-02                atac
#> 193   1.538462e-01                atac
#> 194   8.367347e-01                atac
#> 195   1.000000e+00                 rna
#> 196   1.000000e+00                 rna
#> 197   4.848647e-01                 rna
#> 198   1.000000e+00                 rna
#> 199   3.802351e-01                 rna
#> 1100  9.747094e-01 Joint (Cross-Modal)
#> 1102  5.582026e-01 Joint (Cross-Modal)
#> 1103  1.000000e+00 Joint (Cross-Modal)
#> 1104  6.765538e-01 Joint (Cross-Modal)
#> 1105  5.821401e-01 Joint (Cross-Modal)
#> 
#> $metrics_summary_table
#>                                         Category
#> 1                      Distributional Properties
#> 2                      Distributional Properties
#> 3                      Distributional Properties
#> 4                      Distributional Properties
#> 5                      Distributional Properties
#> 6                      Distributional Properties
#> 7                      Distributional Properties
#> 8                      Distributional Properties
#> 9                      Distributional Properties
#> 10                     Distributional Properties
#> 11                     Distributional Properties
#> 12                     Distributional Properties
#> 13                     Distributional Properties
#> 14                     Distributional Properties
#> 15                     Distributional Properties
#> 16                     Distributional Properties
#> 17                     Distributional Properties
#> 18                     Distributional Properties
#> 19                     Distributional Properties
#> 20                     Distributional Properties
#> 21                     Distributional Properties
#> 22                     Distributional Properties
#> 23                     Distributional Properties
#> 24                     Distributional Properties
#> 25                     Distributional Properties
#> 26                     Distributional Properties
#> 27                     Distributional Properties
#> 28                     Distributional Properties
#> 29                     Distributional Properties
#> 30                     Distributional Properties
#> 31                     Distributional Properties
#> 32                     Distributional Properties
#> 33                     Distributional Properties
#> 34                     Distributional Properties
#> 35                     Distributional Properties
#> 36                     Distributional Properties
#> 37                     Distributional Properties
#> 38                     Distributional Properties
#> 39                     Distributional Properties
#> 40                     Distributional Properties
#> 41                     Distributional Properties
#> 42                     Distributional Properties
#> 43                     Distributional Properties
#> 44                     Distributional Properties
#> 45                     Distributional Properties
#> 46                     Distributional Properties
#> 47                     Distributional Properties
#> 48                     Distributional Properties
#> 49                     Distributional Properties
#> 50                     Distributional Properties
#> 51                     Distributional Properties
#> 52                     Distributional Properties
#> 53                     Distributional Properties
#> 54                     Distributional Properties
#> 55                     Distributional Properties
#> 56                     Distributional Properties
#> 57                     Distributional Properties
#> 58                     Distributional Properties
#> 59                     Distributional Properties
#> 60                     Distributional Properties
#> 61                     Distributional Properties
#> 62                     Distributional Properties
#> 63                     Distributional Properties
#> 64                     Distributional Properties
#> 65                     Distributional Properties
#> 66                     Distributional Properties
#> 67                     Distributional Properties
#> 68                     Distributional Properties
#> 69                     Distributional Properties
#> 70                     Distributional Properties
#> 71                     Distributional Properties
#> 72                     Distributional Properties
#> 73                     Distributional Properties
#> 74                     Distributional Properties
#> 75                     Distributional Properties
#> 76                     Distributional Properties
#> 77                     Distributional Properties
#> 78                     Distributional Properties
#> 79                     Distributional Properties
#> 80                     Distributional Properties
#> 81                     Distributional Properties
#> 82                     Distributional Properties
#> 83                     Distributional Properties
#> 84                     Distributional Properties
#> 85                     Distributional Properties
#> 86                     Distributional Properties
#> 87                     Distributional Properties
#> 88                     Distributional Properties
#> 89                     Distributional Properties
#> 90                     Distributional Properties
#> 91                     Distributional Properties
#> 92                     Distributional Properties
#> 93                     Distributional Properties
#> 94                     Distributional Properties
#> 95                     Distributional Properties
#> 96                     Distributional Properties
#> 97                     Distributional Properties
#> 98                     Distributional Properties
#> 99                     Distributional Properties
#> 100                    Distributional Properties
#> 101                    Distributional Properties
#> 102                    Distributional Properties
#> 103                    Distributional Properties
#> 104                    Distributional Properties
#> 105                    Distributional Properties
#> 106                    Distributional Properties
#> 107                    Distributional Properties
#> 108                    Distributional Properties
#> 109                    Distributional Properties
#> 110                    Distributional Properties
#> 111                    Distributional Properties
#> 112                    Distributional Properties
#> 113                    Distributional Properties
#> 114                    Distributional Properties
#> 115                    Distributional Properties
#> 116                    Distributional Properties
#> 117                    Distributional Properties
#> 118                    Distributional Properties
#> 119                    Distributional Properties
#> 120                    Distributional Properties
#> 121                    Distributional Properties
#> 122                    Distributional Properties
#> 123                    Distributional Properties
#> 124                    Distributional Properties
#> 125                    Distributional Properties
#> 126                    Distributional Properties
#> 127                    Distributional Properties
#> 128                    Distributional Properties
#> 129                    Distributional Properties
#> 130                    Distributional Properties
#> 131                    Distributional Properties
#> 132                    Distributional Properties
#> 133                    Distributional Properties
#> 134                    Distributional Properties
#> 135                    Distributional Properties
#> 136                    Distributional Properties
#> 137                    Distributional Properties
#> 138                    Distributional Properties
#> 139                    Distributional Properties
#> 140                    Distributional Properties
#> 141           Sparsity & Zero-Inflation Dynamics
#> 142           Sparsity & Zero-Inflation Dynamics
#> 143           Sparsity & Zero-Inflation Dynamics
#> 144           Sparsity & Zero-Inflation Dynamics
#> 145   Cellular Manifold & Clustering Concordance
#> 146   Cellular Manifold & Clustering Concordance
#> 147   Cellular Manifold & Clustering Concordance
#> 151  Technical Variation & Batch Mixing Dynamics
#> 152  Technical Variation & Batch Mixing Dynamics
#> 153  Technical Variation & Batch Mixing Dynamics
#> 154  Technical Variation & Batch Mixing Dynamics
#> 155  Technical Variation & Batch Mixing Dynamics
#> 156    Biological Signal & Differential Fidelity
#> 157    Biological Signal & Differential Fidelity
#> 158    Biological Signal & Differential Fidelity
#> 159    Biological Signal & Differential Fidelity
#> 160    Biological Signal & Differential Fidelity
#> 161    Biological Signal & Differential Fidelity
#> 162                    Distributional Properties
#> 210                    Distributional Properties
#> 310                    Distributional Properties
#> 410                    Distributional Properties
#> 510                    Distributional Properties
#> 610                    Distributional Properties
#> 710                    Distributional Properties
#> 810                    Distributional Properties
#> 910                    Distributional Properties
#> 1010                   Distributional Properties
#> 1110                   Distributional Properties
#> 1210                   Distributional Properties
#> 1310                   Distributional Properties
#> 1410                   Distributional Properties
#> 1510                   Distributional Properties
#> 163                    Distributional Properties
#> 171                    Distributional Properties
#> 181                    Distributional Properties
#> 191                    Distributional Properties
#> 201                    Distributional Properties
#> 211                    Distributional Properties
#> 221                    Distributional Properties
#> 231                    Distributional Properties
#> 241                    Distributional Properties
#> 251                    Distributional Properties
#> 261                    Distributional Properties
#> 271                    Distributional Properties
#> 281                    Distributional Properties
#> 291                    Distributional Properties
#> 301                    Distributional Properties
#> 311                    Distributional Properties
#> 321                    Distributional Properties
#> 331                    Distributional Properties
#> 341                    Distributional Properties
#> 351                    Distributional Properties
#> 361                    Distributional Properties
#> 371                    Distributional Properties
#> 381                    Distributional Properties
#> 391                    Distributional Properties
#> 401                    Distributional Properties
#> 411                    Distributional Properties
#> 421                    Distributional Properties
#> 431                    Distributional Properties
#> 441                    Distributional Properties
#> 451                    Distributional Properties
#> 461                    Distributional Properties
#> 471                    Distributional Properties
#> 481                    Distributional Properties
#> 491                    Distributional Properties
#> 501                    Distributional Properties
#> 511                    Distributional Properties
#> 521                    Distributional Properties
#> 531                    Distributional Properties
#> 541                    Distributional Properties
#> 551                    Distributional Properties
#> 561                    Distributional Properties
#> 571                    Distributional Properties
#> 581                    Distributional Properties
#> 591                    Distributional Properties
#> 601                    Distributional Properties
#> 611                    Distributional Properties
#> 621                    Distributional Properties
#> 631                    Distributional Properties
#> 641                    Distributional Properties
#> 651                    Distributional Properties
#> 661                    Distributional Properties
#> 671                    Distributional Properties
#> 681                    Distributional Properties
#> 691                    Distributional Properties
#> 701                    Distributional Properties
#> 711                    Distributional Properties
#> 721                    Distributional Properties
#> 731                    Distributional Properties
#> 741                    Distributional Properties
#> 751                    Distributional Properties
#> 761                    Distributional Properties
#> 771                    Distributional Properties
#> 781                    Distributional Properties
#> 791                    Distributional Properties
#> 801                    Distributional Properties
#> 811                    Distributional Properties
#> 821                    Distributional Properties
#> 831                    Distributional Properties
#> 841                    Distributional Properties
#> 851                    Distributional Properties
#> 861                    Distributional Properties
#> 871                    Distributional Properties
#> 881                    Distributional Properties
#> 891                    Distributional Properties
#> 901                    Distributional Properties
#> 911                    Distributional Properties
#> 921                    Distributional Properties
#> 931                    Distributional Properties
#> 941                    Distributional Properties
#> 951                    Distributional Properties
#> 961                    Distributional Properties
#> 971                    Distributional Properties
#> 981                    Distributional Properties
#> 991                    Distributional Properties
#> 1001                   Distributional Properties
#> 1011                   Distributional Properties
#> 1021                   Distributional Properties
#> 1031                   Distributional Properties
#> 1041                   Distributional Properties
#> 1051                   Distributional Properties
#> 1061                   Distributional Properties
#> 1071                   Distributional Properties
#> 1081                   Distributional Properties
#> 1091                   Distributional Properties
#> 1101                   Distributional Properties
#> 1111                   Distributional Properties
#> 1121                   Distributional Properties
#> 1131                   Distributional Properties
#> 1141                   Distributional Properties
#> 1151                   Distributional Properties
#> 1161                   Distributional Properties
#> 1171                   Distributional Properties
#> 1181                   Distributional Properties
#> 1191                   Distributional Properties
#> 1201                   Distributional Properties
#> 1211                   Distributional Properties
#> 1221                   Distributional Properties
#> 1231                   Distributional Properties
#> 1241                   Distributional Properties
#> 1251                   Distributional Properties
#> 1261                   Distributional Properties
#> 1271                   Distributional Properties
#> 1281                   Distributional Properties
#> 1291                   Distributional Properties
#> 1301                   Distributional Properties
#> 1311          Sparsity & Zero-Inflation Dynamics
#> 1321          Sparsity & Zero-Inflation Dynamics
#> 1331          Sparsity & Zero-Inflation Dynamics
#> 1341  Cellular Manifold & Clustering Concordance
#> 1351  Cellular Manifold & Clustering Concordance
#> 1361  Cellular Manifold & Clustering Concordance
#> 1401 Technical Variation & Batch Mixing Dynamics
#> 1411 Technical Variation & Batch Mixing Dynamics
#> 1421 Technical Variation & Batch Mixing Dynamics
#> 1431 Technical Variation & Batch Mixing Dynamics
#> 1441 Technical Variation & Batch Mixing Dynamics
#> 1451   Biological Signal & Differential Fidelity
#> 1461   Biological Signal & Differential Fidelity
#> 1471   Biological Signal & Differential Fidelity
#> 148    Biological Signal & Differential Fidelity
#> 149    Biological Signal & Differential Fidelity
#> 150    Biological Signal & Differential Fidelity
#> 164   Cellular Manifold & Clustering Concordance
#> 165   Cellular Manifold & Clustering Concordance
#> 166   Cellular Manifold & Clustering Concordance
#> 167   Cellular Manifold & Clustering Concordance
#> 168   Cellular Manifold & Clustering Concordance
#> 169   Cellular Manifold & Clustering Concordance
#> 170  Technical Variation & Batch Mixing Dynamics
#> 172  Technical Variation & Batch Mixing Dynamics
#> 173  Technical Variation & Batch Mixing Dynamics
#> 174  Technical Variation & Batch Mixing Dynamics
#> 175  Technical Variation & Batch Mixing Dynamics
#> 176  Technical Variation & Batch Mixing Dynamics
#> 177  Technical Variation & Batch Mixing Dynamics
#> 178  Technical Variation & Batch Mixing Dynamics
#> 179  Technical Variation & Batch Mixing Dynamics
#> 180  Technical Variation & Batch Mixing Dynamics
#> 182    Biological Signal & Differential Fidelity
#> 183    Biological Signal & Differential Fidelity
#> 184    Biological Signal & Differential Fidelity
#> 185    Biological Signal & Differential Fidelity
#> 186    Biological Signal & Differential Fidelity
#> 187    Biological Signal & Differential Fidelity
#> 188    Biological Signal & Differential Fidelity
#> 189    Biological Signal & Differential Fidelity
#> 190    Biological Signal & Differential Fidelity
#> 192    Biological Signal & Differential Fidelity
#> 193    Biological Signal & Differential Fidelity
#> 194    Biological Signal & Differential Fidelity
#> 195          (VII) Trajectory & Lineage Topology
#> 196          (VII) Trajectory & Lineage Topology
#> 197          (VII) Trajectory & Lineage Topology
#> 198          (VII) Trajectory & Lineage Topology
#> 199          (VII) Trajectory & Lineage Topology
#> 1100     (VI) Multiomics Integration & Alignment
#> 1102     (VI) Multiomics Integration & Alignment
#> 1103     (VI) Multiomics Integration & Alignment
#> 1104     (VI) Multiomics Integration & Alignment
#> 1105     (VI) Multiomics Integration & Alignment
#>                               Property                         Metric
#> 1                         library_size                    Wasserstein
#> 2                         library_size                             KS
#> 3                         library_size                            MAD
#> 4                         library_size                  Bhattacharyya
#> 5                         library_size                 Runs_Statistic
#> 6                         library_size                    NN_Mismatch
#> 7                         library_size    Between_Dataset_Silh_Global
#> 8                         library_size                            MAE
#> 9                         library_size                           RMSE
#> 10                        library_size                        KDE_ISE
#> 11                    log_library_size                    Wasserstein
#> 12                    log_library_size                             KS
#> 13                    log_library_size                            MAD
#> 14                    log_library_size                  Bhattacharyya
#> 15                    log_library_size                 Runs_Statistic
#> 16                    log_library_size                    NN_Mismatch
#> 17                    log_library_size    Between_Dataset_Silh_Global
#> 18                    log_library_size                            MAE
#> 19                    log_library_size                           RMSE
#> 20                    log_library_size                        KDE_ISE
#> 21                          tmm_factor                    Wasserstein
#> 22                          tmm_factor                             KS
#> 23                          tmm_factor                            MAD
#> 24                          tmm_factor                  Bhattacharyya
#> 25                          tmm_factor                 Runs_Statistic
#> 26                          tmm_factor                    NN_Mismatch
#> 27                          tmm_factor    Between_Dataset_Silh_Global
#> 28                          tmm_factor                            MAE
#> 29                          tmm_factor                           RMSE
#> 30                          tmm_factor                        KDE_ISE
#> 31                  zero_fraction_cell                    Wasserstein
#> 32                  zero_fraction_cell                             KS
#> 33                  zero_fraction_cell                            MAD
#> 34                  zero_fraction_cell                  Bhattacharyya
#> 35                  zero_fraction_cell                 Runs_Statistic
#> 36                  zero_fraction_cell                    NN_Mismatch
#> 37                  zero_fraction_cell    Between_Dataset_Silh_Global
#> 38                  zero_fraction_cell                            MAE
#> 39                  zero_fraction_cell                           RMSE
#> 40                  zero_fraction_cell                        KDE_ISE
#> 41                 detection_freq_cell                    Wasserstein
#> 42                 detection_freq_cell                             KS
#> 43                 detection_freq_cell                            MAD
#> 44                 detection_freq_cell                  Bhattacharyya
#> 45                 detection_freq_cell                 Runs_Statistic
#> 46                 detection_freq_cell                    NN_Mismatch
#> 47                 detection_freq_cell    Between_Dataset_Silh_Global
#> 48                 detection_freq_cell                            MAE
#> 49                 detection_freq_cell                           RMSE
#> 50                 detection_freq_cell                        KDE_ISE
#> 51                            cell_cor                    Wasserstein
#> 52                            cell_cor                             KS
#> 53                            cell_cor                            MAD
#> 54                            cell_cor                  Bhattacharyya
#> 55                            cell_cor                 Runs_Statistic
#> 56                            cell_cor                    NN_Mismatch
#> 57                            cell_cor    Between_Dataset_Silh_Global
#> 58                            cell_cor                            MAE
#> 59                            cell_cor                           RMSE
#> 60                            cell_cor                        KDE_ISE
#> 61                            pca_dist                    Wasserstein
#> 62                            pca_dist                             KS
#> 63                            pca_dist                            MAD
#> 64                            pca_dist                  Bhattacharyya
#> 65                            pca_dist                 Runs_Statistic
#> 66                            pca_dist                    NN_Mismatch
#> 67                            pca_dist    Between_Dataset_Silh_Global
#> 68                            pca_dist                            MAE
#> 69                            pca_dist                           RMSE
#> 70                            pca_dist                        KDE_ISE
#> 71                     mean_expression                    Wasserstein
#> 72                     mean_expression                             KS
#> 73                     mean_expression                            MAD
#> 74                     mean_expression                  Bhattacharyya
#> 75                     mean_expression                 Runs_Statistic
#> 76                     mean_expression                    NN_Mismatch
#> 77                     mean_expression    Between_Dataset_Silh_Global
#> 78                     mean_expression                            MAE
#> 79                     mean_expression                           RMSE
#> 80                     mean_expression                        KDE_ISE
#> 81                            variance                    Wasserstein
#> 82                            variance                             KS
#> 83                            variance                            MAD
#> 84                            variance                  Bhattacharyya
#> 85                            variance                 Runs_Statistic
#> 86                            variance                    NN_Mismatch
#> 87                            variance    Between_Dataset_Silh_Global
#> 88                            variance                            MAE
#> 89                            variance                           RMSE
#> 90                            variance                        KDE_ISE
#> 91                     scaled_variance                    Wasserstein
#> 92                     scaled_variance                             KS
#> 93                     scaled_variance                            MAD
#> 94                     scaled_variance                  Bhattacharyya
#> 95                     scaled_variance                 Runs_Statistic
#> 96                     scaled_variance                    NN_Mismatch
#> 97                     scaled_variance    Between_Dataset_Silh_Global
#> 98                     scaled_variance                            MAE
#> 99                     scaled_variance                           RMSE
#> 100                    scaled_variance                        KDE_ISE
#> 101                                 cv                    Wasserstein
#> 102                                 cv                             KS
#> 103                                 cv                            MAD
#> 104                                 cv                  Bhattacharyya
#> 105                                 cv                 Runs_Statistic
#> 106                                 cv                    NN_Mismatch
#> 107                                 cv    Between_Dataset_Silh_Global
#> 108                                 cv                            MAE
#> 109                                 cv                           RMSE
#> 110                                 cv                        KDE_ISE
#> 111              zero_fraction_feature                    Wasserstein
#> 112              zero_fraction_feature                             KS
#> 113              zero_fraction_feature                            MAD
#> 114              zero_fraction_feature                  Bhattacharyya
#> 115              zero_fraction_feature                 Runs_Statistic
#> 116              zero_fraction_feature                    NN_Mismatch
#> 117              zero_fraction_feature    Between_Dataset_Silh_Global
#> 118              zero_fraction_feature                            MAE
#> 119              zero_fraction_feature                           RMSE
#> 120              zero_fraction_feature                        KDE_ISE
#> 121                         dispersion                    Wasserstein
#> 122                         dispersion                             KS
#> 123                         dispersion                            MAD
#> 124                         dispersion                  Bhattacharyya
#> 125                         dispersion                 Runs_Statistic
#> 126                         dispersion                    NN_Mismatch
#> 127                         dispersion    Between_Dataset_Silh_Global
#> 128                         dispersion                            MAE
#> 129                         dispersion                           RMSE
#> 130                         dispersion                        KDE_ISE
#> 131                           gene_cor                    Wasserstein
#> 132                           gene_cor                             KS
#> 133                           gene_cor                            MAD
#> 134                           gene_cor                  Bhattacharyya
#> 135                           gene_cor                 Runs_Statistic
#> 136                           gene_cor                    NN_Mismatch
#> 137                           gene_cor    Between_Dataset_Silh_Global
#> 138                           gene_cor                            MAE
#> 139                           gene_cor                           RMSE
#> 140                           gene_cor                        KDE_ISE
#> 141                      dropout_curve                    slope_error
#> 142                      dropout_curve                 midpoint_error
#> 143                      dropout_curve                 r2_discrepancy
#> 144                excess_zero_weights excess_zero_weight_discrepancy
#> 145                           manifold                            MMD
#> 146                           manifold           Generative_Precision
#> 147                           manifold              Generative_Recall
#> 151                       batch_mixing               batch_silhouette
#> 152                       batch_mixing                         pcr_r2
#> 153                       batch_mixing                            cms
#> 154                       batch_mixing                            isi
#> 155                       batch_mixing                 cross_batch_F1
#> 156      scRNA_differential_expression                pval_dist_score
#> 157      scRNA_differential_expression                       accuracy
#> 158      scRNA_differential_expression                      precision
#> 159      scRNA_differential_expression                         recall
#> 160      scRNA_differential_expression                       f1_score
#> 161      scRNA_differential_expression                            auc
#> 162                       library_size                    Wasserstein
#> 210                       library_size                             KS
#> 310                       library_size                            MAD
#> 410                       library_size                  Bhattacharyya
#> 510                       library_size                 Runs_Statistic
#> 610                       library_size                    NN_Mismatch
#> 710                       library_size    Between_Dataset_Silh_Global
#> 810                       library_size                            MAE
#> 910                       library_size                           RMSE
#> 1010                      library_size                        KDE_ISE
#> 1110                  log_library_size                    Wasserstein
#> 1210                  log_library_size                             KS
#> 1310                  log_library_size                            MAD
#> 1410                  log_library_size                  Bhattacharyya
#> 1510                  log_library_size                 Runs_Statistic
#> 163                   log_library_size                    NN_Mismatch
#> 171                   log_library_size    Between_Dataset_Silh_Global
#> 181                   log_library_size                            MAE
#> 191                   log_library_size                           RMSE
#> 201                   log_library_size                        KDE_ISE
#> 211                 zero_fraction_cell                    Wasserstein
#> 221                 zero_fraction_cell                             KS
#> 231                 zero_fraction_cell                            MAD
#> 241                 zero_fraction_cell                  Bhattacharyya
#> 251                 zero_fraction_cell                 Runs_Statistic
#> 261                 zero_fraction_cell                    NN_Mismatch
#> 271                 zero_fraction_cell    Between_Dataset_Silh_Global
#> 281                 zero_fraction_cell                            MAE
#> 291                 zero_fraction_cell                           RMSE
#> 301                 zero_fraction_cell                        KDE_ISE
#> 311                detection_freq_cell                    Wasserstein
#> 321                detection_freq_cell                             KS
#> 331                detection_freq_cell                            MAD
#> 341                detection_freq_cell                  Bhattacharyya
#> 351                detection_freq_cell                 Runs_Statistic
#> 361                detection_freq_cell                    NN_Mismatch
#> 371                detection_freq_cell    Between_Dataset_Silh_Global
#> 381                detection_freq_cell                            MAE
#> 391                detection_freq_cell                           RMSE
#> 401                detection_freq_cell                        KDE_ISE
#> 411                           cell_cor                    Wasserstein
#> 421                           cell_cor                             KS
#> 431                           cell_cor                            MAD
#> 441                           cell_cor                  Bhattacharyya
#> 451                           cell_cor                 Runs_Statistic
#> 461                           cell_cor                    NN_Mismatch
#> 471                           cell_cor    Between_Dataset_Silh_Global
#> 481                           cell_cor                            MAE
#> 491                           cell_cor                           RMSE
#> 501                           cell_cor                        KDE_ISE
#> 511                           pca_dist                    Wasserstein
#> 521                           pca_dist                             KS
#> 531                           pca_dist                            MAD
#> 541                           pca_dist                  Bhattacharyya
#> 551                           pca_dist                 Runs_Statistic
#> 561                           pca_dist                    NN_Mismatch
#> 571                           pca_dist    Between_Dataset_Silh_Global
#> 581                           pca_dist                            MAE
#> 591                           pca_dist                           RMSE
#> 601                           pca_dist                        KDE_ISE
#> 611                    mean_expression                    Wasserstein
#> 621                    mean_expression                             KS
#> 631                    mean_expression                            MAD
#> 641                    mean_expression                  Bhattacharyya
#> 651                    mean_expression                 Runs_Statistic
#> 661                    mean_expression                    NN_Mismatch
#> 671                    mean_expression    Between_Dataset_Silh_Global
#> 681                    mean_expression                            MAE
#> 691                    mean_expression                           RMSE
#> 701                    mean_expression                        KDE_ISE
#> 711                           variance                    Wasserstein
#> 721                           variance                             KS
#> 731                           variance                            MAD
#> 741                           variance                  Bhattacharyya
#> 751                           variance                 Runs_Statistic
#> 761                           variance                    NN_Mismatch
#> 771                           variance    Between_Dataset_Silh_Global
#> 781                           variance                            MAE
#> 791                           variance                           RMSE
#> 801                           variance                        KDE_ISE
#> 811                    scaled_variance                    Wasserstein
#> 821                    scaled_variance                             KS
#> 831                    scaled_variance                            MAD
#> 841                    scaled_variance                  Bhattacharyya
#> 851                    scaled_variance                 Runs_Statistic
#> 861                    scaled_variance                    NN_Mismatch
#> 871                    scaled_variance    Between_Dataset_Silh_Global
#> 881                    scaled_variance                            MAE
#> 891                    scaled_variance                           RMSE
#> 901                    scaled_variance                        KDE_ISE
#> 911                                 cv                    Wasserstein
#> 921                                 cv                             KS
#> 931                                 cv                            MAD
#> 941                                 cv                  Bhattacharyya
#> 951                                 cv                 Runs_Statistic
#> 961                                 cv                    NN_Mismatch
#> 971                                 cv    Between_Dataset_Silh_Global
#> 981                                 cv                            MAE
#> 991                                 cv                           RMSE
#> 1001                                cv                        KDE_ISE
#> 1011             zero_fraction_feature                    Wasserstein
#> 1021             zero_fraction_feature                             KS
#> 1031             zero_fraction_feature                            MAD
#> 1041             zero_fraction_feature                  Bhattacharyya
#> 1051             zero_fraction_feature                 Runs_Statistic
#> 1061             zero_fraction_feature                    NN_Mismatch
#> 1071             zero_fraction_feature    Between_Dataset_Silh_Global
#> 1081             zero_fraction_feature                            MAE
#> 1091             zero_fraction_feature                           RMSE
#> 1101             zero_fraction_feature                        KDE_ISE
#> 1111                        dispersion                    Wasserstein
#> 1121                        dispersion                             KS
#> 1131                        dispersion                            MAD
#> 1141                        dispersion                  Bhattacharyya
#> 1151                        dispersion                 Runs_Statistic
#> 1161                        dispersion                    NN_Mismatch
#> 1171                        dispersion    Between_Dataset_Silh_Global
#> 1181                        dispersion                            MAE
#> 1191                        dispersion                           RMSE
#> 1201                        dispersion                        KDE_ISE
#> 1211                          gene_cor                    Wasserstein
#> 1221                          gene_cor                             KS
#> 1231                          gene_cor                            MAD
#> 1241                          gene_cor                  Bhattacharyya
#> 1251                          gene_cor                 Runs_Statistic
#> 1261                          gene_cor                    NN_Mismatch
#> 1271                          gene_cor    Between_Dataset_Silh_Global
#> 1281                          gene_cor                            MAE
#> 1291                          gene_cor                           RMSE
#> 1301                          gene_cor                        KDE_ISE
#> 1311      accessibility_sparsity_curve                     curve_rmse
#> 1321      accessibility_sparsity_curve                 r2_discrepancy
#> 1331               excess_zero_weights excess_zero_weight_discrepancy
#> 1341                          manifold                            MMD
#> 1351                          manifold           Generative_Precision
#> 1361                          manifold              Generative_Recall
#> 1401                      batch_mixing               batch_silhouette
#> 1411                      batch_mixing                         pcr_r2
#> 1421                      batch_mixing                            cms
#> 1431                      batch_mixing                            isi
#> 1441                      batch_mixing                 cross_batch_F1
#> 1451 scATAC_differential_accessibility                pval_dist_score
#> 1461 scATAC_differential_accessibility                       accuracy
#> 1471 scATAC_differential_accessibility                      precision
#> 148  scATAC_differential_accessibility                         recall
#> 149  scATAC_differential_accessibility                       f1_score
#> 150  scATAC_differential_accessibility                            auc
#> 164                         clustering                 silhouette_sim
#> 165                         clustering                 silhouette_sim
#> 166                         clustering                            cdi
#> 167                         clustering                            cdi
#> 168                         clustering            neighborhood_purity
#> 169                         clustering            neighborhood_purity
#> 170                       batch_mixing               batch_silhouette
#> 172                       batch_mixing                         pcr_r2
#> 173                       batch_mixing                            cms
#> 174                       batch_mixing                            isi
#> 175                       batch_mixing                 cross_batch_F1
#> 176                       batch_mixing               batch_silhouette
#> 177                       batch_mixing                         pcr_r2
#> 178                       batch_mixing                            cms
#> 179                       batch_mixing                            isi
#> 180                       batch_mixing                 cross_batch_F1
#> 182      scRNA_differential_expression                pval_dist_score
#> 183      scRNA_differential_expression                       accuracy
#> 184      scRNA_differential_expression                      precision
#> 185      scRNA_differential_expression                         recall
#> 186      scRNA_differential_expression                       f1_score
#> 187      scRNA_differential_expression                            auc
#> 188  scATAC_differential_accessibility                pval_dist_score
#> 189  scATAC_differential_accessibility                       accuracy
#> 190  scATAC_differential_accessibility                      precision
#> 192  scATAC_differential_accessibility                         recall
#> 193  scATAC_differential_accessibility                       f1_score
#> 194  scATAC_differential_accessibility                            auc
#> 195                trajectory_topology                            him
#> 196                trajectory_topology                      edge_flip
#> 197                  geodesic_ordering           geodesic_correlation
#> 198                 branch_concordance                    F1_branches
#> 199                   dynamic_features                  wcor_features
#> 1100            multiomics_integration                      batch_asw
#> 1102            multiomics_integration                  cell_type_asw
#> 1103            multiomics_integration             graph_connectivity
#> 1104            multiomics_integration                            ari
#> 1105            multiomics_integration                            nmi
#>              Value            Modality
#> 1     1.506250e+01                 rna
#> 2     2.125000e-01                 rna
#> 3     1.500000e+01                 rna
#> 4     1.188839e-04                 rna
#> 5    -1.744748e+00                 rna
#> 6     1.750000e-01                 rna
#> 7     2.247078e-02                 rna
#> 8     1.506250e+01                 rna
#> 9     1.571504e+01                 rna
#> 10    7.980370e-04                 rna
#> 11    6.194965e-02                 rna
#> 12    2.125000e-01                 rna
#> 13    6.363671e-02                 rna
#> 14    3.834591e-06                 rna
#> 15   -1.744748e+00                 rna
#> 16    1.875000e-01                 rna
#> 17    2.176201e-02                 rna
#> 18    6.194965e-02                 rna
#> 19    6.566950e-02                 rna
#> 20    1.962832e-01                 rna
#> 21    3.268387e-02                 rna
#> 22    1.018380e-01                 rna
#> 23    2.555864e-02                 rna
#> 24    2.186586e-04                 rna
#> 25   -9.516807e-01                 rna
#> 26    5.000000e-02                 rna
#> 27   -6.195226e-04                 rna
#> 28    3.268387e-02                 rna
#> 29    4.613293e-02                 rna
#> 30    7.601569e-02                 rna
#> 31    1.562500e-02                 rna
#> 32    1.250000e-01                 rna
#> 33    1.666667e-02                 rna
#> 34    8.854773e-04                 rna
#> 35   -9.358194e+00                 rna
#> 36    6.875000e-01                 rna
#> 37    7.524680e-03                 rna
#> 38    1.562500e-02                 rna
#> 39    1.980811e-02                 rna
#> 40    3.253205e-01                 rna
#> 41    1.562500e-02                 rna
#> 42    1.250000e-01                 rna
#> 43    1.666667e-02                 rna
#> 44    2.679272e-05                 rna
#> 45   -9.040967e+00                 rna
#> 46    6.875000e-01                 rna
#> 47    7.524680e-03                 rna
#> 48    1.562500e-02                 rna
#> 49    1.980811e-02                 rna
#> 50    3.253205e-01                 rna
#> 51    5.728994e-03                 rna
#> 52    3.323906e-02                 rna
#> 53    5.684572e-03                 rna
#> 54    2.566053e-05                 rna
#> 55   -2.188896e+00                 rna
#> 56    1.666667e-03                 rna
#> 57    1.423443e-03                 rna
#> 58    5.728994e-03                 rna
#> 59    6.328019e-03                 rna
#> 60    6.120211e-03                 rna
#> 61    1.594340e+00                 rna
#> 62    8.747786e-02                 rna
#> 63    1.428096e+00                 rna
#> 64    3.105370e-05                 rna
#> 65    2.515972e-02                 rna
#> 66    1.500000e-02                 rna
#> 67    9.907338e-03                 rna
#> 68    1.594340e+00                 rna
#> 69    1.676818e+00                 rna
#> 70    8.585380e-04                 rna
#> 71    2.473408e-01                 rna
#> 72    2.055826e-01                 rna
#> 73    2.571903e-01                 rna
#> 74    3.690266e-05                 rna
#> 75    3.666923e-01                 rna
#> 76    5.000000e-02                 rna
#> 77   -9.070018e-05                 rna
#> 78    2.473408e-01                 rna
#> 79    2.887192e-01                 rna
#> 80    1.208325e-02                 rna
#> 81    2.036911e+00                 rna
#> 82    1.533879e-01                 rna
#> 83    1.969816e+00                 rna
#> 84    3.196235e-04                 rna
#> 85   -3.666923e-01                 rna
#> 86    4.166667e-02                 rna
#> 87    9.752197e-04                 rna
#> 88    2.036911e+00                 rna
#> 89    2.271370e+00                 rna
#> 90    1.659081e-03                 rna
#> 91    1.119261e-01                 rna
#> 92    6.742514e-02                 rna
#> 93    7.339222e-02                 rna
#> 94    2.010495e-03                 rna
#> 95   -1.283423e+00                 rna
#> 96    1.666667e-02                 rna
#> 97   -1.234822e-02                 rna
#> 98    1.119261e-01                 rna
#> 99    1.504922e-01                 rna
#> 100   6.265612e-03                 rna
#> 101   7.432972e+00                 rna
#> 102   2.906651e-01                 rna
#> 103   7.484040e+00                 rna
#> 104   7.165290e-05                 rna
#> 105  -1.100077e+00                 rna
#> 106   4.166667e-02                 rna
#> 107   5.142707e-02                 rna
#> 108   7.432972e+00                 rna
#> 109   7.844679e+00                 rna
#> 110   4.198168e-03                 rna
#> 111   1.687500e-02                 rna
#> 112   1.333333e-01                 rna
#> 113   1.250000e-02                 rna
#> 114   6.100311e-04                 rna
#> 115  -4.400308e+00                 rna
#> 116   9.166667e-02                 rna
#> 117  -1.887003e-03                 rna
#> 118   1.687500e-02                 rna
#> 119   1.916214e-02                 rna
#> 120   1.732278e-01                 rna
#> 121   4.456931e-01                 rna
#> 122   1.741354e-01                 rna
#> 123   4.388661e-01                 rna
#> 124   5.586790e-04                 rna
#> 125   1.466769e+00                 rna
#> 126   2.500000e-02                 rna
#> 127   1.024133e-02                 rna
#> 128   4.456931e-01                 rna
#> 129   4.890173e-01                 rna
#> 130   1.955467e-02                 rna
#> 131   3.299875e-03                 rna
#> 132   2.672058e-02                 rna
#> 133   2.118671e-03                 rna
#> 134   6.113770e-05                 rna
#> 135  -2.050782e+00                 rna
#> 136   3.000000e-02                 rna
#> 137  -2.482463e-06                 rna
#> 138   3.299875e-03                 rna
#> 139   5.067781e-03                 rna
#> 140   1.403010e-02                 rna
#> 141   2.560909e-01                 rna
#> 142   2.624144e-01                 rna
#> 143   1.646825e-01                 rna
#> 144   3.712431e-02                 rna
#> 145   1.103894e-01                 rna
#> 146   7.625000e-01                 rna
#> 147   9.625000e-01                 rna
#> 151  -2.710457e-03                 rna
#> 152   9.476487e-03                 rna
#> 153   0.000000e+00                 rna
#> 154   1.938397e+00                 rna
#> 155   5.126418e-01                 rna
#> 156   8.719073e-01                 rna
#> 157   9.500000e-01                 rna
#> 158   0.000000e+00                 rna
#> 159   0.000000e+00                 rna
#> 160   0.000000e+00                 rna
#> 161   6.465517e-01                 rna
#> 162   2.337500e+00                atac
#> 210   1.750000e-01                atac
#> 310   2.000000e+00                atac
#> 410   1.055641e-04                atac
#> 510  -6.027311e+00                atac
#> 610   1.937500e-01                atac
#> 710   2.188673e-02                atac
#> 810   2.337500e+00                atac
#> 910   2.487469e+00                atac
#> 1010  4.478217e-03                atac
#> 1110  7.991730e-02                atac
#> 1210  1.750000e-01                atac
#> 1310  8.223810e-02                atac
#> 1410  1.155313e-05                atac
#> 1510 -6.027311e+00                atac
#> 163   1.937500e-01                atac
#> 171   1.990590e-02                atac
#> 181   7.991730e-02                atac
#> 191   8.389317e-02                atac
#> 201   1.268716e-01                atac
#> 211   2.750000e-02                atac
#> 221   2.125000e-01                atac
#> 231   3.333333e-02                atac
#> 241   6.907067e-05                atac
#> 251  -7.296219e+00                atac
#> 261   2.812500e-01                atac
#> 271   2.566981e-02                atac
#> 281   2.750000e-02                atac
#> 291   3.016160e-02                atac
#> 301   4.176638e-01                atac
#> 311   2.750000e-02                atac
#> 321   2.125000e-01                atac
#> 331   3.333333e-02                atac
#> 341   1.795082e-04                atac
#> 351  -6.978992e+00                atac
#> 361   2.812500e-01                atac
#> 371   2.566981e-02                atac
#> 381   2.750000e-02                atac
#> 391   3.016160e-02                atac
#> 401   4.176638e-01                atac
#> 411   3.342856e-03                atac
#> 421   2.056815e-02                atac
#> 431   2.270508e-03                atac
#> 441   5.400110e-05                atac
#> 451  -2.515972e-01                atac
#> 461   4.333333e-02                atac
#> 471  -2.274271e-03                atac
#> 481   3.342856e-03                atac
#> 491   6.619754e-03                atac
#> 501   4.263335e-03                atac
#> 511   1.502554e-01                atac
#> 521   2.948321e-02                atac
#> 531   1.327661e-01                atac
#> 541   2.878871e-06                atac
#> 551   7.547916e-01                atac
#> 561   2.666667e-02                atac
#> 571   6.595987e-04                atac
#> 581   1.502554e-01                atac
#> 591   1.883752e-01                atac
#> 601   3.459236e-04                atac
#> 611   1.501367e-01                atac
#> 621   1.879472e-01                atac
#> 631   1.536693e-01                atac
#> 641   2.959894e-04                atac
#> 651  -1.833462e+00                atac
#> 661   2.500000e-02                atac
#> 671   3.291150e-03                atac
#> 681   1.501367e-01                atac
#> 691   1.677762e-01                atac
#> 701   2.809910e-02                atac
#> 711   2.170317e-01                atac
#> 721   2.753225e-01                atac
#> 731   2.419452e-01                atac
#> 741   5.389967e-05                atac
#> 751  -2.566846e+00                atac
#> 761   1.166667e-01                atac
#> 771   3.168098e-02                atac
#> 781   2.170317e-01                atac
#> 791   2.541541e-01                atac
#> 801   1.429083e-01                atac
#> 811   1.993187e-01                atac
#> 821   1.174387e-01                atac
#> 831   1.419895e-01                atac
#> 841   2.917823e-03                atac
#> 851  -9.167309e-01                atac
#> 861   8.333333e-02                atac
#> 871  -3.112124e-03                atac
#> 881   1.993187e-01                atac
#> 891   2.711405e-01                atac
#> 901   3.669152e-02                atac
#> 911   8.864353e+00                atac
#> 921   1.763639e-01                atac
#> 931   8.752987e+00                atac
#> 941   2.185211e-04                atac
#> 951  -1.833462e+00                atac
#> 961   3.333333e-02                atac
#> 971   6.571696e-03                atac
#> 981   8.864353e+00                atac
#> 991   1.046915e+01                atac
#> 1001  7.968311e-04                atac
#> 1011  2.958333e-02                atac
#> 1021  2.000000e-01                atac
#> 1031  2.500000e-02                atac
#> 1041  1.230865e-04                atac
#> 1051 -1.650116e+00                atac
#> 1061  5.000000e-02                atac
#> 1071  7.445303e-03                atac
#> 1081  2.958333e-02                atac
#> 1091  3.307189e-02                atac
#> 1101  1.835132e-01                atac
#> 1111  2.544810e-02                atac
#> 1121  1.136364e-01                atac
#> 1131  2.119535e-02                atac
#> 1141  1.193835e-04                atac
#> 1151 -5.500385e-01                atac
#> 1161  8.333333e-03                atac
#> 1171 -7.584339e-04                atac
#> 1181  2.544810e-02                atac
#> 1191  2.978845e-02                atac
#> 1201  1.190016e-01                atac
#> 1211  2.718776e-03                atac
#> 1221  1.759557e-02                atac
#> 1231  1.838118e-03                atac
#> 1241  8.878612e-05                atac
#> 1251 -6.723876e-02                atac
#> 1261  4.666667e-02                atac
#> 1271 -2.426294e-04                atac
#> 1281  2.718776e-03                atac
#> 1291  5.260476e-03                atac
#> 1301  6.132365e-03                atac
#> 1311  1.233884e-02                atac
#> 1321  8.535249e-04                atac
#> 1331  8.480348e-04                atac
#> 1341  1.096098e-01                atac
#> 1351  1.000000e+00                atac
#> 1361  9.500000e-01                atac
#> 1401  8.780031e-04                atac
#> 1411  1.375709e-02                atac
#> 1421  0.000000e+00                atac
#> 1431  1.931492e+00                atac
#> 1441  7.346700e-01                atac
#> 1451  8.563379e-01                atac
#> 1461  8.166667e-01                atac
#> 1471  5.000000e-01                atac
#> 148   9.090909e-02                atac
#> 149   1.538462e-01                atac
#> 150   8.367347e-01                atac
#> 164   3.402312e-03                 rna
#> 165   4.149183e-02                atac
#> 166   4.984432e+00                 rna
#> 167   1.687605e+00                atac
#> 168   4.775000e-01                 rna
#> 169   7.275000e-01                atac
#> 170  -2.710457e-03                 rna
#> 172   9.476487e-03                 rna
#> 173   0.000000e+00                 rna
#> 174   1.938397e+00                 rna
#> 175   5.126418e-01                 rna
#> 176   8.780031e-04                atac
#> 177   1.375709e-02                atac
#> 178   0.000000e+00                atac
#> 179   1.931492e+00                atac
#> 180   7.346700e-01                atac
#> 182   8.719073e-01                 rna
#> 183   9.500000e-01                 rna
#> 184   0.000000e+00                 rna
#> 185   0.000000e+00                 rna
#> 186   0.000000e+00                 rna
#> 187   6.465517e-01                 rna
#> 188   8.563379e-01                atac
#> 189   8.166667e-01                atac
#> 190   5.000000e-01                atac
#> 192   9.090909e-02                atac
#> 193   1.538462e-01                atac
#> 194   8.367347e-01                atac
#> 195   1.000000e+00                 rna
#> 196   1.000000e+00                 rna
#> 197   4.848647e-01                 rna
#> 198   1.000000e+00                 rna
#> 199   3.802351e-01                 rna
#> 1100  9.747094e-01 Joint (Cross-Modal)
#> 1102  5.582026e-01 Joint (Cross-Modal)
#> 1103  1.000000e+00 Joint (Cross-Modal)
#> 1104  6.765538e-01 Joint (Cross-Modal)
#> 1105  5.821401e-01 Joint (Cross-Modal)
#> 
#> $summary_table
#>                                         Category
#> 1                      Distributional Properties
#> 2                      Distributional Properties
#> 3                      Distributional Properties
#> 4                      Distributional Properties
#> 5                      Distributional Properties
#> 6                      Distributional Properties
#> 7                      Distributional Properties
#> 8                      Distributional Properties
#> 9                      Distributional Properties
#> 10                     Distributional Properties
#> 11                     Distributional Properties
#> 12                     Distributional Properties
#> 13                     Distributional Properties
#> 14                     Distributional Properties
#> 15                     Distributional Properties
#> 16                     Distributional Properties
#> 17                     Distributional Properties
#> 18                     Distributional Properties
#> 19                     Distributional Properties
#> 20                     Distributional Properties
#> 21                     Distributional Properties
#> 22                     Distributional Properties
#> 23                     Distributional Properties
#> 24                     Distributional Properties
#> 25                     Distributional Properties
#> 26                     Distributional Properties
#> 27                     Distributional Properties
#> 28                     Distributional Properties
#> 29                     Distributional Properties
#> 30                     Distributional Properties
#> 31                     Distributional Properties
#> 32                     Distributional Properties
#> 33                     Distributional Properties
#> 34                     Distributional Properties
#> 35                     Distributional Properties
#> 36                     Distributional Properties
#> 37                     Distributional Properties
#> 38                     Distributional Properties
#> 39                     Distributional Properties
#> 40                     Distributional Properties
#> 41                     Distributional Properties
#> 42                     Distributional Properties
#> 43                     Distributional Properties
#> 44                     Distributional Properties
#> 45                     Distributional Properties
#> 46                     Distributional Properties
#> 47                     Distributional Properties
#> 48                     Distributional Properties
#> 49                     Distributional Properties
#> 50                     Distributional Properties
#> 51                     Distributional Properties
#> 52                     Distributional Properties
#> 53                     Distributional Properties
#> 54                     Distributional Properties
#> 55                     Distributional Properties
#> 56                     Distributional Properties
#> 57                     Distributional Properties
#> 58                     Distributional Properties
#> 59                     Distributional Properties
#> 60                     Distributional Properties
#> 61                     Distributional Properties
#> 62                     Distributional Properties
#> 63                     Distributional Properties
#> 64                     Distributional Properties
#> 65                     Distributional Properties
#> 66                     Distributional Properties
#> 67                     Distributional Properties
#> 68                     Distributional Properties
#> 69                     Distributional Properties
#> 70                     Distributional Properties
#> 71                     Distributional Properties
#> 72                     Distributional Properties
#> 73                     Distributional Properties
#> 74                     Distributional Properties
#> 75                     Distributional Properties
#> 76                     Distributional Properties
#> 77                     Distributional Properties
#> 78                     Distributional Properties
#> 79                     Distributional Properties
#> 80                     Distributional Properties
#> 81                     Distributional Properties
#> 82                     Distributional Properties
#> 83                     Distributional Properties
#> 84                     Distributional Properties
#> 85                     Distributional Properties
#> 86                     Distributional Properties
#> 87                     Distributional Properties
#> 88                     Distributional Properties
#> 89                     Distributional Properties
#> 90                     Distributional Properties
#> 91                     Distributional Properties
#> 92                     Distributional Properties
#> 93                     Distributional Properties
#> 94                     Distributional Properties
#> 95                     Distributional Properties
#> 96                     Distributional Properties
#> 97                     Distributional Properties
#> 98                     Distributional Properties
#> 99                     Distributional Properties
#> 100                    Distributional Properties
#> 101                    Distributional Properties
#> 102                    Distributional Properties
#> 103                    Distributional Properties
#> 104                    Distributional Properties
#> 105                    Distributional Properties
#> 106                    Distributional Properties
#> 107                    Distributional Properties
#> 108                    Distributional Properties
#> 109                    Distributional Properties
#> 110                    Distributional Properties
#> 111                    Distributional Properties
#> 112                    Distributional Properties
#> 113                    Distributional Properties
#> 114                    Distributional Properties
#> 115                    Distributional Properties
#> 116                    Distributional Properties
#> 117                    Distributional Properties
#> 118                    Distributional Properties
#> 119                    Distributional Properties
#> 120                    Distributional Properties
#> 121                    Distributional Properties
#> 122                    Distributional Properties
#> 123                    Distributional Properties
#> 124                    Distributional Properties
#> 125                    Distributional Properties
#> 126                    Distributional Properties
#> 127                    Distributional Properties
#> 128                    Distributional Properties
#> 129                    Distributional Properties
#> 130                    Distributional Properties
#> 131                    Distributional Properties
#> 132                    Distributional Properties
#> 133                    Distributional Properties
#> 134                    Distributional Properties
#> 135                    Distributional Properties
#> 136                    Distributional Properties
#> 137                    Distributional Properties
#> 138                    Distributional Properties
#> 139                    Distributional Properties
#> 140                    Distributional Properties
#> 141           Sparsity & Zero-Inflation Dynamics
#> 142           Sparsity & Zero-Inflation Dynamics
#> 143           Sparsity & Zero-Inflation Dynamics
#> 144           Sparsity & Zero-Inflation Dynamics
#> 145   Cellular Manifold & Clustering Concordance
#> 146   Cellular Manifold & Clustering Concordance
#> 147   Cellular Manifold & Clustering Concordance
#> 151  Technical Variation & Batch Mixing Dynamics
#> 152  Technical Variation & Batch Mixing Dynamics
#> 153  Technical Variation & Batch Mixing Dynamics
#> 154  Technical Variation & Batch Mixing Dynamics
#> 155  Technical Variation & Batch Mixing Dynamics
#> 156    Biological Signal & Differential Fidelity
#> 157    Biological Signal & Differential Fidelity
#> 158    Biological Signal & Differential Fidelity
#> 159    Biological Signal & Differential Fidelity
#> 160    Biological Signal & Differential Fidelity
#> 161    Biological Signal & Differential Fidelity
#> 162                    Distributional Properties
#> 210                    Distributional Properties
#> 310                    Distributional Properties
#> 410                    Distributional Properties
#> 510                    Distributional Properties
#> 610                    Distributional Properties
#> 710                    Distributional Properties
#> 810                    Distributional Properties
#> 910                    Distributional Properties
#> 1010                   Distributional Properties
#> 1110                   Distributional Properties
#> 1210                   Distributional Properties
#> 1310                   Distributional Properties
#> 1410                   Distributional Properties
#> 1510                   Distributional Properties
#> 163                    Distributional Properties
#> 171                    Distributional Properties
#> 181                    Distributional Properties
#> 191                    Distributional Properties
#> 201                    Distributional Properties
#> 211                    Distributional Properties
#> 221                    Distributional Properties
#> 231                    Distributional Properties
#> 241                    Distributional Properties
#> 251                    Distributional Properties
#> 261                    Distributional Properties
#> 271                    Distributional Properties
#> 281                    Distributional Properties
#> 291                    Distributional Properties
#> 301                    Distributional Properties
#> 311                    Distributional Properties
#> 321                    Distributional Properties
#> 331                    Distributional Properties
#> 341                    Distributional Properties
#> 351                    Distributional Properties
#> 361                    Distributional Properties
#> 371                    Distributional Properties
#> 381                    Distributional Properties
#> 391                    Distributional Properties
#> 401                    Distributional Properties
#> 411                    Distributional Properties
#> 421                    Distributional Properties
#> 431                    Distributional Properties
#> 441                    Distributional Properties
#> 451                    Distributional Properties
#> 461                    Distributional Properties
#> 471                    Distributional Properties
#> 481                    Distributional Properties
#> 491                    Distributional Properties
#> 501                    Distributional Properties
#> 511                    Distributional Properties
#> 521                    Distributional Properties
#> 531                    Distributional Properties
#> 541                    Distributional Properties
#> 551                    Distributional Properties
#> 561                    Distributional Properties
#> 571                    Distributional Properties
#> 581                    Distributional Properties
#> 591                    Distributional Properties
#> 601                    Distributional Properties
#> 611                    Distributional Properties
#> 621                    Distributional Properties
#> 631                    Distributional Properties
#> 641                    Distributional Properties
#> 651                    Distributional Properties
#> 661                    Distributional Properties
#> 671                    Distributional Properties
#> 681                    Distributional Properties
#> 691                    Distributional Properties
#> 701                    Distributional Properties
#> 711                    Distributional Properties
#> 721                    Distributional Properties
#> 731                    Distributional Properties
#> 741                    Distributional Properties
#> 751                    Distributional Properties
#> 761                    Distributional Properties
#> 771                    Distributional Properties
#> 781                    Distributional Properties
#> 791                    Distributional Properties
#> 801                    Distributional Properties
#> 811                    Distributional Properties
#> 821                    Distributional Properties
#> 831                    Distributional Properties
#> 841                    Distributional Properties
#> 851                    Distributional Properties
#> 861                    Distributional Properties
#> 871                    Distributional Properties
#> 881                    Distributional Properties
#> 891                    Distributional Properties
#> 901                    Distributional Properties
#> 911                    Distributional Properties
#> 921                    Distributional Properties
#> 931                    Distributional Properties
#> 941                    Distributional Properties
#> 951                    Distributional Properties
#> 961                    Distributional Properties
#> 971                    Distributional Properties
#> 981                    Distributional Properties
#> 991                    Distributional Properties
#> 1001                   Distributional Properties
#> 1011                   Distributional Properties
#> 1021                   Distributional Properties
#> 1031                   Distributional Properties
#> 1041                   Distributional Properties
#> 1051                   Distributional Properties
#> 1061                   Distributional Properties
#> 1071                   Distributional Properties
#> 1081                   Distributional Properties
#> 1091                   Distributional Properties
#> 1101                   Distributional Properties
#> 1111                   Distributional Properties
#> 1121                   Distributional Properties
#> 1131                   Distributional Properties
#> 1141                   Distributional Properties
#> 1151                   Distributional Properties
#> 1161                   Distributional Properties
#> 1171                   Distributional Properties
#> 1181                   Distributional Properties
#> 1191                   Distributional Properties
#> 1201                   Distributional Properties
#> 1211                   Distributional Properties
#> 1221                   Distributional Properties
#> 1231                   Distributional Properties
#> 1241                   Distributional Properties
#> 1251                   Distributional Properties
#> 1261                   Distributional Properties
#> 1271                   Distributional Properties
#> 1281                   Distributional Properties
#> 1291                   Distributional Properties
#> 1301                   Distributional Properties
#> 1311          Sparsity & Zero-Inflation Dynamics
#> 1321          Sparsity & Zero-Inflation Dynamics
#> 1331          Sparsity & Zero-Inflation Dynamics
#> 1341  Cellular Manifold & Clustering Concordance
#> 1351  Cellular Manifold & Clustering Concordance
#> 1361  Cellular Manifold & Clustering Concordance
#> 1401 Technical Variation & Batch Mixing Dynamics
#> 1411 Technical Variation & Batch Mixing Dynamics
#> 1421 Technical Variation & Batch Mixing Dynamics
#> 1431 Technical Variation & Batch Mixing Dynamics
#> 1441 Technical Variation & Batch Mixing Dynamics
#> 1451   Biological Signal & Differential Fidelity
#> 1461   Biological Signal & Differential Fidelity
#> 1471   Biological Signal & Differential Fidelity
#> 148    Biological Signal & Differential Fidelity
#> 149    Biological Signal & Differential Fidelity
#> 150    Biological Signal & Differential Fidelity
#> 164   Cellular Manifold & Clustering Concordance
#> 165   Cellular Manifold & Clustering Concordance
#> 166   Cellular Manifold & Clustering Concordance
#> 167   Cellular Manifold & Clustering Concordance
#> 168   Cellular Manifold & Clustering Concordance
#> 169   Cellular Manifold & Clustering Concordance
#> 170  Technical Variation & Batch Mixing Dynamics
#> 172  Technical Variation & Batch Mixing Dynamics
#> 173  Technical Variation & Batch Mixing Dynamics
#> 174  Technical Variation & Batch Mixing Dynamics
#> 175  Technical Variation & Batch Mixing Dynamics
#> 176  Technical Variation & Batch Mixing Dynamics
#> 177  Technical Variation & Batch Mixing Dynamics
#> 178  Technical Variation & Batch Mixing Dynamics
#> 179  Technical Variation & Batch Mixing Dynamics
#> 180  Technical Variation & Batch Mixing Dynamics
#> 182    Biological Signal & Differential Fidelity
#> 183    Biological Signal & Differential Fidelity
#> 184    Biological Signal & Differential Fidelity
#> 185    Biological Signal & Differential Fidelity
#> 186    Biological Signal & Differential Fidelity
#> 187    Biological Signal & Differential Fidelity
#> 188    Biological Signal & Differential Fidelity
#> 189    Biological Signal & Differential Fidelity
#> 190    Biological Signal & Differential Fidelity
#> 192    Biological Signal & Differential Fidelity
#> 193    Biological Signal & Differential Fidelity
#> 194    Biological Signal & Differential Fidelity
#> 195          (VII) Trajectory & Lineage Topology
#> 196          (VII) Trajectory & Lineage Topology
#> 197          (VII) Trajectory & Lineage Topology
#> 198          (VII) Trajectory & Lineage Topology
#> 199          (VII) Trajectory & Lineage Topology
#> 1100     (VI) Multiomics Integration & Alignment
#> 1102     (VI) Multiomics Integration & Alignment
#> 1103     (VI) Multiomics Integration & Alignment
#> 1104     (VI) Multiomics Integration & Alignment
#> 1105     (VI) Multiomics Integration & Alignment
#>                               Property                         Metric
#> 1                         library_size                    Wasserstein
#> 2                         library_size                             KS
#> 3                         library_size                            MAD
#> 4                         library_size                  Bhattacharyya
#> 5                         library_size                 Runs_Statistic
#> 6                         library_size                    NN_Mismatch
#> 7                         library_size    Between_Dataset_Silh_Global
#> 8                         library_size                            MAE
#> 9                         library_size                           RMSE
#> 10                        library_size                        KDE_ISE
#> 11                    log_library_size                    Wasserstein
#> 12                    log_library_size                             KS
#> 13                    log_library_size                            MAD
#> 14                    log_library_size                  Bhattacharyya
#> 15                    log_library_size                 Runs_Statistic
#> 16                    log_library_size                    NN_Mismatch
#> 17                    log_library_size    Between_Dataset_Silh_Global
#> 18                    log_library_size                            MAE
#> 19                    log_library_size                           RMSE
#> 20                    log_library_size                        KDE_ISE
#> 21                          tmm_factor                    Wasserstein
#> 22                          tmm_factor                             KS
#> 23                          tmm_factor                            MAD
#> 24                          tmm_factor                  Bhattacharyya
#> 25                          tmm_factor                 Runs_Statistic
#> 26                          tmm_factor                    NN_Mismatch
#> 27                          tmm_factor    Between_Dataset_Silh_Global
#> 28                          tmm_factor                            MAE
#> 29                          tmm_factor                           RMSE
#> 30                          tmm_factor                        KDE_ISE
#> 31                  zero_fraction_cell                    Wasserstein
#> 32                  zero_fraction_cell                             KS
#> 33                  zero_fraction_cell                            MAD
#> 34                  zero_fraction_cell                  Bhattacharyya
#> 35                  zero_fraction_cell                 Runs_Statistic
#> 36                  zero_fraction_cell                    NN_Mismatch
#> 37                  zero_fraction_cell    Between_Dataset_Silh_Global
#> 38                  zero_fraction_cell                            MAE
#> 39                  zero_fraction_cell                           RMSE
#> 40                  zero_fraction_cell                        KDE_ISE
#> 41                 detection_freq_cell                    Wasserstein
#> 42                 detection_freq_cell                             KS
#> 43                 detection_freq_cell                            MAD
#> 44                 detection_freq_cell                  Bhattacharyya
#> 45                 detection_freq_cell                 Runs_Statistic
#> 46                 detection_freq_cell                    NN_Mismatch
#> 47                 detection_freq_cell    Between_Dataset_Silh_Global
#> 48                 detection_freq_cell                            MAE
#> 49                 detection_freq_cell                           RMSE
#> 50                 detection_freq_cell                        KDE_ISE
#> 51                            cell_cor                    Wasserstein
#> 52                            cell_cor                             KS
#> 53                            cell_cor                            MAD
#> 54                            cell_cor                  Bhattacharyya
#> 55                            cell_cor                 Runs_Statistic
#> 56                            cell_cor                    NN_Mismatch
#> 57                            cell_cor    Between_Dataset_Silh_Global
#> 58                            cell_cor                            MAE
#> 59                            cell_cor                           RMSE
#> 60                            cell_cor                        KDE_ISE
#> 61                            pca_dist                    Wasserstein
#> 62                            pca_dist                             KS
#> 63                            pca_dist                            MAD
#> 64                            pca_dist                  Bhattacharyya
#> 65                            pca_dist                 Runs_Statistic
#> 66                            pca_dist                    NN_Mismatch
#> 67                            pca_dist    Between_Dataset_Silh_Global
#> 68                            pca_dist                            MAE
#> 69                            pca_dist                           RMSE
#> 70                            pca_dist                        KDE_ISE
#> 71                     mean_expression                    Wasserstein
#> 72                     mean_expression                             KS
#> 73                     mean_expression                            MAD
#> 74                     mean_expression                  Bhattacharyya
#> 75                     mean_expression                 Runs_Statistic
#> 76                     mean_expression                    NN_Mismatch
#> 77                     mean_expression    Between_Dataset_Silh_Global
#> 78                     mean_expression                            MAE
#> 79                     mean_expression                           RMSE
#> 80                     mean_expression                        KDE_ISE
#> 81                            variance                    Wasserstein
#> 82                            variance                             KS
#> 83                            variance                            MAD
#> 84                            variance                  Bhattacharyya
#> 85                            variance                 Runs_Statistic
#> 86                            variance                    NN_Mismatch
#> 87                            variance    Between_Dataset_Silh_Global
#> 88                            variance                            MAE
#> 89                            variance                           RMSE
#> 90                            variance                        KDE_ISE
#> 91                     scaled_variance                    Wasserstein
#> 92                     scaled_variance                             KS
#> 93                     scaled_variance                            MAD
#> 94                     scaled_variance                  Bhattacharyya
#> 95                     scaled_variance                 Runs_Statistic
#> 96                     scaled_variance                    NN_Mismatch
#> 97                     scaled_variance    Between_Dataset_Silh_Global
#> 98                     scaled_variance                            MAE
#> 99                     scaled_variance                           RMSE
#> 100                    scaled_variance                        KDE_ISE
#> 101                                 cv                    Wasserstein
#> 102                                 cv                             KS
#> 103                                 cv                            MAD
#> 104                                 cv                  Bhattacharyya
#> 105                                 cv                 Runs_Statistic
#> 106                                 cv                    NN_Mismatch
#> 107                                 cv    Between_Dataset_Silh_Global
#> 108                                 cv                            MAE
#> 109                                 cv                           RMSE
#> 110                                 cv                        KDE_ISE
#> 111              zero_fraction_feature                    Wasserstein
#> 112              zero_fraction_feature                             KS
#> 113              zero_fraction_feature                            MAD
#> 114              zero_fraction_feature                  Bhattacharyya
#> 115              zero_fraction_feature                 Runs_Statistic
#> 116              zero_fraction_feature                    NN_Mismatch
#> 117              zero_fraction_feature    Between_Dataset_Silh_Global
#> 118              zero_fraction_feature                            MAE
#> 119              zero_fraction_feature                           RMSE
#> 120              zero_fraction_feature                        KDE_ISE
#> 121                         dispersion                    Wasserstein
#> 122                         dispersion                             KS
#> 123                         dispersion                            MAD
#> 124                         dispersion                  Bhattacharyya
#> 125                         dispersion                 Runs_Statistic
#> 126                         dispersion                    NN_Mismatch
#> 127                         dispersion    Between_Dataset_Silh_Global
#> 128                         dispersion                            MAE
#> 129                         dispersion                           RMSE
#> 130                         dispersion                        KDE_ISE
#> 131                           gene_cor                    Wasserstein
#> 132                           gene_cor                             KS
#> 133                           gene_cor                            MAD
#> 134                           gene_cor                  Bhattacharyya
#> 135                           gene_cor                 Runs_Statistic
#> 136                           gene_cor                    NN_Mismatch
#> 137                           gene_cor    Between_Dataset_Silh_Global
#> 138                           gene_cor                            MAE
#> 139                           gene_cor                           RMSE
#> 140                           gene_cor                        KDE_ISE
#> 141                      dropout_curve                    slope_error
#> 142                      dropout_curve                 midpoint_error
#> 143                      dropout_curve                 r2_discrepancy
#> 144                excess_zero_weights excess_zero_weight_discrepancy
#> 145                           manifold                            MMD
#> 146                           manifold           Generative_Precision
#> 147                           manifold              Generative_Recall
#> 151                       batch_mixing               batch_silhouette
#> 152                       batch_mixing                         pcr_r2
#> 153                       batch_mixing                            cms
#> 154                       batch_mixing                            isi
#> 155                       batch_mixing                 cross_batch_F1
#> 156      scRNA_differential_expression                pval_dist_score
#> 157      scRNA_differential_expression                       accuracy
#> 158      scRNA_differential_expression                      precision
#> 159      scRNA_differential_expression                         recall
#> 160      scRNA_differential_expression                       f1_score
#> 161      scRNA_differential_expression                            auc
#> 162                       library_size                    Wasserstein
#> 210                       library_size                             KS
#> 310                       library_size                            MAD
#> 410                       library_size                  Bhattacharyya
#> 510                       library_size                 Runs_Statistic
#> 610                       library_size                    NN_Mismatch
#> 710                       library_size    Between_Dataset_Silh_Global
#> 810                       library_size                            MAE
#> 910                       library_size                           RMSE
#> 1010                      library_size                        KDE_ISE
#> 1110                  log_library_size                    Wasserstein
#> 1210                  log_library_size                             KS
#> 1310                  log_library_size                            MAD
#> 1410                  log_library_size                  Bhattacharyya
#> 1510                  log_library_size                 Runs_Statistic
#> 163                   log_library_size                    NN_Mismatch
#> 171                   log_library_size    Between_Dataset_Silh_Global
#> 181                   log_library_size                            MAE
#> 191                   log_library_size                           RMSE
#> 201                   log_library_size                        KDE_ISE
#> 211                 zero_fraction_cell                    Wasserstein
#> 221                 zero_fraction_cell                             KS
#> 231                 zero_fraction_cell                            MAD
#> 241                 zero_fraction_cell                  Bhattacharyya
#> 251                 zero_fraction_cell                 Runs_Statistic
#> 261                 zero_fraction_cell                    NN_Mismatch
#> 271                 zero_fraction_cell    Between_Dataset_Silh_Global
#> 281                 zero_fraction_cell                            MAE
#> 291                 zero_fraction_cell                           RMSE
#> 301                 zero_fraction_cell                        KDE_ISE
#> 311                detection_freq_cell                    Wasserstein
#> 321                detection_freq_cell                             KS
#> 331                detection_freq_cell                            MAD
#> 341                detection_freq_cell                  Bhattacharyya
#> 351                detection_freq_cell                 Runs_Statistic
#> 361                detection_freq_cell                    NN_Mismatch
#> 371                detection_freq_cell    Between_Dataset_Silh_Global
#> 381                detection_freq_cell                            MAE
#> 391                detection_freq_cell                           RMSE
#> 401                detection_freq_cell                        KDE_ISE
#> 411                           cell_cor                    Wasserstein
#> 421                           cell_cor                             KS
#> 431                           cell_cor                            MAD
#> 441                           cell_cor                  Bhattacharyya
#> 451                           cell_cor                 Runs_Statistic
#> 461                           cell_cor                    NN_Mismatch
#> 471                           cell_cor    Between_Dataset_Silh_Global
#> 481                           cell_cor                            MAE
#> 491                           cell_cor                           RMSE
#> 501                           cell_cor                        KDE_ISE
#> 511                           pca_dist                    Wasserstein
#> 521                           pca_dist                             KS
#> 531                           pca_dist                            MAD
#> 541                           pca_dist                  Bhattacharyya
#> 551                           pca_dist                 Runs_Statistic
#> 561                           pca_dist                    NN_Mismatch
#> 571                           pca_dist    Between_Dataset_Silh_Global
#> 581                           pca_dist                            MAE
#> 591                           pca_dist                           RMSE
#> 601                           pca_dist                        KDE_ISE
#> 611                    mean_expression                    Wasserstein
#> 621                    mean_expression                             KS
#> 631                    mean_expression                            MAD
#> 641                    mean_expression                  Bhattacharyya
#> 651                    mean_expression                 Runs_Statistic
#> 661                    mean_expression                    NN_Mismatch
#> 671                    mean_expression    Between_Dataset_Silh_Global
#> 681                    mean_expression                            MAE
#> 691                    mean_expression                           RMSE
#> 701                    mean_expression                        KDE_ISE
#> 711                           variance                    Wasserstein
#> 721                           variance                             KS
#> 731                           variance                            MAD
#> 741                           variance                  Bhattacharyya
#> 751                           variance                 Runs_Statistic
#> 761                           variance                    NN_Mismatch
#> 771                           variance    Between_Dataset_Silh_Global
#> 781                           variance                            MAE
#> 791                           variance                           RMSE
#> 801                           variance                        KDE_ISE
#> 811                    scaled_variance                    Wasserstein
#> 821                    scaled_variance                             KS
#> 831                    scaled_variance                            MAD
#> 841                    scaled_variance                  Bhattacharyya
#> 851                    scaled_variance                 Runs_Statistic
#> 861                    scaled_variance                    NN_Mismatch
#> 871                    scaled_variance    Between_Dataset_Silh_Global
#> 881                    scaled_variance                            MAE
#> 891                    scaled_variance                           RMSE
#> 901                    scaled_variance                        KDE_ISE
#> 911                                 cv                    Wasserstein
#> 921                                 cv                             KS
#> 931                                 cv                            MAD
#> 941                                 cv                  Bhattacharyya
#> 951                                 cv                 Runs_Statistic
#> 961                                 cv                    NN_Mismatch
#> 971                                 cv    Between_Dataset_Silh_Global
#> 981                                 cv                            MAE
#> 991                                 cv                           RMSE
#> 1001                                cv                        KDE_ISE
#> 1011             zero_fraction_feature                    Wasserstein
#> 1021             zero_fraction_feature                             KS
#> 1031             zero_fraction_feature                            MAD
#> 1041             zero_fraction_feature                  Bhattacharyya
#> 1051             zero_fraction_feature                 Runs_Statistic
#> 1061             zero_fraction_feature                    NN_Mismatch
#> 1071             zero_fraction_feature    Between_Dataset_Silh_Global
#> 1081             zero_fraction_feature                            MAE
#> 1091             zero_fraction_feature                           RMSE
#> 1101             zero_fraction_feature                        KDE_ISE
#> 1111                        dispersion                    Wasserstein
#> 1121                        dispersion                             KS
#> 1131                        dispersion                            MAD
#> 1141                        dispersion                  Bhattacharyya
#> 1151                        dispersion                 Runs_Statistic
#> 1161                        dispersion                    NN_Mismatch
#> 1171                        dispersion    Between_Dataset_Silh_Global
#> 1181                        dispersion                            MAE
#> 1191                        dispersion                           RMSE
#> 1201                        dispersion                        KDE_ISE
#> 1211                          gene_cor                    Wasserstein
#> 1221                          gene_cor                             KS
#> 1231                          gene_cor                            MAD
#> 1241                          gene_cor                  Bhattacharyya
#> 1251                          gene_cor                 Runs_Statistic
#> 1261                          gene_cor                    NN_Mismatch
#> 1271                          gene_cor    Between_Dataset_Silh_Global
#> 1281                          gene_cor                            MAE
#> 1291                          gene_cor                           RMSE
#> 1301                          gene_cor                        KDE_ISE
#> 1311      accessibility_sparsity_curve                     curve_rmse
#> 1321      accessibility_sparsity_curve                 r2_discrepancy
#> 1331               excess_zero_weights excess_zero_weight_discrepancy
#> 1341                          manifold                            MMD
#> 1351                          manifold           Generative_Precision
#> 1361                          manifold              Generative_Recall
#> 1401                      batch_mixing               batch_silhouette
#> 1411                      batch_mixing                         pcr_r2
#> 1421                      batch_mixing                            cms
#> 1431                      batch_mixing                            isi
#> 1441                      batch_mixing                 cross_batch_F1
#> 1451 scATAC_differential_accessibility                pval_dist_score
#> 1461 scATAC_differential_accessibility                       accuracy
#> 1471 scATAC_differential_accessibility                      precision
#> 148  scATAC_differential_accessibility                         recall
#> 149  scATAC_differential_accessibility                       f1_score
#> 150  scATAC_differential_accessibility                            auc
#> 164                         clustering                 silhouette_sim
#> 165                         clustering                 silhouette_sim
#> 166                         clustering                            cdi
#> 167                         clustering                            cdi
#> 168                         clustering            neighborhood_purity
#> 169                         clustering            neighborhood_purity
#> 170                       batch_mixing               batch_silhouette
#> 172                       batch_mixing                         pcr_r2
#> 173                       batch_mixing                            cms
#> 174                       batch_mixing                            isi
#> 175                       batch_mixing                 cross_batch_F1
#> 176                       batch_mixing               batch_silhouette
#> 177                       batch_mixing                         pcr_r2
#> 178                       batch_mixing                            cms
#> 179                       batch_mixing                            isi
#> 180                       batch_mixing                 cross_batch_F1
#> 182      scRNA_differential_expression                pval_dist_score
#> 183      scRNA_differential_expression                       accuracy
#> 184      scRNA_differential_expression                      precision
#> 185      scRNA_differential_expression                         recall
#> 186      scRNA_differential_expression                       f1_score
#> 187      scRNA_differential_expression                            auc
#> 188  scATAC_differential_accessibility                pval_dist_score
#> 189  scATAC_differential_accessibility                       accuracy
#> 190  scATAC_differential_accessibility                      precision
#> 192  scATAC_differential_accessibility                         recall
#> 193  scATAC_differential_accessibility                       f1_score
#> 194  scATAC_differential_accessibility                            auc
#> 195                trajectory_topology                            him
#> 196                trajectory_topology                      edge_flip
#> 197                  geodesic_ordering           geodesic_correlation
#> 198                 branch_concordance                    F1_branches
#> 199                   dynamic_features                  wcor_features
#> 1100            multiomics_integration                      batch_asw
#> 1102            multiomics_integration                  cell_type_asw
#> 1103            multiomics_integration             graph_connectivity
#> 1104            multiomics_integration                            ari
#> 1105            multiomics_integration                            nmi
#>              Value            Modality
#> 1     1.506250e+01                 rna
#> 2     2.125000e-01                 rna
#> 3     1.500000e+01                 rna
#> 4     1.188839e-04                 rna
#> 5    -1.744748e+00                 rna
#> 6     1.750000e-01                 rna
#> 7     2.247078e-02                 rna
#> 8     1.506250e+01                 rna
#> 9     1.571504e+01                 rna
#> 10    7.980370e-04                 rna
#> 11    6.194965e-02                 rna
#> 12    2.125000e-01                 rna
#> 13    6.363671e-02                 rna
#> 14    3.834591e-06                 rna
#> 15   -1.744748e+00                 rna
#> 16    1.875000e-01                 rna
#> 17    2.176201e-02                 rna
#> 18    6.194965e-02                 rna
#> 19    6.566950e-02                 rna
#> 20    1.962832e-01                 rna
#> 21    3.268387e-02                 rna
#> 22    1.018380e-01                 rna
#> 23    2.555864e-02                 rna
#> 24    2.186586e-04                 rna
#> 25   -9.516807e-01                 rna
#> 26    5.000000e-02                 rna
#> 27   -6.195226e-04                 rna
#> 28    3.268387e-02                 rna
#> 29    4.613293e-02                 rna
#> 30    7.601569e-02                 rna
#> 31    1.562500e-02                 rna
#> 32    1.250000e-01                 rna
#> 33    1.666667e-02                 rna
#> 34    8.854773e-04                 rna
#> 35   -9.358194e+00                 rna
#> 36    6.875000e-01                 rna
#> 37    7.524680e-03                 rna
#> 38    1.562500e-02                 rna
#> 39    1.980811e-02                 rna
#> 40    3.253205e-01                 rna
#> 41    1.562500e-02                 rna
#> 42    1.250000e-01                 rna
#> 43    1.666667e-02                 rna
#> 44    2.679272e-05                 rna
#> 45   -9.040967e+00                 rna
#> 46    6.875000e-01                 rna
#> 47    7.524680e-03                 rna
#> 48    1.562500e-02                 rna
#> 49    1.980811e-02                 rna
#> 50    3.253205e-01                 rna
#> 51    5.728994e-03                 rna
#> 52    3.323906e-02                 rna
#> 53    5.684572e-03                 rna
#> 54    2.566053e-05                 rna
#> 55   -2.188896e+00                 rna
#> 56    1.666667e-03                 rna
#> 57    1.423443e-03                 rna
#> 58    5.728994e-03                 rna
#> 59    6.328019e-03                 rna
#> 60    6.120211e-03                 rna
#> 61    1.594340e+00                 rna
#> 62    8.747786e-02                 rna
#> 63    1.428096e+00                 rna
#> 64    3.105370e-05                 rna
#> 65    2.515972e-02                 rna
#> 66    1.500000e-02                 rna
#> 67    9.907338e-03                 rna
#> 68    1.594340e+00                 rna
#> 69    1.676818e+00                 rna
#> 70    8.585380e-04                 rna
#> 71    2.473408e-01                 rna
#> 72    2.055826e-01                 rna
#> 73    2.571903e-01                 rna
#> 74    3.690266e-05                 rna
#> 75    3.666923e-01                 rna
#> 76    5.000000e-02                 rna
#> 77   -9.070018e-05                 rna
#> 78    2.473408e-01                 rna
#> 79    2.887192e-01                 rna
#> 80    1.208325e-02                 rna
#> 81    2.036911e+00                 rna
#> 82    1.533879e-01                 rna
#> 83    1.969816e+00                 rna
#> 84    3.196235e-04                 rna
#> 85   -3.666923e-01                 rna
#> 86    4.166667e-02                 rna
#> 87    9.752197e-04                 rna
#> 88    2.036911e+00                 rna
#> 89    2.271370e+00                 rna
#> 90    1.659081e-03                 rna
#> 91    1.119261e-01                 rna
#> 92    6.742514e-02                 rna
#> 93    7.339222e-02                 rna
#> 94    2.010495e-03                 rna
#> 95   -1.283423e+00                 rna
#> 96    1.666667e-02                 rna
#> 97   -1.234822e-02                 rna
#> 98    1.119261e-01                 rna
#> 99    1.504922e-01                 rna
#> 100   6.265612e-03                 rna
#> 101   7.432972e+00                 rna
#> 102   2.906651e-01                 rna
#> 103   7.484040e+00                 rna
#> 104   7.165290e-05                 rna
#> 105  -1.100077e+00                 rna
#> 106   4.166667e-02                 rna
#> 107   5.142707e-02                 rna
#> 108   7.432972e+00                 rna
#> 109   7.844679e+00                 rna
#> 110   4.198168e-03                 rna
#> 111   1.687500e-02                 rna
#> 112   1.333333e-01                 rna
#> 113   1.250000e-02                 rna
#> 114   6.100311e-04                 rna
#> 115  -4.400308e+00                 rna
#> 116   9.166667e-02                 rna
#> 117  -1.887003e-03                 rna
#> 118   1.687500e-02                 rna
#> 119   1.916214e-02                 rna
#> 120   1.732278e-01                 rna
#> 121   4.456931e-01                 rna
#> 122   1.741354e-01                 rna
#> 123   4.388661e-01                 rna
#> 124   5.586790e-04                 rna
#> 125   1.466769e+00                 rna
#> 126   2.500000e-02                 rna
#> 127   1.024133e-02                 rna
#> 128   4.456931e-01                 rna
#> 129   4.890173e-01                 rna
#> 130   1.955467e-02                 rna
#> 131   3.299875e-03                 rna
#> 132   2.672058e-02                 rna
#> 133   2.118671e-03                 rna
#> 134   6.113770e-05                 rna
#> 135  -2.050782e+00                 rna
#> 136   3.000000e-02                 rna
#> 137  -2.482463e-06                 rna
#> 138   3.299875e-03                 rna
#> 139   5.067781e-03                 rna
#> 140   1.403010e-02                 rna
#> 141   2.560909e-01                 rna
#> 142   2.624144e-01                 rna
#> 143   1.646825e-01                 rna
#> 144   3.712431e-02                 rna
#> 145   1.103894e-01                 rna
#> 146   7.625000e-01                 rna
#> 147   9.625000e-01                 rna
#> 151  -2.710457e-03                 rna
#> 152   9.476487e-03                 rna
#> 153   0.000000e+00                 rna
#> 154   1.938397e+00                 rna
#> 155   5.126418e-01                 rna
#> 156   8.719073e-01                 rna
#> 157   9.500000e-01                 rna
#> 158   0.000000e+00                 rna
#> 159   0.000000e+00                 rna
#> 160   0.000000e+00                 rna
#> 161   6.465517e-01                 rna
#> 162   2.337500e+00                atac
#> 210   1.750000e-01                atac
#> 310   2.000000e+00                atac
#> 410   1.055641e-04                atac
#> 510  -6.027311e+00                atac
#> 610   1.937500e-01                atac
#> 710   2.188673e-02                atac
#> 810   2.337500e+00                atac
#> 910   2.487469e+00                atac
#> 1010  4.478217e-03                atac
#> 1110  7.991730e-02                atac
#> 1210  1.750000e-01                atac
#> 1310  8.223810e-02                atac
#> 1410  1.155313e-05                atac
#> 1510 -6.027311e+00                atac
#> 163   1.937500e-01                atac
#> 171   1.990590e-02                atac
#> 181   7.991730e-02                atac
#> 191   8.389317e-02                atac
#> 201   1.268716e-01                atac
#> 211   2.750000e-02                atac
#> 221   2.125000e-01                atac
#> 231   3.333333e-02                atac
#> 241   6.907067e-05                atac
#> 251  -7.296219e+00                atac
#> 261   2.812500e-01                atac
#> 271   2.566981e-02                atac
#> 281   2.750000e-02                atac
#> 291   3.016160e-02                atac
#> 301   4.176638e-01                atac
#> 311   2.750000e-02                atac
#> 321   2.125000e-01                atac
#> 331   3.333333e-02                atac
#> 341   1.795082e-04                atac
#> 351  -6.978992e+00                atac
#> 361   2.812500e-01                atac
#> 371   2.566981e-02                atac
#> 381   2.750000e-02                atac
#> 391   3.016160e-02                atac
#> 401   4.176638e-01                atac
#> 411   3.342856e-03                atac
#> 421   2.056815e-02                atac
#> 431   2.270508e-03                atac
#> 441   5.400110e-05                atac
#> 451  -2.515972e-01                atac
#> 461   4.333333e-02                atac
#> 471  -2.274271e-03                atac
#> 481   3.342856e-03                atac
#> 491   6.619754e-03                atac
#> 501   4.263335e-03                atac
#> 511   1.502554e-01                atac
#> 521   2.948321e-02                atac
#> 531   1.327661e-01                atac
#> 541   2.878871e-06                atac
#> 551   7.547916e-01                atac
#> 561   2.666667e-02                atac
#> 571   6.595987e-04                atac
#> 581   1.502554e-01                atac
#> 591   1.883752e-01                atac
#> 601   3.459236e-04                atac
#> 611   1.501367e-01                atac
#> 621   1.879472e-01                atac
#> 631   1.536693e-01                atac
#> 641   2.959894e-04                atac
#> 651  -1.833462e+00                atac
#> 661   2.500000e-02                atac
#> 671   3.291150e-03                atac
#> 681   1.501367e-01                atac
#> 691   1.677762e-01                atac
#> 701   2.809910e-02                atac
#> 711   2.170317e-01                atac
#> 721   2.753225e-01                atac
#> 731   2.419452e-01                atac
#> 741   5.389967e-05                atac
#> 751  -2.566846e+00                atac
#> 761   1.166667e-01                atac
#> 771   3.168098e-02                atac
#> 781   2.170317e-01                atac
#> 791   2.541541e-01                atac
#> 801   1.429083e-01                atac
#> 811   1.993187e-01                atac
#> 821   1.174387e-01                atac
#> 831   1.419895e-01                atac
#> 841   2.917823e-03                atac
#> 851  -9.167309e-01                atac
#> 861   8.333333e-02                atac
#> 871  -3.112124e-03                atac
#> 881   1.993187e-01                atac
#> 891   2.711405e-01                atac
#> 901   3.669152e-02                atac
#> 911   8.864353e+00                atac
#> 921   1.763639e-01                atac
#> 931   8.752987e+00                atac
#> 941   2.185211e-04                atac
#> 951  -1.833462e+00                atac
#> 961   3.333333e-02                atac
#> 971   6.571696e-03                atac
#> 981   8.864353e+00                atac
#> 991   1.046915e+01                atac
#> 1001  7.968311e-04                atac
#> 1011  2.958333e-02                atac
#> 1021  2.000000e-01                atac
#> 1031  2.500000e-02                atac
#> 1041  1.230865e-04                atac
#> 1051 -1.650116e+00                atac
#> 1061  5.000000e-02                atac
#> 1071  7.445303e-03                atac
#> 1081  2.958333e-02                atac
#> 1091  3.307189e-02                atac
#> 1101  1.835132e-01                atac
#> 1111  2.544810e-02                atac
#> 1121  1.136364e-01                atac
#> 1131  2.119535e-02                atac
#> 1141  1.193835e-04                atac
#> 1151 -5.500385e-01                atac
#> 1161  8.333333e-03                atac
#> 1171 -7.584339e-04                atac
#> 1181  2.544810e-02                atac
#> 1191  2.978845e-02                atac
#> 1201  1.190016e-01                atac
#> 1211  2.718776e-03                atac
#> 1221  1.759557e-02                atac
#> 1231  1.838118e-03                atac
#> 1241  8.878612e-05                atac
#> 1251 -6.723876e-02                atac
#> 1261  4.666667e-02                atac
#> 1271 -2.426294e-04                atac
#> 1281  2.718776e-03                atac
#> 1291  5.260476e-03                atac
#> 1301  6.132365e-03                atac
#> 1311  1.233884e-02                atac
#> 1321  8.535249e-04                atac
#> 1331  8.480348e-04                atac
#> 1341  1.096098e-01                atac
#> 1351  1.000000e+00                atac
#> 1361  9.500000e-01                atac
#> 1401  8.780031e-04                atac
#> 1411  1.375709e-02                atac
#> 1421  0.000000e+00                atac
#> 1431  1.931492e+00                atac
#> 1441  7.346700e-01                atac
#> 1451  8.563379e-01                atac
#> 1461  8.166667e-01                atac
#> 1471  5.000000e-01                atac
#> 148   9.090909e-02                atac
#> 149   1.538462e-01                atac
#> 150   8.367347e-01                atac
#> 164   3.402312e-03                 rna
#> 165   4.149183e-02                atac
#> 166   4.984432e+00                 rna
#> 167   1.687605e+00                atac
#> 168   4.775000e-01                 rna
#> 169   7.275000e-01                atac
#> 170  -2.710457e-03                 rna
#> 172   9.476487e-03                 rna
#> 173   0.000000e+00                 rna
#> 174   1.938397e+00                 rna
#> 175   5.126418e-01                 rna
#> 176   8.780031e-04                atac
#> 177   1.375709e-02                atac
#> 178   0.000000e+00                atac
#> 179   1.931492e+00                atac
#> 180   7.346700e-01                atac
#> 182   8.719073e-01                 rna
#> 183   9.500000e-01                 rna
#> 184   0.000000e+00                 rna
#> 185   0.000000e+00                 rna
#> 186   0.000000e+00                 rna
#> 187   6.465517e-01                 rna
#> 188   8.563379e-01                atac
#> 189   8.166667e-01                atac
#> 190   5.000000e-01                atac
#> 192   9.090909e-02                atac
#> 193   1.538462e-01                atac
#> 194   8.367347e-01                atac
#> 195   1.000000e+00                 rna
#> 196   1.000000e+00                 rna
#> 197   4.848647e-01                 rna
#> 198   1.000000e+00                 rna
#> 199   3.802351e-01                 rna
#> 1100  9.747094e-01 Joint (Cross-Modal)
#> 1102  5.582026e-01 Joint (Cross-Modal)
#> 1103  1.000000e+00 Joint (Cross-Modal)
#> 1104  6.765538e-01 Joint (Cross-Modal)
#> 1105  5.821401e-01 Joint (Cross-Modal)
#> 
#> $pairing
#> [1] "paired"
#> 
#> $mod1_unimodal
#> $mod1_unimodal$metrics_summary_table
#>                                        Category                      Property
#> 1                     Distributional Properties                  library_size
#> 2                     Distributional Properties                  library_size
#> 3                     Distributional Properties                  library_size
#> 4                     Distributional Properties                  library_size
#> 5                     Distributional Properties                  library_size
#> 6                     Distributional Properties                  library_size
#> 7                     Distributional Properties                  library_size
#> 8                     Distributional Properties                  library_size
#> 9                     Distributional Properties                  library_size
#> 10                    Distributional Properties                  library_size
#> 11                    Distributional Properties              log_library_size
#> 12                    Distributional Properties              log_library_size
#> 13                    Distributional Properties              log_library_size
#> 14                    Distributional Properties              log_library_size
#> 15                    Distributional Properties              log_library_size
#> 16                    Distributional Properties              log_library_size
#> 17                    Distributional Properties              log_library_size
#> 18                    Distributional Properties              log_library_size
#> 19                    Distributional Properties              log_library_size
#> 20                    Distributional Properties              log_library_size
#> 21                    Distributional Properties                    tmm_factor
#> 22                    Distributional Properties                    tmm_factor
#> 23                    Distributional Properties                    tmm_factor
#> 24                    Distributional Properties                    tmm_factor
#> 25                    Distributional Properties                    tmm_factor
#> 26                    Distributional Properties                    tmm_factor
#> 27                    Distributional Properties                    tmm_factor
#> 28                    Distributional Properties                    tmm_factor
#> 29                    Distributional Properties                    tmm_factor
#> 30                    Distributional Properties                    tmm_factor
#> 31                    Distributional Properties            zero_fraction_cell
#> 32                    Distributional Properties            zero_fraction_cell
#> 33                    Distributional Properties            zero_fraction_cell
#> 34                    Distributional Properties            zero_fraction_cell
#> 35                    Distributional Properties            zero_fraction_cell
#> 36                    Distributional Properties            zero_fraction_cell
#> 37                    Distributional Properties            zero_fraction_cell
#> 38                    Distributional Properties            zero_fraction_cell
#> 39                    Distributional Properties            zero_fraction_cell
#> 40                    Distributional Properties            zero_fraction_cell
#> 41                    Distributional Properties           detection_freq_cell
#> 42                    Distributional Properties           detection_freq_cell
#> 43                    Distributional Properties           detection_freq_cell
#> 44                    Distributional Properties           detection_freq_cell
#> 45                    Distributional Properties           detection_freq_cell
#> 46                    Distributional Properties           detection_freq_cell
#> 47                    Distributional Properties           detection_freq_cell
#> 48                    Distributional Properties           detection_freq_cell
#> 49                    Distributional Properties           detection_freq_cell
#> 50                    Distributional Properties           detection_freq_cell
#> 51                    Distributional Properties                      cell_cor
#> 52                    Distributional Properties                      cell_cor
#> 53                    Distributional Properties                      cell_cor
#> 54                    Distributional Properties                      cell_cor
#> 55                    Distributional Properties                      cell_cor
#> 56                    Distributional Properties                      cell_cor
#> 57                    Distributional Properties                      cell_cor
#> 58                    Distributional Properties                      cell_cor
#> 59                    Distributional Properties                      cell_cor
#> 60                    Distributional Properties                      cell_cor
#> 61                    Distributional Properties                      pca_dist
#> 62                    Distributional Properties                      pca_dist
#> 63                    Distributional Properties                      pca_dist
#> 64                    Distributional Properties                      pca_dist
#> 65                    Distributional Properties                      pca_dist
#> 66                    Distributional Properties                      pca_dist
#> 67                    Distributional Properties                      pca_dist
#> 68                    Distributional Properties                      pca_dist
#> 69                    Distributional Properties                      pca_dist
#> 70                    Distributional Properties                      pca_dist
#> 71                    Distributional Properties               mean_expression
#> 72                    Distributional Properties               mean_expression
#> 73                    Distributional Properties               mean_expression
#> 74                    Distributional Properties               mean_expression
#> 75                    Distributional Properties               mean_expression
#> 76                    Distributional Properties               mean_expression
#> 77                    Distributional Properties               mean_expression
#> 78                    Distributional Properties               mean_expression
#> 79                    Distributional Properties               mean_expression
#> 80                    Distributional Properties               mean_expression
#> 81                    Distributional Properties                      variance
#> 82                    Distributional Properties                      variance
#> 83                    Distributional Properties                      variance
#> 84                    Distributional Properties                      variance
#> 85                    Distributional Properties                      variance
#> 86                    Distributional Properties                      variance
#> 87                    Distributional Properties                      variance
#> 88                    Distributional Properties                      variance
#> 89                    Distributional Properties                      variance
#> 90                    Distributional Properties                      variance
#> 91                    Distributional Properties               scaled_variance
#> 92                    Distributional Properties               scaled_variance
#> 93                    Distributional Properties               scaled_variance
#> 94                    Distributional Properties               scaled_variance
#> 95                    Distributional Properties               scaled_variance
#> 96                    Distributional Properties               scaled_variance
#> 97                    Distributional Properties               scaled_variance
#> 98                    Distributional Properties               scaled_variance
#> 99                    Distributional Properties               scaled_variance
#> 100                   Distributional Properties               scaled_variance
#> 101                   Distributional Properties                            cv
#> 102                   Distributional Properties                            cv
#> 103                   Distributional Properties                            cv
#> 104                   Distributional Properties                            cv
#> 105                   Distributional Properties                            cv
#> 106                   Distributional Properties                            cv
#> 107                   Distributional Properties                            cv
#> 108                   Distributional Properties                            cv
#> 109                   Distributional Properties                            cv
#> 110                   Distributional Properties                            cv
#> 111                   Distributional Properties         zero_fraction_feature
#> 112                   Distributional Properties         zero_fraction_feature
#> 113                   Distributional Properties         zero_fraction_feature
#> 114                   Distributional Properties         zero_fraction_feature
#> 115                   Distributional Properties         zero_fraction_feature
#> 116                   Distributional Properties         zero_fraction_feature
#> 117                   Distributional Properties         zero_fraction_feature
#> 118                   Distributional Properties         zero_fraction_feature
#> 119                   Distributional Properties         zero_fraction_feature
#> 120                   Distributional Properties         zero_fraction_feature
#> 121                   Distributional Properties                    dispersion
#> 122                   Distributional Properties                    dispersion
#> 123                   Distributional Properties                    dispersion
#> 124                   Distributional Properties                    dispersion
#> 125                   Distributional Properties                    dispersion
#> 126                   Distributional Properties                    dispersion
#> 127                   Distributional Properties                    dispersion
#> 128                   Distributional Properties                    dispersion
#> 129                   Distributional Properties                    dispersion
#> 130                   Distributional Properties                    dispersion
#> 131                   Distributional Properties                      gene_cor
#> 132                   Distributional Properties                      gene_cor
#> 133                   Distributional Properties                      gene_cor
#> 134                   Distributional Properties                      gene_cor
#> 135                   Distributional Properties                      gene_cor
#> 136                   Distributional Properties                      gene_cor
#> 137                   Distributional Properties                      gene_cor
#> 138                   Distributional Properties                      gene_cor
#> 139                   Distributional Properties                      gene_cor
#> 140                   Distributional Properties                      gene_cor
#> 141          Sparsity & Zero-Inflation Dynamics                 dropout_curve
#> 142          Sparsity & Zero-Inflation Dynamics                 dropout_curve
#> 143          Sparsity & Zero-Inflation Dynamics                 dropout_curve
#> 144          Sparsity & Zero-Inflation Dynamics           excess_zero_weights
#> 145  Cellular Manifold & Clustering Concordance                      manifold
#> 146  Cellular Manifold & Clustering Concordance                      manifold
#> 147  Cellular Manifold & Clustering Concordance                      manifold
#> 148  Cellular Manifold & Clustering Concordance                    clustering
#> 149  Cellular Manifold & Clustering Concordance                    clustering
#> 150  Cellular Manifold & Clustering Concordance                    clustering
#> 151 Technical Variation & Batch Mixing Dynamics                  batch_mixing
#> 152 Technical Variation & Batch Mixing Dynamics                  batch_mixing
#> 153 Technical Variation & Batch Mixing Dynamics                  batch_mixing
#> 154 Technical Variation & Batch Mixing Dynamics                  batch_mixing
#> 155 Technical Variation & Batch Mixing Dynamics                  batch_mixing
#> 156   Biological Signal & Differential Fidelity scRNA_differential_expression
#> 157   Biological Signal & Differential Fidelity scRNA_differential_expression
#> 158   Biological Signal & Differential Fidelity scRNA_differential_expression
#> 159   Biological Signal & Differential Fidelity scRNA_differential_expression
#> 160   Biological Signal & Differential Fidelity scRNA_differential_expression
#> 161   Biological Signal & Differential Fidelity scRNA_differential_expression
#>                             Metric         Value
#> 1                      Wasserstein  1.506250e+01
#> 2                               KS  2.125000e-01
#> 3                              MAD  1.500000e+01
#> 4                    Bhattacharyya  1.188839e-04
#> 5                   Runs_Statistic -1.744748e+00
#> 6                      NN_Mismatch  1.750000e-01
#> 7      Between_Dataset_Silh_Global  2.247078e-02
#> 8                              MAE  1.506250e+01
#> 9                             RMSE  1.571504e+01
#> 10                         KDE_ISE  7.980370e-04
#> 11                     Wasserstein  6.194965e-02
#> 12                              KS  2.125000e-01
#> 13                             MAD  6.363671e-02
#> 14                   Bhattacharyya  3.834591e-06
#> 15                  Runs_Statistic -1.744748e+00
#> 16                     NN_Mismatch  1.875000e-01
#> 17     Between_Dataset_Silh_Global  2.176201e-02
#> 18                             MAE  6.194965e-02
#> 19                            RMSE  6.566950e-02
#> 20                         KDE_ISE  1.962832e-01
#> 21                     Wasserstein  3.268387e-02
#> 22                              KS  1.018380e-01
#> 23                             MAD  2.555864e-02
#> 24                   Bhattacharyya  2.186586e-04
#> 25                  Runs_Statistic -9.516807e-01
#> 26                     NN_Mismatch  5.000000e-02
#> 27     Between_Dataset_Silh_Global -6.195226e-04
#> 28                             MAE  3.268387e-02
#> 29                            RMSE  4.613293e-02
#> 30                         KDE_ISE  7.601569e-02
#> 31                     Wasserstein  1.562500e-02
#> 32                              KS  1.250000e-01
#> 33                             MAD  1.666667e-02
#> 34                   Bhattacharyya  8.854773e-04
#> 35                  Runs_Statistic -9.358194e+00
#> 36                     NN_Mismatch  6.875000e-01
#> 37     Between_Dataset_Silh_Global  7.524680e-03
#> 38                             MAE  1.562500e-02
#> 39                            RMSE  1.980811e-02
#> 40                         KDE_ISE  3.253205e-01
#> 41                     Wasserstein  1.562500e-02
#> 42                              KS  1.250000e-01
#> 43                             MAD  1.666667e-02
#> 44                   Bhattacharyya  2.679272e-05
#> 45                  Runs_Statistic -9.040967e+00
#> 46                     NN_Mismatch  6.875000e-01
#> 47     Between_Dataset_Silh_Global  7.524680e-03
#> 48                             MAE  1.562500e-02
#> 49                            RMSE  1.980811e-02
#> 50                         KDE_ISE  3.253205e-01
#> 51                     Wasserstein  5.728994e-03
#> 52                              KS  3.323906e-02
#> 53                             MAD  5.684572e-03
#> 54                   Bhattacharyya  2.566053e-05
#> 55                  Runs_Statistic -2.188896e+00
#> 56                     NN_Mismatch  1.666667e-03
#> 57     Between_Dataset_Silh_Global  1.423443e-03
#> 58                             MAE  5.728994e-03
#> 59                            RMSE  6.328019e-03
#> 60                         KDE_ISE  6.120211e-03
#> 61                     Wasserstein  1.594340e+00
#> 62                              KS  8.747786e-02
#> 63                             MAD  1.428096e+00
#> 64                   Bhattacharyya  3.105370e-05
#> 65                  Runs_Statistic  2.515972e-02
#> 66                     NN_Mismatch  1.500000e-02
#> 67     Between_Dataset_Silh_Global  9.907338e-03
#> 68                             MAE  1.594340e+00
#> 69                            RMSE  1.676818e+00
#> 70                         KDE_ISE  8.585380e-04
#> 71                     Wasserstein  2.473408e-01
#> 72                              KS  2.055826e-01
#> 73                             MAD  2.571903e-01
#> 74                   Bhattacharyya  3.690266e-05
#> 75                  Runs_Statistic  3.666923e-01
#> 76                     NN_Mismatch  5.000000e-02
#> 77     Between_Dataset_Silh_Global -9.070018e-05
#> 78                             MAE  2.473408e-01
#> 79                            RMSE  2.887192e-01
#> 80                         KDE_ISE  1.208325e-02
#> 81                     Wasserstein  2.036911e+00
#> 82                              KS  1.533879e-01
#> 83                             MAD  1.969816e+00
#> 84                   Bhattacharyya  3.196235e-04
#> 85                  Runs_Statistic -3.666923e-01
#> 86                     NN_Mismatch  4.166667e-02
#> 87     Between_Dataset_Silh_Global  9.752197e-04
#> 88                             MAE  2.036911e+00
#> 89                            RMSE  2.271370e+00
#> 90                         KDE_ISE  1.659081e-03
#> 91                     Wasserstein  1.119261e-01
#> 92                              KS  6.742514e-02
#> 93                             MAD  7.339222e-02
#> 94                   Bhattacharyya  2.010495e-03
#> 95                  Runs_Statistic -1.283423e+00
#> 96                     NN_Mismatch  1.666667e-02
#> 97     Between_Dataset_Silh_Global -1.234822e-02
#> 98                             MAE  1.119261e-01
#> 99                            RMSE  1.504922e-01
#> 100                        KDE_ISE  6.265612e-03
#> 101                    Wasserstein  7.432972e+00
#> 102                             KS  2.906651e-01
#> 103                            MAD  7.484040e+00
#> 104                  Bhattacharyya  7.165290e-05
#> 105                 Runs_Statistic -1.100077e+00
#> 106                    NN_Mismatch  4.166667e-02
#> 107    Between_Dataset_Silh_Global  5.142707e-02
#> 108                            MAE  7.432972e+00
#> 109                           RMSE  7.844679e+00
#> 110                        KDE_ISE  4.198168e-03
#> 111                    Wasserstein  1.687500e-02
#> 112                             KS  1.333333e-01
#> 113                            MAD  1.250000e-02
#> 114                  Bhattacharyya  6.100311e-04
#> 115                 Runs_Statistic -4.400308e+00
#> 116                    NN_Mismatch  9.166667e-02
#> 117    Between_Dataset_Silh_Global -1.887003e-03
#> 118                            MAE  1.687500e-02
#> 119                           RMSE  1.916214e-02
#> 120                        KDE_ISE  1.732278e-01
#> 121                    Wasserstein  4.456931e-01
#> 122                             KS  1.741354e-01
#> 123                            MAD  4.388661e-01
#> 124                  Bhattacharyya  5.586790e-04
#> 125                 Runs_Statistic  1.466769e+00
#> 126                    NN_Mismatch  2.500000e-02
#> 127    Between_Dataset_Silh_Global  1.024133e-02
#> 128                            MAE  4.456931e-01
#> 129                           RMSE  4.890173e-01
#> 130                        KDE_ISE  1.955467e-02
#> 131                    Wasserstein  3.299875e-03
#> 132                             KS  2.672058e-02
#> 133                            MAD  2.118671e-03
#> 134                  Bhattacharyya  6.113770e-05
#> 135                 Runs_Statistic -2.050782e+00
#> 136                    NN_Mismatch  3.000000e-02
#> 137    Between_Dataset_Silh_Global -2.482463e-06
#> 138                            MAE  3.299875e-03
#> 139                           RMSE  5.067781e-03
#> 140                        KDE_ISE  1.403010e-02
#> 141                    slope_error  2.560909e-01
#> 142                 midpoint_error  2.624144e-01
#> 143                 r2_discrepancy  1.646825e-01
#> 144 excess_zero_weight_discrepancy  3.712431e-02
#> 145                            MMD  1.103894e-01
#> 146           Generative_Precision  7.625000e-01
#> 147              Generative_Recall  9.625000e-01
#> 148                 silhouette_sim  3.402312e-03
#> 149                            cdi  4.984432e+00
#> 150            neighborhood_purity  4.775000e-01
#> 151               batch_silhouette -2.710457e-03
#> 152                         pcr_r2  9.476487e-03
#> 153                            cms  0.000000e+00
#> 154                            isi  1.938397e+00
#> 155                 cross_batch_F1  5.126418e-01
#> 156                pval_dist_score  8.719073e-01
#> 157                       accuracy  9.500000e-01
#> 158                      precision  0.000000e+00
#> 159                         recall  0.000000e+00
#> 160                       f1_score  0.000000e+00
#> 161                            auc  6.465517e-01
#> 
#> $mod1_unimodal$cell_metrics
#> $mod1_unimodal$cell_metrics$library_size
#> $mod1_unimodal$cell_metrics$library_size$library_size_MAD
#> [1] 15
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_KS
#> [1] 0.2125
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_MAE
#> [1] 15.0625
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_RMSE
#> [1] 15.71504
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_Bhattacharyya
#> [1] 0.0001188839
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_Wasserstein
#> [1] 15.0625
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_Runs_Statistic
#> [1] -1.744748
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_Runs_PValue
#> [1] 0.04051437
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_NN_Mismatch
#> [1] 0.175
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_Between_Dataset_Silh_Global
#> [1] 0.02247078
#> 
#> $mod1_unimodal$cell_metrics$library_size$library_size_KDE_ISE
#> [1] 0.000798037
#> 
#> 
#> $mod1_unimodal$cell_metrics$log_library_size
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_MAD
#> [1] 0.06363671
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_KS
#> [1] 0.2125
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_MAE
#> [1] 0.06194965
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_RMSE
#> [1] 0.0656695
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_Bhattacharyya
#> [1] 3.834591e-06
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_Wasserstein
#> [1] 0.06194965
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_Runs_Statistic
#> [1] -1.744748
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_Runs_PValue
#> [1] 0.04051437
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_NN_Mismatch
#> [1] 0.1875
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_Between_Dataset_Silh_Global
#> [1] 0.02176201
#> 
#> $mod1_unimodal$cell_metrics$log_library_size$log_library_size_KDE_ISE
#> [1] 0.1962832
#> 
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_MAD
#> [1] 0.02555864
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_KS
#> [1] 0.101838
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_MAE
#> [1] 0.03268387
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_RMSE
#> [1] 0.04613293
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_Bhattacharyya
#> [1] 0.0002186586
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_Wasserstein
#> [1] 0.03268387
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_Runs_Statistic
#> [1] -0.9516807
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_Runs_PValue
#> [1] 0.1706295
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_NN_Mismatch
#> [1] 0.05
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_Between_Dataset_Silh_Global
#> [1] -0.0006195226
#> 
#> $mod1_unimodal$cell_metrics$tmm_factor$tmm_factor_KDE_ISE
#> [1] 0.07601569
#> 
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_MAD
#> [1] 0.01666667
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_KS
#> [1] 0.125
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_MAE
#> [1] 0.015625
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_RMSE
#> [1] 0.01980811
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Bhattacharyya
#> [1] 0.0008854773
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Wasserstein
#> [1] 0.015625
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Runs_Statistic
#> [1] -9.358194
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Runs_PValue
#> [1] 4.055568e-21
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_NN_Mismatch
#> [1] 0.6875
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Between_Dataset_Silh_Global
#> [1] 0.00752468
#> 
#> $mod1_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_KDE_ISE
#> [1] 0.3253205
#> 
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_MAD
#> [1] 0.01666667
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_KS
#> [1] 0.125
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_MAE
#> [1] 0.015625
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_RMSE
#> [1] 0.01980811
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Bhattacharyya
#> [1] 2.679272e-05
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Wasserstein
#> [1] 0.015625
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Runs_Statistic
#> [1] -9.040967
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Runs_PValue
#> [1] 7.764603e-20
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_NN_Mismatch
#> [1] 0.6875
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Between_Dataset_Silh_Global
#> [1] 0.00752468
#> 
#> $mod1_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_KDE_ISE
#> [1] 0.3253205
#> 
#> 
#> $mod1_unimodal$cell_metrics$cell_cor
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_MAD
#> [1] 0.005684572
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_KS
#> [1] 0.03323906
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_MAE
#> [1] 0.005728994
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_RMSE
#> [1] 0.006328019
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_Bhattacharyya
#> [1] 2.566053e-05
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_Wasserstein
#> [1] 0.005728994
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_Runs_Statistic
#> [1] -2.188896
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_Runs_PValue
#> [1] 0.01430221
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_NN_Mismatch
#> [1] 0.001666667
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_Between_Dataset_Silh_Global
#> [1] 0.001423443
#> 
#> $mod1_unimodal$cell_metrics$cell_cor$cell_cor_KDE_ISE
#> [1] 0.006120211
#> 
#> 
#> $mod1_unimodal$cell_metrics$pca_dist
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_MAD
#> [1] 1.428096
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_KS
#> [1] 0.08747786
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_MAE
#> [1] 1.59434
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_RMSE
#> [1] 1.676818
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_Bhattacharyya
#> [1] 3.10537e-05
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_Wasserstein
#> [1] 1.59434
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_Runs_Statistic
#> [1] 0.02515972
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_Runs_PValue
#> [1] 0.5100362
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_NN_Mismatch
#> [1] 0.015
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_Between_Dataset_Silh_Global
#> [1] 0.009907338
#> 
#> $mod1_unimodal$cell_metrics$pca_dist$pca_dist_KDE_ISE
#> [1] 0.000858538
#> 
#> 
#> 
#> $mod1_unimodal$feature_metrics
#> $mod1_unimodal$feature_metrics$mean_expression
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_MAD
#> [1] 0.2571903
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_KS
#> [1] 0.2055826
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_MAE
#> [1] 0.2473408
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_RMSE
#> [1] 0.2887192
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_Bhattacharyya
#> [1] 3.690266e-05
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_Wasserstein
#> [1] 0.2473408
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_Runs_Statistic
#> [1] 0.3666923
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_Runs_PValue
#> [1] 0.6430757
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_NN_Mismatch
#> [1] 0.05
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_Between_Dataset_Silh_Global
#> [1] -9.070018e-05
#> 
#> $mod1_unimodal$feature_metrics$mean_expression$mean_expression_KDE_ISE
#> [1] 0.01208325
#> 
#> 
#> $mod1_unimodal$feature_metrics$variance
#> $mod1_unimodal$feature_metrics$variance$variance_MAD
#> [1] 1.969816
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_KS
#> [1] 0.1533879
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_MAE
#> [1] 2.036911
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_RMSE
#> [1] 2.27137
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_Bhattacharyya
#> [1] 0.0003196235
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_Wasserstein
#> [1] 2.036911
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_Runs_Statistic
#> [1] -0.3666923
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_Runs_PValue
#> [1] 0.3569243
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_NN_Mismatch
#> [1] 0.04166667
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_Between_Dataset_Silh_Global
#> [1] 0.0009752197
#> 
#> $mod1_unimodal$feature_metrics$variance$variance_KDE_ISE
#> [1] 0.001659081
#> 
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_MAD
#> [1] 0.07339222
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_KS
#> [1] 0.06742514
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_MAE
#> [1] 0.1119261
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_RMSE
#> [1] 0.1504922
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_Bhattacharyya
#> [1] 0.002010495
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_Wasserstein
#> [1] 0.1119261
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_Runs_Statistic
#> [1] -1.283423
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_Runs_PValue
#> [1] 0.09967192
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_NN_Mismatch
#> [1] 0.01666667
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_Between_Dataset_Silh_Global
#> [1] -0.01234822
#> 
#> $mod1_unimodal$feature_metrics$scaled_variance$scaled_variance_KDE_ISE
#> [1] 0.006265612
#> 
#> 
#> $mod1_unimodal$feature_metrics$cv
#> $mod1_unimodal$feature_metrics$cv$cv_MAD
#> [1] 7.48404
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_KS
#> [1] 0.2906651
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_MAE
#> [1] 7.432972
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_RMSE
#> [1] 7.844679
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_Bhattacharyya
#> [1] 7.16529e-05
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_Wasserstein
#> [1] 7.432972
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_Runs_Statistic
#> [1] -1.100077
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_Runs_PValue
#> [1] 0.1356493
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_NN_Mismatch
#> [1] 0.04166667
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_Between_Dataset_Silh_Global
#> [1] 0.05142707
#> 
#> $mod1_unimodal$feature_metrics$cv$cv_KDE_ISE
#> [1] 0.004198168
#> 
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_MAD
#> [1] 0.0125
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_KS
#> [1] 0.1333333
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_MAE
#> [1] 0.016875
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_RMSE
#> [1] 0.01916214
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Bhattacharyya
#> [1] 0.0006100311
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Wasserstein
#> [1] 0.016875
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Runs_Statistic
#> [1] -4.400308
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Runs_PValue
#> [1] 5.404863e-06
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_NN_Mismatch
#> [1] 0.09166667
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Between_Dataset_Silh_Global
#> [1] -0.001887003
#> 
#> $mod1_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_KDE_ISE
#> [1] 0.1732278
#> 
#> 
#> $mod1_unimodal$feature_metrics$dispersion
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_MAD
#> [1] 0.4388661
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_KS
#> [1] 0.1741354
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_MAE
#> [1] 0.4456931
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_RMSE
#> [1] 0.4890173
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_Bhattacharyya
#> [1] 0.000558679
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_Wasserstein
#> [1] 0.4456931
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_Runs_Statistic
#> [1] 1.466769
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_Runs_PValue
#> [1] 0.9287806
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_NN_Mismatch
#> [1] 0.025
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_Between_Dataset_Silh_Global
#> [1] 0.01024133
#> 
#> $mod1_unimodal$feature_metrics$dispersion$dispersion_KDE_ISE
#> [1] 0.01955467
#> 
#> 
#> $mod1_unimodal$feature_metrics$gene_cor
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_MAD
#> [1] 0.002118671
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_KS
#> [1] 0.02672058
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_MAE
#> [1] 0.003299875
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_RMSE
#> [1] 0.005067781
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_Bhattacharyya
#> [1] 6.11377e-05
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_Wasserstein
#> [1] 0.003299875
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_Runs_Statistic
#> [1] -2.050782
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_Runs_PValue
#> [1] 0.02014407
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_NN_Mismatch
#> [1] 0.03
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_Between_Dataset_Silh_Global
#> [1] -2.482463e-06
#> 
#> $mod1_unimodal$feature_metrics$gene_cor$gene_cor_KDE_ISE
#> [1] 0.0140301
#> 
#> 
#> 
#> $mod1_unimodal$bivariate_metrics
#> list()
#> 
#> $mod1_unimodal$zero_probability_curve
#>    slope_error midpoint_error r2_discrepancy 
#>      0.2560909      0.2624144      0.1646825 
#> 
#> $mod1_unimodal$accessibility_sparsity_curve
#> NULL
#> 
#> $mod1_unimodal$excess_zero_weights
#> $mod1_unimodal$excess_zero_weights$ref
#> $mod1_unimodal$excess_zero_weights$ref$mean_excess_zero_weight
#> [1] 0.125175
#> 
#> $mod1_unimodal$excess_zero_weights$ref$median_excess_zero_weight
#> [1] 0.05941929
#> 
#> $mod1_unimodal$excess_zero_weights$ref$mean_zero_inflation
#> [1] 0.02386441
#> 
#> $mod1_unimodal$excess_zero_weights$ref$gene_excess_weights
#>  [1] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.09777775
#>  [7] 0.40423218 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
#> [13] 0.00000000 0.44714706 0.08562163 0.18806978 0.27267149 0.06694220
#> [19] 0.00000000 0.00000000 0.01533320 0.00000000 0.52388338 0.00000000
#> [25] 0.00000000 0.00000000 0.00000000 0.00000000 0.49862464 0.28463070
#> [31] 0.00000000 0.02307153 0.31073986 0.00000000 0.06367513 0.15088815
#> [37] 0.24689382 0.00000000 0.09061268 0.31942074 0.21672462 0.18998773
#> [43] 0.07329427 0.00000000 0.22786326 0.00000000 0.34310829 0.06448393
#> [49] 0.05516346 0.00000000 0.30971601 0.19499308 0.28160864 0.37442912
#> [55] 0.00000000 0.27711280 0.26124999 0.45666487 0.09386449 0.00000000
#> 
#> 
#> $mod1_unimodal$excess_zero_weights$sim
#> $mod1_unimodal$excess_zero_weights$sim$mean_excess_zero_weight
#> [1] 0.08805069
#> 
#> $mod1_unimodal$excess_zero_weights$sim$median_excess_zero_weight
#> [1] 0
#> 
#> $mod1_unimodal$excess_zero_weights$sim$mean_zero_inflation
#> [1] 0.0174943
#> 
#> $mod1_unimodal$excess_zero_weights$sim$gene_excess_weights
#>  [1] 0.00000000 0.17239956 0.35154270 0.00000000 0.05919581 0.00000000
#>  [7] 0.38441563 0.00000000 0.25623390 0.01706470 0.00000000 0.00000000
#> [13] 0.00000000 0.20065580 0.00000000 0.00000000 0.20684773 0.00000000
#> [19] 0.00000000 0.13521241 0.00000000 0.21730222 0.30087690 0.13531351
#> [25] 0.37070733 0.00000000 0.00000000 0.00000000 0.00000000 0.28294619
#> [31] 0.00000000 0.00000000 0.00000000 0.00000000 0.39060049 0.00000000
#> [37] 0.00000000 0.26814034 0.00000000 0.05538841 0.00000000 0.45866800
#> [43] 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.07892515
#> [49] 0.00000000 0.19200186 0.00000000 0.00000000 0.00000000 0.00000000
#> [55] 0.00000000 0.28262082 0.00000000 0.22806525 0.00000000 0.23791688
#> 
#> 
#> 
#> $mod1_unimodal$generative_manifold
#> $mod1_unimodal$generative_manifold$mmd
#> [1] 0.1103894
#> 
#> $mod1_unimodal$generative_manifold$precision_recall
#> $mod1_unimodal$generative_manifold$precision_recall$generative_precision
#> [1] 0.7625
#> 
#> $mod1_unimodal$generative_manifold$precision_recall$generative_recall
#> [1] 0.9625
#> 
#> $mod1_unimodal$generative_manifold$precision_recall$generative_f1
#> [1] 0.8509058
#> 
#> 
#> 
#> $mod1_unimodal$ref_properties
#> $mod1_unimodal$ref_properties$cell
#> $mod1_unimodal$ref_properties$cell$library_size
#>  [1] 225 208 262 227 255 256 287 258 258 209 220 245 273 238 285 224 301 229 213
#> [20] 243 261 230 184 219 343 267 270 172 264 194 183 209 271 222 237 196 265 247
#> [39] 302 227 307 221 301 272 273 303 268 320 315 276 269 283 301 227 303 204 310
#> [58] 229 250 300 220 278 227 317 256 279 321 251 301 256 221 307 226 251 294 253
#> [77] 273 272 209 282
#> 
#> $mod1_unimodal$ref_properties$cell$log_library_size
#>  [1] 5.420535 5.342334 5.572154 5.429346 5.545177 5.549076 5.662960 5.556828
#>  [9] 5.556828 5.347108 5.398163 5.505332 5.613128 5.476464 5.655992 5.416100
#> [17] 5.710427 5.438079 5.365976 5.497168 5.568345 5.442418 5.220356 5.393628
#> [25] 5.840642 5.590987 5.602119 5.153292 5.579730 5.273000 5.214936 5.347108
#> [33] 5.605802 5.407172 5.472271 5.283204 5.583496 5.513429 5.713733 5.429346
#> [41] 5.730100 5.402677 5.710427 5.609472 5.613128 5.717028 5.594711 5.771441
#> [49] 5.755742 5.624018 5.598422 5.648974 5.710427 5.429346 5.717028 5.323010
#> [57] 5.739793 5.438079 5.525453 5.707110 5.398163 5.631212 5.429346 5.762051
#> [65] 5.549076 5.634790 5.774552 5.529429 5.710427 5.549076 5.402677 5.730100
#> [73] 5.424950 5.529429 5.686975 5.537334 5.613128 5.609472 5.347108 5.645447
#> 
#> $mod1_unimodal$ref_properties$cell$tmm_factor
#>  [1] 1.2677665 1.0608912 0.9049099 0.8659920 0.8541085 1.1025580 0.8813490
#>  [8] 1.0676245 1.0128242 0.9703222 0.9763756 1.1185974 0.8624278 1.5618503
#> [15] 0.9140547 1.0148082 0.9259660 1.0562338 0.8866491 0.8154473 0.8652361
#> [22] 0.7406297 0.9973933 1.0481011 0.9407201 1.1001742 0.9158304 1.2243780
#> [29] 0.8246904 1.1039272 1.0855407 1.0795148 1.1816197 1.1433863 1.3097731
#> [36] 1.2233884 0.7238186 1.0035751 0.9711287 1.0911025 0.8331748 0.9351113
#> [43] 1.0097258 0.9119055 1.0008139 1.1945043 1.0330545 0.8741315 0.9652980
#> [50] 0.8000045 1.3516526 1.2151585 0.9005631 0.9530173 1.1226875 0.9858696
#> [57] 0.7789178 0.8780626 0.8889611 1.0143402 1.0904013 0.9795808 0.7353033
#> [64] 1.1582880 1.0355863 1.0079076 0.9840788 1.1396867 0.7843954 1.0134781
#> [71] 1.2229822 0.8139370 0.8860124 0.9547651 0.9345196 1.2740610 1.2514075
#> [78] 1.0676455 1.3134852 0.9530813
#> 
#> $mod1_unimodal$ref_properties$cell$effective_library_size
#>  [1] 285.2475 220.6654 237.0864 196.5802 217.7977 282.2549 252.9472 275.4471
#>  [9] 261.3086 202.7973 214.8026 274.0564 235.4428 371.7204 260.5056 227.3170
#> [17] 278.7158 241.8775 188.8563 198.1537 225.8266 170.3448 183.5204 229.5341
#> [25] 322.6670 293.7465 247.2742 210.5930 217.7183 214.1619 198.6540 225.6186
#> [33] 320.2189 253.8318 310.4162 239.7841 191.8119 247.8831 293.2809 247.6803
#> [41] 255.7847 206.6596 303.9275 248.0383 273.2222 361.9348 276.8586 279.7221
#> [49] 304.0689 220.8012 363.5946 343.8898 271.0695 216.3349 340.1743 201.1174
#> [57] 241.4645 201.0763 222.2403 304.3021 239.8883 272.3235 166.9138 367.1773
#> [65] 265.1101 281.2062 315.8893 286.0614 236.1030 259.4504 270.2791 249.8787
#> [73] 200.2388 239.6460 274.7488 322.3374 341.6343 290.3996 274.5184 268.7689
#> 
#> $mod1_unimodal$ref_properties$cell$zero_fraction_cell
#>  [1] 0.21666667 0.18333333 0.20000000 0.21666667 0.15000000 0.23333333
#>  [7] 0.08333333 0.13333333 0.15000000 0.16666667 0.10000000 0.23333333
#> [13] 0.18333333 0.25000000 0.16666667 0.20000000 0.18333333 0.18333333
#> [19] 0.21666667 0.11666667 0.13333333 0.25000000 0.20000000 0.26666667
#> [25] 0.13333333 0.21666667 0.10000000 0.28333333 0.06666667 0.25000000
#> [31] 0.21666667 0.21666667 0.15000000 0.23333333 0.21666667 0.18333333
#> [37] 0.10000000 0.20000000 0.11666667 0.18333333 0.18333333 0.21666667
#> [43] 0.15000000 0.16666667 0.21666667 0.18333333 0.21666667 0.10000000
#> [49] 0.16666667 0.16666667 0.23333333 0.16666667 0.20000000 0.20000000
#> [55] 0.23333333 0.21666667 0.16666667 0.18333333 0.18333333 0.18333333
#> [61] 0.13333333 0.13333333 0.20000000 0.18333333 0.13333333 0.10000000
#> [67] 0.18333333 0.16666667 0.15000000 0.16666667 0.28333333 0.10000000
#> [73] 0.18333333 0.13333333 0.10000000 0.23333333 0.15000000 0.23333333
#> [79] 0.23333333 0.11666667
#> 
#> $mod1_unimodal$ref_properties$cell$detection_freq_cell
#>  [1] 0.7833333 0.8166667 0.8000000 0.7833333 0.8500000 0.7666667 0.9166667
#>  [8] 0.8666667 0.8500000 0.8333333 0.9000000 0.7666667 0.8166667 0.7500000
#> [15] 0.8333333 0.8000000 0.8166667 0.8166667 0.7833333 0.8833333 0.8666667
#> [22] 0.7500000 0.8000000 0.7333333 0.8666667 0.7833333 0.9000000 0.7166667
#> [29] 0.9333333 0.7500000 0.7833333 0.7833333 0.8500000 0.7666667 0.7833333
#> [36] 0.8166667 0.9000000 0.8000000 0.8833333 0.8166667 0.8166667 0.7833333
#> [43] 0.8500000 0.8333333 0.7833333 0.8166667 0.7833333 0.9000000 0.8333333
#> [50] 0.8333333 0.7666667 0.8333333 0.8000000 0.8000000 0.7666667 0.7833333
#> [57] 0.8333333 0.8166667 0.8166667 0.8166667 0.8666667 0.8666667 0.8000000
#> [64] 0.8166667 0.8666667 0.9000000 0.8166667 0.8333333 0.8500000 0.8333333
#> [71] 0.7166667 0.9000000 0.8166667 0.8666667 0.9000000 0.7666667 0.8500000
#> [78] 0.7666667 0.7666667 0.8833333
#> 
#> $mod1_unimodal$ref_properties$cell$cell_cor
#>    [1] -0.0370302447  0.2154919205  0.0030832110  0.1234317559  0.1048672306
#>    [6]  0.1596921410 -0.1281623967  0.0642524208  0.0356629562 -0.1477820293
#>   [11] -0.0933382455  0.2173552646  0.1139414435  0.1886852279  0.1063877047
#>   [16] -0.0924150166  0.0019122672  0.0226358914  0.1282094443  0.0086084125
#>   [21]  0.1786275089  0.1421893740  0.0378730079  0.0472881570  0.3416191274
#>   [26] -0.0551530725 -0.0933281991 -0.2392610396  0.0507451996 -0.0369129581
#>   [31] -0.0633330482  0.0401035619 -0.0011681429  0.2054474700  0.0667712188
#>   [36]  0.1244669181  0.0166834265  0.0044672003 -0.0998619063 -0.1951267540
#>   [41]  0.1271586244  0.1351757256  0.2419450083 -0.1532736630  0.1908821758
#>   [46]  0.0000000000  0.0513536521  0.2951821814  0.2189161641  0.0206349748
#>   [51]  0.1822016112  0.1444826589  0.1550796671  0.2181696080  0.2018369954
#>   [56]  0.0439943388  0.2308064179  0.1252248265  0.2050004392  0.0157979413
#>   [61]  0.1606230612 -0.1216955461  0.0962446718 -0.0073393013 -0.0301514660
#>   [66]  0.1337413462 -0.0206716883  0.2553700394 -0.0208484033  0.1451659037
#>   [71]  0.2972027462  0.2560401048 -0.0462596387 -0.0307030017  0.0225170769
#>   [76]  0.0185256008  0.0023015416  0.2069239213  0.2265370137  0.0548828097
#>   [81] -0.1432157715  0.1613474651 -0.0123330192  0.1787533984  0.0013740189
#>   [86]  0.0433345798  0.3853807094  0.1969992405  0.1383580916 -0.0624929340
#>   [91]  0.0001134238  0.0016510355  0.0704260779  0.1447758000 -0.0921421407
#>   [96]  0.0589577062 -0.0408239873  0.0188932532 -0.0782262727 -0.0072837721
#>  [101]  0.0889591083 -0.0773822665 -0.1721949566 -0.0704340599  0.0168578458
#>  [106]  0.2518817221 -0.2785489299  0.2389652995  0.0578343173 -0.1824054653
#>  [111] -0.2500073558  0.0685569142  0.1677094849 -0.0435708454 -0.1446760700
#>  [116]  0.1736098872 -0.0817770775 -0.1585282048 -0.2375491504 -0.0203573305
#>  [121] -0.0545338760  0.0396232112  0.1016393112  0.1107306235  0.0726558412
#>  [126] -0.0163275625  0.0095913991  0.0343598672  0.1734305264  0.1085666655
#>  [131] -0.0693102711 -0.2017382956 -0.1012103286 -0.0894680213  0.0701561895
#>  [136] -0.1756140661 -0.0732605434  0.2445687432 -0.0447303109  0.0098472661
#>  [141]  0.1376099900 -0.0607763875  0.0919614024  0.1358436798  0.0718857204
#>  [146]  0.0399076330  0.0136988911  0.0459910177  0.1894146528  0.1196483695
#>  [151]  0.0790897749  0.1213633938  0.0043634429 -0.1396955978  0.0564760000
#>  [156]  0.0634790524 -0.0337835982 -0.0334651154  0.3453808728 -0.0129193442
#>  [161] -0.1193855534 -0.0492914097 -0.0305209446  0.1758525084  0.0466950572
#>  [166] -0.0908150796  0.1206485322  0.0932109522 -0.0920329906  0.0470903624
#>  [171] -0.2022522578 -0.1555878616  0.1715009430 -0.1508137035  0.0176655013
#>  [176]  0.0708337363  0.2176438025 -0.0878595751  0.0663944716  0.1908001280
#>  [181]  0.0668381188  0.2020735247  0.0646884672  0.2682030777  0.1588978510
#>  [186]  0.0977595628 -0.1556200951 -0.1242039931  0.2140235619 -0.0775987422
#>  [191] -0.1051761936  0.2936556988  0.0314113395 -0.0392664772  0.1803859995
#>  [196]  0.1750855982 -0.0274336446 -0.0121903072  0.1004355672  0.0881680335
#>  [201]  0.1639707201  0.0830588253  0.1640298597  0.0163521911  0.0281982424
#>  [206] -0.1655348378  0.0486932078  0.1379902357  0.2169797054  0.2528398607
#>  [211]  0.0952802756  0.0554853980  0.0021702920 -0.0609083738  0.0430512168
#>  [216] -0.1120756900 -0.0613163386  0.0683641157  0.0918251062  0.0991800945
#>  [221] -0.0390926105 -0.0086489585 -0.0055883106  0.1261929852  0.2192146678
#>  [226]  0.0348324199  0.0623678124  0.3142723758 -0.1763071510  0.1669615270
#>  [231] -0.0300741760  0.0505013780 -0.0455650981  0.1967365078  0.0337449838
#>  [236]  0.0468683300  0.0624905555  0.0071354068 -0.0068149466  0.1516172535
#>  [241] -0.0017042844  0.0909906993  0.1341726976  0.0706661775  0.1064148521
#>  [246] -0.0394675290  0.1971291285 -0.0939156551  0.1083752164 -0.0576535379
#>  [251]  0.0664778602 -0.1186297257  0.2805284554 -0.0387402738 -0.1561226870
#>  [256]  0.1009124181  0.1870544021  0.1556966832  0.0211377791  0.0867511309
#>  [261] -0.0782551160 -0.0770357770 -0.1708792425 -0.0942175675 -0.0470975916
#>  [266]  0.0035078407 -0.0067973134  0.2389822924 -0.0964900340  0.1027587008
#>  [271] -0.0125147849 -0.1664867514  0.1403741307 -0.0685250939  0.1092761429
#>  [276] -0.1758900644  0.1628403228  0.3866422817  0.0998201286 -0.0662673807
#>  [281] -0.0200427904  0.0770786733 -0.1036786124  0.1360237145  0.0601565145
#>  [286]  0.3626810409  0.0917708157  0.2532997639  0.2212074908 -0.0335865574
#>  [291] -0.0439841666 -0.2063684672  0.2176035784  0.0049501259  0.1183432163
#>  [296] -0.0311627924  0.2095156562  0.1974245188  0.0196948172 -0.1408931646
#>  [301] -0.1400181662  0.3043646286 -0.0992141254  0.3024842984 -0.1350982843
#>  [306]  0.1131032527  0.0388083328 -0.0952079103  0.0786472047 -0.0263035228
#>  [311] -0.0538198042  0.2608045980  0.2375308828  0.0733550556 -0.1024942571
#>  [316] -0.1751656930 -0.0420857472  0.2450020957  0.0204680585  0.2998116830
#>  [321]  0.2484356478  0.0520065064 -0.0160852192  0.0531897960  0.0857107150
#>  [326]  0.0432992016 -0.1682420226  0.1218134824 -0.0037362147 -0.0150418313
#>  [331] -0.0435096885  0.2286690091 -0.0605381973 -0.0396145562  0.3018232868
#>  [336]  0.0059213825 -0.2675117945  0.0555848098 -0.1284063512 -0.0589311404
#>  [341] -0.0189390069  0.0381626736  0.0956622224 -0.3031709975  0.0086184905
#>  [346] -0.0371916373  0.1871565083  0.0530906750 -0.0374126739  0.1012372636
#>  [351]  0.0188860552 -0.0184722168  0.1871833535  0.0947145715  0.1651965820
#>  [356]  0.2372058973  0.2112582318  0.2357624263  0.1185401657  0.0103253423
#>  [361]  0.0119481179  0.1268152399  0.2452100691  0.1253644662 -0.0080069655
#>  [366] -0.1288668588  0.0038384208  0.0225261227  0.0732264857  0.1101484164
#>  [371] -0.0017288097  0.1022265751  0.0423139273  0.2462463182  0.0725983141
#>  [376]  0.0432765782  0.0148522718 -0.1789257580 -0.0502551211  0.0212866086
#>  [381] -0.0626868412  0.1239736170  0.0294521905  0.0519975663 -0.0969062654
#>  [386]  0.0319979177  0.0667879005 -0.2711450335 -0.1754036142  0.0375201884
#>  [391]  0.3650838729 -0.1138145852 -0.1881459292  0.0348867579 -0.1556960716
#>  [396] -0.0546308808 -0.0163109374  0.1986164236  0.1260541250 -0.1139064551
#>  [401]  0.1429731024  0.1901377357 -0.0721758491  0.1425380223 -0.0185639363
#>  [406] -0.0928812706 -0.1098756688  0.1261296453  0.0148268384 -0.0444083810
#>  [411]  0.1368860218 -0.1286566310 -0.0314739033  0.1403031230  0.2079846811
#>  [416]  0.0558582773  0.0808668389 -0.1089162090  0.1335975102 -0.0995677144
#>  [421]  0.0211439143 -0.0197721826 -0.0719289709  0.0801323269 -0.0245020493
#>  [426]  0.0414613586  0.0829642535  0.1589800342  0.0776410135 -0.2953710207
#>  [431]  0.1156705384 -0.0202710523  0.0670590879 -0.0063511448  0.0545682660
#>  [436]  0.0695180971 -0.0473462528  0.0954728116  0.3065883687 -0.0242364956
#>  [441]  0.0726161700  0.2001205728 -0.0533198234  0.1058454400 -0.0556922514
#>  [446]  0.0852282554 -0.1918864754  0.3106907906  0.0757977348  0.1452455566
#>  [451]  0.0258284108  0.1570538697  0.0376242270  0.0225799360 -0.0192433624
#>  [456]  0.0648428451  0.0463024152  0.2248754102 -0.0161648676  0.0386011806
#>  [461]  0.0440891820  0.1025516301  0.0499882138  0.1708418198  0.0831481043
#>  [466] -0.0645135614 -0.1412798605  0.2291665895  0.1115704995  0.2777266317
#>  [471]  0.1691285946  0.3356368082 -0.1184354462 -0.0334291312  0.0036341514
#>  [476]  0.1523653054  0.0279411018 -0.1477821737 -0.0427312987  0.2464854881
#>  [481]  0.0301192212 -0.0733974541  0.0523955968  0.0617201761  0.0527626039
#>  [486]  0.2723967732  0.0330511994 -0.1219901077  0.2637505598 -0.3641952949
#>  [491]  0.0155249521 -0.0271700496  0.0094864819  0.0194292261 -0.0238556592
#>  [496] -0.0970795166  0.0801505042 -0.0147107623 -0.0262961347  0.0965722694
#>  [501] -0.1788165345  0.3165912244  0.0648676760  0.0965477398 -0.0380609519
#>  [506]  0.0292078016  0.1808907631  0.1148958441  0.2041453182  0.1732358308
#>  [511] -0.0342555619 -0.0896587079 -0.1561615593 -0.0758505257  0.0718190791
#>  [516]  0.1024479481  0.1180860684 -0.1444916834 -0.0415621402 -0.2141915647
#>  [521]  0.1697874835  0.0033870588 -0.0123674186 -0.1025481765 -0.0159535223
#>  [526] -0.0760574357  0.1350389122 -0.0407644589  0.0528344274  0.1365416507
#>  [531] -0.1049363452  0.1752211281  0.1103499524  0.0611155548  0.0695084342
#>  [536]  0.2472792173 -0.0161512700  0.0255966915  0.1374633442  0.0218010698
#>  [541] -0.2097343757  0.1337060249  0.0547651758  0.0855035519  0.1824718269
#>  [546]  0.0955242594  0.0750436467  0.0120297287  0.0289066246  0.0278755823
#>  [551]  0.1411463475  0.0804095049 -0.1178269347  0.0306886981 -0.3279054739
#>  [556]  0.4466331454 -0.0449779442 -0.1540950614  0.0653358520  0.1368641076
#>  [561] -0.0569702465 -0.0237667536  0.0838907152 -0.0821949805  0.0456387842
#>  [566]  0.3631454401  0.0694224522  0.1609538218  0.0859192638 -0.0097548974
#>  [571]  0.2960773908  0.2446257626 -0.0343358834  0.1203825950  0.1473703985
#>  [576] -0.1130432213 -0.1858999495  0.1083058340 -0.0569864196  0.1402952934
#>  [581]  0.1679461731  0.2745253648  0.1810494270 -0.1042013266 -0.0135088181
#>  [586]  0.2040471680 -0.1080802241  0.1258184589  0.0750989512 -0.0244198345
#>  [591]  0.2323014896  0.0087539430  0.1187878949  0.0045308955 -0.0525384990
#>  [596]  0.0272104057 -0.2011546361  0.0581912041  0.1270403046 -0.1223597952
#>  [601]  0.1103465262  0.0614226817 -0.0006219541  0.0256177105  0.0858944396
#>  [606]  0.1952309457 -0.0861463028 -0.1888970905  0.0351195847  0.0913730709
#>  [611]  0.1187198018  0.1171143489  0.0130495773  0.0736900661  0.1563159614
#>  [616]  0.1071081830  0.1670661965  0.1205258952 -0.0706412858 -0.0454620936
#>  [621] -0.0058021115  0.2416431622 -0.0720563850 -0.2810755190  0.0746883295
#>  [626]  0.0614542804  0.1993699549  0.1452138079 -0.0089162949  0.1925609782
#>  [631]  0.0711662535  0.0224710721  0.0054417556 -0.0744557062 -0.0413206885
#>  [636] -0.1519729377  0.0968175130  0.0849233565  0.0690549595 -0.0319334675
#>  [641]  0.2784906306  0.0466551048 -0.1112777635  0.0880742683  0.0123103918
#>  [646]  0.2065086657 -0.0986205766  0.3099536531 -0.0848810190  0.0899801417
#>  [651] -0.0132795207  0.0606721509  0.0908521115  0.0134516146  0.0297069516
#>  [656]  0.0736129079 -0.1345593534 -0.0449705705 -0.1611290001 -0.0804272194
#>  [661] -0.0722230411  0.0437126759  0.0218543810  0.1164462132 -0.0003812538
#>  [666]  0.1566242119  0.1772659656  0.3153879913  0.1402754954  0.1152361558
#>  [671] -0.1659237756  0.1971925185  0.0409054318  0.1515383102 -0.0400198442
#>  [676] -0.0367944945  0.0350634375  0.0450758205  0.1038765264  0.1956475551
#>  [681]  0.2542925804 -0.0715662849 -0.0915142830  0.2702005868  0.0088227670
#>  [686]  0.1573943592  0.1429600375  0.0664678706  0.0122811711  0.0865693617
#>  [691]  0.2298419308  0.2099351550 -0.1106835214 -0.0225513071 -0.1419404532
#>  [696]  0.0810528346  0.0043682724  0.0649216561  0.2335164348  0.1058026244
#>  [701] -0.1701616665  0.0265708695  0.3672384324  0.0077156889  0.1541403886
#>  [706]  0.1301785446 -0.0749146465  0.0122565682  0.2906953763 -0.1044750946
#>  [711]  0.0600209128  0.2042672238  0.0068715990  0.1855400496  0.1436728366
#>  [716]  0.0368173465  0.0871351193  0.0637301823 -0.0213231136 -0.1737131380
#>  [721] -0.0270341082  0.2186623843  0.1854683045  0.2190252776  0.1236706136
#>  [726]  0.1530507666 -0.2748467394  0.1231739662  0.1609340907 -0.2859374329
#>  [731]  0.1438082201 -0.0394157468  0.3588702870 -0.0097610029 -0.0074228344
#>  [736]  0.2242068361 -0.0633528789 -0.0274701674  0.0739341248 -0.0532474472
#>  [741]  0.1295371707 -0.0418137510  0.0345605775  0.0873308624 -0.1371945958
#>  [746]  0.1256370972  0.0304636481  0.0221637297 -0.1466810083 -0.0820737013
#>  [751] -0.0595650156 -0.0897443662  0.1024546823  0.2349640647 -0.1126960497
#>  [756] -0.1485943644 -0.1028635487  0.1219694514 -0.0332340968  0.0116676843
#>  [761] -0.0205737508  0.0200855580 -0.1917759209  0.0740694385  0.1236762380
#>  [766]  0.0593532612  0.0863391706  0.1004219872 -0.0291385791  0.1077909148
#>  [771]  0.1267304849 -0.0887032786  0.0024152891 -0.1250166359 -0.2312548056
#>  [776]  0.1245647157 -0.1091510091 -0.0274890543 -0.0373769247 -0.1688527845
#>  [781] -0.0854925141 -0.0527622022  0.2508208996  0.0491972396  0.1679561457
#>  [786] -0.0730072600  0.1123077997 -0.0292654049 -0.1145520555 -0.1154355380
#>  [791] -0.0669323195  0.1222763672 -0.0740967668 -0.2701519687 -0.0156428287
#>  [796]  0.0627628647 -0.0360166288  0.0082098356 -0.1072436058 -0.0924110789
#>  [801] -0.1121926862  0.0291970807 -0.0278930338  0.1977753613 -0.0391492864
#>  [806]  0.2033280476  0.2230987436  0.1138286982 -0.0019836200  0.1391897396
#>  [811] -0.2440317701  0.1869072682 -0.2298238197 -0.1311342193  0.0739220204
#>  [816]  0.1392254546  0.1354304820  0.0821351954  0.0054044223  0.1950682203
#>  [821]  0.0149890237  0.0114909943  0.2096993688  0.1058507897  0.0687172832
#>  [826]  0.0356040406 -0.1472151168  0.0560427375  0.0261024108 -0.1495989429
#>  [831]  0.1897512295 -0.1891202401 -0.1290432509  0.0054618464  0.2137467638
#>  [836] -0.0601367321  0.1938139135 -0.1592171754  0.1194580745  0.1339347353
#>  [841]  0.1128861175  0.0333084161 -0.1435069437  0.3773440021 -0.0470364728
#>  [846] -0.0794186519 -0.0885349993 -0.0636546425  0.0259897827  0.1479255000
#>  [851] -0.0387747658  0.1357732097 -0.1052234344  0.1970453135  0.1521943829
#>  [856]  0.1233292249 -0.0683289093  0.2317808687  0.1030667048 -0.0816975018
#>  [861]  0.1036573341  0.0977474882 -0.0527237484  0.0276000205  0.0018162139
#>  [866] -0.0161492475 -0.0759065228 -0.0533970979  0.1950612329 -0.0142335270
#>  [871] -0.1932312592 -0.1302280045 -0.0911754794  0.0716930721  0.1240451808
#>  [876]  0.0006746218  0.1701345302 -0.0435443107  0.1875414285  0.1696062925
#>  [881]  0.0165237785  0.1209607889  0.0060894827  0.2044448333  0.0090944795
#>  [886] -0.0161763675  0.0702696940 -0.0631336662  0.1544654601  0.3369905191
#>  [891] -0.1021595849  0.1171790443 -0.1402533144 -0.0218633395  0.0772604621
#>  [896] -0.1375882621 -0.1592900380 -0.0259395402 -0.1498337208  0.0540593726
#>  [901]  0.1772261458 -0.2301343143 -0.1576967077  0.0493691479  0.0938246017
#>  [906]  0.0574201152  0.0416660211  0.0001126713  0.0654400213 -0.1267810820
#>  [911]  0.2236937317 -0.0046422654 -0.0631402264 -0.0058221387 -0.0500553566
#>  [916] -0.0175208986 -0.0219968589  0.2235777718  0.1338290534  0.0328129619
#>  [921]  0.0086707673  0.2365362131 -0.0032265816  0.1456769441  0.0438251987
#>  [926] -0.1003753565  0.0892680456  0.1327832931 -0.1863896236 -0.1087069754
#>  [931]  0.0236072365  0.2158257315 -0.1035412879 -0.2146998518 -0.0197829138
#>  [936] -0.0033383109 -0.0438821175  0.0815784174 -0.0005790993 -0.0471211795
#>  [941] -0.0272721260  0.2071255527 -0.0626126869  0.1006733346  0.1171973855
#>  [946]  0.2276716999 -0.3194349800 -0.1162477720  0.0126707424  0.0697805259
#>  [951] -0.0908963506  0.1423199034  0.1618783847  0.0075935118 -0.0314566593
#>  [956]  0.1331778201  0.1910099073  0.0115889243  0.1160413106 -0.1624630998
#>  [961]  0.1285669305  0.0988334807  0.1090711035 -0.1072834162  0.2610519958
#>  [966]  0.0768864454  0.0740960502 -0.0169318153 -0.0491599919  0.0713876084
#>  [971]  0.0696623247  0.0874394158  0.0356397981  0.0045422243  0.1120703515
#>  [976]  0.2020858417  0.1482018112  0.0438637519 -0.0162143983  0.0884643967
#>  [981]  0.1779728133  0.1782820460 -0.1028729117 -0.0263618900  0.1978842177
#>  [986] -0.0060179687  0.1122903430  0.3522423390  0.0333151734 -0.0306325888
#>  [991]  0.1936818748 -0.0132709690  0.1365846868  0.1450048218  0.3533018755
#>  [996] -0.0835158875 -0.0761195408  0.1001508679  0.0051951403 -0.0791346741
#> [1001]  0.1069115310  0.1866171104  0.1751120090  0.0665538642 -0.0244757899
#> [1006] -0.0659012855  0.0466744794  0.0680095129  0.0029854725  0.1827624792
#> [1011]  0.3298577052  0.1514347663  0.0608270340  0.2691993492  0.1669976996
#> [1016]  0.2977613956 -0.0313622824  0.0923156009  0.0622979612  0.0045517138
#> [1021] -0.0086862701  0.2312017750 -0.1144693838  0.2652329903  0.2420269872
#> [1026]  0.0261231779  0.1874964770  0.0987482220  0.0560256420  0.1390345821
#> [1031]  0.1120760165  0.1381192385  0.0408517322  0.0021696871 -0.0023995383
#> [1036] -0.1459203316  0.0568634021  0.1237302233  0.0065325844  0.0215035763
#> [1041]  0.1435216908  0.2584118850 -0.0133380458  0.1296177848  0.0668930982
#> [1046]  0.0401633136  0.2399572386  0.0863803369  0.0391038810 -0.0753075157
#> [1051] -0.0605319671 -0.1206918159  0.2335391469 -0.0401117276  0.1880265608
#> [1056]  0.1013418147 -0.0623437395  0.2195179441  0.2438698309 -0.0250282958
#> [1061]  0.2601845447  0.1252681785  0.1093266410  0.1794260513 -0.1010926013
#> [1066] -0.1753655237  0.2099619150 -0.0375372493  0.0167586953 -0.0318267907
#> [1071] -0.0045423079  0.2110962878  0.1569422484 -0.1010205766  0.3815816943
#> [1076]  0.1297329551 -0.0537534717  0.0969732152 -0.0062478070  0.0645293469
#> [1081]  0.2786569718  0.0770106414 -0.0474116167  0.1026797229  0.0880282869
#> [1086]  0.1742209483  0.0558957226  0.0365226211  0.0873273321  0.0874791697
#> [1091]  0.2562689085  0.1182930305 -0.1759852998  0.0748024952  0.0963485460
#> [1096]  0.1743449569  0.1070046006  0.1201210317 -0.0981349610 -0.0580974650
#> [1101]  0.1030466013  0.1367686471  0.1448482697  0.1180385654  0.2917784481
#> [1106]  0.0989307190 -0.2303095208  0.0637103083  0.0579519635  0.1524057824
#> [1111]  0.1895755907  0.1355587997  0.2296558406 -0.1607456425  0.1301068463
#> [1116]  0.2769945186  0.0716339543 -0.0016034496  0.1823759467  0.0237125964
#> [1121]  0.2062136633  0.1143093600  0.1663494368 -0.0687569447  0.2379477664
#> [1126]  0.1704499459  0.2137008491  0.0860273884  0.2516170531  0.2212670200
#> [1131]  0.1750588531  0.1071437357  0.0165965252  0.0717543263  0.2873969292
#> [1136] -0.0430319854  0.0181659687  0.2475068668  0.1203906942  0.1428480971
#> [1141]  0.2174962587  0.0516297043  0.0737099702  0.0931483679 -0.1031023678
#> [1146]  0.1412308834 -0.0376152375  0.1313883264  0.1966692624 -0.0506334970
#> [1151] -0.0342575005  0.1859611195  0.0563604111  0.1671511624  0.0520043407
#> [1156]  0.1004184509  0.2147257437  0.1758517510 -0.0516099870  0.1377738038
#> [1161] -0.1220330174  0.0230313875  0.1498066319 -0.1437556479  0.0515479012
#> [1166]  0.2379253475  0.0332756969  0.2133262701  0.1457876692  0.0306823715
#> [1171]  0.0250997865  0.0455312981  0.1283323498  0.1906104692  0.3374221723
#> [1176]  0.3677516861  0.0651435375  0.1333967997  0.2549396914 -0.0102130121
#> [1181]  0.1804334148  0.0473474940 -0.1311524403 -0.1042507299 -0.2161911507
#> [1186]  0.1052944721  0.0873535509  0.0113011975 -0.0931200477 -0.0656291521
#> [1191] -0.0262851003 -0.0809183696 -0.0049683313  0.1847556947  0.0464249848
#> [1196] -0.1559019871  0.0636830587  0.0238693345  0.1315339575  0.0174074444
#> [1201]  0.2054822689  0.1047376914  0.2236016144  0.0529988471 -0.0753485809
#> [1206]  0.0834191570 -0.0584187869 -0.0429922765 -0.0283741651 -0.0609568050
#> [1211]  0.0187838618  0.1136241567 -0.0072341011 -0.0393920301  0.2179770239
#> [1216] -0.0029052078  0.1013956564 -0.0622689567  0.0231497568 -0.0181238830
#> [1221] -0.0976879361  0.1752277443  0.0971505676  0.0073827488  0.0617649368
#> [1226] -0.0574279860  0.1082921184  0.2120631011  0.1128685672  0.2434835179
#> [1231] -0.0232347534  0.1235273592 -0.1295770205 -0.0789177510  0.1752570248
#> [1236]  0.1208444703  0.1007038569  0.1443354468 -0.1071455162  0.1663570755
#> [1241] -0.0431227870  0.3012996243  0.1562216991  0.1105962153 -0.0070379804
#> [1246]  0.2933478645  0.0472262103  0.0931910436  0.1736886249  0.2386334704
#> [1251]  0.1837055571 -0.0018490876  0.0393602431 -0.0398247180  0.1949923609
#> [1256]  0.2191489793  0.0570242620  0.0104088919  0.0026490265  0.3115417806
#> [1261]  0.1575392285  0.1526910617  0.1493510856  0.0024962078  0.1782922051
#> [1266]  0.2111779353  0.2435763932 -0.0368388907 -0.0428671758  0.2514141917
#> [1271]  0.2295699590  0.1276656349  0.2435286380  0.2544346438  0.2280878045
#> [1276]  0.0779829706  0.0854738231 -0.0021708182  0.0442048652 -0.0030729764
#> [1281]  0.0600873711  0.0575305314 -0.0227442576  0.1123424465 -0.0712944831
#> [1286] -0.0313621066 -0.0320611635  0.0099480197  0.1440666813  0.0958947376
#> [1291]  0.0318986437 -0.0432957337  0.1946691306  0.1011438932  0.0783271481
#> [1296]  0.2433786730 -0.0869465095  0.0412479580  0.1169690733 -0.0498671523
#> [1301]  0.1972106443  0.0466656506  0.0526612684  0.0588717761  0.0384010471
#> [1306]  0.1748121893 -0.0219955044 -0.0821772998  0.2267715850  0.1311082876
#> [1311]  0.2205327716  0.1373124233  0.1312294924 -0.0739777276  0.0006082841
#> [1316]  0.0562103500  0.1688800530  0.0598385681 -0.0148927833  0.0162189708
#> [1321]  0.3289666653  0.2368619053  0.0540048400  0.2386849196  0.0078548953
#> [1326]  0.2739567876 -0.0384496238  0.0642775931  0.1226749448  0.1722144944
#> [1331]  0.0729083187 -0.0773101014 -0.1339228039  0.1019622365  0.2112577119
#> [1336] -0.1038053844 -0.2293296536 -0.0717014598  0.3779868631  0.1975916449
#> [1341]  0.1514479203 -0.0503660010 -0.0039888967  0.0852506376 -0.0053651040
#> [1346] -0.0142125474  0.1593613714 -0.0272719032  0.0393375262  0.0887458321
#> [1351] -0.0040541252  0.1708705561  0.0192971981 -0.1503801169  0.4127022368
#> [1356]  0.0616833216  0.2807308997 -0.0092666771 -0.0826775191 -0.1532632914
#> [1361] -0.0115117907 -0.1275813021 -0.1765104355  0.0219985052 -0.1147527486
#> [1366]  0.0030427079 -0.0395345648 -0.0016974368  0.3866515221  0.0421600490
#> [1371]  0.0804795852  0.0513359538  0.0605114972  0.0656774249  0.0681601559
#> [1376] -0.1611500542  0.1557488567  0.1910922549 -0.2119614957  0.1267200264
#> [1381]  0.1242227006  0.0363955583  0.2698006032 -0.0154530145  0.0102401798
#> [1386]  0.1568989219  0.0382049477  0.1668185716  0.1242184209  0.1960240373
#> [1391]  0.0811722794 -0.1018891805  0.0263728438 -0.1497557406 -0.0382726773
#> [1396]  0.1653797117 -0.1449245518  0.0736514097  0.1765153591  0.0220400905
#> [1401]  0.0907913295  0.1039539906  0.0104322222  0.2265412222  0.1344007152
#> [1406]  0.1968119788  0.0514408601  0.1790156809 -0.2403513605  0.1426132293
#> [1411] -0.0688891762  0.0632015203  0.0747642697 -0.1106048302  0.0305193007
#> [1416] -0.0259898627  0.1451450430  0.1830573582  0.3077155495  0.0383097327
#> [1421]  0.0496243312  0.0729391402  0.0897223644  0.1777016368  0.3148011328
#> [1426]  0.0895192100  0.3486657891  0.2864797811  0.0662080389 -0.1534718125
#> [1431]  0.0131526016  0.0203199690  0.0851921321  0.0055697213  0.0208891844
#> [1436] -0.0672466093 -0.1325713489 -0.1239684112  0.1436981243  0.0597410632
#> [1441] -0.1275869097  0.1201038995 -0.2231956978  0.0019555889  0.1422206632
#> [1446] -0.1311236658  0.2038671121 -0.1724608088 -0.0234443695  0.1144726375
#> [1451]  0.0724422396  0.1663416169  0.0940107139  0.2559935506 -0.1594159291
#> [1456]  0.0328784463  0.0414939017  0.0242470434  0.0890457456  0.2424252240
#> [1461]  0.2536023803  0.1504827812 -0.3126485971 -0.0994618249 -0.0224765310
#> [1466]  0.2070061977 -0.0362788350  0.1337006240  0.0156825318  0.1929178495
#> [1471] -0.0202581935  0.0463744111  0.1282683934  0.4139937022  0.1069906925
#> [1476]  0.1501234183 -0.0454926609 -0.0780839958  0.1503182409  0.1183617156
#> [1481]  0.0784750289  0.0651281101  0.1607252891  0.2239222184  0.0868924146
#> [1486]  0.0248617590  0.1353972975 -0.0192986822 -0.0298981234  0.0503831433
#> [1491] -0.0555611896 -0.0505457257  0.2061806414 -0.0286739824  0.1961782121
#> [1496]  0.0957024261  0.1754449627  0.1206707502  0.1587894561 -0.0413897752
#> [1501]  0.0742712864 -0.1292101773  0.1386317018  0.0001278228  0.1914558308
#> [1506]  0.1386271486  0.0139443257 -0.0058596006 -0.0728207065  0.0911225298
#> [1511]  0.1564701765 -0.0808793493  0.0932206493 -0.0469566088  0.0239688417
#> [1516] -0.0779165893 -0.0801774295  0.1364238964  0.1446272888  0.0978533511
#> [1521]  0.0064844820  0.2235511747  0.2125067222  0.0781241690  0.0411336309
#> [1526]  0.2529594770 -0.1083831845 -0.0582104700 -0.0928847854 -0.0097618704
#> [1531]  0.1233059946  0.1383947028  0.0575227482  0.2093231148 -0.1242876808
#> [1536]  0.1475975082  0.2312967072  0.0615692929  0.0738997808  0.1433803714
#> [1541]  0.2344882529  0.1560739531  0.1279322353  0.2786140904 -0.0356518593
#> [1546] -0.1409200481 -0.0278848618  0.3117579084 -0.0658557055 -0.0204682180
#> [1551] -0.1601678745  0.0377607612 -0.1743569917  0.2038612297  0.1523240888
#> [1556]  0.1336017892  0.1770840066  0.2265575293 -0.0030707572 -0.1068610622
#> [1561] -0.1784263070  0.2325880249  0.0601403910  0.2057343782  0.0239643761
#> [1566]  0.0699533036 -0.2028098048  0.3610296895 -0.2182547928 -0.0318660404
#> [1571]  0.0064907286 -0.0742772434 -0.3208369368  0.3395662084 -0.0129620822
#> [1576]  0.0897462022  0.0158253290  0.1353461324  0.1507534477 -0.1930094177
#> [1581]  0.0596263690  0.2401026051  0.1520699253  0.1472333647  0.0923798339
#> [1586]  0.1255376468  0.0210903838  0.1039385184 -0.0015914389  0.1541843835
#> [1591]  0.1663749293  0.1955942295  0.0632927505 -0.0450119386  0.0636299223
#> [1596]  0.1200517534 -0.0444583456 -0.0981502775  0.2568787482 -0.0174467151
#> [1601]  0.2512215605  0.0516843840 -0.0991956806 -0.0048949849  0.0457979761
#> [1606] -0.1630434068  0.1464340693  0.0602927644  0.0160348331 -0.0426450129
#> [1611] -0.1489326491 -0.0041016762  0.0904691230  0.0146631387  0.0651352965
#> [1616] -0.0289478341 -0.0741986149  0.0806000459  0.2002703355 -0.1062085614
#> [1621] -0.0834769547 -0.0524475881 -0.0200032840  0.2134559278 -0.0405961288
#> [1626]  0.1551095175  0.2296189574  0.1094896270  0.0127727210  0.2262609449
#> [1631]  0.1163432393  0.0056734390  0.0529647536 -0.1722346301  0.0938557997
#> [1636]  0.0924677479 -0.0136925165  0.1780511383 -0.0390268805 -0.2670082747
#> [1641] -0.0430966068  0.1711336489  0.1178464887 -0.0837170026 -0.2066825045
#> [1646]  0.1058116007  0.1431873626  0.2203243632 -0.0480128828 -0.0628687291
#> [1651] -0.0400556254 -0.0040615772  0.0148147878  0.0457012835  0.0552720450
#> [1656]  0.0504701688  0.2301073792  0.0044267131  0.0899293290 -0.0650768437
#> [1661]  0.1983754328 -0.0664358276  0.0248512775  0.0706436009  0.0065175768
#> [1666]  0.3152493104 -0.0335655791 -0.1680917796 -0.0787582203 -0.2242529614
#> [1671] -0.0099374022  0.0083635682  0.0162479884 -0.0384273484 -0.0832663179
#> [1676]  0.1587042267 -0.2221429976 -0.0146775701  0.1098405108  0.1681252116
#> [1681]  0.1333610022  0.1292802261  0.4184456197  0.0367351104  0.1177805898
#> [1686] -0.0341210336  0.0287792135  0.0943250250  0.1799382872 -0.2362474611
#> [1691]  0.0290211099  0.1292468444  0.0469280458  0.1174670844 -0.0404147398
#> [1696] -0.1453346130 -0.0629738491  0.1036917163  0.1534946713  0.0156984600
#> [1701]  0.1843967030  0.1293626973  0.1231351380 -0.1268351680  0.1062422437
#> [1706]  0.0902164448  0.0722743515  0.0490221478 -0.0121033596  0.0568132132
#> [1711]  0.0155653831 -0.0052640019  0.0991898800  0.1458580201  0.0836719547
#> [1716]  0.1691023822  0.1343609090  0.1313614945 -0.0035067481 -0.1281517707
#> [1721] -0.0405163741 -0.0176443626  0.1525840747  0.2884497723 -0.0214945077
#> [1726]  0.0535925282  0.0310401318 -0.1255065916  0.2781433799  0.1928962500
#> [1731]  0.0739164849  0.0815436286  0.0448869569  0.0999296319  0.0416031887
#> [1736]  0.1927315565  0.0859360691  0.0367435994  0.1515673775  0.1558318931
#> [1741]  0.2861680140  0.2637937411 -0.0346715782  0.0807249124  0.0257607302
#> [1746]  0.2492180949  0.1656912063 -0.1352285303 -0.0149704292  0.1867776736
#> [1751]  0.1161530114  0.2529817999  0.0519972260  0.0897603794  0.1016732166
#> [1756]  0.2063588806  0.2120897749  0.0576062298  0.0251152093  0.1505548775
#> [1761]  0.2972597669  0.3280435996  0.3653716176  0.1597178402  0.0150977464
#> [1766]  0.0562308672  0.2309710456  0.2072852967  0.1790027602  0.2568496898
#> [1771] -0.0178493075  0.0271431331  0.0442132767  0.1193034277 -0.0681203599
#> [1776] -0.0265876163 -0.1846776507  0.2068964243  0.1862942185  0.0592008468
#> [1781]  0.2327766769  0.0344402824 -0.2303271462  0.1938600010  0.1838609148
#> [1786] -0.1132809077  0.0911198444 -0.0497408951  0.0422130488  0.0238944025
#> [1791]  0.0965312194  0.0029284539 -0.0798027750  0.0719268919  0.1518110776
#> [1796] -0.0563404330 -0.0885843161 -0.2150216863 -0.0960735254  0.0372046145
#> [1801] -0.1267129897  0.0111377902  0.1530925294  0.0584894172 -0.1075358719
#> [1806] -0.0431082199 -0.0364374111  0.1690985919  0.0817238111 -0.2385853751
#> [1811]  0.0345174067  0.2514306919  0.0355972316  0.2086134100  0.0764727749
#> [1816]  0.0031783754 -0.1387058892  0.2140624470  0.0866054931 -0.0297202055
#> [1821]  0.0747523352 -0.0251741531 -0.0181042322  0.1895497403  0.1380444514
#> [1826] -0.0147040334 -0.0515527287 -0.2908313536 -0.0798825601 -0.2226073113
#> [1831]  0.0713737802 -0.0414672674 -0.0915434871  0.2698732698 -0.0015050942
#> [1836]  0.1615983619  0.0748352159  0.0832139038  0.0641376984  0.0688864302
#> [1841] -0.0077955880  0.0205387865  0.1044440007  0.3171719167  0.0571842340
#> [1846] -0.0753779125  0.1550445403 -0.0024118030  0.2329787140 -0.0944008279
#> [1851]  0.1055903686  0.0924261922 -0.0649290882  0.2211656061  0.1520640275
#> [1856]  0.1284149064 -0.0007029981 -0.0053798639 -0.0014000805 -0.0700873938
#> [1861]  0.0586561124  0.0133841140  0.0842115578  0.0859330005  0.2416255009
#> [1866]  0.1578822895  0.0148476904  0.1640337297 -0.0038912149  0.2251520272
#> [1871]  0.1977747163  0.0626030032  0.1338735276  0.1938700991  0.2592485051
#> [1876]  0.2740518337  0.0119306430  0.3271444377  0.0459240442 -0.0247129564
#> [1881]  0.2168168224  0.1484038643  0.2166324097  0.0508507653  0.1189368630
#> [1886]  0.2662860802  0.2118274687 -0.1957560439  0.0231089172  0.2867183170
#> [1891]  0.1042467503  0.0310497675  0.1629967209  0.2296448073  0.1350345257
#> [1896]  0.1305082925  0.3173042382  0.0449202196 -0.1557049081  0.1659152522
#> [1901]  0.0777631688  0.1030702551  0.0608085389  0.1920799193  0.0831960930
#> [1906]  0.0931544264  0.0649652989  0.0592764604  0.1582449122 -0.0083626761
#> [1911]  0.1417811452  0.1737961425 -0.1252746176  0.0631416429 -0.0342302563
#> [1916] -0.0801379260  0.1267744800 -0.0427102958 -0.0241221749  0.0599252748
#> [1921] -0.1410931018  0.1519879946  0.1685189611  0.0052623001  0.0630258799
#> [1926] -0.1368668918 -0.0249794430  0.0625357157  0.1497462895  0.0939344222
#> [1931] -0.1151864458 -0.0410244634 -0.0400936483 -0.0289761146 -0.1350423709
#> [1936]  0.0418538684  0.0163988752  0.0790053731  0.1321692633  0.1357794206
#> [1941] -0.0508335864  0.2489761504  0.2597649120  0.2756304470 -0.1509233905
#> [1946] -0.0689645492  0.2731430328 -0.0285283251  0.1481160742 -0.0188923242
#> [1951] -0.0274314324  0.0413767258 -0.0279574457  0.0368939847  0.0816323083
#> [1956]  0.1108827253 -0.0715513958  0.0679376698  0.0656491896 -0.1798571817
#> [1961]  0.1957493394  0.1516229902  0.1084308009 -0.0639140056  0.0070601972
#> [1966]  0.3216479491  0.1013069066  0.0955475855 -0.0697447339 -0.0843266692
#> [1971]  0.1474504165 -0.0315188660  0.1692394399  0.3088420264  0.1447381930
#> [1976] -0.0544352380  0.0842197158  0.3120073789  0.2553236511  0.0138769656
#> [1981]  0.0602748535  0.1230186683  0.1697488361  0.1705953043 -0.1411898934
#> [1986] -0.0361347507 -0.0020763159 -0.0146389328 -0.1948669655 -0.0500958347
#> [1991]  0.2569612874  0.2174787714 -0.0010880261  0.0721447958  0.0503499096
#> [1996]  0.1826029069  0.1081008845  0.1368670448  0.3354979668  0.2188716532
#> [2001]  0.1489501069  0.1878885979 -0.0000282232  0.1536590200  0.3385031448
#> [2006]  0.3269641199  0.1095720374  0.1675084887  0.2996783230  0.1227690457
#> [2011]  0.0959476928 -0.0080158971  0.2103574906  0.0655455915  0.0976203345
#> [2016]  0.1173971754  0.0742385568  0.2133342579  0.1477705850  0.0249627192
#> [2021]  0.2152890087  0.1083369497  0.0149848217  0.0806438055  0.1929662698
#> [2026]  0.0737537067 -0.0274593896 -0.1449524352  0.0782856840  0.0450963770
#> [2031]  0.1207558997 -0.0727651060  0.3682596472  0.0104632279 -0.0444902581
#> [2036]  0.0525827819  0.1427401860 -0.0846656659  0.0134663307  0.0601727696
#> [2041]  0.2999644080  0.1015802853  0.1817158109  0.0491237044  0.0417958871
#> [2046]  0.0383447468  0.1717883552 -0.1293737704 -0.1595797544  0.0188282079
#> [2051]  0.1999153213 -0.1308992943 -0.0660157762  0.1305594124 -0.0394181817
#> [2056]  0.1773328318  0.0925207864  0.1814766027  0.1502728620 -0.0141373249
#> [2061]  0.1054140493  0.1681361162 -0.0667504636  0.1459671023  0.1575731843
#> [2066]  0.0670888370  0.1466973398  0.1774000270  0.2794721617  0.1891031604
#> [2071]  0.2356976187 -0.1729757675  0.0279627553 -0.0001556956 -0.0991153724
#> [2076]  0.0046356369  0.1372952136  0.1408163664  0.0561401094  0.1749174115
#> [2081] -0.1432731844  0.0751826624  0.0804377594 -0.0655214834  0.1311849572
#> [2086]  0.1722113665  0.1318383088 -0.1017653869  0.0400945503  0.2931645333
#> [2091]  0.1008660868  0.0140965178  0.0225483085 -0.0104900954  0.0250481997
#> [2096] -0.2082530696 -0.0251988458 -0.0969937840  0.1524826202  0.0725031061
#> [2101]  0.2611702001 -0.0013611957  0.1956934143 -0.0716280459  0.1213886198
#> [2106]  0.0061075484  0.2815531249  0.0563706035  0.1312841515  0.0821802530
#> [2111]  0.2337845478  0.0241568924 -0.1748133790 -0.0266331132  0.0409396333
#> [2116]  0.1028312205 -0.1773626088 -0.2220250939  0.1092091098  0.0218096651
#> [2121] -0.0532625317 -0.1410375345  0.1643653882 -0.1098659761  0.2718583789
#> [2126] -0.0031983553  0.0666141180  0.2121330911  0.2153008727  0.2703490304
#> [2131]  0.1453641872 -0.0713822980  0.0998168233  0.2006257809  0.0935726709
#> [2136]  0.0173249497 -0.0777848542 -0.1361835851  0.1034920546  0.1461643591
#> [2141]  0.0911711701 -0.0080426384  0.1768854947 -0.0327938088  0.0133346585
#> [2146] -0.1012544937  0.0517349865 -0.0139352693 -0.0608430423  0.2749573770
#> [2151]  0.0169850837  0.0656346928  0.0771399215  0.0621218079 -0.0676756921
#> [2156] -0.0902767429  0.0738168896  0.1148482239 -0.1032592944  0.0310269707
#> [2161] -0.0845557607  0.1276775684 -0.0111188436  0.0108342845  0.1387630756
#> [2166]  0.1839519163  0.1455635461  0.1116524146  0.2698090584  0.0540655443
#> [2171]  0.0262318959  0.0655887842  0.1180182956  0.0262926865  0.0377469545
#> [2176]  0.0271399599 -0.0155163429 -0.1094247636 -0.0796340300  0.1940671791
#> [2181]  0.2592001500  0.0416020586  0.0270163722  0.0877647976  0.2456711504
#> [2186]  0.2448985379  0.1534739746  0.0258553124  0.2274097465  0.0128257713
#> [2191]  0.2270600327  0.2574466694  0.2578062578  0.0587697341  0.0469399875
#> [2196]  0.2849685744  0.2249302448  0.0775750503  0.1592979331  0.2266166081
#> [2201]  0.1108809689  0.1744646871 -0.0669660152 -0.0361744772  0.1713807588
#> [2206]  0.0050835373  0.2471695589 -0.1487521808  0.1135141505  0.1121906870
#> [2211]  0.0701267660  0.1121568954  0.0326114964  0.1305377486 -0.0938500665
#> [2216]  0.3225597494  0.2832449870  0.1040855664  0.0730369556  0.0839106675
#> [2221]  0.3736813089  0.1157156764  0.0698860016  0.1848899086  0.0712086731
#> [2226] -0.0464871015 -0.0236644882  0.0806554164 -0.0048618604  0.1038548269
#> [2231]  0.0115713690  0.0816248963  0.1259213002  0.1732571561 -0.0603448712
#> [2236]  0.2699202855 -0.1847210402  0.1516407212  0.1756983517  0.0852558086
#> [2241]  0.0598529089 -0.0884355818  0.0650905579  0.1087177201  0.0075506140
#> [2246]  0.2132824883  0.1130449424 -0.0439019581 -0.0792018039  0.2695424607
#> [2251]  0.0811942722  0.0959483482 -0.1325769372  0.2747324727  0.3308491059
#> [2256]  0.0533826698  0.1315002820  0.0020466650  0.3442473733  0.1022475337
#> [2261]  0.1459571956  0.1438644060 -0.0427407074  0.1167508528  0.1296779961
#> [2266]  0.1514819347  0.1296892595  0.1103213897 -0.1615979281  0.1345546361
#> [2271]  0.1965833773  0.1084037906  0.3160545408  0.0409193677  0.1200350003
#> [2276]  0.1170139344  0.2210387542  0.3314915792  0.1231822557  0.0214340580
#> [2281]  0.1310529915 -0.1426243957  0.0579392173  0.1494026215  0.0564071103
#> [2286] -0.0875412602  0.0268970661 -0.0113483161 -0.0302140572  0.0598685898
#> [2291]  0.1306539635 -0.0682452011  0.1811853067  0.0151148547 -0.1792307032
#> [2296]  0.0643194169  0.0848662433 -0.0299010675  0.1798560697  0.1418511551
#> [2301]  0.1704219469 -0.0198710511  0.1614911732  0.0915744861  0.0348379158
#> [2306]  0.0145319891 -0.2028062639 -0.0447765676  0.0042170441  0.1253114398
#> [2311]  0.1053765909 -0.1096340761 -0.1378201632  0.1132326266  0.1499234129
#> [2316]  0.4544210162  0.0943149602 -0.0407547602  0.0630869492  0.0950815388
#> [2321] -0.0149247471  0.0015054946 -0.1250743548  0.0639442943  0.1181249964
#> [2326]  0.0978597561  0.0440591205 -0.0071526685  0.1151464574  0.0610187747
#> [2331]  0.0424198725  0.0305582274  0.1188543360 -0.0521286239 -0.0882884084
#> [2336] -0.0850003050  0.0776165269 -0.0371318467  0.2040328444 -0.0423400242
#> [2341]  0.1459230261  0.0723860993  0.0017891384  0.0172789853  0.3088218120
#> [2346]  0.0681969208  0.2103377103  0.1473113614  0.1103404083  0.2451687341
#> [2351]  0.1305308367  0.1299039290  0.2247023219  0.1418152745 -0.1796093036
#> [2356] -0.0402920641  0.0398049973  0.1058900452  0.1305525489  0.0637506772
#> [2361] -0.0037717014 -0.0597341053  0.0292669347  0.0519869152  0.0884782484
#> [2366] -0.0158691970  0.0640862396 -0.1007748418 -0.0962338049  0.1634399916
#> [2371]  0.0628898608  0.0878881428  0.1280694815  0.2489945349 -0.2092204061
#> [2376] -0.2584709447  0.1266881089  0.1267190536 -0.0171527283  0.1234947931
#> [2381]  0.0730178381 -0.0005104646 -0.0019097549  0.2531893579 -0.1904259003
#> [2386]  0.1082906485  0.2280967080  0.0958999375 -0.0363197056 -0.0446018090
#> [2391] -0.0033980448  0.1977605508  0.1555455005  0.1623222927  0.2032902955
#> [2396]  0.0587077097  0.1733574598 -0.0026592105  0.0145044925  0.1401959881
#> [2401] -0.0894931714  0.0514509745  0.2071784195 -0.0453281126  0.0549585019
#> [2406]  0.1025085202  0.1036941298  0.1612605911  0.0344034768  0.0492247599
#> [2411]  0.2464951817  0.1389998399  0.3134703644  0.0011598030  0.3570525170
#> [2416]  0.1582648452  0.1671397870 -0.1345299490  0.2298636390 -0.0516055918
#> [2421]  0.1476103171 -0.2064553692  0.2306333473  0.0688691239 -0.0852736127
#> [2426] -0.1758069520  0.1566828014  0.2283897045  0.2248920377 -0.1949445961
#> [2431] -0.1430059597 -0.1024910359  0.1322687118  0.0157404825 -0.0372833645
#> [2436] -0.0962717759 -0.1379219972 -0.0324966143 -0.0263434073  0.1844951783
#> [2441]  0.1311948214 -0.1555569870 -0.0932371188  0.2041776957  0.1953299873
#> [2446] -0.1257243359 -0.1850206673  0.1294426717  0.0289215034 -0.0196360052
#> [2451] -0.0734818968 -0.0144798630  0.1575712336  0.0797196593  0.1320388786
#> [2456]  0.1245266779  0.0481698855  0.0508394596  0.2274280194 -0.0496127673
#> [2461]  0.1058358598 -0.0248704460  0.1019380945  0.0981258656  0.0767538187
#> [2466]  0.0399507447  0.1178738355  0.0693058050  0.0511219748  0.1149432145
#> [2471]  0.0118483076  0.1321338262 -0.1736014499  0.4348794992  0.3014403232
#> [2476]  0.2557402175  0.3079530925 -0.0693266221  0.0949685229  0.0056626037
#> [2481] -0.1843542547  0.0797956212  0.2018589349  0.0367842550  0.0920257384
#> [2486]  0.2764056081  0.0184415120  0.2999198546 -0.1138189551 -0.0565684634
#> [2491]  0.0742050020  0.0339679464 -0.0323207781 -0.0034593320  0.0292786682
#> [2496] -0.0388283802  0.0377536852  0.1214300288 -0.0558199283  0.0840167107
#> [2501]  0.0293019244 -0.0798142587 -0.0263557014  0.0263428494 -0.0862794832
#> [2506] -0.0230703495  0.1164464950  0.1238311825  0.0456368811  0.2946443759
#> [2511]  0.0092209022  0.0956361505 -0.1575892108  0.0329522703  0.1914699807
#> [2516]  0.1316544083  0.0449408626  0.1311180550 -0.1153708578  0.0353804902
#> [2521]  0.1451171932  0.1420302178  0.2355483553  0.1859173770  0.2535521980
#> [2526]  0.0773466258  0.0685472527 -0.1515616280 -0.0360248821 -0.0471586217
#> [2531]  0.2221425700  0.1128141377  0.2117907154  0.1473318533  0.2154457461
#> [2536]  0.1535006287  0.1348170050 -0.0845992433 -0.1144695880 -0.1299605454
#> [2541]  0.0519246368  0.0360870793  0.1740967066  0.2416189949  0.3227031666
#> [2546] -0.1851295838  0.1292041377 -0.0756932948  0.0388716941 -0.0291096996
#> [2551]  0.0460328224  0.1187526793  0.0846854291  0.2456012069  0.0668499854
#> [2556]  0.2484102642 -0.1019885412 -0.1187501661  0.0369790687 -0.2928698636
#> [2561]  0.0875976107 -0.0163554046  0.3104575126 -0.0799978629  0.1085088937
#> [2566]  0.1655340707  0.0241285500 -0.1303258541 -0.0043736115 -0.0215422861
#> [2571]  0.0194533073  0.0866442755  0.0029071504  0.2123765171 -0.0754022831
#> [2576] -0.0302193089  0.1200956107  0.0469636725 -0.1246582058  0.0805927249
#> [2581]  0.0114021597 -0.0472137076  0.0800075913 -0.0935564039 -0.0400604107
#> [2586]  0.0147580064 -0.0846692164  0.2765264093  0.0251500859  0.0574205029
#> [2591]  0.1021637780 -0.0346639215  0.1926627229  0.1720765111 -0.1925615246
#> [2596]  0.1349365253  0.1344953230 -0.0515522574  0.0420211230 -0.0879520709
#> [2601]  0.1174250193  0.2060766580  0.2289964614  0.1279524231  0.0017914404
#> [2606] -0.0990844395  0.0132450332  0.0540010625  0.1137710179  0.0510716004
#> [2611] -0.0592260385  0.1453507105 -0.0934031399 -0.1286689883 -0.1131706316
#> [2616]  0.1379340381  0.0091552286  0.2839503810 -0.0521012768  0.1765791654
#> [2621]  0.0480019606 -0.0659254266  0.0989906827  0.2077857242  0.1151583758
#> [2626]  0.0276845305 -0.0893208722 -0.0652345321  0.1926146461  0.0811715022
#> [2631]  0.2422265760  0.1415531979 -0.2945941781  0.1322654317  0.1541980407
#> [2636]  0.1868878513  0.0195744667 -0.0877424533  0.0851509907  0.0030800667
#> [2641]  0.1391694567 -0.0126001455 -0.0630463430  0.2186184569 -0.1337041567
#> [2646]  0.2590101915  0.1076235160  0.0419031132  0.0651417069  0.0581224111
#> [2651]  0.0225177384 -0.0341950991  0.2271931402  0.1870593227  0.1836898605
#> [2656] -0.0046633299  0.2243833549  0.1793593216  0.1548616269 -0.0214316566
#> [2661]  0.2211216078 -0.0542846710 -0.0941933261 -0.0433199638  0.2428463198
#> [2666]  0.5007269361  0.0167011502  0.1017411281  0.0342126042 -0.0172407355
#> [2671]  0.1074323283  0.0459359008 -0.0424514811  0.0634777534  0.2700303580
#> [2676]  0.1321566748  0.4657557116  0.0631020343  0.0890673527  0.1582523130
#> [2681]  0.0407078585  0.0801284920  0.1108997268  0.0597702304 -0.1128057153
#> [2686]  0.0892166790  0.0706351816  0.1548207641 -0.0395765256 -0.0338814228
#> [2691]  0.0371944608  0.2155427661  0.1246281702 -0.0389956077 -0.1353902516
#> [2696] -0.0933703633  0.2232263064  0.1604567975  0.1151810169  0.3478139323
#> [2701]  0.0240001704  0.0244702412 -0.1187333456  0.1217360750 -0.0342543421
#> [2706]  0.5419132636  0.0893073289  0.2565223750 -0.1114630733 -0.1975111663
#> [2711]  0.1677942042  0.1771382341 -0.1630631777  0.0964486592 -0.1552504289
#> [2716]  0.1268669031  0.1306447601 -0.0414904403  0.1079351463 -0.0345553808
#> [2721] -0.1178850070  0.2909600746  0.0961928694  0.0818300742  0.0262979137
#> [2726] -0.0159178005 -0.1415442940  0.1041817278  0.1290077218 -0.0602281319
#> [2731]  0.0299619110  0.0978114797  0.2814921919  0.0960285723  0.1152488792
#> [2736]  0.3138515470  0.0237204911  0.0613429847 -0.1098773851  0.1917087257
#> [2741] -0.0438176164  0.1032181544  0.0651375828  0.0620103930  0.2145028286
#> [2746]  0.0903617116  0.2777410412 -0.0761815272  0.1073417839  0.0659717701
#> [2751]  0.2457601311  0.2600612178 -0.0596057116 -0.1303239426  0.2954750030
#> [2756]  0.0912784389 -0.0305560509 -0.1285334392  0.0829935998 -0.0860630820
#> [2761]  0.1928595887 -0.0497615655  0.0946866048  0.0964953363  0.0381653964
#> [2766]  0.0882393654  0.0972282039  0.1059568442  0.2376919759  0.1346530076
#> [2771]  0.0725564593 -0.1120556992  0.0291266431  0.1011746239 -0.1126235794
#> [2776] -0.1679982087  0.1366222234  0.3318603623  0.0738915748  0.1085145513
#> [2781] -0.1285737017 -0.0061849015  0.0198182815 -0.0097392404  0.0041800689
#> [2786] -0.0177100068 -0.1294002610 -0.0601909627 -0.1061135136  0.3303331791
#> [2791] -0.1152913696  0.2958308688  0.0624442097  0.0849693950 -0.1158097118
#> [2796]  0.1582903513 -0.0453953801  0.1125319040  0.0772074874  0.1005187905
#> [2801]  0.0069061236  0.1470974526  0.0698195252 -0.0059671990  0.0401315437
#> [2806]  0.1130245292  0.0588670604 -0.1900482163  0.0297832536 -0.1473158007
#> [2811]  0.0236529577  0.1048809655  0.1789552093 -0.0818675348 -0.0321208319
#> [2816]  0.2121725476  0.2509178049  0.1726549234  0.1873086094  0.1208061724
#> [2821]  0.0135567303  0.1478060760  0.1527463606 -0.0731857411  0.0918451379
#> [2826]  0.2782223899  0.1746640783  0.2208303862  0.1200809381  0.2324094567
#> [2831]  0.0117915915  0.2459890408  0.0050636879 -0.0419212810 -0.0727875908
#> [2836]  0.2997409118  0.0590946799  0.1404051872  0.1621678016  0.1350750429
#> [2841]  0.1409449178  0.2406391792  0.1476821292  0.1828206055  0.0130313964
#> [2846] -0.0828797776 -0.1936221517  0.1511350164 -0.1055194583 -0.0107746545
#> [2851]  0.0233659684  0.0694280323 -0.0134108035  0.0935083943 -0.0337455865
#> [2856]  0.0018510944 -0.0179650916  0.1327672222 -0.0081131489 -0.0401437363
#> [2861] -0.2412375465 -0.1242790223  0.0762228902 -0.0665131407  0.1961060657
#> [2866] -0.2060574993  0.3669540426 -0.0216761550  0.0042891309 -0.1062722798
#> [2871] -0.2566734011  0.0331200178  0.0315967102  0.0321001586  0.2729841325
#> [2876]  0.0973732729  0.0535203748 -0.0532228411 -0.0759935945  0.1906624750
#> [2881]  0.2819138012 -0.1284546223  0.0146677231  0.1077367804 -0.0201980319
#> [2886]  0.0465564562  0.0650713949  0.3744158495 -0.1430707250  0.1186481082
#> [2891]  0.2523483649  0.0353764138 -0.1136846998 -0.2277931526  0.0736849661
#> [2896]  0.0280301349 -0.1272753161  0.1333208465 -0.0750050490 -0.0074481390
#> [2901]  0.1493422607  0.0829947585  0.1086600616 -0.1625320880 -0.0142184544
#> [2906]  0.0825101525  0.0723909800  0.1071200590  0.1610559351  0.1441435717
#> [2911]  0.0522195547  0.2564237042 -0.1276122614  0.1390988590  0.3431277080
#> [2916] -0.1739309445  0.0439405624  0.0077305052  0.0953469258  0.1504179541
#> [2921]  0.2086656223  0.2751088610  0.0876680502  0.1461042340 -0.2673676432
#> [2926]  0.1871503789  0.0761625602  0.1758340632  0.1240396190  0.0230331227
#> [2931]  0.3090392422  0.0706134130 -0.0879011220 -0.0816020508  0.1176421989
#> [2936] -0.0496486063  0.1391725633  0.1003018408  0.1892910422  0.1774818801
#> [2941]  0.1069721796 -0.2433636366 -0.0061155083  0.1764005806 -0.0361281883
#> [2946]  0.0928028563  0.3730203206 -0.1655843809 -0.1460636673  0.1435111657
#> [2951]  0.0547004346  0.2047974839 -0.0093785823 -0.1121970321  0.0397228070
#> [2956] -0.0833586140  0.1070391680  0.2162484511  0.0735880859 -0.0356930244
#> [2961]  0.0680075258  0.0261177638  0.3128905320  0.2582346313 -0.1009435580
#> [2966]  0.1287614808  0.1467987315 -0.0443502673 -0.0427175095  0.0052541920
#> [2971] -0.2123109700  0.2668965000  0.2292064177  0.0752984704  0.1475369719
#> [2976]  0.1094730505  0.1699003363  0.3437435180  0.2197885952  0.1501977539
#> [2981]  0.0954237507  0.1095961258 -0.2307242066  0.0606623871 -0.0799792933
#> [2986] -0.0492901409  0.1793554101  0.0259700836  0.2847133869  0.2886108308
#> [2991]  0.1926126970 -0.0588304896  0.2262864515 -0.0257270076  0.2259397214
#> [2996]  0.1223963340  0.0313512614 -0.0333349716  0.1034238948  0.0188480146
#> [3001]  0.1583553031  0.2932415008  0.0469093490  0.1325295946 -0.0397556278
#> [3006]  0.2194972695  0.0238732878  0.0836598002  0.1353763376  0.0802725852
#> [3011] -0.0289289305  0.0361254220  0.0605661687  0.1091621749 -0.0069985772
#> [3016] -0.0985514094 -0.0199349313  0.0116537384  0.1029462533  0.0324424176
#> [3021]  0.1733474715  0.1624725222 -0.1793186859 -0.0391517584  0.2729973194
#> [3026]  0.2812300378 -0.1114985921  0.1450038314 -0.1566047433  0.3274845341
#> [3031]  0.0179991526 -0.1330269094  0.2104867853 -0.0309943005  0.0586269188
#> [3036] -0.0109295360 -0.0673795856  0.0221183475  0.3028638509  0.0186951868
#> [3041]  0.0682157064  0.0365863555  0.0145893201  0.1612859894  0.0352701061
#> [3046]  0.0041973095  0.1043068481  0.0579908339  0.0990110808  0.1733771423
#> [3051]  0.0583846315 -0.0804107973  0.3282174940  0.0667775180  0.0880132404
#> [3056] -0.2089552386 -0.1653720433  0.0089328910 -0.1220701472  0.1013037978
#> [3061]  0.0791833340  0.1997410758  0.1968714506  0.1127539177  0.1648847444
#> [3066] -0.1306714461 -0.0879706971 -0.1225004711  0.1048504165  0.1716771169
#> [3071]  0.2235254448  0.1639587889  0.0361216711  0.0992325370  0.2568811687
#> [3076]  0.1636931206  0.0205372215  0.1065897073  0.2823440046  0.1282681761
#> [3081] -0.0292297662  0.0403998440  0.0326545120  0.1415057765  0.2952469776
#> [3086]  0.1814292128 -0.0518537512 -0.2081724159  0.3479879879  0.1148778840
#> [3091] -0.0892671495  0.0798296333  0.0851280135  0.2009946094 -0.0210964425
#> [3096]  0.1084810127 -0.0598836202  0.1035147548  0.1630025039 -0.2083621634
#> [3101]  0.0921836420  0.1413254462  0.1582558402 -0.1456252838  0.2565689504
#> [3106]  0.0073211202  0.0853624855 -0.0274947969 -0.0639557137  0.2451362243
#> [3111]  0.3642938720 -0.0415218236  0.3402411129 -0.0227378448  0.0824702502
#> [3116]  0.0442318336  0.1117707437  0.0524655829  0.1987093817  0.0668353770
#> [3121]  0.0703677229  0.2027033926  0.1973893085 -0.0506233177  0.1544836875
#> [3126]  0.0120951013  0.2387378564  0.0668895751  0.3890453123  0.2216299840
#> [3131]  0.0814349605  0.2068576738  0.0155618599  0.2329490629  0.3170452962
#> [3136]  0.0406748745  0.1699222993  0.1021662375  0.0261967139  0.2906663716
#> [3141]  0.1435001556  0.2082485130  0.2077911771  0.1148288192  0.1808488231
#> [3146] -0.0995070503 -0.0148877808  0.1303687357  0.2038418257  0.0595449049
#> [3151] -0.1078082961  0.2809821831  0.0857584926  0.0671606897  0.1908931574
#> [3156]  0.0281684369  0.1858380701  0.1143720282  0.2870252353 -0.0805980424
#> 
#> $mod1_unimodal$ref_properties$cell$pca_dist
#>    [1] 51.87341 46.67670 56.53049 59.05324 60.20071 53.52597 46.84286 54.55789
#>    [9] 46.80347 49.49131 57.61624 55.02936 45.96594 49.83714 48.37199 53.23457
#>   [17] 62.88841 59.22314 46.68545 52.90125 57.68060 63.24569 58.93584 51.39566
#>   [25] 64.07696 48.62531 60.25136 50.06205 63.59018 52.83642 62.97357 53.86443
#>   [33] 50.27888 59.17583 59.77984 48.55978 47.02272 49.28927 55.11271 63.70889
#>   [41] 61.53022 51.54498 47.84367 63.25134 46.50528 60.38734 43.48069 45.20812
#>   [49] 56.53444 62.27628 50.94805 49.74244 67.13517 60.94970 50.89631 47.60359
#>   [57] 57.04428 50.71333 56.50634 51.07065 51.20517 47.82880 50.86423 45.96653
#>   [65] 51.43282 59.47009 48.97518 46.34724 50.76814 54.90716 39.12867 57.75884
#>   [73] 45.50837 50.89627 62.35921 50.63855 59.15457 64.36343 51.58677 47.00574
#>   [81] 49.35145 39.93622 51.81713 47.10969 38.22890 50.89090 50.31117 39.96351
#>   [89] 47.82031 35.50842 50.28307 42.95358 52.03541 46.04271 37.87734 42.57306
#>   [97] 45.57890 40.79828 50.13896 56.12532 53.36836 27.14061 39.75644 41.99336
#>  [105] 45.17279 33.68247 52.49517 41.52167 57.16413 45.14079 52.02518 44.90203
#>  [113] 49.40434 47.52292 49.12756 46.14054 51.95753 48.45548 56.04846 43.48865
#>  [121] 46.66413 44.04889 47.63347 52.45537 43.35559 43.13742 34.44133 46.23753
#>  [129] 43.21218 47.37554 49.32853 51.43591 43.22983 39.75985 44.78882 47.97911
#>  [137] 37.84726 49.12531 51.37126 45.26094 41.11230 42.81881 44.25010 52.57595
#>  [145] 48.83060 47.64899 48.74919 55.64861 47.99249 53.44954 46.27110 37.33959
#>  [153] 51.29525 49.23865 47.60548 56.23435 44.83587 46.02722 51.08549 52.97473
#>  [161] 50.48329 45.23383 54.67293 51.28621 41.42298 41.97941 52.18856 61.80813
#>  [169] 52.13561 46.18426 51.75582 56.07914 54.10126 46.39830 45.95230 55.92464
#>  [177] 51.30056 54.05159 50.38046 55.72184 32.53156 60.31201 44.84726 58.38018
#>  [185] 54.14238 50.29196 52.84645 58.08179 57.99389 48.36792 51.43692 43.72323
#>  [193] 52.21585 47.95913 49.45161 51.33755 48.07408 51.78832 57.60513 42.15306
#>  [201] 49.26820 41.91904 44.88335 37.94173 53.03088 49.51933 46.18228 56.53431
#>  [209] 53.52800 50.48751 46.16198 52.86121 54.51321 51.40998 49.53866 54.51744
#>  [217] 38.29769 48.23198 48.14760 41.11374 51.77504 41.07050 45.38491 50.05394
#>  [225] 66.07320 46.46607 51.51978 42.53887 44.97790 45.99562 54.94642 52.44256
#>  [233] 47.46916 47.90643 50.95968 43.52642 43.09921 39.88423 40.19561 60.33334
#>  [241] 40.06643 43.07012 42.50484 51.06845 57.71282 53.49763 43.22815 55.14394
#>  [249] 55.58107 44.94297 35.67282 59.13292 54.31495 51.77199 57.38395 46.38130
#>  [257] 45.73629 55.52593 43.58450 61.63556 42.47835 50.36170 49.25027 54.77501
#>  [265] 53.40558 49.16701 49.19336 49.53642 47.65384 46.46349 43.67690 45.03912
#>  [273] 55.77302 54.35234 54.27427 51.07292 57.83448 39.85708 45.16917 52.94411
#>  [281] 52.40236 47.06670 54.02877 53.94969 57.55514 54.10519 39.57347 51.65695
#>  [289] 45.63209 47.92775 47.95395 42.37536 36.38506 52.31382 47.94866 46.98113
#>  [297] 44.93491 51.08039 43.78018 32.82889 52.27417 48.70456 59.50167 51.30435
#>  [305] 51.94476 50.76996 51.19930 53.08550 55.27476 43.27676 50.33667 33.45897
#>  [313] 40.54422 46.45981 50.42745 37.52759 51.97078 44.63378 57.26819 47.25406
#>  [321] 45.83337 52.42717 41.12381 51.40308 37.24406 27.80972 53.87463 50.99619
#>  [329] 54.68590 36.19606 45.31052 42.66714 58.38828 36.67663 53.69003 43.09274
#>  [337] 43.51928 44.90638 48.56959 30.94030 42.36952 38.43557 54.92098 46.70331
#>  [345] 48.72811 47.66629 56.78314 42.84671 49.88688 52.47029 30.83070 46.31410
#>  [353] 37.91466 46.78948 37.55813 46.57870 36.80753 47.98813 44.45888 44.19737
#>  [361] 54.08172 50.06568 49.06363 47.99380 42.82512 50.37003 46.32573 46.55215
#>  [369] 45.33364 37.23134 35.37002 34.11773 41.89460 40.74772 41.30575 55.70780
#>  [377] 40.13357 46.85674 51.60192 23.96375 53.10881 53.08730 41.79202 51.30491
#>  [385] 48.46511 45.78229 42.88498 44.85525 54.85190 51.42834 42.31601 45.99488
#>  [393] 57.70209 57.33792 65.47473 57.07866 59.06485 45.02678 51.28208 50.03422
#>  [401] 68.17957 51.62102 58.95153 53.17134 59.86945 53.05467 59.12482 42.17193
#>  [409] 65.40630 49.47749 49.76284 43.30199 49.46250 53.77476 50.92591 58.09581
#>  [417] 56.00292 46.65948 44.29270 54.30428 55.53309 61.87503 51.99227 53.55677
#>  [425] 54.43771 51.03652 48.15041 46.73399 51.30796 64.76594 54.40849 58.95348
#>  [433] 57.90257 63.11043 59.39458 60.60820 51.98763 50.04305 52.96029 53.71663
#>  [441] 46.93327 41.85637 51.63459 55.94525 50.53272 56.59823 49.01471 54.09486
#>  [449] 50.30016 53.92598 46.85382 49.61863 53.22925 50.51264 63.71969 59.22601
#>  [457] 57.05111 49.93181 52.57360 32.88173 34.58191 39.44685 28.53030 53.23744
#>  [465] 40.97538 54.84038 41.28653 42.70167 44.72638 47.50701 51.32849 32.26617
#>  [473] 34.16558 51.67907 42.60614 53.17276 45.15121 50.09764 36.17703 52.58050
#>  [481] 34.74582 51.69971 39.26030 43.20730 35.17867 43.99305 31.28642 39.22629
#>  [489] 29.87431 45.65854 40.76122 40.44936 39.37471 48.09010 47.57830 46.08551
#>  [497] 44.10996 42.79065 46.85192 26.44542 36.52718 46.61573 47.65330 43.84135
#>  [505] 49.19273 40.64781 50.86663 45.74570 47.07855 46.61657 43.66382 46.69641
#>  [513] 38.86606 37.68129 40.41854 41.00621 36.34418 29.34650 27.30225 41.26288
#>  [521] 33.30555 30.98654 58.91742 33.54456 39.90780 41.64251 31.23060 50.44466
#>  [529] 42.58630 49.94467 49.62470 40.37117 40.12928 48.85363 31.01011 41.05683
#>  [537] 40.33963 54.17199 42.27067 40.27465 46.51702 41.71008 47.64789 35.61606
#>  [545] 36.37583 48.09499 42.09220 54.07939 42.09828 43.31140 34.48368 48.88110
#>  [553] 27.21887 56.53109 43.77549 51.20273 34.67719 41.02775 46.62693 41.02214
#>  [561] 30.52250 39.02233 34.14170 45.00853 45.72595 49.01026 44.83487 34.24776
#>  [569] 45.09504 45.05601 51.33078 35.25781 41.89999 39.42177 57.58376 46.59254
#>  [577] 50.90741 42.06616 47.03994 47.37117 43.97933 49.18180 39.88668 47.32007
#>  [585] 38.30893 38.57182 41.52633 47.24278 36.51799 38.35867 42.54017 40.31506
#>  [593] 37.60402 40.41506 49.22617 32.71109 52.12228 30.56220 33.45713 53.71577
#>  [601] 41.89659 50.88219 47.43193 39.25825 44.00506 31.86480 56.59368 44.31960
#>  [609] 47.14621 50.17097 53.19035 39.53824 50.00370 54.65357 42.28453 32.92892
#>  [617] 52.04592 44.47966 58.04283 48.64898 51.36169 41.89097 55.94473 37.26239
#>  [625] 54.58120 27.04186 52.76933 43.64355 50.63334 47.36406 50.36152 47.74084
#>  [633] 52.56194 42.19173 48.05802 57.42928 51.69533 47.31665 48.61779 52.85018
#>  [641] 45.49573 47.79669 33.61448 50.49147 50.40559 60.00697 47.48149 48.82681
#>  [649] 51.54841 53.28292 56.77379 51.47744 44.59591 45.87188 53.36216 39.70601
#>  [657] 44.81587 36.41393 43.32327 40.00504 36.79923 46.78162 43.22746 46.96520
#>  [665] 43.96748 54.25878 42.73297 46.00901 48.01132 42.08993 53.04311 48.63106
#>  [673] 52.60234 49.62571 42.66863 43.07920 57.78193 50.99993 51.66723 47.44025
#>  [681] 58.80680 49.09962 54.30822 53.36732 47.96216 48.79510 54.45607 53.95910
#>  [689] 62.51305 39.82260 54.92322 38.51399 62.51451 44.52313 56.09883 53.73518
#>  [697] 56.24436 46.66072 54.48878 45.41760 55.21168 46.75010 51.69601 44.35408
#>  [705] 49.54539 53.18972 60.21740 50.48885 44.22477 51.79424 48.36121 45.95666
#>  [713] 36.26523 43.44562 44.03762 55.92646 54.78107 48.38653 48.69267 61.30332
#>  [721] 46.66429 53.86419 57.64391 50.95703 53.36235 36.00100 47.70221 42.90720
#>  [729] 42.54165 45.91888 32.16199 48.25489 36.73355 45.00201 50.06714 62.51237
#>  [737] 43.09923 39.25881 50.01552 40.34207 54.34348 53.30024 55.50759 53.47200
#>  [745] 49.16667 44.08032 44.26625 50.42973 38.38311 35.53796 40.41073 41.15516
#>  [753] 47.21606 31.79314 24.71336 45.74937 43.42080 52.69044 45.24951 43.13366
#>  [761] 35.67694 46.73978 32.04290 52.19178 34.08913 47.13054 44.43437 46.12063
#>  [769] 44.31116 48.27248 33.82478 47.30541 40.10319 47.78112 49.79760 48.24261
#>  [777] 38.83465 41.98553 45.40757 37.77118 45.08454 29.51191 38.84511 39.73566
#>  [785] 48.25235 40.73899 44.88347 41.06241 43.99948 44.56385 41.68425 43.24947
#>  [793] 50.89390 45.29893 36.74539 44.35535 34.68157 41.31276 37.38996 30.21451
#>  [801] 38.01190 39.72316 39.64990 39.53301 60.73349 39.27714 42.94866 38.94486
#>  [809] 31.94710 49.51251 49.62502 47.38486 53.98378 38.10407 55.18500 61.19613
#>  [817] 59.95376 55.06040 61.46545 55.19671 42.43278 52.70780 44.20755 64.40187
#>  [825] 49.03733 63.48094 53.09116 50.66532 50.80058 52.87545 45.63482 65.96187
#>  [833] 59.81574 58.06498 52.83594 53.94884 63.52864 60.67694 51.14945 54.12444
#>  [841] 43.03021 47.90037 53.39813 56.67118 55.24456 48.38882 53.03928 51.07676
#>  [849] 46.53321 52.68529 47.32762 40.74316 57.32983 50.25189 57.25871 53.28971
#>  [857] 55.01817 50.14869 49.73551 53.89059 55.85788 45.77190 52.51836 52.66596
#>  [865] 43.14098 55.22147 56.29754 49.12507 54.96622 49.20966 46.70952 50.44443
#>  [873] 56.41489 47.24232 57.76268 50.67002 49.72612 61.69235 64.18311 56.05887
#>  [881] 55.65070 55.91122 46.92346 51.87662 55.25851 50.68187 46.18324 59.50472
#>  [889] 38.20988 39.13201 56.42229 57.98341 55.38924 41.40050 46.45366 39.07009
#>  [897] 62.56826 35.45456 55.30073 35.82405 60.87515 32.50725 61.68044 50.16590
#>  [905] 46.68550 47.89093 47.27430 45.21642 49.47689 49.26345 54.69812 47.33977
#>  [913] 57.47756 49.74222 52.73487 56.49419 44.38052 43.80712 49.38725 55.68199
#>  [921] 47.53740 43.63819 51.33292 55.33580 45.46307 53.32026 51.36166 48.63422
#>  [929] 47.86852 49.73289 50.39843 39.17013 35.05245 42.17018 43.41568 52.64041
#>  [937] 47.39418 44.86022 40.07119 56.34174 48.15735 52.70924 48.57026 37.58474
#>  [945] 57.06162 49.73544 50.20705 61.81309 41.56126 51.47883 65.22466 54.12193
#>  [953] 55.30015 55.39913 43.74413 45.16280 52.66932 63.88046 59.65625 53.39835
#>  [961] 56.21511 57.93186 62.22095 46.74025 64.53215 43.57856 57.52578 51.90172
#>  [969] 54.33315 56.86871 58.13085 52.08844 54.30979 45.40003 62.87605 66.70093
#>  [977] 55.60573 47.17510 52.01756 60.49558 56.55993 62.64222 47.51341 56.71492
#>  [985] 62.00229 58.17362 45.11782 42.86670 68.16479 59.82395 50.56233 45.06776
#>  [993] 46.98976 58.77836 51.67992 50.26909 47.22156 42.07198 49.70202 49.61869
#> [1001] 51.39556 59.66812 49.06799 50.12915 47.37957 61.76109 53.99561 56.28243
#> [1009] 61.19874 52.05142 61.04323 56.38375 52.55397 70.36389 49.48101 55.94052
#> [1017] 47.53074 43.71690 47.59918 36.12969 48.16305 48.61333 55.53359 37.59274
#> [1025] 44.84056 54.26191 50.26491 54.72733 41.41648 62.46595 49.67292 46.93564
#> [1033] 47.25269 44.06271 44.85165 47.85108 44.89000 45.95917 47.08085 54.00931
#> [1041] 56.90987 41.01176 47.82673 38.32761 44.35953 49.53250 48.91308 33.43226
#> [1049] 44.87125 50.72150 49.64075 44.59049 51.86239 52.05284 61.59900 53.43186
#> [1057] 44.18878 49.45902 58.30713 49.47272 37.14971 48.31451 51.67017 44.82856
#> [1065] 48.06297 45.88578 47.28375 53.73601 48.55849 49.61136 64.67574 45.90268
#> [1073] 47.93862 45.19977 39.06099 47.51803 51.06713 52.66263 57.03327 44.93830
#> [1081] 54.95381 54.34692 60.91004 44.38437 40.24127 56.95060 57.19298 67.13381
#> [1089] 56.67058 59.94615 43.92972 57.86057 46.61032 55.72799 49.27390 60.77899
#> [1097] 51.99351 51.49042 56.00011 53.21152 34.30917 48.42485 50.85021 58.19949
#> [1105] 53.95712 58.74142 47.46598 56.68584 53.97580 40.67253 61.17966 43.60620
#> [1113] 52.67233 53.37886 59.69695 46.67329 52.81929 60.11572 41.57301 47.72028
#> [1121] 51.44920 55.16646 53.55625 52.42829 53.04767 52.42322 48.73746 50.20104
#> [1129] 43.94846 47.48604 49.84282 48.56420 47.39823 52.19802 61.47132 42.06884
#> [1137] 51.83157 45.50115 41.28915 60.03191 50.62733 54.52703 62.08001 51.99334
#> [1145] 52.38695 50.47680 50.11223 43.21490 59.64589 61.69830 53.19022 46.55392
#> [1153] 48.66670 43.53873 60.50693 38.25615 56.76052 42.09941 53.58660 52.87154
#> [1161] 54.55922 47.30124 52.97943 50.56139 53.65255 51.85263 44.28455 50.34516
#> [1169] 50.22243 50.76778 50.54298 46.78948 48.07181 57.88351 31.62781 45.27948
#> [1177] 49.80529 44.39653 52.51703 47.21047 56.06481 56.49700 55.63625 42.25089
#> [1185] 54.38382 51.91773 49.92578 44.27615 35.42221 48.72327 50.61910 41.54882
#> [1193] 41.16534 48.64810 51.08281 56.92759 47.74017 57.58675 48.90063 54.15843
#> [1201] 45.21025 48.67545 48.98919 39.14274 56.77437 57.68261 39.07605 56.66203
#> [1209] 44.06232 38.89881 44.36919 50.16387 52.97847 40.60126 31.89131 46.24754
#> [1217] 59.41560 40.12890 59.81771 50.40673 51.92570 46.76063 47.29581 54.54358
#> [1225] 43.17355 40.17614 49.77053 51.06402 53.72652 52.28160 53.35911 39.66623
#> [1233] 52.43329 45.86815 42.19777 43.84319 44.54092 53.64107 36.35840 51.98735
#> [1241] 34.10541 48.49893 44.44771 54.27880 44.17670 49.13795 46.12239 51.49558
#> [1249] 38.12175 49.20588 49.73446 46.69289 46.28536 47.34945 45.97446 54.38075
#> [1257] 51.18681 51.61434 51.94336 57.33942 46.80854 45.93095 49.14709 40.17695
#> [1265] 51.43442 50.13468 48.63924 52.47776 38.87467 53.24468 51.55966 61.56125
#> [1273] 54.02104 62.30617 46.30042 59.56370 58.21606 48.21911 39.54965 56.38324
#> [1281] 54.09007 51.73998 52.10269 47.48556 49.50197 56.81056 54.24860 56.87457
#> [1289] 39.56711 50.14186 55.16106 49.93859 55.55728 40.21203 38.34179 55.69875
#> [1297] 52.99772 46.78128 50.57942 48.35781 47.97225 51.68579 55.60898 62.84378
#> [1305] 54.18797 53.11079 45.29704 44.45102 57.22870 40.41680 48.75502 44.97219
#> [1313] 52.64999 53.53072 57.09102 49.24959 55.28706 51.98513 55.86193 59.25113
#> [1321] 60.73383 49.32035 50.48011 53.23746 48.08094 60.21383 57.62946 58.55634
#> [1329] 55.17376 55.63817 33.11770 39.59843 47.77505 45.67804 46.31390 51.92142
#> [1337] 40.90335 61.30037 33.23658 52.25894 38.84116 45.79419 39.55847 51.30459
#> [1345] 41.60825 38.04724 35.77389 47.33434 37.14575 47.81938 51.42535 41.80049
#> [1353] 38.33170 42.99253 47.60920 38.20160 51.17542 30.27362 39.44832 48.58618
#> [1361] 50.25172 39.24405 41.19830 49.33194 47.58711 49.91840 45.83484 46.96871
#> [1369] 48.08435 43.17731 40.29646 43.84617 39.84464 41.43807 40.88512 34.29815
#> [1377] 39.28556 39.67476 36.52259 39.32379 58.20025 34.82735 44.94445 44.39254
#> [1385] 28.54898 54.97932 52.84201 54.02670 56.40578 36.98651 51.50056 48.22261
#> [1393] 54.77434 42.84580 36.19207 38.49977 56.88559 32.98227 51.56606 31.91369
#> [1401] 46.56544 41.33472 50.28144 44.40522 46.72996 34.56330 47.67899 38.29826
#> [1409] 46.97290 49.80375 49.92389 37.68136 48.45027 50.79524 31.77858 44.44015
#> [1417] 35.19811 45.96383 40.07146 44.92652 31.46976 40.37491 45.36971 40.69611
#> [1425] 46.73943 42.57308 43.31937 47.90507 40.85207 42.49380 41.91623 34.06959
#> [1433] 41.08835 34.34177 32.83193 37.92628 40.03332 36.68476 33.99251 51.87827
#> [1441] 39.18651 47.65010 45.74202 28.96649 50.64754 49.44925 37.75763 54.27593
#> [1449] 39.97216 47.46335 53.33600 48.72610 58.15939 48.88401 55.50850 46.59990
#> [1457] 57.44489 54.88559 58.90182 61.22435 55.65346 51.33487 46.85903 49.12951
#> [1465] 59.67159 52.66059 62.97259 57.46893 51.74484 37.43419 46.10831 54.87667
#> [1473] 50.31494 58.93907 40.27281 57.50480 50.61401 64.10039 50.07998 50.31195
#> [1481] 53.63771 52.10973 58.58419 46.13779 51.40879 49.12835 50.59163 47.17863
#> [1489] 57.57021 52.20583 51.42138 50.58698 45.54921 51.11765 43.66455 50.00361
#> [1497] 59.07254 65.44996 51.54252 51.13716 59.49742 46.93104 50.58968 59.48070
#> [1505] 66.41253 54.40089 47.53659 64.20103 49.93856 59.15827 42.55704 52.37218
#> [1513] 44.14115 61.05710 52.76874 56.45463 48.33761 50.20328 56.47084 46.43803
#> [1521] 46.50913 53.15944 47.00899 42.75989 56.75751 55.24666 45.40590 51.90568
#> [1529] 54.41701 51.11251 37.62801 46.50194 56.19773 41.99260 60.65959 49.45778
#> [1537] 53.34881 52.86522 54.01316 54.63326 52.84541 38.83743 50.44676 50.49829
#> [1545] 53.42178 54.01018 47.31681 53.48509 46.13616 38.58094 48.90437 47.10318
#> [1553] 41.10498 50.25382 66.18743 38.73159 51.36447 46.67922 46.80466 54.71311
#> [1561] 51.80116 64.89292 37.91161 52.61588 53.91714 57.41387 58.83804 60.71034
#> [1569] 50.94440 68.40618 59.81471 50.79673 63.02176 56.45038 53.43728 55.95424
#> [1577] 58.84191 58.60962 61.24435 52.69830 53.23773 41.53682 52.93327 50.49070
#> [1585] 58.00964 51.59496 53.74999 38.56423 45.71881 58.59934 55.55363 47.59424
#> [1593] 54.33359 50.62802 65.58584 67.62727 49.49205 59.45453 57.63112 56.74870
#> [1601] 49.32233 55.54789 59.88104 46.97755 55.76512 57.18568 47.82996 59.95312
#> [1609] 54.70604 51.15709 60.09622 57.09917 56.75604 58.46371 52.11830 42.88140
#> [1617] 59.49695 56.50512 61.77556 46.14232 45.59737 37.37106 53.30394 33.67664
#> [1625] 52.88766 47.41630 59.83729 45.11155 53.09909 39.83484 49.92588 47.33192
#> [1633] 52.83713 44.44235 43.87693 50.53660 58.06031 37.32576 43.76928 45.67746
#> [1641] 42.10565 42.23462 39.15570 45.72643 29.29373 47.96753 45.26707 41.77694
#> [1649] 46.08260 49.18239 49.98554 45.83431 49.04374 45.00128 41.39201 46.13513
#> [1657] 48.47004 49.78027 38.37091 34.45772 34.91971 46.61352 41.66984 44.71216
#> [1665] 49.36357 51.34669 41.59794 49.92870 47.47245 30.35320 48.12957 46.31513
#> [1673] 51.04375 48.96480 48.43424 47.01277 63.41662 41.21421 57.13679 53.56835
#> [1681] 52.47859 49.73147 58.66641 57.63155 53.98973 43.35954 50.11553 50.02266
#> [1689] 51.99364 46.80596 60.27769 51.04837 52.18168 57.72064 48.22170 51.13306
#> [1697] 49.31967 49.60844 38.30919 47.05467 42.13564 52.60535 38.18919 58.53285
#> [1705] 47.59579 47.80222 54.05125 50.25259 45.67507 50.29709 50.27354 50.22905
#> [1713] 51.73452 50.54192 49.23010 52.87865 54.65738 47.67660 45.72317 55.52056
#> [1721] 51.61532 59.29252 46.86175 46.00358 55.94069 52.32091 42.46352 59.43136
#> [1729] 36.15869 60.36694 35.51595 54.99917 45.03940 54.87747 39.01123 55.65261
#> [1737] 47.17459 39.38193 39.53120 40.49152 46.12331 42.59924 44.78215 54.46262
#> [1745] 40.37027 49.95026 51.13244 42.11958 46.96731 34.47514 44.74341 31.27599
#> [1753] 52.81856 50.22797 41.36370 45.94648 49.85699 44.76196 46.47713 52.29126
#> [1761] 45.56141 50.07128 42.48242 45.52528 35.12305 43.45737 30.15465 23.71389
#> [1769] 42.75093 32.73743 37.95080 40.34656 61.75038 38.34188 48.93294 36.24158
#> [1777] 33.65097 44.09732 40.49844 52.54464 42.71669 42.25926 52.39561 69.10687
#> [1785] 55.76999 65.38039 66.27192 44.25307 55.48549 66.58670 55.62644 64.01781
#> [1793] 60.13875 61.06896 58.35658 63.36868 54.42090 52.98905 55.64137 62.45551
#> [1801] 59.35483 50.75592 53.98659 55.75975 63.36812 55.24069 63.21625 57.91271
#> [1809] 58.04629 56.37995 44.42427 46.88605 58.57014 56.40833 59.27382 59.60020
#> [1817] 57.11335 55.21453 53.16050 57.42260 54.11507 59.45842 54.84160 55.85670
#> [1825] 67.78629 58.66663 59.97913 60.63325 58.98454 54.95836 57.58882 62.65806
#> [1833] 63.72526 61.80938 43.61011 37.06805 46.78052 35.51530 48.42583 42.69740
#> [1841] 42.96364 35.54594 44.63720 31.47421 38.22034 45.03051 47.54092 33.21706
#> [1849] 34.37542 42.42855 41.37712 47.83672 32.75697 39.05735 37.11326 47.75434
#> [1857] 42.40115 32.84209 41.15147 34.17617 46.73506 45.50966 40.40222 44.11267
#> [1865] 43.58340 39.74724 35.15431 43.30605 36.16966 33.62788 32.32006 41.96710
#> [1873] 33.01213 41.90236 43.96681 52.52076 37.11507 43.74982 31.87783 25.92800
#> [1881] 51.83426 40.27633 44.94698 51.33373 29.32418 51.65671 58.84686 52.78566
#> [1889] 71.71428 54.48959 56.51047 54.07510 56.16710 44.33371 53.87072 48.64848
#> [1897] 59.74407 57.43028 54.45462 58.01375 52.02588 62.27784 51.02890 53.39453
#> [1905] 57.02397 52.46061 51.73806 53.66749 59.29073 50.73708 57.59829 59.01753
#> [1913] 49.63053 43.60388 49.34466 54.52119 62.06952 62.60499 44.59226 61.76694
#> [1921] 54.81153 58.31791 55.71630 50.98086 64.61444 59.46272 49.83463 52.41504
#> [1929] 51.23045 46.24697 67.08098 59.30487 57.97306 59.85678 43.48372 55.64819
#> [1937] 42.30405 55.22878 46.86837 49.41355 49.99960 53.22097 46.37225 52.53927
#> [1945] 60.59132 55.52937 46.33942 57.09008 53.98229 46.14202 55.78932 40.67159
#> [1953] 49.82155 52.92400 54.42825 43.66399 46.25225 60.72430 51.83698 53.18439
#> [1961] 49.86771 38.25665 50.89287 48.96042 51.84872 50.81293 38.01909 42.57885
#> [1969] 40.06489 42.23054 52.80180 50.26537 48.97541 44.27757 59.75472 45.42810
#> [1977] 51.72016 48.97536 41.22886 58.93140 49.45436 52.10621 59.41427 44.21116
#> [1985] 53.88322 45.61068 45.28269 38.99257 47.48257 51.73268 54.29581 53.06346
#> [1993] 47.12469 49.00506 58.68606 50.40074 60.31107 45.28382 54.99335 40.99772
#> [2001] 50.16217 52.32981 48.29560 42.96557 56.47489 55.44960 61.69509 58.21306
#> [2009] 50.52877 48.49374 54.55197 50.90401 52.00341 44.09449 48.81704 57.97556
#> [2017] 57.63489 49.50824 45.73151 50.30776 50.32954 48.55904 68.06290 49.51713
#> [2025] 45.67070 53.60962 49.59605 53.87548 56.29464 46.59455 50.04954 42.56230
#> [2033] 52.90912 51.01221 38.25300 38.04713 30.62202 31.69040 44.52841 49.45503
#> [2041] 48.95037 53.33220 50.91886 44.35949 52.74286 49.58676 44.58098 49.86185
#> [2049] 49.02123 53.34677 45.24753 49.69898 56.56505 58.00566 40.85593 53.42790
#> [2057] 46.06678 46.32446 44.99290 43.43046 43.40354 41.21562 41.90936 44.65513
#> [2065] 42.38776 52.54228 49.19259 38.32335 41.74063 55.10887 35.06341 49.35047
#> [2073] 37.03100 35.08561 59.28395 41.68615 48.99610 49.88713 44.01675 50.03114
#> [2081] 47.81969 39.96168 50.38571 53.38134 54.64006 56.20338 54.33333 52.28119
#> [2089] 47.35202 50.35468 47.87382 51.88297 40.93511 53.68677 52.32727 60.30718
#> [2097] 42.32464 58.42515 59.91266 59.26609 49.66529 45.38554 46.25594 52.33689
#> [2105] 49.13853 51.55751 38.05722 47.36069 56.80271 47.17900 51.42954 50.89565
#> [2113] 51.35001 51.85978 48.95699 57.68634 42.26122 46.01159 52.48121 51.80302
#> [2121] 53.87620 44.95603 53.70186 51.54310 53.35249 42.60985 46.54323 59.63291
#> [2129] 50.16242 52.09216 44.12958 53.73790 52.96593 47.08950 49.32153 46.28256
#> [2137] 57.93085 30.84637 45.71077 49.68688 47.53431 50.83049 55.83515 50.91171
#> [2145] 55.92380 59.85378 45.77438 54.84400 45.72043 47.78681 46.21134 43.92260
#> [2153] 52.25922 50.52130 42.80802 36.97904 29.70537 45.47371 45.00583 43.94869
#> [2161] 60.82277 45.85955 49.25587 58.53907 36.78673 50.63432 52.32073 55.95259
#> [2169] 52.15990 52.62921 41.99144 41.28691 49.25099 50.22277 43.02194 44.12204
#> [2177] 52.38580 54.07503 48.33189 46.65701 55.89925 38.61405 54.01092 46.64624
#> [2185] 52.66545 44.42203 53.92244 58.98960 59.64438 50.69508 47.15511 47.52741
#> [2193] 42.33268 42.00485 49.68644 42.98563 42.56270 54.16950 48.04492 42.11984
#> [2201] 47.70471 46.84553 44.94679 46.36713 61.98356 39.80656 47.59751 48.48986
#> [2209] 42.78714 50.13942 44.38763 57.26698 36.85958 43.44945 37.61238 38.57090
#> [2217] 47.31442 42.88973 51.26688 40.60285 44.03598 47.58302 41.68529 51.06491
#> [2225] 36.12251 45.82745 45.99769 48.29008 37.11357 44.52033 45.11178 43.35463
#> [2233] 36.05325 42.71472 46.91142 46.29573 44.42507 42.46075 34.77186 45.32292
#> [2241] 46.12509 35.68115 37.78271 37.23447 41.56794 32.05561 36.32192 56.83276
#> [2249] 32.54739 46.07649 37.21385 33.43015 53.39543 37.16320 45.00341 55.54110
#> [2257] 38.63376 40.03011 51.48125 47.43092 43.75409 54.89520 49.36464 51.40022
#> [2265] 54.39936 53.15901 45.19449 54.69901 51.96521 53.15438 43.17669 51.45426
#> [2273] 60.89917 60.78125 39.12479 43.96575 47.61526 48.40514 49.38907 41.51592
#> [2281] 46.72395 44.59082 45.97658 48.59168 50.79833 55.23695 52.97807 37.25569
#> [2289] 44.80221 58.09841 41.18263 54.16626 34.76824 45.21402 48.74250 37.24973
#> [2297] 43.69437 48.41520 44.36841 46.67544 51.09947 43.43628 48.64497 31.93512
#> [2305] 42.85921 49.70674 49.71834 40.96380 49.16269 49.00963 52.57121 45.12066
#> [2313] 47.78122 54.03767 47.61988 47.79262 46.15945 46.19988 42.94650 41.94879
#> [2321] 32.18892 42.12775 45.67491 42.67727 48.00439 41.68960 48.12747 42.81907
#> [2329] 35.51533 47.97244 46.59533 34.69265 50.23596 41.30367 33.38648 59.54630
#> [2337] 51.23275 50.72875 51.62631 44.78856 42.93665 52.66614 50.32543 52.78181
#> [2345] 49.34644 43.86711 37.61360 39.83445 34.95096 42.53652 42.72323 46.52729
#> [2353] 45.57158 48.50029 55.82832 51.39214 51.39746 47.34777 44.35792 43.82171
#> [2361] 53.19426 41.22248 51.16500 45.53216 44.28645 39.58967 45.18210 50.75958
#> [2369] 43.01396 41.69544 52.12571 32.36719 50.28290 37.83507 45.47988 55.72638
#> [2377] 42.04628 59.78622 47.99537 40.25255 48.02988 58.32295 52.96591 46.40325
#> [2385] 52.25124 56.98235 39.43666 41.76473 48.59208 43.93391 49.11859 56.18412
#> [2393] 43.75953 57.02956 47.44584 44.49372 56.53264 38.68180 43.40201 47.71388
#> [2401] 44.08856 52.99320 47.73922 53.51736 47.64370 38.51789 50.85707 42.39189
#> [2409] 43.87275 56.01582 48.19001 49.72180 50.93332 47.83318 46.63240 46.69590
#> [2417] 52.86069 49.68697 42.03639 54.71718 47.14191 38.29499 56.21675 54.32562
#> [2425] 36.58866 56.28212 60.46978 51.04347 45.96150 54.14437 62.20538 60.07610
#> [2433] 57.94240 41.18293 51.65522 54.24499 45.14589 35.05929 45.45621 50.61349
#> [2441] 47.52541 57.18818 51.88693 47.76932 55.96707 49.57300 51.46634 58.38989
#> [2449] 51.40221 50.98073 55.59202 48.83795 40.83787 53.55098 57.44733 50.53685
#> [2457] 48.80928 48.80893 52.46714 39.22863 46.22522 39.01511 52.40323 41.22641
#> [2465] 54.25502 39.56655 22.13689 49.67448 30.50376 49.37854 47.07851 43.26103
#> [2473] 56.78562 48.85314 49.82665 47.70935 46.64924 38.00181 29.91701 34.81809
#> [2481] 46.45486 33.96964 46.67156 49.80070 62.75509 46.44836 45.88094 49.11488
#> [2489] 33.72787 46.93236 47.57488 50.56997 59.26320 43.27038 49.96008 50.43385
#> [2497] 53.77637 36.41097 47.31416 47.22965 58.84107 52.76274 53.16482 46.42062
#> [2505] 51.14174 59.47503 45.06439 55.14864 46.65003 52.97501 29.84736 44.11448
#> [2513] 54.88682 50.00589 49.66093 45.35317 45.37485 41.36232 44.72390 54.89941
#> [2521] 51.92020 43.31066 52.94838 44.00367 39.16064 53.95137 54.84236 53.86386
#> [2529] 52.01918 45.81897 56.27310 52.28195 39.85009 51.40004 51.25126 49.29494
#> [2537] 50.41451 54.82926 58.50361 56.43730 45.32648 46.40690 52.11966 53.90905
#> [2545] 34.71573 40.53808 42.47140 50.18398 45.33285 53.75102 45.95720 51.98246
#> [2553] 54.02976 55.81867 56.17873 57.62545 48.10685 43.42926 54.07479 44.47870
#> [2561] 52.01449 50.39140 63.17343 53.99598 51.99399 40.99145 34.76003 43.15094
#> [2569] 36.84650 49.15535 32.15069 43.17095 49.69636 45.29951 53.28721 49.94197
#> [2577] 49.22172 44.34953 41.53673 50.67102 47.41089 43.61547 44.91062 40.99816
#> [2585] 37.95750 44.40244 40.63895 44.99462 47.82881 49.20499 32.78742 41.71558
#> [2593] 47.33551 34.54761 53.95440 53.10947 48.12801 49.19665 42.10817 45.27633
#> [2601] 50.04311 36.68335 51.10428 39.04863 49.19805 48.00881 59.53984 51.82416
#> [2609] 55.08070 44.79553 55.31865 45.26661 49.66990 55.47384 47.53904 42.06024
#> [2617] 49.14356 41.61611 50.89366 51.74740 44.99526 50.14537 58.34618 40.94800
#> [2625] 43.62237 49.68490 43.07008 49.72970 55.97799 52.59693 44.01976 51.90465
#> [2633] 35.48656 43.01944 44.87869 40.35965 41.02469 43.09298 47.35120 49.28030
#> [2641] 35.76807 47.37318 39.88251 42.50333 28.52739 31.78674 38.51479 37.44597
#> [2649] 33.57896 27.89961 28.05328 34.57342 38.58205 36.84496 51.11514 34.84282
#> [2657] 36.85645 43.12652 32.42460 36.63583 40.28123 49.49933 44.71364 35.56322
#> [2665] 45.26623 45.29306 48.67123 49.82018 40.85492 58.74460 47.66598 49.15565
#> [2673] 53.92695 46.36143 46.39499 49.38355 46.89385 47.43144 43.26698 48.79185
#> [2681] 41.79855 43.84169 48.72291 43.36857 40.59946 56.98725 38.75866 47.67120
#> [2689] 41.16745 43.77366 59.63356 53.32364 57.88984 61.22552 36.03875 48.87864
#> [2697] 45.52658 47.19757 41.27243 51.28301 47.73000 47.72437 49.53879 48.11521
#> [2705] 41.48543 50.05601 50.71766 42.92667 49.05978 41.56642 32.59878 48.90774
#> [2713] 39.72040 45.67406 50.10332 58.82368 41.00900 49.43828 42.65553 36.70478
#> [2721] 50.01184 51.89859 52.08447 41.02402 45.90616 40.56127 49.34351 58.36423
#> [2729] 60.27831 46.37488 37.96367 46.76249 57.68353 34.90374 55.70671 47.25850
#> [2737] 55.61573 51.88239 52.28331 46.27918 45.89606 60.83567 46.24867 43.66777
#> [2745] 66.31693 48.96048 56.74711 52.37055 48.42392 52.79810 49.51373 53.25470
#> [2753] 61.28820 46.75996 41.26143 54.03421 48.75712 40.80113 40.01645 32.15720
#> [2761] 46.59458 31.39747 50.10964 46.13280 42.52971 38.47100 44.82823 47.14984
#> [2769] 47.44365 50.93134 38.58892 41.11052 49.76759 38.13704 42.65379 49.09881
#> [2777] 39.53148 48.47492 46.76274 37.04549 52.52235 39.98169 55.22862 35.51701
#> [2785] 46.73381 49.94167 43.88154 58.68303 50.50540 51.38424 44.34093 45.71243
#> [2793] 34.15297 34.27546 37.16225 50.05211 35.06466 47.63766 47.94570 63.96989
#> [2801] 48.28705 45.04670 47.64450 37.51617 49.60577 44.87008 46.68881 63.10142
#> [2809] 39.36987 55.18879 58.79670 58.47201 62.32491 51.94342 58.30282 46.95421
#> [2817] 55.05629 57.53441 47.19074 49.82307 43.85927 39.71236 47.81865 46.90644
#> [2825] 48.82536 61.86075 52.27297 52.79181 49.52188 41.42550 54.52017 59.46654
#> [2833] 53.14671 59.58736 41.56295 58.03449 55.74310 50.56781 58.37787 55.52464
#> [2841] 55.72361 50.53083 56.83226 45.29331 40.94201 45.11668 48.14332 37.49091
#> [2849] 52.58646 58.58440 62.71496 52.31353 51.78355 53.20045 37.08012 58.41452
#> [2857] 54.37864 51.22218 65.13418 50.46348 45.09633 43.76106 54.91289 38.22636
#> [2865] 52.89455 47.45088 39.35001 44.68878 49.93498 46.79762 58.01817 51.57172
#> [2873] 43.95569 45.41703 65.33856 44.70930 45.13722 46.41286 48.68046 58.12435
#> [2881] 44.49877 51.16670 61.74792 45.66801 44.03632 47.17302 36.94026 43.39994
#> [2889] 41.18765 44.82210 51.95163 43.44564 44.68396 42.39345 51.53879 38.91264
#> [2897] 38.75000 56.79908 46.17908 58.18442 51.38979 49.06543 38.69681 42.94130
#> [2905] 52.56756 52.46093 48.00346 52.16346 39.97246 56.99619 51.55471 45.25873
#> [2913] 43.96045 48.28225 48.42761 55.76425 53.67093 44.87602 49.28204 64.02670
#> [2921] 43.76997 46.72799 47.48356 47.93994 54.81411 45.85119 50.28049 51.80880
#> [2929] 44.42642 42.63647 47.53120 48.36894 50.94171 49.19918 50.73490 49.44466
#> [2937] 49.06013 49.73549 40.12680 49.07933 33.94980 33.76148 53.36663 48.17964
#> [2945] 47.17727 57.20669 47.42063 58.89095 40.70729 42.79656 50.17093 45.04409
#> [2953] 43.82298 46.88774 52.61639 45.20612 51.60148 53.59733 45.01196 48.20398
#> [2961] 49.16426 38.03997 46.00343 52.80801 43.96148 54.83635 51.24711 56.19287
#> [2969] 49.10908 45.62441 42.42071 44.59834 41.08978 47.91233 39.58089 41.67808
#> [2977] 41.43012 42.48460 48.67334 52.64786 45.39203 45.14211 47.96141 36.92775
#> [2985] 42.48665 51.82274 50.07003 49.07869 45.32666 44.34415 52.27591 38.37588
#> [2993] 38.65589 40.64775 41.20644 47.56646 38.72902 52.81909 39.81129 47.43439
#> [3001] 46.14688 44.31601 51.12896 34.18592 52.05578 53.84399 41.27622 45.59755
#> [3009] 43.41919 34.88408 51.37258 37.31956 43.75387 38.88710 60.81687 43.25855
#> [3017] 41.98926 50.80813 44.15476 50.48179 52.13803 50.10457 49.33050 46.40504
#> [3025] 43.24311 43.28046 48.08497 43.70575 43.03841 48.29101 53.54830 46.27980
#> [3033] 35.92992 47.91006 32.96322 45.88885 49.04245 41.58298 57.70725 40.92135
#> [3041] 28.42765 37.01824 36.21916 39.38319 34.61121 58.43487 39.54350 53.12481
#> [3049] 42.81789 32.67838 44.82019 33.24146 48.47701 54.54294 46.67206 32.90134
#> [3057] 27.09014 37.76973 37.81476 62.83268 36.39103 41.99075 44.00029 28.49776
#> [3065] 46.16448 45.77581 54.08582 46.26120 43.73624 42.12351 34.18080 32.79991
#> [3073] 60.82134 44.43401 51.28960 53.44483 35.72946 42.09058 49.57084 51.29191
#> [3081] 53.22749 49.19512 44.59756 49.06609 61.04851 44.28379 39.58879 50.03289
#> [3089] 33.21310 50.00596 52.30769 50.61078 50.57953 45.65635 28.27631 54.53972
#> [3097] 32.23751 54.01004 42.33992 35.90207 49.46735 46.15859 46.75351 49.20620
#> [3105] 45.66187 58.75835 40.55545 56.79794 47.63150 43.48884 48.53416 42.02605
#> [3113] 47.52971 57.54929 45.07532 44.67578 63.07572 58.98897 54.13228 65.95341
#> [3121] 58.22146 60.22696 57.18432 55.35126 46.05523 32.86601 35.68510 58.88219
#> [3129] 41.30428 56.27879 42.75054 40.94217 55.16754 41.82229 52.94196 55.31150
#> [3137] 51.91496 51.92707 44.55309 41.05701 59.53726 37.34630 53.12994 52.46660
#> [3145] 33.77392 49.60246 48.51378 41.83374 50.61970 42.26726 47.40429 49.17758
#> [3153] 46.94779 54.58645 52.34249 52.40445 42.40125 62.07382 48.59207 56.20984
#> 
#> $mod1_unimodal$ref_properties$cell$knn_hubness
#>  [1]  0  5  2  1  5  1 11 10  1  1  9  1  4  0  1  0  0  1  1  5 11  0  3  1  2
#> [26]  2  6  0 18  0  3  0  4  1  2  2  8  3  7  4  1  4  4  3  3  2  2 27  4  8
#> [51]  2  9  4  1  1  0  5  1  4  7  6  3  5  4  5 11  7  5  4  3  1  9  1  5  9
#> [76]  3  2  0  2  7
#> 
#> $mod1_unimodal$ref_properties$cell$prop_outliers_cell
#> [1] 0
#> 
#> 
#> $mod1_unimodal$ref_properties$feature
#> $mod1_unimodal$ref_properties$feature$mean_expression
#>  [1] 12.862039 12.220038 11.290137 11.942478 13.316057 12.340675 12.231143
#>  [8] 12.639862 13.452811 12.638611 13.079575 12.410667  9.306157 11.662829
#> [15] 11.914306 12.276579 11.181936  8.945510  9.645397 12.577307  9.609804
#> [22]  9.996135 10.585018 12.255251 12.132763 12.173704 11.179566 10.325577
#> [29] 11.351419 10.640077 11.390373 11.800145 11.515839 10.069515 11.138959
#> [36]  9.605871 10.896118 11.107756  8.877162 11.887138 11.568087 10.816459
#> [43] 12.435335  9.965387  9.459759 11.142859 10.461643 12.358673 11.217901
#> [50] 12.251147 11.714935 10.467915 12.358032 11.977190 11.451387 11.305978
#> [57] 11.932945 11.332935  8.062381 11.256280
#> 
#> $mod1_unimodal$ref_properties$feature$variance
#>  [1] 20.17357 20.55967 33.62462 24.57170 13.88677 20.81701 25.41493 16.92489
#>  [9] 11.91160 19.45229 15.63600 24.07219 40.93384 32.59825 26.52592 20.28092
#> [17] 30.16775 39.64102 34.56522 19.22258 38.68067 34.67953 36.33711 22.93984
#> [25] 25.33686 20.27472 30.35738 30.10723 31.04307 34.35202 28.99541 26.15281
#> [33] 29.41980 35.23599 32.45479 36.26078 30.85527 29.90967 39.10881 26.29925
#> [41] 27.57009 30.55304 18.53221 34.62153 39.56805 30.18403 32.99389 20.90096
#> [49] 30.66678 23.37569 28.17614 35.57424 23.77448 26.80421 26.91426 31.03658
#> [57] 26.81005 33.67446 40.06448 26.16351
#> 
#> $mod1_unimodal$ref_properties$feature$scaled_variance
#>  [1] -1.16596178 -1.11053554  0.76499234 -0.53459240 -2.06845861 -1.07359360
#>  [7] -0.41354277 -1.63232376 -2.35200347 -1.26950512 -1.81734986 -0.60629857
#> [13]  1.81426162  0.61765251 -0.25405618 -1.15055112  0.26874440  1.62867237
#> [19]  0.90002026 -1.30248139  1.49080968  0.91642933  1.15438230 -0.76885300
#> [25] -0.42475031 -1.15144208  0.29596620  0.26005710  0.39440046  0.86941386
#> [31]  0.10045020 -0.30761759  0.16137248  0.99631193  0.59705834  1.14342538
#> [37]  0.36744092  0.23169614  1.55227092 -0.28659479 -0.10416074  0.32405442
#> [43] -1.40158701  0.90810298  1.61819743  0.27108220  0.67444852 -1.06154299
#> [49]  0.34038253 -0.70628371 -0.01716017  1.04486925 -0.64903621 -0.21410590
#> [55] -0.19830722  0.39346926 -0.21326744  0.77214704  1.68946205 -0.30608108
#> 
#> $mod1_unimodal$ref_properties$feature$sd
#>  [1] 4.491500 4.534278 5.798674 4.956985 3.726496 4.562566 5.041322 4.113987
#>  [9] 3.451318 4.410475 3.954238 4.906342 6.397956 5.709487 5.150332 4.503435
#> [17] 5.492518 6.296111 5.879220 4.384356 6.219378 5.888933 6.028027 4.789555
#> [25] 5.033573 4.502746 5.509753 5.487006 5.571631 5.861059 5.384739 5.113981
#> [33] 5.424002 5.935991 5.696910 6.021693 5.554752 5.468973 6.253703 5.128280
#> [41] 5.250723 5.527480 4.304905 5.884006 6.290314 5.494000 5.744031 4.571756
#> [49] 5.537760 4.834842 5.308120 5.964415 4.875908 5.177278 5.187896 5.571049
#> [57] 5.177842 5.802970 6.329651 5.115028
#> 
#> $mod1_unimodal$ref_properties$feature$cv
#>  [1] 106.80192  94.47199 134.40497 106.94000 102.43863  88.44417  88.44609
#>  [8]  85.71559 105.69305  93.87292  92.19436 113.22554 160.03210  96.71457
#> [15] 101.15520  82.28975  93.43986 118.47808 128.53677  94.28025 122.32675
#> [22] 110.93329  90.03075 102.75258 109.39111  95.24258 113.28536 115.51572
#> [29]  86.23895 100.38950 109.14360 101.04325  93.49135 115.03561 111.67178
#> [36] 104.82284  93.40967 112.37991 114.20099  89.39993  94.94627  94.39124
#> [43]  83.16354 115.72165 111.54921 111.91075  90.03886  90.64175 106.99082
#> [50] 103.72842  93.82192 106.72478  91.46087  89.93028 111.33267  98.49431
#> [57]  94.82855  95.98240 114.70480 105.67711
#> 
#> $mod1_unimodal$ref_properties$feature$zero_fraction_feature
#>  [1] 0.1000 0.1125 0.2000 0.1375 0.0625 0.1125 0.1375 0.0875 0.0500 0.1000
#> [11] 0.0750 0.1250 0.3125 0.1875 0.1500 0.1125 0.1875 0.3250 0.2625 0.1000
#> [21] 0.2875 0.2500 0.2375 0.1250 0.1375 0.1125 0.1875 0.2125 0.1875 0.2250
#> [31] 0.1750 0.1500 0.1750 0.2500 0.2000 0.2750 0.2000 0.1875 0.3250 0.1500
#> [41] 0.1625 0.2000 0.1000 0.2500 0.3000 0.1875 0.2250 0.1125 0.1875 0.1250
#> [51] 0.1625 0.2375 0.1250 0.1500 0.1625 0.1875 0.1500 0.2000 0.3750 0.1625
#> 
#> $mod1_unimodal$ref_properties$feature$detection_freq_feature
#>  [1] 0.9000 0.8875 0.8000 0.8625 0.9375 0.8875 0.8625 0.9125 0.9500 0.9000
#> [11] 0.9250 0.8750 0.6875 0.8125 0.8500 0.8875 0.8125 0.6750 0.7375 0.9000
#> [21] 0.7125 0.7500 0.7625 0.8750 0.8625 0.8875 0.8125 0.7875 0.8125 0.7750
#> [31] 0.8250 0.8500 0.8250 0.7500 0.8000 0.7250 0.8000 0.8125 0.6750 0.8500
#> [41] 0.8375 0.8000 0.9000 0.7500 0.7000 0.8125 0.7750 0.8875 0.8125 0.8750
#> [51] 0.8375 0.7625 0.8750 0.8500 0.8375 0.8125 0.8500 0.8000 0.6250 0.8375
#> 
#> $mod1_unimodal$ref_properties$feature$dispersion
#>  [1]  7.799297  3.993918 10.071068  5.603720  7.516093  3.637403  4.370940
#>  [8]  3.499086  8.420157  4.967835  5.323000  8.413140  8.355352  5.179615
#> [15]  5.052235  2.818679  3.208645  2.982874  3.841296  4.688824  4.189873
#> [22]  3.199611  3.049708  5.424220  6.686235  3.764526  5.101345  3.202532
#> [29]  3.170088  3.602903  4.988286  4.645431  3.736620  3.837625  5.424705
#> [36]  2.595874  2.824837  4.814899  2.592071  3.656499  3.876361  2.739734
#> [43]  2.965310  3.632445  3.390779  4.790436  2.340895  3.974455  4.721903
#> [50]  5.890873  4.181213  3.901136  4.830840  4.104384  4.896010  4.219991
#> [57]  4.552430  4.640858  1.957132  3.810961
#> 
#> $mod1_unimodal$ref_properties$feature$bcv
#>  [1] 0.9972024 0.8179438 1.2755773 0.9692960 0.9538088 0.7531158 0.7767240
#>  [8] 0.7243912 0.9921855 0.8389451 0.8308425 1.0628365 1.5015039 0.8687830
#> [15] 0.9059282 0.6609992 0.7752367 0.9659799 1.1054694 0.8362434 1.0673521
#> [22] 0.9197848 0.7380873 0.9279882 1.0087966 0.8161807 1.0157682 0.9579778
#> [29] 0.7135205 0.8532789 0.9759235 0.8950941 0.8000906 0.9891873 1.0085503
#> [36] 0.8218895 0.7507705 1.0003146 0.8950094 0.7620075 0.8178761 0.7521751
#> [43] 0.6770386 0.9851328 0.9366697 0.9954716 0.6814539 0.7841392 0.9498839
#> [50] 0.9451509 0.8183697 0.9203516 0.8144624 0.7821133 0.9931423 0.8603650
#> [57] 0.8376841 0.8501476 0.8021535 0.9075929
#> 
#> $mod1_unimodal$ref_properties$feature$gene_cor
#>    [1]  1.165479e-01 -4.466283e-02 -4.305952e-02  5.510029e-03  6.136618e-02
#>    [6] -2.708937e-01 -1.151011e-01  1.573380e-01 -1.160835e-01 -2.345642e-04
#>   [11] -1.578053e-01 -3.922623e-02 -8.842931e-02 -9.493410e-02 -5.876795e-03
#>   [16] -7.699883e-02 -1.020776e-01 -1.807130e-01  2.980437e-02 -1.084929e-02
#>   [21] -9.984687e-02 -3.074126e-03 -4.559209e-02 -8.048548e-02 -1.377085e-01
#>   [26]  2.698579e-02 -6.942772e-03  1.471521e-01  9.996405e-03  1.155636e-01
#>   [31]  2.329461e-02 -4.099895e-02 -1.153111e-01 -5.521158e-02  1.527621e-01
#>   [36]  8.163119e-02  2.005559e-01  7.856645e-02  1.341374e-01  5.597908e-02
#>   [41] -7.441323e-02 -1.417602e-01 -2.952198e-02 -1.070693e-01 -9.771086e-02
#>   [46] -1.639368e-01 -1.421822e-01 -1.499075e-01  9.801820e-02  7.924901e-02
#>   [51] -9.115419e-02 -7.524190e-02  7.022695e-03 -1.313507e-01 -1.539094e-01
#>   [56] -4.558309e-03 -9.441596e-02  5.703230e-02 -2.110884e-02 -1.676752e-02
#>   [61]  2.131429e-01 -3.952735e-02  1.256085e-01  5.651109e-02 -1.059875e-01
#>   [66] -1.386723e-01  4.494135e-02  9.569126e-02 -2.203890e-01  1.249238e-01
#>   [71] -5.866814e-02 -2.655910e-01 -6.553552e-02 -2.182685e-02  1.080447e-01
#>   [76] -1.900018e-01  2.075360e-01 -5.485986e-03 -1.201882e-01 -2.739235e-02
#>   [81]  1.814787e-01  1.189319e-01  1.591090e-01 -9.989263e-02 -6.699430e-02
#>   [86] -1.771297e-01 -2.072055e-02  1.149498e-01  3.063352e-02  7.604076e-02
#>   [91] -9.306308e-02 -2.769807e-02 -1.602871e-01 -5.794838e-02  1.275944e-01
#>   [96] -3.257220e-02  1.113319e-01 -7.872929e-02 -1.127296e-01 -1.552363e-01
#>  [101]  1.593381e-02 -1.545007e-01 -2.649819e-01 -1.692872e-01 -1.706113e-01
#>  [106] -2.913737e-02 -1.592336e-01  4.381229e-02 -8.045553e-02 -1.539292e-01
#>  [111] -2.860901e-02  2.934001e-02  6.086181e-02 -2.647766e-02  5.724275e-02
#>  [116]  5.232794e-02 -7.610821e-02 -1.277513e-01 -1.244528e-02  7.132878e-02
#>  [121] -8.636157e-02 -2.071484e-02  1.705737e-02 -1.492090e-01  4.536463e-02
#>  [126] -5.585391e-02  6.899835e-02 -1.232959e-01  6.603253e-02 -1.184759e-02
#>  [131] -1.743025e-01  2.314804e-02 -4.022468e-02 -2.785456e-02  4.050408e-02
#>  [136] -9.934190e-02  2.171804e-01 -2.045169e-01 -1.941796e-02 -6.913283e-02
#>  [141]  2.086559e-01 -8.676395e-02  1.737260e-02 -1.486531e-01 -1.354318e-01
#>  [146]  2.289960e-03  1.861098e-02 -2.247132e-01 -7.914037e-02  1.248409e-01
#>  [151]  9.530066e-02 -4.147217e-02  1.561834e-01 -6.065093e-02  1.694048e-01
#>  [156] -1.178343e-01  8.508711e-02 -2.754444e-02  6.139004e-03 -9.686867e-03
#>  [161] -4.199913e-02  5.175418e-02 -4.930831e-02  8.754921e-02  7.258283e-02
#>  [166]  6.757617e-02  1.288341e-01 -1.130419e-02 -3.259298e-01  4.841093e-02
#>  [171] -1.466936e-01 -5.688797e-03 -2.015105e-02  1.250412e-02  2.702452e-02
#>  [176] -9.710107e-02 -6.371636e-03 -1.223198e-01 -7.802146e-02  8.830279e-02
#>  [181] -7.303914e-02  2.627142e-02  1.857780e-01  1.878457e-01  2.069919e-02
#>  [186] -2.481390e-03 -7.946395e-02  9.623287e-02 -1.217706e-01 -1.146523e-01
#>  [191]  5.247283e-02  1.441803e-01 -8.932047e-02  2.852735e-02 -1.008476e-01
#>  [196] -1.361606e-01 -1.759075e-02 -1.406764e-01  1.179209e-01  5.597716e-02
#>  [201] -6.280777e-02 -7.682682e-02 -2.816082e-02  6.724811e-03  1.345411e-01
#>  [206]  9.314789e-02 -1.038794e-01 -7.072447e-02 -6.976462e-02  9.673360e-02
#>  [211] -2.043259e-01  3.256659e-03 -1.826088e-01  1.177857e-01  4.612112e-02
#>  [216] -6.671027e-02  1.389399e-01  1.065443e-01  2.421870e-02 -9.698133e-02
#>  [221]  1.811137e-01 -1.499441e-01  1.379864e-01 -4.893527e-04 -1.527641e-01
#>  [226]  2.304586e-02  7.796594e-02  1.381460e-01 -8.299122e-02  1.340719e-01
#>  [231]  5.615236e-03  1.074472e-01  7.697950e-02 -1.645455e-01 -1.093193e-01
#>  [236] -5.407425e-02  6.198238e-02 -2.010970e-02  5.291333e-02 -7.473904e-02
#>  [241]  1.991720e-01 -5.684991e-02  7.201492e-02 -1.039388e-01 -7.411644e-02
#>  [246] -8.011595e-02  1.070326e-01 -2.084889e-02  9.577328e-03  3.055940e-02
#>  [251]  4.865815e-02 -3.086232e-01 -1.227354e-01  7.644898e-02 -1.923541e-01
#>  [256] -3.072205e-02  1.352989e-01  4.468437e-02 -1.430626e-01  8.803109e-02
#>  [261] -1.099579e-01 -1.928763e-01  4.226686e-02  1.058786e-01 -2.289812e-01
#>  [266]  1.523058e-02 -5.051754e-03  7.047414e-02  1.550893e-01 -3.480602e-01
#>  [271]  5.227305e-02 -1.464489e-01 -2.768819e-01  1.141339e-01 -9.161603e-02
#>  [276] -1.862292e-01 -1.915579e-01 -7.844770e-02 -4.981603e-02  1.346384e-01
#>  [281] -2.885049e-02 -1.036983e-01 -6.702394e-03 -1.340224e-01  9.540819e-02
#>  [286] -1.199506e-01  5.906287e-02 -4.503586e-02 -1.013343e-01  9.272386e-02
#>  [291]  2.797845e-02  7.645673e-02 -1.716074e-01 -1.106055e-01  6.335215e-02
#>  [296]  2.823762e-02 -1.314983e-01 -7.484573e-02 -3.575589e-02  1.124782e-01
#>  [301] -1.086844e-01  6.415738e-02 -4.781339e-02 -6.258043e-02 -5.033756e-02
#>  [306]  8.742717e-02 -1.059502e-01 -1.172617e-01 -2.408594e-01 -2.781553e-01
#>  [311] -1.603187e-01 -1.439320e-01  1.162045e-01  5.350321e-03  1.337748e-01
#>  [316] -4.094302e-02 -3.643093e-02 -3.791050e-02  1.872112e-01  3.602080e-02
#>  [321] -5.045592e-02  2.266229e-02 -5.302997e-02 -6.333251e-02  5.558800e-02
#>  [326] -8.500159e-02  5.975009e-02 -1.191030e-02 -6.825351e-03  4.926352e-02
#>  [331]  1.366187e-01  7.155332e-02  3.972843e-02 -1.688744e-01 -2.653688e-01
#>  [336] -3.056626e-01  1.795919e-02  2.407528e-02 -7.701386e-02  2.252996e-01
#>  [341]  1.120299e-01 -1.208586e-01 -9.864221e-02  9.587789e-02  3.053211e-02
#>  [346] -1.096561e-02 -5.942350e-02 -1.105156e-01 -9.625725e-02  1.712702e-02
#>  [351]  3.477305e-01  1.590693e-01  9.706704e-02  2.684638e-01 -8.208107e-02
#>  [356] -2.020079e-01 -1.794800e-01 -2.008283e-01 -1.164623e-01  5.566230e-03
#>  [361]  8.630038e-02  4.678658e-03  8.024655e-03  6.494141e-02  7.029320e-02
#>  [366]  6.380819e-02 -6.227227e-02 -5.132362e-03  1.330122e-02 -5.145671e-02
#>  [371]  1.938589e-03  1.932115e-02 -8.234314e-02 -2.803317e-02  1.589671e-01
#>  [376] -1.347644e-02 -3.952847e-02 -1.885654e-01 -1.418549e-01 -1.007108e-01
#>  [381]  2.075798e-01 -2.339739e-01  1.814199e-01  3.597761e-02  8.402621e-02
#>  [386]  1.350181e-01 -6.099647e-02 -9.281649e-02  1.050782e-01 -3.930234e-02
#>  [391] -1.532831e-01  1.672692e-01 -1.358656e-01  7.839177e-02 -7.901414e-02
#>  [396]  5.167682e-02 -1.246776e-02 -1.005826e-01 -2.717917e-01  4.268930e-02
#>  [401] -1.486182e-01  7.580316e-02  1.623683e-01 -1.068913e-01  3.421092e-02
#>  [406] -2.358587e-01  1.303657e-01  3.699048e-02 -1.929321e-01  3.714155e-02
#>  [411] -5.945009e-02 -1.725711e-01  1.191476e-02 -8.605167e-02 -1.377207e-01
#>  [416]  1.424846e-01  9.749248e-02 -2.170999e-01 -1.169235e-01 -1.637948e-01
#>  [421]  3.672734e-02  6.736814e-02  1.277916e-01  2.088986e-01 -5.501902e-02
#>  [426] -1.488017e-01  1.439831e-01  7.707229e-02  1.174393e-02  5.033124e-02
#>  [431] -8.030080e-02 -1.011226e-01 -2.852852e-01  1.132317e-01 -1.275624e-01
#>  [436] -9.290244e-02 -9.787777e-02  7.444549e-02 -2.254602e-02 -2.729374e-01
#>  [441] -3.480034e-02 -9.632960e-02 -2.808553e-02 -5.299763e-02 -2.291601e-01
#>  [446]  1.399560e-02 -7.156623e-03  3.596604e-02 -9.605555e-02 -8.734228e-03
#>  [451] -5.729116e-02 -1.289043e-01 -7.688404e-02  1.057936e-01  6.814309e-03
#>  [456]  1.696264e-01 -4.436978e-02 -2.082104e-01  1.102898e-01  6.217596e-02
#>  [461]  9.977760e-02 -1.237516e-01  5.681590e-02 -5.768743e-02  1.079190e-01
#>  [466]  5.932622e-02  6.806175e-02 -1.289075e-01 -3.878821e-02  5.642819e-02
#>  [471]  2.584972e-02 -2.445815e-01  2.380435e-02 -1.102773e-01 -7.985777e-02
#>  [476]  6.318830e-02  2.813032e-02  3.004252e-02 -8.052568e-03  5.671187e-02
#>  [481]  8.512636e-02  1.478637e-02 -3.940316e-02  6.635966e-02 -6.584152e-02
#>  [486] -4.564469e-02  3.545877e-03  7.699135e-02  8.359937e-02  1.090516e-01
#>  [491] -1.459576e-01  2.323274e-02 -9.593286e-02  5.310866e-02  4.556569e-02
#>  [496] -8.063432e-02 -5.900294e-02 -1.042200e-01 -1.250433e-01 -6.946955e-02
#>  [501]  1.493684e-01  3.384735e-01 -6.329152e-02 -1.313774e-02 -8.308060e-03
#>  [506]  5.080092e-02 -5.176989e-02 -1.068200e-01  7.845257e-03 -1.223237e-01
#>  [511]  9.910068e-03  7.813133e-02  1.214723e-01  1.116253e-01 -1.179269e-01
#>  [516] -4.000495e-03 -7.054576e-02  2.045823e-01 -3.424800e-03 -1.273275e-01
#>  [521] -1.356091e-01 -3.636247e-03  9.893633e-02 -1.220539e-01 -2.199400e-02
#>  [526] -3.176145e-02 -2.474781e-01 -8.391826e-02 -2.906059e-03 -1.185023e-01
#>  [531]  6.764485e-02 -7.910722e-02 -5.150798e-02 -7.958714e-03 -9.206481e-02
#>  [536] -1.286733e-01  2.557280e-02 -1.531153e-02  8.591292e-02 -4.097870e-02
#>  [541] -4.761121e-04 -2.504664e-01  8.653750e-02  1.754913e-01  7.615443e-02
#>  [546] -1.029760e-01 -3.706639e-01  1.498798e-01 -6.340554e-02 -1.234013e-01
#>  [551]  5.485026e-02  2.983690e-02 -4.831603e-02 -1.617058e-01 -9.433242e-02
#>  [556] -3.829029e-02 -1.987693e-01  1.198195e-01 -4.580986e-02  9.476230e-02
#>  [561] -4.530397e-02  2.889738e-01  6.960246e-02  7.016765e-02 -3.735690e-02
#>  [566] -7.296016e-02 -3.701222e-02 -1.228628e-01  9.553993e-02 -7.143788e-02
#>  [571]  9.672845e-02 -1.969456e-01  1.543238e-01 -8.251597e-02  8.145050e-02
#>  [576] -1.009830e-01  1.176341e-02  1.156117e-02 -2.459393e-02 -5.484906e-02
#>  [581]  7.568815e-02  4.192061e-02 -6.790841e-02  1.431580e-01 -2.739788e-01
#>  [586] -1.307053e-01 -1.219170e-01 -6.384632e-02  1.449193e-01 -1.308555e-01
#>  [591]  7.335919e-02  6.528673e-02 -1.047195e-01 -1.432544e-01 -1.070794e-01
#>  [596]  6.153738e-02  9.647086e-03 -3.471587e-02 -1.533816e-01 -1.069363e-01
#>  [601]  1.623943e-01 -7.242706e-02  4.107693e-04  1.289038e-01 -1.899316e-01
#>  [606] -2.165680e-01 -2.582520e-02 -2.434301e-02 -5.500093e-02  2.913787e-02
#>  [611]  5.759734e-02  1.681117e-01 -2.376119e-02 -7.150809e-02  1.452800e-02
#>  [616] -1.243035e-01 -9.469701e-02  1.186727e-01 -1.716757e-01 -1.343840e-01
#>  [621] -1.704523e-02  1.265945e-01 -1.996421e-02 -1.251042e-01 -5.541790e-02
#>  [626]  1.341099e-01  6.966062e-02  1.192324e-01  1.268066e-01  1.681523e-01
#>  [631] -8.235467e-02  3.102049e-02  3.279108e-02 -1.383064e-01 -1.654344e-02
#>  [636]  6.130824e-02 -1.312273e-01 -6.989352e-02 -1.804944e-01 -2.362187e-02
#>  [641]  1.489677e-01  7.310958e-02  3.619657e-02  9.581651e-02  4.714007e-02
#>  [646]  5.244048e-03 -1.194040e-01  7.356175e-02  1.727179e-01 -1.289790e-01
#>  [651]  1.861599e-01  1.183496e-01  7.284460e-02 -1.732725e-01 -5.425221e-02
#>  [656]  1.048679e-01  4.624390e-02  1.351769e-01 -1.058760e-01 -6.575799e-02
#>  [661]  7.501387e-02 -4.856559e-02  6.040197e-02 -2.769830e-01  7.575125e-02
#>  [666] -9.208014e-02 -2.299468e-01 -6.568368e-02  9.130295e-02 -1.231641e-01
#>  [671] -2.339888e-01 -2.763521e-02  5.658829e-02  1.147025e-01  4.445976e-02
#>  [676] -3.923697e-02  1.322015e-02  1.594970e-02 -1.316460e-01  8.834104e-02
#>  [681]  1.247570e-01  2.104714e-01  1.892067e-02 -6.262656e-02 -9.715593e-02
#>  [686]  2.251139e-02 -1.335343e-01  1.648470e-01 -6.979002e-03 -8.385925e-02
#>  [691]  1.018532e-01  1.053007e-01  3.359471e-02 -2.964998e-02  1.266244e-02
#>  [696] -1.559278e-01  9.301535e-02 -3.763176e-02 -5.898233e-02  2.340447e-02
#>  [701]  2.289186e-02  2.205887e-01 -2.549738e-02 -1.876961e-02  2.167689e-02
#>  [706]  1.965729e-01 -5.687600e-02 -1.921487e-01  8.925172e-03  6.200777e-02
#>  [711] -2.205854e-02 -1.072800e-01  1.185936e-02 -1.476877e-01  1.294769e-03
#>  [716]  5.874365e-02  7.643271e-04  1.189678e-02 -2.457743e-02  1.351139e-01
#>  [721]  1.331878e-01  1.749726e-01 -1.671223e-01 -7.723532e-02  6.752663e-02
#>  [726] -5.903702e-02 -2.219220e-02 -1.073689e-01  1.021360e-01  5.511661e-02
#>  [731] -3.223217e-02  1.432760e-01  1.882986e-02  4.604821e-02 -1.156428e-01
#>  [736]  6.839400e-02 -2.453726e-01 -9.006870e-02 -2.447958e-02  7.483612e-02
#>  [741]  5.352085e-02  1.573544e-01 -1.467114e-01  1.050333e-01 -2.310743e-01
#>  [746] -1.143035e-01  6.799687e-03  6.657682e-02 -2.737304e-02 -9.416584e-02
#>  [751] -1.281524e-02  1.211612e-01 -1.189392e-02 -1.806554e-02  9.599043e-03
#>  [756] -1.980381e-01 -4.876373e-02 -1.561210e-01  4.164856e-02 -2.209514e-02
#>  [761]  2.429556e-02  2.680711e-02 -3.932238e-02  3.222073e-02  1.831463e-02
#>  [766] -1.536331e-01  1.851589e-01 -1.792712e-02  5.222296e-02  2.442937e-02
#>  [771] -9.317193e-02 -8.340486e-03 -2.108875e-01 -1.825568e-01  5.336663e-03
#>  [776]  1.067100e-01 -1.501229e-01  1.126889e-01  5.977518e-02  1.053766e-01
#>  [781]  8.073201e-02 -1.147128e-01  1.032546e-01 -3.725183e-02 -2.516074e-02
#>  [786]  1.399421e-01  9.357615e-03 -1.954688e-02  1.551580e-01 -7.649333e-02
#>  [791] -3.050295e-01  1.557727e-01 -1.223413e-01 -1.329499e-01  6.605251e-02
#>  [796]  2.268110e-02  2.430536e-01  1.948208e-01 -1.500858e-01 -6.802994e-02
#>  [801] -1.502486e-01 -2.478595e-02 -7.343045e-02 -8.732498e-03 -1.794639e-01
#>  [806] -3.377509e-02 -8.133957e-02  1.447985e-01  5.353811e-02  1.027636e-02
#>  [811] -8.037245e-03 -2.232777e-02  6.880447e-02  3.289741e-02  5.256994e-02
#>  [816]  1.594602e-01 -1.505394e-01  1.609217e-02  4.122280e-02 -1.567272e-01
#>  [821] -1.432181e-01 -5.944845e-02  7.326839e-02 -5.256941e-02  1.113740e-01
#>  [826]  9.031230e-02 -1.604066e-02  7.352208e-02 -1.128329e-01  2.631196e-02
#>  [831]  1.035333e-01  1.260744e-02  1.484254e-02 -1.902527e-01 -1.186544e-01
#>  [836] -1.547926e-01 -1.937793e-01 -7.698692e-03 -5.230557e-02 -9.907863e-02
#>  [841] -1.126982e-01  3.243284e-02 -2.360955e-01  1.434050e-01  4.471275e-02
#>  [846] -1.892276e-02 -4.022400e-02 -9.617035e-02 -1.756269e-02  2.659255e-02
#>  [851] -5.913398e-05  4.305695e-03  2.280222e-02  1.927091e-01 -1.135237e-01
#>  [856] -1.288150e-01 -6.835754e-02  7.237781e-02 -9.697858e-04  1.660793e-01
#>  [861] -1.826393e-01 -2.697484e-01  1.186406e-01  5.680765e-02 -7.027123e-02
#>  [866] -9.412545e-02  1.963870e-01 -9.799196e-02  5.607400e-02  1.370094e-02
#>  [871]  3.972123e-02  7.459665e-03 -8.606392e-02 -4.817793e-02 -1.109383e-01
#>  [876] -7.476042e-02  6.149099e-02 -1.349757e-02 -1.813633e-01 -1.182102e-01
#>  [881] -1.751772e-01 -4.194690e-02  1.006693e-01  9.661587e-03 -8.081730e-02
#>  [886]  4.226020e-02  1.062543e-01 -7.325066e-02  1.658207e-02 -1.743296e-01
#>  [891] -1.401398e-01 -6.764033e-02  3.764706e-02  2.198361e-01 -7.212461e-03
#>  [896] -1.087054e-01 -6.312980e-02  1.296874e-01  1.729990e-01 -6.591001e-02
#>  [901] -1.267497e-01 -6.830344e-02  1.543127e-01 -1.945217e-02  4.451189e-02
#>  [906]  8.498151e-02 -1.327844e-01  9.631781e-02 -3.984630e-02  7.056021e-02
#>  [911]  1.826229e-01  1.763502e-01 -2.840195e-02 -4.539501e-02 -1.332525e-01
#>  [916]  4.702900e-02 -1.260333e-01 -2.300610e-01 -1.269315e-01 -4.785324e-02
#>  [921]  2.811165e-02  1.359257e-01 -1.300035e-01 -1.283595e-01  2.954876e-02
#>  [926]  1.084943e-01 -1.262012e-01  3.825352e-02 -9.149633e-03 -7.215886e-02
#>  [931]  8.561320e-03  8.609268e-02 -1.047583e-01  5.810729e-02 -1.036067e-01
#>  [936]  2.323003e-02 -1.856505e-01 -5.769223e-02  2.279795e-02 -4.401015e-02
#>  [941]  1.262148e-02  5.359966e-02  1.297674e-01 -3.337500e-02  1.288959e-01
#>  [946]  2.235369e-01  1.623519e-01 -1.725753e-01  3.643696e-02 -4.745091e-02
#>  [951] -4.716278e-02  6.063632e-02 -9.971787e-05 -2.720650e-01 -2.164214e-01
#>  [956]  8.413425e-02 -7.181685e-02  3.505682e-02  2.360624e-02  1.195363e-01
#>  [961]  6.893577e-02  8.425697e-03 -1.328291e-01  1.474127e-01 -8.939313e-02
#>  [966] -6.884847e-03  2.055749e-01 -1.155719e-01 -1.200442e-01  1.896262e-01
#>  [971] -1.101435e-01 -8.768181e-02 -5.390543e-02  9.296211e-02 -1.513809e-01
#>  [976] -1.215620e-02  2.206180e-01  4.672755e-02 -1.230501e-01  5.781073e-02
#>  [981]  2.574642e-01  5.241205e-02  1.971061e-01 -2.493926e-01 -9.117630e-02
#>  [986]  1.208996e-01 -2.933314e-02 -4.341033e-02  3.747059e-02 -5.923422e-02
#>  [991] -1.829948e-01 -1.627337e-01 -1.541309e-01 -3.463124e-02  1.186832e-04
#>  [996]  2.355263e-01  1.241062e-02 -1.013211e-01 -1.193248e-01  3.799072e-02
#> [1001]  1.369424e-01  3.604893e-02 -1.760419e-01 -1.484230e-01  1.131430e-01
#> [1006]  2.094759e-02  8.821695e-02 -1.402152e-01 -5.339269e-02  5.327479e-02
#> [1011] -3.032670e-01 -5.998893e-02  2.435419e-01 -1.505246e-01  1.688909e-01
#> [1016]  5.496127e-02 -8.040999e-02  4.119844e-02  6.042407e-02  1.751916e-02
#> [1021] -1.044723e-01  1.135346e-01  1.340937e-01  5.645961e-02 -9.638623e-02
#> [1026] -3.126826e-02  1.403512e-02  1.137021e-01 -2.612966e-01 -2.739617e-02
#> [1031] -9.145442e-02  1.060787e-01  2.376569e-01 -6.639982e-02 -5.071841e-02
#> [1036]  1.917083e-02 -2.842784e-02 -4.342506e-02  9.037678e-02 -5.048177e-02
#> [1041] -3.058991e-02  3.717159e-02 -1.411122e-01  1.066801e-01  4.649757e-02
#> [1046] -1.055353e-02 -1.348852e-02 -9.731703e-02  4.793558e-02 -2.119511e-01
#> [1051]  1.199473e-02 -1.158612e-01 -6.124406e-02 -1.532452e-01  1.703300e-01
#> [1056]  1.642595e-01 -8.931148e-02 -4.240963e-02  9.530293e-02 -1.129314e-03
#> [1061]  1.266759e-02 -9.603365e-02  5.830467e-02 -2.262409e-02  1.065031e-01
#> [1066] -8.611647e-02 -1.729361e-01  9.009811e-02  7.471964e-02  5.825098e-02
#> [1071]  4.249454e-04 -1.334125e-01 -4.671513e-02 -2.630345e-01  1.331710e-01
#> [1076] -5.392390e-03 -5.428213e-02  1.938564e-02 -5.818909e-03  2.809219e-02
#> [1081]  8.410579e-02  4.196451e-02  4.830873e-02  9.672690e-02  1.162585e-02
#> [1086] -1.507528e-02  1.373717e-01 -1.201120e-01  6.938204e-02  6.911602e-02
#> [1091]  1.317012e-01  2.263254e-02  1.188963e-01  2.502346e-02 -1.057328e-01
#> [1096] -2.581370e-02  1.235097e-03 -1.882112e-01 -1.328556e-01 -4.078286e-02
#> [1101]  5.969230e-02 -4.078269e-02 -4.949430e-02 -2.027018e-02  1.480931e-02
#> [1106] -1.081276e-01 -2.224181e-01  2.276572e-02 -6.186949e-02  9.841463e-02
#> [1111] -1.364965e-01 -1.746667e-01 -1.634958e-02  7.091826e-02  4.129792e-03
#> [1116] -3.576413e-02 -7.193198e-02 -1.135566e-01 -3.446101e-02  1.189597e-01
#> [1121]  9.667834e-02  3.869484e-02  2.204974e-01  7.336888e-02 -7.289224e-02
#> [1126] -1.442422e-01 -2.368490e-02 -1.175513e-01 -1.992322e-01 -5.800522e-02
#> [1131] -1.138499e-01  1.203649e-01  1.048894e-01 -3.328003e-02 -1.023005e-01
#> [1136] -5.495019e-02 -1.826433e-02 -6.048503e-02  3.626076e-02 -1.891947e-04
#> [1141]  1.242013e-04  4.992136e-02  1.123551e-01 -2.310651e-02  1.758631e-01
#> [1146]  5.100052e-04  1.806033e-01 -4.664653e-02 -1.585405e-01 -6.190229e-02
#> [1151] -5.556878e-02 -8.069609e-02  1.932067e-01  6.791343e-02  6.230716e-02
#> [1156]  1.741955e-02  5.292033e-02 -1.389093e-02 -1.035419e-01 -9.084531e-02
#> [1161]  2.114202e-01 -1.915454e-01 -1.539548e-01  1.406150e-01  4.037342e-02
#> [1166]  1.019768e-01  2.481848e-02 -2.636103e-01  2.765955e-01 -1.293064e-01
#> [1171]  2.474995e-02 -5.216863e-02 -3.727103e-01  5.933815e-02 -1.077223e-01
#> [1176] -5.853252e-02  5.547056e-02  9.710274e-02  1.784053e-04  1.086528e-01
#> [1181]  6.383836e-02 -1.885104e-01  6.342151e-02  5.066710e-02 -3.902414e-02
#> [1186]  1.398924e-01 -1.016785e-01  4.646781e-02  5.681070e-02 -7.736319e-02
#> [1191]  2.814808e-02 -3.331661e-01  2.509784e-02 -5.644492e-02  1.005737e-01
#> [1196] -4.826291e-02 -3.069944e-01  9.659863e-02  9.438345e-03  3.904202e-03
#> [1201]  1.083645e-01 -1.144915e-01 -6.525369e-03  2.454237e-01 -5.434291e-02
#> [1206]  3.713484e-03 -2.051894e-01 -4.549375e-02  5.350745e-02  2.540513e-02
#> [1211]  1.896860e-02 -4.834691e-02 -7.044659e-02  3.732529e-02 -5.482540e-02
#> [1216] -1.740574e-01  1.078776e-01  2.935899e-01  5.491179e-02  2.705570e-02
#> [1221] -2.149929e-01  5.679127e-02 -1.719590e-01  7.996255e-02  2.367993e-01
#> [1226] -5.875025e-02  1.398907e-02  3.774680e-02  4.365871e-02  1.286859e-01
#> [1231] -1.474539e-01  1.530243e-01 -1.544915e-01 -1.742534e-01  6.270679e-02
#> [1236]  2.794036e-02 -1.533109e-01 -1.790680e-01  2.880243e-02  8.744543e-02
#> [1241] -5.178637e-02  1.008399e-01  4.019424e-02 -1.048702e-01 -1.274536e-01
#> [1246]  2.227847e-02 -1.143069e-01 -6.349991e-02  1.142051e-01  1.143393e-01
#> [1251]  8.493069e-02 -1.129779e-01  3.382753e-02  1.249078e-01  3.249723e-01
#> [1256] -4.235048e-02 -2.001977e-04  9.674737e-02  4.736293e-02 -2.330761e-01
#> [1261] -3.046024e-01 -8.738548e-02 -1.253783e-01 -1.197653e-01  2.957077e-02
#> [1266]  1.761310e-02 -7.783404e-02  2.719775e-03 -5.438801e-02 -7.596826e-02
#> [1271]  7.287751e-02  1.520783e-01 -1.577467e-02 -2.729153e-02 -1.705954e-02
#> [1276]  7.949768e-02 -1.315161e-01  1.639787e-01 -8.377784e-02 -1.567462e-01
#> [1281] -9.830056e-03 -2.296527e-01 -6.302227e-02 -1.834880e-01 -4.722353e-02
#> [1286]  1.053838e-01  1.475604e-01  7.370173e-02 -4.091423e-02  1.484533e-02
#> [1291]  7.505071e-02 -1.186713e-02 -2.434451e-03 -1.507537e-01  4.577462e-02
#> [1296]  2.194230e-02 -1.755077e-01 -1.011816e-02  7.981765e-02 -1.315254e-01
#> [1301] -2.598287e-01  5.112403e-02  2.283285e-01  1.022226e-01  1.028878e-02
#> [1306] -6.373820e-02  1.932921e-01  4.447560e-02  1.631771e-01  1.412740e-01
#> [1311]  1.248465e-01  4.289056e-02 -8.425410e-02  3.236120e-02  2.891242e-02
#> [1316] -8.068707e-02  1.424518e-02  1.955601e-02 -2.504403e-01  1.558274e-01
#> [1321]  3.066972e-02  1.163172e-01  1.588345e-01  5.584006e-02  4.570301e-02
#> [1326] -1.523806e-01  9.021489e-02 -1.260323e-01 -1.138006e-01  1.806030e-02
#> [1331]  7.565788e-02 -1.602042e-01 -7.123950e-02  2.097319e-02 -8.215662e-02
#> [1336]  4.561752e-02  1.726670e-01 -6.141661e-02  1.082934e-01 -1.073225e-02
#> [1341] -2.098773e-02 -1.589392e-01 -1.530813e-01  1.531159e-01 -1.230542e-02
#> [1346] -7.538227e-03  4.422715e-02 -1.381012e-01 -4.752796e-02  1.925135e-01
#> [1351]  1.415703e-01 -1.625304e-01 -6.794634e-02 -1.605065e-01 -8.437628e-02
#> [1356] -1.108120e-01 -2.230291e-02  2.914308e-01  7.443191e-02  5.173058e-02
#> [1361] -2.603429e-01 -1.146820e-01  6.444216e-02 -9.746858e-02 -1.113125e-01
#> [1366] -1.404501e-02 -1.851067e-01  1.123579e-01 -1.331226e-01 -9.358122e-02
#> [1371]  5.640833e-02 -1.347628e-01 -4.994273e-02 -1.951089e-03 -3.806722e-02
#> [1376]  3.447120e-02  5.724413e-02 -4.149922e-02 -9.952513e-02 -1.652219e-01
#> [1381]  1.044143e-02 -2.114595e-01 -1.518242e-01 -5.271489e-02 -9.346881e-02
#> [1386]  1.772796e-01  7.544006e-02 -1.943983e-01  2.333521e-01  2.692020e-01
#> [1391] -5.460158e-02 -9.322069e-02 -3.109810e-01  1.527048e-02 -6.698630e-02
#> [1396] -7.434852e-02  1.664663e-02  6.387016e-02 -8.249527e-02  7.179927e-02
#> [1401]  1.012397e-01 -1.654466e-01  1.267453e-01 -1.182566e-01 -1.307949e-01
#> [1406]  1.464682e-01  1.734470e-01  8.898780e-02  2.437623e-01  7.855693e-03
#> [1411]  2.562971e-02  9.573185e-03  1.856418e-01  5.831913e-02  7.161483e-02
#> [1416] -1.143788e-02 -2.476529e-01 -1.103685e-01  1.440730e-01 -3.342727e-02
#> [1421] -7.251283e-02  2.127530e-02 -8.115681e-02  1.926080e-01  1.774946e-01
#> [1426] -1.554586e-01  2.601013e-02  1.115441e-01 -9.985012e-02  2.089246e-01
#> [1431] -1.142794e-02 -2.155463e-03 -4.238577e-02  6.429090e-02 -1.314541e-01
#> [1436]  4.932470e-02 -8.276731e-02  5.825110e-02  6.038288e-02 -1.688657e-01
#> [1441] -1.865300e-02 -1.279916e-01 -1.008466e-01 -1.547207e-02 -4.671955e-02
#> [1446]  3.037518e-01 -5.102149e-03 -1.439303e-02  1.549427e-01 -7.145647e-02
#> [1451] -2.126099e-01 -4.656933e-02  1.975654e-02 -1.924720e-02 -7.834157e-02
#> [1456]  6.785937e-02  4.893689e-02  1.550524e-01  7.456381e-02  4.164038e-02
#> [1461] -6.847542e-02  5.102276e-02 -3.710217e-02  1.092604e-01 -1.256414e-01
#> [1466] -1.007469e-01 -1.250013e-02  2.025288e-01 -2.066691e-02  1.201939e-01
#> [1471] -9.605301e-02  2.510454e-02  1.807522e-02 -3.893278e-02  9.537044e-02
#> [1476] -9.904207e-02 -6.381242e-02 -3.447305e-01 -4.889392e-02  2.934721e-01
#> [1481]  1.107268e-01  1.345706e-01  4.008285e-02 -2.908885e-03 -5.799915e-02
#> [1486]  5.286230e-02  2.262926e-02 -3.817742e-02  1.448358e-02 -1.153617e-01
#> [1491] -1.336612e-01  5.825051e-02  1.624620e-01 -3.678155e-03  8.473352e-02
#> [1496] -5.273775e-02 -1.160406e-01  1.418822e-01 -2.669698e-02  8.991584e-03
#> [1501]  7.510334e-02  1.052328e-01 -7.599959e-02 -5.081494e-02 -1.509163e-01
#> [1506] -6.159537e-02  2.476224e-02 -3.083624e-02  1.764982e-01  3.041921e-02
#> [1511]  2.202163e-02 -9.378032e-02  1.265123e-01  7.391792e-04  9.448539e-02
#> [1516] -1.433898e-01 -1.558090e-01 -2.427115e-01  5.792383e-02  9.856947e-02
#> [1521] -1.237695e-01 -1.963700e-01 -4.364780e-03  6.180684e-02 -5.583711e-02
#> [1526] -5.175461e-02  3.212842e-03  5.007793e-02  1.237644e-01  4.101028e-02
#> [1531] -5.326368e-02 -8.253973e-03 -1.059788e-01 -1.412484e-01  3.657370e-02
#> [1536]  9.709103e-02 -1.912148e-01 -2.125636e-01 -1.084275e-01  9.552975e-03
#> [1541] -1.398063e-01 -3.640459e-02 -8.914200e-02  2.078792e-01  2.032519e-01
#> [1546]  4.231926e-02 -1.107075e-02 -8.640799e-02 -1.091232e-01  1.260510e-01
#> [1551]  5.131321e-02  9.051936e-02 -1.358313e-01  4.799234e-02 -5.475406e-02
#> [1556]  1.326552e-01 -1.090257e-01 -8.869870e-02 -2.011761e-01  4.247283e-03
#> [1561] -1.447114e-01  4.227616e-02 -1.675616e-02  1.698835e-01  1.455855e-01
#> [1566] -2.144949e-01 -1.168097e-01 -1.055732e-02  8.067266e-02 -5.857673e-03
#> [1571]  6.786284e-02  1.102699e-01  7.191901e-02 -8.428307e-02  3.605551e-02
#> [1576]  8.866196e-02 -4.287362e-02 -4.281463e-02 -2.145115e-01 -2.826756e-01
#> [1581]  3.136460e-02  1.030898e-01  7.869522e-02 -9.749341e-02  9.876256e-02
#> [1586]  2.217943e-01  1.446468e-01  3.143589e-03  5.710788e-02 -2.542822e-02
#> [1591] -5.349887e-02  7.405415e-02  1.092856e-01  6.225944e-02 -7.616404e-03
#> [1596] -2.541887e-01 -1.147255e-02  3.515335e-02 -1.204825e-01  6.487002e-02
#> [1601] -6.895794e-02 -6.168394e-02  2.527559e-02  1.990653e-02  3.443641e-01
#> [1606]  1.069038e-01  1.004806e-02 -4.093868e-02  2.627586e-02 -7.382326e-02
#> [1611] -1.499246e-01  1.188662e-01  1.416086e-02 -1.306692e-01 -4.592146e-02
#> [1616] -4.847261e-02  1.865333e-01  5.111020e-02  4.055895e-02 -7.423709e-02
#> [1621] -1.573546e-01 -1.839706e-01 -8.974954e-02 -1.124631e-01  1.049444e-01
#> [1626]  5.083257e-02 -1.316121e-02 -3.091905e-02 -5.673464e-02  7.817003e-02
#> [1631] -2.981131e-02  3.571782e-02 -3.291512e-02 -1.803129e-01 -7.138460e-02
#> [1636]  8.109310e-02  1.141195e-01 -1.233709e-01  1.100343e-02  2.221686e-02
#> [1641]  7.322816e-02  3.774946e-02  1.344445e-01  1.994253e-02 -4.026732e-02
#> [1646] -2.605049e-01 -1.778825e-01 -9.441685e-02 -1.581544e-01  8.074403e-02
#> [1651] -2.395007e-01  7.180293e-02 -3.805822e-02  2.689779e-02 -1.029842e-02
#> [1656] -3.080348e-01 -4.913997e-02  5.718839e-02  1.003110e-01  1.693378e-02
#> [1661] -9.001733e-03  5.796468e-02 -1.480821e-01  7.670530e-02  8.959790e-03
#> [1666]  1.349502e-01 -1.095124e-01  6.502110e-02 -2.595404e-02 -9.312161e-02
#> [1671] -5.171271e-02  1.981835e-02  9.873647e-02 -1.831473e-02 -7.678288e-02
#> [1676]  1.013947e-01 -1.291342e-01  2.527365e-02  1.257956e-01  2.509152e-01
#> [1681] -1.584036e-01 -1.163422e-01 -1.893735e-01 -1.953208e-01 -7.006808e-02
#> [1686] -5.140613e-02  8.570273e-02 -2.390494e-04 -1.489117e-01 -5.290608e-02
#> [1691]  5.038683e-02 -1.768487e-01  2.513832e-01 -1.539446e-01 -8.453659e-02
#> [1696]  1.740095e-01 -2.484244e-02  9.864678e-03  9.936955e-02  1.461923e-02
#> [1701]  4.578229e-02 -9.137279e-02 -5.849678e-02 -5.739538e-02 -7.757292e-03
#> [1706] -3.352973e-02 -7.247597e-02 -8.047665e-02  1.073828e-01 -2.626529e-01
#> [1711]  2.653464e-02 -1.217606e-01 -5.380440e-04  1.269976e-01  1.156238e-01
#> [1716] -7.415223e-02 -7.057637e-02 -9.382997e-02  1.982323e-02 -9.811327e-02
#> [1721] -5.530877e-02  1.807827e-01 -1.256755e-01 -9.322158e-02 -6.724150e-02
#> [1726] -1.727246e-01  1.873234e-01 -3.293805e-01  3.816417e-02 -8.418827e-02
#> [1731] -1.347772e-01  3.160425e-03  1.001668e-01  6.310030e-02  2.140357e-01
#> [1736] -7.247012e-03 -1.651633e-02 -2.583234e-02 -2.492794e-02 -8.934765e-03
#> [1741] -2.790394e-02  1.029951e-01 -1.537856e-01 -4.576998e-02  1.631377e-01
#> [1746] -1.117880e-01  2.841798e-02  1.451499e-01  8.462463e-02 -2.863458e-02
#> [1751]  3.250996e-02  3.806386e-02  1.624082e-01  8.428094e-03  1.061800e-01
#> [1756] -9.742088e-02 -2.180195e-01  2.458781e-01  1.792709e-02 -3.330550e-02
#> [1761] -4.974843e-02 -5.055770e-02  3.633658e-02  5.816920e-02  8.575853e-02
#> [1766] -5.422917e-02  1.223269e-03  1.245501e-01 -1.563808e-03 -2.476664e-01
#> 
#> $mod1_unimodal$ref_properties$feature$prop_outliers_feature
#> [1] 0
#> 
#> 
#> 
#> $mod1_unimodal$sim_properties
#> $mod1_unimodal$sim_properties$cell
#> $mod1_unimodal$sim_properties$cell$library_size
#>  [1] 223 201 188 186 209 207 256 279 202 199 217 217 192 210 146 201 259 218 217
#> [20] 265 205 239 211 222 228 252 160 227 251 196 239 215 253 261 177 267 255 194
#> [39] 237 244 232 221 288 322 235 296 237 286 286 211 206 255 330 332 246 258 302
#> [58] 256 231 250 321 294 258 250 246 293 244 325 231 268 263 328 264 251 287 183
#> [77] 261 262 213 239
#> 
#> $mod1_unimodal$sim_properties$cell$log_library_size
#>  [1] 5.411646 5.308268 5.241747 5.231109 5.347108 5.337538 5.549076 5.634790
#>  [9] 5.313206 5.298317 5.384495 5.384495 5.262690 5.351858 4.990433 5.308268
#> [17] 5.560682 5.389072 5.384495 5.583496 5.327876 5.480639 5.356586 5.407172
#> [25] 5.433722 5.533389 5.081404 5.429346 5.529429 5.283204 5.480639 5.375278
#> [33] 5.537334 5.568345 5.181784 5.590987 5.545177 5.273000 5.472271 5.501258
#> [41] 5.451038 5.402677 5.666427 5.777652 5.463832 5.693732 5.472271 5.659482
#> [49] 5.659482 5.356586 5.332719 5.545177 5.802118 5.808142 5.509388 5.556828
#> [57] 5.713733 5.549076 5.446737 5.525453 5.774552 5.686975 5.556828 5.525453
#> [65] 5.509388 5.683580 5.501258 5.786897 5.446737 5.594711 5.575949 5.796058
#> [73] 5.579730 5.529429 5.662960 5.214936 5.568345 5.572154 5.365976 5.480639
#> 
#> $mod1_unimodal$sim_properties$cell$tmm_factor
#>  [1] 1.2241491 1.1990518 0.9425501 1.0253628 1.2142413 1.0050382 0.9533712
#>  [8] 1.2313072 1.0474625 1.0952495 0.7116900 1.2718076 0.8977008 0.9640806
#> [15] 1.4798378 0.9402670 1.0063943 0.9897486 0.9612711 0.8468148 1.0066578
#> [22] 0.9163983 1.3748273 1.2242553 0.9583008 0.9534466 1.2550834 0.8104542
#> [29] 1.0836651 1.0800008 1.1847363 0.9071323 0.9402444 0.8052025 1.1425884
#> [36] 1.2491647 1.0572792 0.8133728 0.7413009 0.9002491 0.7044330 1.1903803
#> [43] 0.8756284 0.9446347 1.2420449 0.9878424 0.8612045 1.3137172 0.8311057
#> [50] 0.7868001 1.2071196 0.8432905 0.8579475 0.7449676 0.8136617 0.9587171
#> [57] 1.0213516 1.0823555 1.3359564 1.0951720 0.9878042 1.0495024 0.9611525
#> [64] 0.7510275 1.0280683 0.7722798 0.9277693 0.7036136 1.3899001 1.2250927
#> [71] 1.4115454 0.7664577 1.7306634 1.0351209 1.0131255 0.9348562 0.8023686
#> [78] 0.8341818 1.0276488 1.0004524
#> 
#> $mod1_unimodal$sim_properties$cell$effective_library_size
#>  [1] 272.9852 241.0094 177.1994 190.7175 253.7764 208.0429 244.0630 343.5347
#>  [9] 211.5874 217.9547 154.4367 275.9823 172.3586 202.4569 216.0563 188.9937
#> [17] 260.6561 215.7652 208.5958 224.4059 206.3649 219.0192 290.0886 271.7847
#> [25] 218.4926 240.2685 200.8134 183.9731 271.9999 211.6802 283.1520 195.0334
#> [33] 237.8818 210.1578 202.2381 333.5270 269.6062 157.7943 175.6883 219.6608
#> [41] 163.4285 263.0740 252.1810 304.1724 291.8806 292.4013 204.1055 375.7231
#> [49] 237.6962 166.0148 248.6666 215.0391 283.1227 247.3292 200.1608 247.3490
#> [57] 308.4482 277.0830 308.6059 273.7930 317.0852 308.5537 247.9773 187.7569
#> [65] 252.9048 226.2780 226.3757 228.6744 321.0669 328.3248 371.2364 251.3981
#> [73] 456.8951 259.8154 290.7670 171.0787 209.4182 218.5556 218.8892 239.1081
#> 
#> $mod1_unimodal$sim_properties$cell$zero_fraction_cell
#>  [1] 0.1666667 0.2166667 0.1500000 0.2166667 0.2666667 0.1833333 0.1833333
#>  [8] 0.1833333 0.2166667 0.2833333 0.1833333 0.2166667 0.1833333 0.1666667
#> [15] 0.2666667 0.2166667 0.2666667 0.1833333 0.2666667 0.1666667 0.2166667
#> [22] 0.1333333 0.1833333 0.2333333 0.2333333 0.2166667 0.3333333 0.1166667
#> [29] 0.1333333 0.2500000 0.1666667 0.1833333 0.2166667 0.1166667 0.2000000
#> [36] 0.1500000 0.1666667 0.1833333 0.1666667 0.2166667 0.2166667 0.2666667
#> [43] 0.2166667 0.1166667 0.2333333 0.2166667 0.1500000 0.2000000 0.2166667
#> [50] 0.1666667 0.1833333 0.1833333 0.1333333 0.1666667 0.1833333 0.2000000
#> [57] 0.2333333 0.2000000 0.2500000 0.2000000 0.2000000 0.1166667 0.1333333
#> [64] 0.1833333 0.1500000 0.1333333 0.1833333 0.1666667 0.2500000 0.1833333
#> [71] 0.2000000 0.1166667 0.2500000 0.2166667 0.2000000 0.1833333 0.1500000
#> [78] 0.2166667 0.2333333 0.1833333
#> 
#> $mod1_unimodal$sim_properties$cell$detection_freq_cell
#>  [1] 0.8333333 0.7833333 0.8500000 0.7833333 0.7333333 0.8166667 0.8166667
#>  [8] 0.8166667 0.7833333 0.7166667 0.8166667 0.7833333 0.8166667 0.8333333
#> [15] 0.7333333 0.7833333 0.7333333 0.8166667 0.7333333 0.8333333 0.7833333
#> [22] 0.8666667 0.8166667 0.7666667 0.7666667 0.7833333 0.6666667 0.8833333
#> [29] 0.8666667 0.7500000 0.8333333 0.8166667 0.7833333 0.8833333 0.8000000
#> [36] 0.8500000 0.8333333 0.8166667 0.8333333 0.7833333 0.7833333 0.7333333
#> [43] 0.7833333 0.8833333 0.7666667 0.7833333 0.8500000 0.8000000 0.7833333
#> [50] 0.8333333 0.8166667 0.8166667 0.8666667 0.8333333 0.8166667 0.8000000
#> [57] 0.7666667 0.8000000 0.7500000 0.8000000 0.8000000 0.8833333 0.8666667
#> [64] 0.8166667 0.8500000 0.8666667 0.8166667 0.8333333 0.7500000 0.8166667
#> [71] 0.8000000 0.8833333 0.7500000 0.7833333 0.8000000 0.8166667 0.8500000
#> [78] 0.7833333 0.7666667 0.8166667
#> 
#> $mod1_unimodal$sim_properties$cell$cell_cor
#>    [1]  0.1549533566 -0.0198304203 -0.0258113817  0.1128361658  0.1158512527
#>    [6] -0.1183800807  0.2092394714 -0.0404089783 -0.1973127648  0.1682957942
#>   [11]  0.1522498625 -0.1159468775 -0.1549242927  0.2202216389  0.1638476210
#>   [16]  0.0039824945 -0.1366187188  0.0936969298  0.0954442689  0.1502587042
#>   [21] -0.0180348086 -0.1469047349 -0.0941205967  0.0313927915  0.0159499616
#>   [26]  0.1782355100 -0.0425277267  0.2969591374  0.0100829907  0.0805502342
#>   [31]  0.0155274650  0.0139531054 -0.0232672358  0.0792147739  0.0806785161
#>   [36] -0.1686707080  0.1165996963  0.1828066172 -0.1735558844 -0.1135997999
#>   [41]  0.0977671936  0.1134432623 -0.1957434826 -0.1852552893  0.0880754631
#>   [46]  0.1058930794  0.0641066449  0.0652647038  0.1836586714 -0.0849809748
#>   [51]  0.1147719244 -0.0197079953 -0.0260440895  0.3452542897 -0.0963281783
#>   [56]  0.1037358119 -0.0446622865  0.0888395511 -0.1018596144  0.0120186796
#>   [61]  0.1594479140 -0.1492267380 -0.0662460119  0.1353202106  0.0017785670
#>   [66]  0.1185023647  0.1086227827 -0.0588433720 -0.0646723636  0.1269078757
#>   [71]  0.1730598995 -0.0247255437 -0.0367377847  0.2281443002  0.0652575250
#>   [76]  0.0946062362  0.2413679763  0.0914000909 -0.0084578963 -0.1098226733
#>   [81]  0.0359958796  0.1143246060  0.0980784869  0.1933100095  0.1066347731
#>   [86]  0.4734626459  0.0388843320 -0.1195719060 -0.0562205740  0.0777573253
#>   [91]  0.1144801812  0.0377803090  0.0084881830 -0.0340245114  0.0540148919
#>   [96] -0.1070349248 -0.0515033123  0.1464625265  0.2058278763 -0.1789211034
#>  [101]  0.0477964115  0.0749585101 -0.1022657043  0.1350751288  0.0023967498
#>  [106] -0.1772321241 -0.1705197217  0.1623641208  0.0045024762 -0.1073818749
#>  [111]  0.0777070180  0.1248953086  0.0961040814  0.1069360838  0.1649716761
#>  [116] -0.0292701639 -0.0127821378  0.0831651276  0.0900150756 -0.1200002505
#>  [121]  0.0446236342  0.0007838572  0.2484746777  0.0148558020  0.0595137712
#>  [126] -0.2456005129 -0.0655864685  0.1279718259 -0.1696288957 -0.0820338414
#>  [131] -0.1635097913  0.2693459388  0.2023643548  0.1105540560 -0.1288505548
#>  [136]  0.0392622926 -0.0959045671 -0.0192952187  0.1125350758  0.0480843143
#>  [141]  0.0381369169  0.3905428049  0.1918408330  0.0414291057 -0.0358192394
#>  [146]  0.1970863976  0.1839007741  0.0397485117 -0.0371419454 -0.1143659575
#>  [151] -0.0419104400  0.2060661511 -0.0464394891 -0.0452981549 -0.0563388725
#>  [156]  0.0248998918  0.2701685629  0.1046340172  0.0230800639  0.0874828729
#>  [161]  0.2053174102  0.0185455791  0.0923902478 -0.0393211816  0.1391708978
#>  [166]  0.0928616885  0.2102404172 -0.2810405623  0.1035141256  0.3288785784
#>  [171]  0.1404637350  0.0723397670  0.0332944474  0.0288713702 -0.0640915134
#>  [176]  0.0320953138 -0.0541541296  0.2039899337  0.1437262673 -0.1930908741
#>  [181] -0.0319554360 -0.2133909304  0.1583856478 -0.1735968778 -0.1368846526
#>  [186]  0.1953084371  0.0429214696  0.0994006239  0.1068811570  0.0182103678
#>  [191]  0.0441416264  0.2587860478 -0.0293356108  0.0879208302  0.0039964779
#>  [196]  0.0183183118  0.0189373500  0.1732899529  0.2026967723  0.1522415430
#>  [201]  0.1383079339  0.0452771931 -0.0855148977  0.1122527123 -0.0331054691
#>  [206]  0.2288179135 -0.1107791731  0.1206430910 -0.0401080844  0.0374611923
#>  [211]  0.1530150646  0.2165440790  0.1944903182  0.1677731699  0.0112494721
#>  [216] -0.1321109923  0.0161644147 -0.0538948531 -0.0085926075 -0.0105698797
#>  [221] -0.1250517663 -0.0345930458 -0.0617404615 -0.0066085967  0.2761577802
#>  [226]  0.1056753146  0.2028629032 -0.1071950035  0.0879412947  0.1946213744
#>  [231]  0.0048109771  0.0176155981  0.0386265336  0.1520121746  0.0356225165
#>  [236]  0.0601774920  0.0488484408  0.1383053325  0.0305936011  0.0462237926
#>  [241] -0.1061207377  0.0522812173  0.2953802518  0.0984622460  0.2116615291
#>  [246]  0.2438969988  0.0668032670  0.1538924604 -0.0580858842 -0.0338040657
#>  [251]  0.1659877394  0.2251101581 -0.0258467687  0.0320770566 -0.0289544547
#>  [256] -0.1354585993  0.2049701870  0.0215164322  0.1133824677  0.0825639182
#>  [261] -0.0318933507 -0.0785330735 -0.0753814855 -0.1000823385  0.2742967273
#>  [266]  0.2365836755  0.0023847380  0.0681281793  0.0561441922  0.1487550798
#>  [271]  0.0607853120 -0.0284038082  0.0181944887  0.0006113081 -0.0110341770
#>  [276]  0.1046940270 -0.1292049307  0.1184128271 -0.1032311220 -0.0149565235
#>  [281] -0.0816206586  0.1001486611 -0.1157926667  0.1093410473  0.1427316369
#>  [286]  0.1391728276  0.1401459928  0.1699038080  0.1941767696  0.4019793608
#>  [291] -0.0844006365  0.2321532211  0.1249026249 -0.1754284276  0.0760531420
#>  [296] -0.0857701405  0.0880525790  0.0585997975  0.0360652558 -0.0119157470
#>  [301]  0.1119395597  0.1674029735  0.1378297524  0.0316875487  0.0761843247
#>  [306]  0.0355698935  0.1226327860  0.0338665697  0.1250606217  0.1263489744
#>  [311] -0.1723982541 -0.0438445056  0.0777766743  0.2119157528 -0.0670605183
#>  [316]  0.2049945871  0.1793015377 -0.0423017604  0.1764091581  0.0058554804
#>  [321]  0.0314873669  0.1269680584  0.2225382411 -0.0369698992  0.0125730983
#>  [326] -0.1984756847  0.1691102881 -0.0046836733  0.0480261758 -0.1051549521
#>  [331] -0.0676170385  0.0794436301  0.2261767486 -0.1385888136  0.1289637304
#>  [336]  0.0219426222 -0.0384720671  0.0709519924  0.0841300267  0.2037709546
#>  [341]  0.2625002027  0.1164671772  0.3109112274  0.2684359673  0.1880111822
#>  [346]  0.1301710943  0.2139239676 -0.0684822910 -0.0390119310  0.1238886601
#>  [351] -0.1385503992 -0.0459653322 -0.1178903713 -0.1355006376  0.0357565860
#>  [356]  0.1539142515  0.0053247084  0.2270682201  0.3402724879 -0.1132795481
#>  [361] -0.1417698443 -0.1822507538  0.1337691083  0.3132457891  0.2550539088
#>  [366]  0.0470125526  0.0981162484  0.3294618160  0.0640491364 -0.0156689084
#>  [371]  0.0636418979 -0.0047632225 -0.2250504005  0.2716748927  0.2296502726
#>  [376]  0.0940944778 -0.0988070116  0.0300126704 -0.1590362471  0.0756970297
#>  [381]  0.0953843292  0.1091238766 -0.1583011870  0.0794090878 -0.1180127906
#>  [386]  0.2302442418 -0.2220227186  0.0666931126 -0.1811912306  0.1419228122
#>  [391]  0.0492340789  0.0203171017  0.0662332566  0.1845237451  0.2376846565
#>  [396]  0.2257146974  0.1204318818  0.1030069090  0.1231239372  0.0338207693
#>  [401] -0.0604907179  0.1890427911  0.0849663236 -0.0363626074  0.2602796070
#>  [406]  0.1576931939 -0.0663008811 -0.0864062428 -0.1036563728  0.0307347774
#>  [411] -0.0476287190 -0.2665297179  0.0553367657  0.0353329679 -0.0797856797
#>  [416]  0.0133165720 -0.2463047611 -0.2215480942 -0.0632228667 -0.1970023502
#>  [421] -0.2023754279  0.0780784970  0.0443262582 -0.1377501945  0.3529309907
#>  [426]  0.1282121613 -0.0949607528 -0.0272106261 -0.2369140597 -0.1784987143
#>  [431] -0.0622554597 -0.1918283363  0.0382654337  0.0449544191 -0.0856012925
#>  [436]  0.1734136767 -0.0101718855 -0.1395883978  0.2332209746  0.3093308617
#>  [441]  0.0363159499  0.1497272722  0.0720192748 -0.0959190595 -0.1354162050
#>  [446] -0.0190112267 -0.0590281735 -0.0038595715 -0.0857281108  0.0120000234
#>  [451] -0.1926366708  0.0560170105  0.0099713547  0.0418857510  0.2024506268
#>  [456]  0.0331793801  0.0622644571 -0.0506249797  0.2170964941 -0.0409773888
#>  [461]  0.0534454649 -0.2253272059  0.0669366145  0.0675507230 -0.0632815155
#>  [466]  0.0190891609  0.2275374628  0.1262479876  0.0933021161  0.0490492574
#>  [471]  0.0424075102  0.1733592264  0.2654841202 -0.0612364754 -0.0291514580
#>  [476] -0.0970713516 -0.0690717557 -0.0969767522  0.1869263238 -0.0711551023
#>  [481]  0.0441882776  0.1047820710  0.2427664120  0.2433066595  0.1049009282
#>  [486]  0.2173536347  0.0364750423  0.1005305861 -0.1897257325 -0.0547385285
#>  [491]  0.2206204033  0.1926403778  0.1123611927  0.0788040613  0.1974874338
#>  [496] -0.3079604991  0.0743129534  0.2944640354  0.1715151833  0.0502954521
#>  [501] -0.0183115060 -0.0330440729  0.0530463036  0.2313819730 -0.1236924761
#>  [506] -0.0056581505 -0.1409437082  0.0861215211  0.1199562021  0.1343202946
#>  [511] -0.0001424295 -0.1157362502  0.2232835573  0.0333055126  0.1884959629
#>  [516] -0.0007635423  0.0785188273 -0.0284625398  0.0444894003  0.0446262585
#>  [521] -0.0467535999  0.2586850845 -0.0721511338  0.2838732211  0.1702335380
#>  [526]  0.0197477394 -0.2027529927  0.3167503676 -0.0102685197  0.1109099392
#>  [531] -0.1521731405  0.1566092082  0.3451499107  0.2136708269  0.0318091455
#>  [536]  0.1237545150 -0.0352695834  0.0370962332 -0.0721560252  0.1783083848
#>  [541]  0.1260733160  0.1058592080 -0.1196033032 -0.1068753856  0.0749433133
#>  [546]  0.0849712610  0.1947286038  0.2166590795 -0.0605465879  0.0425214136
#>  [551]  0.1042395166  0.0524330360  0.1035553101  0.2083512796  0.0526788517
#>  [556]  0.0691343911 -0.0687761374  0.0350268995  0.1401082730  0.1592764256
#>  [561]  0.0197842424  0.0316893994  0.0880342714 -0.0344313795  0.1385871526
#>  [566]  0.0652314639 -0.0519835184  0.0248349302  0.1350661364  0.0130836185
#>  [571] -0.1320465686 -0.1310753666  0.0055891009  0.1442432269  0.1692759221
#>  [576] -0.1290355006 -0.1323470116  0.1421340991 -0.0825121655  0.3549582722
#>  [581] -0.1150235139 -0.0582070196  0.1771572730 -0.0189821281  0.1847525435
#>  [586]  0.0284624180  0.2800235800 -0.0187468144  0.1457212124 -0.0819354131
#>  [591]  0.1193065570  0.2101563648  0.0690521970  0.1562480217  0.0520248722
#>  [596] -0.3684735904 -0.1473323151 -0.0947682847  0.1622470914  0.0474106821
#>  [601]  0.1912753405  0.2761994740  0.1030031582  0.0938392669  0.0250417391
#>  [606]  0.0961319116  0.0789029457  0.1351696410  0.0400273465  0.1687208892
#>  [611]  0.0776793712 -0.0555104343  0.2825339930  0.0803876174  0.0365125637
#>  [616]  0.0431621544 -0.0911303919  0.0866188729  0.0293495231  0.1724693726
#>  [621] -0.1043251358  0.1861520559  0.1665820482  0.3549066611 -0.2664029883
#>  [626] -0.0581742860 -0.0630859018  0.1203085210  0.0565598737 -0.1537613142
#>  [631]  0.0446600536  0.0489227961 -0.0886764086 -0.0839342170 -0.1550470023
#>  [636] -0.1013493322  0.3014910428 -0.0209286192 -0.1184721819 -0.0360592376
#>  [641] -0.1481200403 -0.2300883264 -0.0315492935 -0.0651887524  0.0242652950
#>  [646] -0.0788342169  0.0430699674 -0.0080205870 -0.0761327227 -0.0950380133
#>  [651] -0.2242022294  0.0636353283 -0.2130984826 -0.0276841894 -0.1531845420
#>  [656] -0.0282001284 -0.0490351620  0.1787506060 -0.3103531371  0.0442625027
#>  [661] -0.2042941216  0.0380982368  0.2582696214 -0.0683071970  0.0643774116
#>  [666] -0.0829387632  0.1349823354  0.2360823471 -0.0551264514 -0.1821257558
#>  [671] -0.0487557425 -0.0976020770  0.0566575409 -0.1671154670  0.0384523449
#>  [676]  0.1385598137 -0.0390415563  0.1826140023 -0.0398109667 -0.0932446950
#>  [681]  0.1677240712 -0.1793287082 -0.1539682952 -0.0956007785  0.0723408848
#>  [686]  0.1399007587 -0.0912476203  0.1632201710 -0.1072020036  0.0738475457
#>  [691] -0.0315730867  0.1650785527  0.1132891742 -0.1709144598 -0.0125264554
#>  [696]  0.0767216817 -0.0419591513  0.0753529656  0.1554796267 -0.0461762446
#>  [701]  0.0767404073 -0.0886567956 -0.2193796913  0.1149651379  0.1010497425
#>  [706]  0.1020742242 -0.0635104662 -0.0275222857  0.0200054320  0.1566583637
#>  [711]  0.3305646911 -0.0461175920 -0.2201982196 -0.0003546049  0.0080671491
#>  [716]  0.0343860535  0.3240064150  0.1044687210 -0.0282737674  0.1621753133
#>  [721] -0.0121180429  0.1453813219  0.0244759622 -0.0242006781  0.2782708753
#>  [726]  0.0159279147 -0.0128334362  0.1689638065  0.0380174126  0.0100044784
#>  [731] -0.0825868418  0.0794060998  0.0852561860  0.1013765252  0.1571150753
#>  [736]  0.0277648347  0.2845496357  0.0665725881 -0.1376200477  0.0536826627
#>  [741] -0.1183438136 -0.0560394341 -0.0197580665 -0.0254822516 -0.0795984079
#>  [746] -0.1027324720 -0.0031070293 -0.1854980630 -0.1187315401  0.1117674610
#>  [751]  0.0057938817 -0.1349711426  0.1171937833  0.0750292429  0.0008815354
#>  [756] -0.0587187375  0.0176488194  0.0351611599 -0.0646777422 -0.1526655239
#>  [761]  0.0967433386  0.0155500210 -0.0212359907  0.0504383674  0.1033772888
#>  [766] -0.0819281753  0.1030136449 -0.1574017142  0.1416122943 -0.1124401808
#>  [771] -0.0087163983 -0.1628291287 -0.0332269031  0.2233436557  0.0173931408
#>  [776]  0.1410133421  0.0545294530  0.1600558478 -0.0503260902 -0.2089730713
#>  [781]  0.0605972857  0.1828734282  0.0561925641 -0.0183604480 -0.0997635372
#>  [786] -0.0048770438  0.1068617598  0.0768576475 -0.1350997859 -0.1075826538
#>  [791]  0.0352816425  0.1077561415 -0.2610095287 -0.0740399480  0.0364039422
#>  [796] -0.1155554369  0.2190562977  0.1146044130 -0.0624137616  0.1900894822
#>  [801]  0.0379087572  0.0252348428  0.2367048361 -0.0233705416 -0.0191685560
#>  [806]  0.0638650616 -0.1370250267  0.1688012657  0.1761549702 -0.1030506523
#>  [811]  0.0833048157  0.2200477019  0.0481299399  0.0550186155 -0.0462428663
#>  [816]  0.0954131026  0.0935213938 -0.1552789539 -0.0227214808 -0.0376767524
#>  [821]  0.0817451280 -0.0812031012  0.1331279316  0.0160131374  0.0233898180
#>  [826] -0.0113818327 -0.0013536232  0.0677788265 -0.1863226100 -0.0757074751
#>  [831] -0.0741467243 -0.0212812738  0.0873214726 -0.0389656548  0.0344962057
#>  [836]  0.0104395063  0.0526363837  0.0836124546  0.2480303870  0.0828326953
#>  [841] -0.3881736312  0.1814606235 -0.2801169668 -0.1899366766 -0.0867643889
#>  [846]  0.0181573527 -0.0069172195 -0.0531932836  0.1167506966  0.1409249712
#>  [851]  0.1082099785 -0.0734571836  0.2513054520 -0.1907907146  0.1100074983
#>  [856]  0.0549625247  0.0724074999  0.0914053423  0.0954510343  0.0175481285
#>  [861] -0.0130982355  0.0115027095  0.0541487747  0.1199589717 -0.2120019691
#>  [866]  0.1059181374  0.0134309837 -0.0895020105  0.2053792407 -0.0435980075
#>  [871] -0.0395198353 -0.1822187236  0.2620551389  0.0773507348 -0.1308340870
#>  [876] -0.0083632119 -0.1052539731  0.1079386083 -0.0771898060 -0.0354109243
#>  [881]  0.1744388486  0.1432084269  0.1091167066  0.1611314722  0.0561793928
#>  [886]  0.1062371917 -0.0681863203 -0.1180133798  0.0682197614 -0.0345100814
#>  [891] -0.0353145280 -0.0201361664  0.0766340947  0.1205760095  0.1912648245
#>  [896] -0.0174137957 -0.0159896285 -0.0115670635 -0.0991779061 -0.0128948284
#>  [901]  0.1847908357  0.1220011861 -0.0881056226 -0.2322972439  0.1267730529
#>  [906] -0.1763629869 -0.1100482613 -0.0137264543 -0.0377579648  0.2581891581
#>  [911]  0.1408194955 -0.1357377909 -0.0968303557 -0.2820587635 -0.1277777586
#>  [916]  0.0230823884  0.1440763324  0.0792017537 -0.0875681373 -0.0552774951
#>  [921]  0.0810785425  0.1154971199  0.0997670704 -0.0426865482  0.1274333328
#>  [926]  0.1192802161 -0.0802201341 -0.0370093835 -0.0450500649  0.3630426808
#>  [931]  0.1378994282 -0.0178700300  0.1705071491 -0.2780832629  0.4024083891
#>  [936]  0.0328635460  0.2071333886  0.1638274810  0.2004794830  0.2078295871
#>  [941]  0.1627875207  0.2239806640  0.1353292361  0.0355157933 -0.1058929615
#>  [946]  0.0946516075 -0.0726969730  0.0567407855  0.0683016569 -0.0838708798
#>  [951] -0.0312329067  0.0282357110 -0.0240909035  0.0180706370 -0.0607421042
#>  [956]  0.2134099753  0.0156476124  0.1830068221  0.1281602212  0.0504356727
#>  [961] -0.0425959874  0.0711575164  0.2942657086  0.0523779643  0.1669865388
#>  [966]  0.0921108526 -0.0603903324  0.0580709984  0.0257239135 -0.0256850782
#>  [971]  0.1412506096  0.2249088846  0.2033889683 -0.0220870052  0.0576648584
#>  [976] -0.2879649058 -0.1249831532 -0.1019380060  0.1735573103  0.1526606979
#>  [981] -0.0275323843  0.3272045854  0.1460618163 -0.1734347224  0.0790388746
#>  [986]  0.2239956282 -0.0012080221  0.1702430604  0.1133363106  0.0144183888
#>  [991]  0.2402177686  0.1239561567  0.2152467589  0.1883580933  0.2106866917
#>  [996]  0.2455044991  0.1055402412  0.1479774720 -0.0990524951  0.1564227497
#> [1001] -0.1296154386  0.0615602839  0.1907148070  0.1683918039 -0.0049965194
#> [1006]  0.0365766688  0.2108438147  0.1118244940  0.0736333226  0.1863005758
#> [1011] -0.1454527134  0.1592229785  0.0690372119 -0.0268828087  0.1557028669
#> [1016]  0.1637945010  0.0634367545  0.1422679841  0.0727069912  0.0141364351
#> [1021]  0.1138419814 -0.0024836712  0.1850046214  0.4276320058  0.0892435199
#> [1026]  0.0513569787 -0.0574734911 -0.1511745314  0.2391415283  0.0117895198
#> [1031] -0.0109590163  0.1755173621  0.1355237805  0.0287944524  0.2146959898
#> [1036] -0.1062716271 -0.0603062074  0.1259329924 -0.0726380164 -0.0563396120
#> [1041] -0.0658121416 -0.0241045939  0.0657307434 -0.1126272118  0.2107664473
#> [1046] -0.1869696086 -0.1571183350  0.0593396088 -0.0429103034  0.0770447457
#> [1051]  0.0040257784  0.0948508779  0.1284024834  0.1603348469  0.0057407098
#> [1056] -0.1205222210 -0.0621008712 -0.0053462405 -0.1748401847 -0.1108085177
#> [1061]  0.0695387495  0.0192167530  0.1165184966  0.0774075576  0.0506847822
#> [1066] -0.1406771063  0.0214098934  0.3812775626 -0.0475937843  0.0336174512
#> [1071]  0.2710229619  0.1199312042 -0.1854888099 -0.0774223069  0.2597549069
#> [1076]  0.0349245695  0.2736852966  0.1544388159  0.0547190787  0.3415663544
#> [1081]  0.3448394161  0.1356088052  0.1188456748 -0.0992940887  0.1446862164
#> [1086] -0.1152553900 -0.0314068101 -0.0588958685 -0.0627929973 -0.1651618343
#> [1091] -0.0262312210 -0.0646143880 -0.1276698523  0.0362690566  0.1409874963
#> [1096] -0.0471332429 -0.0603698901  0.2101949579  0.0191159428  0.2890675020
#> [1101]  0.0819876140  0.0247986659 -0.0202345223 -0.2152335768  0.0057682349
#> [1106] -0.0247444964 -0.0056301369  0.1843353987  0.0553004013  0.0585922039
#> [1111]  0.1725803541  0.0371991864 -0.0314131336  0.2404293539 -0.2928868073
#> [1116]  0.1115098235 -0.0249954390  0.1943650739  0.0117705462  0.0402901906
#> [1121]  0.0891649103 -0.0491547405  0.4237957464 -0.2060433790  0.0498566913
#> [1126]  0.1372545756  0.0366925866  0.0926041474  0.0303921723  0.0345444547
#> [1131] -0.1648168918 -0.1299757515 -0.1338993009  0.1152878767  0.0471814465
#> [1136]  0.0105694111  0.0619384803 -0.0183655870 -0.2132888082  0.2040785503
#> [1141]  0.0608985614 -0.0070530360  0.0220028634 -0.0118050109  0.0685585195
#> [1146] -0.0007780782  0.1081936635  0.1695720406  0.0091772241  0.0656011453
#> [1151]  0.1348109738 -0.1103698353  0.0986342671  0.0502093273  0.0067243400
#> [1156]  0.1041796357 -0.0422883844  0.1166481789 -0.1978604730  0.1959694290
#> [1161]  0.1550393087  0.0010731137  0.1240384642  0.1397228238  0.1839021270
#> [1166] -0.0842121585  0.2050399416  0.2894452775  0.0060847064  0.0756684880
#> [1171]  0.3662555844  0.2963597882  0.3053413176  0.0991478693  0.4191651085
#> [1176]  0.0895374667  0.0983628226 -0.1861878199 -0.0391934669 -0.2949952092
#> [1181]  0.2214177607  0.0399410732  0.1775526838  0.1925316985 -0.1755864995
#> [1186]  0.0462662483 -0.0045141051  0.2428683068  0.1727511059 -0.1601699601
#> [1191] -0.0478278226  0.1410980688  0.2045874435  0.1285019892 -0.0070115197
#> [1196]  0.2095621069 -0.0464663205 -0.1113872839 -0.0883185551  0.1668791071
#> [1201] -0.0892704304  0.0884830135  0.0401900247  0.0905724237  0.0557918760
#> [1206] -0.0661047830  0.0844701808  0.0730617904  0.0250867983  0.1230592323
#> [1211] -0.0607745980  0.0398340078  0.1070320559  0.0047604077  0.0479959262
#> [1216]  0.2052852310  0.0698453043  0.0876609006  0.2205024955  0.0793667838
#> [1221]  0.1232785956  0.0242201521 -0.0270163772  0.0027449515  0.0627383418
#> [1226] -0.0102734850 -0.0632183226  0.0871878834  0.0478147298  0.0060236917
#> [1231]  0.2498812544  0.0707545209  0.2553376959  0.0625444622 -0.1234563535
#> [1236]  0.1324651600  0.1057252965  0.0734120707  0.2833090611  0.0393092906
#> [1241] -0.0150762358  0.2212753419 -0.1008222183  0.1721772745  0.0827722877
#> [1246] -0.0928911978 -0.0216936415  0.1133081636  0.1424733543  0.2133520541
#> [1251]  0.0201510783 -0.0607075842  0.0776772420 -0.0588346398 -0.0679177194
#> [1256]  0.0352777283  0.0961763539  0.0742683490  0.0913788682 -0.0516797350
#> [1261]  0.0782413990 -0.1247206531 -0.0366927386  0.0025438812 -0.1244232614
#> [1266]  0.2023651983 -0.1293956458  0.1974861559 -0.0524233215 -0.0255707082
#> [1271]  0.1214579010 -0.0206644541  0.0232681385 -0.0680751659  0.0797131737
#> [1276] -0.1312332541  0.1229062267  0.0989798886 -0.1756982567  0.0486075809
#> [1281]  0.0750681648 -0.0107538187  0.1043841336  0.3346352820 -0.0676951741
#> [1286] -0.0270890999 -0.0512667696  0.1836507092  0.1377118521 -0.1837142993
#> [1291]  0.0260358966  0.0334716151 -0.0634335188 -0.0578867390 -0.1442491584
#> [1296] -0.0319003785 -0.0440677931  0.0427982717 -0.1531752364  0.1185615212
#> [1301]  0.1233685541 -0.1471607282  0.0079090168 -0.0425672341 -0.0521636373
#> [1306] -0.0340843847 -0.0599700477 -0.0076870402  0.1219570033  0.1126097022
#> [1311]  0.0526393208 -0.0158958882 -0.0745565232  0.0407226895  0.2860322875
#> [1316]  0.0023945302  0.0188129094  0.0880802494  0.1481204100 -0.0282298851
#> [1321]  0.1890026659  0.0868743884 -0.1408206196  0.0422810427  0.0518148246
#> [1326]  0.0995995557 -0.3045145554 -0.0001699823  0.0726849559 -0.1206358722
#> [1331] -0.1941036202 -0.2477399160  0.2053143513  0.0387402025  0.0400928394
#> [1336]  0.1628138353 -0.1021377958  0.0947418037  0.0493019047 -0.1361689125
#> [1341]  0.0269732513  0.0663211846  0.1173765310 -0.0212098380  0.3276848437
#> [1346]  0.2495537880  0.0373284501  0.0036204499  0.1124043306 -0.1902677856
#> [1351] -0.1045735244  0.0693242589  0.1196396880  0.0155924801  0.0729748398
#> [1356]  0.2010393719 -0.2369993731  0.1999092963  0.1440217431 -0.0179180005
#> [1361] -0.1128560370  0.2565521536  0.0650813904  0.2314521113 -0.0220263974
#> [1366] -0.0014007142  0.0362235022  0.0421600003  0.0178215056  0.2465376409
#> [1371]  0.2314544994 -0.1622556605  0.2336207090  0.1334036171  0.3307559184
#> [1376] -0.0004961289  0.0591641475 -0.0023509848 -0.1569092650  0.0551421066
#> [1381]  0.2364874903 -0.1290054598 -0.1223331127 -0.2103761235  0.1130186762
#> [1386]  0.1893493923  0.0569367457 -0.2091439211  0.1871209446 -0.0805798246
#> [1391]  0.0538483470  0.1540822539  0.1314539533  0.1092699476  0.3352383294
#> [1396]  0.0008326317  0.0759980507  0.0210002889 -0.0149845389  0.2651192724
#> [1401]  0.2641816104 -0.1026939075  0.1270268262  0.1214644330  0.2807020493
#> [1406]  0.0901122170 -0.1272614339 -0.0630600842 -0.0976472062  0.1605633834
#> [1411] -0.1089330320  0.1476020447  0.1128517277 -0.0239721835  0.0204922744
#> [1416] -0.2244863952  0.3657167714 -0.0084527894  0.1466958041 -0.0694060767
#> [1421] -0.0098342555  0.2398260770  0.2606755333  0.1730295561  0.1410591417
#> [1426] -0.0414434649  0.1079475747  0.1270291793  0.1197999376  0.2637542730
#> [1431]  0.0792082195  0.0137445063  0.0404935760  0.0904064265 -0.0575428510
#> [1436] -0.0240997877  0.1278566241 -0.0668114410  0.1077550022 -0.0196260413
#> [1441] -0.1097236867  0.2461437925  0.0132868127  0.2041252843  0.1623127299
#> [1446]  0.0138974865 -0.0120746246 -0.0120810578  0.0724178096 -0.0207465765
#> [1451]  0.1243116643  0.1490862536 -0.0810551479  0.0498997590 -0.1109906784
#> [1456]  0.1813457097 -0.0688827891  0.0873857795 -0.0297650185 -0.0698141771
#> [1461]  0.0203062155  0.0615527833  0.0004680058  0.0820823478  0.0163599183
#> [1466] -0.1553359676 -0.2193374633 -0.0071983991 -0.0298975629  0.1352221692
#> [1471]  0.0112988605 -0.1325875539  0.1092009845  0.0071201616 -0.2209199002
#> [1476] -0.0178852802 -0.0361629547  0.0332141706  0.1743737656  0.1142131689
#> [1481] -0.0076271589 -0.0114199441  0.1548488898  0.0457214544  0.0950583961
#> [1486]  0.1906984729  0.3238140386  0.0578013606  0.1356319281  0.1606121221
#> [1491] -0.0506090428  0.1226330257  0.1123829601 -0.0488238101  0.1624468243
#> [1496] -0.0514115569 -0.0787414290  0.1009729871  0.1581471967  0.0129286314
#> [1501] -0.1468219125  0.0635532092  0.0168970321  0.0908204916  0.1210247637
#> [1506]  0.2135051950  0.0984434917  0.0501297422 -0.0514979429 -0.0035170835
#> [1511]  0.1264421593  0.0173629996 -0.1350280570  0.0756204680 -0.1335826885
#> [1516]  0.1396496158  0.1295245787  0.2113857200  0.0599041144  0.1410731396
#> [1521]  0.1587638299 -0.0067806799 -0.1137195154  0.1364339562  0.3230372507
#> [1526]  0.0128334588  0.1836257348  0.1029338031  0.0964291655  0.2854019986
#> [1531]  0.2934121886  0.3182735430  0.2166290086  0.1924214121  0.2180173129
#> [1536]  0.0278164958  0.2629379618  0.0143335154  0.0563294305  0.1215321768
#> [1541] -0.0365548204  0.1676636600 -0.0070807449 -0.0702197623  0.0069329852
#> [1546] -0.0758124736  0.0350922583  0.0619062417  0.0981985366  0.1764976819
#> [1551]  0.0266349387  0.0029879915  0.0314620770 -0.0446570554  0.2111999658
#> [1556]  0.1779053695 -0.0547989380  0.3359617212 -0.1170794383 -0.0526576521
#> [1561]  0.4179254159  0.1865735090  0.0929781623 -0.0779026562  0.0992111347
#> [1566] -0.0776862576  0.2593505158  0.0495730263  0.2749055417 -0.1552751441
#> [1571] -0.0279944254  0.1913143657  0.0423481872  0.0219996462 -0.0329264463
#> [1576]  0.2816558835 -0.0421631064 -0.0888775122  0.1504118513  0.0777330516
#> [1581]  0.0959453541 -0.0857292517  0.0481501417  0.2748704759 -0.0656630676
#> [1586]  0.0837729435  0.1475847323 -0.1421516779  0.1847518959  0.0221127862
#> [1591] -0.3190449698  0.0434755578 -0.0477523680  0.1563155916 -0.1703480528
#> [1596]  0.1452044044  0.0489163582  0.1583875132  0.0663410017  0.0879452047
#> [1601] -0.0898095325  0.0011019128  0.2382828792  0.1794779873 -0.0066891847
#> [1606] -0.0176914872  0.1530799459  0.1111787949  0.0460470289  0.0646592200
#> [1611]  0.0129125108  0.1780930870 -0.0180161999  0.0024465795  0.0729437350
#> [1616]  0.2652816221  0.3528984430 -0.0056729150  0.2014055488 -0.0163182604
#> [1621]  0.1229503654  0.0894129865  0.1176647691  0.1374936472  0.1695486027
#> [1626] -0.0857818723  0.2357491545  0.0111989163  0.0557330300 -0.0667363302
#> [1631]  0.1473062402  0.0365984653  0.0690799558 -0.0936462619  0.0340248607
#> [1636]  0.1718233500  0.2094470862  0.0434914156 -0.0360546439 -0.0167984250
#> [1641]  0.0703120305  0.0209643015  0.1249757357  0.1596059192  0.1968596416
#> [1646]  0.2011153752 -0.1549216428  0.0401221229  0.1983524358  0.0611457876
#> [1651]  0.3661120145  0.3428791302  0.1642389266  0.0470784344  0.1611383746
#> [1656] -0.1702106823 -0.0073664693 -0.0801135422 -0.0393088860  0.2119139005
#> [1661]  0.3413564355  0.0146651330 -0.0196629736  0.0924147157  0.0070783224
#> [1666]  0.1523855031  0.1322667443  0.0978699434  0.0656787384  0.2124593507
#> [1671]  0.0332530276  0.1046995456  0.0046624470  0.0629994289 -0.1067691214
#> [1676]  0.1164012763  0.2710602954  0.0635853429  0.1721794743 -0.0326499967
#> [1681]  0.2684807233  0.1287393729 -0.0759538102  0.0152800421 -0.0575268272
#> [1686]  0.3020685987 -0.1391276630  0.0827121475  0.0216392877  0.2581697136
#> [1691] -0.0497209800  0.1486498192  0.0799699350  0.1878501636 -0.0562672375
#> [1696] -0.1152054346 -0.0755429473  0.2332192754 -0.0949490958  0.0724748890
#> [1701]  0.0458811616  0.1404808015  0.1360080936  0.1476916350  0.1091715271
#> [1706]  0.0173334659  0.1386399571  0.1201456962  0.0949405794 -0.0494491394
#> [1711]  0.1504489209 -0.0139565188  0.1066810820  0.0304286809 -0.0720145213
#> [1716] -0.1177400582 -0.1058967511  0.0996300581  0.1537489529 -0.0364232951
#> [1721]  0.1009699101 -0.0667696464 -0.1242581766  0.3180196516 -0.1078572162
#> [1726]  0.2682818868  0.1532000355  0.1842389375 -0.0148528874 -0.0438986330
#> [1731]  0.0821565317 -0.1451310303  0.1749785883 -0.1300008596  0.2364172460
#> [1736]  0.1497152771 -0.1352434217  0.1569936139  0.1689696396  0.2510192239
#> [1741]  0.0359435955 -0.0379449390 -0.1154776628  0.0461769373 -0.2898213238
#> [1746]  0.0542887401  0.0499541932  0.3298131079 -0.1220245880  0.1264814913
#> [1751]  0.0394063916 -0.0743713382  0.1709188150  0.0621735313 -0.0100911808
#> [1756]  0.0497100811  0.0266329325  0.0297493664  0.2375810789  0.1064510592
#> [1761]  0.1986674434  0.0856679012 -0.0077923452  0.0641658334 -0.0081588159
#> [1766]  0.0456229568  0.1201350592  0.0529572528  0.2242310129  0.2793036176
#> [1771] -0.0223330123  0.1631907373 -0.1910169338 -0.1674295256  0.0268820574
#> [1776]  0.0973622033  0.0624145703  0.1212373054  0.0025330691  0.0269169055
#> [1781] -0.0156727778  0.3358725722 -0.0727783589  0.0949027737 -0.0511620650
#> [1786] -0.0230575397  0.1241450941  0.0513875713  0.1395813653  0.1344362154
#> [1791]  0.1678085929 -0.0462656033  0.2847791395  0.2524207648  0.2078473885
#> [1796]  0.0595918545  0.0143853438  0.1427645669 -0.1598270171 -0.0745701941
#> [1801] -0.0605979495  0.1678650181  0.1082116824  0.3177449363  0.1884945140
#> [1806] -0.1215435492  0.1913021855 -0.1202947236  0.1359564602  0.0516131088
#> [1811]  0.1822400660 -0.3214557840  0.4440784732  0.2041328764  0.0704069161
#> [1816]  0.1671864002  0.0617878549 -0.1073300924  0.3430714036  0.2215015677
#> [1821]  0.1889609376 -0.0683841521 -0.0660486788  0.1978517424 -0.0487151755
#> [1826] -0.0184819922  0.1633626466  0.1301846670  0.1670052593  0.0485178073
#> [1831] -0.0505738872  0.3019660852  0.1109275670  0.2196527354 -0.1126966823
#> [1836] -0.1301676541  0.0114290596  0.0443698912 -0.2011856182  0.1455224995
#> [1841] -0.1139083297  0.0130970547  0.0070431753  0.1027788770 -0.0119845983
#> [1846]  0.1447586079  0.3736258594 -0.0064332823  0.2301741117  0.1190476381
#> [1851] -0.1464428079  0.3155599416 -0.1319235041  0.1043750870  0.2938867226
#> [1856]  0.0912458179  0.2528854492  0.0352487608  0.2668778726  0.0175403079
#> [1861] -0.0676146942  0.1059769180  0.2295126490 -0.0006209341 -0.1088511259
#> [1866]  0.0654625188  0.0847613400  0.1106726530  0.0285816230  0.0501600866
#> [1871]  0.0195623981  0.2032497007  0.0966059746 -0.0187566194  0.2788240619
#> [1876]  0.3243834941  0.2133204012  0.2885469560  0.0615789047 -0.0256244133
#> [1881]  0.1720293786 -0.0764217761  0.1426725036  0.1493997929  0.0969458510
#> [1886]  0.2848414399 -0.0963760751  0.0918143957  0.1021901126  0.2567814759
#> [1891]  0.0545283156 -0.0246710992  0.0117690230  0.1303939476  0.0933694931
#> [1896] -0.2762379643  0.0478068616 -0.0111701236 -0.1165629535  0.0843179612
#> [1901] -0.1807335216  0.1022050437  0.0912591344  0.0995271990  0.1642467357
#> [1906] -0.1204713401 -0.0003690057  0.0846609665  0.0222307141  0.2543127619
#> [1911] -0.0429771481 -0.0984653017  0.0729696663  0.2344651998 -0.0881183283
#> [1916]  0.0785152312  0.2025040690  0.0640708623  0.0600257442 -0.0804402574
#> [1921]  0.0229023855 -0.3176811881  0.2654436234  0.0998939267  0.1476634096
#> [1926]  0.2398654913  0.0738063561  0.0334470258 -0.0553035153  0.0679777569
#> [1931]  0.0619990782  0.0578698998 -0.0536372985  0.0270927949  0.2421574098
#> [1936]  0.0651401706  0.1350220851  0.2737385201  0.0566902407  0.4363837076
#> [1941] -0.0812460812  0.0913777619  0.0499380075  0.1905363816  0.4064992481
#> [1946] -0.0195491474 -0.0464497453  0.0338498578  0.0813577203 -0.0651632994
#> [1951] -0.0986101063  0.2874048056  0.0423451069 -0.0727136436 -0.1633448312
#> [1956] -0.1368884254  0.0006129589  0.0140111087  0.0289392041  0.0996738595
#> [1961]  0.0214398508  0.0194534651  0.0130749136  0.1919371116 -0.0420101415
#> [1966]  0.5147725107 -0.1031809431 -0.0024374058 -0.0545324868 -0.0134398971
#> [1971]  0.0379327445 -0.0271397922 -0.0843368694 -0.2343878749 -0.3311619274
#> [1976] -0.0513187183  0.0650546650  0.2079243031 -0.1221936100 -0.1696653741
#> [1981]  0.1954322497 -0.0310434193  0.0220601769  0.0810614831 -0.1634701413
#> [1986]  0.0829441809  0.0639063091 -0.3095699574  0.2079169599  0.1079135147
#> [1991] -0.1477175644 -0.0722437579 -0.0290509056 -0.0352484535  0.1347604261
#> [1996]  0.0348952471 -0.1683441976  0.0730490488  0.1249064036  0.3067335054
#> [2001] -0.0099070180  0.1498373295  0.0531865285  0.1743197528  0.1794966542
#> [2006]  0.1823300305 -0.0378585712  0.4104004209  0.0260322309 -0.0778652560
#> [2011]  0.0764138965  0.2227283579  0.1490319931 -0.0477981325  0.1333494975
#> [2016]  0.0124414398  0.1419778440  0.0074356035 -0.1207210928 -0.0583521593
#> [2021]  0.1272064994  0.0966405986  0.0236027392  0.3873142087 -0.0386940349
#> [2026]  0.0033137472 -0.0534449764 -0.0734836118  0.3505056527  0.3094289092
#> [2031]  0.0010828928  0.0714765681  0.1533598936  0.0833392315  0.0413494653
#> [2036] -0.2581229540  0.1681124723 -0.1900620136 -0.0635127295  0.1874364924
#> [2041]  0.1699826112  0.0873662767  0.0025444753  0.2809179409  0.0835302615
#> [2046] -0.0627328279 -0.0014263152  0.1819912416  0.2320054505  0.1401126736
#> [2051]  0.1455073243  0.1604466613  0.0101344596 -0.0648870554  0.2517255877
#> [2056] -0.0564615605  0.0994767309 -0.1573672622  0.1249257073  0.1886128735
#> [2061] -0.0237123606  0.1127009993  0.1674706486  0.0335158401  0.1140375841
#> [2066]  0.1602549034  0.2116670463  0.1320818415 -0.0811483175  0.0255119709
#> [2071] -0.0831038629  0.0343141811  0.2097865504 -0.1041007162  0.1514982267
#> [2076] -0.0221136844  0.1270617887 -0.1467308381  0.0704048013  0.0891638648
#> [2081] -0.0450295770  0.0416680964 -0.0117689775  0.0765987719 -0.0220943389
#> [2086]  0.0033168210  0.0386339928  0.2306001515  0.0375355216  0.0390701391
#> [2091]  0.0349543817 -0.0602520104  0.0878493648 -0.0629449651  0.1299454379
#> [2096]  0.1094701667 -0.0335776836  0.0463010890  0.1519345002  0.2832535470
#> [2101] -0.0229921614  0.0215798154 -0.0280273909  0.1910683945 -0.0124251162
#> [2106]  0.0614185857  0.0077491095  0.0567449085  0.1408030819  0.2385137677
#> [2111]  0.0071753815  0.1397091768  0.0334402189  0.0555061979  0.2218450364
#> [2116]  0.0285835155 -0.0420390607 -0.0383553558 -0.1066209719  0.2953473511
#> [2121]  0.3095478218  0.1216329335  0.1272056154  0.0832099290  0.0982192317
#> [2126]  0.0019759505  0.1551137563 -0.0093209904  0.1014106786  0.1006216914
#> [2131]  0.1411504656  0.1738639017  0.0331058378  0.1333212191 -0.0538615595
#> [2136]  0.0841414248  0.0453139758  0.1815090043  0.2867372404  0.1234033298
#> [2141]  0.1404389280 -0.0122212844  0.0978080853  0.0486864970  0.0511549855
#> [2146] -0.1136409192  0.0413792180 -0.1357809981 -0.0758449975  0.0583022759
#> [2151]  0.2034288531  0.1764905528 -0.0054855215  0.2436298738  0.0257820994
#> [2156]  0.3065894392  0.1879572100  0.0644060187 -0.0592299778  0.0473054396
#> [2161] -0.1623199999 -0.2185836267  0.1868074484  0.0471429199 -0.0212772023
#> [2166] -0.0041672901 -0.1488189178  0.2453639913 -0.0002407227 -0.0313278388
#> [2171]  0.0162537374 -0.0552363977 -0.0793832344 -0.0365350580 -0.2458430686
#> [2176] -0.1143490362  0.1386813265  0.0199793771 -0.1053886815  0.0416496762
#> [2181]  0.2680212180  0.0336700835  0.0806598967 -0.1103188294  0.0837821822
#> [2186] -0.0144854696 -0.1101944929  0.1236591066  0.2176815336  0.1435048211
#> [2191] -0.1862819921  0.0320827791 -0.1905805615  0.2110595157  0.1533837299
#> [2196]  0.0786880048  0.0907930611  0.1534826153 -0.0525411632 -0.0429131908
#> [2201]  0.0861850006  0.0061923064  0.1197231909  0.1682677031  0.1597852612
#> [2206]  0.3119775167 -0.1620305015  0.1940166521  0.0773343336 -0.1264966635
#> [2211]  0.1146606645  0.1657176863 -0.0503961970  0.0337722841 -0.1580671488
#> [2216]  0.1258385414  0.2204625626  0.2440513264  0.1259605762  0.0245171635
#> [2221]  0.1556069639  0.0139790210  0.1052467652  0.1066054091 -0.0071613904
#> [2226]  0.0825929672  0.0427641806  0.0750368903  0.2141780396 -0.0332727031
#> [2231]  0.4950807389 -0.0194939075  0.0652567627  0.1982009494  0.0171989486
#> [2236]  0.0820258732  0.1061604849 -0.0629230666  0.1274928474 -0.0130759204
#> [2241] -0.1376374047  0.2028521420 -0.0808442203  0.0011422491  0.0810016137
#> [2246]  0.0228295035  0.1514188498  0.0943875692 -0.1097443679 -0.0077977912
#> [2251]  0.0356154647  0.2580461686  0.0542059589  0.0268018376 -0.1443421056
#> [2256]  0.1903551524  0.2000230464 -0.0901636681 -0.0653948612  0.0464885113
#> [2261]  0.0951050792  0.1144959509  0.0373225374  0.0288425082 -0.0063999777
#> [2266]  0.1250071812  0.0254315501 -0.0751348154  0.1209249870  0.1862250443
#> [2271]  0.1744154236  0.1114864800 -0.0881766344 -0.1068662746  0.0919437315
#> [2276] -0.1210840550  0.1407972104  0.1480628497  0.0885212844  0.1362749900
#> [2281] -0.1512167599  0.1536515593  0.0439452198  0.1666703539 -0.0416116931
#> [2286]  0.1120474267  0.0903865877  0.1502306508 -0.2019725880 -0.0182779366
#> [2291]  0.0512662835  0.1743112978 -0.0066887229  0.1187275731 -0.0581162148
#> [2296] -0.0623973365  0.1965855732  0.0422250102  0.3383914728  0.0588159920
#> [2301]  0.1130342682  0.1046061085  0.0134651851 -0.0021858653 -0.0032732284
#> [2306]  0.0871095461  0.0324610605  0.2017077704 -0.0989290532  0.2177005092
#> [2311]  0.1695553054 -0.0627743934  0.1213789666  0.0055686106 -0.0792448537
#> [2316]  0.0451447023  0.0150675415  0.3506730512 -0.1142853143 -0.1545547088
#> [2321]  0.1286148469  0.1837617250 -0.1721664303 -0.0635925213  0.0420916574
#> [2326]  0.1788339654  0.1956851841 -0.0032965233  0.2629226339  0.1447189620
#> [2331]  0.1637262037 -0.1801487416 -0.0249705195  0.2804421969  0.0812548631
#> [2336]  0.1217369731  0.1713315783  0.2149113190  0.0035690015 -0.0080119793
#> [2341]  0.0577761356 -0.1174119169  0.1196541878  0.1044648656  0.1093445309
#> [2346] -0.0425931561  0.0114780492 -0.0331799566  0.2808094276 -0.2060735909
#> [2351] -0.2192074115 -0.1860082826  0.1194700919  0.1146547973 -0.0245261439
#> [2356] -0.0435361344 -0.0486900664  0.0180968640  0.2661170334  0.0753008937
#> [2361]  0.1266862183  0.0773229417  0.2532655098  0.0098151638 -0.0042464499
#> [2366]  0.0801237570  0.1803684228  0.0621742172  0.3066753084 -0.0997405169
#> [2371]  0.0084208480  0.2244519760 -0.0128518022  0.1697182843 -0.0209948648
#> [2376] -0.1113381567  0.0255460054 -0.0644642552  0.1128505979 -0.1900231362
#> [2381]  0.0053466039  0.0234960835  0.0694943502 -0.0869617509  0.0149457286
#> [2386]  0.2399086008  0.1196582494  0.1176306554  0.2137415970  0.0535679121
#> [2391]  0.1788252083  0.0375820300  0.1554236702  0.1749554990  0.1752197251
#> [2396]  0.2431115998  0.0070096344  0.2166295796  0.2434933147  0.2719811477
#> [2401]  0.2469797517  0.2660953040  0.0317064440  0.2005419677  0.2266772379
#> [2406]  0.1592394398 -0.0273872250 -0.0294105402  0.0954191204  0.1536902016
#> [2411]  0.0521761303 -0.0286568121  0.0105778373  0.1754853475  0.1085521431
#> [2416] -0.0635273659 -0.0087538293  0.1371347357 -0.0010665044 -0.0449988164
#> [2421]  0.0228436930  0.0993459640  0.2013356724  0.1293731712 -0.0659671843
#> [2426]  0.0244689609  0.3488813701  0.1609455041  0.1583175924 -0.1710408207
#> [2431]  0.1358425767  0.1508474363 -0.0652421144  0.2269717114  0.0019338834
#> [2436]  0.1766292739 -0.1605165030  0.2797068346  0.1043115628  0.1203548322
#> [2441]  0.1731274958 -0.1933090179  0.2783602542  0.2286634482 -0.0875555624
#> [2446]  0.0823151399 -0.0386873843  0.2990192686 -0.2011757662  0.2778531749
#> [2451]  0.1631870580 -0.1104688108 -0.0688536177 -0.0780494445  0.1548480780
#> [2456]  0.0685263316  0.0591236249  0.1538279919 -0.2396689728  0.1392351314
#> [2461]  0.0471400250  0.1815727247  0.0208951330  0.3293766306  0.0881959224
#> [2466]  0.1148908460  0.0922173111  0.1722039128 -0.0664639208  0.1154583741
#> [2471]  0.2438170658  0.0170633950  0.4103746725  0.4410725162  0.1048857325
#> [2476]  0.1380986640  0.1004345768  0.1520955911  0.1665141098 -0.0419168500
#> [2481]  0.1062285056  0.1877645586  0.1303216238  0.2658409841  0.2412738703
#> [2486] -0.1011831594  0.3378753549 -0.0806191796  0.0440543057 -0.1812458612
#> [2491]  0.0357571433 -0.0446551915 -0.0520296346  0.1279322039  0.1391152724
#> [2496]  0.0150781571  0.1598944252  0.0094971134  0.0925353501  0.2137596412
#> [2501] -0.0016711290 -0.2050917639 -0.0657725884  0.0074833858 -0.0435906176
#> [2506]  0.2350315774  0.1602118674  0.1273298053  0.1664677330  0.2443362732
#> [2511]  0.0489575922  0.1191762590 -0.1990129347  0.0051243813 -0.0957960290
#> [2516] -0.2173680343  0.1680232599  0.2512953590 -0.0721020497 -0.1674393566
#> [2521]  0.0803421816 -0.1049262981  0.3954593172 -0.0986661098  0.1178248681
#> [2526]  0.0012174908 -0.0502323852  0.0532788874  0.1295155984 -0.0716886774
#> [2531] -0.1048438232 -0.0161568035 -0.0865185600  0.0948989183 -0.1186983115
#> [2536]  0.1524890199  0.0517595555  0.0693778790 -0.1123483384  0.1624006718
#> [2541]  0.1356146163  0.1674744586 -0.0039216286  0.2103755963 -0.0539378331
#> [2546]  0.1665516323  0.2148963865 -0.0086102448  0.1665896181 -0.0175352403
#> [2551]  0.0922146920  0.1683448737 -0.0917049290  0.2601930449  0.1273895236
#> [2556]  0.0896641161 -0.1063298133  0.0385076551 -0.0196120335 -0.0742032331
#> [2561] -0.1937751884 -0.1025482288  0.1591802387  0.3942587397 -0.0511810115
#> [2566] -0.1570795131 -0.0325547187  0.1709803686 -0.1037454682  0.1682163677
#> [2571]  0.2175481019 -0.0705771779 -0.0312384447 -0.1963404739  0.1533831606
#> [2576]  0.0874950832 -0.0938456570  0.0030329352  0.0505106139 -0.1033365625
#> [2581]  0.0155831019  0.0538855326  0.1058353200 -0.0259383432  0.1917480180
#> [2586] -0.0129782370 -0.1447534129  0.1415949660  0.1244863510  0.1281405527
#> [2591]  0.0342542687  0.0137174363  0.0539790133  0.0703129290  0.2347785411
#> [2596] -0.0137468308  0.1339843825  0.1561168196  0.1805557126  0.0479185878
#> [2601]  0.1068505122 -0.0018580399  0.1332158065 -0.0627443086  0.2661141084
#> [2606]  0.0983917620  0.1711893931 -0.0251108864  0.2780058388  0.0959577358
#> [2611]  0.1181186464  0.2267818588 -0.0826545009  0.2351433037  0.1563572505
#> [2616]  0.0874557838  0.2067252788  0.2849003475  0.0642825358  0.0434160632
#> [2621] -0.0641134392  0.0907993916  0.0849943535 -0.1445481280 -0.0534013075
#> [2626] -0.0117005491  0.2461321442  0.1244938148  0.0379438729  0.0583037857
#> [2631]  0.1663430705  0.0708508640  0.1010820997  0.0643473774  0.0756702012
#> [2636]  0.0143476970  0.0337656390 -0.0407697456 -0.0087701919  0.1506684455
#> [2641]  0.2220314861  0.1794603803  0.0229026698  0.0358102755  0.2194919670
#> [2646]  0.2064274095  0.0499153401  0.0969029242 -0.1296544945  0.0328957010
#> [2651]  0.2991209926  0.0563812346  0.0817670688  0.2783643936  0.0857121620
#> [2656]  0.0611001648  0.0319165519 -0.3141489367 -0.1206286237  0.0779083605
#> [2661]  0.1082890180  0.0652396396 -0.0471259723  0.0743918841 -0.0420803114
#> [2666] -0.0713592145  0.0573644210 -0.1130055149  0.0950991536  0.1336703820
#> [2671] -0.0521302219 -0.0846759802  0.2035642021  0.1752062414  0.1457774330
#> [2676] -0.0048112988  0.2814658374 -0.0489810300  0.0432543637  0.0957142022
#> [2681]  0.0933240738  0.1474750919  0.2076512689  0.1256486889  0.0767640731
#> [2686]  0.0268867546  0.2618244090  0.0448462882  0.0145040577  0.0824814382
#> [2691]  0.2504493558  0.3105865391  0.0025621426 -0.0501774371  0.0035204832
#> [2696]  0.1150580246 -0.1757237701  0.1630467340  0.2351435087  0.0984991993
#> [2701]  0.0875866404  0.0272131036  0.2258886733 -0.0600570449  0.1054356776
#> [2706]  0.1395753443  0.1311607408  0.0385114929  0.0793438944 -0.0099622134
#> [2711]  0.0023302803  0.0292948508 -0.0352793395  0.1398709827  0.0743662998
#> [2716] -0.0537959554 -0.0153317162 -0.0111547977  0.0271932951  0.1188291669
#> [2721]  0.0551724188  0.1090040017 -0.0846197323  0.0477745914 -0.1639217783
#> [2726]  0.1781168886  0.1977861047  0.1087087073  0.0303526790  0.1517940929
#> [2731] -0.0008550505  0.0441469356  0.1811569969  0.2428165922 -0.0853340424
#> [2736]  0.0876940599  0.1572366761 -0.1319531073  0.0349070875 -0.0426307196
#> [2741]  0.2217230093 -0.0095982712  0.1486220204 -0.0391624742  0.0909591636
#> [2746]  0.2021274989  0.0197056165  0.2074564025  0.1862243784  0.3911335669
#> [2751]  0.1351761467  0.0295074752  0.1954128904  0.0768756863  0.0694731744
#> [2756]  0.2967789242  0.4236690935  0.1254366985  0.4029667312  0.2894874744
#> [2761]  0.1640453946 -0.0919607606  0.1058192101  0.1683441625  0.2178533326
#> [2766] -0.0824178503  0.1701200328  0.2102245529 -0.0126145493  0.3086253238
#> [2771]  0.1786607635  0.3630067900  0.2043757665  0.0542818279  0.3539155754
#> [2776]  0.0887256337  0.2213259670 -0.0372467580  0.0079306760 -0.0626465740
#> [2781] -0.1541078928  0.0650573853  0.0847785762 -0.0330417518  0.0710813534
#> [2786] -0.0147019396 -0.0276441530  0.0634474195 -0.1111333482  0.3254201633
#> [2791] -0.0951086221 -0.0232604473 -0.0246118911 -0.0877264424  0.3168050290
#> [2796]  0.1362782490  0.0739032664  0.2197522546 -0.0425937476 -0.2095515067
#> [2801]  0.0678437972  0.1707985629 -0.0103206092  0.0425277104  0.0339860176
#> [2806] -0.2190354585  0.2275881835  0.0337066023  0.1420066149 -0.0895496110
#> [2811]  0.0798104691 -0.0111307220  0.0395289268  0.0926984830  0.1479287690
#> [2816]  0.0969945359 -0.0137579007  0.0631769781  0.1806968906  0.2485184790
#> [2821]  0.0760674548  0.1202554464  0.0020571561  0.1806134774  0.0314302173
#> [2826] -0.1182414643  0.0162927482  0.2502915539  0.2123648356 -0.0284494648
#> [2831]  0.2941856150  0.0822276472  0.0388386510  0.1314764840 -0.0498737964
#> [2836] -0.0505054242  0.0396952113 -0.0504258851  0.0400926831 -0.0987614379
#> [2841]  0.2557543468 -0.0344874332  0.0542958974 -0.0014377634  0.0542480056
#> [2846]  0.0102394644  0.0119494186  0.2896741403  0.2359370145  0.0053454316
#> [2851] -0.1108913506  0.0785354496 -0.0037407278  0.0214753280  0.0023882202
#> [2856] -0.0787110303  0.3639878352  0.0704607101  0.0861434373 -0.1450797663
#> [2861]  0.0718955955  0.0482339019  0.0336974505  0.1506472188 -0.0552142450
#> [2866]  0.0297452238 -0.0324749326  0.0176135653  0.1222558789  0.1366532870
#> [2871] -0.0216712136 -0.1079949027  0.0571283628 -0.0481411323  0.0443367235
#> [2876]  0.0105995153  0.0848515904  0.3366419177 -0.0853036922  0.0571826857
#> [2881] -0.0225663678  0.2410223127  0.1464579951  0.2973785586  0.1123437084
#> [2886]  0.2138121108  0.1272851960  0.0475063058 -0.0186682403  0.0599779651
#> [2891]  0.1531921614 -0.1063655400 -0.0490559450  0.2654579044  0.1566451870
#> [2896]  0.0703079806  0.1099377282 -0.0314102977  0.0360805914 -0.0385918849
#> [2901]  0.1742892117  0.1941230160  0.2012820865  0.2202260739  0.0499513228
#> [2906]  0.1513793268 -0.0584552051  0.1869086455  0.0385875380 -0.2061367165
#> [2911]  0.1009392207  0.1179103450  0.1554538959  0.1855383979  0.0386714094
#> [2916]  0.2084656339  0.0647137716  0.0026347353 -0.0161765176 -0.1144441210
#> [2921]  0.1641959993 -0.0253528831  0.1907771617 -0.0779682500  0.0802796068
#> [2926]  0.2169309120  0.0039958109  0.2513431139 -0.2225443857  0.0681088886
#> [2931]  0.0699427872 -0.2143237191  0.2612393972  0.1597777623 -0.0186781879
#> [2936]  0.0287119914  0.0311502656  0.1939486760 -0.0090043406  0.1694124011
#> [2941] -0.0361138135 -0.1876495305  0.1746580661 -0.0890407800  0.1944371753
#> [2946]  0.1939065196  0.1556861296  0.0022742807  0.0895500057  0.0513082645
#> [2951]  0.1776092490  0.0008502327  0.2293818823  0.1676746797 -0.1336546839
#> [2956]  0.0542076717  0.0912143543  0.1104775861  0.1141964477  0.2198435175
#> [2961]  0.0777509215  0.0977215805  0.2062713861  0.0045426326  0.0622427839
#> [2966]  0.0387793900  0.1574251736 -0.0635327956  0.1477321062  0.1412691944
#> [2971]  0.3079554341  0.0762022833 -0.0458851199  0.2533155773  0.0813350496
#> [2976]  0.1805541410  0.1423090962 -0.0493938323  0.1894518912  0.0900248741
#> [2981] -0.0064940886  0.3396920037 -0.1640421381  0.1668346506  0.1021581818
#> [2986]  0.0113275553  0.2864777672  0.1928625819  0.0119584168  0.0049904330
#> [2991]  0.0085458057 -0.0389187932  0.1942083987  0.0572576845  0.0444944170
#> [2996]  0.1120639762  0.0328536242  0.0260570229  0.2499186351 -0.0206109718
#> [3001]  0.0038512833  0.2388155349  0.3129195497  0.0602096662 -0.0163733196
#> [3006]  0.0329060233 -0.0266118310 -0.0107328798 -0.0110732341 -0.2020401511
#> [3011]  0.1649075884 -0.1345687719  0.0491368048 -0.0025316221 -0.0129793427
#> [3016]  0.2439621367  0.0596271922  0.2756845154  0.2405428879  0.0161389323
#> [3021] -0.0409860186 -0.0159351962  0.1851602705  0.1050870058  0.1241200169
#> [3026]  0.0551597446  0.0184197150  0.0191373745  0.0669343566  0.0380910674
#> [3031]  0.1688571657  0.1084204730 -0.1242418502  0.1191623255 -0.2350910811
#> [3036]  0.1047124818 -0.0294894445  0.0444446357  0.0405510181 -0.1910701215
#> [3041] -0.0114366518 -0.0276466577  0.0267056446 -0.2150207104  0.1029768711
#> [3046] -0.0117078053 -0.2328146798  0.1060214913  0.0949729356  0.0227684166
#> [3051]  0.1561849091 -0.0498222786  0.0329124298  0.1494407936  0.1143207900
#> [3056]  0.0893493008  0.1792659795  0.1146357521  0.0458382850 -0.0195132432
#> [3061] -0.0373398389  0.1393574330  0.0980816933 -0.1492254391  0.0319805453
#> [3066]  0.0076287780  0.0288847570  0.0999084569  0.1068670554 -0.2105478930
#> [3071]  0.1316606591  0.1445838804  0.3105416192  0.1746012738 -0.0292872928
#> [3076] -0.0214958820  0.0168991698  0.0796290352  0.0468523267  0.0847099202
#> [3081] -0.1152460067  0.1039152369 -0.0450763468  0.0122304570 -0.0625589078
#> [3086]  0.0984752935  0.0133045212  0.1207777615  0.2379349308  0.0658131420
#> [3091] -0.0197187860 -0.1236259588 -0.0190792957  0.1484300149  0.0877155643
#> [3096] -0.0301172136 -0.0086784645  0.0231844425 -0.0007917662  0.0646937626
#> [3101]  0.1125869186  0.3220374936 -0.0310233948  0.0907831992  0.0574117574
#> [3106] -0.0983083863 -0.2449926955  0.0607999291  0.1797127646 -0.0382212598
#> [3111]  0.1287487780  0.1048259132  0.0260393883  0.0424786366  0.0885793926
#> [3116]  0.1258884803 -0.0106537628  0.0036013336 -0.1698608917  0.1851736754
#> [3121]  0.2160823015 -0.1901893349 -0.1706867078  0.2605206067  0.1750003850
#> [3126] -0.1254694078  0.1399355983  0.1652910176  0.0141181500  0.2683792292
#> [3131]  0.2346777519 -0.0093001554  0.3051531037  0.0401285894  0.2051422940
#> [3136]  0.1587162975  0.3239098525  0.2563027834  0.1716315713  0.0566951445
#> [3141] -0.0681214647  0.1239193610 -0.1660170593  0.1503378441 -0.0044547572
#> [3146]  0.2716684534  0.0625300022 -0.1166021865  0.0424876209  0.3845641428
#> [3151]  0.1987729628  0.2417895203 -0.0695851431 -0.0173723465 -0.0909772990
#> [3156]  0.0486991900  0.1241432671  0.0835274697  0.0436168856  0.1986059262
#> 
#> $mod1_unimodal$sim_properties$cell$pca_dist
#>    [1] 51.00288 48.62139 44.28552 51.40935 47.04041 42.78467 50.39172 44.90481
#>    [9] 55.81469 35.85980 47.95673 51.18518 48.13931 51.72267 56.22688 53.29813
#>   [17] 44.37666 50.51072 48.81291 52.27351 44.25791 50.10263 47.76405 57.19027
#>   [25] 50.78360 68.18656 45.22848 48.43569 58.53570 37.91282 42.84128 51.18851
#>   [33] 42.74205 48.60520 49.20288 43.42606 42.49998 42.17740 58.25431 46.34688
#>   [41] 47.96679 46.72602 49.80950 53.43760 50.73492 49.48387 35.13628 56.18952
#>   [49] 44.77022 38.73228 55.07283 50.90880 46.47453 47.94873 42.98074 44.46531
#>   [57] 49.04402 54.07419 49.35498 48.91807 38.20346 49.94309 51.51481 44.06213
#>   [65] 50.33058 50.36666 49.43276 53.30322 46.21368 50.75245 45.96298 51.32926
#>   [73] 51.92761 53.92155 46.26525 48.83304 43.27266 52.74985 47.17994 43.09600
#>   [81] 55.14518 53.95094 50.95990 59.76847 56.56597 51.07644 49.13121 51.64166
#>   [89] 49.06886 49.53599 48.30412 50.11794 56.23708 50.47793 55.40225 60.79367
#>   [97] 45.35392 48.14118 35.02544 40.43025 55.72271 47.04478 46.19507 56.83227
#>  [105] 47.79484 47.51024 60.74947 56.49351 33.15947 42.88448 51.34813 47.03571
#>  [113] 52.18500 55.31705 34.46953 52.08830 59.45464 45.58691 64.55514 54.93769
#>  [121] 43.08868 52.61425 51.84856 52.51860 53.55385 54.83501 61.16755 45.88949
#>  [129] 48.46784 51.63171 42.89132 56.30171 49.89367 50.03067 57.06204 54.59549
#>  [137] 52.87704 46.95980 35.81500 46.13605 61.94031 46.74167 43.98442 53.90361
#>  [145] 46.29826 57.28746 56.36131 57.76482 33.70627 55.92717 46.26003 46.43624
#>  [153] 40.79054 44.08063 54.32212 49.81446 56.45620 51.04875 60.61439 52.41485
#>  [161] 50.72683 51.86575 47.71767 65.59043 43.07085 49.40421 50.06768 43.09106
#>  [169] 59.75144 52.83111 45.46059 45.76497 55.91510 33.57358 54.70534 35.84357
#>  [177] 41.89886 58.04673 51.01426 49.52878 65.49439 42.96436 45.26387 52.91399
#>  [185] 53.46782 41.73867 49.79257 46.40529 49.46561 49.51157 48.40657 50.13863
#>  [193] 44.39991 51.17816 41.99376 48.17276 44.89766 46.04287 51.78595 43.28466
#>  [201] 45.34477 51.25856 52.61250 49.94091 41.61003 37.22040 42.37713 33.14605
#>  [209] 42.57852 45.94118 53.67496 54.21167 55.36491 48.10055 51.77233 34.84860
#>  [217] 37.78796 45.12865 50.18644 41.58911 54.08376 33.71394 60.73830 43.92843
#>  [225] 55.54949 35.99654 45.95054 46.23425 48.48039 43.16121 40.01715 58.00357
#>  [233] 51.51245 50.22623 48.06786 46.55300 46.65260 49.83812 53.61272 66.01003
#>  [241] 51.30021 62.14999 55.48412 44.70844 61.30259 53.22985 60.50258 53.14552
#>  [249] 45.10775 49.17443 50.44028 39.71947 53.41115 54.79199 57.00409 55.13204
#>  [257] 58.06089 54.25851 42.07006 45.58617 44.45467 44.64483 54.79856 46.90301
#>  [265] 40.45172 54.15599 49.59752 56.58242 44.83783 61.82130 56.45257 47.73721
#>  [273] 59.17089 47.18045 65.12764 48.86039 53.53897 44.30243 58.77475 53.63239
#>  [281] 49.99998 58.47813 51.74175 53.70794 57.19059 46.35355 46.50085 43.79660
#>  [289] 67.02152 50.56680 56.58117 31.82700 41.24039 57.62045 55.50168 50.45027
#>  [297] 57.44122 56.94991 46.49460 57.34117 57.64983 44.66657 50.14158 56.04424
#>  [305] 45.83164 58.86921 49.30574 56.03753 60.69409 50.31063 34.56068 53.80598
#>  [313] 50.87501 53.84595 58.60549 58.43391 58.78175 50.29267 48.37893 61.25719
#>  [321] 58.84081 65.63547 59.80501 59.33242 55.73810 55.07690 50.39297 57.35858
#>  [329] 58.00328 59.66795 49.38180 67.14031 53.45346 54.40586 59.10747 45.29703
#>  [337] 51.36002 48.84326 40.58976 48.49640 56.26825 62.11656 52.76047 56.41317
#>  [345] 59.66817 60.92836 55.65051 46.40303 51.18011 57.73485 43.72854 55.19271
#>  [353] 59.90207 63.16013 51.11651 44.25936 50.88246 57.93386 59.67600 57.80355
#>  [361] 49.83653 61.50255 56.87645 62.67973 56.32013 56.10380 51.14702 58.66581
#>  [369] 50.49506 51.14147 47.17136 56.04780 54.66055 56.51979 62.71006 54.19657
#>  [377] 54.80331 54.71791 61.38474 48.46189 58.28147 42.31534 53.32574 54.39865
#>  [385] 50.76430 52.13511 49.81840 47.76348 52.56856 47.57210 47.54825 50.03619
#>  [393] 44.68213 47.97127 51.29547 59.99255 43.79514 58.24262 49.88951 50.02017
#>  [401] 41.71603 43.53534 53.74264 52.69478 54.93383 61.36421 42.97673 42.60885
#>  [409] 58.95570 41.84619 51.40547 50.90930 39.42413 52.18312 38.36909 51.38761
#>  [417] 52.23286 47.31609 52.18338 52.62434 49.68510 45.54857 44.84826 48.05605
#>  [425] 44.63263 45.00323 55.88639 49.45942 41.85179 49.31180 46.84925 51.02234
#>  [433] 50.34231 52.83359 49.84512 48.84770 51.26869 55.18187 49.77219 40.46302
#>  [441] 43.70677 45.07506 47.09004 46.44258 43.25044 39.29213 47.12444 56.81224
#>  [449] 50.91446 54.59844 45.67641 55.19163 56.43695 42.45077 54.70424 42.10100
#>  [457] 55.31905 53.36719 45.64092 42.84727 44.33323 64.73395 49.23288 56.78952
#>  [465] 46.61830 50.89583 57.53901 50.90448 62.79775 51.31343 53.34028 47.16478
#>  [473] 48.81284 45.68527 51.58163 54.56261 58.84705 49.64560 69.73366 45.90618
#>  [481] 39.65823 50.82349 37.32835 49.53476 53.73783 44.55443 51.70617 41.58658
#>  [489] 40.05250 46.59576 40.48223 61.97073 47.62619 59.13026 44.75916 41.82325
#>  [497] 61.78446 47.14513 48.13148 53.23899 48.92473 44.97068 41.36850 48.34357
#>  [505] 38.28363 48.43032 53.04878 50.01807 47.80711 38.39414 49.78401 48.91430
#>  [513] 57.23092 44.88712 50.18483 41.99864 46.80547 38.47215 47.30188 47.75575
#>  [521] 54.86212 42.67067 43.82268 42.89812 41.48077 56.35455 52.59512 48.33612
#>  [529] 41.44438 51.32591 56.03556 42.03953 62.86916 67.44527 58.98994 58.44236
#>  [537] 45.07322 37.65395 59.07117 52.34227 52.12678 53.19230 54.84062 44.15489
#>  [545] 47.06277 47.92562 54.08019 55.39852 47.94023 55.85491 61.07706 33.23155
#>  [553] 38.46821 54.27463 45.88155 51.14724 51.99879 40.57096 51.99991 49.42428
#>  [561] 51.93484 53.93071 27.59903 61.98847 50.15854 57.23002 41.56285 45.68888
#>  [569] 57.50674 33.54317 48.74410 57.08290 54.82214 49.64857 45.52898 49.04332
#>  [577] 48.43055 43.13467 52.71308 50.60346 50.56841 45.99752 51.50903 51.72733
#>  [585] 46.45827 45.04328 55.41958 52.06416 34.87512 41.80103 58.83524 47.20806
#>  [593] 57.25280 53.60196 45.72911 48.73105 44.70178 53.83146 53.97663 48.48007
#>  [601] 51.18795 51.38069 53.89006 35.29889 49.50134 39.97385 43.39386 52.99843
#>  [609] 50.57576 62.45335 41.95874 65.42458 55.22858 52.67721 52.74112 49.70852
#>  [617] 47.24936 43.65752 59.77316 52.55980 48.83730 72.85383 56.26462 55.64750
#>  [625] 60.17270 52.41925 47.76992 59.47508 52.77131 55.79212 48.67440 51.08200
#>  [633] 48.76063 56.50087 57.05572 57.22089 61.76425 49.55919 51.57918 60.91958
#>  [641] 54.11904 54.90881 56.09639 47.94732 47.38543 45.91222 46.50472 48.79263
#>  [649] 52.50915 57.11370 52.55137 49.67223 54.67194 55.26288 56.04158 54.90570
#>  [657] 48.23155 47.55755 51.89278 54.43309 50.34454 49.16942 50.35024 51.02223
#>  [665] 51.36147 51.40085 41.18050 52.95752 55.28302 54.93485 50.62321 48.03004
#>  [673] 58.18661 60.99106 53.95382 63.59430 58.15326 58.58739 55.72032 56.23061
#>  [681] 46.01148 59.33195 55.23341 60.84968 57.27084 56.13480 53.41526 57.49019
#>  [689] 62.44352 46.47298 58.73086 59.06200 57.01784 57.81952 63.80277 56.56168
#>  [697] 59.27034 59.49745 60.23306 66.16064 55.50556 62.96327 53.43337 67.42305
#>  [705] 64.05772 59.43806 65.10841 64.43786 57.40238 56.59514 63.83818 58.67725
#>  [713] 56.81726 59.43500 58.70011 60.22320 62.05485 56.80841 64.81265 60.64896
#>  [721] 57.22525 53.42634 64.41071 56.76478 60.92094 60.41210 52.94530 64.15998
#>  [729] 61.68086 55.34894 52.48826 64.46895 49.91656 64.54206 60.06925 62.96102
#>  [737] 50.32890 68.31659 55.76367 57.46374 54.71484 57.37816 56.06252 61.11506
#>  [745] 63.49204 50.57177 50.29833 50.03483 58.83770 54.82726 62.88126 51.59036
#>  [753] 57.61485 53.75346 56.16685 40.15281 44.56242 59.56812 53.82826 55.28747
#>  [761] 67.40002 53.45499 53.55801 64.23825 51.32412 46.60864 60.82973 52.14874
#>  [769] 55.23580 46.36802 50.00985 45.79233 45.44442 64.77259 54.97727 53.44608
#>  [777] 48.30227 53.52451 57.72256 55.10974 58.46926 41.98778 56.24834 50.67552
#>  [785] 39.56118 55.78803 55.18251 41.31731 44.93949 57.58654 56.08998 52.53897
#>  [793] 56.56284 47.87462 53.03847 45.12024 41.07685 48.20254 53.53669 55.87836
#>  [801] 41.25333 52.38989 56.35887 48.85816 52.17139 45.30440 52.84037 52.83238
#>  [809] 51.93485 53.58508 48.58761 53.30923 53.40718 55.39683 55.30227 50.10693
#>  [817] 61.40145 57.54987 55.63899 54.73847 55.81184 50.12835 57.58527 53.91562
#>  [825] 42.05735 51.25643 49.02840 55.62687 70.34099 47.14449 51.42415 68.26632
#>  [833] 57.77595 48.90717 55.16107 44.89274 60.12599 47.05536 53.47682 45.93172
#>  [841] 49.51181 49.76446 53.92202 63.26054 49.24178 48.44092 44.39795 57.54251
#>  [849] 52.83465 60.17011 46.86412 45.27149 56.10000 43.17415 47.69179 56.46298
#>  [857] 55.02911 48.66110 51.95812 54.46023 50.12851 57.84277 39.90137 43.22224
#>  [865] 53.33461 54.08287 52.03346 51.71495 43.00753 50.77448 61.85550 56.83732
#>  [873] 48.67206 43.14729 48.50391 48.19128 50.36135 49.68325 47.74665 53.09271
#>  [881] 66.74461 56.60705 42.37248 58.09937 49.88774 56.13170 54.20512 60.64089
#>  [889] 55.39397 55.18403 44.55720 49.43464 43.74721 46.55083 41.88466 59.08714
#>  [897] 41.28898 50.05597 60.91772 45.50684 51.40037 50.52962 46.78496 44.57532
#>  [905] 49.85036 46.65846 45.30377 45.18332 57.15874 55.88506 57.36675 49.49938
#>  [913] 47.58423 57.13760 43.08582 49.47277 54.75031 44.59107 50.84754 48.69928
#>  [921] 38.18305 48.83967 48.31545 43.85540 51.98818 49.92144 44.06834 44.76297
#>  [929] 44.36182 51.68176 47.12052 45.58524 39.51565 41.34956 36.48398 54.24216
#>  [937] 43.48084 60.33615 45.89909 45.40899 50.44972 55.95595 32.71166 45.91302
#>  [945] 51.85590 45.17906 52.10152 42.72865 43.69561 61.97837 47.75220 44.26891
#>  [953] 54.32325 47.53509 40.22732 43.35760 45.64046 40.07141 50.32854 30.78661
#>  [961] 54.22612 55.92153 34.84459 43.33449 46.77643 42.75575 49.57687 54.44613
#>  [969] 40.42136 50.47323 48.47600 51.16490 49.22546 41.29591 47.86889 51.25293
#>  [977] 52.68215 48.32377 42.83148 54.73661 44.33905 52.80193 45.91071 54.86683
#>  [985] 48.40912 43.23894 38.30927 48.31117 46.65720 41.27031 46.56848 50.52460
#>  [993] 44.45569 51.22218 46.32135 49.32325 40.36423 45.59097 47.05638 40.98552
#> [1001] 42.56118 54.13877 39.85584 41.71182 44.75160 48.82992 38.80147 54.13731
#> [1009] 51.15204 45.60284 51.44002 44.74562 45.96948 50.80035 32.79157 62.28234
#> [1017] 71.78938 53.64568 71.52979 46.06364 57.47320 37.15789 57.32040 53.53042
#> [1025] 67.64841 69.29011 64.28359 52.76285 49.55163 70.19153 51.81913 59.09716
#> [1033] 52.92798 59.01856 66.98967 48.59015 55.39305 49.88875 49.94316 63.37217
#> [1041] 59.39212 56.62282 56.92640 51.06426 56.31127 57.02290 48.80862 52.50222
#> [1049] 54.40039 64.67431 58.83052 60.12543 51.11569 57.03012 59.88828 53.76012
#> [1057] 49.12401 65.00179 53.32643 48.92216 58.50813 46.23110 60.14743 61.22948
#> [1065] 56.12309 45.90995 51.54022 57.04490 65.24453 54.94531 66.74189 46.46852
#> [1073] 63.00839 56.73302 54.11948 45.49778 57.58210 63.79046 43.78765 58.29623
#> [1081] 60.95886 49.40156 52.29864 50.01753 47.18195 42.61923 53.94887 53.71547
#> [1089] 43.91136 51.15554 52.83865 52.32210 46.90753 53.16538 59.58881 54.60115
#> [1097] 60.44501 55.20998 55.55240 46.17003 51.64595 58.96619 53.10347 52.20564
#> [1105] 55.25980 58.06893 48.64466 51.52492 58.16082 52.65666 49.39483 56.89359
#> [1113] 51.44325 42.20318 48.95962 52.71604 49.00617 51.89077 54.96932 57.25978
#> [1121] 45.55256 52.54264 44.03924 47.25831 53.40939 46.58778 48.36794 53.01201
#> [1129] 48.33220 41.99097 57.79084 48.84433 51.79684 55.21801 49.18796 43.57061
#> [1137] 51.72825 50.51842 51.33588 50.63294 55.32354 57.47660 54.47248 51.04841
#> [1145] 54.81015 51.79513 50.27309 60.92755 57.37595 51.63529 63.32774 43.03580
#> [1153] 50.80134 66.72179 39.19031 46.88315 54.29747 56.04201 53.63352 52.18329
#> [1161] 56.56911 48.30319 55.42175 54.56850 56.83318 52.72864 57.55762 37.45124
#> [1169] 54.30897 57.56261 53.30733 46.45899 54.81441 51.17369 54.22664 59.34421
#> [1177] 54.73880 52.95901 51.26556 51.09416 49.14826 52.77264 50.07378 62.44533
#> [1185] 53.15378 55.53813 56.22949 53.32598 43.76914 53.49406 57.80805 48.27693
#> [1193] 49.28873 64.67754 36.11933 65.25539 49.81879 51.16669 52.23007 58.77949
#> [1201] 53.83715 52.94493 54.47030 54.01259 53.18536 65.01674 55.34267 64.40678
#> [1209] 48.43772 50.56717 44.43857 52.65446 53.25959 57.26750 57.51979 58.41060
#> [1217] 41.08896 42.81034 62.73756 50.90104 43.04404 45.06741 43.16234 61.25405
#> [1225] 42.10052 47.49547 53.91770 53.78589 54.72426 41.10586 49.74905 48.60845
#> [1233] 50.07467 54.95187 56.16406 35.89065 54.26828 54.82065 42.82811 54.40746
#> [1241] 56.26860 50.42136 50.49472 44.22163 53.56497 37.86865 53.02479 57.09575
#> [1249] 58.60670 47.34497 42.37148 49.11937 46.16717 44.13496 48.02860 53.59186
#> [1257] 42.90461 65.71992 49.92029 59.80740 43.99266 55.95353 47.10095 49.86453
#> [1265] 49.36196 48.35232 57.52536 58.57378 50.56637 51.02217 63.59024 54.13135
#> [1273] 63.17018 67.48674 52.53823 50.85563 58.02067 54.21519 48.93572 43.40290
#> [1281] 58.82792 54.91677 59.43219 53.70027 48.22491 53.16800 52.24492 53.79795
#> [1289] 50.33504 60.51384 64.88459 49.21216 61.38171 41.62212 52.24935 51.28801
#> [1297] 53.06952 43.70074 51.09509 47.22135 51.81107 55.05120 38.40833 58.18447
#> [1305] 60.51846 48.76028 59.58706 55.54181 62.19450 58.11803 61.10998 43.22817
#> [1313] 45.93803 61.88448 50.59765 52.68740 58.47024 56.68662 52.97985 55.42342
#> [1321] 52.97888 52.29869 53.53821 59.76277 50.24244 58.17421 52.38590 50.56152
#> [1329] 66.17962 51.27511 44.91429 34.31523 44.54552 55.83872 49.26351 57.60037
#> [1337] 53.78362 41.71751 39.28651 46.12847 48.68174 46.87214 46.48729 42.65916
#> [1345] 55.59786 46.56403 48.69754 48.50174 40.26507 47.03159 45.27940 46.43093
#> [1353] 40.62530 39.91717 44.36993 44.23753 41.59032 43.32665 50.63834 50.56503
#> [1361] 44.50963 42.15995 33.24430 41.30780 40.12623 36.42192 51.31517 50.11103
#> [1369] 50.56138 44.46243 51.28065 31.68423 48.58270 47.75826 50.30149 36.96607
#> [1377] 51.77442 34.02192 52.63442 44.85389 51.97734 31.07805 38.05217 52.32944
#> [1385] 46.72719 29.17133 39.86942 45.61490 44.81215 45.04926 46.94074 34.53486
#> [1393] 59.84796 46.60038 55.83157 58.07263 50.84652 45.09860 64.12711 49.44032
#> [1401] 41.24128 51.16504 52.36792 57.68760 47.86113 57.13942 54.85058 52.36221
#> [1409] 53.72475 47.48821 66.83551 46.59542 52.47492 60.80064 59.48111 51.81096
#> [1417] 54.36770 54.55515 56.74198 51.20659 53.59336 54.77177 47.63611 50.86165
#> [1425] 47.80727 41.71280 40.32454 57.39085 59.42621 49.59551 47.90746 54.35857
#> [1433] 61.70511 50.85006 51.91712 53.74269 49.33419 39.21757 50.70799 47.24839
#> [1441] 36.20153 51.04116 59.31971 50.95350 39.25858 56.11646 50.22561 52.91820
#> [1449] 46.10069 43.90893 48.00575 49.59414 48.22590 48.61183 46.79365 34.11077
#> [1457] 51.99014 46.32673 39.42405 47.02242 45.77759 46.03038 39.49134 46.97708
#> [1465] 42.00271 37.16976 54.00681 46.72801 46.45505 46.76623 38.34357 49.52075
#> [1473] 42.72233 44.44141 41.21886 49.93420 49.37718 37.96766 47.09122 40.57200
#> [1481] 36.08098 46.76465 44.51908 41.93424 49.77768 47.87159 34.84459 50.02874
#> [1489] 25.82958 36.23150 48.25377 46.81814 33.29911 44.51897 44.43705 53.09153
#> [1497] 46.50387 53.39618 28.16908 43.91755 42.03745 36.86351 40.04842 40.31632
#> [1505] 53.30591 36.66655 50.11821 56.61923 40.06095 51.60437 64.57530 47.68970
#> [1513] 44.13881 64.52256 45.13337 45.19749 56.39385 50.11786 52.12030 43.58907
#> [1521] 53.06938 51.48708 47.03629 49.62477 42.46882 64.71623 49.08995 50.82220
#> [1529] 50.31347 56.98089 54.94757 54.40793 47.86145 53.34219 50.72245 41.00322
#> [1537] 51.64095 39.04889 47.20265 44.50197 47.88027 37.68502 53.41943 51.37761
#> [1545] 41.75834 42.54874 43.90897 54.48834 54.96430 48.44419 42.43912 40.22150
#> [1553] 51.00275 41.27397 46.33083 34.06730 49.77496 52.99422 46.17312 43.63751
#> [1561] 47.82084 50.59038 53.54263 50.93095 59.30242 59.30887 65.32226 46.97424
#> [1569] 50.67942 63.16325 51.94835 57.13528 56.38167 48.37906 55.21016 53.61070
#> [1577] 46.77294 50.40384 50.90011 46.85619 53.52705 61.33873 54.61546 51.74774
#> [1585] 61.38358 60.57745 53.29471 56.97446 61.22175 50.09859 53.82992 51.33729
#> [1593] 57.80607 63.27628 54.89745 55.96160 44.81812 52.66801 39.45399 39.82441
#> [1601] 52.71989 43.12261 56.31107 51.31650 50.09716 40.80785 55.70105 57.08939
#> [1609] 57.25301 57.50931 56.80822 48.59488 60.51763 45.89329 51.78337 57.22408
#> [1617] 53.63713 56.72526 52.60328 50.63076 53.54925 52.66669 46.19395 44.41154
#> [1625] 59.06870 52.54026 52.43642 61.20497 52.52138 57.93792 48.17755 61.39560
#> [1633] 53.76346 47.92708 58.04516 52.29914 62.61315 56.82194 51.95956 51.33723
#> [1641] 55.54571 61.14662 52.94068 55.76837 53.27359 52.57334 49.22842 51.93407
#> [1649] 49.95073 47.40334 52.24080 53.65187 46.01980 51.32317 52.75414 51.75712
#> [1657] 44.63544 52.98525 53.53298 49.75255 49.75382 56.76962 42.58890 55.27915
#> [1665] 52.02344 45.25062 40.91524 53.25819 45.82382 45.56837 54.49528 53.07665
#> [1673] 49.46381 59.86146 54.01372 63.49573 55.20743 47.91229 62.06881 51.97599
#> [1681] 42.93489 50.85237 49.84550 34.79445 49.93775 56.15360 49.21121 55.11008
#> [1689] 53.32541 51.15310 61.71696 52.98442 46.81984 48.15684 49.90772 51.14942
#> [1697] 57.30089 51.21952 42.68902 41.63616 41.29507 46.52601 44.90307 57.68590
#> [1705] 43.82985 53.14767 49.65609 54.84264 56.98685 55.76750 48.31043 45.23240
#> [1713] 56.05674 45.66431 41.95304 56.01490 49.04721 62.39811 47.08979 45.52600
#> [1721] 50.13304 51.38403 48.38699 44.19472 50.40970 45.88563 51.40427 51.03886
#> [1729] 54.20495 60.61481 52.97039 59.23790 69.38710 57.46986 60.77700 55.95212
#> [1737] 63.27618 57.03919 62.43450 58.54004 61.61574 63.35551 69.10633 64.78055
#> [1745] 67.75765 48.64506 57.09226 66.43839 58.08021 51.79599 61.45677 59.82854
#> [1753] 65.70965 66.67376 55.37292 60.75068 52.43390 61.87204 57.29613 57.21605
#> [1761] 64.73845 60.78180 59.55080 50.30069 51.46261 66.47728 56.81154 56.05871
#> [1769] 63.10464 59.96384 64.50380 64.77362 66.70600 53.35526 62.93417 55.87314
#> [1777] 49.73759 56.08348 56.00356 49.30261 54.66637 59.24950 39.51891 47.98553
#> [1785] 47.20608 47.38785 39.69132 39.24946 50.61966 41.81961 35.38549 41.02902
#> [1793] 39.09478 52.37968 38.22649 50.88337 40.35806 36.35244 50.77266 40.11916
#> [1801] 33.93806 53.22074 46.09848 45.53211 45.95846 40.99915 42.38290 45.38855
#> [1809] 44.29487 53.06274 48.48585 47.81634 41.44231 46.74808 38.47983 39.33068
#> [1817] 45.75558 41.97661 26.08253 34.35971 49.73706 33.47063 55.59827 45.26916
#> [1825] 48.52050 42.89529 54.42266 47.67010 49.62465 45.46372 40.73709 49.95663
#> [1833] 53.69766 29.65973 52.30990 41.05404 43.90167 48.42096 43.86379 49.78842
#> [1841] 25.41970 48.62179 47.59639 33.51444 50.76721 38.22109 51.61286 50.33414
#> [1849] 32.69461 44.32406 48.75182 42.64303 48.78099 49.62884 41.03461 43.82018
#> [1857] 46.57174 34.21367 43.32654 52.51285 43.41016 39.43027 40.49048 45.76445
#> [1865] 44.14492 49.53652 30.55433 43.33370 49.95012 41.48150 32.75279 40.86713
#> [1873] 45.41529 54.32565 41.16152 45.66049 29.66415 42.55940 48.83983 31.25020
#> [1881] 48.32052 43.78387 52.95727 52.59530 45.95394 58.80397 61.36332 59.28990
#> [1889] 51.81963 47.52333 60.91186 41.25303 56.53996 55.41666 64.48067 55.75589
#> [1897] 55.93853 58.63109 43.23072 65.73547 49.15569 53.25690 54.35180 61.61563
#> [1905] 55.47391 51.16572 57.72825 48.00205 57.78926 58.75264 61.23001 66.38823
#> [1913] 54.23768 59.81777 48.28262 60.47656 43.91489 47.13319 54.12277 54.29937
#> [1921] 43.50659 65.55143 46.84765 56.47697 59.78288 62.55869 53.36899 57.20391
#> [1929] 67.62410 60.99741 58.92295 42.02909 52.79042 64.62420 46.36043 54.99850
#> [1937] 58.32806 41.53517 52.52776 46.74456 53.47414 53.43124 42.31245 58.73325
#> [1945] 49.37065 54.22964 54.48786 50.48946 55.98252 48.37193 56.01687 47.40955
#> [1953] 57.11927 49.24820 46.86045 47.68236 49.25020 49.72966 49.84903 38.03480
#> [1961] 46.57976 40.93014 58.10590 49.52071 57.59638 44.82231 55.58409 45.63833
#> [1969] 51.55057 43.36307 52.79017 44.89194 59.22540 37.05621 52.67251 46.57663
#> [1977] 53.27239 55.64751 49.22563 54.37645 46.05378 48.82401 50.65568 46.07845
#> [1985] 33.37594 44.31109 47.26112 48.32923 52.33113 36.21263 51.02569 58.89356
#> [1993] 44.03793 58.46251 48.41887 41.95006 57.12867 52.90551 44.97911 50.07703
#> [2001] 54.44899 53.37824 43.08660 52.09465 48.74261 45.41775 50.73521 47.52407
#> [2009] 42.40088 50.23466 60.68892 60.85815 48.44444 34.73168 44.09512 57.27937
#> [2017] 42.42400 50.01731 51.69347 50.47858 51.13138 56.45578 51.56358 33.23357
#> [2025] 47.51991 44.66103 42.90597 41.75618 43.49542 52.63816 57.19498 51.97286
#> [2033] 50.42741 50.02417 50.28659 48.18080 38.84709 57.06035 56.98682 45.48227
#> [2041] 49.59481 44.00645 42.18284 53.64701 51.48382 28.74533 54.34983 50.25622
#> [2049] 58.07757 50.73280 51.55912 44.82250 57.06980 49.74977 50.45739 53.85985
#> [2057] 56.05916 55.13738 58.42929 53.01775 40.70635 52.35545 51.53055 42.39448
#> [2065] 43.77813 54.26277 45.99967 55.27218 58.27610 51.10321 45.09493 54.39016
#> [2073] 47.18533 43.49577 41.28125 44.10073 56.07551 55.30947 50.44407 52.15610
#> [2081] 48.42277 49.00935 44.99536 43.39707 51.04884 51.42049 57.58346 44.78035
#> [2089] 37.44217 51.02772 41.59174 48.28850 54.50100 56.03545 35.13967 44.91948
#> [2097] 41.53949 46.52596 47.25800 45.19904 44.54540 46.22852 47.57409 57.40587
#> [2105] 55.39214 41.00858 41.14585 48.26148 43.47444 38.86703 41.37484 51.13292
#> [2113] 44.23544 62.09258 51.25051 56.49160 44.24906 45.86142 51.77070 49.92414
#> [2121] 48.74643 29.11304 40.43751 53.74174 37.08969 58.19813 43.47485 52.98389
#> [2129] 48.67325 59.13486 49.77495 55.07881 53.61751 48.32018 55.99338 43.26148
#> [2137] 49.63647 53.56609 50.97237 52.38457 46.05600 49.14107 52.07288 44.90930
#> [2145] 59.64745 47.61747 57.72234 43.75543 57.30672 46.78553 48.37209 39.65752
#> [2153] 37.05297 58.68938 50.15844 42.88980 59.03472 48.99494 55.90575 54.64414
#> [2161] 48.40038 54.39463 51.80266 53.30623 50.09531 52.89193 47.02149 51.61753
#> [2169] 54.57670 50.81615 47.99548 44.42003 43.17913 47.77147 44.85621 53.19946
#> [2177] 46.76087 34.10222 41.24945 56.09404 40.85051 50.17419 44.39066 36.35044
#> [2185] 44.51376 45.49279 34.86028 47.69752 49.79515 51.45745 43.87541 46.21054
#> [2193] 40.23715 46.57571 50.68276 42.30082 42.90946 45.80706 41.44505 37.34787
#> [2201] 24.42589 47.94272 54.00900 36.73896 44.17015 33.52340 50.25696 51.84646
#> [2209] 31.39118 49.74120 44.73320 54.45674 51.21376 46.30742 44.26451 46.75988
#> [2217] 56.76183 45.56520 53.75892 47.48408 39.85072 58.97035 50.07772 35.84288
#> [2225] 50.50318 40.58892 47.62930 51.17231 50.47251 42.93840 49.87744 50.42442
#> [2233] 57.68375 50.59993 43.26582 44.38868 42.36351 43.13946 39.84063 36.09204
#> [2241] 46.79828 43.47196 38.87585 48.56264 43.08015 54.63724 47.49554 51.88455
#> [2249] 48.99366 54.94833 53.46454 56.60124 46.93911 42.50514 49.93048 56.54148
#> [2257] 38.49252 49.95900 60.74595 50.26777 61.36821 50.42274 30.67587 54.09799
#> [2265] 53.90586 47.80741 46.78455 51.78476 53.87196 40.43660 49.12986 43.08016
#> [2273] 52.15794 50.70160 55.45561 51.63844 55.08635 46.98499 51.29185 51.90739
#> [2281] 39.95368 47.39519 50.82582 37.00909 41.98148 42.87846 50.52642 53.79471
#> [2289] 52.47552 52.83561 37.12924 55.83526 47.58134 46.05344 44.33383 35.69886
#> [2297] 46.93131 51.25176 47.06581 60.86012 49.65670 51.60849 45.25808 43.62092
#> [2305] 50.18084 35.96167 50.91824 47.55316 48.58531 47.89305 43.77291 47.65846
#> [2313] 42.45257 35.99237 51.56330 45.15264 47.88020 43.33746 47.69473 40.30239
#> [2321] 44.36521 34.64027 46.47306 50.90400 43.54052 40.88846 47.73802 48.04740
#> [2329] 58.39967 47.24899 44.78773 42.73002 38.24720 48.12033 48.69375 46.75083
#> [2337] 49.34220 50.91001 50.68589 44.13805 55.21933 57.76008 50.30443 49.54937
#> [2345] 44.67766 64.48554 46.34814 57.97746 56.24104 40.30035 54.59729 36.86846
#> [2353] 49.07608 59.30522 50.68614 45.60115 52.66756 56.30770 47.75790 52.99324
#> [2361] 56.08252 51.25502 52.31358 56.36648 54.15591 44.06641 51.80644 50.88951
#> [2369] 52.55285 45.19364 55.59377 46.74100 60.31236 59.34427 45.63085 50.79763
#> [2377] 53.80448 54.93715 49.81285 47.81997 59.60943 44.21228 50.29263 51.77756
#> [2385] 56.79492 45.07622 57.01171 59.86256 50.29790 43.21629 52.76294 51.37391
#> [2393] 42.55859 51.23108 52.43595 52.60336 45.01944 46.14291 48.52431 46.80169
#> [2401] 39.61785 52.07960 50.76223 47.39333 41.23763 52.64572 35.33569 60.33182
#> [2409] 47.38109 47.32470 40.27798 49.03276 55.88152 51.78688 44.66971 46.48624
#> [2417] 53.14916 59.07319 52.92165 51.78486 55.36815 50.65411 49.51510 41.83420
#> [2425] 39.89155 56.36357 53.79529 53.22884 57.38509 46.11737 58.92276 45.83996
#> [2433] 47.64174 62.52242 61.39952 60.91906 48.27152 65.30149 41.95024 51.64720
#> [2441] 46.58267 57.68268 50.99117 58.60119 49.27478 60.53077 53.40775 56.65181
#> [2449] 55.79396 56.32479 52.57922 45.77795 57.05888 55.20638 61.06983 55.31906
#> [2457] 57.31975 50.75402 48.69192 43.86474 40.12599 56.35266 52.42753 43.54792
#> [2465] 36.61243 45.15024 48.41669 43.68887 42.31677 51.30690 59.34049 50.52855
#> [2473] 41.90470 46.49468 44.63737 48.07677 54.12346 40.97688 44.57060 41.83596
#> [2481] 47.07785 40.41181 53.23061 52.59353 36.99488 47.44457 38.59341 56.34897
#> [2489] 53.23576 36.01189 44.63145 47.90550 47.81607 43.60556 47.56255 46.40172
#> [2497] 38.15175 48.51645 44.02903 41.71626 42.99459 39.86611 28.12505 48.25565
#> [2505] 52.13322 49.62898 46.09858 49.55147 47.36051 48.79609 48.66719 35.51658
#> [2513] 35.83301 49.86686 32.91236 30.87384 39.45782 47.18380 51.54691 44.59357
#> [2521] 53.70060 33.49908 51.50692 51.87452 37.84331 45.93934 30.33855 47.87469
#> [2529] 50.14599 35.11340 56.61269 46.94181 50.31069 49.20428 41.69000 52.33970
#> [2537] 46.06426 40.20933 49.12203 48.88597 39.83720 60.93658 55.05655 48.36851
#> [2545] 51.71659 49.84034 46.01092 53.09048 54.69213 53.09330 45.52793 43.93826
#> [2553] 44.41826 67.58057 45.71041 46.45073 49.98118 48.13571 53.80604 40.73728
#> [2561] 44.40514 49.11072 44.68407 50.38577 58.72093 46.43352 59.17762 48.76096
#> [2569] 51.01091 45.30132 42.39796 45.54226 40.55653 56.57199 48.27980 56.87084
#> [2577] 56.75337 57.85788 50.53576 48.21742 43.52519 48.34475 49.07734 41.91860
#> [2585] 38.37487 61.92792 42.68974 65.32079 55.50328 54.75727 52.99469 47.41655
#> [2593] 51.51429 56.34016 49.35332 43.88719 58.03502 52.17748 42.79763 52.93688
#> [2601] 36.61708 44.58034 52.25382 44.82523 34.44489 51.06482 47.83002 48.15036
#> [2609] 46.48952 52.65501 47.47645 52.20728 46.45365 40.31045 42.83289 47.38036
#> [2617] 38.31813 36.12958 49.14928 40.24477 56.17462 46.98871 50.53954 45.47165
#> [2625] 52.97196 48.24028 44.89627 39.00043 45.26078 55.93959 51.39920 40.94717
#> [2633] 57.58569 54.85752 45.52942 59.48859 46.17033 51.54188 39.59051 43.66455
#> [2641] 54.22308 50.33884 57.31326 45.66560 62.46910 38.01046 47.61924 54.71225
#> [2649] 53.92249 53.10139 51.96049 51.49921 47.43357 44.40021 53.32043 45.17359
#> [2657] 55.21186 53.00846 45.59784 47.58513 52.70371 40.60636 45.63727 52.48287
#> [2665] 51.70564 62.89524 42.18516 34.19528 51.28177 54.97558 48.72138 49.09649
#> [2673] 47.27402 51.78567 56.09145 45.70519 46.98909 40.66416 54.47623 49.81448
#> [2681] 44.58706 46.91542 46.20573 60.93921 46.08678 47.46400 50.28039 53.57264
#> [2689] 46.20705 49.58030 43.98877 51.69947 57.34628 57.03222 48.65922 42.87927
#> [2697] 42.53765 42.63886 48.89728 50.28683 48.60422 47.64294 48.90570 46.03461
#> [2705] 51.10387 46.02805 47.06278 45.07809 44.94729 40.47421 40.24984 44.07328
#> [2713] 47.65662 60.74556 43.80738 48.01919 47.97099 48.04271 55.21641 46.41022
#> [2721] 53.22750 42.49586 46.41498 55.07301 42.62141 46.77021 48.26614 38.92060
#> [2729] 49.09089 51.16011 57.43382 53.43422 45.12069 41.21726 56.37291 42.99629
#> [2737] 48.20362 46.54455 40.45258 38.75335 47.81892 46.76611 48.93443 46.73765
#> [2745] 45.34739 40.16718 46.41907 57.81605 47.62635 46.61818 40.24860 47.42798
#> [2753] 43.24345 44.76998 37.54787 45.29035 47.91317 41.10948 52.29031 50.77702
#> [2761] 47.22280 50.16716 50.58468 46.57694 43.93221 44.10484 44.51340 32.26253
#> [2769] 49.31737 37.74319 58.02866 37.23848 49.45437 42.35596 52.05408 48.14042
#> [2777] 44.57451 46.03454 38.87147 54.21091 46.31416 39.19652 48.64637 48.00069
#> [2785] 39.27661 49.28748 48.20993 46.66761 48.40402 55.89589 36.17786 42.53694
#> [2793] 44.76707 45.05317 31.92185 43.65498 40.55461 57.13403 38.07862 46.08860
#> [2801] 36.81674 42.24814 46.80142 38.16038 40.05784 38.32817 50.07806 51.00864
#> [2809] 45.20371 48.72439 47.69189 52.32837 47.73017 54.50743 48.84943 42.61066
#> [2817] 42.89992 38.93672 55.57315 45.20840 44.04152 52.35903 38.73601 62.86235
#> [2825] 41.35847 51.46621 43.43729 45.16010 54.14278 54.27886 38.51276 46.63974
#> [2833] 49.75161 43.67040 47.08640 48.08526 58.00419 44.89125 52.74618 48.43792
#> [2841] 55.48720 45.66617 49.84063 33.47641 52.55005 48.98547 51.51339 34.76389
#> [2849] 57.12727 42.83173 48.07818 46.52614 50.07222 45.46389 46.73556 43.53572
#> [2857] 43.10930 39.13217 44.55730 47.93081 46.85870 45.67304 59.80508 54.33305
#> [2865] 54.85289 37.54174 53.11367 55.34305 55.47819 45.20720 56.17112 43.57151
#> [2873] 57.42905 41.95115 48.40761 43.46608 42.10643 48.60394 43.16892 39.04761
#> [2881] 49.14007 44.86513 48.20436 52.83793 49.24332 56.38192 59.10834 48.56553
#> [2889] 40.61660 48.62026 61.41393 43.34562 46.97308 55.31878 55.56531 56.24965
#> [2897] 48.98936 60.16271 37.42263 58.22685 44.76475 49.18973 49.16058 56.50338
#> [2905] 59.17900 55.86169 47.23067 52.50007 48.02866 46.67684 42.18978 45.27088
#> [2913] 48.58413 54.36305 46.48653 47.78120 44.11858 51.36661 44.11423 35.77425
#> [2921] 46.07019 41.92783 52.05471 46.61536 49.08605 47.17590 39.87352 57.28556
#> [2929] 47.18175 32.49589 52.18442 48.33802 55.30064 46.98946 44.82437 32.55696
#> [2937] 45.23570 47.59436 56.13543 49.27479 40.27939 46.16499 53.92901 51.25383
#> [2945] 50.02583 47.82021 52.56204 53.68971 48.11641 49.61349 55.21895 37.77543
#> [2953] 46.42007 43.31172 54.37661 32.61894 47.17319 45.87124 53.55726 48.55963
#> [2961] 46.93395 45.46135 49.79448 51.60134 46.55028 51.40427 47.58695 52.42217
#> [2969] 41.81750 51.78666 43.84658 44.65021 59.06216 42.98776 49.01125 49.70265
#> [2977] 45.91367 62.89129 58.39429 52.33382 50.26546 51.13256 53.94471 57.79487
#> [2985] 45.12323 48.38869 48.23204 60.56513 46.11202 36.43490 47.59479 45.48779
#> [2993] 35.46909 46.50984 40.40553 49.87549 48.77272 47.76002 27.86367 38.54471
#> [3001] 37.73518 35.58919 40.78426 38.03031 47.34073 51.32742 48.77241 52.56116
#> [3009] 46.08641 42.47633 43.58211 46.39183 53.20900 43.94489 55.33844 42.13765
#> [3017] 54.78733 49.63047 42.63994 51.24011 40.70976 54.03227 51.57396 43.89941
#> [3025] 52.19357 38.87292 48.92488 36.16123 64.63815 47.14056 47.24717 49.70980
#> [3033] 48.61994 46.69523 47.33232 55.26746 34.76186 52.76606 53.87159 49.63850
#> [3041] 37.72217 51.98899 46.81799 55.19852 48.23011 51.78673 44.94201 58.62719
#> [3049] 49.74302 50.03740 47.45113 45.01376 50.31964 50.42505 28.72101 46.49410
#> [3057] 35.12485 55.83093 40.29387 46.63918 36.70755 46.52812 45.31234 39.91667
#> [3065] 42.01959 33.33251 50.63729 39.69091 38.67024 54.91531 54.15714 45.38419
#> [3073] 45.80196 40.64433 51.39845 58.48816 37.43151 52.44286 43.88289 53.01320
#> [3081] 54.31094 52.23756 61.75112 38.78079 46.85406 42.97698 46.94210 47.16980
#> [3089] 50.47628 37.66252 36.43676 44.76951 49.31734 44.99602 57.64114 50.14912
#> [3097] 44.10139 61.85990 66.91366 51.13425 54.59206 60.69069 56.00576 57.94065
#> [3105] 46.20851 50.61125 43.55268 57.62610 55.36779 46.06372 46.93454 48.07073
#> [3113] 49.48863 42.17885 41.86650 49.30889 36.96630 50.67472 44.66740 46.40341
#> [3121] 53.42595 46.01035 57.14556 53.75590 44.08101 45.68012 35.07530 38.84604
#> [3129] 39.61714 51.52057 49.17716 46.94640 51.97726 49.39800 40.98758 45.00454
#> [3137] 47.07854 59.21981 58.96830 42.21678 50.48735 51.57769 58.55969 54.44841
#> [3145] 59.58149 53.03177 44.01575 54.52230 49.33246 54.40629 46.71020 41.38894
#> [3153] 40.93932 47.85256 42.82200 51.62242 43.76954 49.73675 46.94392 47.77946
#> 
#> $mod1_unimodal$sim_properties$cell$knn_hubness
#>  [1]  4  3  3  1  1  3  2  5  3  0  2  2  2  3  0  2  1  1  1  8  2 14  4  1  3
#> [26]  1  0 13  8  1  4  3  2  6  2  8  4  3  4  2  2  1  3 12  1  3  5  4  1  3
#> [51]  3  5  9  3  4  2  1  4  4  3  1 21  4  3  3 18  3  8  1  2  3 15  2  1  7
#> [76]  7  6  3  1  6
#> 
#> $mod1_unimodal$sim_properties$cell$prop_outliers_cell
#> [1] 0
#> 
#> 
#> $mod1_unimodal$sim_properties$feature
#> $mod1_unimodal$sim_properties$feature$mean_expression
#>  [1] 12.351221 12.148863 11.000753 11.955018 11.804892 12.410481 12.973868
#>  [8] 11.697096 12.289612 12.449871 12.450975 12.357034 11.315901  9.880961
#> [15] 12.239692 11.343633  9.773410  8.637613 10.474719 10.934556 10.253266
#> [22]  8.931317  9.697346 11.664570 12.194100 11.073026 10.967787  8.992583
#> [29] 13.115604 11.035537 10.856682 12.291518 12.402456 10.755332 11.334035
#> [36]  8.962758 11.031285 10.869918  8.763965 10.585361 10.838214 10.984303
#> [43] 11.087264  9.972900 10.733358  9.909553 10.487023 11.945403  9.934252
#> [50] 12.263342 11.165369 11.037627 12.244533 12.391237 11.116308 11.488850
#> [57] 11.443752 11.128507  8.688491 11.635368
#> 
#> $mod1_unimodal$sim_properties$feature$variance
#>  [1] 21.55092 27.81482 36.74076 26.97810 26.19284 23.92838 20.23390 28.23542
#>  [9] 28.50145 21.48720 21.45778 26.45556 26.28789 40.99848 23.04773 28.96139
#> [17] 40.05241 39.18586 28.55080 36.69644 34.25149 41.95426 41.75456 30.66618
#> [25] 25.11464 27.48012 31.67581 40.38019 10.72921 32.14445 33.42716 18.31830
#> [33] 23.46884 32.74777 33.74322 40.17714 27.40481 33.15085 38.16178 33.75451
#> [41] 30.92148 33.56867 30.18345 39.13627 30.20021 39.02026 33.58222 26.83123
#> [49] 36.62681 22.99377 27.84502 32.11458 23.23923 21.06392 30.10819 29.39716
#> [57] 29.23077 32.22107 37.40067 28.04191
#> 
#> $mod1_unimodal$sim_properties$feature$scaled_variance
#>  [1] -1.32319058 -0.37511351  0.97587929 -0.50175485 -0.62060940 -0.96334828
#>  [7] -1.52252833 -0.31145319 -0.27118665 -1.33283521 -1.33728800 -0.58084483
#> [13] -0.60622289  1.62031090 -1.09663978 -0.20157245  1.47711738  1.34595929
#> [19] -0.26371756  0.96917105  0.59911459  1.76497286  1.73474741  0.05645661
#> [25] -0.78380065 -0.42577106  0.20926990  1.52672824 -2.96111815  0.28020163
#> [31]  0.47434729 -1.81246511 -1.03290189  0.37151792  0.52218496  1.49599545
#> [37] -0.43717072  0.43252552  1.19095869  0.52389373  0.09509795  0.49576509
#> [43] -0.01660712  1.33845467 -0.01406953  1.32089525  0.49781721 -0.52398426
#> [49]  0.95863298 -1.10480700 -0.37054187  0.27568012 -1.06765469 -1.39689983
#> [55] -0.02799863 -0.13561723 -0.16080132  0.29179857  1.07576091 -0.34074089
#> 
#> $mod1_unimodal$sim_properties$feature$sd
#>  [1] 4.642297 5.273975 6.061415 5.194045 5.117894 4.891664 4.498211 5.313701
#>  [9] 5.338675 4.635428 4.632254 5.143497 5.127171 6.403006 4.800805 5.381579
#> [17] 6.328698 6.259861 5.343295 6.057758 5.852477 6.477211 6.461777 5.537705
#> [25] 5.011451 5.242149 5.628126 6.354541 3.275548 5.669608 5.781623 4.279989
#> [33] 4.844465 5.722567 5.808892 6.338544 5.234960 5.757677 6.177522 5.809863
#> [41] 5.560708 5.793847 5.493947 6.255899 5.495472 6.246620 5.795017 5.179887
#> [49] 6.052009 4.795182 5.276838 5.666972 4.820708 4.589545 5.487093 5.421914
#> [57] 5.406548 5.676361 6.115609 5.295462
#> 
#> $mod1_unimodal$sim_properties$feature$cv
#>  [1] 130.55767 103.95416 106.25960 111.95181  98.14919 106.41408  83.64384
#>  [8] 110.48257 104.52145  94.26704 118.44160 112.88781 104.81625 120.68695
#> [15] 105.83086 107.11749 115.32781 119.41742 100.01927 119.97535 120.98654
#> [22] 115.28308 114.52381 107.58913  86.06260 116.12402 115.17596 158.46266
#> [29]  99.13179  96.30449 113.85643  96.33931 103.01956 122.53741  98.64338
#> [36] 136.77820 100.91043  97.63758 125.49460 108.22896 109.76487  88.48904
#> [43] 118.68170 136.07973 107.51032 133.33971 114.64982 101.41153 118.91615
#> [50]  88.78895 123.30338 111.62377 100.75538 106.43936 125.32453  92.86932
#> [57] 131.71075 100.59417 115.23526  93.86104
#> 
#> $mod1_unimodal$sim_properties$feature$zero_fraction_feature
#>  [1] 0.1125 0.1500 0.2250 0.1500 0.1500 0.1250 0.1000 0.1625 0.1500 0.1125
#> [11] 0.1125 0.1375 0.1625 0.2875 0.1250 0.1750 0.2875 0.3375 0.2000 0.2250
#> [21] 0.2375 0.3375 0.3000 0.1750 0.1375 0.1750 0.2000 0.3250 0.0500 0.2000
#> [31] 0.2125 0.1000 0.1250 0.2125 0.2000 0.3250 0.1750 0.2125 0.3250 0.2250
#> [41] 0.2000 0.2125 0.1875 0.2750 0.2000 0.2750 0.2250 0.1500 0.2625 0.1250
#> [51] 0.1750 0.2000 0.1250 0.1125 0.1875 0.1750 0.1750 0.2000 0.3250 0.1625
#> 
#> $mod1_unimodal$sim_properties$feature$detection_freq_feature
#>  [1] 0.8875 0.8500 0.7750 0.8500 0.8500 0.8750 0.9000 0.8375 0.8500 0.8875
#> [11] 0.8875 0.8625 0.8375 0.7125 0.8750 0.8250 0.7125 0.6625 0.8000 0.7750
#> [21] 0.7625 0.6625 0.7000 0.8250 0.8625 0.8250 0.8000 0.6750 0.9500 0.8000
#> [31] 0.7875 0.9000 0.8750 0.7875 0.8000 0.6750 0.8250 0.7875 0.6750 0.7750
#> [41] 0.8000 0.7875 0.8125 0.7250 0.8000 0.7250 0.7750 0.8500 0.7375 0.8750
#> [51] 0.8250 0.8000 0.8750 0.8875 0.8125 0.8250 0.8250 0.8000 0.6750 0.8375
#> 
#> $mod1_unimodal$sim_properties$feature$dispersion
#>  [1] 9.460144 6.618962 5.476184 6.767932 4.166386 6.879304 4.433900 5.706492
#>  [9] 7.770216 4.631970 7.873446 8.601969 3.652993 5.298143 5.586085 4.403208
#> [17] 4.206284 2.442114 2.088304 7.376968 4.336431 3.089969 4.442945 6.062627
#> [25] 3.601542 4.096005 4.924818 5.869559 4.704729 3.466365 4.926049 3.631294
#> [33] 5.757569 5.292935 4.719301 4.396445 2.965774 3.324668 2.775743 3.733681
#> [41] 3.885585 2.848275 5.352432 6.133986 3.323060 6.133920 4.091251 5.077872
#> [49] 3.994847 3.675671 5.074243 4.921647 5.151961 5.537214 5.791676 3.460665
#> [57] 7.502890 3.971781 2.008474 3.766230
#> 
#> $mod1_unimodal$sim_properties$feature$bcv
#>  [1] 1.2346458 0.9578003 0.9606896 1.0335065 0.8556357 0.9837615 0.7360966
#>  [8] 1.0033624 0.9756412 0.8347340 1.1066469 1.0612340 0.8932484 1.0870234
#> [15] 0.9589139 0.9417162 1.0068986 0.9176657 0.7220412 1.1154759 1.0612358
#> [22] 0.9481093 1.0081505 0.9831639 0.7314517 1.0095845 1.0281974 1.4433412
#> [29] 0.8796783 0.8123404 1.0164501 0.8200824 0.9364674 1.1035641 0.8757089
#> [36] 1.2022045 0.8215496 0.8164382 1.0037498 0.9260806 0.9459150 0.7128232
#> [43] 1.0702228 1.2449428 0.8988998 1.2198741 0.9965807 0.9087894 1.0296227
#> [50] 0.7575424 1.1048741 0.9964046 0.9045001 0.9634989 1.1399287 0.7831027
#> [57] 1.2261967 0.8701389 0.8165534 0.8044074
#> 
#> $mod1_unimodal$sim_properties$feature$gene_cor
#>    [1] -0.0040909414 -0.1878354751 -0.1934749560 -0.1221652169  0.0753980251
#>    [6]  0.1034648701 -0.0783372826  0.0335689647 -0.0067389482  0.1833412641
#>   [11] -0.1731617129  0.0857255419  0.0389468680  0.1216251710 -0.1112377988
#>   [16]  0.1056257276 -0.0885740884 -0.0094839910 -0.0288655109  0.1796024918
#>   [21] -0.3641422014  0.3086874975  0.0999788318 -0.1136567121  0.0430599161
#>   [26]  0.0684901898 -0.0749600151  0.0277243758  0.1265777094  0.0132317088
#>   [31]  0.1255081553 -0.2234196400 -0.0244391999 -0.0239680391 -0.0264462351
#>   [36] -0.1679999438 -0.0511414943 -0.1265833864 -0.0598724001 -0.1272097856
#>   [41] -0.1589567955  0.0238384053  0.1310425387  0.0450222128 -0.2223623198
#>   [46]  0.0878932454 -0.0907408503 -0.0982072661  0.1103962464  0.2785733998
#>   [51] -0.0870518363  0.0667952394  0.0985040151 -0.0088345277 -0.1564737392
#>   [56]  0.0856493993 -0.0650144347 -0.0388532948 -0.1818672600 -0.1046487867
#>   [61] -0.1085172391  0.0515848451 -0.0291540604  0.1112897112  0.0807301273
#>   [66]  0.0457010441  0.1302716044 -0.0110524859  0.0993108504  0.0362195491
#>   [71] -0.0994028431  0.1046502331 -0.0159171132  0.0784201017 -0.0227215583
#>   [76] -0.1464864291 -0.0970065998  0.0659193691  0.1962689060  0.0096514822
#>   [81] -0.0407748872 -0.1691246809 -0.0485113195 -0.0368323221 -0.0825607495
#>   [86] -0.0034423829 -0.1330239015  0.0199928654  0.0200759606 -0.0493384586
#>   [91] -0.1237665362 -0.0655439089  0.0065533274 -0.1606997585  0.0416334045
#>   [96] -0.1361716049  0.2571475891 -0.0779299393 -0.0313316265 -0.0376169700
#>  [101] -0.2142293842  0.0273000009 -0.0252303622 -0.0228675718 -0.0829950010
#>  [106] -0.0552029867 -0.1114216337  0.0324508200  0.1501785123 -0.0389847671
#>  [111] -0.0055318867 -0.0403492549  0.2178453459  0.0786472617 -0.0567689468
#>  [116] -0.0813377464 -0.0835865545  0.1417496216 -0.1058344113  0.1107108370
#>  [121] -0.2829864119 -0.1998437519  0.0784014493  0.1054006489  0.0647771378
#>  [126] -0.1121420381  0.0742700290 -0.0739073909 -0.0693116343  0.0107166542
#>  [131]  0.1318669704 -0.0379482408 -0.0722925523 -0.2553286500 -0.1469050434
#>  [136]  0.0347073303 -0.0306704856 -0.0306114649 -0.1035287428 -0.0809791887
#>  [141] -0.1540961449  0.2209625731 -0.1400025517  0.0787816423 -0.0906860036
#>  [146]  0.1150122717 -0.0920164522 -0.1974406609  0.0619553942  0.1629365406
#>  [151] -0.0128949533 -0.0993608394 -0.0605745598 -0.0808582747 -0.1013898998
#>  [156]  0.0785211268  0.1500152586  0.0181035856 -0.0997815233 -0.2453569813
#>  [161] -0.1647282159 -0.0561917135 -0.1388907966 -0.0549209314  0.0511629010
#>  [166] -0.1008957854  0.1395900211 -0.0473448382  0.0564321057 -0.1948178570
#>  [171] -0.0144977836 -0.0142515117 -0.0034890737 -0.0832557725 -0.0865550626
#>  [176]  0.0344422988  0.0839214609 -0.0962602735  0.0544879490  0.0043008886
#>  [181]  0.2513084694 -0.1571273416 -0.1520303722 -0.0056363007 -0.1766045470
#>  [186] -0.1971170046  0.0585775157 -0.0123680579  0.0289382646 -0.0676997788
#>  [191] -0.0829685704  0.0913361642 -0.1677263016  0.0583460528 -0.1460134214
#>  [196] -0.0310856586 -0.0746703597 -0.1022509758 -0.1225990697  0.1488375201
#>  [201]  0.0052740285 -0.0341386212  0.0898286143 -0.0763505772 -0.1002586377
#>  [206] -0.2061156565  0.0855342996  0.1489050065  0.0312485396  0.2016515089
#>  [211] -0.1486783017 -0.1848724104 -0.1002480443 -0.0850933942 -0.1200507018
#>  [216] -0.2850331440  0.0806932749 -0.2220850854  0.0204280193  0.2136469345
#>  [221] -0.1123283773  0.1252394707 -0.0156727690 -0.1009782270 -0.0556109145
#>  [226] -0.1308796902  0.1191671898 -0.0833079204  0.0754655133 -0.0483594969
#>  [231]  0.3452444182  0.0278048100  0.0196231432 -0.0071351986 -0.0901817823
#>  [236] -0.1130624543  0.1785738894 -0.1424971587  0.0119116490 -0.0605445604
#>  [241]  0.1388965740  0.0058211439 -0.1906307267 -0.1557778425  0.0950278912
#>  [246] -0.0320476440  0.0235857366  0.0019621167  0.1166777124  0.0783923647
#>  [251]  0.1052237476  0.0738255959 -0.1124752787 -0.0135879275  0.0043464746
#>  [256] -0.1420935714 -0.1266574750 -0.0500319157 -0.1137399801  0.0837622922
#>  [261]  0.1175446625  0.0164485350  0.0670318068 -0.0935187004 -0.0030530430
#>  [266] -0.0315397555 -0.1062266001  0.1355449754 -0.0817281073 -0.0750798024
#>  [271] -0.1121844657 -0.1429173022 -0.1718556731 -0.1671506793  0.0019230996
#>  [276]  0.0262264068 -0.0934476437 -0.0767544890  0.0786057095 -0.0461546771
#>  [281] -0.2216849264 -0.1156165533 -0.0121561839 -0.1589623796 -0.2346428311
#>  [286]  0.1798155724 -0.0117489891  0.0893192704 -0.1284030536  0.0901150055
#>  [291] -0.1271697993 -0.2771324052 -0.1323272871 -0.0758473350  0.1810586276
#>  [296] -0.0482149160  0.1278459136  0.2355860962 -0.0187973525  0.0980549142
#>  [301] -0.0640526460 -0.0939521159 -0.0922413463 -0.0790069858 -0.0184698507
#>  [306]  0.0174869335  0.0025504655  0.0546238285 -0.0330621706  0.0565723208
#>  [311]  0.0328403326  0.0068404064 -0.0858642770 -0.0503014455 -0.0668774835
#>  [316]  0.0799994053 -0.0271686611  0.0520656090 -0.0847635913  0.2296118197
#>  [321] -0.0434210299 -0.1040868758  0.0621583068  0.0221727714  0.0774868045
#>  [326] -0.1156297056  0.0865558816 -0.1088298293 -0.0530629098  0.0223644792
#>  [331] -0.1266657352 -0.0943381214 -0.0521515865 -0.0764922750  0.0908481658
#>  [336] -0.0570774729  0.0566880336 -0.0834432995  0.0576440611 -0.1071196395
#>  [341] -0.0798190580  0.2280783519 -0.0582290606 -0.2183362838 -0.0131247027
#>  [346]  0.0260104607  0.1675606504 -0.1826809372  0.1020431151 -0.1269638207
#>  [351] -0.1383044556  0.0687264328 -0.1300784381  0.1159447149 -0.2275695320
#>  [356] -0.1549676579  0.0446045761 -0.0749997343  0.0830039137 -0.1569463368
#>  [361]  0.0664719705 -0.0116117490  0.1937522116 -0.1328395876  0.0498996548
#>  [366]  0.0065974896 -0.0087294384 -0.0266421856  0.0973603285  0.0259632249
#>  [371] -0.1318064538 -0.0696268796 -0.0388983205  0.0521513652  0.0602498592
#>  [376]  0.1208931007  0.0374480582 -0.1314598746  0.3049959870  0.1630930076
#>  [381] -0.1260809006  0.0196232592 -0.1956219642 -0.1401436650 -0.0085281314
#>  [386]  0.0802851730 -0.0120215429  0.0509436450  0.0364799695  0.0358555945
#>  [391]  0.1857301403 -0.0329641380 -0.0008025616 -0.0396558111 -0.1702162803
#>  [396]  0.0329196176 -0.0189510271 -0.0400494394  0.0493450965  0.0307218201
#>  [401]  0.1087901528  0.1476183456  0.0449394614  0.0562197178  0.1148924861
#>  [406] -0.1384337742 -0.0672643212  0.0408794665  0.0363293649 -0.2235574356
#>  [411] -0.1063206903  0.0702845114 -0.1105530400 -0.1870010044  0.1572477594
#>  [416] -0.1322553654 -0.1134583288  0.1045087521 -0.0677958276  0.0666179655
#>  [421]  0.0704751479 -0.0542823593 -0.0600097027  0.0006718983  0.0281685060
#>  [426] -0.0060847733 -0.2131180405  0.0149059929 -0.2150355583 -0.0676112224
#>  [431] -0.0317395219 -0.0589753890  0.0711807402  0.0229751932  0.0812827650
#>  [436] -0.1390062942 -0.0530326640  0.1096040627  0.0474190569  0.1486039022
#>  [441]  0.0211980266  0.0330557268 -0.0970581490  0.0110717452 -0.1501475200
#>  [446]  0.0710168532 -0.2996296228 -0.0645180414  0.0660400956  0.1098672358
#>  [451]  0.0682899789 -0.0711990079 -0.0199748330  0.0451034405 -0.1680560859
#>  [456] -0.1343941609 -0.0624124719 -0.0789329351 -0.0690736501 -0.0018508249
#>  [461]  0.0389642052 -0.1299103777 -0.1831028027  0.1164321760  0.0660648822
#>  [466] -0.0480630813  0.0558282281 -0.0875774558 -0.1683823869  0.0632841539
#>  [471]  0.0527370765 -0.0025022962 -0.1001252075  0.0442114849 -0.0603270806
#>  [476]  0.1163354921  0.0780832882 -0.0921356422 -0.0207288610  0.0417404330
#>  [481] -0.0965567657  0.0493982373  0.0066949595 -0.0671305803  0.1195050679
#>  [486] -0.0322067860 -0.2313403092 -0.0670315873  0.0282360697 -0.1317090044
#>  [491]  0.0732034727  0.0648654292 -0.0446470249 -0.0298356158 -0.0552405529
#>  [496] -0.0493751725 -0.0646457323  0.0284569396 -0.0960124014  0.0799124176
#>  [501]  0.0673200005  0.0430002753 -0.0041651401 -0.0344918852 -0.1264074235
#>  [506] -0.0315384399  0.0539672474 -0.0680971270 -0.0106630966 -0.2219522012
#>  [511] -0.0182588933 -0.1094204715  0.2053703019  0.2587816327 -0.1206436093
#>  [516]  0.0054184191  0.2287439958  0.2123869828  0.1538818356 -0.1364710105
#>  [521]  0.0751180947 -0.1480591765 -0.1582730688  0.0321594287  0.0771925167
#>  [526] -0.1810196331 -0.1850233165 -0.1700871768 -0.1086655694 -0.1871417905
#>  [531] -0.1122285566  0.0831766844 -0.0668639984  0.1176367294 -0.1840329942
#>  [536] -0.0909556879  0.1125231993 -0.0625773305  0.2809471320  0.0759790177
#>  [541]  0.0979250680 -0.0202983994  0.0123660364  0.0332537060  0.0791035041
#>  [546] -0.1288075177  0.1277051133 -0.0602113767 -0.1053139245  0.1075679801
#>  [551]  0.0720393904 -0.0709986400 -0.1059748598 -0.0339229337  0.0361283126
#>  [556] -0.0737571568  0.1295717566  0.1243245963  0.0454516244  0.0040245037
#>  [561] -0.0351435952  0.0094374682 -0.1686257822  0.0293303646  0.0832594990
#>  [566] -0.0509893928 -0.1391036328 -0.1472225350 -0.0220371515 -0.0791697108
#>  [571] -0.0626896450 -0.1107655563  0.0048148547  0.0719811426  0.0064682419
#>  [576] -0.0306326193 -0.0660113887 -0.0506244235 -0.0169927570  0.1550487277
#>  [581]  0.1137037724  0.1139002573  0.2803174643 -0.0553065311  0.0700145160
#>  [586]  0.0974307024 -0.0179173142  0.0409635035 -0.0418798107  0.0378709566
#>  [591]  0.1055523398 -0.2018553732 -0.0474941310  0.1088840282  0.0679149458
#>  [596]  0.0794594796 -0.1133217772 -0.0046076769  0.1146186546  0.0411278097
#>  [601] -0.0168369797 -0.0477374795 -0.0919797967  0.0080265246  0.0264702974
#>  [606]  0.0819788362 -0.1817062605 -0.0706799245 -0.1419807562  0.1078590053
#>  [611] -0.1503949150 -0.0621146205 -0.0087067028  0.0227092649 -0.0575624537
#>  [616]  0.1105814247  0.1393908146  0.1453581474 -0.0214837212  0.0106650576
#>  [621] -0.0761725230 -0.0534246839 -0.0422653052  0.0022307149 -0.1896198144
#>  [626] -0.0319916262  0.1218376641  0.1200854442  0.1000677339 -0.0048229246
#>  [631] -0.0048723189  0.0522954676  0.0141548261 -0.1315593614  0.0784502121
#>  [636]  0.0416046486  0.0100954058  0.0594865730  0.1634238999 -0.0607909727
#>  [641] -0.1496847172  0.0349969504  0.1472346349 -0.1489503904 -0.1134718652
#>  [646] -0.1453063837  0.0869484745  0.0450051323 -0.1032054046  0.0100073253
#>  [651] -0.3385671684 -0.1588762839 -0.2582670905  0.2642093769  0.0159707479
#>  [656] -0.0364672870  0.0139585601  0.0101690939 -0.1260561242  0.0522824077
#>  [661]  0.0001944980  0.0866899204 -0.0612170799  0.0863240519 -0.0089380610
#>  [666] -0.1291426135  0.0896522169  0.0644124946 -0.0870837603 -0.1725545052
#>  [671] -0.1312682518 -0.0021006909 -0.0179546936  0.0700012088 -0.0496989039
#>  [676]  0.2360609566 -0.1613280553 -0.0494170857 -0.0420907507 -0.1065486286
#>  [681] -0.0471408174  0.0487555074  0.0418469359  0.1253444523 -0.0181682121
#>  [686]  0.2433922583  0.0946595073 -0.0436788493  0.1972826019  0.0854384042
#>  [691]  0.2051232361 -0.0848313968 -0.0585579880  0.1436873934  0.0717903954
#>  [696] -0.0945312798 -0.1171216435 -0.0943890120  0.1192320580 -0.1492966094
#>  [701] -0.1824773564 -0.1192963319 -0.0315301232 -0.1087523537  0.1216318644
#>  [706] -0.1249633314 -0.0542510135  0.1729944381 -0.1557622700 -0.0268069795
#>  [711]  0.0348807406 -0.2379499534 -0.0857400226  0.1368758514 -0.0227957998
#>  [716] -0.1736800905 -0.0894779786 -0.0619191143 -0.1120220392  0.1808871716
#>  [721] -0.1601603305 -0.0378562970  0.0640188065  0.0636573720  0.0652558829
#>  [726]  0.0061265147 -0.0311308151  0.0464524906 -0.3147624629  0.2723882181
#>  [731]  0.1040796281 -0.1915345143  0.0265734266 -0.2375533947  0.0635476261
#>  [736]  0.1774490038  0.0239391976  0.0026396766  0.0896554990 -0.0350326822
#>  [741]  0.0262542492 -0.1405970949  0.0132668684 -0.0371862922 -0.3189075790
#>  [746] -0.0125870689  0.0044027215 -0.1372410457 -0.0975308761  0.1049763145
#>  [751] -0.0142819266 -0.0808449751 -0.1227288238  0.0169395073  0.1010104170
#>  [756]  0.1659649640 -0.1204633973  0.0251328177 -0.0186055091  0.0003880922
#>  [761] -0.1955907833  0.0161371431  0.0521984275  0.1126557341  0.1570669194
#>  [766] -0.0617195710 -0.0330853797  0.0297592003  0.0012281123 -0.1238334227
#>  [771] -0.0009321597 -0.0128609667  0.0421457459  0.0122635695 -0.0363878464
#>  [776] -0.0673850795  0.0111828251  0.0710974083 -0.0456655050  0.1215662192
#>  [781] -0.1179015928  0.0176874129 -0.1336631383 -0.0934792637 -0.2153419166
#>  [786]  0.0940651692  0.0338804108  0.0563303833 -0.1684019257 -0.1971504164
#>  [791] -0.0479688837  0.0065988128 -0.0178361156 -0.0293430444  0.3404880654
#>  [796]  0.0905696732 -0.0218713278  0.0692865614 -0.0809450311  0.0505454767
#>  [801]  0.0190788306 -0.0957302694  0.1098656865 -0.0381485124  0.0052971279
#>  [806]  0.0630654645 -0.1132385978 -0.0133026793  0.0130091386  0.1063198219
#>  [811]  0.0510082684 -0.0001731682  0.0462380731 -0.0199414301  0.0604046517
#>  [816] -0.1924214002 -0.1065766543  0.0629213642  0.0295100461  0.1557593321
#>  [821] -0.0745208806 -0.0974162118  0.0811712885 -0.2730187611 -0.1358135985
#>  [826]  0.0959023592 -0.1869046068  0.0849111955  0.1325263455  0.0818977541
#>  [831]  0.0375153475 -0.0236590799 -0.0992176749  0.0893744103 -0.0298261209
#>  [836] -0.0537813786 -0.0345825129 -0.0236397747  0.0500109152  0.0528203778
#>  [841]  0.0920598290  0.0344309228  0.1681822598 -0.0597096231 -0.0614590482
#>  [846] -0.0948807273 -0.0892381204 -0.0394176993 -0.0357467129 -0.1977315566
#>  [851]  0.2152617217  0.0108089252 -0.0703646953  0.1226578515 -0.0374547249
#>  [856]  0.1355351181  0.0880389731 -0.0226860344 -0.0015310043  0.1033365035
#>  [861] -0.0657812322 -0.1759101548 -0.0187068426  0.2557996120 -0.1367998319
#>  [866] -0.1713183010  0.0740090059 -0.0290376064 -0.1993592917 -0.0301891731
#>  [871]  0.0748087668 -0.1266621888  0.0521863693  0.0599445623 -0.1818882282
#>  [876]  0.0384450510 -0.0488911155 -0.0632250812  0.0646197316 -0.0168868725
#>  [881]  0.0374831636  0.0738438438 -0.0264705256 -0.0723024178 -0.0278046779
#>  [886] -0.0500468298  0.0861852016 -0.0958460006  0.1052737708 -0.0598561138
#>  [891]  0.0035918309  0.1594471140  0.1429317296 -0.1355033055 -0.2393911327
#>  [896] -0.0380531312 -0.1967890252 -0.0136642367 -0.1125385613 -0.1489198052
#>  [901] -0.0438584461 -0.0436844533  0.0077236572  0.0644755381  0.0316057279
#>  [906] -0.0504088809  0.1511711507  0.0019252553 -0.0218370490  0.1450196279
#>  [911] -0.1145218770 -0.0132810672 -0.0676650121  0.1231860871  0.0638802640
#>  [916]  0.0158448761 -0.0568682584  0.1552276760 -0.1024749990 -0.0486448919
#>  [921] -0.1264000156 -0.0695616578  0.0299852436 -0.0227639817 -0.0114481223
#>  [926] -0.0640271126 -0.2190296905  0.0871660663  0.0326510548 -0.1628474010
#>  [931] -0.1295406155  0.0035926080 -0.2935134666  0.0638744826 -0.0109410667
#>  [936]  0.1546287097  0.1467517201 -0.1101990224  0.1310965688 -0.0177887694
#>  [941]  0.1381688422 -0.0735617436 -0.1960723856  0.0518424182 -0.0883025582
#>  [946]  0.0156740206  0.0775180805 -0.1428793997 -0.0983096849 -0.0504654073
#>  [951] -0.0668499170 -0.1515598614  0.0615534176  0.0226179669 -0.1764748757
#>  [956] -0.0312217031  0.0212035294 -0.0533666152  0.0259468898  0.1900558895
#>  [961] -0.0085409047  0.0251216664 -0.0314100140  0.0515764545 -0.0225897584
#>  [966] -0.3124860391 -0.1563696343 -0.0280119332 -0.1896243561  0.2795930366
#>  [971]  0.0994787824 -0.0506847467 -0.0073551042  0.0451107971  0.0932366393
#>  [976] -0.0186526806 -0.0676181411  0.0494448364  0.0019115147 -0.2026638740
#>  [981]  0.0237869798 -0.1205567042 -0.0122732966 -0.2124965320  0.0818317357
#>  [986]  0.0469652493 -0.0133688757 -0.0409034462  0.1638313489 -0.2299586991
#>  [991] -0.0345598824 -0.0461364725 -0.0053816602 -0.0743086936 -0.1402268093
#>  [996]  0.0598755667 -0.1302342860 -0.0627930198  0.0226433548 -0.0071681902
#> [1001]  0.2087281675 -0.0669033375 -0.0794156493  0.1023816812  0.0595898290
#> [1006] -0.0347421616 -0.0459743814  0.2850659028 -0.0436891889 -0.0669084808
#> [1011]  0.0984996495 -0.1448263617  0.2419098577  0.1500449042  0.0617721891
#> [1016] -0.0855649092 -0.0680062955  0.0510937845  0.0322693272  0.0380450515
#> [1021] -0.1463391267 -0.0488477729  0.1312038570 -0.1311674854  0.1669971986
#> [1026]  0.0147298410 -0.1680943323  0.0286550723  0.0417142641 -0.0105545817
#> [1031]  0.1046127699 -0.0619440881  0.0509056323 -0.1160683140  0.1209399052
#> [1036] -0.1651517560 -0.1766974635  0.0343406872  0.0954270126  0.1662042048
#> [1041]  0.1274214415 -0.0609415812 -0.1062139649  0.0606756932  0.0830883263
#> [1046]  0.2951332207  0.1372279154  0.1564570156 -0.0765143226 -0.1283852819
#> [1051]  0.0902776953  0.0677353004 -0.0727100800  0.0317533751 -0.0718525690
#> [1056]  0.0893606296 -0.0577880071  0.0307770162 -0.2105589344 -0.0622012819
#> [1061]  0.0911451594 -0.2093219253 -0.0943354886 -0.2366342710  0.0058762154
#> [1066] -0.1694547972 -0.0283544245  0.2015542869  0.1608602151 -0.0143998492
#> [1071]  0.0595423257  0.0211979088 -0.0988140259 -0.0129692543  0.1520939879
#> [1076] -0.1113451073 -0.0736376418 -0.1104556771  0.0149138399 -0.0835361296
#> [1081]  0.0679182648  0.2365955169  0.1546326417 -0.1264723085  0.0885091858
#> [1086]  0.1737005086  0.1041014764  0.0939128419  0.0735581409 -0.0763645691
#> [1091] -0.0531039194 -0.2056200323 -0.1246711286 -0.0779261006  0.0543975166
#> [1096] -0.0261311738 -0.1029969897 -0.2056553036 -0.0865938473 -0.0355673983
#> [1101]  0.0328076283 -0.2243514612 -0.2408038625  0.0772466916 -0.0305807000
#> [1106] -0.1791700691 -0.1271627411  0.1488690761 -0.0362330826 -0.0979833327
#> [1111]  0.0389759645 -0.0850888680 -0.1054742753 -0.1374766513 -0.2238134824
#> [1116] -0.0394655545 -0.1569564527  0.0269351216  0.1518412993  0.0687702898
#> [1121]  0.0356732810  0.0961929284 -0.2317497915 -0.0366525210 -0.0437645790
#> [1126] -0.0478668682 -0.1106061759 -0.1453564947  0.0886470027  0.1865823829
#> [1131] -0.1935017270  0.2207438930 -0.1144823556  0.0164480645 -0.0253227585
#> [1136] -0.0529042019 -0.1574127188  0.0411790570 -0.0435095293 -0.0909692964
#> [1141] -0.0889010300  0.1231038694  0.1212716018 -0.1334246186 -0.0992725219
#> [1146] -0.0163180844 -0.2005675362  0.0377656031  0.1300365333 -0.0960087946
#> [1151] -0.1009838799 -0.0151203799  0.1565294808  0.0292504668  0.0214206888
#> [1156] -0.0526490554  0.0576300326 -0.0555132249  0.1332678910  0.1075279182
#> [1161] -0.0544595014 -0.1812081725  0.1125798528  0.0477315729 -0.2021744084
#> [1166] -0.0280881912 -0.0692405530 -0.0497183360  0.0033826432 -0.1482388020
#> [1171]  0.0369553397 -0.0012282759  0.0259930409  0.0282564847 -0.2440322646
#> [1176]  0.0943853593 -0.0348587058 -0.0216887157  0.0827339617  0.1004145158
#> [1181]  0.1995888673 -0.0109772028 -0.1706988417  0.0896815642 -0.0520539894
#> [1186] -0.1523231518 -0.0721465223 -0.0314958936  0.2441193111  0.0746260000
#> [1191] -0.1895461724 -0.0252645894  0.0296720256  0.1610530801  0.2230399180
#> [1196]  0.1526774021 -0.1205801230 -0.0578101942 -0.0230469457 -0.1384726961
#> [1201] -0.0417796556  0.0660104269 -0.0066391861 -0.1068311588 -0.0056093347
#> [1206] -0.0778299621  0.1077553570 -0.0162857650 -0.0445908208  0.1246593072
#> [1211]  0.0964251002 -0.0320916071  0.0307321295 -0.0494402174 -0.1261617991
#> [1216] -0.1547371625 -0.0500711018 -0.0756891054 -0.0497756442 -0.0701981376
#> [1221] -0.1709861534 -0.1650536475  0.0300184304  0.0579260656 -0.0914701319
#> [1226]  0.2089868751 -0.0134827466 -0.0906890315 -0.0599636412 -0.0654604561
#> [1231] -0.2066983858  0.1510268658 -0.0834091847  0.1743513116 -0.0738269859
#> [1236] -0.0924539158 -0.0231999864 -0.1369582403  0.1608313485 -0.0138445304
#> [1241] -0.1869040734  0.1097572066 -0.0043836956 -0.0312594273 -0.1588211992
#> [1246] -0.0402185707 -0.0470643659 -0.0624642013 -0.1596936893 -0.1031209140
#> [1251] -0.1764452907 -0.0933607808  0.1143329534 -0.0934873527  0.0539091978
#> [1256]  0.0903348108 -0.0769354341  0.0572201817  0.0366936551 -0.1127974821
#> [1261]  0.1920866584 -0.0458394700 -0.0047662650 -0.0898899078 -0.0680034447
#> [1266] -0.0968301643  0.0053030935 -0.1088746791  0.1874889934 -0.1141405598
#> [1271] -0.0489330405  0.0076140087  0.0021217217  0.1606980096  0.0414786246
#> [1276]  0.1817099525  0.0001651565  0.0330344619 -0.0865539719 -0.1751792930
#> [1281]  0.0975850446 -0.0377486028  0.0595160239 -0.0719339837  0.0314528502
#> [1286]  0.0820131349  0.0516503676  0.1930143824 -0.0589980331 -0.0346987827
#> [1291] -0.0649259386  0.0380181630 -0.1069394618 -0.2253721065  0.0810207070
#> [1296] -0.1617089137 -0.3409038414  0.1586796330  0.0663823708 -0.1676234374
#> [1301]  0.1505489982  0.0260559693  0.1415807723  0.0879855626  0.0315066404
#> [1306] -0.1974126860  0.0598005978 -0.1562585147  0.0566109205 -0.1678085578
#> [1311] -0.0547643162 -0.0398486586  0.0041667553 -0.0263468865  0.1163813547
#> [1316]  0.1074880457 -0.0364020725  0.0387962756  0.0006046928 -0.0359962485
#> [1321]  0.1635450436  0.0937981831  0.1118268591  0.0419178869 -0.0625188604
#> [1326]  0.0903967197 -0.0478055281 -0.2321286270 -0.0609176407  0.0463039841
#> [1331]  0.1768937590 -0.1298406375 -0.0283866143  0.1656519469 -0.1376670644
#> [1336]  0.0277808024  0.0922248947 -0.1187189664  0.0007939101 -0.0485231429
#> [1341]  0.0423473590  0.3326685691 -0.0503565833 -0.0549442014  0.0124365559
#> [1346]  0.0411093631  0.0379191013  0.0298977598 -0.0839465204  0.0721447243
#> [1351]  0.0087923886  0.2107750448  0.1391129113  0.1000097827 -0.0531331330
#> [1356] -0.1013442195  0.1286287674 -0.0991835896 -0.1969505529 -0.1736728146
#> [1361]  0.0823732412 -0.0263857031 -0.0367330210 -0.1637599839  0.0264074937
#> [1366] -0.0963711910 -0.2063378401 -0.0043154033 -0.0232612529 -0.1181178617
#> [1371]  0.0424097442 -0.0022692445 -0.0116182559 -0.0317455787 -0.0203618330
#> [1376] -0.0097553559 -0.1831287694 -0.1258721088  0.0578960037  0.0723515573
#> [1381]  0.1967297244  0.0130127768  0.0728722778  0.0919964760  0.0192083785
#> [1386]  0.1089293914 -0.0259529779 -0.0136401299  0.1300067255  0.0668087543
#> [1391] -0.1055296053  0.1840716131  0.0354966012 -0.1433231449  0.0752100156
#> [1396] -0.1560141199 -0.2060680861 -0.1252385144 -0.1086967427 -0.0485427847
#> [1401] -0.0652204388  0.1107162747 -0.1377126979  0.0209598525  0.1318730751
#> [1406]  0.0549423437 -0.1159513204 -0.1551330978 -0.0721936998  0.0910963457
#> [1411] -0.1613520506 -0.1284469762  0.0054218017 -0.1128553892  0.0838044602
#> [1416] -0.2083982099 -0.0531697725  0.1799990966  0.1676656415  0.0451506731
#> [1421]  0.0096411165 -0.1780465808  0.0011302728  0.0510530202  0.0291624392
#> [1426]  0.0622310846  0.0642113726 -0.1118245602 -0.1398095577  0.0665007513
#> [1431] -0.0739026895 -0.2514461603  0.0981830665  0.1187797443 -0.1603152569
#> [1436] -0.1704077390  0.0660525946  0.0365751553 -0.2257042230 -0.2008967500
#> [1441] -0.1033259685 -0.2284892505 -0.1079197652 -0.0243270147 -0.0636945416
#> [1446]  0.1886126099  0.0287861660  0.0817742953 -0.0301537562 -0.0872619508
#> [1451]  0.0311029050  0.0981039571 -0.0175920609 -0.0472902312 -0.2081098362
#> [1456]  0.1102152162  0.0646486431 -0.0863219008 -0.1161865100 -0.2145843743
#> [1461]  0.0458764180  0.0877194444 -0.0435990474  0.0949853160 -0.3871039964
#> [1466] -0.1350942107 -0.0771785704 -0.1127021651 -0.0182102072  0.1032530954
#> [1471]  0.1485099513  0.1754138649  0.0468986917  0.2262302896  0.1184396944
#> [1476]  0.0316174763 -0.0613687861 -0.0824279418 -0.0282065064  0.0410960785
#> [1481] -0.1334753423  0.0748405542  0.0189530794  0.1206119808 -0.1176513270
#> [1486] -0.0090596786 -0.1451770096 -0.0785579124 -0.0622688935  0.0612317811
#> [1491] -0.2714767777 -0.0113235440 -0.2246522716  0.0320826982 -0.0812669012
#> [1496] -0.1738958015 -0.2009093685 -0.0055273382  0.1562436102 -0.2727121214
#> [1501] -0.0008174108  0.0682422877 -0.0683387138  0.1955814529  0.1635877965
#> [1506]  0.0963445008  0.2683880421  0.0455955979 -0.1346437128  0.0657663657
#> [1511] -0.1319558419  0.0510242116 -0.2329926673 -0.0898828577  0.0803672894
#> [1516]  0.1354578313 -0.0246422814 -0.1360068166  0.0373192963  0.0984640800
#> [1521]  0.0034617060 -0.0349812251  0.0158729087  0.1053886696  0.0243044248
#> [1526] -0.0251645567  0.0422056601 -0.1158280566 -0.2111326557  0.1322666352
#> [1531] -0.0273009562  0.0349829645  0.0700705336 -0.0535163672  0.1767752249
#> [1536]  0.0991714067 -0.2193613953 -0.0388475563 -0.0537615653  0.0264616223
#> [1541] -0.1335648940  0.2079433523 -0.0293728886  0.0384983650 -0.0701992046
#> [1546]  0.0968553953  0.0889328696 -0.1815725049 -0.0686156952  0.0211197566
#> [1551] -0.0054917423 -0.0895319942 -0.0071332957 -0.0009070917  0.2028334866
#> [1556]  0.0197322053  0.0313089002 -0.0793583545  0.0976192019 -0.0627421541
#> [1561] -0.0661649989  0.0470158139 -0.1065759671  0.0095730850 -0.0221046450
#> [1566] -0.0115198604 -0.0636390055 -0.0042125407 -0.0764254078  0.0847354424
#> [1571]  0.0300173445 -0.0771737803  0.1018485399  0.0024814824 -0.0722982653
#> [1576]  0.0311610171 -0.0425231456 -0.0355053242  0.0180366783  0.1671273274
#> [1581] -0.0255177896 -0.3187623707 -0.0307098302  0.1514696359 -0.1384457738
#> [1586]  0.0627848078  0.0792533703  0.0161956025  0.1317199677 -0.0415144445
#> [1591]  0.1903094386 -0.0702110087 -0.1340761104 -0.0199561824  0.1400879795
#> [1596] -0.0752508734 -0.0125901862 -0.1738417673 -0.1061252410 -0.1678471558
#> [1601] -0.1414097947 -0.1049982203 -0.0980065205 -0.1715275445  0.0316816303
#> [1606]  0.1803600395 -0.0117438349  0.3006894315 -0.2681234112 -0.0713276029
#> [1611]  0.1143556617 -0.0994799913  0.0201836119  0.1117400816  0.1092223269
#> [1616] -0.1221548113 -0.0452863865  0.1272868203  0.0487027636 -0.0160687656
#> [1621] -0.0437914607  0.1077943713  0.0822502106  0.1911128011  0.1121083489
#> [1626]  0.1685019554 -0.2310092306  0.1154582258  0.0665939420  0.0804121020
#> [1631] -0.0627312486  0.0101421845 -0.1327012258  0.0565880651 -0.0922877019
#> [1636]  0.1679342604  0.1777185149 -0.1659177069  0.0368335549 -0.0645142666
#> [1641]  0.0912088735  0.0488495052  0.0512864827  0.0563464321 -0.1246473461
#> [1646] -0.0757361011  0.0522281861  0.0823809983 -0.1730842540  0.0100664004
#> [1651] -0.0503006325 -0.0244207110  0.0515040746  0.1543634856 -0.1946252873
#> [1656] -0.0826016376 -0.1040498860 -0.1829797331  0.0294886639  0.0796203905
#> [1661] -0.3110377584  0.0222316513  0.0864632363 -0.0888512494 -0.0318808943
#> [1666] -0.0143808257  0.2146865759 -0.1246671801 -0.0926390186  0.0608792275
#> [1671] -0.1074001234 -0.0447216222 -0.0108123827  0.0866154791  0.1325227829
#> [1676] -0.0164568899 -0.1613753588  0.0944293322 -0.0826594141 -0.0675129436
#> [1681] -0.0587773202 -0.0355380302  0.1776814381 -0.0011074776 -0.1607692086
#> [1686] -0.0683033409  0.0699295916  0.0108879877 -0.0685042837  0.1184155905
#> [1691]  0.0869324973 -0.0634276996  0.0543500542  0.0495877951  0.1811030780
#> [1696] -0.0152647207  0.0441619908  0.0780109780  0.0198729382  0.0956440394
#> [1701]  0.0215604270 -0.0353113757 -0.2021668340  0.1736063554  0.1148283528
#> [1706] -0.0633439385 -0.0385215911  0.0458529895  0.2285874481  0.0082950899
#> [1711] -0.0500120657 -0.1630509207  0.0129266443  0.2156597823 -0.0075134772
#> [1716] -0.0077481775  0.0239838665 -0.0876313837 -0.0392343799 -0.0513107870
#> [1721]  0.1377209046 -0.0851518515  0.0272784405 -0.0322818937 -0.0098109127
#> [1726]  0.0503438212 -0.1378495464 -0.1986697457 -0.0133030527  0.1544553236
#> [1731] -0.0185370464  0.1230559242  0.0396265889 -0.0903168798  0.0086516756
#> [1736]  0.2174065360 -0.0669576919 -0.1055197964 -0.0419276630 -0.0337424707
#> [1741] -0.0606296700  0.1559946136 -0.0181170883  0.0126191427 -0.2039419739
#> [1746]  0.0847057503 -0.1126374550  0.0112864164  0.0098541935 -0.0122193573
#> [1751] -0.0002152624  0.1025627022  0.1578179738  0.2199863268 -0.0209412220
#> [1756]  0.1692981354 -0.0643774798 -0.0899030046 -0.0961696072 -0.0255439096
#> [1761] -0.1960490880 -0.2426546145 -0.3965361015  0.0486992299  0.1491295031
#> [1766]  0.0855070047 -0.0414992694 -0.1738315148 -0.0115378228 -0.1003924810
#> 
#> $mod1_unimodal$sim_properties$feature$prop_outliers_feature
#> [1] 0
#> 
#> 
#> 
#> $mod1_unimodal$batch_effects
#> $mod1_unimodal$batch_effects$batch_silhouette
#> [1] -0.002710457
#> 
#> $mod1_unimodal$batch_effects$pcr_r2
#> [1] 0.009476487
#> 
#> $mod1_unimodal$batch_effects$cms
#> [1] 0
#> 
#> $mod1_unimodal$batch_effects$isi
#> [1] 1.938397
#> 
#> $mod1_unimodal$batch_effects$cross_batch_F1
#> [1] 0.5126418
#> 
#> 
#> $mod1_unimodal$biological_signal
#> $mod1_unimodal$biological_signal$deg_summary_table
#>                         SubCategory          Metric     Value
#> 1 scRNA-seq Differential Expression pval_dist_score 0.8719073
#> 2 scRNA-seq Differential Expression        accuracy 0.9500000
#> 3 scRNA-seq Differential Expression       precision 0.0000000
#> 4 scRNA-seq Differential Expression          recall 0.0000000
#> 5 scRNA-seq Differential Expression        f1_score 0.0000000
#> 6 scRNA-seq Differential Expression             auc 0.6465517
#>                                                         Direction
#> 1    Higher is better (1.0 = uniform null p-values, no inflation)
#> 2             Higher is better (1.0 = perfect DEG classification)
#> 3                 Higher is better (1.0 = zero false discoveries)
#> 4                 Higher is better (1.0 = 100% true DEG recovery)
#> 5 Higher is better (1.0 = harmonic balance of precision & recall)
#> 6         Higher is better (1.0 = perfect ranking discrimination)
#>                                            Framework
#> 1 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 2 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 3 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 4 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 5 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 6 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 
#> $mod1_unimodal$biological_signal$pairwise_results
#>             Pair                       SubCategory N_total N_true_DEG N_sim_DEG
#> 1 TypeA_vs_TypeB scRNA-seq Differential Expression      60          2         1
#>   TP FP FN TN pval_dist_score accuracy precision recall f1_score       auc
#> 1  0  1  2 57       0.8719073     0.95         0      0        0 0.6465517
#>   LogFC_Pearson_r LogFC_Spearman_rho Top_DEG_Jaccard
#> 1       0.4413862          0.3846624               1
#> 
#> $mod1_unimodal$biological_signal$pseudobulks_ref
#> $mod1_unimodal$biological_signal$pseudobulks_ref$TypeA
#>         Batch1 Batch2
#> Gene_01    164     87
#> Gene_02     65     67
#> Gene_03     75     81
#> Gene_04     72     95
#> Gene_05    112     99
#> Gene_06    110     72
#> Gene_07     97     81
#> Gene_08    104     76
#> Gene_09    156    115
#> Gene_10    105     68
#> Gene_11    102     81
#> Gene_12    123    108
#> Gene_13     45     30
#> Gene_14     97     55
#> Gene_15     72     75
#> Gene_16     67     97
#> Gene_17     75     63
#> Gene_18     62     32
#> Gene_19     56     44
#> Gene_20    112     63
#> Gene_21     72     56
#> Gene_22     75     51
#> Gene_23     94     71
#> Gene_24    132    104
#> Gene_25    142    113
#> Gene_26    103     64
#> Gene_27    110     55
#> Gene_28     43     49
#> Gene_29     89     72
#> Gene_30     96     79
#> Gene_31     98     82
#> Gene_32     73     91
#> Gene_33     92     80
#> Gene_34     90     48
#> Gene_35     83     63
#> Gene_36     41     42
#> Gene_37     58     55
#> Gene_38     93     75
#> Gene_39     34     35
#> Gene_40     87    102
#> Gene_41     76    104
#> Gene_42     46     60
#> Gene_43     95     63
#> Gene_44     63     40
#> Gene_45     59     42
#> Gene_46     63    118
#> Gene_47     62     48
#> Gene_48     86    113
#> Gene_49     58    113
#> Gene_50    124     97
#> Gene_51     96     62
#> Gene_52     51    117
#> Gene_53    127    137
#> Gene_54    124     87
#> Gene_55     58     64
#> Gene_56    120     52
#> Gene_57     81    140
#> Gene_58     64     93
#> Gene_59     45     32
#> Gene_60     85     62
#> 
#> $mod1_unimodal$biological_signal$pseudobulks_ref$TypeB
#>         Batch1 Batch2
#> Gene_01    137    159
#> Gene_02     93    133
#> Gene_03    175    115
#> Gene_04    119    106
#> Gene_05    188    174
#> Gene_06     81    109
#> Gene_07    129    140
#> Gene_08     94    107
#> Gene_09    147    185
#> Gene_10    184     94
#> Gene_11    179    139
#> Gene_12    197     97
#> Gene_13     87     99
#> Gene_14    131    160
#> Gene_15    111    137
#> Gene_16     73     96
#> Gene_17     72     84
#> Gene_18     40     36
#> Gene_19     46     40
#> Gene_20     91    156
#> Gene_21     46     50
#> Gene_22     41     41
#> Gene_23     67     69
#> Gene_24     85     90
#> Gene_25     81    111
#> Gene_26     84     81
#> Gene_27    100     53
#> Gene_28     62     38
#> Gene_29     86     94
#> Gene_30     45     66
#> Gene_31     57     98
#> Gene_32    113     87
#> Gene_33     79     91
#> Gene_34     43     51
#> Gene_35    103     99
#> Gene_36     48     58
#> Gene_37     83     63
#> Gene_38     56     81
#> Gene_39     51     39
#> Gene_40     90     87
#> Gene_41     58    106
#> Gene_42     90     50
#> Gene_43     91     94
#> Gene_44     64     50
#> Gene_45     69     48
#> Gene_46     74     51
#> Gene_47     56     65
#> Gene_48    100     88
#> Gene_49     69     90
#> Gene_50     98    119
#> Gene_51    108    114
#> Gene_52     54     52
#> Gene_53    103     95
#> Gene_54     91    104
#> Gene_55    108     86
#> Gene_56     94     82
#> Gene_57    118     66
#> Gene_58    113    133
#> Gene_59     23     19
#> Gene_60     70     56
#> 
#> 
#> $mod1_unimodal$biological_signal$pseudobulks_sim
#> $mod1_unimodal$biological_signal$pseudobulks_sim$TypeA
#>         Batch1 Batch2
#> Gene_01     82     96
#> Gene_02     67     77
#> Gene_03     81     67
#> Gene_04     62     78
#> Gene_05     53     96
#> Gene_06     86     61
#> Gene_07    108    132
#> Gene_08     68     86
#> Gene_09     61    135
#> Gene_10    122     74
#> Gene_11     70     52
#> Gene_12    109    106
#> Gene_13     46     30
#> Gene_14     57     82
#> Gene_15     57     84
#> Gene_16     89     87
#> Gene_17     64     55
#> Gene_18     35     32
#> Gene_19     44     49
#> Gene_20    106     91
#> Gene_21     47     57
#> Gene_22     50     36
#> Gene_23     79     54
#> Gene_24    107    126
#> Gene_25    114     81
#> Gene_26     48     50
#> Gene_27     62     86
#> Gene_28     48     66
#> Gene_29     79    114
#> Gene_30     57    100
#> Gene_31     75     72
#> Gene_32     76     68
#> Gene_33    120     74
#> Gene_34     70     94
#> Gene_35    112     87
#> Gene_36     46     50
#> Gene_37     60     75
#> Gene_38     54     58
#> Gene_39     38     29
#> Gene_40     63     47
#> Gene_41     80     61
#> Gene_42     53     93
#> Gene_43     89     84
#> Gene_44     90     78
#> Gene_45     54     42
#> Gene_46     84     39
#> Gene_47     53     93
#> Gene_48     86    102
#> Gene_49     60     69
#> Gene_50     80     80
#> Gene_51     65     59
#> Gene_52     88     85
#> Gene_53     97    113
#> Gene_54     77     97
#> Gene_55     57     87
#> Gene_56     67     95
#> Gene_57     67     85
#> Gene_58     67     57
#> Gene_59     36     24
#> Gene_60    103     63
#> 
#> $mod1_unimodal$biological_signal$pseudobulks_sim$TypeB
#>         Batch1 Batch2
#> Gene_01    108    158
#> Gene_02    184    162
#> Gene_03    111    129
#> Gene_04    186    106
#> Gene_05     84    113
#> Gene_06    137    202
#> Gene_07    156    111
#> Gene_08    110    110
#> Gene_09    166    207
#> Gene_10     93    128
#> Gene_11    187    140
#> Gene_12    156    169
#> Gene_13     94     96
#> Gene_14     63     89
#> Gene_15    154    104
#> Gene_16     53     78
#> Gene_17     62     72
#> Gene_18     42     28
#> Gene_19     46     28
#> Gene_20     88    125
#> Gene_21     66     67
#> Gene_22     35     65
#> Gene_23     64     74
#> Gene_24    105     81
#> Gene_25    110     84
#> Gene_26     74     71
#> Gene_27     60     89
#> Gene_28     46     27
#> Gene_29    109     81
#> Gene_30     86     56
#> Gene_31     74     83
#> Gene_32     58    111
#> Gene_33    125    115
#> Gene_34     59     59
#> Gene_35    104     85
#> Gene_36     50     42
#> Gene_37     39     59
#> Gene_38    100     67
#> Gene_39     34     40
#> Gene_40     77     68
#> Gene_41     69     48
#> Gene_42     96     49
#> Gene_43     90     41
#> Gene_44     32     65
#> Gene_45     64     70
#> Gene_46     71     82
#> Gene_47     49     54
#> Gene_48    101    106
#> Gene_49     42     55
#> Gene_50    103    110
#> Gene_51     59     84
#> Gene_52     75     68
#> Gene_53    116     80
#> Gene_54     61    156
#> Gene_55     35    116
#> Gene_56     73     86
#> Gene_57     78    116
#> Gene_58    107     83
#> Gene_59     30     31
#> Gene_60     75    101
#> 
#> 
#> $mod1_unimodal$biological_signal$de_results_ref
#> $mod1_unimodal$biological_signal$de_results_ref$TypeA_vs_TypeB
#>         feature         logFC   logCPM       F_stat      PValue        FDR
#> Gene_01 Gene_01  1.238685e-01 14.70815 1.574821e-01 0.696937479 0.85799160
#> Gene_02 Gene_02  6.239397e-01 14.11652 3.907120e+00 0.066296970 0.28099865
#> Gene_03 Gene_03  7.316681e-01 14.42108 5.266377e+00 0.036196721 0.27517630
#> Gene_04 Gene_04  2.591277e-01 14.25533 6.691704e-01 0.425833651 0.68681783
#> Gene_05 Gene_05  6.309193e-01 14.77632 4.718310e+00 0.045862717 0.27517630
#> Gene_06 Gene_06 -6.371856e-02 14.17162 4.291334e-02 0.838606347 0.93178483
#> Gene_07 Gene_07  4.531216e-01 14.42722 2.387011e+00 0.142672602 0.40763601
#> Gene_08 Gene_08  2.531663e-02 14.20603 7.375600e-03 0.932668838 0.97218763
#> Gene_09 Gene_09  1.599954e-01 14.85542 2.928419e-01 0.596156983 0.77759606
#> Gene_10 Gene_10  5.487471e-01 14.42792 2.586612e+00 0.128101479 0.39057322
#> Gene_11 Gene_11  6.538953e-01 14.58179 4.828320e+00 0.043702412 0.27517630
#> Gene_12 Gene_12  1.945438e-01 14.65742 3.309318e-01 0.573431405 0.77759606
#> Gene_13 Gene_13  1.176655e+00 13.65571 1.655902e+01 0.000960716 0.02882148
#> Gene_14 Gene_14  8.166494e-01 14.40094 6.784185e+00 0.019614268 0.23537122
#> Gene_15 Gene_15  6.004514e-01 14.25517 3.848478e+00 0.068150555 0.28099865
#> Gene_16 Gene_16 -1.284917e-01 14.03806 1.476561e-01 0.706046162 0.85799160
#> Gene_17 Gene_17  3.504218e-02 13.84443 1.478805e-02 0.904785635 0.96941318
#> Gene_18 Gene_18 -4.258217e-01 13.07585 2.556344e+00 0.130191074 0.39057322
#> Gene_19 Gene_19 -3.574471e-01 13.20907 2.079732e+00 0.169321041 0.42368029
#> Gene_20 Gene_20  3.817123e-01 14.34034 1.257287e+00 0.279367850 0.55873570
#> Gene_21 Gene_21 -5.521256e-01 13.47104 4.306360e+00 0.055140368 0.28099865
#> Gene_22 Gene_22 -7.497055e-01 13.36652 8.146085e+00 0.011833590 0.23537122
#> Gene_23 Gene_23 -4.150511e-01 13.88186 2.093072e+00 0.168041553 0.42368029
#> Gene_24 Gene_24 -5.696138e-01 14.32594 3.784301e+00 0.070249662 0.28099865
#> Gene_25 Gene_25 -5.457576e-01 14.44601 3.205418e+00 0.093096359 0.32857538
#> Gene_26 Gene_26 -1.435495e-01 14.00958 2.341711e-01 0.635266118 0.81097802
#> Gene_27 Gene_27 -2.285939e-01 13.94197 4.306683e-01 0.521359349 0.74479907
#> Gene_28 Gene_28 -4.184076e-02 13.25524 2.101028e-02 0.886633502 0.96723655
#> Gene_29 Gene_29  2.066557e-02 14.05210 5.046172e-03 0.944284272 0.97218763
#> Gene_30 Gene_30 -7.942410e-01 13.82018 6.838140e+00 0.019211489 0.23537122
#> Gene_31 Gene_31 -3.533610e-01 14.03673 1.187289e+00 0.292655965 0.56643090
#> Gene_32 Gene_32  1.174641e-01 14.15167 1.383590e-01 0.714993001 0.85799160
#> Gene_33 Gene_33 -1.606162e-01 14.06180 2.998297e-01 0.591847526 0.77759606
#> Gene_34 Gene_34 -6.697999e-01 13.51131 5.192875e+00 0.037344978 0.27517630
#> Gene_35 Gene_35  3.299970e-01 14.07297 1.289504e+00 0.273526213 0.55873570
#> Gene_36 Gene_36  2.013795e-01 13.22965 5.948548e-01 0.452235707 0.68681783
#> Gene_37 Gene_37  2.171497e-01 13.66426 5.798783e-01 0.457878556 0.68681783
#> Gene_38 Gene_38 -4.311995e-01 13.90402 2.006781e+00 0.176533454 0.42368029
#> Gene_39 Gene_39  2.283944e-01 12.98891 8.440496e-01 0.372426167 0.68681783
#> Gene_40 Gene_40 -2.578351e-01 14.16743 7.123154e-01 0.411611858 0.68681783
#> Gene_41 Gene_41 -3.003584e-01 14.08788 7.024743e-01 0.414789659 0.68681783
#> Gene_42 Gene_42  2.297639e-01 13.59758 4.846403e-01 0.496714734 0.72689961
#> Gene_43 Gene_43  9.793853e-02 14.05362 1.106754e-01 0.743867872 0.87513867
#> Gene_44 Gene_44  1.529073e-02 13.41191 3.147207e-03 0.955984507 0.97218763
#> Gene_45 Gene_45  7.416194e-02 13.41949 7.239424e-02 0.791458225 0.91322103
#> Gene_46 Gene_46 -7.255692e-01 13.93750 3.810639e+00 0.069379263 0.28099865
#> Gene_47 Gene_47 -5.884604e-05 13.50544 5.181635e-08 0.999821302 0.99982130
#> Gene_48 Gene_48 -2.531384e-01 14.24921 6.374833e-01 0.436779547 0.68681783
#> Gene_49 Gene_49 -2.930161e-01 14.03696 6.157780e-01 0.444540258 0.68681783
#> Gene_50 Gene_50 -1.632461e-01 14.40896 3.027851e-01 0.590044619 0.77759606
#> Gene_51 Gene_51  3.623413e-01 14.19317 1.480713e+00 0.242006955 0.52467350
#> Gene_52 Gene_52 -8.621055e-01 13.79145 4.836415e+00 0.043548525 0.27517630
#> Gene_53 Gene_53 -5.744205e-01 14.50217 3.645005e+00 0.075074344 0.28152879
#> Gene_54 Gene_54 -2.446537e-01 14.29948 6.812593e-01 0.421771937 0.68681783
#> Gene_55 Gene_55  5.082222e-01 13.94117 2.810622e+00 0.113845696 0.37948565
#> Gene_56 Gene_56 -7.467621e-02 14.06516 5.081873e-02 0.824612789 0.93178483
#> Gene_57 Gene_57 -4.559174e-01 14.32490 1.462159e+00 0.244847631 0.52467350
#> Gene_58 Gene_58  4.744774e-01 14.29319 2.058091e+00 0.171420964 0.42368029
#> Gene_59 Gene_59 -1.006405e+00 12.60231 2.025826e+01 0.000398335 0.02390010
#> Gene_60 Gene_60 -3.590200e-01 13.74176 1.590569e+00 0.226040127 0.52163106
#>         is_da
#> Gene_01 FALSE
#> Gene_02 FALSE
#> Gene_03 FALSE
#> Gene_04 FALSE
#> Gene_05 FALSE
#> Gene_06 FALSE
#> Gene_07 FALSE
#> Gene_08 FALSE
#> Gene_09 FALSE
#> Gene_10 FALSE
#> Gene_11 FALSE
#> Gene_12 FALSE
#> Gene_13  TRUE
#> Gene_14 FALSE
#> Gene_15 FALSE
#> Gene_16 FALSE
#> Gene_17 FALSE
#> Gene_18 FALSE
#> Gene_19 FALSE
#> Gene_20 FALSE
#> Gene_21 FALSE
#> Gene_22 FALSE
#> Gene_23 FALSE
#> Gene_24 FALSE
#> Gene_25 FALSE
#> Gene_26 FALSE
#> Gene_27 FALSE
#> Gene_28 FALSE
#> Gene_29 FALSE
#> Gene_30 FALSE
#> Gene_31 FALSE
#> Gene_32 FALSE
#> Gene_33 FALSE
#> Gene_34 FALSE
#> Gene_35 FALSE
#> Gene_36 FALSE
#> Gene_37 FALSE
#> Gene_38 FALSE
#> Gene_39 FALSE
#> Gene_40 FALSE
#> Gene_41 FALSE
#> Gene_42 FALSE
#> Gene_43 FALSE
#> Gene_44 FALSE
#> Gene_45 FALSE
#> Gene_46 FALSE
#> Gene_47 FALSE
#> Gene_48 FALSE
#> Gene_49 FALSE
#> Gene_50 FALSE
#> Gene_51 FALSE
#> Gene_52 FALSE
#> Gene_53 FALSE
#> Gene_54 FALSE
#> Gene_55 FALSE
#> Gene_56 FALSE
#> Gene_57 FALSE
#> Gene_58 FALSE
#> Gene_59  TRUE
#> Gene_60 FALSE
#> 
#> 
#> $mod1_unimodal$biological_signal$de_results_sim
#> $mod1_unimodal$biological_signal$de_results_sim$TypeA_vs_TypeB
#>         feature        logFC   logCPM       F_stat       PValue        FDR
#> Gene_01 Gene_01  0.375749645 14.49426 1.192262e+00 0.2770610558 0.69754319
#> Gene_02 Gene_02  1.075103099 14.62498 9.710716e+00 0.0022925114 0.05920757
#> Gene_03 Gene_03  0.498288251 14.30499 2.014196e+00 0.1584261280 0.55915104
#> Gene_04 Gene_04  0.883949172 14.45905 6.417335e+00 0.0125905964 0.12590596
#> Gene_05 Gene_05  0.204096815 14.14743 3.309457e-01 0.5661799052 0.93260381
#> Gene_06 Gene_06  0.997846175 14.60584 8.372152e+00 0.0045264471 0.06789671
#> Gene_07 Gene_07 -0.028119282 14.70488 6.986109e-03 0.9335272965 0.98927371
#> Gene_08 Gene_08  0.322259554 14.25937 8.393349e-01 0.3614245040 0.82247059
#> Gene_09 Gene_09  0.732128103 14.83749 4.793690e+00 0.0304990871 0.23961319
#> Gene_10 Gene_10 -0.031803078 14.41976 8.461955e-03 0.9268602616 0.98927371
#> Gene_11 Gene_11  1.235027561 14.50173 1.236687e+01 0.0006177416 0.03706450
#> Gene_12 Gene_12  0.400141390 14.77652 1.428736e+00 0.2343276286 0.64360885
#> Gene_13 Gene_13  1.123895479 13.76307 9.204300e+00 0.0029603783 0.05920757
#> Gene_14 Gene_14 -0.071604928 13.91126 3.961578e-02 0.8425705368 0.97219677
#> Gene_15 Gene_15  0.691117584 14.34789 3.875465e+00 0.0513021017 0.34201401
#> Gene_16 Gene_16 -0.628533862 14.00201 3.052444e+00 0.0831741068 0.38388049
#> Gene_17 Gene_17 -0.026742474 13.71780 5.403187e-03 0.9415255584 0.98927371
#> Gene_18 Gene_18 -0.121627951 12.87738 9.661010e-02 0.7564754619 0.93260381
#> Gene_19 Gene_19 -0.510395126 13.15749 1.795192e+00 0.1828265768 0.60942192
#> Gene_20 Gene_20 -0.091238584 14.39624 6.936641e-02 0.7927132365 0.93260381
#> Gene_21 Gene_21  0.161695977 13.62302 1.934948e-01 0.6608149241 0.93260381
#> Gene_22 Gene_22  0.007704124 13.28375 4.218245e-04 0.9836480175 0.99457218
#> Gene_23 Gene_23 -0.145636522 13.81768 1.628195e-01 0.6872919150 0.93260381
#> Gene_24 Gene_24 -0.509428601 14.44576 2.169023e+00 0.1434337163 0.53787644
#> Gene_25 Gene_25 -0.194645496 14.33385 3.115047e-01 0.5777981403 0.93260381
#> Gene_26 Gene_26  0.372428007 13.65403 1.026974e+00 0.3129094871 0.75098277
#> Gene_27 Gene_27 -0.192095382 13.94261 2.858004e-01 0.5939135120 0.93260381
#> Gene_28 Gene_28 -0.821253386 13.32018 4.710515e+00 0.0319484260 0.23961319
#> Gene_29 Gene_29 -0.205193498 14.31062 3.446941e-01 0.5582351552 0.93260381
#> Gene_30 Gene_30 -0.323417193 13.96512 8.093585e-01 0.3701117645 0.82247059
#> Gene_31 Gene_31 -0.101507753 13.97750 8.008337e-02 0.7776720551 0.93260381
#> Gene_32 Gene_32  0.019655378 14.00818 3.015600e-03 0.9562979172 0.98927371
#> Gene_33 Gene_33  0.113548811 14.47825 1.089370e-01 0.7419331525 0.93260381
#> Gene_34 Gene_34 -0.665544984 13.88667 3.378597e+00 0.0685214872 0.37375357
#> Gene_35 Gene_35 -0.262969111 14.33085 5.680127e-01 0.4525259164 0.90505183
#> Gene_36 Gene_36 -0.249706324 13.31289 4.431156e-01 0.5069000588 0.93260381
#> Gene_37 Gene_37 -0.663062718 13.61473 3.224073e+00 0.0750794932 0.37539747
#> Gene_38 Gene_38  0.393653528 13.85198 1.182548e+00 0.2790172760 0.69754319
#> Gene_39 Gene_39 -0.054950891 12.90967 1.989400e-02 0.8880700229 0.98674447
#> Gene_40 Gene_40  0.206734173 13.72712 3.225255e-01 0.5711549217 0.93260381
#> Gene_41 Gene_41 -0.453939197 13.76217 1.564491e+00 0.2134421384 0.64032642
#> Gene_42 Gene_42 -0.182627595 13.92660 2.573636e-01 0.6128682463 0.93260381
#> Gene_43 Gene_43 -0.574297282 13.99965 2.544824e+00 0.1132866091 0.45314644
#> Gene_44 Gene_44 -1.001420190 13.80208 7.501064e+00 0.0071064840 0.08527781
#> Gene_45 Gene_45  0.283591497 13.57816 5.894405e-01 0.4441436660 0.90505183
#> Gene_46 Gene_46  0.113745898 13.83736 9.926623e-02 0.7532591832 0.93260381
#> Gene_47 Gene_47 -0.694913552 13.71083 3.600351e+00 0.0601703439 0.36102206
#> Gene_48 Gene_48 -0.054627545 14.34600 2.462523e-02 0.8755682381 0.98674447
#> Gene_49 Gene_49 -0.609329660 13.57247 2.709499e+00 0.1023689929 0.43872426
#> Gene_50 Gene_50  0.217457517 14.25783 3.828951e-01 0.5372307517 0.93260381
#> Gene_51 Gene_51  0.002459299 13.78958 4.647149e-05 0.9945721830 0.99457218
#> Gene_52 Gene_52 -0.464900292 14.04410 1.685465e+00 0.1966894015 0.62112443
#> Gene_53 Gene_53 -0.281102723 14.39668 6.569469e-01 0.4192443936 0.89838084
#> Gene_54 Gene_54  0.101191808 14.31455 8.397152e-02 0.7724858793 0.93260381
#> Gene_55 Gene_55 -0.151643149 13.92262 1.778540e-01 0.6739785682 0.93260381
#> Gene_56 Gene_56 -0.222506687 14.05566 3.884317e-01 0.5343076833 0.93260381
#> Gene_57 Gene_57  0.148717706 14.14822 1.759297e-01 0.6756446789 0.93260381
#> Gene_58 Gene_58  0.428087481 14.01531 1.418567e+00 0.2359899100 0.64360885
#> Gene_59 Gene_59 -0.171614401 12.70550 1.839392e-01 0.6687784834 0.93260381
#> Gene_60 Gene_60 -0.119522083 14.14245 1.136218e-01 0.7366467302 0.93260381
#>         is_da
#> Gene_01 FALSE
#> Gene_02 FALSE
#> Gene_03 FALSE
#> Gene_04 FALSE
#> Gene_05 FALSE
#> Gene_06 FALSE
#> Gene_07 FALSE
#> Gene_08 FALSE
#> Gene_09 FALSE
#> Gene_10 FALSE
#> Gene_11  TRUE
#> Gene_12 FALSE
#> Gene_13 FALSE
#> Gene_14 FALSE
#> Gene_15 FALSE
#> Gene_16 FALSE
#> Gene_17 FALSE
#> Gene_18 FALSE
#> Gene_19 FALSE
#> Gene_20 FALSE
#> Gene_21 FALSE
#> Gene_22 FALSE
#> Gene_23 FALSE
#> Gene_24 FALSE
#> Gene_25 FALSE
#> Gene_26 FALSE
#> Gene_27 FALSE
#> Gene_28 FALSE
#> Gene_29 FALSE
#> Gene_30 FALSE
#> Gene_31 FALSE
#> Gene_32 FALSE
#> Gene_33 FALSE
#> Gene_34 FALSE
#> Gene_35 FALSE
#> Gene_36 FALSE
#> Gene_37 FALSE
#> Gene_38 FALSE
#> Gene_39 FALSE
#> Gene_40 FALSE
#> Gene_41 FALSE
#> Gene_42 FALSE
#> Gene_43 FALSE
#> Gene_44 FALSE
#> Gene_45 FALSE
#> Gene_46 FALSE
#> Gene_47 FALSE
#> Gene_48 FALSE
#> Gene_49 FALSE
#> Gene_50 FALSE
#> Gene_51 FALSE
#> Gene_52 FALSE
#> Gene_53 FALSE
#> Gene_54 FALSE
#> Gene_55 FALSE
#> Gene_56 FALSE
#> Gene_57 FALSE
#> Gene_58 FALSE
#> Gene_59 FALSE
#> Gene_60 FALSE
#> 
#> 
#> $mod1_unimodal$biological_signal$ml_classification
#> $mod1_unimodal$biological_signal$ml_classification$accuracy
#> [1] 0.6875
#> 
#> $mod1_unimodal$biological_signal$ml_classification$precision
#> [1] 0.6875
#> 
#> $mod1_unimodal$biological_signal$ml_classification$recall
#> [1] 0.6904762
#> 
#> $mod1_unimodal$biological_signal$ml_classification$F1
#> [1] 0.6862745
#> 
#> 
#> 
#> $mod1_unimodal$trajectory
#> NULL
#> 
#> $mod1_unimodal$resource_usage
#> $mod1_unimodal$resource_usage$core_seconds
#> NULL
#> 
#> $mod1_unimodal$resource_usage$memory_per_core_mib
#> NULL
#> 
#> $mod1_unimodal$resource_usage$memory_per_core_gib
#> NULL
#> 
#> $mod1_unimodal$resource_usage$memory_mb
#> NULL
#> 
#> $mod1_unimodal$resource_usage$memory_input
#> NULL
#> 
#> $mod1_unimodal$resource_usage$memory_unit
#> [1] "MiB"
#> 
#> $mod1_unimodal$resource_usage$elapsed_time
#> NULL
#> 
#> $mod1_unimodal$resource_usage$elapsed_time_input
#> NULL
#> 
#> $mod1_unimodal$resource_usage$time_unit
#> [1] "seconds"
#> 
#> $mod1_unimodal$resource_usage$sim_cores
#> NULL
#> 
#> 
#> 
#> $mod2_unimodal
#> $mod2_unimodal$metrics_summary_table
#>                                        Category
#> 1                     Distributional Properties
#> 2                     Distributional Properties
#> 3                     Distributional Properties
#> 4                     Distributional Properties
#> 5                     Distributional Properties
#> 6                     Distributional Properties
#> 7                     Distributional Properties
#> 8                     Distributional Properties
#> 9                     Distributional Properties
#> 10                    Distributional Properties
#> 11                    Distributional Properties
#> 12                    Distributional Properties
#> 13                    Distributional Properties
#> 14                    Distributional Properties
#> 15                    Distributional Properties
#> 16                    Distributional Properties
#> 17                    Distributional Properties
#> 18                    Distributional Properties
#> 19                    Distributional Properties
#> 20                    Distributional Properties
#> 21                    Distributional Properties
#> 22                    Distributional Properties
#> 23                    Distributional Properties
#> 24                    Distributional Properties
#> 25                    Distributional Properties
#> 26                    Distributional Properties
#> 27                    Distributional Properties
#> 28                    Distributional Properties
#> 29                    Distributional Properties
#> 30                    Distributional Properties
#> 31                    Distributional Properties
#> 32                    Distributional Properties
#> 33                    Distributional Properties
#> 34                    Distributional Properties
#> 35                    Distributional Properties
#> 36                    Distributional Properties
#> 37                    Distributional Properties
#> 38                    Distributional Properties
#> 39                    Distributional Properties
#> 40                    Distributional Properties
#> 41                    Distributional Properties
#> 42                    Distributional Properties
#> 43                    Distributional Properties
#> 44                    Distributional Properties
#> 45                    Distributional Properties
#> 46                    Distributional Properties
#> 47                    Distributional Properties
#> 48                    Distributional Properties
#> 49                    Distributional Properties
#> 50                    Distributional Properties
#> 51                    Distributional Properties
#> 52                    Distributional Properties
#> 53                    Distributional Properties
#> 54                    Distributional Properties
#> 55                    Distributional Properties
#> 56                    Distributional Properties
#> 57                    Distributional Properties
#> 58                    Distributional Properties
#> 59                    Distributional Properties
#> 60                    Distributional Properties
#> 61                    Distributional Properties
#> 62                    Distributional Properties
#> 63                    Distributional Properties
#> 64                    Distributional Properties
#> 65                    Distributional Properties
#> 66                    Distributional Properties
#> 67                    Distributional Properties
#> 68                    Distributional Properties
#> 69                    Distributional Properties
#> 70                    Distributional Properties
#> 71                    Distributional Properties
#> 72                    Distributional Properties
#> 73                    Distributional Properties
#> 74                    Distributional Properties
#> 75                    Distributional Properties
#> 76                    Distributional Properties
#> 77                    Distributional Properties
#> 78                    Distributional Properties
#> 79                    Distributional Properties
#> 80                    Distributional Properties
#> 81                    Distributional Properties
#> 82                    Distributional Properties
#> 83                    Distributional Properties
#> 84                    Distributional Properties
#> 85                    Distributional Properties
#> 86                    Distributional Properties
#> 87                    Distributional Properties
#> 88                    Distributional Properties
#> 89                    Distributional Properties
#> 90                    Distributional Properties
#> 91                    Distributional Properties
#> 92                    Distributional Properties
#> 93                    Distributional Properties
#> 94                    Distributional Properties
#> 95                    Distributional Properties
#> 96                    Distributional Properties
#> 97                    Distributional Properties
#> 98                    Distributional Properties
#> 99                    Distributional Properties
#> 100                   Distributional Properties
#> 101                   Distributional Properties
#> 102                   Distributional Properties
#> 103                   Distributional Properties
#> 104                   Distributional Properties
#> 105                   Distributional Properties
#> 106                   Distributional Properties
#> 107                   Distributional Properties
#> 108                   Distributional Properties
#> 109                   Distributional Properties
#> 110                   Distributional Properties
#> 111                   Distributional Properties
#> 112                   Distributional Properties
#> 113                   Distributional Properties
#> 114                   Distributional Properties
#> 115                   Distributional Properties
#> 116                   Distributional Properties
#> 117                   Distributional Properties
#> 118                   Distributional Properties
#> 119                   Distributional Properties
#> 120                   Distributional Properties
#> 121                   Distributional Properties
#> 122                   Distributional Properties
#> 123                   Distributional Properties
#> 124                   Distributional Properties
#> 125                   Distributional Properties
#> 126                   Distributional Properties
#> 127                   Distributional Properties
#> 128                   Distributional Properties
#> 129                   Distributional Properties
#> 130                   Distributional Properties
#> 131          Sparsity & Zero-Inflation Dynamics
#> 132          Sparsity & Zero-Inflation Dynamics
#> 133          Sparsity & Zero-Inflation Dynamics
#> 134  Cellular Manifold & Clustering Concordance
#> 135  Cellular Manifold & Clustering Concordance
#> 136  Cellular Manifold & Clustering Concordance
#> 137  Cellular Manifold & Clustering Concordance
#> 138  Cellular Manifold & Clustering Concordance
#> 139  Cellular Manifold & Clustering Concordance
#> 140 Technical Variation & Batch Mixing Dynamics
#> 141 Technical Variation & Batch Mixing Dynamics
#> 142 Technical Variation & Batch Mixing Dynamics
#> 143 Technical Variation & Batch Mixing Dynamics
#> 144 Technical Variation & Batch Mixing Dynamics
#> 145   Biological Signal & Differential Fidelity
#> 146   Biological Signal & Differential Fidelity
#> 147   Biological Signal & Differential Fidelity
#> 148   Biological Signal & Differential Fidelity
#> 149   Biological Signal & Differential Fidelity
#> 150   Biological Signal & Differential Fidelity
#>                              Property                         Metric
#> 1                        library_size                    Wasserstein
#> 2                        library_size                             KS
#> 3                        library_size                            MAD
#> 4                        library_size                  Bhattacharyya
#> 5                        library_size                 Runs_Statistic
#> 6                        library_size                    NN_Mismatch
#> 7                        library_size    Between_Dataset_Silh_Global
#> 8                        library_size                            MAE
#> 9                        library_size                           RMSE
#> 10                       library_size                        KDE_ISE
#> 11                   log_library_size                    Wasserstein
#> 12                   log_library_size                             KS
#> 13                   log_library_size                            MAD
#> 14                   log_library_size                  Bhattacharyya
#> 15                   log_library_size                 Runs_Statistic
#> 16                   log_library_size                    NN_Mismatch
#> 17                   log_library_size    Between_Dataset_Silh_Global
#> 18                   log_library_size                            MAE
#> 19                   log_library_size                           RMSE
#> 20                   log_library_size                        KDE_ISE
#> 21                 zero_fraction_cell                    Wasserstein
#> 22                 zero_fraction_cell                             KS
#> 23                 zero_fraction_cell                            MAD
#> 24                 zero_fraction_cell                  Bhattacharyya
#> 25                 zero_fraction_cell                 Runs_Statistic
#> 26                 zero_fraction_cell                    NN_Mismatch
#> 27                 zero_fraction_cell    Between_Dataset_Silh_Global
#> 28                 zero_fraction_cell                            MAE
#> 29                 zero_fraction_cell                           RMSE
#> 30                 zero_fraction_cell                        KDE_ISE
#> 31                detection_freq_cell                    Wasserstein
#> 32                detection_freq_cell                             KS
#> 33                detection_freq_cell                            MAD
#> 34                detection_freq_cell                  Bhattacharyya
#> 35                detection_freq_cell                 Runs_Statistic
#> 36                detection_freq_cell                    NN_Mismatch
#> 37                detection_freq_cell    Between_Dataset_Silh_Global
#> 38                detection_freq_cell                            MAE
#> 39                detection_freq_cell                           RMSE
#> 40                detection_freq_cell                        KDE_ISE
#> 41                           cell_cor                    Wasserstein
#> 42                           cell_cor                             KS
#> 43                           cell_cor                            MAD
#> 44                           cell_cor                  Bhattacharyya
#> 45                           cell_cor                 Runs_Statistic
#> 46                           cell_cor                    NN_Mismatch
#> 47                           cell_cor    Between_Dataset_Silh_Global
#> 48                           cell_cor                            MAE
#> 49                           cell_cor                           RMSE
#> 50                           cell_cor                        KDE_ISE
#> 51                           pca_dist                    Wasserstein
#> 52                           pca_dist                             KS
#> 53                           pca_dist                            MAD
#> 54                           pca_dist                  Bhattacharyya
#> 55                           pca_dist                 Runs_Statistic
#> 56                           pca_dist                    NN_Mismatch
#> 57                           pca_dist    Between_Dataset_Silh_Global
#> 58                           pca_dist                            MAE
#> 59                           pca_dist                           RMSE
#> 60                           pca_dist                        KDE_ISE
#> 61                    mean_expression                    Wasserstein
#> 62                    mean_expression                             KS
#> 63                    mean_expression                            MAD
#> 64                    mean_expression                  Bhattacharyya
#> 65                    mean_expression                 Runs_Statistic
#> 66                    mean_expression                    NN_Mismatch
#> 67                    mean_expression    Between_Dataset_Silh_Global
#> 68                    mean_expression                            MAE
#> 69                    mean_expression                           RMSE
#> 70                    mean_expression                        KDE_ISE
#> 71                           variance                    Wasserstein
#> 72                           variance                             KS
#> 73                           variance                            MAD
#> 74                           variance                  Bhattacharyya
#> 75                           variance                 Runs_Statistic
#> 76                           variance                    NN_Mismatch
#> 77                           variance    Between_Dataset_Silh_Global
#> 78                           variance                            MAE
#> 79                           variance                           RMSE
#> 80                           variance                        KDE_ISE
#> 81                    scaled_variance                    Wasserstein
#> 82                    scaled_variance                             KS
#> 83                    scaled_variance                            MAD
#> 84                    scaled_variance                  Bhattacharyya
#> 85                    scaled_variance                 Runs_Statistic
#> 86                    scaled_variance                    NN_Mismatch
#> 87                    scaled_variance    Between_Dataset_Silh_Global
#> 88                    scaled_variance                            MAE
#> 89                    scaled_variance                           RMSE
#> 90                    scaled_variance                        KDE_ISE
#> 91                                 cv                    Wasserstein
#> 92                                 cv                             KS
#> 93                                 cv                            MAD
#> 94                                 cv                  Bhattacharyya
#> 95                                 cv                 Runs_Statistic
#> 96                                 cv                    NN_Mismatch
#> 97                                 cv    Between_Dataset_Silh_Global
#> 98                                 cv                            MAE
#> 99                                 cv                           RMSE
#> 100                                cv                        KDE_ISE
#> 101             zero_fraction_feature                    Wasserstein
#> 102             zero_fraction_feature                             KS
#> 103             zero_fraction_feature                            MAD
#> 104             zero_fraction_feature                  Bhattacharyya
#> 105             zero_fraction_feature                 Runs_Statistic
#> 106             zero_fraction_feature                    NN_Mismatch
#> 107             zero_fraction_feature    Between_Dataset_Silh_Global
#> 108             zero_fraction_feature                            MAE
#> 109             zero_fraction_feature                           RMSE
#> 110             zero_fraction_feature                        KDE_ISE
#> 111                        dispersion                    Wasserstein
#> 112                        dispersion                             KS
#> 113                        dispersion                            MAD
#> 114                        dispersion                  Bhattacharyya
#> 115                        dispersion                 Runs_Statistic
#> 116                        dispersion                    NN_Mismatch
#> 117                        dispersion    Between_Dataset_Silh_Global
#> 118                        dispersion                            MAE
#> 119                        dispersion                           RMSE
#> 120                        dispersion                        KDE_ISE
#> 121                          gene_cor                    Wasserstein
#> 122                          gene_cor                             KS
#> 123                          gene_cor                            MAD
#> 124                          gene_cor                  Bhattacharyya
#> 125                          gene_cor                 Runs_Statistic
#> 126                          gene_cor                    NN_Mismatch
#> 127                          gene_cor    Between_Dataset_Silh_Global
#> 128                          gene_cor                            MAE
#> 129                          gene_cor                           RMSE
#> 130                          gene_cor                        KDE_ISE
#> 131      accessibility_sparsity_curve                     curve_rmse
#> 132      accessibility_sparsity_curve                 r2_discrepancy
#> 133               excess_zero_weights excess_zero_weight_discrepancy
#> 134                          manifold                            MMD
#> 135                          manifold           Generative_Precision
#> 136                          manifold              Generative_Recall
#> 137                        clustering                 silhouette_sim
#> 138                        clustering                            cdi
#> 139                        clustering            neighborhood_purity
#> 140                      batch_mixing               batch_silhouette
#> 141                      batch_mixing                         pcr_r2
#> 142                      batch_mixing                            cms
#> 143                      batch_mixing                            isi
#> 144                      batch_mixing                 cross_batch_F1
#> 145 scATAC_differential_accessibility                pval_dist_score
#> 146 scATAC_differential_accessibility                       accuracy
#> 147 scATAC_differential_accessibility                      precision
#> 148 scATAC_differential_accessibility                         recall
#> 149 scATAC_differential_accessibility                       f1_score
#> 150 scATAC_differential_accessibility                            auc
#>             Value
#> 1    2.337500e+00
#> 2    1.750000e-01
#> 3    2.000000e+00
#> 4    1.055641e-04
#> 5   -6.027311e+00
#> 6    1.937500e-01
#> 7    2.188673e-02
#> 8    2.337500e+00
#> 9    2.487469e+00
#> 10   4.478217e-03
#> 11   7.991730e-02
#> 12   1.750000e-01
#> 13   8.223810e-02
#> 14   1.155313e-05
#> 15  -6.027311e+00
#> 16   1.937500e-01
#> 17   1.990590e-02
#> 18   7.991730e-02
#> 19   8.389317e-02
#> 20   1.268716e-01
#> 21   2.750000e-02
#> 22   2.125000e-01
#> 23   3.333333e-02
#> 24   6.907067e-05
#> 25  -7.296219e+00
#> 26   2.812500e-01
#> 27   2.566981e-02
#> 28   2.750000e-02
#> 29   3.016160e-02
#> 30   4.176638e-01
#> 31   2.750000e-02
#> 32   2.125000e-01
#> 33   3.333333e-02
#> 34   1.795082e-04
#> 35  -6.978992e+00
#> 36   2.812500e-01
#> 37   2.566981e-02
#> 38   2.750000e-02
#> 39   3.016160e-02
#> 40   4.176638e-01
#> 41   3.342856e-03
#> 42   2.056815e-02
#> 43   2.270508e-03
#> 44   5.400110e-05
#> 45  -2.515972e-01
#> 46   4.333333e-02
#> 47  -2.274271e-03
#> 48   3.342856e-03
#> 49   6.619754e-03
#> 50   4.263335e-03
#> 51   1.502554e-01
#> 52   2.948321e-02
#> 53   1.327661e-01
#> 54   2.878871e-06
#> 55   7.547916e-01
#> 56   2.666667e-02
#> 57   6.595987e-04
#> 58   1.502554e-01
#> 59   1.883752e-01
#> 60   3.459236e-04
#> 61   1.501367e-01
#> 62   1.879472e-01
#> 63   1.536693e-01
#> 64   2.959894e-04
#> 65  -1.833462e+00
#> 66   2.500000e-02
#> 67   3.291150e-03
#> 68   1.501367e-01
#> 69   1.677762e-01
#> 70   2.809910e-02
#> 71   2.170317e-01
#> 72   2.753225e-01
#> 73   2.419452e-01
#> 74   5.389967e-05
#> 75  -2.566846e+00
#> 76   1.166667e-01
#> 77   3.168098e-02
#> 78   2.170317e-01
#> 79   2.541541e-01
#> 80   1.429083e-01
#> 81   1.993187e-01
#> 82   1.174387e-01
#> 83   1.419895e-01
#> 84   2.917823e-03
#> 85  -9.167309e-01
#> 86   8.333333e-02
#> 87  -3.112124e-03
#> 88   1.993187e-01
#> 89   2.711405e-01
#> 90   3.669152e-02
#> 91   8.864353e+00
#> 92   1.763639e-01
#> 93   8.752987e+00
#> 94   2.185211e-04
#> 95  -1.833462e+00
#> 96   3.333333e-02
#> 97   6.571696e-03
#> 98   8.864353e+00
#> 99   1.046915e+01
#> 100  7.968311e-04
#> 101  2.958333e-02
#> 102  2.000000e-01
#> 103  2.500000e-02
#> 104  1.230865e-04
#> 105 -1.650116e+00
#> 106  5.000000e-02
#> 107  7.445303e-03
#> 108  2.958333e-02
#> 109  3.307189e-02
#> 110  1.835132e-01
#> 111  2.544810e-02
#> 112  1.136364e-01
#> 113  2.119535e-02
#> 114  1.193835e-04
#> 115 -5.500385e-01
#> 116  8.333333e-03
#> 117 -7.584339e-04
#> 118  2.544810e-02
#> 119  2.978845e-02
#> 120  1.190016e-01
#> 121  2.718776e-03
#> 122  1.759557e-02
#> 123  1.838118e-03
#> 124  8.878612e-05
#> 125 -6.723876e-02
#> 126  4.666667e-02
#> 127 -2.426294e-04
#> 128  2.718776e-03
#> 129  5.260476e-03
#> 130  6.132365e-03
#> 131  1.233884e-02
#> 132  8.535249e-04
#> 133  8.480348e-04
#> 134  1.096098e-01
#> 135  1.000000e+00
#> 136  9.500000e-01
#> 137  4.149183e-02
#> 138  1.687605e+00
#> 139  7.275000e-01
#> 140  8.780031e-04
#> 141  1.375709e-02
#> 142  0.000000e+00
#> 143  1.931492e+00
#> 144  7.346700e-01
#> 145  8.563379e-01
#> 146  8.166667e-01
#> 147  5.000000e-01
#> 148  9.090909e-02
#> 149  1.538462e-01
#> 150  8.367347e-01
#> 
#> $mod2_unimodal$cell_metrics
#> $mod2_unimodal$cell_metrics$library_size
#> $mod2_unimodal$cell_metrics$library_size$library_size_MAD
#> [1] 2
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_KS
#> [1] 0.175
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_MAE
#> [1] 2.3375
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_RMSE
#> [1] 2.487469
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_Bhattacharyya
#> [1] 0.0001055641
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_Wasserstein
#> [1] 2.3375
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_Runs_Statistic
#> [1] -6.027311
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_Runs_PValue
#> [1] 8.335488e-10
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_NN_Mismatch
#> [1] 0.19375
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_Between_Dataset_Silh_Global
#> [1] 0.02188673
#> 
#> $mod2_unimodal$cell_metrics$library_size$library_size_KDE_ISE
#> [1] 0.004478217
#> 
#> 
#> $mod2_unimodal$cell_metrics$log_library_size
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_MAD
#> [1] 0.0822381
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_KS
#> [1] 0.175
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_MAE
#> [1] 0.0799173
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_RMSE
#> [1] 0.08389317
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_Bhattacharyya
#> [1] 1.155313e-05
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_Wasserstein
#> [1] 0.0799173
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_Runs_Statistic
#> [1] -6.027311
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_Runs_PValue
#> [1] 8.335488e-10
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_NN_Mismatch
#> [1] 0.19375
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_Between_Dataset_Silh_Global
#> [1] 0.0199059
#> 
#> $mod2_unimodal$cell_metrics$log_library_size$log_library_size_KDE_ISE
#> [1] 0.1268716
#> 
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_MAD
#> [1] 0.03333333
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_KS
#> [1] 0.2125
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_MAE
#> [1] 0.0275
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_RMSE
#> [1] 0.0301616
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Bhattacharyya
#> [1] 6.907067e-05
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Wasserstein
#> [1] 0.0275
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Runs_Statistic
#> [1] -7.296219
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Runs_PValue
#> [1] 1.479836e-13
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_NN_Mismatch
#> [1] 0.28125
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_Between_Dataset_Silh_Global
#> [1] 0.02566981
#> 
#> $mod2_unimodal$cell_metrics$zero_fraction_cell$zero_fraction_cell_KDE_ISE
#> [1] 0.4176638
#> 
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_MAD
#> [1] 0.03333333
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_KS
#> [1] 0.2125
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_MAE
#> [1] 0.0275
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_RMSE
#> [1] 0.0301616
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Bhattacharyya
#> [1] 0.0001795082
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Wasserstein
#> [1] 0.0275
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Runs_Statistic
#> [1] -6.978992
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Runs_PValue
#> [1] 1.486527e-12
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_NN_Mismatch
#> [1] 0.28125
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_Between_Dataset_Silh_Global
#> [1] 0.02566981
#> 
#> $mod2_unimodal$cell_metrics$detection_freq_cell$detection_freq_cell_KDE_ISE
#> [1] 0.4176638
#> 
#> 
#> $mod2_unimodal$cell_metrics$cell_cor
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_MAD
#> [1] 0.002270508
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_KS
#> [1] 0.02056815
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_MAE
#> [1] 0.003342856
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_RMSE
#> [1] 0.006619754
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_Bhattacharyya
#> [1] 5.40011e-05
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_Wasserstein
#> [1] 0.003342856
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_Runs_Statistic
#> [1] -0.2515972
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_Runs_PValue
#> [1] 0.4006762
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_NN_Mismatch
#> [1] 0.04333333
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_Between_Dataset_Silh_Global
#> [1] -0.002274271
#> 
#> $mod2_unimodal$cell_metrics$cell_cor$cell_cor_KDE_ISE
#> [1] 0.004263335
#> 
#> 
#> $mod2_unimodal$cell_metrics$pca_dist
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_MAD
#> [1] 0.1327661
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_KS
#> [1] 0.02948321
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_MAE
#> [1] 0.1502554
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_RMSE
#> [1] 0.1883752
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_Bhattacharyya
#> [1] 2.878871e-06
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_Wasserstein
#> [1] 0.1502554
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_Runs_Statistic
#> [1] 0.7547916
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_Runs_PValue
#> [1] 0.774813
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_NN_Mismatch
#> [1] 0.02666667
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_Between_Dataset_Silh_Global
#> [1] 0.0006595987
#> 
#> $mod2_unimodal$cell_metrics$pca_dist$pca_dist_KDE_ISE
#> [1] 0.0003459236
#> 
#> 
#> 
#> $mod2_unimodal$feature_metrics
#> $mod2_unimodal$feature_metrics$mean_expression
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_MAD
#> [1] 0.1536693
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_KS
#> [1] 0.1879472
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_MAE
#> [1] 0.1501367
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_RMSE
#> [1] 0.1677762
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_Bhattacharyya
#> [1] 0.0002959894
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_Wasserstein
#> [1] 0.1501367
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_Runs_Statistic
#> [1] -1.833462
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_Runs_PValue
#> [1] 0.03336697
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_NN_Mismatch
#> [1] 0.025
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_Between_Dataset_Silh_Global
#> [1] 0.00329115
#> 
#> $mod2_unimodal$feature_metrics$mean_expression$mean_expression_KDE_ISE
#> [1] 0.0280991
#> 
#> 
#> $mod2_unimodal$feature_metrics$variance
#> $mod2_unimodal$feature_metrics$variance$variance_MAD
#> [1] 0.2419452
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_KS
#> [1] 0.2753225
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_MAE
#> [1] 0.2170317
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_RMSE
#> [1] 0.2541541
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_Bhattacharyya
#> [1] 5.389967e-05
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_Wasserstein
#> [1] 0.2170317
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_Runs_Statistic
#> [1] -2.566846
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_Runs_PValue
#> [1] 0.005131403
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_NN_Mismatch
#> [1] 0.1166667
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_Between_Dataset_Silh_Global
#> [1] 0.03168098
#> 
#> $mod2_unimodal$feature_metrics$variance$variance_KDE_ISE
#> [1] 0.1429083
#> 
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_MAD
#> [1] 0.1419895
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_KS
#> [1] 0.1174387
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_MAE
#> [1] 0.1993187
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_RMSE
#> [1] 0.2711405
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_Bhattacharyya
#> [1] 0.002917823
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_Wasserstein
#> [1] 0.1993187
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_Runs_Statistic
#> [1] -0.9167309
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_Runs_PValue
#> [1] 0.1796418
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_NN_Mismatch
#> [1] 0.08333333
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_Between_Dataset_Silh_Global
#> [1] -0.003112124
#> 
#> $mod2_unimodal$feature_metrics$scaled_variance$scaled_variance_KDE_ISE
#> [1] 0.03669152
#> 
#> 
#> $mod2_unimodal$feature_metrics$cv
#> $mod2_unimodal$feature_metrics$cv$cv_MAD
#> [1] 8.752987
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_KS
#> [1] 0.1763639
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_MAE
#> [1] 8.864353
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_RMSE
#> [1] 10.46915
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_Bhattacharyya
#> [1] 0.0002185211
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_Wasserstein
#> [1] 8.864353
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_Runs_Statistic
#> [1] -1.833462
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_Runs_PValue
#> [1] 0.03336697
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_NN_Mismatch
#> [1] 0.03333333
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_Between_Dataset_Silh_Global
#> [1] 0.006571696
#> 
#> $mod2_unimodal$feature_metrics$cv$cv_KDE_ISE
#> [1] 0.0007968311
#> 
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_MAD
#> [1] 0.025
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_KS
#> [1] 0.2
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_MAE
#> [1] 0.02958333
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_RMSE
#> [1] 0.03307189
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Bhattacharyya
#> [1] 0.0001230865
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Wasserstein
#> [1] 0.02958333
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Runs_Statistic
#> [1] -1.650116
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Runs_PValue
#> [1] 0.04945965
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_NN_Mismatch
#> [1] 0.05
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_Between_Dataset_Silh_Global
#> [1] 0.007445303
#> 
#> $mod2_unimodal$feature_metrics$zero_fraction_feature$zero_fraction_feature_KDE_ISE
#> [1] 0.1835132
#> 
#> 
#> $mod2_unimodal$feature_metrics$dispersion
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_MAD
#> [1] 0.02119535
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_KS
#> [1] 0.1136364
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_MAE
#> [1] 0.0254481
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_RMSE
#> [1] 0.02978845
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_Bhattacharyya
#> [1] 0.0001193835
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_Wasserstein
#> [1] 0.0254481
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_Runs_Statistic
#> [1] -0.5500385
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_Runs_PValue
#> [1] 0.2911465
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_NN_Mismatch
#> [1] 0.008333333
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_Between_Dataset_Silh_Global
#> [1] -0.0007584339
#> 
#> $mod2_unimodal$feature_metrics$dispersion$dispersion_KDE_ISE
#> [1] 0.1190016
#> 
#> 
#> $mod2_unimodal$feature_metrics$gene_cor
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_MAD
#> [1] 0.001838118
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_KS
#> [1] 0.01759557
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_MAE
#> [1] 0.002718776
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_RMSE
#> [1] 0.005260476
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_Bhattacharyya
#> [1] 8.878612e-05
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_Wasserstein
#> [1] 0.002718776
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_Runs_Statistic
#> [1] -0.06723876
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_Runs_PValue
#> [1] 0.4731958
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_NN_Mismatch
#> [1] 0.04666667
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_Between_Dataset_Silh_Global
#> [1] -0.0002426294
#> 
#> $mod2_unimodal$feature_metrics$gene_cor$gene_cor_KDE_ISE
#> [1] 0.006132365
#> 
#> 
#> 
#> $mod2_unimodal$bivariate_metrics
#> list()
#> 
#> $mod2_unimodal$zero_probability_curve
#> NULL
#> 
#> $mod2_unimodal$accessibility_sparsity_curve
#> $mod2_unimodal$accessibility_sparsity_curve$curve_rmse
#> [1] 0.01233884
#> 
#> $mod2_unimodal$accessibility_sparsity_curve$r2_discrepancy
#> [1] 0.0008535249
#> 
#> $mod2_unimodal$accessibility_sparsity_curve$ref_r_squared
#> [1] 0.9692997
#> 
#> $mod2_unimodal$accessibility_sparsity_curve$sim_r_squared
#> [1] 0.9684462
#> 
#> 
#> $mod2_unimodal$excess_zero_weights
#> $mod2_unimodal$excess_zero_weights$ref
#> $mod2_unimodal$excess_zero_weights$ref$mean_excess_zero_weight
#> [1] 0.001210918
#> 
#> $mod2_unimodal$excess_zero_weights$ref$median_excess_zero_weight
#> [1] 0
#> 
#> $mod2_unimodal$excess_zero_weights$ref$mean_zero_inflation
#> [1] 0.0009000037
#> 
#> $mod2_unimodal$excess_zero_weights$ref$gene_excess_weights
#>  [1] 0.010238908 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [7] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [13] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [19] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [25] 0.039992389 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [31] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [37] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [43] 0.000000000 0.000000000 0.000000000 0.000000000 0.018656244 0.000000000
#> [49] 0.000000000 0.000000000 0.003767533 0.000000000 0.000000000 0.000000000
#> [55] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> 
#> 
#> $mod2_unimodal$excess_zero_weights$sim
#> $mod2_unimodal$excess_zero_weights$sim$mean_excess_zero_weight
#> [1] 0.002058953
#> 
#> $mod2_unimodal$excess_zero_weights$sim$median_excess_zero_weight
#> [1] 0
#> 
#> $mod2_unimodal$excess_zero_weights$sim$mean_zero_inflation
#> [1] 0.001469021
#> 
#> $mod2_unimodal$excess_zero_weights$sim$gene_excess_weights
#>  [1] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#>  [7] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [13] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [19] 0.000000000 0.010152284 0.000000000 0.000000000 0.000000000 0.000000000
#> [25] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.008808498
#> [31] 0.000000000 0.000000000 0.051350526 0.000000000 0.000000000 0.000000000
#> [37] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [43] 0.000000000 0.000000000 0.000000000 0.000000000 0.053225852 0.000000000
#> [49] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> [55] 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000 0.000000000
#> 
#> 
#> 
#> $mod2_unimodal$generative_manifold
#> $mod2_unimodal$generative_manifold$mmd
#> [1] 0.1096098
#> 
#> $mod2_unimodal$generative_manifold$precision_recall
#> $mod2_unimodal$generative_manifold$precision_recall$generative_precision
#> [1] 1
#> 
#> $mod2_unimodal$generative_manifold$precision_recall$generative_recall
#> [1] 0.95
#> 
#> $mod2_unimodal$generative_manifold$precision_recall$generative_f1
#> [1] 0.974359
#> 
#> 
#> 
#> $mod2_unimodal$ref_properties
#> $mod2_unimodal$ref_properties$cell
#> $mod2_unimodal$ref_properties$cell$library_size
#>  [1] 25 22 23 28 24 25 20 17 32 30 28 28 16 18 27 21 17 23 14 30 23 27 30 25 34
#> [26] 22 18 17 30 28 25 29 31 26 24 23 22 22 34 30 36 34 38 36 33 28 34 29 25 33
#> [51] 36 26 37 30 32 38 34 31 35 27 29 32 31 35 36 35 27 26 26 38 38 41 34 34 27
#> [76] 39 30 37 30 35
#> 
#> $mod2_unimodal$ref_properties$cell$log_library_size
#>  [1] 3.258097 3.135494 3.178054 3.367296 3.218876 3.258097 3.044522 2.890372
#>  [9] 3.496508 3.433987 3.367296 3.367296 2.833213 2.944439 3.332205 3.091042
#> [17] 2.890372 3.178054 2.708050 3.433987 3.178054 3.332205 3.433987 3.258097
#> [25] 3.555348 3.135494 2.944439 2.890372 3.433987 3.367296 3.258097 3.401197
#> [33] 3.465736 3.295837 3.218876 3.178054 3.135494 3.135494 3.555348 3.433987
#> [41] 3.610918 3.555348 3.663562 3.610918 3.526361 3.367296 3.555348 3.401197
#> [49] 3.258097 3.526361 3.610918 3.295837 3.637586 3.433987 3.496508 3.663562
#> [57] 3.555348 3.465736 3.583519 3.332205 3.401197 3.496508 3.465736 3.583519
#> [65] 3.610918 3.583519 3.332205 3.295837 3.295837 3.663562 3.663562 3.737670
#> [73] 3.555348 3.555348 3.332205 3.688879 3.433987 3.637586 3.433987 3.583519
#> 
#> $mod2_unimodal$ref_properties$cell$tmm_factor
#>  [1] 1.1183765 1.2708824 1.2156266 0.9985504 1.1649755 1.1183765 1.3979706
#>  [8] 1.6446713 0.8737316 0.9319804 0.9985504 0.9985504 1.7474633 1.5533007
#> [15] 1.0355338 1.3314006 1.6446713 1.2156266 1.9971009 0.9319804 1.2156266
#> [22] 1.0355338 0.9319804 1.1183765 0.8223357 1.2708824 1.5533007 1.6446713
#> [29] 0.9319804 0.9985504 1.1183765 0.9641177 0.9019165 1.0753620 1.1649755
#> [36] 1.2156266 1.2708824 1.2708824 0.8223357 0.9319804 0.7766503 0.8223357
#> [43] 0.7357740 0.7766503 0.8472549 0.9985504 0.9985504 0.9641177 1.1183765
#> [50] 0.8472549 0.7766503 1.0753620 0.7556598 0.9319804 0.8737316 0.7357740
#> [57] 0.8223357 0.9019165 0.7988403 1.0355338 0.9641177 0.8737316 0.9019165
#> [64] 0.7988403 0.7766503 0.7988403 1.0355338 1.0753620 1.0753620 0.7357740
#> [71] 0.7357740 0.6819369 0.8223357 0.8223357 1.0355338 0.7169080 0.9319804
#> [78] 0.7556598 0.9319804 0.7988403
#> 
#> $mod2_unimodal$ref_properties$cell$effective_library_size
#>  [1] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#>  [9] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [17] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [25] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [33] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [41] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 33.95071 27.95941
#> [49] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [57] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [65] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> [73] 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941 27.95941
#> 
#> $mod2_unimodal$ref_properties$cell$zero_fraction_cell
#>  [1] 0.6166667 0.6666667 0.6833333 0.5666667 0.6333333 0.6500000 0.6833333
#>  [8] 0.7166667 0.5333333 0.5666667 0.5833333 0.5666667 0.7500000 0.7333333
#> [15] 0.5500000 0.6833333 0.7166667 0.6500000 0.8000000 0.6166667 0.6666667
#> [22] 0.6166667 0.5166667 0.5833333 0.5000000 0.6833333 0.7166667 0.7166667
#> [29] 0.5666667 0.6166667 0.6333333 0.5833333 0.5666667 0.6166667 0.6333333
#> [36] 0.6333333 0.6666667 0.6833333 0.5500000 0.5500000 0.5500000 0.5666667
#> [43] 0.5000000 0.5000000 0.5500000 0.6000000 0.5833333 0.5833333 0.6500000
#> [50] 0.5500000 0.4833333 0.6166667 0.4500000 0.5666667 0.5500000 0.5166667
#> [57] 0.5500000 0.5333333 0.5666667 0.6500000 0.5500000 0.5666667 0.5333333
#> [64] 0.5000000 0.5166667 0.4833333 0.6000000 0.6166667 0.6000000 0.5500000
#> [71] 0.4833333 0.4166667 0.5166667 0.5833333 0.6000000 0.5000000 0.5500000
#> [78] 0.5000000 0.5666667 0.5000000
#> 
#> $mod2_unimodal$ref_properties$cell$detection_freq_cell
#>  [1] 0.3833333 0.3333333 0.3166667 0.4333333 0.3666667 0.3500000 0.3166667
#>  [8] 0.2833333 0.4666667 0.4333333 0.4166667 0.4333333 0.2500000 0.2666667
#> [15] 0.4500000 0.3166667 0.2833333 0.3500000 0.2000000 0.3833333 0.3333333
#> [22] 0.3833333 0.4833333 0.4166667 0.5000000 0.3166667 0.2833333 0.2833333
#> [29] 0.4333333 0.3833333 0.3666667 0.4166667 0.4333333 0.3833333 0.3666667
#> [36] 0.3666667 0.3333333 0.3166667 0.4500000 0.4500000 0.4500000 0.4333333
#> [43] 0.5000000 0.5000000 0.4500000 0.4000000 0.4166667 0.4166667 0.3500000
#> [50] 0.4500000 0.5166667 0.3833333 0.5500000 0.4333333 0.4500000 0.4833333
#> [57] 0.4500000 0.4666667 0.4333333 0.3500000 0.4500000 0.4333333 0.4666667
#> [64] 0.5000000 0.4833333 0.5166667 0.4000000 0.3833333 0.4000000 0.4500000
#> [71] 0.5166667 0.5833333 0.4833333 0.4166667 0.4000000 0.5000000 0.4500000
#> [78] 0.5000000 0.4333333 0.5000000
#> 
#> $mod2_unimodal$ref_properties$cell$cell_cor
#>    [1] -9.251908e-02 -5.438642e-02 -1.552964e-01  6.262229e-02  1.201706e-01
#>    [6] -4.354333e-02 -7.535929e-02  1.398382e-01 -9.825150e-02  2.391122e-01
#>   [11] -1.484335e-01 -7.943744e-02  5.550852e-02  3.333449e-01  1.158085e-01
#>   [16]  3.303659e-02  1.803289e-01 -1.394343e-01  3.523855e-01  8.650869e-02
#>   [21] -1.803503e-02  8.202538e-02  8.745531e-02 -5.468010e-02 -5.604029e-02
#>   [26] -5.576976e-02  8.679883e-02 -1.159254e-01 -1.447554e-01 -6.216982e-02
#>   [31]  1.146061e-01 -9.928915e-02  5.057136e-02  1.647377e-02 -2.045344e-01
#>   [36] -2.928104e-02 -1.695951e-01  1.274398e-02  1.649145e-02 -1.139664e-02
#>   [41]  4.380624e-02  6.981498e-02 -6.073913e-02 -1.344761e-01 -2.044857e-02
#>   [46] -4.557153e-02  4.302824e-02  7.804074e-03  1.129605e-01  6.502919e-02
#>   [51]  1.808987e-01 -7.916050e-02  2.537626e-02 -3.909387e-02  1.342005e-01
#>   [56]  7.324278e-02  4.460075e-02 -1.574482e-01  2.271607e-01  1.392028e-01
#>   [61]  1.668299e-01  3.005951e-01  1.119083e-02  6.653883e-02 -1.559474e-02
#>   [66] -2.091609e-01  7.777573e-02 -1.248303e-02 -2.134651e-01  2.637553e-01
#>   [71]  3.255870e-01  3.948327e-02  4.351511e-01 -2.757210e-01 -1.426827e-02
#>   [76] -7.360663e-02 -3.566144e-02  2.693224e-01 -1.175918e-01 -1.489760e-04
#>   [81]  2.490568e-01  2.109561e-01  1.989526e-01 -1.584783e-02 -5.388979e-02
#>   [86]  9.050313e-03  2.610172e-01  1.062515e-01 -6.008387e-02 -3.188824e-02
#>   [91] -9.393778e-05  1.889149e-01 -7.330815e-03 -3.290752e-02  4.510912e-02
#>   [96]  7.242406e-02  1.518144e-01  1.359595e-01 -3.521218e-02 -1.174554e-01
#>  [101]  7.568086e-02 -2.637079e-02  2.071704e-01 -2.502897e-02 -7.845821e-02
#>  [106]  3.017619e-02  3.307625e-01 -5.300310e-02  1.897551e-01  1.867168e-01
#>  [111] -2.716122e-02 -1.118821e-01 -9.614550e-02  4.570351e-02 -1.978229e-01
#>  [116]  2.527392e-02 -4.243429e-02  1.931129e-02 -1.560066e-01 -1.092736e-01
#>  [121]  2.115839e-02 -6.296446e-02  3.251163e-02 -1.507489e-03 -7.103690e-02
#>  [126]  8.229533e-02  4.822388e-02 -2.766467e-01 -1.913862e-01  6.686047e-02
#>  [131]  3.795227e-02 -6.475684e-02 -1.721587e-01 -2.115762e-02  2.730788e-01
#>  [136] -4.733778e-02  8.372531e-02 -1.639277e-02  2.099320e-02  6.623479e-02
#>  [141]  1.082697e-01  9.262771e-02 -9.365381e-02 -1.694494e-01  3.758781e-02
#>  [146] -5.686448e-02 -1.212771e-01  3.471121e-02  1.974125e-02  1.305744e-01
#>  [151]  1.240372e-01 -1.148498e-01  7.687591e-02 -2.413515e-01 -7.182990e-03
#>  [156] -6.265639e-02  6.913593e-02  1.254975e-01  7.119505e-02 -7.246704e-02
#>  [161] -9.590805e-02 -1.047324e-01 -5.892366e-03 -4.535852e-02 -5.942755e-02
#>  [166]  2.361654e-01  2.176527e-01 -1.569218e-01  4.875149e-02 -8.005901e-03
#>  [171]  1.005445e-01 -8.428027e-02 -1.686628e-01 -1.174203e-01  5.371305e-02
#>  [176] -2.338862e-02  1.011512e-01  1.178870e-01 -2.380994e-02  1.372733e-01
#>  [181] -4.090100e-01 -1.336312e-01  1.958197e-02  2.086953e-01  3.005794e-02
#>  [186] -2.091287e-02 -2.241031e-01 -1.980161e-01  1.766057e-01  2.356654e-01
#>  [191]  1.533185e-01  2.214984e-02 -8.868110e-02 -7.050129e-02  6.255147e-02
#>  [196]  7.884939e-02 -6.865824e-03 -1.376954e-02 -1.049829e-01  2.782270e-01
#>  [201] -1.585624e-01  2.637424e-02  1.354754e-01 -2.363045e-01 -9.741786e-02
#>  [206] -1.485042e-03  1.562560e-03 -5.325877e-02 -4.281093e-02 -2.856309e-01
#>  [211] -2.996860e-02  3.486771e-02 -4.491217e-02 -6.796974e-02 -1.728373e-02
#>  [216] -2.181471e-02 -9.771119e-02 -1.797298e-01  9.345658e-02 -8.934258e-02
#>  [221] -1.628589e-01  2.583017e-01 -2.506554e-02 -1.604081e-01  1.350439e-01
#>  [226] -2.414943e-02  2.451983e-01  4.418792e-02 -1.929748e-01  8.452366e-02
#>  [231] -8.896228e-03 -8.571598e-02 -1.906907e-02  4.119564e-01 -5.807879e-02
#>  [236] -6.390261e-02  2.078578e-02 -3.058428e-02 -1.265480e-01  1.225343e-01
#>  [241] -3.694305e-02  3.393035e-03 -4.885488e-02 -2.132173e-01  1.711805e-01
#>  [246]  1.912347e-01 -9.059487e-02  1.467646e-01  3.816217e-02  1.639455e-01
#>  [251] -2.073653e-01 -1.003698e-01 -1.408977e-01  1.234700e-02 -8.178751e-02
#>  [256]  8.584928e-02 -2.088753e-01 -1.487178e-01 -1.923880e-01 -2.754206e-02
#>  [261] -3.534901e-01  1.205232e-01 -1.284768e-01 -1.101423e-01 -1.084734e-02
#>  [266] -1.414896e-01 -1.280712e-01 -5.066992e-03  1.938896e-02  1.163227e-01
#>  [271]  9.475263e-04 -1.766009e-01 -3.584851e-02 -1.129344e-01  1.616790e-01
#>  [276]  9.063051e-02  2.363579e-01 -1.690898e-02  1.231658e-01  1.031577e-01
#>  [281]  9.627708e-02  1.644931e-01 -2.596380e-02  1.592941e-01  9.756062e-02
#>  [286]  2.154469e-02  1.477202e-01  1.794796e-01 -7.670551e-02  1.035831e-01
#>  [291] -3.524846e-02  6.927738e-02 -3.744635e-02  9.881257e-02 -2.054402e-01
#>  [296] -2.710155e-01  9.218835e-02 -1.445197e-01  1.069939e-01 -7.857582e-02
#>  [301]  1.551117e-02  4.290734e-02 -7.075096e-02  1.489267e-03 -1.581349e-01
#>  [306]  2.882114e-02  1.376355e-01 -8.386088e-02 -1.278941e-01  6.464292e-02
#>  [311]  8.219497e-02 -4.880123e-02  8.089760e-02 -2.535516e-02  2.023107e-01
#>  [316] -1.071873e-01  7.387509e-02  1.780235e-01  5.277488e-02  6.155685e-02
#>  [321]  1.258447e-01 -3.946571e-02 -1.311387e-01  1.406282e-01 -2.693645e-01
#>  [326] -5.222550e-02  2.581444e-02 -1.284128e-01 -5.355335e-02 -1.555897e-01
#>  [331] -1.590871e-01  1.675510e-01 -6.657140e-02  1.846635e-02 -8.028876e-02
#>  [336] -1.261143e-02  3.371148e-02  9.987191e-02 -8.905319e-02  7.686636e-02
#>  [341]  9.975662e-03  2.240595e-02  1.725463e-01 -1.235884e-02  2.768010e-01
#>  [346]  9.756407e-02  1.341030e-01 -2.383648e-01 -2.077285e-02 -2.212164e-01
#>  [351]  3.348874e-01  3.127536e-01 -9.888245e-02 -1.632747e-01  9.489244e-03
#>  [356]  1.254755e-02 -1.750328e-01  1.501349e-01 -2.466747e-02 -1.060724e-01
#>  [361] -4.108157e-02  1.834319e-01 -7.616010e-03  1.067691e-01 -1.622722e-01
#>  [366]  1.052460e-01 -9.249197e-02  9.808923e-02 -8.750445e-02 -1.847571e-01
#>  [371]  4.742534e-02  1.299813e-01  3.985334e-02 -1.047435e-01 -2.010516e-01
#>  [376] -3.772581e-02 -9.742486e-03  2.728672e-02 -6.454961e-02  1.119770e-01
#>  [381] -1.168833e-01  9.581805e-02 -1.084310e-01 -2.237480e-01  4.529565e-01
#>  [386] -1.221209e-02 -3.120090e-01  8.189399e-02 -7.206466e-04  1.768573e-01
#>  [391]  1.825874e-01  1.636447e-02  8.280664e-02 -1.082459e-01  6.530078e-02
#>  [396] -1.228291e-02  2.170246e-02 -6.993283e-02 -4.447688e-02  1.878654e-02
#>  [401] -1.332909e-01 -2.227429e-01 -1.260931e-01  2.184792e-01  3.401402e-01
#>  [406]  3.014708e-02 -4.588605e-02 -8.285039e-02 -1.608621e-02 -9.834273e-02
#>  [411]  9.728181e-02  7.715042e-02 -2.646345e-01 -1.726126e-01  1.473698e-01
#>  [416] -1.608424e-01 -3.653640e-02 -1.160981e-01 -6.470853e-02  5.164461e-02
#>  [421]  1.371267e-01  1.542814e-01  1.534781e-01 -2.265248e-02  1.558289e-01
#>  [426]  5.328590e-02  4.731865e-04  1.149183e-01  1.022404e-01  2.175093e-01
#>  [431]  1.011754e-01 -4.848854e-03 -5.994764e-02 -1.323802e-01 -2.015178e-01
#>  [436]  9.619076e-02  2.586873e-01  1.267085e-01 -4.884524e-02 -2.464330e-02
#>  [441] -2.623018e-01  3.677847e-02  7.796549e-02 -2.238507e-01 -5.043843e-02
#>  [446] -3.358968e-01  2.435602e-02 -5.734485e-03  1.244235e-01  1.251203e-01
#>  [451] -1.821484e-01  6.504390e-02 -4.249786e-02 -4.549978e-02 -5.951741e-02
#>  [456] -2.688347e-02  1.218825e-01 -6.630861e-02  7.403267e-02 -1.516505e-01
#>  [461]  8.828480e-02 -6.386871e-02 -1.948239e-01  1.428617e-01 -8.922530e-02
#>  [466]  1.784747e-01  1.895055e-01  1.249418e-01 -2.301796e-02 -2.133566e-02
#>  [471] -6.069507e-02  1.243855e-01 -5.928182e-02  3.354863e-01  1.920616e-01
#>  [476] -1.737351e-02  1.795256e-01  1.527366e-01  1.968845e-02  4.490047e-02
#>  [481] -5.417530e-03  1.165092e-01  1.117643e-01 -2.673477e-01 -5.574731e-02
#>  [486]  7.103491e-02  4.485253e-02  9.597094e-02  4.896742e-02  4.167129e-02
#>  [491] -2.186531e-01 -8.955742e-02  2.000103e-01 -1.009735e-01 -2.625410e-01
#>  [496]  7.615486e-02  5.724087e-03  7.232695e-03  2.166261e-01  4.689517e-02
#>  [501] -1.250891e-01  1.646659e-01 -2.070603e-01  3.543584e-02  3.405762e-01
#>  [506]  1.086864e-03 -1.675561e-01  8.782809e-03 -2.329670e-01  2.854212e-02
#>  [511] -2.861347e-01 -5.826084e-03 -2.776496e-01  5.082361e-02 -9.740158e-02
#>  [516]  6.720680e-02  1.838277e-01  1.110964e-01  1.315281e-01 -3.031815e-02
#>  [521]  1.348537e-01 -1.235179e-01  7.014517e-03 -2.011215e-01 -1.071841e-01
#>  [526]  6.109377e-03 -1.103404e-01  7.203077e-02 -3.826188e-02  1.746130e-01
#>  [531] -7.130909e-02  1.676940e-01 -6.251222e-02 -4.703890e-02  1.197263e-02
#>  [536]  4.342202e-02  2.688023e-02 -5.370094e-02 -9.611027e-02 -5.480797e-02
#>  [541] -9.837690e-02  2.527908e-01 -1.554956e-01 -5.222150e-02 -7.273249e-02
#>  [546]  8.464230e-03  1.502576e-01  2.794918e-03 -1.273795e-02 -5.123319e-02
#>  [551]  5.014441e-02  9.776478e-02 -3.393389e-02  2.954668e-01 -2.486953e-02
#>  [556] -5.896775e-02  7.451002e-02  4.631327e-02  1.694386e-01 -3.050444e-01
#>  [561]  1.611407e-01  1.639496e-01  5.210665e-02 -6.390603e-02 -1.617238e-01
#>  [566] -1.037797e-01  7.542800e-02 -1.262598e-01  1.005774e-02 -5.427772e-02
#>  [571] -1.052440e-01 -7.394586e-02 -1.600281e-01  5.999550e-02 -1.900534e-01
#>  [576]  2.375022e-01  1.116114e-01  2.610904e-03  2.187049e-01 -5.464627e-02
#>  [581]  5.474200e-02  2.780647e-03  1.636530e-03 -1.065854e-01  1.636957e-01
#>  [586] -1.513551e-03  1.314745e-01  9.977886e-02  7.700985e-02 -4.929118e-02
#>  [591]  1.450575e-01  1.161261e-01  1.262595e-01 -9.686811e-02  5.975767e-02
#>  [596] -8.277777e-02 -2.652017e-01 -1.664600e-02 -1.989200e-01 -1.344732e-02
#>  [601] -9.828556e-02 -1.747887e-01  5.333452e-02  4.890548e-04 -9.126443e-02
#>  [606] -8.825138e-02  5.461361e-02 -4.913002e-02 -3.780517e-02 -6.642266e-02
#>  [611]  9.420388e-02 -7.155377e-02  4.508131e-02  3.202646e-02  1.850990e-01
#>  [616] -1.267547e-01  1.668398e-01 -1.134782e-01 -5.889720e-02 -1.221105e-01
#>  [621] -8.166899e-02  2.693527e-01 -6.483406e-02  8.709273e-02 -6.248682e-02
#>  [626]  1.834504e-02 -9.978107e-02 -1.779438e-01 -1.951965e-01 -7.746610e-02
#>  [631]  1.817476e-03  4.983120e-02  2.306011e-02  2.625602e-01  7.575040e-03
#>  [636]  2.881584e-01  1.097414e-01  1.376646e-01  9.306773e-02  4.954867e-02
#>  [641]  9.118246e-02  1.601263e-01 -3.516232e-02  2.651825e-02  1.413569e-02
#>  [646] -9.087570e-02 -1.871052e-02 -4.213979e-02 -4.850392e-02  7.055205e-02
#>  [651] -6.206815e-02  5.156434e-02 -1.330937e-01  4.473334e-02  1.350621e-01
#>  [656] -1.226436e-01 -1.849968e-01 -3.471163e-01 -5.028537e-02  1.934350e-01
#>  [661]  3.129622e-02  1.097344e-02  1.665384e-01  6.854050e-02 -3.728279e-02
#>  [666] -2.947745e-01  1.585340e-01  1.385530e-01 -4.259660e-02  2.386502e-02
#>  [671]  1.811730e-01  1.043312e-01  5.044898e-02  6.598616e-02  8.463945e-02
#>  [676] -5.647066e-02  1.324570e-01  1.308166e-01  8.401841e-02 -3.011405e-02
#>  [681] -3.962021e-02  2.173895e-01 -9.217254e-02  2.275159e-01  5.084599e-02
#>  [686] -4.048546e-02  1.693280e-01 -3.232620e-02  1.312810e-02 -1.190370e-01
#>  [691]  2.945923e-01  8.947551e-02 -1.812913e-01  6.768772e-02  6.738338e-03
#>  [696] -6.184086e-02 -2.480657e-01  1.216637e-01  2.799384e-02 -1.231407e-02
#>  [701]  7.910924e-02  1.627499e-01 -1.284008e-01 -2.784950e-02 -6.067876e-02
#>  [706]  1.678812e-01 -1.908485e-01  3.273875e-02 -2.099820e-01 -1.217512e-03
#>  [711] -6.410561e-02 -2.082034e-02  1.209319e-01 -2.338834e-02 -8.987143e-03
#>  [716]  1.981700e-01 -5.078405e-02 -3.114475e-02  2.581654e-04 -1.431989e-02
#>  [721]  7.039766e-03  1.210867e-01 -2.003694e-02 -9.269061e-03  2.975918e-01
#>  [726] -2.347102e-01  3.469338e-02 -1.082515e-01  2.245993e-01  7.339741e-02
#>  [731] -9.246466e-03  3.109028e-01  1.683567e-02  1.121851e-01 -1.649257e-02
#>  [736] -8.092548e-03  2.413997e-02  4.602713e-02  1.573123e-01  4.611314e-02
#>  [741] -8.980100e-02  5.891415e-02 -5.539721e-02  1.317131e-01  4.090369e-02
#>  [746] -1.345076e-01  2.049867e-01 -4.829224e-02  1.564566e-01  1.101936e-01
#>  [751]  9.685172e-02  2.208354e-01  1.072008e-02 -3.166558e-02 -5.479006e-02
#>  [756]  1.724546e-01  3.551440e-02 -2.103081e-01  1.124877e-01  1.380916e-02
#>  [761] -6.494286e-02 -2.642729e-02  2.944464e-02  1.220366e-01 -2.171245e-01
#>  [766]  4.432983e-02 -2.035675e-02 -5.015430e-02 -9.188575e-03 -1.342277e-01
#>  [771] -1.433743e-01 -6.437942e-02  6.695309e-02  5.047075e-02 -2.251792e-01
#>  [776] -7.585638e-02  2.565071e-02  4.063268e-02 -6.601057e-03  1.843364e-02
#>  [781] -3.016479e-01  1.982717e-01  1.059400e-01 -5.401374e-02 -1.726802e-02
#>  [786]  1.491244e-01 -3.240797e-01 -3.326830e-02 -2.372858e-02  1.082297e-01
#>  [791]  2.311299e-01 -2.140589e-01 -2.690307e-01 -3.766272e-02  1.136417e-01
#>  [796]  1.713112e-01  1.573465e-01 -3.330488e-02 -3.933848e-02 -2.365645e-01
#>  [801] -1.300934e-01 -8.600768e-02  1.824741e-01 -5.240554e-02 -7.670665e-02
#>  [806]  4.029141e-02 -9.667630e-02 -1.400179e-01 -1.030566e-01  1.418162e-01
#>  [811] -2.074904e-01 -4.075535e-02 -3.112518e-03 -1.880413e-01  2.203205e-01
#>  [816] -2.052688e-01 -2.044434e-02 -1.367114e-01 -9.577305e-02 -9.366823e-03
#>  [821]  1.723664e-01  6.361542e-02 -3.309552e-01 -2.136839e-01  4.633574e-02
#>  [826] -1.146426e-01 -1.232018e-01  5.711857e-02 -1.224516e-01  7.411559e-02
#>  [831]  1.736587e-01 -5.643805e-02 -4.168635e-02 -1.660727e-01 -5.762803e-02
#>  [836] -1.273821e-01 -2.682376e-02  1.075882e-01 -1.016739e-01 -1.426574e-01
#>  [841]  7.364595e-02 -7.409269e-03 -1.686818e-01 -6.213326e-02  6.459029e-02
#>  [846]  4.837927e-02  1.783702e-02  2.990551e-01 -2.272573e-02 -2.832450e-01
#>  [851] -5.163442e-02  3.261086e-02 -5.135870e-02 -1.957029e-01 -3.269476e-02
#>  [856] -1.084346e-01 -1.795331e-01 -7.095730e-02  1.766685e-02  1.896651e-01
#>  [861]  1.307347e-01  1.094435e-01 -8.284517e-02  6.767548e-02 -2.414787e-01
#>  [866] -1.585700e-01 -1.827701e-02 -1.306113e-01  3.828216e-02 -1.795412e-01
#>  [871]  1.729815e-01 -1.180069e-01  3.303437e-02 -1.709560e-01 -1.224025e-01
#>  [876] -4.561469e-02 -2.518032e-02 -4.262232e-02 -9.766156e-03  1.593049e-01
#>  [881] -1.029805e-01  9.402963e-02  1.993025e-02 -7.582997e-03 -8.891213e-02
#>  [886] -5.872078e-02  3.590687e-02 -1.016112e-01  1.270607e-03  4.381564e-02
#>  [891] -1.499467e-01 -9.983345e-02 -9.203861e-02 -2.781808e-02  3.412062e-02
#>  [896]  1.110052e-01  1.877829e-02 -1.606297e-01  3.129504e-01  1.324593e-01
#>  [901] -1.274910e-01  1.277251e-01 -1.009500e-01 -1.061140e-01  9.405142e-02
#>  [906] -4.163056e-02 -1.195467e-01  5.935101e-02 -1.377730e-02 -2.607142e-02
#>  [911] -1.233828e-01  1.123163e-01 -6.382842e-02 -1.338154e-01  7.952106e-02
#>  [916] -2.610212e-02  1.038677e-01  3.580119e-03  7.408188e-02 -1.006215e-01
#>  [921]  4.902184e-02  1.108717e-01 -2.666992e-02 -2.328839e-01 -2.089757e-01
#>  [926]  5.464791e-02  7.412821e-02 -2.309891e-01 -1.311345e-01 -9.615351e-02
#>  [931] -2.112350e-01 -4.677097e-02 -3.984144e-02 -1.852178e-02 -2.001385e-02
#>  [936] -7.890423e-02 -1.094898e-02  2.345966e-02 -6.056418e-02  1.437712e-01
#>  [941] -8.606538e-02  3.104008e-02 -8.237898e-02  6.821333e-02 -6.439220e-02
#>  [946]  2.506947e-01  3.400962e-02  7.294635e-02 -6.100485e-02  1.048374e-01
#>  [951] -8.474228e-02 -2.657114e-01  9.880956e-02  5.484467e-02 -2.100681e-01
#>  [956]  8.136696e-02  5.825378e-02 -9.963939e-02 -1.855115e-01  1.580727e-01
#>  [961] -2.806216e-02 -4.584585e-02  2.149084e-01 -5.985987e-02 -8.859371e-02
#>  [966] -2.250827e-01 -2.769838e-01 -1.055817e-01  4.155238e-02 -2.213317e-02
#>  [971] -1.012140e-01 -9.980397e-02 -1.220309e-01  2.501881e-01 -1.532911e-02
#>  [976] -1.668374e-01  2.604023e-01  1.187544e-01 -2.554454e-01  4.668300e-02
#>  [981] -6.270116e-02  1.607057e-01 -1.252866e-01 -1.797984e-01  1.783571e-02
#>  [986]  3.273606e-02 -8.696160e-02  1.938939e-01 -2.170219e-01  8.845710e-02
#>  [991]  5.111748e-02 -9.441877e-02 -2.112050e-01 -1.052447e-01  2.910565e-02
#>  [996] -2.478494e-01 -1.822189e-01  1.613296e-01 -2.040798e-01  7.084627e-02
#> [1001]  8.545612e-02 -2.881320e-01 -9.023628e-02  5.462140e-02 -2.243625e-01
#> [1006] -6.476092e-02 -1.349788e-01 -2.108329e-01  9.518295e-02  1.846485e-01
#> [1011] -1.292570e-01 -1.571528e-01 -2.321895e-01 -7.504966e-02 -6.617543e-02
#> [1016] -2.175512e-01 -1.904427e-01  2.687279e-01 -8.474863e-02 -9.698181e-02
#> [1021] -9.174661e-03 -2.110773e-02 -9.430139e-02 -1.545799e-01 -9.651553e-02
#> [1026]  1.043089e-01 -1.571183e-03 -6.076554e-02 -3.897797e-02 -7.607885e-02
#> [1031] -1.442425e-01  3.042797e-01 -1.776420e-01 -7.394085e-02  3.128414e-01
#> [1036] -4.183004e-02  1.711570e-01 -7.689630e-02 -2.861947e-01 -9.493261e-02
#> [1041] -3.479979e-01 -8.495752e-02 -1.877908e-01  2.992807e-02  1.305163e-01
#> [1046]  9.786791e-02 -1.649214e-01  3.049802e-02 -2.125244e-01 -9.088122e-02
#> [1051]  5.399102e-02 -1.414084e-01 -6.337384e-02  1.088464e-02 -9.170389e-03
#> [1056] -1.460076e-01 -5.356687e-02 -6.445319e-02 -7.272002e-02 -5.045290e-02
#> [1061] -4.496339e-02 -2.821077e-02  1.823841e-01  1.023354e-01 -6.371034e-02
#> [1066] -7.117435e-02  5.817381e-02  4.745815e-02 -5.054877e-02  6.357947e-02
#> [1071] -1.500074e-01 -8.299020e-02 -2.177229e-01  2.603573e-01 -6.239629e-02
#> [1076]  2.745506e-01  1.621339e-01  7.872378e-02  5.126504e-03 -8.149352e-02
#> [1081]  1.023897e-01 -9.150418e-02  1.179196e-01 -1.733377e-01  4.123694e-02
#> [1086] -1.166078e-01 -2.066594e-01  3.400676e-01  1.007703e-01 -3.939823e-01
#> [1091]  8.064253e-03 -4.090972e-02  1.149117e-01 -9.574342e-03 -5.346698e-02
#> [1096] -1.497679e-01  3.059785e-02  4.482740e-02 -2.533948e-01  4.167409e-02
#> [1101]  2.586246e-02  4.260320e-02 -9.908441e-02 -3.432385e-01 -2.993673e-01
#> [1106] -7.819558e-02 -1.124218e-01  7.071484e-02  1.068142e-01  4.651238e-01
#> [1111] -2.448641e-01  1.629286e-01 -1.535981e-01 -2.043379e-01 -9.432549e-02
#> [1116] -2.037202e-01  1.043207e-01  6.834071e-02 -7.533497e-02 -3.188923e-02
#> [1121] -2.274968e-01 -9.376385e-02 -5.434395e-02  1.058204e-01 -5.378868e-02
#> [1126]  1.058569e-01  1.183779e-01  1.764033e-01  2.058596e-02  2.543456e-02
#> [1131]  1.692471e-01 -9.935029e-02  4.128978e-02 -1.365447e-01 -2.829335e-01
#> [1136]  2.954741e-02 -3.939382e-02 -1.071441e-01 -9.245012e-02 -1.192318e-01
#> [1141] -2.507950e-01 -2.469316e-02 -6.720040e-02  9.508040e-02 -2.334553e-01
#> [1146]  7.649729e-02 -6.369717e-02 -6.859376e-02  1.525789e-01  2.669869e-02
#> [1151]  1.367803e-01 -1.119882e-01  2.727496e-02 -1.059897e-01  4.576247e-02
#> [1156] -8.516334e-03 -2.141261e-01 -5.488693e-02  1.678161e-01 -1.218958e-01
#> [1161]  3.449555e-01  3.767815e-02 -1.165850e-01 -3.822489e-02 -7.142431e-02
#> [1166] -1.619661e-01 -6.753758e-03  1.171129e-01  1.176726e-01  2.383487e-01
#> [1171] -1.998536e-02  9.366462e-02  3.573655e-05  5.364676e-02  1.486833e-01
#> [1176] -1.135754e-01 -5.517623e-02 -1.375544e-01  4.237529e-02  2.268203e-02
#> [1181]  5.987321e-02 -1.829986e-01  5.113467e-02  2.040473e-02  9.170713e-03
#> [1186] -1.621988e-01 -1.110078e-01 -4.459837e-03 -1.057289e-01  5.414659e-02
#> [1191]  6.862402e-02  9.736360e-03 -1.266945e-01  2.310220e-02 -2.412110e-01
#> [1196]  1.586231e-02  8.476243e-02  7.701429e-02  1.809468e-01  5.408746e-02
#> [1201] -1.722821e-02 -5.613044e-02 -6.845716e-02  1.065527e-01  1.542879e-01
#> [1206]  3.217508e-03 -8.911931e-02 -6.304477e-02  1.432576e-01 -3.207015e-02
#> [1211] -2.068719e-01  1.232720e-01  2.394505e-03 -3.798500e-02  4.271880e-02
#> [1216] -7.608333e-02 -1.732253e-01 -1.847416e-01 -1.409817e-01  4.880434e-03
#> [1221] -2.850333e-02  1.189599e-01 -8.447832e-03  1.636760e-01  2.226202e-02
#> [1226]  1.866737e-01 -1.877577e-01  3.968913e-02 -2.686023e-01 -1.037624e-01
#> [1231] -1.707859e-01 -2.898729e-01 -1.990153e-01 -7.591940e-02  2.699476e-01
#> [1236] -4.862941e-02 -5.121774e-02 -8.195443e-02  5.185151e-02  1.234306e-01
#> [1241] -1.131140e-01 -2.598453e-04  9.659030e-02 -2.403416e-01 -2.342400e-01
#> [1246]  8.928340e-02  6.429637e-02 -1.306835e-01  1.852496e-01 -1.387821e-01
#> [1251]  5.168421e-02 -1.131787e-01  2.650921e-01 -3.947606e-02 -2.086985e-01
#> [1256]  2.322524e-01  1.388951e-01 -2.081976e-01  1.298335e-03  2.827401e-01
#> [1261]  1.834992e-01 -1.875504e-01 -1.024367e-01  1.754431e-01  4.198014e-02
#> [1266] -5.309993e-02  2.322563e-01  2.089696e-01  2.077676e-01  2.140138e-01
#> [1271]  1.009879e-01  1.442453e-01 -3.118035e-02  1.142917e-01  1.112236e-01
#> [1276]  5.523488e-02  1.150689e-01  1.364274e-02 -1.028570e-01 -1.290005e-01
#> [1281] -2.708144e-01  1.402504e-01 -2.167185e-01  6.293548e-02 -1.546657e-02
#> [1286]  3.942810e-02  8.687756e-02 -2.013771e-01 -3.124035e-02  1.810063e-01
#> [1291] -6.984927e-02  1.373997e-01 -1.609381e-01 -1.588093e-01  1.277723e-02
#> [1296] -2.593109e-02 -1.157916e-02  2.485590e-01  2.457157e-01  5.374618e-02
#> [1301]  1.593608e-01 -1.404330e-02  2.868989e-01  7.918293e-03 -7.905483e-02
#> [1306]  8.711253e-02  5.971949e-02  3.544556e-02  1.152953e-01 -8.664257e-02
#> [1311] -2.316044e-01 -2.011048e-01  2.909306e-02 -2.106411e-01 -1.286984e-01
#> [1316]  3.431548e-02  1.495170e-01  1.673839e-02 -1.113951e-01 -2.243787e-03
#> [1321] -1.294788e-02  1.869663e-01  8.091924e-02 -2.796600e-04  1.481994e-01
#> [1326]  9.780021e-02  1.445744e-01 -7.496647e-02  2.673877e-02  4.633218e-02
#> [1331]  6.001004e-02 -9.121017e-03  7.226009e-02  1.661310e-01 -3.835834e-01
#> [1336] -1.076837e-01  3.921602e-02 -6.980946e-02  1.086621e-01 -1.889510e-01
#> [1341]  2.202962e-01 -1.080854e-01  2.387233e-02 -1.154736e-01  3.131311e-02
#> [1346]  7.934236e-02  6.149840e-02 -4.042839e-02 -1.391911e-01 -1.789477e-01
#> [1351]  7.115871e-02 -7.175524e-02  1.822471e-02  3.941965e-03  2.281483e-01
#> [1356] -1.982202e-02  1.715521e-01 -2.202580e-01 -2.509816e-01 -7.985296e-02
#> [1361]  8.283063e-02  9.157255e-02  8.426293e-03  8.512261e-02  1.791405e-01
#> [1366]  3.735232e-03 -1.426802e-01 -2.755571e-02  1.789017e-01  8.574122e-02
#> [1371] -4.222582e-02  1.996561e-02 -3.965693e-02  5.312194e-02 -8.301005e-02
#> [1376]  3.764090e-02  1.041567e-02 -3.339821e-02  7.037705e-02 -5.235184e-02
#> [1381]  2.833286e-01 -3.604138e-01 -2.244673e-01 -1.473566e-02 -2.647184e-01
#> [1386] -9.680457e-03 -1.535222e-01  1.232953e-01 -2.275793e-01 -1.534580e-01
#> [1391] -2.962196e-01  2.050917e-02 -1.226237e-01 -9.699922e-02 -9.432444e-02
#> [1396]  8.498544e-02 -2.015819e-01 -2.736821e-01  2.875244e-01 -3.583366e-02
#> [1401]  1.226314e-01  2.251203e-02  1.398013e-01 -9.818853e-02 -7.177080e-02
#> [1406]  7.760676e-02  1.949794e-02 -7.320967e-03 -8.181545e-04 -9.590699e-02
#> [1411]  2.502123e-01  5.075427e-02  1.173983e-01 -3.421701e-02 -7.328304e-02
#> [1416] -4.984711e-02  4.948231e-02 -4.852086e-02  7.899335e-02  1.023302e-01
#> [1421]  1.885749e-01  1.119312e-01 -8.078219e-02 -3.301355e-02  7.257946e-02
#> [1426]  9.197149e-03  1.812038e-01  1.832716e-01  3.382429e-01  7.285987e-02
#> [1431]  7.650880e-02  1.719908e-01 -8.870823e-02 -9.611448e-02  9.300433e-02
#> [1436]  2.409321e-01  9.788290e-02 -2.371277e-01 -2.184143e-02 -9.797455e-02
#> [1441] -1.010491e-01 -1.752485e-02 -6.314772e-02  1.143018e-01  1.671217e-01
#> [1446] -3.706483e-02 -1.585926e-02 -1.522018e-01  1.603834e-01  9.693358e-02
#> [1451] -4.797197e-02 -3.916354e-02 -9.881212e-02 -6.602264e-02 -9.483627e-02
#> [1456]  1.277619e-01 -1.616874e-01 -2.248900e-01  4.309418e-02 -2.548990e-01
#> [1461] -6.273585e-02  8.369670e-02 -1.935986e-01 -9.704651e-02 -8.995396e-02
#> [1466]  1.277947e-01 -5.470612e-02 -1.559046e-01 -2.046867e-02 -3.064473e-01
#> [1471]  6.168642e-02 -4.075408e-02  1.063748e-01  9.240897e-03  1.352273e-01
#> [1476]  4.005668e-02  1.054915e-01  1.919131e-02 -4.785244e-02  3.812963e-02
#> [1481] -2.421868e-02  2.692415e-01 -2.270324e-01  4.677411e-02  1.007514e-01
#> [1486]  9.865853e-02  3.251700e-02 -8.739774e-03 -1.642316e-01 -1.053944e-01
#> [1491] -2.902442e-01 -1.858519e-01 -1.301725e-01 -1.045723e-01 -1.208471e-01
#> [1496] -1.107596e-01 -3.653951e-02 -1.141009e-01 -7.460542e-03 -2.549632e-01
#> [1501]  6.349936e-02 -6.633723e-02  2.494436e-01 -9.235529e-02 -1.022735e-02
#> [1506] -7.330607e-02 -6.854597e-02 -7.574781e-02  1.034111e-01 -3.799253e-02
#> [1511] -4.094192e-02  1.607946e-01 -1.386387e-01 -5.575488e-02 -1.863398e-01
#> [1516]  1.169292e-01 -1.412722e-01 -1.235457e-01  4.309542e-03 -3.378769e-03
#> [1521]  1.868718e-01 -1.993690e-01 -7.560865e-02 -2.153721e-01 -6.171910e-02
#> [1526] -6.817958e-02  1.425637e-01  3.132184e-02  2.469235e-01  1.763936e-01
#> [1531] -5.387381e-02 -6.373033e-02  1.569495e-01  1.033736e-01 -1.032510e-02
#> [1536]  2.295076e-01  1.156094e-02 -4.632078e-02  2.013725e-01  3.367626e-01
#> [1541]  3.352614e-01  5.906232e-02  8.869528e-03  3.972654e-02 -7.250564e-02
#> [1546] -1.031293e-01 -5.718428e-02  1.772464e-02  2.247531e-01 -1.133068e-01
#> [1551] -2.872455e-03  5.691667e-02 -1.797482e-02 -2.429223e-01 -1.527176e-01
#> [1556]  2.842121e-01 -1.779015e-01 -7.221773e-02 -2.215384e-01 -1.471587e-01
#> [1561]  8.035337e-02 -1.287029e-03 -1.252183e-02  1.247118e-01  2.521556e-01
#> [1566] -3.280552e-01 -1.820325e-01  9.161927e-02 -1.332036e-01  2.376221e-01
#> [1571] -2.161169e-01  2.219452e-01  1.118450e-01 -1.634785e-01  3.323287e-02
#> [1576] -6.883174e-02  1.385525e-01  2.289388e-01 -8.448615e-02  1.647306e-01
#> [1581] -6.310210e-02 -6.910030e-02  8.854727e-02  3.250315e-02 -1.200787e-01
#> [1586] -8.578991e-02  1.170949e-01 -9.285983e-02 -9.815849e-02 -3.626727e-02
#> [1591] -8.467267e-03 -6.022836e-02  2.217232e-02  1.155334e-01  2.333489e-02
#> [1596]  8.527079e-02 -5.041890e-02  1.667262e-01  7.492689e-02  6.851367e-01
#> [1601] -5.977042e-02  1.507797e-01  2.429481e-01 -1.644332e-01 -2.905767e-02
#> [1606]  5.441765e-02  4.773286e-02  4.137347e-02  1.415948e-01  1.455163e-01
#> [1611] -4.982935e-02  1.332927e-01 -1.064611e-02  7.957713e-02 -1.595421e-01
#> [1616] -4.534448e-02 -9.130559e-02  4.206788e-02 -1.855431e-01 -2.257193e-01
#> [1621]  6.610268e-02  8.053699e-02 -1.000067e-02  1.408444e-01  2.567170e-01
#> [1626] -8.868415e-02  1.826929e-02 -6.709631e-02  5.698690e-02  2.243958e-01
#> [1631] -1.551463e-01 -1.100247e-01  1.210371e-01  1.835168e-02  1.300296e-02
#> [1636]  4.754848e-02 -8.107438e-02 -2.330667e-01 -2.260296e-01 -1.285910e-01
#> [1641]  1.082275e-01 -4.402086e-02 -3.179493e-02  7.658345e-02 -1.002444e-01
#> [1646]  1.145629e-01 -5.976426e-02 -3.556740e-02  2.579294e-02 -2.451671e-03
#> [1651] -3.625089e-02 -9.387033e-02  9.894809e-02  5.771887e-02 -1.749699e-01
#> [1656] -8.261316e-02 -2.351372e-01  1.851829e-01 -1.344652e-01 -3.359371e-02
#> [1661]  1.837610e-02 -1.888234e-02  3.004365e-02 -1.423105e-01 -1.055028e-02
#> [1666]  9.899245e-02 -3.476216e-02  2.386582e-02 -1.847720e-01 -2.398291e-02
#> [1671]  6.884560e-02  7.585966e-03  2.740777e-02  2.295989e-01  1.547870e-01
#> [1676] -1.099451e-01 -1.634580e-01  6.777629e-02 -1.071082e-01 -4.487277e-02
#> [1681]  3.353862e-01 -9.985393e-02 -4.894856e-02  1.033400e-01 -9.014327e-02
#> [1686] -5.123981e-02 -1.735114e-01  6.581644e-04  8.606120e-02 -2.534229e-01
#> [1691]  2.082556e-01 -3.161524e-02  4.147295e-02 -2.070049e-01  1.994191e-01
#> [1696]  7.288373e-02 -2.354478e-01 -5.815304e-02  2.254643e-01  1.595773e-03
#> [1701]  1.166670e-01  1.983985e-01  6.841161e-02  1.587035e-01  6.719137e-02
#> [1706] -6.990447e-02 -9.069602e-03  1.159488e-01 -2.114384e-02 -1.374098e-01
#> [1711] -5.248755e-02  1.215136e-01  2.158644e-03  5.258588e-02 -7.358181e-02
#> [1716] -1.753189e-01 -1.538167e-01 -1.904862e-01 -3.312089e-02  1.918507e-02
#> [1721] -3.969817e-02 -1.079295e-01 -1.791806e-01 -1.128127e-01  6.186459e-02
#> [1726] -1.054218e-01  3.894518e-02 -1.157742e-01  1.571467e-01  5.849585e-02
#> [1731] -4.449045e-02  1.411787e-01 -3.141871e-02  1.095050e-01 -1.037985e-01
#> [1736]  1.110985e-01 -8.921625e-02  1.932059e-01  1.436217e-01 -6.542366e-02
#> [1741] -1.760503e-01 -1.928606e-02 -1.281050e-01  1.067774e-01  1.710316e-01
#> [1746]  1.320302e-01 -2.676377e-02 -3.156545e-01  4.107857e-02 -2.522180e-01
#> [1751]  6.271634e-02 -1.325865e-02  2.177542e-01  1.506022e-01  6.828004e-02
#> [1756] -1.200863e-02 -8.843719e-04  1.468753e-02 -3.491598e-02  4.225870e-01
#> [1761] -1.503813e-01  2.022558e-01  8.340259e-02 -4.581307e-02  2.904192e-01
#> [1766]  3.483817e-01  3.020506e-01  2.599823e-02 -4.394525e-02  1.711441e-01
#> [1771] -1.021892e-01 -1.522400e-01 -9.514603e-03 -1.751039e-01 -1.723163e-01
#> [1776] -6.673167e-02 -2.883467e-01  2.605588e-01 -1.885272e-01  7.349976e-02
#> [1781]  1.141635e-01 -3.649311e-01 -1.471964e-01  2.563475e-02 -1.560968e-01
#> [1786] -6.483206e-02 -7.950074e-02 -2.894330e-02  1.974416e-01 -4.838458e-02
#> [1791] -1.615901e-01 -2.912414e-01 -1.276234e-01 -1.526339e-02 -1.778857e-01
#> [1796]  7.991366e-02  3.828509e-05 -2.404304e-02 -1.166271e-01  5.902528e-02
#> [1801]  8.656016e-02 -1.743675e-01 -2.684494e-01  5.627647e-02  5.467559e-02
#> [1806]  1.399491e-01 -1.703642e-01 -2.135904e-01 -2.188333e-02  6.329684e-02
#> [1811]  9.363333e-02  1.625073e-01  2.926534e-02  3.155856e-02  2.835689e-01
#> [1816]  2.390421e-01 -9.100087e-02 -1.504269e-02  1.476616e-01 -3.239779e-01
#> [1821]  1.310605e-01 -3.186449e-01  3.087102e-02 -6.293811e-02  9.790326e-02
#> [1826]  2.470198e-01 -1.282957e-01 -1.888195e-01 -1.773297e-02  7.479379e-02
#> [1831] -1.439874e-03 -2.977275e-01  8.458490e-02 -7.117015e-02 -2.609165e-01
#> [1836]  1.700852e-02 -1.112075e-01 -9.228473e-02 -2.932551e-01  1.285622e-02
#> [1841]  1.417248e-02 -1.435731e-01 -5.395088e-02 -2.256442e-01  9.206191e-02
#> [1846] -1.947649e-01  1.820429e-01  7.358796e-02  2.405480e-02  4.392894e-02
#> [1851]  3.973372e-02  1.411249e-01 -8.497645e-02  5.892350e-02 -1.276675e-01
#> [1856]  1.225190e-01  1.069817e-01  4.444686e-04  1.486998e-01  1.543223e-02
#> [1861]  2.254587e-01 -1.462829e-01 -7.116048e-02 -1.440491e-01  1.833579e-01
#> [1866]  1.375462e-01  1.120689e-02 -2.056089e-01  1.008507e-01  5.589079e-02
#> [1871]  5.269678e-02 -4.736994e-02 -1.487867e-02 -1.438906e-01 -1.086144e-01
#> [1876]  1.262967e-01  4.279344e-02  1.728270e-01  1.494276e-01  2.281470e-02
#> [1881]  1.269737e-01 -4.447161e-02  7.156305e-02 -9.103415e-03 -1.251001e-01
#> [1886] -4.089724e-02 -1.034216e-01  5.700648e-02  2.830477e-01  5.265189e-02
#> [1891]  1.811121e-01  1.512494e-01  1.293507e-01 -3.359592e-01 -2.693953e-02
#> [1896]  2.365296e-01  5.147356e-02  1.671086e-01 -1.804595e-01 -1.911478e-01
#> [1901]  6.637954e-03 -8.581694e-02  1.790354e-01  3.386512e-01 -2.093812e-01
#> [1906] -4.378328e-02  1.481818e-01  9.861106e-02  3.223678e-02  1.797207e-01
#> [1911] -5.881375e-02  1.954622e-01 -1.301776e-02 -2.449428e-01  6.426916e-02
#> [1916]  1.188604e-02  1.006139e-01  2.011652e-03 -8.280035e-02 -7.289025e-03
#> [1921]  4.774154e-02  2.576384e-02 -1.194699e-01 -3.115888e-01 -9.246415e-02
#> [1926] -4.269335e-02  1.947041e-01 -6.486547e-02  1.903289e-01  7.418949e-02
#> [1931] -1.732709e-01 -9.089451e-02  1.507677e-01  6.242992e-02  1.589994e-01
#> [1936] -1.608152e-01 -1.219853e-01 -7.501871e-02  1.096336e-01 -2.660369e-02
#> [1941] -2.960590e-01  5.943177e-03 -7.250156e-02  6.013472e-04 -1.280957e-01
#> [1946]  2.564497e-02  1.419683e-01 -2.033947e-02 -8.525240e-02  2.003731e-01
#> [1951]  6.334482e-02 -2.868249e-02  1.108180e-01 -4.224121e-02 -1.329482e-01
#> [1956] -1.047857e-01  1.607015e-01  1.795464e-02 -3.672726e-02 -6.186463e-02
#> [1961] -2.057718e-01 -1.362483e-01 -2.476247e-02 -4.426156e-02 -1.987256e-01
#> [1966]  4.747541e-02  4.565374e-02 -7.317972e-02 -9.251089e-02 -1.763323e-01
#> [1971]  2.596843e-01  2.956178e-01  8.562759e-02  4.247221e-02 -1.849685e-01
#> [1976]  1.986554e-01 -5.816683e-02  1.891137e-02  2.050316e-01 -1.117351e-01
#> [1981]  1.056065e-01 -1.028190e-01 -4.131765e-02 -1.310434e-02 -2.939799e-01
#> [1986]  7.134931e-02  3.400751e-01 -1.324315e-02 -2.820093e-02 -2.900343e-01
#> [1991]  1.011459e-01 -1.965840e-01  1.385975e-01 -1.310072e-01  1.427247e-01
#> [1996] -4.033154e-03 -1.055749e-01  5.381855e-02  1.043173e-01 -2.468646e-02
#> [2001] -1.811067e-01  1.147988e-01  1.434825e-01  8.553145e-02  2.616945e-01
#> [2006]  7.061113e-02  1.168888e-01  2.171705e-01  1.113506e-01 -8.947171e-02
#> [2011]  1.814339e-01  5.232091e-02  2.248584e-01  9.062369e-02  3.208830e-02
#> [2016] -9.363385e-02  1.899928e-01 -4.121030e-03  7.848929e-02 -1.547428e-01
#> [2021] -5.094681e-02 -6.593051e-02 -1.899094e-01  3.162594e-02  1.833073e-01
#> [2026] -1.616187e-01 -9.730427e-02 -1.698335e-01 -6.358555e-03  6.346329e-02
#> [2031]  4.692788e-02  1.469679e-01 -1.833793e-02 -1.047277e-02 -1.101073e-01
#> [2036] -6.169901e-02 -8.983670e-02  1.401587e-01  1.742480e-01 -1.225185e-01
#> [2041]  1.666021e-01 -1.928280e-01  8.455220e-02  6.433197e-02 -1.474185e-01
#> [2046]  1.968206e-01 -4.333783e-02  1.394890e-01  1.990486e-02 -2.319831e-01
#> [2051]  5.564604e-02  1.970899e-01 -1.798303e-01 -7.253798e-02  8.966404e-03
#> [2056]  1.006580e-01  7.258361e-02  7.621613e-02 -1.354032e-01 -6.606061e-02
#> [2061]  1.355050e-01  5.060357e-02  3.253061e-02 -2.770213e-01  2.537314e-01
#> [2066] -3.139009e-02  2.329773e-02 -1.150806e-01 -1.053728e-01  1.389416e-01
#> [2071]  3.494883e-02  1.067180e-01  2.575911e-01 -2.447067e-02  7.313836e-03
#> [2076]  2.259726e-01  1.079175e-01 -7.082056e-02 -4.883079e-02  3.254181e-02
#> [2081]  1.571555e-02 -5.864404e-02  1.015327e-02 -2.108226e-01  9.208168e-02
#> [2086] -6.836360e-02 -1.615463e-01  2.312364e-01  3.033625e-01 -4.698547e-02
#> [2091]  1.252719e-01 -1.732576e-01  3.560625e-02  7.978510e-02  1.514381e-01
#> [2096]  0.000000e+00 -1.976216e-01 -3.409281e-01 -1.347547e-01 -9.225978e-02
#> [2101] -1.009286e-01 -9.571441e-02 -1.336521e-01 -1.177482e-01  2.712484e-03
#> [2106]  4.543851e-02 -1.169856e-01  9.631268e-02 -1.266318e-01  1.185892e-01
#> [2111]  5.491218e-02 -1.933313e-02  1.162833e-02 -7.829373e-02  6.140830e-02
#> [2116] -1.271487e-01 -1.610237e-01  7.757637e-03  7.983256e-02 -1.786893e-02
#> [2121]  6.458019e-02 -1.508686e-02  7.889860e-02  1.423243e-01 -1.270038e-02
#> [2126]  1.907520e-02  9.317861e-02 -1.083414e-01  6.451600e-02 -5.822946e-02
#> [2131]  7.894861e-02 -6.413440e-02 -3.842268e-02 -2.151884e-01  1.836970e-01
#> [2136] -1.823696e-01  3.063297e-02 -1.490773e-01  2.377200e-01 -1.829139e-02
#> [2141]  1.614155e-01 -1.247703e-01 -1.252056e-01 -1.635639e-01  1.800314e-01
#> [2146] -2.464400e-02 -3.505166e-02 -5.437940e-03 -1.971993e-01  3.791729e-02
#> [2151] -5.812569e-02 -1.622090e-01 -7.973714e-02 -2.194139e-01 -7.004829e-02
#> [2156]  3.571771e-02 -1.114247e-01 -2.646408e-01 -2.399591e-01  1.572136e-01
#> [2161]  1.317845e-01 -3.774196e-02  1.004606e-01 -7.676624e-02 -8.534792e-02
#> [2166]  7.722504e-02  5.702910e-03 -7.412352e-02  1.076319e-01  2.920442e-02
#> [2171]  4.658280e-03 -1.397527e-02  1.179805e-01 -1.691503e-01  1.050696e-01
#> [2176] -2.530753e-02 -2.211929e-01 -1.611589e-02 -1.418297e-01  1.342532e-01
#> [2181] -8.592497e-03 -1.053628e-01  1.260942e-01 -9.918071e-02  2.050194e-02
#> [2186]  4.827562e-02  7.439143e-02  6.714835e-02  6.882471e-02 -1.146902e-01
#> [2191]  1.293947e-01  8.298785e-02  1.465264e-02  3.335052e-01 -5.660991e-02
#> [2196]  2.317006e-01  1.225027e-01 -9.707929e-02  1.582846e-01  1.441442e-01
#> [2201]  2.906617e-03  1.245775e-01 -4.524745e-02  3.606458e-01  1.809278e-01
#> [2206] -7.836317e-03  2.570347e-01  1.637261e-01  7.799898e-02 -7.916738e-02
#> [2211]  1.792236e-01  4.425566e-02  1.607910e-01 -9.621636e-04  3.099224e-02
#> [2216]  1.007641e-02  1.209175e-01  1.845723e-01  2.451349e-01 -4.669550e-02
#> [2221] -3.349810e-02  5.665327e-02 -6.559640e-02  2.744777e-02  2.981361e-02
#> [2226]  9.826791e-02 -6.712053e-02  1.669343e-01 -2.288160e-01 -3.200453e-02
#> [2231]  1.858608e-02 -2.100722e-01 -6.468709e-02 -1.565017e-01 -1.333922e-01
#> [2236] -1.458662e-01  1.935449e-01 -1.556634e-01  2.433240e-02 -1.234865e-01
#> [2241] -4.265128e-02  2.825887e-01  1.411983e-01 -2.453233e-01 -2.540881e-04
#> [2246]  2.513728e-02 -1.548035e-01  2.381624e-01 -1.089510e-01 -4.675370e-02
#> [2251]  8.104960e-02  1.175800e-02 -9.035423e-02 -2.001925e-01  8.311795e-02
#> [2256]  2.390373e-01  1.239441e-01 -4.611885e-04  4.133102e-02 -7.004931e-02
#> [2261] -3.279314e-01 -1.154948e-01 -1.449081e-01 -4.877570e-02 -2.412056e-01
#> [2266]  7.691657e-02 -1.753732e-01 -1.566579e-01  2.310369e-02 -1.029083e-01
#> [2271] -1.417530e-01  1.753495e-01  4.075177e-02  1.117154e-01 -2.527961e-01
#> [2276]  4.314049e-02  2.183791e-01 -2.828244e-02 -8.514048e-02  5.483972e-02
#> [2281]  1.085667e-02 -6.166055e-02 -1.393304e-01 -5.482694e-02 -5.835466e-02
#> [2286]  3.601923e-02  5.699073e-02 -3.002413e-01 -1.705938e-01 -1.585878e-01
#> [2291]  1.556610e-01  5.290661e-02 -1.352752e-01  5.601359e-02 -1.808108e-01
#> [2296]  1.935302e-01  1.856012e-01  2.684135e-01 -2.072265e-01 -1.208529e-01
#> [2301] -1.085300e-01  3.748911e-02 -9.618331e-02  1.194020e-01  7.060546e-02
#> [2306] -1.427505e-01 -2.923212e-02 -1.150383e-01  4.369242e-02  1.306061e-02
#> [2311] -1.660480e-01 -1.178696e-01  2.622152e-01  1.644353e-01 -2.856697e-01
#> [2316]  8.472042e-02  1.103841e-01 -3.815045e-02  1.377953e-01  4.125459e-02
#> [2321]  1.298827e-01  2.051789e-01 -5.877918e-03 -5.665125e-02  2.153965e-02
#> [2326] -1.456650e-01 -1.995138e-02 -2.214003e-01  4.631930e-03 -1.761203e-01
#> [2331]  5.041516e-02 -1.652730e-01  1.269382e-01  1.381971e-01 -1.390624e-01
#> [2336] -7.403686e-02 -1.360812e-01  7.959206e-02  1.189346e-01 -1.279599e-01
#> [2341]  1.467035e-01  1.222801e-02  1.168560e-01  1.803199e-01 -6.650364e-02
#> [2346]  5.437592e-02  3.722420e-03 -1.099800e-01  1.395153e-01 -2.348930e-01
#> [2351] -1.777137e-01 -1.430446e-01  1.984848e-02 -1.663743e-01 -1.114512e-01
#> [2356] -9.091461e-02  3.586320e-03 -1.884881e-01 -1.984236e-01 -2.150184e-01
#> [2361]  6.999000e-02 -1.305904e-01  5.153458e-02  1.291772e-01 -1.752965e-01
#> [2366] -2.313215e-02 -2.223775e-02  1.030666e-01  2.175914e-01  6.183749e-02
#> [2371]  1.486742e-01 -3.348593e-02 -2.023434e-01  2.568249e-01  1.160925e-03
#> [2376] -3.300934e-02  3.661266e-02 -6.127551e-02 -5.390366e-02 -1.380872e-01
#> [2381]  3.210514e-02 -6.451544e-02 -1.355594e-01  3.558730e-02  1.408380e-02
#> [2386]  4.953251e-02 -6.336878e-02  7.894293e-02  4.354970e-02 -1.874543e-01
#> [2391]  1.930193e-02  6.446273e-02  3.149442e-01  6.654700e-02  1.435036e-01
#> [2396]  2.260704e-01  1.793315e-01  4.007108e-01 -1.388763e-02  1.289761e-01
#> [2401]  1.180079e-02  3.916917e-02 -3.378333e-02  2.337812e-02  2.572880e-01
#> [2406]  1.554533e-02 -1.380121e-01  1.574794e-01 -1.537439e-01  1.865729e-01
#> [2411]  8.441579e-02 -5.903632e-02  3.091152e-01 -1.436132e-01  1.276515e-01
#> [2416] -1.848687e-01 -1.138821e-01 -1.508760e-01  4.020907e-02  1.755801e-01
#> [2421]  1.658333e-01  1.028265e-01 -6.826822e-02  2.883353e-04  3.664745e-02
#> [2426] -3.447074e-02 -7.770856e-02  1.461186e-01  6.212118e-02 -8.808713e-03
#> [2431] -1.357527e-01 -1.419179e-01 -3.667869e-01 -2.147478e-02 -2.212180e-02
#> [2436]  1.168998e-01 -1.480886e-01  3.478104e-02  1.708467e-01 -2.199108e-01
#> [2441] -9.506427e-02 -1.872154e-01 -1.756351e-01 -8.879648e-02  1.045165e-01
#> [2446] -1.014258e-01 -1.739569e-01  3.823012e-02  1.045314e-02 -1.552789e-01
#> [2451] -1.703139e-01  3.341263e-02 -9.016442e-02 -1.999705e-01 -1.560293e-01
#> [2456] -4.349714e-02 -1.023019e-01 -3.009107e-02  2.104922e-01 -5.309848e-02
#> [2461]  4.573265e-02 -3.330530e-01 -4.727697e-02 -1.569588e-01  1.634024e-01
#> [2466] -1.129613e-01 -2.000163e-01  8.065582e-02 -4.945293e-02  1.511140e-01
#> [2471] -7.591280e-02 -1.070265e-01 -1.233453e-01 -1.239627e-03 -1.360486e-01
#> [2476]  3.249142e-02 -1.224764e-01  3.262834e-02 -4.873084e-02 -2.371505e-01
#> [2481]  2.423839e-01 -1.187883e-01 -1.396127e-02 -9.734468e-02 -3.192700e-01
#> [2486] -9.024882e-02 -2.102226e-01  5.735432e-02 -8.334422e-02 -2.160926e-01
#> [2491]  1.394161e-01 -2.071530e-01 -7.824374e-02  1.625314e-01  2.487983e-02
#> [2496] -1.537967e-01  6.359276e-02 -1.215189e-01 -2.387969e-01 -1.050409e-01
#> [2501]  2.375983e-01 -7.446114e-02 -1.615170e-01 -1.065672e-01 -9.137667e-02
#> [2506]  2.011827e-03  2.027991e-03  6.441815e-02  2.345179e-01 -3.332699e-02
#> [2511] -6.838327e-02 -1.495179e-01 -2.622688e-01 -3.381871e-01  4.919761e-02
#> [2516] -9.131440e-02  4.081035e-02  1.657408e-01 -1.407316e-01 -3.205573e-02
#> [2521]  1.188463e-01  3.956906e-02 -1.253318e-01  1.815464e-02  1.603975e-01
#> [2526]  1.609111e-01 -3.314402e-02  8.884626e-02  1.215544e-01  1.036574e-01
#> [2531] -9.116829e-02 -1.661126e-01 -1.013727e-01  1.585499e-01 -8.852808e-02
#> [2536]  4.229439e-02 -1.975359e-01 -3.919796e-01 -2.052690e-02 -5.158998e-02
#> [2541]  1.784202e-01  2.295877e-01 -1.631549e-01 -1.228466e-01  7.497158e-02
#> [2546]  1.416726e-01  2.888489e-02  5.176753e-02 -1.193164e-01  2.097599e-01
#> [2551]  8.903392e-02  1.355469e-02 -3.732958e-02  1.389656e-02 -1.569242e-01
#> [2556]  1.973604e-01 -5.287500e-02  2.090530e-01 -2.608041e-01 -1.884726e-02
#> [2561] -9.001756e-02  9.375866e-02 -3.120205e-03 -6.066894e-02  1.235294e-01
#> [2566] -3.259300e-03  3.841106e-02 -1.527592e-02  4.830615e-02 -1.495329e-01
#> [2571]  2.355588e-01  1.528339e-01  2.099690e-01  2.874555e-02 -8.508133e-02
#> [2576] -1.639785e-01 -1.227237e-01  9.696224e-02 -1.577963e-02 -1.665903e-01
#> [2581] -1.364184e-01  2.058356e-01 -1.839283e-02  1.686793e-02  1.106767e-01
#> [2586]  2.287740e-01 -4.766681e-02 -9.279719e-02 -1.900518e-01  1.544402e-02
#> [2591]  1.794733e-01 -1.876349e-01  4.185067e-02 -3.528745e-02 -3.950208e-02
#> [2596]  1.609015e-01  3.348006e-01 -2.205592e-02  4.055595e-02  2.337649e-01
#> [2601] -7.045281e-02 -2.260789e-01  1.955388e-01 -1.294165e-01 -9.052214e-02
#> [2606] -1.118169e-01 -2.633141e-02  4.850590e-02  5.651337e-04 -1.387341e-01
#> [2611] -9.064473e-03 -1.788995e-03  1.047180e-01  2.179796e-01 -3.347271e-02
#> [2616] -2.869503e-02  2.048727e-02  1.507653e-01  1.337601e-01  6.879094e-03
#> [2621]  1.830872e-01  1.488778e-01 -5.599159e-02  2.065547e-01 -4.039165e-03
#> [2626]  9.142543e-02 -1.253366e-01 -4.125318e-02 -1.839695e-02 -3.799986e-02
#> [2631] -6.085824e-02 -1.912583e-01 -1.319026e-01 -1.205428e-01 -2.112384e-01
#> [2636]  1.742719e-01 -1.264763e-01 -1.178602e-01 -4.852046e-02 -1.187413e-01
#> [2641] -1.612764e-01  1.231098e-01 -5.718418e-02 -1.007283e-01 -1.494598e-01
#> [2646]  1.569122e-01  8.521565e-02  1.490931e-01 -1.579472e-01 -1.123012e-01
#> [2651] -6.977220e-02 -1.236738e-01 -9.085645e-02 -1.634263e-01  1.128489e-02
#> [2656] -1.847169e-01 -2.107044e-01 -1.908572e-01  2.950627e-02 -5.580121e-02
#> [2661] -1.315548e-01 -2.665744e-01 -8.144634e-02  1.491168e-01 -7.983088e-02
#> [2666] -1.600379e-01 -8.558904e-02  1.584041e-01  1.125662e-01  1.665353e-01
#> [2671]  8.053638e-02  1.773796e-01  7.673778e-02  1.507005e-01  2.129843e-01
#> [2676]  1.559038e-01  2.742402e-01  4.200338e-02  1.588380e-01 -1.714730e-01
#> [2681]  6.676582e-02  5.056804e-02  2.295767e-01  2.081999e-01 -4.932892e-02
#> [2686] -2.037543e-01  6.095140e-02  1.428008e-01  4.886834e-02  4.116552e-02
#> [2691]  6.856851e-02 -1.090044e-02  9.373516e-02 -2.463471e-02  2.055787e-01
#> [2696]  6.824171e-02  2.628711e-01 -2.277766e-02 -1.109026e-02  2.647563e-02
#> [2701] -1.795346e-01 -6.527126e-02  4.756645e-02  4.183179e-04  2.202478e-01
#> [2706]  2.000040e-01 -3.518986e-02  1.206720e-01 -7.030539e-02 -3.000479e-02
#> [2711] -4.911187e-02  1.517093e-01 -8.623095e-02  7.052798e-02  1.347134e-01
#> [2716]  7.277297e-02  1.227573e-01 -8.770580e-02 -3.178081e-02  1.190391e-02
#> [2721] -2.003479e-01 -1.853541e-01 -2.627198e-01  1.907236e-01  5.656091e-02
#> [2726] -2.004289e-01 -6.209992e-02 -1.793360e-01  5.858769e-02 -1.712209e-01
#> [2731] -9.411586e-02 -1.617523e-01 -2.827905e-02 -2.624954e-01 -1.307520e-01
#> [2736] -1.842073e-01 -1.770153e-01 -8.351883e-02 -2.120025e-01 -1.925937e-01
#> [2741]  1.666679e-01  2.318996e-01  5.553531e-02 -1.017913e-01  3.064259e-01
#> [2746]  2.717392e-01 -4.717349e-02  1.020953e-01 -3.996316e-02  5.878876e-02
#> [2751]  1.541618e-01  4.832126e-03 -6.529765e-02 -3.166958e-02 -1.884771e-01
#> [2756]  2.641162e-01  6.692216e-02 -2.373399e-02  2.671329e-02 -7.066395e-02
#> [2761]  4.011137e-02  1.334971e-01 -3.358780e-02 -6.303933e-03  9.024718e-03
#> [2766] -9.752604e-03  1.196699e-01 -7.986886e-02  1.061391e-01  3.629211e-02
#> [2771] -3.872026e-02  3.276373e-01  2.686595e-02  2.129716e-01  2.039033e-01
#> [2776] -2.773818e-01  2.249939e-01  1.021203e-01 -5.041857e-03  1.561891e-01
#> [2781] -2.241792e-01 -2.995958e-02 -8.575244e-02  7.815445e-02  2.300016e-01
#> [2786]  5.328556e-02 -2.551788e-01 -3.456077e-02 -1.261565e-02 -1.887802e-01
#> [2791]  1.112633e-01 -1.556827e-01 -3.081424e-02  4.096702e-02 -1.537843e-01
#> [2796]  4.091360e-02  1.107824e-01 -3.343399e-02  1.642448e-03 -4.985789e-03
#> [2801] -2.153590e-01 -1.645170e-01 -4.366568e-02 -4.085761e-02 -1.292471e-01
#> [2806]  1.316432e-01 -3.046788e-02 -6.548597e-02  5.859641e-02  4.523237e-03
#> [2811]  3.556012e-02 -1.255833e-01  1.509606e-01  4.293090e-02  6.230282e-03
#> [2816] -3.660270e-02 -6.857773e-02  4.947681e-02 -4.056843e-02  1.235139e-01
#> [2821]  1.232029e-01  6.650972e-02  4.165919e-02  7.686846e-02 -1.808957e-02
#> [2826]  1.045666e-01  8.504596e-03 -7.182527e-03 -9.140103e-02  1.420720e-01
#> [2831]  8.530716e-02 -1.196299e-01  1.243297e-01  2.933110e-01  2.389813e-01
#> [2836]  9.057345e-02 -1.524745e-01 -1.320787e-02  1.212747e-01 -5.335556e-02
#> [2841]  1.466880e-01  1.282518e-01 -1.283059e-01  1.179674e-01  2.773416e-01
#> [2846] -3.434760e-02 -9.022995e-02 -4.858483e-02 -1.352827e-01  1.288266e-02
#> [2851] -5.571733e-02  3.297932e-02  1.821756e-01 -3.435230e-01  8.278214e-02
#> [2856] -1.225035e-01 -8.504281e-02 -8.284673e-02 -4.428528e-03  2.605766e-01
#> [2861] -7.327677e-02 -1.319819e-01 -2.237156e-01 -6.887347e-02  3.077847e-02
#> [2866] -1.095746e-01  2.094032e-01  4.106003e-02  3.045189e-02 -2.473350e-01
#> [2871] -6.094827e-03  5.334478e-02  1.275684e-01 -4.181634e-03 -2.781984e-03
#> [2876] -1.356285e-01  4.339977e-02 -5.643692e-02  1.263621e-02 -1.165646e-01
#> [2881]  1.025664e-01  1.632282e-03 -1.008473e-01 -1.106666e-01 -7.376446e-02
#> [2886]  1.536150e-01 -1.761109e-01 -3.068406e-02  1.136761e-02 -1.963399e-02
#> [2891]  6.122177e-02  1.825268e-01  9.364387e-02 -4.353492e-02  1.330200e-01
#> [2896]  3.992387e-02  2.165573e-02  1.101523e-02  4.595491e-02 -2.176145e-01
#> [2901] -3.172152e-02 -6.002266e-02 -1.591730e-01  1.211420e-01 -5.673995e-02
#> [2906]  1.547014e-01 -1.221894e-01 -1.675708e-01  3.281040e-01  2.073530e-01
#> [2911] -2.345039e-02  1.369937e-01  1.975390e-01 -1.139276e-01  1.562508e-01
#> [2916]  4.154077e-02  6.764543e-02  2.355118e-02 -5.337982e-02  1.868427e-02
#> [2921]  9.694917e-03  8.802541e-02 -1.499793e-02  1.611515e-02  4.814064e-02
#> [2926]  2.123643e-01 -1.252635e-01 -2.163573e-02 -3.385467e-02 -1.334389e-01
#> [2931] -7.972961e-02 -9.305634e-02  3.402915e-02 -9.640297e-03 -1.314744e-01
#> [2936]  7.395900e-02  1.544509e-02 -4.069784e-01 -1.108541e-01  3.435482e-03
#> [2941] -2.917750e-02 -3.518833e-01  3.800572e-02  4.180131e-02 -1.558345e-02
#> [2946]  1.246048e-02 -9.376109e-02 -1.612533e-01  1.658128e-01  1.550894e-01
#> [2951] -1.243291e-01  1.126306e-01 -1.260103e-01 -3.127835e-02 -1.340046e-01
#> [2956] -1.413186e-01  7.711266e-02 -6.176647e-02  8.686717e-02  2.475279e-01
#> [2961] -6.673171e-02 -2.005017e-01 -5.880007e-02 -1.858061e-01  6.791943e-02
#> [2966] -1.260229e-01  5.054507e-02  7.575529e-02  7.699514e-02  9.001413e-02
#> [2971]  2.837513e-01  3.530921e-02  7.946785e-02 -1.381859e-01  5.119995e-02
#> [2976]  2.584962e-02 -5.882949e-02  5.006078e-02 -4.035367e-03  9.756057e-02
#> [2981] -5.682646e-02 -1.430073e-02 -2.097380e-01 -1.714680e-01 -1.583258e-01
#> [2986] -6.501761e-02  1.565865e-01 -5.979196e-02 -2.354554e-01  2.653064e-01
#> [2991] -2.623824e-02  1.284604e-01 -1.296125e-01  1.707829e-02  1.266916e-01
#> [2996]  1.183629e-01  2.458463e-01  1.148287e-01 -1.487964e-01 -4.330310e-03
#> [3001]  1.508034e-01 -1.778314e-03  9.180256e-02  2.114085e-02  1.299527e-01
#> [3006] -1.063232e-01 -5.284076e-02 -1.037121e-01 -1.694843e-01  7.589856e-02
#> [3011] -2.278057e-02 -2.735131e-01 -1.093846e-01 -1.279426e-01 -4.473023e-02
#> [3016]  1.493196e-02 -5.021603e-02 -7.423607e-02  3.255589e-02  1.144783e-01
#> [3021]  4.236575e-02 -1.024603e-01  8.943502e-02  8.464642e-02  7.896966e-02
#> [3026] -4.069996e-02 -1.398056e-01  2.515822e-02 -1.666292e-01  1.348240e-01
#> [3031]  6.945069e-02  7.921330e-02 -9.048351e-02 -1.152511e-02 -1.155052e-02
#> [3036] -7.569247e-02 -3.019574e-02 -8.173274e-02 -6.958883e-02 -8.873652e-02
#> [3041] -9.570930e-02 -2.414346e-01 -1.526261e-01 -1.090241e-01  1.482702e-02
#> [3046] -8.724911e-02 -1.484180e-01  9.459482e-02  2.266151e-01 -5.065034e-02
#> [3051]  2.639331e-01 -2.959594e-02  1.684294e-01 -2.367954e-01 -5.291468e-02
#> [3056] -3.504714e-02 -9.680542e-03  6.906774e-02  3.935564e-01 -9.801749e-02
#> [3061] -1.407602e-02  1.041288e-01  1.018732e-01  3.467341e-02  7.438859e-03
#> [3066] -1.657543e-01 -6.974124e-02  4.758715e-02 -1.569117e-01 -1.189716e-01
#> [3071] -1.838596e-01 -7.213114e-02  1.752740e-01  4.922375e-02 -1.524739e-02
#> [3076]  1.846448e-02 -6.873841e-02  5.261794e-02  9.351779e-02  1.660627e-01
#> [3081]  8.637601e-02  9.770515e-02  1.543312e-01  1.844760e-02  1.560159e-03
#> [3086] -5.499724e-02 -2.381092e-01  1.442862e-01  5.259026e-02 -1.883676e-02
#> [3091]  1.220679e-01 -1.802921e-01  4.451391e-02 -5.261202e-02 -1.318802e-01
#> [3096]  2.200041e-02  2.196434e-01  3.333336e-02 -1.015616e-01 -1.339326e-02
#> [3101] -1.726078e-01  7.332819e-02  1.418624e-01  3.799256e-02  4.563943e-02
#> [3106]  6.243252e-02  8.689091e-02  1.103489e-01  5.026051e-02  4.330379e-02
#> [3111]  4.854754e-02  3.800481e-03 -2.309946e-02  3.827617e-02  1.216213e-01
#> [3116]  2.206548e-02 -8.769524e-03 -6.074575e-02  1.036623e-01  4.060395e-02
#> [3121] -1.991056e-02 -7.170762e-02  1.543631e-03  1.064624e-01 -1.715279e-01
#> [3126]  8.855050e-02 -4.758625e-02 -8.694734e-02  4.578954e-03  6.914215e-02
#> [3131] -1.472357e-01 -1.050381e-01  8.927138e-02 -1.736399e-01 -2.184328e-01
#> [3136] -2.228574e-01  1.127067e-01  6.067616e-03 -8.014381e-02  9.471985e-02
#> [3141]  5.073236e-02  2.573364e-01 -8.407721e-02  4.550918e-02 -8.386523e-02
#> [3146]  4.829922e-02  5.052470e-02 -4.585859e-02  8.742587e-03  3.235539e-02
#> [3151]  1.132839e-01 -1.983235e-01  6.509222e-02  9.366469e-02 -4.626151e-01
#> [3156] -1.426967e-01  2.517981e-01  4.745341e-02 -5.414245e-02  2.219924e-01
#> 
#> $mod2_unimodal$ref_properties$cell$pca_dist
#>    [1] 27.85083 30.95346 32.04039 31.70145 32.34540 29.31993 24.69458 31.94323
#>    [9] 32.64213 28.63502 27.46478 28.49470 32.41435 26.20480 29.13963 30.71185
#>   [17] 27.63892 32.74534 32.50445 24.52701 28.78376 31.22369 29.85001 21.51332
#>   [25] 29.65532 29.91148 20.40388 29.85599 29.92296 27.49250 24.04842 31.29883
#>   [33] 30.14816 20.22824 29.86534 31.21312 25.12790 30.80586 25.56300 34.69573
#>   [41] 23.34028 28.73964 30.80233 28.41558 26.63596 31.08258 32.40651 26.55531
#>   [49] 30.34597 28.20163 28.29985 23.74393 25.50614 27.55908 28.06095 19.88950
#>   [57] 30.40416 29.24428 24.07443 30.26786 28.66672 27.19802 28.90864 23.26488
#>   [65] 30.50263 27.69470 29.48846 29.27938 27.79020 35.28621 31.79393 29.08740
#>   [73] 31.04509 33.18887 33.91782 31.87610 32.07151 28.98260 27.31891 32.28069
#>   [81] 26.96625 25.57770 33.00077 26.64027 25.57041 29.44491 28.20696 26.80874
#>   [89] 28.68380 29.42933 28.00828 28.22407 21.25122 28.43514 29.29022 26.74153
#>   [97] 33.82174 30.26289 27.72899 31.02510 32.13830 27.54852 27.02948 28.75870
#>  [105] 28.34036 26.60541 31.02050 27.15461 29.39720 30.75323 23.00428 24.92689
#>  [113] 34.81243 27.21169 26.60185 28.85853 31.11717 23.66612 26.06648 27.52543
#>  [121] 24.85814 29.27131 29.95606 22.80086 25.82464 26.52425 32.49653 33.07203
#>  [129] 25.30021 29.23661 30.98280 29.49016 28.56498 28.36804 26.34738 33.60543
#>  [137] 26.18739 31.40191 34.29010 25.82307 30.29415 30.23471 28.69298 29.06355
#>  [145] 23.94689 27.69167 31.34236 34.37650 32.76265 24.18980 30.20111 27.69769
#>  [153] 25.42892 28.91677 29.12964 28.09362 24.78196 32.19353 31.98566 29.68713
#>  [161] 33.00823 28.16608 24.80920 29.44503 26.46195 31.40665 35.95759 19.86942
#>  [169] 26.30174 32.72162 27.97739 27.80687 28.01365 31.63174 31.61645 27.32412
#>  [177] 17.83611 28.75609 25.52685 30.97320 30.68697 32.86207 30.70863 27.81131
#>  [185] 28.11662 28.34188 26.86179 28.81707 31.46785 26.95751 27.54823 30.83756
#>  [193] 28.33104 23.95076 29.90047 34.65783 28.20839 30.48398 29.17878 32.74988
#>  [201] 32.06652 33.60384 24.28860 27.14992 30.25121 28.67415 29.77522 25.61024
#>  [209] 31.81761 29.19130 29.22684 29.99001 30.92685 27.69225 29.90585 27.63042
#>  [217] 34.56935 29.05902 23.30604 30.36944 29.34382 30.14369 30.12705 26.24797
#>  [225] 34.71871 29.24507 31.77606 29.24031 30.85001 28.07552 23.88360 29.42523
#>  [233] 30.92393 29.58453 29.24365 23.08462 25.89822 32.09022 34.48288 30.55850
#>  [241] 27.89087 27.18128 26.46726 26.47156 31.15873 27.59289 31.92985 29.36145
#>  [249] 27.84645 31.85671 32.98301 32.58301 32.49997 37.10861 29.94661 30.78108
#>  [257] 33.51383 32.45825 28.12247 33.45976 33.04019 33.27724 31.26908 27.64729
#>  [265] 32.26250 34.63531 25.74557 30.76146 35.02198 26.44255 31.98431 34.21600
#>  [273] 35.13411 30.68895 28.69980 31.54111 35.94714 31.71752 32.14062 30.70415
#>  [281] 36.58395 34.47115 30.94241 36.22406 27.33167 31.58356 30.88799 14.96526
#>  [289] 36.34120 33.51447 33.05042 32.44770 31.81574 27.99194 33.64263 33.67280
#>  [297] 33.97551 29.47441 31.16215 33.26742 32.54516 32.61626 27.97203 35.66509
#>  [305] 27.71925 31.62072 36.34254 33.48845 29.02309 29.46802 28.25513 30.38784
#>  [313] 30.65850 27.51572 29.58817 27.78128 28.23222 23.23204 22.07851 29.11660
#>  [321] 27.17821 29.54965 29.06903 26.23242 32.84220 27.37442 26.55453 33.13226
#>  [329] 33.08294 25.09653 32.07053 32.74779 30.23335 31.77164 27.19858 31.71114
#>  [337] 32.16652 33.76952 30.73774 29.82348 32.11025 29.59380 25.80634 30.43517
#>  [345] 33.83923 30.88655 29.02327 33.43756 30.13632 33.61556 29.70442 30.36521
#>  [353] 32.23010 28.78122 29.58184 31.69030 32.79525 31.25305 32.27649 23.07124
#>  [361] 32.70866 32.15023 30.72798 26.13598 30.61879 35.56481 35.23804 24.11369
#>  [369] 31.67483 30.12762 26.44671 28.64226 29.14049 31.59614 32.82110 27.43437
#>  [377] 36.63125 31.34938 30.91377 28.59775 25.02972 27.94959 34.67762 31.00295
#>  [385] 33.48571 31.61365 28.65261 32.61834 29.06883 24.78548 28.27879 29.18293
#>  [393] 29.81294 28.33016 32.62834 29.84834 27.10758 25.66138 30.15854 27.91115
#>  [401] 31.79050 31.47216 35.50200 27.72735 28.91041 34.12995 35.10624 35.31786
#>  [409] 28.89412 36.89070 33.87968 29.44715 32.68553 28.38658 35.07512 19.99156
#>  [417] 27.91889 35.02616 24.91739 28.47222 33.09804 31.78094 32.97705 35.83359
#>  [425] 33.89710 37.30162 36.82956 31.11017 35.19646 36.40458 37.66377 31.39243
#>  [433] 34.47437 27.26369 35.11520 31.15095 30.24647 33.45712 34.25078 32.56927
#>  [441] 27.71101 29.79945 31.55786 34.23450 33.10315 30.42095 28.64128 31.84266
#>  [449] 33.83231 30.87312 32.14478 30.70426 33.04282 30.81130 34.63824 33.45281
#>  [457] 33.23717 33.52425 33.86388 30.93407 34.57008 32.26629 31.37865 23.11772
#>  [465] 22.77493 31.29402 25.38057 33.55096 28.25416 33.12431 29.38452 28.83766
#>  [473] 30.22098 30.14659 31.97262 33.60555 32.26998 28.38303 27.65466 27.23566
#>  [481] 22.32219 36.17066 30.07831 29.03760 35.48371 30.93498 32.29076 34.27213
#>  [489] 30.19200 30.19386 30.93624 31.55724 37.03301 32.96353 32.40341 32.17265
#>  [497] 30.38918 33.98781 33.88652 24.98241 35.89477 29.03809 37.88380 26.79080
#>  [505] 27.84335 36.13678 34.82483 32.45402 33.11985 28.78372 32.23206 34.03892
#>  [513] 36.94439 32.50577 30.55324 33.63174 35.56019 33.33207 35.36429 28.59791
#>  [521] 31.97988 30.80065 29.78557 36.29877 31.68656 35.94575 30.08910 33.91069
#>  [529] 32.08480 30.07497 25.09456 28.42031 33.19324 31.77601 25.92580 32.33304
#>  [537] 32.50686 27.78718 32.08810 30.87480 34.53097 33.44973 25.41453 30.04118
#>  [545] 28.75509 33.30149 32.48117 38.08896 26.89164 31.78820 32.25666 28.16210
#>  [553] 32.73941 32.10789 29.39436 31.52443 28.98081 29.90276 28.60520 31.41414
#>  [561] 24.68697 26.74536 31.84006 26.20833 31.84833 28.00399 27.33911 30.20368
#>  [569] 29.79619 24.75513 31.26102 27.05424 24.45754 34.02603 35.47423 34.10505
#>  [577] 26.84668 30.24027 29.01626 33.51987 27.84728 34.38863 31.37692 27.73392
#>  [585] 26.23430 30.92692 30.64659 34.23213 29.89405 26.97445 28.37809 22.06496
#>  [593] 28.87238 33.95869 33.78804 33.43545 34.66833 25.87686 32.13848 31.45770
#>  [601] 32.04079 30.05498 31.40568 29.08420 31.09104 30.55745 27.54137 30.96438
#>  [609] 28.02894 30.93750 29.53023 34.58565 29.21788 32.26927 30.82753 32.80589
#>  [617] 22.77037 26.42741 27.51177 28.09136 34.91320 31.89456 32.28175 36.64483
#>  [625] 29.95342 34.28257 21.06777 23.60594 32.38712 30.09355 32.79244 31.31103
#>  [633] 30.55391 28.09307 26.41615 29.92329 31.97922 31.37197 27.67481 34.31499
#>  [641] 34.59548 27.05706 38.59549 25.77472 30.35621 31.56380 30.77595 36.52371
#>  [649] 32.44751 30.44353 30.40654 25.27884 33.72620 30.44090 28.83413 34.77864
#>  [657] 35.29048 33.09729 33.01667 24.04132 25.75598 32.62975 31.56376 29.44223
#>  [665] 31.63825 32.81266 26.39975 29.37050 32.44306 30.17295 27.19627 28.19725
#>  [673] 32.17692 33.14291 28.98064 24.10999 27.91655 31.31372 30.88044 27.73342
#>  [681] 32.98364 30.66951 31.42916 30.24958 39.57451 21.82026 29.31849 30.12635
#>  [689] 32.71911 29.41163 26.49074 33.50721 31.22873 26.34687 36.10208 29.02247
#>  [697] 28.93917 30.10685 30.27810 29.79978 32.62208 29.53999 31.13115 23.96538
#>  [705] 25.15276 25.22648 23.31687 23.98098 28.03657 29.59610 30.38728 27.68693
#>  [713] 30.04285 27.02583 32.48344 22.49820 31.28995 28.84172 26.49426 29.74245
#>  [721] 31.41983 32.93823 28.03014 28.62156 29.16419 28.45022 26.59164 25.45579
#>  [729] 28.95616 31.95110 28.40014 30.40081 30.74339 32.34181 29.39928 29.34767
#>  [737] 27.57813 28.58201 31.86942 29.47790 23.59270 24.92493 28.94475 33.38546
#>  [745] 27.33674 34.55907 31.00673 27.66829 25.89551 26.37166 26.76129 29.34392
#>  [753] 26.12282 32.74773 29.11637 31.74979 26.72045 32.62771 25.68477 26.30949
#>  [761] 32.70809 21.67713 30.37879 28.09122 35.25973 29.34704 32.71863 31.67188
#>  [769] 26.60189 29.86914 28.35068 25.87094 28.44018 22.04807 19.44531 22.81572
#>  [777] 26.58744 30.75830 24.91134 22.68378 26.47492 30.91340 27.22625 29.57020
#>  [785] 29.59358 29.16499 27.87304 30.53955 28.75116 33.93779 28.69028 26.48356
#>  [793] 30.13558 30.49675 24.83673 26.74770 29.19841 27.97639 26.14062 26.02048
#>  [801] 25.75392 27.27823 28.81161 25.17574 31.48507 31.15969 26.41876 29.43083
#>  [809] 23.51958 27.63761 27.83622 28.06828 30.57581 29.00535 23.50988 31.79497
#>  [817] 24.55270 29.79519 31.97996 29.83887 30.38182 32.78334 25.27798 24.75301
#>  [825] 31.60544 30.90162 30.24335 30.57938 29.46511 32.29781 27.69308 34.79022
#>  [833] 28.38844 24.64919 30.29740 32.65228 30.59388 32.07993 27.42931 28.30305
#>  [841] 30.46785 26.29444 36.20389 33.24477 30.33775 27.02865 35.68412 37.22508
#>  [849] 32.77278 30.32449 30.12853 30.45668 29.76910 29.64969 31.53273 32.41999
#>  [857] 29.87823 30.46396 26.14453 30.18150 30.00675 29.60583 38.09937 30.67989
#>  [865] 24.72509 33.26624 33.63509 32.83299 30.79495 29.54239 33.85190 33.56512
#>  [873] 30.80426 29.63683 29.93091 31.85279 31.06214 32.71754 31.83354 37.66683
#>  [881] 33.25381 30.80157 28.62846 27.69655 30.94841 33.40016 28.89996 23.47919
#>  [889] 26.34167 27.95784 29.09340 36.16637 35.61566 33.24160 25.71168 25.76041
#>  [897] 27.44292 25.69447 32.42056 31.00758 29.67496 35.30512 31.25687 26.72988
#>  [905] 31.27622 33.21154 29.32635 26.20762 29.21465 35.87892 29.55801 34.01303
#>  [913] 31.51578 33.73534 32.23280 31.68221 33.34055 32.50495 32.59618 33.94645
#>  [921] 35.15923 28.91847 37.30108 26.99675 33.08503 32.81602 29.27718 29.49612
#>  [929] 32.08707 33.36493 31.68172 26.42975 30.79120 31.73030 28.60597 35.26107
#>  [937] 29.99040 27.32674 34.91548 29.54548 36.31399 29.40897 34.32181 29.01891
#>  [945] 32.80395 33.81619 34.28577 31.10084 31.25914 30.52164 31.59804 29.58214
#>  [953] 26.43713 22.46426 30.50938 32.86498 30.48748 25.88472 33.25572 27.89684
#>  [961] 30.34696 30.98959 32.10943 30.77338 30.57072 26.75970 30.71750 31.32865
#>  [969] 25.12006 32.35997 28.95374 29.19913 30.67079 30.56478 29.40973 31.96034
#>  [977] 32.64502 34.88486 30.55309 27.47182 29.50407 32.58607 33.58606 25.77670
#>  [985] 28.86254 32.56316 32.02585 32.18024 31.02529 25.92605 30.41101 34.97309
#>  [993] 27.76368 31.84445 28.83420 31.02344 33.77721 32.35396 27.79188 27.68502
#> [1001] 30.61872 32.38156 27.25258 28.80843 32.28289 33.37412 35.56818 33.94268
#> [1009] 28.85412 28.14316 27.81003 28.03262 30.91439 31.11788 33.12073 30.54357
#> [1017] 23.81371 29.24105 29.40490 32.36702 28.95704 27.14035 22.93408 25.32870
#> [1025] 27.35051 24.71247 30.36228 28.52887 24.65594 24.46147 24.89950 28.88109
#> [1033] 35.68926 31.30097 22.23324 32.67042 26.84619 31.51105 26.70277 25.40306
#> [1041] 26.68383 30.26554 26.74101 25.20156 30.39889 33.62905 26.92090 30.76321
#> [1049] 29.19370 28.14784 25.67938 24.06450 21.28024 28.02404 28.68272 29.63008
#> [1057] 27.65840 28.83691 29.45910 30.16259 31.67403 22.96154 27.03126 28.85098
#> [1065] 27.13397 26.84586 26.48369 26.19941 30.57274 24.09905 28.35724 31.52635
#> [1073] 19.95506 28.51281 24.06168 31.52016 25.83189 29.76167 30.00163 32.14045
#> [1081] 31.01828 30.94756 27.59225 34.70481 31.07211 27.91227 32.27019 30.67213
#> [1089] 29.05806 30.95983 29.08974 29.55808 30.82848 24.90117 35.40901 33.84906
#> [1097] 32.67465 29.84927 28.04836 27.81309 31.38399 25.15166 30.70705 30.17294
#> [1105] 25.42719 30.77467 28.15022 25.88203 30.96004 30.62986 29.34422 29.99030
#> [1113] 26.47543 30.02708 32.14098 31.20584 31.52624 32.08771 29.53389 29.04929
#> [1121] 21.64802 26.47287 34.92119 27.00123 30.24660 33.74197 24.69431 31.05682
#> [1129] 26.59216 29.68706 25.22239 32.01718 30.12149 33.18435 33.42257 27.88453
#> [1137] 25.11710 31.51398 26.59319 27.62584 31.80352 36.29774 31.07324 26.73686
#> [1145] 27.97710 30.94654 35.12334 30.61942 23.89276 31.39871 26.08192 30.47834
#> [1153] 26.82226 27.06598 30.32055 30.01404 28.21594 28.66406 30.81300 38.56430
#> [1161] 32.92813 27.71366 30.02146 27.48704 32.43558 32.08956 34.15064 27.05934
#> [1169] 30.53058 34.34835 34.43229 27.85432 32.55997 34.17233 30.96028 32.31803
#> [1177] 34.10209 32.90098 28.32427 32.33698 32.83213 34.51046 30.22773 34.97275
#> [1185] 30.86296 31.84549 32.32532 33.07140 25.52392 25.91060 34.74238 30.06076
#> [1193] 36.41979 28.83044 25.77286 32.93916 28.95452 36.08872 33.76778 28.01976
#> [1201] 34.30455 32.03029 32.73903 22.43546 32.34808 27.35350 29.85983 25.80555
#> [1209] 25.65717 29.84224 25.84493 30.15399 28.83899 26.64973 26.65659 26.58046
#> [1217] 31.12718 31.39842 32.41329 31.58843 28.23619 29.42840 29.19375 21.24786
#> [1225] 26.61439 30.78633 25.61043 30.29612 27.80981 28.25923 28.04951 28.68882
#> [1233] 29.26970 33.06473 31.75183 29.23279 34.80976 27.20318 30.90460 28.62339
#> [1241] 30.85627 30.68991 27.39042 24.64905 22.81137 30.54600 28.50452 30.34695
#> [1249] 25.15078 32.63155 27.77309 26.85834 22.35764 28.53502 36.35848 27.04232
#> [1257] 32.83537 21.00494 25.73435 38.17811 31.58039 29.89451 26.48397 32.19907
#> [1265] 30.80472 28.18730 30.24644 28.28301 32.34318 23.88585 28.77555 31.82577
#> [1273] 29.50493 34.43680 33.10748 20.22093 24.35742 30.80447 29.94566 29.10273
#> [1281] 30.51425 36.37200 32.05788 23.82650 27.24335 27.44841 30.80699 26.16455
#> [1289] 28.39775 26.77462 28.86037 30.06370 27.59680 29.11296 31.77083 29.97386
#> [1297] 31.36952 30.89464 25.07327 34.27484 33.93103 32.14745 29.95229 34.52916
#> [1305] 25.98159 30.32704 35.33669 31.04851 29.98440 25.96175 24.40482 28.55952
#> [1313] 25.83258 25.04004 32.28336 29.26524 27.96081 26.87793 23.77700 32.99831
#> [1321] 31.91941 33.14564 30.88000 27.91278 26.90191 28.27199 29.37497 30.40256
#> [1329] 30.88964 29.13043 35.24909 28.15591 35.05306 33.00569 36.13521 28.97365
#> [1337] 26.26133 29.73923 33.99798 31.26362 33.82041 33.68536 29.63297 31.30263
#> [1345] 27.30222 29.46367 30.72251 31.62433 30.55809 29.97936 35.80281 34.50303
#> [1353] 32.70590 32.23453 34.07645 28.85261 31.53229 32.38591 29.99113 29.52133
#> [1361] 37.31188 33.49837 31.05015 37.01685 29.81077 30.92933 33.70087 33.32116
#> [1369] 29.86851 32.41277 31.51744 28.62425 35.03389 29.36350 33.63714 31.14591
#> [1377] 30.02966 29.89574 22.73302 29.48434 33.98530 34.05501 34.38277 29.01553
#> [1385] 32.43734 34.25739 35.80850 29.33579 27.35879 32.35228 26.43658 32.88170
#> [1393] 31.37005 26.32619 25.26672 26.24086 26.94531 27.38687 30.64541 31.03287
#> [1401] 31.97089 25.84053 28.18655 26.91540 31.43108 29.83519 26.59162 27.26547
#> [1409] 29.19987 32.65358 22.49898 25.79024 33.11089 34.98045 30.21401 31.86014
#> [1417] 29.96225 24.17529 29.70471 27.63136 30.61771 27.13153 22.14722 29.44939
#> [1425] 30.90390 28.81043 32.32652 24.75593 23.95427 31.67555 25.97812 24.18860
#> [1433] 28.94357 31.16973 31.68251 26.87356 34.47748 35.15462 30.32930 28.97921
#> [1441] 31.03096 32.69058 33.57120 34.53893 29.29620 28.10927 32.47711 30.11257
#> [1449] 27.24536 34.04805 25.34184 29.09132 29.62346 23.51848 29.51885 27.78641
#> [1457] 25.93388 25.99522 25.43099 26.45311 31.38799 25.05510 26.57462 25.37914
#> [1465] 30.09458 23.39962 30.93135 29.91430 30.55140 29.90184 30.43017 33.03790
#> [1473] 30.68417 26.17489 28.49392 25.59414 27.42465 27.83452 28.78093 30.76656
#> [1481] 28.55700 31.17949 27.30845 26.70122 28.57231 23.71530 28.39148 34.21652
#> [1489] 24.08965 26.40946 33.18740 26.29594 29.39434 25.40731 28.39857 29.51478
#> [1497] 26.32648 34.13935 28.15214 26.81784 31.32864 34.75861 25.61526 24.79992
#> [1505] 33.59283 26.87299 24.99759 26.96458 24.84435 31.60016 34.47306 32.22545
#> [1513] 33.23194 29.65933 32.34091 31.10840 30.13332 28.24379 31.56413 34.64658
#> [1521] 33.64743 32.28733 33.71948 23.18752 28.41899 32.32353 27.41236 25.91306
#> [1529] 29.17849 35.40386 31.40869 36.43024 25.59139 28.43674 30.43217 23.26866
#> [1537] 29.60999 26.60128 29.17311 26.89410 29.16306 31.09364 33.04797 25.11322
#> [1545] 30.42345 32.06407 34.26208 24.16095 25.05465 29.19181 30.11674 33.20116
#> [1553] 31.68130 25.47525 30.42265 28.29563 27.64502 30.47074 23.29749 29.17413
#> [1561] 24.75984 26.54910 29.66049 31.26780 32.63195 28.15189 31.84869 32.05707
#> [1569] 36.16610 27.29840 30.37383 28.48404 33.26402 29.16034 25.90996 30.03512
#> [1577] 32.55484 33.02877 29.19938 33.99995 30.79716 31.38052 30.29607 26.04643
#> [1585] 28.30963 32.92209 33.31242 37.17973 29.98661 27.11154 24.96922 25.22546
#> [1593] 32.39618 28.87798 33.01019 26.04469 29.60881 35.36601 33.96737 30.73248
#> [1601] 32.25699 30.97902 29.63656 29.29457 29.89482 33.12616 28.42620 32.97457
#> [1609] 31.21839 28.13358 28.71349 24.18269 30.81548 32.96277 29.90376 32.44604
#> [1617] 28.95022 26.81765 31.14463 30.59555 35.54318 36.09157 28.64232 31.61599
#> [1625] 28.00271 31.99379 28.69619 28.58086 31.19941 26.66507 32.36426 29.75538
#> [1633] 21.91290 33.61429 26.02663 30.54484 27.57647 25.88288 29.88782 31.42145
#> [1641] 30.22889 31.65222 33.05184 26.71670 29.17797 30.14710 29.21193 24.33390
#> [1649] 22.30602 23.77875 26.95288 21.20940 27.89032 26.61813 24.58228 33.18266
#> [1657] 32.15949 27.61245 27.40433 21.70405 27.79511 26.47741 32.74297 29.98430
#> [1665] 24.25772 32.73449 31.80733 28.24524 31.21219 32.08581 25.23265 25.27183
#> [1673] 30.95284 26.83426 28.04233 23.89388 28.22245 26.57712 28.94320 28.06456
#> [1681] 33.47959 33.10489 20.73253 20.86527 31.73834 30.23500 29.16479 22.73389
#> [1689] 29.76075 28.00610 27.78296 27.49925 32.54572 32.15459 31.77332 30.00192
#> [1697] 33.72366 29.49244 33.52343 30.18770 27.75804 28.27893 32.76816 33.42861
#> [1705] 33.40488 35.93595 30.48883 33.01181 30.44441 27.72472 25.13789 28.35056
#> [1713] 25.10378 34.71244 33.18923 30.46117 29.49290 28.42489 30.15195 33.94980
#> [1721] 34.22813 28.64570 34.74189 31.55145 32.41946 32.37587 28.32436 33.64192
#> [1729] 27.87018 28.57676 24.55737 30.86675 31.13122 32.33141 30.07255 29.77744
#> [1737] 27.08934 27.61612 32.43973 33.98048 26.17135 32.03628 31.26020 29.09540
#> [1745] 32.49654 33.54004 34.14041 33.47909 28.98780 28.84175 25.74600 31.82786
#> [1753] 35.72598 32.00811 32.21507 30.92056 34.56218 27.85349 32.86643 32.19728
#> [1761] 31.91817 25.04378 31.36540 25.84454 29.65227 32.61918 28.85167 35.77940
#> [1769] 31.79972 32.29235 29.41943 33.26972 36.78799 32.91437 30.34420 32.09650
#> [1777] 32.64973 34.42516 28.00490 32.84895 25.68809 28.27742 29.35738 32.43489
#> [1785] 32.54686 26.73114 33.29570 31.98939 26.82023 30.57026 35.31996 28.21833
#> [1793] 30.22680 26.28633 32.41334 18.74972 31.32967 34.08811 25.18960 18.86805
#> [1801] 25.49651 28.41528 26.70581 26.64360 29.39333 23.41137 28.63711 30.90978
#> [1809] 30.00589 32.57530 28.34876 28.87636 24.71332 26.84227 29.44637 28.08817
#> [1817] 30.45126 27.97053 26.95430 28.21781 25.95860 29.95706 31.56825 23.40816
#> [1825] 34.09328 34.62386 30.12281 31.56434 28.68120 31.06022 31.56299 31.39368
#> [1833] 26.79791 27.91091 33.63933 25.74505 32.82990 34.28796 30.14712 29.15549
#> [1841] 29.89177 31.30521 32.95018 21.40340 30.58492 32.33306 31.43758 27.43514
#> [1849] 30.95449 31.59297 33.25840 27.43194 21.74223 31.42993 26.73992 31.09871
#> [1857] 30.68005 20.57854 27.55358 33.89093 30.10690 31.62410 22.09749 31.64987
#> [1865] 30.85045 33.46810 26.49787 30.05458 30.74894 30.26020 31.82828 33.33240
#> [1873] 33.35218 30.57067 28.04960 33.59459 35.76997 27.10807 33.31890 31.67741
#> [1881] 29.17049 29.09392 31.61485 26.86664 28.01156 32.33933 37.60139 31.99561
#> [1889] 28.92834 23.68965 32.58663 27.30717 31.50730 29.20207 32.31651 27.95311
#> [1897] 34.78000 31.40339 31.83205 34.64772 32.32761 32.61000 36.93182 27.04292
#> [1905] 29.87234 33.43684 32.89141 26.87195 30.62689 30.02743 33.06436 26.82641
#> [1913] 31.38362 31.87793 31.48498 30.28708 28.14634 31.73353 31.49889 24.61064
#> [1921] 25.84449 26.35652 31.35845 32.50972 30.31358 30.45021 32.55469 24.99773
#> [1929] 34.68387 31.03913 30.62945 30.97201 33.23275 31.54265 30.56142 30.69915
#> [1937] 34.88877 25.94626 26.26049 31.31610 28.85403 35.74380 26.44698 31.18080
#> [1945] 35.03316 31.62548 31.14144 27.68204 27.25515 30.31999 30.12596 28.26597
#> [1953] 26.73634 33.02646 25.50106 28.83493 25.05478 29.26685 29.03324 27.46292
#> [1961] 34.26382 29.34183 28.71609 28.35784 29.36776 26.15716 28.43952 30.53405
#> [1969] 30.29746 28.85560 30.34141 22.34460 29.42882 29.87733 34.84072 33.21023
#> [1977] 28.84883 29.58442 33.79771 28.15279 27.72610 30.74025 30.64963 28.22629
#> [1985] 30.90435 35.74200 28.43635 31.16691 30.63850 30.30694 28.33958 25.72811
#> [1993] 33.23678 27.77399 32.98950 30.18608 29.01337 31.99534 28.93946 35.69729
#> [2001] 31.36092 31.36251 29.02855 30.09330 34.55553 32.38647 33.44937 31.60402
#> [2009] 25.77700 32.75239 32.86082 31.95576 35.99059 33.38257 31.57868 35.46878
#> [2017] 26.87370 32.31175 34.99887 28.34198 29.93867 31.14301 36.59842 29.61794
#> [2025] 31.82908 32.22527 32.42425 34.04447 30.62285 32.34656 33.75060 30.49415
#> [2033] 26.21343 30.70860 36.24019 26.68419 32.24027 30.76860 26.77228 32.01887
#> [2041] 31.27012 29.26742 31.61032 37.00390 32.55787 29.55766 34.75162 21.04063
#> [2049] 27.92842 34.47263 32.58328 33.36341 25.53508 31.27800 30.44230 26.38579
#> [2057] 30.64624 31.67849 26.00173 33.73150 31.56322 37.28853 29.16358 30.33800
#> [2065] 31.33091 31.14285 35.49228 33.12181 30.99751 33.98998 27.68250 33.49617
#> [2073] 33.33140 35.47813 29.45315 33.72367 29.95150 30.24500 27.68123 26.39435
#> [2081] 35.86663 30.52140 29.95716 28.82305 32.11833 31.96489 31.04518 29.57027
#> [2089] 30.44047 30.84705 31.50894 32.18289 34.50206 25.71668 30.99902 33.07804
#> [2097] 26.24948 29.41149 29.99495 30.55516 29.82703 34.23473 27.43080 34.51213
#> [2105] 26.61251 29.18101 33.52996 33.01317 21.28935 33.52889 31.71407 32.08878
#> [2113] 31.30759 29.96019 31.24951 33.21857 33.42039 30.97260 36.28079 32.80456
#> [2121] 27.57506 33.90721 26.00550 30.86190 24.66798 30.44853 27.50630 27.04095
#> [2129] 23.76020 27.08157 24.53181 23.54899 24.66095 26.47043 31.03104 28.18621
#> [2137] 24.65232 34.00039 25.92994 32.57194 26.05438 29.70258 21.93113 26.97741
#> [2145] 25.46557 26.43990 25.57713 29.72335 29.84382 26.27413 26.22323 23.66572
#> [2153] 25.01977 25.32699 25.18856 27.30684 25.82166 26.92288 21.03133 26.06875
#> [2161] 33.62739 29.99210 21.40490 28.47508 30.64969 30.50519 29.53697 26.37366
#> [2169] 29.29967 27.45282 35.07701 27.00262 25.68328 29.87017 34.50514 32.34559
#> [2177] 26.97773 29.91582 27.19673 27.28998 34.84852 28.62553 29.35090 28.95735
#> [2185] 26.56997 35.85217 29.51480 30.25543 30.81365 26.26270 29.47053 32.07279
#> [2193] 29.42602 29.24521 26.90658 27.16027 24.74776 31.49042 26.62547 32.96633
#> [2201] 26.60955 33.46220 23.84308 30.01540 35.00486 27.36366 35.67123 27.01968
#> [2209] 33.20201 29.00617 27.50110 34.15470 30.63460 29.17601 32.43120 30.68953
#> [2217] 29.52137 27.34869 34.46077 32.10526 30.06495 33.53245 31.95952 32.25491
#> [2225] 30.19406 28.99267 31.04701 34.86105 35.15252 29.40493 31.07870 31.72205
#> [2233] 34.65826 28.64792 29.13723 35.17083 35.33939 33.94909 26.90600 32.43748
#> [2241] 36.42972 33.52442 31.82696 29.72724 21.93759 33.22105 33.60301 32.05989
#> [2249] 31.62671 30.01537 30.17409 32.66473 33.88590 33.35547 31.64108 32.53439
#> [2257] 32.49363 30.77872 29.12251 33.54517 30.40360 22.47296 32.44514 34.92591
#> [2265] 31.18280 34.23283 32.52124 32.01273 32.09897 32.40320 29.85263 29.29998
#> [2273] 31.89082 29.90213 32.59771 25.84374 32.10957 28.93898 29.36009 34.19938
#> [2281] 34.79246 23.54357 28.48740 32.10072 32.16399 26.11152 33.81278 27.24773
#> [2289] 28.87001 35.46481 33.95131 34.39654 33.42674 35.54917 27.21369 32.24185
#> [2297] 34.68948 33.58336 26.97526 30.10134 29.20477 29.75305 23.57250 28.88962
#> [2305] 31.84888 30.63259 24.78556 30.44270 29.38224 28.07361 25.58092 33.05690
#> [2313] 24.16335 28.48802 34.08931 32.05744 30.51171 30.19070 30.47759 32.43933
#> [2321] 28.63302 26.66056 29.66031 31.21715 28.07231 25.60228 32.55769 31.98817
#> [2329] 25.29938 29.31709 31.09785 30.48952 29.06496 32.40383 32.99965 26.28913
#> [2337] 29.10705 27.71397 33.24198 26.23017 27.79155 26.71431 25.19416 26.23726
#> [2345] 27.78541 28.25321 27.75116 32.72528 22.37814 31.37764 26.97054 30.05344
#> [2353] 27.51123 29.58550 24.28045 28.50226 23.73669 25.57334 27.81738 24.42844
#> [2361] 24.13241 23.35212 29.94806 24.88575 23.83527 24.95397 26.96026 27.65092
#> [2369] 26.72207 25.92629 32.88987 24.32862 24.17117 26.75829 24.66504 28.77891
#> [2377] 28.88691 30.83586 30.72718 27.25970 25.49690 27.08191 25.32555 30.39726
#> [2385] 30.71697 22.03193 32.82246 25.35901 34.15933 30.18242 31.91005 31.26316
#> [2393] 29.44548 29.57430 29.81689 30.59594 29.38545 35.40488 29.22924 26.66040
#> [2401] 27.69960 29.15548 30.58644 27.50165 28.89182 27.06733 28.27239 28.41710
#> [2409] 29.63893 33.71777 28.05403 20.72121 27.37820 23.27249 28.24099 25.50053
#> [2417] 27.59200 30.97598 31.44458 28.11869 28.71480 25.43543 20.37314 23.34476
#> [2425] 30.11653 22.19957 32.37248 24.93172 28.60400 29.03225 25.68851 24.57525
#> [2433] 27.10077 30.26599 31.60094 25.94802 21.40846 24.64962 27.60463 24.65699
#> [2441] 26.60271 26.70539 27.87554 26.33347 29.32239 28.77794 27.53977 31.52412
#> [2449] 29.63410 28.34369 27.42890 26.80886 28.33921 25.96064 26.52904 26.64520
#> [2457] 28.26310 23.37463 32.48039 30.88009 27.11269 27.00649 26.72819 29.99780
#> [2465] 24.90757 28.41983 20.68215 23.42963 29.50833 26.15110 25.25952 31.67656
#> [2473] 27.98509 25.28030 26.17278 26.99829 24.60872 25.44754 29.45328 26.11322
#> [2481] 23.56300 33.77039 23.86930 23.73376 30.40030 25.19185 28.54664 27.25823
#> [2489] 31.19871 22.24244 26.08822 26.04759 30.14671 24.08468 26.25589 30.32912
#> [2497] 25.82817 28.73244 24.41375 28.71217 23.58985 30.35982 25.92150 27.08839
#> [2505] 21.73662 19.46649 25.16708 29.98558 33.14664 25.27980 25.95875 31.06528
#> [2513] 24.75636 27.23088 27.53650 24.07616 26.76929 27.42564 23.27010 30.38292
#> [2521] 24.03951 19.61840 24.41128 19.87958 18.64480 26.18361 27.00485 25.01198
#> [2529] 28.30827 30.50828 19.17117 32.84297 28.84175 28.46918 29.59287 27.49412
#> [2537] 29.85356 28.85293 31.95159 28.01525 25.35666 32.08801 26.48227 32.74715
#> [2545] 29.50077 20.52689 29.64219 29.85995 28.58707 26.25370 28.23152 30.58889
#> [2553] 25.05655 27.14433 28.46775 29.71023 25.84133 30.06209 29.47165 23.48123
#> [2561] 28.75868 27.79769 23.76052 24.93076 26.50228 28.06057 26.41341 24.04330
#> [2569] 27.30443 26.89693 31.37284 27.67890 30.01317 24.68498 29.98221 31.01952
#> [2577] 28.46937 25.33432 28.34946 21.66875 27.01450 29.93160 27.72241 28.22009
#> [2585] 25.00259 23.47730 26.46895 27.74207 25.96633 30.78102 30.50804 32.92266
#> [2593] 26.15728 28.29991 26.92316 31.47986 27.32885 26.20681 28.54223 27.78519
#> [2601] 22.61277 28.74696 26.78410 26.21187 29.12832 27.65158 29.24882 30.07646
#> [2609] 29.53248 29.51683 28.37633 26.32679 31.62319 28.48071 31.94469 29.56091
#> [2617] 27.40154 25.03678 26.32133 29.63015 28.44247 24.91206 35.80138 32.98811
#> [2625] 23.42943 24.36938 26.78973 26.25172 28.23724 30.09385 29.94628 31.78340
#> [2633] 28.77054 27.28835 32.24916 29.48112 25.79738 29.41865 32.44103 28.81777
#> [2641] 32.23389 28.36321 27.90136 29.36972 31.97812 27.45942 27.73601 34.34685
#> [2649] 34.34384 31.29609 25.33786 28.77034 31.90348 28.09016 32.65548 33.25802
#> [2657] 32.66663 26.05845 31.00121 27.14788 28.27050 32.26327 23.33746 28.76414
#> [2665] 28.09793 24.88292 27.71885 28.92892 21.02453 21.36955 23.09606 25.99923
#> [2673] 28.99030 24.45597 13.21266 22.74624 24.05351 28.19997 24.11403 22.01786
#> [2681] 24.64577 20.34329 27.65135 29.02386 27.21110 31.53782 24.64870 26.55449
#> [2689] 21.59712 25.91713 24.34107 24.27461 29.40141 26.38931 27.23113 27.76932
#> [2697] 26.71785 26.31210 23.41939 30.94528 29.62611 26.53583 26.23408 30.39078
#> [2705] 31.38246 35.69611 31.01466 35.62365 27.42046 29.42245 30.62492 27.23188
#> [2713] 36.99661 33.72039 22.99340 26.81183 31.41989 33.87580 29.44431 28.30982
#> [2721] 29.74636 33.06683 29.45473 26.83524 33.53925 28.84378 26.99711 23.13976
#> [2729] 25.21560 26.15585 29.36244 30.66021 26.12581 26.76482 28.00993 25.19185
#> [2737] 24.84225 25.07527 28.22207 28.03341 23.05552 32.90440 29.98094 24.38508
#> [2745] 31.00994 25.58628 29.25601 25.28736 29.24698 25.84544 27.29718 31.65546
#> [2753] 34.20184 32.22251 31.08959 29.59454 33.23367 29.50439 31.65828 30.33723
#> [2761] 27.68964 25.89701 35.20987 31.30656 31.80432 24.93995 32.41997 32.08062
#> [2769] 24.83146 31.63866 34.41868 19.64117 35.11050 33.75511 28.78751 32.86032
#> [2777] 29.64700 28.05796 28.57125 30.32345 28.14315 28.14264 22.83808 26.22667
#> [2785] 26.89741 26.09364 26.39829 28.30444 28.88171 26.50474 24.10664 28.20806
#> [2793] 26.05444 26.31944 23.83822 27.37062 31.04328 25.34518 23.30295 27.07480
#> [2801] 31.76052 25.71475 28.56624 29.72711 27.49647 28.90377 24.81264 25.40431
#> [2809] 27.84753 27.81672 23.41189 25.21628 30.79000 29.47257 21.44209 31.20066
#> [2817] 27.19878 30.61345 26.28332 24.55540 31.80614 25.30622 36.97316 32.33546
#> [2825] 24.53350 31.63158 28.73660 31.48480 26.05848 32.26723 28.27527 25.21665
#> [2833] 27.79789 27.31637 31.73818 21.93887 28.61785 27.38073 24.49218 21.91423
#> [2841] 25.57689 28.79325 23.64613 21.55929 25.75687 23.62818 23.67837 28.17500
#> [2849] 24.38035 28.32165 27.25867 28.25892 26.34752 21.97145 23.15251 25.05593
#> [2857] 27.58810 29.91989 26.37059 32.30850 27.52760 28.34386 29.22558 17.66219
#> [2865] 26.07341 28.17582 24.36987 23.36740 23.60439 30.91441 26.50341 32.34342
#> [2873] 21.73110 26.21194 30.99582 20.99836 26.42138 24.53622 26.88477 25.68570
#> [2881] 21.43280 25.06312 18.24747 26.25057 28.80943 32.48641 26.59329 32.29935
#> [2889] 31.15710 29.39771 32.02244 20.93543 27.37178 26.54328 33.24758 30.12110
#> [2897] 29.00168 32.19724 24.16712 26.20050 29.88930 30.22445 30.99659 32.21389
#> [2905] 33.84773 30.13910 27.56439 32.11109 30.43160 31.10440 27.95656 32.03860
#> [2913] 24.16485 28.33422 30.58927 30.36454 31.86974 29.88920 25.31799 34.57660
#> [2921] 31.62665 22.80686 34.34915 28.35731 25.86908 32.57844 32.28330 25.93055
#> [2929] 26.20287 24.71356 31.01247 24.20271 25.45677 27.97908 28.99792 24.25260
#> [2937] 21.04251 31.14816 31.03407 21.47653 32.08197 32.95943 31.19358 29.63081
#> [2945] 33.88481 21.17855 23.68079 34.21596 26.34741 27.63775 25.24188 27.64668
#> [2953] 24.72436 22.04981 22.44937 29.03974 22.77264 31.88907 27.59201 26.93972
#> [2961] 34.07685 26.12381 27.17484 24.91007 28.05208 24.99006 22.08789 30.13440
#> [2969] 23.60120 25.39244 24.40875 28.07872 26.31935 25.53006 23.60450 27.99267
#> [2977] 26.78723 22.85847 31.01748 30.21300 23.30108 27.89469 27.80850 25.83220
#> [2985] 27.11655 27.52683 23.97473 28.50146 24.53595 27.90363 28.52175 28.03481
#> [2993] 29.59482 22.27668 27.51571 29.26808 23.16430 33.70967 28.05712 25.03687
#> [3001] 27.21475 31.11418 30.48275 25.84613 30.98513 27.20501 28.70281 30.17593
#> [3009] 30.07184 29.56434 23.37826 27.71533 26.35148 31.45115 30.07173 27.36858
#> [3017] 27.51532 27.73742 29.28495 25.31669 22.91757 33.88187 29.31275 26.90255
#> [3025] 29.66063 30.04976 25.53501 35.69114 25.33683 20.96962 32.02690 29.98882
#> [3033] 28.02978 30.15729 27.43023 25.01178 30.71883 25.98159 27.30870 28.44865
#> [3041] 23.90882 29.04403 31.36541 26.20788 27.07060 33.59923 26.45626 23.19385
#> [3049] 28.95324 27.34040 27.36881 22.41250 30.13849 25.66741 26.87106 27.98529
#> [3057] 26.33208 26.85653 28.28551 23.84855 28.09052 24.15574 29.17051 24.86656
#> [3065] 20.78904 27.72508 27.45796 29.50065 26.27605 29.16263 30.20675 20.83210
#> [3073] 32.63449 28.19857 28.27421 21.59269 28.91043 24.75733 26.55688 33.66034
#> [3081] 28.01820 30.16694 28.37932 33.70979 33.35138 31.91221 27.46571 26.73972
#> [3089] 27.60107 31.47419 29.15388 30.57798 31.01421 29.84514 28.14328 33.66793
#> [3097] 27.42566 27.78774 25.51670 29.00369 26.39291 27.93629 23.17189 26.45978
#> [3105] 26.35731 34.04452 31.32409 26.65451 29.21894 30.71066 23.74956 25.80044
#> [3113] 27.18732 23.63033 27.76332 27.14575 32.89265 32.54904 23.16807 30.14139
#> [3121] 31.44363 25.38339 29.22835 33.25300 27.96575 28.97116 26.62301 27.77176
#> [3129] 27.13834 27.18630 28.70987 24.70833 31.83584 23.69744 27.24146 26.43238
#> [3137] 30.80695 27.12562 26.39001 24.14516 30.76877 27.23903 30.62145 29.80491
#> [3145] 37.44835 30.35344 28.08976 27.02186 27.60946 34.10496 21.89942 28.72648
#> [3153] 26.57600 20.96252 27.94837 23.01798 26.65418 23.44085 27.36328 24.47141
#> 
#> $mod2_unimodal$ref_properties$cell$knn_hubness
#>  [1]  7  4  2  2  2  2  4  2  3  2  4  3  4  3 10  3  0  3  6  1  3  6  4  0  5
#> [26]  5  2  5  6  1  1  2  2  4  8  0  4  3  4  4  5  6  4  8  3  5  2  3 12  0
#> [51]  3  3  6  6  6  8  5  4  4  8  6  2  4  5  6  2  8  3  7  6  1  3  3  2  4
#> [76]  6  7  3  4  1
#> 
#> $mod2_unimodal$ref_properties$cell$prop_outliers_cell
#> [1] 0
#> 
#> 
#> $mod2_unimodal$ref_properties$feature
#> $mod2_unimodal$ref_properties$feature$mean_expression
#>  [1] 2.515016 3.079690 3.009241 3.012685 3.825412 3.106977 3.255403 3.422643
#>  [9] 3.429865 4.009745 3.690031 3.964169 3.293251 3.086738 3.615171 2.260877
#> [17] 2.708356 1.585968 2.714355 2.727232 2.038118 2.828071 2.344227 2.907054
#> [25] 1.481160 2.540602 2.487980 2.674601 3.315180 2.040646 2.798478 2.629667
#> [33] 2.202299 2.924119 3.113251 2.979615 2.019863 2.041088 2.401802 2.711177
#> [41] 1.651634 3.529243 2.755073 2.636220 3.334218 2.258513 2.251378 3.250704
#> [49] 1.407267 2.339233 1.478937 1.976841 1.897893 2.733963 3.236528 2.898465
#> [57] 3.524952 3.090174 2.907645 2.187672
#> 
#> $mod2_unimodal$ref_properties$feature$variance
#>  [1] 10.705826 10.657424 10.686240 10.727286  9.456150 10.852786 10.800694
#>  [8] 10.252386 10.307432  8.824670  9.746971  9.117682 10.491397 10.707373
#> [15]  9.842093 10.772938 11.169928  9.458253 11.227135 11.326249 10.434543
#> [22] 10.997915 10.970942 11.057358  8.914374 10.936613 11.066427 10.895355
#> [29] 10.637033 10.469323 11.320174 11.094742 10.825762 11.164130 10.879057
#> [36] 11.018792 10.256142 10.473951 10.867837 11.194229  9.533124 10.399472
#> [43] 10.989574 11.152432 10.774184 10.754517 10.697305 10.739998  8.695127
#> [50] 10.910872  8.879274 10.457783 10.261985 10.827565 10.675327 10.971812
#> [57] 10.365650 10.731177 11.051602 10.704844
#> 
#> $mod2_unimodal$ref_properties$feature$scaled_variance
#>  [1]  0.25811397  0.18440969  0.22828983  0.29079300 -1.64482949  0.48189724
#>  [7]  0.40257480 -0.43236110 -0.34853958 -2.60641425 -1.20198204 -2.16023161
#> [13] -0.06840779  0.26046954 -1.05713391  0.36030898  0.96482494 -1.64162668
#> [19]  1.05193711  1.20286238 -0.15498122  0.70289222  0.66181950  0.79340921
#> [25] -2.46981846  0.60954492  0.80721957  0.54671878  0.15336024 -0.10202047
#> [31]  1.19361178  0.85033638  0.44074710  0.95599644  0.52190197  0.73468220
#> [37] -0.42664200 -0.09497365  0.50481678  1.00183014 -1.52761654 -0.20838562
#> [43]  0.69019176  0.93818247  0.36220680  0.33225860  0.24513948  0.31014897
#> [49] -2.80367552  0.57034799 -2.52326672 -0.11959291 -0.41774408  0.44349205
#> [55]  0.21167233  0.66314335 -0.25988938  0.29671788  0.78464373  0.25661885
#> 
#> $mod2_unimodal$ref_properties$feature$sd
#>  [1] 3.271976 3.264571 3.268982 3.275254 3.075085 3.294357 3.286441 3.201935
#>  [9] 3.210519 2.970635 3.122014 3.019550 3.239043 3.272212 3.137211 3.282215
#> [17] 3.342144 3.075427 3.350692 3.365449 3.230254 3.316310 3.312241 3.325261
#> [25] 2.985695 3.307055 3.326624 3.300811 3.261447 3.235633 3.364547 3.330877
#> [33] 3.290253 3.341277 3.298342 3.319457 3.202521 3.236348 3.296640 3.345778
#> [41] 3.087576 3.224821 3.315053 3.339526 3.282405 3.279408 3.270674 3.277194
#> [49] 2.948750 3.303161 2.979811 3.233850 3.203433 3.290527 3.267312 3.312372
#> [57] 3.219573 3.275848 3.324395 3.271826
#> 
#> $mod2_unimodal$ref_properties$feature$cv
#>  [1] 142.01764 117.40274 118.20020 120.18708  92.04002 117.69947 111.71406
#>  [8] 104.14235 104.82284  86.01514  96.17832  87.90590 106.36466 116.66442
#> [15]  98.26143 151.45893 128.46741 193.72125 128.46741 128.46741 158.41816
#> [22] 127.28917 150.83290 124.16029 215.82928 140.98974 144.76127 130.33748
#> [29] 109.12773 162.60883 122.89131 131.90292 153.06520 118.89404 111.91153
#> [36] 117.65058 162.60883 165.72052 140.48273 128.46741 196.16578 102.13807
#> [43] 131.20667 133.81136 107.28943 151.45893 157.10998 106.50084 209.47997
#> [50] 150.83290 211.85375 167.80695 168.67387 131.69474 112.05785 120.65950
#> [57] 102.66407 115.35524 120.65950 158.13430
#> 
#> $mod2_unimodal$ref_properties$feature$zero_fraction_feature
#>  [1] 0.6250 0.5250 0.5375 0.5375 0.3875 0.5250 0.5000 0.4625 0.4625 0.3500
#> [11] 0.4125 0.3625 0.4875 0.5250 0.4250 0.6750 0.6000 0.7875 0.6000 0.6000
#> [21] 0.7125 0.5750 0.6625 0.5625 0.8000 0.6250 0.6375 0.6000 0.4875 0.7125
#> [31] 0.5875 0.6125 0.6875 0.5625 0.5250 0.5500 0.7125 0.7125 0.6500 0.6000
#> [41] 0.7750 0.4500 0.5875 0.6125 0.4875 0.6750 0.6750 0.5000 0.8125 0.6625
#> [51] 0.8000 0.7250 0.7375 0.5875 0.5000 0.5625 0.4500 0.5250 0.5625 0.6875
#> 
#> $mod2_unimodal$ref_properties$feature$detection_freq_feature
#>  [1] 0.3750 0.4750 0.4625 0.4625 0.6125 0.4750 0.5000 0.5375 0.5375 0.6500
#> [11] 0.5875 0.6375 0.5125 0.4750 0.5750 0.3250 0.4000 0.2125 0.4000 0.4000
#> [21] 0.2875 0.4250 0.3375 0.4375 0.2000 0.3750 0.3625 0.4000 0.5125 0.2875
#> [31] 0.4125 0.3875 0.3125 0.4375 0.4750 0.4500 0.2875 0.2875 0.3500 0.4000
#> [41] 0.2250 0.5500 0.4125 0.3875 0.5125 0.3250 0.3250 0.5000 0.1875 0.3375
#> [51] 0.2000 0.2750 0.2625 0.4125 0.5000 0.4375 0.5500 0.4750 0.4375 0.3125
#> 
#> $mod2_unimodal$ref_properties$feature$dispersion
#>  [1] 0.9580280 0.8270042 0.7509567 0.8486399 0.6988876 0.8658228 0.7644020
#>  [8] 0.7049659 0.7416784 0.6473779 0.7053331 0.6375144 0.6505228 0.7826087
#> [15] 0.7120790 0.8028933 0.7014147 0.7974684 0.7014147 0.7014147 0.7215190
#> [22] 0.8101266 0.8815843 0.7900587 1.1645570 0.8945148 0.9168174 0.7432188
#> [29] 0.7443038 0.7932489 0.6418466 0.7176832 0.7614411 0.6537804 0.6418648
#> [36] 0.6747809 0.7932489 0.8582278 0.7154081 0.7014147 0.9620253 0.7041725
#> [43] 0.8607595 0.7609829 0.6762726 0.8028933 0.9873418 0.6096556 0.8227848
#> [50] 0.8815843 1.0098453 0.8095762 0.7468354 0.8888546 0.7848101 0.6915390
#> [57] 0.7377939 0.7318757 0.6915390 0.8752260
#> 
#> $mod2_unimodal$ref_properties$feature$bcv
#>  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>  [8] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [15] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [22] 0.0000000 0.0000000 0.0000000 0.8113124 0.0000000 0.0000000 0.0000000
#> [29] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [36] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [43] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [50] 0.0000000 0.2091813 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [57] 0.0000000 0.0000000 0.0000000 0.0000000
#> 
#> $mod2_unimodal$ref_properties$feature$gene_cor
#>    [1]  0.1448735818 -0.0230042315 -0.0262185314  0.0339289720 -0.0601816539
#>    [6] -0.0453220091 -0.2100195129  0.0228149938 -0.0355853263  0.1574744873
#>   [11]  0.1796128492 -0.1369680287  0.0806149097 -0.1457595400 -0.2726254517
#>   [16]  0.2364056413  0.3423199873 -0.2282143838  0.1007130220  0.0033415375
#>   [21] -0.1045733375  0.0580089627 -0.0452551715  0.0236781130 -0.1577830624
#>   [26]  0.1124456805  0.1969515631 -0.0662069701 -0.0420504511  0.0541217251
#>   [31] -0.0420645737 -0.0592861416 -0.0063750038  0.0324153203 -0.0221132530
#>   [36] -0.1649298382 -0.0293894885 -0.0135275419 -0.0109011663 -0.0126688913
#>   [41]  0.0208967881  0.0395130926 -0.0267496445 -0.0843321873 -0.0063765059
#>   [46]  0.1808425387  0.0994331355 -0.0083825013  0.0123700405  0.0078567107
#>   [51] -0.0535635332 -0.0661291884 -0.0217857035 -0.0223074592  0.0857770425
#>   [56]  0.3075960602  0.0879584776 -0.2704554335  0.1048023105 -0.0578755799
#>   [61] -0.0525307285  0.1285909017 -0.1618890848  0.1843506002  0.1884703085
#>   [66]  0.2806708831 -0.0203243878 -0.0749972906  0.1692692097  0.1810197226
#>   [71] -0.1381021056  0.0230389285 -0.1847181876  0.0180870465  0.0792951825
#>   [76]  0.0246556690 -0.0751268539 -0.0316828329 -0.0834961876 -0.0088914057
#>   [81] -0.0606107840 -0.1529301399  0.2192026673 -0.0703989088 -0.0230208510
#>   [86]  0.1211446639 -0.1164400667  0.0069813307 -0.0577360477 -0.1382111531
#>   [91] -0.1004847533 -0.0249427312  0.0766174480 -0.0138290750  0.0432859872
#>   [96]  0.1618880453 -0.1911831822  0.0644245776  0.0161783422 -0.2437417501
#>  [101] -0.1515438922 -0.1151015519  0.0022023847 -0.1513464388  0.0955732589
#>  [106]  0.2680409362  0.1017824170  0.1358321672  0.1400771945  0.0590061024
#>  [111]  0.0721070616  0.2066189304 -0.0384014011  0.1421232529 -0.0497609538
#>  [116]  0.0944439136 -0.1765446015  0.0808718211  0.0145710579 -0.0247269411
#>  [121]  0.1460690700  0.1407680763 -0.0298844847 -0.0495388757 -0.0492309895
#>  [126]  0.0827044195  0.0732342627  0.0844291330 -0.1112053777 -0.0740052258
#>  [131]  0.0470063353  0.1698970101  0.1160266641  0.0638887654 -0.0987603903
#>  [136] -0.0498724463 -0.1706430583 -0.0332868385 -0.1181889882  0.0197373346
#>  [141] -0.1124189972  0.1150917177  0.0124507088 -0.1127330365  0.0351614796
#>  [146]  0.0034620893 -0.2190037914  0.0276339813 -0.0947963873 -0.0852398951
#>  [151]  0.2306319538 -0.1809700106  0.0590920043 -0.1190385406 -0.0235571670
#>  [156]  0.0215312615  0.1119318735 -0.0783044259  0.0092364956  0.0916574177
#>  [161] -0.2283831240 -0.0236705049 -0.1869449564 -0.1926412092 -0.2039990302
#>  [166] -0.0279375490 -0.0003781771 -0.0820270922  0.0394512434  0.0161722469
#>  [171] -0.1005065165  0.0482313186  0.1222140193  0.2420942234  0.0437946026
#>  [176] -0.0478358346 -0.0760595502 -0.0094104018 -0.1377503390  0.0128702946
#>  [181]  0.1611744430  0.1849862609  0.2865721237  0.0684483650 -0.0373864319
#>  [186]  0.0500406273 -0.0726840061  0.0585309353  0.1366174080  0.0540146645
#>  [191]  0.1401767918  0.2397453506 -0.1958763439 -0.1733549480 -0.1382518196
#>  [196] -0.1590573115  0.0641316920 -0.0439552875  0.0041972825 -0.1759189542
#>  [201] -0.0506486078  0.0594611791  0.0127555557  0.0260099838  0.2965078866
#>  [206] -0.0153219756 -0.0636575493  0.1859398133 -0.0520441797  0.0831000707
#>  [211]  0.0218440229  0.0401123275 -0.0008376513 -0.1034239096 -0.1784398809
#>  [216]  0.0854593183  0.0309099516 -0.1770771078 -0.3802625265  0.0779787195
#>  [221] -0.0347744589  0.0846032315 -0.1843335193  0.0783677075  0.0413108209
#>  [226] -0.2517909198  0.0376240892  0.0084453734  0.2139638626  0.0889822686
#>  [231] -0.1230577103  0.0636165856 -0.1908326255 -0.2419374416 -0.1329679042
#>  [236] -0.1069824764  0.0122072195 -0.1864300084  0.0830815833 -0.0776756789
#>  [241]  0.0613627880 -0.0061908470  0.1530858886  0.0087238261 -0.0467065704
#>  [246] -0.1437663149 -0.2263900633  0.1040341760 -0.0638963352 -0.0587157197
#>  [251] -0.0993865088  0.0977895331 -0.1393211929 -0.0463736540  0.1268629623
#>  [256]  0.1188792870  0.1841949593  0.1943495852 -0.2428058295 -0.0375460382
#>  [261] -0.0985616166  0.1217531999 -0.0462265334  0.0825031390 -0.1624782443
#>  [266]  0.0813837059  0.0632225263  0.1384199446  0.2518758949 -0.2157411892
#>  [271] -0.1893903789 -0.0039960320  0.1181410423 -0.0085096361 -0.0882619975
#>  [276] -0.1909998441  0.0474466974 -0.1212918420  0.0734396861 -0.0731185099
#>  [281] -0.1697797916 -0.0034697839 -0.1159569775  0.0124250106  0.0037855088
#>  [286] -0.1301049281 -0.1936346652 -0.1056510171 -0.0470788977 -0.0357644216
#>  [291] -0.0326385523  0.0415127849  0.0574311067  0.1498758697 -0.0753734411
#>  [296] -0.0364619492  0.0700741950 -0.0592217474 -0.1034658621  0.0160234248
#>  [301] -0.1322254641  0.1354206106  0.1282775425  0.1796725507  0.1042608182
#>  [306] -0.1426365883 -0.0840281996  0.0451524204  0.0699309655  0.0582467983
#>  [311]  0.0842234666 -0.0397751405  0.0403296672 -0.0792706918  0.0437951541
#>  [316] -0.1027042164  0.1301275719 -0.0166082220 -0.0091135864 -0.0766629383
#>  [321] -0.0628513693 -0.1716969988  0.0694548067 -0.0875984996  0.0360956859
#>  [326]  0.1973483059 -0.0569813207 -0.0051549110  0.0387993902 -0.0268403392
#>  [331]  0.1275420014  0.2161128001 -0.0104186005 -0.1022046196 -0.0325230527
#>  [336]  0.0801534084 -0.0214491752  0.0749056300 -0.0646926702  0.1044063023
#>  [341]  0.0560804467 -0.0703484988 -0.0923846411 -0.1489070816  0.0436357943
#>  [346] -0.0574461176 -0.0857664818 -0.0713442444 -0.0604373737 -0.2330916058
#>  [351] -0.2647019934 -0.1712684839  0.0484572882  0.0066421230 -0.1649731557
#>  [356]  0.0866333179 -0.0706737843 -0.1079336208  0.1284188997  0.2283307389
#>  [361] -0.1637709817 -0.0676166908  0.0276508235 -0.1264397863 -0.0795217453
#>  [366]  0.1948859144 -0.1184435838 -0.0790371174 -0.0588680933  0.1071036520
#>  [371] -0.0246536389  0.1573348162 -0.0066427377  0.0137805186 -0.2143385187
#>  [376] -0.0599863568  0.1172245391 -0.1692130593  0.0521399521 -0.3085485754
#>  [381] -0.0044876108 -0.0403272049  0.0480821960  0.2330989560 -0.0086486571
#>  [386]  0.0779182523  0.0393098307  0.1069661493 -0.1042883349 -0.1339651577
#>  [391]  0.1156016961  0.0889860474 -0.2082082649  0.1046426304  0.0970957328
#>  [396] -0.0295064993 -0.0603929875 -0.0717355915 -0.2012846364  0.1290430498
#>  [401] -0.0833160872  0.1643397884 -0.1037731436 -0.2678520651  0.2342196578
#>  [406] -0.1415010516  0.0725694849  0.0671232487  0.1759628515  0.0879652088
#>  [411] -0.1142411193  0.0566298702  0.0195109154  0.1314000941 -0.1275069283
#>  [416] -0.1539915321  0.1455790259  0.0720066888  0.0431667709 -0.0535490060
#>  [421]  0.0212006686 -0.1537532456  0.0722654634 -0.0025591865  0.0789051217
#>  [426]  0.1111431023 -0.2234494764 -0.0044565419  0.0389193611 -0.1495583121
#>  [431] -0.1235393237 -0.0039927692  0.1741676908 -0.0307453045 -0.0936313193
#>  [436] -0.1638785619 -0.3567024527  0.0410836570 -0.0164013118 -0.0123451721
#>  [441]  0.0988549316 -0.3554396318 -0.1044468417 -0.0330435239 -0.0923282916
#>  [446] -0.2006968792 -0.2214113942 -0.1299706419  0.1511442695  0.1279511844
#>  [451] -0.1145915051  0.0589185625  0.0398842532  0.1617669262 -0.1415161462
#>  [456] -0.0233777232  0.0730665720 -0.0441112216 -0.1012210345  0.0653405543
#>  [461] -0.0673303112 -0.0708481373  0.1225590959  0.1106521175 -0.0311218836
#>  [466]  0.0369168831 -0.0064252842 -0.1492737526 -0.0763200817  0.0478104800
#>  [471]  0.0091795068 -0.0318588852  0.1177238328 -0.0294534938  0.0595872449
#>  [476] -0.2238025725  0.0287154132  0.0237356518  0.1004305593  0.0515052485
#>  [481] -0.1247784976 -0.0095213153 -0.0407833247 -0.1283240447 -0.0863787150
#>  [486]  0.2278617675  0.0414949995 -0.0106548618 -0.0526564395 -0.1373058266
#>  [491] -0.0131028497 -0.0360139948  0.0293396126  0.2248795528 -0.1456599545
#>  [496]  0.0536903031 -0.0716988997 -0.1358755624 -0.2553245125  0.0463818032
#>  [501] -0.0502086738 -0.0297314137  0.1173756930 -0.0037492241 -0.0164715415
#>  [506] -0.0782784702 -0.1087214442  0.0502256379 -0.0261166093 -0.0309343521
#>  [511] -0.0949663751  0.0296291131  0.1522969790  0.1760381727  0.2324029055
#>  [516]  0.0201579371 -0.0134336687 -0.0051043598 -0.0532095451 -0.2342332877
#>  [521]  0.0789438382 -0.1700615078 -0.0376192956 -0.0356077984 -0.0139715795
#>  [526]  0.0534058208  0.0769099267  0.1495726666 -0.0866572886 -0.0978428205
#>  [531]  0.1095912836  0.1078736647  0.0700015838 -0.1158922230 -0.1879350376
#>  [536]  0.0284997961 -0.0997271119 -0.0860961905 -0.0794743861 -0.1678008362
#>  [541] -0.0033451078  0.2051873262  0.1731365908  0.0799528253 -0.1938520418
#>  [546]  0.0044044657  0.0414430792 -0.0062866144  0.0682535695 -0.0849143338
#>  [551] -0.0833334729  0.0520649749  0.1889575531  0.0639272715 -0.1080595333
#>  [556]  0.0005914848 -0.2757225298 -0.1618113858  0.1070355175 -0.1442394747
#>  [561] -0.0606061123 -0.0527417218 -0.1065021177  0.0092195818  0.0490361532
#>  [566]  0.0608137845 -0.1158830397 -0.0793463592  0.2116138504 -0.1304137940
#>  [571] -0.0646160490 -0.0358406383 -0.1265524932  0.1178402754  0.0311733043
#>  [576]  0.0570722010 -0.0004934536 -0.0801402624 -0.0422010836 -0.3347204882
#>  [581]  0.0324074789 -0.0158288534 -0.1680290660 -0.0279953568  0.1443834900
#>  [586]  0.0326594569 -0.1806396056  0.0500342817 -0.0795710392  0.0916368128
#>  [591]  0.0829741143 -0.0947488589 -0.0740205561 -0.1667477829 -0.0592952129
#>  [596]  0.0302354815 -0.1729413632 -0.0129957535 -0.0334999270 -0.0602814038
#>  [601]  0.0084990735 -0.0756575983 -0.0190801931 -0.2829490715  0.0335324187
#>  [606]  0.0401688565 -0.0294986944 -0.1503923455  0.1270841074  0.1332572410
#>  [611] -0.2849454859  0.2238516187  0.0349606341 -0.0133778578  0.0264168848
#>  [616] -0.0330730108  0.1503391898  0.3209333917 -0.1784038220 -0.0902774031
#>  [621]  0.1314468504  0.0520568141 -0.0322470978  0.0939676295  0.1677644709
#>  [626]  0.2220664512  0.0382070898 -0.0581987217 -0.0412379500 -0.0993644413
#>  [631] -0.0223572662  0.1836000842 -0.0636657020  0.0968985396  0.0782156601
#>  [636] -0.0358173301 -0.0156093512  0.0487000028 -0.1603953003  0.1093059315
#>  [641] -0.0397731910 -0.1008798464  0.1862532050  0.0177522308 -0.0486671976
#>  [646] -0.0731361229 -0.0210968528  0.0409973941 -0.0607583715 -0.0476150159
#>  [651]  0.0083214388 -0.1361250662 -0.0456094508 -0.0238328734  0.0893657504
#>  [656]  0.0550965552  0.0163340755 -0.1443925039 -0.2230702253  0.1748366570
#>  [661] -0.2350553249  0.0477809149 -0.0381958177 -0.0367143835  0.2470414953
#>  [666] -0.1636002713 -0.0348408901 -0.2836283760 -0.0829370209 -0.1108112320
#>  [671] -0.2663184756  0.1174389610  0.0354196266  0.0035969365  0.0021648318
#>  [676] -0.0193520467 -0.0754118884  0.0258827516 -0.2126066429 -0.0825028702
#>  [681]  0.0073758704  0.0443750541 -0.0516819330 -0.1051606193  0.0750644836
#>  [686] -0.0231385021 -0.1742117746  0.1783400648  0.0248944705 -0.1739606457
#>  [691]  0.0736472077 -0.2184122992 -0.0268200668  0.1346027088  0.0287233136
#>  [696] -0.0220738448  0.1790159534 -0.0907987991  0.0649035280 -0.0440633158
#>  [701] -0.0479099916  0.0823532126 -0.2079699727  0.1764392588  0.0442046952
#>  [706] -0.1584089608  0.1074454134 -0.0754228866 -0.0690894842  0.0651311950
#>  [711]  0.0661912304 -0.1153102978 -0.0116354210 -0.0741686671  0.0531772213
#>  [716]  0.1262374547 -0.1411400178 -0.0097374445 -0.0530862525 -0.0533372808
#>  [721] -0.0625149740 -0.0305446477  0.0286646725  0.0307077173  0.0014564640
#>  [726]  0.2632822321 -0.0588268935 -0.2170727634 -0.0489121799  0.1103389994
#>  [731]  0.1159738243 -0.0087597130 -0.0104244350 -0.1981276619 -0.0240281973
#>  [736]  0.0223316210 -0.0184532758 -0.0049372121  0.0988198084 -0.1166107579
#>  [741] -0.0058722109 -0.0074424370 -0.0992944799 -0.1275949008  0.0251867487
#>  [746] -0.0941992835  0.0700645716  0.0103932202 -0.0050866135 -0.1283952197
#>  [751] -0.1284203097 -0.0474229906 -0.0729934519 -0.0829789349 -0.1147254867
#>  [756] -0.0933772094 -0.0796166905  0.0371589528  0.0115182196  0.0692908940
#>  [761] -0.0921403925 -0.0601185851  0.1894646707  0.0055227565 -0.0568226288
#>  [766] -0.0813298482 -0.1289308432 -0.0460766855  0.0007941202  0.2703457642
#>  [771] -0.1273427965  0.0483862258  0.0076114745  0.1024408669 -0.1411966191
#>  [776]  0.0619723921  0.0594738900 -0.1619320951  0.2684578510  0.1101053707
#>  [781] -0.0543174926  0.0643067385  0.1309761381  0.0264387283  0.0274974099
#>  [786] -0.0334096849  0.0206389444 -0.0719829924  0.0132434372  0.0578548448
#>  [791]  0.2354590643  0.0398395053  0.0229060633 -0.0686155074  0.1090745886
#>  [796]  0.1567669148 -0.1505885217 -0.0322631917  0.0149048869 -0.0573989879
#>  [801] -0.0790352665 -0.0640917032 -0.1159305227  0.3178575094 -0.1220598544
#>  [806] -0.1937169450  0.1052706736 -0.1220401020  0.1046724754 -0.0445100260
#>  [811] -0.1696567809 -0.0037840235 -0.0678437176 -0.0323020744  0.0558904467
#>  [816] -0.0495526875  0.0203758515 -0.1106013737 -0.1599595702  0.0244928342
#>  [821]  0.1564363678  0.1243193771 -0.2901952420  0.0842731785  0.0633311211
#>  [826] -0.2380819176  0.1965727607 -0.0833290248  0.0256353644  0.0617485989
#>  [831]  0.2529769753  0.1379796265 -0.0643832993 -0.0950782381 -0.0655598412
#>  [836]  0.1666930219 -0.1768179186 -0.1048495506 -0.0413769069 -0.0466121813
#>  [841]  0.1035083870  0.0164478535  0.0344736640  0.0260974971 -0.0115372995
#>  [846]  0.0493203503 -0.0087720631 -0.0671133262 -0.2966037158 -0.2181289376
#>  [851] -0.1913570472 -0.2464363585  0.0253999024  0.0650434201 -0.0001410523
#>  [856] -0.1788480558 -0.0228862096  0.0462187170  0.1273625377  0.0366973471
#>  [861] -0.0823776776 -0.1091946334  0.0692828253 -0.1259524996  0.0992348635
#>  [866]  0.3456368561 -0.1070117359  0.1936331972  0.0019781465  0.2304007947
#>  [871]  0.1439033354  0.0645788816  0.1464132164 -0.1145379803 -0.0249369088
#>  [876] -0.0376379005 -0.0070242671 -0.2803666426 -0.0859242915 -0.1448540926
#>  [881]  0.1126744434 -0.0682495330 -0.1006437066 -0.1041729036  0.1696969874
#>  [886]  0.1264023974  0.1029439659 -0.0828198915  0.0893596065  0.0181319678
#>  [891] -0.2054922506 -0.1762365058 -0.1282800980 -0.1236219619 -0.0545501437
#>  [896]  0.0798887152 -0.2229527309  0.0193158838  0.0772645484 -0.0435382775
#>  [901] -0.0917186020  0.1048841457  0.1364892895  0.0277806432 -0.0176012731
#>  [906]  0.0388155149 -0.1271253077  0.0279287348 -0.1083953398 -0.1715847518
#>  [911] -0.0678964813 -0.1391097776  0.0578639990  0.1423963585 -0.2187998809
#>  [916] -0.0711486803 -0.0242829398 -0.0043209866 -0.0974225355 -0.0511505789
#>  [921]  0.0696145148  0.0920242463  0.2031688606  0.0039449452  0.0166857300
#>  [926]  0.1838523759  0.0920948379 -0.1430502195 -0.1603618630  0.1528000017
#>  [931] -0.0861983552  0.0782704690  0.2236808725  0.0600734124 -0.0633458748
#>  [936] -0.0862322987 -0.1439066877  0.0581923681  0.0611216919 -0.0072234426
#>  [941] -0.1185127143  0.2208497001 -0.0622983051 -0.0053987550 -0.0863181864
#>  [946] -0.1232467292 -0.1230316986  0.0699355704  0.0349683432 -0.0409354077
#>  [951] -0.2660448878  0.1488618324 -0.0575890299 -0.0011630167  0.0393909787
#>  [956]  0.1744714584  0.0235401075 -0.1184675744 -0.0497913675 -0.1309653856
#>  [961] -0.2635958766  0.0729398402 -0.0331736786 -0.0483422745  0.0904156095
#>  [966] -0.1420863624  0.0058905959 -0.0136052068  0.0344899518  0.0596747111
#>  [971] -0.0481425865  0.1325081272 -0.1077937998 -0.0904767058 -0.0556709797
#>  [976] -0.0888455139 -0.1075716205 -0.1232935663 -0.2039696061  0.0506135737
#>  [981] -0.0700493360 -0.1371508521 -0.1801653895  0.0282436944 -0.1211348896
#>  [986]  0.0209593261  0.0009772678  0.0195689802 -0.1187352385 -0.0773286969
#>  [991] -0.1234735920 -0.0246584430  0.0750444758 -0.1639395526  0.0218826000
#>  [996] -0.0677713140 -0.1443302715  0.0630014128  0.0896590254 -0.0862821439
#> [1001]  0.1275125719 -0.2398598552  0.0715139638  0.2337764875 -0.1073015892
#> [1006]  0.0350186543 -0.2271577725 -0.0926466523 -0.0140216727  0.0309056427
#> [1011]  0.2070461775 -0.1951898614 -0.0004751766  0.1721056980 -0.0608051400
#> [1016] -0.0016421225 -0.0888845149 -0.0404272006  0.0193838448 -0.2289336786
#> [1021] -0.0858961044  0.1168949077 -0.1310674928  0.2241502158  0.2206279220
#> [1026] -0.0971668375 -0.0747420268 -0.1185573303  0.0842700142  0.0046045413
#> [1031] -0.0889676531 -0.0134609612  0.0028138635  0.0906085998  0.0866009182
#> [1036] -0.1851405193 -0.1234363850 -0.0629670971 -0.0610317960  0.1054349236
#> [1041] -0.1666029848 -0.2104549089 -0.0158763163 -0.2316390282  0.0915723226
#> [1046] -0.0465681160 -0.0123859721 -0.1559146066  0.1550148267  0.0801079011
#> [1051] -0.0165141803 -0.1687108003  0.1573401976  0.0043383266 -0.0047613439
#> [1056] -0.0893186307  0.1751585091 -0.0989654549 -0.0948907342 -0.0411538999
#> [1061] -0.0606985242 -0.1027751450  0.1344502442 -0.1169233530 -0.0186819322
#> [1066]  0.2786154364  0.0555906630  0.1071601329  0.1919628033 -0.1368718789
#> [1071]  0.0100005236  0.0499878211  0.0011086076 -0.0148784919  0.1242739907
#> [1076]  0.0204827626  0.0693341287 -0.1403812980  0.0270454725 -0.1029750625
#> [1081] -0.0899480937 -0.1023498380  0.1640459784  0.0469823540 -0.0520574748
#> [1086]  0.0304403692  0.0615814933 -0.0932197549 -0.0932736637 -0.0140769016
#> [1091] -0.1290957638  0.0001587242 -0.1457609819 -0.0078516905  0.0425954703
#> [1096]  0.1278176707 -0.0305662986  0.0734679015  0.0368057343 -0.0658092735
#> [1101] -0.0093029653  0.2568308622  0.0157921909 -0.0720935348 -0.0259277723
#> [1106] -0.0292620717  0.1102756845  0.0209967313  0.0838090343 -0.1318858747
#> [1111] -0.0077486049  0.0788275233  0.0289672918 -0.2266199385 -0.0126275707
#> [1116]  0.0362398984 -0.1478911277 -0.0385813054 -0.0446998809 -0.1344869620
#> [1121] -0.0658284271 -0.1263501781 -0.1285148463 -0.0913273446  0.1114802500
#> [1126]  0.1864173495  0.1299378528 -0.0802267002 -0.1433435218 -0.0358701346
#> [1131] -0.1508227765 -0.0009696863  0.1428160600  0.0011169615  0.0959596280
#> [1136] -0.0503055336 -0.0791778186 -0.0657802585 -0.2178765821 -0.2438357359
#> [1141] -0.1216101713  0.1281320651  0.2086967129 -0.0464596636 -0.0551236485
#> [1146] -0.0605815254  0.0673699571 -0.1628543291 -0.0226769701  0.0777143097
#> [1151]  0.1097445593  0.0954790380 -0.0171341983 -0.0844743230  0.0335634899
#> [1156]  0.0235618317  0.2686811567 -0.0244363479  0.2071629810 -0.0547699094
#> [1161] -0.0884261474  0.0029114617  0.0938873761  0.0553653810 -0.1004195371
#> [1166]  0.0990810039  0.1063624294  0.1113009518 -0.0930907531 -0.1132347582
#> [1171]  0.1674196463  0.0343310601 -0.1259903004 -0.0013857294 -0.1651101824
#> [1176]  0.0306631911 -0.0577357490 -0.0415363570 -0.0318977509 -0.0338946955
#> [1181]  0.0394097847 -0.0329991468  0.0304174102  0.1546312304  0.0240864126
#> [1186] -0.2042856172  0.0049324412 -0.0317556919  0.0558237985 -0.1046134192
#> [1191]  0.0228626582  0.0226759171  0.0437814251 -0.1135949037  0.1024754799
#> [1196] -0.1021271109 -0.0516827193 -0.1176194308  0.0602945597  0.0407256628
#> [1201] -0.1773867957  0.0470489501 -0.1002825290 -0.0367529909 -0.1035664264
#> [1206]  0.1094358203  0.0030173628 -0.0397729943 -0.0308755279 -0.0027361148
#> [1211]  0.0178646904  0.0017609938 -0.0670346443  0.0719520517 -0.0787214358
#> [1216]  0.0867935341  0.1499233164 -0.0544361472  0.0864984620 -0.0867012527
#> [1221]  0.1237411983 -0.2293633933 -0.0284702836 -0.0763790241  0.0197793293
#> [1226] -0.1265029472 -0.2830406512 -0.2358060736 -0.0265238264 -0.0480070317
#> [1231] -0.0400355663 -0.1337869276  0.0507961134 -0.1615569704 -0.0801732229
#> [1236]  0.1319654010 -0.0850298574 -0.2024863287 -0.0436249716 -0.0405610690
#> [1241] -0.0966534116 -0.0564933119  0.0207048974  0.0786864775 -0.1742600833
#> [1246] -0.0703196299  0.0922132252  0.2200446648 -0.1340101236 -0.2001613248
#> [1251] -0.0703995527 -0.0016791005 -0.0314317791 -0.0475607286  0.0888471699
#> [1256]  0.1676344230 -0.0611367358  0.2120560644 -0.0468189264 -0.0981902290
#> [1261]  0.0856278425 -0.0265292701  0.1244178100  0.0747364617  0.1053207631
#> [1266] -0.0451228805  0.1955709597 -0.2066139161  0.1724383949 -0.0270658746
#> [1271] -0.0935892327  0.1588252702 -0.0848976356  0.0131995215 -0.0628385624
#> [1276] -0.0241777464 -0.0502395160 -0.0451679805 -0.1888543784 -0.2477834541
#> [1281] -0.0677057982 -0.0481157948 -0.0463778956  0.0018929255  0.0023472036
#> [1286] -0.0911424616  0.0014451608 -0.0023308389  0.0216576807 -0.0179058102
#> [1291] -0.1694113036 -0.1291153684  0.0956824615 -0.1793850393 -0.1083917774
#> [1296] -0.0027912150  0.1057581996 -0.1040905470 -0.0292283433 -0.0647422654
#> [1301] -0.0349178695 -0.0199616493 -0.1761153175  0.0002470467 -0.1922743614
#> [1306]  0.0119673420  0.2953939198  0.0725379237 -0.0502435883 -0.0623227168
#> [1311]  0.1103912553 -0.1848522297  0.0528879249  0.1203774392  0.0882539073
#> [1316]  0.1366061273 -0.0880141482 -0.1646203622  0.0520044239 -0.0163267168
#> [1321]  0.0249976157  0.0486733198 -0.0431287704 -0.2084879682  0.0716195150
#> [1326]  0.1820396661 -0.2187599351 -0.1171372059  0.0564583296 -0.1569266033
#> [1331]  0.1778659480 -0.0451926195 -0.1551368389  0.0383306270 -0.0491906130
#> [1336] -0.1098927579  0.0633404818 -0.0105195792  0.1392403276 -0.1946386334
#> [1341] -0.0388798243 -0.0873647996  0.0056799105 -0.0172951931  0.0172384655
#> [1346] -0.1240797262  0.0280678595 -0.1521023641  0.1620418840 -0.1476542482
#> [1351] -0.0771701161  0.0205549779 -0.0570517859  0.0418085242 -0.1590700312
#> [1356]  0.0094505182 -0.0401412474 -0.0983152595  0.0998841038  0.1017705058
#> [1361] -0.1304467304 -0.0613215534  0.0881999493 -0.1738396878  0.1209879988
#> [1366] -0.2039931849 -0.0118793904  0.0438666601 -0.0009272782  0.0536990414
#> [1371] -0.0711064475  0.0509744703  0.1440466133 -0.1372136706 -0.0099961650
#> [1376]  0.0467299918  0.1113532754 -0.1248331796  0.0363964026  0.2556637406
#> [1381]  0.2392300185  0.0988183092 -0.0689220121 -0.1212235972  0.1222006238
#> [1386] -0.3180201925  0.2447290489 -0.1464049100 -0.0645423914  0.1642899103
#> [1391]  0.0159609890 -0.0752082543 -0.1332547684  0.0511768775  0.1437431599
#> [1396] -0.0630441946 -0.0448820831 -0.0450625514 -0.1859465732 -0.1122638065
#> [1401] -0.1487315744 -0.0183488837  0.0169772384  0.1082035093 -0.1821877303
#> [1406] -0.0338320155 -0.1286543818 -0.0209141748 -0.1404363590  0.0201670068
#> [1411] -0.0190953127 -0.1547103551  0.0717317827 -0.1664627271  0.1030628854
#> [1416] -0.0925591183 -0.2548660786  0.1086621432  0.0501271784 -0.0530513577
#> [1421] -0.0239029208 -0.1906057376 -0.0482348884  0.0652195998 -0.1375973989
#> [1426]  0.0189192539 -0.1635576916 -0.0998006582 -0.1160078590  0.1075946023
#> [1431] -0.0525256158 -0.1311349140 -0.1796726361 -0.0155062296 -0.1571523819
#> [1436]  0.0025104487 -0.0297287095  0.0036962358 -0.0659211521  0.1716566594
#> [1441] -0.0133387999 -0.1079742767 -0.1773412701 -0.2203874826  0.1175857221
#> [1446] -0.1337477871  0.0843276130 -0.2389200286 -0.1040119515  0.0875135834
#> [1451] -0.1255518868 -0.0557028490 -0.0662452759 -0.0999982026  0.0344418657
#> [1456]  0.1269450473  0.0675891060 -0.1212277846  0.0022864964  0.0073890989
#> [1461] -0.0735030068  0.1662639629 -0.0854075245  0.0977473389  0.0173497351
#> [1466] -0.2210680196  0.0247365245 -0.1092364597  0.1460191053 -0.2477364218
#> [1471] -0.1108681635 -0.1336020166  0.1795345450  0.1116164799 -0.0583952539
#> [1476]  0.1529938806 -0.0144613532 -0.0076882693  0.1030474682  0.0014595408
#> [1481] -0.0003500750  0.0662361581  0.0723323628 -0.0179523445  0.0789283999
#> [1486] -0.1635146136 -0.1374827659  0.0944178593 -0.1791855150  0.0514046685
#> [1491]  0.1763298155 -0.0229999245  0.1158579195 -0.1650028297  0.0522132408
#> [1496] -0.1083272315 -0.0914032560 -0.1160418357  0.0353277394 -0.0440764565
#> [1501] -0.0806248210  0.0716930984  0.0013124263  0.0249330614 -0.1080625793
#> [1506] -0.2589241997  0.2249161694  0.0505848921 -0.3174085539  0.0949429776
#> [1511] -0.0083981174 -0.0031331210  0.2320465640  0.1217818673  0.1898370933
#> [1516]  0.0542053957 -0.0654487322  0.0870376150 -0.1058793179 -0.1453013628
#> [1521]  0.0587975313 -0.0703381265  0.3179411712  0.0064599367  0.0483878611
#> [1526] -0.3046781881 -0.0395406852  0.0323834825  0.0314728634 -0.0298996998
#> [1531] -0.2892791870  0.0768770437 -0.1640782021  0.0032564812 -0.1021590296
#> [1536]  0.2410490460 -0.1230430855  0.0506437463 -0.1405195058 -0.0059927120
#> [1541] -0.0504173415 -0.0471274011 -0.1736461164 -0.0255931802 -0.0524809651
#> [1546]  0.1965384269  0.0065867553 -0.0893868909  0.0170753540  0.0064498911
#> [1551] -0.1030731150  0.1168908291 -0.1792427408 -0.0857070516 -0.0557818416
#> [1556] -0.0473785815  0.0717260133 -0.0779482931  0.0796002479 -0.3056094777
#> [1561] -0.0993707355  0.0458856440 -0.1185352707 -0.1687477616  0.1794336945
#> [1566] -0.0709608610 -0.0851940621  0.0125716806 -0.1243984275 -0.0514837569
#> [1571]  0.0757020590 -0.0427037048  0.2423472518 -0.0530946922  0.1819950438
#> [1576] -0.1988334164  0.0202708039  0.1186214022 -0.1998539516  0.0578597626
#> [1581] -0.0336843825  0.0267445403  0.1309220228 -0.2263342259  0.0284160579
#> [1586]  0.0010297838 -0.0136570548  0.1159669449  0.0136858026 -0.0334869258
#> [1591] -0.0501891611 -0.0677227336  0.0147216003  0.2797036597  0.0035794275
#> [1596]  0.0281079965  0.1013694725 -0.1763521647 -0.0898613121 -0.0177610390
#> [1601]  0.1345635132 -0.0834129928 -0.1707209997  0.0382682335 -0.0748484536
#> [1606] -0.0875052908 -0.1172680311  0.1222946100 -0.0394149888  0.0000000000
#> [1611]  0.0918370526 -0.2316958792 -0.1283023687  0.0105814672 -0.0536033517
#> [1616] -0.0463128349  0.0898462088 -0.0344243564  0.1753812699  0.0268084936
#> [1621]  0.1108427528 -0.0103117611  0.1025184030 -0.0706776424 -0.0173086706
#> [1626]  0.1259253021 -0.0704601619  0.0680946969 -0.0783298215  0.0928654281
#> [1631]  0.0855325883  0.1612180313  0.0430897754 -0.2513714285  0.1066975979
#> [1636]  0.0012674216 -0.0619438377 -0.0588997305 -0.1108989936  0.0211734605
#> [1641] -0.0775216289 -0.1901485281  0.0278508358 -0.1007379694  0.0756620462
#> [1646]  0.2433131066 -0.0685093843 -0.0016799628  0.0302533090 -0.2057110505
#> [1651] -0.0119600918 -0.1040637450 -0.1082671755 -0.0848357572 -0.0596453762
#> [1656] -0.1369648753 -0.2298857206 -0.1636013511  0.0270402567 -0.0470546421
#> [1661] -0.1112989586  0.0061191685  0.0029479654  0.0482336862  0.0858228525
#> [1666] -0.0447540008  0.1582222823  0.0238811637 -0.2214286238 -0.0466074826
#> [1671]  0.0550293460  0.0531886104  0.0572746764  0.0551405632  0.0870547516
#> [1676] -0.1126794015 -0.1303932584  0.1263334192 -0.2444907238 -0.0771276338
#> [1681] -0.0067617700 -0.1069632741 -0.0596390956  0.1743644044  0.0634801658
#> [1686]  0.1107565057 -0.0361667757 -0.0544911396  0.0668506347 -0.0838501300
#> [1691]  0.0443894585 -0.0738998398  0.1736490584  0.0320775655 -0.0168004971
#> [1696] -0.0885485252 -0.0508210370 -0.1283003929  0.0714995277  0.1639147014
#> [1701]  0.0078434835 -0.0205337643 -0.0789049496  0.0185327263  0.2494233121
#> [1706] -0.0938839559 -0.0126494849  0.0201634262  0.0177738462  0.1820950154
#> [1711] -0.0063400344 -0.1587921246  0.1513260163  0.0854500450 -0.0356500634
#> [1716]  0.1704516369 -0.1143968697 -0.0641082725  0.0587836061 -0.1921535135
#> [1721] -0.0025540042  0.0967239936 -0.0610213150  0.2075580154 -0.0457421026
#> [1726]  0.1517211429  0.0724867555  0.0911549860 -0.0609141585 -0.1060129284
#> [1731]  0.0693608110  0.2331380175 -0.1303910328 -0.0195884172  0.1241924738
#> [1736] -0.1136441174 -0.0137919745  0.0211024320 -0.1205456725 -0.0796432189
#> [1741] -0.0050923466 -0.1088963754  0.0707711206 -0.1199633246 -0.0305272046
#> [1746]  0.2796884860 -0.0315401841  0.2378053253 -0.2256149442 -0.1043567774
#> [1751] -0.1836781908  0.2472121841 -0.0562806477  0.0964697199  0.0377714661
#> [1756]  0.0767565640  0.0246621900 -0.1539620419  0.0485139917  0.0250013561
#> [1761]  0.1267101608  0.0166839217 -0.1705578920  0.3421397954 -0.0884433207
#> [1766] -0.0264527355 -0.1625885777 -0.0947696212  0.0577472442 -0.2277959291
#> 
#> $mod2_unimodal$ref_properties$feature$prop_outliers_feature
#> [1] 0
#> 
#> 
#> 
#> $mod2_unimodal$sim_properties
#> $mod2_unimodal$sim_properties$cell
#> $mod2_unimodal$sim_properties$cell$library_size
#>  [1] 24 26 21 24 24 17 21 20 22 27 14 23 27 20 18 17 20 22 23 19 25 23 27 20 23
#> [26] 28 26 15 20 18 26 20 32 21 25 19 28 30 23 31 27 31 34 36 29 33 33 34 29 35
#> [51] 34 34 26 35 26 32 24 29 25 28 31 25 37 27 36 35 27 30 30 32 29 31 29 21 28
#> [76] 31 31 27 31 26
#> 
#> $mod2_unimodal$sim_properties$cell$log_library_size
#>  [1] 3.218876 3.295837 3.091042 3.218876 3.218876 2.890372 3.091042 3.044522
#>  [9] 3.135494 3.332205 2.708050 3.178054 3.332205 3.044522 2.944439 2.890372
#> [17] 3.044522 3.135494 3.178054 2.995732 3.258097 3.178054 3.332205 3.044522
#> [25] 3.178054 3.367296 3.295837 2.772589 3.044522 2.944439 3.295837 3.044522
#> [33] 3.496508 3.091042 3.258097 2.995732 3.367296 3.433987 3.178054 3.465736
#> [41] 3.332205 3.465736 3.555348 3.610918 3.401197 3.526361 3.526361 3.555348
#> [49] 3.401197 3.583519 3.555348 3.555348 3.295837 3.583519 3.295837 3.496508
#> [57] 3.218876 3.401197 3.258097 3.367296 3.465736 3.258097 3.637586 3.332205
#> [65] 3.610918 3.583519 3.332205 3.433987 3.433987 3.496508 3.401197 3.465736
#> [73] 3.401197 3.091042 3.367296 3.465736 3.465736 3.332205 3.465736 3.295837
#> 
#> $mod2_unimodal$sim_properties$cell$tmm_factor
#>  [1] 1.0711448 0.9887491 1.2241655 1.0711448 1.0711448 1.5122045 1.2241655
#>  [8] 1.2853738 1.1685216 0.9521287 1.8362483 1.1177164 0.9521287 1.2853738
#> [15] 1.4281931 1.5122045 1.2853738 1.1685216 1.1177164 2.2447861 1.0282990
#> [22] 1.1177164 0.9521287 1.2853738 1.1177164 0.9181242 0.9887491 1.7138317
#> [29] 1.2853738 1.4281931 0.9887491 1.2853738 0.8033586 1.2241655 1.0282990
#> [36] 1.3530251 0.9181242 0.8569159 1.1177164 0.8292734 0.9521287 0.8292734
#> [43] 0.7561022 0.7140966 0.8864647 0.7790144 0.7790144 0.7561022 0.8864647
#> [50] 0.7344993 0.7561022 0.7561022 0.9887491 0.7344993 0.9887491 0.8033586
#> [57] 1.0711448 0.8864647 1.0282990 0.9181242 0.8292734 1.0282990 0.6947967
#> [64] 0.9521287 0.7140966 0.7344993 0.9521287 0.8569159 0.8569159 0.8033586
#> [71] 0.8864647 0.8292734 0.8864647 1.2241655 0.9181242 0.8292734 0.8292734
#> [78] 0.9521287 0.8292734 0.9887491
#> 
#> $mod2_unimodal$sim_properties$cell$effective_library_size
#>  [1] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#>  [9] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [17] 25.70748 25.70748 25.70748 42.65094 25.70748 25.70748 25.70748 25.70748
#> [25] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [33] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [41] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [49] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [57] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [65] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> [73] 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748 25.70748
#> 
#> $mod2_unimodal$sim_properties$cell$zero_fraction_cell
#>  [1] 0.6166667 0.6500000 0.6833333 0.6666667 0.6666667 0.7166667 0.7000000
#>  [8] 0.7000000 0.6833333 0.6000000 0.7833333 0.6333333 0.6166667 0.6833333
#> [15] 0.7000000 0.7500000 0.6833333 0.6666667 0.6666667 0.7333333 0.6500000
#> [22] 0.6500000 0.5833333 0.7000000 0.6500000 0.5833333 0.6666667 0.7500000
#> [29] 0.7000000 0.7000000 0.6166667 0.6666667 0.5166667 0.6833333 0.6333333
#> [36] 0.7500000 0.5833333 0.5666667 0.6500000 0.5500000 0.6000000 0.5833333
#> [43] 0.5166667 0.5000000 0.6166667 0.5000000 0.5166667 0.5166667 0.6000000
#> [50] 0.5500000 0.5000000 0.5666667 0.6333333 0.5000000 0.6333333 0.5833333
#> [57] 0.6833333 0.5833333 0.6333333 0.5833333 0.5166667 0.6333333 0.4500000
#> [64] 0.5833333 0.4833333 0.5500000 0.6000000 0.5333333 0.5833333 0.5166667
#> [71] 0.5500000 0.5833333 0.5666667 0.7000000 0.6333333 0.5500000 0.5666667
#> [78] 0.6000000 0.5333333 0.6666667
#> 
#> $mod2_unimodal$sim_properties$cell$detection_freq_cell
#>  [1] 0.3833333 0.3500000 0.3166667 0.3333333 0.3333333 0.2833333 0.3000000
#>  [8] 0.3000000 0.3166667 0.4000000 0.2166667 0.3666667 0.3833333 0.3166667
#> [15] 0.3000000 0.2500000 0.3166667 0.3333333 0.3333333 0.2666667 0.3500000
#> [22] 0.3500000 0.4166667 0.3000000 0.3500000 0.4166667 0.3333333 0.2500000
#> [29] 0.3000000 0.3000000 0.3833333 0.3333333 0.4833333 0.3166667 0.3666667
#> [36] 0.2500000 0.4166667 0.4333333 0.3500000 0.4500000 0.4000000 0.4166667
#> [43] 0.4833333 0.5000000 0.3833333 0.5000000 0.4833333 0.4833333 0.4000000
#> [50] 0.4500000 0.5000000 0.4333333 0.3666667 0.5000000 0.3666667 0.4166667
#> [57] 0.3166667 0.4166667 0.3666667 0.4166667 0.4833333 0.3666667 0.5500000
#> [64] 0.4166667 0.5166667 0.4500000 0.4000000 0.4666667 0.4166667 0.4833333
#> [71] 0.4500000 0.4166667 0.4333333 0.3000000 0.3666667 0.4500000 0.4333333
#> [78] 0.4000000 0.4666667 0.3333333
#> 
#> $mod2_unimodal$sim_properties$cell$cell_cor
#>    [1] -1.260383e-02  1.294424e-01 -1.195173e-01 -8.686826e-02  1.196957e-01
#>    [6]  1.459098e-01 -2.047838e-01  5.127215e-02 -3.958609e-02 -5.553142e-02
#>   [11]  7.003365e-02 -3.123322e-02 -5.625094e-02 -1.703454e-01  4.340123e-01
#>   [16]  1.028444e-02  2.197612e-02  1.296646e-01  1.125203e-01  5.848860e-02
#>   [21]  1.250534e-02 -1.947774e-01  3.597886e-01 -3.030284e-01  2.111671e-01
#>   [26]  3.897196e-02 -6.485530e-02  6.241807e-02  2.790486e-01 -6.191714e-02
#>   [31]  2.863934e-01  3.020738e-01 -1.258320e-01 -5.760991e-02  1.446275e-01
#>   [36] -1.587310e-01  2.317374e-01 -1.197440e-01 -1.110503e-02 -1.630840e-01
#>   [41]  8.404660e-02  9.798920e-02  4.229842e-02 -5.536544e-02  7.572420e-02
#>   [46] -2.123905e-01  2.199140e-01 -2.079002e-02  1.674333e-01  1.611226e-01
#>   [51]  2.939051e-02  1.688454e-02  9.750227e-03 -1.938963e-02  4.398059e-03
#>   [56]  1.781798e-01  1.775335e-01 -2.431992e-01  2.745888e-02 -1.988860e-02
#>   [61]  1.528916e-01 -1.403530e-01 -1.984348e-04 -1.226499e-01 -1.142376e-01
#>   [66] -7.378323e-02  3.923775e-02 -7.327372e-02  9.159299e-03  4.416292e-02
#>   [71] -1.660460e-01 -1.243467e-01  1.773820e-01  8.460663e-02 -1.061081e-01
#>   [76]  1.599797e-02  9.341657e-03  6.270857e-02  1.643454e-02  1.705163e-01
#>   [81]  1.262166e-01  1.229219e-01 -1.270819e-01 -8.502045e-02  2.417828e-02
#>   [86] -8.065657e-02  1.673604e-01 -4.699199e-02  8.457029e-02 -9.731301e-02
#>   [91] -1.458132e-01 -2.200375e-01  2.600215e-02 -7.280051e-02  3.083354e-02
#>   [96]  8.424186e-02 -6.472867e-02  2.393005e-01  4.296710e-02 -1.635360e-01
#>  [101] -4.721749e-02  2.263100e-01  7.963527e-02 -1.730993e-01  3.036953e-02
#>  [106]  2.215905e-02 -4.780334e-02 -1.148795e-01  9.149128e-04  1.583112e-01
#>  [111]  1.038093e-01 -7.551850e-03 -8.985575e-02 -6.439670e-02  5.758264e-02
#>  [116] -9.634497e-03 -4.466021e-02 -2.567610e-01 -1.083804e-01  3.806744e-01
#>  [121]  3.021152e-02 -7.176074e-03 -8.744440e-02  2.691111e-01 -2.108668e-01
#>  [126] -4.274060e-02  1.341355e-01  1.116411e-01  6.180245e-02 -1.351674e-01
#>  [131] -2.957443e-01 -8.888625e-03  1.130719e-01  1.308331e-01 -5.693419e-02
#>  [136]  6.399679e-02 -7.281669e-02  9.991027e-02 -8.910429e-03 -5.283108e-02
#>  [141] -3.892027e-02  2.866018e-01 -5.565301e-02  2.010577e-01  1.103837e-02
#>  [146]  1.704036e-02 -1.273459e-01  6.739549e-02 -4.292910e-02  1.401127e-01
#>  [151] -1.654950e-02 -1.431072e-02 -1.085752e-02  8.888496e-02  1.387594e-01
#>  [156] -3.389508e-02  1.731191e-01 -1.491815e-01  5.979000e-03  5.411964e-02
#>  [161]  1.378427e-01  8.744222e-02  4.909995e-02  1.141541e-01 -1.614584e-01
#>  [166] -1.435106e-01 -1.264422e-02  9.917246e-02  1.375937e-01  3.023187e-01
#>  [171]  1.609847e-01  9.854578e-03  5.772476e-02  1.661342e-02  1.261945e-02
#>  [176] -1.396714e-01 -4.345286e-02  4.217280e-02 -5.285364e-03 -2.746141e-01
#>  [181] -2.423724e-01 -1.641164e-01 -2.085349e-02  2.317764e-02 -6.178745e-02
#>  [186]  5.957298e-02  1.183338e-02  2.189037e-01  3.638709e-02  1.040948e-01
#>  [191]  1.165272e-01  1.669827e-01  1.898086e-01 -9.468513e-02 -6.648135e-02
#>  [196] -8.413519e-03  1.063769e-01  2.113250e-01  1.352187e-01  1.574017e-01
#>  [201]  1.442677e-01 -9.037037e-02 -5.216690e-02 -4.282102e-02 -1.485234e-02
#>  [206] -6.247785e-02 -1.798531e-01  1.110309e-01 -2.628911e-02 -1.311882e-01
#>  [211]  3.449272e-03 -6.397487e-03 -2.189107e-01 -1.367700e-01 -1.280831e-01
#>  [216] -7.198827e-02 -1.248332e-01  1.073449e-01 -7.907680e-04 -3.360718e-02
#>  [221]  1.703451e-01  1.839549e-01  2.129226e-02 -8.493195e-02 -1.734215e-01
#>  [226] -1.203666e-02  4.691136e-02 -1.245087e-01 -2.943970e-02  3.825261e-02
#>  [231] -1.149227e-03 -2.245905e-01 -3.861244e-02  8.857614e-02  2.760646e-01
#>  [236]  2.175349e-01  3.314921e-02 -1.544641e-01 -7.374763e-02  1.023108e-02
#>  [241]  3.786119e-02  3.272972e-01 -2.871194e-02  2.369711e-02  2.159899e-01
#>  [246]  1.895463e-02 -1.725186e-01 -1.451881e-01  1.193743e-01  7.763896e-02
#>  [251]  1.097912e-01  2.314151e-02  2.354284e-02 -1.090346e-01 -1.465477e-01
#>  [256]  1.213360e-01 -1.924716e-02  6.212625e-03 -1.285894e-01 -4.387449e-02
#>  [261] -8.601514e-02  5.248798e-02 -6.858733e-02  3.458239e-01 -1.190212e-01
#>  [266]  1.844342e-01  9.617304e-03  1.239322e-01 -2.887286e-01 -3.165416e-02
#>  [271] -5.691971e-02 -1.294867e-01  3.063309e-02  9.860892e-02  1.146502e-01
#>  [276]  2.920613e-01  2.411507e-01  8.067729e-02  1.432784e-01  1.449183e-01
#>  [281] -1.075237e-01 -2.935348e-01 -1.935388e-01 -1.554026e-01  2.615465e-02
#>  [286] -6.526131e-02 -1.171094e-01  2.925525e-02 -7.144654e-02 -1.048601e-01
#>  [291]  9.545528e-02  2.721443e-01 -9.629381e-02  1.095746e-01  1.376899e-01
#>  [296]  3.020384e-02 -1.072612e-01 -1.152697e-01  8.721054e-02 -2.457220e-02
#>  [301]  2.224265e-01 -1.126120e-01  1.071188e-01 -6.711984e-02 -5.149708e-02
#>  [306]  4.089964e-02 -2.577146e-01 -9.015513e-02  3.477814e-01  1.735305e-01
#>  [311] -1.860857e-01 -5.113528e-02  1.618889e-01 -8.094917e-02 -2.142831e-01
#>  [316] -2.653504e-02 -1.010674e-01  1.812722e-01 -5.815985e-02 -7.314099e-02
#>  [321]  1.087104e-01  1.257683e-01  1.487709e-01  1.503862e-01  1.052503e-01
#>  [326] -1.349552e-01 -5.535985e-02  2.743716e-01  6.226680e-02 -1.297514e-01
#>  [331] -1.685161e-01 -8.505099e-02 -5.624438e-02  9.053290e-02 -3.021938e-01
#>  [336]  1.828579e-01 -1.013693e-01  1.483750e-01 -1.738881e-01 -1.284298e-01
#>  [341] -1.782161e-01 -3.024392e-02 -6.492686e-02  1.258636e-01 -5.006850e-02
#>  [346] -8.034511e-02  2.893413e-02 -5.934377e-02  6.482508e-02  7.616598e-02
#>  [351] -3.026384e-02 -2.696197e-02 -5.745703e-02 -1.773194e-01  1.436161e-01
#>  [356] -3.645035e-02 -2.567586e-02 -8.593123e-02  2.015281e-01  1.247158e-01
#>  [361]  1.947856e-01 -1.333585e-01  2.822136e-02  2.152781e-01 -9.398356e-02
#>  [366] -4.558245e-02  4.527867e-02 -5.039903e-02  1.382247e-01 -3.369515e-02
#>  [371] -2.604690e-01 -1.962255e-01 -4.316787e-02 -1.755353e-02 -2.441330e-01
#>  [376]  2.498244e-01  9.487069e-02 -3.894437e-02  2.880034e-01  3.501690e-02
#>  [381]  2.205731e-02 -1.364873e-02 -1.423999e-01 -2.081924e-01  9.715820e-02
#>  [386] -1.043896e-01  1.870209e-01 -3.298696e-02 -3.043973e-03  1.367216e-02
#>  [391]  1.222379e-01  3.094985e-02 -1.011165e-01  8.600092e-02  3.987925e-03
#>  [396] -6.580702e-02  2.478185e-02 -4.723794e-02  1.665686e-01 -2.499482e-01
#>  [401] -1.901892e-01 -3.266799e-02  1.941627e-01  3.405946e-02  9.514352e-02
#>  [406]  8.275276e-02 -4.419369e-02  1.584414e-01 -9.001543e-02  7.592618e-02
#>  [411]  6.004176e-02  1.710214e-01  5.017550e-02  3.114771e-01  3.946465e-03
#>  [416] -1.120666e-01  2.016632e-01  6.091950e-03  1.253722e-03  2.956854e-01
#>  [421] -1.302969e-01 -8.437784e-02  8.621811e-02  4.688827e-02  9.316178e-03
#>  [426]  1.011707e-01  2.016049e-01  2.094775e-01  2.136231e-01 -1.424705e-01
#>  [431] -1.551611e-01 -1.167671e-01  1.635234e-01 -6.253232e-02 -2.120776e-01
#>  [436] -1.340109e-01  7.140925e-02  2.475320e-01 -6.301687e-02 -7.480590e-02
#>  [441] -8.271802e-02 -4.783732e-02 -2.103902e-02 -8.490189e-02 -1.334360e-02
#>  [446]  7.975660e-02  1.186885e-01  1.501996e-01 -1.000120e-01  1.035344e-01
#>  [451]  1.057202e-01 -2.815060e-01 -2.432271e-02 -2.718222e-01 -1.222539e-01
#>  [456]  7.345876e-02  2.968238e-02 -2.005386e-01 -2.178341e-01 -6.958206e-02
#>  [461]  3.139513e-02  2.052152e-01 -7.309596e-04  1.514261e-02 -9.571382e-02
#>  [466] -1.717353e-01  7.723698e-02 -1.886643e-01  6.489379e-02  2.225228e-01
#>  [471]  2.374564e-01  2.204464e-02  4.218923e-02 -2.208677e-01 -7.429671e-02
#>  [476]  2.423295e-01  1.428247e-01 -5.831771e-02  1.204235e-02  1.069180e-01
#>  [481]  1.577011e-02 -1.984457e-01 -4.244315e-02 -1.510551e-01 -2.772654e-02
#>  [486]  3.109254e-02  5.030515e-02 -4.569264e-03  5.834985e-02 -1.779056e-01
#>  [491] -2.799091e-01 -1.335452e-01 -1.407949e-01 -1.461606e-01  3.869559e-02
#>  [496]  5.682238e-02  2.055890e-01  1.529577e-01 -1.088061e-02  7.903031e-02
#>  [501] -2.352999e-01 -1.094314e-01  2.672863e-01  2.079878e-01  1.186785e-01
#>  [506]  1.676712e-01  6.935080e-02  8.493083e-02  2.880298e-01 -1.056615e-01
#>  [511] -1.081523e-01 -2.101315e-01 -1.675683e-02 -6.539554e-02 -1.329647e-01
#>  [516] -1.920124e-01  2.444816e-01 -6.040995e-02 -2.773207e-01 -3.112617e-02
#>  [521] -6.241307e-02  9.188860e-02  7.269001e-02  9.043942e-02  2.398155e-01
#>  [526]  1.523816e-01  1.171618e-02 -6.459281e-02  1.222124e-01  5.910871e-02
#>  [531]  6.687885e-02  2.651571e-01  9.105975e-02  9.835845e-02  1.258481e-01
#>  [536] -9.666564e-02  1.666565e-01 -6.037655e-02  5.699969e-03  8.586490e-02
#>  [541]  4.887724e-02 -1.431574e-01  1.371323e-01  1.952899e-01  2.407574e-03
#>  [546] -1.615181e-01  1.436937e-02 -9.222121e-02  3.585991e-02 -1.415446e-02
#>  [551] -1.788521e-02  2.550610e-02  9.067948e-02 -1.673800e-02 -9.871919e-02
#>  [556]  1.283227e-01 -3.587056e-02 -5.936193e-02 -8.178279e-03  1.485734e-01
#>  [561] -1.620190e-01 -2.143680e-02  1.185175e-01 -1.133296e-01 -1.685179e-01
#>  [566] -8.003687e-02 -1.513594e-01  1.667708e-01 -2.052290e-01 -3.567351e-02
#>  [571] -1.418940e-02  2.879126e-01  8.160757e-02  1.407364e-01  1.001823e-01
#>  [576]  4.098005e-02 -1.983656e-01 -2.145435e-01  6.362064e-02 -1.129350e-01
#>  [581] -2.409712e-01  1.424638e-01  3.106952e-02 -4.959193e-02  1.303964e-01
#>  [586] -8.988439e-02 -5.133483e-02  1.417152e-01 -2.361247e-01  1.199980e-01
#>  [591]  3.800169e-02  2.087238e-01  2.646642e-01  1.996867e-01  1.994167e-02
#>  [596] -2.004870e-03  1.530577e-01  1.388762e-01  5.750604e-03 -7.407544e-03
#>  [601] -3.411107e-02  2.306923e-01  1.498576e-01  1.896378e-01  1.006719e-01
#>  [606]  3.194443e-02 -7.818863e-02  1.456454e-01 -1.134242e-01 -8.205545e-03
#>  [611] -7.738524e-02 -7.572959e-02  2.587568e-02 -5.013858e-02 -9.602421e-02
#>  [616]  2.563357e-01 -1.148255e-01 -7.960129e-02 -1.705509e-01 -1.121222e-01
#>  [621]  4.629121e-02 -1.468002e-02  1.261205e-01  1.094628e-01  7.078110e-02
#>  [626]  1.202922e-01 -1.307288e-01  2.461783e-01 -5.584504e-02 -9.008106e-02
#>  [631]  1.254520e-01  2.937286e-02  1.726329e-02  5.638759e-02  8.252570e-03
#>  [636] -5.251371e-03  1.242869e-01 -3.185664e-02 -4.993119e-02 -8.308141e-02
#>  [641]  5.809323e-02 -2.985543e-01  1.110320e-02  1.892009e-02 -1.834323e-01
#>  [646] -1.183013e-01 -8.583100e-02 -9.295449e-02 -1.464490e-01  3.082662e-02
#>  [651] -5.418934e-02 -4.358553e-02 -2.242754e-01  1.109431e-01  8.899321e-02
#>  [656]  1.048745e-02 -3.226058e-02  2.690087e-02 -8.825669e-03  4.659800e-02
#>  [661] -1.599145e-01  2.485076e-01  3.550384e-01 -3.197035e-03  2.359308e-01
#>  [666] -5.371880e-02  8.118308e-02  1.487840e-01  1.052472e-01  1.144517e-02
#>  [671]  3.756590e-02 -5.410684e-03  1.227802e-01  2.574993e-01  1.569682e-01
#>  [676] -1.001094e-01  6.232104e-02 -1.140883e-01 -7.122832e-02  1.386932e-02
#>  [681] -2.446612e-02 -2.102417e-01 -9.590652e-02 -9.420922e-02 -6.804170e-02
#>  [686] -2.052099e-01  7.019715e-02  2.207862e-01 -2.072275e-01  5.503222e-02
#>  [691] -9.234745e-03  3.047412e-02  3.587021e-02  7.685376e-02  4.779314e-02
#>  [696] -1.125470e-01  9.958270e-02  1.197528e-01  1.174770e-01 -3.464325e-02
#>  [701]  5.186186e-02  4.097153e-02  2.304354e-01 -7.124008e-02 -1.772108e-01
#>  [706] -3.570031e-02  2.643172e-01 -2.022060e-01 -1.477701e-01  3.425151e-02
#>  [711]  1.026740e-01  1.998628e-01 -6.282801e-03  4.006179e-02 -8.309919e-04
#>  [716]  1.891470e-01  1.066176e-01 -2.012512e-03 -1.970026e-01  2.826686e-01
#>  [721] -2.971194e-02  8.771461e-02 -1.660236e-01  4.397625e-02  1.066866e-01
#>  [726]  1.322164e-01  1.154731e-01 -1.471968e-01 -3.174213e-02  1.983350e-01
#>  [731]  2.043569e-01  2.443086e-02  1.244887e-01 -2.376596e-01 -2.879220e-01
#>  [736]  1.566974e-01  1.791019e-02  1.369299e-02  4.737104e-02 -6.063958e-02
#>  [741] -4.822897e-02  3.022572e-02 -2.694235e-02  6.516288e-02  1.464752e-01
#>  [746]  6.052479e-02  1.224666e-01 -9.505902e-02  1.317010e-01  8.729484e-02
#>  [751] -1.916599e-01 -9.268321e-02  9.689164e-03 -1.639284e-01  8.534368e-02
#>  [756] -1.755024e-01 -5.015281e-02 -7.442222e-02  6.293854e-03  5.858429e-02
#>  [761]  2.684355e-02 -4.963157e-02 -5.227673e-02  2.663400e-01 -5.010121e-02
#>  [766]  4.923947e-02 -1.905343e-02 -1.369022e-01 -9.403927e-02 -3.743043e-02
#>  [771]  2.709858e-01 -1.202768e-01  7.220702e-02 -1.779933e-01 -2.502253e-02
#>  [776] -8.704717e-02 -1.119078e-01 -4.099489e-02 -8.046861e-02  2.429768e-03
#>  [781] -1.407534e-01  1.509928e-01 -1.836513e-01  1.992252e-02 -6.971452e-02
#>  [786] -2.057694e-01 -2.043101e-01  2.784920e-01 -2.456948e-01  1.192974e-01
#>  [791]  1.898568e-02 -1.716835e-02  1.921370e-01  2.679361e-01 -2.311508e-02
#>  [796]  2.175344e-02  1.581390e-01 -9.670543e-03 -2.001407e-03  3.603023e-02
#>  [801] -4.832782e-02  2.450155e-01  1.582107e-01 -8.858471e-02 -8.919537e-02
#>  [806] -2.047309e-02 -2.549456e-02  1.898324e-01 -2.880552e-02  1.722523e-01
#>  [811]  4.553682e-02  2.431665e-04 -1.283733e-01 -9.962606e-02  5.555848e-02
#>  [816] -9.357245e-02 -8.455373e-02 -1.680347e-01  2.043476e-01  7.264583e-02
#>  [821] -1.524402e-01 -1.953285e-01 -1.534997e-01 -4.912055e-02 -1.914322e-01
#>  [826]  1.327576e-03 -1.008991e-02 -1.139184e-01  6.486803e-02  3.558374e-02
#>  [831]  6.890706e-03 -2.109628e-01 -9.367425e-02 -1.604369e-02  1.256219e-01
#>  [836] -1.833022e-02  5.414649e-02  2.137084e-01  1.319002e-01 -1.572546e-01
#>  [841] -8.823234e-02 -1.181192e-01 -7.059150e-02 -4.866899e-02  1.091172e-01
#>  [846]  4.940711e-02  9.553050e-02  2.163241e-01  4.480151e-02 -3.909602e-02
#>  [851]  5.869496e-02 -2.096385e-01  1.481770e-01 -2.032139e-01  6.652017e-02
#>  [856]  5.841713e-02  9.700794e-02 -1.268426e-01  2.272320e-02 -1.072536e-01
#>  [861]  5.954580e-02  2.935314e-02 -5.254311e-02  1.448713e-01  3.459472e-02
#>  [866] -1.905316e-01 -5.911714e-02 -1.407829e-01 -2.213200e-01  8.016652e-02
#>  [871]  2.302998e-01 -1.108474e-01  5.161840e-02  7.109071e-02  2.039070e-01
#>  [876] -1.377238e-01 -1.680664e-01  1.052300e-01  7.009381e-02 -1.122835e-01
#>  [881] -1.746855e-01  5.363161e-02 -4.713064e-02 -3.410422e-02 -7.293016e-03
#>  [886]  4.217507e-04  1.737309e-01 -9.119687e-02  1.157422e-02 -6.360618e-02
#>  [891] -1.299449e-01  3.995859e-02 -2.299591e-01  1.457312e-02 -1.157131e-01
#>  [896]  1.062038e-01 -2.550990e-02 -9.849758e-02 -1.192755e-01 -4.878156e-02
#>  [901] -1.697164e-01  1.361784e-01  1.490607e-01 -8.430039e-02 -1.838979e-01
#>  [906]  3.470058e-02 -1.740665e-01  2.460187e-01  1.294052e-01 -1.493871e-01
#>  [911] -1.391629e-01 -2.560225e-01 -2.082170e-01 -1.738066e-01  1.154175e-01
#>  [916]  2.568395e-02 -8.368230e-02 -1.400483e-01 -1.187524e-01 -1.058699e-01
#>  [921]  8.027212e-02 -1.715176e-01  1.141594e-01 -7.347055e-02  9.344000e-02
#>  [926]  4.486112e-02 -1.491489e-02 -2.113474e-01 -1.919538e-02  3.255790e-02
#>  [931] -9.993793e-02 -2.203886e-01 -6.478999e-02 -1.075908e-01  1.442207e-01
#>  [936] -2.362447e-01  2.379463e-02  4.527794e-02 -9.282518e-02  1.164955e-01
#>  [941]  5.274699e-02 -1.981072e-02 -1.054900e-01 -3.281280e-02 -1.355337e-01
#>  [946] -3.511934e-02 -1.428234e-01  1.990401e-01 -1.373820e-02 -1.247019e-01
#>  [951] -4.926205e-02 -1.226022e-01 -8.458243e-02  9.167341e-03 -3.125728e-01
#>  [956] -5.231338e-02  1.947541e-01 -4.020409e-02  4.210297e-03 -8.976925e-02
#>  [961]  6.875797e-03 -2.297232e-02 -2.291718e-01 -3.009416e-02  4.016261e-02
#>  [966]  3.844096e-01  1.128295e-01  5.608024e-02  9.790490e-03  1.066389e-01
#>  [971] -2.721980e-03 -3.256402e-02 -6.951803e-02 -1.880970e-01 -1.966669e-02
#>  [976] -1.337492e-01 -6.502368e-02 -6.102112e-02  3.742131e-02 -1.854925e-01
#>  [981] -1.429348e-01 -1.445098e-02 -1.429504e-02 -1.743656e-01 -1.220579e-01
#>  [986] -3.060666e-02  5.197263e-02  1.828686e-02  1.443631e-02  2.218579e-02
#>  [991] -1.190002e-02  2.321214e-02 -1.790492e-01 -5.586294e-02  8.330623e-02
#>  [996]  4.301120e-02 -1.214832e-01  7.486881e-02  2.321950e-03 -2.204498e-01
#> [1001] -1.977148e-01  8.021917e-02  1.775667e-02 -1.911097e-01 -2.810906e-02
#> [1006] -6.884193e-02 -1.216389e-01 -3.211644e-02 -6.708202e-02 -1.164132e-02
#> [1011] -5.741013e-02 -1.778318e-01 -7.069389e-02 -5.881764e-03  6.806684e-02
#> [1016] -6.646321e-02 -8.205996e-02  6.452311e-02 -2.565255e-04 -2.104425e-01
#> [1021] -1.251115e-02 -1.591658e-02 -7.002502e-02  3.018474e-02 -1.550274e-01
#> [1026] -1.636968e-02 -6.269490e-02 -4.266465e-02  1.889721e-01  2.405568e-01
#> [1031] -2.079351e-01 -1.371950e-01 -1.049117e-01 -1.064958e-01  5.289365e-02
#> [1036] -1.958676e-01  4.821395e-02  2.364644e-02 -9.124763e-02 -5.236846e-02
#> [1041] -2.262727e-01 -2.110087e-02 -8.597265e-02 -7.159652e-02 -2.442310e-02
#> [1046]  4.717363e-02 -1.418760e-01 -2.511115e-01  1.220527e-01 -4.816666e-02
#> [1051]  6.483243e-02  5.179217e-02 -1.539403e-02  1.576653e-01  1.345544e-01
#> [1056] -7.186271e-02  1.394409e-01 -5.101600e-02  5.438596e-02 -3.840951e-02
#> [1061] -1.818657e-02 -8.950839e-02 -2.562138e-01  6.725325e-02 -7.482224e-02
#> [1066] -8.825331e-03 -2.690834e-01 -1.092199e-01 -1.924797e-01  3.135977e-03
#> [1071]  7.045809e-03  4.275701e-02 -2.001181e-01  1.771091e-02  1.675312e-01
#> [1076]  2.040047e-01  1.244352e-01  3.430377e-02 -4.565725e-02  3.235437e-01
#> [1081] -8.013342e-02  6.379140e-02 -1.384129e-01 -1.256045e-01 -8.009295e-03
#> [1086]  1.357567e-01  1.031338e-01 -1.090418e-01  8.565329e-03  4.551562e-02
#> [1091] -2.276525e-02  2.678260e-02  2.210552e-01 -1.191001e-01 -1.502077e-01
#> [1096]  2.326028e-02 -7.096863e-02  1.819318e-01  4.512676e-02  1.799327e-01
#> [1101] -2.051254e-01  5.049391e-02  1.875725e-01 -8.034669e-02  1.294398e-01
#> [1106] -1.286091e-02 -1.287018e-01  1.814488e-01 -6.222890e-02 -6.748889e-02
#> [1111] -9.654455e-03 -1.618728e-01 -1.043980e-01 -6.397073e-02  1.525921e-01
#> [1116]  1.148629e-01 -1.683407e-01 -5.724524e-02 -3.970670e-03  2.828467e-01
#> [1121] -6.409433e-02 -8.259375e-02  6.363438e-02  8.069248e-02  1.598061e-01
#> [1126] -1.408807e-01  4.794258e-02  9.335010e-02  2.905049e-03  5.924636e-02
#> [1131]  7.740096e-02 -2.511378e-01 -2.663967e-02 -1.376273e-01  7.461662e-02
#> [1136]  2.874337e-02 -2.168158e-01 -2.959419e-02  1.146816e-01 -3.963419e-01
#> [1141] -8.315253e-02  9.442702e-02 -4.202322e-02  1.986216e-01  3.700035e-02
#> [1146] -2.250093e-01  1.436882e-01  2.440919e-01  6.705617e-02  6.602580e-02
#> [1151] -4.456823e-02  4.537000e-02 -4.612963e-02 -1.880698e-01  1.543297e-01
#> [1156] -2.198095e-01  7.450048e-02  2.277597e-01 -1.762575e-01 -5.492913e-02
#> [1161] -2.047935e-01 -3.500393e-02 -1.011101e-01  8.914329e-03 -1.179449e-02
#> [1166]  9.638288e-02 -7.750146e-02 -3.988338e-02  6.200603e-02 -2.763746e-01
#> [1171] -1.942161e-01 -5.392546e-02  4.581479e-02 -1.753486e-01  5.411188e-02
#> [1176] -1.963753e-01 -2.053746e-01 -1.128547e-01 -8.506690e-02 -2.845440e-01
#> [1181]  4.412800e-02 -1.121606e-01  1.100497e-01 -2.973895e-01 -2.547553e-01
#> [1186] -3.372258e-02 -2.069217e-02 -2.577223e-01 -4.379515e-02 -5.174985e-02
#> [1191]  1.214503e-01  1.633032e-01  5.370148e-02 -7.734778e-02  9.650142e-03
#> [1196]  1.091251e-01 -1.136387e-01 -2.115579e-02 -1.625228e-01 -3.433656e-02
#> [1201]  8.649931e-02 -1.490429e-01 -7.466162e-02 -1.066379e-01  2.162509e-01
#> [1206] -2.026026e-01  2.130622e-02  7.328448e-03  2.065322e-02 -4.449920e-02
#> [1211]  5.560581e-02  7.885569e-02  7.023541e-02 -1.137744e-01 -1.671168e-01
#> [1216] -3.078833e-01 -2.921937e-02  2.664779e-01  1.803889e-01 -7.195512e-02
#> [1221]  2.023024e-01  5.311716e-02  2.860542e-01  7.719623e-03  1.140253e-01
#> [1226]  5.646414e-02  2.816459e-01 -2.651307e-02  6.857923e-02  8.447909e-02
#> [1231]  6.571572e-02  1.762964e-01  2.706473e-02  1.253913e-01 -9.111331e-02
#> [1236] -1.226245e-01  1.045392e-01 -1.296332e-01 -6.169014e-02  4.934895e-02
#> [1241]  2.283770e-02  1.413908e-01  1.015948e-01 -2.142336e-01 -1.014560e-01
#> [1246]  8.986893e-02 -1.639402e-01 -1.272418e-01  3.930525e-02  9.714332e-03
#> [1251]  6.802599e-02 -6.936920e-02 -1.124275e-01  9.381952e-02 -1.035496e-01
#> [1256]  9.438984e-03 -1.120234e-01  1.171619e-01 -3.420020e-04  3.841904e-02
#> [1261]  4.188338e-02 -1.287804e-01  1.741977e-01 -2.354590e-01 -1.301424e-04
#> [1266] -2.283285e-01 -1.991074e-02  2.023450e-01 -1.338923e-01 -1.580095e-01
#> [1271]  1.061691e-01 -9.399207e-02 -1.510015e-01 -7.523205e-02  7.843492e-02
#> [1276] -1.483032e-01  2.077163e-01  2.049715e-01  1.647568e-01  4.031807e-02
#> [1281]  2.552794e-01  5.414799e-02 -3.383065e-02  2.316450e-01 -7.844863e-02
#> [1286]  1.651576e-01  1.021492e-01 -4.761996e-02  2.216376e-01 -1.627730e-01
#> [1291] -1.593873e-01  1.028064e-01  3.128960e-01  1.448713e-01 -8.454538e-02
#> [1296] -8.739610e-02  1.116779e-01  1.009232e-01 -3.504362e-02  2.025752e-02
#> [1301]  7.274503e-02  6.121917e-02  8.338938e-02 -9.692433e-02  1.355879e-01
#> [1306] -2.718130e-02 -2.464198e-01 -3.918573e-02 -1.149079e-01 -5.332660e-02
#> [1311]  2.084381e-01 -6.249781e-02 -7.292215e-02  1.385170e-01  6.775782e-02
#> [1316]  5.226438e-02  1.459928e-01  3.453080e-01 -2.646696e-01  9.562101e-02
#> [1321] -3.706912e-02  3.143693e-01  7.851280e-02 -2.111272e-01  9.591710e-04
#> [1326]  2.450587e-02 -2.681254e-01  4.568824e-02 -1.366621e-01  1.738871e-01
#> [1331]  2.276630e-02 -1.349284e-01 -2.484100e-01 -2.122538e-01  1.384054e-03
#> [1336] -4.534793e-02  2.031458e-01  6.119636e-02 -6.835507e-02  2.060571e-02
#> [1341]  8.096785e-03  1.701002e-02 -5.111150e-02 -3.334427e-02  1.996399e-01
#> [1346] -1.929382e-01 -1.291863e-01 -3.188167e-02  8.258134e-02 -2.024111e-02
#> [1351] -7.237592e-02 -1.218776e-01  1.736782e-01 -3.027456e-02  3.194920e-02
#> [1356] -1.737362e-01 -3.855314e-02 -2.205055e-02 -6.499131e-02  4.998088e-02
#> [1361]  1.687712e-01 -1.603181e-01 -1.347740e-01 -3.335628e-01  1.061254e-01
#> [1366] -2.306459e-01  7.841413e-02  1.406591e-01  2.470636e-01  4.130559e-02
#> [1371]  6.114486e-02  2.584650e-02  1.183528e-01  3.222571e-01 -4.092724e-02
#> [1376]  2.262114e-01 -1.426698e-01  9.119595e-02 -3.594058e-02 -1.051169e-01
#> [1381]  2.082953e-01 -1.602889e-01 -1.491449e-01  5.374024e-02 -2.030561e-01
#> [1386] -2.520563e-01 -3.412998e-02 -9.297846e-02 -6.034552e-02 -1.375786e-02
#> [1391] -2.393775e-02 -1.618360e-02 -5.420677e-02  1.877674e-01  1.153843e-02
#> [1396]  1.338586e-02 -9.104332e-02  1.410429e-01 -4.197970e-02 -7.544979e-02
#> [1401] -8.027322e-02  1.073278e-01  3.084616e-01 -5.327766e-02 -7.418738e-02
#> [1406]  6.800568e-02  2.107524e-01 -2.021582e-01 -1.269889e-02  1.813112e-02
#> [1411] -1.645727e-01 -1.185705e-01 -1.965556e-01 -2.491203e-01  3.816081e-02
#> [1416] -6.689258e-02 -1.348862e-01  1.072082e-01 -2.164241e-03 -7.889880e-02
#> [1421]  7.157347e-02 -1.151269e-02  1.609684e-01 -2.159656e-03  1.458803e-01
#> [1426] -7.840682e-02  5.336026e-02  1.139710e-01 -1.813562e-01  1.678673e-02
#> [1431] -1.155704e-02 -1.641129e-01  6.245663e-02 -3.326409e-01 -1.309258e-01
#> [1436]  2.148006e-03 -1.162043e-01  6.016432e-02  6.595851e-02 -1.298209e-01
#> [1441]  3.137480e-02  7.034495e-02 -1.404868e-01 -1.684147e-01 -1.791251e-01
#> [1446]  1.594065e-01  1.298105e-01 -1.418450e-01 -1.448161e-01  5.744214e-02
#> [1451] -1.048855e-01 -4.957530e-02  2.343776e-02 -1.299200e-01 -1.672775e-02
#> [1456] -2.054440e-01 -3.580840e-02 -7.949103e-02 -1.778082e-01 -2.377800e-01
#> [1461] -2.820103e-01  9.908909e-03  9.945203e-02 -1.124420e-01 -4.056334e-02
#> [1466]  1.274373e-01 -9.011834e-03  5.937916e-02  7.416046e-03 -2.061766e-01
#> [1471] -6.572654e-02 -2.234569e-03  4.332015e-02 -9.625206e-02 -3.511295e-02
#> [1476]  8.086626e-02  1.614671e-01 -2.428947e-02 -2.011807e-01  9.007850e-02
#> [1481] -7.130043e-02  1.953690e-01 -2.361144e-01  1.289271e-01 -3.075341e-01
#> [1486]  2.248193e-01 -7.543258e-02 -2.893237e-01 -2.057448e-01 -2.330952e-01
#> [1491] -1.099022e-01 -7.460270e-02 -1.489182e-01 -4.599824e-02  5.768074e-02
#> [1496]  7.751103e-02  8.253518e-02  1.935688e-01 -1.050706e-01 -4.173486e-03
#> [1501]  1.392561e-01 -4.910015e-02 -2.034719e-04  1.796621e-01  3.331734e-02
#> [1506] -6.493637e-02  2.065286e-01 -1.419979e-01  1.024567e-01  1.543999e-01
#> [1511]  8.313255e-02 -7.526461e-02  5.422753e-02  3.910851e-01 -2.467230e-01
#> [1516] -6.752028e-02 -9.951624e-02  5.666279e-02 -2.753821e-01  8.898824e-02
#> [1521] -7.575753e-02 -6.489522e-02 -1.166529e-01 -1.207080e-01 -7.440559e-02
#> [1526] -1.272446e-02  1.464577e-01  3.085830e-02 -1.163733e-01  1.317700e-01
#> [1531] -7.999446e-02  1.025246e-01  1.988855e-02  3.024699e-02  6.893998e-02
#> [1536] -9.405023e-02  2.062561e-02  2.397813e-01  2.066881e-01  6.075942e-02
#> [1541] -2.738553e-01  8.496283e-02 -2.433083e-01 -9.856435e-02 -1.398561e-01
#> [1546] -1.627380e-01  1.002548e-01  2.397019e-01 -7.524534e-02 -1.157091e-01
#> [1551] -8.423463e-02 -1.111423e-01 -5.676456e-02  1.277264e-01  2.494119e-01
#> [1556]  1.703042e-01  9.731913e-02  1.121831e-02 -8.980491e-02  2.445504e-01
#> [1561] -1.556106e-01  1.195568e-01 -3.966180e-02 -6.860197e-02 -1.618405e-01
#> [1566] -1.758359e-01 -2.378074e-02 -3.556131e-02 -1.371832e-01  1.021449e-01
#> [1571]  2.195153e-01  1.587421e-02 -1.421734e-01 -6.385931e-02 -8.908086e-02
#> [1576]  8.126206e-02 -9.679100e-02  8.860567e-02 -1.100265e-01 -9.274739e-03
#> [1581]  2.066980e-01  5.137970e-02 -2.004766e-01 -1.007901e-02 -1.154855e-01
#> [1586] -4.330442e-02  3.891508e-02 -2.604749e-01  1.453146e-01 -2.672056e-02
#> [1591]  9.494668e-02 -2.198650e-01 -2.893276e-01 -9.050371e-02  3.373773e-01
#> [1596] -7.296086e-02 -6.612373e-02  1.459193e-01  6.416487e-02 -7.921780e-02
#> [1601]  2.923607e-03  2.846331e-02  2.378186e-01 -5.707110e-02 -2.841331e-01
#> [1606]  3.505605e-03  6.502765e-03 -2.772325e-02 -7.148032e-02 -8.749280e-02
#> [1611]  1.110454e-02 -1.730871e-01  2.309804e-01 -8.298528e-02  1.153003e-01
#> [1616]  1.799780e-01  7.474395e-02  9.599188e-02  4.680431e-03  2.381632e-02
#> [1621] -4.770439e-02 -2.951792e-01 -5.714144e-02 -2.804025e-01 -5.787518e-02
#> [1626] -8.672662e-03 -2.426784e-01 -3.089287e-03 -1.283777e-02  1.164117e-02
#> [1631]  1.530855e-02  1.861821e-01 -7.100140e-02 -2.738506e-02 -9.041431e-02
#> [1636]  7.658630e-02  2.713338e-02 -1.792494e-02  7.342310e-02 -4.494330e-02
#> [1641]  1.857323e-01 -4.849441e-02  1.944659e-01  5.468896e-02  1.615680e-01
#> [1646]  2.512166e-01  2.049538e-01  1.490741e-01  4.164931e-02 -1.281396e-01
#> [1651]  1.577662e-01 -1.417874e-02 -7.708436e-02 -4.337079e-02 -4.935329e-02
#> [1656]  3.853745e-02 -6.870856e-02  9.902662e-02  2.505166e-02  1.816971e-01
#> [1661] -1.097208e-01  2.300107e-01  3.410242e-01  1.912758e-01 -2.523512e-01
#> [1666] -2.899353e-02  1.293875e-01 -5.745252e-02 -4.804785e-02  9.298489e-03
#> [1671]  3.819611e-02  5.761899e-02 -2.905598e-01  1.106792e-02 -1.146099e-01
#> [1676] -4.855068e-02  2.554198e-01 -1.915391e-01 -4.740070e-02  1.562711e-02
#> [1681]  1.260723e-01 -1.369277e-03 -1.225404e-01 -2.036901e-01  5.297984e-02
#> [1686]  1.142177e-01 -1.353252e-01  1.331458e-01 -7.446952e-03  1.888416e-01
#> [1691]  1.067136e-03  4.969231e-02 -1.668039e-01 -2.126814e-01 -5.016958e-02
#> [1696]  3.777590e-02 -2.534958e-01 -9.157192e-02 -1.292349e-01 -3.039008e-02
#> [1701]  1.507279e-02  1.319215e-01 -7.706729e-02  3.037121e-02  9.849061e-02
#> [1706]  1.676941e-01  1.057695e-02  1.799762e-01  1.328829e-01 -1.823365e-01
#> [1711]  2.197500e-02 -3.414501e-01  1.002592e-01  2.192627e-01  6.508591e-02
#> [1716]  1.924466e-03 -1.449961e-01 -1.717744e-01 -1.633121e-01 -2.439327e-01
#> [1721] -8.529433e-02  2.724366e-01 -1.401124e-01  1.627076e-01  1.436948e-02
#> [1726] -9.515772e-03 -1.256558e-01 -2.093979e-01 -7.066758e-03 -4.837489e-02
#> [1731]  2.270441e-01 -2.274103e-01 -5.699943e-02  2.358672e-01  1.752398e-01
#> [1736]  5.856655e-02 -2.815290e-02  2.798529e-01  5.396782e-02  2.676557e-02
#> [1741] -1.621251e-01  2.037526e-01  2.009031e-02 -1.677159e-01 -1.065608e-01
#> [1746] -1.092838e-01  3.004490e-02 -2.125262e-02 -8.003028e-02 -9.183876e-02
#> [1751] -9.283819e-02  8.697708e-02  7.544684e-02  1.814911e-01  1.260911e-01
#> [1756]  4.322434e-01 -8.068302e-02  7.077506e-02 -2.468109e-01  1.373700e-01
#> [1761]  2.865999e-01 -5.863611e-02  1.330850e-01  2.086491e-01  1.113358e-01
#> [1766] -1.033125e-01  5.092737e-02  3.757499e-02  9.281281e-02 -5.049215e-02
#> [1771]  6.328445e-02 -1.904723e-01  7.688968e-03  1.844635e-01 -6.780462e-03
#> [1776] -2.243222e-01  3.427940e-01 -2.880001e-02  1.273110e-01  1.940727e-01
#> [1781]  4.705057e-02 -4.567068e-02  1.340433e-01  3.026933e-02  4.746664e-02
#> [1786] -1.467929e-01 -5.941058e-02 -1.881058e-01 -9.778812e-02 -1.910703e-01
#> [1791] -5.186933e-02  7.731349e-02  1.012998e-02  6.118157e-02 -4.716123e-02
#> [1796]  9.177159e-02 -7.226189e-02  2.119372e-01  7.556311e-02 -2.332616e-01
#> [1801]  4.903690e-02 -6.204192e-02  2.640962e-01 -1.283784e-01  1.734026e-01
#> [1806]  4.085913e-02  2.554137e-01  3.033885e-01  2.151734e-01 -2.693173e-01
#> [1811]  3.329761e-03  1.183878e-01  1.225497e-01 -4.266707e-02 -4.275744e-04
#> [1816] -8.741528e-02 -8.700838e-02  8.397266e-02 -2.061455e-01  4.271610e-02
#> [1821] -1.374168e-01  1.403273e-02  3.153404e-02 -1.404291e-01  1.207929e-01
#> [1826]  3.732390e-02 -1.129408e-01 -1.708359e-01  1.505345e-01 -7.729675e-02
#> [1831] -1.215564e-01 -1.204961e-02  1.848798e-01 -6.845545e-02  1.018404e-02
#> [1836] -3.065521e-02 -5.615706e-02 -2.944773e-02 -2.052612e-02 -2.632480e-01
#> [1841] -1.830967e-01 -1.657417e-03  7.631289e-02  8.171797e-02  1.057702e-02
#> [1846] -7.633067e-02 -1.559213e-01  2.529392e-01  2.160831e-02 -4.321571e-02
#> [1851]  1.925068e-01 -6.708499e-02  1.497045e-01  3.230519e-02 -5.416660e-02
#> [1856]  1.369836e-01  2.227171e-01 -6.207006e-02  3.214644e-03  4.048071e-02
#> [1861] -4.656722e-02 -2.157039e-02 -7.470864e-02  7.299978e-02  9.617104e-02
#> [1866] -1.937792e-02 -1.388654e-01 -2.133086e-01  1.695412e-01  4.624374e-02
#> [1871]  1.224494e-01  7.025505e-03 -3.815499e-02  2.846918e-01  1.254838e-01
#> [1876]  2.701486e-02  3.089082e-02  4.136989e-02  3.700208e-02  3.238266e-02
#> [1881] -3.638839e-02  2.708919e-02  2.540226e-01 -1.433569e-01 -7.251532e-02
#> [1886] -1.140832e-01 -1.249230e-01  4.105484e-02 -6.237080e-02  1.120827e-01
#> [1891] -1.016675e-01  1.653470e-01 -7.339819e-02  1.007768e-01  3.283814e-03
#> [1896]  1.445337e-01  1.348317e-01 -3.394664e-02 -8.359001e-02  1.875772e-02
#> [1901] -2.986310e-01 -3.963582e-02  5.899161e-02 -1.431086e-01  1.007871e-01
#> [1906] -2.652194e-02 -1.784291e-01  5.298682e-03  1.012690e-01  4.529131e-02
#> [1911]  1.853021e-01  5.373418e-02 -5.462509e-02 -5.450647e-02 -4.027568e-03
#> [1916] -4.257745e-02 -4.763618e-02  1.615685e-02 -1.736326e-01  1.108660e-01
#> [1921]  1.251603e-01 -3.656424e-02 -6.949321e-05 -3.833669e-02  1.011518e-01
#> [1926] -1.214975e-03 -1.227169e-01  1.060463e-01 -6.756880e-02 -2.086050e-02
#> [1931]  6.822366e-02 -2.109891e-01 -4.824023e-02 -5.099493e-02  1.388583e-01
#> [1936]  1.787303e-02 -7.541344e-02  1.771406e-01  2.361356e-01 -8.609873e-02
#> [1941]  1.503338e-01 -3.979651e-02  8.071185e-02  4.144416e-02  1.633049e-02
#> [1946] -3.248021e-01  2.690902e-02 -1.945067e-01 -8.312581e-02 -1.326373e-02
#> [1951] -6.326638e-02 -2.070213e-01  1.909465e-01  1.185293e-01 -2.184241e-01
#> [1956] -1.885242e-01 -2.459813e-01  2.455487e-01  2.424668e-01 -2.791576e-02
#> [1961]  1.954869e-02 -1.823289e-01  4.022680e-01 -3.518689e-02  1.160776e-02
#> [1966]  5.142417e-02 -2.445196e-02  5.349789e-02  1.325415e-01 -1.015635e-03
#> [1971]  1.630434e-01 -1.824096e-02 -2.595808e-01  4.245941e-02  1.077178e-01
#> [1976] -6.962311e-02 -7.394421e-02 -1.836906e-01  1.217884e-01 -1.664442e-01
#> [1981]  5.290918e-02 -4.036105e-02  2.086980e-02  5.598299e-02 -1.342569e-01
#> [1986] -2.626287e-02 -6.063714e-02  2.420470e-02  2.082024e-02 -1.372651e-01
#> [1991] -1.026857e-01 -1.083353e-01 -2.416434e-01  7.607823e-02  2.574723e-01
#> [1996]  1.612706e-01  1.541876e-01 -1.171747e-01 -3.196793e-01 -1.751916e-02
#> [2001]  1.901910e-01 -1.059861e-01  8.659642e-02  8.027624e-02 -2.708177e-02
#> [2006]  6.630685e-02 -2.534639e-01  8.468475e-02  2.394300e-01 -1.392491e-01
#> [2011]  9.430855e-03  1.411697e-01 -1.369255e-01  1.220864e-01 -3.054673e-03
#> [2016]  1.429185e-01  1.427615e-01 -3.562119e-01  7.146276e-02  1.343967e-01
#> [2021] -7.210982e-03 -1.835156e-02 -1.401179e-01 -1.390443e-01  7.493503e-02
#> [2026]  1.259515e-02 -1.528016e-01  1.286627e-01 -9.434260e-02 -1.312847e-01
#> [2031] -1.899133e-01 -5.276992e-02 -2.813581e-02 -1.529888e-01  1.614782e-01
#> [2036] -1.287841e-01 -1.693476e-01 -1.825484e-02 -4.481668e-02 -1.035235e-01
#> [2041] -4.960328e-02  2.260606e-01  2.920462e-02  5.238065e-03  1.285445e-01
#> [2046] -2.744673e-01 -7.616867e-02 -5.109591e-02 -9.494840e-02 -2.513275e-02
#> [2051] -6.137919e-02 -1.466623e-02 -7.898665e-02 -2.995978e-02  2.003837e-01
#> [2056]  2.332791e-02 -4.356657e-02  4.857505e-02  1.258474e-01  2.175494e-01
#> [2061]  4.895568e-02  2.788006e-02  6.425016e-02  8.858452e-02 -1.384604e-01
#> [2066] -3.106929e-02 -2.646367e-01  8.571440e-02  1.759400e-01 -4.748503e-02
#> [2071] -9.235395e-02  8.306842e-02 -3.579573e-01 -1.255504e-01 -1.838677e-01
#> [2076]  3.938004e-02  3.184013e-01  1.648159e-01 -1.853790e-02  6.897633e-02
#> [2081]  3.544423e-02 -1.672340e-02  2.227118e-01 -1.379724e-01 -1.536031e-01
#> [2086] -9.869627e-02 -1.881533e-01  8.561980e-02 -1.264656e-01 -2.340143e-01
#> [2091] -2.118906e-01  5.165377e-02  1.518293e-01 -2.308856e-02 -2.223093e-01
#> [2096] -2.923610e-01  9.148189e-02 -1.262803e-01 -1.100623e-01  1.858582e-01
#> [2101] -7.253812e-03  9.001332e-02 -1.147833e-01 -1.567475e-02 -6.759039e-02
#> [2106]  1.550161e-02  2.694060e-01 -1.691322e-01 -1.121261e-01  1.867829e-01
#> [2111]  1.507510e-01 -2.487878e-01 -3.807478e-02 -1.926265e-01 -1.701537e-01
#> [2116] -4.510885e-02 -1.002131e-01 -5.266050e-03 -2.615660e-02 -1.959869e-02
#> [2121]  3.796546e-02 -9.506963e-02  2.142754e-01 -8.594021e-02  1.312652e-01
#> [2126]  3.569105e-02 -1.997264e-01  4.875658e-02  2.036509e-01 -1.331711e-02
#> [2131] -8.732848e-02  1.311453e-01  6.837509e-02 -5.697356e-02 -1.106343e-01
#> [2136] -4.541306e-02 -8.194496e-02  2.430583e-02 -1.628982e-01  1.803623e-01
#> [2141] -4.449980e-02  1.378550e-01  8.683881e-02 -7.889030e-02  8.076184e-02
#> [2146]  2.986009e-02 -8.874845e-02  4.198010e-01  2.117286e-01  1.801705e-02
#> [2151] -6.170386e-02  6.980321e-02 -3.457446e-01  6.180074e-02 -1.874169e-01
#> [2156] -1.470876e-01 -2.259925e-02 -8.789447e-02  1.771567e-01  3.801977e-02
#> [2161] -4.430944e-02 -5.991697e-02 -1.391399e-01 -5.016345e-02 -1.219559e-01
#> [2166] -1.549798e-01 -3.148847e-01 -1.407450e-02 -7.570725e-02  2.154135e-01
#> [2171] -1.189585e-01  1.081077e-02  3.500094e-02  9.278058e-02 -1.610739e-01
#> [2176]  6.047341e-02 -7.094472e-02 -1.845638e-01  2.659469e-01 -6.985317e-02
#> [2181] -6.756504e-02  1.945299e-01  3.448995e-02  5.901097e-02  6.129971e-02
#> [2186] -1.473638e-01 -6.920057e-02  1.705013e-01  5.267214e-02 -8.294018e-02
#> [2191]  4.601337e-02  4.720421e-02 -1.686601e-01 -2.629667e-02  4.423920e-02
#> [2196] -1.084209e-01  1.759578e-01  1.231560e-01  2.315614e-01 -1.722161e-01
#> [2201] -1.425029e-01 -2.071054e-01 -1.384640e-02 -1.678687e-02  2.063361e-01
#> [2206]  1.288059e-01  1.038403e-01  1.445092e-01 -1.778804e-01  2.477582e-01
#> [2211]  1.503586e-01 -5.590950e-02  2.392772e-01 -1.673233e-01 -1.312672e-02
#> [2216] -6.798466e-03 -7.313022e-02  5.476483e-02  8.435235e-02  2.309656e-01
#> [2221] -4.504553e-02 -1.895330e-01 -7.132161e-02 -2.442916e-02 -4.473257e-02
#> [2226] -5.439496e-02 -1.174684e-01  1.506410e-01  1.307354e-01 -3.425196e-02
#> [2231]  1.091551e-01 -8.368343e-02 -6.625729e-02  1.099916e-01  1.107206e-02
#> [2236]  6.032407e-02  1.343218e-01 -5.674498e-02  1.064551e-01 -8.820431e-02
#> [2241] -8.498504e-02 -1.463681e-01 -1.970471e-01 -6.614343e-02 -1.662882e-01
#> [2246] -1.973468e-01  2.092547e-01 -2.703529e-01 -3.232152e-02 -1.072565e-01
#> [2251] -3.304054e-02 -1.082353e-01  1.568041e-01 -4.671598e-02 -1.809120e-01
#> [2256]  4.364945e-02  2.512715e-01 -1.814621e-01 -2.964713e-01 -4.147894e-02
#> [2261] -3.315950e-02  4.111135e-01  8.566721e-02 -3.794996e-02 -1.360799e-01
#> [2266]  2.559602e-01  5.167202e-02  1.258124e-01  1.014136e-01  7.050476e-02
#> [2271]  6.355628e-02 -9.252931e-02 -7.465717e-02 -2.693537e-01 -6.217786e-02
#> [2276] -2.175765e-01  1.699992e-02 -2.027623e-01 -1.606622e-01 -1.059956e-01
#> [2281]  1.179490e-01 -2.421674e-02  7.327521e-03  3.006413e-02  2.425665e-02
#> [2286] -7.943240e-02 -1.590679e-01 -3.381156e-02 -1.116236e-01 -4.114468e-02
#> [2291]  1.380475e-01  7.652799e-02  4.903229e-02 -7.103103e-02  1.465772e-01
#> [2296]  2.771109e-02 -2.806201e-02  5.188374e-02 -6.658079e-02 -3.126843e-02
#> [2301]  2.516685e-01  9.321676e-02  3.520431e-02 -3.903826e-02 -2.166898e-01
#> [2306] -1.688171e-01 -1.122851e-01 -4.679026e-02 -9.414828e-02 -7.775014e-02
#> [2311] -3.371586e-01 -5.694583e-02 -3.487546e-02 -6.004691e-02 -1.145799e-01
#> [2316] -1.743652e-01  4.846878e-02  2.959094e-01  2.536305e-01  7.163041e-02
#> [2321]  2.358500e-01  5.196102e-02  1.444250e-01 -6.858566e-02  1.553186e-01
#> [2326] -4.664942e-02 -5.064882e-02  7.655997e-02 -1.325323e-01  1.977416e-01
#> [2331] -5.782031e-02  2.495429e-01 -3.128015e-02 -9.989253e-02 -1.402149e-01
#> [2336]  1.884967e-01 -1.616591e-01  1.150019e-01  7.585651e-02  1.894272e-01
#> [2341] -8.623530e-02  9.673561e-02  4.382731e-02  1.600445e-01  2.738959e-01
#> [2346]  2.795234e-02 -1.115067e-02  6.166911e-02 -6.926655e-03 -4.628980e-02
#> [2351] -2.619860e-01 -2.080935e-01  3.427239e-01 -4.260513e-02 -3.441779e-02
#> [2356] -7.490393e-02 -5.854800e-03 -2.188031e-01  2.631182e-01  1.435342e-01
#> [2361] -2.789388e-01 -2.031836e-01 -1.214934e-02 -2.852261e-01 -3.709953e-03
#> [2366]  8.146185e-02  1.102292e-01 -5.973646e-04 -8.825266e-02 -1.254156e-01
#> [2371] -1.622021e-01 -1.981841e-01 -9.858130e-02 -6.289803e-02  1.055714e-01
#> [2376]  1.566495e-01 -1.164663e-02  3.321187e-02  1.126942e-01 -3.597655e-02
#> [2381]  1.527001e-01  1.483067e-01  1.707830e-01  7.248128e-02 -2.185610e-02
#> [2386]  4.536259e-02  5.590754e-02 -7.789236e-03  7.150456e-02 -5.498033e-02
#> [2391]  6.790935e-02 -1.721262e-02  6.005222e-02 -2.627236e-01  3.409572e-01
#> [2396]  1.686189e-01 -1.310845e-01 -7.404526e-02 -3.840890e-02  8.525402e-02
#> [2401]  5.746289e-02 -4.128129e-03  1.136552e-01  1.584509e-02 -2.259186e-02
#> [2406] -7.931177e-02  2.373477e-01 -2.385525e-02 -1.030532e-01 -2.689386e-01
#> [2411] -1.642690e-02  1.138200e-01  8.568203e-02 -1.045586e-01  1.154285e-01
#> [2416] -1.264438e-01 -3.355567e-02  4.256040e-02  1.372356e-01 -1.275191e-01
#> [2421] -2.482462e-01 -1.790953e-02 -7.494219e-02  2.402007e-01 -2.985482e-02
#> [2426] -6.169542e-02 -3.286134e-01 -2.226245e-01  9.549404e-02 -7.058981e-02
#> [2431] -9.807233e-03  1.935581e-01 -1.167997e-01  5.897009e-02  1.948030e-01
#> [2436] -1.220345e-01 -4.568788e-02 -7.938172e-02 -4.509297e-02 -1.364025e-01
#> [2441]  2.006811e-01  1.659361e-01 -7.041024e-02 -3.058137e-02  1.057488e-01
#> [2446] -7.176163e-02 -9.860589e-02  4.442007e-02 -1.194986e-01 -9.811902e-02
#> [2451]  1.004292e-01  2.236209e-01 -7.825438e-02 -1.715273e-02 -6.249896e-02
#> [2456] -2.029853e-02  2.707232e-01 -9.638940e-03  2.667973e-03 -3.415709e-02
#> [2461] -2.694695e-01  2.017688e-01 -3.078128e-01  1.533774e-01 -6.235520e-02
#> [2466] -3.652426e-02 -7.218904e-03 -3.682326e-03 -1.512592e-01  1.223380e-01
#> [2471] -6.589642e-02  2.684182e-01  4.283479e-03  1.913612e-01  1.334894e-01
#> [2476] -9.556052e-02 -4.945409e-02 -2.407343e-02 -1.388002e-01  1.130625e-01
#> [2481] -4.210070e-02  1.136726e-01  1.473731e-01 -2.090066e-01  2.265041e-02
#> [2486] -1.649546e-03 -8.064921e-02 -1.847279e-01 -2.063511e-01 -8.963260e-02
#> [2491]  2.346564e-02  1.573523e-01 -1.786607e-01  8.797821e-02 -2.624296e-01
#> [2496]  5.826650e-02  1.119351e-02  9.187793e-02 -7.372566e-03  2.572041e-01
#> [2501]  2.788172e-02  5.829989e-04  8.324427e-03  2.217398e-01 -5.686873e-02
#> [2506] -3.170514e-02 -1.008389e-01 -1.553556e-01  1.158269e-02 -2.272245e-02
#> [2511] -7.353477e-02  1.463035e-01 -2.326822e-03  1.982842e-01 -5.854526e-02
#> [2516] -6.066072e-03  2.016361e-02 -1.075163e-01  9.100425e-02  2.877001e-01
#> [2521]  8.400852e-02 -2.397337e-02 -1.450241e-01  1.310844e-01 -1.041886e-01
#> [2526]  4.676477e-02  1.725752e-01  4.383319e-02 -1.363728e-02 -1.418887e-01
#> [2531]  1.472750e-01  3.476784e-02  6.080637e-02 -7.217018e-02  2.944601e-01
#> [2536] -9.729208e-02  8.250815e-02  9.981101e-02 -2.671303e-02  1.155084e-01
#> [2541]  1.148859e-01  3.176274e-02 -2.385245e-02 -6.369696e-02 -6.122803e-02
#> [2546] -1.718389e-02  1.873071e-01  2.236876e-01 -5.363799e-02  3.387555e-02
#> [2551] -8.991457e-02  8.186906e-02 -1.175548e-01  1.395091e-01  7.885774e-02
#> [2556] -7.484190e-02 -5.955197e-02 -1.935476e-01  1.040540e-02 -1.281838e-03
#> [2561] -1.013767e-01 -2.039793e-01 -1.584992e-01 -8.989960e-02  1.571145e-02
#> [2566] -1.097368e-01 -9.595232e-02 -2.221260e-02 -2.408765e-02 -2.375430e-02
#> [2571] -1.069024e-01 -1.154835e-01  3.566649e-02 -9.775986e-02 -1.221471e-01
#> [2576]  1.045670e-01 -5.895955e-02  1.472484e-01  3.052508e-02 -5.459544e-02
#> [2581]  5.632738e-02  1.015101e-01 -1.735071e-01 -3.505232e-02  7.904115e-02
#> [2586] -9.912833e-02  1.070952e-03  1.534971e-01 -1.466618e-01 -1.340286e-02
#> [2591]  3.005661e-02 -1.726916e-01  8.793458e-02  1.462831e-01 -2.760223e-01
#> [2596]  1.104469e-01  1.228358e-02 -8.155747e-02  1.681214e-01  1.063133e-01
#> [2601] -1.321237e-01 -4.683507e-02  1.530435e-02 -1.631256e-01 -4.924880e-02
#> [2606]  1.422409e-01  2.974136e-02 -1.191059e-01  8.127733e-03  1.317628e-01
#> [2611]  3.805772e-02  5.473107e-02  8.171686e-02  1.246205e-01 -1.626644e-01
#> [2616]  4.092388e-02 -1.544218e-01 -1.125391e-01  5.720609e-02 -2.325833e-01
#> [2621]  2.801852e-02  1.240108e-01  1.024792e-01 -2.078015e-02  4.812893e-02
#> [2626]  8.243964e-02  1.922256e-01 -7.892802e-02 -2.491772e-02  1.531486e-01
#> [2631] -2.267624e-01  1.071466e-01  3.429206e-02  3.292355e-02 -1.021060e-01
#> [2636]  1.290650e-01  2.380341e-02  2.888295e-02  9.079601e-02 -2.936836e-02
#> [2641]  1.535809e-01  7.238943e-02 -3.833629e-02 -8.782689e-02  1.339486e-01
#> [2646] -7.763635e-02  5.971752e-02 -7.912631e-02 -2.535658e-02 -1.193392e-02
#> [2651]  1.050693e-01  1.106694e-01  3.892108e-02  1.183178e-01  7.316060e-02
#> [2656] -1.627554e-02  6.366558e-02  8.058105e-02 -4.415451e-02 -5.528526e-02
#> [2661]  1.054893e-01 -2.662148e-01  3.540209e-02 -1.036921e-01 -4.858904e-02
#> [2666] -8.156291e-02 -4.079189e-02  3.457716e-02  1.201004e-01  2.303832e-01
#> [2671]  2.721680e-01 -1.182564e-01  1.119082e-01  6.480934e-02 -2.721144e-02
#> [2676]  1.602676e-01 -2.532320e-01 -1.794738e-02  9.354466e-02  8.449133e-02
#> [2681]  1.444659e-01  1.051655e-01 -5.540198e-02  2.713991e-01 -6.704425e-02
#> [2686] -4.906851e-02 -3.022331e-02  6.597578e-02  1.877909e-01 -2.825712e-02
#> [2691] -1.538207e-01  2.350821e-01  8.336106e-02 -4.127997e-02 -1.942914e-01
#> [2696]  1.919467e-01  1.875056e-01 -6.584597e-02 -7.745778e-02 -5.866013e-02
#> [2701] -2.426759e-01  7.417755e-02  6.691252e-02  1.534338e-01  3.990210e-02
#> [2706] -2.966878e-02  9.024233e-03 -1.117147e-02 -3.539946e-02  8.486934e-02
#> [2711]  1.849389e-02 -1.593444e-02  6.577281e-02 -5.586756e-02 -1.368413e-02
#> [2716] -1.888315e-01 -1.372547e-01 -4.678049e-04  4.218053e-02  1.985130e-01
#> [2721] -9.479823e-02  2.830933e-02  6.457068e-02  6.547084e-02 -4.543989e-02
#> [2726] -1.076493e-03  1.890282e-01 -1.712112e-02 -1.601830e-01  2.222189e-01
#> [2731] -1.377285e-01  4.540733e-03 -1.681853e-01  2.538994e-01 -2.290295e-02
#> [2736]  5.139856e-02  2.005413e-02 -1.798409e-02 -8.939992e-02 -8.609888e-02
#> [2741]  9.377000e-02 -4.071200e-02  1.581129e-01  1.639460e-01 -1.215518e-01
#> [2746]  1.014716e-02  1.631148e-01  1.669047e-01  1.028319e-01 -1.710897e-01
#> [2751]  2.218429e-01  7.323837e-02  2.243712e-01  3.380945e-01 -1.480635e-01
#> [2756]  2.068863e-02  9.674257e-02 -2.996948e-01  2.439663e-01 -1.164125e-01
#> [2761]  2.278224e-02  1.690930e-02  1.156983e-01  9.651585e-03  4.711209e-02
#> [2766]  1.900154e-01  7.224565e-02 -3.164504e-02  1.167927e-01 -2.497335e-03
#> [2771] -5.460345e-02 -8.921488e-02 -1.255135e-01  2.029648e-01  1.069326e-01
#> [2776] -4.365455e-03  1.550182e-01  1.016786e-01  9.144652e-02 -2.365451e-01
#> [2781]  7.963347e-03  6.044835e-02 -2.290655e-01  3.050175e-02 -1.141580e-01
#> [2786]  8.880446e-02  2.166641e-01  1.236589e-01  1.271388e-01  1.017528e-01
#> [2791]  5.560008e-02  1.481905e-01  8.382086e-02 -6.105507e-02  3.646508e-02
#> [2796] -1.707025e-01 -2.558622e-02 -1.838158e-01 -5.556743e-03 -2.506091e-01
#> [2801] -6.907965e-02  3.413715e-02 -1.521260e-01 -8.494683e-02 -2.737259e-01
#> [2806]  2.837894e-01 -1.058962e-01  6.228541e-03  1.260615e-01  1.018349e-01
#> [2811]  6.888635e-02 -1.375256e-01 -1.423059e-01 -5.049251e-02 -2.140929e-01
#> [2816]  4.410258e-02  4.115345e-02  3.306692e-01  5.449731e-02  1.953739e-01
#> [2821] -7.927451e-02  6.251287e-02 -4.664398e-02 -2.288163e-01 -2.103474e-02
#> [2826]  1.482096e-01  1.723085e-01  1.990185e-01  2.246805e-02  1.934050e-01
#> [2831]  4.867789e-02  2.270342e-01  1.141562e-01  4.316611e-02  2.107102e-01
#> [2836] -2.093118e-02 -1.005978e-01 -8.511764e-02 -4.110388e-02 -5.657880e-02
#> [2841] -2.216667e-02  1.436093e-01  3.386881e-03 -6.312927e-02 -3.100858e-02
#> [2846]  1.352068e-01  7.247798e-02 -1.086204e-01 -2.951081e-02 -2.621436e-02
#> [2851] -1.986596e-01  9.848157e-02 -3.775112e-02 -1.383333e-01 -8.352189e-02
#> [2856] -4.338283e-02 -6.513394e-02 -8.971315e-02 -1.818792e-01 -1.308385e-01
#> [2861]  3.201361e-02 -7.446544e-02  1.225316e-01 -2.384819e-02 -1.214988e-01
#> [2866] -6.596683e-02 -7.217311e-02  1.710444e-01 -1.599447e-01  1.224135e-01
#> [2871] -4.509143e-02  1.624946e-01  3.382067e-02  9.702816e-03 -2.146393e-02
#> [2876] -4.681015e-02 -9.878666e-02  7.018783e-02 -1.656301e-01 -5.186837e-02
#> [2881] -7.683027e-02  2.743558e-01  3.160351e-02 -7.986140e-03  7.487728e-03
#> [2886]  2.829669e-01  1.632464e-01 -1.013299e-01 -2.044085e-01 -7.592746e-02
#> [2891] -2.758695e-02  8.065078e-04 -3.779550e-02  2.322418e-01  3.238248e-01
#> [2896] -4.925075e-02  6.353590e-03 -5.022983e-02 -8.200317e-02  3.288444e-01
#> [2901] -7.868343e-02  2.125021e-02  8.926581e-02  9.420632e-02 -4.243402e-02
#> [2906]  7.063808e-02 -1.435213e-02  2.724139e-01  4.293518e-02  3.230780e-01
#> [2911] -1.408094e-01  8.583000e-02 -7.886652e-02 -5.934044e-02 -1.766303e-01
#> [2916] -2.686909e-02 -1.361925e-01  1.715338e-01 -3.112615e-02 -4.649072e-03
#> [2921] -3.600671e-03 -1.519588e-01  1.705730e-01 -6.306830e-02  1.966938e-02
#> [2926]  1.850782e-01 -1.077295e-01 -5.033439e-02  2.006024e-02 -3.357467e-01
#> [2931]  2.415369e-02 -4.567586e-02 -1.036595e-01 -7.282298e-02 -2.875541e-01
#> [2936]  8.864865e-02 -1.867987e-02 -1.157576e-01 -6.276896e-03 -4.991802e-02
#> [2941] -2.671450e-02 -9.637871e-02 -1.710990e-01  9.982388e-02  1.515004e-03
#> [2946] -1.305850e-01 -7.384201e-02 -1.790814e-02  5.468548e-03  2.243536e-02
#> [2951]  1.223538e-01  1.787795e-03 -7.528259e-02  1.102241e-01 -1.450652e-01
#> [2956] -1.423939e-01 -1.524030e-01 -2.068675e-02  3.569585e-02 -1.886691e-01
#> [2961]  7.157906e-02 -5.179849e-02  1.870109e-01 -2.619709e-02 -1.730166e-01
#> [2966] -8.006598e-02  2.287993e-02  1.233095e-01  1.599031e-01  1.873384e-01
#> [2971]  1.046929e-01 -8.907437e-02  1.071517e-02  6.124491e-02  1.209044e-01
#> [2976]  1.679852e-01 -1.598550e-01 -8.703327e-02  2.272544e-01  2.220661e-01
#> [2981]  1.635872e-01  1.643226e-01 -2.221268e-01  8.805612e-02  1.338226e-01
#> [2986]  1.022391e-01  1.010706e-01 -6.096710e-02 -2.330847e-02  2.081873e-01
#> [2991] -4.579661e-02  1.269171e-01  1.685334e-02 -3.897833e-02  9.990444e-02
#> [2996]  4.042434e-02 -2.203502e-01 -7.999366e-02  5.884704e-02 -5.749183e-03
#> [3001]  1.392418e-01 -1.049744e-01  2.938288e-01 -2.545964e-01 -2.965872e-01
#> [3006]  1.316497e-01  1.366684e-01 -3.539692e-02 -3.512472e-01 -1.512015e-01
#> [3011]  7.536814e-03  6.710855e-02  1.094674e-01  4.123215e-02 -1.469116e-01
#> [3016] -1.264965e-01 -1.912242e-02  5.525370e-02 -9.213818e-02  3.045507e-02
#> [3021] -1.396605e-01  1.207059e-01 -1.600192e-01 -2.011689e-01  8.184709e-02
#> [3026]  2.576369e-01  9.635957e-02  2.839232e-02  1.036277e-02  1.990667e-01
#> [3031]  9.302875e-02 -2.007770e-01 -7.732213e-02 -1.254965e-01 -2.869991e-01
#> [3036] -1.014379e-01 -9.249859e-02 -1.574486e-01 -1.373674e-01 -3.108217e-01
#> [3041] -5.076962e-02  8.343957e-02 -1.101956e-01  2.359725e-02 -6.977838e-02
#> [3046]  1.807632e-01  1.896452e-02 -4.572039e-02 -2.142251e-01  3.017205e-02
#> [3051]  9.780508e-03  4.344004e-03 -2.334662e-02 -5.125261e-02 -6.738951e-03
#> [3056]  2.355179e-01 -1.183137e-01  1.006636e-01 -7.789093e-02  3.184308e-02
#> [3061]  7.878936e-02  1.160233e-01  1.200922e-01  1.075635e-01  8.265960e-02
#> [3066] -1.248564e-01  8.712068e-03  1.738968e-01  1.486451e-01  9.605358e-02
#> [3071]  2.954519e-02  3.594049e-02 -2.233758e-01  1.453150e-01 -6.633741e-02
#> [3076]  2.183023e-01 -5.214359e-02  3.978271e-02 -2.856686e-02 -7.026809e-02
#> [3081]  7.066872e-02 -1.762472e-01  2.918864e-02 -1.678197e-01 -3.840079e-02
#> [3086] -3.015403e-02 -3.956374e-02 -1.068972e-01 -1.414263e-01  1.280100e-01
#> [3091]  5.684708e-02  7.164926e-02  1.224096e-01  2.106221e-01 -7.992421e-02
#> [3096] -3.021936e-02  9.083419e-03 -1.188755e-01 -7.657836e-02 -2.268285e-02
#> [3101] -1.146132e-01 -3.245302e-01 -3.616980e-02 -2.538169e-02  7.314450e-02
#> [3106] -1.574942e-01 -1.016563e-01  1.779094e-01  7.181610e-02 -7.187776e-02
#> [3111] -1.142321e-01  1.029355e-01 -5.061792e-02 -9.379796e-02  8.768646e-02
#> [3116]  1.232771e-01 -1.960356e-02 -1.301229e-01 -3.228144e-01  6.800181e-02
#> [3121]  4.271646e-02  1.512929e-01 -1.265429e-02  1.482215e-01 -1.321161e-01
#> [3126]  3.713652e-02  2.209642e-01  6.261320e-02  5.320887e-02 -1.173626e-01
#> [3131]  8.380883e-03 -4.438971e-02  1.290202e-01  4.419489e-01  3.895511e-02
#> [3136]  2.044532e-01  1.688143e-01  7.106548e-02 -1.253708e-01  2.175270e-01
#> [3141]  2.501768e-01 -9.317359e-02 -6.806926e-02 -2.375546e-01 -6.207106e-02
#> [3146] -1.040245e-01  7.671080e-02  8.643514e-02  1.354336e-01  6.358263e-02
#> [3151] -4.776848e-02 -4.584013e-02  2.027355e-01 -1.638673e-01  1.128428e-01
#> [3156] -2.571196e-02  3.856347e-01  5.751756e-02  3.441694e-02  4.804785e-02
#> 
#> $mod2_unimodal$sim_properties$cell$pca_dist
#>    [1] 28.49896 30.33067 32.82978 35.43662 29.42893 30.18457 34.42269 20.99348
#>    [9] 27.33423 35.24269 26.10286 30.19539 30.06146 35.50999 31.57713 28.90776
#>   [17] 29.05613 28.33414 32.02926 25.88401 28.72518 36.60415 31.29330 28.47522
#>   [25] 21.55328 34.59361 32.43416 20.25224 33.25403 35.11347 36.14950 24.80630
#>   [33] 26.64831 31.10012 29.32461 28.23026 29.14772 29.87525 29.63456 35.30291
#>   [41] 35.18635 29.75696 31.06058 34.09193 30.58663 30.63648 26.38132 33.26107
#>   [49] 33.72351 28.89414 28.81235 34.45660 27.21405 35.80504 23.62233 36.39285
#>   [57] 32.62286 30.63052 38.36538 30.59179 31.76709 24.42631 32.42940 27.94332
#>   [65] 31.95670 30.59086 33.66586 34.16050 30.41629 32.43912 30.22467 31.87308
#>   [73] 29.42990 26.36823 31.37007 35.24767 32.88315 37.49443 33.21320 31.58591
#>   [81] 28.17251 29.28279 27.94090 27.53587 22.64630 29.67478 31.03550 24.43461
#>   [89] 23.89982 28.73166 24.28927 30.12479 30.61451 27.59145 24.11875 25.22512
#>   [97] 28.73173 25.75604 26.52871 29.34037 30.79959 29.00267 31.66407 29.59968
#>  [105] 31.32354 23.79545 23.24236 29.29133 27.23008 24.28183 27.40107 26.94566
#>  [113] 23.72706 29.29051 28.44184 32.86791 23.56873 27.02717 31.82718 30.13235
#>  [121] 34.68896 25.79315 28.96480 25.49883 29.22988 29.18445 29.70773 22.65764
#>  [129] 22.19433 28.00686 29.78240 30.21993 27.61153 28.72606 22.46880 28.40384
#>  [137] 28.91301 34.60935 29.13162 26.36402 33.07319 35.71893 29.03976 31.96479
#>  [145] 25.87631 31.35087 27.85715 30.77377 32.08392 30.98028 24.13565 22.79700
#>  [153] 26.49976 26.25374 29.71661 34.83825 26.86486 27.30125 30.56520 31.75307
#>  [161] 28.32695 37.03270 24.87847 31.26367 30.19735 35.31738 30.31389 28.00816
#>  [169] 31.91307 34.34816 32.37130 30.41881 31.52247 31.76187 26.64206 35.23894
#>  [177] 29.29932 29.41855 28.72527 28.32957 27.12132 33.39836 26.15070 32.87560
#>  [185] 27.78808 35.58491 30.91335 26.16286 33.48411 24.65714 29.90889 29.23545
#>  [193] 30.54290 30.93966 35.07583 33.60246 27.97248 31.43205 31.18246 33.96377
#>  [201] 29.48825 33.79005 31.02366 30.16263 30.53603 24.23088 33.10957 26.53903
#>  [209] 38.55998 36.06990 36.92564 29.08361 30.26496 25.52681 31.78559 26.49879
#>  [217] 26.93280 35.95393 29.90640 27.45408 18.53881 33.57046 29.15416 32.06464
#>  [225] 28.64548 34.24590 31.01179 33.56307 27.60008 30.20771 31.93749 32.47913
#>  [233] 29.40840 33.98762 31.22615 32.57504 27.98310 27.57698 23.18954 32.47489
#>  [241] 26.85306 29.81341 31.89380 27.49312 29.81597 32.67267 25.32521 30.05411
#>  [249] 24.16911 30.75013 32.26136 29.00279 23.71319 31.18977 27.77741 31.78173
#>  [257] 30.54177 27.38959 29.41399 29.40546 33.41660 30.60024 31.31818 24.24524
#>  [265] 34.58125 27.77408 29.82058 29.11433 24.40021 25.37856 29.49607 31.34446
#>  [273] 31.24895 35.52626 32.73563 33.07189 30.13452 29.52597 37.68578 36.00023
#>  [281] 32.53868 24.90505 30.60906 33.24831 34.74594 35.49934 32.28488 29.14488
#>  [289] 32.96294 31.05163 28.50061 32.32265 31.74572 37.68320 30.07875 35.08895
#>  [297] 28.13715 31.75433 30.50790 33.18173 28.97710 34.94855 30.16299 27.38200
#>  [305] 26.38510 31.10962 32.88455 36.54009 28.48788 31.96133 16.02706 28.96089
#>  [313] 29.83610 34.92716 28.38037 25.42926 30.10195 33.60530 30.71691 25.38119
#>  [321] 27.09458 34.84458 26.69805 28.39124 32.10221 31.35235 30.34937 26.39897
#>  [329] 29.81681 32.34490 32.58683 32.97380 31.64294 34.37482 29.16010 31.08180
#>  [337] 23.21760 35.48714 26.36782 31.82597 29.02733 29.18929 28.12595 35.64810
#>  [345] 27.38427 32.45252 31.50566 31.98078 25.79272 28.86736 29.24295 30.29790
#>  [353] 29.68547 31.62449 29.71901 28.36965 26.98726 29.58352 30.98575 27.81913
#>  [361] 33.98856 33.20139 28.17888 28.36706 28.73566 30.70239 29.01842 28.72703
#>  [369] 24.07542 30.54721 33.14132 28.96308 31.63655 27.35273 34.15647 34.91890
#>  [377] 32.65940 33.82047 26.24614 28.88152 33.50319 30.78071 25.81646 31.91042
#>  [385] 30.61648 30.17400 30.79650 32.62894 27.99665 29.47883 24.83560 31.42054
#>  [393] 28.69068 30.94066 28.82787 31.97421 20.76356 29.15718 31.01526 30.55694
#>  [401] 30.12346 28.50135 32.40837 34.30895 29.84031 35.20310 30.44417 33.42096
#>  [409] 25.69044 33.10861 26.90999 32.13782 26.84654 32.77269 27.55304 29.93043
#>  [417] 31.85512 33.94445 27.67002 33.11114 32.39568 32.08604 24.45951 30.95670
#>  [425] 31.84770 31.66304 28.43836 35.40136 32.93097 30.83831 24.97976 32.45467
#>  [433] 30.34155 34.70707 32.56330 34.83921 31.10564 29.27521 33.38856 34.01399
#>  [441] 29.24121 25.49456 23.91121 32.04236 34.37075 31.63689 33.11578 29.41562
#>  [449] 35.34051 35.82936 32.17260 36.58746 27.22990 30.06217 31.44118 29.64681
#>  [457] 28.72375 37.47660 31.24540 28.76275 28.16214 28.54520 29.27352 33.67100
#>  [465] 27.12897 27.37106 24.24223 31.36720 25.48816 31.19246 28.05866 30.37988
#>  [473] 26.14531 33.31704 34.50000 30.02383 34.47540 35.54498 33.36752 31.31206
#>  [481] 25.54801 29.86543 32.38885 29.90691 24.71257 24.59812 27.05226 19.95398
#>  [489] 24.94849 27.39032 26.34274 31.30283 32.21998 30.89774 31.00447 31.85196
#>  [497] 32.52501 30.39948 29.64993 31.65935 30.54537 26.93218 26.64951 29.61683
#>  [505] 34.88189 33.70448 27.88164 33.19430 28.37966 23.75208 26.72283 33.72969
#>  [513] 25.67939 29.22007 30.22056 31.93864 33.75055 34.23325 29.54634 30.35737
#>  [521] 28.22635 23.85973 30.77772 25.94234 33.30611 30.16670 30.29914 30.10181
#>  [529] 31.09186 32.55771 33.48010 32.38426 34.11870 30.87236 29.15594 29.56065
#>  [537] 28.56190 31.06310 30.50393 32.28470 28.46395 27.78619 26.31032 29.49692
#>  [545] 27.49688 24.24213 28.80882 33.80255 32.10699 32.98341 30.98523 26.16114
#>  [553] 31.47051 18.36252 32.83996 30.09050 28.57941 31.20169 33.94084 24.51343
#>  [561] 33.22862 27.89727 27.96597 25.79381 24.75987 32.28760 35.29631 34.95971
#>  [569] 30.89596 32.05616 33.16399 31.04437 30.56577 36.18548 31.20302 31.12047
#>  [577] 35.28393 37.04272 31.19897 33.23624 25.66364 29.27947 35.06549 34.10344
#>  [585] 33.10565 30.40646 33.95150 31.64013 36.05797 30.07229 38.73655 28.73821
#>  [593] 32.46353 29.94180 33.35871 35.58970 33.89807 23.30106 30.08709 38.28963
#>  [601] 32.39836 32.55360 33.21784 33.64960 29.28402 31.34206 30.88668 31.26272
#>  [609] 24.19938 32.10573 33.32280 28.45695 26.44831 26.92394 34.77521 27.95635
#>  [617] 31.74703 29.70190 28.25150 28.19021 21.23555 30.48996 28.80099 23.68593
#>  [625] 32.61316 32.81991 36.38597 27.85448 26.69704 30.26415 26.08684 30.12026
#>  [633] 28.39350 24.82419 28.57694 34.53464 31.67262 28.23203 34.72235 37.45500
#>  [641] 31.23099 31.39995 30.61120 36.14592 35.76495 27.65151 25.05541 32.64081
#>  [649] 31.90802 33.34981 29.03226 33.34335 34.07181 25.71795 35.52507 27.93755
#>  [657] 29.75271 28.16766 35.09008 28.98866 34.50494 27.95734 28.16254 33.17069
#>  [665] 32.69154 26.97117 29.14360 31.21906 29.42958 26.39832 29.66903 35.47997
#>  [673] 35.70408 29.98286 30.22633 29.68114 33.76503 28.48162 30.10503 32.68036
#>  [681] 26.63346 32.00641 29.15118 28.87942 36.42863 25.83489 27.96800 31.49089
#>  [689] 31.75630 32.77030 25.09038 36.99555 25.73352 28.56112 32.87844 33.86635
#>  [697] 34.83733 25.47601 29.88697 32.36031 24.93365 30.11390 32.40075 30.67889
#>  [705] 33.83466 28.83984 30.55729 25.03199 32.65452 29.14524 35.82727 30.05742
#>  [713] 31.16708 31.49121 30.58485 31.92193 28.15600 30.22222 28.82653 29.81993
#>  [721] 26.52586 35.27134 30.72118 23.73053 33.36747 26.17387 33.39506 35.01188
#>  [729] 20.66932 30.15495 34.77418 33.11368 32.31767 30.07767 30.28040 32.97070
#>  [737] 36.68407 34.94686 25.39124 27.74910 32.50868 31.31956 24.92614 29.94433
#>  [745] 30.46395 30.33595 29.33371 26.89837 23.06150 29.99157 34.57331 29.73114
#>  [753] 26.36646 31.32399 26.09617 25.46206 23.04143 21.57641 32.13765 32.97519
#>  [761] 26.05564 33.05702 28.97171 26.32201 28.91934 22.09764 29.59347 27.97132
#>  [769] 22.17852 27.74361 27.37645 29.74761 29.96123 30.38418 28.23408 31.14573
#>  [777] 32.88578 33.69233 25.63149 35.45431 27.45353 28.97684 30.39356 31.67219
#>  [785] 33.68020 28.11938 25.50133 30.95057 29.15702 29.37442 31.33021 29.27887
#>  [793] 25.67319 26.23870 30.72385 29.66649 30.40729 31.55672 34.59633 34.57496
#>  [801] 33.54378 34.94391 32.04730 32.16882 31.79537 31.94409 34.57681 24.85861
#>  [809] 29.08219 29.70266 29.11758 29.67055 30.70765 27.45900 25.78953 30.19555
#>  [817] 31.67248 32.30407 29.10321 25.70891 30.10056 30.64572 31.59663 21.99779
#>  [825] 30.48941 32.54122 31.61625 29.71898 31.14517 30.48786 28.39308 29.45273
#>  [833] 27.02834 27.48906 30.22868 25.12069 27.16161 30.70092 34.50046 30.98237
#>  [841] 29.54167 26.73725 26.32687 33.60429 27.78657 26.19267 30.47782 27.33660
#>  [849] 30.91765 22.61343 38.72983 34.12948 29.27773 26.17659 27.39811 29.10693
#>  [857] 33.84697 26.72354 31.74825 30.43661 35.96244 32.44211 32.23445 28.75347
#>  [865] 26.08718 30.12840 28.90815 29.00186 32.58361 32.44677 29.80306 31.99469
#>  [873] 37.68781 28.81846 30.85243 24.32880 26.66376 24.06227 29.54666 32.31622
#>  [881] 36.05812 25.37744 31.09036 30.39307 35.19749 29.12561 28.16882 31.69513
#>  [889] 29.05175 29.33301 27.37379 30.68027 26.89969 30.72657 29.54996 26.59086
#>  [897] 22.34022 26.31093 29.47912 27.12485 31.53211 22.62563 28.75832 26.42867
#>  [905] 21.40462 28.97065 30.83805 21.65346 32.09921 25.99807 29.68078 27.05398
#>  [913] 27.78937 26.99611 27.09540 32.06489 28.88673 32.77512 28.78996 33.55801
#>  [921] 28.43140 30.31530 25.57720 33.20354 25.01018 30.57631 30.24766 29.53088
#>  [929] 25.80761 27.07783 26.35380 31.92533 31.24050 31.28884 24.75385 31.10334
#>  [937] 28.77568 27.16414 24.97390 34.52763 27.62235 32.68034 24.56199 31.30074
#>  [945] 27.78544 23.11530 26.54642 32.59705 22.38169 30.01571 33.04201 25.66927
#>  [953] 23.41484 26.77317 30.80851 29.13956 28.22713 23.83134 28.58120 32.52520
#>  [961] 29.92817 32.78161 31.35646 25.81054 22.83411 31.63316 30.32848 30.77034
#>  [969] 31.28030 25.86239 28.31718 30.19081 28.89207 27.98056 23.12234 22.86422
#>  [977] 28.50649 24.79982 31.94623 32.24215 35.35796 23.25161 30.92221 29.69235
#>  [985] 30.11900 29.10248 19.99883 29.73723 26.94442 33.73484 28.54224 28.08793
#>  [993] 28.82197 24.71264 30.92301 29.53926 27.20839 23.40224 30.77554 32.21397
#> [1001] 30.42538 26.38119 29.44432 23.62186 25.53197 26.40067 28.00597 29.32133
#> [1009] 25.76643 28.32879 28.23671 31.20927 29.75054 30.17591 29.17352 21.68659
#> [1017] 29.76561 29.15173 26.43643 28.62738 29.94286 30.44666 29.67382 25.70948
#> [1025] 27.19415 35.00358 29.46758 31.27202 26.58419 32.85957 27.81984 27.37743
#> [1033] 35.27583 25.71492 27.10096 27.11394 32.87926 30.37575 30.94807 33.26155
#> [1041] 29.54624 27.73475 33.34337 33.04502 30.92989 31.82517 29.33004 29.93333
#> [1049] 31.74639 25.95691 28.98324 32.66376 31.24311 30.97163 27.60978 28.82629
#> [1057] 24.40585 29.12244 31.70330 29.54923 31.13241 29.12961 31.38245 29.53339
#> [1065] 35.33448 35.59724 32.28801 31.66649 28.57933 34.12791 34.47854 23.66641
#> [1073] 33.17377 25.48960 32.07869 31.55392 32.67288 32.31239 30.39076 29.87890
#> [1081] 31.48290 29.70977 25.26867 30.12036 32.27137 31.69444 34.32096 36.27087
#> [1089] 26.52909 33.04535 35.02551 28.74391 26.81081 34.44160 31.15922 31.23130
#> [1097] 36.72653 24.56019 35.56760 29.06188 35.26283 36.39874 36.93443 34.17327
#> [1105] 30.49589 31.92736 35.51029 35.35351 31.73201 36.36480 31.19365 34.05034
#> [1113] 30.22679 30.46293 33.33713 34.29445 31.36500 30.03331 29.85327 28.81209
#> [1121] 27.48276 33.36410 33.21368 33.41791 37.52282 34.14391 34.81473 29.24209
#> [1129] 35.46772 38.15037 32.67046 35.16668 34.22108 36.85252 33.50477 33.02579
#> [1137] 33.75119 30.18198 32.08488 32.58511 34.71370 32.39412 34.20164 30.89155
#> [1145] 27.32951 22.52426 25.19311 33.84177 26.15989 31.75307 31.66906 30.13382
#> [1153] 33.70608 32.91663 29.45701 28.43412 28.45376 37.61675 36.71272 30.30503
#> [1161] 26.62488 35.05417 28.40182 32.62143 32.39187 24.53923 30.98378 27.39942
#> [1169] 30.28035 27.06981 30.98288 33.05477 33.94450 27.48031 24.46109 31.37863
#> [1177] 29.14657 29.46055 24.90370 31.95890 29.77593 35.28391 29.03418 28.46606
#> [1185] 22.71867 31.22270 33.59364 32.30269 31.83024 26.14934 31.64763 33.52504
#> [1193] 29.51572 30.32201 28.82645 27.55365 30.93843 28.79487 28.66119 30.18534
#> [1201] 27.89218 29.20437 28.35982 30.46356 31.83305 27.53639 31.76082 24.36663
#> [1209] 30.59865 27.29047 28.61495 23.95401 30.38596 27.49410 23.37514 30.48948
#> [1217] 24.47801 26.94101 24.41923 30.42689 30.66820 27.65751 29.01500 29.17585
#> [1225] 21.94257 32.25471 31.99839 30.26438 26.83868 28.15724 22.39992 26.34959
#> [1233] 27.86997 31.06035 30.98052 28.02694 27.67858 35.20266 28.23899 26.72421
#> [1241] 21.54320 28.56127 29.46829 34.62506 26.05125 30.92699 29.44332 27.97258
#> [1249] 29.63628 33.90756 22.19389 22.77526 24.10917 32.03011 31.81681 32.05460
#> [1257] 23.94959 28.32859 33.85634 30.50694 26.64539 30.99850 23.42536 24.02311
#> [1265] 30.20776 23.98482 26.09279 30.97448 29.36593 27.31916 29.40585 26.44733
#> [1273] 27.18254 30.12337 25.36246 29.69925 26.25811 27.98836 23.42852 26.04063
#> [1281] 35.35041 31.57094 28.99473 24.70526 30.48526 26.41986 28.41534 33.48870
#> [1289] 24.69653 25.67193 26.92643 24.64558 29.37905 31.72945 27.38443 29.22657
#> [1297] 20.77175 21.24360 26.76021 26.77193 31.61556 23.50852 21.46368 30.80984
#> [1305] 28.64177 25.02314 31.10298 24.71167 27.23819 29.66533 31.02165 24.30435
#> [1313] 23.41952 29.34502 26.38217 30.94320 28.74931 29.55391 29.39029 29.80181
#> [1321] 24.96183 24.96151 30.26398 24.32419 20.74871 31.60415 31.14353 27.02298
#> [1329] 26.12410 27.77966 32.82705 26.12064 28.68221 28.61953 30.60242 31.46663
#> [1337] 31.44089 35.19330 29.78038 28.24663 33.96529 31.20991 34.37056 29.53850
#> [1345] 36.41003 27.71246 30.40485 33.03503 31.93128 27.01634 28.61442 32.60997
#> [1353] 34.96189 26.93092 22.04390 31.66098 25.44469 33.05831 26.60470 29.33623
#> [1361] 33.89716 31.13006 34.99737 27.10560 33.47929 29.33202 23.72104 26.29943
#> [1369] 37.30162 25.21814 35.74925 30.24589 27.09840 36.63671 32.50381 27.08467
#> [1377] 31.56524 29.38100 28.12910 29.48757 25.81722 31.81480 28.25372 29.24623
#> [1385] 31.24903 29.83132 26.85592 33.02309 35.35059 32.24583 27.78366 28.47482
#> [1393] 26.90863 32.74455 25.04129 33.45743 34.38893 23.29130 27.50739 31.59350
#> [1401] 30.65797 24.91015 26.22418 26.13982 22.42561 31.82867 29.81489 29.44845
#> [1409] 28.91512 30.51222 34.51384 30.52281 32.36355 29.41272 34.11889 30.29019
#> [1417] 31.16140 29.92076 32.15437 28.84782 31.21876 33.43392 31.93476 33.39442
#> [1425] 30.40201 34.83774 27.27337 31.34385 33.89733 32.42965 23.56698 27.98960
#> [1433] 29.39883 32.94141 32.50814 35.10979 33.01162 30.57315 28.73539 34.34239
#> [1441] 32.75460 31.72152 29.22543 24.78605 35.15303 30.98081 30.80638 33.39097
#> [1449] 36.88607 24.92932 25.76449 32.54609 25.98480 29.68520 31.53021 30.36573
#> [1457] 25.48065 30.68213 27.67347 29.11980 28.99521 27.26706 30.78951 31.63267
#> [1465] 26.92272 27.42645 28.10681 20.16443 30.91642 25.29018 23.94721 25.74737
#> [1473] 33.57403 26.51429 22.07255 29.57326 32.03873 30.90085 28.40000 26.35773
#> [1481] 29.01871 28.50541 23.29511 26.83609 26.43200 32.46594 29.12429 27.85976
#> [1489] 29.53667 26.77553 24.96343 27.22880 25.37597 35.08652 31.15134 28.43166
#> [1497] 28.52987 29.58638 32.29727 25.78264 21.60994 25.86786 28.43499 28.07223
#> [1505] 26.61975 25.85246 29.47965 22.91960 30.85095 26.40992 32.80451 30.73229
#> [1513] 33.04238 26.35228 34.20750 27.27067 35.93290 29.64373 31.73114 29.82665
#> [1521] 33.82526 32.19808 28.07129 21.22939 23.96584 30.61652 28.52418 29.52543
#> [1529] 26.51221 32.76471 27.47338 31.12702 33.83625 33.78192 33.63486 24.84083
#> [1537] 27.59530 28.37794 33.71983 31.40044 32.85295 29.00642 31.41405 25.13609
#> [1545] 31.49455 25.14382 28.71569 31.32987 30.07929 31.63470 31.33556 28.07427
#> [1553] 22.29943 32.61324 31.34039 34.36928 28.17639 23.90995 26.40019 33.43548
#> [1561] 27.06920 28.46529 26.37983 30.75253 30.75117 27.00625 30.27831 35.54600
#> [1569] 30.09516 32.58800 35.63852 30.31916 32.02236 30.06390 26.33059 30.68416
#> [1577] 28.83273 27.92015 27.55753 30.65062 32.69129 33.47500 30.06006 29.66684
#> [1585] 27.69968 30.41916 29.29696 27.68348 30.59183 30.58959 30.24671 29.20842
#> [1593] 31.17883 26.35698 31.21593 25.34208 34.04444 28.16454 25.90798 27.91062
#> [1601] 28.11238 29.28732 27.60624 32.02318 32.83962 30.94752 33.06881 29.05926
#> [1609] 27.12089 32.33811 33.33030 29.26765 32.05372 24.83049 30.23199 32.17727
#> [1617] 29.08979 28.20009 28.49931 30.00645 27.66223 29.36738 24.89107 20.10984
#> [1625] 35.68567 31.79451 33.58387 32.33474 27.86111 35.14279 30.22960 30.04843
#> [1633] 28.50686 33.81671 29.20568 32.90429 27.72853 31.85846 36.16088 31.24928
#> [1641] 29.11568 29.93856 29.27697 32.98393 27.41714 29.50510 28.99528 31.42626
#> [1649] 22.15957 37.12092 25.17110 34.35905 31.83247 34.12555 28.50736 33.15121
#> [1657] 31.61010 29.21921 35.79731 31.62664 31.36701 26.21158 30.08042 31.83769
#> [1665] 35.20313 33.18555 31.71591 27.99155 25.11880 26.32699 36.22291 31.66919
#> [1673] 29.16377 31.64003 32.08815 32.69290 29.49074 26.51875 33.35659 31.72077
#> [1681] 36.06694 28.12100 31.68401 32.53847 26.33665 31.84002 29.37048 30.70130
#> [1689] 28.05079 31.89581 29.08690 25.33816 27.96441 29.78145 30.26299 29.31562
#> [1697] 30.65119 34.88607 32.95182 28.67200 27.26550 32.04426 27.84510 33.19182
#> [1705] 23.27703 34.42257 33.78389 30.84350 30.77646 28.41623 27.40226 27.59733
#> [1713] 27.93040 24.22621 30.08839 32.25236 25.86761 30.71202 32.79264 28.54182
#> [1721] 33.22955 28.22021 25.13373 22.91442 32.50425 30.21523 28.96031 30.45039
#> [1729] 32.14129 30.35673 28.25432 28.38180 23.54513 34.05751 29.15394 31.09793
#> [1737] 28.26652 26.51313 30.61718 30.06128 27.27076 34.67646 32.18240 28.52333
#> [1745] 32.34595 31.87732 31.35663 31.79323 30.24232 27.46584 30.03248 31.28752
#> [1753] 34.65358 27.74538 27.20228 33.79059 33.53676 31.00843 31.02245 33.24477
#> [1761] 31.05519 24.20960 32.75169 25.63083 28.29943 35.56981 30.81514 23.71845
#> [1769] 29.21749 33.65687 36.54051 33.36132 27.30465 27.70647 34.48826 29.39016
#> [1777] 29.95574 30.31225 32.55127 31.70054 28.17389 26.17432 28.36850 29.79438
#> [1785] 30.76585 32.73589 27.81080 28.87431 34.07011 24.12681 30.94104 31.51349
#> [1793] 26.58757 32.17425 27.92518 26.45756 30.73840 33.91306 35.37138 32.08175
#> [1801] 36.16169 32.04250 35.96383 33.12003 35.07162 29.68612 31.65247 30.05582
#> [1809] 34.52088 28.88454 30.98805 36.16637 29.38932 31.83250 30.18519 31.24204
#> [1817] 36.10454 30.57852 33.10499 33.20114 31.49809 29.48569 32.50373 32.57642
#> [1825] 33.55025 32.54618 33.00270 25.65630 32.16372 35.17130 29.06572 27.85762
#> [1833] 30.64771 26.82534 31.61734 27.22442 31.39902 22.75955 25.48865 24.17591
#> [1841] 23.74226 26.14508 28.37780 28.49332 27.39713 28.24852 26.21106 27.18864
#> [1849] 33.67764 28.44124 28.53325 23.23420 27.48651 28.86551 23.25231 27.52950
#> [1857] 26.90813 28.30708 23.17347 32.76209 20.82276 31.19674 27.98475 29.99222
#> [1865] 29.75185 28.58012 25.09298 24.10564 33.15423 27.99654 29.52061 25.05997
#> [1873] 31.10916 29.28041 24.01074 29.25749 23.73621 27.06328 25.69761 22.22674
#> [1881] 29.64839 30.37471 28.89486 33.51740 29.37721 32.86645 28.62649 29.21906
#> [1889] 31.75695 30.84644 26.10925 30.83450 29.89536 26.93670 24.90434 25.00080
#> [1897] 32.70101 34.98743 31.53200 31.97257 36.91349 30.80633 30.49122 27.76277
#> [1905] 36.62291 35.52059 27.11993 34.15777 34.67837 36.23939 34.17617 27.95511
#> [1913] 31.39099 31.24300 33.39315 35.68948 26.90776 26.87022 31.98514 36.06449
#> [1921] 28.74438 35.21557 33.85617 31.85898 29.30739 29.71441 33.18393 35.19028
#> [1929] 27.43279 32.28137 35.93164 31.94967 31.87139 34.63351 33.36152 30.66737
#> [1937] 30.54914 31.91956 26.35963 24.72498 34.57121 28.44305 34.69612 34.23770
#> [1945] 28.10860 28.38978 28.75036 31.07023 30.92097 33.16645 31.15822 34.66653
#> [1953] 36.07820 30.08983 30.16063 30.08768 31.12621 32.22015 32.26464 32.16651
#> [1961] 26.50798 34.79235 36.28032 24.84350 30.67932 28.27333 32.90481 32.09765
#> [1969] 31.15224 27.72436 29.46913 33.71358 33.27082 31.30624 32.55808 31.00539
#> [1977] 31.97893 28.60688 30.44782 25.99582 31.73287 33.90860 33.94955 27.72738
#> [1985] 32.95205 25.01380 23.61664 31.66141 25.26007 26.80295 36.35479 27.22530
#> [1993] 29.75061 33.90989 35.93864 27.20209 30.80122 33.31317 32.94301 32.94560
#> [2001] 32.95184 32.10590 33.35952 34.37440 29.64734 31.24411 28.60223 32.63652
#> [2009] 30.48050 30.40384 32.33407 29.96564 32.49994 30.19998 29.95612 33.80478
#> [2017] 34.10057 36.04640 32.56153 35.99096 31.88341 28.85458 34.73626 32.09567
#> [2025] 26.79594 30.75135 31.37915 32.94764 25.94884 28.14578 36.05807 32.11937
#> [2033] 31.53146 25.00276 19.15695 21.47600 27.35450 28.43919 34.02830 32.55868
#> [2041] 26.74879 27.70836 32.45312 29.82796 33.00272 29.40650 28.34499 34.27646
#> [2049] 28.35874 29.72733 27.13790 31.66853 32.51514 33.50854 26.34631 34.76316
#> [2057] 28.50161 26.63561 32.08677 26.02248 30.30875 29.09788 30.90225 31.67458
#> [2065] 31.62160 34.04853 31.42361 35.60812 25.54723 29.55069 31.14271 33.71690
#> [2073] 26.99885 25.09885 29.88613 25.90767 28.39721 35.23274 31.56859 29.58969
#> [2081] 26.86970 30.10250 31.27591 27.90384 28.46773 31.69377 36.20929 31.38172
#> [2089] 27.73665 32.34100 28.67072 32.63775 27.27896 30.74320 31.47013 29.45912
#> [2097] 30.06968 26.93833 28.49634 29.87414 31.81170 32.87083 25.36737 30.01626
#> [2105] 31.64518 33.57102 26.35983 27.73698 32.85923 30.29177 32.65690 25.48040
#> [2113] 33.05066 29.90505 31.25249 33.52792 28.99061 26.82045 32.77020 26.61893
#> [2121] 27.80604 28.95288 30.28707 29.38264 29.19586 28.56056 27.40674 27.80162
#> [2129] 29.75423 32.49967 28.26640 28.57421 27.13612 29.05548 32.13598 32.19767
#> [2137] 28.22840 27.15631 33.25460 27.21051 29.96368 30.21321 24.99718 31.66271
#> [2145] 27.60093 24.77604 31.86368 30.10983 26.24452 30.87979 26.34071 24.68556
#> [2153] 26.10499 29.51307 31.09730 32.53589 32.57452 34.34117 31.11686 24.65051
#> [2161] 33.47149 21.27709 30.25548 27.78071 27.96073 27.53218 28.74181 27.19958
#> [2169] 32.39420 27.45466 27.55822 29.99536 26.42625 30.10578 29.02860 24.52215
#> [2177] 28.27115 30.36195 26.41899 29.06387 28.05584 32.77225 30.95006 26.52261
#> [2185] 26.87326 25.57186 31.45340 32.15900 29.08458 29.55948 26.39142 26.46862
#> [2193] 27.34226 26.03198 29.02741 23.54403 30.19658 28.92785 31.59733 29.75262
#> [2201] 30.75509 23.75706 30.01239 26.85972 26.20333 28.56434 32.36703 26.33664
#> [2209] 25.18015 27.29374 23.69033 29.61691 31.63938 27.90694 22.39809 32.81504
#> [2217] 30.80183 36.09421 28.70166 33.23147 29.20487 29.78309 32.42993 28.64520
#> [2225] 30.39608 30.06552 27.24874 32.75286 29.54758 32.70196 28.46853 31.90959
#> [2233] 28.91077 34.82282 29.60778 24.93458 29.80264 24.53939 33.28366 26.71323
#> [2241] 35.11649 29.84571 33.41224 27.47530 34.53010 34.24666 25.13286 26.10322
#> [2249] 31.24774 30.22934 29.76987 29.19001 33.33344 27.39285 26.43830 35.89529
#> [2257] 33.13308 33.51593 28.86443 33.32573 31.36066 30.32917 29.22408 34.72568
#> [2265] 31.41324 32.64457 30.03782 31.88195 31.39575 23.66709 30.14455 36.64125
#> [2273] 30.46562 31.42662 31.43787 29.84500 30.27174 30.77045 32.07782 22.87964
#> [2281] 33.47710 27.81120 32.59496 30.85240 28.32739 29.05034 31.16205 32.22985
#> [2289] 28.07612 30.44185 33.04829 26.73115 26.35536 29.73485 35.28664 32.90394
#> [2297] 29.95094 31.54829 35.85722 29.67409 26.45913 31.69582 28.77706 31.94608
#> [2305] 33.66851 30.52462 31.25218 24.53625 34.29779 35.42995 37.06838 27.31688
#> [2313] 30.49305 33.33887 34.75203 30.49277 33.18572 31.90857 30.43677 34.00096
#> [2321] 26.03054 25.18742 30.51158 33.24653 28.10581 29.38821 31.89769 32.85746
#> [2329] 27.32557 28.11666 32.71567 26.08547 35.39789 26.37537 30.77687 32.81209
#> [2337] 34.80287 33.71384 28.33278 27.99610 25.13609 33.46496 30.78008 30.71102
#> [2345] 27.17419 24.03887 25.11452 31.17237 30.92253 33.66481 29.39176 23.89035
#> [2353] 30.50106 27.98768 32.18257 31.13228 31.80207 27.74397 34.21655 30.10912
#> [2361] 30.82739 25.06800 25.09524 35.67692 26.38388 28.16647 26.95379 27.02893
#> [2369] 22.20004 24.67546 31.02791 32.18956 26.24756 24.54006 23.33585 34.33125
#> [2377] 31.92272 30.97561 34.73350 30.70495 27.82597 23.86793 30.13027 25.58788
#> [2385] 33.86697 23.84409 29.56912 29.93207 29.56432 34.34955 27.00563 25.39003
#> [2393] 27.91790 30.66970 26.48827 24.68309 28.13227 34.46505 26.95487 29.32736
#> [2401] 26.14910 31.02857 27.81557 29.03900 26.03428 31.89741 30.89622 21.18032
#> [2409] 23.24921 31.32901 30.90564 27.10073 21.80346 27.86762 29.55739 29.12462
#> [2417] 27.41173 27.86284 26.62482 24.64657 30.67913 29.25283 33.25762 22.47766
#> [2425] 28.56628 36.74914 19.96037 31.24100 25.66859 27.55748 31.83870 32.31130
#> [2433] 24.65892 30.10325 30.88256 31.18980 26.00934 25.29450 27.13881 27.93665
#> [2441] 25.92347 27.46151 33.05660 31.46797 28.47360 29.64216 30.75383 26.53455
#> [2449] 24.85798 31.04224 21.87816 24.84690 31.74036 27.43003 26.27377 29.36432
#> [2457] 30.75778 25.75158 28.15013 30.51751 22.09541 24.75845 34.73491 22.86148
#> [2465] 25.35238 19.04551 22.63706 25.43527 31.37444 24.37816 34.63171 23.96526
#> [2473] 28.23927 25.98919 22.27830 26.02688 25.71098 22.39962 22.32086 24.91287
#> [2481] 25.39076 27.03877 21.03549 26.41914 30.42817 28.92189 26.08654 23.65220
#> [2489] 20.15856 23.61573 27.25635 22.82890 22.24764 27.25196 29.23155 30.10704
#> [2497] 29.57729 26.89071 32.10003 27.60741 31.87137 29.96250 29.10049 28.51756
#> [2505] 30.72063 29.40422 32.19269 27.82505 33.07210 27.24176 28.72693 24.20859
#> [2513] 22.31925 25.12637 23.06302 25.12696 28.59983 32.06629 26.06536 28.26225
#> [2521] 30.36823 27.76963 24.65731 30.90437 28.62047 28.41278 24.60523 23.56690
#> [2529] 29.73496 31.66784 28.84118 18.18767 30.75359 28.54006 23.86739 32.72146
#> [2537] 25.95529 25.95043 24.83608 30.27015 25.74258 32.61211 22.42327 33.92739
#> [2545] 17.52539 30.83474 28.39906 29.46904 31.27746 27.87216 26.48307 31.17304
#> [2553] 29.34600 24.98420 26.74792 30.27254 34.52873 30.83380 23.05285 23.09923
#> [2561] 27.87534 23.49438 25.80983 33.71378 27.55604 30.17480 29.50927 35.09255
#> [2569] 29.65677 25.70572 28.00486 29.59980 30.82552 27.18184 28.54203 35.14363
#> [2577] 28.35559 33.80560 29.40886 28.74136 25.92812 29.93414 35.46590 25.74901
#> [2585] 27.46702 28.56524 22.51057 26.65263 28.19836 36.54775 26.41102 31.24282
#> [2593] 26.07002 24.83239 32.20793 31.80386 30.11665 35.15140 24.58016 26.03225
#> [2601] 28.07930 19.22169 29.07220 19.69448 22.49929 25.79361 29.81025 23.41127
#> [2609] 31.21190 21.10718 30.58690 22.72523 26.98564 25.55474 21.27684 29.13403
#> [2617] 23.16533 27.00625 25.84417 29.20666 22.64331 24.06541 23.93369 27.05463
#> [2625] 27.96531 22.41106 18.27688 25.98938 28.86405 26.33845 29.66823 28.52205
#> [2633] 34.54935 28.46204 30.88626 23.84900 23.54835 29.51126 34.35714 24.30430
#> [2641] 37.53830 26.07195 30.83451 33.05045 27.18313 27.42665 20.58134 25.28683
#> [2649] 24.96071 28.26424 31.11040 34.96527 29.04168 32.63945 35.34772 24.80478
#> [2657] 32.96313 22.24807 23.59619 30.96040 32.95403 26.69538 27.94118 29.97787
#> [2665] 29.16872 34.35800 33.14230 31.99725 28.60280 29.35724 30.21200 29.56339
#> [2673] 27.39374 28.37050 29.15226 36.64765 29.64075 29.47089 34.42685 34.93603
#> [2681] 24.61960 29.55918 33.59207 31.76048 25.68818 26.50101 33.72508 30.05733
#> [2689] 33.65909 32.63540 36.22685 33.13364 25.64688 30.19021 32.61040 26.65727
#> [2697] 26.79305 25.70605 24.31217 29.90518 23.46653 32.40478 20.78817 29.81428
#> [2705] 21.98545 27.89409 26.31846 24.83003 27.02626 28.65915 28.53407 25.64453
#> [2713] 28.72821 24.78693 25.94097 29.55767 24.17712 27.00225 26.66369 23.69985
#> [2721] 28.42624 23.64314 20.17531 29.41387 30.30819 25.53923 32.94387 31.97482
#> [2729] 27.40461 29.46048 30.79552 22.56011 29.97041 31.85438 29.94030 30.99927
#> [2737] 27.41151 29.44592 33.20075 30.50070 30.54919 20.19771 29.70374 31.74466
#> [2745] 32.95482 30.57615 28.31460 26.74601 22.78417 28.44585 30.55533 30.41600
#> [2753] 30.99592 32.29158 23.68276 25.94320 34.13174 28.00372 34.66756 23.75820
#> [2761] 25.01328 24.93816 27.45563 26.12277 21.63824 28.20751 26.09704 25.40939
#> [2769] 22.70393 25.69732 23.53379 29.38751 27.89114 29.71372 32.06220 21.20418
#> [2777] 19.72511 23.93407 28.35099 26.84992 27.96017 24.48724 29.23637 26.81980
#> [2785] 24.77032 36.21457 25.89828 27.21682 23.79688 30.68092 24.34212 27.29582
#> [2793] 28.46759 24.59417 29.66595 27.58331 31.89129 28.73299 29.79715 30.81415
#> [2801] 27.72460 28.36841 27.24088 20.14273 24.65152 27.03007 22.07023 23.75459
#> [2809] 20.79318 37.26946 22.14782 35.33223 28.48740 30.50794 24.29119 31.01111
#> [2817] 30.65762 25.88313 32.28898 29.51156 25.99843 21.56162 32.28192 22.64318
#> [2825] 28.26813 33.70472 32.03960 26.33358 26.05653 28.70519 29.75160 26.38134
#> [2833] 22.72452 32.61575 28.33949 28.58518 25.83281 28.51816 27.60572 30.67076
#> [2841] 27.90129 31.34046 35.66016 29.87084 30.27047 34.12930 34.94094 25.96296
#> [2849] 32.59512 27.47466 29.30131 29.71418 29.66217 29.76934 28.52547 29.43134
#> [2857] 31.90296 28.85247 28.94815 25.37789 30.74277 29.82941 27.28720 28.80812
#> [2865] 27.59789 29.23678 26.34240 26.37364 27.71261 29.40311 32.56934 27.55397
#> [2873] 29.03139 27.84335 30.68704 24.38932 28.45576 20.21335 24.96502 28.22598
#> [2881] 27.06276 22.41128 31.10815 24.04386 32.73113 35.32270 30.63283 33.58162
#> [2889] 32.91851 34.78696 33.85885 36.79124 31.33178 35.84366 27.66927 32.61216
#> [2897] 29.39656 26.02357 30.20852 29.67353 29.04017 35.43806 29.60152 32.05978
#> [2905] 35.71394 33.09191 29.77704 30.21280 25.30753 30.43194 26.85149 24.20291
#> [2913] 29.62767 30.47629 26.61923 28.13179 26.21518 23.26223 25.63877 30.14061
#> [2921] 30.75289 25.82808 27.62026 19.88953 26.82231 23.65628 24.76358 26.73820
#> [2929] 31.83330 31.82208 29.04788 31.75934 30.24423 28.64228 34.64828 34.44427
#> [2937] 29.78648 29.34086 33.00954 30.63024 28.72676 30.12927 35.44323 29.74139
#> [2945] 30.42362 29.65063 29.67384 23.68535 28.02893 25.09021 30.86006 25.98893
#> [2953] 29.42122 33.34465 28.31568 22.96436 24.55873 27.13245 26.00140 28.69676
#> [2961] 26.38533 32.54645 26.45400 26.99990 25.48737 25.51147 20.27159 23.96719
#> [2969] 27.03702 24.74123 31.59818 31.30827 28.22977 20.80812 31.17379 28.75991
#> [2977] 31.30906 25.85888 23.61880 31.00326 28.88115 30.48394 22.20113 26.90033
#> [2985] 31.79471 33.45456 28.62553 29.08869 30.42550 21.36848 28.94025 24.64195
#> [2993] 24.73892 27.37003 28.68353 23.78288 25.75043 30.48512 21.21396 27.46945
#> [3001] 28.69871 22.36090 31.58445 27.50081 25.74839 27.99921 29.18050 29.51247
#> [3009] 25.37456 25.42178 24.74972 31.95528 26.69242 28.41457 27.82812 22.86690
#> [3017] 26.51072 27.91821 21.89646 29.44059 28.93664 26.16550 31.30246 33.14916
#> [3025] 29.49751 32.73200 35.99293 32.17638 27.95726 33.96824 35.84471 30.55659
#> [3033] 34.83027 23.58508 27.98229 32.63018 31.07877 22.40766 27.91294 32.22901
#> [3041] 26.64713 24.20387 32.89114 25.61716 26.02944 28.68945 28.60504 26.12726
#> [3049] 26.87276 20.31248 30.97986 33.73308 27.95623 28.70024 30.18791 25.35443
#> [3057] 29.15172 26.05104 25.02768 29.58714 31.62575 27.44045 27.95541 27.79045
#> [3065] 30.51484 30.18001 26.02651 27.88352 28.41442 32.95427 24.80348 26.86049
#> [3073] 28.17397 29.84021 25.63598 32.28306 26.64483 28.42708 31.92891 27.87542
#> [3081] 29.07614 28.76925 27.54955 30.68531 28.86757 30.98139 28.60423 25.18334
#> [3089] 25.18906 29.76948 25.66583 28.91403 29.10152 26.61474 23.93527 34.67493
#> [3097] 28.26754 27.33306 23.14433 26.41729 30.91060 29.08319 25.68806 28.67451
#> [3105] 29.45270 27.86125 28.08867 24.85957 28.85903 27.09238 30.82917 28.69273
#> [3113] 26.20379 33.45370 29.43098 32.44294 26.28107 32.88150 27.51037 28.99804
#> [3121] 29.95386 31.40516 27.85880 32.12471 31.93596 28.35684 29.35325 29.60531
#> [3129] 32.45978 30.05341 32.45491 27.06031 32.31238 24.07630 30.96874 25.23005
#> [3137] 25.85673 26.20515 34.14868 23.40167 30.85719 29.46787 25.94265 29.92992
#> [3145] 25.30765 27.44674 26.41703 24.95547 26.78429 29.58710 25.87373 32.68705
#> [3153] 31.97572 20.84641 21.20346 31.34619 28.48593 25.23893 27.98316 30.81035
#> 
#> $mod2_unimodal$sim_properties$cell$knn_hubness
#>  [1]  3  4  1  3  2  2  4  3  4  2  5  2  6  5  5  2  2  7  7  3  0  7  5  2  4
#> [26]  2  2  2  7  2  1  3  2  3  3 10  3  2  2  4  6  1  7  3  5  1 11  6  1  5
#> [51]  3 11  6  3  0  7  3  7  2  6  4  6  4  4  3  2  2  3  6  3  1  3  0 11 15
#> [76]  3  3 10  0  5
#> 
#> $mod2_unimodal$sim_properties$cell$prop_outliers_cell
#> [1] 0
#> 
#> 
#> $mod2_unimodal$sim_properties$feature
#> $mod2_unimodal$sim_properties$feature$mean_expression
#>  [1] 3.130320 2.851670 2.894024 2.851446 3.539986 3.486060 3.038775 3.475203
#>  [9] 2.174359 3.660603 4.074889 3.956205 3.537117 3.685633 3.772017 2.569579
#> [17] 2.601790 1.175417 1.963698 1.747854 1.659279 2.674022 2.006106 3.221571
#> [25] 2.459400 2.290853 2.592022 2.447732 3.000076 1.331818 2.121959 2.353598
#> [33] 1.888293 2.599530 2.495010 2.338081 2.297425 1.761964 1.666367 2.797380
#> [41] 1.598760 3.173896 2.955787 2.605011 3.008812 1.978093 2.353964 3.167935
#> [49] 1.894754 2.813672 1.903438 1.657412 2.062132 2.576450 3.702718 2.509255
#> [57] 3.303397 2.650407 2.358220 2.272985
#> 
#> $mod2_unimodal$sim_properties$feature$variance
#>  [1] 11.000569 11.177434 10.952893 11.184623 10.421311 10.630753 10.916725
#>  [8] 10.580200 10.558518 10.116044  9.118071  9.579023 10.424933 10.216084
#> [15] 10.178985 11.165124 11.455502  7.937451 10.314840  9.966189  9.609825
#> [22] 11.483076 10.756192 11.118283 11.399087 11.060653 11.361659 11.283451
#> [29] 11.186716  8.487769 10.657831 11.035831 10.176719 11.426864 11.107601
#> [36] 10.902490 11.129979 10.107373  9.706134 11.326387  9.604087 11.346069
#> [43] 11.409223 11.472230 11.254990 10.466100 11.051578 11.286139 10.226605
#> [50] 11.456169 10.324347  9.593934 10.697448 11.240964 10.330629 11.269996
#> [57] 11.127542 11.262475 11.077426 10.890455
#> 
#> $mod2_unimodal$sim_properties$feature$scaled_variance
#>  [1]  0.43146888  0.66716178  0.36793417  0.67674208 -0.34045997 -0.06135394
#>  [7]  0.31973710 -0.12872210 -0.15761569 -0.74726395 -2.07717739 -1.46290549
#> [13] -0.33563364 -0.61394812 -0.66338685  0.65075759  1.03771940 -3.65048974
#> [19] -0.48234446 -0.94696192 -1.42185825  1.07446531  0.10580850  0.58833607
#> [25]  0.96253924  0.51153700  0.91266278  0.80844156  0.67953075 -2.91712772
#> [31] -0.02527009  0.47845967 -0.66640699  0.99955566  0.57410122  0.30076730
#> [37]  0.60392197 -0.75881860 -1.29351631  0.86565813 -1.42950484  0.89188773
#> [43]  0.97604761  1.06001069  0.77051418 -0.28077283  0.49944393  0.81202367
#> [49] -0.59992824  1.03860754 -0.46967603 -1.44303455  0.02752430  0.75182203
#> [55] -0.46130370  0.79051120  0.60067478  0.78048875  0.53388937  0.28472949
#> 
#> $mod2_unimodal$sim_properties$feature$sd
#>  [1] 3.316711 3.343267 3.309515 3.344342 3.228206 3.260484 3.304047 3.252722
#>  [9] 3.249387 3.180573 3.019614 3.095000 3.228766 3.196261 3.190452 3.341426
#> [17] 3.384598 2.817348 3.211672 3.156927 3.099972 3.388669 3.279663 3.334409
#> [25] 3.376253 3.325756 3.370706 3.359085 3.344655 2.913378 3.264633 3.322022
#> [33] 3.190097 3.380364 3.332807 3.301892 3.336162 3.179209 3.115467 3.365470
#> [41] 3.099046 3.368393 3.377754 3.387068 3.354846 3.235135 3.324391 3.359485
#> [49] 3.197906 3.384696 3.213152 3.097408 3.270695 3.352755 3.214130 3.357081
#> [57] 3.335797 3.355961 3.328277 3.300069
#> 
#> $mod2_unimodal$sim_properties$feature$cv
#>  [1] 114.43339 126.38105 126.00588 126.38105 101.73542 102.25243 120.36045
#>  [8] 105.02591 160.86965  98.41598  84.98863  88.74813 102.66407  97.07512
#> [15]  94.81522 135.47852 137.42713 239.54951 173.36579 191.67745 192.35357
#> [22] 135.31635 163.39325 113.33521 140.48273 153.58566 133.01331 137.13683
#> [29] 121.13303 230.57441 164.53908 149.25978 182.24578 133.01331 139.20685
#> [36] 150.83290 151.45893 180.30990 192.35357 129.62574 199.70286 117.59966
#> [43] 122.09388 129.91396 121.13303 171.03968 153.62661 114.43339 176.69758
#> [50] 126.99417 173.33489 186.76309 167.99585 135.47852  95.12104 145.35577
#> [57] 110.65487 129.50734 147.18008 153.58566
#> 
#> $mod2_unimodal$sim_properties$feature$zero_fraction_feature
#>  [1] 0.5250 0.5750 0.5625 0.5750 0.4500 0.4625 0.5375 0.4625 0.6875 0.4250
#> [11] 0.3500 0.3750 0.4500 0.4250 0.4125 0.6250 0.6250 0.8500 0.7250 0.7625
#> [21] 0.7750 0.6125 0.7250 0.5125 0.6500 0.6750 0.6250 0.6500 0.5500 0.8250
#> [31] 0.7000 0.6625 0.7375 0.6250 0.6375 0.6625 0.6750 0.7625 0.7750 0.5875
#> [41] 0.7875 0.5250 0.5625 0.6250 0.5500 0.7250 0.6625 0.5250 0.7375 0.5875
#> [51] 0.7375 0.7750 0.7125 0.6250 0.4250 0.6375 0.5000 0.6125 0.6625 0.6750
#> 
#> $mod2_unimodal$sim_properties$feature$detection_freq_feature
#>  [1] 0.4750 0.4250 0.4375 0.4250 0.5500 0.5375 0.4625 0.5375 0.3125 0.5750
#> [11] 0.6500 0.6250 0.5500 0.5750 0.5875 0.3750 0.3750 0.1500 0.2750 0.2375
#> [21] 0.2250 0.3875 0.2750 0.4875 0.3500 0.3250 0.3750 0.3500 0.4500 0.1750
#> [31] 0.3000 0.3375 0.2625 0.3750 0.3625 0.3375 0.3250 0.2375 0.2250 0.4125
#> [41] 0.2125 0.4750 0.4375 0.3750 0.4500 0.2750 0.3375 0.4750 0.2625 0.4125
#> [51] 0.2625 0.2250 0.2875 0.3750 0.5750 0.3625 0.5000 0.3875 0.3375 0.3250
#> 
#> $mod2_unimodal$sim_properties$feature$dispersion
#>  [1] 0.7038563 0.7786433 0.8931083 0.7786433 0.6856938 0.6404030 0.8691983
#>  [8] 0.7583429 0.9704641 0.7385350 0.5688166 0.5907173 0.7377939 0.6478711
#> [15] 0.6292948 0.7341772 0.7790564 0.8607595 0.9392405 1.0103567 0.8787475
#> [22] 0.8010850 0.7341772 0.7385801 0.7154081 0.8550851 0.6855860 0.6582278
#> [29] 0.7703436 1.0632911 0.9475588 0.8354430 1.1209564 0.6855860 0.7509187
#> [36] 0.8815843 0.8028933 0.7721519 0.8787475 0.7981346 0.8973277 0.8470679
#> [43] 0.7267121 0.6329114 0.7703436 0.8776371 1.0325497 0.7038563 0.8976335
#> [50] 0.7257384 0.8262371 0.7848101 0.9172347 0.7341772 0.5881207 0.9507736
#> [57] 0.7193644 0.6708861 0.7852466 0.8550851
#> 
#> $mod2_unimodal$sim_properties$feature$bcv
#>  [1] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#>  [8] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [15] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.1940640 0.0000000
#> [22] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [29] 0.0000000 0.5625440 0.0000000 0.0000000 0.5986563 0.0000000 0.0000000
#> [36] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [43] 0.0000000 0.0000000 0.0000000 0.0000000 0.2727625 0.0000000 0.0000000
#> [50] 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000 0.0000000
#> [57] 0.0000000 0.0000000 0.0000000 0.0000000
#> 
#> $mod2_unimodal$sim_properties$feature$gene_cor
#>    [1] -0.1434627339  0.0939460109  0.0247455345  0.0367083391  0.1331559886
#>    [6] -0.0440181000 -0.1446023208  0.0647105366 -0.1018869203  0.2909876448
#>   [11]  0.1141215497 -0.0129094282 -0.0178844576  0.0774000801 -0.0470479137
#>   [16]  0.1236417872 -0.0047974835  0.0768909638 -0.0446119931  0.0668484821
#>   [21]  0.1547654142  0.2256244830 -0.0800676279  0.1168209719  0.0282620835
#>   [26] -0.0800077569  0.0915505002  0.0578798975  0.1103975911 -0.1300277264
#>   [31]  0.0201037761  0.0224048611  0.2491067320  0.0292701023  0.1148489617
#>   [36]  0.0096967143  0.0505220712 -0.0925930768  0.1336181226 -0.0153975375
#>   [41] -0.2197149491  0.1196367291 -0.0088105990  0.0911378115 -0.0786185323
#>   [46]  0.0387625394 -0.0171261899 -0.1060382751 -0.1253563223 -0.0384785527
#>   [51]  0.1388122925 -0.0349574948  0.0834820931  0.1540016862  0.1269085028
#>   [56]  0.0215820245 -0.1211641665  0.0508443907  0.1285949609  0.0492792790
#>   [61]  0.0814894409  0.1737486060 -0.0637858480 -0.1307977256 -0.0848453750
#>   [66] -0.0784664399 -0.0134855370  0.1510158794 -0.0540428532  0.1539434736
#>   [71] -0.0166780451 -0.1870136720  0.0001039444 -0.1301445651  0.0336066494
#>   [76] -0.1496959249 -0.1137098303 -0.0414578094 -0.2593334015 -0.0748140409
#>   [81] -0.0979980920 -0.1587622119  0.0054175254  0.0446613965  0.0070903615
#>   [86] -0.1613491422  0.1323338997 -0.0190398741 -0.0239525272 -0.0548404950
#>   [91] -0.1478944887 -0.0835639786  0.1202934076  0.1200626813 -0.0759913352
#>   [96] -0.0966573488  0.2086244638 -0.1426989567 -0.0318022756 -0.0789055407
#>  [101]  0.0961371697  0.0835184582 -0.2003726157 -0.1546734512  0.0872698777
#>  [106] -0.1313505103 -0.0449740579 -0.0149518762  0.2907101198  0.0607407410
#>  [111] -0.1360057752 -0.2924553816  0.0865144488 -0.1449184890 -0.0239756858
#>  [116] -0.0483394129  0.0736106261 -0.0499664268  0.0081216269  0.0215067874
#>  [121]  0.0914754466 -0.0577509908  0.0557283166 -0.0140542708 -0.0283768342
#>  [126] -0.0178235548 -0.0955149396  0.1628143180  0.0669883046  0.0280940558
#>  [131] -0.0545144776 -0.0531174636  0.0937646718 -0.1699436289 -0.1011245439
#>  [136] -0.0792876040  0.0514734202 -0.0133237986 -0.0309446881 -0.1430902562
#>  [141]  0.0869601921 -0.0193720069 -0.0862109085 -0.0920997157 -0.0377946234
#>  [146] -0.0837565345  0.1837804527  0.0376665939 -0.0668757605 -0.0991663456
#>  [151] -0.0631845686 -0.0344193649 -0.0015770389 -0.2043824013 -0.0611758729
#>  [156] -0.1002963307  0.0349011961  0.1701822023 -0.1904980001 -0.2274839501
#>  [161]  0.0952312843  0.1184010880 -0.0201058085  0.1136224512  0.0417139315
#>  [166]  0.0246012787 -0.1063650503 -0.1273256036  0.0631628424  0.2305518078
#>  [171]  0.2179226394 -0.0245710659  0.0227020656 -0.1550660033 -0.0148767510
#>  [176]  0.1351976109  0.0062991868 -0.1270748751  0.0212521718 -0.0229458097
#>  [181] -0.2067279776  0.0672356410 -0.0509487116  0.0223178228  0.0193643210
#>  [186] -0.1699625932  0.0548212958  0.1207796907  0.1701016590  0.2344699648
#>  [191]  0.1408230116 -0.1616177273  0.0359632550  0.0281944710 -0.0556645238
#>  [196]  0.0306615013 -0.0800378634  0.1620670830  0.0148361362  0.0342827719
#>  [201]  0.0154002543 -0.0024396745 -0.0101299018 -0.1443086881  0.0426877968
#>  [206]  0.0407031106  0.0986407556 -0.2291645171 -0.2409294136 -0.1933547720
#>  [211] -0.0206140267  0.1744329519 -0.0044566489 -0.0093923533  0.0325264115
#>  [216] -0.1495300533 -0.0400050208 -0.0671726653 -0.0444712731 -0.2567993475
#>  [221]  0.0217112360 -0.1293669811  0.0878547386 -0.0755683610  0.0968466501
#>  [226]  0.0462113915 -0.2751300601 -0.1323178660 -0.0535436718 -0.0302477053
#>  [231]  0.0027333094 -0.0630818520  0.1906768614 -0.0291034243 -0.0178055103
#>  [236] -0.0066196327 -0.0103783295  0.0662403749  0.0678860799 -0.0727240965
#>  [241] -0.0285102757 -0.0014925614 -0.2569765412 -0.1091261386 -0.0536954634
#>  [246]  0.0569489512 -0.0709671973 -0.1510602016 -0.1828578977 -0.1259495123
#>  [251]  0.0298492812 -0.0884547313  0.2365199389 -0.0761395617  0.1674857876
#>  [256]  0.0722458942 -0.0909883944 -0.2716487977 -0.0108085930  0.0339365215
#>  [261]  0.0198581147 -0.1528619683  0.1765557077 -0.0345236278 -0.1324755777
#>  [266]  0.0183470996  0.0304302745  0.1526193339 -0.0409069329 -0.1187536139
#>  [271]  0.0549596329  0.0200322168 -0.0361155856 -0.2586438917  0.0812727909
#>  [276]  0.1632115503  0.0120229883  0.0669613744 -0.0259349056  0.0106882952
#>  [281] -0.0638950365  0.0915587760 -0.0640066382 -0.1054562106 -0.2257626798
#>  [286] -0.0334493122  0.0829262410 -0.1324264000 -0.0175706079 -0.0321870711
#>  [291]  0.1726389104  0.1117409908  0.0283524334 -0.2307912082 -0.2087402316
#>  [296]  0.0867687214 -0.0913443123  0.2376990960  0.2794965273  0.1854237988
#>  [301] -0.1551738744 -0.1115349182 -0.0715800258 -0.0645723665  0.0062471195
#>  [306]  0.0893362678  0.0647635682 -0.1487409898 -0.0059025169 -0.1346395083
#>  [311]  0.1108721066 -0.0526064433 -0.2698929793  0.1265470845  0.1054281005
#>  [316]  0.1221427726 -0.1822778709 -0.0206198144 -0.0628510027  0.0886162856
#>  [321] -0.1396277834  0.1457852643 -0.1599999079  0.0347294576  0.1531339053
#>  [326]  0.0062014729  0.0632429767  0.0146776682  0.2530857914 -0.0324180497
#>  [331]  0.1425615178 -0.1282918904  0.0130251541  0.0205302185 -0.1422702234
#>  [336]  0.1254809616  0.2094151209  0.1078997689 -0.2071813782  0.0477466415
#>  [341]  0.0835547222 -0.0537361499 -0.0359077165  0.0210058423  0.0948580181
#>  [346]  0.0719037690  0.0039657627  0.0175559194  0.0862457842  0.0009916946
#>  [351] -0.0239932146  0.1788861081 -0.0235924218  0.0476272093 -0.0578825329
#>  [356] -0.2826477912  0.0071926745  0.1438864034 -0.0431510392 -0.2358304834
#>  [361]  0.1461388469 -0.0697746790 -0.0760660330 -0.0591538221 -0.0308453295
#>  [366]  0.0241977239  0.1956687474 -0.1965148770  0.0059022061 -0.1974009481
#>  [371]  0.0175740033 -0.1486353581  0.1349572094  0.1909855178  0.2035372693
#>  [376]  0.2346718678  0.1173587820 -0.1308156351 -0.1695013959 -0.1858956844
#>  [381] -0.1810264127 -0.1072802763  0.1280005490 -0.1739632776 -0.2447773203
#>  [386] -0.3091624934 -0.1763810747 -0.1740317215 -0.0602763958  0.0173357573
#>  [391] -0.1944544003  0.1025548959 -0.0622386162  0.0955021741 -0.1764787232
#>  [396]  0.2380274148 -0.0855778229  0.1460537454 -0.0101766957  0.0655009332
#>  [401] -0.0128825990 -0.1178647165  0.0435683710 -0.0947522372 -0.0991823107
#>  [406] -0.0258328387 -0.0186668002 -0.1420121592 -0.2305735748 -0.1804353896
#>  [411] -0.0173511801 -0.0618497893 -0.0608499316 -0.0057781858 -0.1560477387
#>  [416] -0.1544072771 -0.1282428913 -0.0096710096  0.0761804498  0.0340230617
#>  [421] -0.0836203406 -0.1125667942  0.1332137880  0.0344042397  0.1328595805
#>  [426]  0.1462507570 -0.0310346012 -0.1203782532 -0.2228332896 -0.1313715357
#>  [431]  0.0137872671  0.0960022148 -0.0820356345 -0.1690789511  0.1402637433
#>  [436] -0.1110498032  0.0786991122 -0.1724377559 -0.1240217564 -0.1397418297
#>  [441] -0.3047111143 -0.1346702138  0.0046074571 -0.1811891344 -0.1549706770
#>  [446] -0.0926495776 -0.0442460028  0.0894811139 -0.1370348728 -0.2444895748
#>  [451] -0.1481136052  0.0172630736 -0.0057917666  0.1811450830  0.3227744698
#>  [456] -0.2290259804 -0.0674171500  0.0769406871 -0.0201637442 -0.1071363088
#>  [461] -0.0299591786  0.0283994003 -0.0631385904  0.1478934016  0.2734364713
#>  [466] -0.0253852211 -0.0320703230 -0.0706691128  0.0805952226  0.0089593473
#>  [471] -0.0169185213 -0.0133532336 -0.0843726905 -0.0772680487 -0.2113354180
#>  [476] -0.3233222100 -0.0399371190  0.1379393391 -0.0287672949 -0.0229626794
#>  [481] -0.1230164876  0.0337257280 -0.0489271996 -0.0060868450 -0.0540991072
#>  [486]  0.0595946414  0.0459232493  0.0094740671  0.0029999733  0.1518887963
#>  [491] -0.0687526934  0.0294686390  0.0306673939  0.0441386266  0.0894369772
#>  [496]  0.1747690116 -0.3055774565 -0.0914548876 -0.0673467645 -0.1448848085
#>  [501]  0.0433751569 -0.0764750815  0.0171447527  0.0901282636  0.0107239941
#>  [506] -0.1956568083  0.0167644318 -0.0459922671  0.0090830065  0.1121681971
#>  [511]  0.0169432789  0.0732352908 -0.1117769312 -0.1199741946  0.0935525497
#>  [516] -0.0306933067 -0.0259326969  0.1332838376 -0.0302260925 -0.1318367671
#>  [521] -0.0663413185  0.2165703212 -0.1651160972 -0.0232921000  0.0366486297
#>  [526]  0.0394730455  0.0221043560 -0.0453547786  0.0545080557 -0.1211297852
#>  [531] -0.0033626504  0.0270120897  0.0597367799 -0.0110436549  0.2754408957
#>  [536] -0.0573655216 -0.0630059281  0.0238804830 -0.0215263533  0.3026551945
#>  [541] -0.0560071223 -0.0106126594 -0.1310535644  0.1183474534 -0.1292273541
#>  [546] -0.0279712722 -0.1440719462 -0.1507697241  0.1086066803 -0.0019513516
#>  [551] -0.0709980210 -0.1200902943 -0.1392511936  0.1526253680  0.0247065679
#>  [556] -0.0937148966  0.0696550481 -0.0184801763 -0.0113575841  0.0427882721
#>  [561]  0.0539032337 -0.1348969714  0.0126572938 -0.0443757100 -0.0915579755
#>  [566] -0.1872330476  0.0350830465 -0.0837599782  0.1098256744 -0.0677996262
#>  [571]  0.0518520730  0.0781225694 -0.0103192569  0.0744865328  0.0870782423
#>  [576] -0.0228455526 -0.0491873881 -0.0918035839 -0.1958891756 -0.1964435184
#>  [581] -0.1766883212  0.0429639428  0.1615433030  0.0424130982  0.1750912941
#>  [586]  0.1159938228 -0.0380455369  0.0076476934 -0.0930194390  0.0736337365
#>  [591] -0.0700128612 -0.1410359029 -0.1148218227 -0.0149773176 -0.0008484911
#>  [596] -0.1347119746 -0.0016303025 -0.0371584052  0.0064549173 -0.0768001050
#>  [601]  0.2230253906 -0.0833467510  0.1452999315 -0.0559559724  0.1608619916
#>  [606]  0.0296601840 -0.0227216693 -0.0859586088  0.0908304801  0.0493407399
#>  [611]  0.0189152762  0.0517303755 -0.1267641155 -0.0465120808  0.0367686920
#>  [616]  0.0318720833 -0.1842258067 -0.0300130261 -0.1150350265 -0.1052465288
#>  [621] -0.1541916318 -0.1060389526 -0.0847381149 -0.0915410742 -0.0723723071
#>  [626]  0.0082857759 -0.0602757692 -0.0869401275 -0.0899810561  0.1734218033
#>  [631]  0.0924884035  0.0465830862 -0.0828443109 -0.1604685462 -0.1374577613
#>  [636]  0.1352269609 -0.0036617630 -0.0197972068 -0.0412208665 -0.0263719348
#>  [641] -0.1601144640  0.0169673769  0.0067796204  0.1092143386 -0.0653918719
#>  [646]  0.0161176398 -0.0140084721 -0.0389236836  0.1614882061  0.1231827027
#>  [651] -0.0958415626  0.1716370684  0.0828803132 -0.0156122976  0.0086453409
#>  [656]  0.0013282584 -0.1698098970  0.2864503250 -0.0976353251  0.0637248145
#>  [661] -0.0046701878  0.0252017040 -0.0666383508 -0.2279421736 -0.1211708756
#>  [666] -0.0100869087 -0.2096315118  0.0480815583 -0.0073520562  0.0437431996
#>  [671]  0.1754637014 -0.1678032114 -0.1637413728 -0.1913504436 -0.1802474663
#>  [676] -0.1272041618 -0.0184106729  0.1498702699  0.2307064710 -0.1032036765
#>  [681] -0.0193326527  0.0971243083  0.1621414294 -0.0154156277  0.2602652498
#>  [686] -0.0401158130  0.1572070405 -0.1151963890 -0.0507026081 -0.3401033681
#>  [691] -0.1137525880 -0.1859178578  0.0321888077 -0.1656394497  0.1040528400
#>  [696]  0.2335993741  0.0841300293 -0.0593266396  0.1053517273 -0.0883396740
#>  [701] -0.1128231237 -0.0009523077 -0.1370287239  0.0470956678  0.0151717785
#>  [706] -0.0019494448 -0.2964642734 -0.0837546720 -0.0717236401 -0.0321306998
#>  [711] -0.1008466717 -0.1100669960  0.0830185740 -0.3164351374 -0.1263952711
#>  [716] -0.0320886655  0.0337126350 -0.0501113413 -0.1917071433  0.0291605051
#>  [721]  0.0070009621 -0.0481979188 -0.2250803141 -0.0656080370 -0.1784325410
#>  [726] -0.0325065511  0.2693439974  0.1021716436 -0.2820723733 -0.1520458529
#>  [731] -0.0503806852  0.0963903481  0.1267289282  0.0811940038  0.2199877532
#>  [736] -0.0298311423 -0.0677383465 -0.0865561909  0.0342519794  0.0791381581
#>  [741] -0.0669250439  0.1722443649  0.0518177622  0.1223691381 -0.0326412287
#>  [746] -0.1188948292  0.1626658895  0.1057827229  0.1456364079  0.0532964369
#>  [751]  0.2600220769 -0.0179809135 -0.0265134733 -0.1008366972  0.0409123751
#>  [756]  0.0011995855 -0.1582590868 -0.1741435029 -0.0366560486  0.0688293756
#>  [761]  0.0444870393 -0.0914241688 -0.2031638442 -0.0163169910  0.1904891953
#>  [766] -0.1431912648 -0.0225227803  0.1055805318  0.0338171116 -0.1686804921
#>  [771]  0.0389638404  0.0933621282 -0.2334645636 -0.1404331507 -0.0562825490
#>  [776] -0.0345537224 -0.0335646957  0.0245663710 -0.1506387514 -0.0105080672
#>  [781] -0.1970496057 -0.1694394741 -0.0925198735 -0.1484916700  0.0990005351
#>  [786] -0.0390370116 -0.0085716718  0.1411038069 -0.0380792292  0.1790527418
#>  [791]  0.2203762945 -0.1484336260 -0.2565551596  0.0492004421  0.0710386339
#>  [796]  0.1348147623 -0.0223041572  0.2390510865  0.2732924411 -0.1176899577
#>  [801] -0.0718525360 -0.1838662160  0.0663194327  0.1652077857 -0.2471092322
#>  [806] -0.0458846212 -0.1885816763 -0.1640998053 -0.0251285385  0.0731906792
#>  [811] -0.0529817315 -0.1676779872  0.1153039389  0.0068866021 -0.0551184791
#>  [816]  0.0939381051 -0.0643852055  0.0652400540  0.0225289021  0.0317662782
#>  [821] -0.0729445420 -0.1019656360 -0.0712097164 -0.2224656010  0.0051714277
#>  [826] -0.0747648415 -0.1148096290 -0.1896824242 -0.0571543197 -0.0027274244
#>  [831] -0.0735616931  0.0133655738  0.0528244953  0.1735688040 -0.0648156295
#>  [836] -0.0882732415 -0.1790930286  0.1851935178  0.0023713398 -0.0139788008
#>  [841] -0.0226927958 -0.0507014651  0.0500481525  0.0804004703 -0.1122267198
#>  [846] -0.0920724940 -0.1837990484 -0.0594735417  0.2297906067  0.2073258829
#>  [851]  0.2053125713  0.0391159930  0.0497035154  0.0576444477 -0.1278224866
#>  [856] -0.0274334778 -0.0001101403  0.0581601569  0.2201449660  0.0904612297
#>  [861]  0.0247186807 -0.0034416486 -0.0303138883  0.1517507180  0.2433206766
#>  [866]  0.1340462731 -0.1386856853 -0.0962797784 -0.2033864271  0.1281295488
#>  [871]  0.0369652600 -0.1642930619 -0.0680578669  0.1125991326 -0.0052152939
#>  [876] -0.0217094923  0.0230851351  0.0908772893 -0.1010106653 -0.0656927884
#>  [881] -0.1182634612  0.1219325607 -0.1717280049 -0.1749420060 -0.1460203871
#>  [886] -0.0593481029 -0.1334303980  0.1305278441 -0.1687413262  0.0548986962
#>  [891] -0.0805791988 -0.0966679878 -0.0758285410 -0.2322775594 -0.1678619406
#>  [896]  0.0687384830  0.1172765764 -0.4300090759  0.1869231956  0.0220257299
#>  [901] -0.0288496936 -0.0811804857 -0.0691581141  0.0322592611 -0.1705003396
#>  [906] -0.0072611855 -0.0183827943 -0.1197510345  0.1144184635 -0.0083951069
#>  [911]  0.1637611492  0.0626523925  0.0975417734  0.0817866971  0.0698376846
#>  [916] -0.1133035851  0.1086323163 -0.0348778852  0.0039727129  0.1766871149
#>  [921]  0.0446491491  0.0385228375  0.1061824721 -0.1786129854 -0.1056097534
#>  [926] -0.1254050744  0.0971036721  0.1335339401  0.1907733154 -0.0270145980
#>  [931]  0.0786260100 -0.0887264994  0.0685747145  0.0135242597 -0.1609841251
#>  [936]  0.2128156290 -0.1264826093  0.0488548164 -0.1779628169 -0.0197170965
#>  [941] -0.1746759006 -0.0229699997 -0.0538486124 -0.0439133998 -0.1828868182
#>  [946] -0.0904189576 -0.0453573118 -0.2247285154 -0.0572657273  0.0492084036
#>  [951] -0.0540635833 -0.2973680146 -0.0536680502 -0.0411349135 -0.0146217673
#>  [956] -0.0653590476 -0.0519876023 -0.1049041985  0.0887501886 -0.0743317706
#>  [961] -0.0362177736 -0.0039592354  0.1248709079 -0.1244642315 -0.0632607841
#>  [966] -0.0496828312  0.0725394217 -0.0625712778  0.1851749739  0.0320872810
#>  [971]  0.0450109130  0.0736053527  0.0237779747  0.0288029633 -0.0552894698
#>  [976]  0.0730668225  0.0162704405 -0.0511040291  0.1468409252 -0.0474219971
#>  [981] -0.0392733863 -0.1575188909 -0.1281043850  0.1764377020 -0.0715652298
#>  [986] -0.1353199066  0.1202079120 -0.0822349800  0.1341186484  0.1207752373
#>  [991] -0.1044738972  0.0438272256  0.1278699262 -0.0193439789 -0.0739316757
#>  [996]  0.1409449437 -0.1008272297  0.0127737987 -0.1283866390  0.1040971642
#> [1001] -0.0422836891 -0.0723472730  0.0011088257  0.0357211606 -0.0082809801
#> [1006] -0.1207672921 -0.1131567603  0.0465247786 -0.2066413512 -0.0948610906
#> [1011] -0.1288685391 -0.1581099115 -0.0471427460  0.0459054936  0.0387312481
#> [1016] -0.0536850598 -0.0259032843 -0.0134015734  0.1033219087 -0.2694432037
#> [1021] -0.0234751659 -0.0070515661 -0.0113418812 -0.1277046701  0.1065786746
#> [1026]  0.1195626688  0.0689680966 -0.0747295828  0.1847641772 -0.0967568852
#> [1031] -0.1245110990 -0.0401494555  0.0605930663  0.0150929362  0.0438074300
#> [1036] -0.0895922783  0.0149718038  0.1957938808 -0.2228589927 -0.0227658629
#> [1041] -0.0923345319  0.0770540210 -0.0867093274  0.1454448215 -0.0374594516
#> [1046] -0.1881371661  0.0665929023 -0.0836460038  0.0006397793 -0.1541274920
#> [1051] -0.0784546387 -0.0080408641  0.0439971565 -0.1682151609 -0.1806344594
#> [1056] -0.0717852086 -0.2516407619 -0.0246247707 -0.1974369195 -0.1507010640
#> [1061]  0.0470897356 -0.1328611191 -0.0619284917  0.0095693786  0.0813174261
#> [1066] -0.0188048213 -0.1999941768  0.2005394874 -0.0776844079 -0.0620675026
#> [1071]  0.0605124216 -0.0794491850  0.1584987428  0.1212152790  0.0738191084
#> [1076]  0.0336457290  0.0786243968  0.1599108115 -0.0227637300  0.0218799710
#> [1081]  0.0804536446  0.0838656250  0.1200774090 -0.1677880089 -0.0330262399
#> [1086]  0.0370720049 -0.0545413024 -0.2450865418 -0.0140719534  0.0742817063
#> [1091]  0.0248557867 -0.0490099853 -0.1422014547 -0.0003531856 -0.0166367139
#> [1096]  0.0317894451 -0.0166317620 -0.1688493873 -0.0593346335 -0.1225684559
#> [1101]  0.0224315383  0.1497384157  0.0481214117  0.0987915667 -0.1429879228
#> [1106] -0.0815842310 -0.2723223813  0.0989133530 -0.1935415657  0.1190812106
#> [1111] -0.0450188582 -0.0334780008  0.0459885025 -0.1320205335 -0.1311767232
#> [1116] -0.0100067443  0.0961501210 -0.0324838804  0.0168599915  0.0767798026
#> [1121] -0.0255017551 -0.0225389612  0.1480579637  0.1015836620 -0.0936103188
#> [1126]  0.0173265582  0.1700869894  0.0323524071  0.0204870219  0.0856257850
#> [1131]  0.1452265076 -0.1263658740 -0.1391058045 -0.1315965478 -0.1485893538
#> [1136] -0.0991222759 -0.0637800958 -0.0479589181 -0.1221265259 -0.1659480478
#> [1141]  0.0971449806  0.2263616905 -0.0738429839 -0.1574648565  0.3382955444
#> [1146]  0.0464992921 -0.0326635962  0.0085904504 -0.0020294067 -0.0360304383
#> [1151] -0.2028520114  0.0625187661 -0.1338434696 -0.0436007622 -0.1584932266
#> [1156] -0.1410950255  0.0882725309  0.0960375967  0.0656819811  0.0824088473
#> [1161] -0.0976803976 -0.0972527424  0.0242274006  0.0414289742  0.0662822780
#> [1166] -0.0035267336  0.1027157532 -0.1330953914 -0.0237994000  0.1914054620
#> [1171] -0.0261201849 -0.1148689645 -0.0526265245  0.1856104015  0.0588269973
#> [1176]  0.0114621997 -0.0002224509 -0.0643292434 -0.0436346666 -0.2232206745
#> [1181]  0.0553888292 -0.1951495637 -0.0536543618 -0.0777204044  0.0165056454
#> [1186]  0.1411590900 -0.2167116417 -0.1743014595 -0.0864330308  0.0381761156
#> [1191] -0.0078923496 -0.1029418160  0.0720834752  0.0673262094 -0.0425083138
#> [1196] -0.1944810095 -0.1080117199 -0.0386180710 -0.0502809797 -0.0097178402
#> [1201] -0.1769767244  0.0436624982 -0.1400592313 -0.0181388214 -0.0668643423
#> [1206]  0.0416854068  0.1066086395  0.1152365698 -0.0243174338 -0.0716256256
#> [1211] -0.0871733937 -0.1015930260  0.0303761630 -0.1530431888  0.1537146339
#> [1216] -0.0717469755  0.1957330839 -0.1041891705  0.1658749137  0.0668719595
#> [1221]  0.1265809012  0.1360342098  0.0491911577  0.0517668766  0.0078512159
#> [1226]  0.0340947934  0.1408434773  0.1001679472 -0.0449556343  0.0938500699
#> [1231] -0.0649561294  0.1339960955 -0.1009430353  0.0711330046 -0.0474485522
#> [1236] -0.0014950845  0.1198926997 -0.0411136206  0.1433616889 -0.1201770243
#> [1241] -0.2310277132 -0.0690773061  0.0310130587  0.0204798766  0.0034972684
#> [1246] -0.0523035729 -0.0347435369  0.1006925295  0.1023221461 -0.0370598999
#> [1251] -0.0433365341 -0.2024791766  0.1550705706  0.0380868808 -0.1514007601
#> [1256] -0.0943789854 -0.1568241506  0.0524478269 -0.1469932152  0.0359148962
#> [1261]  0.0565009849  0.0264177040  0.1086578258 -0.0002294490  0.0451486309
#> [1266] -0.1290084011 -0.0248980453 -0.0631485788 -0.0074098673 -0.0954436523
#> [1271] -0.0742516812  0.0452418947 -0.2196254414 -0.0220155797 -0.0673236133
#> [1276]  0.1380232722 -0.0376520442 -0.0820154539  0.0269637765 -0.0461333191
#> [1281] -0.0177440779  0.0999115540 -0.0076986532  0.0680480656  0.0249159041
#> [1286] -0.2279701064  0.1523300697  0.3087517700 -0.3323908724 -0.2524731367
#> [1291] -0.0648602499  0.1127602187  0.0433230735 -0.2018962545 -0.0225194037
#> [1296] -0.0466106488 -0.2042410020 -0.2298454341 -0.1103789118 -0.0343447466
#> [1301] -0.1112793822  0.0272074499 -0.0338956143 -0.0748043552  0.0310523577
#> [1306]  0.0751186541  0.2265988043 -0.0749791788  0.1311049632 -0.0052696357
#> [1311]  0.0644562844 -0.1955630334 -0.1095079654  0.0667778128 -0.1098340625
#> [1316] -0.1201041937 -0.1050953951  0.1056706493 -0.1253357747  0.0268183128
#> [1321]  0.1193037168  0.1208469634  0.0702897436 -0.0243292915  0.1137298692
#> [1326] -0.1699408334 -0.1321356567 -0.0550205235 -0.1289541657  0.1152644903
#> [1331] -0.1318738282 -0.0813950651 -0.0269295156 -0.0973093839  0.1285245113
#> [1336] -0.1497017301  0.0357665278 -0.0292638462 -0.0337592987  0.1586088854
#> [1341] -0.0288719219  0.1355489380  0.0816153598  0.0093463407 -0.2103868224
#> [1346]  0.0221967075  0.1511623582 -0.2232581809 -0.0523630967 -0.1066785106
#> [1351] -0.0431580356  0.0056254103 -0.0560051151  0.0026651413  0.0491828769
#> [1356] -0.2139663938 -0.1384674668 -0.0174376388 -0.0216741591 -0.0806345188
#> [1361]  0.0948394613  0.1469407643 -0.1069034586 -0.0330724617 -0.1101739015
#> [1366] -0.0924908311 -0.0395901210  0.0513445879  0.1557934935 -0.0102545651
#> [1371]  0.2316048273 -0.0556521333 -0.0217150274 -0.0229007688  0.0250864489
#> [1376] -0.1426831005  0.0266652586  0.0680691059 -0.0159307819  0.0297660531
#> [1381] -0.0822440967  0.1246663187  0.1453257986  0.0882157483  0.1431722809
#> [1386] -0.2353767605 -0.1656406600 -0.0347559561 -0.0263796212 -0.0246245027
#> [1391]  0.0032782881  0.0204185235  0.0320959558 -0.0355485764  0.0947640227
#> [1396]  0.0294124154 -0.1365431784  0.0129294345 -0.1424993840  0.0923205480
#> [1401] -0.0802936027  0.0636243116 -0.0709654314  0.0803293734 -0.0313638592
#> [1406] -0.0981195762  0.0328247026  0.1089419650 -0.0530550628 -0.1853496689
#> [1411] -0.0766817211  0.1159235725 -0.0148190022  0.0028681310 -0.0418156293
#> [1416]  0.0031165152  0.0161272865  0.0787602423  0.0916323970 -0.0492610057
#> [1421]  0.0643847831 -0.0785492697  0.0404009245 -0.0085367403 -0.1943671993
#> [1426] -0.2004611993 -0.0083590636  0.1928613426 -0.0506038984 -0.0351104157
#> [1431]  0.0479443044 -0.0622781860 -0.0820589388 -0.0923265802 -0.0610045999
#> [1436]  0.1441252954 -0.2369907912  0.1157772298 -0.0414545842  0.1348697067
#> [1441] -0.0819277656 -0.0254224065  0.0445193420  0.0405920659 -0.0883213349
#> [1446] -0.1893754714 -0.0790917271  0.1536973140 -0.0878702716  0.1681714327
#> [1451]  0.0178538137  0.0861065227  0.1198678908  0.1557052644 -0.1639043849
#> [1456] -0.0242474875 -0.1627520593 -0.3018695947 -0.0181167561  0.1494515046
#> [1461] -0.0764608493 -0.0763797082  0.0567632412  0.0036138020  0.1542684832
#> [1466]  0.0496367497  0.0083603797  0.1164083040  0.0526268470  0.0447170925
#> [1471] -0.1320584850 -0.0753414537 -0.1722153542 -0.0228869271 -0.1504366147
#> [1476]  0.0316413779 -0.0752293516 -0.0877567201 -0.1076709265 -0.1909182751
#> [1481] -0.0719151155  0.2078898038  0.0254338430  0.0921749602  0.0758949794
#> [1486] -0.1130500902  0.0037345884  0.1201607801 -0.1956354724 -0.2190400459
#> [1491] -0.3115359002 -0.0930855152 -0.1186704615 -0.0021584951  0.0523179943
#> [1496] -0.1199089453  0.0152773578 -0.0504036792 -0.0787323650 -0.0609076742
#> [1501] -0.1005411487 -0.0895481378 -0.0894846676 -0.0859864786 -0.1196387684
#> [1506] -0.0208937446  0.1302532202  0.0458832901 -0.1418251196 -0.1877085636
#> [1511]  0.1592330665 -0.0427957241  0.1575738712 -0.0941304269 -0.0841225432
#> [1516]  0.2376486477  0.0006663145  0.1272858769  0.1302497168  0.0029898867
#> [1521] -0.0732263516  0.0250336347 -0.0827185637  0.0580504262 -0.1001060001
#> [1526] -0.2157491264 -0.0739377670 -0.1694784412 -0.0508503713 -0.0124820006
#> [1531]  0.0881439157  0.2709990694 -0.0845968548  0.0163245158  0.2714026516
#> [1536] -0.0122158455  0.1915172185 -0.0697620876 -0.1235006469  0.0749199553
#> [1541] -0.0865028367 -0.1512975620  0.0468638433  0.0368246967 -0.0603242078
#> [1546]  0.0491065531  0.0298451333  0.0577892814 -0.0170836770 -0.0349098325
#> [1551]  0.0466609441 -0.0954690650 -0.1123430222  0.1318508415 -0.0480710314
#> [1556] -0.0198496932 -0.1039156659  0.1031953285  0.0067201336  0.0734391683
#> [1561]  0.0391837794 -0.0722158702 -0.0898397445  0.0345278311 -0.1039121344
#> [1566]  0.2501298261 -0.0762865041  0.0432681145 -0.0290787490 -0.0880894936
#> [1571] -0.0247219910  0.0809793158 -0.0255785942 -0.0141856996 -0.1448321545
#> [1576]  0.1304467196  0.0333274258 -0.2053316504 -0.1102545447 -0.0208834986
#> [1581]  0.0439784435  0.0242339809 -0.0733169851  0.0776717302  0.0069206659
#> [1586]  0.1556153535  0.0244801056 -0.0879521833  0.1453602995  0.0539325981
#> [1591] -0.0401997871 -0.1738563823  0.0339215141 -0.1069670175 -0.1235426062
#> [1596]  0.0056510877 -0.2242934870  0.2107306519  0.0672523280 -0.0831967180
#> [1601] -0.1891027173 -0.0045804666 -0.0914752114 -0.0022885808 -0.5172549606
#> [1606]  0.0294512021  0.0466389157 -0.0565387172 -0.2302407111 -0.1359170277
#> [1611]  0.1116862475  0.1273126481 -0.1363353217  0.1182062542 -0.1184662729
#> [1616]  0.1013189194 -0.2334893381  0.0815569117  0.3838064601  0.1437699675
#> [1621]  0.1749371658  0.0375242325  0.1396129016  0.1462754117  0.1562081224
#> [1626] -0.1043071565  0.1540443615 -0.0491901504  0.0633639790  0.0353268637
#> [1631] -0.0837539546  0.0001848295 -0.0524027577  0.0543533019 -0.0574048648
#> [1636]  0.0072814150  0.1993830605 -0.0384589010 -0.1529172199 -0.1264728769
#> [1641]  0.0144187982  0.1307431950  0.0960653547  0.0811855409 -0.0005807319
#> [1646] -0.0256230516 -0.0344147195 -0.0574155541 -0.0556408419  0.0406526262
#> [1651] -0.1924532482  0.0616756073 -0.1013541924  0.1805566435 -0.1542095962
#> [1656] -0.0994673176  0.0738854187 -0.0484160318  0.0783896533  0.1560250257
#> [1661]  0.1815698109 -0.0999384149 -0.1652291270 -0.1551960288  0.1425457438
#> [1666] -0.0208838912 -0.0722098126 -0.0577741235 -0.1031651280  0.2080115068
#> [1671] -0.1914858218 -0.1192984758 -0.0148553831  0.2236837068  0.0329384456
#> [1676] -0.0696065470 -0.1484405816  0.0408622714 -0.0257485183  0.1972163134
#> [1681] -0.0271142757 -0.1098615244  0.0937363410 -0.0851777569  0.1088879887
#> [1686] -0.1165592540  0.0408364210  0.1078157574 -0.2935129713  0.0047183839
#> [1691] -0.1262672353 -0.1359823308 -0.0193600051 -0.2016408674 -0.1117839463
#> [1696]  0.1564984014  0.0341612196  0.1489474747  0.0065427140 -0.1332869475
#> [1701]  0.0768012283 -0.1013932096  0.0783504406 -0.1951546110 -0.0068798499
#> [1706]  0.0031230166  0.1335935008  0.0933705565 -0.0367031075  0.0188631289
#> [1711] -0.0645198737 -0.1203220770 -0.0914550572 -0.0322719568 -0.0469828285
#> [1716]  0.0234445758 -0.1672696392 -0.2081519108  0.0400471120  0.0103043381
#> [1721]  0.0260477370  0.2092260423  0.0555802936 -0.0348122385  0.1302295906
#> [1726] -0.0508584713  0.2050216784 -0.0719176344 -0.0393559805 -0.0437215864
#> [1731] -0.1046623998  0.1103444410 -0.1065962286  0.0686602597 -0.0082214246
#> [1736]  0.1877898182 -0.0264157735 -0.1829954201 -0.0817471472  0.0984270467
#> [1741] -0.1446059152 -0.0807264210 -0.2034660392  0.0348294809  0.0806340938
#> [1746]  0.3430768425  0.1897298208 -0.3488394459  0.1395189664 -0.0180045185
#> [1751] -0.1456611006  0.0630987479 -0.0675581443  0.0842127387  0.1451124035
#> [1756]  0.0773730070 -0.0817421149 -0.0732851838 -0.0347342701 -0.1667697198
#> [1761] -0.0271526498  0.2285158307 -0.0774094986  0.2627779771 -0.0472285576
#> [1766]  0.1925336801 -0.0144395508 -0.0402155198 -0.0705808674 -0.1914784545
#> 
#> $mod2_unimodal$sim_properties$feature$prop_outliers_feature
#> [1] 0
#> 
#> 
#> 
#> $mod2_unimodal$batch_effects
#> $mod2_unimodal$batch_effects$batch_silhouette
#> [1] 0.0008780031
#> 
#> $mod2_unimodal$batch_effects$pcr_r2
#> [1] 0.01375709
#> 
#> $mod2_unimodal$batch_effects$cms
#> [1] 0
#> 
#> $mod2_unimodal$batch_effects$isi
#> [1] 1.931492
#> 
#> $mod2_unimodal$batch_effects$cross_batch_F1
#> [1] 0.73467
#> 
#> 
#> $mod2_unimodal$biological_signal
#> $mod2_unimodal$biological_signal$da_summary_table
#>                             SubCategory          Metric      Value
#> 1 scATAC-seq Differential Accessibility pval_dist_score 0.85633794
#> 2 scATAC-seq Differential Accessibility        accuracy 0.81666667
#> 3 scATAC-seq Differential Accessibility       precision 0.50000000
#> 4 scATAC-seq Differential Accessibility          recall 0.09090909
#> 5 scATAC-seq Differential Accessibility        f1_score 0.15384615
#> 6 scATAC-seq Differential Accessibility             auc 0.83673469
#>                                                         Direction
#> 1    Higher is better (1.0 = uniform null p-values, no inflation)
#> 2             Higher is better (1.0 = perfect DAR classification)
#> 3                 Higher is better (1.0 = zero false discoveries)
#> 4                 Higher is better (1.0 = 100% true DAR recovery)
#> 5 Higher is better (1.0 = harmonic balance of precision & recall)
#> 6         Higher is better (1.0 = perfect ranking discrimination)
#>                                                      Framework
#> 1 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 2 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 3 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 4 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 5 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 6 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 
#> $mod2_unimodal$biological_signal$pairwise_results
#>             Pair                           SubCategory N_total N_true_DAR
#> 1 TypeA_vs_TypeB scATAC-seq Differential Accessibility      60         11
#>   N_sim_DAR TP FP FN TN pval_dist_score  accuracy precision     recall
#> 1         2  1  1 10 48       0.8563379 0.8166667       0.5 0.09090909
#>    f1_score       auc LogFC_Pearson_r LogFC_Spearman_rho Top_DAR_Jaccard
#> 1 0.1538462 0.8367347       0.5558168          0.4911988               1
#> 
#> $mod2_unimodal$biological_signal$pseudobulks_ref
#> $mod2_unimodal$biological_signal$pseudobulks_ref$TypeA
#>         Batch1 Batch2
#> Peak_01      3      5
#> Peak_02      4      6
#> Peak_03      5      7
#> Peak_04      4      6
#> Peak_05      7      9
#> Peak_06      4      8
#> Peak_07      6     10
#> Peak_08      9      7
#> Peak_09      7      5
#> Peak_10     10      8
#> Peak_11     10      6
#> Peak_12      7     12
#> Peak_13      6      6
#> Peak_14     10      4
#> Peak_15      8      7
#> Peak_16      7      7
#> Peak_17      9      8
#> Peak_18      6      3
#> Peak_19      8     11
#> Peak_20     10      9
#> Peak_21      7      6
#> Peak_22      8     11
#> Peak_23      7      7
#> Peak_24     10      6
#> Peak_25      3      3
#> Peak_26      6      9
#> Peak_27      8      5
#> Peak_28      6      5
#> Peak_29     13     10
#> Peak_30      5      4
#> Peak_31     11      9
#> Peak_32     10      6
#> Peak_33      7      8
#> Peak_34      7      8
#> Peak_35      8     12
#> Peak_36      9      8
#> Peak_37      7      4
#> Peak_38      7      5
#> Peak_39      6      7
#> Peak_40     11      9
#> Peak_41      5      5
#> Peak_42      9     11
#> Peak_43      6      7
#> Peak_44     10      6
#> Peak_45     11     10
#> Peak_46      5      5
#> Peak_47      7      9
#> Peak_48      7     10
#> Peak_49      3      8
#> Peak_50      7      5
#> Peak_51      5      3
#> Peak_52      5      8
#> Peak_53      4      8
#> Peak_54     10      4
#> Peak_55     11      7
#> Peak_56     10      7
#> Peak_57     10     14
#> Peak_58     11      8
#> Peak_59      7     10
#> Peak_60      7      7
#> 
#> $mod2_unimodal$biological_signal$pseudobulks_ref$TypeB
#>         Batch1 Batch2
#> Peak_01     11     11
#> Peak_02     12     16
#> Peak_03     13     12
#> Peak_04     17     10
#> Peak_05     15     18
#> Peak_06     13     13
#> Peak_07     12     12
#> Peak_08     13     14
#> Peak_09     15     16
#> Peak_10     17     17
#> Peak_11     17     14
#> Peak_12     14     18
#> Peak_13     15     14
#> Peak_14     13     11
#> Peak_15     16     15
#> Peak_16      7      5
#> Peak_17      6      9
#> Peak_18      5      3
#> Peak_19      7      6
#> Peak_20      6      7
#> Peak_21      5      5
#> Peak_22      8      7
#> Peak_23      6      7
#> Peak_24     10      9
#> Peak_25      4      6
#> Peak_26      6      9
#> Peak_27      8      8
#> Peak_28     11     10
#> Peak_29      9      9
#> Peak_30      5      9
#> Peak_31      8      5
#> Peak_32      6      9
#> Peak_33      6      4
#> Peak_34     11      9
#> Peak_35      9      9
#> Peak_36      7     12
#> Peak_37      5      7
#> Peak_38      6      5
#> Peak_39      8      7
#> Peak_40      7      5
#> Peak_41      1      7
#> Peak_42     13     11
#> Peak_43     12      8
#> Peak_44      7      8
#> Peak_45      9     11
#> Peak_46      9      7
#> Peak_47      4      6
#> Peak_48     10     13
#> Peak_49      2      2
#> Peak_50      8      7
#> Peak_51      5      3
#> Peak_52      3      6
#> Peak_53      5      4
#> Peak_54      8     11
#> Peak_55     10     12
#> Peak_56     12      6
#> Peak_57     12      8
#> Peak_58      9     10
#> Peak_59      8     10
#> Peak_60      7      4
#> 
#> 
#> $mod2_unimodal$biological_signal$pseudobulks_sim
#> $mod2_unimodal$biological_signal$pseudobulks_sim$TypeA
#>         Batch1 Batch2
#> Peak_01      6      4
#> Peak_02      7      3
#> Peak_03      4      3
#> Peak_04      3      6
#> Peak_05     10      4
#> Peak_06     11      4
#> Peak_07      7      1
#> Peak_08      8      6
#> Peak_09      3      4
#> Peak_10      6      9
#> Peak_11     12      8
#> Peak_12      7     13
#> Peak_13      6      9
#> Peak_14      6     10
#> Peak_15     12      9
#> Peak_16      7      4
#> Peak_17      9      8
#> Peak_18      2      5
#> Peak_19      5      6
#> Peak_20      2      6
#> Peak_21      3      4
#> Peak_22      8     10
#> Peak_23      5      5
#> Peak_24     11      9
#> Peak_25      5     10
#> Peak_26      7      5
#> Peak_27     10      6
#> Peak_28      6      6
#> Peak_29      6     10
#> Peak_30      2      4
#> Peak_31      5      6
#> Peak_32     10      4
#> Peak_33      5      4
#> Peak_34      6      6
#> Peak_35      9      6
#> Peak_36      6      4
#> Peak_37      5      9
#> Peak_38      5      5
#> Peak_39      2      7
#> Peak_40      8      9
#> Peak_41      2      6
#> Peak_42     10     13
#> Peak_43      8     10
#> Peak_44     10      5
#> Peak_45     10     10
#> Peak_46      6      4
#> Peak_47      9      6
#> Peak_48     10     11
#> Peak_49      5      4
#> Peak_50     12      9
#> Peak_51      5      7
#> Peak_52      4      3
#> Peak_53      7      8
#> Peak_54      9      5
#> Peak_55      5     13
#> Peak_56      8      7
#> Peak_57     16      6
#> Peak_58      7      6
#> Peak_59     10      7
#> Peak_60      7      5
#> 
#> $mod2_unimodal$biological_signal$pseudobulks_sim$TypeB
#>         Batch1 Batch2
#> Peak_01     14     14
#> Peak_02      9     15
#> Peak_03     14     14
#> Peak_04     11     14
#> Peak_05     14     16
#> Peak_06     12     16
#> Peak_07     16     13
#> Peak_08     15     14
#> Peak_09     11      7
#> Peak_10     16     15
#> Peak_11     17     15
#> Peak_12     13     17
#> Peak_13     13     16
#> Peak_14     16     14
#> Peak_15     14     12
#> Peak_16      9     10
#> Peak_17      8      5
#> Peak_18      2      3
#> Peak_19      5      6
#> Peak_20      7      4
#> Peak_21      8      3
#> Peak_22      6      7
#> Peak_23      6      6
#> Peak_24     10      9
#> Peak_25      8      5
#> Peak_26      7      7
#> Peak_27      5      9
#> Peak_28      9      7
#> Peak_29      9     11
#> Peak_30      4      4
#> Peak_31      7      6
#> Peak_32      6      7
#> Peak_33      5      7
#> Peak_34     11      7
#> Peak_35      9      5
#> Peak_36      8      9
#> Peak_37      4      8
#> Peak_38      5      4
#> Peak_39      4      5
#> Peak_40      7      9
#> Peak_41      4      5
#> Peak_42      6      9
#> Peak_43     11      6
#> Peak_44      9      6
#> Peak_45      8      8
#> Peak_46      5      7
#> Peak_47      8      4
#> Peak_48      6     11
#> Peak_49      7      5
#> Peak_50      6      6
#> Peak_51      6      3
#> Peak_52      5      6
#> Peak_53      3      5
#> Peak_54      5     11
#> Peak_55     16     12
#> Peak_56      7      7
#> Peak_57      8     10
#> Peak_58      7     11
#> Peak_59      7      3
#> Peak_60      7      7
#> 
#> 
#> $mod2_unimodal$biological_signal$da_results_ref
#> $mod2_unimodal$biological_signal$da_results_ref$TypeA_vs_TypeB
#>         feature       logFC   logCPM       F_stat       PValue         FDR
#> Peak_01 Peak_01  1.28531261 14.22614 11.540251884 0.0009244519 0.011093423
#> Peak_02 Peak_02  1.31598148 14.50177 15.402310051 0.0001453586 0.004568652
#> Peak_03 Peak_03  0.89663746 14.47007  7.365534953 0.0076281536 0.041608111
#> Peak_04 Peak_04  1.26401102 14.47006 13.885816108 0.0002976166 0.005952331
#> Peak_05 Peak_05  0.88549848 14.81110 10.444579056 0.0015881324 0.015881324
#> Peak_06 Peak_06  0.95264340 14.50178  8.514421397 0.0042077395 0.028051596
#> Peak_07 Peak_07  0.43033667 14.56319  1.949167601 0.1652540108 0.367231135
#> Peak_08 Peak_08  0.59852365 14.65064  4.077300889 0.0456913356 0.152304452
#> Peak_09 Peak_09  1.20406871 14.65062 15.302990463 0.0001522884 0.004568652
#> Peak_10 Peak_10  0.76064935 14.88510  8.607271703 0.0040123221 0.028051596
#> Peak_11 Peak_11  0.79603617 14.75956  8.006089246 0.0054664471 0.032798682
#> Peak_12 Peak_12  0.59686319 14.86087  5.215576039 0.0241439580 0.103474106
#> Peak_13 Peak_13  1.10868988 14.59292 12.378623328 0.0006142302 0.009213453
#> Peak_14 Peak_14  0.61997547 14.50177  3.743741204 0.0553589620 0.158168463
#> Peak_15 Peak_15  0.88774186 14.73309  9.545246130 0.0024917097 0.021357512
#> Peak_16 Peak_16 -0.36447582 14.06570  0.847826649 0.3590152736 0.566866221
#> Peak_17 Peak_17 -0.32460904 14.30017  0.880196597 0.3500323240 0.566866221
#> Peak_18 Peak_18 -0.30991147 13.62208  0.308556118 0.5796016620 0.739917015
#> Peak_19 Peak_19 -0.68561609 14.30018  3.875542327 0.0512998099 0.153899430
#> Peak_20 Peak_20 -0.68561609 14.30017  3.876416004 0.0512739925 0.153899430
#> Peak_21 Peak_21 -0.51626906 13.93243  1.375398233 0.2432093521 0.461824530
#> Peak_22 Peak_22 -0.48298338 14.37057  2.070471061 0.1527783764 0.364775027
#> Peak_23 Peak_23 -0.25137983 14.10752  0.433535058 0.5115207278 0.682027637
#> Peak_24 Peak_24  0.09740973 14.40449  0.086984027 0.7685565786 0.853951754
#> Peak_25 Peak_25  0.56937183 13.56319  0.908256599 0.3424931616 0.566866221
#> Peak_26 Peak_26 -0.14668022 14.22617  0.170279736 0.6805996416 0.787464022
#> Peak_27 Peak_27  0.14787656 14.18768  0.168176410 0.6824688189 0.787464022
#> Peak_28 Peak_28  0.77094771 14.30014  4.704472039 0.0320565347 0.120212005
#> Peak_29 Peak_29 -0.49663556 14.59295  2.726646077 0.1013010009 0.264263481
#> Peak_30 Peak_30  0.47763510 13.93242  1.152464806 0.2851877882 0.518523251
#> Peak_31 Peak_31 -0.75872735 14.33580  4.878430508 0.0290937213 0.116374885
#> Peak_32 Peak_32 -0.23838590 14.26363  0.461216281 0.4983638845 0.679587115
#> Peak_33 Peak_33 -0.71927045 14.02264  3.034056376 0.0840963051 0.229353559
#> Peak_34 Peak_34  0.26213113 14.40449  0.625291992 0.4306475528 0.630215931
#> Peak_35 Peak_35 -0.29720802 14.50181  0.894639487 0.3461238190 0.566866221
#> Peak_36 Peak_36  0.01118660 14.43766  0.001187637 0.9725659140 0.972565914
#> Peak_37 Peak_37 -0.02309002 13.93242  0.002761925 0.9581746416 0.972565914
#> Peak_38 Peak_38 -0.26880382 13.93243  0.374948045 0.5414780071 0.706275661
#> Peak_39 Peak_39  0.05632117 14.14816  0.023188601 0.8792233825 0.909541430
#> Peak_40 Peak_40 -0.87182334 14.30017  6.189300867 0.0142251628 0.071125814
#> Peak_41 Peak_41 -0.45820680 13.67868  0.738617126 0.3918178407 0.602796678
#> Peak_42 Peak_42  0.11262793 14.67867  0.153694146 0.6957251917 0.787613425
#> Peak_43 Peak_43  0.46513252 14.33578  1.819333952 0.1799324714 0.385569582
#> Peak_44 Peak_44 -0.23838590 14.26363  0.461281312 0.4983336603 0.679587115
#> Peak_45 Peak_45 -0.21686753 14.59295  0.522636240 0.4711259920 0.673037131
#> Peak_46 Peak_46  0.51865685 14.06569  1.668316503 0.1989664410 0.411654705
#> Peak_47 Peak_47 -0.81097614 14.06571  4.112556493 0.0447804460 0.152304452
#> Peak_48 Peak_48  0.28339252 14.56320  0.853505208 0.3574165881 0.566866221
#> Peak_49 Peak_49 -1.54748112 13.50182  5.602992618 0.0195298585 0.090137808
#> Peak_50 Peak_50  0.16964664 14.10751  0.195889915 0.6588557942 0.787464022
#> Peak_51 Peak_51 -0.14460782 13.56319  0.060624347 0.8059325121 0.879199104
#> Peak_52 Peak_52 -0.66417154 13.88514  2.108434046 0.1490987469 0.364775027
#> Peak_53 Peak_53 -0.55084607 13.83623  1.357359258 0.2463064160 0.461824530
#> Peak_54 Peak_54  0.28704852 14.33577  0.701030106 0.4041021870 0.606153281
#> Peak_55 Peak_55  0.13864382 14.56318  0.205762503 0.6509284399 0.787464022
#> Peak_56 Peak_56 -0.06572241 14.40449  0.039746560 0.8423143530 0.886646687
#> Peak_57 Peak_57 -0.40749793 14.67869  2.017684112 0.1580691785 0.364775027
#> Peak_58 Peak_58 -0.14718774 14.50179  0.219570926 0.6402171365 0.787464022
#> Peak_59 Peak_59 -0.06572241 14.40451  0.039756379 0.8422951343 0.886646687
#> Peak_60 Peak_60 -0.48719757 14.02263  1.407326552 0.2378441759 0.461824530
#>         is_da
#> Peak_01  TRUE
#> Peak_02  TRUE
#> Peak_03  TRUE
#> Peak_04  TRUE
#> Peak_05  TRUE
#> Peak_06  TRUE
#> Peak_07 FALSE
#> Peak_08 FALSE
#> Peak_09  TRUE
#> Peak_10  TRUE
#> Peak_11  TRUE
#> Peak_12 FALSE
#> Peak_13  TRUE
#> Peak_14 FALSE
#> Peak_15  TRUE
#> Peak_16 FALSE
#> Peak_17 FALSE
#> Peak_18 FALSE
#> Peak_19 FALSE
#> Peak_20 FALSE
#> Peak_21 FALSE
#> Peak_22 FALSE
#> Peak_23 FALSE
#> Peak_24 FALSE
#> Peak_25 FALSE
#> Peak_26 FALSE
#> Peak_27 FALSE
#> Peak_28 FALSE
#> Peak_29 FALSE
#> Peak_30 FALSE
#> Peak_31 FALSE
#> Peak_32 FALSE
#> Peak_33 FALSE
#> Peak_34 FALSE
#> Peak_35 FALSE
#> Peak_36 FALSE
#> Peak_37 FALSE
#> Peak_38 FALSE
#> Peak_39 FALSE
#> Peak_40 FALSE
#> Peak_41 FALSE
#> Peak_42 FALSE
#> Peak_43 FALSE
#> Peak_44 FALSE
#> Peak_45 FALSE
#> Peak_46 FALSE
#> Peak_47 FALSE
#> Peak_48 FALSE
#> Peak_49 FALSE
#> Peak_50 FALSE
#> Peak_51 FALSE
#> Peak_52 FALSE
#> Peak_53 FALSE
#> Peak_54 FALSE
#> Peak_55 FALSE
#> Peak_56 FALSE
#> Peak_57 FALSE
#> Peak_58 FALSE
#> Peak_59 FALSE
#> Peak_60 FALSE
#> 
#> 
#> $mod2_unimodal$biological_signal$da_results_sim
#> $mod2_unimodal$biological_signal$da_results_sim$TypeA_vs_TypeB
#>         feature       logFC   logCPM       F_stat       PValue        FDR is_da
#> Peak_01 Peak_01  1.25446850 14.59948 7.586161e+00 0.0067979617 0.10691384 FALSE
#> Peak_02 Peak_02  1.03435648 14.46824 5.376410e+00 0.0221029709 0.16577228 FALSE
#> Peak_03 Peak_03  1.75510805 14.50218 1.359780e+01 0.0003414759 0.02048856  TRUE
#> Peak_04 Peak_04  1.24098975 14.46824 7.495383e+00 0.0071275892 0.10691384 FALSE
#> Peak_05 Peak_05  0.87702858 14.77637 3.872594e+00 0.0513870476 0.30506181 FALSE
#> Peak_06 Peak_06  0.68043706 14.74836 2.392691e+00 0.1245368453 0.43954181 FALSE
#> Peak_07 Peak_07  1.61840612 14.56777 1.183112e+01 0.0008018078 0.02405423  TRUE
#> Peak_08 Peak_08  0.82855990 14.74836 3.487400e+00 0.0642774646 0.30506181 FALSE
#> Peak_09 Peak_09  1.12526552 14.12034 5.865343e+00 0.0169379568 0.16577228 FALSE
#> Peak_10 Peak_10  0.82595269 14.83082 3.439836e+00 0.0660967264 0.30506181 FALSE
#> Peak_11 Peak_11  0.46182873 14.98283 1.097282e+00 0.2969686806 0.66541089 FALSE
#> Peak_12 Peak_12  0.36951881 14.93392 7.124508e-01 0.4003109531 0.71256709 FALSE
#> Peak_13 Peak_13  0.73059088 14.77637 2.737091e+00 0.1006564795 0.43138491 FALSE
#> Peak_14 Peak_14  0.68732239 14.83082 2.421468e+00 0.1223162273 0.43954181 FALSE
#> Peak_15 Peak_15  0.09543170 14.85732 4.877506e-02 0.8255839592 0.93462335 FALSE
#> Peak_16 Peak_16  0.56698231 14.32387 1.720615e+00 0.1921178489 0.57943439 FALSE
#> Peak_17 Peak_17 -0.58875225 14.32391 1.924845e+00 0.1678947611 0.55964920 FALSE
#> Peak_18 Peak_18 -0.66866607 13.39790 1.086087e+00 0.2994348999 0.66541089 FALSE
#> Peak_19 Peak_19 -0.20697673 13.98286 1.952197e-01 0.6594025811 0.82902037 FALSE
#> Peak_20 Peak_20  0.24133761 13.83085 2.263297e-01 0.6351236936 0.82902037 FALSE
#> Peak_21 Peak_21  0.42819336 13.77640 6.671542e-01 0.4156641353 0.71256709 FALSE
#> Peak_22 Peak_22 -0.67013634 14.36138 2.484602e+00 0.1175980884 0.43954181 FALSE
#> Peak_23 Peak_23  0.05021920 13.98285 1.142159e-02 0.9150691344 0.98063157 FALSE
#> Peak_24 Peak_24 -0.28205934 14.63056 4.379152e-01 0.5093991605 0.72771309 FALSE
#> Peak_25 Peak_25 -0.41076410 14.24590 9.363288e-01 0.3351699576 0.67244446 FALSE
#> Peak_26 Peak_26  0.01051796 14.16342 5.918404e-04 0.9806315710 0.98063157 FALSE
#> Peak_27 Peak_27 -0.39767143 14.32390 8.822018e-01 0.3494860786 0.67642467 FALSE
#> Peak_28 Peak_28  0.19976792 14.24588 2.187320e-01 0.6408560372 0.82902037 FALSE
#> Peak_29 Peak_29  0.10873999 14.53539 6.484054e-02 0.7994385242 0.92242907 FALSE
#> Peak_30 Peak_30  0.19629573 13.53539 1.102660e-01 0.7404208351 0.87108334 FALSE
#> Peak_31 Peak_31  0.02873680 14.07596 4.090317e-03 0.9491120090 0.98063157 FALSE
#> Peak_32 Peak_32 -0.31279510 14.20525 5.388756e-01 0.4643310019 0.72596357 FALSE
#> Peak_33 Peak_33  0.19859693 13.93394 1.693535e-01 0.6814210812 0.82902037 FALSE
#> Peak_34 Peak_34  0.36704674 14.32388 7.337291e-01 0.3933844286 0.71256709 FALSE
#> Peak_35 Peak_35 -0.30593424 14.28543 5.219352e-01 0.4714229369 0.72596357 FALSE
#> Peak_36 Peak_36  0.54340967 14.20523 1.563302e+00 0.2136150042 0.57943439 FALSE
#> Peak_37 Peak_37 -0.42584323 14.16343 9.735706e-01 0.3257766046 0.67244446 FALSE
#> Peak_38 Peak_38 -0.35433923 13.83086 4.959861e-01 0.4826336156 0.72596357 FALSE
#> Peak_39 Peak_39 -0.20596150 13.77640 1.589319e-01 0.6908503055 0.82902037 FALSE
#> Peak_40 Peak_40 -0.29467244 14.43352 4.832431e-01 0.4883008949 0.72596357 FALSE
#> Peak_41 Peak_41 -0.04055543 13.71982 5.807827e-03 0.9393796537 0.98063157 FALSE
#> Peak_42 Peak_42 -0.81722631 14.59955 3.621140e+00 0.0594457272 0.30506181 FALSE
#> Peak_43 Peak_43 -0.28999399 14.50224 4.661599e-01 0.4960751045 0.72596357 FALSE
#> Peak_44 Peak_44 -0.20820830 14.32390 2.419183e-01 0.6237226056 0.82902037 FALSE
#> Peak_45 Peak_45 -0.52622504 14.53541 1.527236e+00 0.2189414068 0.57943439 FALSE
#> Peak_46 Peak_46  0.05021920 13.98285 1.142027e-02 0.9150740112 0.98063157 FALSE
#> Peak_47 Peak_47 -0.52381222 14.20526 1.506236e+00 0.2221165162 0.57943439 FALSE
#> Peak_48 Peak_48 -0.50976488 14.59954 1.428429e+00 0.2343775912 0.58594398 FALSE
#> Peak_49 Peak_49  0.19859694 13.93394 1.693463e-01 0.6814275035 0.82902037 FALSE
#> Peak_50 Peak_50 -1.00295534 14.43354 5.430891e+00 0.0214533938 0.16577228 FALSE
#> Peak_51 Peak_51 -0.61191845 13.93396 1.619299e+00 0.2056496425 0.57943439 FALSE
#> Peak_52 Peak_52  0.42819335 13.77639 6.674835e-01 0.4155494049 0.71256709 FALSE
#> Peak_53 Peak_53 -1.09308454 14.03018 5.485528e+00 0.0208220688 0.16577228 FALSE
#> Peak_54 Peak_54 -0.01871528 14.32388 1.946987e-03 0.9648783627 0.98063157 FALSE
#> Peak_55 Peak_55  0.42106484 14.83084 9.322443e-01 0.3362222318 0.67244446 FALSE
#> Peak_56 Peak_56 -0.30593424 14.28542 5.221181e-01 0.4713454294 0.72596357 FALSE
#> Peak_57 Peak_57 -0.49494699 14.66094 1.341463e+00 0.2490760021 0.59778241 FALSE
#> Peak_58 Peak_58  0.25367316 14.36135 3.530982e-01 0.5534831479 0.77230207 FALSE
#> Peak_59 Peak_59 -0.95853821 14.20527 4.927859e+00 0.0283056924 0.18870462 FALSE
#> Peak_60 Peak_60  0.01051796 14.16342 5.918404e-04 0.9806315710 0.98063157 FALSE
#> 
#> 
#> 
#> $mod2_unimodal$trajectory
#> NULL
#> 
#> $mod2_unimodal$resource_usage
#> $mod2_unimodal$resource_usage$core_seconds
#> NULL
#> 
#> $mod2_unimodal$resource_usage$memory_per_core_mib
#> NULL
#> 
#> $mod2_unimodal$resource_usage$memory_per_core_gib
#> NULL
#> 
#> $mod2_unimodal$resource_usage$memory_mb
#> NULL
#> 
#> $mod2_unimodal$resource_usage$memory_input
#> NULL
#> 
#> $mod2_unimodal$resource_usage$memory_unit
#> [1] "MiB"
#> 
#> $mod2_unimodal$resource_usage$elapsed_time
#> NULL
#> 
#> $mod2_unimodal$resource_usage$elapsed_time_input
#> NULL
#> 
#> $mod2_unimodal$resource_usage$time_unit
#> [1] "seconds"
#> 
#> $mod2_unimodal$resource_usage$sim_cores
#> NULL
#> 
#> 
#> 
#> $clustering
#> $clustering$m1
#> $clustering$m1$silhouette_sim
#> [1] 0.003402312
#> 
#> $clustering$m1$silhouette
#> [1] 0.003402312
#> 
#> $clustering$m1$davies_bouldin
#> [1] 7.916816
#> 
#> $clustering$m1$davies_bouldin_sim
#> [1] 7.916816
#> 
#> $clustering$m1$neighborhood_purity
#> [1] 0.4775
#> 
#> $clustering$m1$silhouette_ref
#> [1] 0.001852738
#> 
#> $clustering$m1$davies_bouldin_ref
#> [1] 8.436759
#> 
#> $clustering$m1$cdi
#> [1] 4.984432
#> 
#> $clustering$m1$CDI_AIC
#> [1] 4.984432
#> 
#> $clustering$m1$CDI_BIC
#> [1] 5.104526
#> 
#> 
#> $clustering$m2
#> $clustering$m2$silhouette_sim
#> [1] 0.04149183
#> 
#> $clustering$m2$silhouette
#> [1] 0.04149183
#> 
#> $clustering$m2$davies_bouldin
#> [1] 4.163203
#> 
#> $clustering$m2$davies_bouldin_sim
#> [1] 4.163203
#> 
#> $clustering$m2$neighborhood_purity
#> [1] 0.7275
#> 
#> $clustering$m2$silhouette_ref
#> [1] 0.04023238
#> 
#> $clustering$m2$davies_bouldin_ref
#> [1] 4.208554
#> 
#> $clustering$m2$cdi
#> [1] 1.687605
#> 
#> $clustering$m2$CDI_AIC
#> [1] 1.687605
#> 
#> $clustering$m2$CDI_BIC
#> [1] 1.806707
#> 
#> 
#> 
#> $batch
#> $batch$rna
#> $batch$rna$batch_silhouette
#> [1] -0.002710457
#> 
#> $batch$rna$pcr_r2
#> [1] 0.009476487
#> 
#> $batch$rna$cms
#> [1] 0
#> 
#> $batch$rna$isi
#> [1] 1.938397
#> 
#> $batch$rna$cross_batch_F1
#> [1] 0.5126418
#> 
#> 
#> $batch$atac
#> $batch$atac$batch_silhouette
#> [1] 0.0008780031
#> 
#> $batch$atac$pcr_r2
#> [1] 0.01375709
#> 
#> $batch$atac$cms
#> [1] 0
#> 
#> $batch$atac$isi
#> [1] 1.931492
#> 
#> $batch$atac$cross_batch_F1
#> [1] 0.73467
#> 
#> 
#> 
#> $signal
#> $signal$rna
#> $signal$rna$deg_summary_table
#>                         SubCategory          Metric     Value
#> 1 scRNA-seq Differential Expression pval_dist_score 0.8719073
#> 2 scRNA-seq Differential Expression        accuracy 0.9500000
#> 3 scRNA-seq Differential Expression       precision 0.0000000
#> 4 scRNA-seq Differential Expression          recall 0.0000000
#> 5 scRNA-seq Differential Expression        f1_score 0.0000000
#> 6 scRNA-seq Differential Expression             auc 0.6465517
#>                                                         Direction
#> 1    Higher is better (1.0 = uniform null p-values, no inflation)
#> 2             Higher is better (1.0 = perfect DEG classification)
#> 3                 Higher is better (1.0 = zero false discoveries)
#> 4                 Higher is better (1.0 = 100% true DEG recovery)
#> 5 Higher is better (1.0 = harmonic balance of precision & recall)
#> 6         Higher is better (1.0 = perfect ranking discrimination)
#>                                            Framework
#> 1 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 2 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 3 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 4 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 5 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 6 Libra Pseudobulk edgeR-QLF (Squair 2021; Duo 2024)
#> 
#> $signal$rna$pairwise_results
#>             Pair                       SubCategory N_total N_true_DEG N_sim_DEG
#> 1 TypeA_vs_TypeB scRNA-seq Differential Expression      60          2         1
#>   TP FP FN TN pval_dist_score accuracy precision recall f1_score       auc
#> 1  0  1  2 57       0.8719073     0.95         0      0        0 0.6465517
#>   LogFC_Pearson_r LogFC_Spearman_rho Top_DEG_Jaccard
#> 1       0.4413862          0.3846624               1
#> 
#> $signal$rna$pseudobulks_ref
#> $signal$rna$pseudobulks_ref$TypeA
#>         Batch1 Batch2
#> Gene_01    164     87
#> Gene_02     65     67
#> Gene_03     75     81
#> Gene_04     72     95
#> Gene_05    112     99
#> Gene_06    110     72
#> Gene_07     97     81
#> Gene_08    104     76
#> Gene_09    156    115
#> Gene_10    105     68
#> Gene_11    102     81
#> Gene_12    123    108
#> Gene_13     45     30
#> Gene_14     97     55
#> Gene_15     72     75
#> Gene_16     67     97
#> Gene_17     75     63
#> Gene_18     62     32
#> Gene_19     56     44
#> Gene_20    112     63
#> Gene_21     72     56
#> Gene_22     75     51
#> Gene_23     94     71
#> Gene_24    132    104
#> Gene_25    142    113
#> Gene_26    103     64
#> Gene_27    110     55
#> Gene_28     43     49
#> Gene_29     89     72
#> Gene_30     96     79
#> Gene_31     98     82
#> Gene_32     73     91
#> Gene_33     92     80
#> Gene_34     90     48
#> Gene_35     83     63
#> Gene_36     41     42
#> Gene_37     58     55
#> Gene_38     93     75
#> Gene_39     34     35
#> Gene_40     87    102
#> Gene_41     76    104
#> Gene_42     46     60
#> Gene_43     95     63
#> Gene_44     63     40
#> Gene_45     59     42
#> Gene_46     63    118
#> Gene_47     62     48
#> Gene_48     86    113
#> Gene_49     58    113
#> Gene_50    124     97
#> Gene_51     96     62
#> Gene_52     51    117
#> Gene_53    127    137
#> Gene_54    124     87
#> Gene_55     58     64
#> Gene_56    120     52
#> Gene_57     81    140
#> Gene_58     64     93
#> Gene_59     45     32
#> Gene_60     85     62
#> 
#> $signal$rna$pseudobulks_ref$TypeB
#>         Batch1 Batch2
#> Gene_01    137    159
#> Gene_02     93    133
#> Gene_03    175    115
#> Gene_04    119    106
#> Gene_05    188    174
#> Gene_06     81    109
#> Gene_07    129    140
#> Gene_08     94    107
#> Gene_09    147    185
#> Gene_10    184     94
#> Gene_11    179    139
#> Gene_12    197     97
#> Gene_13     87     99
#> Gene_14    131    160
#> Gene_15    111    137
#> Gene_16     73     96
#> Gene_17     72     84
#> Gene_18     40     36
#> Gene_19     46     40
#> Gene_20     91    156
#> Gene_21     46     50
#> Gene_22     41     41
#> Gene_23     67     69
#> Gene_24     85     90
#> Gene_25     81    111
#> Gene_26     84     81
#> Gene_27    100     53
#> Gene_28     62     38
#> Gene_29     86     94
#> Gene_30     45     66
#> Gene_31     57     98
#> Gene_32    113     87
#> Gene_33     79     91
#> Gene_34     43     51
#> Gene_35    103     99
#> Gene_36     48     58
#> Gene_37     83     63
#> Gene_38     56     81
#> Gene_39     51     39
#> Gene_40     90     87
#> Gene_41     58    106
#> Gene_42     90     50
#> Gene_43     91     94
#> Gene_44     64     50
#> Gene_45     69     48
#> Gene_46     74     51
#> Gene_47     56     65
#> Gene_48    100     88
#> Gene_49     69     90
#> Gene_50     98    119
#> Gene_51    108    114
#> Gene_52     54     52
#> Gene_53    103     95
#> Gene_54     91    104
#> Gene_55    108     86
#> Gene_56     94     82
#> Gene_57    118     66
#> Gene_58    113    133
#> Gene_59     23     19
#> Gene_60     70     56
#> 
#> 
#> $signal$rna$pseudobulks_sim
#> $signal$rna$pseudobulks_sim$TypeA
#>         Batch1 Batch2
#> Gene_01     82     96
#> Gene_02     67     77
#> Gene_03     81     67
#> Gene_04     62     78
#> Gene_05     53     96
#> Gene_06     86     61
#> Gene_07    108    132
#> Gene_08     68     86
#> Gene_09     61    135
#> Gene_10    122     74
#> Gene_11     70     52
#> Gene_12    109    106
#> Gene_13     46     30
#> Gene_14     57     82
#> Gene_15     57     84
#> Gene_16     89     87
#> Gene_17     64     55
#> Gene_18     35     32
#> Gene_19     44     49
#> Gene_20    106     91
#> Gene_21     47     57
#> Gene_22     50     36
#> Gene_23     79     54
#> Gene_24    107    126
#> Gene_25    114     81
#> Gene_26     48     50
#> Gene_27     62     86
#> Gene_28     48     66
#> Gene_29     79    114
#> Gene_30     57    100
#> Gene_31     75     72
#> Gene_32     76     68
#> Gene_33    120     74
#> Gene_34     70     94
#> Gene_35    112     87
#> Gene_36     46     50
#> Gene_37     60     75
#> Gene_38     54     58
#> Gene_39     38     29
#> Gene_40     63     47
#> Gene_41     80     61
#> Gene_42     53     93
#> Gene_43     89     84
#> Gene_44     90     78
#> Gene_45     54     42
#> Gene_46     84     39
#> Gene_47     53     93
#> Gene_48     86    102
#> Gene_49     60     69
#> Gene_50     80     80
#> Gene_51     65     59
#> Gene_52     88     85
#> Gene_53     97    113
#> Gene_54     77     97
#> Gene_55     57     87
#> Gene_56     67     95
#> Gene_57     67     85
#> Gene_58     67     57
#> Gene_59     36     24
#> Gene_60    103     63
#> 
#> $signal$rna$pseudobulks_sim$TypeB
#>         Batch1 Batch2
#> Gene_01    108    158
#> Gene_02    184    162
#> Gene_03    111    129
#> Gene_04    186    106
#> Gene_05     84    113
#> Gene_06    137    202
#> Gene_07    156    111
#> Gene_08    110    110
#> Gene_09    166    207
#> Gene_10     93    128
#> Gene_11    187    140
#> Gene_12    156    169
#> Gene_13     94     96
#> Gene_14     63     89
#> Gene_15    154    104
#> Gene_16     53     78
#> Gene_17     62     72
#> Gene_18     42     28
#> Gene_19     46     28
#> Gene_20     88    125
#> Gene_21     66     67
#> Gene_22     35     65
#> Gene_23     64     74
#> Gene_24    105     81
#> Gene_25    110     84
#> Gene_26     74     71
#> Gene_27     60     89
#> Gene_28     46     27
#> Gene_29    109     81
#> Gene_30     86     56
#> Gene_31     74     83
#> Gene_32     58    111
#> Gene_33    125    115
#> Gene_34     59     59
#> Gene_35    104     85
#> Gene_36     50     42
#> Gene_37     39     59
#> Gene_38    100     67
#> Gene_39     34     40
#> Gene_40     77     68
#> Gene_41     69     48
#> Gene_42     96     49
#> Gene_43     90     41
#> Gene_44     32     65
#> Gene_45     64     70
#> Gene_46     71     82
#> Gene_47     49     54
#> Gene_48    101    106
#> Gene_49     42     55
#> Gene_50    103    110
#> Gene_51     59     84
#> Gene_52     75     68
#> Gene_53    116     80
#> Gene_54     61    156
#> Gene_55     35    116
#> Gene_56     73     86
#> Gene_57     78    116
#> Gene_58    107     83
#> Gene_59     30     31
#> Gene_60     75    101
#> 
#> 
#> $signal$rna$de_results_ref
#> $signal$rna$de_results_ref$TypeA_vs_TypeB
#>         feature         logFC   logCPM       F_stat      PValue        FDR
#> Gene_01 Gene_01  1.238685e-01 14.70815 1.574821e-01 0.696937479 0.85799160
#> Gene_02 Gene_02  6.239397e-01 14.11652 3.907120e+00 0.066296970 0.28099865
#> Gene_03 Gene_03  7.316681e-01 14.42108 5.266377e+00 0.036196721 0.27517630
#> Gene_04 Gene_04  2.591277e-01 14.25533 6.691704e-01 0.425833651 0.68681783
#> Gene_05 Gene_05  6.309193e-01 14.77632 4.718310e+00 0.045862717 0.27517630
#> Gene_06 Gene_06 -6.371856e-02 14.17162 4.291334e-02 0.838606347 0.93178483
#> Gene_07 Gene_07  4.531216e-01 14.42722 2.387011e+00 0.142672602 0.40763601
#> Gene_08 Gene_08  2.531663e-02 14.20603 7.375600e-03 0.932668838 0.97218763
#> Gene_09 Gene_09  1.599954e-01 14.85542 2.928419e-01 0.596156983 0.77759606
#> Gene_10 Gene_10  5.487471e-01 14.42792 2.586612e+00 0.128101479 0.39057322
#> Gene_11 Gene_11  6.538953e-01 14.58179 4.828320e+00 0.043702412 0.27517630
#> Gene_12 Gene_12  1.945438e-01 14.65742 3.309318e-01 0.573431405 0.77759606
#> Gene_13 Gene_13  1.176655e+00 13.65571 1.655902e+01 0.000960716 0.02882148
#> Gene_14 Gene_14  8.166494e-01 14.40094 6.784185e+00 0.019614268 0.23537122
#> Gene_15 Gene_15  6.004514e-01 14.25517 3.848478e+00 0.068150555 0.28099865
#> Gene_16 Gene_16 -1.284917e-01 14.03806 1.476561e-01 0.706046162 0.85799160
#> Gene_17 Gene_17  3.504218e-02 13.84443 1.478805e-02 0.904785635 0.96941318
#> Gene_18 Gene_18 -4.258217e-01 13.07585 2.556344e+00 0.130191074 0.39057322
#> Gene_19 Gene_19 -3.574471e-01 13.20907 2.079732e+00 0.169321041 0.42368029
#> Gene_20 Gene_20  3.817123e-01 14.34034 1.257287e+00 0.279367850 0.55873570
#> Gene_21 Gene_21 -5.521256e-01 13.47104 4.306360e+00 0.055140368 0.28099865
#> Gene_22 Gene_22 -7.497055e-01 13.36652 8.146085e+00 0.011833590 0.23537122
#> Gene_23 Gene_23 -4.150511e-01 13.88186 2.093072e+00 0.168041553 0.42368029
#> Gene_24 Gene_24 -5.696138e-01 14.32594 3.784301e+00 0.070249662 0.28099865
#> Gene_25 Gene_25 -5.457576e-01 14.44601 3.205418e+00 0.093096359 0.32857538
#> Gene_26 Gene_26 -1.435495e-01 14.00958 2.341711e-01 0.635266118 0.81097802
#> Gene_27 Gene_27 -2.285939e-01 13.94197 4.306683e-01 0.521359349 0.74479907
#> Gene_28 Gene_28 -4.184076e-02 13.25524 2.101028e-02 0.886633502 0.96723655
#> Gene_29 Gene_29  2.066557e-02 14.05210 5.046172e-03 0.944284272 0.97218763
#> Gene_30 Gene_30 -7.942410e-01 13.82018 6.838140e+00 0.019211489 0.23537122
#> Gene_31 Gene_31 -3.533610e-01 14.03673 1.187289e+00 0.292655965 0.56643090
#> Gene_32 Gene_32  1.174641e-01 14.15167 1.383590e-01 0.714993001 0.85799160
#> Gene_33 Gene_33 -1.606162e-01 14.06180 2.998297e-01 0.591847526 0.77759606
#> Gene_34 Gene_34 -6.697999e-01 13.51131 5.192875e+00 0.037344978 0.27517630
#> Gene_35 Gene_35  3.299970e-01 14.07297 1.289504e+00 0.273526213 0.55873570
#> Gene_36 Gene_36  2.013795e-01 13.22965 5.948548e-01 0.452235707 0.68681783
#> Gene_37 Gene_37  2.171497e-01 13.66426 5.798783e-01 0.457878556 0.68681783
#> Gene_38 Gene_38 -4.311995e-01 13.90402 2.006781e+00 0.176533454 0.42368029
#> Gene_39 Gene_39  2.283944e-01 12.98891 8.440496e-01 0.372426167 0.68681783
#> Gene_40 Gene_40 -2.578351e-01 14.16743 7.123154e-01 0.411611858 0.68681783
#> Gene_41 Gene_41 -3.003584e-01 14.08788 7.024743e-01 0.414789659 0.68681783
#> Gene_42 Gene_42  2.297639e-01 13.59758 4.846403e-01 0.496714734 0.72689961
#> Gene_43 Gene_43  9.793853e-02 14.05362 1.106754e-01 0.743867872 0.87513867
#> Gene_44 Gene_44  1.529073e-02 13.41191 3.147207e-03 0.955984507 0.97218763
#> Gene_45 Gene_45  7.416194e-02 13.41949 7.239424e-02 0.791458225 0.91322103
#> Gene_46 Gene_46 -7.255692e-01 13.93750 3.810639e+00 0.069379263 0.28099865
#> Gene_47 Gene_47 -5.884604e-05 13.50544 5.181635e-08 0.999821302 0.99982130
#> Gene_48 Gene_48 -2.531384e-01 14.24921 6.374833e-01 0.436779547 0.68681783
#> Gene_49 Gene_49 -2.930161e-01 14.03696 6.157780e-01 0.444540258 0.68681783
#> Gene_50 Gene_50 -1.632461e-01 14.40896 3.027851e-01 0.590044619 0.77759606
#> Gene_51 Gene_51  3.623413e-01 14.19317 1.480713e+00 0.242006955 0.52467350
#> Gene_52 Gene_52 -8.621055e-01 13.79145 4.836415e+00 0.043548525 0.27517630
#> Gene_53 Gene_53 -5.744205e-01 14.50217 3.645005e+00 0.075074344 0.28152879
#> Gene_54 Gene_54 -2.446537e-01 14.29948 6.812593e-01 0.421771937 0.68681783
#> Gene_55 Gene_55  5.082222e-01 13.94117 2.810622e+00 0.113845696 0.37948565
#> Gene_56 Gene_56 -7.467621e-02 14.06516 5.081873e-02 0.824612789 0.93178483
#> Gene_57 Gene_57 -4.559174e-01 14.32490 1.462159e+00 0.244847631 0.52467350
#> Gene_58 Gene_58  4.744774e-01 14.29319 2.058091e+00 0.171420964 0.42368029
#> Gene_59 Gene_59 -1.006405e+00 12.60231 2.025826e+01 0.000398335 0.02390010
#> Gene_60 Gene_60 -3.590200e-01 13.74176 1.590569e+00 0.226040127 0.52163106
#>         is_da
#> Gene_01 FALSE
#> Gene_02 FALSE
#> Gene_03 FALSE
#> Gene_04 FALSE
#> Gene_05 FALSE
#> Gene_06 FALSE
#> Gene_07 FALSE
#> Gene_08 FALSE
#> Gene_09 FALSE
#> Gene_10 FALSE
#> Gene_11 FALSE
#> Gene_12 FALSE
#> Gene_13  TRUE
#> Gene_14 FALSE
#> Gene_15 FALSE
#> Gene_16 FALSE
#> Gene_17 FALSE
#> Gene_18 FALSE
#> Gene_19 FALSE
#> Gene_20 FALSE
#> Gene_21 FALSE
#> Gene_22 FALSE
#> Gene_23 FALSE
#> Gene_24 FALSE
#> Gene_25 FALSE
#> Gene_26 FALSE
#> Gene_27 FALSE
#> Gene_28 FALSE
#> Gene_29 FALSE
#> Gene_30 FALSE
#> Gene_31 FALSE
#> Gene_32 FALSE
#> Gene_33 FALSE
#> Gene_34 FALSE
#> Gene_35 FALSE
#> Gene_36 FALSE
#> Gene_37 FALSE
#> Gene_38 FALSE
#> Gene_39 FALSE
#> Gene_40 FALSE
#> Gene_41 FALSE
#> Gene_42 FALSE
#> Gene_43 FALSE
#> Gene_44 FALSE
#> Gene_45 FALSE
#> Gene_46 FALSE
#> Gene_47 FALSE
#> Gene_48 FALSE
#> Gene_49 FALSE
#> Gene_50 FALSE
#> Gene_51 FALSE
#> Gene_52 FALSE
#> Gene_53 FALSE
#> Gene_54 FALSE
#> Gene_55 FALSE
#> Gene_56 FALSE
#> Gene_57 FALSE
#> Gene_58 FALSE
#> Gene_59  TRUE
#> Gene_60 FALSE
#> 
#> 
#> $signal$rna$de_results_sim
#> $signal$rna$de_results_sim$TypeA_vs_TypeB
#>         feature        logFC   logCPM       F_stat       PValue        FDR
#> Gene_01 Gene_01  0.375749645 14.49426 1.192262e+00 0.2770610558 0.69754319
#> Gene_02 Gene_02  1.075103099 14.62498 9.710716e+00 0.0022925114 0.05920757
#> Gene_03 Gene_03  0.498288251 14.30499 2.014196e+00 0.1584261280 0.55915104
#> Gene_04 Gene_04  0.883949172 14.45905 6.417335e+00 0.0125905964 0.12590596
#> Gene_05 Gene_05  0.204096815 14.14743 3.309457e-01 0.5661799052 0.93260381
#> Gene_06 Gene_06  0.997846175 14.60584 8.372152e+00 0.0045264471 0.06789671
#> Gene_07 Gene_07 -0.028119282 14.70488 6.986109e-03 0.9335272965 0.98927371
#> Gene_08 Gene_08  0.322259554 14.25937 8.393349e-01 0.3614245040 0.82247059
#> Gene_09 Gene_09  0.732128103 14.83749 4.793690e+00 0.0304990871 0.23961319
#> Gene_10 Gene_10 -0.031803078 14.41976 8.461955e-03 0.9268602616 0.98927371
#> Gene_11 Gene_11  1.235027561 14.50173 1.236687e+01 0.0006177416 0.03706450
#> Gene_12 Gene_12  0.400141390 14.77652 1.428736e+00 0.2343276286 0.64360885
#> Gene_13 Gene_13  1.123895479 13.76307 9.204300e+00 0.0029603783 0.05920757
#> Gene_14 Gene_14 -0.071604928 13.91126 3.961578e-02 0.8425705368 0.97219677
#> Gene_15 Gene_15  0.691117584 14.34789 3.875465e+00 0.0513021017 0.34201401
#> Gene_16 Gene_16 -0.628533862 14.00201 3.052444e+00 0.0831741068 0.38388049
#> Gene_17 Gene_17 -0.026742474 13.71780 5.403187e-03 0.9415255584 0.98927371
#> Gene_18 Gene_18 -0.121627951 12.87738 9.661010e-02 0.7564754619 0.93260381
#> Gene_19 Gene_19 -0.510395126 13.15749 1.795192e+00 0.1828265768 0.60942192
#> Gene_20 Gene_20 -0.091238584 14.39624 6.936641e-02 0.7927132365 0.93260381
#> Gene_21 Gene_21  0.161695977 13.62302 1.934948e-01 0.6608149241 0.93260381
#> Gene_22 Gene_22  0.007704124 13.28375 4.218245e-04 0.9836480175 0.99457218
#> Gene_23 Gene_23 -0.145636522 13.81768 1.628195e-01 0.6872919150 0.93260381
#> Gene_24 Gene_24 -0.509428601 14.44576 2.169023e+00 0.1434337163 0.53787644
#> Gene_25 Gene_25 -0.194645496 14.33385 3.115047e-01 0.5777981403 0.93260381
#> Gene_26 Gene_26  0.372428007 13.65403 1.026974e+00 0.3129094871 0.75098277
#> Gene_27 Gene_27 -0.192095382 13.94261 2.858004e-01 0.5939135120 0.93260381
#> Gene_28 Gene_28 -0.821253386 13.32018 4.710515e+00 0.0319484260 0.23961319
#> Gene_29 Gene_29 -0.205193498 14.31062 3.446941e-01 0.5582351552 0.93260381
#> Gene_30 Gene_30 -0.323417193 13.96512 8.093585e-01 0.3701117645 0.82247059
#> Gene_31 Gene_31 -0.101507753 13.97750 8.008337e-02 0.7776720551 0.93260381
#> Gene_32 Gene_32  0.019655378 14.00818 3.015600e-03 0.9562979172 0.98927371
#> Gene_33 Gene_33  0.113548811 14.47825 1.089370e-01 0.7419331525 0.93260381
#> Gene_34 Gene_34 -0.665544984 13.88667 3.378597e+00 0.0685214872 0.37375357
#> Gene_35 Gene_35 -0.262969111 14.33085 5.680127e-01 0.4525259164 0.90505183
#> Gene_36 Gene_36 -0.249706324 13.31289 4.431156e-01 0.5069000588 0.93260381
#> Gene_37 Gene_37 -0.663062718 13.61473 3.224073e+00 0.0750794932 0.37539747
#> Gene_38 Gene_38  0.393653528 13.85198 1.182548e+00 0.2790172760 0.69754319
#> Gene_39 Gene_39 -0.054950891 12.90967 1.989400e-02 0.8880700229 0.98674447
#> Gene_40 Gene_40  0.206734173 13.72712 3.225255e-01 0.5711549217 0.93260381
#> Gene_41 Gene_41 -0.453939197 13.76217 1.564491e+00 0.2134421384 0.64032642
#> Gene_42 Gene_42 -0.182627595 13.92660 2.573636e-01 0.6128682463 0.93260381
#> Gene_43 Gene_43 -0.574297282 13.99965 2.544824e+00 0.1132866091 0.45314644
#> Gene_44 Gene_44 -1.001420190 13.80208 7.501064e+00 0.0071064840 0.08527781
#> Gene_45 Gene_45  0.283591497 13.57816 5.894405e-01 0.4441436660 0.90505183
#> Gene_46 Gene_46  0.113745898 13.83736 9.926623e-02 0.7532591832 0.93260381
#> Gene_47 Gene_47 -0.694913552 13.71083 3.600351e+00 0.0601703439 0.36102206
#> Gene_48 Gene_48 -0.054627545 14.34600 2.462523e-02 0.8755682381 0.98674447
#> Gene_49 Gene_49 -0.609329660 13.57247 2.709499e+00 0.1023689929 0.43872426
#> Gene_50 Gene_50  0.217457517 14.25783 3.828951e-01 0.5372307517 0.93260381
#> Gene_51 Gene_51  0.002459299 13.78958 4.647149e-05 0.9945721830 0.99457218
#> Gene_52 Gene_52 -0.464900292 14.04410 1.685465e+00 0.1966894015 0.62112443
#> Gene_53 Gene_53 -0.281102723 14.39668 6.569469e-01 0.4192443936 0.89838084
#> Gene_54 Gene_54  0.101191808 14.31455 8.397152e-02 0.7724858793 0.93260381
#> Gene_55 Gene_55 -0.151643149 13.92262 1.778540e-01 0.6739785682 0.93260381
#> Gene_56 Gene_56 -0.222506687 14.05566 3.884317e-01 0.5343076833 0.93260381
#> Gene_57 Gene_57  0.148717706 14.14822 1.759297e-01 0.6756446789 0.93260381
#> Gene_58 Gene_58  0.428087481 14.01531 1.418567e+00 0.2359899100 0.64360885
#> Gene_59 Gene_59 -0.171614401 12.70550 1.839392e-01 0.6687784834 0.93260381
#> Gene_60 Gene_60 -0.119522083 14.14245 1.136218e-01 0.7366467302 0.93260381
#>         is_da
#> Gene_01 FALSE
#> Gene_02 FALSE
#> Gene_03 FALSE
#> Gene_04 FALSE
#> Gene_05 FALSE
#> Gene_06 FALSE
#> Gene_07 FALSE
#> Gene_08 FALSE
#> Gene_09 FALSE
#> Gene_10 FALSE
#> Gene_11  TRUE
#> Gene_12 FALSE
#> Gene_13 FALSE
#> Gene_14 FALSE
#> Gene_15 FALSE
#> Gene_16 FALSE
#> Gene_17 FALSE
#> Gene_18 FALSE
#> Gene_19 FALSE
#> Gene_20 FALSE
#> Gene_21 FALSE
#> Gene_22 FALSE
#> Gene_23 FALSE
#> Gene_24 FALSE
#> Gene_25 FALSE
#> Gene_26 FALSE
#> Gene_27 FALSE
#> Gene_28 FALSE
#> Gene_29 FALSE
#> Gene_30 FALSE
#> Gene_31 FALSE
#> Gene_32 FALSE
#> Gene_33 FALSE
#> Gene_34 FALSE
#> Gene_35 FALSE
#> Gene_36 FALSE
#> Gene_37 FALSE
#> Gene_38 FALSE
#> Gene_39 FALSE
#> Gene_40 FALSE
#> Gene_41 FALSE
#> Gene_42 FALSE
#> Gene_43 FALSE
#> Gene_44 FALSE
#> Gene_45 FALSE
#> Gene_46 FALSE
#> Gene_47 FALSE
#> Gene_48 FALSE
#> Gene_49 FALSE
#> Gene_50 FALSE
#> Gene_51 FALSE
#> Gene_52 FALSE
#> Gene_53 FALSE
#> Gene_54 FALSE
#> Gene_55 FALSE
#> Gene_56 FALSE
#> Gene_57 FALSE
#> Gene_58 FALSE
#> Gene_59 FALSE
#> Gene_60 FALSE
#> 
#> 
#> $signal$rna$ml_classification
#> $signal$rna$ml_classification$accuracy
#> [1] 0.6875
#> 
#> $signal$rna$ml_classification$precision
#> [1] 0.6875
#> 
#> $signal$rna$ml_classification$recall
#> [1] 0.6904762
#> 
#> $signal$rna$ml_classification$F1
#> [1] 0.6862745
#> 
#> 
#> 
#> $signal$atac
#> $signal$atac$da_summary_table
#>                             SubCategory          Metric      Value
#> 1 scATAC-seq Differential Accessibility pval_dist_score 0.85633794
#> 2 scATAC-seq Differential Accessibility        accuracy 0.81666667
#> 3 scATAC-seq Differential Accessibility       precision 0.50000000
#> 4 scATAC-seq Differential Accessibility          recall 0.09090909
#> 5 scATAC-seq Differential Accessibility        f1_score 0.15384615
#> 6 scATAC-seq Differential Accessibility             auc 0.83673469
#>                                                         Direction
#> 1    Higher is better (1.0 = uniform null p-values, no inflation)
#> 2             Higher is better (1.0 = perfect DAR classification)
#> 3                 Higher is better (1.0 = zero false discoveries)
#> 4                 Higher is better (1.0 = 100% true DAR recovery)
#> 5 Higher is better (1.0 = harmonic balance of precision & recall)
#> 6         Higher is better (1.0 = perfect ranking discrimination)
#>                                                      Framework
#> 1 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 2 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 3 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 4 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 5 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 6 Libra Pseudobulk edgeR-QLF (Squair 2021; Teo 2024; Duo 2024)
#> 
#> $signal$atac$pairwise_results
#>             Pair                           SubCategory N_total N_true_DAR
#> 1 TypeA_vs_TypeB scATAC-seq Differential Accessibility      60         11
#>   N_sim_DAR TP FP FN TN pval_dist_score  accuracy precision     recall
#> 1         2  1  1 10 48       0.8563379 0.8166667       0.5 0.09090909
#>    f1_score       auc LogFC_Pearson_r LogFC_Spearman_rho Top_DAR_Jaccard
#> 1 0.1538462 0.8367347       0.5558168          0.4911988               1
#> 
#> $signal$atac$pseudobulks_ref
#> $signal$atac$pseudobulks_ref$TypeA
#>         Batch1 Batch2
#> Peak_01      3      5
#> Peak_02      4      6
#> Peak_03      5      7
#> Peak_04      4      6
#> Peak_05      7      9
#> Peak_06      4      8
#> Peak_07      6     10
#> Peak_08      9      7
#> Peak_09      7      5
#> Peak_10     10      8
#> Peak_11     10      6
#> Peak_12      7     12
#> Peak_13      6      6
#> Peak_14     10      4
#> Peak_15      8      7
#> Peak_16      7      7
#> Peak_17      9      8
#> Peak_18      6      3
#> Peak_19      8     11
#> Peak_20     10      9
#> Peak_21      7      6
#> Peak_22      8     11
#> Peak_23      7      7
#> Peak_24     10      6
#> Peak_25      3      3
#> Peak_26      6      9
#> Peak_27      8      5
#> Peak_28      6      5
#> Peak_29     13     10
#> Peak_30      5      4
#> Peak_31     11      9
#> Peak_32     10      6
#> Peak_33      7      8
#> Peak_34      7      8
#> Peak_35      8     12
#> Peak_36      9      8
#> Peak_37      7      4
#> Peak_38      7      5
#> Peak_39      6      7
#> Peak_40     11      9
#> Peak_41      5      5
#> Peak_42      9     11
#> Peak_43      6      7
#> Peak_44     10      6
#> Peak_45     11     10
#> Peak_46      5      5
#> Peak_47      7      9
#> Peak_48      7     10
#> Peak_49      3      8
#> Peak_50      7      5
#> Peak_51      5      3
#> Peak_52      5      8
#> Peak_53      4      8
#> Peak_54     10      4
#> Peak_55     11      7
#> Peak_56     10      7
#> Peak_57     10     14
#> Peak_58     11      8
#> Peak_59      7     10
#> Peak_60      7      7
#> 
#> $signal$atac$pseudobulks_ref$TypeB
#>         Batch1 Batch2
#> Peak_01     11     11
#> Peak_02     12     16
#> Peak_03     13     12
#> Peak_04     17     10
#> Peak_05     15     18
#> Peak_06     13     13
#> Peak_07     12     12
#> Peak_08     13     14
#> Peak_09     15     16
#> Peak_10     17     17
#> Peak_11     17     14
#> Peak_12     14     18
#> Peak_13     15     14
#> Peak_14     13     11
#> Peak_15     16     15
#> Peak_16      7      5
#> Peak_17      6      9
#> Peak_18      5      3
#> Peak_19      7      6
#> Peak_20      6      7
#> Peak_21      5      5
#> Peak_22      8      7
#> Peak_23      6      7
#> Peak_24     10      9
#> Peak_25      4      6
#> Peak_26      6      9
#> Peak_27      8      8
#> Peak_28     11     10
#> Peak_29      9      9
#> Peak_30      5      9
#> Peak_31      8      5
#> Peak_32      6      9
#> Peak_33      6      4
#> Peak_34     11      9
#> Peak_35      9      9
#> Peak_36      7     12
#> Peak_37      5      7
#> Peak_38      6      5
#> Peak_39      8      7
#> Peak_40      7      5
#> Peak_41      1      7
#> Peak_42     13     11
#> Peak_43     12      8
#> Peak_44      7      8
#> Peak_45      9     11
#> Peak_46      9      7
#> Peak_47      4      6
#> Peak_48     10     13
#> Peak_49      2      2
#> Peak_50      8      7
#> Peak_51      5      3
#> Peak_52      3      6
#> Peak_53      5      4
#> Peak_54      8     11
#> Peak_55     10     12
#> Peak_56     12      6
#> Peak_57     12      8
#> Peak_58      9     10
#> Peak_59      8     10
#> Peak_60      7      4
#> 
#> 
#> $signal$atac$pseudobulks_sim
#> $signal$atac$pseudobulks_sim$TypeA
#>         Batch1 Batch2
#> Peak_01      6      4
#> Peak_02      7      3
#> Peak_03      4      3
#> Peak_04      3      6
#> Peak_05     10      4
#> Peak_06     11      4
#> Peak_07      7      1
#> Peak_08      8      6
#> Peak_09      3      4
#> Peak_10      6      9
#> Peak_11     12      8
#> Peak_12      7     13
#> Peak_13      6      9
#> Peak_14      6     10
#> Peak_15     12      9
#> Peak_16      7      4
#> Peak_17      9      8
#> Peak_18      2      5
#> Peak_19      5      6
#> Peak_20      2      6
#> Peak_21      3      4
#> Peak_22      8     10
#> Peak_23      5      5
#> Peak_24     11      9
#> Peak_25      5     10
#> Peak_26      7      5
#> Peak_27     10      6
#> Peak_28      6      6
#> Peak_29      6     10
#> Peak_30      2      4
#> Peak_31      5      6
#> Peak_32     10      4
#> Peak_33      5      4
#> Peak_34      6      6
#> Peak_35      9      6
#> Peak_36      6      4
#> Peak_37      5      9
#> Peak_38      5      5
#> Peak_39      2      7
#> Peak_40      8      9
#> Peak_41      2      6
#> Peak_42     10     13
#> Peak_43      8     10
#> Peak_44     10      5
#> Peak_45     10     10
#> Peak_46      6      4
#> Peak_47      9      6
#> Peak_48     10     11
#> Peak_49      5      4
#> Peak_50     12      9
#> Peak_51      5      7
#> Peak_52      4      3
#> Peak_53      7      8
#> Peak_54      9      5
#> Peak_55      5     13
#> Peak_56      8      7
#> Peak_57     16      6
#> Peak_58      7      6
#> Peak_59     10      7
#> Peak_60      7      5
#> 
#> $signal$atac$pseudobulks_sim$TypeB
#>         Batch1 Batch2
#> Peak_01     14     14
#> Peak_02      9     15
#> Peak_03     14     14
#> Peak_04     11     14
#> Peak_05     14     16
#> Peak_06     12     16
#> Peak_07     16     13
#> Peak_08     15     14
#> Peak_09     11      7
#> Peak_10     16     15
#> Peak_11     17     15
#> Peak_12     13     17
#> Peak_13     13     16
#> Peak_14     16     14
#> Peak_15     14     12
#> Peak_16      9     10
#> Peak_17      8      5
#> Peak_18      2      3
#> Peak_19      5      6
#> Peak_20      7      4
#> Peak_21      8      3
#> Peak_22      6      7
#> Peak_23      6      6
#> Peak_24     10      9
#> Peak_25      8      5
#> Peak_26      7      7
#> Peak_27      5      9
#> Peak_28      9      7
#> Peak_29      9     11
#> Peak_30      4      4
#> Peak_31      7      6
#> Peak_32      6      7
#> Peak_33      5      7
#> Peak_34     11      7
#> Peak_35      9      5
#> Peak_36      8      9
#> Peak_37      4      8
#> Peak_38      5      4
#> Peak_39      4      5
#> Peak_40      7      9
#> Peak_41      4      5
#> Peak_42      6      9
#> Peak_43     11      6
#> Peak_44      9      6
#> Peak_45      8      8
#> Peak_46      5      7
#> Peak_47      8      4
#> Peak_48      6     11
#> Peak_49      7      5
#> Peak_50      6      6
#> Peak_51      6      3
#> Peak_52      5      6
#> Peak_53      3      5
#> Peak_54      5     11
#> Peak_55     16     12
#> Peak_56      7      7
#> Peak_57      8     10
#> Peak_58      7     11
#> Peak_59      7      3
#> Peak_60      7      7
#> 
#> 
#> $signal$atac$da_results_ref
#> $signal$atac$da_results_ref$TypeA_vs_TypeB
#>         feature       logFC   logCPM       F_stat       PValue         FDR
#> Peak_01 Peak_01  1.28531261 14.22614 11.540251884 0.0009244519 0.011093423
#> Peak_02 Peak_02  1.31598148 14.50177 15.402310051 0.0001453586 0.004568652
#> Peak_03 Peak_03  0.89663746 14.47007  7.365534953 0.0076281536 0.041608111
#> Peak_04 Peak_04  1.26401102 14.47006 13.885816108 0.0002976166 0.005952331
#> Peak_05 Peak_05  0.88549848 14.81110 10.444579056 0.0015881324 0.015881324
#> Peak_06 Peak_06  0.95264340 14.50178  8.514421397 0.0042077395 0.028051596
#> Peak_07 Peak_07  0.43033667 14.56319  1.949167601 0.1652540108 0.367231135
#> Peak_08 Peak_08  0.59852365 14.65064  4.077300889 0.0456913356 0.152304452
#> Peak_09 Peak_09  1.20406871 14.65062 15.302990463 0.0001522884 0.004568652
#> Peak_10 Peak_10  0.76064935 14.88510  8.607271703 0.0040123221 0.028051596
#> Peak_11 Peak_11  0.79603617 14.75956  8.006089246 0.0054664471 0.032798682
#> Peak_12 Peak_12  0.59686319 14.86087  5.215576039 0.0241439580 0.103474106
#> Peak_13 Peak_13  1.10868988 14.59292 12.378623328 0.0006142302 0.009213453
#> Peak_14 Peak_14  0.61997547 14.50177  3.743741204 0.0553589620 0.158168463
#> Peak_15 Peak_15  0.88774186 14.73309  9.545246130 0.0024917097 0.021357512
#> Peak_16 Peak_16 -0.36447582 14.06570  0.847826649 0.3590152736 0.566866221
#> Peak_17 Peak_17 -0.32460904 14.30017  0.880196597 0.3500323240 0.566866221
#> Peak_18 Peak_18 -0.30991147 13.62208  0.308556118 0.5796016620 0.739917015
#> Peak_19 Peak_19 -0.68561609 14.30018  3.875542327 0.0512998099 0.153899430
#> Peak_20 Peak_20 -0.68561609 14.30017  3.876416004 0.0512739925 0.153899430
#> Peak_21 Peak_21 -0.51626906 13.93243  1.375398233 0.2432093521 0.461824530
#> Peak_22 Peak_22 -0.48298338 14.37057  2.070471061 0.1527783764 0.364775027
#> Peak_23 Peak_23 -0.25137983 14.10752  0.433535058 0.5115207278 0.682027637
#> Peak_24 Peak_24  0.09740973 14.40449  0.086984027 0.7685565786 0.853951754
#> Peak_25 Peak_25  0.56937183 13.56319  0.908256599 0.3424931616 0.566866221
#> Peak_26 Peak_26 -0.14668022 14.22617  0.170279736 0.6805996416 0.787464022
#> Peak_27 Peak_27  0.14787656 14.18768  0.168176410 0.6824688189 0.787464022
#> Peak_28 Peak_28  0.77094771 14.30014  4.704472039 0.0320565347 0.120212005
#> Peak_29 Peak_29 -0.49663556 14.59295  2.726646077 0.1013010009 0.264263481
#> Peak_30 Peak_30  0.47763510 13.93242  1.152464806 0.2851877882 0.518523251
#> Peak_31 Peak_31 -0.75872735 14.33580  4.878430508 0.0290937213 0.116374885
#> Peak_32 Peak_32 -0.23838590 14.26363  0.461216281 0.4983638845 0.679587115
#> Peak_33 Peak_33 -0.71927045 14.02264  3.034056376 0.0840963051 0.229353559
#> Peak_34 Peak_34  0.26213113 14.40449  0.625291992 0.4306475528 0.630215931
#> Peak_35 Peak_35 -0.29720802 14.50181  0.894639487 0.3461238190 0.566866221
#> Peak_36 Peak_36  0.01118660 14.43766  0.001187637 0.9725659140 0.972565914
#> Peak_37 Peak_37 -0.02309002 13.93242  0.002761925 0.9581746416 0.972565914
#> Peak_38 Peak_38 -0.26880382 13.93243  0.374948045 0.5414780071 0.706275661
#> Peak_39 Peak_39  0.05632117 14.14816  0.023188601 0.8792233825 0.909541430
#> Peak_40 Peak_40 -0.87182334 14.30017  6.189300867 0.0142251628 0.071125814
#> Peak_41 Peak_41 -0.45820680 13.67868  0.738617126 0.3918178407 0.602796678
#> Peak_42 Peak_42  0.11262793 14.67867  0.153694146 0.6957251917 0.787613425
#> Peak_43 Peak_43  0.46513252 14.33578  1.819333952 0.1799324714 0.385569582
#> Peak_44 Peak_44 -0.23838590 14.26363  0.461281312 0.4983336603 0.679587115
#> Peak_45 Peak_45 -0.21686753 14.59295  0.522636240 0.4711259920 0.673037131
#> Peak_46 Peak_46  0.51865685 14.06569  1.668316503 0.1989664410 0.411654705
#> Peak_47 Peak_47 -0.81097614 14.06571  4.112556493 0.0447804460 0.152304452
#> Peak_48 Peak_48  0.28339252 14.56320  0.853505208 0.3574165881 0.566866221
#> Peak_49 Peak_49 -1.54748112 13.50182  5.602992618 0.0195298585 0.090137808
#> Peak_50 Peak_50  0.16964664 14.10751  0.195889915 0.6588557942 0.787464022
#> Peak_51 Peak_51 -0.14460782 13.56319  0.060624347 0.8059325121 0.879199104
#> Peak_52 Peak_52 -0.66417154 13.88514  2.108434046 0.1490987469 0.364775027
#> Peak_53 Peak_53 -0.55084607 13.83623  1.357359258 0.2463064160 0.461824530
#> Peak_54 Peak_54  0.28704852 14.33577  0.701030106 0.4041021870 0.606153281
#> Peak_55 Peak_55  0.13864382 14.56318  0.205762503 0.6509284399 0.787464022
#> Peak_56 Peak_56 -0.06572241 14.40449  0.039746560 0.8423143530 0.886646687
#> Peak_57 Peak_57 -0.40749793 14.67869  2.017684112 0.1580691785 0.364775027
#> Peak_58 Peak_58 -0.14718774 14.50179  0.219570926 0.6402171365 0.787464022
#> Peak_59 Peak_59 -0.06572241 14.40451  0.039756379 0.8422951343 0.886646687
#> Peak_60 Peak_60 -0.48719757 14.02263  1.407326552 0.2378441759 0.461824530
#>         is_da
#> Peak_01  TRUE
#> Peak_02  TRUE
#> Peak_03  TRUE
#> Peak_04  TRUE
#> Peak_05  TRUE
#> Peak_06  TRUE
#> Peak_07 FALSE
#> Peak_08 FALSE
#> Peak_09  TRUE
#> Peak_10  TRUE
#> Peak_11  TRUE
#> Peak_12 FALSE
#> Peak_13  TRUE
#> Peak_14 FALSE
#> Peak_15  TRUE
#> Peak_16 FALSE
#> Peak_17 FALSE
#> Peak_18 FALSE
#> Peak_19 FALSE
#> Peak_20 FALSE
#> Peak_21 FALSE
#> Peak_22 FALSE
#> Peak_23 FALSE
#> Peak_24 FALSE
#> Peak_25 FALSE
#> Peak_26 FALSE
#> Peak_27 FALSE
#> Peak_28 FALSE
#> Peak_29 FALSE
#> Peak_30 FALSE
#> Peak_31 FALSE
#> Peak_32 FALSE
#> Peak_33 FALSE
#> Peak_34 FALSE
#> Peak_35 FALSE
#> Peak_36 FALSE
#> Peak_37 FALSE
#> Peak_38 FALSE
#> Peak_39 FALSE
#> Peak_40 FALSE
#> Peak_41 FALSE
#> Peak_42 FALSE
#> Peak_43 FALSE
#> Peak_44 FALSE
#> Peak_45 FALSE
#> Peak_46 FALSE
#> Peak_47 FALSE
#> Peak_48 FALSE
#> Peak_49 FALSE
#> Peak_50 FALSE
#> Peak_51 FALSE
#> Peak_52 FALSE
#> Peak_53 FALSE
#> Peak_54 FALSE
#> Peak_55 FALSE
#> Peak_56 FALSE
#> Peak_57 FALSE
#> Peak_58 FALSE
#> Peak_59 FALSE
#> Peak_60 FALSE
#> 
#> 
#> $signal$atac$da_results_sim
#> $signal$atac$da_results_sim$TypeA_vs_TypeB
#>         feature       logFC   logCPM       F_stat       PValue        FDR is_da
#> Peak_01 Peak_01  1.25446850 14.59948 7.586161e+00 0.0067979617 0.10691384 FALSE
#> Peak_02 Peak_02  1.03435648 14.46824 5.376410e+00 0.0221029709 0.16577228 FALSE
#> Peak_03 Peak_03  1.75510805 14.50218 1.359780e+01 0.0003414759 0.02048856  TRUE
#> Peak_04 Peak_04  1.24098975 14.46824 7.495383e+00 0.0071275892 0.10691384 FALSE
#> Peak_05 Peak_05  0.87702858 14.77637 3.872594e+00 0.0513870476 0.30506181 FALSE
#> Peak_06 Peak_06  0.68043706 14.74836 2.392691e+00 0.1245368453 0.43954181 FALSE
#> Peak_07 Peak_07  1.61840612 14.56777 1.183112e+01 0.0008018078 0.02405423  TRUE
#> Peak_08 Peak_08  0.82855990 14.74836 3.487400e+00 0.0642774646 0.30506181 FALSE
#> Peak_09 Peak_09  1.12526552 14.12034 5.865343e+00 0.0169379568 0.16577228 FALSE
#> Peak_10 Peak_10  0.82595269 14.83082 3.439836e+00 0.0660967264 0.30506181 FALSE
#> Peak_11 Peak_11  0.46182873 14.98283 1.097282e+00 0.2969686806 0.66541089 FALSE
#> Peak_12 Peak_12  0.36951881 14.93392 7.124508e-01 0.4003109531 0.71256709 FALSE
#> Peak_13 Peak_13  0.73059088 14.77637 2.737091e+00 0.1006564795 0.43138491 FALSE
#> Peak_14 Peak_14  0.68732239 14.83082 2.421468e+00 0.1223162273 0.43954181 FALSE
#> Peak_15 Peak_15  0.09543170 14.85732 4.877506e-02 0.8255839592 0.93462335 FALSE
#> Peak_16 Peak_16  0.56698231 14.32387 1.720615e+00 0.1921178489 0.57943439 FALSE
#> Peak_17 Peak_17 -0.58875225 14.32391 1.924845e+00 0.1678947611 0.55964920 FALSE
#> Peak_18 Peak_18 -0.66866607 13.39790 1.086087e+00 0.2994348999 0.66541089 FALSE
#> Peak_19 Peak_19 -0.20697673 13.98286 1.952197e-01 0.6594025811 0.82902037 FALSE
#> Peak_20 Peak_20  0.24133761 13.83085 2.263297e-01 0.6351236936 0.82902037 FALSE
#> Peak_21 Peak_21  0.42819336 13.77640 6.671542e-01 0.4156641353 0.71256709 FALSE
#> Peak_22 Peak_22 -0.67013634 14.36138 2.484602e+00 0.1175980884 0.43954181 FALSE
#> Peak_23 Peak_23  0.05021920 13.98285 1.142159e-02 0.9150691344 0.98063157 FALSE
#> Peak_24 Peak_24 -0.28205934 14.63056 4.379152e-01 0.5093991605 0.72771309 FALSE
#> Peak_25 Peak_25 -0.41076410 14.24590 9.363288e-01 0.3351699576 0.67244446 FALSE
#> Peak_26 Peak_26  0.01051796 14.16342 5.918404e-04 0.9806315710 0.98063157 FALSE
#> Peak_27 Peak_27 -0.39767143 14.32390 8.822018e-01 0.3494860786 0.67642467 FALSE
#> Peak_28 Peak_28  0.19976792 14.24588 2.187320e-01 0.6408560372 0.82902037 FALSE
#> Peak_29 Peak_29  0.10873999 14.53539 6.484054e-02 0.7994385242 0.92242907 FALSE
#> Peak_30 Peak_30  0.19629573 13.53539 1.102660e-01 0.7404208351 0.87108334 FALSE
#> Peak_31 Peak_31  0.02873680 14.07596 4.090317e-03 0.9491120090 0.98063157 FALSE
#> Peak_32 Peak_32 -0.31279510 14.20525 5.388756e-01 0.4643310019 0.72596357 FALSE
#> Peak_33 Peak_33  0.19859693 13.93394 1.693535e-01 0.6814210812 0.82902037 FALSE
#> Peak_34 Peak_34  0.36704674 14.32388 7.337291e-01 0.3933844286 0.71256709 FALSE
#> Peak_35 Peak_35 -0.30593424 14.28543 5.219352e-01 0.4714229369 0.72596357 FALSE
#> Peak_36 Peak_36  0.54340967 14.20523 1.563302e+00 0.2136150042 0.57943439 FALSE
#> Peak_37 Peak_37 -0.42584323 14.16343 9.735706e-01 0.3257766046 0.67244446 FALSE
#> Peak_38 Peak_38 -0.35433923 13.83086 4.959861e-01 0.4826336156 0.72596357 FALSE
#> Peak_39 Peak_39 -0.20596150 13.77640 1.589319e-01 0.6908503055 0.82902037 FALSE
#> Peak_40 Peak_40 -0.29467244 14.43352 4.832431e-01 0.4883008949 0.72596357 FALSE
#> Peak_41 Peak_41 -0.04055543 13.71982 5.807827e-03 0.9393796537 0.98063157 FALSE
#> Peak_42 Peak_42 -0.81722631 14.59955 3.621140e+00 0.0594457272 0.30506181 FALSE
#> Peak_43 Peak_43 -0.28999399 14.50224 4.661599e-01 0.4960751045 0.72596357 FALSE
#> Peak_44 Peak_44 -0.20820830 14.32390 2.419183e-01 0.6237226056 0.82902037 FALSE
#> Peak_45 Peak_45 -0.52622504 14.53541 1.527236e+00 0.2189414068 0.57943439 FALSE
#> Peak_46 Peak_46  0.05021920 13.98285 1.142027e-02 0.9150740112 0.98063157 FALSE
#> Peak_47 Peak_47 -0.52381222 14.20526 1.506236e+00 0.2221165162 0.57943439 FALSE
#> Peak_48 Peak_48 -0.50976488 14.59954 1.428429e+00 0.2343775912 0.58594398 FALSE
#> Peak_49 Peak_49  0.19859694 13.93394 1.693463e-01 0.6814275035 0.82902037 FALSE
#> Peak_50 Peak_50 -1.00295534 14.43354 5.430891e+00 0.0214533938 0.16577228 FALSE
#> Peak_51 Peak_51 -0.61191845 13.93396 1.619299e+00 0.2056496425 0.57943439 FALSE
#> Peak_52 Peak_52  0.42819335 13.77639 6.674835e-01 0.4155494049 0.71256709 FALSE
#> Peak_53 Peak_53 -1.09308454 14.03018 5.485528e+00 0.0208220688 0.16577228 FALSE
#> Peak_54 Peak_54 -0.01871528 14.32388 1.946987e-03 0.9648783627 0.98063157 FALSE
#> Peak_55 Peak_55  0.42106484 14.83084 9.322443e-01 0.3362222318 0.67244446 FALSE
#> Peak_56 Peak_56 -0.30593424 14.28542 5.221181e-01 0.4713454294 0.72596357 FALSE
#> Peak_57 Peak_57 -0.49494699 14.66094 1.341463e+00 0.2490760021 0.59778241 FALSE
#> Peak_58 Peak_58  0.25367316 14.36135 3.530982e-01 0.5534831479 0.77230207 FALSE
#> Peak_59 Peak_59 -0.95853821 14.20527 4.927859e+00 0.0283056924 0.18870462 FALSE
#> Peak_60 Peak_60  0.01051796 14.16342 5.918404e-04 0.9806315710 0.98063157 FALSE
#> 
#> 
#> 
#> 
#> $trajectory
#> $trajectory$him
#> [1] 1
#> 
#> $trajectory$edge_flip
#> [1] 1
#> 
#> $trajectory$geodesic_correlation
#> [1] 0.4848647
#> 
#> $trajectory$F1_branches
#> [1] 1
#> 
#> $trajectory$wcor_features
#> [1] 0.3802351
#> 
#> $trajectory$cor_features
#> [1] 0.3578252
#> 
#> $trajectory$n_hvgs
#> [1] 60
#> 
#> $trajectory$ref_trajectory
#> $id
#> [1] "20261113_053551__data_wrapper__IcLTIqicPk"
#> 
#> $cell_ids
#>  [1] "Cell_01" "Cell_02" "Cell_03" "Cell_04" "Cell_05" "Cell_06" "Cell_07"
#>  [8] "Cell_08" "Cell_09" "Cell_10" "Cell_11" "Cell_12" "Cell_13" "Cell_14"
#> [15] "Cell_15" "Cell_16" "Cell_17" "Cell_18" "Cell_19" "Cell_20" "Cell_21"
#> [22] "Cell_22" "Cell_23" "Cell_24" "Cell_25" "Cell_26" "Cell_27" "Cell_28"
#> [29] "Cell_29" "Cell_30" "Cell_31" "Cell_32" "Cell_33" "Cell_34" "Cell_35"
#> [36] "Cell_36" "Cell_37" "Cell_38" "Cell_39" "Cell_40" "Cell_41" "Cell_42"
#> [43] "Cell_43" "Cell_44" "Cell_45" "Cell_46" "Cell_47" "Cell_48" "Cell_49"
#> [50] "Cell_50" "Cell_51" "Cell_52" "Cell_53" "Cell_54" "Cell_55" "Cell_56"
#> [57] "Cell_57" "Cell_58" "Cell_59" "Cell_60" "Cell_61" "Cell_62" "Cell_63"
#> [64] "Cell_64" "Cell_65" "Cell_66" "Cell_67" "Cell_68" "Cell_69" "Cell_70"
#> [71] "Cell_71" "Cell_72" "Cell_73" "Cell_74" "Cell_75" "Cell_76" "Cell_77"
#> [78] "Cell_78" "Cell_79" "Cell_80"
#> 
#> $cell_info
#> # A tibble: 80 × 1
#>    cell_id
#>    <chr>  
#>  1 Cell_01
#>  2 Cell_02
#>  3 Cell_03
#>  4 Cell_04
#>  5 Cell_05
#>  6 Cell_06
#>  7 Cell_07
#>  8 Cell_08
#>  9 Cell_09
#> 10 Cell_10
#> # ℹ 70 more rows
#> 
#> $feature_ids
#> NULL
#> 
#> $feature_info
#> NULL
#> 
#> $milestone_ids
#> [1] "TypeA" "TypeB"
#> 
#> $milestone_network
#>    from    to   length directed
#> 1 TypeA TypeB 1.543713    FALSE
#> 
#> $divergence_regions
#> # A tibble: 0 × 3
#> # ℹ 3 variables: divergence_id <chr>, milestone_id <chr>, is_start <lgl>
#> 
#> $milestone_percentages
#>     cell_id milestone_id percentage
#> 1   Cell_01        TypeA 0.86808122
#> 2   Cell_02        TypeA 1.00000000
#> 3   Cell_03        TypeA 0.32642222
#> 4   Cell_04        TypeA 1.00000000
#> 5   Cell_05        TypeA 0.56074964
#> 6   Cell_06        TypeA 1.00000000
#> 7   Cell_07        TypeA 1.00000000
#> 8   Cell_08        TypeA 1.00000000
#> 9   Cell_09        TypeA 1.00000000
#> 10  Cell_10        TypeA 1.00000000
#> 11  Cell_11        TypeA 1.00000000
#> 12  Cell_12        TypeA 1.00000000
#> 13  Cell_13        TypeA 0.98272020
#> 14  Cell_14        TypeA 1.00000000
#> 15  Cell_15        TypeA 0.60710083
#> 16  Cell_16        TypeA 0.49702075
#> 17  Cell_17        TypeA 0.51749571
#> 18  Cell_18        TypeA 1.00000000
#> 19  Cell_19        TypeA 0.63126146
#> 20  Cell_20        TypeA 1.00000000
#> 21  Cell_21        TypeA 0.61462635
#> 22  Cell_22        TypeA 1.00000000
#> 23  Cell_23        TypeA 1.00000000
#> 24  Cell_24        TypeA 0.00000000
#> 25  Cell_25        TypeA 0.45533598
#> 26  Cell_26        TypeA 0.50025326
#> 27  Cell_27        TypeA 0.49310324
#> 28  Cell_28        TypeA 1.00000000
#> 29  Cell_29        TypeA 0.22345545
#> 30  Cell_30        TypeA 0.44135267
#> 31  Cell_31        TypeA 1.00000000
#> 32  Cell_32        TypeA 1.00000000
#> 33  Cell_33        TypeA 1.00000000
#> 34  Cell_34        TypeA 1.00000000
#> 35  Cell_35        TypeA 0.76460593
#> 36  Cell_36        TypeA 0.94914199
#> 37  Cell_37        TypeA 0.51675572
#> 38  Cell_38        TypeA 0.49801588
#> 39  Cell_39        TypeA 1.00000000
#> 40  Cell_40        TypeA 0.00000000
#> 41  Cell_41        TypeA 0.00000000
#> 42  Cell_42        TypeA 0.00000000
#> 43  Cell_43        TypeA 0.00000000
#> 44  Cell_44        TypeA 0.00000000
#> 45  Cell_45        TypeA 0.34513730
#> 46  Cell_46        TypeA 0.00000000
#> 47  Cell_47        TypeA 0.00000000
#> 48  Cell_48        TypeA 0.13974858
#> 49  Cell_49        TypeA 0.30843782
#> 50  Cell_50        TypeA 0.00000000
#> 51  Cell_51        TypeA 0.00000000
#> 52  Cell_52        TypeA 0.00000000
#> 53  Cell_53        TypeA 0.00000000
#> 54  Cell_54        TypeA 0.25235120
#> 55  Cell_55        TypeA 0.00000000
#> 56  Cell_56        TypeA 0.57646785
#> 57  Cell_57        TypeA 0.00000000
#> 58  Cell_58        TypeA 1.00000000
#> 59  Cell_59        TypeA 0.77137947
#> 60  Cell_60        TypeA 0.00000000
#> 61  Cell_61        TypeA 0.37379461
#> 62  Cell_62        TypeA 0.00000000
#> 63  Cell_63        TypeA 1.00000000
#> 64  Cell_64        TypeA 0.00000000
#> 65  Cell_65        TypeA 0.00000000
#> 66  Cell_66        TypeA 0.90733529
#> 67  Cell_67        TypeA 0.00000000
#> 68  Cell_68        TypeA 0.20823078
#> 69  Cell_69        TypeA 0.00000000
#> 70  Cell_70        TypeA 0.03375102
#> 71  Cell_71        TypeA 0.00000000
#> 72  Cell_72        TypeA 0.29721141
#> 73  Cell_73        TypeA 0.04925066
#> 74  Cell_74        TypeA 0.94301457
#> 75  Cell_75        TypeA 0.64495507
#> 76  Cell_76        TypeA 0.00000000
#> 77  Cell_77        TypeA 0.00000000
#> 78  Cell_78        TypeA 0.00000000
#> 79  Cell_79        TypeA 0.18322522
#> 80  Cell_80        TypeA 0.27768589
#> 81  Cell_01        TypeB 0.13191878
#> 82  Cell_02        TypeB 0.00000000
#> 83  Cell_03        TypeB 0.67357778
#> 84  Cell_04        TypeB 0.00000000
#> 85  Cell_05        TypeB 0.43925036
#> 86  Cell_06        TypeB 0.00000000
#> 87  Cell_07        TypeB 0.00000000
#> 88  Cell_08        TypeB 0.00000000
#> 89  Cell_09        TypeB 0.00000000
#> 90  Cell_10        TypeB 0.00000000
#> 91  Cell_11        TypeB 0.00000000
#> 92  Cell_12        TypeB 0.00000000
#> 93  Cell_13        TypeB 0.01727980
#> 94  Cell_14        TypeB 0.00000000
#> 95  Cell_15        TypeB 0.39289917
#> 96  Cell_16        TypeB 0.50297925
#> 97  Cell_17        TypeB 0.48250429
#> 98  Cell_18        TypeB 0.00000000
#> 99  Cell_19        TypeB 0.36873854
#> 100 Cell_20        TypeB 0.00000000
#> 101 Cell_21        TypeB 0.38537365
#> 102 Cell_22        TypeB 0.00000000
#> 103 Cell_23        TypeB 0.00000000
#> 104 Cell_24        TypeB 1.00000000
#> 105 Cell_25        TypeB 0.54466402
#> 106 Cell_26        TypeB 0.49974674
#> 107 Cell_27        TypeB 0.50689676
#> 108 Cell_28        TypeB 0.00000000
#> 109 Cell_29        TypeB 0.77654455
#> 110 Cell_30        TypeB 0.55864733
#> 111 Cell_31        TypeB 0.00000000
#> 112 Cell_32        TypeB 0.00000000
#> 113 Cell_33        TypeB 0.00000000
#> 114 Cell_34        TypeB 0.00000000
#> 115 Cell_35        TypeB 0.23539407
#> 116 Cell_36        TypeB 0.05085801
#> 117 Cell_37        TypeB 0.48324428
#> 118 Cell_38        TypeB 0.50198412
#> 119 Cell_39        TypeB 0.00000000
#> 120 Cell_40        TypeB 1.00000000
#> 121 Cell_41        TypeB 1.00000000
#> 122 Cell_42        TypeB 1.00000000
#> 123 Cell_43        TypeB 1.00000000
#> 124 Cell_44        TypeB 1.00000000
#> 125 Cell_45        TypeB 0.65486270
#> 126 Cell_46        TypeB 1.00000000
#> 127 Cell_47        TypeB 1.00000000
#> 128 Cell_48        TypeB 0.86025142
#> 129 Cell_49        TypeB 0.69156218
#> 130 Cell_50        TypeB 1.00000000
#> 131 Cell_51        TypeB 1.00000000
#> 132 Cell_52        TypeB 1.00000000
#> 133 Cell_53        TypeB 1.00000000
#> 134 Cell_54        TypeB 0.74764880
#> 135 Cell_55        TypeB 1.00000000
#> 136 Cell_56        TypeB 0.42353215
#> 137 Cell_57        TypeB 1.00000000
#> 138 Cell_58        TypeB 0.00000000
#> 139 Cell_59        TypeB 0.22862053
#> 140 Cell_60        TypeB 1.00000000
#> 141 Cell_61        TypeB 0.62620539
#> 142 Cell_62        TypeB 1.00000000
#> 143 Cell_63        TypeB 0.00000000
#> 144 Cell_64        TypeB 1.00000000
#> 145 Cell_65        TypeB 1.00000000
#> 146 Cell_66        TypeB 0.09266471
#> 147 Cell_67        TypeB 1.00000000
#> 148 Cell_68        TypeB 0.79176922
#> 149 Cell_69        TypeB 1.00000000
#> 150 Cell_70        TypeB 0.96624898
#> 151 Cell_71        TypeB 1.00000000
#> 152 Cell_72        TypeB 0.70278859
#> 153 Cell_73        TypeB 0.95074934
#> 154 Cell_74        TypeB 0.05698543
#> 155 Cell_75        TypeB 0.35504493
#> 156 Cell_76        TypeB 1.00000000
#> 157 Cell_77        TypeB 1.00000000
#> 158 Cell_78        TypeB 1.00000000
#> 159 Cell_79        TypeB 0.81677478
#> 160 Cell_80        TypeB 0.72231411
#> 
#> $progressions
#>    cell_id  from    to percentage
#> 1  Cell_01 TypeA TypeB 0.13191878
#> 2  Cell_02 TypeA TypeB 0.00000000
#> 3  Cell_03 TypeA TypeB 0.67357778
#> 4  Cell_04 TypeA TypeB 0.00000000
#> 5  Cell_05 TypeA TypeB 0.43925036
#> 6  Cell_06 TypeA TypeB 0.00000000
#> 7  Cell_07 TypeA TypeB 0.00000000
#> 8  Cell_08 TypeA TypeB 0.00000000
#> 9  Cell_09 TypeA TypeB 0.00000000
#> 10 Cell_10 TypeA TypeB 0.00000000
#> 11 Cell_11 TypeA TypeB 0.00000000
#> 12 Cell_12 TypeA TypeB 0.00000000
#> 13 Cell_13 TypeA TypeB 0.01727980
#> 14 Cell_14 TypeA TypeB 0.00000000
#> 15 Cell_15 TypeA TypeB 0.39289917
#> 16 Cell_16 TypeA TypeB 0.50297925
#> 17 Cell_17 TypeA TypeB 0.48250429
#> 18 Cell_18 TypeA TypeB 0.00000000
#> 19 Cell_19 TypeA TypeB 0.36873854
#> 20 Cell_20 TypeA TypeB 0.00000000
#> 21 Cell_21 TypeA TypeB 0.38537365
#> 22 Cell_22 TypeA TypeB 0.00000000
#> 23 Cell_23 TypeA TypeB 0.00000000
#> 24 Cell_24 TypeA TypeB 1.00000000
#> 25 Cell_25 TypeA TypeB 0.54466402
#> 26 Cell_26 TypeA TypeB 0.49974674
#> 27 Cell_27 TypeA TypeB 0.50689676
#> 28 Cell_28 TypeA TypeB 0.00000000
#> 29 Cell_29 TypeA TypeB 0.77654455
#> 30 Cell_30 TypeA TypeB 0.55864733
#> 31 Cell_31 TypeA TypeB 0.00000000
#> 32 Cell_32 TypeA TypeB 0.00000000
#> 33 Cell_33 TypeA TypeB 0.00000000
#> 34 Cell_34 TypeA TypeB 0.00000000
#> 35 Cell_35 TypeA TypeB 0.23539407
#> 36 Cell_36 TypeA TypeB 0.05085801
#> 37 Cell_37 TypeA TypeB 0.48324428
#> 38 Cell_38 TypeA TypeB 0.50198412
#> 39 Cell_39 TypeA TypeB 0.00000000
#> 40 Cell_40 TypeA TypeB 1.00000000
#> 41 Cell_41 TypeA TypeB 1.00000000
#> 42 Cell_42 TypeA TypeB 1.00000000
#> 43 Cell_43 TypeA TypeB 1.00000000
#> 44 Cell_44 TypeA TypeB 1.00000000
#> 45 Cell_45 TypeA TypeB 0.65486270
#> 46 Cell_46 TypeA TypeB 1.00000000
#> 47 Cell_47 TypeA TypeB 1.00000000
#> 48 Cell_48 TypeA TypeB 0.86025142
#> 49 Cell_49 TypeA TypeB 0.69156218
#> 50 Cell_50 TypeA TypeB 1.00000000
#> 51 Cell_51 TypeA TypeB 1.00000000
#> 52 Cell_52 TypeA TypeB 1.00000000
#> 53 Cell_53 TypeA TypeB 1.00000000
#> 54 Cell_54 TypeA TypeB 0.74764880
#> 55 Cell_55 TypeA TypeB 1.00000000
#> 56 Cell_56 TypeA TypeB 0.42353215
#> 57 Cell_57 TypeA TypeB 1.00000000
#> 58 Cell_58 TypeA TypeB 0.00000000
#> 59 Cell_59 TypeA TypeB 0.22862053
#> 60 Cell_60 TypeA TypeB 1.00000000
#> 61 Cell_61 TypeA TypeB 0.62620539
#> 62 Cell_62 TypeA TypeB 1.00000000
#> 63 Cell_63 TypeA TypeB 0.00000000
#> 64 Cell_64 TypeA TypeB 1.00000000
#> 65 Cell_65 TypeA TypeB 1.00000000
#> 66 Cell_66 TypeA TypeB 0.09266471
#> 67 Cell_67 TypeA TypeB 1.00000000
#> 68 Cell_68 TypeA TypeB 0.79176922
#> 69 Cell_69 TypeA TypeB 1.00000000
#> 70 Cell_70 TypeA TypeB 0.96624898
#> 71 Cell_71 TypeA TypeB 1.00000000
#> 72 Cell_72 TypeA TypeB 0.70278859
#> 73 Cell_73 TypeA TypeB 0.95074934
#> 74 Cell_74 TypeA TypeB 0.05698543
#> 75 Cell_75 TypeA TypeB 0.35504493
#> 76 Cell_76 TypeA TypeB 1.00000000
#> 77 Cell_77 TypeA TypeB 1.00000000
#> 78 Cell_78 TypeA TypeB 1.00000000
#> 79 Cell_79 TypeA TypeB 0.81677478
#> 80 Cell_80 TypeA TypeB 0.72231411
#> 
#> $trajectory_type
#> [1] "linear"
#> 
#> $directed
#> [1] FALSE
#> 
#> $pca_coords
#>                 PC1         PC2         PC3          PC4          PC5
#> Cell_01  1.26333798  2.55737125  0.49844133 -0.612744604 -2.131707651
#> Cell_02  2.07548594 -0.67754117 -0.37388152 -0.007869354  0.401937625
#> Cell_03 -2.00935863  0.31472326  0.08733576 -0.288322194 -1.320997078
#> Cell_04  1.62520668  0.55546654  1.92874025  0.313211444 -2.007344856
#> Cell_05 -2.28085949 -2.04887867 -1.02529427  0.108014514  1.294016337
#> Cell_06  0.30595498 -0.95725026 -1.68115214  2.283402245 -1.753950342
#> Cell_07 -2.00806260 -0.77294985 -1.05692738  0.969466364 -1.641765135
#> Cell_08  1.87507393  2.07878477  0.75790196 -0.578317882 -0.177283452
#> Cell_09  0.89082775  0.60972858 -0.21909732  2.129704153  0.459420248
#> Cell_10 -0.97505796  0.31935666 -1.90952288  1.111065963  0.895561425
#> Cell_11 -0.82258029 -0.09208471 -1.34655219  2.003874321 -1.410597901
#> Cell_12  1.56157466 -1.47638699 -1.24061334 -1.217769114 -2.164829847
#> Cell_13  2.80960575 -2.17257530 -1.80719993 -0.139069067  1.189757206
#> Cell_14  2.57696710  1.04137862  0.56451287  2.736648878 -0.969174694
#> Cell_15 -1.41809077  0.63716450  2.15423135  1.469577784  0.249461833
#> Cell_16 -0.85834464  2.62843617 -0.47725774  0.194132738 -1.961834898
#> Cell_17 -1.35682245  1.11416033  1.86300027  0.914044404  1.590417623
#> Cell_18  0.60245344 -0.56139066 -0.34242317 -0.884642049 -1.173812972
#> Cell_19 -0.09936044  1.05463374 -0.81890868  2.259866412 -0.059628860
#> Cell_20  1.37991827 -1.52846406 -0.25479862  1.367461279 -0.088888167
#> Cell_21  0.25367810 -2.52861490 -0.49885199  2.001756618  1.117063482
#> Cell_22 -0.62407533  2.62030129 -0.17650233 -0.978955934  0.399022372
#> Cell_23 -0.54251117  1.84189594 -2.42042762  0.094245484 -0.909344071
#> Cell_24 -1.64550886 -1.80090011  4.22841603 -0.515088180 -0.024590650
#> Cell_25  1.21801972  0.89025131 -1.22996833 -1.222266815  1.772415772
#> Cell_26  2.43460511 -2.53141781  0.16704672 -0.702149426 -0.445431279
#> Cell_27 -1.60248088  0.21006920 -1.65668231 -1.014156140  0.630519561
#> Cell_28 -0.99402176  0.76660028 -0.75603780  0.662678734 -1.554526321
#> Cell_29  1.70338700 -0.79329429 -0.62687640  0.446403561  1.518871048
#> Cell_30  0.42128293  1.34645854 -2.30398699 -1.756341244  2.284855499
#> Cell_31  0.67399661  0.89544555 -0.76124033  2.747388152  0.300418546
#> Cell_32 -3.17169613 -2.26819072  0.84474885  0.576119145 -2.372313821
#> Cell_33  2.02082974 -0.23690276 -1.75531138  1.187299867 -1.692753764
#> Cell_34 -0.29105366  1.32235074  1.69252872  1.694625396 -2.510389183
#> Cell_35 -1.37514922 -0.41302562 -1.30043646  0.415541229  1.865628368
#> Cell_36 -1.88489646  1.45774703 -0.17275580  0.016258750 -1.282020104
#> Cell_37 -0.05267390  0.13412801  0.42404379 -0.266103389 -1.847475463
#> Cell_38  2.27084725 -0.21529358  1.83093926 -0.736434967 -2.061227294
#> Cell_39  0.58160239  1.41789321 -1.81350704  0.762940993  0.130547826
#> Cell_40 -0.37490656 -1.96504004 -1.09336800 -2.144231426  1.111216776
#> Cell_41 -1.62189291 -0.96588691  0.95464956 -3.699960506 -0.037863610
#> Cell_42 -1.44815337  1.00546980  3.32699166  0.671831079  0.429591685
#> Cell_43  0.83501032  1.37717491 -0.17214298  1.126233665  1.810786974
#> Cell_44 -0.30981519  1.70789279  0.73049383 -0.167714461  1.170026352
#> Cell_45 -1.11325356  0.74531577  0.04284898  1.543545082  1.421047681
#> Cell_46 -0.20862793 -1.65065409  0.76552284 -1.387975946 -0.380506389
#> Cell_47 -0.74303493 -2.95911990 -0.24379227 -1.246196681 -1.417073315
#> Cell_48 -1.11727523  0.49195914  1.09209289  0.175526764  1.061645320
#> Cell_49  0.51679288 -2.06090262 -0.20916268 -0.532231119 -0.221844759
#> Cell_50 -1.88675189 -0.15417855 -1.63616704 -2.070905669  0.154748549
#> Cell_51 -1.44289716 -1.39069961  0.53151050 -0.396137909  0.983918633
#> Cell_52  0.64048081 -0.70458559  1.07023279  0.030496313 -0.838042414
#> Cell_53  2.27145731 -0.80891929  1.37864837  1.998234578  2.692695244
#> Cell_54 -0.87526409 -2.14076002 -0.57148145 -1.442027472  1.290399233
#> Cell_55  1.09551896  1.98184439 -0.25822617  1.027908050  2.882612108
#> Cell_56  1.79008250 -0.92773573 -0.27938087 -0.137770661 -0.912800199
#> Cell_57 -0.23879139  3.18968754  2.94729425 -0.841696603 -0.592652267
#> Cell_58 -1.01896285  0.30231810 -1.39233155  0.762074864 -1.720655188
#> Cell_59  1.40534709  0.92230005 -1.46012404 -2.401476947 -0.300095329
#> Cell_60  0.17012581  0.31486764 -1.63836076 -1.476250173  0.457401216
#> Cell_61  0.29835711  0.88156410  1.90710468  0.739287862  1.041345467
#> Cell_62  0.52491868  0.60290015  1.46311038 -0.291444000  0.632608941
#> Cell_63  0.40203530 -1.74900725  0.15626803  3.032547555 -1.379335061
#> Cell_64  1.76629304 -0.72713049  0.34048969  0.201790605  1.658926190
#> Cell_65  0.29078403  0.15278760  0.78681463  0.715589839  2.115031166
#> Cell_66 -1.50978120 -0.21660024 -1.69341643  1.728559496  1.436829755
#> Cell_67 -1.18648261 -0.83927691  1.27102190 -1.945691522  1.774922806
#> Cell_68 -0.90248053  1.05710650 -1.37105361 -0.232033040  1.735294098
#> Cell_69 -0.22985159 -0.45924109  0.09165239 -0.570732371 -1.039284421
#> Cell_70 -0.34165125 -1.02360191  1.34554373 -0.700496598 -1.574916376
#> Cell_71  3.99371998  1.11402795  0.42274452 -3.014426431  0.558503368
#> Cell_72 -0.12146840  1.14276145 -1.17926346 -2.128972323 -1.052120186
#> Cell_73 -1.29577954 -1.00751343  0.15520710  0.068519056  0.172319029
#> Cell_74  1.82901209  0.50012392 -0.80470810 -0.890888518 -1.725236859
#> Cell_75 -2.60845696 -0.36928463 -1.61270975  0.737788333  0.553690304
#> Cell_76 -1.90813194  0.36996920  2.85254207  0.256213314  1.796402490
#> Cell_77  1.35318928  1.48061773  1.52668953 -1.581419677  0.219491632
#> Cell_78  0.59937769 -3.84344031  1.34532791  0.411231844 -0.007144562
#> Cell_79 -2.59694132  2.58966810 -1.54146009 -1.914312477 -1.053057024
#> Cell_80  0.78609885 -0.70596228  1.63667379 -1.557883836  0.553145976
#>                  PC6          PC7         PC8         PC9        PC10
#> Cell_01 -1.748466612 -0.277479821 -0.53276802  0.63406196 -0.75536385
#> Cell_02  1.323714393 -0.372093051 -0.41089887  1.53053913  1.55909556
#> Cell_03 -0.729742480 -0.405268010  0.89618449  1.24084298 -1.44066123
#> Cell_04  2.648431907  0.793463069  1.20201087 -0.33182272 -2.83288713
#> Cell_05  1.321271772 -0.130006033  0.18066313 -0.58797432  1.75784711
#> Cell_06  1.406700779  2.494243047 -1.73718818  1.48091223 -0.87971843
#> Cell_07 -0.577203645  0.327860378  0.22236479  0.20796481 -1.82688140
#> Cell_08  1.092350922 -1.225013113 -0.94212321 -0.85790770 -1.11172117
#> Cell_09 -1.036894168 -0.137462282 -0.80730389 -0.50859306 -0.53017348
#> Cell_10 -1.488630176  0.072926048 -1.87856273  0.86556903  0.24108411
#> Cell_11  0.522378349 -0.504710386 -1.02163243 -0.84177993 -0.65088159
#> Cell_12  2.611135438 -1.180395775 -1.06638952  1.62135495  0.17581977
#> Cell_13  0.338736137  0.915809398  1.09200621  0.19047823 -0.77545873
#> Cell_14 -0.291326114 -0.379840328 -1.35216523 -0.46186449  1.21529123
#> Cell_15 -1.235227216 -0.204330459 -0.85947448 -0.44158386  0.47424272
#> Cell_16 -1.373564688 -2.336388166  1.35738472 -1.63390043 -1.96721562
#> Cell_17  0.075412935  2.298909795  0.69290804  1.54721772 -0.01944474
#> Cell_18 -1.209223539 -2.634370196  0.72876653 -0.56561011  1.87116233
#> Cell_19  1.852288448  1.797713074 -1.16378624  1.07999074  1.03395527
#> Cell_20  0.155328531 -0.527810374 -0.78051569 -1.92850348  0.33975362
#> Cell_21  0.068693760 -0.067421499 -0.98709247 -0.55508129  1.10012360
#> Cell_22 -1.161109192 -1.857565338  0.11468619 -1.05682454  1.47339872
#> Cell_23 -0.018864875 -1.743210873  1.53413163  1.27584692 -0.07189720
#> Cell_24  0.903087178 -0.161977269  0.36301902 -0.32667569 -1.29748826
#> Cell_25 -0.456034003  0.816577216 -1.12694440  2.50335190  1.18883630
#> Cell_26  0.577324688 -0.427848602  1.14898364  0.52056527  0.55754323
#> Cell_27 -2.387404774  0.115910333  0.58364625  0.78511976 -2.02579264
#> Cell_28  4.094575143 -1.580769027  1.11134043  1.93183265  0.91922591
#> Cell_29  0.770943248 -0.654145721  0.96141258 -1.03947630 -2.44145267
#> Cell_30 -0.371893371  0.576992787  0.88951415 -3.45831139  0.38058417
#> Cell_31 -0.976378026  1.866900508  3.58820788  0.26218205 -0.83051361
#> Cell_32  0.306981868  0.254516028 -1.21754580 -2.64645498 -0.83110291
#> Cell_33 -0.721530390  1.988291486 -2.12006793 -0.12016703  0.24492406
#> Cell_34  2.444539605 -0.677518396  0.77358594 -0.14469265  1.50574589
#> Cell_35  1.146064546  1.760775975 -0.61167924 -2.11672827  0.91533663
#> Cell_36 -0.608850426  1.968799194 -0.40575174 -1.67402402  0.40294856
#> Cell_37 -2.295485712 -1.320611586 -0.42753133 -0.84124135  1.47463628
#> Cell_38 -2.461346401  0.858200921 -1.11330100  0.45051542  0.68128968
#> Cell_39  1.040205729 -0.102786787 -1.71922622 -0.88797199  1.07458840
#> Cell_40 -0.132557584  0.948952135  0.89346929  0.96896040 -1.37894594
#> Cell_41  0.265271405  0.217337693  0.12949595 -0.37186568 -0.45580412
#> Cell_42  0.948660761  1.671842761 -0.29387939 -1.61434728  0.23794875
#> Cell_43 -0.056258324 -3.137767653  1.01679976  1.18852537 -0.51173925
#> Cell_44  1.378088165 -1.404019596 -2.86954110 -0.24309196 -0.49057147
#> Cell_45  1.480174156  1.534687144  0.71912474 -0.84025882 -1.42764076
#> Cell_46  0.631213612 -0.165230770  0.34841251 -0.91018198  1.34998827
#> Cell_47 -0.277072644 -1.866291613  0.19847343  0.92830038 -0.83050499
#> Cell_48 -0.172829530  0.360432527 -0.63540146 -0.84093612 -1.47039202
#> Cell_49 -0.201073000 -1.170051444 -0.21154768 -0.14458318 -2.32522733
#> Cell_50  0.375193450 -0.535168909  0.18888010  1.46499409  1.20341319
#> Cell_51  0.007274488  1.682265505  1.65044903  0.48551058  1.63923159
#> Cell_52 -0.615035987  0.355092410  2.11806324 -0.66621627  1.79876260
#> Cell_53 -0.965209086 -0.792613408  2.12609662 -0.06335297 -0.97479298
#> Cell_54  1.610616242 -2.500449982 -1.17706228 -0.05484682 -1.12253607
#> Cell_55 -0.428313054 -1.882935354  0.85181143 -0.48942725  0.27850803
#> Cell_56 -0.372749321 -0.885089890  0.18361787 -0.98091979  1.93171331
#> Cell_57  2.006710271 -1.654121695  1.20632814  0.87766130  1.39669214
#> Cell_58  0.657058684  0.755738177  3.39313200 -0.57016647  1.89879685
#> Cell_59  1.579904661  1.129447752  0.58569772 -1.58121482 -1.35436886
#> Cell_60  1.644455622  1.242997308  1.91350656 -1.09042533  1.34507352
#> Cell_61 -0.429508393 -0.281289978 -3.32457433  0.37990143 -0.52088206
#> Cell_62  0.470225737  1.473899021 -0.90993099 -0.17051409 -0.11147742
#> Cell_63 -1.124100210  0.003188456  1.28437822  0.51205857 -0.02139704
#> Cell_64 -1.304586046 -0.814112477  0.60891593  0.08787362  1.13577058
#> Cell_65 -0.769355970  0.502167904  0.72512642  1.75440778 -0.67861989
#> Cell_66  0.261664739 -0.624333352  0.02832781  1.42528186 -1.56377346
#> Cell_67  0.321430940  0.251628059 -0.78411404  0.91464426  0.63611889
#> Cell_68  0.229007243 -0.720678812 -2.10890030  1.08031887 -0.27888225
#> Cell_69 -2.130696065  0.509176730 -1.13737203  1.46256928  0.53938453
#> Cell_70  0.691516095  0.886444799 -0.14787150  2.58863049 -1.16166051
#> Cell_71  1.519498945  1.405522972 -1.57634005 -0.25255707 -0.27217242
#> Cell_72 -1.508721045  1.891782104  0.15916273 -0.03616256  0.24594489
#> Cell_73 -2.254418509 -0.429638679 -0.78565875 -0.86698887  0.33627447
#> Cell_74 -2.132472201 -0.113385419  0.44294193 -0.02974535 -2.08331868
#> Cell_75  0.019117807 -0.431567579 -0.28390921 -1.07171306  0.47178274
#> Cell_76 -0.746741142 -0.904278017  0.48262186  2.13616370  0.62899156
#> Cell_77 -1.389632270  3.059568032  1.45829323  1.08403707  0.35337535
#> Cell_78 -1.904955233 -0.013434252 -0.51940106  0.19988384  1.68976052
#> Cell_79 -0.998452273  0.815663040 -0.78261716  0.93219137  0.46489883
#> Cell_80  0.266665296 -0.370241813 -0.32387308 -2.68355468 -0.84650158
#>                  PC11        PC12        PC13        PC14        PC15
#> Cell_01 -1.4438539700 -2.68243728  1.01416088  2.62723921  0.05646503
#> Cell_02  0.8896202818  1.29268138  1.82251497  0.17627367  0.46085367
#> Cell_03  2.3431039489 -1.71268690  1.43202331  0.56181154  0.50430647
#> Cell_04 -0.2297583807 -0.16477451  1.67589582 -0.69931988 -1.45438522
#> Cell_05 -1.3195114746 -2.01823722  1.21916377  0.21628358  0.46780945
#> Cell_06  0.5310019135 -0.64319528 -0.90321945 -1.28985269  0.23528383
#> Cell_07 -1.3475349252  2.02877008 -0.99777481  0.50129149  0.22476460
#> Cell_08 -0.6821097359  0.08096278  1.43625459  0.21778657 -0.20262891
#> Cell_09 -1.3475059410  0.24664655  0.73796645 -1.99616308 -1.57884707
#> Cell_10 -2.3888043595  2.09395981  1.01232850  0.44603921  0.43486741
#> Cell_11  0.4109490242  0.69145990  1.57004962  0.85901769 -1.17886939
#> Cell_12  0.5587692023 -0.11581534  0.30077540 -0.44152465  0.86718300
#> Cell_13 -0.6511029242 -1.20789033  0.41092057 -0.32413941  1.87660208
#> Cell_14 -2.5699478479 -0.53223477 -0.23649381 -0.55261601 -2.03406591
#> Cell_15  1.0944903881  0.17487361  0.33580125 -0.73326909  2.73836101
#> Cell_16  0.7513052296 -0.13181648 -1.22335589  2.49362467  0.27119829
#> Cell_17 -1.3718748884  0.70000424  1.93763329 -0.70588244 -0.09859127
#> Cell_18 -0.7846226486  1.17541640  0.25381165 -1.65532853  0.96943784
#> Cell_19  2.0637880532  0.48823626 -2.62530666  0.51707391  0.11314081
#> Cell_20 -0.1707461757  0.72385146  0.31857756 -1.55767389  0.23574516
#> Cell_21  1.1620272093 -0.34752849  0.14692725  0.70397039 -0.51781770
#> Cell_22 -1.5783974899  0.09325992  1.21654698 -2.19643403  2.10872066
#> Cell_23 -0.0425854245 -0.53405231 -0.34138172 -1.91893122 -0.96499400
#> Cell_24 -1.1849862200 -1.04675526 -0.77743074 -0.94372512  1.57435605
#> Cell_25 -0.4655079458  0.64458758  1.30400632  0.45676134  1.66745994
#> Cell_26  0.8631056613  2.18365667  0.22311269 -1.18125869 -0.75472470
#> Cell_27 -0.9229388165  0.59770262  1.05860418 -0.84787566 -1.24002750
#> Cell_28 -0.6946204142  0.70879641  0.56961799  0.79534459  0.39085871
#> Cell_29  0.3452283098 -1.31892967 -1.50559690 -1.10280139  0.62810114
#> Cell_30  1.5568788958  1.22021470  1.58970456 -0.11148537 -0.04486637
#> Cell_31 -0.4879280287 -0.96512086  0.59631229 -0.34421723  0.64257220
#> Cell_32  0.2370508861 -0.88843798  0.06357409 -0.98692383  0.59030475
#> Cell_33  0.3211260556 -0.69298664 -0.42568672 -0.30372343  0.44309679
#> Cell_34 -1.4884741147  0.91031849  0.04251710  0.60847571  0.08511698
#> Cell_35 -2.3062115361  0.31746794  0.78086429  1.87267370 -1.02849886
#> Cell_36  0.3202528639  0.58326023 -0.31208267 -1.43480239 -1.78154477
#> Cell_37  0.0468065840  0.82922034 -0.12679731  0.92037154 -1.13188579
#> Cell_38  1.6474097241  1.18282253  0.24131975  0.20437286  1.21647905
#> Cell_39  2.6088186459 -0.38816934  0.23171024 -0.21865280  0.61084926
#> Cell_40 -0.9940226827 -0.14746869 -1.15825574  0.51754343 -0.51042075
#> Cell_41  1.1545650732  1.47931177 -0.51878466  0.67695534 -1.05312830
#> Cell_42  2.1877501500 -0.14674451  0.70635069  0.07544156 -0.05264995
#> Cell_43 -0.1680247983 -1.35350359 -2.02984137 -0.10927847  0.41617106
#> Cell_44  0.7686879108 -1.62699220 -1.35905212  0.04552050  0.95672282
#> Cell_45  1.5780779599  2.99502001 -1.58845338  0.16129637  1.17410443
#> Cell_46 -1.5027114779 -1.68455466 -0.07356876  0.65247364  0.06961510
#> Cell_47 -0.4261745239  0.78457999 -2.03046109 -1.88555353 -0.81661138
#> Cell_48 -1.1555357469 -0.26709109  0.14726832  0.80873446  0.30016082
#> Cell_49  0.1763944753  1.31819628  0.38562679  2.32012541  0.48822028
#> Cell_50  1.0584436605 -0.90410112  1.16920047 -0.63810374 -1.33349773
#> Cell_51  0.9089138737  0.67818647  0.34967936  1.20317057  0.87174263
#> Cell_52  0.1341646680 -0.58451987 -1.91981843  0.39881915 -1.74427004
#> Cell_53 -0.1778595672 -1.03603379 -0.94373891 -0.84154095  0.21567854
#> Cell_54  0.4212588066  1.12541774  1.77764121 -0.06876413 -0.42087532
#> Cell_55  2.4083947976  0.05941423 -0.83303696  1.83058078 -3.07885940
#> Cell_56 -0.7889481812  2.20023768 -1.05541234  1.87600694 -0.17945293
#> Cell_57 -0.6236428313  0.71212177 -0.12967644  0.07045073  0.65337218
#> Cell_58  0.2518685608 -1.56028398  0.88177333 -0.42055954 -0.20932040
#> Cell_59  0.1517670245  0.02271253  0.68662467 -0.61225083 -0.89424591
#> Cell_60 -0.0005712374 -0.18813266 -2.32656973  0.31576516  1.27463349
#> Cell_61  0.8652615641  0.35084916  0.85713441 -0.36909679 -0.52122098
#> Cell_62 -2.1296454825  0.65424741 -2.42634867  0.41831624 -0.61822765
#> Cell_63  0.8860139102 -0.39766997  0.46971667  0.04350766  0.29914782
#> Cell_64  0.5439913225  0.07027404 -0.33122312  0.46858245  1.41818813
#> Cell_65 -0.7511170520 -0.08772402  1.62393307  1.05061872 -0.85852278
#> Cell_66  0.1155035278  0.93912773 -0.12769593 -0.57861071 -0.24581652
#> Cell_67 -0.1872935007 -0.99755338 -1.33472528  0.10293941 -1.13795176
#> Cell_68 -1.3968849885 -1.21006366 -1.03059714  0.09333635  0.62442824
#> Cell_69  1.6052529776 -1.24010706 -0.32453345  0.39901515  0.62477477
#> Cell_70 -0.7252060868 -0.02640811  0.16242518  1.98454790 -0.34309840
#> Cell_71 -0.2152554288 -1.23876021 -1.04320857 -0.16576105 -0.93702467
#> Cell_72  0.1058808287 -1.00921325 -0.36518580  0.15476130  1.05906535
#> Cell_73 -1.8733288307  1.45884797 -2.26959139 -0.01767734  0.72413672
#> Cell_74  1.2627298537  0.65408750  0.24927765  0.48085462  0.91287979
#> Cell_75 -0.0298000571 -1.67044656  0.41790956  1.85000374  0.86071139
#> Cell_76  2.0953201327  0.57691513  0.01054593 -1.40713663 -1.04222383
#> Cell_77 -0.7368690224  1.34726610  0.34280234 -0.47715056  0.34545565
#> Cell_78  0.5716115052 -2.23432213  0.76467532  0.24548965 -2.12636373
#> Cell_79  0.3164408595 -0.61948406 -1.31722184 -2.21334014 -1.21090077
#> Cell_80  0.0418884384 -0.77273621  1.46684746 -1.07184367  0.64328728
#> 
#> $centroids
#>              PC1        PC2        PC3        PC4        PC5        PC6
#> TypeA  0.1039286  0.0693536 -0.3019424  0.3850317 -0.3587689  0.0854608
#> TypeB -0.1039286 -0.0693536  0.3019424 -0.3850317  0.3587689 -0.0854608
#>               PC7         PC8         PC9        PC10       PC11        PC12
#> TypeA  0.02369546 -0.09869157 -0.09824709 -0.02700419 -0.1680328  0.08574081
#> TypeB -0.02369546  0.09869157  0.09824709  0.02700419  0.1680328 -0.08574081
#>             PC13       PC14       PC15
#> TypeA  0.3227328 -0.1712661  0.1222941
#> TypeB -0.3227328  0.1712661 -0.1222941
#> 
#> $cell_milestone_geodesic
#>              TypeA      TypeB
#> Cell_01 0.20364474 1.34006830
#> Cell_02 0.00000000 1.54371305
#> Cell_03 1.03981081 0.50390224
#> Cell_04 0.00000000 1.54371305
#> Cell_05 0.67807651 0.86563654
#> Cell_06 0.00000000 1.54371305
#> Cell_07 0.00000000 1.54371305
#> Cell_08 0.00000000 1.54371305
#> Cell_09 0.00000000 1.54371305
#> Cell_10 0.00000000 1.54371305
#> Cell_11 0.00000000 1.54371305
#> Cell_12 0.00000000 1.54371305
#> Cell_13 0.02667505 1.51703800
#> Cell_14 0.00000000 1.54371305
#> Cell_15 0.60652357 0.93718948
#> Cell_16 0.77645563 0.76725741
#> Cell_17 0.74484816 0.79886489
#> Cell_18 0.00000000 1.54371305
#> Cell_19 0.56922649 0.97448656
#> Cell_20 0.00000000 1.54371305
#> Cell_21 0.59490633 0.94880672
#> Cell_22 0.00000000 1.54371305
#> Cell_23 0.00000000 1.54371305
#> Cell_24 1.54371305 0.00000000
#> Cell_25 0.84080495 0.70290810
#> Cell_26 0.77146556 0.77224748
#> Cell_27 0.78250314 0.76120991
#> Cell_28 0.00000000 1.54371305
#> Cell_29 1.19876196 0.34495109
#> Cell_30 0.86239117 0.68132188
#> Cell_31 0.00000000 1.54371305
#> Cell_32 0.00000000 1.54371305
#> Cell_33 0.00000000 1.54371305
#> Cell_34 0.00000000 1.54371305
#> Cell_35 0.36338090 1.18033215
#> Cell_36 0.07851017 1.46520287
#> Cell_37 0.74599050 0.79772255
#> Cell_38 0.77491944 0.76879361
#> Cell_39 0.00000000 1.54371305
#> Cell_40 1.54371305 0.00000000
#> Cell_41 1.54371305 0.00000000
#> Cell_42 1.54371305 0.00000000
#> Cell_43 1.54371305 0.00000000
#> Cell_44 1.54371305 0.00000000
#> Cell_45 1.01092009 0.53279296
#> Cell_46 1.54371305 0.00000000
#> Cell_47 1.54371305 0.00000000
#> Cell_48 1.32798135 0.21573170
#> Cell_49 1.06757356 0.47613949
#> Cell_50 1.54371305 0.00000000
#> Cell_51 1.54371305 0.00000000
#> Cell_52 1.54371305 0.00000000
#> Cell_53 1.54371305 0.00000000
#> Cell_54 1.15415521 0.38955784
#> Cell_55 1.54371305 0.00000000
#> Cell_56 0.65381211 0.88990094
#> Cell_57 1.54371305 0.00000000
#> Cell_58 0.00000000 1.54371305
#> Cell_59 0.35292450 1.19078855
#> Cell_60 1.54371305 0.00000000
#> Cell_61 0.96668143 0.57703162
#> Cell_62 1.54371305 0.00000000
#> Cell_63 0.00000000 1.54371305
#> Cell_64 1.54371305 0.00000000
#> Cell_65 1.54371305 0.00000000
#> Cell_66 0.14304772 1.40066533
#> Cell_67 1.54371305 0.00000000
#> Cell_68 1.22226447 0.32144858
#> Cell_69 1.54371305 0.00000000
#> Cell_70 1.49161115 0.05210189
#> Cell_71 1.54371305 0.00000000
#> Cell_72 1.08490391 0.45880914
#> Cell_73 1.46768417 0.07602888
#> Cell_74 0.08796916 1.45574389
#> Cell_75 0.54808749 0.99562556
#> Cell_76 1.54371305 0.00000000
#> Cell_77 1.54371305 0.00000000
#> Cell_78 1.54371305 0.00000000
#> Cell_79 1.26086589 0.28284716
#> Cell_80 1.11504572 0.42866733
#> 
#> $branch_assignment
#>  [1] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#>  [6] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [11] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [16] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [21] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [26] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [31] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [36] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [41] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [46] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [51] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [56] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [61] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [66] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [71] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [76] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> 
#> attr(,"class")
#> [1] "dynwrap::with_trajectory" "dynwrap::data_wrapper"   
#> [3] "list"                    
#> 
#> $trajectory$sim_trajectory
#> $id
#> [1] "20261113_053551__data_wrapper__b3tEb1BRXw"
#> 
#> $cell_ids
#>  [1] "Cell_01" "Cell_02" "Cell_03" "Cell_04" "Cell_05" "Cell_06" "Cell_07"
#>  [8] "Cell_08" "Cell_09" "Cell_10" "Cell_11" "Cell_12" "Cell_13" "Cell_14"
#> [15] "Cell_16" "Cell_17" "Cell_18" "Cell_19" "Cell_20" "Cell_21" "Cell_22"
#> [22] "Cell_23" "Cell_24" "Cell_25" "Cell_26" "Cell_28" "Cell_29" "Cell_30"
#> [29] "Cell_31" "Cell_32" "Cell_33" "Cell_34" "Cell_35" "Cell_36" "Cell_37"
#> [36] "Cell_38" "Cell_39" "Cell_40" "Cell_41" "Cell_42" "Cell_43" "Cell_44"
#> [43] "Cell_45" "Cell_46" "Cell_47" "Cell_48" "Cell_49" "Cell_50" "Cell_51"
#> [50] "Cell_52" "Cell_53" "Cell_54" "Cell_55" "Cell_56" "Cell_57" "Cell_58"
#> [57] "Cell_59" "Cell_60" "Cell_61" "Cell_62" "Cell_63" "Cell_64" "Cell_65"
#> [64] "Cell_66" "Cell_67" "Cell_68" "Cell_69" "Cell_70" "Cell_71" "Cell_72"
#> [71] "Cell_73" "Cell_74" "Cell_75" "Cell_76" "Cell_77" "Cell_78" "Cell_79"
#> [78] "Cell_80"
#> 
#> $cell_info
#> # A tibble: 78 × 1
#>    cell_id
#>    <chr>  
#>  1 Cell_01
#>  2 Cell_02
#>  3 Cell_03
#>  4 Cell_04
#>  5 Cell_05
#>  6 Cell_06
#>  7 Cell_07
#>  8 Cell_08
#>  9 Cell_09
#> 10 Cell_10
#> # ℹ 68 more rows
#> 
#> $feature_ids
#> NULL
#> 
#> $feature_info
#> NULL
#> 
#> $milestone_ids
#> [1] "TypeA" "TypeB"
#> 
#> $milestone_network
#>    from    to    length directed
#> 1 TypeA TypeB 0.7686201    FALSE
#> 
#> $divergence_regions
#> # A tibble: 0 × 3
#> # ℹ 3 variables: divergence_id <chr>, milestone_id <chr>, is_start <lgl>
#> 
#> $milestone_percentages
#>     cell_id milestone_id percentage
#> 1   Cell_01        TypeA 1.00000000
#> 2   Cell_02        TypeA 0.08138571
#> 3   Cell_03        TypeA 0.53717799
#> 4   Cell_04        TypeA 1.00000000
#> 5   Cell_05        TypeA 1.00000000
#> 6   Cell_06        TypeA 0.64823395
#> 7   Cell_07        TypeA 1.00000000
#> 8   Cell_08        TypeA 1.00000000
#> 9   Cell_09        TypeA 0.00000000
#> 10  Cell_10        TypeA 0.00000000
#> 11  Cell_11        TypeA 1.00000000
#> 12  Cell_12        TypeA 0.00000000
#> 13  Cell_13        TypeA 1.00000000
#> 14  Cell_14        TypeA 0.18854886
#> 15  Cell_16        TypeA 0.00000000
#> 16  Cell_17        TypeA 1.00000000
#> 17  Cell_18        TypeA 1.00000000
#> 18  Cell_19        TypeA 1.00000000
#> 19  Cell_20        TypeA 0.00000000
#> 20  Cell_21        TypeA 1.00000000
#> 21  Cell_22        TypeA 0.00000000
#> 22  Cell_23        TypeA 1.00000000
#> 23  Cell_24        TypeA 1.00000000
#> 24  Cell_25        TypeA 0.00000000
#> 25  Cell_26        TypeA 1.00000000
#> 26  Cell_28        TypeA 1.00000000
#> 27  Cell_29        TypeA 1.00000000
#> 28  Cell_30        TypeA 0.98469417
#> 29  Cell_31        TypeA 1.00000000
#> 30  Cell_32        TypeA 1.00000000
#> 31  Cell_33        TypeA 0.11171760
#> 32  Cell_34        TypeA 1.00000000
#> 33  Cell_35        TypeA 1.00000000
#> 34  Cell_36        TypeA 0.01245950
#> 35  Cell_37        TypeA 0.73523955
#> 36  Cell_38        TypeA 1.00000000
#> 37  Cell_39        TypeA 1.00000000
#> 38  Cell_40        TypeA 0.00000000
#> 39  Cell_41        TypeA 1.00000000
#> 40  Cell_42        TypeA 0.00000000
#> 41  Cell_43        TypeA 0.00000000
#> 42  Cell_44        TypeA 0.37951751
#> 43  Cell_45        TypeA 0.00000000
#> 44  Cell_46        TypeA 0.00000000
#> 45  Cell_47        TypeA 0.00000000
#> 46  Cell_48        TypeA 0.68366976
#> 47  Cell_49        TypeA 0.00000000
#> 48  Cell_50        TypeA 1.00000000
#> 49  Cell_51        TypeA 0.00000000
#> 50  Cell_52        TypeA 0.00000000
#> 51  Cell_53        TypeA 0.55601318
#> 52  Cell_54        TypeA 0.00000000
#> 53  Cell_55        TypeA 0.00000000
#> 54  Cell_56        TypeA 0.00000000
#> 55  Cell_57        TypeA 0.00000000
#> 56  Cell_58        TypeA 0.96840390
#> 57  Cell_59        TypeA 1.00000000
#> 58  Cell_60        TypeA 0.00000000
#> 59  Cell_61        TypeA 0.00000000
#> 60  Cell_62        TypeA 0.00000000
#> 61  Cell_63        TypeA 0.03673746
#> 62  Cell_64        TypeA 0.08964335
#> 63  Cell_65        TypeA 1.00000000
#> 64  Cell_66        TypeA 0.00000000
#> 65  Cell_67        TypeA 0.00000000
#> 66  Cell_68        TypeA 0.28163256
#> 67  Cell_69        TypeA 0.00000000
#> 68  Cell_70        TypeA 1.00000000
#> 69  Cell_71        TypeA 0.45798690
#> 70  Cell_72        TypeA 0.00000000
#> 71  Cell_73        TypeA 0.04366513
#> 72  Cell_74        TypeA 0.68120379
#> 73  Cell_75        TypeA 0.00000000
#> 74  Cell_76        TypeA 0.93135170
#> 75  Cell_77        TypeA 0.45938046
#> 76  Cell_78        TypeA 1.00000000
#> 77  Cell_79        TypeA 0.00000000
#> 78  Cell_80        TypeA 0.65596735
#> 79  Cell_01        TypeB 0.00000000
#> 80  Cell_02        TypeB 0.91861429
#> 81  Cell_03        TypeB 0.46282201
#> 82  Cell_04        TypeB 0.00000000
#> 83  Cell_05        TypeB 0.00000000
#> 84  Cell_06        TypeB 0.35176605
#> 85  Cell_07        TypeB 0.00000000
#> 86  Cell_08        TypeB 0.00000000
#> 87  Cell_09        TypeB 1.00000000
#> 88  Cell_10        TypeB 1.00000000
#> 89  Cell_11        TypeB 0.00000000
#> 90  Cell_12        TypeB 1.00000000
#> 91  Cell_13        TypeB 0.00000000
#> 92  Cell_14        TypeB 0.81145114
#> 93  Cell_16        TypeB 1.00000000
#> 94  Cell_17        TypeB 0.00000000
#> 95  Cell_18        TypeB 0.00000000
#> 96  Cell_19        TypeB 0.00000000
#> 97  Cell_20        TypeB 1.00000000
#> 98  Cell_21        TypeB 0.00000000
#> 99  Cell_22        TypeB 1.00000000
#> 100 Cell_23        TypeB 0.00000000
#> 101 Cell_24        TypeB 0.00000000
#> 102 Cell_25        TypeB 1.00000000
#> 103 Cell_26        TypeB 0.00000000
#> 104 Cell_28        TypeB 0.00000000
#> 105 Cell_29        TypeB 0.00000000
#> 106 Cell_30        TypeB 0.01530583
#> 107 Cell_31        TypeB 0.00000000
#> 108 Cell_32        TypeB 0.00000000
#> 109 Cell_33        TypeB 0.88828240
#> 110 Cell_34        TypeB 0.00000000
#> 111 Cell_35        TypeB 0.00000000
#> 112 Cell_36        TypeB 0.98754050
#> 113 Cell_37        TypeB 0.26476045
#> 114 Cell_38        TypeB 0.00000000
#> 115 Cell_39        TypeB 0.00000000
#> 116 Cell_40        TypeB 1.00000000
#> 117 Cell_41        TypeB 0.00000000
#> 118 Cell_42        TypeB 1.00000000
#> 119 Cell_43        TypeB 1.00000000
#> 120 Cell_44        TypeB 0.62048249
#> 121 Cell_45        TypeB 1.00000000
#> 122 Cell_46        TypeB 1.00000000
#> 123 Cell_47        TypeB 1.00000000
#> 124 Cell_48        TypeB 0.31633024
#> 125 Cell_49        TypeB 1.00000000
#> 126 Cell_50        TypeB 0.00000000
#> 127 Cell_51        TypeB 1.00000000
#> 128 Cell_52        TypeB 1.00000000
#> 129 Cell_53        TypeB 0.44398682
#> 130 Cell_54        TypeB 1.00000000
#> 131 Cell_55        TypeB 1.00000000
#> 132 Cell_56        TypeB 1.00000000
#> 133 Cell_57        TypeB 1.00000000
#> 134 Cell_58        TypeB 0.03159610
#> 135 Cell_59        TypeB 0.00000000
#> 136 Cell_60        TypeB 1.00000000
#> 137 Cell_61        TypeB 1.00000000
#> 138 Cell_62        TypeB 1.00000000
#> 139 Cell_63        TypeB 0.96326254
#> 140 Cell_64        TypeB 0.91035665
#> 141 Cell_65        TypeB 0.00000000
#> 142 Cell_66        TypeB 1.00000000
#> 143 Cell_67        TypeB 1.00000000
#> 144 Cell_68        TypeB 0.71836744
#> 145 Cell_69        TypeB 1.00000000
#> 146 Cell_70        TypeB 0.00000000
#> 147 Cell_71        TypeB 0.54201310
#> 148 Cell_72        TypeB 1.00000000
#> 149 Cell_73        TypeB 0.95633487
#> 150 Cell_74        TypeB 0.31879621
#> 151 Cell_75        TypeB 1.00000000
#> 152 Cell_76        TypeB 0.06864830
#> 153 Cell_77        TypeB 0.54061954
#> 154 Cell_78        TypeB 0.00000000
#> 155 Cell_79        TypeB 1.00000000
#> 156 Cell_80        TypeB 0.34403265
#> 
#> $progressions
#>    cell_id  from    to percentage
#> 1  Cell_01 TypeA TypeB 0.00000000
#> 2  Cell_02 TypeA TypeB 0.91861429
#> 3  Cell_03 TypeA TypeB 0.46282201
#> 4  Cell_04 TypeA TypeB 0.00000000
#> 5  Cell_05 TypeA TypeB 0.00000000
#> 6  Cell_06 TypeA TypeB 0.35176605
#> 7  Cell_07 TypeA TypeB 0.00000000
#> 8  Cell_08 TypeA TypeB 0.00000000
#> 9  Cell_09 TypeA TypeB 1.00000000
#> 10 Cell_10 TypeA TypeB 1.00000000
#> 11 Cell_11 TypeA TypeB 0.00000000
#> 12 Cell_12 TypeA TypeB 1.00000000
#> 13 Cell_13 TypeA TypeB 0.00000000
#> 14 Cell_14 TypeA TypeB 0.81145114
#> 15 Cell_16 TypeA TypeB 1.00000000
#> 16 Cell_17 TypeA TypeB 0.00000000
#> 17 Cell_18 TypeA TypeB 0.00000000
#> 18 Cell_19 TypeA TypeB 0.00000000
#> 19 Cell_20 TypeA TypeB 1.00000000
#> 20 Cell_21 TypeA TypeB 0.00000000
#> 21 Cell_22 TypeA TypeB 1.00000000
#> 22 Cell_23 TypeA TypeB 0.00000000
#> 23 Cell_24 TypeA TypeB 0.00000000
#> 24 Cell_25 TypeA TypeB 1.00000000
#> 25 Cell_26 TypeA TypeB 0.00000000
#> 26 Cell_28 TypeA TypeB 0.00000000
#> 27 Cell_29 TypeA TypeB 0.00000000
#> 28 Cell_30 TypeA TypeB 0.01530583
#> 29 Cell_31 TypeA TypeB 0.00000000
#> 30 Cell_32 TypeA TypeB 0.00000000
#> 31 Cell_33 TypeA TypeB 0.88828240
#> 32 Cell_34 TypeA TypeB 0.00000000
#> 33 Cell_35 TypeA TypeB 0.00000000
#> 34 Cell_36 TypeA TypeB 0.98754050
#> 35 Cell_37 TypeA TypeB 0.26476045
#> 36 Cell_38 TypeA TypeB 0.00000000
#> 37 Cell_39 TypeA TypeB 0.00000000
#> 38 Cell_40 TypeA TypeB 1.00000000
#> 39 Cell_41 TypeA TypeB 0.00000000
#> 40 Cell_42 TypeA TypeB 1.00000000
#> 41 Cell_43 TypeA TypeB 1.00000000
#> 42 Cell_44 TypeA TypeB 0.62048249
#> 43 Cell_45 TypeA TypeB 1.00000000
#> 44 Cell_46 TypeA TypeB 1.00000000
#> 45 Cell_47 TypeA TypeB 1.00000000
#> 46 Cell_48 TypeA TypeB 0.31633024
#> 47 Cell_49 TypeA TypeB 1.00000000
#> 48 Cell_50 TypeA TypeB 0.00000000
#> 49 Cell_51 TypeA TypeB 1.00000000
#> 50 Cell_52 TypeA TypeB 1.00000000
#> 51 Cell_53 TypeA TypeB 0.44398682
#> 52 Cell_54 TypeA TypeB 1.00000000
#> 53 Cell_55 TypeA TypeB 1.00000000
#> 54 Cell_56 TypeA TypeB 1.00000000
#> 55 Cell_57 TypeA TypeB 1.00000000
#> 56 Cell_58 TypeA TypeB 0.03159610
#> 57 Cell_59 TypeA TypeB 0.00000000
#> 58 Cell_60 TypeA TypeB 1.00000000
#> 59 Cell_61 TypeA TypeB 1.00000000
#> 60 Cell_62 TypeA TypeB 1.00000000
#> 61 Cell_63 TypeA TypeB 0.96326254
#> 62 Cell_64 TypeA TypeB 0.91035665
#> 63 Cell_65 TypeA TypeB 0.00000000
#> 64 Cell_66 TypeA TypeB 1.00000000
#> 65 Cell_67 TypeA TypeB 1.00000000
#> 66 Cell_68 TypeA TypeB 0.71836744
#> 67 Cell_69 TypeA TypeB 1.00000000
#> 68 Cell_70 TypeA TypeB 0.00000000
#> 69 Cell_71 TypeA TypeB 0.54201310
#> 70 Cell_72 TypeA TypeB 1.00000000
#> 71 Cell_73 TypeA TypeB 0.95633487
#> 72 Cell_74 TypeA TypeB 0.31879621
#> 73 Cell_75 TypeA TypeB 1.00000000
#> 74 Cell_76 TypeA TypeB 0.06864830
#> 75 Cell_77 TypeA TypeB 0.54061954
#> 76 Cell_78 TypeA TypeB 0.00000000
#> 77 Cell_79 TypeA TypeB 1.00000000
#> 78 Cell_80 TypeA TypeB 0.34403265
#> 
#> $trajectory_type
#> [1] "linear"
#> 
#> $directed
#> [1] FALSE
#> 
#> $pca_coords
#>                   PC1         PC2          PC3         PC4          PC5
#> Cell_01  0.2172365628  0.12246973 -0.400246147  2.56407681  1.847211623
#> Cell_02  0.9326033245 -0.46316506 -0.148794921 -0.74014695  1.484271267
#> Cell_03  1.0226666847  0.25362931 -0.592907505  0.41850392 -0.434076901
#> Cell_04 -1.2678144970  1.17264906 -0.182050433  2.35879804  0.298434472
#> Cell_05  1.6301382876 -1.53939538  0.376108461  1.54192429  0.578242272
#> Cell_06 -1.0513578799  0.92676906  0.533410532  0.85638363  0.061190050
#> Cell_07  0.0378690058  0.88240147 -0.402596427  1.23840458 -0.632304974
#> Cell_08 -1.1509935725 -0.28202818 -0.823858530  0.77329373 -0.791394852
#> Cell_09  0.7690456228 -0.10014474  0.427862531 -1.03413694  1.604288309
#> Cell_10 -0.1101049338  1.00801999 -0.178267108 -1.20149889  0.167539243
#> Cell_11 -0.9908253329  1.71165999 -0.258867684  1.65428960  0.662918902
#> Cell_12  0.5854889818 -0.39201252  0.037874288  0.20670999  1.190776090
#> Cell_13  0.5900033442  0.84010142  0.415769355  1.67782625  0.064509181
#> Cell_14  0.2897674458 -0.31823342  1.021142962 -0.24981514 -0.416658831
#> Cell_16 -1.2410528087  1.18137318  0.979851459 -0.51472158 -1.098841558
#> Cell_17  0.4466767608 -2.82630407 -0.068716542  1.16652917 -0.380845558
#> Cell_18 -1.1694692312  0.72529907 -0.832243848  1.39452198 -0.616893719
#> Cell_19 -0.7458663391 -0.22651039  0.002747098 -0.24695826 -0.545723131
#> Cell_20 -0.7884055223 -1.75682178 -0.301806436  0.37465690  0.807071768
#> Cell_21 -0.6589055863 -1.02312458  0.102443639  1.16716842 -0.983018557
#> Cell_22 -0.3460336658 -1.38011031  0.922875319 -0.80367061  0.098703556
#> Cell_23  1.1810528304 -2.01057464  1.354340436  1.41419386 -0.616756713
#> Cell_24 -0.4351721445 -0.16060879  0.354980354  1.95421471 -0.005079907
#> Cell_25 -0.6591593476  0.08305122  0.155884103 -1.94464609 -0.717805240
#> Cell_26  2.4046938092  0.08874302  1.113758416  0.06558808 -0.936578480
#> Cell_28 -0.6134344690 -0.83465646  0.300495582  1.55774451 -0.021911729
#> Cell_29 -1.7403820214  0.24128317 -1.120936513  0.13674282 -0.729625628
#> Cell_30 -1.5587855265  0.06661074 -1.356370740 -0.62748611 -0.581444213
#> Cell_31 -0.8286243360 -0.40315486  0.307485771  1.07241138 -0.365635264
#> Cell_32  0.0952947299 -0.11310789 -0.594720541  1.08935977 -0.380529121
#> Cell_33  0.5458102048  0.90931430  1.129603832  0.69112491  0.622055536
#> Cell_34  1.6325472229 -0.18911110  0.775073794  1.13078435 -0.175947576
#> Cell_35  0.1084072862 -0.77236671 -0.163216786 -0.47386589 -0.352959048
#> Cell_36  0.1882086890  1.19399814 -0.032970163  0.16754496  0.175933484
#> Cell_37  1.2514244417  1.54746467 -0.142992347 -0.13147558  0.499864123
#> Cell_38  0.2302934107  1.02249541 -0.006082775 -0.83729410 -0.027754869
#> Cell_39 -0.4728037398 -1.20485317 -0.070640231  0.13639810 -0.018825087
#> Cell_40  1.1512506061  0.53936862  1.803599525  0.59608958  2.301265686
#> Cell_41 -0.0734849912 -1.26275028 -1.332272440  0.75622176 -0.565474380
#> Cell_42 -1.5066432841  1.17714667 -0.106058879 -1.63727492  0.806727510
#> Cell_43  1.0741839302 -0.48614399  0.443116428  0.06036634  0.066575089
#> Cell_44  1.1359431831 -0.86430367 -0.017109432 -0.34048994 -0.481870249
#> Cell_45  0.9802556732 -0.44635287  0.781283381 -0.61949555  1.558692322
#> Cell_46 -0.3899931829  0.25640059 -0.043548389  0.17053180  1.395975363
#> Cell_47 -0.1365681120  0.96558070  0.284081849 -0.34785277  0.882875271
#> Cell_48 -1.8077357615 -0.62541103  0.063580358  0.35651176  0.941108976
#> Cell_49 -0.0006331437 -0.46648834  0.829303484 -0.88888873  1.299045310
#> Cell_50  0.3908456605  0.21328039  0.580708952  1.50636669 -0.172101194
#> Cell_51 -1.0996032670  0.19131474  1.070091510  0.44017236 -0.349859329
#> Cell_52  1.4343246920  0.22831489  0.756235131 -1.74423277  0.177060423
#> Cell_53 -0.1357570423  0.76325523 -0.087272020 -0.84916237 -0.708430680
#> Cell_54  0.6210847826 -1.92821596  0.951487016 -0.05110422 -0.282044071
#> Cell_55 -0.3599644383  0.18207484  1.118774917  0.16323050 -0.815821761
#> Cell_56  0.7481913114  0.55681791  1.236885192  0.33540616  1.718502555
#> Cell_57  0.5074951895 -0.42285299 -0.118029111  0.25995386  0.692695485
#> Cell_58 -1.1199202510  0.33157068  0.373447946  0.31936358  0.603860587
#> Cell_59  0.0807634557 -0.46263309 -0.139302829  0.42708186  0.056665618
#> Cell_60 -1.3947131922 -0.33154981 -0.882863254 -0.51906024  2.113883836
#> Cell_61  0.6001302527 -0.97651016  0.686893009  0.08109355  0.039966374
#> Cell_62 -0.6401115307  1.24083660  0.477188156 -0.71228631  0.538003557
#> Cell_63  0.7623439488  0.13919456  0.582677165  0.72000268 -0.081140312
#> Cell_64 -0.2318004442  2.23365350  0.207620878  0.91772465 -0.584059734
#> Cell_65  0.5037399392  0.03574757 -0.283640538  1.73206523 -1.262531335
#> Cell_66 -0.4972703621 -0.64094194 -1.310099518  0.09599443  0.735975793
#> Cell_67  0.0557837119  0.95083653 -0.197282406 -0.36133038  0.878352825
#> Cell_68  0.1553627119 -1.23340343 -1.011140797  0.35840273  0.159040739
#> Cell_69 -0.8699573963 -0.06874587  1.263522850  0.94358565  0.930055593
#> Cell_70  1.3267081872 -0.35383376  0.795616532  0.22401658 -1.445939777
#> Cell_71  0.1087469856  0.21559947  0.678483064 -0.05326305  0.180018739
#> Cell_72  0.6521594606  1.25107269  1.006738802 -0.77063595 -0.558670557
#> Cell_73 -0.5629671317  2.16078833  0.552522518 -0.88450926 -0.145610378
#> Cell_74  0.8205898928 -1.31040165  0.291452354  0.37960710  0.236828417
#> Cell_75 -0.3122614670 -0.20818423  0.785004141 -0.21848645  1.661034667
#> Cell_76 -0.2137395594 -0.75910973 -0.263794631  1.39981733  1.163806542
#> Cell_77  0.3139479927  0.80869838  1.151149708  0.98910130  0.516926152
#> Cell_78 -0.9426449634 -0.26310623  0.947979298  0.10793558  0.563667164
#> Cell_79 -0.3868297396  0.28748134  1.466726389  0.72951548 -0.507028531
#> Cell_80  0.2450432194 -0.86093240  1.149611784  1.79838991 -0.116373608
#>                 PC6          PC7         PC8         PC9        PC10
#> Cell_01 -0.14735202 -0.803367320 -0.30600552 -0.79652258 -1.76042181
#> Cell_02  0.48262679  1.904596690 -1.14772923  0.89247694 -0.68443042
#> Cell_03  0.06084429 -0.033253043  0.56569870  0.88968180  0.81226906
#> Cell_04  0.27753146  0.854101672 -0.05502387 -0.54767235  1.15718900
#> Cell_05 -0.91526150 -1.747240529 -0.10784397  0.51537263  0.26786175
#> Cell_06 -0.55786381 -0.399591188 -0.49055900 -0.57026624  0.08559347
#> Cell_07  0.33556282  0.165569264  0.30730020 -0.98960150  0.49884379
#> Cell_08 -0.68932089 -0.485923065 -1.08614205 -0.74162673  0.52978562
#> Cell_09 -0.75697318  0.378374196 -0.32122001 -0.53674360  0.28270377
#> Cell_10  0.45886032  0.358383451 -0.95883913  0.59041839 -0.48247733
#> Cell_11 -0.79675584 -0.431534934 -0.41316303  0.67163518  0.49853107
#> Cell_12  0.61731347  0.211454468  0.39446764  0.96209376  0.64400237
#> Cell_13  0.30709149 -1.338677329  0.06213987 -0.27783585  1.50263077
#> Cell_14 -1.70650499 -0.316813167 -1.33188873 -1.83687347  1.25043178
#> Cell_16 -0.57273810  0.062600790  0.44700580  0.93089319  0.04336814
#> Cell_17  0.64748861 -0.370512460  1.20521748  0.26515278 -0.48905766
#> Cell_18  1.27855582 -0.884480161 -0.50597455  0.13906295  0.86844701
#> Cell_19 -0.23497096 -0.001709605  0.39376276 -0.69858009 -0.28988311
#> Cell_20  1.37924328  0.707997775 -0.05331507  0.18885128 -0.35315731
#> Cell_21  0.39813818  1.300323628 -2.37192895  0.78758063 -0.33564357
#> Cell_22 -0.06004185 -0.134748949  1.15172448  0.05850950 -0.50405254
#> Cell_23  0.32420263  0.529650038 -1.05500577  1.18911019  0.18319956
#> Cell_24  1.44556533 -0.273251702  0.15195168  0.18915343  1.01876185
#> Cell_25 -1.87419540  0.137382899  0.11175319 -0.50272720  1.00104014
#> Cell_26 -1.12856933  2.285104091  0.12743981  0.28935770  0.17954009
#> Cell_28  0.33047051 -1.200827552 -0.49551804 -0.82179389  1.04867253
#> Cell_29  1.14470265  0.075491552  0.02302924  0.21507439  0.62166368
#> Cell_30 -0.10349718 -0.670765661  0.36463602  0.55185732 -1.33179382
#> Cell_31 -0.38243943 -0.177720738 -0.84196092 -0.55775136 -0.94060551
#> Cell_32 -0.79384497 -0.729561132 -0.03402734  0.77026956  1.17846697
#> Cell_33  1.28341522  0.898406948 -1.53103121 -1.11509328  0.69835690
#> Cell_34 -0.36962812  0.272009735  0.51637994  0.43719972  0.63047080
#> Cell_35 -0.20664273 -0.413539155  0.76050996  0.39712021 -0.18778335
#> Cell_36  0.62046694 -0.251466466 -0.04905605 -0.69069335  1.39818432
#> Cell_37  0.57903808 -0.025046350  0.01755450 -3.05606464 -0.18004381
#> Cell_38  0.52907189  0.668913336 -0.29984024  0.79923046 -0.65076855
#> Cell_39 -1.29600993 -0.797426203  0.26251823 -1.19083447 -0.38953891
#> Cell_40  0.40710840  0.902619484 -1.58476944  0.60301068  0.60164730
#> Cell_41 -0.70326056 -0.062644265 -0.87749296  0.69893110 -0.73321172
#> Cell_42 -0.37945739 -0.971856926  0.54314757 -1.11339429 -0.37240021
#> Cell_43  0.09341664 -1.139877343 -0.60220526  1.42473584 -0.48128103
#> Cell_44  0.24362371 -0.434908089 -0.80875511  0.35304658  0.89921906
#> Cell_45  0.77107474  1.160464272  0.36379858 -0.38859631  0.21315058
#> Cell_46  0.07835793 -0.397572253  0.86787794  0.85576176  0.91676268
#> Cell_47  0.15521209 -0.075305561 -1.38661742  0.35682922  0.58515498
#> Cell_48 -0.37763087  0.658718621 -0.68619582 -2.08361169 -0.63211912
#> Cell_49  0.40732384  0.197635941  0.37997336 -0.83491488  0.17096820
#> Cell_50 -0.24245055 -1.024132329  1.15167893  1.43304759  0.40333715
#> Cell_51 -0.45920806 -1.023832356 -0.12067672  0.37702342 -0.86720292
#> Cell_52 -1.54703159 -0.634214462 -1.41374290  0.83675683  0.09762355
#> Cell_53  0.48177603  0.551912513 -0.53187779 -0.37463595  0.26355237
#> Cell_54 -1.81560144 -0.652464987  0.51316854  1.12857435  0.65179186
#> Cell_55 -0.68712326  0.011634604 -1.41727817 -0.65029605  0.10857764
#> Cell_56  0.42188258 -0.281681439 -0.52912509  1.12038054 -0.09596401
#> Cell_57  0.26766643 -0.225769102 -0.88526671  0.21453217  1.27434463
#> Cell_58 -0.28922843  1.980332741 -1.60650489  0.51665125  1.25026549
#> Cell_59  1.17879843  0.203829949 -1.28268660 -1.27704060 -0.41037990
#> Cell_60  0.65756170 -0.348875347  0.91586227 -0.64221527  0.23806974
#> Cell_61  0.18488922  0.640931417  0.08523879  1.02019538  0.13203000
#> Cell_62  0.51644615  0.165896885  0.44118067  0.54975216  0.41082015
#> Cell_63 -0.01944476  0.184213112  0.63329220  0.17093542  1.88716199
#> Cell_64 -0.29486020 -0.816106235 -0.57575207 -0.94757737  1.32821853
#> Cell_65 -0.73604268 -1.087628847 -1.05764792 -0.42196288  0.29790710
#> Cell_66 -0.61221590 -0.082868444 -0.80039319  0.85887580  0.74674767
#> Cell_67  0.69111932 -0.036110358 -1.08529752  0.87504736  1.54945893
#> Cell_68  0.98874053 -0.115479755 -0.38562324 -1.73777405 -1.38350386
#> Cell_69 -0.09559517 -0.270517131 -2.30720248  0.12827484 -1.13341703
#> Cell_70 -0.04032409  0.275462918 -0.70199472 -0.00401383 -0.06055673
#> Cell_71  0.90167595 -0.488811535 -0.05094759 -0.15292164  0.33264947
#> Cell_72  0.60962804  1.669751109 -1.27290758  0.69363238  0.40601542
#> Cell_73 -0.32822327  0.991284034  1.01291460 -0.21355429 -0.41671802
#> Cell_74  1.13603788  0.284084021  0.72820919 -1.01752837  1.96770150
#> Cell_75 -0.20675967  1.309549565 -0.95394259  0.85912960  1.75660324
#> Cell_76  0.12829818 -0.404809809  0.30981569  0.08233848 -0.38406837
#> Cell_77 -0.88255280 -0.448603826  0.41714188  0.28318003  0.76011372
#> Cell_78  0.06108564  0.923434244  0.14622253 -0.39955342  0.47335868
#> Cell_79 -0.69899460 -0.252863811  0.56302608 -0.94236027 -0.32600517
#> Cell_80  0.70078506 -1.817548651 -1.37687340 -0.87638169  0.05459009
#>                PC11        PC12        PC13        PC14         PC15
#> Cell_01  0.28856562  0.06762215  0.93955542 -0.92928828  0.869029479
#> Cell_02  0.03585864 -0.72138349  0.75112706  0.28758565  1.633183349
#> Cell_03 -0.61239649 -2.64478850  0.62875131  0.09463638  0.599526324
#> Cell_04  0.88519768  0.73086756  0.60087791 -0.84443708 -0.059967919
#> Cell_05 -0.01974335  0.86642724 -0.15316792 -1.40507678  0.513880036
#> Cell_06 -0.92909539  0.72159276 -1.60832526 -1.18005192 -0.045957244
#> Cell_07 -0.09291311  0.67649373  1.56640187  0.20476073  0.280988888
#> Cell_08 -1.23892933  0.37778034  0.60816361 -0.15527913 -0.438935247
#> Cell_09  0.08023453 -1.19343961  1.66510325 -0.47520952 -1.391410281
#> Cell_10  0.25469409 -0.28251058 -0.11569052 -0.37745128  0.684121771
#> Cell_11  0.26022537 -1.72104449  1.74432655 -0.78903563  0.219455973
#> Cell_12 -1.03603964  0.59240744 -0.65835167 -0.52668553 -0.188368781
#> Cell_13  0.44719988 -0.09059321  0.89435727 -0.50891078 -0.945517932
#> Cell_14 -1.38160107  0.80451516  0.93804606  0.11553531 -0.414047358
#> Cell_16 -0.46146627 -0.55198847  0.12556306  0.40655327 -0.616236305
#> Cell_17  0.19905253  0.27698783  1.28674537  0.94750300  1.100073683
#> Cell_18 -1.17036354 -0.37445748 -0.83288346 -0.40903419  0.518944313
#> Cell_19  0.75505176  1.12863719  1.58904392 -0.56378799 -0.707101045
#> Cell_20 -0.40076811  0.41513852 -1.56984880  1.04095588 -0.308183128
#> Cell_21 -1.29131521  0.05648395 -0.08676888 -1.38189546 -0.011345280
#> Cell_22  0.26671067 -0.70770013 -0.09662983  0.54921409 -0.651129162
#> Cell_23 -0.27527779  0.49341038 -0.16866345 -0.69419753  0.269112231
#> Cell_24 -0.41047799  0.55498568 -0.85334248  0.52763845  0.557305600
#> Cell_25  0.43514973  0.01073579  1.05422819  0.53188666  0.078680641
#> Cell_26 -0.66489178  1.30882670  1.46338772 -0.13411697  0.411802619
#> Cell_28 -0.91179481  1.44093469 -0.23059616 -0.05233258  0.511077763
#> Cell_29 -0.12515699  0.51922278 -0.36043923 -0.45249331  0.317808808
#> Cell_30  3.05265793  0.43876439 -0.17549626  0.94164694 -0.197201330
#> Cell_31 -0.02543208  0.34534002  0.18017639 -1.68010485  1.330668511
#> Cell_32 -0.83908574  1.68923985  0.08646673  0.35514686  0.715830935
#> Cell_33 -0.04856511  0.51107503  0.11060696  1.21522129  1.324696245
#> Cell_34  0.58131187  1.31781869 -0.46334146  0.05140398 -0.061260742
#> Cell_35 -0.70314219  2.36231918  1.33654810  0.12635648 -1.314697369
#> Cell_36  0.02165396  0.37979218 -0.01694713 -0.20707906 -1.260455056
#> Cell_37  0.93204192  1.03236401 -0.02503731  1.31325577  0.747111792
#> Cell_38 -0.09439241  0.35612099  0.73105476 -0.51215027  1.054838463
#> Cell_39  0.79570170 -0.35594625  1.57851156  0.62462347  0.359120830
#> Cell_40  0.91091281 -0.17710629 -1.42232874  0.34937700 -1.519520822
#> Cell_41  0.79029030  0.75796277 -0.83393897  0.37973035  1.361556943
#> Cell_42  1.14791794 -0.45756533 -0.29689689 -0.14115753 -0.349090340
#> Cell_43 -0.03856234  0.06540419 -1.78005938 -0.57204618 -0.743354713
#> Cell_44  1.06390591  2.58191906 -0.11768544  0.50501310 -1.186715641
#> Cell_45 -0.64343022 -0.30432940  1.17810773 -0.57395488 -0.226512650
#> Cell_46 -0.15127880 -1.32562955  0.04650710  0.30216439  0.067403493
#> Cell_47 -0.74947890 -0.73674266 -0.14171305 -0.25937698 -0.493122844
#> Cell_48  0.74838372  0.23656859  0.34864918  1.13256307  0.599573269
#> Cell_49 -0.45084449  0.56757383 -1.05692967 -0.78533557 -1.186169421
#> Cell_50  0.18284223  1.16108633  0.32713170  0.68440633  0.438077320
#> Cell_51 -0.87784745 -0.30999401  0.05411039  2.22928645  0.221539188
#> Cell_52  1.00675077 -0.53264333  0.09481815 -0.06989180 -0.958945202
#> Cell_53 -0.02788241  0.40471695  0.87728993  0.15811585 -0.315000279
#> Cell_54 -0.65719466 -0.59239784  0.98812955  1.20753619 -0.005185379
#> Cell_55 -0.90751227 -1.29206417 -0.90896385  0.10956212  1.370921217
#> Cell_56  0.37451710 -0.02282563  0.57791738  0.12283989  0.099103177
#> Cell_57 -0.43388460  0.03155908 -0.53156838 -0.51453175 -0.700646033
#> Cell_58 -0.71899292  0.89346052  0.35313444 -1.42416496 -0.125635805
#> Cell_59 -0.03965631 -0.26496510  1.32506723  0.16513574  0.988478265
#> Cell_60  0.35582203 -0.27691131  0.42907565  1.57797373 -0.300910694
#> Cell_61 -0.81891839  1.33944799 -1.27783901  0.39044212 -0.454902635
#> Cell_62  0.05669897 -1.25535383 -0.10402778  0.92937860  1.087491958
#> Cell_63 -0.53102183  0.83380813  0.32746421  0.93866594 -0.921045624
#> Cell_64  0.54058050 -1.01821311 -0.60992989  0.20901282 -0.139139272
#> Cell_65 -0.20461886  0.40688168  0.71026151  0.10981046 -0.029346566
#> Cell_66  1.43738823  0.15914288 -1.27992331  0.50409316 -1.039270116
#> Cell_67 -0.07923754  0.94106456  0.70841553 -0.58852970 -1.187486788
#> Cell_68 -0.78227143 -0.36704850 -1.14472485  0.37896992 -0.575503194
#> Cell_69 -0.49757525  0.69171445 -0.08529687  1.15775658 -1.200060853
#> Cell_70  0.07448409 -0.47803765  0.31886193 -0.40801761  0.102958114
#> Cell_71 -1.96761252 -0.22572806  0.79197339 -0.99518277 -0.241954880
#> Cell_72  0.36370229 -1.13379989 -1.26915839  0.35377269 -0.595387192
#> Cell_73 -1.04983377  0.53507351  0.59515563  0.55775983 -0.018241704
#> Cell_74 -1.07691447 -0.59380752  0.32043951 -0.26926591  0.148293280
#> Cell_75 -0.31876708  1.15689246  0.74019110 -0.17951689  0.290012534
#> Cell_76  0.35877301 -0.34808276  0.26530565 -0.52871143 -1.093677461
#> Cell_77 -0.28628683  1.56157762  0.87808570  0.77148495 -0.613832958
#> Cell_78  0.87061975  1.08415112  1.23303911  1.01786848  1.225228502
#> Cell_79 -1.39881182 -0.65195737 -1.30642433  0.64224987  0.220100928
#> Cell_80 -0.64978182 -0.39811707 -0.08598820 -0.22970563 -0.177035635
#> 
#> $centroids
#>                PC1          PC2       PC3       PC4        PC5         PC6
#> TypeA -0.013650308 -0.038936381 0.1167637 0.4894623 0.04299117 0.008281261
#> TypeB -0.004123877  0.005444855 0.3684942 0.1318596 0.29600973 0.006484870
#>                PC7        PC8         PC9      PC10        PC11       PC12
#> TypeA  0.005934824 -0.2152040 -0.06835705 0.2216317 -0.06659547 0.28039331
#> TypeB -0.047133673 -0.3411115  0.01898248 0.2969842 -0.14963850 0.07059479
#>             PC13       PC14       PC15
#> TypeA 0.29055749 -0.0946136  0.1043664
#> TypeB 0.01645159  0.2249051 -0.1664359
#> 
#> $cell_milestone_geodesic
#>              TypeA       TypeB
#> Cell_01 0.00000000 0.768620055
#> Cell_02 0.70606537 0.062554687
#> Cell_03 0.35573428 0.412885775
#> Cell_04 0.00000000 0.768620055
#> Cell_05 0.00000000 0.768620055
#> Cell_06 0.27037444 0.498245616
#> Cell_07 0.00000000 0.768620055
#> Cell_08 0.00000000 0.768620055
#> Cell_09 0.76862005 0.000000000
#> Cell_10 0.76862005 0.000000000
#> Cell_11 0.00000000 0.768620055
#> Cell_12 0.76862005 0.000000000
#> Cell_13 0.00000000 0.768620055
#> Cell_14 0.62369762 0.144922438
#> Cell_16 0.76862005 0.000000000
#> Cell_17 0.00000000 0.768620055
#> Cell_18 0.00000000 0.768620055
#> Cell_19 0.00000000 0.768620055
#> Cell_20 0.76862005 0.000000000
#> Cell_21 0.00000000 0.768620055
#> Cell_22 0.76862005 0.000000000
#> Cell_23 0.00000000 0.768620055
#> Cell_24 0.00000000 0.768620055
#> Cell_25 0.76862005 0.000000000
#> Cell_26 0.00000000 0.768620055
#> Cell_28 0.00000000 0.768620055
#> Cell_29 0.00000000 0.768620055
#> Cell_30 0.01176437 0.756855685
#> Cell_31 0.00000000 0.768620055
#> Cell_32 0.00000000 0.768620055
#> Cell_33 0.68275167 0.085868389
#> Cell_34 0.00000000 0.768620055
#> Cell_35 0.00000000 0.768620055
#> Cell_36 0.75904343 0.009576623
#> Cell_37 0.20350019 0.565119866
#> Cell_38 0.00000000 0.768620055
#> Cell_39 0.00000000 0.768620055
#> Cell_40 0.76862005 0.000000000
#> Cell_41 0.00000000 0.768620055
#> Cell_42 0.76862005 0.000000000
#> Cell_43 0.76862005 0.000000000
#> Cell_44 0.47691528 0.291704772
#> Cell_45 0.76862005 0.000000000
#> Cell_46 0.76862005 0.000000000
#> Cell_47 0.76862005 0.000000000
#> Cell_48 0.24313777 0.525482288
#> Cell_49 0.76862005 0.000000000
#> Cell_50 0.00000000 0.768620055
#> Cell_51 0.76862005 0.000000000
#> Cell_52 0.76862005 0.000000000
#> Cell_53 0.34125717 0.427362882
#> Cell_54 0.76862005 0.000000000
#> Cell_55 0.76862005 0.000000000
#> Cell_56 0.76862005 0.000000000
#> Cell_57 0.76862005 0.000000000
#> Cell_58 0.02428539 0.744334660
#> Cell_59 0.00000000 0.768620055
#> Cell_60 0.76862005 0.000000000
#> Cell_61 0.76862005 0.000000000
#> Cell_62 0.76862005 0.000000000
#> Cell_63 0.74038291 0.028237146
#> Cell_64 0.69971838 0.068901679
#> Cell_65 0.00000000 0.768620055
#> Cell_66 0.76862005 0.000000000
#> Cell_67 0.76862005 0.000000000
#> Cell_68 0.55215162 0.216468431
#> Cell_69 0.76862005 0.000000000
#> Cell_70 0.00000000 0.768620055
#> Cell_71 0.41660214 0.352017913
#> Cell_72 0.76862005 0.000000000
#> Cell_73 0.73505816 0.033561892
#> Cell_74 0.24503316 0.523586894
#> Cell_75 0.76862005 0.000000000
#> Cell_76 0.05276446 0.715855594
#> Cell_77 0.41553102 0.353089038
#> Cell_78 0.00000000 0.768620055
#> Cell_79 0.76862005 0.000000000
#> Cell_80 0.26443040 0.504189657
#> 
#> $branch_assignment
#>  [1] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#>  [6] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [11] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [16] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [21] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [26] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [31] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [36] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [41] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [46] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [51] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [56] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [61] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [66] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [71] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> [76] "TypeA->TypeB" "TypeA->TypeB" "TypeA->TypeB"
#> 
#> attr(,"class")
#> [1] "dynwrap::with_trajectory" "dynwrap::data_wrapper"   
#> [3] "list"                    
#> 
#> $trajectory$preprocessing_summary
#> $trajectory$preprocessing_summary$n_hvgs
#> [1] 60
#> 
#> $trajectory$preprocessing_summary$n_strict_hvgs
#> [1] 0
#> 
#> $trajectory$preprocessing_summary$n_ref_cells
#> [1] 80
#> 
#> $trajectory$preprocessing_summary$n_sim_cells
#> [1] 78
#> 
#> $trajectory$preprocessing_summary$hvg_method
#> [1] "scran_modelGeneVar"
#> 
#> 
#> 
#> $cross_modality
#> $cross_modality$correlation
#> $metrics
#> $metrics$batch_asw
#> [1] 0.9747094
#> 
#> $metrics$cell_type_asw
#> [1] 0.5582026
#> 
#> $metrics$graph_connectivity
#> [1] 1
#> 
#> $metrics$ari
#> [1] 0.6765538
#> 
#> $metrics$nmi
#> [1] 0.5821401
#> 
#> 
#> $integration_summary_table
#>                                  Category               Property
#> 1 (VI) Multiomics Integration & Alignment multiomics_integration
#> 2 (VI) Multiomics Integration & Alignment multiomics_integration
#> 3 (VI) Multiomics Integration & Alignment multiomics_integration
#> 4 (VI) Multiomics Integration & Alignment multiomics_integration
#> 5 (VI) Multiomics Integration & Alignment multiomics_integration
#>               Metric     Value            Modality
#> 1          batch_asw 0.9747094 Joint (Cross-Modal)
#> 2      cell_type_asw 0.5582026 Joint (Cross-Modal)
#> 3 graph_connectivity 1.0000000 Joint (Cross-Modal)
#> 4                ari 0.6765538 Joint (Cross-Modal)
#> 5                nmi 0.5821401 Joint (Cross-Modal)
#> 
#> $embedding
#>                 PC1        PC2          PC3         PC4         PC5         PC1
#> Cell_01 -1.17846364 -3.3620143  2.832131471 -0.49628064 -1.21988018 -0.78014965
#> Cell_02  1.42838960 -1.2228432 -0.989926755 -0.69393958 -3.39190107 -2.57283173
#> Cell_03 -1.42490520 -1.0468425  0.202934552 -1.21119303 -2.31916863 -0.60515437
#> Cell_04 -1.06933980 -2.7504402 -1.228661162  0.05345014 -2.09346790 -2.76780368
#> Cell_05  0.98358953 -2.2988903 -0.037091539 -2.67446836  0.72489341 -1.21484367
#> Cell_06 -0.32524365 -1.6600872 -0.351170597  0.17951809  0.11368101 -1.67404688
#> Cell_07  1.58213201 -1.1048090  0.720067029  1.52175549  1.06914614 -3.03141842
#> Cell_08  0.73070844  0.2403091 -0.973331994 -0.70108525  2.15213063 -4.33788275
#> Cell_09  0.70502400 -2.3210341 -2.184280026  0.17430096 -0.46025961 -1.80168481
#> Cell_10 -0.59232754 -1.1825113 -0.615439494 -0.98171086 -0.87372370 -0.67878696
#> Cell_11  0.58712932 -2.9950572 -1.152764703  0.22484763  1.32808848 -2.71167285
#> Cell_12  0.43126754 -1.1020843  0.151172037  0.26639312 -1.75269790 -0.54419822
#> Cell_13 -1.42113109 -1.5759778 -1.602256798 -0.57246799  2.36613284 -0.74331447
#> Cell_14 -0.05192821  0.1663499 -0.500022248  0.22140916 -0.72672626 -1.56948214
#> Cell_15  0.11691228 -1.3867138 -0.161213253 -0.45435700  0.26365585 -2.46461787
#> Cell_16 -0.78571376 -0.8376623 -1.772675099 -2.37024973  1.36700316 -2.57750024
#> Cell_17 -1.69517230  0.6603628  1.428849595 -1.38929978 -0.99904504 -0.98543065
#> Cell_18  0.40064894 -1.9457076 -1.506406725 -1.61636128  0.50071219 -0.73806625
#> Cell_19 -0.41980888  0.4538118 -1.349502609 -0.19922586 -1.58040528 -0.33303823
#> Cell_20  1.46434486 -1.2343204  3.625451572  1.74450263  0.60380046 -1.85861445
#> Cell_21  1.76104326 -1.4770661 -3.322837531 -2.06714684 -0.31450330 -2.42402398
#> Cell_22  0.63928316 -2.0352659  2.776254623 -1.86796824 -2.96652688 -0.96670014
#> Cell_23  1.82753342 -1.6076166 -0.002012967  1.09813762  0.48419065 -0.41494719
#> Cell_24  0.21364763 -1.1728183 -1.159937359 -0.23879431 -0.33261871 -0.96333301
#> Cell_25 -0.16555814 -0.8151257  0.602545147 -0.29354402 -0.19021145 -0.25734473
#> Cell_26 -0.53083948 -2.6068450 -0.623662682  2.02999571 -2.96372915  0.43131756
#> Cell_27  0.36484841  0.1900586 -0.180105477 -2.28183196  0.16538689  0.03820620
#> Cell_28  0.08425883  0.2777805 -0.633757543 -0.95290488  1.98340989 -1.89532966
#> Cell_29 -2.05559864  0.9670976 -0.435364354 -2.79498216  0.65675849 -0.81525798
#> Cell_30 -0.91529542  0.5404075 -0.698822727 -1.80112008 -0.13398394 -3.24310098
#> Cell_31 -1.19222993 -3.1898022  2.102352135 -0.69015468  0.98116089 -1.35241288
#> Cell_32  1.44325898  0.6880820 -1.233151877 -0.97699905 -1.31830683 -2.51773455
#> Cell_33 -0.93381381  1.0919875 -1.259411154  1.04168088 -2.12151800 -1.31670438
#> Cell_34  1.52776477  0.1609286  2.280730344 -1.44829845 -0.29505126 -2.05074992
#> Cell_35  1.05594101 -0.2607200 -0.367692481  0.26765782 -1.27080005 -0.43566205
#> Cell_36  0.83513551  0.1710158 -0.605634983  0.48034248  2.44228630 -2.44029502
#> Cell_37 -0.67996363  1.2806676  1.460735617  3.02102688  1.60645816 -0.54579848
#> Cell_38  0.52722865 -2.2127014 -1.423009012 -0.51930066 -1.21970962 -2.56594900
#> Cell_39 -0.04636575  1.5506284  1.282945487 -0.89850288 -0.55987114 -0.76356831
#> Cell_40  0.87222301  0.7021129 -0.139275593  1.38522542 -0.83672673  0.24414215
#> Cell_41  0.77314823 -0.3842292  2.647102629  0.60616368 -0.23056131 -0.70516407
#> Cell_42 -3.13880075  0.6721441  1.880416367 -0.74058236 -1.32398594  1.95422632
#> Cell_43  2.77788628  2.6329101  0.659576916  0.10232690 -2.69710922  4.30851814
#> Cell_44  5.26375951  4.3232000  1.405115461 -0.86216699  2.17953512  2.52275785
#> Cell_45 -0.85015764  0.8426889  2.270578152  1.87636293  0.48424770  0.89576750
#> Cell_46 -1.86192312  0.5600133  4.684151714 -2.33008158 -1.53822145  0.82175137
#> Cell_47 -2.07413549  2.1078234  0.199719238 -0.17616694  0.42043908  2.52691283
#> Cell_48 -3.95839996  2.2011691  1.609539134 -1.01805856  0.82672594  2.88847038
#> Cell_49 -0.48860323  5.8756782 -0.283282530  2.57406197 -0.96842795 -0.85277685
#> Cell_50  1.43474033 -0.8368804  0.513037754 -0.16589704  0.08727049  4.18577692
#> Cell_51 -0.82716733 -0.7769833  0.021868755  0.19387276 -0.46067533 -0.65446694
#> Cell_52  0.82836039  0.2334939 -1.082878307 -0.53488708  0.55097451  2.19051905
#> Cell_53 -0.70510429  2.4316846 -1.718947646  0.51915339  1.13462839  3.09380811
#> Cell_54  0.74003196  3.3347574  3.378151865 -1.51498769  1.24798349  2.56057546
#> Cell_55 -3.28541452 -0.1242079 -0.381109882  1.35695894  2.32123568 -0.54968957
#> Cell_56 -0.40085090  0.7764992  0.071089893  0.78223475 -0.75962643  2.48048702
#> Cell_57  4.37883099  2.0442258 -1.066291015 -4.78467392  1.86550836 -4.03756054
#> Cell_58 -0.52108522  0.1190477 -0.620382467  0.68407896  1.72752732  1.44914106
#> Cell_59 -0.49140667 -1.0737989 -1.151688671  1.95356152  1.40204104 -0.24968072
#> Cell_60 -1.45703689  0.8395974  0.059775165 -0.48696564  2.62600778  1.52694968
#> Cell_61  5.00865983  2.6022940  0.627327655  1.40323500 -1.07943896  0.64999926
#> Cell_62 -3.88920095  2.1167509  2.428906340 -1.62319543 -2.33770331  1.36915880
#> Cell_63 -0.68317593  1.3906733 -0.505104626  1.70922885 -2.68037270  1.99828146
#> Cell_64 -4.07386716  0.9591172 -1.305755151  1.91871141  4.46517914  1.20869779
#> Cell_65  1.16216539  0.4232923 -1.678342020 -2.61219965  1.98078305  4.05273420
#> Cell_66  1.99252469  1.2130451  0.026579570 -0.02747427  0.66688601  3.30894482
#> Cell_67  2.50501102 -1.9950947 -1.779449457  3.25911035  0.69469875  1.63749492
#> Cell_68  2.90678374 -4.3700472  5.654472532  3.72765703  3.72327716 -0.03503347
#> Cell_69 -0.05623038  1.2134511 -2.403465870 -0.31678245 -0.44250132  2.77216664
#> Cell_70 -1.34483513 -0.1563568  0.196560575  0.72983676  0.20374742  0.58763771
#> Cell_71 -1.72216327  0.1381487 -1.882667889  1.77403418  0.10204593 -0.14104038
#> Cell_72  0.66606488  0.3806816 -4.240226320  3.39525461 -2.14135289  1.38825787
#> Cell_73 -0.22480183  2.4832289 -0.341905231  1.72109551 -1.85591999  1.43065026
#> Cell_74 -2.31540849  0.4219395  0.204966366  2.69084157 -0.24581173 -0.06908903
#> Cell_75 -2.30726646  1.5282493 -1.415171706  2.10024444  1.27243270  2.39790483
#> Cell_76  0.10832176 -0.4362481 -0.099493586  0.08859031 -0.31664358  1.14876094
#> Cell_77  1.91432571  1.3323137 -0.196262428  0.27496846  0.46421740  2.10020071
#> Cell_78  2.09375687  0.0623231  1.402852565  2.43246117 -1.57448230  2.04955703
#> Cell_79 -2.33673376 -1.9826386  0.087219202 -2.33371795  0.76635274  2.29481553
#> Cell_80  0.36078353  1.1393994 -0.833333914 -1.86426436  0.57102642  1.73936172
#>                 PC2         PC3          PC4          PC5
#> Cell_01  2.98807650 -1.31382697 -0.014214359  1.527894315
#> Cell_02 -2.10899603 -0.60552198  0.445326676  0.348037409
#> Cell_03  1.86537414 -3.40887143 -1.020815830  0.024592707
#> Cell_04  1.79997795 -1.81519616  1.420911235 -1.558321682
#> Cell_05 -1.39826635 -0.07499597  1.376136341  4.370142012
#> Cell_06  0.45516478  0.49221753  2.657124632  3.132163568
#> Cell_07  1.20115539 -0.09694159 -3.216476025  1.076417169
#> Cell_08 -1.58259514 -0.37651344  1.616679344 -0.155276767
#> Cell_09  4.90423488 -1.10604214 -0.720742651 -1.189834475
#> Cell_10  2.35922485  3.20527021  0.689635533  1.278711186
#> Cell_11 -0.07666269  1.44945632  0.241518505  0.076573716
#> Cell_12  0.95822459  0.88124969  2.574232164  0.978772144
#> Cell_13  0.52725966  1.08483641  0.094569888 -1.909893843
#> Cell_14  0.78621418 -1.02711870  0.258095048 -0.808102857
#> Cell_15 -0.72597527  2.02706971  0.008155624  0.002692197
#> Cell_16 -0.76220234  1.85834153  0.506867978 -0.267753926
#> Cell_17  0.56316934 -0.50505486  1.146588130 -1.965886755
#> Cell_18  0.22226714  0.85461994  2.653126100  0.272102765
#> Cell_19  0.78434458 -0.22374246  1.741423998 -0.942366551
#> Cell_20 -3.78920646 -1.35798121  1.253710239 -1.680674087
#> Cell_21  0.34435679 -0.74752870 -0.266288861  2.333752304
#> Cell_22 -1.44576788  1.00471665  1.395392019 -0.642638371
#> Cell_23 -0.32540152 -0.65009861  3.692231662 -0.013512087
#> Cell_24  0.19869079  0.82904731 -0.378837508  0.849725696
#> Cell_25  0.70868079 -1.19266228  1.092976704 -0.993740709
#> Cell_26  2.01190132 -0.21827838  1.077926648 -0.078424715
#> Cell_27  1.69339121 -1.13904102 -1.459175596 -3.144007694
#> Cell_28  1.96835874  1.37878618  1.228777839 -1.228238857
#> Cell_29  1.34109686 -1.07397810 -1.963715892 -0.275497076
#> Cell_30 -0.14390660 -1.35289087  1.712099369 -0.477402663
#> Cell_31 -0.27856396  0.81215440 -1.573450088 -1.323572873
#> Cell_32 -1.85206113  0.58551329  0.294639735  2.947949003
#> Cell_33  2.43209120  0.54387053 -3.069232504  0.486529096
#> Cell_34  0.88869737 -0.57799743  0.027403539  1.762139009
#> Cell_35  0.94132320  1.30443963 -2.551774347  1.506886676
#> Cell_36  0.90931674 -0.31223888 -1.263956436 -0.527060492
#> Cell_37 -0.10118316 -1.11546094 -2.538749301  1.841532608
#> Cell_38 -0.31738035 -2.46471019 -2.287574298  1.814872370
#> Cell_39  3.03667715 -0.06558983  1.059401750 -2.089820731
#> Cell_40 -0.85613238 -4.19360311  2.980851913 -0.038041283
#> Cell_41 -1.59281157  1.12240269  2.454708212 -2.749815384
#> Cell_42  1.48308544  2.52726658 -0.142540028 -1.910479371
#> Cell_43  1.15070047  0.47804470  0.195013447 -0.585691006
#> Cell_44 -2.17293229  0.11384776  1.284630146  3.106018546
#> Cell_45 -4.03157741  0.47251713  0.250861105 -0.325231011
#> Cell_46 -0.34926908 -0.62940204  0.148240063  0.777902396
#> Cell_47 -1.58395496 -0.90356222 -0.056930026 -0.969483664
#> Cell_48  2.19109931  1.22424038  1.646394964  1.832379017
#> Cell_49 -3.37817506 -2.30902327 -2.587436384 -0.390925367
#> Cell_50 -2.47288281  1.07020363 -3.110544223  0.925264861
#> Cell_51 -0.39552248  0.12034344 -1.585827081  1.657766624
#> Cell_52  2.22004799 -1.41212655  1.784034416 -1.911269104
#> Cell_53  0.40104647  1.92493980 -0.159399235 -0.136229691
#> Cell_54 -0.75636799 -1.19285460  0.349358471 -0.800966953
#> Cell_55 -2.59538945  2.96765583 -2.075568300  1.004566528
#> Cell_56  0.28676454  2.93594753 -0.162082357 -1.259273275
#> Cell_57 -4.09959773  1.38039308 -1.150281859 -3.318334653
#> Cell_58 -2.55802742 -0.39803066 -0.515220106  0.919139858
#> Cell_59  1.99106978  1.84994286 -2.303668504  0.852786275
#> Cell_60 -2.30175195  0.14798511  0.008742470 -1.729781251
#> Cell_61  2.25510477  0.69201959 -1.917540711  0.897184154
#> Cell_62 -0.04676941 -1.17390887  0.670107930  0.344029843
#> Cell_63  0.36954789 -2.83006278  0.472466043  2.732774078
#> Cell_64  0.92947195  3.33816806  1.654751985  2.769209357
#> Cell_65  1.94936994 -2.25494787  0.426692003  0.984175189
#> Cell_66 -0.91940726 -3.37440973 -0.109896042 -2.058512411
#> Cell_67  1.26601950 -3.52323586 -1.392358266  0.025322848
#> Cell_68 -0.87607563  0.86761169  0.253781122 -1.060650512
#> Cell_69 -1.20477722 -1.46589142  2.792562802 -0.343601585
#> Cell_70 -0.94333263 -2.21613811 -3.110601647 -0.455059362
#> Cell_71 -0.99372892 -1.09033241 -2.120912246 -2.291648131
#> Cell_72  1.30952225  1.32440126 -2.013802797 -0.413959369
#> Cell_73 -2.41974283 -2.00794332 -0.828039136  0.150334636
#> Cell_74  0.77474481  1.76080712  1.844460001 -0.156436138
#> Cell_75  0.61199855 -0.74101576 -0.338327237  1.529014855
#> Cell_76 -0.45505621  2.59022937 -1.557341178 -1.838106902
#> Cell_77 -2.75699938  1.75873655 -0.037612347  0.214193530
#> Cell_78 -1.60938835  1.16367154  0.377513580  2.608206954
#> Cell_79  0.39655272  0.52825626 -0.161998885 -1.622183330
#> Cell_80  0.80242276  3.86754949 -0.097189080 -2.324048534
#> 
#> attr(,"class")
#> [1] "scSimEval_integration" "list"                 
#> 
#> $cross_modality$foscttm
#> NULL
#> 
#> $cross_modality$generation
#> NULL
#> 
#> $cross_modality$coupling
#> NULL
#> 
#> $cross_modality$coaccessibility
#> NULL
#> 
#> $cross_modality$coexpression
#> NULL
#> 
#> $cross_modality$profile
#> NULL
#> 
#> 
#> $resource_usage
#> $resource_usage$mode
#> [1] "combined"
#> 
#> $resource_usage$core_seconds
#> NULL
#> 
#> $resource_usage$memory_per_core_mib
#> NULL
#> 
#> $resource_usage$memory_per_core_gib
#> NULL
#> 
#> $resource_usage$memory_mb
#> NULL
#> 
#> $resource_usage$memory_input
#> NULL
#> 
#> $resource_usage$memory_unit
#> [1] "MiB"
#> 
#> $resource_usage$elapsed_time
#> NULL
#> 
#> $resource_usage$elapsed_time_input
#> NULL
#> 
#> $resource_usage$time_unit
#> [1] "seconds"
#> 
#> $resource_usage$sim_cores
#> NULL
#> 
#> $resource_usage$rna
#> $resource_usage$rna$core_seconds
#> NULL
#> 
#> $resource_usage$rna$memory_per_core_mib
#> NULL
#> 
#> $resource_usage$rna$elapsed_time
#> NULL
#> 
#> $resource_usage$rna$elapsed_time_input
#> NULL
#> 
#> $resource_usage$rna$time_unit
#> [1] "seconds"
#> 
#> $resource_usage$rna$memory_mb
#> NULL
#> 
#> $resource_usage$rna$memory_input
#> NULL
#> 
#> $resource_usage$rna$memory_unit
#> [1] "MiB"
#> 
#> $resource_usage$rna$sim_cores
#> NULL
#> 
#> 
#> $resource_usage$atac
#> $resource_usage$atac$core_seconds
#> NULL
#> 
#> $resource_usage$atac$memory_per_core_mib
#> NULL
#> 
#> $resource_usage$atac$elapsed_time
#> NULL
#> 
#> $resource_usage$atac$elapsed_time_input
#> NULL
#> 
#> $resource_usage$atac$time_unit
#> [1] "seconds"
#> 
#> $resource_usage$atac$memory_mb
#> NULL
#> 
#> $resource_usage$atac$memory_input
#> NULL
#> 
#> $resource_usage$atac$memory_unit
#> [1] "MiB"
#> 
#> $resource_usage$atac$sim_cores
#> NULL
#> 
#> 
#> 
```
