# Demonstration 3: Multiomics (scRNA-seq + scATAC-seq) Simulation Evaluation

## Introduction

Single-cell multiomics technologies allow simultaneous measurement of
multiple cellular features, such as gene expression (scRNA-seq) and
chromatin accessibility (scATAC-seq). Simulating multiomics data
presents unique challenges because algorithms must preserve not only the
distinct statistical properties of each individual modality, but also
the biological coupling and regulatory connections between them.

------------------------------------------------------------------------

## Paired vs. Unpaired Multiomics: How `scSimEval` Works

A common question is whether `scSimEval` supports both **paired** and
**unpaired** single-cell multiomics data.

**The answer is yes — `scSimEval` provides comprehensive evaluation for
both data types:**

#### 1. Paired Multiomics (e.g., 10x Chromium Multiome, SHARE-seq, SNARE-seq)

- **Definition:** Both scRNA-seq expression counts and scATAC-seq
  chromatin accessibility peaks are measured simultaneously from the
  **exact same individual cells** (each column in the RNA matrix
  corresponds to the same cell barcode in the ATAC matrix).
- **Evaluation in `scSimEval`:**
  - Runs all unimodal evaluations across both modalities (Categories 1
    through 6 and Category 8).
  - **Unlocks full Category 7 (Cross-Modal Coupling):**
    - **FOSCTTM (Fraction of Samples Closer Than The True Match):**
      Evaluates whether the paired cell is its nearest neighbor in joint
      embedding space.
    - **<Match@1> Rate:** Percentage of cells whose true multiomics
      profile is rank-1 in joint space.
    - **Cross-Modal Generation Fidelity:** Predicts scATAC from scRNA
      (or vice versa) and computes cell-wise Pearson correlation, RMSE,
      and cosine similarity.
    - **Direct Peak-to-Gene Linkage:** Correlates promoter/enhancer peak
      accessibility with target gene transcription across matched single
      cells.

#### 2. Unpaired Multiomics (Independent scRNA-seq and scATAC-seq)

- **Definition:** scRNA-seq and scATAC-seq are profiled from separate
  single cells derived from the same biological tissue or cell line.
  Cell barcodes do not match one-to-one.
- **Evaluation in `scSimEval`:**
  - **Unimodal Fidelity:** Evaluates distributions, zero patterns,
    clustering, batch mixing, marker genes, and trajectories for RNA and
    ATAC independently.
  - **Cell-Type / Cluster-Level Cross-Modal Transfer
    (`evaluate_cross_modal_prediction`):** Tests whether a model trained
    on ATAC accessibility can accurately predict cell types defined in
    RNA.
  - **Regulatory Modularity (`calc_network_jaccard`,
    `calc_coexpression_module_fidelity`):** Compares co-accessibility
    and co-expression network topologies without requiring cell-by-cell
    pairing.
  - **Multi-Dataset Benchmarking (`evaluate_multiple_datasets`):**
    Consolidates and ranks simulation accuracy across unpaired
    modalities side-by-side.

------------------------------------------------------------------------

## 1. Loading the Packaged Multiomics Dataset

`scSimEval` includes paired scRNA-seq and scATAC-seq matrices
(`example_multiomics`):

``` r
library(scSimEval)

# Load packaged multiomics example dataset
data(example_multiomics)

cat("Reference RNA dimensions:", dim(example_multiomics$ref_multi$rna), "\n")
#> Reference RNA dimensions: 60 80
cat("Reference ATAC dimensions:", dim(example_multiomics$ref_multi$atac), "\n")
#> Reference ATAC dimensions: 60 80
cat("Simulated RNA dimensions:", dim(example_multiomics$sim_multi$rna), "\n")
#> Simulated RNA dimensions: 60 80
cat("Simulated ATAC dimensions:", dim(example_multiomics$sim_multi$atac), "\n")
#> Simulated ATAC dimensions: 60 80
cat("Cell types:", levels(example_multiomics$cell_types), "\n")
#> Cell types: TypeA TypeB
```

The dataset contains: \* `ref_multi`: Named list containing `rna` (count
matrix) and `atac` (peak accessibility matrix). \* `sim_multi`: Named
list containing simulated `rna` and `atac` matrices. \* `cell_types`:
Cell type annotations across the cells. \* `batch_info`: Batch labels
across the cells. \* `resource_stats`: Hardware logs (runtime in seconds
and peak RAM in MiB).

------------------------------------------------------------------------

## 2. Category 7: Cross-Modality Coupling Analysis

For multiomics data, `scSimEval` provides specialized functions to
evaluate inter-modality coupling.

#### 2.1 Alignment in Joint Space (FOSCTTM & <Match@1>)

Fraction of Samples Closer Than The True Match (FOSCTTM) quantifies how
well the two modalities are aligned in joint latent space:

``` r
# Compute reduced representations (e.g. PCA) for RNA and ATAC
pca_rna  <- prcomp(t(example_multiomics$ref_multi$rna), rank. = 5)$x
pca_atac <- prcomp(t(example_multiomics$ref_multi$atac), rank. = 5)$x

foscttm_res <- calc_foscttm(x = pca_rna, y = pca_atac)

cat("Bidirectional FOSCTTM score:", round(foscttm_res$foscttm, 4), "\n")
#> Bidirectional FOSCTTM score: 0.455
cat("Top-1 Match Rate (Match@1):", round(foscttm_res$match_at_1, 4), "\n")
#> Top-1 Match Rate (Match@1): 0.0125
```

A score of $`0.00`$ indicates perfect alignment (where each cell’s true
match is its closest neighbor), while $`0.50`$ represents random chance.

#### 2.2 Cross-Modal Cell Type Label Transfer

We test whether chromatin accessibility profiles in ATAC can accurately
predict cell types annotated in RNA:

``` r
transfer_res <- evaluate_cross_modal_prediction(
  mod1_data  = example_multiomics$sim_multi$atac,
  cell_types = example_multiomics$cell_types
)

cat("Cross-Modal Prediction Accuracy:", round(transfer_res$cross_modal_accuracy, 4), "\n")
#> Cross-Modal Prediction Accuracy: 1
cat("Cross-Modal Macro F1 Score:", round(transfer_res$cross_modal_F1, 4), "\n")
#> Cross-Modal Macro F1 Score: 1
```

------------------------------------------------------------------------

## 3. Master Multiomics Pipeline: `evaluate_multiomics_accuracy`

The master function
[`evaluate_multiomics_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiomics_accuracy.md)
coordinates evaluation across all 8 categories:

``` r
master_eval <- evaluate_multiomics_accuracy(
  ref_multi  = example_multiomics$ref_multi,
  sim_multi  = example_multiomics$sim_multi,
  cell_types = example_multiomics$cell_types,
  batch_info = example_multiomics$batch_info,
  cpu_time   = example_multiomics$resource_stats$cpu_time,
  memory_mb  = example_multiomics$resource_stats$memory_mb,
  verbose    = FALSE
)

# Inspect tidy master summary table
head(master_eval$benchmark_summary_table)
#>                    Category     Property        Metric        Value Modality
#> 1 Distributional Properties library_size           MAD 1.500000e+01      rna
#> 2 Distributional Properties library_size            KS 2.125000e-01      rna
#> 3 Distributional Properties library_size           MAE 1.506250e+01      rna
#> 4 Distributional Properties library_size          RMSE 1.571504e+01      rna
#> 5 Distributional Properties library_size            OV 8.269394e-01      rna
#> 6 Distributional Properties library_size Bhattacharyya 1.188839e-04      rna
```

The output contains: \* `benchmark_summary_table`: Master table
containing metric names, categories, raw values, and standardized
scores. \* `category_scores`: Mean normalized score for each of the 8
categories. \* `composite_score`: Overall unified simulation accuracy
index.

------------------------------------------------------------------------

## 4. Multi-Method / Multi-Dataset Consolidated Benchmarking

To benchmark multiple simulators (or multiple datasets) simultaneously,
use
[`evaluate_multiple_datasets()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_multiple_datasets.md):

``` r
# Create benchmark collection comparing multiple simulated datasets
dataset_collection <- list(
  "Method 1-RNA"  = list(ref = example_multiomics$ref_multi$rna,  sim = example_multiomics$sim_multi$rna),
  "Method 1-ATAC" = list(ref = example_multiomics$ref_multi$atac, sim = example_multiomics$sim_multi$atac),
  "Method 2-RNA"  = list(ref = example_multiomics$ref_multi$rna,  sim = example_multiomics$sim_multi$rna),
  "Method 2-ATAC" = list(ref = example_multiomics$ref_multi$atac, sim = example_multiomics$sim_multi$atac)
)

consolidated <- evaluate_multiple_datasets(
  datasets          = dataset_collection,
  pair_by_prefix    = TRUE,
  compute_bivariate = FALSE,
  verbose           = FALSE
)

# View dataset overview and ranking
print(consolidated$dataset_overview)
#>    Dataset      Data_Name            Modality Total_Metrics Mean_KS_Distance
#> 1 Method 1   Method 1-RNA           scRNA-seq           153           0.1616
#> 2 Method 1  Method 1-ATAC          scATAC-seq           155           0.1488
#> 3 Method 1 Method 1-Joint Joint (Cross-Modal)            10           0.3167
#> 4 Method 2   Method 2-RNA           scRNA-seq           153           0.1616
#> 5 Method 2  Method 2-ATAC          scATAC-seq           155           0.1488
#> 6 Method 2 Method 2-Joint Joint (Cross-Modal)            10           0.3167
#>   Mean_Wasserstein Mean_RMSE
#> 1           2.3207    2.4504
#> 2           1.2090    1.4124
#> 3           0.0635    0.0663
#> 4           2.3207    2.4504
#> 5           1.2090    1.4124
#> 6           0.0635    0.0663
```

------------------------------------------------------------------------

## 5. Visualizing Benchmark Results

`scSimEval` provides built-in functions to generate high-resolution
figures:

``` r
# 1. Master Benchmark Bubble Matrix across the 8 categories
plot_benchmark_bubble_matrix(
  data  = consolidated,
  title = "Benchmarking Single-Cell Multiomics Simulators"
)

# 2. Evaluation Summary Horizontal Bar Matrix
plot_evaluation_summary(
  data  = consolidated,
  title = "Multiomics Simulator Leaderboard"
)
```

------------------------------------------------------------------------

## 6. Summary

Whether your study employs paired single-cell multiomics (e.g. 10x
Multiome) or independent unpaired datasets, `scSimEval` provides a
rigorous, automated, and ground-truth-free framework to evaluate
simulation fidelity across all biological and technical dimensions.
