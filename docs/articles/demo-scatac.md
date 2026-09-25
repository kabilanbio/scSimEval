# Demonstration 2: scATAC-seq Simulation Evaluation

## Introduction

Single-cell assay for transposase-accessible chromatin sequencing
(scATAC-seq) profiles epigenomic landscapes at single-cell resolution.
Unlike scRNA-seq expression counts, scATAC-seq data have unique
statistical characteristics: \* **High Sparsity:** Over 90–98% of values
in peak-by-cell matrices are zeroes due to the diploid nature of
eukaryotic genomes (each locus having at most 2 copies per cell). \*
**Near-Binary Nature:** Values primarily represent open (accessible)
vs. closed chromatin states. \* **Peak Co-Accessibility:** Genomic loci
in spatial proximity or under shared transcription factor control
exhibit correlated accessibility.

This tutorial demonstrates how to evaluate synthetic scATAC-seq data
using **`scSimEval`**.

------------------------------------------------------------------------

## 1. Loading the Packaged scATAC-seq Dataset

`scSimEval` includes an empirical reference and simulated scATAC-seq
peak accessibility dataset (`example_scatac`):

``` r
library(scSimEval)

# Load packaged scATAC-seq example data
data(example_scatac)

cat("Reference scATAC-seq dimensions:", dim(example_scatac$ref), "\n")
#> Reference scATAC-seq dimensions: 60 80
cat("Simulated scATAC-seq dimensions:", dim(example_scatac$sim), "\n")
#> Simulated scATAC-seq dimensions: 60 80
cat("Cell types:", levels(example_scatac$cell_types), "\n")
#> Cell types: TypeA TypeB
cat("Batches:", levels(example_scatac$batch_info), "\n")
#> Batches: Batch1 Batch2
```

The matrices contain integer insertion counts across genomic peaks
($`60\text{ peaks} \times 80\text{ cells}`$).

------------------------------------------------------------------------

## 2. Epigenomic Distribution & Summary Properties

We evaluate whether simulated peak insertion frequencies and cell
library sizes match empirical distributions using
[`evaluate_simulation_accuracy()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_simulation_accuracy.md):

``` r
atac_res <- evaluate_simulation_accuracy(
  ref_data          = example_scatac$ref,
  sim_data          = example_scatac$sim,
  compute_bivariate = FALSE,
  verbose           = FALSE
)

# Inspect peak accessibility and library size statistical distances
head(atac_res$metrics_summary_table, 8)
#>                    Category     Property        Metric        Value
#> 1 Distributional Properties library_size           MAD 2.0000000000
#> 2 Distributional Properties library_size            KS 0.1750000000
#> 3 Distributional Properties library_size           MAE 2.3375000000
#> 4 Distributional Properties library_size          RMSE 2.4874685928
#> 5 Distributional Properties library_size            OV 0.8346458491
#> 6 Distributional Properties library_size Bhattacharyya 0.0001055641
#> 7 Distributional Properties library_size   Wasserstein 2.3375000000
#> 8 Distributional Properties library_size ECDF_DiffArea 0.0863425926
```

------------------------------------------------------------------------

## 3. Epigenomic Cluster & Cell-Type Concordance

In scATAC-seq, identifying cell types requires clustering cells based on
chromatin accessibility profiles. We benchmark whether the simulator
preserves genuine epigenomic separation using
[`evaluate_clustering_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_clustering_metrics.md):

``` r
atac_clust <- evaluate_clustering_metrics(
  data       = example_scatac$sim,
  cell_types = example_scatac$cell_types,
  ref_data   = example_scatac$ref
)

cat("Average Silhouette Width (Simulated):", round(atac_clust$silhouette, 4), "\n")
#> Average Silhouette Width (Simulated): 0.0465
cat("Davies-Bouldin Index:", round(atac_clust$davies_bouldin, 4), "\n")
#> Davies-Bouldin Index: NA
cat("Adjusted Rand Index (ARI):", round(atac_clust$ARI, 4), "\n")
#> Adjusted Rand Index (ARI): 0.7627
cat("Normalized Mutual Information (NMI):", round(atac_clust$NMI, 4), "\n")
#> Normalized Mutual Information (NMI): 0.721
```

------------------------------------------------------------------------

## 4. Epigenomic Batch Integration & Confounder Handling

When simulating multi-sample or multi-donor epigenomic experiments,
batch variations can confound biological signals. We evaluate batch
mixing metrics using
[`evaluate_batch_metrics()`](https://kabilanbio.github.io/scSimEval/reference/evaluate_batch_metrics.md):

``` r
atac_batch <- evaluate_batch_metrics(
  data       = example_scatac$sim,
  batch_info = example_scatac$batch_info,
  cell_types = example_scatac$cell_types,
  verbose    = FALSE
)

cat("Batch Shannon Entropy:", round(atac_batch$shannon_entropy, 4), "\n")
#> Batch Shannon Entropy: 0.9671
cat("Principal Component Regression (PCR R2):", round(atac_batch$pcr_r2, 4), "\n")
#> Principal Component Regression (PCR R2): 0.0125
cat("Cross-Batch Transfer Accuracy:", round(atac_batch$cross_batch_accuracy, 4), "\n")
#> Cross-Batch Transfer Accuracy: 0.825
```

------------------------------------------------------------------------

## 5. Chromatin Peak Co-Accessibility & Regulatory Coupling

In real biological cells, peaks within the same topologically
associating domain (TAD) or co-regulated by the same transcription
factor complexes exhibit correlated accessibility:

``` r
# Evaluate peak co-accessibility fidelity
coacc_res <- calc_peak_coaccessibility_fidelity(
  ref_atac = example_scatac$ref,
  sim_atac = example_scatac$sim
)

cat("Co-Accessibility RV Coefficient:", round(coacc_res$rv_coefficient, 4), "\n")
#> Co-Accessibility RV Coefficient: 0.6014
cat("Peak Co-Accessibility Pearson Correlation:", round(coacc_res$coaccessibility_pearson, 4), "\n")
#> Peak Co-Accessibility Pearson Correlation: 0.1128
cat("Peak Co-Accessibility Spearman Correlation:", round(coacc_res$coaccessibility_spearman, 4), "\n")
#> Peak Co-Accessibility Spearman Correlation: 0.0902
cat("Frobenius Matrix Distance:", round(coacc_res$frobenius_distance, 4), "\n")
#> Frobenius Matrix Distance: 0.0053
```

------------------------------------------------------------------------

## 6. Master Unimodal scATAC-seq Pipeline

Execute the comprehensive evaluation workflow in a single coordinated
call:

``` r
master_atac <- evaluate_simulation_accuracy(
  ref_data          = example_scatac$ref,
  sim_data          = example_scatac$sim,
  cpu_time          = 18.2,
  memory_mb         = 512.4,
  compute_bivariate = FALSE,
  verbose           = FALSE
)

# Inspect tidy summary table
head(master_atac$metrics_summary_table)
#>                    Category     Property        Metric        Value
#> 1 Distributional Properties library_size           MAD 2.0000000000
#> 2 Distributional Properties library_size            KS 0.1750000000
#> 3 Distributional Properties library_size           MAE 2.3375000000
#> 4 Distributional Properties library_size          RMSE 2.4874685928
#> 5 Distributional Properties library_size            OV 0.8346458491
#> 6 Distributional Properties library_size Bhattacharyya 0.0001055641
```

------------------------------------------------------------------------

## 7. Conclusion

`scSimEval` effectively benchmarks scATAC-seq simulations, validating
peak insertion distributions, extreme sparsity, cell type clustering,
batch mixing, and chromatin peak co-accessibility without requiring
synthetic ground-truth annotations.
