# Evaluate scATAC-seq Differential Accessibility (DA) Fidelity (Sub-category V-B)

Evaluates the preservation of Differentially Accessible Regions (DARs)
between reference and simulated scATAC-seq datasets using biological
replicate pseudobulks and edgeR-QLF (Squair et al., Nature
Communications 2021; Teo et al., Nature Communications 2024 / Libra)
across standard benchmarking metrics (Duo et al., Genome Biology 2024 /
Simpipe).

## Usage

``` r
evaluate_atac_da_fidelity(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  batch_ref = NULL,
  batch_sim = NULL,
  binarize = TRUE,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  min_cells = 3,
  min_reps = 2
)
```

## Arguments

- ref_data:

  Reference scATAC-seq count matrix (peaks x cells).

- sim_data:

  Simulated scATAC-seq count matrix (peaks x cells).

- ref_celltypes:

  Factor or character vector of reference cell types.

- sim_celltypes:

  Factor or character vector of simulated cell types.

- batch_ref:

  Factor or character vector of biological replicate / batch labels for
  reference cells.

- batch_sim:

  Factor or character vector of biological replicate / batch labels for
  simulated cells.

- binarize:

  Logical. If TRUE (default for scATAC-seq), binarizes single-cell
  counts prior to pseudobulk aggregation.

- fdr_cutoff:

  Adjusted p-value significance threshold (default 0.05).

- logfc_cutoff:

  Absolute log2 fold-change cutoff (default 0).

- min_cells:

  Minimum number of cells per cell type (default 3).

- min_reps:

  Minimum number of replicates per group (default 2).

## Value

A list containing:

- da_summary_table:

  Tidy summary data.frame with 6 scATAC-seq DA metrics: `da_score`,
  `pval_dist_score`, `true_proportion`, `accuracy`, `precision`,
  `recall`, `f1_score`, and `auc`

- pairwise_results:

  Detailed contrast-level data.frame with per-contrast metrics and
  confusion matrix counts

- pseudobulks_ref:

  List of reference pseudobulk matrices per cell type

- pseudobulks_sim:

  List of simulated pseudobulk matrices per cell type

- da_results_ref:

  List of edgeR-QLF DA results for reference

- da_results_sim:

  List of edgeR-QLF DA results for simulation

## Details

Category V is divided into two distinct sub-categories:

- **Sub-category V-A (scRNA-seq Differential Expression)**: Computes 8
  metrics: `de_score`, `pval_dist_score`, `true_proportion`, `accuracy`,
  `precision`, `recall`, `f1_score`, and `auc`.

- **Sub-category V-B (scATAC-seq Differential Accessibility)**: Computes
  8 metrics: `da_score`, `pval_dist_score`, `true_proportion`,
  `accuracy`, `precision`, `recall`, `f1_score`, and `auc`.

If `batch_ref` or `batch_sim` is NULL, the analysis is gracefully
skipped, adhering to the strict rule that single-cell ATAC-seq
biological signal evaluation requires genuine biological replicates to
prevent pseudoreplication bias.

## Examples

``` r
data(example_scatac, package = "scSimEval")
rep_id <- factor(rep(c("Rep1", "Rep2"), length.out = ncol(example_scatac$ref)))
res <- evaluate_atac_da_fidelity(
  ref_data = example_scatac$ref,
  sim_data = example_scatac$sim,
  ref_celltypes = example_scatac$cell_types,
  sim_celltypes = example_scatac$cell_types,
  batch_ref = rep_id,
  batch_sim = rep_id
)
#> calcNormFactors has been renamed to normLibSizes
#> calcNormFactors has been renamed to normLibSizes
res$da_summary_table
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
```
