# Evaluate scRNA-seq Differential Expression (DE) Fidelity (Sub-category V-A)

Evaluates the preservation of Differentially Expressed Genes (DEGs)
between reference and simulated scRNA-seq datasets using biological
replicate pseudobulks and edgeR-QLF (Squair et al., Nature
Communications 2021 / Libra) across standard benchmarking metrics (Duo
et al., Genome Biology 2024 / Simpipe).

## Usage

``` r
evaluate_deg_fidelity(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  batch_ref = NULL,
  batch_sim = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  classifier = c("svm", "knn", "rf"),
  min_cells = 3,
  min_reps = 2
)

evaluate_rna_de_fidelity(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  batch_ref = NULL,
  batch_sim = NULL,
  fdr_cutoff = 0.05,
  logfc_cutoff = 0,
  classifier = c("svm", "knn", "rf"),
  min_cells = 3,
  min_reps = 2
)
```

## Arguments

- ref_data:

  Reference scRNA-seq count matrix (genes x cells).

- sim_data:

  Simulated scRNA-seq count matrix (genes x cells).

- ref_celltypes:

  Factor or character vector of reference cell types.

- sim_celltypes:

  Factor or character vector of simulated cell types.

- batch_ref:

  Factor or character vector of biological replicate / batch labels for
  reference cells. Required.

- batch_sim:

  Factor or character vector of biological replicate / batch labels for
  simulated cells. Required.

- fdr_cutoff:

  Adjusted p-value significance threshold for DEGs (default 0.05).

- logfc_cutoff:

  Absolute log2 fold-change cutoff (default 0).

- classifier:

  Classifier algorithm for marker gene predictability: "svm" (default),
  "knn", or "rf".

- min_cells:

  Minimum number of cells per cell type (default 3).

- min_reps:

  Minimum number of replicates per group (default 2).

## Value

A list containing:

- deg_summary_table:

  Tidy summary data.frame with 6 scRNA-seq DE metrics: `de_score`,
  `pval_dist_score`, `true_proportion`, `accuracy`, `precision`,
  `recall`, `f1_score`, and `auc`

- pairwise_results:

  Detailed contrast-level data.frame with per-contrast metrics and
  confusion matrix counts

- pseudobulks_ref:

  List of reference pseudobulk matrices per cell type

- pseudobulks_sim:

  List of simulated pseudobulk matrices per cell type

- de_results_ref:

  List of edgeR-QLF DE results for reference

- de_results_sim:

  List of edgeR-QLF DE results for simulation

- ml_classification:

  Detailed classification results using simulated pseudobulk DEGs

## Details

Category V is divided into two distinct sub-categories:

- **Sub-category V-A (scRNA-seq Differential Expression)**: Computes 8
  metrics: `de_score`, `pval_dist_score`, `true_proportion`, `accuracy`,
  `precision`, `recall`, `f1_score`, and `auc`.

- **Sub-category V-B (scATAC-seq Differential Accessibility)**: Computes
  8 metrics: `da_score`, `pval_dist_score`, `true_proportion`,
  `accuracy`, `precision`, `recall`, `f1_score`, and `auc`.

If `batch_ref` or `batch_sim` is NULL, the analysis is gracefully
skipped, adhering to the strict rule that single-cell biological signal
evaluation requires genuine biological replicates to prevent false
positive inflation and pseudoreplication bias.
