# Evaluate Epigenomic Supervised Cell-Type Annotation (EpiAnno / SCAN-ATAC-Sim)

Trains a supervised classification model on reference single-cell
chromatin accessibility profiles (scATAC-seq / scCAS) and evaluates
cell-type prediction / projection accuracy on simulated cells against
ground-truth cell-type labels (Chen et al., Nat Mach Intell 2022).

## Usage

``` r
evaluate_epigenomic_annotation(
  ref_data,
  sim_data,
  ref_celltypes,
  sim_celltypes,
  method = c("knn", "centroid"),
  n_pcs = 20,
  k = 5
)
```

## Arguments

- ref_data:

  Reference count/accessibility matrix (features x cells).

- sim_data:

  Simulated count/accessibility matrix (features x cells).

- ref_celltypes:

  Factor or character vector of cell types for reference cells.

- sim_celltypes:

  Factor or character vector of true cell types for simulated cells.

- method:

  Classification method: "knn" (k-nearest neighbors on PCA, default) or
  "centroid" (nearest centroid classifier).

- n_pcs:

  Number of principal components for dimensionality reduction (default
  20).

- k:

  Number of nearest neighbors for k-NN (default 5).

## Value

A list containing overall accuracy, balanced accuracy, macro F1, macro
precision, macro recall, Cohen's kappa, and per-class performance
metrics.
