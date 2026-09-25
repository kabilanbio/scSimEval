# Generative Precision and Recall for Single-Cell Manifolds

Evaluates the realism (Precision) and coverage/diversity (Recall) of
generated single-cell distributions relative to empirical reference data
using k-NN hyperspheres.

## Usage

``` r
calc_generative_precision_recall(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  k = 5,
  max_cells = 500
)
```

## Arguments

- ref_mat:

  Matrix of reference cells (cells x features or features x cells).

- sim_mat:

  Matrix of simulated cells (cells x features or features x cells).

- cells_as_cols:

  Logical, whether cells are columns (default TRUE).

- k:

  Number of nearest neighbors to define the local manifold radius
  (default 5).

- max_cells:

  Maximum number of cells to subsample (default 500).

## Value

A list containing Generative Precision, Generative Recall, and
Generative F1.
