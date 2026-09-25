# Maximum Mean Discrepancy (MMD) with Gaussian RBF Kernel

Evaluates the distributional discrepancy between empirical reference and
generated/simulated single-cell high-dimensional profiles or latent
embeddings.

## Usage

``` r
calc_mmd(ref_mat, sim_mat, cells_as_cols = TRUE, sigma = NULL, max_cells = 500)
```

## Arguments

- ref_mat:

  Matrix of reference cells (cells x features or features x cells).

- sim_mat:

  Matrix of simulated cells (cells x features or features x cells).

- cells_as_cols:

  Logical, whether cells are columns. If TRUE (default for scRNA-seq),
  the matrix is transposed so rows represent cells.

- sigma:

  Optional kernel bandwidth. If NULL, uses the median pairwise distance
  heuristic.

- max_cells:

  Maximum number of cells to subsample for speed (default 500).

## Value

A list with MMD, squared MMD, and the kernel bandwidth sigma used.
