# Fréchet Single-Cell Distance (FSD)

Computes the single-cell analogue of Fréchet Inception Distance (FID) on
low-dimensional PCA embeddings between reference and simulated cell
populations.

## Usage

``` r
calc_frechet_singlecell_distance(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  n_pcs = 15
)
```

## Arguments

- ref_mat:

  Matrix of reference cells (cells x features or features x cells).

- sim_mat:

  Matrix of simulated cells (cells x features or features x cells).

- cells_as_cols:

  Logical, whether cells are columns (default TRUE).

- n_pcs:

  Number of principal components to evaluate (default 15).

## Value

A list with Fréchet distance (FSD), mean discrepancy, and covariance
trace discrepancy.
