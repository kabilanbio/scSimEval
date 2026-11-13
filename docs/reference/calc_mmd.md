# Calculate Maximum Mean Discrepancy (MMD) on Single-Cell PCA Manifolds

Integrated from ACTIVA (Heydari et al., Bioinformatics 2022;
SindiLab/ACTIVA) and scGAN (Marouf et al., Nature Communications 2020).
Evaluates how closely the multivariate distribution of simulated
single-cell profiles matches the reference biological manifold in a
shared Principal Component Analysis (PCA) latent space using a
multi-scale Radial Basis Function (RBF) Gaussian kernel in a Reproducing
Kernel Hilbert Space (RKHS).

## Usage

``` r
calc_mmd(
  ref_mat,
  sim_mat,
  cells_as_cols = TRUE,
  n_pcs = 50,
  scales = c(0.01, 0.1, 1, 10, 100),
  max_cells = 1000,
  max_features = 2000
)
```

## Arguments

- ref_mat:

  Matrix of reference cells (features x cells by default, or
  SingleCellExperiment).

- sim_mat:

  Matrix of simulated cells (features x cells by default, or
  SingleCellExperiment).

- cells_as_cols:

  Logical, whether cells are columns (default TRUE).

- n_pcs:

  Number of principal components for shared latent projection (default
  50, as in ACTIVA).

- scales:

  Numeric vector of bandwidth multipliers for the multi-scale RBF kernel
  (default `c(0.01, 0.1, 1, 10, 100)`).

- max_cells:

  Maximum number of cells per dataset to subsample for speed (default
  1000).

- max_features:

  Maximum number of highly variable features to retain prior to PCA
  (default 2000).

## Value

Numeric scalar representing the multi-scale RBF Maximum Mean Discrepancy
(`MMD >= 0`; lower values indicate higher manifold fidelity).

## Examples

``` r
ref <- matrix(stats::rnorm(400), 20, 20)
sim <- matrix(stats::rnorm(400), 20, 20)
calc_mmd(ref, sim)
#> [1] 0.2215625
```
