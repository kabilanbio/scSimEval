# Fast Truncated PCA on Top Highly Variable Features

Fast Truncated PCA on Top Highly Variable Features

## Usage

``` r
fast_pca(mat, n_pcs = 20, n_hvg = 2000, scale = FALSE, normalize = FALSE)
```

## Arguments

- mat:

  Matrix or dgCMatrix (cells x features).

- n_pcs:

  Number of principal components. Default 20.

- n_hvg:

  Number of top highly variable genes to retain. Default 2000.

- scale:

  Logical, whether to scale features. Default FALSE.

- normalize:

  Logical, whether to apply library-size log-normalization if raw
  counts. Default FALSE.

## Value

Matrix of cell principal component embeddings (cells x n_pcs).
