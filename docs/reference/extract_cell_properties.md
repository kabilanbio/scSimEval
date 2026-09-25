# Extract Comprehensive Cell-Level Properties

Computes library size, log-library size, zero fraction / detection
frequency, pairwise cell correlation (on top HVGs), cell k-NN hubness,
pairwise PCA distance, and outlier proportions.

## Usage

``` r
extract_cell_properties(
  data,
  max_cells_cor = 1000,
  n_top_hvgs = 500,
  verbose = FALSE
)
```

## Arguments

- data:

  Count matrix (features x cells).

- max_cells_cor:

  Maximum number of cells for computing pairwise correlation (default
  1000).

- n_top_hvgs:

  Number of highly variable genes to use for correlation and PCA
  (default 500).

- verbose:

  Logical, whether to print progress messages.

## Value

A named list of cell-level summary vectors.
