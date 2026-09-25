# Cell-Specific Mixing Score (CMS)

Tests the hypothesis that group-specific distance distributions of
k-nearest neighbor cells have the same underlying unspecified
distribution. Computes cell-specific p-values via Anderson-Darling tests
(if kSamples installed) or two-sample Kolmogorov-Smirnov /
Kruskal-Wallis tests, with optional delegation to CellMixS::cms.

## Usage

``` r
calc_cms(
  data = NULL,
  batch_info,
  k = NULL,
  n_pcs = 30,
  coords = NULL,
  cell_min = 4,
  unbalanced = FALSE
)
```

## Arguments

- data:

  Count or log-normalized expression matrix (features x cells) or NULL
  if coords provided.

- batch_info:

  Factor or character vector of batch assignments for each cell.

- k:

  Neighborhood size k for k-NN graph. Default is
  min(table(batch_info))/2 or 30.

- n_pcs:

  Number of principal components to compute if coords is NULL. Default
  is 30.

- coords:

  Optional precomputed embedding matrix (cells x dimensions, e.g. PCA or
  UMAP).

- cell_min:

  Minimum number of cells per batch in neighborhood to perform test.
  Default is 4.

- unbalanced:

  Logical, whether to set pure neighborhoods to NA. Default is FALSE.

## Value

A named list:

- cms_scores:

  Numeric vector of CMS p-values for each cell (higher means better
  mixed)

- mean_cms:

  Mean CMS score across cells

- median_cms:

  Median CMS score across cells
