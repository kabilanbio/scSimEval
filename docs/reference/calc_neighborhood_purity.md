# Calculate Neighborhood Purity

Integrated from simPIC (Chugh et al., 2024; bluster::neighborPurity).
For each cell, evaluates the proportion of its k-nearest neighbors that
share the same cluster or cell-type identity.

## Usage

``` r
calc_neighborhood_purity(data, cluster_labels, k = NULL, is_distance = FALSE)
```

## Arguments

- data:

  Matrix of coordinates (cells x dimensions) or expression matrix
  (features x cells) or dist matrix.

- cluster_labels:

  Vector of cluster or cell-type labels.

- k:

  Number of nearest neighbors (default 5% of cells).

- is_distance:

  Logical, whether data is already a distance matrix. Default FALSE.

## Value

A numeric vector of neighborhood purity scores (values from 0 to 1),
with an attribute "mean_purity" containing the overall average.
