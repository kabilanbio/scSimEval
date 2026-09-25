# Local Structure Preservation Metric

Calculates the proportion of overlapping k-nearest neighbors between an
original/reference embedding and an integrated/simulated embedding
(adapted from Seurat LocalStruct and CellMixS locStructure).

## Usage

``` r
calc_local_structure_metric(coords_pre, coords_post, k = 30)
```

## Arguments

- coords_pre:

  Matrix of coordinates before integration / reference (cells x
  dimensions).

- coords_post:

  Matrix of coordinates after integration / simulated (cells x
  dimensions).

- k:

  Number of nearest neighbors. Default is 30.

## Value

A named list:

- cell_overlaps:

  Numeric vector of overlap fractions per cell

- mean_local_structure:

  Mean overlap fraction across cells (higher indicates better
  preservation)

- median_local_structure:

  Median overlap fraction
