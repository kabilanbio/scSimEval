# Seurat Mixing Metric

Calculates the mixing metric originally proposed in Seurat (Stuart et
al., Cell 2019) and adapted by CellMixS. For each cell, finds the rank
of the k_pos-th neighbor from each batch in the sorted neighborhood, and
computes the median rank across batches.

## Usage

``` r
calc_seurat_mixing_metric(coords, batch_info, k = 300, k_pos = 5)
```

## Arguments

- coords:

  Matrix of cell coordinates (cells x dimensions).

- batch_info:

  Factor or character vector of batch assignments.

- k:

  Maximum neighborhood size to search. Default is 300.

- k_pos:

  The rank position to extract per batch. Default is 5.

## Value

A named list:

- mixing_metrics:

  Numeric vector of median ranks per cell

- mean_mixing_metric:

  Mean mixing rank (lower indicates better mixing)

- median_mixing_metric:

  Median mixing rank
