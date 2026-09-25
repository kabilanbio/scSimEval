# Inverse Simpson Index (ISI) for Batch Mixing

Evaluates the diversity of batch labels within each cell's k-nearest
neighborhood using the Inverse Simpson Index (1 / sum(p_b^2)).
Implements a pure-R, lightweight formulation of LISI with optional
distance weighting.

## Usage

``` r
calc_isi(coords, batch_info, k = 30, weighted = TRUE)
```

## Arguments

- coords:

  Matrix of cell coordinates (cells x dimensions).

- batch_info:

  Factor or character vector of batch assignments.

- k:

  Number of nearest neighbors. Default is 30.

- weighted:

  Logical, whether to weight neighbor contributions by inverse distance.
  Default is TRUE.

## Value

A named list:

- isi_scores:

  Numeric vector of cell-level ISI values

- mean_isi:

  Mean ISI score (ranges from 1 to number of batches; higher indicates
  better mixing)

- median_isi:

  Median ISI score
