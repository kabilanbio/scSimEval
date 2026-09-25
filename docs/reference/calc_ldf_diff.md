# Local Density Differences (ldfDiff)

Quantifies cell-specific changes in the Local Density Factor (LDF)
before and after data integration or between reference and simulated
datasets.

## Usage

``` r
calc_ldf_diff(coords_pre, coords_post, k = 15, h = 1, c = 1)
```

## Arguments

- coords_pre:

  Matrix of coordinates before integration / reference (cells x
  dimensions).

- coords_post:

  Matrix of coordinates after integration / simulated (cells x
  dimensions).

- k:

  Number of nearest neighbors. Default is 15.

- h:

  Bandwidth parameter for Gaussian kernel. Default is 1.

- c:

  Scaling constant. Default is 1.

## Value

A named list:

- diff:

  Numeric vector of absolute LDF differences per cell

- mean_ldf_diff:

  Mean absolute LDF difference (lower indicates better structure
  preservation)

- median_ldf_diff:

  Median absolute LDF difference

- ldf_pre:

  LDF in pre-integration / reference space

- ldf_post:

  LDF in post-integration / simulated space
