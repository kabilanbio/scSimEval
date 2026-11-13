# Subsample Large Vectors for Accelerated Distributional Distance Estimation

Subsample Large Vectors for Accelerated Distributional Distance
Estimation

## Usage

``` r
subsample_dist_vector(x, max_n = 10000)
```

## Arguments

- x:

  Numeric vector.

- max_n:

  Maximum sample size (default 10000, matching SimBench framework).

## Value

Subsampled numeric vector if length(x) \> max_n, else x.
