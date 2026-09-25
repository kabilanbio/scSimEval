# Align Two Empirical Distributions

Sorts both distributions and aligns their lengths using either quantile
interpolation or nearest-neighbor resampling, ensuring element-wise
distance metrics (MAD, MAE, RMSE) are mathematically well-defined even
when the number of cells or features differs.

## Usage

``` r
align_distributions(
  ref,
  sim,
  method = c("quantile", "resample"),
  n_points = NULL
)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- method:

  Alignment method: "quantile" (default, interpolated quantiles) or
  "resample" (resampling).

- n_points:

  Number of points to evaluate when method = "quantile". Default is
  max(length(ref), length(sim)).

## Value

A list with aligned numeric vectors: `list(ref = ..., sim = ...)`.
