# Calculate Root Mean Squared Error (RMSE)

Calculate Root Mean Squared Error (RMSE)

## Usage

``` r
calc_rmse(ref, sim, align = TRUE)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- align:

  Logical, whether to sort and quantile-align vectors. Default TRUE.

## Value

Quadratic error penalizing large discrepancies.
