# Calculate Mean Absolute Error (MAE)

Calculate Mean Absolute Error (MAE)

## Usage

``` r
calc_mae(ref, sim, align = TRUE)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- align:

  Logical, whether to sort and quantile-align vectors. Default TRUE.

## Value

Mean absolute difference.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_mae(ref, sim)
#> [1] 0.2081537
```
