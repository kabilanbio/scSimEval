# Calculate Median Absolute Deviation (MAD)

Calculate Median Absolute Deviation (MAD)

## Usage

``` r
calc_mad(ref, sim, align = TRUE)
```

## Arguments

- ref:

  Numeric vector of reference distribution.

- sim:

  Numeric vector of simulated distribution.

- align:

  Logical, whether to sort and quantile-align vectors. Default TRUE.

## Value

Numeric MAD value.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_mad(ref, sim)
#> [1] 0.3038233
```
