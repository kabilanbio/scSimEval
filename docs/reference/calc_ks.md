# Calculate Kolmogorov-Smirnov Distance (1D KS)

Calculate Kolmogorov-Smirnov Distance (1D KS)

## Usage

``` r
calc_ks(ref, sim)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

## Value

Maximum vertical distance between empirical CDFs (between 0 and 1).

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_ks(ref, sim)
#> [1] 0.0646207
```
