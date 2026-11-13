# Calculate Bhattacharyya Distance (BH)

Calculate Bhattacharyya Distance (BH)

## Usage

``` r
calc_bhattacharyya(ref, sim, align = TRUE)
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- align:

  Logical, whether to align lengths. Default TRUE.

## Value

Statistical divergence between discrete probability measures.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_bhattacharyya(ref, sim)
#> [1] 0.004823546
```
