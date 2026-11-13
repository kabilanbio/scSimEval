# Calculate Fasano-Franceschini 2D Kolmogorov-Smirnov Test Statistic

Integrated from simpipe.

## Usage

``` r
calc_fasano_franceschini(ref_mat, sim_mat, threads = 1)
```

## Arguments

- ref_mat:

  2-column numeric matrix for reference.

- sim_mat:

  2-column numeric matrix for simulation.

- threads:

  CPU threads. Default 1.

## Value

Estimated 2D KS statistic.

## Examples

``` r
ref <- matrix(stats::rnorm(100), 50, 2)
sim <- matrix(stats::rnorm(100), 50, 2)
calc_fasano_franceschini(ref, sim)
#> [1] 1100
```
