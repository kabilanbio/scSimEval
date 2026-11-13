# Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint Distributions

Compute All Bivariate (2D) Distance & Accuracy Metrics for Joint
Distributions

## Usage

``` r
calc_all_bivariate_metrics(ref_mat, sim_mat, metric_prefix = "", threads = 1)
```

## Arguments

- ref_mat:

  2-column matrix for reference.

- sim_mat:

  2-column matrix for simulation.

- metric_prefix:

  Optional prefix string for metric names.

- threads:

  CPU threads.

## Value

Named list of bivariate metrics.

## Examples

``` r
ref <- matrix(stats::rnorm(100), 50, 2)
sim <- matrix(stats::rnorm(100), 50, 2)
calc_all_bivariate_metrics(ref, sim)
#> $Fasano_Franceschini_2D_KS
#> [1] 850
#> 
#> $EMD_2D
#> [1] 0.03219795
#> 
#> $NN_Mismatch_2D
#> [1] 0.04
#> 
#> $Between_Dataset_Silh_2D_Global
#> [1] -0.01104689
#> 
```
