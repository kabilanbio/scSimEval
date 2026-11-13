# Compute All Univariate Distance & Accuracy Metrics for a Feature

Compute All Univariate Distance & Accuracy Metrics for a Feature

## Usage

``` r
calc_all_univariate_metrics(ref, sim, metric_prefix = "")
```

## Arguments

- ref:

  Numeric vector of reference values.

- sim:

  Numeric vector of simulated values.

- metric_prefix:

  Optional prefix string for metric names.

## Value

Named list of univariate distance metrics.

## Examples

``` r
ref <- stats::rnorm(50)
sim <- stats::rnorm(50)
calc_all_univariate_metrics(ref, sim)
#> $MAD
#> [1] 0.2552427
#> 
#> $KS
#> [1] 0.2205651
#> 
#> $MAE
#> [1] 0.3120819
#> 
#> $RMSE
#> [1] 0.4446398
#> 
#> $Bhattacharyya
#> [1] 0.006985793
#> 
#> $Wasserstein
#> [1] 0.3120819
#> 
#> $Runs_Statistic
#> [1] -0.4020356
#> 
#> $Runs_PValue
#> [1] 0.3438289
#> 
#> $NN_Mismatch
#> [1] 0.05
#> 
#> $Between_Dataset_Silh_Global
#> [1] 0.01116562
#> 
#> $KDE_ISE
#> [1] 0.02391995
#> 
```
