# Evaluate Accessibility-Sparsity Curve Concordance (simATAC)

Compares the non-linear relationship between peak mean accessibility and
non-zero proportion (NZP) between reference and simulated scATAC-seq
datasets.

## Usage

``` r
evaluate_accessibility_sparsity_curve(ref_data, sim_data, poly_degree = 2)
```

## Arguments

- ref_data:

  Reference count matrix (peaks x cells).

- sim_data:

  Simulated count matrix (peaks x cells).

- poly_degree:

  Degree of polynomial (default: 2).

## Value

A list containing `curve_rmse` and `r2_discrepancy` (along with fitted
reference and simulated curve parameters).

## Examples

``` r
data(example_scatac, package = "scSimEval")
evaluate_accessibility_sparsity_curve(example_scatac$ref, example_scatac$sim)
#> $curve_rmse
#> [1] 0.01233884
#> 
#> $r2_discrepancy
#> [1] 0.0008535249
#> 
#> $ref_r_squared
#> [1] 0.9692997
#> 
#> $sim_r_squared
#> [1] 0.9684462
#> 
```
