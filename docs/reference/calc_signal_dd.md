# Detect Differential Distribution (DD) Genes (Kolmogorov-Smirnov Test)

Detect Differential Distribution (DD) Genes (Kolmogorov-Smirnov Test)

## Usage

``` r
calc_signal_dd(exprs_mat, cell_types)
```

## Arguments

- exprs_mat:

  Log-normalized matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

## Value

Vector of adjusted p-values.

## Examples

``` r
data(example_scrna, package = "scSimEval")
calc_signal_dd(example_scrna$ref, example_scrna$cell_types)
#>  [1] 0.9895096 0.5642625 0.5642625 0.8388210 0.4794990 0.9895096 0.6921413
#>  [8] 0.9895096 0.8388210 0.6921413 0.5642625 0.9895096 0.4794990 0.2070355
#> [15] 0.5642625 0.9895096 0.9895096 0.9895096 0.9895096 0.9176670 0.9176670
#> [22] 0.2096603 0.8388210 0.8388210 0.9895096 0.9895096 0.9895096 0.9895096
#> [29] 0.9895096 0.4794990 0.9895096 0.9895096 0.9895096 0.9176670 0.9895096
#> [36] 0.9895096 0.9895096 0.9895096 0.9895096 0.9895096 0.9895096 0.9176670
#> [43] 0.9895096 0.9895096 0.8388210 0.9895096 0.9895096 0.9895096 0.9895096
#> [50] 0.9895096 0.8388210 0.6921413 0.5642625 0.9895096 0.9176670 0.9895096
#> [57] 0.9895096 0.2095607 0.5642625 0.9895096
```
