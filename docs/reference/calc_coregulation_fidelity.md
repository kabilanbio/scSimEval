# Multi-Omics Co-Regulation & Modularity Fidelity

Evaluates the preservation of co-regulated gene/feature modules between
reference and simulated multi-omics datasets (Monzo et al., 2025;
Arzalluz-Luque et al., 2022).

## Usage

``` r
calc_coregulation_fidelity(
  ref_data,
  sim_data,
  modules,
  method = c("pearson", "spearman")
)
```

## Arguments

- ref_data:

  Reference feature-by-cell matrix or data frame.

- sim_data:

  Simulated feature-by-cell matrix or data frame.

- modules:

  A list of character vectors representing feature clusters/modules, or
  a named vector/factor of module assignments per feature.

- method:

  Correlation method: "pearson" (default) or "spearman".

## Value

A list containing module correlation r, RMSE, MAE, and modularity
fidelity.

## Examples

``` r
data(example_scrna, package = "scSimEval")
modules <- list(Module_1 = rownames(example_scrna$ref)[seq_len(30)],
                Module_2 = rownames(example_scrna$ref)[31:60])
calc_coregulation_fidelity(example_scrna$ref, example_scrna$sim, modules)
#> $module_correlation_r
#> [1] -0.003588998
#> 
#> $module_correlation_rmse
#> [1] 0.1635834
#> 
#> $module_correlation_mae
#> [1] 0.1316224
#> 
#> $ref_modularity_ratio
#> [1] -0.05300986
#> 
#> $sim_modularity_ratio
#> [1] -0.08805737
#> 
#> $modularity_fidelity
#> [1] 0
#> 
```
