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
