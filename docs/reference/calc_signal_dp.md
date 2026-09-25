# Detect Differential Proportion / Zero-Inflation (DP) Genes (Chisq Test)

Detect Differential Proportion / Zero-Inflation (DP) Genes (Chisq Test)

## Usage

``` r
calc_signal_dp(exprs_mat, cell_types, threshold = 0)
```

## Arguments

- exprs_mat:

  Count or expression matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

- threshold:

  Value threshold for zero/detection status (default 0).

## Value

Vector of adjusted p-values.
