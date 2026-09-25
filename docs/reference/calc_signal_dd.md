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
