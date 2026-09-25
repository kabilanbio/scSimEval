# Detect Differential Variability (DV) Genes (Bartlett Test)

Detect Differential Variability (DV) Genes (Bartlett Test)

## Usage

``` r
calc_signal_dv(exprs_mat, cell_types)
```

## Arguments

- exprs_mat:

  Log-normalized matrix (genes x cells).

- cell_types:

  Factor of 2 cell types.

## Value

Vector of adjusted p-values.
